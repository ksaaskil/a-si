%clear
path='/wrk/kisaaski/lammps_sisu/a-si/';
file_forces=strcat(path,folder,'.Fij.dat');
%path='/proj/quantum-data/Kimmo/lammps/'
%file_forces=strcat(path,'060214d_Si.Fij.dat');
%file_forces=strcat(path,'060514c_fcc.Fij.dat');
%file_forces=strcat(path,'140514a5_cnt.Fij.dat');
%file_forces=strcat(path,'090614a_Si.Fij.dat');
%file_forces=strcat(path,'010814a2_cnt.Fij.dat');
%file_forces=strcat(path,'270115a.Fij.dat');

%file_forces=strcat(path,'090614a2_cnt.Fij_type7.dat');
%file_forces=strcat(path,'180214a2_acnt.Fij.dat');
%file_forces=strcat(path,'040314b_SiGe.Fij_ids7.dat');

k_tol=1e-5;
file_forces
fid1=fopen(file_forces,'r');

% File format
% NL $NL
% NR $NR
% "atom_id atom_type" of interface atoms, NL+NR lines
% HSTEP $hstep

ss=textscan(fid1,'%s%d',1);
NL=ss{2};
ss=textscan(fid1,'%s%d',1);
NR=ss{2};

ss=textscan(fid1,'%d%d',NL+NR);
atom_ids=ss{1};
atom_types=ss{2};

ids_L=find(atom_types==6);
ids_R=find(atom_types==7);

ss=textscan(fid1,'%s%f',1);
HSTEP=ss{2};

%return

if 1 % Read forces from file
    
    Kii=zeros(3*NL,3*NL);
    Kij=zeros(3*NL,3*NR);

    for i=1:3*NL % Loop over the particles
        if (NL>1e2 && mod(i,1e1)==0)
           fprintf('i=%d/%d\n',i,3*NL); 
        end
        % ITEM: TIMESTEP
        %ss=textscan(fid1,'%s',9);
        %ss=textscan(fid1,'%s',1);
        %ss
        % ITEM: ENTRIES index c_atomids[1] c_atomids[2] c_forces
        ss=textscan(fid1,'%d%f%f%f',NL+NR,'headerlines',10);
        
        % Positive direction shift for the given composite index i
        % (particle number and component)
        inds=ss{1};
        fxs1=ss{2};
        fys1=ss{3};
        fzs1=ss{4};
        
        ss=textscan(fid1,'%d%f%f%f',NL+NR,'headerlines',10);      
        
        % Negative direction shift for the given composite index i
        
        inds=ss{1};
        fxs2=ss{2};
        fys2=ss{3};
        fzs2=ss{4};
        
        % Forces on ids_L and ids_7 due to atom displacements of type ids_L
        
        Fii=[fxs1(ids_L)-fxs2(ids_L),fys1(ids_L)-fys2(ids_L),fzs1(ids_L)-fzs2(ids_L)];
        Fij=[fxs1(ids_R)-fxs2(ids_R),fys1(ids_R)-fys2(ids_R),fzs1(ids_R)-fzs2(ids_R)];
        
        Fii=Fii';
        Fij=Fij';
        Fii=Fii(:);
        Fij=Fij(:);

        Kij(i,:)=Fij;
        Kii(i,:)=Fii;
        %forces(i,:)=ss1{4};

    end
end

% The spring constant matrices Kij=d^2V/du_idu_j=-dF_j/du_i come with minus

Kii=-Kii/(2*HSTEP);
Kij=-Kij/(2*HSTEP);

fprintf('Sparsifying with tolerance %d.\n',k_tol);

% Sparse spring constant matrices
Kii_s=Kii;
Kii_s(abs(Kii_s)<k_tol)=0;
Kii_s=sparse(Kii_s);

Kij_s=Kij;
Kij_s(abs(Kij_s)<k_tol)=0;
Kij_s=sparse(Kij_s);

fclose(fid1);
