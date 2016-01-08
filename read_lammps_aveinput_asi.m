%file='150914a_restart';
% file='130815n';
%path='~/Documents/koodit/lammps/250214a_avepos_tar/';
%path='/wrk/ksaaskil/lammps/';
%path=strcat('/wrk/ksaaskil/lammps/liquid-solid/');
% path=strcat('/proj/quantum-data/Kimmo/lammps/a-si/',file,'_tar/');
path=strcat('/proj/quantum-data/Kimmo/lammps/a-si/',file,'_tar/');
%path=strcat('/proj/quantum-data/Kimmo/lammps/fpu_chain/')
%return
%path=strcat('/wrk/ksaaskil/lammps/',file,'_anharm_tar/');

filehot=strcat(path,file,'.aveinput_hot.dat');
filecold=strcat(path,file,'.aveinput_cold.dat');
%filehot=strcat(path,file,'.aveinput_hot_start.dat');
%filecold=strcat(path,file,'.aveinput_cold_start.dat');

fid1=fopen(filehot,'r');
fid2=fopen(filecold,'r');

dt_md=2.5e-15;
%dt_md=1.0e-15;
%dt_md=0.002*2.14e-12;
% =0.005;
%dt_md=1.0e-15;

% Pairs Np
ss=textscan(fid1,'%f%f','headerlines',2);

ts=ss{1};
ts=ts-ts(1);
curr_hot=-ss{2};

ss=textscan(fid2,'%f%f','headerlines',2);

ts2=ss{1};
ts2=ts2-ts2(1);
curr_cold=ss{2};

% Energy is in units of eV
curr_hot=curr_hot*1.602e-19;
curr_cold=curr_cold*1.602e-19;

% Energy is in units of kcal/mol
%curr_hot=curr_hot*4.2e3/6.022e23;
%curr_cold=curr_cold*4.2e3/6.022e23;

% Energy is in units of epsilon for LJ
%curr_hot=curr_hot*1.67e-21;
%curr_cold=curr_cold*1.67e-21;
if 1
    figure(23223);clf
    set(gca,'fontsize',24);
    if 1
        Qdiff=diff(curr_hot)/(ts(2)-ts(1))/dt_md;
        %g1=ones(10000,1);
        %g1=g1/sum(g1);
        %Qsmooth=conv(Qdiff,g1,'same');
        %Nave=100000;
        Nave=10;
        Nmax=floor((length(ts)-1)/Nave)*Nave;
        Qdiff=reshape(Qdiff(1:Nmax),Nave,Nmax/Nave);
        Qdiff=mean(Qdiff,1);
        %plot(ts(2:Nave:Nmax)*dt_md*1e12,Qdiff,'ro-');
        plot(ts(2:Nave:Nmax)*dt_md,Qdiff,'ro-');
        hold on
        %plot(ts(2:end)*dt_md,Qsmooth,'k-');
        %return
        QH=mean(Qdiff);
        Qdiff_H=Qdiff;
        %plot(ts(1:Nmax)*dt_md*1e12,ones(Nmax,1)*QH,'r-','linewidth',3);
        plot(ts(1:Nmax)*dt_md,ones(Nmax,1)*QH,'r-','linewidth',3);
        
        Qdiff=diff(curr_cold)/(ts(2)-ts(1))/dt_md;
        Qdiff=reshape(Qdiff(1:Nmax),Nave,Nmax/Nave);
        Qdiff=mean(Qdiff,1);
        hold on
        %plot(ts(2:Nave:Nmax)*dt_md*1e12,Qdiff,'bo-');
        plot(ts(2:Nave:Nmax)*dt_md,Qdiff,'bo-');
        QC=mean(Qdiff);
        Qdiff_C=Qdiff;
        %plot(ts(1:Nmax)*dt_md*1e12,ones(Nmax,1)*QC,'b-','linewidth',3);
        plot(ts(1:Nmax)*dt_md,ones(Nmax,1)*QC,'b-','linewidth',3);
        %return
    else
        %plot(ts*dt_md/1e-12,curr_hot./(ts*dt_md),'b-','linewidth',3);
        plot(ts(2:end)*dt_md/1e-12,diff(curr_hot),'b-','linewidth',3);
        hold on;
        %plot(ts2*dt_md/1e-12,curr_cold./(ts2*dt_md),'r--','linewidth',3);
        plot(ts2(2:end)*dt_md/1e-12,diff(curr_cold),'r--','linewidth',3);
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        QH=curr_hot(end)/(dt_md*ts(end));

        QC=curr_cold(end)/(dt_md*ts(end));
    end
    %xlabel('t (ps)')
    %ylabel('P (W)')
    xlabel('t');
    ylabel('P');
    
end



fclose(fid1);
fclose(fid2);



Qs=diff(curr_hot+curr_cold)/2/((ts(2)-ts(1))*dt_md);

varQ=var(Qs);
dQ=sqrt(varQ)/sqrt(length(Qs));

%R=sqrt(3*10^2)/pi*sqrt(3)*1.44e-10/2;
%R=sqrt(3*5^2)/pi*0.246e-9/2;
%A=2*pi*R*0.34e-9;
fprintf('The current is %.3e.\n',mean([QH,QC]));
Q=mean([QH,QC]);
%dT=100
%G=mean([QH,QC])/A/dT
%dG=dQ/A/dT
return
kappa=-mean([QH,QC])/(A*Tgrad)
dkappa=-dQ/(A*Tgrad)

kappas(1)=kappa;
dkappas(1)=dkappa;
%%
%figure(5675);
%plot(ts(2:end)*dt_md,Qs);

%QH/A/dT

return
%%
load(strcat(path,file,'.mat'));

T1=mean(vel_ave(1:end/2));
T2=mean(vel_ave(end/2+1:end))*1;
dT=(T1-T2);%*40/0.3305;
f_THz=oms_fft/(2*pi)/1e12;

inds=find(f_THz<18);
%dom=(oms_fft(2)-oms_fft(1))*om0;
df=f_THz(2)-f_THz(1);
a_latt=1.5496*3.4e-10;
%A=64*a_latt^2;
% Current per THz
J_THz=Jom_ave2(inds)*1e12;
%J=sum(-J_Thz)*dom/(2*pi)
J=sum(-J_THz)*df
%dT=0.3305/1/5;
A=(5*alatt_Si*1e-10)^2;

% Flux GW/m^2
JA=J/A/1e9

return
kB=1.38e-23;
GJ=kB*J/dT/A;
GQ=kB*QH/dT/A;
GJ=GJ/1e6;
GQ=GQ/1e6*om0;
ind=2;
GJs(ind)=GJ;
GQs(ind)=GQ;

return
Np=ss{2}

% dt_timestep
ss=textscan(fid1,'%s%d',1);
% Read the number of entries
d_timesteps=ss{2}
dt=double(d_timesteps)*dt_md

ss=textscan(fid1,'%s%d',1);
% Read N_left
N_left=ss{2};
ss=textscan(fid1,'%s%d',1);
% Read N_right
N_right=ss{2};

% Read the pair ids
ss=textscan(fid1,'%s%s',1); % "Pair ids:"
ss=textscan(fid1,'%d%d%d%d%d',Np);
ids1=ss{2};
ids2=ss{3};
ss=textscan(fid1,'%s',1);


if 1
    disp('Reading forces...')
    data_f=textscan(fid1,'%.12f');
    data_f=data_f{1};
    data_f=reshape(data_f,3*Np,length(data_f)/(3*Np));
end

fclose(fid1);

%return
% Pairs
% Pairs Np
ss=textscan(fid2,'%s%d',1);
Np=ss{2}

% dt_timestep
ss=textscan(fid2,'%s%d',1);
% Read the number of entries
d_timesteps=ss{2}
dt=double(d_timesteps)*dt_md

ss=textscan(fid2,'%s%d',1);
% Read N_left
N_left=ss{2};
ss=textscan(fid2,'%s%d',1);
% Read N_right
N_right=ss{2};

% Read the pair ids
ss=textscan(fid2,'%s%s',1); % "Pair ids:"
ss=textscan(fid2,'%d%d%d%d%d',Np);
ids1=ss{2};
ids2=ss{3};
ss=textscan(fid2,'%s',1);


if 1
    disp('Reading velocities...')
    data_v=textscan(fid2,'%.12f');
    data_v=data_v{1};
    data_v=reshape(data_v,3*Np,length(data_v)/(3*Np));
end

fclose(fid2);
vels=data_v';
forces=data_f';

%% Now time runs in columns
dims=2;
heatcurr=sum(mean(forces.*vels,1))*.5

%dt=0.1;
dT=0.1;

force_fft=fft(forces)*dt;
vel_fft=fft(vels)/2*dt;

oms_fft=(0:size(force_fft,1)-1)/(size(force_fft,1)*dt)*2*pi;

Tom_fft=2*sum(real(force_fft.*conj(vel_fft)),2)/(dt*length(force_fft));

%vel1_fft=fft(vels(:,dims))*dt;
%vel2_fft=fft(vels(:,3+dims))*dt;

%Tom_fft2=2*imag(vel1_fft.*conj(vel2_fft))./(oms_fft'*length(force_fft)*dt);

%vec=vels(:,4)*dt;
%Tom_fft=2*real(fft(vec).*conj(fft(vec)))/(dt*length(f_fft)*mean(vec.^2/dt^2));
Tom_fft=Tom_fft/dT;
%Tom_fft2=Tom_fft2/dT;
%%
dom=oms_fft(2)-oms_fft(1);
win=round(.05/dom);

g = gausswin(win); % <-- this value determines the width of the smoothing window
g = g/sum(g);
Tom_ave = conv(Tom_fft, g, 'same');
%Tom_ave2 = conv(Tom_fft2, g, 'same');

% Need to multiply by 4/3 for some reason?
figure(2227);clf;
set(gca,'fontsize',24);
%plot(oms_fft,Tom_fft,'b-');
hold on
plot(oms_fft,Tom_ave,'r-','linewidth',3);
%plot(oms_fft,Tom_ave2,'g-.','linewidth',3);
set(gca,'xlim',[0,3]);
set(gca,'ylim',[0,5]);
xlabel('\omega');
