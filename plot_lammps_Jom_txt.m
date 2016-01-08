
%clear
%file='220514c';
%file='220714b';
% file='210714c4';
%file='101014';
file='181215b';
%file='020215c';
%load(strcat('../koodit/lammps/',file,'_avepos_tar/',file,'_avepos.mat'));
%load(strcat('/wrk/ksaaskil/lammps/',file,'_anharm_tar/',file,'_anharm.mat'),...
%    'Jom_ave1','Jom_ave2','oms_fft','DoS_ave1');
%load(strcat('/wrk/ksaaskil/lammps/',file,'_tar/',file,'_2.mat'));
%load(strcat('/wrk/ksaaskil/lammps/',file,'.mat'));
%load(strcat('/proj/quantum-data/Kimmo/lammps/',file,'_tar/',file,'_Tom.mat'));
% load(strcat(file,'_tar/',file,'_Tom.mat'));

fid=fopen(strcat('/proj/quantum-data/Kimmo/lammps/a-si/',file,'_tar/',file,'_SHC.txt'),'r');

A=textscan(fid,'%f%f');
fclose(fid);

freqs=A{1}/(2*pi)/1e12;
ys=A{2};

figure(242);%clf;
hold on
set(gca,'fontsize',18);

plot(freqs,ys,'-','linewidth',2);
xlim([0,20])

return

%load(strcat(file,'_Tom.mat'));
%load(strcat('/proj/quantum-data/Kimmo/lammps/',file,'_tar/',file,'_Tomega_alpha.mat'));
%return
%load(strcat('/proj/quantum-data/Kimmo/',file,'_Tom_R.mat'));
%load(strcat('/wrk/ksaaskil/lammps/',file,'.mat'));
if 0
load(strcat('../koodit/lammps/',file,'_avepos_tar/',file,'_func.mat'),...
    'Jom_ave1','Jom_ave2','Jom2_ave1','Jom2_ave2','oms_fft','oms_short');
end

run color_palette.m
%dT=30;
%alatt_Si=5.5;
%A=(5*(alatt_Si)*1e-10)^2;

%R=1.36e-9/2;
%R=sqrt(3*10^2)/pi*sqrt(3)*1.44e-10/2;
%R=sqrt(3*5^2)/pi*0.246e-9/2;
A=(5e-9)^2;
%A=pi*R^2;
%A=(7*mean([alatt_Si,alatt_Ge])*1e-10)^2;
% A=(10*1.579*3.4*1e-10)^2;
%A=(abs(0.138814+15.6512)*3.4e-10)^2;
%A=1; 
%dT=4.91;
% Current in GW/m^2/THz
%Jom=-Jom_ave2/A/1e9*1e12;
% Current in W/mK/THz
% Jom=-Jom_ave2/A*1e12/(-Tgrad);
% dT=30;
dT=100;
%dT=20;

k_B=1.38e-23;

%dT=50/5;
%Jom=-Jom_ave1/A*1e12/(dT)/1e6;
%Jom=-Jom_ave1/(k_B*dT);

Jom=-Jom_ave_smooth1/(k_B*dT);
% Standard error
dJom=sqrt(Jom_ave_smooth1_2-Jom_ave_smooth1.^2)/sqrt(k)/(k_B*dT);
% 95 % confidence interval, tinv(0.5+p/2,dof)
dJom=dJom*tinv(.975,k-1);
Qom=-Jom_ave_smooth1;
Qom=-Jom_ave1;
% Standard error
dQom=sqrt(Jom_ave_smooth1_2-Jom_ave_smooth1.^2)/sqrt(k);
% 95 % confidence interval, tinv(0.5+p/2,dof)
dQom=dQom*tinv(.975,k-1);
%Qom=Qom/A*1e12/abs(Tgrad);
%dQom=dQom/A*1e12/abs(Tgrad);

%Current in MW/m^2K/THz
%Jom=-Jom_ave1/A*1e12/1e6;%/(dT);
%Jom=-Jom_ave2/(1.38e-23*100);%*1e12/(dT);
% Current in MW/m^2K/THz
%Jom=-Jom_ave2/(5*1.58*3.4e-10)^2*1e12/1e6;
%Jom=-Jom_ave1/(10*1.58*3.4e-10)^2*1e12/(-Tgrad);
%Jom=-Jom_ave2/(5*1.579*3red).4e-10)^2*1e12/(dT)/1e6;
%Jom_anharm=real(Jom2_ave2)/A/1e6/dT;

f_THz=oms_fft/(2*pi*1e12);
%oms_lj=oms_fft*2.14e-12;
inds=find(f_THz<25);
%inds=inds(1:10:end);
figure(499);%clf;
hold on
set(gca,'fontsize',20); 
color=gr;
%plot(f_THz,real(Tom),'-','linewidth',2,'color',color);
%plot(oms_lj/(2*pi),real(Jom),'-','linewidth',2,'color','g');
pp4=boundedline(f_THz(inds),Jom(inds),[dJom(inds)',dJom(inds)'],'alpha','transparency',0.5,...
    'cmap',repmat(color,3,1));
%pp2=boundedline(f_THz(inds),Qom(inds),[dQom(inds)',dQom(inds)'],'alpha','transparency',0.5,...
%    'cmap',repmat(color,3,1));
set(pp4,'linewidth',1,'color',color);
%inds=find(f_THz<55);
%return
inds=find(f_THz<20);
Jtot1=sum(Jom(inds))*(f_THz(2)-f_THz(1));
fprintf('Total sum is %.2f.\n',Jtot1);
fprintf('Current Q=%g.\n',real(Jtot1)*k_B*dT*1e12) % Metal units
%fprintf('Current Q=%g.\n',real(Jtot1)*dT*1e12)
set(gca,'xlim',[0,20]);
%set(gca,'xlim',[0,65]);
%set(gca,'ylim',[0,10]);
xlabel('Frequency (THz)')
%ylabel('g^{el}(\omega) [MW/(m^2\cdotK\cdotTHz)]')
ylabel('Transmission')
return
if 1
    %%
    hold on;
    k_B=1.38e-23;
    %load('280714a2_airebo_Tom.mat','om_interp','om_mult','Tom');
    %load('290714a_rebo_Tom.mat','om_interp','om_mult','Tom')
    %load('/wrk/ksaaskil/lammps/290714e_Tom_landauer.mat','om_interp','om_mult','Tom')
    %load('290714e_Tom.mat','om_interp','om_mult','Tom','f_landauer')
    %load('/wrk/ksaaskil/lammps/290714e_Tom2.mat','f_landauer','om_interp','om_mult','Tom')
    %load('/proj/quantum-data/Kimmo/260814b_4.Tom.mat','f_landauer','Tom')
    load('/proj/quantum-data/Kimmo/lammps/110914a_2.Tom.mat','f_landauer','Tom')
    %load('/proj/quantum-data/Kimmo/lammps/101014a.Tom.mat','f_landauer','Tom')
    %load('/proj/quantum-data/Kimmo/lammps/011014a2_Tom.mat','om_interp','Tom','om_unit')
    %load('/proj/quantum-data/Kimmo/lammps/011014a2_Tom_tau1e13.mat','om_interp','Tom','om_unit')
    %f_landauer=om_interp*om_unit/1e12/(2*pi);
    %load('/wrk/ksaaskil/lammps/010814a2_Tom.mat','om_interp','om_mult','Tom')
    %plot(f_landauer,k_B*real(Tom)/(A)*1e12/1e6,'-','linewidth',3,'color','k');
    plot(f_landauer,real(Tom),'-','linewidth',3,'color','k');
    %set(gca,'xscale','log');
    return
    if 1 % Find mean-free paths
        %%
        Mom=real(Tom);
        Tom_md=-real(Jom_ave1)/(k_B*dT);
        %figure(2342);%clf;
        %plot(om_interp*om_mult*0.2418,Mom,'b-','linewidth',3);
        %hold on
        %plot(f_THz,Tom_md,'r--','linewidth',3);
        Tom_landauer=@(x)interp1(f_landauer,Mom,x);
   
        Tom_landauer=Tom_landauer(f_THz);

        %G1=Gom_all(:,j);
        %G2=Gom_landauer;
        %L=176;
        L=(860*sqrt(3)-500)*1.42e-10*1e9;
        L=(1100*sqrt(3)-500)*1.42e-10*1e9;
        L=(2320*sqrt(3)-500)*1.42e-10*1e9;
        fprintf('Using L=%.2f nm.\n',L);
        Lambda=L*Tom_md./(Tom_landauer-Tom_md);
        
        figure(2422);clf;
        set(gca,'fontsize',24);
        hold on
        plot(f_THz,Lambda,'-','linewidth',3,'color',bl);
        set(gca,'xlim',[.0,52]);
        set(gca,'ylim',[1e0,1e4]);
        set(gca,'yscale','log');
        ylabel('\Lambda (nm)');
        
        if 1
            %%
            l_klemens=2.5e23/300./(2*pi*f_THz*1e12).^2*1e9;
            l_klemens=l_klemens*1;
            %l_mingo=2e4*1/(32/27*1.8^4*(1.38e-23*500/(12*1.66e-27*4e8))^2*(2*pi*30e12))*1e9;
            l_sec=3e3;
            l_mingo=1e3;
            hold on;
            plot(f_THz,1./(1./l_klemens+1./l_mingo),'r--','linewidth',3);
            %plot(f_THz,1./(1./l_klemens),'k--','linewidth',3);
            
            
            
        end
        return
        if 1 % Determine scattering times
            %%
            fs_1=f_THz(f_THz<55);
            fs_1=fs_1(15:20:end);
            Lambda_1=Lambda(f_THz<55);
            Lambda_1=Lambda_1(15:20:end);
            load /wrk/ksaaskil/lammps/290714e_cnt.Kij.mat
            K0=Kii;
            KV=Kij;
            vel_ave=bandstructure_vels(K0,KV,linspace(.01,pi,500),fs_1);
            %%
            figure(433);%clf;
            hold on
            set(gca,'fontsize',24);
            taus=Lambda_1*1e-9./vel_ave_max'*1e12;
            plot(fs_1,taus,'r--','linewidth',3);
            xlabel('Frequency (THz)');
            ylabel('\tau (ps)');
            set(gca,'xscale','log');
            set(gca,'yscale','log');
            set(gca,'xlim',[.25,55]);
            xdata=log(fs_1(2:100));
            ydata=real(log(taus(2:100)));
            a=polyfit(xdata,ydata,1)
            hold on
            %plot(fs_1,exp(a(2))*fs_1.^(a(1)),'k-','linewidth',3);
            plot(fs_1,exp(a(2))*fs_1.^(-1),'k-','linewidth',3);
            %%
            figure(2422);hold on;
           
            plot(fs_1,exp(a(2))*fs_1.^(-1).*vel_ave'*1e-3,'m-','linewidth',3);
        end
        
    end
    
end

%%
if 0 % Cumulative distribution
    figure(75675);%clf,
    hold on;
    fs=f_THz(inds);
    Jom_fs=real(Jom(inds));
    plot(fs,cumsum(Jom_fs)/sum(Jom_fs),'m-','linewidth',3);
end
%return
%ylabel('q(\omega)/A (GW/m^2/THz)')
%ylabel('q(\omega)/(-A\partial_xT) (W/mK/THz)')
%%
if 0
    %Jom_anharm=real(Jom2_ave2)/A/1e6/dT;
    figure(232342);clf;
    set(gca,'fontsize',24);
    pcolor(oms_short/(2*pi*2.14e-12)/1e12,f_THz_2,real(Jom2_ave2)'/A/1e6/dT);
    shading flat;
    axis equal;
    colorbar
    
    
    axis tight
    Jtot2=sum(sum(Jom_anharm,1))*(f_THz(2)-f_THz(1))^2
    fprintf('Total harmonic+anharmonic sum=%.2f\n',abs(Jtot1)+Jtot2);
    
end