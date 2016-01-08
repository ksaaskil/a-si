load 181215_Lambdas.mat

freqs_fit=freqs2';
Lambdas_mod=Lambdas2;
Mom_mod=Mom2;
Mom_mod(freqs_fit<2)=Mom_mod(freqs_fit<2)*0.95;
area=(70e-10)^2;

Ts=linspace(1,500,50);
Ls_est=linspace(1,1000,1);

Ls=[1:10,20,50,75,100]*1e-9;
N_files=length(Ls);
dT_MD=100*ones(1,N_files);

% T=300;
hconst=6.626e-34;
kB=1.38e-23;
kappa_est=zeros(length(Ls_est),length(Ts));
kappa_qm=zeros(length(Ls_est),length(Ts));


for k=1:length(Ts)
    T=Ts(k);
    fs=hconst*freqs_fit*1e12/(kB*T);
    weights=fs.^2.*exp(fs)./(exp(fs)-1).^2;

    for j=1:length(Ls_est)
       L=Ls_est(j);
       Tomega=Mom_mod./(1+L./(2*Lambdas_mod));
       kappa_qm(j,k)=L*1e-9*sum((Tomega(1:end-1).*weights(1:end-1)'+...
           Tomega(2:end).*weights(2:end)')/2.*diff(freqs_fit'))*1e12;
       kappa_qm(j,k)=kappa_qm(j,k)*1.38e-23/area;
       kappa_est(j,k)=L*1e-9*sum((Tomega(1:end-1)+Tomega(2:end))/2.*diff(freqs_fit'))*1e12;
       kappa_est(j,k)=kappa_est(j,k)*1.38e-23/area;
       kappa_spectrum=L*1e-9*(Tomega(1:end-1).*weights(1:end-1)'+...
           Tomega(2:end).*weights(2:end)')/2*1e12*1.38e-23/area;
    end
end
figure(5454);% clf;
hold on
plot(freqs_fit(1:end-1),kappa_spectrum,'-','linewidth',2);
xlabel('Frequency (THz)');
ylabel('$\kappa(\omega)$ (W$^{-1}$m$^{-1}$K$^{-1}$THz$^{-1}$)','interpreter','latex');
set(gca,'fontsize',18);
return
figure(3358);clf;
hold on
p1=errorbar(Ls*1e9,currs.*Ls./dT_MD,1.96*dcurrs.*Ls./dT_MD,'bo','linewidth',2);
set(gca,'fontsize',18);
plot(Ls_est,kappa_est(:,1),'r-','linewidth',2);
xlabel('Length (nm)');
ylabel('\kappa (W/mK)');
% plot(Ls_est,kappa_qm,'b-','linewidth',2);
xlabel('Length (nm)','interpreter','latex');
ylabel('Thermal conductivity (Wm$^{-1}$K$^{-1}$)','interpreter','latex')
legend('NEMD','\kappa_{est}');
box on
% return
figure(3379);%clf;
hold on
set(gca,'fontsize',18);
plot(Ts,kappa_qm(end,:),'-','linewidth',2,'color',gr);
xlabel('Temperature (K)','interpreter','latex');
set(gca,'fontsize',18);
ylabel('Thermal conductivity (Wm$^{-1}$K$^{-1}$)','interpreter','latex');
return
fid=fopen('~/matlab/a-si/experiment_data/cahill_data.csv','r');
A=textscan(fid,'%f%f','delimiter',',');
fclose(fid);
hold on
plot(A{1},A{2},'ks','markerfacecolor','k','markersize',5);