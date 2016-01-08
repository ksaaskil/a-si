clear

files={'181215a','181215b','181215c','181215d','181215e','181215f',...
    '181215g','181215h','181215i','181215j','170415b','170415c','150515a','150515b'};
N_files=length(files);
Ls=[1:10,20,50,75,100]*1e-9;
dT_MD=100*ones(1,N_files);


run color_palette.m

currs_int=zeros(1,N_files);

flag_aveinput=1;
if flag_aveinput==1
    currs=zeros(1,N_files);
    dcurrs=zeros(1,N_files);
end

area=(7e-9)^2;

for iter=1:length(files)
    file=files{iter};
    %folder=strcat('/wrk/ksaaskil/lammps/',file,'_tar/');
    folder=strcat('/proj/quantum-data/Kimmo/lammps/a-si/',file,'_tar/');
    
    fid=fopen(strcat(folder,file,'_SHC.txt'),'r');

    A=textscan(fid,'%f%f');
    fclose(fid);
    
    freqs=A{1}/(2*pi*1e12);
    SHC=A{2};% /area;
    Nfreqs=length(freqs);
    %return
    
    if iter==1
        Toms=zeros(Nfreqs,length(files));
        kappaoms=zeros(Nfreqs,length(files));
        dkappaoms=zeros(Nfreqs,length(files));
        dToms=zeros(Nfreqs,length(files));
    end

    
    Toms(:,iter)=SHC/(1.38e-23*dT_MD(iter));
    %kappaoms(:,iter)=Toms(:,iter)*Ls(iter);
    %dToms(:,iter)=sqrt(real(Jom2_ave_smooth1_2-Jom2_ave_smooth1.^2))/dT_MD(iter)/sqrt(k);
    %dToms(:,iter)=dToms(:,iter)*tinv(0.975,k-1);
    %dkappaoms(:,iter)=dToms(:,iter)*Ls(iter);
    % return
    
    %inds=find(oms_fft<3);
    currs_int(iter)=sum(Toms(freqs<18,iter))*dT_MD(iter)*(freqs(2)-freqs(1))*1e12/area*1.38e-23;
    
    if flag_aveinput
       run read_lammps_aveinput_asi;
       currs(iter)=mean([QH,QC])/area;
       dcurrs(iter)=dQ/area;
       
    end
    
    
end

clear curr_cold curr_hot Jom* ts ts2 Jsmooth Jsmooth2
inds=find(freqs<18);
Toms=Toms(inds,:);
dToms=dToms(inds,:);
kappaoms=kappaoms(inds,:);
dkappaoms=dkappaoms(inds,:);
freqs=freqs(inds);
% oms=oms_fft(inds);
% inds=inds(1:2:end); % Necessary for 230615 files, where only the half read
% return
% clear oms_fft
%%
figure(3349);clf;
%plot(Ls,currs.*Ls,'bo-','linewidth',2);
p2=errorbar(Ls*1e9,currs.*Ls./dT_MD,1.96*dcurrs.*Ls./dT_MD,'bo','linewidth',2);
ylabel('Thermal conductivity (W/mK)')
xlabel('Length (nm)');
set(gca,'fontsize',18);
hold on;
plot(Ls*1e9,currs_int.*Ls./dT_MD,'rs','linewidth',2);
xdata=log(Ls(end-5:end));
ydata=log(currs(end-5:end).*Ls(end-5:end)./dT_MD(end-5:end));
%xdata=log(Ls(1:3));
%ydata=log(currs(1:3).*Ls(1:3)./dT_MD(1:3));
a=polyfit(xdata,ydata,1)
% plot(Ls,exp(a(2))*Ls.^(a(1)),'k--','linewidth',2);
%plot(Ls,exp(a(2))*Ls.^(0.33),'m--','linewidth',2);
% set(gca,'xscale','log');set(gca,'yscale','log');
% set(gca,'xlim',[min(Ls),max(Ls)])
return
%%
figure(233);clf;
set(gca,'fontsize',18);

% plot(freqs,Toms(:,[1,5,10]),'-','linewidth',2);
% plot(freqs,Toms(:,[1,5,10])*1.38e-23*100/area*1e12,'-','linewidth',1);

plot(freqs,Toms(:,1)*1.38e-23*100/area*1e12,'-','linewidth',1,'color',bl);
hold on;
plot(freqs,Toms(:,5)*1.38e-23*100/area*1e12,'-','linewidth',1,'color',red);
plot(freqs,Toms(:,10)*1.38e-23*100/area*1e12,'-','linewidth',1,'color',gr);

xlabel('Frequency (THz)');
ylabel('q [W/(m^2K THz)]')
set(gca,'fontsize',18);
xlim([0,17])

%%
if 1 % Mean free paths from curve fitting
    %ind_Ls=[1,2,3,4,5,6,7,8,9];
    %ind_Ls=3:10;
    ind_Ls=[1,2:10]; % Low frequencies, up to 0.5
    % ind_Ls=[1,2:7];
    % ind_Ls=[3:10] % Intermediate frequencies, 0.5-1.5
    % ind_Ls=1:5; % High frequencies, 1.5-2.4
    %ind_Ls=3:7;
    % Figure 181215_kappaest: Ls: 1-7 for 10-15, 1-10 for 0.01-10
    
    Ls_fit=Ls(ind_Ls)*1e9;
    % freqs_fit=linspace(.1,10,5000);
    freqs_fit=linspace(2,10,5000);
    % freqs_fit=linspace(10,15,5000);

    ind_max=length(freqs_fit);
    
    Toms_fit=zeros(length(freqs_fit),length(ind_Ls));
    dToms_fit=Toms_fit;
    
    iter=0;
    for i=ind_Ls
        %i
        iter=iter+1;
        %Tom_func=@(x)interp1(fs',Tom_alpha_mat(alpha_ind,:,i),x);
        Tom_func=@(x)interp1(freqs',Toms(inds,i),x);
        Toms_fit(:,iter)=Tom_func(freqs_fit);
        %dTom_func=@(x)interp1(oms',dToms(inds,i),x);
        %dToms_fit(:,iter)=dTom_func(oms_fit);
    end
    % return
    
    Lambdas=zeros(ind_max,1);
    T0=zeros(ind_max,1);
    dLambdas=Lambdas;
    unity_fit=zeros(ind_max,1);
    
    err_as=zeros(ind_max,1);
    err_bs=err_as;
    as=err_as;
    bs=as;
    
    plot_or_not=(ind_max<6);
    
    if(plot_or_not);
        fit_window=7457;
        figure(fit_window);
        clf;
    end
    
    % %%
    for j=1:ind_max
       xdata=Ls_fit;
       ydata=1./Toms_fit(j,:);
       % yerrdata=1./Toms_fit(j,:).^2.*dToms_fit(j,:);
       if 1
           a=polyfit(xdata,ydata,1);
       else
           X=[xdata',ones(length(Ls_fit),1)];
           Y=ydata';
           Ws=1./yerrdata.^2;
           [a,da]=lscov(X,Y,Ws);
           
           %[a,da]=lscov(X'*W*X,X'*Y,ones(size(X'*X,1),1));
           % da=da*3;
           %[a,da]=lscov(X,Y,ones(length(Ls_fit),1));
       end
       %return
       as(j)=a(1);
       bs(j)=a(2); % Should be one
       T0(j)=1./a(2);
       Lambdas(j)=bs(j)./as(j)/2;
       
       x_mean=mean(xdata);
       y_est=a(1)*xdata+a(2);
       n=length(xdata);
       dof=n-2;
       s_a1=sqrt(sum((y_est-ydata).^2)/dof)/sqrt(sum((xdata-x_mean).^2));
       s_a2=s_a1*sqrt(sum(xdata.^2)/n);
       t_factor=tinv(0.975,dof);

       as(j)=a(1);
       err_as(j)=t_factor*s_a1;
       err_bs(j)=t_factor*s_a2;
       %return
       dLambdas(j)=(abs(err_bs(j)/as(j))+abs(bs(j)/as(j)^2*err_as(j)))/2;
       % dLambdas(j)=2*(abs(da(2)/a(1))+abs(a(2)/a(1)^2*da(1)))/2; % Division by two from MFP def.
       %dLambdas(j)=1;
       if (plot_or_not)
           colors=[bl;red;gr;[0,0,0]];
          color=colors(j,:);
          %figure(fit_window);
          set(gca,'fontsize',16);
          plot(xdata,ydata,'o','color',color,'markersize',5,'markerfacecolor',color);
          %dydata=1./Toms_fit(j,:).^2.*dToms_fit(j,:)*1.96;
          %p3=errorbar(xdata,ydata,dydata,'o','markersize',6);% ,'markerfacecolor','g','linewidth',2); %,dydata,'o','color',color);
          hold on;
          %plot(xdata,a(1)*xdata+a(2),'k--','linewidth',2);
          xplot=linspace(1,max(xdata),100);
          plot(xplot,a(1)*xplot+a(2),'k--','linewidth',2);
          set(gca,'xlim',[0,max(Ls_fit)]);
          xlabel('1/N');
          ylabel('[q(\omega)/\DeltaT]^{-1}')
          %plot(Ls,a(1)*Ls+a(2),'k-','linewidth',3);
          %return
       end
    
    end

    if ind_max>10
        figure(7999);%clf;
        color=bl;
        [~,ind_fit1]=min(abs(freqs_fit-2));
        [~,ind_fit2]=min(abs(freqs_fit-5));
        spacing=10;
        inds_fit=ind_fit1:spacing:ind_fit2;

        xdata=log(freqs_fit(inds_fit));
        xdata_mean=mean(xdata);
        ydata=real(log(Lambdas(inds_fit)))';
        a=polyfit(xdata,ydata,1);    

        fprintf('Fitted the slope %f pm %f.\n',-a(1),0); 
        
        xdata=freqs_fit(inds_fit).^(-2);
        ydata=Lambdas(inds_fit)';
        a2=polyfit(xdata,ydata,1);

        if 0
            pp1=boundedline(freqs_fit,Lambdas,[dLambdas,dLambdas],'alpha','transparency',0.5,...
                'cmap',repmat(color,3,1));
            % Fit a line
            hold on;
            plot(freqs_fit,exp(a(2))*freqs_fit.^(a(1)),'r--','linewidth',2);
            % set(gca,'xscale','log');%set(gca,'yscale','log');
        else % Log-log plot
            ub=log10(Lambdas+dLambdas)-log10(Lambdas);
            lb=log10(Lambdas)-log10(Lambdas-dLambdas);

            [hh1,h1]=boundedline(freqs_fit,log10(Lambdas),...
                [lb,ub],...
                'transparency',.5,'cmap',repmat(color,3,1));
            hold on;
            %set(gca,'xlim',[min(oms_fit),max(oms_fit)])
            fit1=a(2)+a(1)*log(freqs_fit);
            %p2=plot(oms_fit,log10(exp(a(2))*oms_fit.^(a(1))),'k-','linewidth',2);
            % p2=plot(freqs_fit,log10(exp(fit1)),'k-','linewidth',2);
            p3=plot(freqs_fit,log10(a2(1)*freqs_fit.^(-2)),'g-','linewidth',2);
            set(gca,'xscale','log')
            set(gca,'ylim',[-1,1.5])
            set(gca,'fontsize',18);
            set(gca,'ytick',-2:2);
            set(gca,'yticklabel',{'10^{-2}','10^{-1}','10^0','10^1','10^2'});
            xlabel('Frequency (THz)');
            ylabel('\Lambda(\omega) (nm)')
            % set(gca,'ylim',[2,5])
            set(gca,'xlim',[min(freqs_fit),max(freqs_fit)]);

        end
    
    end
    
    
end


%% Estimate kappa from Lambdas
Ls_est=linspace(1,100);
kappa_est=zeros(length(Ls_est),1);
Lambdas_mod=Lambdas;

% Fit a line for T0
[~,ind_fit1]=min(abs(freqs_fit-.65));
[~,ind_fit2]=min(abs(freqs_fit-1.5));
%xdata=log(freqs_fit(ind_fit1:ind_fit2))';
%ydata=log(T0(ind_fit1:ind_fit2));
xdata=freqs_fit(ind_fit1:ind_fit2)'.^2;
ydata=T0(ind_fit1:ind_fit2);
a_T0=polyfit(xdata,ydata,1);
figure(232);hold on;
%plot(freqs_fit(ind_fit1:ind_fit2),exp(a(2))*freqs_fit(ind_fit1:ind_fit2).^(a(1)),'k--',...
%    'linewidth',2);
plot(freqs_fit,a_T0(1)*freqs_fit.^2,'k-');
Mom_mod=T0;
Mom_mod(freqs_fit<2)=a_T0(1)*freqs_fit(freqs_fit<2).^2;
%Mom_mod(oms_fit<.01)=0;
figure(2347);clf;
set(gca,'fontsize',18);
plot(freqs_fit,Mom_mod,'b-','linewidth',2);
title('Extrapolated T0');
% Lambdas_mod(freqs_fit<1.5)=exp(a(2))*freqs_fit(freqs_fit<1.5).^(a(1));
% Lambdas_mod(freqs_fit<1.5)=exp(a(2))*freqs_fit(freqs_fit<1.5).^(-2);
Lambdas_mod(freqs_fit<2)=a2(1)*freqs_fit(freqs_fit<2).^(-2);
%Lambdas_mod(oms_fit>1)=1e3./oms_fit(oms_fit>1).^2;
figure(2348);clf
set(gca,'fontsize',18);
hold on;
plot(freqs_fit,Lambdas_mod,'r-','linewidth',2);
title('Extrapolated MFP');
set(gca,'xscale','log');
set(gca,'yscale','log');
%Mom_mod(fs_fit<.25)=4;
% return
for j=1:length(Ls_est)
   L=Ls_est(j);
   Tomega=Mom_mod./(1+L./(2*Lambdas_mod));
   kappa_est(j)=L*1e-9*sum((Tomega(1:end-1)+Tomega(2:end))/2.*diff(freqs_fit'))*1e12;
   kappa_est(j)=kappa_est(j)*1.38e-23/area;
end
% return
%% 
figure(3348);%clf;
hold on
p1=errorbar(Ls*1e9,currs.*Ls./dT_MD,1.96*dcurrs.*Ls./dT_MD,'bo','linewidth',2);
set(gca,'fontsize',18);
plot(Ls_est,kappa_est+kappa_est2,'b-','linewidth',2);
% +kappa_est2
return

