clear
close all

%% Load the data
% load('HotWireData_FlowInjection.mat')
load('HotWireData_Baseline.mat')

LL=393217;
sig=transpose(Vel);

%% Extract the relevant signals
sig1=sig(4,1:LL);    sig1=(sig1(:)-mean(sig1(:)));
sig2=sig(20,1:LL);    sig2=(sig2(:)-mean(sig2(:)));
sig3=sig(70,1:LL);    sig3=(sig3(:)-mean(sig3(:)));
%% Compute PSD
[F1a,PSD1a] = myPSD(Vel(1:LL,4),Fs,1000);   [F1b,PSD1b] = myPSD(Vel(1:LL,4),Fs,2000);   [F1c,PSD1c] = myPSD(Vel(1:LL,4),Fs,3000);   [F1d,PSD1d] = myPSD(Vel(1:LL,4),Fs,4000);
[F2a,PSD2a] = myPSD(Vel(1:LL,20),Fs,1000);  [F2b,PSD2b] = myPSD(Vel(1:LL,20),Fs,2000);  [F2c,PSD2c] = myPSD(Vel(1:LL,20),Fs,3000);  [F2d,PSD2d] = myPSD(Vel(1:LL,20),Fs,4000);
[F3a,PSD3a] = myPSD(Vel(1:LL,70),Fs,1000);  [F3b,PSD3b] = myPSD(Vel(1:LL,70),Fs,2000);  [F3c,PSD3c] = myPSD(Vel(1:LL,70),Fs,3000);  [F3d,PSD3d] = myPSD(Vel(1:LL,70),Fs,4000);

%% Compute spectral measures
nu1a=MomentMeas(ErgodicMoments(sig1,1000),'filt','cosine'); nu1b=MomentMeas(ErgodicMoments(sig1,2000),'filt','cosine'); nu1c=MomentMeas(ErgodicMoments(sig1,3000),'filt','cosine'); nu1d=MomentMeas(ErgodicMoments(sig1,4000),'filt','cosine');
nu2a=MomentMeas(ErgodicMoments(sig2,1000),'filt','cosine'); nu2b=MomentMeas(ErgodicMoments(sig2,2000),'filt','cosine'); nu2c=MomentMeas(ErgodicMoments(sig2,3000),'filt','cosine'); nu2d=MomentMeas(ErgodicMoments(sig2,4000),'filt','cosine');
nu3a=MomentMeas(ErgodicMoments(sig3,1000),'filt','cosine'); nu3b=MomentMeas(ErgodicMoments(sig3,2000),'filt','cosine'); nu3c=MomentMeas(ErgodicMoments(sig3,3000),'filt','cosine'); nu3d=MomentMeas(ErgodicMoments(sig3,4000),'filt','cosine');

%% Plot results
close all
figure
h=plot(F1a,10*log10(PSD1a(:)),F1b,10*log10(PSD1b(:)),F1c,10*log10(PSD1c(:)),F1d,10*log10(PSD1d(:)),'linewidth',1);
hold on
c = get(h,'Color');
plot(F2a,10*log10(PSD2a(:)),'linewidth',1,'color',c{1})
plot(F2b,10*log10(PSD2b(:)),'linewidth',1,'color',c{2})
plot(F2c,10*log10(PSD2c(:)),'linewidth',1,'color',c{3})
plot(F2d,10*log10(PSD2d(:)),'linewidth',1,'color',c{4})

plot(F3a,10*log10(PSD3a(:)),'linewidth',1,'color',c{1})
plot(F3b,10*log10(PSD3b(:)),'linewidth',1,'color',c{2})
plot(F3c,10*log10(PSD3c(:)),'linewidth',1,'color',c{3})
plot(F3d,10*log10(PSD3d(:)),'linewidth',1,'color',c{4})

set(gca, 'xscale','log')
ax=gca; ax.FontSize=14;
xlim([2,1000])
grid on
set(gca,'XMinorGrid','on');
set(gca,'YMinorGrid','on');
legend({'$N_{\mathrm{ac}}=1000$','$N_{\mathrm{ac}}=2000$','$N_{\mathrm{ac}}=3000$','$N_{\mathrm{ac}}=4000$'},'interpreter','latex','fontsize',14,'location','southwest')
ylim([-60,-10])

figure
h=plot(10.^(-8:0.001:0)*pi/(2*pi*(t(2)-t(1))),10*log10(4*pi*(t(2)-t(1))*nu1a(10.^(-8:0.001:0)*pi)),10.^(-8:0.001:0)*pi/(2*pi*(t(2)-t(1))),10*log10(4*pi*(t(2)-t(1))*nu1b(10.^(-8:0.001:0)*pi)),...
    10.^(-8:0.001:0)*pi/(2*pi*(t(2)-t(1))),10*log10(4*pi*(t(2)-t(1))*nu1c(10.^(-8:0.001:0)*pi)),10.^(-8:0.001:0)*pi/(2*pi*(t(2)-t(1))),10*log10(4*pi*(t(2)-t(1))*nu1d(10.^(-8:0.001:0)*pi)),'linewidth',1);
c = get(h,'Color');
hold on

plot(10.^(-8:0.001:0)*pi/(2*pi*(t(2)-t(1))),10*log10(4*pi*(t(2)-t(1))*nu2a(10.^(-8:0.001:0)*pi)),'linewidth',1,'color',c{1})
plot(10.^(-8:0.001:0)*pi/(2*pi*(t(2)-t(1))),10*log10(4*pi*(t(2)-t(1))*nu2b(10.^(-8:0.001:0)*pi)),'linewidth',1,'color',c{2})
plot(10.^(-8:0.001:0)*pi/(2*pi*(t(2)-t(1))),10*log10(4*pi*(t(2)-t(1))*nu2c(10.^(-8:0.001:0)*pi)),'linewidth',1,'color',c{3})
plot(10.^(-8:0.001:0)*pi/(2*pi*(t(2)-t(1))),10*log10(4*pi*(t(2)-t(1))*nu2d(10.^(-8:0.001:0)*pi)),'linewidth',1,'color',c{4})

plot(10.^(-8:0.001:0)*pi/(2*pi*(t(2)-t(1))),10*log10(4*pi*(t(2)-t(1))*nu3a(10.^(-8:0.001:0)*pi)),'linewidth',1,'color',c{1})
plot(10.^(-8:0.001:0)*pi/(2*pi*(t(2)-t(1))),10*log10(4*pi*(t(2)-t(1))*nu3b(10.^(-8:0.001:0)*pi)),'linewidth',1,'color',c{2})
plot(10.^(-8:0.001:0)*pi/(2*pi*(t(2)-t(1))),10*log10(4*pi*(t(2)-t(1))*nu3c(10.^(-8:0.001:0)*pi)),'linewidth',1,'color',c{3})
plot(10.^(-8:0.001:0)*pi/(2*pi*(t(2)-t(1))),10*log10(4*pi*(t(2)-t(1))*nu3d(10.^(-8:0.001:0)*pi)),'linewidth',1,'color',c{4})
legend({'$N_{\mathrm{ac}}=1000$','$N_{\mathrm{ac}}=2000$','$N_{\mathrm{ac}}=3000$','$N_{\mathrm{ac}}=4000$'},'interpreter','latex','fontsize',14,'location','southwest')
set(gca, 'xscale','log')
ax=gca; ax.FontSize=14;
xlim([2,1000])
grid on
set(gca,'XMinorGrid','on');
set(gca,'YMinorGrid','on');
ylim([-60,-10])

%% Compute spectral measures using different filters for zoom in plot
nu1a1=MomentMeas(ErgodicMoments(sig1,1000),'filt','fejer'); nu1b1=MomentMeas(ErgodicMoments(sig1,2000),'filt','fejer'); nu1c1=MomentMeas(ErgodicMoments(sig1,3000),'filt','fejer'); nu1d1=MomentMeas(ErgodicMoments(sig1,4000),'filt','fejer');
nu2a1=MomentMeas(ErgodicMoments(sig2,1000),'filt','fejer'); nu2b1=MomentMeas(ErgodicMoments(sig2,2000),'filt','fejer'); nu2c1=MomentMeas(ErgodicMoments(sig2,3000),'filt','fejer'); nu2d1=MomentMeas(ErgodicMoments(sig2,4000),'filt','fejer');
nu3a1=MomentMeas(ErgodicMoments(sig3,1000),'filt','fejer'); nu3b1=MomentMeas(ErgodicMoments(sig3,2000),'filt','fejer'); nu3c1=MomentMeas(ErgodicMoments(sig3,3000),'filt','fejer'); nu3d1=MomentMeas(ErgodicMoments(sig3,4000),'filt','fejer');
nu1a3=MomentMeas(ErgodicMoments(sig1,1000),'filt','vand');  nu1b3=MomentMeas(ErgodicMoments(sig1,2000),'filt','vand');  nu1c3=MomentMeas(ErgodicMoments(sig1,3000),'filt','vand');  nu1d3=MomentMeas(ErgodicMoments(sig1,4000),'filt','vand');
nu2a3=MomentMeas(ErgodicMoments(sig2,1000),'filt','vand');  nu2b3=MomentMeas(ErgodicMoments(sig2,2000),'filt','vand');  nu2c3=MomentMeas(ErgodicMoments(sig2,3000),'filt','vand');  nu2d3=MomentMeas(ErgodicMoments(sig2,4000),'filt','vand');
nu3a3=MomentMeas(ErgodicMoments(sig3,1000),'filt','vand');  nu3b3=MomentMeas(ErgodicMoments(sig3,2000),'filt','vand');  nu3c3=MomentMeas(ErgodicMoments(sig3,3000),'filt','vand');  nu3d3=MomentMeas(ErgodicMoments(sig3,4000),'filt','vand');
nu1a4=MomentMeas(ErgodicMoments(sig1,1000));    nu1b4=MomentMeas(ErgodicMoments(sig1,2000));    nu1c4=MomentMeas(ErgodicMoments(sig1,3000));    nu1d4=MomentMeas(ErgodicMoments(sig1,4000));
nu2a4=MomentMeas(ErgodicMoments(sig2,1000));    nu2b4=MomentMeas(ErgodicMoments(sig2,2000));    nu2c4=MomentMeas(ErgodicMoments(sig2,3000));    nu2d4=MomentMeas(ErgodicMoments(sig2,4000));
nu3a4=MomentMeas(ErgodicMoments(sig3,1000));    nu3b4=MomentMeas(ErgodicMoments(sig3,2000));    nu3c4=MomentMeas(ErgodicMoments(sig3,3000));    nu3d4=MomentMeas(ErgodicMoments(sig3,4000));

%%
figure
h=plot(F1d,10*log10(PSD1d(:)),10.^(-8:0.001:0)*pi/(2*pi*(t(2)-t(1))),10*log10(4*pi*(t(2)-t(1))*nu1d1(10.^(-8:0.001:0)*pi)),...
    10.^(-8:0.001:0)*pi/(2*pi*(t(2)-t(1))),10*log10(4*pi*(t(2)-t(1))*nu1d(10.^(-8:0.001:0)*pi)),...
    10.^(-8:0.001:0)*pi/(2*pi*(t(2)-t(1))),10*log10(4*pi*(t(2)-t(1))*nu1d3(10.^(-8:0.001:0)*pi)),...
    10.^(-8:0.001:0)*pi/(2*pi*(t(2)-t(1))),10*log10(4*pi*(t(2)-t(1))*nu1d4(10.^(-8:0.001:0)*pi)),'linewidth',1);
c = get(h,'Color');
hold on
plot(F2d,10*log10(PSD2d(:)),'linewidth',1,'color',c{1})
plot(10.^(-8:0.001:0)*pi/(2*pi*(t(2)-t(1))),10*log10(4*pi*(t(2)-t(1))*nu2d1(10.^(-8:0.001:0)*pi)),'linewidth',1,'color',c{2})
plot(10.^(-8:0.001:0)*pi/(2*pi*(t(2)-t(1))),10*log10(4*pi*(t(2)-t(1))*nu2d(10.^(-8:0.001:0)*pi)),'linewidth',1,'color',c{3})
plot(10.^(-8:0.001:0)*pi/(2*pi*(t(2)-t(1))),10*log10(4*pi*(t(2)-t(1))*nu2d3(10.^(-8:0.001:0)*pi)),'linewidth',1,'color',c{4})
plot(10.^(-8:0.001:0)*pi/(2*pi*(t(2)-t(1))),10*log10(4*pi*(t(2)-t(1))*nu2d4(10.^(-8:0.001:0)*pi)),'linewidth',1,'color',c{5})

plot(F3d,10*log10(PSD3d(:)),'linewidth',1,'color',c{1})
plot(10.^(-8:0.001:0)*pi/(2*pi*(t(2)-t(1))),10*log10(4*pi*(t(2)-t(1))*nu3d1(10.^(-8:0.001:0)*pi)),'linewidth',1,'color',c{2})
plot(10.^(-8:0.001:0)*pi/(2*pi*(t(2)-t(1))),10*log10(4*pi*(t(2)-t(1))*nu3d(10.^(-8:0.001:0)*pi)),'linewidth',1,'color',c{3})
plot(10.^(-8:0.001:0)*pi/(2*pi*(t(2)-t(1))),10*log10(4*pi*(t(2)-t(1))*nu3d3(10.^(-8:0.001:0)*pi)),'linewidth',1,'color',c{4})
plot(10.^(-8:0.001:0)*pi/(2*pi*(t(2)-t(1))),10*log10(4*pi*(t(2)-t(1))*nu3d4(10.^(-8:0.001:0)*pi)),'linewidth',1,'color',c{5})


set(gca, 'xscale','log')
ax=gca; ax.FontSize=14;
xlim([2,10])
grid on
set(gca,'XMinorGrid','on');
set(gca,'YMinorGrid','on');
legend({'PSD','$\varphi_{\mathrm{hat}}$','$\varphi_{\mathrm{cos}}$','$\varphi_{\mathrm{four}}$','$\varphi_{\mathrm{bump}}$'},'interpreter','latex','fontsize',14,'location','west')
ylim([-30,-15])


%% Weak convergence plot
Nvec=[10:10:100,200:100:1000];
MU=ErgodicMoments(sig2,2000);

phi0=chebfun(@(t) exp(sin(t)),[-pi pi],'trig'); % test function
  
for j=1:length(Nvec)
    [F,PSD] = myPSD(sig2',Fs,Nvec(j)); % Power spectral density
    F=F(:);
    integral1(j)=(sum(phi0(2*pi*(t(2)-t(1))*F(1:end-1)).*PSD(1:end-1))+sum(phi0(-2*pi*(t(2)-t(1))*F(2:end)).*PSD(2:end)))*(F(2)-F(1))*2*pi*(t(2)-t(1));
    nu=MomentMeas(MU((2000-Nvec(j)+1):(2000-Nvec(j)+1+2*Nvec(j))),'filt','fejer');
    integral2(j)=4*pi*(t(2)-t(1))*sum(nu*phi0);
    nu=MomentMeas(MU((2000-Nvec(j)+1):(2000-Nvec(j)+1+2*Nvec(j))),'filt','cosine');
    integral3(j)=4*pi*(t(2)-t(1))*sum(nu*phi0);
    nu=MomentMeas(MU((2000-Nvec(j)+1):(2000-Nvec(j)+1+2*Nvec(j))),'filt','vand');
    integral4(j)=4*pi*(t(2)-t(1))*sum(nu*phi0);
    nu=MomentMeas(MU((2000-Nvec(j)+1):(2000-Nvec(j)+1+2*Nvec(j))));
    integral5(j)=4*pi*(t(2)-t(1))*sum(nu*phi0);
end
%%
figure
loglog(Nvec,abs(integral1-integral2(end))/abs(integral5(end)),'-*','linewidth',2)
hold on
loglog(Nvec,abs(integral2-integral5(end))/abs(integral5(end)),'-*','linewidth',2)
loglog(Nvec,abs(integral3-integral5(end))/abs(integral5(end)),'-*','linewidth',2)
loglog(Nvec,abs(integral4-integral5(end))/abs(integral5(end)),'-*','linewidth',2)
loglog(Nvec,max(abs(integral5-integral5(end))/abs(integral5(end)),10^(-16)*4),'-*','linewidth',2)
ax=gca; ax.FontSize=14;

legend({'PSD','$\varphi_{\mathrm{hat}}(x)$','$\varphi_{\mathrm{cos}}(x)$','$\varphi_{\mathrm{four}}(x)$','$\varphi_{\mathrm{bump}}(x)$'},'interpreter','latex','fontsize',14,'location','southeast')
