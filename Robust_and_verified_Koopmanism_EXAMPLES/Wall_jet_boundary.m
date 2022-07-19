clear
close all

%% The following data files are available upon request
% load('dataset_1s.mat') % this loads DATA from experimental realisation #1
% load('dataset_1s_b.mat') % this loads DATA from experimental realisation #2
load('EDMD_canopy_final.mat') % this is saved data from running the commented code below

%%%%% UNCOMMENT THE FOLLOWING TO PERFORM COMPUTATIONS FROM DATA %%%%%
% 
% %% Algorithmic parameters
% N=2000;		% number of basis functions used
% M1=N;		% number of snapshots to compute the basis
% M2=11999;	% number of snapshots used for ResDMD matrices
% ind1=1:M1;  ind2=1:M2;
% use_DMD=1;	% basis choice: set to 1 to use DMD, set to 0 to use kEDMD
% comp_pseudospec=1;	% set to 1 to compute and show pseudospectra
% %% Apply ResDMD
% if use_DMD~=1
%     tic
%     [PSI_x,PSI_y] = kernel_ResDMD(DATA(:,ind1),DATA(:,ind1+1),DATA2(:,ind2),DATA2(:,ind2+1),'N',N,'sketch','off','Parallel','on','cut_off',10^(-12));
%     toc
% else
%     [~,S,V]=svd(transpose(DATA(:,ind1))/sqrt(M1),'econ');
%     PSI_x=transpose(DATA2(:,ind2))*V(:,1:N)*diag(1./diag(S(1:N,1:N)));
%     PSI_y=transpose(DATA2(:,ind2+1))*V(:,1:N)*diag(1./diag(S(1:N,1:N)));
%     clear S V
% end
%% Compute the ResDMD matrices
G_matrix=(PSI_x(1:M2,1:N)'*PSI_x(1:M2,1:N))/M2;    
A_matrix=(PSI_x(1:M2,1:N)'*PSI_y(1:M2,1:N))/M2;
L_matrix=(PSI_y(1:M2,1:N)'*PSI_y(1:M2,1:N))/M2;

%% (uncomment to) Compute pseudospectra
% if comp_pseudospec==1
%     x_pts=-1.5:0.05:1.5;
%     y_pts=x_pts;
%     z_pts=kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:)); z_pts=z_pts(:);
%     RES = reshape(KoopPseudoSpec(G_matrix,A_matrix,L_matrix,z_pts,'Parallel','on'),length(y_pts),length(x_pts));
% end
[V,D]=eig(A_matrix,G_matrix);   E=diag(D);	% EDMD eigenvalues
RES2=real(sqrt(dot(V,L_matrix*V-A_matrix'*V*D-A_matrix*V*D'+G_matrix*V*abs(D).^2)./dot(V,G_matrix*V))); % residual of the eigenpairs
%%
if comp_pseudospec==1
    figure
    v=(10.^(-4:0.2:0));
	contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RES)),log10(v));
    hold on
	cbh=colorbar;
	cbh.Ticks=log10([0.005,0.01,0.1,1]);
	cbh.TickLabels=[0,0.01,0.1,1];
	caxis([log10(0.01),log10(1)]);
	reset(gcf)
	set(gca,'YDir','normal')
	colormap bone
	th=-pi:0.001:pi;
	c1=exp(1i*th-th.^2/18); c2=exp(1i*th-th.^2/30-0.005-abs(th)/20);    c3=exp(1i*th-abs(th)/12);
	plot(real(c3),imag(c3),'-b','linewidth',2)

	ax=gca; ax.FontSize=14; axis equal tight;   axis([x_pts(1),x_pts(end),y_pts(1),y_pts(end)])
	plot(real(E),imag(E),'.r');
end


%%%%% UNCOMMENT THE FOLLOWING TO COMPUTE KOOPMAN MODES (REQUIRES RAW DATA) %%%%%

% %% Koopman mode decomposition
% sigu=DATA2(1:51150,ind2);       sigu=sigu(:,1:M2);
% sigv=DATA2(51151:end,ind2);     sigv=sigv(:,1:M2);
% sig_vel=sqrt(sigu.^2+sigv.^2);
% 
% I_p=find(RES2<0.5);
% K_modes_res_u=pinv(PSI_x(1:M2,1:N)*V(:,I_p))*transpose(sigu);
% K_modes_res_v=pinv(PSI_x(1:M2,1:N)*V(:,I_p))*transpose(sigv);
% K_modes_res_vel=pinv(PSI_x(1:M2,1:N)*V(:,I_p))*transpose(sig_vel);
% 
% K_modes_u=pinv(PSI_x(1:M2,1:N)*V)*transpose(sigu);  
% K_modes_v=pinv(PSI_x(1:M2,1:N)*V)*transpose(sigv);
% 
% %% Plot specific modes
% th=0.01;0.05;
% close all
% 
% figure
% h=tiledlayout(1,3);
% 
% lam=0.9439+0.2458i; % spectral parameter
% mode=find(abs(E-lam)==min(abs(E-lam)));
% lam=E(mode);
% clc
% RES2(mode)
% tt=norm(PSI_x(1:M2,1:N)*V(:,I_p==mode))/sqrt(M2);
% 
% XIu=reshape((K_modes_res_u(I_p==mode,:)),[310,165])*tt;
% myContours = linspace(min(real(XIu(:))),max(real(XIu(:))), 21);
% nexttile
% contourf(Xgrid,Ygrid,(real(XIu)),myContours,'edgecolor','none')
% set(gca,'ydir','normal')
% colorbar
% colormap(brighten(redblueTecplot(21),-0.6))
% 
% XIv=reshape((K_modes_res_v(I_p==mode,:)),[310,165])*tt;
% myContours = linspace(min(real(XIv(:))),max(real(XIv(:))), 21);
% nexttile
% contourf(Xgrid,Ygrid,(real(XIv)),myContours,'edgecolor','none')
% set(gca,'ydir','normal')
% colorbar
% colormap(brighten(redblueTecplot(21),-0.6))
% 
% XIvel=reshape((K_modes_res_vel(I_p==mode,:)),[310,165])*tt;
% myContours = linspace(min(real(XIvel(:))),max(real(XIvel(:))), 21);
% nexttile
% contourf(Xgrid,Ygrid,(real(XIvel)),myContours,'edgecolor','none')
% set(gca,'ydir','normal')
% colorbar
% colormap(brighten(redblueTecplot(21),-0.6))
% 
% 
% %% Energy test
% sige=sum(sigu.^2+sigv.^2);
% K_modes_energy=pinv(PSI_x(1:M2,1:N)*V)*transpose(sige(1:M2));
% 
% %% error of decomposition
% en=mean(sige);
% en_k=zeros(M2,1);
% for j=1:M2
%     en_k(j)=sum(abs(PSI_x(1:M3,1:N0)*V*diag(E.^(j-1))*K_modes_energy),1)/M3;
%     figure(100)
%     semilogy(en_k*0+en,'--k')
%     hold on
%     semilogy(abs(en_k),'b')
%     hold off
%     ylim([10^(-6),1000000])
%     pause(0.1)
% 
% end

%%%%% UNCOMMENT THE FOLLOWING TO COMPUTE WAVENUMBER SPECTRA AND SPECTRAL MEASURES (REQUIRES RAW DATA) %%%%%

% %% Wavenumber spectra and spectral measures
% 
% % compute wavenumber spectra
% for j=1:12000
%     u_field=reshape(transpose(DATA2(1:51150,j)),[310,165]);
%     v_field=reshape(transpose(DATA2(51151:end,j)),[310,165]);
%     
%     sig1=u_field(:,105);sig1=sig1(:);
%     sig2=u_field(:,160);sig2=sig2(:);
%     if j==1
%         MU1a=ErgodicMoments(sig1-mean(sig1),309);
%         MU2a=ErgodicMoments(sig2-mean(sig2),309);
%         
%     else
%         MU1a=MU1a+ErgodicMoments(sig1-mean(sig1),309);
%         MU2a=MU2a+ErgodicMoments(sig2-mean(sig2),309);
%     end
%     
%     sig1=v_field(:,105);sig1=sig1(:);
%     sig2=v_field(:,160);sig2=sig2(:);
%     if j==1
%         MU1b=ErgodicMoments(sig1-mean(sig1),309);
%         MU2b=ErgodicMoments(sig2-mean(sig2),309);
%         
%     else
%         MU1b=MU1b+ErgodicMoments(sig1-mean(sig1),309);
%         MU2b=MU2b+ErgodicMoments(sig2-mean(sig2),309);
%     end
% end
% %% plot wavenumber spectra
% nu1a=MomentMeas(MU1a/12000,'filt','cosine');    nu2a=MomentMeas(MU2a/12000,'filt','cosine');
% nu1b=MomentMeas(MU1b/12000,'filt','cosine');    nu2b=MomentMeas(MU2b/12000,'filt','cosine');
% 
% figure
% theta=-pi/10*0:0.0001:pi;
% lambda=1./(theta/(2*pi*0.25*0.001));
% kx=2*pi./lambda;
% semilogx(kx,4*pi*(2*pi*0.25*0.001)*abs(nu1a(theta).*kx),'linewidth',2)
% hold on
% plot(kx,4*pi*(2*pi*0.25*0.001)*abs(nu2a(theta).*kx),'linewidth',2)
% xlim([10,10000])
% set(gca,'fontsize',14)
% set(gca,'xminorgrid','on')
% set(gca,'yminorgrid','on')
% 
% figure
% theta=-pi/10*0:0.0001:pi;
% lambda=1./(theta/(2*pi*0.25*0.001));
% kx=2*pi./lambda;
% semilogx(kx,4*pi*(2*pi*0.25*0.001)*abs(nu1b(theta).*kx),'linewidth',2)
% hold on
% plot(kx,4*pi*(2*pi*0.25*0.001)*abs(nu2b(theta).*kx),'linewidth',2)
% xlim([10,10000])
% set(gca,'fontsize',14)
% set(gca,'xminorgrid','on')
% set(gca,'yminorgrid','on')
% 
% %% compute spectral measures
% Y2=Ygrid(:);
% I1=find(Y2==Ygrid(1,105));
% I2=find(Y2==Ygrid(1,160));
% 
% u_data1=DATA2(I1,:);
% v_data1=DATA2(I1+51150,:);
% 
% u_data2=DATA2(I2,:);
% v_data2=DATA2(I2+51150,:);
% 
% for j=1:length(I1)
%     sigu1=u_data1(j,:);sigu1=sigu1(:);
%     sigu2=u_data2(j,:);sigu2=sigu2(:);
%     
%     sigv1=v_data1(j,:);sigv1=sigv1(:);
%     sigv2=v_data2(j,:);sigv2=sigv2(:);
%     if j==1
%         MUu1=ErgodicMoments(sigu1-mean(sigu1),10000);
%         MUu2=ErgodicMoments(sigu2-mean(sigu2),10000);
%         MUv1=ErgodicMoments(sigv1-mean(sigv1),10000);
%         MUv2=ErgodicMoments(sigv2-mean(sigv2),10000); 
%     else
%         MUu1=MUu1+ErgodicMoments(sigu1-mean(sigu1),10000);
%         MUu2=MUu2+ErgodicMoments(sigu2-mean(sigu2),10000);
%         MUv1=MUv1+ErgodicMoments(sigv1-mean(sigv1),10000);
%         MUv2=MUv2+ErgodicMoments(sigv2-mean(sigv2),10000);
%     end
% end
% MUu1=MUu1/length(I1);   MUu2=MUu2/length(I1);   MUv1=MUv1/length(I1);   MUv2=MUv2/length(I1);
% 
% %% plot spectral measures
% nu1=MomentMeas(MUu1((10001-100):10001+100),'filt','cosine');
% nu2=MomentMeas(MUu2((10001-100):10001+100),'filt','cosine');
% 
% figure
% plot(10.^(-8:0.001:0)*pi/(2*pi*(delta_t)),10*log10(4*pi*(delta_t)*nu1(10.^(-8:0.001:0)*pi)),'linewidth',2)
% hold on
% plot(10.^(-8:0.001:0)*pi/(2*pi*(delta_t)),10*log10(4*pi*(delta_t)*nu2(10.^(-8:0.001:0)*pi)),'linewidth',2)
% set(gca, 'xscale','log')
% ax=gca; ax.FontSize=14;
% xlim([10,5000])
% grid on
% set(gca,'XMinorGrid','on');
% set(gca,'YMinorGrid','on');
% 
% nu1=MomentMeas(MUv1((10001-100):10001+100),'filt','cosine');
% nu2=MomentMeas(MUv2((10001-100):10001+100),'filt','cosine');
% 
% figure
% plot(10.^(-8:0.001:0)*pi/(2*pi*(delta_t)),10*log10(4*pi*(delta_t)*nu1(10.^(-8:0.001:0)*pi)),'linewidth',2)
% hold on
% plot(10.^(-8:0.001:0)*pi/(2*pi*(delta_t)),10*log10(4*pi*(delta_t)*nu2(10.^(-8:0.001:0)*pi)),'linewidth',2)
% set(gca, 'xscale','log')
% ax=gca; ax.FontSize=14;
% xlim([10,5000])
% grid on
% set(gca,'XMinorGrid','on');
% set(gca,'YMinorGrid','on');
