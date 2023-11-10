clear
close all
addpath(genpath('./data_from_runs'))

load('circle_data.mat')

%%%% UNCOMMENT TO RUN
% %% Parameters
% c=1/5;
% f = @(theta) (1/2)*sin(2*pi*theta)/(2*pi);
% F = @(theta) theta + c +(rand(size(theta,1),1)-0.5)/(5*pi)+f(theta);
% 
% %% Dictionary
% n = 20; 
% Phi = [chebfun(@(t) exp(2*pi*1i*((-n:n))*t),[0 1],'periodic')];
% 
% %% Form ResDMD matrices
% % M=10^7;
% M1=4*n;
% M2=2*10^4;
% M=M1*M2;
% xM = linspace(0,1,M1+1).'; xM(end) = [];
% xM = repmat(xM,M2,1);
% yM = mod(F(xM),1);
% yMb = mod(F(xM),1);
% 
% PX = Phi(xM,:); 
% PY1 = Phi(yM,:);
% PY2 = Phi(yMb,:);
% 
% G=eye(2*n+1);
% A = (1/M)*PX'*(PY1+PY2)/2;
% L = (1/M)*(PY1'*PY1+PY2'*PY2)/2;
% H = (1/M)*(PY1'*PY2+PY2'*PY1)/2;
% 
% 
% %% Convergence of ResDMD matrices in Frobenius norm
% 
% Mvec=round(10.^(2:0.1:4));
% E1=0*Mvec;
% E2=0*Mvec;
% E3=0*Mvec;
% E4=0*Mvec;
% 
% for jj=1:length(Mvec)
%     jj
%     M2=Mvec(jj);
%     M=M1*M2;
%     xM = linspace(0,1,M1+1).'; xM(end) = [];
%     xM = repmat(xM,M2,1);
%     yM = mod(F(xM),1);
%     yMb = mod(F(xM),1);
% 
%     PXb = PX(1:M,:);%Phi(xM,:); 
%     PY1b = PY1(1:M,:);%Phi(yM,:);
%     PY2b = PY2(1:M,:);%Phi(yMb,:);
% 
%     G2=(1/M)*(PXb'*PXb);
%     A2 = (1/M)*PXb'*(PY1b+PY2b)/2;
%     L2 = (1/M)*(PY1b'*PY1b+PY2b'*PY2b)/2;
%     H2 = (1/M)*(PY1b'*PY2b+PY2b'*PY1b)/2;
% 
%     E1(jj) = norm(G-G2,'fro'); E2(jj) = norm(A-A2,'fro'); E3(jj) = norm(L-L2,'fro'); E4(jj) = norm(H-H2,'fro');
% end
% 
% %% Compute eigenvalues and residuals
% 
% [V,lam] = eig(A,'vector');
% lam_var = 0*lam;
% lam_res = 0*lam;
% for j=1:length(lam)
%     lam_var(j)=sqrt(max(real(V(:,j)'*(L-lam(j)*A'-conj(lam(j))*A+abs(lam(j))^2*G)*V(:,j))/norm(V(:,j))^2,0));
%     lam_res(j)=sqrt(max(real(V(:,j)'*(H-lam(j)*A'-conj(lam(j))*A+abs(lam(j))^2*G)*V(:,j))/norm(V(:,j))^2,0));
% end


% %% Compute pseudospectra
% x_pts=-1.5:0.01:1.5;
% y_pts=x_pts;
% z_pts=kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:));
% z_pts=z_pts(:);
% RES = KoopPseudoSpec(G,A,L,z_pts,'Parallel','on');
% RES = reshape(RES,length(y_pts),length(x_pts));
% %
% RESb = KoopPseudoSpec(G,A,H,z_pts,'Parallel','on');
% RESb = reshape(RESb,length(y_pts),length(x_pts));


%% Plot results

v=(10.^(-3:0.1:0));

figure
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RES)),log10(v));
cbh=colorbar;
cbh.Ticks=log10([0.005,0.01,0.1,1]);
cbh.TickLabels=[0,0.01,0.1,1];
clim([log10(0.01),0]);
reset(gcf)
set(gca,'YDir','normal')
colormap gray
ax=gca; ax.FontSize=16; axis equal tight;   axis([x_pts(1),x_pts(end),y_pts(1),y_pts(end)])
hold on
plot(real(lam),imag(lam),'.r','markersize',15);
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',16)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',16)
title('$\mathrm{Sp}_\epsilon^{\mathrm{var}}(\mathcal{K}_{(1)})$','interpreter','latex','fontsize',16)
exportgraphics(gcf,'circle_var_pseudospec.pdf','Resolution',1000);
%
figure
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(max(real(RESb),0.01)),log10(v));
cbh=colorbar;
cbh.Ticks=log10([0.005,0.01,0.1,1]);
cbh.TickLabels=[0,0.01,0.1,1];
clim([log10(0.01),0]);
reset(gcf)
set(gca,'YDir','normal')
colormap gray
ax=gca; ax.FontSize=16; axis equal tight;   axis([x_pts(1),x_pts(end),y_pts(1),y_pts(end)])
hold on
plot(real(lam),imag(lam),'.r','markersize',15);
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',16)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',16)
title('$\mathrm{Sp}_\epsilon(\mathcal{K}_{(1)})$','interpreter','latex','fontsize',16)
exportgraphics(gcf,'circle_pseudospec.pdf','Resolution',1000);
%%

% Covariance plot
COV = L -H;
COV(1,1)=1;

figure
imagesc([-n,n],[-n,n],min(max(log10(abs(COV)),-10),0))
title('$|\tilde{L}-\tilde{H}|$','interpreter','latex','fontsize',16)
cbh=colorbar;
cbh.Ticks=-10:2:0;
cbh.TickLabels={'1e-10','1e-8','1e-6','1e-4','1e-2','1'};
ax=gca; ax.FontSize=16; axis equal tight;
colormap hot
exportgraphics(gcf,'circle_covariance.pdf','Resolution',1000);


figure
plot(abs(lam),lam_var,'.','markersize',20)
hold on
plot(abs(lam),abs(lam_res),'.','markersize',20)
xlabel('$|\lambda|$','interpreter','latex','fontsize',16)
legend({'$\mathrm{res}^{\mathrm{var}}(\lambda,g)$','$\mathrm{res}(\lambda,g)$'},'interpreter','latex','fontsize',16)
ax = gca; ax.FontSize = 16;
xlim([0,1])
exportgraphics(gcf,'circle_residual.pdf','Resolution',1000);
%%
figure
loglog(Mvec,E2,'.','markersize',20)
hold on
loglog(Mvec,E3,'.','markersize',20)
loglog(Mvec,E4,'.','markersize',20)
loglog(Mvec(6:end-5),2./sqrt(Mvec(6:end-5)),'k','linewidth',2)
text(500,0.04,'$\mathcal{O}(M_2^{-1/2})$','interpreter','latex','fontsize',16)
xlabel('$M_2$','interpreter','latex','fontsize',16)
legend({'$\|\tilde{A}-A\|_{\mathrm{Fr}}$','$\|\tilde{L}-L\|_{\mathrm{Fr}}$','$\|\tilde{H}-H\|_{\mathrm{Fr}}$'},'interpreter','latex','fontsize',16)
ax = gca; ax.FontSize = 16;
ylim([0.01,1])
exportgraphics(gcf,'circle_convergence.pdf','Resolution',1000);
