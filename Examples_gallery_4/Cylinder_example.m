clear
close all

load('Vorticity_data.mat')
%%
r = 50;
M = r;
ind = (1:M);
X = VORT(:,ind);
Y = VORT(:,ind+1);

% Kernel ResDMD
[G,K,L,PX,PY] = kernel_ResDMD(X,Y,'N',r,'type',"Laplacian");

[W,LAM,W2] = eig(K,'vector');
R = abs(sqrt(real(diag(W2'*L*W2)./diag(W2'*W2)-abs(LAM).^2)));

figure
scatter(real(LAM),imag(LAM),300,R,'.','LineWidth',1);
hold on
% scatter(real(LAM),imag(LAM),250,R,'.');
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'-k')
axis equal
axis([-1.15,1.15,-1.15,1.15])
clim([0,1])
load('cmap.mat')
colormap(cmap2); colorbar
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18)
title(sprintf('Residuals ($M=%d$)',M),'interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

exportgraphics(gcf,sprintf('cylinder_res_M%d.pdf',M),'ContentType','vector','BackgroundColor','none')

figure
loglog([0.001,1],[0.001,1],'k','linewidth',2)
hold on
loglog(sqrt(abs(abs(LAM).^2-1)),R,'b.','markersize',20)
xlabel('$\sqrt{|1-|\lambda|^2|}$','interpreter','latex','fontsize',18)
ylabel('residual','interpreter','latex','fontsize',18)
title(sprintf('Residuals ($M=%d$)',M),'interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

exportgraphics(gcf,sprintf('cylinder_res2_M%d.pdf',M),'ContentType','vector','BackgroundColor','none')

%%


x_pts=-1.2:0.02:1.2;    y_pts=-0.02:0.02:1.2;
z_pts=kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:));    z_pts=z_pts(:);		% complex points where we compute pseudospectra
RES0 = KoopPseudoSpecQR(PX,PY,1/M,z_pts);
RES0=reshape(RES0,length(y_pts),length(x_pts));

RES = KoopPseudoSpec(double(G),double(K),double(L),z_pts,'Parallel','off');	% compute pseudospectra
RES=reshape(RES,length(y_pts),length(x_pts));

%% Plot pseudospectra

figure
hold on
v=(10.^(-10:0.2:0));
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RES0)),log10(v));
hold on
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),-reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RES0)),log10(v));
cbh=colorbar;
cbh.Ticks=log10(10.^(-4:1:0));
cbh.TickLabels=10.^(-4:1:0);
clim([-4,0]);
reset(gcf)
set(gca,'YDir','normal')
colormap gray
axis equal;

title(sprintf('Naive Residual ($M=%d$)',M),'interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)

ax=gca; ax.FontSize=18; axis equal tight;   axis([x_pts(1),x_pts(end),-y_pts(end),y_pts(end)])
hold on
plot(real(LAM),imag(LAM),'.r','markersize',12);
box on
exportgraphics(gcf,sprintf('cylinder_pseudoW_M%d.pdf',M),'ContentType','vector','BackgroundColor','none')




figure
hold on
v=(10.^(-10:0.2:0));
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RES)),log10(v));
hold on
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),-reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RES)),log10(v));
cbh=colorbar;
cbh.Ticks=log10(10.^(-4:1:0));
cbh.TickLabels=10.^(-4:1:0);
clim([-4,0]);
reset(gcf)
set(gca,'YDir','normal')
colormap gray
axis equal;

title(sprintf('Pseudospectrum ($M=%d$)',M),'interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)

ax=gca; ax.FontSize=18; axis equal tight;   axis([x_pts(1),x_pts(end),-y_pts(end),y_pts(end)])
hold on
plot(real(LAM),imag(LAM),'.r','markersize',12);
box on
exportgraphics(gcf,sprintf('cylinder_pseudo_M%d.pdf',M),'ContentType','vector','BackgroundColor','none')


