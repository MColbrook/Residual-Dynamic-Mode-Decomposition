clear
close all

%% The following data files are available from the dropbox link

use_DMD=1;  % choice of basis: set to 1 to use DMD, set to 0 to use kEDMD
if use_DMD==1 % linear dictionary
    load('Cylinder_DMD.mat') % precomputed matrices
elseif use_DMD==2 % combined
    load('Cylinder_DMD.mat')
    PSI_x0=PSI_x;   PSI_y0=PSI_y;
    load('Cylinder_EDMD.mat')
    N=2*N; % non-linear dictionary
    PSI_x=[PSI_x,PSI_x0];
    PSI_y=[PSI_y,PSI_y0];    
else
    load('Cylinder_EDMD.mat') % precomputed matrices
end

%%%%% UNCOMMENT THE FOLLOWING TO PERFORM COMPUTATIONS FROM DATA FILE 'Cylinder_data.mat' %%%%%
% load('Cylinder_data.mat');
% %% Algorithmic parameters
% N=200;      % number of basis functions used
% M1=500;     % number of snapshots to compute the basis
% M2=1000;    % number of snapshots used for ResDMD matrices
% 
% ind1=(1:M1)+6000; % 5000 is roughly where we reach post transient regime
% ind2=(1:M2)+ind1(end)+500;
% 
% %% Apply ResDMD
% if use_DMD~=1
%     [~,~,~,~,~,PSI_x,PSI_y] = kernel_ResDMD(DATA(:,ind1),DATA(:,ind1+1),'Xb',DATA(:,ind2),'Yb',DATA(:,ind2+1),'N',N,'type','Laplacian');
% else
%     [~,S,V]=svd(transpose(DATA(:,ind1))/sqrt(M1),'econ');
%     PSI_x=transpose(DATA(:,ind2))*V(:,1:N)*diag(1./(diag(S(1:N,1:N))));
%     PSI_y=transpose(DATA(:,ind2+1))*V(:,1:N)*diag(1./(diag(S(1:N,1:N))));
%     clear S V
% end
% 
% %%
G_matrix=(PSI_x'*PSI_x)/M2;     
A_matrix=(PSI_x'*PSI_y)/M2;
L_matrix=(PSI_y'*PSI_y)/M2;

%% Compute pseudospectra
x_pts=-1.5:0.05:1.5;    y_pts=x_pts;
z_pts=kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:));    z_pts=z_pts(:);		% complex points where we compute pseudospectra
RES = KoopPseudoSpec(G_matrix,A_matrix,L_matrix,z_pts,'Parallel','on');	% compute pseudospectra
RES=reshape(RES,length(y_pts),length(x_pts));
[V,D]=eig((G_matrix)\(A_matrix));   E=diag(D);  % EDMD eigenvalues

%% Plot pseudospectra
figure
hold on
v=(10.^(-2:0.2:0));
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RES)),log10(v));
cbh=colorbar;
cbh.Ticks=log10([0.005,0.01,0.1,1]);
cbh.TickLabels=[0,0.01,0.1,1];
caxis([log10(0.01),log10(1)]);
reset(gcf)
set(gca,'YDir','normal')
colormap bone
ax=gca; ax.FontSize=14; axis equal tight;   axis([x_pts(1),x_pts(end),y_pts(1),y_pts(end)])
hold on
plot(real(E),imag(E),'.r');

%%
RES2 = zeros(N,1);
for j=1:N
    RES2(j)=abs(sqrt(V(:,j)'*(L_matrix-E(j)*A_matrix'-conj(E(j))*A_matrix+abs(E(j))^2*G_matrix)*V(:,j))/sqrt(V(:,j)'*G_matrix*V(:,j)));
end
%% Plot EDMD eigenvalues against residual
[RES_p,I]=sort(RES2,'ascend');
figure
semilogy(angle(E(I)),RES_p,'.r','markersize',10)
set(gca,'yscale','log')
ax=gca; ax.FontSize=14;
ylim([10^(-8),1])

%% Check the lattice structure
evec_x=PSI_x*V(:,I);
lam=E(I);
t1=0.967585567481353 + 0.252543401421919i;
t1=lam(abs(lam-t1)==min(abs(lam-t1)));
lam1=zeros(100,1);  ang1=zeros(100,1);  res1=zeros(100,1);
for j=1:100
    I2=find(abs(lam-t1^j)<max(0.001,0*max(lam1)));
    if abs(length(I2)-1)>0
        break
    end
    b1=evec_x(:,abs(lam-t1)<max(0.001,0*max(lam1))).^j;
    b2=evec_x(:,I2);
    ang1(j)=acos(abs(b1(:)'*b2(:)/(norm(b1)*norm(b2))));
    lam1(j)=abs(lam(I2)-t1^j);
    res1(j)=RES2(I(I2));
end
%%
figure
lam1(1)=0;  ang1(1)=0;  res1(1)=0;
semilogy(lam1,'*-','linewidth',2)
hold on
semilogy(res1,'d-','linewidth',2)
semilogy(real(ang1),'s-','linewidth',2)
% semilogy([35,35],[10^(-14),1],'--k','linewidth',2)
legend({'$|\lambda_j-\lambda_1^j|$','$\mathrm{res}(\lambda_j,g_j)$','eigenspace error'},'interpreter','latex','fontsize',22,'location','southeast')
ax=gca; ax.FontSize=19;
ylim([10^(-14),1])
yticks(10.^(-14:2:0));
xlim([0,100])

%%%%% UNCOMMENT THE FOLLOWING TO COMPUTE KOOPMAN MODES (REQUIRES RAW DATA) %%%%%

% %% Plot Koopman modes
% myContours1=cell(20,1); myContours2=cell(20,1);
% for pow=[1,2,20]
%     lambda=t1^pow;
%     Idd=find(abs(E-lambda)==min(abs(E-lambda)));
%     tt=norm(PSI_x*V(:,Idd))/sqrt(M2);
%     XI=pinv(V)*pinv(PSI_x)*transpose(DATA(1:size(DATA,1)/2,ind2));   XI=XI(Idd,:); XI=-1i*reshape(XI,[400,100])*tt;
% 
%     h=figure;
%     subplot1=subplot(2,1,1);
%     d = 2*obst_r;
%     myContours1{pow} = linspace(min(real(XI(:))),max(real(XI(:))), 21);
%     contourf((x-obst_x)/d,(y-obst_y)/d,real(XI),myContours1{pow},'edgecolor','none')
%     colorbar
%     axis equal tight
%     colormap(viridisCMap)
%     hold on
%     fill(obst_r*cos(0:0.01:2*pi)/d,obst_r*sin(0:0.01:2*pi)/d,[200,200,200]/255,'edgecolor','none')
%     xlim([-2,max((x(:)-obst_x)/d)])
%     T=sprintf('Mode %d (real part)',pow);
%     title(T,'interpreter','latex','fontsize',16)
%     set(subplot1,'BoxStyle','full','DataAspectRatio',[1 1 1],'Layer','top',...
%     'PlotBoxAspectRatio',[8.25 2.25 1]);
%     caxis([min(myContours1{pow}),max(myContours1{pow})])
%     
%     subplot2=subplot(2,1,2);
%     d = 2*obst_r;
%     myContours2{pow} = linspace(0,max(abs(XI(:))), 21);
%     contourf((x-obst_x)/d,(y-obst_y)/d,abs(XI),myContours2{pow},'edgecolor','none')
%     colorbar
%     axis equal tight
%     colormap(viridisCMap)
%     hold on
%     fill(obst_r*cos(0:0.01:2*pi)/d,obst_r*sin(0:0.01:2*pi)/d,[200,200,200]/255,'edgecolor','none')
%     xlim([-2,max((x(:)-obst_x)/d)])
%     T=sprintf('Mode %d (absolute value)',pow);
%     title(T,'interpreter','latex','fontsize',16)
%     set(subplot2,'BoxStyle','full','DataAspectRatio',[1 1 1],'Layer','top',...
%     'PlotBoxAspectRatio',[8.25 2.25 1]);
%     caxis([min(myContours2{pow}),max(myContours2{pow})])
% 
%     h.Position =[  360.0000  262.3333  560.0000  355.6667];
% end

