clear
close all
addpath(genpath('./data_from_runs'))

load('van_der_pol_data') % load the trajectory data


% %% Compute ResDMD matrices
% [PX,PY1,PY2] = kernel_ResDMD(Xa,Ya,Xb,Yb,'Y2',Y2,'parallel','off','N',N);
% %%
% G = (PX'*PX)/M;
% A = PX'*(PY1+PY2)/(2*M);
% L = (PY1'*PY1+PY2'*PY2)/(2*M); L = (L+L')/2;
% H = (PY1'*PY2+PY2'*PY1)/(2*M); H = (H+H')/2; 
% 
% [V,E] = eig(A,G,'vector'); % EDMD eigenvalues

% %% UNCOMMENT TO COMPUTE PSEUDOSPECTRA
% x_pts=-1.5:0.1:1.5;    y_pts=x_pts;
% z_pts=kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:));    z_pts=z_pts(:);
% 
% RES = KoopPseudoSpec(G,A,L,z_pts,'Parallel','on');
% RES = reshape(RES,length(y_pts),length(x_pts));
% %
% RESb = KoopPseudoSpec(G,A,H,z_pts,'Parallel','on');
% RESb = reshape(RESb,length(y_pts),length(x_pts));


%% Pseudospectra plots

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
plot(real(E),imag(E),'.r','markersize',15);
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',16)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',16)
title('$\mathrm{Sp}_\epsilon^{\mathrm{var}}(\mathcal{K}_{(1)})$','interpreter','latex','fontsize',16)
% exportgraphics(gcf,'SDE_var_pseudospec.pdf','Resolution',1000);

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
plot(real(E),imag(E),'.r','markersize',15);
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',16)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',16)
title('$\mathrm{Sp}_\epsilon(\mathcal{K}_{(1)})$','interpreter','latex','fontsize',16)
% exportgraphics(gcf,'SDE_pseudospec.pdf','Resolution',1000);

%% Plot EDMD eigenfunctions
clear LAMBDA
clear RESIDUALS
clear VARIANCE

kk=8;
for m=0:2
    for k=0:4
    close all
    lam=exp(-m*mu*delta_t+k*1i*delta_t*(1-mu^2/16));
    [~,I]=sort(abs(E-lam),'ascend');
    lam=E(I(1));
    v=V(:,I(1));
    v=v(:,1)/sqrt(abs(v(:,1)'*G*v(:,1))); % normalise
    LAMBDA(m+1,k+kk+1)=lam;
    RESIDUALS(m+1,k+kk+1) =sqrt(abs(v'*(H-lam*A'-lam'*A+abs(lam)^2*G)*v)/abs(v'*G*v));
    VARIANCE(m+1,k+kk+1)=sqrt(abs(v'*(L-lam*A'-lam'*A+abs(lam)^2*G)*v)/abs(v'*G*v));
   
    C=(PX*v(:,1));
    CC=1;
    scatter(Xb(1,:)',Xb(2,:)',2,real(C),'filled');
    clim([-CC,CC]);
    colorbar
    xlabel('$X_1$','interpreter','latex','fontsize',16)
    ylabel('$X_2$','interpreter','latex','fontsize',16)
    if k==0
        title(sprintf('$\\lambda= %5.3f$',real(lam)),'fontsize',16,'interpreter','latex');
    else
        if imag(lam)>0
            title(sprintf('$\\lambda= %5.3f+%5.3fi$',real(lam),imag(lam)),'fontsize',16,'interpreter','latex');
        else
            title(sprintf('$\\lambda= %5.3f-%5.3fi$',real(lam),-imag(lam)),'fontsize',16,'interpreter','latex');
        end
    end
    
    T = sprintf('SDE_mode_m%dk%d.pdf',m,k);

    ax = gca; ax.FontSize = 16;
    axis tight
    box on
    colormap(coolwarm)
    grid on

    % exportgraphics(gcf,T,'ContentType','image','BackgroundColor','none','Resolution',200)
    pause(0.1)

    end
end

% %% Subspace error
% 
% ww=1000;
% 
% c1=PX\(Xb(1,:)');
% C1=c1;
% c2=PX\(Xb(2,:)');
% C2=c2;
% K=G\A;
% 
% NT=ceil(ww/delta_t);
% 
% E1=zeros(NT,1);
% E2=zeros(NT,1);
% x1_pred=E1;
% x2_pred=E2;
% 
% XX=PX(100,:);
% 
% 
% for j=1:NT
%     e1 = sqrt(max(c1'*H*c1-2*real(c1'*K'*A*c1)+c1'*K'*G*K*c1,0));
%     e2 = sqrt(max(c2'*H*c2-2*real(c2'*K'*A*c2)+c2'*K'*G*K*c2,0));
%     if j==1
%         E1(j)=e1;
%         E2(j)=e2;
%     else
%         E1(j)=E1(j-1)+e1;
%         E2(j)=E2(j-1)+e2;
%     end
%     c1=K*c1;
%     c2=K*c2;
%     x1_pred(j)=XX*c1;
%     x2_pred(j)=XX*c2;
% end

%%
close all
figure
plot((1:NT)*delta_t,E1/sqrt(C1'*G*C1),'linewidth',1);
hold on
plot((1:NT)*delta_t,E2/sqrt(C2'*G*C2),'linewidth',1);
xlabel('$t$','interpreter','latex','fontsize',16)
legend({'$X_1$','$X_2$'},'interpreter','latex','fontsize',16,'location','northwest')
title('$\delta_n(X_i)$','interpreter','latex','fontsize',16)
xlim([0,1000])
ax = gca; ax.FontSize = 16;
% exportgraphics(gcf,'SDE_subspace.pdf','Resolution',1000);
%%

figure
plot((1:NT)*delta_t,var(y(1:NT,1:2:end),0,2),'linewidth',1)
hold on
plot((1:NT)*delta_t,var(y(1:NT,2:2:end),0,2),'linewidth',1)
xlabel('$t$','interpreter','latex','fontsize',16)
legend({'$X_1$','$X_2$'},'interpreter','latex','fontsize',16,'location','northwest')
title('$\mathrm{Var}(X_i(\textbf{x}_n))$','interpreter','latex','fontsize',16)
xlim([0,1000])
ax = gca; ax.FontSize = 16;
% exportgraphics(gcf,'SDE_variance.pdf','Resolution',1000);
%%

figure
plot((1:NT)*delta_t,x1_pred,'r','linewidth',2);
hold on
plot((1:NT)*delta_t,mean(y(1:NT,1:2:end),2),'k','linewidth',1)
xlabel('$t$','interpreter','latex','fontsize',16)
legend({'$K^nX_1$','$\mathcal{K}^nX_1$'},'interpreter','latex','fontsize',16)
xlim([0,1000])
ax = gca; ax.FontSize = 16;
% exportgraphics(gcf,'SDE_pred_X1.pdf','Resolution',1000);


figure
plot((1:NT)*delta_t,x2_pred,'r','linewidth',2);
hold on
plot((1:NT)*delta_t,mean(y(1:NT,2:2:end),2),'k','linewidth',1)
xlabel('$t$','interpreter','latex','fontsize',16)
legend({'$K^nX_2$','$\mathcal{K}^nX_2$'},'interpreter','latex','fontsize',16)
xlim([0,1000])
ax = gca; ax.FontSize = 16;
% exportgraphics(gcf,'SDE_pred_X2.pdf','Resolution',1000);





