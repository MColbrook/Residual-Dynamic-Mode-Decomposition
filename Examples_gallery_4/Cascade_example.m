clear
close all


%% Load and process the data
load('pressure_data.mat')
DATA = DATA-mean(DATA,1);
%% Run new DMD with residuals
M = 700;    r = M;

X = DATA(1:M,:)'/sqrt(M);
Y = DATA(2:(M+1),:)'/sqrt(M);

[U,S,V] = svd(X,'econ');
r = min(rank(S),r);
U = U(:,1:r); V = V(:,1:r); S = S(1:r,1:r); Sinv = diag(1./diag(S));
K = (U')*Y*V*Sinv; % DMD matrix

PX2 = X*V*diag(1./diag(S));
PY2 = Y*V*diag(1./diag(S));
[W,LAM,W2] = eig(K,'vector');
Phi = Y*V*Sinv*W; % DMD modes

R = vecnorm(PY2*W-PX2*W*diag(LAM))./vecnorm(PX2*W); % residuals
%%
[~,I]=sort(abs(1-LAM),'ascend'); % reorder modes
R = R(I); LAM = LAM(I); W = W(:,I); W2 = W2(:,I); Phi = Phi(:,I);

%%

figure
scatter(real(LAM),imag(LAM),300,R,'.','LineWidth',1);
hold on
plot(cos(-1:0.01:2*pi),sin(-1:0.01:2*pi),'-k')
xlim([0.95,1.01])
clim([0,max(R(real(LAM)>0.95))])
load('cmap.mat')
colormap(cmap2); colorbar
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18)
title('Residuals','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

% exportgraphics(gcf,'cascade_res1.pdf','ContentType','vector','BackgroundColor','none')

figure
loglog([0.01,1],[0.01,1],'k','linewidth',2)
hold on
loglog(sqrt(abs(abs(LAM).^2-1)),R,'b.','markersize',20)
xlabel('$\sqrt{|1-|\lambda|^2|}$','interpreter','latex','fontsize',18)
ylabel('residual','interpreter','latex','fontsize',18)
title('Residuals','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

% exportgraphics(gcf,'cascade_res2.pdf','ContentType','vector','BackgroundColor','none')
return

%% Plot DMD modes

close all
sz2=abs(X_data-0.1);
sz2=max(0,sz2-0.04);
sz2=3-2.1*exp(-(sz2*5));

for ind = 1:11
    u = real(Phi(:,ind));
    res = R(ind);
    L = LAM(ind);
    if imag(L)<0
    else
        figure;
        scatter(X_data,Z_data,sz2,u,'o','filled');
        hold on
        scatter(X_data,Z_data+max(Z_data(:)),sz2,u,'o','filled');
        title(sprintf('$\\lambda=%0.4f+%0.4fi,r=%0.4f$',real(L),imag(L),res),'interpreter','latex','fontsize',18)
        load('cmap.mat')
        colormap(cmap)
        clim([-3*std(u)+mean(u),3*std(u)+mean(u)])
        % xlim([-0.1,0.3])
        axis equal tight
        axis off
        % exportgraphics(gcf,sprintf('cascade_mode%d.png',ind),'ContentType','image','BackgroundColor','none')
        % close all
    end
    pause(1)
    % close all
end
