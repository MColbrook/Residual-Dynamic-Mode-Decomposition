clear
close all
rng(1)

load('Plasma_data.mat')
delay=10;
pG = pG(randperm(65),:);

figure
plot(1:size(pG,2),pG(1,:),':k')
hold on
plot(1:size(pG,2),pG(61:65,:),'r')
plot(1:size(pG,2),pG(1:60,:),':k')
plot(1:size(pG,2),pG(61:65,:),'r')
xlim([1,110])
xlabel('time step','interpreter','latex','fontsize',18)
legend({'Snapshot Data','Test Data'},'interpreter','latex','fontsize',14)
ax=gca; ax.FontSize=18;
exportgraphics(gcf,'LIP_DATA.pdf','ContentType','vector','BackgroundColor','none')

%% Extract the data matrices using delay embedding
X = []; Y = [];
for ii=1:60
    DATA=pG(ii,1:(end-delay+1));
    for j=2:delay
        DATA=[DATA;pG(ii,j:(end-delay+j))];
    end
    X = [X,DATA(:,1:end-1)];    Y = [Y,DATA(:,2:end)];
end

test_DATA=[];
for ii=61:65
    test_DATA=[test_DATA,pG(ii,1:delay)'];
end

%%
[G,K,L,PX2,PY2,PX3] = kernel_ResDMD(X,Y,'type',"Gaussian",'N',200,'Xb',test_DATA);

%%
[W,LAM,W2] = eig(K,'vector');
R = (sqrt(real(diag(W2'*L*W2)./diag(W2'*W2)-abs(LAM).^2)));
[~,I] = sort(R,'ascend');
W = W(:,I); LAM = LAM(I); W2 = W2(:,I); R = R(I);

PV1 = PX2*W; PV2 = PY2*W; PV3 = PX3*W;

%%
close all
figure
scatter(real(LAM),imag(LAM),300,R,'.','LineWidth',1);
hold on
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'-k')
axis equal
axis([-1.15,1.15,-1.15,1.15])
clim([0,max(R(real(LAM)>0))])
load('cmap.mat')
colormap(cmap2); colorbar
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18)
title('Residuals','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
exportgraphics(gcf,'LIP_res.pdf','ContentType','vector','BackgroundColor','none')

figure
loglog([0.001,1],[0.001,1],'k','linewidth',2)
hold on
loglog(sqrt(abs(abs(LAM).^2-1)),R,'b.','markersize',20)
xlabel('$\sqrt{|1-|\lambda|^2|}$','interpreter','latex','fontsize',18)
ylabel('residual','interpreter','latex','fontsize',18)
title('Residuals','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
exportgraphics(gcf,'LIP_res2.pdf','ContentType','vector','BackgroundColor','none')


%% Test prediction
sig = X;
load('LIP_times.mat')
t=t-t(1);

nvec = 1:length(LAM);
Er = zeros(length(nvec),2);

for ct = 1:min(length(nvec),100)
    ct
    nn = nvec(ct);
    K2 = PX2(:,1:nn)\PY2(:,1:nn);
    for ii=61:65
        coeffs=PX2(:,1:nn)\transpose(sig);
        DATA=pG(ii,1:(end-delay+1));
        for j=2:delay
            DATA=[DATA;pG(ii,j:(end-delay+j))];
        end

        y = zeros(delay,110);        y(:,1) =  real(PV3(ii-60,1:nn)*coeffs)';
        
        for jj=1:109
            coeffs = K2*coeffs; 
            y(:,jj+1) = real(PX3(ii-60,1:nn)*coeffs)';
        end

        if ct==40
            if ii ==64
                Y1 = y;
            end
        end
        
        Er(ct,2) = Er(ct,2) + sum((vecnorm((y-DATA(:,1:length(y)))')./vecnorm(DATA(:,1:length(y))')).^2);
    end

    
    K2 = PV1(:,1:nn)\PV2(:,1:nn);

    for ii=61:65
        coeffs=PV1(:,1:nn)\transpose(sig);
        DATA=pG(ii,1:(end-delay+1));
        for j=2:delay
            DATA=[DATA;pG(ii,j:(end-delay+j))];
        end
        
        y = zeros(delay,110);        y(:,1) =  real(PV3(ii-60,1:nn)*coeffs)';%sig(:,1);
        
        for jj=1:109
            coeffs = K2*coeffs; 
            y(:,jj+1) = real(PV3(ii-60,1:nn)*coeffs)';
        end

        if ct==40
            if ii ==64
                Y2 = y;
                Z = DATA;
            end
        end

        Er(ct,1) = Er(ct,1) + sum((vecnorm((y-DATA(:,1:length(y)))')./vecnorm(DATA(:,1:length(y))')).^2);
    end
% plot(t(1:length(y)),y(1,:),'linewidth',2)
end
%%
figure
semilogx(nvec,Er/(5*delay),'linewidth',2)
xlim([0,100])
xlabel('number of modes','interpreter','latex','fontsize',18)
ylabel('MISE','interpreter','latex','fontsize',18)
legend({'Residual Ordering','Principal Component Ordering'},'interpreter','latex','fontsize',14)
ax=gca; ax.FontSize=18;
grid on
exportgraphics(gcf,'LIP_MISE.pdf','ContentType','vector','BackgroundColor','none')

%%

figure
plot(Y2(1,:),'linewidth',2)
hold on
plot(Y1(1,:),'linewidth',2)
plot(DATA(1,:),'k--','linewidth',2)
xlim([1,110])
xlabel('time step','interpreter','latex','fontsize',18)
legend({'Residual Ordering','Principal Component Ordering','True Signal'},'interpreter','latex','fontsize',14)
ax=gca; ax.FontSize=18;
exportgraphics(gcf,'LIP_signal.pdf','ContentType','vector','BackgroundColor','none')

%%

x_pts=-1.2:0.02:1.5;    y_pts=-0.02:0.02:1.5;
z_pts=kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:));    z_pts=z_pts(:);		% complex points where we compute pseudospectra
RES = KoopPseudoSpec(G,K,L,z_pts,'Parallel','off');	% compute pseudospectra
RES = reshape(RES,length(y_pts),length(x_pts));

%% Plot pseudospectra
figure
hold on
v=(10.^(-10:0.2:0));
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RES)),log10(v));
hold on
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),-reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RES)),log10(v));
cbh=colorbar;
cbh.Ticks=log10(10.^(-10:1:0));
cbh.TickLabels=10.^(-10:1:0);
clim([-4,0]);
reset(gcf)
set(gca,'YDir','normal')
colormap gray
axis equal;

title('Pseudospectrum','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(z)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(z)$','interpreter','latex','fontsize',18)

ax=gca; ax.FontSize=18; axis equal tight;   axis([x_pts(1),x_pts(end),-y_pts(end),y_pts(end)])
hold on
plot(real(LAM),imag(LAM),'.r','markersize',12);
box on
exportgraphics(gcf,'LIP_pseudospectra.pdf','ContentType','vector','BackgroundColor','none')
