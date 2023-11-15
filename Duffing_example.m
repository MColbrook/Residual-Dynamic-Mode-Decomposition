clear
close all
% Simple example for the unforced duffing oscillator

rng(1)
%% Set parameters
M1=10^3; % number of data points
M2=50;
delta_t=0.25; % time step
ODEFUN=@(t,y) [y(2);y(1)-y(1).^3];
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
N=100;
PHI = @(r) exp(-r); % radial basis function

x_pts = -1.2:0.04:1.2;    y_pts = -0.05:0.04:1.2;
v=(10.^(-2:0.1:1));

%% Produce the data
X=[];
Y=[];
for jj=1:M1
    Y0=(rand(2,1)-0.5)*4;
    [~,Y1]=ode45(ODEFUN,[0 0.000001 (1:(3+M2))*delta_t],Y0,options);
    Y1=Y1';
    X = [X,Y1(:,[1,3:M2+1])];
    Y = [Y,Y1(:,3:M2+2)];
end
M = M1*M2;

d=mean(vecnorm(X-mean(X')')); % scaling for radial function

[~,C] = kmeans([X';Y'],N); % find centers

PX = zeros(M,N); PY = zeros(M,N);

for j = 1:N
    R = sqrt((X(1,:)-C(j,1)).^2+(X(2,:)-C(j,2)).^2);
    PX(:,j) = PHI(R(:)/d);
    R = sqrt((Y(1,:)-C(j,1)).^2+(Y(2,:)-C(j,2)).^2);
    PY(:,j) = PHI(R(:)/d);
end

%% Apply ResDMD algorithm 1 (residuals computed after EDMD)
K = PX\PY;
[V,LAM] = eig(K,'vector');
res = (vecnorm(PY*V-PX*V*diag(LAM))./vecnorm(PX*V))'; % residuals

%%
figure
scatter(real(LAM),imag(LAM),250,res,'.');
colormap turbo; colorbar;
ax=gca; ax.FontSize=14; axis equal tight;
xlim([-1,1]);
ylim([-1,1]);

%% ResDMD for pseudospectrum
z_pts=kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:));    z_pts=z_pts(:);

RES = KoopPseudoSpecQR(PX,PY,1/M,z_pts,'Parallel','off');
RES = reshape(RES,length(y_pts),length(x_pts));

%%
figure
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(max(min(v),real(RES))),log10(v));
hold on
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),-reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(max(min(v),real(RES))),log10(v));
cbh=colorbar;
cbh.Ticks=log10([0.005,0.01,0.1,1]);
cbh.TickLabels=[0,0.01,0.1,1];
clim([log10(min(v)),0]);
reset(gcf)
set(gca,'YDir','normal')
colormap gray
ax=gca; ax.FontSize=18; axis equal tight;   axis([x_pts(1),x_pts(end),-y_pts(end),y_pts(end)])
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18)
title(sprintf('$\\mathrm{Sp}_\\epsilon(\\mathcal{K})$, $N=%d$',N),'interpreter','latex','fontsize',18)











