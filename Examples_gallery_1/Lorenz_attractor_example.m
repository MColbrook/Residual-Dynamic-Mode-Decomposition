clear
close all
rng(1);

%% Set parameters
M=10^5;
delta_t=0.05;
SIGMA=10;   BETA=8/3;   RHO=28;
ODEFUN=@(t,y) [SIGMA*(y(2)-y(1));y(1).*(RHO-y(3))-y(2);y(1).*y(2)-BETA*y(3)];
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
N = 100; % number of delay embeddings
%% Produce the data
Y0=(rand(3,1)-0.5)*4;
[~,Y]=ode45(ODEFUN,[0.000001 (1:(10000+M+2+N))*delta_t],Y0,options);
Y=Y(10001:end,:); % sample after when on the attractor

%% Reference scalar-valued spectral measure
PX1=zeros(M,N+1);   PX1(:,1)=Y(1:M,1);
X1=Y(1:M,1);
X2=Y(1:M,2);
X3=Y(1:M,3);

PX2=zeros(M,N+1);   PX2(:,1)=Y(1:M,2);
PX3=zeros(M,N+1);   PX3(:,1)=Y(1:M,3);

for j=2:N+1
    PX1(:,j)=Y((1:M)+j,1);
    PX2(:,j)=Y((1:M)+j,2);
    PX3(:,j)=Y((1:M)+j,3);
end

% PX = PX3(:,1:N);
PX = [PX1(:,1:N),PX2(:,1:N),PX3(:,1:N)];
% PY = PX3(:,2:(N+1));
PY = [PX1(:,2:(N+1)),PX2(:,2:(N+1)),PX3(:,2:(N+1))]; clear PX1 PX2 PX3

%%

G = (PX'*PX)/M;
A = (PX'*PY)/M;
L = (PY'*PY)/M;

G = (G+L)/2;
L=G;


E=exp(1i*[0.008,0.019,0.05,0.072,0.099,0.31,0.404,0.499,0.78]);
[~,RES,V] = KoopPseudoSpec(G,A,L,[],'z_pts2',E);

%%

for j=1:length(E)
    C=PX(1:min(2*10^4,M),:)*V(:,j);
    
    figure1=figure;
    axes1 = axes('Parent',figure1);
    hold(axes1,'on');
    scatter3(X1(1:min(2*10^4,M)),X2(1:min(2*10^4,M)),X3(1:min(2*10^4,M)),3,real(C),'filled');
    colormap turbo; colorbar; axis equal; view(axes1,[44.1307171302158 20.3999998682605]);
    grid(axes1,'on'); axis(axes1,'tight'); hold(axes1,'off');
    set(axes1,'DataAspectRatio',[1 1 1]);
    xlabel('$X$','interpreter','latex','fontsize',14)
    ylabel('$Y$','interpreter','latex','fontsize',14)
    zlabel('$Z$','interpreter','latex','fontsize',14,'rotation',0)
    title(sprintf('$\\theta=$%d, $\\mathrm{res}=$%d',imag(log(E(j))),RES(j)),'interpreter','latex','fontsize',20);
    % saveas(gcf,sprintf('LM%d.png',j))
    % close all
end
