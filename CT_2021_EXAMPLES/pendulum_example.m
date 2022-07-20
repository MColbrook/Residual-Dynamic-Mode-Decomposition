clear
close all

%%%%% UNCOMMENT THE FOLLOWING TO PERFORM COMPUTATIONS TO SIMULATE DATA %%%%%
% %%
% delta_t=0.5;
% ODEFUN=@(t,y) [y(2);-sin(y(1))];
% options = odeset('RelTol',1e-16,'AbsTol',1e-16);
% 
% %% set up the computational grid for integration
% M1=300;
% M2=M1;
% N1=50; N2=100;
% %%
% 
% 
% L=15; % cut off
% x1=linspace(-pi,pi,M1+1);
% x1=x1+(x1(2)-x1(1))/2;
% x1=x1(1:end-1);
% x2=linspace(-L,L,M2);
% [X1,X2] = meshgrid(x1(1:end-1),x2);
% X1=X1(:); X2=X2(:);
% M=length(X1); % number of data points
% 
% DATA=zeros(M,2);
% 
% pf = parfor_progress(M);
% pfcleanup = onCleanup(@() delete(pf));
% parfor j=1:M
%     Y0=[X1(j);X2(j)];
%     [~,Y]=ode45(ODEFUN,[0.000001 delta_t 2*delta_t],Y0,options);
%     DATA(j,:)=[Y(2,1),Y(2,2)];
%     parfor_progress(pf);
% end
% W=zeros(M,1)+(x1(2)-x1(1))*(x2(2)-x2(1));
% 
% %% construct the relevant matrices using the bases
% A=zeros((2*N1+1)*N2);
% 
% %compute Hermite function values
% XH=zeros(M,N2);
% XH(:,1)=exp(-0.5*X2(:).^2)/(pi^(1/4));
% XH(:,2)=exp(-0.5*X2(:).^2)/(pi^(1/4))*sqrt(2).*X2(:);
% for j=3:N2
%     XH(:,j)=sqrt(2)/sqrt(j-1)*XH(:,j-1).*X2(:)-sqrt(j-2)/sqrt(j-1)*XH(:,j-2);
% end
% 
% YH=zeros(M,N2);
% YH(:,1)=exp(-0.5*DATA(:,2).^2)/(pi^(1/4));
% YH(:,2)=exp(-0.5*DATA(:,2).^2)/(pi^(1/4))*sqrt(2).*DATA(:,2);
% for j=3:N2
%     YH(:,j)=sqrt(2)/sqrt(j-1)*YH(:,j-1).*DATA(:,2)-sqrt(j-2)/sqrt(j-1)*YH(:,j-2);
% end
% %%
% 
% whos A
% pf = parfor_progress(M);
% pfcleanup = onCleanup(@() delete(pf));
% parfor (j=1:M,15)
%     X=kron(exp(1i*(-N1:N1)*X1(j))/sqrt(2*pi),XH(j,:));
%     Y=kron(exp(1i*(-N1:N1)*DATA(j,1))/sqrt(2*pi),YH(j,:));
%     A=A+(X')*Y;% kron(X',Y) - slightly slower;
%     parfor_progress(pf);
% end
% A=A*W(1);

%% Load data from above computation - file available from the dropbox link
clear
load('pendulum_data.mat','A','N1','N2')

Id1=max(kron(abs(-N1:N1),ones(1,N2)),1); Id2=kron(ones(1,2*N1+1),1:N2);
N_trun=N2;
Id=find((abs(Id1).*abs(Id2)<N_trun+1)); A=A(Id,Id); % sort out indexing

%% Compute pseudospectra (recall that, for this example, basis is orthonormal and Koopman operator is unitary)
x_pts=-1.5:0.05:1.5;    y_pts=x_pts;
z_pts=kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:));    z_pts=z_pts(:);

RES=KoopPseudoSpec(eye(size(A)),A,eye(size(A)),z_pts,'parallel','on');	% compute pseudospectra
RES=reshape(RES,length(y_pts),length(x_pts));

%%
E=eig((A)); % EDMD eigenvalues

%%
figure
hold on
v=[0.00000000000000000000000000000000000000000000000000000000000000000001,0.25];
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(1./real(RES)),log10(1./v),'k','linewidth',1.5)
reset(gcf)
set(gca,'YDir','normal')
colormap bone

ax=gca; ax.FontSize=14; axis equal tight
axis([x_pts(1),x_pts(end),y_pts(1),y_pts(end)])
hold on
plot(real(E),imag(E),'.m')
plot(real(exp(1i*(0:0.001:2*pi))),imag(exp(1i*(0:0.001:2*pi))),'-r')
    
%% Compute approximate eigenvalues and eigenfunctions
E_test=[0.880347701848197 + 0.473226384533343i;
  0.558620606228836 + 0.826563250517901i;
  0.124674806328568 + 0.987209196182015i;
 -0.314745042094750 + 0.936136228237800i];

z=E_test(4);
z=z/abs(z);

xp1=linspace(-pi,pi,500);
xp1=xp1+(xp1(2)-xp1(1))/2;
xp1=xp1(1:end-1);
xp2=linspace(-4,4,length(xp1));
XpH=zeros(length(xp1),N2);
XpH(:,1)=exp(-0.5*xp2(:).^2)/(pi^(1/4));
XpH(:,2)=exp(-0.5*xp2(:).^2)/(pi^(1/4))*sqrt(2).*xp2(:);
for j=3:N2
    XpH(:,j)=sqrt(2)/sqrt(j-1)*XpH(:,j-1).*xp2(:)-sqrt(j-2)/sqrt(j-1)*XpH(:,j-2);
end
XpF=exp(1i*kron((-N1:N1),xp1(:)));

H=eye(size(A))-(z*A')-(z*A')'+(abs(z)^2)*eye(size(A));
[V,D]=eigs(H,1,'smallestabs');
sqrt(D)
V2=zeros(N2*(2*N1+1),1);
V2(Id)=V;
V=reshape(V2,N2,2*N1+1);
PSI=transpose(XpF*transpose(V)*(XpH'));

figure
PhasePlot(xp1+1i*xp2,PSI,'m');
axis on
ax=gca;
ax.FontSize=14;



