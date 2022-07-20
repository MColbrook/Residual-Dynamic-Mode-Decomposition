clear
close all

%%%%% UNCOMMENT THE FOLLOWING TO PERFORM COMPUTATIONS TO SIMULATE DATA %%%%%
% %% y is [theta1,theta2,p1,p2] , ml^2=6, g/l=1/3
% delta_t=1;
% ODEFUN=@(t,y) [(2*y(3)-3*y(4)*cos(y(1)-y(2)))/(16-9*cos(y(1)-y(2))^2);
%     (8*y(4)-3*y(3)*cos(y(1)-y(2)))/(16-9*cos(y(1)-y(2))^2);
%     -3*((2*y(3)-3*y(4)*cos(y(1)-y(2)))/(16-9*cos(y(1)-y(2))^2)*(8*y(4)-3*y(3)*cos(y(1)-y(2)))/(16-9*cos(y(1)-y(2))^2)*sin(y(1)-y(2))+sin(y(1)));
%     -3*(-(2*y(3)-3*y(4)*cos(y(1)-y(2)))/(16-9*cos(y(1)-y(2))^2)*(8*y(4)-3*y(3)*cos(y(1)-y(2)))/(16-9*cos(y(1)-y(2))^2)*sin(y(1)-y(2))+sin(y(2))/3)];
% options = odeset('RelTol',1e-13,'AbsTol',1e-14);
% 
% %% Set up the computational grid for integration
% M1=50;  M2=25;
% L=5; % cut off
% x1=linspace(-pi,pi,M1+1);
% x1=x1+(x1(2)-x1(1))/2;
% 
% theta_grid=x1(1:end-1);
% p_grid=linspace(-L,L,M2);
% 
% X1=kron(theta_grid(:),kron(ones(M1,1),kron(ones(M2,1),ones(M2,1))));
% X2=kron(ones(M1,1),kron(theta_grid(:),kron(ones(M2,1),ones(M2,1))));
% X3=kron(ones(M1,1),kron(ones(M1,1),kron(p_grid(:),ones(M2,1))));
% X4=kron(ones(M1,1),kron(ones(M1,1),kron(ones(M2,1),p_grid(:))));
% 
% M=length(X1); % number of data points
% 
% DATA=zeros(M,4);
% 
% pf = parfor_progress(M);
% pfcleanup = onCleanup(@() delete(pf));
% parfor j=1:M
%     Y0=[X1(j);X2(j);X3(j);X4(j)];
%     [~,Y]=ode45(ODEFUN,[0.000001 delta_t 2*delta_t],Y0,options);
%     DATA(j,:)=Y(2,:);
%     parfor_progress(pf);
% end
% W=zeros(M,1)+(theta_grid(2)-theta_grid(1))^2*(p_grid(2)-p_grid(1))^2;
% 
% %% Construct the relevant matrices using the bases
% N1=30; N2=30;
% 
% % Compute Hermite function values
% XH1=zeros(M,N2);
% XH1(:,1)=exp(-0.5*X3(:).^2)/(pi^(1/4));
% XH1(:,2)=exp(-0.5*X3(:).^2)/(pi^(1/4))*sqrt(2).*X3(:);
% for j=3:N2
%     XH1(:,j)=sqrt(2)/sqrt(j-1)*XH1(:,j-1).*X3(:)-sqrt(j-2)/sqrt(j-1)*XH1(:,j-2);
% end
% XH2=zeros(M,N2);
% XH2(:,1)=exp(-0.5*X4(:).^2)/(pi^(1/4));
% XH2(:,2)=exp(-0.5*X4(:).^2)/(pi^(1/4))*sqrt(2).*X4(:);
% for j=3:N2
%     XH2(:,j)=sqrt(2)/sqrt(j-1)*XH2(:,j-1).*X4(:)-sqrt(j-2)/sqrt(j-1)*XH2(:,j-2);
% end
% 
% YH1=zeros(M,N2);
% YH1(:,1)=exp(-0.5*DATA(:,3).^2)/(pi^(1/4));
% YH1(:,2)=exp(-0.5*DATA(:,3).^2)/(pi^(1/4))*sqrt(2).*DATA(:,3);
% for j=3:N2
%     YH1(:,j)=sqrt(2)/sqrt(j-1)*YH1(:,j-1).*DATA(:,3)-sqrt(j-2)/sqrt(j-1)*YH1(:,j-2);
% end
% YH2=zeros(M,N2);
% YH2(:,1)=exp(-0.5*DATA(:,4).^2)/(pi^(1/4));
% YH2(:,2)=exp(-0.5*DATA(:,4).^2)/(pi^(1/4))*sqrt(2).*DATA(:,4);
% for j=3:N2
%     YH2(:,j)=sqrt(2)/sqrt(j-1)*YH2(:,j-1).*DATA(:,4)-sqrt(j-2)/sqrt(j-1)*YH2(:,j-2);
% end
% 
% %%
% I1=max(abs(-N1:N1),1);
% I2=(1:N2);
% IDX=kron(I1,kron(I1,kron(I2,I2)));
% I=find(IDX<25);
% A=zeros(length(I));
% whos A
% 
% pf = parfor_progress(M);
% pfcleanup = onCleanup(@() delete(pf));
% parfor (j=1:M,40)
%     X=kron(exp(1i*(-N1:N1)*X1(j))/sqrt(2*pi),kron(exp(1i*(-N1:N1)*X2(j))/sqrt(2*pi),kron(XH1(j,:),XH2(j,:))));
%     Y=kron(exp(1i*(-N1:N1)*DATA(j,1))/sqrt(2*pi),kron(exp(1i*(-N1:N1)*DATA(j,2))/sqrt(2*pi),kron(YH1(j,:),YH2(j,:))));
%     X=X(I); Y=Y(I);
%     A=A+(X')*Y;
%     parfor_progress(pf);
% end
% A=A*W(1);

%% Load data from above computation - file available from the dropbox link
load('double_pendulum_data.mat')
I1=max(abs(-N1:N1),1);  I2=(1:N2);  IDX=kron(I1,kron(I1,kron(I2,I2)));  I=find(IDX<25); IDX=IDX(I); I2=find(IDX<N1);
A=A(I2,I2);     % sort out indexing

%% Compute spectral measures
theta=-pi:0.005:pi;                                     % points where we compute spectral measure
theta=[theta(abs(theta)>0.002),0]; theta=sort(theta);
epsilon=0.1;                                            % smoothing parameter

c1=zeros(1,2*N1+1); c2=c1; c3=zeros(1,N2); c4=c3;
CASE=1; % the observables in the paper
if CASE==1
    c1(N1+2)=1; c2(N1+1)=1; c3(1)=1; c4(1)=1;
elseif CASE ==2
    c1(N1+1)=1; c2(N1+2)=1; c3(1)=1; c4(1)=1;
elseif CASE==3
    c1(N1+1)=1; c2(N1+1)=1; c3(2)=1; c4(1)=1;
else
    c1(N1+1)=1; c2(N1+1)=1; c3(1)=1; c4(2)=1;
end
f=transpose(kron(c1,kron(c2,kron(c3,c4)))); f=f(I); f=f(I2); f=f/norm(f);   % coefficient vector for the test functions

nu1 = IsomMeas(eye(size(A)),A,eye(size(A)),f,theta,epsilon,'parallel','on','order',1);
nu6 = IsomMeas(eye(size(A)),A,eye(size(A)),f,theta,epsilon,'parallel','on','order',6);
%% Plot the results
figure
plot(theta,nu1,'linewidth',2)
hold on
plot(theta,nu6,'linewidth',2)
xlim([-pi,pi]); ax = gca; ax.FontSize=14;
legend({'$m=1$','$m=6$'},'interpreter','latex','fontsize',20,'location','northeast')
