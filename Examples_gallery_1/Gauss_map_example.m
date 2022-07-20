clear
close all
%% Compute quadrature errors
test=quad_error(200,1);

Mvec=unique(round(10.^(0.2:0.05:4)));   % vector of sample sizes to track convergence
Mvec(Mvec==200)=201;
E=zeros(length(Mvec),4);
ct=1;   % counter for forloop
for j=Mvec   
    for k=0:3
        A=quad_error(j,k);            
        E(ct,k+1)=max(abs(test(:)-A(:)));
    end
    ct=ct+1;
end
%% Plot quadrature errors
close all
figure
loglog(Mvec,E(:,2),'linewidth',2)
hold on
loglog(Mvec,E(:,3),'linewidth',2)
loglog(Mvec,E(:,4),'linewidth',2)
loglog(Mvec,E(:,1),'linewidth',2)
loglog(Mvec(50:end),20./Mvec(50:end),'k:','linewidth',2)
loglog(Mvec(50:end),5000./Mvec(50:end).^2,'k:','linewidth',2)
loglog(Mvec(50:end),10./Mvec(50:end).^0.5,'k:','linewidth',2)
legend({'Gauss--Legendre','Trapezoidal','Riemann sum','Monte Carlo'},'interpreter','latex','fontsize',14,'location','southeast')
ax = gca; ax.FontSize = 14;
xlim([10,10000])
ylim([10^(-15),10])

%% Construct the three ResDMD matrices
alpha=2; beta=-1-exp(-alpha);
N=40; M=100;
[X, W] = legpts(M,[-1,0]);
DATA=exp(-alpha*X(:).^2)+beta;

chebpoly(0);
XH=zeros(M,N); YH=XH;
for j=1:N
    f=legpoly(j-1,[-1,0]);
    XH(:,j)=f(X)/sqrt(2/(2*j-1))*sqrt(2);
    YH(:,j)=f(DATA)/sqrt(2/(2*j-1))*sqrt(2);
end
G=zeros(N); A=zeros(N); L=zeros(N);

pf = parfor_progress(M);
pfcleanup = onCleanup(@() delete(pf));
for j=1:M % 
    G=G+(XH(j,:).*W(j))'*XH(j,:);
    A=A+(XH(j,:).*W(j))'*YH(j,:);
    L=L+(YH(j,:).*W(j))'*YH(j,:);
    parfor_progress(pf);
end

%% Pseudospectra plot
x_pts=-1.5:0.05:1.5;    y_pts=x_pts;
z_pts=kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:));    z_pts=z_pts(:); % spectral parameter where we compute spectra

RES=KoopPseudoSpec(G,A,L,z_pts);	% compute pseudospectra
RES=reshape(RES,length(y_pts),length(x_pts));

E=eig(A); % EDMD eigenvalues

%%
figure
v=([0.001,0.01,0.1,0.3]);
contour(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),real(RES*0.99),v,'k',...
    'linewidth',1.5,'ShowText','on')
set(gca,'YDir','normal')
colormap bone
ax=gca; ax.FontSize=14;
axis equal tight;   axis([x_pts(1),x_pts(end),y_pts(1),y_pts(end)])

hold on
II1=find((abs(imag(E))<10^(-6))&(real(E)>-0.01));
C = setdiff( 1:length(E),II1 );
plot(real(E(C)),imag(E(C)),'.m')
plot(real(E(II1)),imag(E(II1)),'xb','MarkerSize',7)

function A = quad_error(M,type)
    % This function constructs the Koopman matrix for various quadrature
    % rules. type 1 is Legendre-quad, 2 is trapezoidal, 0 is random, other
    % is Riemann sum
    alpha=2;
    beta=-1-exp(-alpha);
    N=40;

    if type==1
        [X, W] = legpts(M,[-1,0]);
    elseif type==2
        X=linspace(-1,0,M);
        X=X(:); W=0*X+X(2)-X(1);
        W(1)=W(1)/2;
        W(end)=W(end)/2;
    else
        X=linspace(-1,0,M);
        X=X(:); W=0*X+X(2)-X(1);
    end
    if type==0
        X=-rand(M,1);
        X=X(:); W=0*X+1/M;
    end
    DATA=exp(-alpha*X(:).^2)+beta;

    XH=zeros(M,N); YH=XH;
    for j=1:N
        f=legpoly(j-1,[-1,0]);
        XH(:,j)=f(X)/sqrt(2/(2*j-1))*sqrt(2);
        YH(:,j)=f(DATA)/sqrt(2/(2*j-1))*sqrt(2);
    end
    A=zeros(N);
    for j=1:M
        A=A+(XH(j,:).*W(j))'*YH(j,:);
    end
end