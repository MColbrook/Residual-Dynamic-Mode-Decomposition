clear
% close all

SCALE=1/10;

%% Set up the simulation
N=16;
dim=3; % dimension of phase space
delta_t=0.05;
SIGMA=10;
BETA=8/3;
RHO=40;
asq=4/BETA-1;
dvec=1:20;
dvec=((2*dvec-1).^2+asq)/(1+asq);

% these correspond to various generalisations of the Lorenz model
if dim==3
    ODEFUN=@(t,y) [SIGMA*(y(2)-y(1));
                    y(1).*(RHO-y(3)*SCALE)-y(2);
                    y(1).*y(2)*SCALE-BETA*y(3)];
elseif dim==5
    ODEFUN=@(t,y) [SIGMA*(y(2)-y(1));
                    y(1).*(RHO-y(3)*SCALE)-y(2);
                    y(1).*y(2)*SCALE-BETA*y(3)-y(1).*y(4)*SCALE;
                    y(1).*y(3)*SCALE-2*y(1).*y(5)*SCALE-dvec(2)*y(4);
                    2*y(1).*y(4)*SCALE-4*BETA*y(5)];
elseif dim==6
    ODEFUN=@(t,y) [SIGMA*(y(2)-y(1));
                    -y(1).*y(3)*SCALE+y(4).*y(3)*SCALE-2*y(4).*y(6)*SCALE+RHO*y(1)-y(2);
                    y(1).*y(2)*SCALE-y(1).*y(5)*SCALE-y(4).*y(2)*SCALE-BETA*y(3);
                    -dvec(2)*SIGMA*y(4)+SIGMA/dvec(2)*y(5);
                    y(1).*y(3)*SCALE-2*y(1).*y(6)*SCALE+RHO*y(4)-dvec(2)*y(5);
                    2*y(1).*y(5)*SCALE+2*y(4).*y(2)*SCALE-4*BETA*y(6)];
elseif dim==8
    ODEFUN=@(t,y) [SIGMA*(y(2)-y(1));
                    -y(1).*y(3)*SCALE+y(4).*y(3)*SCALE-2*y(4).*y(6)*SCALE+RHO*y(1)-y(2);
                    y(1).*y(2)*SCALE-BETA*y(3)-y(1).*y(5)*SCALE-y(4).*y(2)*SCALE-y(4).*y(7)*SCALE;
                    -dvec(2)*SIGMA*y(4)+SIGMA/dvec(2)*y(5);
                    y(1).*y(3)*SCALE-2*y(1).*y(6)*SCALE-dvec(2)*y(5)+RHO*y(4)-3*y(4).*y(8)*SCALE;
                    2*y(1).*y(5)*SCALE-4*BETA*y(6)+2*y(4).*y(2)*SCALE-2*y(1).*y(7)*SCALE;
                    2*y(1).*y(6)*SCALE-3*y(1).*y(8)*SCALE+y(4).*y(3)*SCALE-dvec(3)*y(7);
                    3*y(1).*y(7)*SCALE+3*y(4).*y(5)*SCALE-9*BETA*y(8)];
                    
elseif dim==9  
    ODEFUN=@(t,y) [SIGMA*(y(2)-y(1));
                    -y(1).*y(3)*SCALE+RHO*y(1)-y(2)+y(4).*y(3)*SCALE-2*y(4).*y(6)*SCALE+2*y(7).*y(6)*SCALE-3*y(7).*y(9)*SCALE;
                    y(1).*y(2)*SCALE-BETA*y(3)-y(1).*y(5)*SCALE-y(4).*y(2)*SCALE-y(4).*y(8)*SCALE-y(7).*y(5)*SCALE;
                    -dvec(2)*SIGMA*y(4)+SIGMA/dvec(2)*y(5);
                    y(1).*y(3)*SCALE-2*y(1).*y(6)*SCALE-dvec(2)*y(5)+RHO*y(4)-3*y(4).*y(9)*SCALE+y(7).*y(3)*SCALE;
                    2*y(1).*y(5)*SCALE-4*BETA*y(6)+2*y(4).*y(2)*SCALE-2*y(1).*y(8)*SCALE-2*y(7).*y(2)*SCALE;
                    -dvec(3)*SIGMA*y(7)+SIGMA/dvec(3)*y(8);
                    2*y(1).*y(6)*SCALE-3*y(1).*y(9)*SCALE+y(4).*y(3)*SCALE-dvec(3)*y(8)+RHO*y(7);
                    3*y(1).*y(8)*SCALE+3*y(4).*y(5)*SCALE-9*BETA*y(9)+3*y(7).*y(2)*SCALE];   
elseif dim==11
    ODEFUN=@(t,y) [SIGMA*(y(2)-y(1));
                    -y(1).*y(3)*SCALE+RHO*y(1)-y(2)+y(4).*y(3)*SCALE-2*y(4).*y(6)*SCALE+2*y(7).*y(6)*SCALE-3*y(7).*y(9)*SCALE;
                    y(1).*y(2)*SCALE-BETA*y(3)-y(1).*y(5)*SCALE-y(4).*y(2)*SCALE-y(4).*y(8)*SCALE-y(7).*y(5)*SCALE-y(7).*y(10)*SCALE;
                    -dvec(2)*SIGMA*y(4)+SIGMA/dvec(2)*y(5);
                    y(1).*y(3)*SCALE-2*y(1).*y(6)*SCALE-dvec(2)*y(5)+RHO*y(4)-3*y(4).*y(9)*SCALE+y(7).*y(3)*SCALE-4*y(7).*y(11)*SCALE;
                    2*y(1).*y(5)*SCALE-4*BETA*y(6)+2*y(4).*y(2)*SCALE-2*y(1).*y(8)*SCALE-2*y(7).*y(2)*SCALE-2*y(4).*y(10)*SCALE;
                    -dvec(3)*SIGMA*y(7)+SIGMA/dvec(3)*y(8);
                    2*y(1).*y(6)*SCALE-3*y(1).*y(9)*SCALE+y(4).*y(3)*SCALE-dvec(3)*y(8)+RHO*y(7)-4*y(4).*y(11)*SCALE;
                    3*y(1).*y(8)*SCALE+3*y(4).*y(5)*SCALE-9*BETA*y(9)+3*y(7).*y(2)*SCALE-3*y(1).*y(10)*SCALE;
                    -4*y(1).*y(11)*SCALE+2*y(4).*y(6)*SCALE+3*y(1).*y(9)*SCALE+y(7).*y(3)*SCALE-dvec(4)*y(10);
                    4*y(1).*y(10)*SCALE+4*y(4).*y(8)*SCALE+4*y(7).*y(5)*SCALE-16*BETA*y(11)];
                    
end

options = odeset('RelTol',1e-13,'AbsTol',1e-14);

%% Find the quadrature/data points
warning('off','SparseGKit:uint16')
[lev2knots,idxset]=define_functions_for_rule('SM',dim);
knots=@(n) hermpts(n);
[S,~]=smolyak_grid(dim,log2(2*N)+4,knots,lev2knots,idxset);     % sparse grids

Sr=reduce_sparse_grid(S);
xQ=transpose(Sr.knots);
wQ=transpose(Sr.weights)/sqrt(pi^dim); % correct so that we have a probability measure (just convention)

M=length(xQ); % number of data points
DATA=zeros(M,dim); % matrix of DATA points

pf = parfor_progress(M);
pfcleanup = onCleanup(@() delete(pf));
parfor j=1:M
    Y0=transpose(xQ(j,:));
    [~,Y]=ode45(ODEFUN,[0.000001 delta_t 2*delta_t],Y0,options);
    DATA(j,:)=Y(2,:);
    parfor_progress(pf);
end

%% Construct cells containing Hermite evaluations
X_hermite=cell(dim,1);
parfor jj=1:dim
    X_hermite{jj}=zeros(M,N);
    X_hermite{jj}(:,1)=0*xQ(:,jj)+1;
    X_hermite{jj}(:,2)=sqrt(2).*xQ(:,jj);
    for j=3:N
        X_hermite{jj}(:,j)=sqrt(2)/sqrt(j-1)*X_hermite{jj}(:,j-1).*xQ(:,jj)-sqrt(j-2)/sqrt(j-1)*X_hermite{jj}(:,j-2);
    end
end

Y_hermite=cell(dim,1);
parfor jj=1:dim
    Y_hermite{jj}=zeros(M,N);
    Y_hermite{jj}(:,1)=0*DATA(:,jj)+1;
    Y_hermite{jj}(:,2)=sqrt(2).*DATA(:,jj);
    for j=3:N
        Y_hermite{jj}(:,j)=sqrt(2)/sqrt(j-1)*Y_hermite{jj}(:,j-1).*DATA(:,jj)-sqrt(j-2)/sqrt(j-1)*Y_hermite{jj}(:,j-2);
    end
end

%% Create list of hyperbolic cross indices
Index=cell(dim,1);
Index0=ones(N,1);
Index{1}=exp(transpose(1:N));
for kk=2:dim
    Index{kk}=kron(Index0,exp(transpose(1:N)));
    II=log(Index{kk});
    for jj=1:kk-1
        Index{jj}=kron(Index{jj},ones(N,1));
        II=II.*log(Index{jj});
    end
    I=find(II<N+1);
    for jj=1:kk
        ff=Index{jj};
        Index{jj}=ff(I);
    end
    Index0=ones(length(I),1);
end

for kk=1:dim
    Index{kk}=log(Index{kk});
end
Ntrunc=length(Index{1});

%% Form the ResDMD matrices
A_matrix=zeros(Ntrunc,Ntrunc);  L_matrix=A_matrix;  G_matrix=A_matrix;

pf = parfor_progress(Ntrunc);
pfcleanup = onCleanup(@() delete(pf));
parfor (ii=1:Ntrunc,10)
    Xprod1=X_hermite{1}(:,Index{1}(ii));    Yprod1=Y_hermite{1}(:,Index{1}(ii));
    for ll=2:dim
        Xprod1=Xprod1.*X_hermite{ll}(:,Index{ll}(ii));  Yprod1=Yprod1.*Y_hermite{ll}(:,Index{ll}(ii));
    end
    for jj=1:Ntrunc
        Xprod2=X_hermite{1}(:,Index{1}(jj));    Yprod2=Y_hermite{1}(:,Index{1}(jj));
        for ll=2:dim
            Xprod2=Xprod2.*X_hermite{ll}(:,Index{ll}(jj));  Yprod2=Yprod2.*Y_hermite{ll}(:,Index{ll}(jj));
        end
        G_matrix(ii,jj)=sum(wQ(:).*conj(Xprod1).*Xprod2);
        A_matrix(ii,jj)=sum(wQ(:).*conj(Xprod1).*Yprod2);
        L_matrix(ii,jj)=sum(wQ(:).*conj(Yprod1).*Yprod2);
    end
    parfor_progress(pf);
end

L_matrix=(L_matrix+L_matrix')/2;
G_matrix=(G_matrix+G_matrix')/2;

%% compute pseudospectra
x_pts=-2:0.05:8;    y_pts=-3:0.05:3;
z_pts=kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:));    z_pts=z_pts(:);
RES=KoopPseudoSpec(G_matrix,A_matrix,L_matrix,z_pts,'parallel','on');	% compute pseudospectra
RES=reshape(RES,length(y_pts),length(x_pts));
[V,D]=eig((A_matrix)); % EDMD eigenvalues
E=diag(D);
%%
figure
v=[0.01,0.05,0.1,1];
contour(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),real(RES),v,'k','linewidth',1.5,...
    'ShowText','on')
set(gca,'YDir','normal')
ax=gca; ax.FontSize=14; axis equal tight
hold on
plot(real(E),imag(E),'.m')
axis([x_pts(1),x_pts(end),y_pts(1),y_pts(end)]); ylim([-3,3])
