clear
close all
%% Load the data from the dropbox link
load('LIP.mat')
pG=transpose(pG);
pG=pG(:,40:150);

%% Algorithmic parameters
N=200;
delay=10;

%% Extract the data matrices using delay embedding
X1=[]; Y1=[];  
for ii=1:10
    DATA=pG(ii,1:(end-delay+1));
    for j=2:delay
        DATA=[DATA;pG(ii,j:(end-delay+j))];
    end
    X1=[X1,DATA(:,1:end-1)];    Y1=[Y1,DATA(:,2:end)];
end

X2=[]; Y2=[];  
for ii=11:60
    DATA=pG(ii,1:(end-delay+1));
    for j=2:delay
        DATA=[DATA;pG(ii,j:(end-delay+j))];
    end
    X2=[X2,DATA(:,1:end-1)];    Y2=[Y2,DATA(:,2:end)];
end

X3={}; Y3={}; ct=1; 
for ii=61:65
    DATA=pG(ii,1:(end-delay+1));
    for j=2:delay
        DATA=[DATA;pG(ii,j:(end-delay+j))];
    end
    X3{ct}=DATA(:,1:end-1);    Y3{ct}=DATA(:,2:end);
    ct=ct+1;
end

%% Apply kernel ResDMD
M1=size(X1,2);
[~,~,~,~,~,PSI_x,PSI_y] = kernel_ResDMD(X1,Y1,'Xb',X2,'Yb',Y2,'N',N,'Parallel','on','type','Lorentzian');
[~,S2,V2]=svd(transpose(X1)/sqrt(M1));
if delay>1
    PSI_x_DMD=transpose(X2)*V2*diag(1./diag(S2));
    PSI_y_DMD=transpose(Y2)*V2*diag(1./diag(S2));
else
    PSI_x_DMD=transpose(X2);
    PSI_y_DMD=transpose(Y2);
end
if delay>1
    PSI_x=[PSI_x,PSI_x_DMD];
    PSI_y=[PSI_y,PSI_y_DMD];
end

%%
M2=size(X2,2);
G_matrix=(PSI_x'*PSI_x)/M2; A_matrix=(PSI_x'*PSI_y)/M2; L_matrix=(PSI_y'*PSI_y)/M2;
G_matrix_DMD=(PSI_x_DMD'*PSI_x_DMD)/M2; A_matrix_DMD=(PSI_x_DMD'*PSI_y_DMD)/M2; L_matrix_DMD=(PSI_y_DMD'*PSI_y_DMD)/M2;

%% Compute pseudospectra
x_pts=-1.5:0.02:1.5;    y_pts=x_pts;
z_pts=kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:));    z_pts=z_pts(:);
RES = KoopPseudoSpec(G_matrix,A_matrix,L_matrix,z_pts,'Parallel','on'); RES=reshape(RES,length(y_pts),length(x_pts));

[V,D]=eig((A_matrix),(G_matrix)); E=diag(D); % EDMD eigenvalues
[V_DMD,D_DMD]=eig((A_matrix_DMD),(G_matrix_DMD)); E_DMD=diag(D_DMD); % DMD eigenvalues

%% Plot results
figure
hold on
v=(10.^(-5:0.25:0));
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RES)),log10(v));
cbh=colorbar;
cbh.Ticks=log10([0.0005,0.001,0.01,0.1,1]);
cbh.TickLabels=[0,0.001,0.01,0.1,1];
caxis([log10(0.001),log10(1)]);
reset(gcf)
set(gca,'YDir','normal')
c = turbo;
c = flipud(c);
colormap(c);
colormap bone

ax=gca; ax.FontSize=14; axis equal tight;   axis([x_pts(1),x_pts(end),y_pts(1),y_pts(end)])
hold on
plot(real(E),imag(E),'.r');


%% Test ResDMD compression
RES2 = KoopPseudoSpec(G_matrix,A_matrix,L_matrix,E,'Parallel','on');

% linear vs non-linear
[PSI_x2,PSI_y2] = kernel_ResDMD(X1,Y1,X3{1},Y3{1},'N',N,'Parallel','on');

if delay==1
    PSI_x2_DMD=transpose(X3{1});
    PSI_y2_DMD=transpose(Y3{1});
else
    PSI_x2_DMD=transpose(X3{1})*V2*diag(1./diag(S2));
    PSI_y2_DMD=transpose(Y3{1})*V2*diag(1./diag(S2));
    PSI_x2=[PSI_x2,PSI_x2_DMD];
    PSI_y2=[PSI_y2,PSI_y2_DMD];
end

xQ2=transpose(X2(1,:));
coeff_mat=G_matrix\(PSI_x'*xQ2/M2);
MAT1=V\coeff_mat;
coeff_mat_DMD=G_matrix_DMD\(PSI_x_DMD'*xQ2/M2);
MAT1_DMD=V_DMD\coeff_mat_DMD;

Dt=diag(D);        Dt_DMD=diag(D_DMD);

y_test=zeros(1,(length(X3{1}(1,:))-1)+1);        y_testb=zeros(1,(length(X3{1}(1,:))-1)+1);
for dd=1:delay
    MAT2=PSI_x2(dd,:)*V;
    MAT2_DMD=PSI_x2_DMD(dd,:)*V_DMD;
    for n=delay:(length(X3{1}(1,:))-1)
        y_test(n+1)=y_test(n+1)+real((MAT2.*transpose(Dt.^(n+1-dd)))*MAT1);
        y_testb(n+1)=y_testb(n+1)+real((MAT2_DMD.*transpose(Dt_DMD.^(n+1-dd)))*MAT1_DMD);
    end
end
y_real=real(X3{1}(1,1:(length(X3{1}(1,:))-1)+1));
y_test=y_test/delay;
y_testb=y_testb/delay;

load('LIP_times.mat')
t=t-t(1);
figure
plot(t(1:end-delay),y_real,'b','linewidth',2)
hold on
y_test(delay)=y_real(delay);
y_testb(delay)=y_real(delay);
plot(t(delay:length(y_real)),y_testb(delay:end),'or','linewidth',1)
plot(t(delay:length(y_real)),y_test(delay:end),'og','linewidth',1)

plot(0*[delay,delay]+t(delay),[-100,100],'--m')
legend({'true','linear dictionary','non-linear dictionary'},'interpreter','latex','fontsize',14,'location','best')
ax=gca; ax.FontSize=14;
xlim([0,t(end-delay)])

% modulus vs residual ordering
[~,IP]=sort(RES2,'ascend');

PSI_x_res=PSI_x*V(:,IP(1:40));
PSI_x2_res=PSI_x2*V(:,IP(1:40));
xQ2=transpose(X2(1,:));
MAT1_res=pinv(PSI_x_res)*xQ2;           
Dt_res=Dt(IP(1:40));

y_residual=zeros(1,(length(X3{1}(1,:))-1)+1);
for dd=1:delay
    MAT2_res=PSI_x2_res(dd,:);
    for n=delay:(length(X3{1}(1,:))-1)
        y_residual(n+1)=y_residual(n+1)+real((MAT2_res.*transpose(Dt_res.^(n+1-dd)))*MAT1_res);
    end
end
y_residual=y_residual/delay;

[~,IP]=sort(abs(E),'descend');

PSI_x_res=PSI_x*V(:,IP(1:40));
PSI_x2_res=PSI_x2*V(:,IP(1:40));
xQ2=transpose(X2(1,:));
MAT1_res=pinv(PSI_x_res)*xQ2;           
Dt_res=Dt(IP(1:40));

y_energy=zeros(1,(length(X3{1}(1,:))-1)+1);
for dd=1:delay
    MAT2_res=PSI_x2_res(dd,:);
    for n=delay:(length(X3{1}(1,:))-1)
        y_energy(n+1)=y_energy(n+1)+real((MAT2_res.*transpose(Dt_res.^(n+1-dd)))*MAT1_res);
    end
end
y_energy=y_energy/delay;

figure
plot(t(1:end-delay),y_real,'b','linewidth',2)
%
hold on
y_energy(delay)=y_real(delay);
y_residual(delay)=y_real(delay);
plot(t(delay:length(y_real)),y_energy(delay:end),'or','linewidth',1)
plot(t(delay:length(y_real)),y_residual(delay:end),'ok','linewidth',1)

plot(0*[delay,delay]+t(delay),[-100,100],'--m')
legend({'true','modulus ordering','residual ordering'},'interpreter','latex','fontsize',14,'location','best')
ax=gca; ax.FontSize=14;
xlim([0,t(end-delay)])

