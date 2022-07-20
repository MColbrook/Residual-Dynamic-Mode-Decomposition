clear
close all
%%
figure(10)
h = plot(1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11);
c = get(h,'Color');
close(10)

%% Choose observable
g = @(x) abs(x-1/3)+(x>0.78)+sin(20*x); gg = chebfun(@(x)  g(x),[0,1],'splitting','on');
f = @(x) g(x)/sqrt(sum(gg*conj(gg))); ff = chebfun(@(x)  f(x),[0,1],'splitting','on');  % function we compute measure wrt
n=15; % trunction size to compute the inner products
%% Setup the koopman matrix wrt Haar system
K=sparse(2^n,2^n);
% set up the indexing
nId=0;    kId=0;
for j=1:n
    nId=[nId(:);    zeros(2^(j-1),1)+j];    kId=[kId(:);    transpose(0:(2^(j-1)-1))];
end
K(1,1)=1;
for j=2:2^(n-1)
    I1 = find((nId==(nId(j)+1)).*(kId==kId(j)));
    if ~isempty(I1)
        K(I1,j)=1/sqrt(2);
    end
    I1 = find((nId==(nId(j)+1)).*(kId==(2^(nId(j))-(kId(j)+1))));
    if ~isempty(I1)
        K(I1,j)=-1/sqrt(2);
    end
end

%% Compute wavelet coefficients
xpts=1/(2^(n+1)):2^(-n):1;
[coeffs,L]=wavedec(f(xpts),round(log2(length(xpts))),'db1');
coeffs=coeffs/sqrt(length(xpts));
at=coeffs(1);
coeffs(1)=0; % subtract off the atomic part, will add on later
coeffs=coeffs(:);

%% Compute the fourier coefficients
MU=zeros(2*(2^n)+1,1);  N=1000;
MU(N+1)=(coeffs'*coeffs)/(2*pi);

Fc=coeffs;
for j=1:N
    Fc=K*Fc;
    MU(N-j+1)=(coeffs'*Fc)/(2*pi);
    MU(N+1+j)=(MU(N-j+1))';
end

MU=MU+abs(at)^2/(2*pi); % add atomic part back on

%% Spectral measures plot
N1=100; N2=1000;

MU01=chebfun(full(MU((N+1-N1):(N+1+N1))),[-pi pi],'trig','coeffs'); % approximation with no filter
MU1=MomentMeas(MU((N+1-N1):(N+1+N1)));  % approximation with filter
MU02=chebfun(full(MU((N+1-N2):(N+1+N2))),[-pi pi],'trig','coeffs'); % approximation with no filter
MU2=MomentMeas(MU((N+1-N2):(N+1+N2)));  % approximation with filter

%% N=100 plot
figure
semilogy(max(real(MU01),0),'color',c{2},'linewidth',1)
hold on
semilogy(max(real(MU1),10^(-16)),'color',c{1},'linewidth',1); hold on
ylim([10^(-2),50]); ax = gca; ax.FontSize = 14;
legend({'no filter','with filter'},'interpreter','latex','fontsize',14,'location','northwest')
%% N=1000 plot
figure
semilogy(max(real(MU02),0),'color',c{2},'linewidth',1)
hold on
semilogy(max(real(MU2),10^(-16)),'color',c{1},'linewidth',1); hold on
ylim([10^(-2),50]); ax = gca; ax.FontSize = 14;
legend({'no filter','with filter'},'interpreter','latex','fontsize',14,'location','northeast')

a2 = axes();    a2.Position = [0.18 0.70 0.3 0.2]; % xlocation, ylocation, xsize, ysize
semilogy(max(real(MU02),0),'color',c{2},'linewidth',1)
hold on
semilogy(max(real(MU2),10^(-16)),'color',c{1},'linewidth',1); hold on
axis([-1,-0.8,10^(-5),1])

%% Convergence in number of data points
% g = chebfun(@(x) abs(x-1/3)+(x>0.78)+sin(20*x),[0,1],'splitting','on');
% g = chebfun(@(x) cos(pi*x-0.1),[0,1],'splitting','on');
g = chebfun(@(x) cos(100*pi*x),[0,1],'splitting','on');
f = @(x) g(x)/sqrt(sum(g*conj(g)));

Nvec=1:12;
LL=zeros(length(Nvec),10);
ct=1;

for n=Nvec

    K=sparse(2^n,2^n);
    % set up the indexing
    nId=0;    kId=0;
    for j=1:n
        nId=[nId(:);    zeros(2^(j-1),1)+j];    kId=[kId(:);    transpose(0:(2^(j-1)-1))];
    end

    K(1,1)=1;

    for j=2:2^(n-1)
        I1 = find((nId==(nId(j)+1)).*(kId==kId(j)));
        if ~isempty(I1)
            K(I1,j)=1/sqrt(2);
        end
        I1 = find((nId==(nId(j)+1)).*(kId==(2^(nId(j))-(kId(j)+1))));
        if ~isempty(I1)
            K(I1,j)=-1/sqrt(2);
        end
    end

    xpts=1/(2^(n+1)):2^(-n):1;
    [coeffs,L]=wavedec(f(xpts),round(log2(length(xpts))),'db1');
    coeffs=coeffs/sqrt(length(xpts));
    at=coeffs(1);
    coeffs(1)=0; % subtract off the atomic part
    coeffs=coeffs(:);

    MU=zeros(2*(2^n)+1,1);  N=2^n;
    MU(N+1)=(coeffs'*coeffs)/(2*pi);
    LL(ct,1)=(coeffs'*coeffs)/(2*pi)+abs(at)^2/(2*pi);
    Fc=coeffs;
    for j=1:N
        Fc=K*Fc;
        MU(N-j+1)=(coeffs'*Fc)/(2*pi);
        MU(N+1+j)=(MU(N-j+1))';
        if j<11
            LL(ct,j+1)=MU(N-j+1)+abs(at)^2/(2*pi);
        end
    end
    MU=MU+abs(at)^2/(2*pi); % add atomic part back on
    
    ct=ct+1;
end

%% Compute the fourier modes via ergodic sampling
N=10;
SVEC=2.^Nvec(1:end-2)*10;
MU_an=LL(end,1:N+1);

ct=1;
for M=SVEC
    MU2=zeros(1,N+1);
    x1=rand(1);
    FT = @(x) max((min(2*(2*x<1).*x+(2*x>=1)*2.*(1-x),1)),0);
    Ya=zeros(M,1);
    Ya(1)=f(x1);
    Yb=Ya;
    xa=x1;
    xb=x1;

    for j=1:M-1
        xa=FT(xa);
        xb=FT(xb)+rand(1,1)*0.0001; % this kick is to stop underflow
        Ya(j+1)=f(xa);
        Yb(j+1)=f(xb);
    end
    for j=0:N
        MU2(j+1)=(Ya(1:(M-j))'*Ya((j+1):M)/(M-j))/(2*pi);
    end
    E(ct,1)=max(abs(MU2-MU_an));
    for j=0:N
        MU2(j+1)=(Yb(1:(M-j))'*Yb((j+1):M)/(M-j))/(2*pi);
    end
    E(ct,2)=max(abs(MU2-MU_an));

    ct=ct+1;
end

%%
figure
loglog(2.^Nvec(1:4)*10,max((abs(LL(1:4,:)-repmat(LL(end,:),[4,1])))'),'o-','color',c{1},'linewidth',2)
hold on
loglog(SVEC,E(1:length(SVEC),2),'o-','color',c{2},'linewidth',2)
ax = gca; ax.FontSize = 14;
loglog(SVEC,SVEC.^(-0.5),'k:','linewidth',2)
ylim([10^(-16),1])
xlim([10,10^4*2]);%[max(min(2.^Nvec)/2,1.0001),max(2.^Nvec)])
    
