clear
close all

%% 1D example for the point spectrum

TRUNC=1000;
omega = pi*1.4;
% omega = pi*sqrt(2);


f = chebfun(@(x)  cos(11*x)/(sqrt(1.001+sin(x))) ,[-pi,pi],'trig');
coeffs = trigcoeffs(f,2*TRUNC+1);

coeffs(abs(coeffs)<10^(-14))=0;

% semilogy(abs(coeffs))
% return
% Compute moments


N=100;
MU = zeros(2*N+1,1)';
Fc=coeffs;
MU(N+1)=(coeffs'*coeffs);
for j=1:N
    K = exp((-TRUNC:TRUNC)'*1i*omega*j);
    Fc=K.*coeffs;
    MU(N+1+j)=(Fc'*coeffs);
    MU(N-j+1)=(MU(N+1+j))';
end
mu1=MomentMeas(MU,'filt','fejer');  % approximation with filter
mu1=real(mu1)/sum(real(mu1));
K0=MomentMeas(MU*0+1,'filt','fejer');   K0=K0(0);
mu1 = mu1/real(K0);

N=10000;
MU = zeros(2*N+1,1)';
Fc=coeffs;
MU(N+1)=(coeffs'*coeffs);
for j=1:N
    K = exp((-TRUNC:TRUNC)'*1i*omega*j);
    Fc=K.*coeffs;
    MU(N+1+j)=(Fc'*coeffs);
    MU(N-j+1)=(MU(N+1+j))';
end
mu2=MomentMeas(MU,'filt','fejer');  % approximation with filter
mu2=real(mu2)/sum(real(mu2));
K0=MomentMeas(MU*0+1,'filt','fejer');   K0=K0(0);
mu2 = mu2/real(K0);



figure
subplot(2,1,1)
semilogy(mu1,'linewidth',0.5)
ax = gca; ax.FontSize = 14;
ylim([10^(-10),1])
title('$N=100$','interpreter','latex','fontsize',14)

subplot(2,1,2)
semilogy(mu2,'linewidth',0.5)
ax = gca; ax.FontSize = 14;
ylim([10^(-10),1])
title('$N=10000$','interpreter','latex','fontsize',14)