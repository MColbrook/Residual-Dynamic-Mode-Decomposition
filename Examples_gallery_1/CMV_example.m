clear
close all
%%
sigma=0.001;                % standard deviation for noisy experiment. NB: results will be slightly different each time due to randomness
m=1;                        % order of kernel
Nvec=[1:200,1000];          % dector of truncation sizes for convergence plot
MU1=zeros(length(Nvec),1);  % vector for measure with epsilon=0.5
MU2=zeros(length(Nvec),1);  % vector for measure with epsilon=0.1
MU3=zeros(length(Nvec),1);  % vector for measure with epsilon=0.01
%% Compute approximations at theta=0.2 for different truncation sizes
for j=1:length(Nvec)
    MU1(j)=CMV_meas(Nvec(j),m,0.5,0.2,sigma);
    MU2(j)=CMV_meas(Nvec(j),m,0.1,0.2,sigma);
    MU3(j)=CMV_meas(Nvec(j),m,0.01,0.2,sigma);
end
%% Plot the convergence (as relative errors)
figure
semilogy(Nvec,abs(MU1-MU1(end))/abs(MU1(end)),'linewidth',2)
hold on
semilogy(Nvec,abs(MU2-MU2(end))/abs(MU2(end)),'linewidth',2)
semilogy(Nvec,abs(MU3-MU3(end))/abs(MU3(end)),'linewidth',2)
ax = gca; ax.FontSize = 14; ylim([10^(-4),1]); xlim([0,200]);
legend({'$\epsilon=0.5$','$\epsilon=0.1$','$\epsilon=0.01$'},'interpreter','latex','fontsize',14,'location','northeast')
%% Compute approximations of measures over the interval [-pi,pi]
theta=-pi:0.005:pi;
mu = CMV_meas(10000,6,0.01,theta,0);                % reference accurate solution
mu1 = CMV_meas(200,1,0.1,theta,sigma);              % smoothed measure with epsilon=0.1
mu2 = CMV_meas(500,1,0.01,-pi:0.001:pi,sigma);      % smoothed measure with epsilon=0.01
%% Plot the results
figure
h2=plot(theta,mu1,'linewidth',2);
xlim([-1,1]); ylim([0,2.5]);
hold on
plot(theta,mu1,'linewidth',2)
plot(-pi:0.001:pi,mu2,'linewidth',2)
plot(theta,mu,'k','linewidth',1)
ax = gca; ax.FontSize = 14;
legend({'$\epsilon=0.1$','$\epsilon=0.1$','$\epsilon=0.01$','$\rho_f$'},'interpreter','latex','fontsize',14,'location','northwest')
delete(h2)

function mu = CMV_meas(n,m,epsilon,theta,sigma)
% For truncation size n, order of kernel m, smoothing parameter epsilon,
% angles theta and noise level sigma, this approximates the spectral measure of the CMV example.
  
    % Set up the matrix.
    q=0.95;
    a_c = @(k) (-1).^k.*q.^((k+1)/2);
    rho_c = @(k) sqrt(1-abs(a_c(k)).^2);
    A=sparse(n+2,n+2);
    A(1:2,1:3)=[conj(a_c(0)) conj(a_c(1))*rho_c(0) rho_c(1)*rho_c(0);
        rho_c(0) -conj(a_c(1))*a_c(0) -rho_c(1)*a_c(0)];
    for j=1:round((n+2)/2)
        A([2*j+1,2*j+2],[2*j:2*j+3])=[conj(a_c(2*j))*rho_c(2*j-1) -conj(a_c(2*j))*a_c(2*j-1) conj(a_c(2*j+1))*rho_c(2*j) rho_c(2*j)*rho_c(2*j+1);
                                      rho_c(2*j)*rho_c(2*j-1)     -rho_c(2*j)*a_c(2*j-1)     -conj(a_c(2*j+1))*a_c(2*j)  -rho_c(2*j+1)*a_c(2*j)];
    end
    A=A(1:n,1:n);
    if n<1000
        A=A+sigma*randn(n,n);
    end
    f=zeros(n,1);
    f(1)=1;
    % Call the routine for computing spectral measures.
    mu = IsomMeas(speye(n),A,speye(n),f,theta,epsilon,'order',m);
end




