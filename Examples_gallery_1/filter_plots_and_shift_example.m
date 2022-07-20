clear
close all
figure(10)
h = plot(1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11,1:10,1:10,1:10,2:11);
c = get(h,'Color'); close(10);  clc

%% Plot the filters
phi_fejer = @(x) 1-abs(x);
phi_cosine = @(x) (1+cos(pi*x))/2;
phi_opt4 = @(x)  1- x.^4.*(-20*abs(x).^3 + 70*x.^2 - 84*abs(x) + 35);
phi_sharp_cosine = @(x) phi_cosine(x).^4.*(35-84*phi_cosine(x)+70*phi_cosine(x).^2-20*phi_cosine(x).^3);
phi_us = @(x) exp(-2./(1-abs(x)).*exp(-0.109550455106347./abs(x).^4));

xp=-1:0.01:1;

figure
plot(xp,phi_fejer(xp),'linewidth',2,'color',c{1})
hold on
plot(xp,phi_cosine(xp),'linewidth',2,'color',c{2})
plot(xp,phi_opt4(xp),'linewidth',2,'color',c{3})
plot(xp,phi_sharp_cosine(xp),'linewidth',2,'color',c{4})
plot(xp,phi_us(xp),'linewidth',2,'color',c{5})
legend({'Fej\''er','Cosine (RC)','Vandeven 4th','Sharp. RC','Bump'},'interpreter','latex','fontsize',14,'location','south')
ax = gca; ax.FontSize = 14;
ylim([0,1])

%% real space plots
N=10;
thetap=-pi:0.001:pi;
X=zeros(2*N+1,length(thetap));
for j=1:2*N+1
    X(j,:)=exp(1i*(j-N-1)*thetap)/(2*pi);
end

figure
plot(thetap,phi_fejer((-N:N)/N)*X,'linewidth',2,'color',c{1})
hold on
plot(thetap,phi_cosine((-N:N)/N)*X,'linewidth',2,'color',c{2})
plot(thetap,phi_opt4((-N:N)/N)*X,'linewidth',2,'color',c{3})
plot(thetap,(phi_sharp_cosine((-N:N)/N)*X),'linewidth',2,'color',c{4})
plot(thetap,(phi_us((-N:N)/N)*X),'linewidth',2,'color',c{5})
ax = gca; ax.FontSize = 14;
xlim([-pi,pi])
legend({'Fej\''er','Cosine (RC)','Vandeven 4th','Sharp. RC','Bump'},'interpreter','latex','fontsize',14,'location','northeast')


%% Shift example
Nvec=10:1:1000;
f = chebfun(@(t) (abs(t)<1),[-pi pi],'splitting','on');
rho = f*conj(f)/(2*pi); % analytic form of density of spectral measure

phi = chebfun(@(t) cos(5*t)/(2+cos(t)),[-pi pi],'trig'); % test function for weak convergence

MU_an = chebfun(@(t) abs(f(t)).^2/(2*pi),[-pi pi],'trig','splitting','on');
I=-0.023176565271898;

Ew=zeros(length(Nvec),5);
X=0;
E=zeros(length(Nvec),5);

for j=1:length(Nvec)
    MU=(sin(-Nvec(j):Nvec(j))/(2*pi^2))./(-Nvec(j):Nvec(j));
    MU(Nvec(j)+1)=1/(2*pi^2);
    
    MU_trunc = MomentMeas(MU,'filt','fejer');
    Ew(j,1)=abs(sum(MU_trunc*phi)-I)/abs(I);    % error of integration against test function
    E(j,1)=abs(MU_trunc(X)-rho(X))/rho(X);      % pointwise error of density
    
    MU_trunc = MomentMeas(MU,'filt','cosine');
    Ew(j,2)=abs(sum(MU_trunc*phi)-I)/abs(I);
    E(j,2)=abs(MU_trunc(X)-rho(X))/rho(X);
    
    MU_trunc = MomentMeas(MU,'filt','vand');
    Ew(j,3)=abs(sum(MU_trunc*phi)-I)/abs(I);
    E(j,3)=abs(MU_trunc(X)-rho(X))/rho(X);
    
    MU_trunc = MomentMeas(MU);
    Ew(j,5)=abs(sum(MU_trunc*phi)-I)/abs(I);
    E(j,5)=abs(MU_trunc(X)-rho(X))/rho(X);
end


%% Pointwise error
figure
loglog(Nvec,E(:,1),'linewidth',2,'color',c{1})
hold on
loglog(Nvec,E(:,2),'linewidth',2,'color',c{2})
ct=3;
for j=[3,5]
    loglog(Nvec,E(:,j),'linewidth',2,'color',c{ct});
    ct=ct+1;
end
ax = gca; ax.FontSize = 14;
ylim([10^(-15),1])

%% Weak convergence error
figure
loglog(Nvec,Ew(:,1),'linewidth',2,'color',c{1})
hold on
loglog(Nvec,Ew(:,2),'linewidth',2,'color',c{2})
ct=3;
for j=[3,5]
    loglog(Nvec,Ew(:,j),'linewidth',2,'color',c{ct});
    ct=ct+1;
end
ax = gca; ax.FontSize = 14;
ylim([10^(-15),1])

