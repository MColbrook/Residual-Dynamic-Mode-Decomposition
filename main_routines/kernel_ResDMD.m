function [PSI_x,PSI_y,PSI_y2] = kernel_ResDMD(Xa,Ya,Xb,Yb,varargin)
% This code applies kernelized ResDMD.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% Xa and Ya: data matrices used in kernel_EDMD to form dictionary (columns
% correspond to instances of the state variable)
% Xb, Yb: data matrices used in ResDMD

% OPTIONAL LABELLED INPUTS
% N: size of computed dictionary, default is number of data points for kernel EDMD
% kernel_f: kernel used, default is normalised Gaussian
% Parallel: parfor (on) or normal for (off) loop, default is "off"
% Sketch: sketching to approximate Gaussian kernel, default is "off"
% s: size of sketching
% cut_off: stability parameter for SVD, default is 10^(-12)
% Y2: additional data matrix for stochastic version

% OUTPUTS
% PSI matrices for ResDMD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Collect the optional inputs
p = inputParser;

M1=size(Xa,2);
M2=size(Xb,2);
d=mean(vecnorm(Xa-mean(Xa')'));

validPar = {'on','off'};
checkPar = @(x) any(validatestring(x,validPar));
validSketch = {'on','off'};
checkSketch = @(x) any(validatestring(x,validSketch));

addParameter(p,'N',size(Xa,2),@(x) x==floor(x))
addParameter(p,'kernel_f',@(x,y) exp(-vecnorm(x-y)/d))
addParameter(p,'Parallel','off',checkPar)
addParameter(p,'Sketch','off',checkSketch)
addParameter(p,'cut_off',10^(-12),@(x) x>=0)
addParameter(p,'s',max(ceil(5*sqrt(M1+M2)*log(M1+M2)),5000),@(x) x==floor(x))
addParameter(p,'Y2',[],@isnumeric)

p.CaseSensitive = false;
parse(p,varargin{:})

% Apply kernel EDMD

if p.Results.Sketch=="off"
    G1=zeros(M1,M1); A1=G1;
    kernel_f=p.Results.kernel_f;
    pf = parfor_progress(M1); pfcleanup = onCleanup(@() delete(pf));
    if p.Results.Parallel=="off"
        for i=1:M1
            G1(i,:)=kernel_f(Xa(:,i),Xa);
            A1(i,:)=kernel_f(Ya(:,i),Xa);
            parfor_progress(pf);
        end
    else
        parfor i=1:M1
            G1(i,:)=kernel_f(Xa(:,i),Xa);
            A1(i,:)=kernel_f(Ya(:,i),Xa);
            parfor_progress(pf);
        end
    end
else
    s=p.Results.s;
    Z=sqrt(2/d^2)*randn([size(Xa,1),s]);
    TH=2*pi*rand(s,1);
    
    psi_xa=sqrt(2/s)*cos(TH+Z'*Xa);
    psi_ya=sqrt(2/s)*cos(TH+Z'*Ya);

    G1=psi_xa'*psi_xa;
    A1=psi_ya'*psi_xa;
    G1=max(G1,zeros(size(G1)));
    A1=max(A1,zeros(size(A1)));
end

% Post processing

[U,D0] = eig(G1+norm(G1)*p.Results.cut_off*eye(size(G1)));
D0(D0<p.Results.cut_off)=0;
SIG=sqrt(D0); SIG_dag = 0*SIG; SIG_dag(SIG>0)=1./(SIG(SIG>0));

K_hat = SIG_dag*U'*A1*U*SIG_dag;
[U1,D1] = eig(K_hat);

I=find(abs(diag(D1))>p.Results.cut_off);
if length(I)>p.Results.N
    [~,I]=sort(abs(diag(D1)),'descend');
    N=p.Results.N;
else
    N=length(I);
end

[P,~,~]=svd(U1(:,I(1:N)),0);
P=U*SIG_dag*P;

% Compute matrices for ResDMD

if p.Results.Sketch=="off"
    PSI_x=zeros(M2,M1); PSI_y=zeros(M2,M1);
    if ~isempty(p.Results.Y2) % stochastic case
        PSI_y2=zeros(M2,M1);
    end
    kernel_f=p.Results.kernel_f;

    pf = parfor_progress(M2); pfcleanup = onCleanup(@() delete(pf));
    if p.Results.Parallel=="off"
        for ii=1:M2
            PSI_x(ii,:)=kernel_f(Xb(:,ii),Xa);
            PSI_y(ii,:)=kernel_f(Yb(:,ii),Xa);
            if ~isempty(p.Results.Y2) % stochastic case
                PSI_y2(ii,:)=kernel_f(p.Results.Y2(:,ii),Xa);
            end
            parfor_progress(pf);
        end
    else
        parfor ii=1:M2
            PSI_x(ii,:)=kernel_f(Xb(:,ii),Xa);
            PSI_y(ii,:)=kernel_f(Yb(:,ii),Xa);
            if ~isempty(p.Results.Y2) % stochastic case
                PSI_y2(ii,:)=kernel_f(p.Results.Y2(:,ii),Xa);
            end
            parfor_progress(pf);
        end
    end

    PSI_x=PSI_x*P;
    PSI_y=PSI_y*P;
    if ~isempty(p.Results.Y2) % stochastic case
        PSI_y2=PSI_y2*P;
    end
else
    psi_xb=sqrt(2/s)*cos(TH+Z'*Xb);
    psi_yb=sqrt(2/s)*cos(TH+Z'*Yb);
    
    PSI_x=psi_xb'*psi_xa;
    PSI_y=psi_yb'*psi_xa;
    
    PSI_x=max(PSI_x,zeros(size(PSI_x)))*P;
    PSI_y=max(PSI_y,zeros(size(PSI_y)))*P;
    
    if ~isempty(p.Results.Y2) % stochastic case
        psi_y2b=sqrt(2/s)*cos(TH+Z'*p.Results.Y2);
        PSI_y2=psi_y2b'*psi_xa;
        PSI_y2=max(PSI_y2,zeros(size(PSI_y2)))*P;
    end

end

if isempty(p.Results.Y2)
    PSI_y2=[];
end

end