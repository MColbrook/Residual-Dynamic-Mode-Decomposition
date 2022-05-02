function [nu] = IsomMeas(G,A,L,f,THETA,epsilon,varargin)
% This code computes smoothed spectral measures of an isometry using the
% resDMD matrices. Currently it is coded for dense matrices. Future
% releases will support sparse matrices.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% G: Gram matrix
% A: 1st Galerkin matrix  i.e. <K psi_j,psi_i>
% L: 2nd Galerkin matrix  i.e. <K psi_j,K psi_i>
% f: vector in \mathbb{C}^N (discretization of the observable)
% THETA: vector of points in periodic interval [-pi,pi] where we compute smoothed measures
% epsilon: the smoothing parameter

% OPTIONAL LABELLED INPUTS
% parallel: parfor (on) or normal for (off) loop, default is "off"
% order: order of kernel used, default is 2

% OUTPUTS
% nu: smoothed measure at points THETA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Collect the optional inputs
p = inputParser;
addRequired(p,'G',@isnumeric);
addRequired(p,'A',@isnumeric);
addRequired(p,'L',@isnumeric);
addRequired(p,'f',@isnumeric);
addRequired(p,'THETA',@isnumeric);
addRequired(p,'epsilon',@isnumeric);

validLQ = {'on','off'};
checkLQ = @(x) any(validatestring(x,validLQ));
validPar = {'on','off'};
checkPar = @(x) any(validatestring(x,validPar));

addParameter(p,'least_squares','off',checkLQ)
addParameter(p,'Parallel','off',checkPar)
addParameter(p,'Order',2,@(x) x==floor(x))

p.CaseSensitive = false;
parse(p,G,A,L,f,THETA,epsilon,varargin{:})

%% compute the poles and residues
delta=(2*(1:p.Results.Order)/(p.Results.Order+1)-1)*1i+1;
[c,d] = unitary_kern(delta,epsilon);

%% perform computation
if norm(G-speye(size(G)))<10^(-14)*norm(G)
    [Q,S]=schur(A);
    Z=Q;
    T=speye(size(A));
else
    [S,T,Q,Z]=qz(A,G);
    Q=Q';
end
v1=T*(Z')*f;
v2=(T')*(Q')*f;
v3=(S')*(Q')*f;
% perform computation
nu=0*THETA;
pf = parfor_progress(length(THETA));
pfcleanup = onCleanup(@() delete(pf));
if p.Results.Parallel=="off"
    for k=1:length(THETA)
        for j=1:p.Results.Order
            lambda=exp(1i*THETA(k))*(1+epsilon*delta(j));
            Ij=(S-lambda*T)\v1;
            nu(k)=nu(k)-real(1/(2*pi)*(c(j)*conj(lambda)*(Ij'*v2)+d(j)*(v3'*Ij)));
        end
        parfor_progress(pf);
    end
else
    parfor k=1:length(THETA)
        for j=1:p.Results.Order
            lambda=exp(1i*THETA(k))*(1+epsilon*delta(j));
            Ij=(S-lambda*T)\v1;
            nu(k)=nu(k)-real(1/(2*pi)*(c(j)*conj(lambda)*(Ij'*v2)+d(j)*(v3'*Ij)));
        end
        parfor_progress(pf);
    end
end


end

function [c,d] = unitary_kern(Z,epsilon)
m=length(Z);
sigma=-conj(Z(:))./(1+epsilon*conj(Z(:)));
V1=zeros(m); V2=V1;
for i=1:m 
    V1(:,i)=(sigma(:)).^(i-1);
    V2(:,i)=(Z(:)).^(i-1);
end

rhs=eye(m,1);
c=transpose(V1)\rhs;
d=transpose(V2)\rhs;
end
