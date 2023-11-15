function [RES,RES2,V2] = KoopPseudoSpecQR(PX,PY,W,z_pts,varargin)
% This code computes pseudospectrum of K (currently written for dense matrices).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% PX: dictionary evaluated at snapshots
% PY: dictionary evaluated at snapshots one time step later
% W: vector of weights for the quadrature
% z_pts: vector of complex points where we want to compute pseudospectra

% OPTIONAL LABELLED INPUTS
% parallel: parfor (on) or normal for (off) loop, default is "off"
% z_pts2: vector of complex points where we want to compute
% pseudoeigenfunctions
% reg_param: regularisation parameter for G

% OUTPUTS
% RES: residual for shifts z_pts.
% RES2: residual for pseudoeigenfunctions corresponding to shifts z_pts2.
% RES2: pseudoeigenfunctions corresponding to shifts z_pts2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Collect the optional inputs
p = inputParser;
validPar = {'on','off'};
checkPar = @(x) any(validatestring(x,validPar));

addParameter(p,'Parallel','off',checkPar)
addParameter(p,'z_pts2',[],@isnumeric)
addParameter(p,'reg_param',10^(-14),@isnumeric)

p.CaseSensitive = false;
parse(p,varargin{:})

%% compute the pseudospectrum
W = W(:);
[Q,R] = qr(sqrt(W).*PX,"econ");

C1 = (sqrt(W).*PY)/R;
L = C1'*C1;
G = eye(size(PX,2));
A = Q'*C1;

z_pts=z_pts(:);
LL=length(z_pts);
RES=zeros(LL,1);

if LL>0
    warning('off','all')
    pf = parfor_progress(LL);
    pfcleanup = onCleanup(@() delete(pf));
    if p.Results.Parallel=="on"
        parfor jj=1:LL
            warning('off','all')
            RES(jj)=sqrt(real(eigs(L-z_pts(jj)*A'-conj(z_pts(jj))*A+abs(z_pts(jj))^2*G,1,'smallestabs')));
            parfor_progress(pf);
        end
    else
        for jj=1:LL
            RES(jj)=sqrt(real(eigs(L-z_pts(jj)*A'-conj(z_pts(jj))*A+abs(z_pts(jj))^2*G,1,'smallestabs')));
            parfor_progress(pf);
        end
    end
end

RES2=[];
V2=[];

if ~isempty(p.Results.z_pts2)
    RES2=zeros(length(p.Results.z_pts2),1);
    V2=zeros(size(G,1),length(p.Results.z_pts2));
    pf = parfor_progress(length(p.Results.z_pts2));
    pfcleanup = onCleanup(@() delete(pf));
    if p.Results.Parallel=="on"
        parfor jj=1:length(p.Results.z_pts2)
            warning('off','all')
            [V,D]=eigs( L-p.Results.z_pts2(jj)*A'-conj(p.Results.z_pts2(jj))*A+abs(p.Results.z_pts2(jj))^2*G,1,'smallestabs');
            V2(:,jj)=V; RES2(jj)=sqrt(real(D(1,1)));
            parfor_progress(pf);
        end
    else
        for jj=1:length(p.Results.z_pts2)
            [V,D]=eigs( L-p.Results.z_pts2(jj)*A'-conj(p.Results.z_pts2(jj))*A+abs(p.Results.z_pts2(jj))^2*G,1,'smallestabs');
            V2(:,jj)=V; RES2(jj)=sqrt(real(D(1,1)));
            parfor_progress(pf);
        end
    end
    V2=R\V2;
end

warning('on','all')



end