function [RES,RES2,V2] = KoopPseudoSpec(G,A,L,z_pts,varargin)
% This code computes pseudospectrum of K (currently written for dense matrices).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% G: Gram matrix
% A: 1st Galerkin matrix  i.e. <K psi_j,psi_i>
% L: 2nd Galerkin matrix  i.e. <K psi_j,K psi_i>
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
% addRequired(p,'G',@isnumeric);
% addRequired(p,'A',@isnumeric);
% addRequired(p,'L',@isnumeric);
% addRequired(p,'z_pts',@isnumeric);

validPar = {'on','off'};
checkPar = @(x) any(validatestring(x,validPar));

addParameter(p,'Parallel','off',checkPar)
addParameter(p,'z_pts2',[],@isnumeric)
addParameter(p,'reg_param',10^(-14),@isnumeric)

p.CaseSensitive = false;
parse(p,varargin{:})

%% compute the pseudospectrum
G=(G+G')/2; L=(L+L')/2; % safeguards
[VG,DG]=eig(G+norm(G)*(p.Results.reg_param)*eye(size(G)));
DG(abs(DG)>0)=sqrt(1./abs(DG(abs(DG)>0)));
SQ=VG*DG*(VG'); % needed to compute pseudospectra according to Gram matrix G

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
            RES(jj)=sqrt(real(eigs( SQ*((L)-z_pts(jj)*A'-conj(z_pts(jj))*A+(abs(z_pts(jj))^2)*G)*SQ,1,'smallestabs')));
            parfor_progress(pf);
        end
    else
        for jj=1:LL
            RES(jj)=sqrt(real(eigs( SQ*((L)-z_pts(jj)*A'-conj(z_pts(jj))*A+(abs(z_pts(jj))^2)*G)*SQ,1,'smallestabs')));
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
            [V,D]=eigs( SQ*((L)-p.Results.z_pts2(jj)*A'-conj(p.Results.z_pts2(jj))*A+(abs(p.Results.z_pts2(jj))^2)*G)*SQ,1,'smallestabs');
            V2(:,jj)=V; RES2(jj)=sqrt(real(D(1,1)));
            parfor_progress(pf);
        end
    else
        for jj=1:length(p.Results.z_pts2)
            [V,D]=eigs( SQ*((L)-p.Results.z_pts2(jj)*A'-conj(p.Results.z_pts2(jj))*A+(abs(p.Results.z_pts2(jj))^2)*G)*SQ,1,'smallestabs');
            V2(:,jj)=V; RES2(jj)=sqrt(real(D(1,1)));
            parfor_progress(pf);
        end
    end
end

warning('on','all')



end
