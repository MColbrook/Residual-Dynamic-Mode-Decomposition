function [nu] = MomentMeas(MU,varargin)
% This code computes smoothed spectral measures of an isometry using the
% computed moments (Fourier coefficients). Requires chebfun.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% MU: vector of Fourier coefficients (-N to N)

% OPTIONAL LABELLED INPUTS
% filt: type of filter

% OUTPUTS
% nu: smoothed measure as a chebfun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addRequired(p,'MU',@isnumeric);
validFILT = {'fejer','cosine','vand','sharp_cosine'};
checkFILT = @(x) any(validatestring(x,validFILT));
addParameter(p,'filt','inft',checkFILT)

p.CaseSensitive = false;
parse(p,MU,varargin{:})

phi_fejer = @(x) 1-abs(x);
phi_cosine = @(x) (1+cos(pi*x))/2;
phi_opt4 = @(x)  1- x.^4.*(-20*abs(x).^3 + 70*x.^2 - 84*abs(x) + 35);
phi_sharp_cosine = @(x) phi_cosine(x).^4.*(35-84*phi_cosine(x)+70*phi_cosine(x).^2-20*phi_cosine(x).^3);
phi_inft = @(x) exp(-2./(1-abs(x)).*exp(-0.109550455106347./abs(x).^4));

N=(length(MU)-1)/2;

if p.Results.filt=="fejer"
    FILTER=phi_fejer(abs((-N:N)/N));
elseif p.Results.filt=="cosine"
    FILTER=phi_cosine(abs((-N:N)/N));
elseif p.Results.filt=="vand"
    FILTER=phi_opt4(abs((-N:N)/N));
elseif p.Results.filt=="sharp_cosine"
    FILTER=phi_sharp_cosine(abs((-N:N)/N));
else
    FILTER=phi_inft(abs((-N:N)/N));
end
FILTER(1)=0;
FILTER(end)=0;

nu = chebfun(FILTER(:).*MU(:),[-pi pi],'trig','coeffs');
end
