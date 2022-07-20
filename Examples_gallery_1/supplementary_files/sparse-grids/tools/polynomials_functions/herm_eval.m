function H = herm_eval(x,k,mi,sigma)

% H = herm_eval(x,k,mi,sigma)
%
% returns the values of the k-th Hermite polynomial orthoNORMAL in (-inf,+inf) w.r.t to rho=1/sqrt(2 pi sigma) * e^( -(x-mi)^2/(2*sigma^2) )  in the points x
% ( x can be a matrix as well)
%
% N.B. the polynomials start from k=0: L_0(x) = 1, L_1(x) = (x - mi)/sigma


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



% this function expresses H as a function of the standard Hermite "probabilistic" polynomial (i.e. orthoGONAL w.r.t. rho=1/sqrt(2 pi) * e^(-x^2/2) ),
% which are recursively calculated through the function standard_herm_eval, coded below in this .m file

% first compute the transformation of x (referred to N(mi,sigma^2)) to z, the standard gaussian

% sigma=sqrt(sigma2);


z = ( x - mi ) / sigma ;

% calculate the standard legendre polynomials in t
H = standard_herm_eval(z,k);
% modify L to take into account normalizations. 
if k>1
    H = H / sqrt ( factorial(k) );
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function H = standard_herm_eval(x,k)

% L = standard_herm_eval(x) 
%
% returns the values of the k-th standard Hermite "probabilistic" polynomial (i.e. orthoGONAL w.r.t. rho=1/sqrt(2 pi) * e^(-x^2/2) ), in the points x
% ( x can be a vector as well)
%
% N.B. the polynomials start from k=0: L_0(x) = 1, L_1(x) = x

% base steps

% read this as L(k-2)
H_2=ones(size(x));

% and this as L(k-1)
H_1=x;

if k==0
      H=H_2;
      return
elseif k==1
      H=H_1;
      return
else
      % recursive step
      for ric=2:k
            H = x .* H_1 - (ric-1) * H_2;
            H_2=H_1;
            H_1=H;
      end
      return
end


