function L = cheb_eval(x,k,a,b)

% L = cheb_eval(x,k,a,b) 
%
% returns the values of the k-th Chebyshev polynomial of the first kind
% on [a,b] (i.e. T_n(a)=(-1)^n, T_n(b)=1), evaluated at x (vector)


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


% first compute the transformation of x in (a,b) to t in (-1,1)
t = ( 2*x - a - b ) / (b - a ) ;


% base steps
% read this as L(k-2)
L_2=ones(size(t));

% and this as L(k-1)
L_1=t;

if k==0
      L=L_2;
      return
elseif k==1
      L=L_1;
      return
else
      % recursive step
      for ric=2:k
            L = 2 * t .* L_1 -  L_2;
            L_2=L_1;
            L_1=L;
      end
      return
end