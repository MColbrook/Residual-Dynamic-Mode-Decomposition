function L = lege_eval_multidim(X,k,a,b)

% L = LEGE_EVAL_MULTIDIM(X,k,a,b)
%
% evaluates the multidim Legendre polynomial order k (multi-index) orthonormal on [a,b]^N 
% on the list of points X (each column is a point in R^N). a,b are scalar values
%
% a,b may differ in each direction, so that actually the domain is not [a,b]^N, but
% [a1,b1] x [a2,b2] x [a3,b3] x ... [aN,bN] 
% in this case, a,b are defined as vectors, each with N components


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


[N, nb_pts] = size(X); % N here is the number of dimensions

L = zeros(1,nb_pts);

% L is a product of N polynomials, one for each dim. I store the evaluation
% of each of these polynomials as columns of a matrix

M = 0*X;

% take care of the fact that a,b may be scalar values

if length(a)==1
    a=a*ones(1,N);
    b=b*ones(1,N);
end

for n=1:N
    M(n,:) = lege_eval(X(n,:),k(n),a(n),b(n));
end

L = prod(M,1);