function H = herm_eval_multidim(X,k,mu,sig)

% H = HERM_EVAL_MULTIDIM(X,k,mu,sig)
%
% evaluates the multidim Hermite polynomial order k (multi-index) orthonormal on [-inf,+inf]^N 
% with respect to rho=prod_i 1/sqrt(2 pi sigma_i) * e^( -(x-mi_i)^2/(2*sigma_i^2) ) 
% on the list of points X (each column is a point in R^N)
% MU, SIGMA, can be scalar values



%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


[N, nb_pts] = size(X); % N here is the number of dimensions

% take care of the fact that mu,sigma may be scalar values

if length(mu)==1
    mu=mu*ones(1,N);
    sig=sig*ones(1,N);
end

% H is a product of N polynomials, one for each dim. I store the evaluation
% of each of these polynomials as rows of a matrix. I do not compute the 
% zeros though, I know already they will be 1

nonzero_n = find(k~=0);
if isempty(nonzero_n)
    H = ones(1,nb_pts);
else
    M = zeros(length(nonzero_n),nb_pts);
    j=0;
    for n=nonzero_n
        j=j+1;
        M(j,:) = herm_eval(X(n,:),k(n),mu(n),sig(n));
    end
    H = prod(M,1);
end