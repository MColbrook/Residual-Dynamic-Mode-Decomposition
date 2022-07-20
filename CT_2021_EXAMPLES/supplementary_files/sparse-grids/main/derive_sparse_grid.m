function grads = derive_sparse_grid(S,Sr,values_on_grid,domain,eval_points,h)

% DERIVE_SPARSE_GRID computes derivatives (gradients) of sparse grid by centered finite differences
% 
% 
% GRADS = DERIVE_SPARSE_GRID(S,SR,VALUES_ON_GRID,DOMAIN,EVAL_POINTS) computes the derivative of the sparse grid.
%        by finite differences (see below for default value of increment size).
%
%        S is a sparse grid struct, SR is the reduced version of S, VALUES_ON_GRID are the values of the interpolated
%        function on SR. 
%        DOMAIN is a 2xN matrix = [a1, a2, a3, ...; b1, b2, b3, ...] defining the lower and upper bound of the 
%        hyper-rectangle on which the sparse grid is defined. The finite differences increment size is chosen according to
%        to the length of each interval [an bn] as h_n = (b_n - a_n)/1E5 
%
%        EVAL_POINTS are the points where the derivative must be evaluated. It is a matrix with points stored as 
%        columns, following the convention of the package
%
%        The output GRADS contains the compute values of the derivatives (gradients). The gradients in each point are stored
%        as columns of GRADS, i.e., size(GRADS) = N x size(EVAL_POINTS,2)
%
%
% GRADS = DERIVE_SPARSE_GRID(S,SR,VALUES_ON_GRID,DOMAIN,EVAL_POINTS,H) uses the input H as finite differences 
%        increment. H can be a scalar or a vector, in which case the n-th entry will be used as increment to approximate
%        the n-th component of the gradient


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



% get dimensions
N=size(domain,2);   % the sparse grid is defined over an N-variate hypercube
M=size(eval_points,2); % number of points where we need to evaluate the derivatives

grads = zeros(N,M);

switch nargin
    case 5
        h = zeros(1,N);
        a = domain(1,:);
        b = domain(2,:);
        
        for n=1:N
            h(n) = ( b(n)-a(n) ) / 1E5;
        end
    case 6
        if length(h) ==1
            h = h*ones(1,N);
        end
end

for k=1:N
    epsi = zeros(N,M);
    epsi(k,:) = h(k)*ones(1,M);
    grads(k,:) = ( interpolate_on_sparse_grid(S,Sr,values_on_grid,eval_points+epsi) - interpolate_on_sparse_grid(S,Sr,values_on_grid,eval_points-epsi) ) / (2*h(k));
end


