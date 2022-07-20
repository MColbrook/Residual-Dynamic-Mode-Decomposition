function map = get_interval_map(a,b,type)

% MAP defines functions to shift sparse grids.
% 
% map = GET_INTERVAL_MAP(a,b,'uniform') where a,b are two vectors, 
%       returns a function that takes column points from (-1,1)^N to the generic box
%
%       ( a(1) b(1) ) x ( a(2) b(2) ) x ( a(3) b(3) ) ...
%
% map = GET_INTERVAL_MAP(a,b,'gaussian') where a,b are two vectors,
%       returns a function that takes points from R^N chosen according to standard gaussian measure
%       to points in R^N chosen according to a product of gaussians with means in a and stdev in b
%
% in both cases the output function must be used as
%
% X = interval_map(T)
%
% where T is a matrix of points in (-1,1)^N, one point per column (as in the output of 
% reduce_sparse_grid or tensor_grid). Similarly X is a matrix of points in 
% ( a(1) b(1) ) x ( a(2) b(2) ) x ( a(3) b(3) ), one point per column


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


switch type
    
    case 'uniform'
    
        % first we build D, resize matrix: it modifies a column vector "c" with components in 0,1
        % into a column vector "d" with components in ( 0 |b(1)-a(1)| ) x ( 0 |b(2)-a(2)| ) x ...
        D=diag(b-a);
        
        % next a translation vector: shifts a column vector "d" with components in ( 0 |b(1)-a(1)| ) x ( 0 |b(2)-a(2)| ) x ...
        % to "e" with components  in ( a(1) b(1) ) x ( a(2) b(2) ) x ...
        
        [r,c]=size(a);
        if r>c
            % a is a column vector
            p=a;
        else
            % a is a row vector
            p=a';
        end
        
        % finally, define the function
        map = @(T) D*(T+1)/2+p*ones(1,size(T,2));
        
    case 'gaussian'

        % D is the stdev matrix
        D=diag(b);
        
        % shift vector, according to means
        [r,c]=size(a);
        if r>c
            % a is a column vector
            p=a;
        else
            % a is a row vector
            p=a';
        end
        
        % finally, define the function
        map = @(T) D*T+p*ones(1,size(T,2));

end