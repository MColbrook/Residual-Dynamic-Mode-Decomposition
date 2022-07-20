function L = lagr_eval_fast(current_knot,other_knots,ok_len,non_grid_points,ng_size)

% L = LAGR_EVAL_FAST(current_knot,other_knots,ok_len,non_grid_points,ng_size)
%
% builds the lagrange function L(x) s.t.
% L(current_knot)=1;
% L(other_knots)=0;
%
% and returns L=L(non_grid_points)
%
% where ok_len = length(other_knots), ng_size = size(non_grid_points);
% 
% this is essentially the same function as LAGR_EVAL, but some quantities are provided as input
% instead of being computed, for speed purposes. 

%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



% each monodim lagrange function is a product of K terms like (x-x_k)/(current_knot-x_k),
% so i compute it with a for loop

% this is the number of iterations
% ok_len=length(other_knots);

% L is the result. It is a column vector, same size as non_grid_points. It is initialized to 1, and then iteratively multiplied by
% each factor
% L=ones(size(non_grid_points));
L=ones(ng_size);

for k=1:ok_len
    % these are the non current-knot, one at a time
    knots_k=other_knots(k);
    % here it comes the computation of the lagrangian factor (x-x_k)/(current_knot-x_k)
    L=L.*(non_grid_points-knots_k)/(current_knot-knots_k);
end

end