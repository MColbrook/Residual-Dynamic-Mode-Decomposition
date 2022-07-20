function L = lagr_eval(current_knot,other_knots,non_grid_points)

% L = LAGR_EVAL(current_knot,other_knots,non_grid_points)
%
% builds the lagrange function L(x) s.t.
% L(current_knot)=1;
% L(other_knots)=0;
%
% and returns L=L(non_grid_points)


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



% each monodim lagrange function is a product of K terms like (x-x_k)/(current_knot-x_k),
% so i compute it with a for loop

% this is the number of iterations
K=length(other_knots);

% L is the result. It is a column vector, same size as non_grid_points. It is initialized to 1, and then iteratively multiplied by
% each factor
L=ones(size(non_grid_points));

for k=1:K
    % these are the non current-knot, one at a time
    knots_k=other_knots(k);
    % here it comes the computation of the lagrangian factor (x-x_k)/(current_knot-x_k)
    L=L.*(non_grid_points-knots_k)/(current_knot-knots_k);
end

end