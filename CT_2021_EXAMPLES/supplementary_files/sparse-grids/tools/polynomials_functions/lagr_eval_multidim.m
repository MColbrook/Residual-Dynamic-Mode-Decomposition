function L= lagr_eval_multidim(current_knot,knots_per_dim,non_grid_points)

% L= LAGR_EVAL_MULTIDIM(current_knot,knots_per_dim,non_grid_points)
%
% evaluates the multidim Lagrange function L centered in current_knot  in non_grid_points.
%
% -> non_grid_points is a matrix; each row is a point to be evaluated
%
% -> knots_per_dim is a cell_array, that contains one array per cell. In the arrays are stored all the
% knots in each direction. E.g.
% knots_per_dim={[0 0.5 1],
%                [0 0.2 0.4 0.6 0.8 1] }
% it has to contain current_knot


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



% we are working on P points in N dimensions,
[P,N]=size(non_grid_points);

% the lagrangian functions in N dim are products of lagrangian functions in 1 D , so I work dimension-wise:
% we evaluate the multidimensional lagrange function on the corresponding column of points non_grid_points
% one dimension at a time and then we multiply
lagrange_monodim_eval=zeros(P,N);

for dim=1:N
    % these are the knots where the lagrange function is zero.
    % I compute it discarding the knot where I'm centering the multidim. lagr. function
    % this is an expensive way, with setdiff:
    % monodim_knots=setdiff(knots_per_dim{dim},current_knot(dim));
    % this is much cheaper, with logical indexing
    monodim_knots_temp=knots_per_dim{dim};
    monodim_knots=monodim_knots_temp(monodim_knots_temp~=current_knot(dim));
    
    % now i compute the lagrange function in the dim-th dimension and evaluate it in the points
    lagrange_monodim_eval(:,dim)=lagr_eval(current_knot(dim),monodim_knots,non_grid_points(:,dim));
end

L=prod(lagrange_monodim_eval,2);



%     switch type
%         case 'lagr'
%             monodim_knots_temp=knots_per_dim{dim};
%             monodim_knots=monodim_knots_temp(monodim_knots_temp~=current_knot(dim));
%             % now i compute the lagrange function in the dim-th dimension and evaluate it in the points
%             lagrange_monodim_eval(:,dim)=lagr_eval(current_knot(dim),monodim_knots,non_grid_points(:,dim));
%         case 'newton'
%             versor=(knots_per_dim{dim}==current_knot(dim));
%             lagrange_monodim_eval(:,dim)=lagr_eval_newton(knots_per_dim{dim},versor,non_grid_points(:,dim)');%lagr_eval(current_knot(dim),monodim_knots,non_grid_points(:,dim));
%     end