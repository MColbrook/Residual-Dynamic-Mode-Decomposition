function f_values = interpolate_on_sparse_grid(S,Sr,function_on_grid,non_grid_points,~) 

% INTERPOLATE_ON_SPARSE_GRID interpolates a function on a sparse grid, i.e. evaluates the sparse  
% grid polynomial approximation (surrogate model) on a generic point of the parameters space.
%   
% F_VALUES = INTERPOLATE_ON_SPARSE_GRID(S,SR,FUNCTION_ON_GRID,NON_GRID_POINTS) evaluates the
%       sparse grid approximation a vector-valued function F: R^N -> R^V based on the sparse grid S. 
%       SR is the reduced version of S. 
%       FUNCTION_ON_GRID is a matrix containing the evaluation of F on the points of SR. 
%       Its dimensions are: V x number_of_points_in_the_sparse_grid
%       NON_GRID_POINTS is the set of points where one wants to evaluate the sparse grid polynomial  
%       approximation. It is a matrix, each column is a different point (i.e. the same convention as 
%       points stored in the knots field of S and SR).
%       F_VALUES is a matrix containing the evaluation of the vector-valued function F 
%       in each of the non_grid_points. Its dimensions are: V x number_of_non_grid_point


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


if nargin == 5
    errmsg=['Too many input arguments. Note that starting from release 14.4 '...
        'interpolate_on_sparse_grid does not accept any longer INTERVAL_MAP as second input argument. '...
        'The new function call is '...
        'F_VALUES = INTERPOLATE_ON_SPARSE_GRID(S,SR,FUNCTION_ON_GRID,NON_GRID_POINTS). '...
        'To fix this, either regenerate your sparse grid using INTERVAL_MAP as input (see help SMOLYAK_GRID) '...
        'or modify the fields KNOTS of your sparse grid and reduced sparse grid applying INTERVAL_MAP '...
        '(e.g.  S(i).knots = interval_map(S(i).knots) ) '...
        'This message will disappear in future releases of SPARSE_GRID_MATLAB_KIT'];
    error('SparseGKit:WrongInput',errmsg)
end

% in previous versions of interpolate_on_sparse_grid we needed
% function_on_grid to row-oriented, i.e. the evaluation of F on each sparse
% grid point was a different row, which is does not follow the convention
% used in other functions, i.e. evaluations being column-vector. We fix
% this as of march 2014 and throw an error if the previous convention is
% detected
[N,nb_points_Sr]=size(Sr.knots);
if size(function_on_grid,2)~=nb_points_Sr
    % if this condition is true, then function_on_grid stores
    % evaluation as rows rather than as columns and we stop the
    % function
    error('SparseGKit:WrongInput',[ 'Incompatible sizes. Starting from release 14.4 function_on_grid needs to store evaluations on grid as columns, '...
            'i.e. the dimensions of function_on_grid needs to be V x number_of_points_in_the_sparse_grid. '...
            'Type help INTERPOLATE_ON_SPARSE_GRID for details. Similarly, the size of the output matrix F_VALUES has been modified to '...
            'V x number_of_non_grid_points. This message will disappear in future releases of sparse-grid-matlab-kit'])
end

% similarly, we have changed the dimensions of non_grid_points to follow
% the convention that points are stored as columns. Again, we throw an
% error if the other convention is detected (this also prevents accidental input errors)
if size(non_grid_points,1)~=N
    error('SparseGKit:WrongInput',[ 'Incompatible sizes. Starting from release 14.4 non_grid_points needs to store points as columns, '...
            'i.e. its dimensions must be N x number_of_queried_evaluations (that is, following the same convention as points'...
             'stored in the ''knots'' sparse grids field). '...
            'Type help INTERPOLATE_ON_SPARSE_GRID for details.'])
end


% the other sizes. Observe that, although we have changed the orientation of the inputs for consistency
% matlab processes faster matrices when working column-wise, so we transpose everything to gain efficiency
V = size(function_on_grid,1);
nb_pts   = size(non_grid_points,2);
f_values = zeros(nb_pts,V); 
function_on_grid = function_on_grid';
non_grid_points = non_grid_points';

nb_grids=length(S);

% I need to go back from each point of each grid to its corresponding value in function_on_grid (F comes from an evaluation over a reduced grid)
% I only have a global mapping from [S.knots] to Sr.knots, so I need a global counter that scrolls [S.knots]
global_knot_counter=1;

% loop over the grids
for i=1:nb_grids
    
    % some of the grids in S structure are empty, so I skip it
    if isempty(S(i).weights)
        continue
    end
    
    % this is the set of points where I build the tensor lagrange function
    knots=S(i).knots;
    
    % I will need the knots in each dimension separately, to collocate the lagrange function.
    % I compute them once for all. As the number of knots is different in each direction, I use a cell array
    
    % we had here 
    % dimtot=size(knots,1);
    % clearly dimtot==N, hence we sobstitute it everywhere
    
    % the following lines are also no longer needed    
    %  knots_per_dim=cell(1,N);
    %  for dim=1:N
    %      knots_per_dim{dim}=unique(knots(dim,:));
    %  end
    
    knots_per_dim=S(i).knots_per_dim;
    
    % I also have to take into account the coefficient of the sparse grid
    
    coeff=S(i).coeff;
    
    
    % we could just loop on the interpolation points of the current tensor grid, 
    % and for each one evaluate the corresponding lagrangian polynomial on
    % the set of non_grid_points, but this is not convenient. Actually doing this we recompute the same thing over and over! 
    % think e.g. of the lagrangian functions for the knot [x1 y1 z1] and
    % for the knot [x1 y1 z2]: the parts on x and y are the same!
    
    % Therefore we evaluate each monodim_lagr just once and save all these
    % evaluations in a cell array. We will combine them suitably afterward.
    %
    % Such cell array stores one matrix per direction, i.e. N matrices. 
    % In turn, the n-th matrix contains the evaluations of each monodimensional lagrange polynomial 
    % in direction n on all the n-th coordinates of the non_grid_points.
    % Its dimension will be therefore (number_of_points_to_evaluate) X (number_of_lagrangian_polynomials_in_the_nth_direction)
    % 
    %
    % This is actually all the information that we will need to combine to
    % get the final interpolant value.
    
    mono_lagr_eval=cell(1,N);
    
    % loop on directions
    for dim=1:N
        
        % this is how many grid points in the current direction for the
        % current tensor grid
        %K=length(knots_per_dim{dim});
        K=S(i).m(dim);
        
        % allocate space for evaluations. Since I will be accessing it one lagrangian polynomial at a time
        % i.e. one knot at a time, it's better to have all information for the same lagrange polynomial
        % on the same column, for speed purposes. Moreover, note that whenever K=1 (one point only in a direction)
        % then the lagrange polynomial is identically one
        mono_lagr_eval{dim}=ones(nb_pts,K);
        
        if K>1
            % loop on each node of the current dimension and evaluate the corresponding monodim lagr polynomial.
            % We will need an auxiliary vector to pick up the current knot (where the lagr pol is centered) and
            % the remaining knots (where the lagr pol is zero). Here we see that mono_lagr_eval it's written
            % one column at a time
            aux=1:K;
            for k=aux
                %mono_lagr_eval{dim}(:,k) = lagr_eval(knots_per_dim{dim}(k),knots_per_dim{dim}(aux~=k),non_grid_points(:,dim));
                mono_lagr_eval{dim}(:,k) = lagr_eval_fast(knots_per_dim{dim}(k),knots_per_dim{dim}(aux~=k),K-1,non_grid_points(:,dim),[nb_pts 1]);
            end
        end

    end
    
    % now put everything together. We have to take the tensor product of
    % each of the monodim lagr pol we have evaluated. That is, we have to
    % pick one column for each matrix in the cell array and dot-multiply them.
    % all the possible combinations have to be generated !
    %
    % once this is done, we have the evaluation of each multidim lagr
    % polynomial on the non_grid_points, which we will then multiply by the
    % corresponding nodal value and eventually sum everything up.

    % We start by generating the combination of column we need to take. We actually don't
    % need to generate them, but only to recover it from the matrix knots,
    % which already contains all the points of the grid, i.e. all the
    % combinations of 1D points!
    %
    % Given a matrix of points like
    %
    % knots=[a1 a1 b1 b1 a1 a1 b1 b1 ...
    %        a2 b2 a2 b2 .....
    %        
    % combi is
    %
    % combi=[1 1 2 2 1 1 2 2 ...
    %        1 2 1 2 ......
    %
    %
    % again, we exploit the fact that the minimum entry of combi is 1 and that for many directions
    % there is only one point, so if we init combi with ones we're good in many cases
  
    combi = ones(N,S(i).size);
    
    % the easiest way to recover combi from knots is to proceed one dimension at a time,
    % and mark with a different label (1,2,...K) all the equal points. We need of course as many labels
    % as the number of different points in each dir!
    
    for dim=1:N
        
        % this is how many points per direction
        % K=length(knots_per_dim{dim});
        K=S(i).m(dim);
        
        % we start from a row of zeroes and we place 1....K in the right
        % positions by summations (each element of the row will be written
        % only once!). Since 1 are already in place, we proceed to place 2 and higher, but only if needed
        if K>1
            for k=2:K
                % here we add to the row of "1" either 0 or (k-1) so we get k where needed
                combi(dim,:) = combi(dim,:)+ (k-1)*( knots(dim,:)==knots_per_dim{dim}(k) ); 
            end
        end
    end


    % Now we can do the dot-multiplications among rows, the
    % multiplication by nodal values and the final sum! We proceed one
    % knot at a time
    
    for kk=1:S(i).size
        
        % dot-multiply all the lagrangian functions according to the
        % combi represented by the current knot. The result F_LOC is a
        % column vector
        
        f_loc=mono_lagr_eval{1}(:,combi(1,kk));
        for dim=2:N
            f_loc=f_loc.*mono_lagr_eval{dim}(:,combi(dim,kk));
        end


        % recover F, the corresponding value for the interpolating function in function_on_grid, with the global counter
        position = Sr.n(global_knot_counter);
        F_value = function_on_grid(position,:);
        
        % add the contribution of this knot to the sparse interpolation.
        % Its contribution is a matrix, since I am evaluating it on a bunch
        % of non_grid_points (one per row of f_values, will be trasposed later), 
        % and my function to be evaluated is vector-valued
        % which gives a bunch of rows in columns in f_values.
        % 
        % we generate this matrix as outer product of f_loc (column vector with lagr pol values) and F_value
        % (row vector with function evaluations)
        
        f_values = f_values + coeff*f_loc*F_value; 
        
        
        % update global counter
        global_knot_counter=global_knot_counter+1;
        
    end
    
    
end % of for loop on tensor grid

% finally, transpose to comply with output orientation
f_values=f_values';



% % -------------------------------------------------
% % old code
% 
% % for each knot in knots I have to build the multidimensional lagrangian function and evaluate it
% 
% for idx_current_knot=1:S(i).size
%     
%     % recover F, the corresponding value for the interpolating function in function_on_grid, with the global counter
%     position = Sr.n(global_knot_counter);
%     F_value = function_on_grid(position);
%     
%     % this is the current knot where the lagrange function is centered
%     current_knot=knots(:,idx_current_knot);
%     
%     % compute the contribute of the current knot to the value of f in non_grid_points
%     f_values = f_values + coeff*F_value*lagr_eval_multidim(current_knot,knots_per_dim,non_grid_points);
%     
%     % update global counter
%     global_knot_counter=global_knot_counter+1;
% end
% 
% % old code
% % -------------------------------------------------

