function [modal_coeffs,K] = convert_to_modal(S,Sr,nodal_values,domain,flags,~)

% CONVERT_TO_MODAL recasts a sparse grid interpolant as a sum of orthogonal polynomials
% i.e. computes the spectral expansion of the interpolant.
%
% 
% [MODAL_COEFFS,K] = CONVERT_TO_MODAL(S,SR,NODAL_VALUES,DOMAIN,'legendre') returns the Legendre expansion
%       of the sparse grid interpolant. S is a sparse grid, SR is its reduced counterpart, NODAL_VALUES
%       are the values of the target function (say F) evaluated on the reduced sparse grid, and has the same
%       size of e.g. F_EVAL, output of EVALUATE_ON_SPARSE_GRID. 
%
%       More precisely, if F: R^N -> R^V and SR.KNOTS contains M points in R^N, then 
%       SR.KNOTS is a matrix NxM and NODAL_VALUES is a matrix VxM, i.e.,
%       evaluations of F on different points of SR must be stored in NODAL_VALUES as columns. 
%       If F is a scalar-valued function then NODAL_VALUES is a row-vector
%
%       DOMAIN is a 2xN matrix = [a1, a2, a3, ...; b1, b2, b3, ...] defining the lower and upper bound
%       of the hyper-rectangle on which the sparse grid is defined
%       The function returns the Legendre expansion as a matrix of coefficients MODAL_COEFFS, where
%       the rows of MODAL_COEFF have V components (V columns), i.e. they are the "function coefficient" of the expansion
%       of F: R^N -> R^V. K is a matrix containing the associated multi-indices, one per row.
%
%       
%       
%
% [MODAL_COEFFS,K] = CONVERT_TO_MODAL(S,SR,NODAL_VALUES,DOMAIN,'chebyshev') returns the Chebyshev expansion
%        of the sparse grid interpolant. See above for inputs and outputs.
%
%
% [MODAL_COEFFS,K] = CONVERT_TO_MODAL(S,SR,NODAL_VALUES,DOMAIN,'hermite') returns the Hermite expansion
%       of the sparse grid interpolant. Here, DOMAIN is a 2XN matrix = [mu1, mu2, mu3, ...; sigma1, sigma2, sigma3,...]
%       i.e. the n-th variable of the sparse grid space has normal distribution with mean mu_n and std sigma_n
%
%
% [MODAL_COEFFS,K] = CONVERT_TO_MODAL(S,SR,NODAL_VALUES,DOMAIN,{<family1>,<family2>,<family3>,...}) returns the expansion
%       of the sparse grid interpolant over polynomial of "mixed" type, according to the families specified in
%       the last argument. For example:
%
%       CONVERT_TO_MODAL(S,SR,NODAL_VALUES,DOMAIN,{'legendre','hermite','legendre','hermite'})
%
%       converts the sparse grid interpolant in a sum of multi-variate polynomials that are products of
%       univariate Legendre polynomials (directions 1 and 3) and Hermite. Here DOMAIN is a 2XN matrix 
%       containing the parameters for each directions as column vectors. E.g. in the case above
%
%       DOMAIN = [a1, mu1, a2, mu2; b1, sigma1, b2 sigma2]


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile, B. Sprungk
% See LICENSE.txt for license
%----------------------------------------------------


% fix input 
errmsg=[' Please note that CONVERT_TO_MODAL does not accept INTERVAL_MAP '...
        'input argument any longer, and FLAG input argument is now mandatory. '...
        'Type help convert_to_modal for help. '...
        'This error message will not be shown in future releases of SPARSE-GRID-MATLAB-KIT'];

if nargin==4
    error('SparseGKit:WrongInput',strcat('not enough input arguments.',errmsg))
end

if any(~ismember(flags,{'legendre','chebyshev','hermite'}));
    error('SparseGKit:WrongInput',strcat('One or more strings in FLAGS unrecognized. ',errmsg));
end

% N is be the number of random variables
N=size(Sr.knots,1);
% V is the output dimensionality, F:R^N -> R^V
V=size(nodal_values,1);

% one tensor grid at a time, compute the modal equivalent, then sum up
% how many grid to process?
nb_tensor_grids = size(S,2);


% each tensor grid will result in vector of coefficients, u, that
% are the coefficients of the modal (Legendre,Hermite) expansion of the 
% lagrangian interpolant. To combine these coefficients properly, I need 
% to store not only the coefficients u themselves, but also the associated
% multi-indices. I store this information into a structure array, U, 
% whose lenght is nb_tensor_grids, with fields U(i).modal_coeffs, U(i).multi_indices, U(i).size (i.e. 
% how many polynomials are associated to each grid)

% I will need this counter
global_values_counter=0; 

% I also do a trick to preallocate U
U(nb_tensor_grids).multi_indices=[];
U(nb_tensor_grids).size=[];
U(nb_tensor_grids).modal_coeffs=[];

for g = 1: nb_tensor_grids
    
    % some of the grids in S structure are empty, so I skip them
    if isempty(S(g).knots)
        continue
    end

    % recover the values of f on the current grid. 
    % To do this, I need to use Sr.n, that is the mapping from [S.knots] to Sr.knots
    % To locate the knots of the g-th grid in the mapping Sr.n, i use global_values_counter    
    grid_knots = global_values_counter + 1 : ( global_values_counter + S(g).size );
    
    % extract the values of the g-th grid. If nodal_values is vector-valued, because F:R^N -> R^V
    % then the evaluation of each component of F, F_1, F_2, ... F_V on each point of the tensor
    % grid must be in a column vector, therefore I need to transpose nodal_values
    S_values = nodal_values(:,Sr.n(grid_knots))';

    % then update global_values_counter
    global_values_counter = global_values_counter + S(g).size ;
        
    % and compute modal
    U(g) = compute_modal_tensor(S(g),S_values,domain,flags);
    
end

% now I need to put together all the U.multi_indices, summing the coefficients that
% share the same multi_index. Same procedure as reduce sparse grid

% First I build a list of all multi_indices and coefficients
% I preallocate them to speed up code

tot_midx=sum([U.size]);

% this is the container of multi_indices
All = zeros(tot_midx,N);

% this is the container of coefficients
all_coeffs = zeros(tot_midx,V);

% this is the index that scrolls them
l=0; 

% fill All and all_coeffs
for g = 1: nb_tensor_grids

    if isempty(S(g).knots)
        continue
    end

    All(l+1:l+U(g).size,:) = U(g).multi_indices;
    all_coeffs(l+1:l+U(g).size,:) = U(g).modal_coeffs*S(g).coeff; % note that this one works even if on the left of = we want
                                                                  % a row vector (row of all_coeffs) and on the right we have
                                                                  % a column vector
    l=l+U(g).size;
end


% get rid of identical ones, by sorting and taking differences between cosecutive rows
% I need to store the sorter vector, to sum up all coefficients with
% the same multi_index

[All_sorted,sorter] = sortrows(All);

% taking differences between consecutive rows, I have a way
% to build the unique version of All. Remember diff has 1 row
% less, so I have to add a ficticious row at the end of All

dAll_sorted=diff( [All_sorted; ...
                   All_sorted(end,:)+1] );

% if a row of dAll_sorted is made of 0s then I have
% to discard the corresponding row of All

selector = sum(abs(dAll_sorted),2); % use abs(dAll_sorted), difference could be negative!
selector = selector>0;

% store uniqued version of All
K = All_sorted(selector,:);

% now I have to sum all the coefficients.

% sort coefficients in the same way as multi_indices
all_coeffs_sorted = all_coeffs(sorter,:);

% initialize modal_coeffs, that is the final output
nb_modal_coeffs = size(K,1);
modal_coeffs = zeros(nb_modal_coeffs,V);

% compute them one at a time

l=1; %scrolls all_coeffs_sorted
for k= 1 : nb_modal_coeffs
    ss=all_coeffs_sorted(l,:);
    while selector(l)~=1
        l=l+1;
        ss = ss+all_coeffs_sorted(l,:);
    end
    modal_coeffs(k,:) = ss; 
    l=l+1;
end