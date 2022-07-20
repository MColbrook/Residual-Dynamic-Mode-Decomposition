function S = asreduced(pts_list,wgs_list)

% S = ASREDUCED(pts_list) returns a reduced sparse grids structure whose "knots" field is set equal to
% pts_list, "size" filed is set to the number of columns of PTS_LIST, while the other fields are []. Thus, isreduced(S)==true and
% S can then be used as input to e.g. PLOT_GRID or EVALUATE_ON_SPARSE_GRID. The knots in pts_list are supposed
% to be all different, i.e. the list is not checked for duplicates
%
% S = ASREDUCED(pts_list,wgs_list) also set the weights field of S as wgs_list. An error is thrown if the number of weights and
% nodes is not identical


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



S.knots = pts_list;
if nargin==2
    if length(S.weights) == size(pts_list,2)
    S.weights=wgs_list;
    S.size=length(S.weights);
    else
        error('SparseGKit:WrongInput','the number of points and of weights are different')
    end
else
    S.weights=[];
    S.size=size(pts_list,2);
end
S.n=[];
S.m=[];
