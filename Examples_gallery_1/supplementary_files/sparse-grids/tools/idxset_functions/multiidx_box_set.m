function [C_with,C_without] = multiidx_box_set(shape,min_idx)

% [C_with,C_without] = MULTIIDX_BOX_SET(ii,min_idx)
%
% given an index ii, generates C_with, the box indeces set up to that ii.
% min_idx is either 0 or 1. C_without is C_with without ii itself


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



N=length(shape);
shape_fun=@(i) i-shape(1:length(i)); 
C_with=multiidx_gen(N,shape_fun,0,min_idx);

C_without=C_with;
C_without(end,:)=[];

