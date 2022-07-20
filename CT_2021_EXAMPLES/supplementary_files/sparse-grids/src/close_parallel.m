function [] = close_parallel()

%CLOSE_PARALLEL closes a parallel toolbox session
%
% CLOSE_PARALLEL() uses the default number of parallel workers
 

%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

if verLessThan('matlab', '8.3')
    eval('matlabpool close')
else
    p=gcp('nocreate');
    delete(p)
end
    