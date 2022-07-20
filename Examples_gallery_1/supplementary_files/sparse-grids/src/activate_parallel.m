function [] = activate_parallel(numworkers)

%ACTIVATE_PARALLEL opens a parallel toolbox session
%
% ACTIVATE_PARALLEL() uses the default number of parallel workers
%
% ACTIVATE_PARALLEL(NUMWORKERS) specifies the number of parallel workers
 

%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

switch nargin
    case 0
        if verLessThan('matlab', '8.3')
            eval('matlabpool open local')
        else
            parpool
        end
    case 1
        if verLessThan('matlab', '8.3')
            eval(['matlabpool open local ',num2str(numworkers)])
        else
            parpool(numworkers)
        end        
end
