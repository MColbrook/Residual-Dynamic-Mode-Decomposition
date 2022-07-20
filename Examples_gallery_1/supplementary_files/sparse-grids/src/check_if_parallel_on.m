function is_on = check_if_parallel_on()


%CHECK_IF_PARALLEL_ON checks if a parallel toolbox session is open
%
% IS_ON = CHECK_IF_PARALLEL_ON() returns TRUE if a parallel session is open and FALSE otherwise
 

%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


if verLessThan('matlab', '8.3')
    is_on = logical( matlabpool('size') );
else
    is_on = ~isempty(gcp('nocreate'));
end
