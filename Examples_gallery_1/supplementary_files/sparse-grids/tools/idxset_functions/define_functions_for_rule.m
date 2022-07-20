function [lev2nodes,idxset] = define_functions_for_rule(rule,input2)

% [lev2nodes,idxset] = DEFINE_FUNCTIONS_FOR_RULE(rule,N)
%
% sets the functions lev2nodes and idxset to use in smolyak_grid.m to build the desired ISOTROPIC sparse grid.
% N is the number of variables.
%
% [lev2nodes,idxset] = DEFINE_FUNCTIONS_FOR_RULE(rule,rates)
%
% sets the functions lev2nodes and idxset to use in smolyak_grid.m to build the desired ANISOTROPIC sparse grid 
% with specified rates. 
%
% rule can be any of the following: 'TP', 'TD', 'HC', 'SM'
% 
% The outputs are anonymous function , lev2nodes =@(i) ... and idxset=@(i) ...


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



if isscalar(input2)
    % here input2 is the number of variables
    N=input2;
    rates=ones(1,N);
else
    % here input2 is the rates vector
    rates=input2;
end

switch rule
    case 'TP'
        lev2nodes=@lev2knots_lin;
        idxset=@(i) max( rates(1:length(i)) .* (i-1) );
    case 'TD'
        lev2nodes=@lev2knots_lin;
        idxset=@(i) sum( rates(1:length(i)) .* (i-1) );
    case 'HC'
        lev2nodes=@lev2knots_lin;
        idxset=@(i) prod ( (i).^rates(1:length(i)) );
    case 'SM'
        lev2nodes=@lev2knots_doubling;
        idxset=@(i) sum( rates(1:length(i)) .* (i-1) );
end
