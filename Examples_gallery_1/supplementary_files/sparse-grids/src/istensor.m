function ist = istensor(S)

% ISTENSOR(S) returns 1 if S is a tensor grid. A tensor grid is a struct with fields 
% 'knots','weights','size','knots_per_dim','m'.


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

if isstruct(S)
    ist=isempty(setxor(fieldnames(S),{'knots','weights','size','knots_per_dim','m'})) && length(S)==1;
else
    ist=0;
end

