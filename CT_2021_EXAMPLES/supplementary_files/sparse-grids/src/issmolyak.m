function iss = issmolyak(S)

% ISSMOLYAK(S) returns 1 if S is a smolyak sparse grid. A smolyak sparse grid is a vector of structs with fields 
% 'knots','weights','size','knots_per_dim','m','coeff','idx'.


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

if isstruct(S)
    iss=isempty(setxor(fieldnames(S),{'knots','weights','size','knots_per_dim','m','coeff','idx'})) && length(S)>=1;
else
    iss=0;
end

