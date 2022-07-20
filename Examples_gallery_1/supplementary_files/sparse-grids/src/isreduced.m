function isr = isreduced(S)

% ISREDUCED(S) returns 1 if S is a reduced sparse grid. A reduced sparse grid is a struct with fields 
% 'knots','m','weights','n','size'.


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


if isstruct(S)
    isr=isempty(setxor(fieldnames(S),{'knots','m','weights','n','size'})) && length(S)==1;
else
    isr=0;
end


