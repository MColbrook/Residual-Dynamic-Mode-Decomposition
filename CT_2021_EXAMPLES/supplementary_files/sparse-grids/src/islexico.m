function isl = islexico(a,b)

% ISLEXICO checks lexicographic order of vectors
% ISLEXICO(A,B) for A and B vectors returns logical 1 (true) if A<=B in lexicographic sense,
% 
% examples:
%
% islexico([1 2 3],[1 2 5]) -> true  
% islexico([1 2 3],[1 2 3]) -> true
% islexico([1 2 3],[1 2 1]) -> false
% islexico([2 2 3],[1 7 1]) -> false
% islexico([1 7 3],[1 5 3]) -> false


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

[~,~,v]=find(a-b,1); % first nonzero element of a-b. if a==b, v is empty

% return true if a==b (i.e., v is empty) or if the first component for which a and b differ is smaller in a
% than in b
isl = isempty(v) || v<0 ;

