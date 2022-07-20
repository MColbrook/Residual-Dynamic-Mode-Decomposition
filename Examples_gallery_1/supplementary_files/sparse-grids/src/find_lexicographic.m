function [found,pos,iter] = find_lexicographic(lookfor,I,nocheck)

% FIND_LEXICOGRAPHIC finds specific rows of a matrix that is sorted lexicographically
%
% [FOUND,POS] = FIND_LEXICOGRAPHIC(LOOKFOR,I) looks for the row vector LOOKFOR among
%               the rows of I. If found, FOUND=TRUE and POS is the rownumber of LOOKFOR,
%               i.e. LOOKFOR == I(POS,:). If not found, FOUND=FALSE and POS=[]. The function
%               performs a preliminar check whether I is actually lexicographically sorted.
%   
% [FOUND,POS] = FIND_LEXICOGRAPHIC(LOOKFOR,I,'nocheck') is the same but does *not*  
%               check whether I is actually lexicographically sorted. This could be useful for speed purposes
%
%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


% first a check on I being sorted. Observe that the function sortrows used is very efficient
% so the cost of this preliminary analysis is negligible 

if nargin == 2 && ~issorted(I,'rows'),
    error('SparseGKit:SetNotSorted','I is not lexicographically sorted')
end

if nargin == 3 && ~strcmp(nocheck,'nocheck') 
    error('SparseGKit:WrongInput','unknown 3rd input')
end

nb_idx = size(I,1);

% we exploit the fact that I is sorted lexicographically and we proceed by binary search
% which guarantees cost ~ log(nb_idx)
%
% Basically, we start from the middle row, compare it with the index to be found, and
% if our index is larger we make a step in the increasing direction (i.e. we look in the upper half
% of the sorting), and viceversa. Of course, the step halves at each iteration: 
% therefore we necessarily terminate in ceil(log2(nb_idx)) steps at most

% the position to compare against -- if found, this is the position to be returned
idx = ceil(nb_idx/2);
% the associated vector 
jj = I(idx,:); 

found = isequal(jj,lookfor);

iter=1;
itermax = ceil(log2(nb_idx));

while ~found && iter <= itermax
    if islexico(jj,lookfor) % look in the upper half, unless we are already at the end
        if idx < nb_idx
            idx = min(idx+ceil(nb_idx/2^iter),nb_idx);
            jj = I(idx,:); % row entry of C2
            found = isequal(lookfor,jj);
            iter=iter+1;
        else
            break
        end
    else % look in the bottom half, unless we are already at the beginning
        if idx > 1
            idx = max(idx-ceil(nb_idx/2^iter),1);
            jj = I(idx,:); % row entry of C2
            found = isequal(lookfor,jj);
            iter=iter+1;
        else
            break
        end
    end
end

pos = [];
if found, 
    pos = idx;
end