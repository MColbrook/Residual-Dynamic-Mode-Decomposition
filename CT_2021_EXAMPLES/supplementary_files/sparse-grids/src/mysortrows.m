function [A,i_ord]=mysortrows(A,Tol,i,i_ord,n)

% Similar to Matlab builtin function sortrows. Given a matrix A
% of real numbers, sorts the rows in lexicographic order; entries 
% that differ less than Tol are treated as equal (default Tol is 1e-14).
%
% usage: 
%   [B,i]=mysortrows(A,Tol)
%   input
%     A: input matrix
%     Tol: (optional) tolerance used to identify coincident entries 
%          (default 1e-14)
%   output
%     B: sorted matrix
%     i: index vector such that A(i,:)=B
%
% mysortrows has a recursive implementation.  
% recursive call: [A,i_ord]=mysortrows(A,Tol,i,i_ord,n) (with i index)
% sorts the submatrix A(i,n:end); modifies the matrix itself and 
% stores the rows permutation in the global vector i_ord.



%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


% main ides: what points are to be considered equal even though numerically different? Let us take consecutive
% differences between the n-th column of A sorted increasingly: any time such difference is >Tol, then we have a match. Observe that this is NOT stable:
% if all the points are separated by less then Tol (say e.g. a:Tol/10:b) they are considered the same
% even if a-b >> Tol. In other words, points are supposed to cluster in groups whose size is ``a few deltas'',
% separated by ``many deltas''.

% more specifically, we need a recursive call like
%
% [A,i_ord]=mysortrows(A,Tol,i,i_ord,n) 
%
% and we will need three auxiliary vectors.
%
% -> i_ord is the GLOBAL final output, that returns the overall sorting. It keeps been partially rewritten at each
% recursive call (i.e. any time a sort a submatrix I need to keep track of the swaps I am doing)
%
% -> i is the GLOBAL partial order: it is as 1:<submatrix_size> and says the global position of the submatrix
% before the ordering. 
%
% -> ii is the LOCAL final order of each submatrix I am considering
%
% E.g: suppose that after sorting the entries A([3 4 5 10 15],:) are ``Tol-similar'' and need to be sorted as A([3 5 10 15 4],:);
% then 
%
% i = [3 4 5 10 15],
% ii = [1 3 4 5 2],
% the entries of i_ord in position  [3 4 5 10 15] must be swapped to positions [3 5 10 15 4],
% AND we look for submatrices that need further ordering


% this if happens only at user calls, not at recursive calls
if nargin==1
    Tol=1e-14;
    n=1;
    i=1:size(A,1); 
    i_ord=i;
elseif nargin==2
    n=1;
    i=1:size(A,1);
    i_ord=i;
end

%disp(['mysortrows level ',num2str(n)])

% Note that although the sorting is according to the rightward part of A, we always sort the WHOLE line
[~,ii]=sort(A(i,n));   
A(i,:)=A(i(ii),:);      
i_ord(i)=i_ord(i(ii));  



if n<size(A,2)

    % main ides: what points are to be considered equal even though numerically different? Let us take consecutive
    % differences between the n-th column of A sorted increasingly: any time such difference is >Tol, then we have a match. Observe that this is NOT stable:
    % if all the points are separated by less then Tol (say e.g. a:Tol/10:b) they are considered the same
    % even if a-b >> Tol. In other words, points are supposed to cluster in groups whose size is ``a few deltas'',
    % separated by ``many deltas''.

    % to begin with, we therefore take differences between to points: whenever the difference is smaller than Tol,
    % we will ignore such point. The first and last points of the list have a special treatment. Indeed, say 
    % the 1st and 2nd points are equal: then we would get a 0 in their difference and lose track of the 1st point!
    % to fix this, we put a "1" in front of the diff vector. For reasons that will be clearer in a second, we also
    % append "1" at the end.    
    j=[1;diff(A(i,n))>Tol;1];
    
    % now j contains a sequence of "1" and "0", e.g. [1 1 0 0 0 1 1 0 0 0 1 0 0 0 1 1]. Every time we run into 
    % a "1 0" fragment, it means that a sequence of points to be considered as equal starts, and ends where we find the 
    % opposite fragement "0 1". We can find the positions of these fragments by taking diff of j and looking for
    % +1 and -1, and we bookkeep the positions in two additional vectors, jm (beginning of a fragment) and jp (end
    % of a fragment). Now it's clear why we appended "1" at the end of j. Indeed, if the point list ends like [b a a]
    % (one points and two occurrences of the same), then diff would be [1 0], meaning that a fragment begins but
    % does not end; hence we need to append 1. Observe that a point list ending in [.. c c b a] implies diff ends in [.. 0 1 1]
    % and appending 1 does not matter (we do NOT focus on consecutive "1" fragments)
    
    jm=find(diff(j)==-1);
    jp=find(diff(j)==1);
    
    % now we need to consider each repeated entry and sort the rightward part of the matrix according to the same
    % criterion --> recursive call. 
    
    for k=1:length(jm)
        % let us identify the rows of A (in the current ordering that need to be sorted)
        v=i(jm(k):jp(k))';
        % we need to sort the submatrix and track all swaps in i_ord, whose lenght is always that of the full matrix
        [A,i_ord]=mysortrows(A,Tol,v,i_ord,n+1);
    end
end