function [is_adm, completed_set, missing_set] = check_index_admissibility(idx,idx_set,sort_option)

% [IS_ADM, COMPLETED_SET, MISSING_SET] = CHECK_INDEX_ADMISSIBILITY(IDX,IDX_SET) given a multiindex IDX 
%       as row vector checks if it is admissible w.r.t. the index set IDX_SET (matrix with indices as
%       rows). If it is so, IS_ADM = true. Otherwise IS_ADM = false and the function will return the 
%       COMPLETED_SET and a list of the indices added stored in MISSING_SET. 
%
%       Note that the starting IDX_SET is supposed to be admissible, so its multi-indices will not 
%       be checked for admissibility. Anyway, being admissible is not necessary for this function to work, 
%       since the core is SETDIFF, that does not assume any ordering
%
% [IS_ADM, COMPLETED_SET, MISSING_SET] = CHECK_INDEX_ADMISSIBILITY(IDX,IDX_SET,'sorting') returns
%       COMPLETED_SET and MISSING_SET sorted in lexicographic order.


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


is_adm=true;
completed_set=idx_set;
missing_set=[];

% now check if idx is admissible. Note that if this is not the case, the indices needed
% may not be included in idx_set, too ! 
% So we add everything we need to check in a queue and keep on checking while the queue is empty

% here's the queue
the_queue=idx;


while ~isempty(the_queue)
    
    % consider the first element of the stack
    i=the_queue(1,:);
    
    % build its needed set
    S = needed_set(i);

    % take the setdiff with idx_set. needed_rows is s.t. missing=needed_set(needed_rows,:)
    missing = setdiff(S,idx_set,'rows');

    % if missing is not empty, 
    if ~isempty(missing)
        
        % the initial set was not admissible
        is_adm=false;

        % missing had to be added to completed_set, missing_set, and the_queue. 
        % I have to make sure I don't add duplicates, 
        % and that the order of completed_set, missing_set, the_queue stays the same!
        % So I cannot use unique, since it automatically reorders rows
        % (and old versions of matlab like 2011b do not provide the 'stable'
        % option that deactivates the ordering).
        %
        % setdiff(A,B,'rows') returns the rows from A that are not in B.

        % add missing to completed_set
        completed_set=[completed_set; setdiff(missing,completed_set,'rows')]; %#ok<AGROW> 

        
        % add missing to missing set
        if isempty(missing_set)
            missing_set=missing;
        else
            missing_set=[missing_set; setdiff(missing,missing_set,'rows')]; %#ok<AGROW>
        end
        
        % add missing to the queue
        the_queue=[the_queue; setdiff(missing,the_queue,'rows')]; %#ok<AGROW>         
        

        % version with unique, does not work on R2011b
        %----------------------------------
        
        % add missing to completed_set
        % completed_set=unique([completed_set; missing],'rows','stable');
        
        % add missing to missing set
        % missing_set=unique([missing_set; missing],'rows','stable');
        
        % add missing to the queue
        % the_queue=unique([the_queue; missing],'rows','stable');
        
    end
    
    % delete current i from the queue
    the_queue(1,:)=[];
end
    

% if requested, sort in ascending lexicographic order
if nargin==3 
    
    if strcmp(sort_option,'sorting')
        missing_set=sortrows(missing_set);
        completed_set=sortrows(completed_set);
    else
        error('SparseGKit:WrongInput','unknown sorting option')
    end
end




function S = needed_set(idx)

% S = needed_set(idx)
%
% computes the indices of the form idx-e_j where e_j is the j-th N-dimensional unit vector,
% and store them as rows of S

N=length(idx);

% I can build the set quickly with matrices operation: [idx; idx; ... ; idx] - eye(N)

S=ones(N,1)*idx - eye(N);

% if idx has 1 inside, like [2 1 1], [3 1 2] etc, care has to be taken: the minimum value for indices is
% 1, so I don't have to check for [2 0 1], [2 1 0] etc to be in the set. I handle this deleting all rows
% that contain 0. Note that 0 can be only in the main diagonal of needed_set

D=diag(S);
S(D==0,:)=[];