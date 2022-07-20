function [tocomp_list,recycle_list,recycle_list_old,discard_list] = lookup_merge_and_diff(pts_list,pts_list_old,Tol)

% [tocomp_list,recycle_list,recycle_list_old,discard_list] = lookup_merge_and_diff(pts_list,pts_list_old,Tol)
%
% looks for points of pts_list in pts_list_old using the same algorithm as reduce_sparse_grid. tol is
% the tolerance for two points to be considered equal. discard_list is the list of points that no longer belong to the grid
% (may happen for non-nested grids)


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

% declare a global variable controlling verbosity
global MATLAB_SPARSE_KIT_VERBOSE
if isempty(MATLAB_SPARSE_KIT_VERBOSE)
    MATLAB_SPARSE_KIT_VERBOSE=1;
end

if nargin==2
    Tol=1e-14;
end

N=size(pts_list,1);
N_old=size(pts_list_old,1);

% the list of indices of points to be evaluated. Init to its max length, will be cut after the search is over
tocomp_list=zeros(N,1);
% if grid are not nested, not all of the old grid will be recycled. Thus,
% we also need 2 list of indices of points to be recycled, storing the
% positions  in the old grid and in the new one
recycle_list_old=zeros(N_old,1);
recycle_list=zeros(N,1);
discard_list=zeros(N_old,1);


% first, merge the two lists 
Merged = [pts_list_old; pts_list];


% and a safety check: are we using a sufficiently fine tol when detecting identical points?
detect_unsufficient_tolerance(Merged,Tol);


% next, let's order the rows of Merged in lexicographic order, obtaining Sorted. If I use mysortrows then two rows like
%
% [a b c d]
% [a-t b c+t d]
% 
% are considered equal (if t < Tol ) and therefore placed one after the other
%
% sorter is an index vector that maps Merged into Sorted, i.e. Merged(sorter,:)==Sorted

[Sorted,sorter] = mysortrows(Merged,Tol/size(Merged,2));

% I also need to remember which points come from pts_list and which from
% pts_list_old, and from which original position. 
% Thus I create a flag vector [-1 -2 ... -N_old 1 2 .. N], i.e. positive
% flags identify the new grid and negative ones the old grid, then I sort
% it according to sorter

flags=[-1:-1:-N_old, 1:N];
flags_sorted=flags(sorter);


% next I take the difference of two consecutive rows. if the difference is small, then the rows are the same, i.e. the knot is the same
dSorted=diff(Sorted,1,1);

% I measure the difference with infty norm instead of L2 norm:
% i take  the maximum component of each row (2 means "operate on columns"):
%       max(abs(dSorted),[],2)    
% then I want to see which ones have this max bigger than Tol
%       diff_eq=max(abs(dSorted),[],2)>Tol
% this command returns a vector of True and false
%       diff_eq=[1 1 0 1 1]
% this means that the 2nd point is different from the 1st, the 3rd from the 2nd, but
% the 4th is equal to the 3rd ( diff(3)=v(4)-v(3) ), hence in common between the
% grids.

diff_eq=max(abs(dSorted),[],2)>Tol;

% now I scroll diff_eq and sort out everything according to these rules:
%
% --> if diff_eq(k)==0, 
%
% in this case either Sorted(k+1) is in the new grid and Sorted(k)
% is in the old grid and both are equal or viceversa. 
%
% Therefore, the point in the new grid goes into recycle_list and the
% old one in recycle_list_old. 
%
% Then I can skip the following (since it's equal and I have sorted it
% already)
%
% --> else, diff_eq(k)==1. 
%
% In this case, 4 cases are possible but actually only two matters
%
% -----> both Sorted(k) and Sorted(k+1) comes from the old_grid. 
%   then, Sorted(k) is to discard
%
% -----> both Sorted(k) and Sorted(k+1) comes from the new_grid. 
%   then, Sorted(k) is to compute
%
% -----> Sorted(k) is new, Sorted(k+1) is old. 
%   then, Sorted(k) is to compute
%
% -----> Sorted(k) is old, Sorted(k+1) is new. 
%   then, Sorted(k) is to discard



i=0; % scrolls recycle lists
j=0; % scroll compute_list
k=1; % scrolls diff_eq
L=length(diff_eq);
discard=0;

while k<=L
    
    if diff_eq(k) % short-hand for diff_eq(k)==1, i.e. two consecutive entries are different
        
       % compute or discard        
       if flags_sorted(k)>0
           j=j+1; 
           tocomp_list(j)=flags_sorted(k);
       else
           discard=discard+1;
           discard_list(discard)=-flags_sorted(k);
       end
       
       % then move to the following
       k=k+1;
       
    else % in this case diff_eq(k)==0, i.e. two consecutive entries are equal
         
        % recycling case
        i=i+1;
        
        if flags_sorted(k)>0    
            recycle_list(i)=flags_sorted(k);
            recycle_list_old(i)=-flags_sorted(k+1);
        else
            recycle_list_old(i)=-flags_sorted(k);
            recycle_list(i)=flags_sorted(k+1);
        end
        
        % then I can skip the k+1 because I have already sorted it
        k=k+2;
    end
    
end

% need to handle the case k==L, since diff is 1 element shorter than sorted. Note
% that 
% --> the node in Sorted(L,:) has been already taken care of inside the while loop
% --> we only need to do something if the diff_eq(L)==1. Indeed, if 
% diff_eq(L)==0, then Sorted(L+1,:) is already taken care of (and in this
% case the final value of k is L+2)

if diff_eq(L) % short-hand for diff_eq(L)==1,
    if flags_sorted(L+1)>0 % then it's a new point and has to be computed
        j=j+1;
        tocomp_list(j)=flags_sorted(L+1);
    else % then it's an old point and has to be discarded
        discard=discard+1;
        discard_list(discard)=-flags_sorted(k);
    end
end

% show some statistics
if MATLAB_SPARSE_KIT_VERBOSE
    disp(strcat('new evaluation needed:',num2str(j),' recycled evaluations:',num2str(i),' discarded evaluations:',num2str(discard)))
end


% remove the extra entries of tocomp_list, recycle_lists, discard_list. Pay attention to special cases

if j~=N % in this case there are no points to recycle and we have completely filled  tocomp_list
    if tocomp_list(j+1)~=0,
        error('SparseGKit:FailedSanityChk','tocomp_list(j+1)~=0'),
    end
    tocomp_list(j+1:end)=[];
end

if i==N % in this case the two grids are the same one
    warning('SparseGKit:GridsAreEqual','the two grids are the same!') 
else
    if i>N
        error('SparseGKit:FailedSanityChk',['The code has detected more points to recycle than points in the new sparse grid! ',...
               'Double check the values of tolerances used to detect identical points (both here and in reduce_sparse_grid) ',...
               'and rerun the code. i>N'])
    end
    if recycle_list(i+1)~=0,
        error('SparseGKit:FailedSanityChk','recycle_list(i+1)~=0'),
    end
    recycle_list(i+1:end)=[];
end

if i~=length(recycle_list_old) % in this case we haven't exhausted the old_grid (non-nested case)
    if recycle_list_old(i+1)~=0,error('SparseGKit:FailedSanityChk','recycle_list_old(j+1)~=0'),end
    recycle_list_old(i+1:end)=[];
end

if discard~=N_old % in this case we are discarding completely the old grid (non-nested case)
    if discard_list(discard+1)~=0,error('SparseGKit:FailedSanityChk','discard_list(discard+1)~=0'),end
    discard_list(discard+1:end)=[];
end



% safety checks
if i+discard~=N_old
    error('SparseGKit:FailedSanityChk','The code has lost track of some points of the old grid, i+discard~=N_old')
end

if length(recycle_list)~=length(recycle_list_old),
    error('SparseGKit:FailedSanityChk','mismatch between the two sets of recycling points. length(recycle_list)~=length(recycle_list_old)'),
end
if ~isempty( setxor( [tocomp_list; recycle_list], 1:N ) ),
    error('SparseGKit:FailedSanityChk',[ 'The code has lost track of some points of the new grid, ',...
            'or some points from the old grid have been mistaken as points of the new grid. ',...
            'Double check the values of tolerances use to detect identical points (both here and in reduce_sparse_grid) ',...
            'and try to rerun the code. ~isempty(setxor([tocomp_list recycle_list],1:N))']),
end


end