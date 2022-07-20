function [pts_in_S_only,pts_in_both_grids_S,pts_in_both_grids_S_old,pts_in_S_old_only] = compare_sparse_grids(S,Sr,S_old,Sr_old,Tol)

% COMPARE_GRIDS compares the points of two sparse grids 
% 
% [PTS_IN_S_ONLY, PTS_IN_BOTH_GRIDS_S, PTS_IN_BOTH_GRIDS_S_OLD, PTS_IN_S_OLD_ONLY] = COMPARE_SPARSE_GRIDS(S,SR,S_OLD,SR_OLD,TOL)
%
% compares the points in S and S_OLD, and determines and determines those in common, and those belonging
% exclusively to each of the two. Points are considered equal is coordinate-wise they are closer than TOL.
% SR is the reduced version of S and SR_OLD is the reduced version of S_OLD.
%
% The function returns 4 set of scalar indices, which refer to the columns of SR.KNOTS and S.KNOTS. 
% In this way, each point of SR.KNOTS and SR_OLD.KNOTS is catalogued as follows:
%
% 1) IN_S_ONLY contains the column index of points of SR.KNOTS that are not in the intersection
% of S and S_OLD
% 
% 2) IN_BOTH_GRIDS_S are the column index of the points of SR.KNOTS that are also in SR_OLD.KNOTS
%
% 3) IN_BOTH_GRIDS_S_OLD are the column index of the points of SR_OLD.KNOTS that are also in SR.KNOTS, i.e.
% 
% SR.KNOTS(:,IN_BOTH_GRIDS_S) == SR_OLD.KNOTS(:,IN_BOTH_GRIDS_S_OLD)
%
% 4) IN_S_OLD_ONLY are column index of the points of SR_OLD.KNOTS that are in SR_OLD but not in SR.
%
% Example: 
%
% Suppose
%
% SR.KNOTS =    [a b c d e;
%                x y z w t]
%
% SR_OLD.KNOTS =[a f b h;
%                x s y r];
%
% Then: 
%
% IN_S_ONLY             = [3 4 5];
% IN_BOTH_GRIDS_S       = [1 2];
% IN_BOTH_GRIDS_S_OLD   = [1 3];
% IN_S_OLD_ONLY         = [2 4];
%
% crucially, the same results would happen if SR_OLD.KNOTS are perturbed by numerical noise below toll
%
% SR_OLD.KNOTS =[a+n f b    h;
%                x   s y-n  r];
%
% COMPARE GRIDS is very efficient because it works with comparing integer indices coming from the multiindices
% in S.IDX and S_OLD.IDX as much as possible. See also LOOKUP_MERGE_AND_DIFF for the same kind of analysis
% based however on comparing the actual coordinates of the points.

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

if nargin==4
    Tol=1e-14;
end

% ----------------------------------------------------------------------------------------
%
%                            PART I: looking for MULTI-IDX 
%
% ----------------------------------------------------------------------------------------

% to get good performance, we proceed by tensor grids: we identify which tensor grids are in common between S
% and S_old, and we classify the points of the two sparse grids depending whether they belong to tensor grids in
% common between the two sparse grids or not. Thus, the first thing to do is to recover the set of multiindices
% corresponding to the tensor grids stored in S and S_old. Note that these are subsets of the index sets that
% generated the grids (i.e, they are only those who "survived" the combination techinque selection)

N = size(Sr.knots,1);

N_old = size(Sr_old.knots,1);

if N~=N_old, error('SparseGKit:FailedSanityChk','grids with different N'), end

nb_idx_S = length(S);
I_S = reshape([S.idx],N,nb_idx_S)'; 

nb_idx_S_old = length(S_old);
I_S_old = reshape([S_old.idx],N_old,nb_idx_S_old)'; 


% we now need to identify which multiindices (hence tensor grids) belongs to both sparse grids and which are
% in one grid only. To do so, we exploit the fact that I_S and I_S_old are lexicographically ordered: we pick
% the smallest of the two sets, and we look for all of its multiidx in the other set. We then generate 4
% vectors of indices that refer to the elements of I_S and I_S_old:
% idx_in_both
% idx_in_both_old
% idx_in_S_only
% idx_in_S_old_only



if nb_idx_S < nb_idx_S_old

    % these are convenience counters
    ctr_S_only = 1; % scrolls idx_in_S_only
    ctr_both = 1;   % scrolls idx_in_both*
                    % note that we do not need a counter for idx_in_S_old_only. This will be clearer in a moment
    
    % these are the vector addressed with the increasing counters: every time the we find a multiidx in common, we
    % append its location in both sets in the indices below. However, instead of actually appending and letting
    % the vectors grow, we preallocate them to full length instead. We'll remove the extra components in one-shot at the end of the procedure
    idx_in_both=zeros(1,nb_idx_S);
    idx_in_both_old=zeros(1,nb_idx_S_old);
    idx_in_S_only=zeros(1,nb_idx_S);
    % idx_in_S_old_only is handled differently. Instead of appending to a vector if a multiidx from S_old is not found,
    % we work on the complementary information: if we find a multiidx from S_old, it is marked as NaN in the right
    % position of the vector below, and at the end those unmarked are those that haven't been found.
    idx_in_S_old_only=1:nb_idx_S_old;

    for i = 1:nb_idx_S
        row_of_I_S = I_S(i,:);

        % this exploits that I_S_old is sorted lexicographically, so it's efficient
        [found,pos] = find_lexicographic(row_of_I_S,I_S_old,'nocheck'); 
        
        if found
            % note that both this one and the next one are sorted because i is growing and the set is sorted so 
            % we look for higher and higher idx in a sorted searchee set so pos keeps increasing
            idx_in_both(ctr_both) = i; 
            idx_in_both_old(ctr_both) = pos;
            ctr_both = ctr_both+1;
            % if found, we need to mark the row of S_old as "found", and at the end we will keep only the unmarked ones,
            % which are those that we have not found during the process, i.e. the indices in I_S_old only
            idx_in_S_old_only(pos)=NaN;
        else
            idx_in_S_only(ctr_S_only) = i;
            ctr_S_only = ctr_S_only+1;
        end
    end
    
    % we remove the extra components as promised
    idx_in_both(ctr_both:end)=[];
    idx_in_both_old(ctr_both:end)=[];
    idx_in_S_only(ctr_S_only:end)=[];
    idx_in_S_old_only = idx_in_S_old_only(~isnan(idx_in_S_old_only));
    
else % same procedure for the opposite case

    ctr_S_old_only = 1;
    ctr_both = 1;
    idx_in_both=zeros(1,nb_idx_S);
    idx_in_both_old=zeros(1,nb_idx_S_old);
    idx_in_S_only=1:nb_idx_S;
    idx_in_S_old_only=zeros(1,nb_idx_S_old);
    
    for i = 1:nb_idx_S_old
        row_of_I_S_old = I_S_old(i,:);

        [found,pos] = find_lexicographic(row_of_I_S_old,I_S,'nocheck'); 
        
        if found
            idx_in_both(ctr_both) = pos; 
            idx_in_both_old(ctr_both) = i;
            ctr_both = ctr_both+1;
            idx_in_S_only(pos)=NaN;
        else
            idx_in_S_old_only(ctr_S_old_only) = i;
            ctr_S_old_only = ctr_S_old_only+1;
        end
    end
    
    idx_in_both(ctr_both:end)=[];
    idx_in_both_old(ctr_both:end)=[];
    idx_in_S_old_only(ctr_S_old_only:end)=[];
    idx_in_S_only = idx_in_S_only(~isnan(idx_in_S_only));
   
    
end

% ----------------------------------------------------------------------------------------
%
%                       PART II: classifying points based on their origin
%
% ----------------------------------------------------------------------------------------

% now, we derive information about sparse grid points (i.e. whether points are in S, S_old, or both) from the multiindex information.
% At the end of this part, we will have split points in 4 categories, two for points in S and two for points
% in S_old
% pts_from_new_idx_in_Sr
% pts_from_common_idx_in_Sr 
% pts_from_common_idx_in_Sr_old
% pts_from_old_idx_in_Sr_old


% --------------------------------- part II.a -------------------------------------------

% We begin by noting that points that belong to tensor grids associated to multiidx in both sparse grids will be necessarily
% present in both grids and we will return them as "in_both_grids". For now, we collect their position in both grids,
% Sr and Sr_old, i.e. in which column of Sr.knots and Sr_old.knots they are.

% we start from Sr. This is the cointainer. Note that it will be not sorted at first, hence the name. we will
% sort it afterwards
if ~isempty(idx_in_both)
    pts_from_common_idx_in_Sr_nosort=[];
    pos_end=0;
    i_prev= 0;
    % for each multiidx, we recover where the points of the associated tensor grid have been placed in Sr,
    % thanks to the vector Sr.n, that says precisely this information. Of course, grids in common and not in
    % common are interwined in Sr, so we need to use the S(i).size information to move to the correct fragment of
    % Sr.n
    for i=idx_in_both
        pos = pos_end +sum([S(i_prev+1:i-1).size])+1;
        pos_end = pos+S(i).size-1;
        pts_from_common_idx_in_Sr_nosort = [pts_from_common_idx_in_Sr_nosort; Sr.n(pos:pos_end)]; %#ok<AGROW> % this has been profiled, not too expensive
        i_prev = i;
    end
    % the vector pts_from_common_idx_in_red_nosort is not uniqued, because points of different tensor grids
    % might be mapped to the same point in Sr. We need them uniqued though, but the official unique is slow, so
    % we reimplement it as follows:
    vs = sort(pts_from_common_idx_in_Sr_nosort);
    vd = diff(vs);
    pts_from_common_idx_in_Sr = vs([logical(vd); true]);
else
    pts_from_common_idx_in_Sr=[];
end

% same procedure to locate the common points in SR_old starting from the common indices
if ~isempty(idx_in_both_old)
    pts_from_common_idx_in_Sr_old_nosort=[];
    pos_end = 0;
    i_prev = 0;
    for i=idx_in_both_old
        pos = pos_end + sum([S_old(i_prev+1:i-1).size])+1;
        pos_end = pos+S_old(i).size-1;
        pts_from_common_idx_in_Sr_old_nosort = [pts_from_common_idx_in_Sr_old_nosort; Sr_old.n(pos:pos_end)]; %#ok<AGROW>
        i_prev = i;
    end
    vs = sort(pts_from_common_idx_in_Sr_old_nosort);
    vd = diff(vs);
    pts_from_common_idx_in_Sr_old = vs([logical(vd); true]);
else
    pts_from_common_idx_in_Sr_old=[];
end

% --------------------------------- part II.b -------------------------------------------

% next, we perform the complementary operation and identify the column index of the points in Sr
% corresponding to multiidx in S only, and same for S_old. Note that it's the multiidx that belongs to S only,
% not the point!  If the 1D knots have some degree of nestedness, it will actually happen that some of the 
% points of multiidx in S_only will also be present in grid of midx in common. 
% Pay attention however, it might happen that idx_in_S_only is empty, for instance if the "old" grid is 
% actually larger. In this case, we throw an error

if isempty(idx_in_S_only)
    error('SparseGKit:ToBeCoded','idx_in_S_only is empty: this case is not handled yet')
else
    pts_from_new_idx_in_Sr_nosort=[];
    pos_end=0;
    i_prev = 0;
    for i=idx_in_S_only
        pos = pos_end + sum([S(i_prev+1:i-1).size])+1;
        pos_end = pos+S(i).size-1;
        pts_from_new_idx_in_Sr_nosort = [pts_from_new_idx_in_Sr_nosort; Sr.n(pos:pos_end)]; %#ok<AGROW>
        i_prev = i;
    end
    vs = sort(pts_from_new_idx_in_Sr_nosort);
    vd = diff(vs);
    pts_from_new_idx_in_Sr = vs([logical(vd); true]);
end

% now, the same for S_old, i.e. the column index of the points in Sr_old corresponding to multiidx in S_old
% only. Pay attention however, it might happen that idx_in_S_old_only is empty ...
if ~isempty(idx_in_S_old_only)
    pts_from_old_idx_in_Sr_old_nosort=[]; 
    for i=idx_in_S_old_only
        pos = sum([S_old(1:i-1).size])+1;
        pos_end = pos+S_old(i).size-1;
        pts_from_old_idx_in_Sr_old_nosort = [pts_from_old_idx_in_Sr_old_nosort; Sr_old.n(pos:pos_end)]; %#ok<AGROW> 
    end
    vs = sort(pts_from_old_idx_in_Sr_old_nosort);
    vd = diff(vs);
    pts_from_old_idx_in_Sr_old = vs([logical(vd); true]);
end 




% ----------------------------------------------------------------------------------------
%
%                       PART III: preparing the final outputs
%
% ----------------------------------------------------------------------------------------

% here, we start from these vectors:
%
% -- pts_from_new_idx_in_Sr
% -- pts_from_common_idx_in_Sr 
% -- pts_from_common_idx_in_Sr_old
% -- pts_from_old_idx_in_Sr_old
%
% and at the end of this section, we'll have points categorized in the final lists:
%
% -- pts_in_S_only
% -- pts_in_both_grids_S,
% -- pts_in_both_grids_S_old,
% -- pts_in_S_old_only

% we'd like to proceed differently whether the points are nested or not, but that's not as easy as it seems.
% so we do not take this into account
    
% As we have said already points from new idx might need be new, but not necessarily. Indeed:
% --> they might again be in the common points (imagine adding [5x1] grid to a [3x1]. the new index [5x1] is new,
% but some of its points might be in the [3x1] - actually, all of them if nested points!
% --> they might be among the points of the idx we are trashing
% at the same time, points from old idx might also be in points from common idx!
%
% So, this is the strategy:
% before trashing points, we check if they are also in common. Otherwise, before trashing them, we make
% sure that they are not needed as points from new idx


if isempty(idx_in_S_old_only)
    % if we are not trashing multiidx, then we init the set of new points as the set of points coming from new
    % idx, we'll figure out soon which points are not really new ...  the rest of the outputs are empty for now
    pts_in_S_only = pts_from_new_idx_in_Sr;
    pts_to_recycle_from_S_old_only_list=[];
    pts_to_recycle_from_S_old_only_list_old=[];
    pts_to_discard_from_S_old_list=[];
    
else
    % first, we make sure that points coming from old idx in Sr_old are really going to be discarded,
    % or if they are also in points in common
    [~,pos_in_arg1,~] = intersect(pts_from_old_idx_in_Sr_old,pts_from_common_idx_in_Sr_old);
    pts_from_old_idx_in_Sr_old(pos_in_arg1)=[];
    
    % now, depending whether pts_from_old_idx_in_Sr_old is still not empty or not, we behave differently
    
    if ~isempty(pts_from_old_idx_in_Sr_old)
        % if they really are going to be discarded, let's see if perhaps we can use them in the points of the new grid
        % coming from new idx. Here we need to resort to comparing the coordinates ...

        if MATLAB_SPARSE_KIT_VERBOSE
            disp('(hopefully) small local comparison between coordinates of points...')
        end
        
        pts_list = Sr.knots(:,pts_from_new_idx_in_Sr)';
        pts_list_old = Sr_old.knots(:,pts_from_old_idx_in_Sr_old)';
        
        % do the comparison. Observe that the output give the position of the common points in the subparts of
        % Sr.knots and SR_old.knots, so we need them to take them back to the "global" counting
        [not_found_in_S_old_only_loc,to_recycle_from_S_old_only_list_loc,to_recycle_from_S_old_only_list_old_loc,to_discard_from_S_old_list_loc] = ...
            lookup_merge_and_diff(pts_list,pts_list_old,Tol);

        if MATLAB_SPARSE_KIT_VERBOSE
            disp('...done. Overall statistics:')
        end

        
        pts_in_S_only = pts_from_new_idx_in_Sr(not_found_in_S_old_only_loc);
        pts_to_recycle_from_S_old_only_list = pts_from_new_idx_in_Sr(to_recycle_from_S_old_only_list_loc);
        pts_to_recycle_from_S_old_only_list_old = pts_from_old_idx_in_Sr_old(to_recycle_from_S_old_only_list_old_loc);
        pts_to_discard_from_S_old_list=pts_from_old_idx_in_Sr_old(to_discard_from_S_old_list_loc);
        
        %         % some of these to_discard might also be already in common, so although we discard it from here we're not
        %         % really losing them... let's find out
        %         [~,in_to_disc,~] = intersect(pts_to_discard_from_S_old_list,pts_from_common_idx_in_Sr_old);
        %         pts_to_discard_from_S_old_list(in_to_disc)=[];
        
    else
        % if there are actually no points going to be lost, the init is like in the previous case, i.e.
        % the set of new points is for now the set of points coming from new idx, and the rest of the outputs are empty for now
        
        pts_in_S_only = pts_from_new_idx_in_Sr;
        pts_to_recycle_from_S_old_only_list=[];
        pts_to_recycle_from_S_old_only_list_old=[];
        pts_to_discard_from_S_old_list=[];
        
    end
end


% in any case, we look we look whether some of the points from the new idx are already in the list of points
% from common idx, and remove them from the in_S_only list. We use this call
%
% [pts_to_comp_already_in_common,in_to_comp,in_com] = intersect(to_comp_list,pts_from_common_idx_in_red);
%
% but we actually care only about the second output:
[~,in_to_comp,~] = intersect(pts_in_S_only,pts_from_common_idx_in_Sr);
pts_in_S_only(in_to_comp)=[];

% note: the points to be recycled are alredy accounted for, no need to add them anywhere
% note (@@): the label of points in common doesn't really matter: they will in any case appear in both
% sets of indices!!

% we are done!! let's set the outputs

pts_in_both_grids_S = sort([pts_to_recycle_from_S_old_only_list; pts_from_common_idx_in_Sr]);
pts_in_both_grids_S_old = sort([pts_to_recycle_from_S_old_only_list_old; pts_from_common_idx_in_Sr_old]);
pts_in_S_old_only = pts_to_discard_from_S_old_list;



% safety checks

L1 = length(pts_in_S_only);
L2 = length(pts_in_both_grids_S);
L3 = length(pts_in_both_grids_S_old);
L4 = length(pts_in_S_old_only);

if L3+L4~=size(Sr_old.knots,2)
    error('SparseGKit:FailedSanityChk','The code has lost track of some points of the old grid, i+discard~=N_old')
end

if L2~=L3
    error('SparseGKit:FailedSanityChk','mismatch between the two sets of recycling points. length(pts_in_both_grids_S)~=length(pts_in_both_grids_S_old)'),
end
if ~isempty( setxor( [pts_in_S_only; pts_in_both_grids_S], 1:size(Sr.knots,2) ) ),
    error('SparseGKit:FailedSanityChk',[ 'The code has lost track of some points of the new grid, ',...
            'or some points from the old grid have been mistaken as points of the new grid. ',...
            'Double check the values of tolerances use to detect identical points (both here and in reduce_sparse_grid) ',...
            'and try to rerun the code. ']),
end



% show some statistics
if MATLAB_SPARSE_KIT_VERBOSE
    disp(strcat('new evaluation needed:',num2str(L1),' recycled evaluations:',num2str(L2),' discarded evaluations:',num2str(L4)))
end



% if using_nested_points
% 
%     % if we know that points are nested, then it's very easy! there are no points in S_old only,
%     % they are all in common with S. The pts_rom_new_idx_in_Sr conversely are a superset
%     % of the pts_from_common_idx_in_Sr, so we only need to take the setdiff.  However, note that
%     % reduce_sparse_grid might have mixed labels, i.e. some points from the new idx in common
%     % with common_idx might have received the label of the point in common. So we cannot rely
%     % here on the labels in pts_from_new_idx and we generate the entire set of labels, i.e.
%     % instead of
%     % pts_in_S_only = setdiff(pts_from_new_idx_in_Sr,pts_from_common_idx_in_Sr);    
%     % we do      
%     
%     pts_in_S_only = setdiff((1:size(Sr.knots,2))',pts_from_common_idx_in_Sr);    
%     pts_in_both_grids_S = sort(pts_from_common_idx_in_Sr);
%     pts_in_both_grids_S_old = (1:size(Sr_old.knots,2))';
%     pts_in_S_old_only = []; 
%     % this trick is not needed below, where we compare more closely the sets (see comment row marked with @@
%     % below): indeed, regardless of the label they will in any case appear in both sets of indices, and we'll take 
%     % intersect, not setdiff !!!
%     
% else
