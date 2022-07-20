function is_unsuf = detect_unsufficient_tolerance(pts,tol)

% is_unsuf = detect_unsufficient_tolerance(pts,tol)
%
% for a matrix of points (one per row, as in mysortrows or lookup_merge_and_diff), 
% computes an approximation of the support of the set of point in each direction 
% (as max - min coord in each dir) and compares it with tol. It they are "same-sized"
% then tol is too large and a warning is thrown


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


N=size(pts,2);
is_unsuf=false;

for n=1:N
    % check one column at a time. exploit the fact that unique returns a sorted sequence
    uu = unique(pts(:,n));    
    
    % if the tolerance is not at least 2 order of magnitude smaller than the support on the n-th direction 
    % throw a warning
    
    if abs(log10(uu(end)-uu(1))-log10(tol)) < 2 
        warning('SparseGKit:TolNotEnough','tol does not seem small enough to detect identical points')
        is_unsuf=true;
        return
    end
end

