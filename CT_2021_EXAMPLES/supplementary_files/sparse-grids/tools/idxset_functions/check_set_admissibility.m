function [adm,C] = check_set_admissibility(I)

% [adm,C] = CHECK_SET_ADMISSIBILITY(I) checks whether I is an admissible set. 
%       If that is the case, adm=true and C contains I ordered in lexicographic order. 
%       If not, adm=false and C contains I plus the multiindeces needed, again in lexicographic order


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



C=I;

% now check amdissibility condition and add what's missing. Print a warning if needed
row = 0;
row_max=size(C,1);
adm=true;

while (row<row_max)
    
    % current index
    row=row+1;
    idx=C(row,:);
    [is_adm, ~, missing_set] = check_index_admissibility(idx,C);
    
    % if it is not admissible, add what's needed to the bottom of C
    % and increase row_max, so that the new indices will also be checked.
    % So we don't sort the C after the update, otherwise I am not sure I am
    % checking everything
    
    if ~is_adm
        % update C and counter
        C = [C; missing_set];
        row_max=size(C,1);
        
        % print a warning
        disp(strcat('the set is not admissible. Adding: ',num2str(missing_set) ) )

        % record non admissibility has been found
        adm=false;
        
    end
    
end

% finally, sort C
C=sortrows(C);