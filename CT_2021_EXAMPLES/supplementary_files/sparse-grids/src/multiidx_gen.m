function MULTI_IDX = multiidx_gen(N,rule,w,base,multiidx,MULTI_IDX)

% MULTI_IDX = multiidx_gen(N,rule,w,base,[],[])
%
% calculates all multi indexes M_I of length N with elements such that rule(M_I) <= w.
% M_I's are stored as rows of the matrix MULTI_IDX
% indices will start from base (either 0 or 1)


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


% multiidx_gen works recursively, exploring in depth the tree of all possible multiindexes. 
% the current multi-index is passed as 4-th input argument, and eventually stored in MULTI_IDX.
% The starting point is the empty multiidx: [], and MULTI_IDX is empty at the first call of the function.
% That's why the call from keyboard comes with [], [] as input argument: multiidx_gen(L,rule,w,[],[])

% disp('multiidx_fabio')

if nargin==3
      base = 0;   % default choice
      multiidx=[];
      MULTI_IDX=[];
elseif nargin==4
      multiidx=[];
      MULTI_IDX=[];
end


if length(multiidx)~=N 
      % recursive step: generates all possible leaves from the current node (i.e. all multiindexes with length le+1 starting from
      % the current multi_idx, which is of length le that are feasible w.r.t. rule)
      
      i=base;
      
      while rule([multiidx, i]) <= w
            % if [multiidx, i] is feasible further explore the branch of the tree that comes out from it.
            MULTI_IDX = multiidx_gen(N,rule,w,base,[multiidx, i],MULTI_IDX);
            i=i+1;
      end

%       for i=0:w
%             if rule([multiidx, i]) <= w
%                   % if [multiidx, i] is feasible further explore the branch of the tree that comes out from it.  
%                   MULTI_IDX = multiidx_gen(L,rule,w,[multiidx, i],MULTI_IDX);
%             end
%       end
      
else
      % base step: if the length of the current multi-index is L then I store it in MULTI_IDX  (the check for feasibility was performed in the previous call  
      MULTI_IDX=[MULTI_IDX; multiidx];
end