function I = fast_TD_set(N,w)

% I = fast_TD_set(N,w) returns the multiindex set TD(w) in N dimensions 
% (one row per multiindex) i.e. {ii in N_+ : sum(ii-1) <= w} 


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



% initialize I
TDsize=nchoosek(N+w,N);

I=NaN(TDsize,N);

r1=1;

% compute one level at a time

for i=0:w
    
    % compute number of rows
    r = setsize(N,i);
    r2= r1 + r-1;
    
    % now we do the recursive call
    I(r1:r2,:)=fast_TD_rec(N,i,r);
    r1=r2+1;

end

% sort rows to get lexicographic order. Note that this is not much time consuming
I = sortrows(I);

% if nargin==3
%     I=I+base;
% end

I=I+1;



function I = fast_TD_rec(N,w,rows)

% I = fast_TD_rec(N,w,rows)
%
% actually returns the sub matrix of the multindices in N dimensions st sum(i)=w. 
% rows is the number of rows and is needed for the recursive base step


if N==1
    
    I=w*ones(rows,1);

else
    
    
    switch w
        
        case 0
            
            I=zeros(rows,N);
            
        case 1
            
            I = eye(N);
            
        otherwise
            
            % initialize with NaN is faster than 0
            I = NaN(rows,N);
            
            % loop on levels and fill the submatrices
            
            r1=1; % first row of the submatrix
            
            for k=0:w
                
                % size of the submatrix
                inner_rows = setsize(N-1,w-k);
                
                % second row of the submatrix
                r2=r1+inner_rows-1;
                
                % rows of the submatrix
                rows_idx=r1:r2;
                
                % the first column
                I(rows_idx,1)=k*ones(inner_rows,1);
                
                % the submatrix
                I(rows_idx,2:end)= fast_TD_rec(N-1,w-k,inner_rows);
                
                % update row indices
                r1=r2+1;
            end
    end
end






function r = setsize(N,w)

% r = setsize(N,w)
%
% returns the number of multiindices i in N dimensions s.t. sum(i)=w

r=nchoosek(N+w-1,N-1);