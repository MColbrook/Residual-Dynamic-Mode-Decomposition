function m = lev2knots_doubling(i)

% m = lev2knots_doubling(i)
%
% relation level / number of points:
%    m = 2^{i-1}+1, for i>1
%    m=1            for i=1
%    m=0            for i=0
%
% i.e. m(i)=2*m(i-1)-1


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


m = 2.^(i-1)+1;
for k=1:length(m(:))
    if i(k)==1 
        m(k) = 1;
    end
    if i(k)==0
        m(k) = 0;
    end
end

