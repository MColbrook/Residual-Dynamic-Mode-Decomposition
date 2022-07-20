function [m] = lev2knots_2step(i)


%   relation level / number of points:
%    m = 2(i-1)+1, 
%
%   i.e. m=1,3,5,7,9,...
%  
%   [m] = lev2knots_2step(i)
%   i: level in each direction
%   m: number of points to be used in each direction


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

m = 2*(i-1)+1;
m(i==0)=0;