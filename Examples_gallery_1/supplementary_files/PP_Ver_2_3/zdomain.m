function  z = zdomain(zll,zur,m,n)
% discrete field of complex numbers, domain of function
%
% Usage: z = zdomain(zll,zur,m,n);
% zll - complex number, lower left corner of the rectangular domain
% zur - complex number, upper right corner of the rectangular domain
% m,n - number of discretization points in x and y direction

% part of the phase plot package - visualization of complex functions
% Version 2.3, January 15, 2014
% Copyright (c) 2012-2014 by Elias Wegert (elias@wegert.com, www.wegert.com)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1, zll=-1-1i; end;
if nargin<2, zur=1+1i; end;
if nargin<3, m=1000; end;
if nargin<4, n=1000; end;

z = fldrect(zll,zur,m,n);
