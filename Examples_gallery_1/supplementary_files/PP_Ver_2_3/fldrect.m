function  z = fldrect(zll,zur,m,n)
% discrete field of complex numbers covering the domain of a function

% Usage: z = fldrect(zll,zur,m,n);
% zll - complex number, lower left corner of the rectangular domain
% zur - complex number, upper right corner of the rectangular domain
% m,n - number of discretization points in x and y direction

% part of the phase plot package - visualization of complex functions
% Version 2.3, January 15, 2014
% Copyright (c) 2012-2014 by Elias Wegert (elias@wegert.com, www.wegert.com)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1, zll=-1-1i; end;
if nargin<2, zur=1+1i; end;
if nargin<3, m=1000; end;
if nargin<4, n=1000; end;

a = min(real(zll),real(zur));
b = max(real(zll),real(zur));
c = min(imag(zll),imag(zur));
d = max(imag(zll),imag(zur));

x = linspace(a,b,m);

y = linspace(c,d,n);

z = ones(n,1)*x + 1i*y'*ones(1,m);

