function t=unitcirc(n)
%  uniformly distributed points on complex unit circle
%
% Usage: t = unitcirc(n)
% n+1 uniformly distributed points on the complex unit circle, t(n+1)=t(1)

% Part of the phase plot package
% Version 2.3, January 15, 2014
% Copyright (c) 2012-2014 by Elias Wegert (elias@wegert.com, www.wegert.com)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1, n=1024; end;

t = exp(1i*linspace(0,2*pi,n+1));
