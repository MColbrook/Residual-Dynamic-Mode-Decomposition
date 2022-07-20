function  z = zplanePP(n)
% discrete field of points covering the entire complex plane
%
% Usage: z = zplanePP(n);

% Part of the phase plot package
% Version 2.3, January 15, 2014
% Copyright (c) 2012-2014 by Elias Wegert (elias@wegert.com, www.wegert.com)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1, n=1200; end;

[p,q,r] = sphere(n);

z = StereoS2P(p,q,r);

