function zz = StereoS2P(x,y,z)
%  stereographic projection from sphere to plane
%
% Usage: zz = StereoS2P(x,y,z)

% Part of the phase plot package
% Version 2.3, January 15, 2014
% Copyright (c) 2012-2014 by Elias Wegert (elias@wegert.com, www.wegert.com)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zz = (x+1i*y)./(1-z);
