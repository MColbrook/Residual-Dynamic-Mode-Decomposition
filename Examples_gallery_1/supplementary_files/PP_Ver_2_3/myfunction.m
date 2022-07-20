function w = myfunction(z)
% example of complex function written as script for use with PPGUI
%
% this is an atomic singular inner function with atoms sitting at the
% fifth root of unity

% Part of the phase plot package
% Version 2.3, January 15, 2014
% Copyright (c) 2012-2014 by Elias Wegert (elias@wegert.com, www.wegert.com)

omega = exp(2i*pi/5);

w = ones(size(z));

for k=0:4
  w = w.*exp((z+omega.^k)./(z-omega.^k));
end

end
