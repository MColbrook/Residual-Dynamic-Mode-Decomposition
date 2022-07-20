function PP = PhasePlot(z,f,cs,pres,t,h)
% phase plot of complex function f(z)
%
% Usage: PP = PhasePlot(z,f,cs,pres,t,h)
%
% z    - complex field of arguments
% f    - complex field of values f(z)
% cs   - color scheme (optional)
%        call 'help colscheme' to see a list of available color schemes 
%        call 'PPDemo' to see a demonstration of all color schemes 
% pres - number of jumps in phase (optional)
% t    - positions of jumps at unit circle (optional)
% h    - height (z-axis) in which the plot is displayed (otional)
% 
% PP   - handle to colored surface

% Part of the phase plot package
% Version 2.3, January 15, 2014
% Copyright (c) 2012-2014 by Elias Wegert (elias@wegert.com, www.wegert.com)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==6
  PP = pltphase(z,f,cs,t,pres,h); 
elseif nargin==5
  PP = pltphase(z,f,cs,t,pres); 
elseif nargin==4
  if cs=='j'  
    PP = pltphase(z,f,cs,pres); 
  else
    PP = pltphase(z,f,cs,[],pres);
  end
elseif nargin==3
  PP = pltphase(z,f,cs);
elseif nargin==2
  PP = pltphase(z,f); 
else
  disp('PhasePlot(z,f) - phase plot of complex function f(z)')
  disp('  call as PhasePlot(z,f,colorscheme) for changing color options');
  disp('  call PPColorScheme to see a list of implemented color schemes');
  disp('  call PPDemo for a demonstration');
  return
end

axis equal
view(0,90)

