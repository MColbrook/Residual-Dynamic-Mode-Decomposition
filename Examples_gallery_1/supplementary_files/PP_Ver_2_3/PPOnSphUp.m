function z = PPOnSphUp(z,w,c,pres,tjmp)
% phase plot of function on upper hemishphere
%
% Usage: function z = PPOnSphUp(z,w,c,pres,tjmp)
%
% z - variable on domain covering exterior of unit circle, e.g. z=0.8./zdomain; 
% w - values of function at points z
% c - color scheme (optional)
% pres - resolution of phase (optional)
% tjmp -  jumps of phase (optional)


% Part of the phase plot package
% Version 2.3, January 15, 2014
% Copyright (c) 2012-2014 by Elias Wegert (elias@wegert.com, www.wegert.com)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2
 z = fldplane(1200);
 c = 'p';
 w = (z-1)./(z.^2+z+1);
end

xmax = max(max(real(1./z))); xmin = min(min(real(1./z))); 
ymax = max(max(imag(1./z))); ymin = min(min(imag(1./z))); 

if xmax<1 || xmin >-1 || ymax<1 || ymin >-1
  disp(' '), disp('Warning: z does not cover exterior of unit circle'); disp(' ')
end

if nargin==5
  RGB = colscheme(w,c,tjmp,pres); 
elseif nargin==4
  RGB = colscheme(w,c,[],pres); 
elseif nargin==3
  RGB = colscheme(w,c,[],30); 
elseif nargin==2 
  RGB = colscheme(w,'m'); 
elseif nargin<=1 
  RGB = colscheme(w,c,[],30); 
end

%[p,q,r] = StereoP2S(z);

%% plane plot of upper hemisphere

zext = abs(z)<=1.0;

% brighten colormap outside the unit circle 

RGB = BrightenRGB(RGB,0.7*zext);

h = surf(-real(1./z),-imag(1./z),zeros(size(z)),RGB);

set(h,'EdgeColor','none');
axis equal
axis off
axis tight
view(180,90)

ax = 1.2*[-1,1,-1,1];

axis(ax);

title('upper hemisphere'); 


