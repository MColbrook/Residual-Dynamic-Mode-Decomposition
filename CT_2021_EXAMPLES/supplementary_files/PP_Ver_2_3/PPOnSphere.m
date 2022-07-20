function PPOnSphere(z,w,c,pres,tjmp)
% phase plot of function on Riemann shphere
%
% Usage: PPOnSphere(z,w,c,pres,tjmp)
%
% z - values on domain
% w - values of function at points z
% c - color scheme (optional)
% pres - resolution of phase (optional)
% tjmp -  jumps of phase (optional)

% Part of the phase plot package
% Version 2.3, January 15, 2014
% Copyright (c) 2012-2014 by Elias Wegert (elias@wegert.com, www.wegert.com)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  set lght = 1 for more realistic representation with light
lght = 0;

if nargin<2
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
  disp(' ')
  disp('The best way to make a phase plot on the entire sphere is')
  disp('the usage of two different charts (no singularities)')
  disp(' ')
  disp('First phase plot on lower hemishpere from transformed unit square:');
  disp(' ')
  disp('z = zdomain(-1-1i,1+1i,1000,1000);')
  disp('w = (z-1)./(z.^2+z+1);')
  disp('figure(1)')
  disp('hold on')
  disp('PPOnSphere(z,w);')
  z = zdomain(-1-1i,1+1i,1000,1000);
  w = (z-1)./(z.^2+z+1);
  figure(1)
  clf
  hold on
  PPOnSphere(z,w);
  disp('paused, press key to continue');
  pause
  disp(' ')
  disp('Next phase plot on upper hemishpere from inverted unit square:');
  disp(' ')
  disp('z = 1./z;')
  disp('w = (z-1)./(z.^2+z+1);')
  disp('PPOnSphere(z,w);')
  z = 1./z;
  w = (z-1)./(z.^2+z+1);
  PPOnSphere(z,w);
  disp('Phase Plot on upper hemishpere added');
  disp(' ')
  return
end

[p,q,r] = StereoP2S(z);

h=surf(p,q,r,RGB);

set(h,'EdgeColor','none');

axis equal
axis off
axis(0.8*[-1,1,-1,1,-1,1])

view(30,20)

if lght==1
    
lightangle(140,60)

set(gcf,'Renderer','zbuffer')
set(findobj(gca,'type','surface'),...
    'FaceLighting','phong',...
    'AmbientStrength',.7,'DiffuseStrength',.3,...
    'SpecularStrength',.5,'SpecularExponent',25,...
    'BackFaceLighting','unlit')

end

end

