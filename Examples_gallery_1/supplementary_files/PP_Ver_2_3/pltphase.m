function PP = pltphase(z,f,cs,t,pres,height)
% phase plot of complex function f(z), internal use 

% Usage: PP = pltphase(z,f,cs,t,pres,height)

% Part of the phase plot package
% Version 2.3, January 15, 2014
% Copyright (c) 2012-2014 by Elias Wegert (elias@wegert.com, www.wegert.com)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%rep='a'; % standard analytic landscae 
%rep='l'; % logarithmic analytic landscae 
rep='c'; % compressed analytic landscae 

if nargin<6, height=0; end

if nargin<5 || isempty(pres), pres=20; end

if nargin<4 || isempty(t), t=exp(1i*linspace(-pi,pi,pres+1)); t=t(1:pres); end

if nargin<3 || isempty(cs), cs = 'p'; end

if nargin<2, z=zdomain; f=sin(5*z); end

% plot

if nargin<6, 
   
   % the phase plot is drawn on a surface
   
   if rep=='a'
   % analytic landscape
   height=abs(f); hmin=0; hmax=5;
   height = height.*(height<=hmax) + hmax*(height>hmax);
   f = f.*(height<hmax);
   
   elseif rep=='l'
   % logarithmic analytic landscape
   height=log(abs(f)); hmin =-2.5; hmax=2.5;
   height = height.*(height<=hmax).*(height>=hmin) ...
          + hmax*(height>hmax)+hmin*(height<hmin);
   f = f.*(height<hmax).*(height>hmin);
   
   elseif rep=='c' 
   % height such that zeros are at level 0 and poles are at level 1
   height=(2/pi)*atan(abs(f)); hmin=0; hmax=1;
   end
   
end

pres = length(t);
RGB = colscheme(f,cs,t,pres); 

PP = surf(real(z),imag(z),height,RGB);
set(PP,'EdgeColor','none');

axis equal
ax = axis;
axis([ax,hmin,hmax]);
axis off
view(30,20)
axis tight



