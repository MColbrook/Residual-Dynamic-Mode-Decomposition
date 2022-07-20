function PPDemo(cs)
% demonstration of phase plots with various color schemes
%
% Usage: PPDemo, PPDemo(cs)
%
% cs - color scheme (for example 'p' or 'pn')
%
% call PPDemo without argument to see a demonstration of available color schemes
% call 'help colscheme' to get a list of available color schemes 

% Part of the phase plot package
% Version 2.3, January 15, 2014
% Copyright (c) 2012-2014 by Elias Wegert (elias@wegert.com, www.wegert.com)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% if there are memory problems reduce the resolution

xres = 800;
yres = 800;

%%

z = zdomain(-2-2i,2+2i,xres,yres);

clf

w = (z-1)./(z.*z+z+1);

if nargin == 1;
    %if strcmp(cs,'f') || strcmp(cs,'fn'), w=z+1./z; end
    PhasePlot(z,w,cs);
else
    
%% proper phase plot

PhasePlot(z,w,'p');
title('p - proper phase plot, hsv coloring')

disp(' ')
disp('DEMO: PHASE PLOTS OF COMPLEX FUNCTIONS')
disp(' ')
disp('To create a phase plot of w=f(z), call PhasePlot(z,w,cs), where') 
disp('z  is a 2D field of complex numbers covering the domain of f')
disp('w  is a 2D field of corresponding values f(z)')
disp('cs is a colorscheme (optional), for example cs=''p''')
disp(' ')
disp('The function f(z)=(z-1)/(z^2+z+1) is depicted in |Re(z)|<2, |Im(z)|<2')
pause

PhasePlot(z,w,'pn');
title('pn - proper phase plot, hsv coloring modified according to NIST standard')
disp(' ')
disp('Note that all color schemes can be modified to the NIST standard')
disp('by appending ''n'' to cs, for example cs=''pn''');
disp('For details see http://dlmf.nist.gov/help/vrml/aboutcolor');
pause

%% colored modifications


PhasePlot(z,w,'m');
title('m - enhanced phase plot with modulus jumps')
disp(' ')
disp('There are various modifications of the standard color scheme.')
disp('Here we use shading to generate contour lines of |f(z)|.')
pause

PhasePlot(z,w,'c');
title('c - phase plot with conformal polar grid')
disp(' ')
disp('Here we further added contour lines of arg f(z)')
disp('which generates tiles of almost square shape.')
pause

PhasePlot(z,w,'j');
title('j - phase plot with some enhanced isochromatic lines')
disp(' ')
disp('The contour lines of arg f(z) are isochromatic.')
pause

%PhasePlot(z,w,'q');
%title('q - phase plot colored in steps ')
%pause

PhasePlot(z,w,'d');
title('d - standard domain coloring')
disp(' ')
disp('This is standard domain coloring, where color corresponds to phase')
disp('and the gray values represent the modulus of f(z).')

pause

PhasePlot(z,w,'e');
title('e - enhanced domain coloring')
disp(' ')
disp('Here we see enhanced domain coloring with contour lines of |f(z)|.')
pause

%% black and white modifications

PhasePlot(z,w,'u');
title('u - polar chessboard')
disp(' ')
disp('The color scheme ''u'' is good to visualize conformal mappings.')
pause

PhasePlot(z,w,'a');
title('a - alternating black and white phase')
disp(' ')
disp('The color scheme ''a'' is appropriate to depict electric field lines.')
pause
hold on
plot3(cos(2*pi/3),sin(2*pi/3),20,'r.','Markersize',20);
plot3(cos(4*pi/3),sin(4*pi/3),20,'r.','Markersize',20);
plot3(1,0,1,'b.','Markersize',20);
disp('Here we have two positive and one negative point charges.')
pause

PhasePlot(z,w,'b');
title('b - alternating black and white modulus')
disp(' ')
disp('Using the scheme ''b'' highlights the corresponding potential lines.')
pause
hold off


%% color schemes based on cartesian coordinates

PhasePlot(z,w,'v');
title('v - Cartesian chessboard')

%% color scheme for depicting stream lines of potential flows

disp(' ')
disp('The scheme ''f'' is designed for depicting stream lines of potential flows.')
disp('Here f(z)=(z+1/z)/2 is the Joukovski function.')
PhasePlot(z,z+1./z,'f');
pause
disp('The function in the interior of the unit circle is irrelevant.')
t = exp(1i*linspace(0,2*pi,1000));
patch(real(t),imag(t),20*ones(size(t)),'k');
title('Plane potential flow around a circular cylinder')
pause

disp(' ')
disp('Call PPDemo(cs) to see the standard example represented by color scheme cs.');
pause

clf
disp(' ')
disp('Phase plots can also be depicted on the Riemann sphere using the routine PPOnSphere.');
disp('This code uses one chart to (almost) cover the sphere:')
disp(' ')
disp('z = zplanePP; w = (z-1)./(z.^2+z+1); PPOnSphere(z,w);')

z = zplanePP;
w = (z-1)./(z.^2+z+1);
PPOnSphere(z,w);
pause

disp(' ')
disp('Using two charts to cover the sphere yields a more uniform distribution of points:');
disp(' ')
disp('z1= zdomain(-1-1i,1+1i); w1=(z1-1)./(z1.^2+z1+1); PPOnSphere(z1,w1);')
disp('z2=1./z1; w2=(z2-1)./(z2.^2+z2+1); PPOnSphere(z2,w2);')

clf, hold on
z = 2*zdomain(-1-1i,1+1i,1000,1000);
w=(z-1)./(z.^2+z+1);
PPOnSphere(z,w);
zz=1./z;
ww=(zz-1)./(zz.^2+zz+1);
PPOnSphere(zz,ww);
hold off
pause

disp(' ')
disp('Visit a gallery of phase plots at www.visual.wegert.com')
disp(' ')
disp('For detailed information read')
disp('E. Wegert: Visual Complex Functions, Springer/Birkhaeuser 2012')
pause

disp(' ')
disp('Now invoking the phase plot user interface PPGUI ...')
PPGUI

end

