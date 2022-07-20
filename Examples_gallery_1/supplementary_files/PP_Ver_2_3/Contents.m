% Contents of the Phase Plot Package
%
% Phase Plot Package - Visualization of complex functions
% Version 2.3, January 15, 2014
% Copyright (c) 2012-2014 by Elias Wegert (elias@wegert.com, www.wegert.com)
%
% For background information on phase portraits see 
% E. Wegert, Visual Complex Functions, Birkhaeuser-Springer 2012
%
% -----------------------------------------------------------------------------
%
% Demonstrations
%
% PPDemo        - basic demonstration of phase plots 
% PPGUI         - graphical user interface
%
% -----------------------------------------------------------------------------
%
% Main routines 
%
% PhasePlot    - phase plot of complex function
% PPOnSphere   - phase plot of function on Riemann sphere
% PPOnSphLow   - stereographic projection of phase plot from the lower hemisphere
% PPOnSphUp    - stereographic projection of phase plot from the upper hemisphere
%
% -----------------------------------------------------------------------------
%
% auxiliary functions
%
% colscheme  - algorithms for converting complex numbers to colors 
% myfunction - example of matlab routine (singular inner function) 
% fldrect    - discrete field of complex numbers, domain of function 
% unitcircle - vector with discrete points on complex unit circle
% zdomain    - the same as fldrect
% zplane     - discrete field of points in the entire complex plane
% zplanePP   - the same as zplane (to avoid confusion with existing m-file)
% 
% -----------------------------------------------------------------------------
%
% auxiliary functions, internal use
%
% BrightenRGB   - modification of color scheme 
% GUI_PhasePlot - graphical user interface 
% pltphase      - phase plot
% StereoP2S     - stereographic projection from plane to sphere
% StereoS2P     - stereographic projection from sphere to plane
% unitcirc      - uniformly distributed points on the complex unit circle
%
% -----------------------------------------------------------------------------
%
% Contents     - this file