function rgb = colscheme(f,cs,t,pres)
% various color schemes for use with phase plots
%
% Usage: rgb = colscheme(f,cs,t,pres)
%
% to get schemes bases on standard hsv coloring, specify string cs as follows: 
%
% a - alternating black and white phase
% b - alternating black and white modulus
% c - phase plot with conformal polar grid
% d - standard domain coloring
% e - enhanced domain coloring
% f - like 'y' but white and blue, especially for stream lines
% i - like 'y' with spacing depending on size of Im f 
% j - colored phase plot with specific phase jumps
% l - like 'y' but spacing in integer fractions of 2 pi
% m - colored phase plot with module jumps
% n - like 'c' - with brighter color for background 
% p - proper phase plot
% q - phase plot colored in steps 
% r - conformal cartesian grid
% s - conformal polar grid
% t - polar chessboard - light gray
% u - polar chessboard
% v - cartesian chessboard
% w - cartesian chessboard - light gray
% x - black and white stripes corresponding to real part
% y - black and white stripes corresponding to imaginary part
%
% appending 'n' to cs yields color schemes based on the NIST scaling of phase
% see  http://dlmf.nist.gov/help/vrml/aboutcolor

% Part of the phase plot package
% Version 2.3, January 15, 2014
% Copyright (c) 2012-2014 by Elias Wegert (elias@wegert.com, www.wegert.com)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if nargin<1, eval('help colscheme'); return, end
if nargin<2, cs='p'; end
if nargin<3 || isempty(t); 
    t=exp(1i*linspace(-pi,pi,21)); t=t(1:20); 
end
if nargin<4 || isempty(pres)
  phaseres = 20; 
  % number of enhanced isochromatic lines 
else
  phaseres=pres;
end

[m,n]=size(f);

if length(cs)==2 
   % used for modified color schemes (e.g. according to NIST)
   ms=cs(end); cs=cs(1:end-1);
else
   ms='s';   
end

if  length(cs)>=3 && isreal(cs)
   rgb(:,:,1) =  cs(1)*ones(m,n);
   rgb(:,:,2) =  cs(2)*ones(m,n);
   rgb(:,:,3) =  cs(3)*ones(m,n);
   return
end
    
%% choose basic color map

% continuous coloring 
if ms=='s' % standard hsv scheme
  colmap=hsv(600); 

elseif ms=='n' % NIST scheme, see http://dlmf.nist.gov/help/vrml/aboutcolor
  cmap=hsv(900);
  colmap(1:150,:)=cmap(1:150,:);
  colmap(151:300,:)=cmap(151:2:450,:);
  colmap(301:450,:)=cmap(451:600,:);
  colmap(451:600,:)=cmap(601:2:900,:);
else
  disp('color scheme not implemented')
end

% stepwise coloring (hsv scheme only)
if cs == 'q'
  colmap=hsv(pres);
end

palette = colormap(colmap);

% number of colors in palette encoding phase

p = size(palette);
p = p(1);

% encoding phase

myangle = (angle(-f)+pi)/(2*pi); 
nphase  = stepfct(myangle,p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% alternating black and white stripes corresponding to phase
%  also appropriate for electric field lines

if  cs == 'a' 
    
   black = sawfct(myangle,1/phaseres,0,1);   
      
   bmax = max(max(black));
   bmin = min(min(black));
   
   black = floor(2*(black-bmin)/(bmax-bmin));

   rgb(:,:,1) =  black.*ones(m,n);
   rgb(:,:,2) =  black.*ones(m,n);
   rgb(:,:,3) =  black.*ones(m,n);
   
end

%%%%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% alternating black and white stripes corresponding to modulus

if cs == 'b'  
        
   black = sawfct(log(abs(f)),2*pi/phaseres,0,1);
   
   bmax = max(max(black));
   bmin = min(min(black));
   black = floor(2*(black-bmin)/(bmax-bmin));
   
   rgb(:,:,1) =  black.*ones(m,n);
   rgb(:,:,2) =  black.*ones(m,n);
   rgb(:,:,3) =  black.*ones(m,n);
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% phase plot with conformal polar grid

if  cs=='c' 
    
   blackp = sawfct(myangle,1/phaseres,0.75,1);
   blackm = sawfct(log(abs(f)),2*pi/phaseres,0.75,1);
 
   black = blackp.*blackm;
   
   rgb(:,:,1) =  black.*reshape(palette(nphase,1),m,n);
   rgb(:,:,2) =  black.*reshape(palette(nphase,2),m,n);
   rgb(:,:,3) =  black.*reshape(palette(nphase,3),m,n);
   
   rgb = BrightenRGB(rgb,.1);
        
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  proper and enhanced domain coloring

if cs =='d' || cs=='e'   

   if cs =='d'
   
     rgb(:,:,1) =  reshape(palette(nphase,1),m,n);
     rgb(:,:,2) =  reshape(palette(nphase,2),m,n);
     rgb(:,:,3) =  reshape(palette(nphase,3),m,n);
   
   elseif cs=='e'
       
     black = sawfct(log(abs(f)),2*pi/phaseres,0.7,1);
 
     rgb(:,:,1) =  black.*reshape(palette(nphase,1),m,n);
     rgb(:,:,2) =  black.*reshape(palette(nphase,2),m,n);
     rgb(:,:,3) =  black.*reshape(palette(nphase,3),m,n);
  
   end

   % modified domain coloring
     
     % version 1
     logf   = log(abs(f));  
     bright = (logf>=0).*((1-(1./(1+logf))).^2)-(logf<0).*(1-(1./(1-logf))).^2;
      
     % version 2
     %bright = (3.2/pi)*atan(abs(f))-0.8;
   
     rgb = BrightenRGB(rgb,bright);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% alternating black and white stripes corresponding to imaginary part
%  in particular recommended for stream lines of potential flow

if cs == 'f' 
   
   imres = 8/phaseres;
   
   f = f-0.25i*imres;
   
   blue = sawfct(imag(f),imres,0,1);
   
   bmax = max(max(blue));
   bmin = min(min(blue));
   blue = floor(2*(blue-bmin)/(bmax-bmin));
   
   rgb(:,:,1) =  (1-blue).*ones(m,n);
   rgb(:,:,2) =  (1-blue).*ones(m,n);
   rgb(:,:,3) =  ones(m,n);
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% alternating black and white stripes corresponding to imaginary part

if cs == 'i' || cs== 'l' 
    
   impart = imag(f);
   
   if cs=='l'
     black = sawfct(imag(f),2*pi/pres,0,1);
   elseif cs=='i'
     fmax = mean(max(impart));
     fmin = mean(min(impart));
     imres = (fmax-fmin)/phaseres;
     black = sawfct(impart,imres,0,1);
   end
   
   bmax = max(max(black));
   bmin = min(min(black));
   black = floor(2*(black-bmin)/(bmax-bmin));
   
   rgb(:,:,1) =  black.*ones(m,n);
   rgb(:,:,2) =  black.*ones(m,n);
   rgb(:,:,3) =  black.*ones(m,n);
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% phase plot with specified jumps at points t

if cs == 'j' 
  
   t = exp(1i*sort(angle(t)));
  
   black = sawfctt(f,t);   
   
   rgb(:,:,1) =  black.*reshape(palette(nphase,1),m,n);
   rgb(:,:,2) =  black.*reshape(palette(nphase,2),m,n);
   rgb(:,:,3) =  black.*reshape(palette(nphase,3),m,n);
   
   %rgb = BrightenRGB(rgb,.1);
      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% contour lines and conformal grid

if cs=='m'  
   
   black = sawfct(log(abs(f)),2*pi/phaseres,0.75,1);
  
   rgb(:,:,1) =  black.*reshape(palette(nphase,1),m,n);
   rgb(:,:,2) =  black.*reshape(palette(nphase,2),m,n);
   rgb(:,:,3) =  black.*reshape(palette(nphase,3),m,n);
   
   rgb = BrightenRGB(rgb,.1);
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% phase plot with conformal polar grid for background

if  cs=='n'

   blackp = sawfct(myangle,1/phaseres,0.7,1);
   blackm = sawfct(log(abs(f)),2*pi/phaseres,0.7,1);
 
   black = blackp.*blackm;
   
   rgb(:,:,1) =  black.*reshape(palette(nphase,1),m,n);
   rgb(:,:,2) =  black.*reshape(palette(nphase,2),m,n);
   rgb(:,:,3) =  black.*reshape(palette(nphase,3),m,n);
  
   rgb = BrightenRGB(rgb,0.6);
   
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% proper phase portrait

if cs=='p'  
   
   rgb(:,:,1) =  reshape(palette(nphase,1),m,n);
   rgb(:,:,2) =  reshape(palette(nphase,2),m,n);
   rgb(:,:,3) =  reshape(palette(nphase,3),m,n);
      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  phase plot colored in steps

if cs=='q'  
   
   rgb(:,:,1) =  reshape(palette(nphase,1),m,n);
   rgb(:,:,2) =  reshape(palette(nphase,2),m,n);
   rgb(:,:,3) =  reshape(palette(nphase,3),m,n);
      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% special conformal grid 

if  cs=='s' 

   phase  = sawfct(myangle,1/phaseres,0,1);
   modul  = sawfct(log(abs(f)),2*pi/phaseres,0,1);
   
   rgb(:,:,1) =  phase.*(1-modul).*ones(m,n);
   rgb(:,:,2) =  modul.*(1-phase).*ones(m,n);
   rgb(:,:,3) =  phase.*modul.*ones(m,n);
   
   rgb = BrightenRGB(rgb,.1);
   
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% polar chessboard - light gray

if  cs=='t' 
   
   blackp = rectfct(myangle,1/phaseres);
   blackm = rectfct(log(abs(f)),2*pi/phaseres);
 
   black = mod(blackp+blackm,2);
  
   rgb(:,:,1) =  0.85*ones(m,n)+0.15*black;
   rgb(:,:,2) =  rgb(:,:,1);
   rgb(:,:,3) =  rgb(:,:,1);
  
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% polar chessboard

if  cs=='u' 

   blackp = rectfct(myangle,1/phaseres);
   blackm = rectfct(log(abs(f)),2*pi/phaseres);
 
   black = mod(blackp+blackm,2);
   
   rgb(:,:,1) =  black.*ones(m,n);
   rgb(:,:,2) =  black.*ones(m,n);
   rgb(:,:,3) =  black.*ones(m,n);
   
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% cartesian chessboard

if  cs=='v' 

   sc = 4; 
    
   blackx = rectfct(real(f),sc/phaseres);
   blacky = rectfct(imag(f),sc/phaseres);
 
   white = mod(blackx+blacky,2);
   
   rgb(:,:,1) =  white.*ones(m,n);
   rgb(:,:,2) =  white.*ones(m,n);
   rgb(:,:,3) =  white.*ones(m,n);
   
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% cartesian chessboard, light gray

if  cs=='w' 
    
   sc = 4; 

   blackx = rectfct(real(f),sc/phaseres);
   blacky = rectfct(imag(f),sc/phaseres);
 
   black = mod(blackx+blacky,2);
   
   rgb(:,:,1) =  0.85*ones(m,n)+0.15*black;
   rgb(:,:,2) =  rgb(:,:,1);
   rgb(:,:,3) =  rgb(:,:,1);
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% alternating black and white stripes corresponding to real part

if cs == 'x'  
    
   reres = 10/phaseres;
   
   black = sawfct(real(f),reres,0,1);
   
   bmax = max(max(black));
   bmin = min(min(black));
   black = floor(2*(black-bmin)/(bmax-bmin));
   
   rgb(:,:,1) =  black.*ones(m,n);
   rgb(:,:,2) =  black.*ones(m,n);
   rgb(:,:,3) =  black.*ones(m,n);
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% alternating black and white stripes corresponding to imaginary part

if cs == 'y' 
    
   
   imres = 10/phaseres;
   
   black = sawfct(imag(f),imres,0,1);
   
   bmax = max(max(black));
   bmin = min(min(black));
   black = floor(2*(black-bmin)/(bmax-bmin));
   
   rgb(:,:,1) =  black.*ones(m,n);
   rgb(:,:,2) =  black.*ones(m,n);
   rgb(:,:,3) =  black.*ones(m,n);
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% coloring zeros black, poles white %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

rgb(:,:,1) =  rgb(:,:,1) .* (abs(f)>0) + (1-rgb(:,:,1)).* (f==Inf);
rgb(:,:,2) =  rgb(:,:,2) .* (abs(f)>0) + (1-rgb(:,:,2)).* (f==Inf);
rgb(:,:,3) =  rgb(:,:,3) .* (abs(f)>0) + (1-rgb(:,:,3)).* (f==Inf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% auxiliary function: stepfct 

function y = stepfct(x,nmax)

%y = stepfct(x,nmax)
% integer step function of x with period 1 such that [0,1] --> [1,nmax]

x = x-floor(x);
y = floor(nmax*x)+1;

y = int16(y);
y = y+int16(y==0);

end

%% auxiliary function: sawfct

function y = sawfct(x,dx,a,b)

%y = sawfct(x,dx,a,b)
% saw tooth function on R with period dx onto [a,b]

x = x/dx-floor(x/dx);
y = a + (b-a)*x;

end

%% auxiliary function: sawfctxk

function y = sawfctx(x,xk)

%y = sawfctxk(x,xk)
% saw tooth function on [0,1] with values in [0,1] and jumps at points xk

xk = sort(xk);
kk = length(xk);
xk = [xk(kk)-1,xk,xk(1)+1];

y = zeros(size(x));

for j = 1 : kk+1
  y = y + ((x>=xk(j)).*(x<xk(j+1))) .* (x-xk(j)) ./ (xk(j+1)-xk(j));
end

y = y-floor(y);

end

%% auxiliary function: rectfct rectangular impulses of length dx

function y = rectfct(x,dx)
    
  y = mod(floor(x/dx),2);

end

%% auxiliary function: sawfctt

function  grayval = sawfctt(f,t,min,max)

%grayval = sawfctt(f,t,min,max)
% sawtooth function on unit circle with jumps min to max at points t
%
% f - function to be represented
% t - vector containing jump points on unit circle
% min (optional) min value of brightness
% max (optional) max value of brightness

if nargin<4, max = 1; end
if nargin<3, min = 0.8; end    

% eliminate multiple values of t

tol = 1e-4;
tau = sort(angle(t));
tt = find(abs(tau(2:end)-tau(1:end-1))>tol);
tau = [tau(1),tau(tt+1)];

t = exp(1i*tau);

ang    = (-angle(f)+pi)/(2*pi);
angjmp = (angle(conj(t))+pi)/(2*pi);

saw = sawfctx(ang,angjmp);

grayval = min+(max-min)*saw;

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RGB = BrightenRGB(RGB,bright)
% modification of color scheme
% bright - scalar value or field  
% between 0 and 1 for brightening
% between -1 and 0 for darkening 

if size(size(bright))==1
    bright = bright * ones(size(RGB(:,:,1)));
end
    
RGB(:,:,1) = (bright>=0).* ...
  ((1-bright).*RGB(:,:,1) + bright.*ones(size(RGB(:,:,1)))) ...
  + (bright<0).*((1+bright).*RGB(:,:,1));
 % 

RGB(:,:,2) = (bright>=0).* ...
    ((1-bright).*RGB(:,:,2) + bright.*ones(size(RGB(:,:,2)))) ...
   + (bright<0).*((1+bright).*RGB(:,:,2));
 %  

RGB(:,:,3) = (bright>=0).* ...
    ((1-bright).*RGB(:,:,3) + bright.*ones(size(RGB(:,:,3)))) ...
    + (bright<0).*((1+bright).*RGB(:,:,3));

end

end


