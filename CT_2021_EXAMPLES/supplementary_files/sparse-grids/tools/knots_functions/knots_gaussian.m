function [x,w]=knots_gaussian(n,mi,sigma)

% [x,w]=KNOTS_GAUSSIAN(n,mi,sigma) 
% 
% calculates the collocation points (x) 
% and the weights (w) for the gaussian integration 
% w.r.t to the weight function 
% rho(x)=1/sqrt(2*pi*sigma) *exp( -(x-mi)^2 / (2*sigma^2) ) 
% i.e. the density of a gaussian random variable 
% with mean mi and standard deviation sigma


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



if (n==1) 
      % the point (traslated if needed) 
      x=mi;
      % the weight is 1:
      w=1;
      return
end

% calculates the values of the recursive relation
[a,b]=coefherm(n); 

% builds the matrix
JacM=diag(a)+diag(sqrt(b(2:n)),1)+diag(sqrt(b(2:n)),-1);

% calculates points and weights from eigenvalues / eigenvectors of JacM
[W,X]=eig(JacM); 
x=diag(X)'; 
w=W(1,:).^2;
[x,ind]=sort(x);  %#ok<UDIM>
w=w(ind);

% modifies points according to mi, sigma (the weigths are unaffected)
x=mi + sqrt(2)*sigma*x;



%----------------------------------------------------------------------
function [a, b] = coefherm(n)

if (n <= 1),  
    disp(' n must be > 1 '); 
    return; 
end
a=zeros(n,1); 
b=zeros(n,1); 

b(1)=sqrt(pi);
k=2:n;
b(k)=0.5*(k-1); 
