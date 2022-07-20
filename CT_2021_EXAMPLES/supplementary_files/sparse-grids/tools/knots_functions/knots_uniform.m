function [x,w]=knots_uniform(n,x_a,x_b,whichrho)

% [x,w]=KNOTS_UNIFORM(n,x_a,x_b) 
%
% calculates the collocation points (x) 
% and the weights (w) for the gaussian integration 
% w.r.t. to the weight function rho(x)=1/(b-a) 
% i.e. the density of a uniform random variable 
% with range going from x=a to x=b.
%
%
% [x,w]=KNOTS_UNIFORM(n,x_a,x_b,'prob') 
% 
% is the same as [x,w]=KNOTS_UNIFORM(n,x_a,x_b) above 
%
%
% [x,w]=KNOTS_UNIFORM(n,x_a,x_b,'nonprob') 
%
% calculates the collocation points (x) 
% and the weights (w) for the gaussian integration 
% w.r.t to the weight function rho(x)=1


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


if nargin==3
    whichrho='prob';
end

if n==1

    x=(x_a+x_b)/2;
    wt=1;

else
    
    % calculates the values of the recursive relation
    [a,b]=coeflege(n);
    
    % builds the matrix
    JacM=diag(a)+diag(sqrt(b(2:n)),1)+diag(sqrt(b(2:n)),-1);
    
    % calculates points and weights from eigenvalues / eigenvectors of JacM
    [W,X]=eig(JacM);
    x=diag(X)';
    wt=W(1,:).^2;
    [x,ind]=sort(x);  %#ok<TRSRT>
    wt=wt(ind);
    
    % modifies points according to the distribution and its interval x_a, x_b
    x = (x_b-x_a)/2*x + (x_a+x_b)/2;
    
end

% finally, fix weights

switch whichrho
    
    case 'nonprob'
        w=(x_b-x_a)*wt;

    case 'prob'
        w=wt;

    otherwise
    error('SparseGKit:WrongInput','4th input not recognized')
    
end



%----------------------------------------------------------------------
function [a, b] = coeflege(n)

if (n <= 1), disp(' n must be > 1 '); 
    return; 
end

a=zeros(1,n);
b=zeros(1,n);

b(1)=2;

k=2:n;
b(k)=1./(4-1./(k-1).^2); 