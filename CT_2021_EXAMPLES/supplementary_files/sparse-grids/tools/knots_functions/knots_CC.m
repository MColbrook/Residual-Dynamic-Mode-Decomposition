function [x,w] = knots_CC(nn,x_a,x_b,whichrho)

% [x,w] = KNOTS_CC(nn,x_a,x_b)
% 
% calculates the collocation points (x) 
% and the weights (w) for the Clenshaw-Curtis integration formula
% w.r.t to the weight function rho(x)=1/(b-a) 
% i.e. the density of a uniform random variable 
% with range going from x=a to x=b. Note that for efficiency reasons
% nn must an odd number
% 
% [x,w] = KNOTS_CC(nn,x_a,x_b,'prob') 
% 
% is the same as [x,w] = KNOTS_CC(nn,x_a,x_b) above 
%
% [x,w]=[x,w] = KNOTS_CC(nn,x_a,x_b,'nonprob') 
%
% calculates the collocation points (x) 
% and the weights (w) for the Clenshaw-Curtis integration formula
% w.r.t to the weight function rho(x)=1


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



if nargin==3
    whichrho='prob';
end


if nn==1
    x=(x_a+x_b)/2; wt=1;

elseif mod(nn,2)==0
    error('SparseGKit:WrongInput','error in knots_CC: Clenshaw-Curtis formula \n use only odd number of points')

else
    n=nn-1;    
    
    N=(1:2:n-1)'; 
    end_N=length(N);   
    l=length(N);     
    m=n-l; 
    v0=[2./N./(N-2); 1/N(end_N); zeros(m,1)];
    end_v0=length(v0);
    v2=-v0(1:end_v0-1)-v0(end_v0:-1:2);
        
    g0=-ones(n,1); 
    g0(1+l)=g0(1+l)+n; 
    g0(1+m)=g0(1+m)+n;
    g=g0/(n^2-1+mod(n,2));
    
    
    wcc=ifft(v2+g);
    wt=[wcc;wcc(1,1)]'/2;
    
    x=cos( (0:n)*pi/n );
    x = (x_b-x_a)/2*x + (x_a+x_b)/2;
end

switch whichrho
    
    case 'nonprob'
        w=(x_b-x_a)*wt;

    case 'prob'
        w=wt;

    otherwise
    error('SparseGKit:WrongInput','4th input not recognized')
    
end
