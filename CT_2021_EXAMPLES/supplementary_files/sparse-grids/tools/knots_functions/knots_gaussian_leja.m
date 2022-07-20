function [x,w]=knots_gaussian_leja(n)

% [x,w] = knots_gaussian_leja(n)
%
% returns the collocation points (x) and the weights (w) 
% for the weighted Leja sequence for integration 
% w.r.t to the weight function 
%
% rho(x)=1/sqrt(2*pi)*exp(-x^2/2) 
%
% i.e. the density of a gaussian random variable 
% with mean 0 and standard deviation 1.
%
% Knots and weights have been precomputed (up to 50) following the work
% 
% A. Narayan, J. Jakeman, "Adaptive Leja sparse grid constructions for stochastic collocation and high-dimensional approximation"
% SIAM Journal on Scientific Computing,  Vol. 36, No. 6, pp. A2952–A2983, 2014
%
% an error is raised if more than 50 points are requested.
%
% knots are sorted increasingly before returning (weights are returned in the corresponding order)

%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile, B. Sprungk
% See LICENSE.txt for license
%----------------------------------------------------


if n>50
   error('SparseGKit:OutOfTable',strcat('this number of points is not available:',num2str(n)))

else
    
    % load points and weights
    SS = load('GaussianLejaPrecomputedKnotsAndWeights.mat');
    x = SS.X(1:n);
    w = SS.W(1:n,n);
    
    % sort knots increasingly and weights accordingly. Weights need to be row vectors
    [x,sorter]=sort(x);
    w=w(sorter)';

end






