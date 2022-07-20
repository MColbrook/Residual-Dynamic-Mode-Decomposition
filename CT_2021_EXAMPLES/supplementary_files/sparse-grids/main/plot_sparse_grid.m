function h= plot_sparse_grid(S,dims,varargin)

% h = plot_sparse_grid(S,dims,varargin)
%
% PLOT_SPARSE_GRID(S) plots S, which is a sparse grid in 2D. S can be either reduced or not. S can also be a tensor grid
%       (use istensor(S) to verify, type HELP ISTENSOR for details)
%
% PLOT_SPARSE_GRID(S,[d1 d2]) plots the d1 and d2 components of the points in S if S is more than 2D
%
% PLOT_SPARSE_GRID(S,[d1,d2],varargin) defines the style of the plot, e.g.
%
% plot_sparse_grid(S,[],'color','g','marker','o')
%
% h=PLOT_SPARSE_GRID(S)
%
% returns the handle to the plot


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



if ~exist('dims','var') || isempty(dims) 
    dims=[1 2];
end

x=[S.knots];

if nargin==1 || nargin==2
    % use a default plot style
    h=plot(x(dims(1),:),x(dims(2),:),'ok','MarkerFaceColor','k');
else
    % use style provided
    h=plot(x(dims(1),:),x(dims(2),:),varargin{:},'LineStyle','none');
end
grid on


