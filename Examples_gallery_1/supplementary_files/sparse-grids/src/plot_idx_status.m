function plot_idx_status(G,I,idx_bin,idx)


% plot_idx_status(G,I,idx_bin,idx)
%
% plots status of a two-dimensional adaptive sparse grid

%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


% plot the indices used 
figure
plot(G(:,1),G(:,2),'xr','LineWidth',2,'MarkerSize',14,'DisplayName','G -- set used to build the sparse grid '),
hold on
plot(I(:,1),I(:,2),'o','MarkerFaceColor','b','MarkerSize',8,'DisplayName',sprintf('I -- set of idx whose neighbour \n has been explored ')), 
if ~isempty(idx_bin)
    plot(idx_bin(:,1),idx_bin(:,2),'sk','MarkerFaceColor','none','MarkerSize',14,...
        'DisplayName',sprintf('idx-bin -- set of idx whose \n neighbour has not been explored ')),
end
plot(idx(1),idx(2),'og','MarkerFaceColor','g','MarkerSize',8,'DisplayName','next idx to be considered'), 

set(gca,'FontSize',16)    
axis([0 12 0 12])
legend show
axis square
