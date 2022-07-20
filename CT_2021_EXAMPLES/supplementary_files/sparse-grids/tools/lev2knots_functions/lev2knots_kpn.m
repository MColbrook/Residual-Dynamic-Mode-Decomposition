function nb_knots = lev2knots_kpn(I)

% nb_knots = lev2knots_kpn(I)
%
% returns the number of knots corresponding to the i-level i, 
% via the i2l map tabulated as in kpn_lev_table



%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



% are there non tabulated levels?
non_tab=I(I>5);

if isempty(non_tab)
    
    [r,c]=size(I);
    % nb knots is a vector after this instruction; I is understood coloumnwise as I(:)
    % Note that I want to access rows I+1, because minimum value for level is 0, 
    % whose data are stored in row 1
    vect_nb_knots=kpn_lev_table(I+1,3)'; 
    nb_knots=reshape(vect_nb_knots,r,c);
else
    %error('SparseGKit:OutOfTable',strcat('levels:',num2str(non_tab),' are not tabulated'))
    warning('SparseGKit:KpnNonTab','asking for non tabulated levels') 
    % put Inf for I>5
    [r,c]=size(I);
    pos=find(I>5); % pos is the linear position 
    I(I>5)=0; % to have it pass through the table. I will change it to Inf after
    vect_nb_knots=kpn_lev_table(I+1,3)';
    vect_nb_knots(pos)=Inf; %#ok<FNDSB>
    nb_knots=reshape(vect_nb_knots,r,c);    
end