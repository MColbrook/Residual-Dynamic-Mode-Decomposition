function h = plot_sparse_grids_interpolant(S,Sr,domain,f_values,varargin) %with_f_values,nb_plot_pts,couples)

% PLOTS_SPARSE_GRIDS_INTERPOLANT plots the sparse grid interpolant of a function. Different plots
% are produces depending on the number of dimensions of the sparse grid:
%
% if N==2, a surf plot will be generated
%
% if N==3, a number of contourfs (i.e. flat surfaces colored according to the value of the interpolant) 
%       will be stacked over the same axes 
%
% if N>3, a number of bidimensional cuts will be considered, and for each of them a surf will be generated.
%       In other words, all variables but two will be frozen to their average value and the resulting
%       two-dimensional plot will be produced
%
%
% PLOTS_SPARSE_GRIDS_INTERPOLANT(S,SR,DOMAIN,F_VALUES) produces the plots discussed above
%
% Additional inputs can be passed to control the behavior of the plots. Any combination of these optional
% inputs is allowed
% 
% PLOTS_SPARSE_GRIDS_INTERPOLANT(S,SR,DOMAIN,F_VALUES,'with_f_values') 
%       adds dots with the values of the sparse grids interpolant to the plots above (case N=2 and N>3).
%       For N==3, adds the sparse grid points in the 3D plot
%
% PLOTS_SPARSE_GRIDS_INTERPOLANT(S,SR,DOMAIN,F_VALUES,'nb_plot_pts',NP) 
%       sets the number of points used in each direction for the surf/contourf plots (default 20)
%
% PLOTS_SPARSE_GRIDS_INTERPOLANT(S,SR,DOMAIN,F_VALUES,'nb_contourfs',NC,'nb_countourf_lines',NL) 
%       sets the number of contourfs in the vertical direction for the case N=3 (default NC=5)
%       as well as the number of contourf lines NL (default 10)
%
% PLOTS_SPARSE_GRIDS_INTERPOLANT(S,SR,DOMAIN,F_VALUES,'two_dim_cuts',C) 
%      specifies the couples of variables to consider for the two-dimensional cuts when  N>3.
%      C is a vector with 2*k components denoting the directions of the cuts. 
%      For instance, the default value is C = [1 2 3 4 ...] and produces cut plots for (y1,y2)
%      (y3,y4),  (y5,y6)

%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



N = size(domain,2);


% look for 'with_f_values' in varargin
if any(cellfun(@(in) strcmp(in,'with_f_values'), varargin))
    with_f_values = true;
else
    with_f_values = false;
end

% value of 'plot_grid_size'
NP = value_of('nb_plot_pts',varargin);
if isempty(NP), NP = 20; end %#ok<*NASGU>

% value of 'nb_contourfs'
NC = value_of('nb_contourfs',varargin);
if isempty(NC) && N==3
        NC = 5; 
end
if ~isempty(NC) && N~=3
        warning('SparseGKit:IgnoredInput','ignoring nb_contourfs input')
end


% value of 'nb_contour_lines'
NL = value_of('nb_contourf_lines',varargin);
if isempty(NL) && N==3
        NL = 10; 
end
if ~isempty(NL) && N~=3
        warning('SparseGKit:IgnoredInput','ignoring nb_contourf_lines input')
end


% value of 'two_dim_cuts'
couples = value_of('two_dim_cuts',varargin);
if isempty(couples) && N > 3
    couples = 1:N;
end
if ~isempty(couples) && N <= 3
    warning('SparseGKit:IgnoredInput','ignoring two_dim_cuts input')
end 


% extract info on lower and upper ends of each direction
aa_vec = domain(1,:);
bb_vec = domain(2,:);
avg_vec = (aa_vec+bb_vec)/2;

% wrap interpolate on sparse grid into a @-function for ease of plotting
f_interp = @(x) interpolate_on_sparse_grid(S,Sr,f_values,x);



switch N
   
    case 2
    
        
        % generate a mesh grid over the cut
        xp = linspace(aa_vec(1),bb_vec(1),NP);
        yp = linspace(aa_vec(2),bb_vec(2),NP);
        
        [XP,YP] = meshgrid(xp,yp);
        
        nb_pts = length(xp)*length(yp);
        PTS = zeros(2,nb_pts);
        PTS(1,:)= XP(:)';
        PTS(2,:)= YP(:)';
        
        % interpolate on sparse grid
        f_interp_eval = f_interp(PTS);
        
        % reshape to use surf
        FIP = reshape(f_interp_eval,size(XP));
        
        h=figure;
        surf(XP,YP,FIP)
        xlabel('y_1')
        ylabel('y_2')

        if with_f_values
            hold on
            plot3(Sr.knots(1,:),Sr.knots(2,:),f_values,'ok','MarkerSize',12,'MarkerFaceColor','r')
        end
        
    case 3

        % generate a mesh grid over the cut
        xp = linspace(aa_vec(1),bb_vec(1),NP);
        yp = linspace(aa_vec(2),bb_vec(2),NP);
        zp = linspace(aa_vec(3),bb_vec(3),NC);

        
        [XP,YP] = meshgrid(xp,yp);
        
        XP_vect = XP(:);
        YP_vect = YP(:);
        
        PTS_XY = [XP_vect YP_vect]';
        nb_pts = size(PTS_XY,2);
        
        
        h=figure;
        
        for z_lev = 1:length(zp)
            
            PTS_Z = zp(z_lev)*ones(1,nb_pts);
            PTS= [PTS_XY; PTS_Z];
            f_interp_eval = f_interp(PTS);
            
            FIP = reshape(f_interp_eval,size(XP));
            
            [~,o2]=contourf(XP,YP,FIP,NL);
            o2.ContourZLevel = zp(z_lev);
            
            hold on
        end
                
        view([-30 20])
        xlabel('y_1')
        ylabel('y_2')
        zlabel('y_3')
        
        
        if with_f_values
            for pp = 1:size(Sr.knots,2)
                plot3(Sr.knots(1,pp),Sr.knots(2,pp),Sr.knots(3,pp),'ok','MarkerSize',8,'MarkerFaceColor','r','Color','k')
            end
        end
    
        colorbar
        
    otherwise
        
        CUTS = floor(length(couples)/2);
        
        h = zeros(1,CUTS);
        
        for ii = 1:CUTS
            
            couple_loc=couples([2*ii-1 2*ii]);
            v1 = couple_loc(1);
            v2 = couple_loc(2);
            
            % generate a mesh grid over the cut
            xp = linspace(aa_vec(v1),bb_vec(v1),NP);
            yp = linspace(aa_vec(v2),bb_vec(v2),NP);
            
            [XP,YP] = meshgrid(xp,yp);
            
            nb_pts = length(xp)*length(yp);
            
            % we need to generate the matrix of points where we want to evaluate our interpolant.
            % As usual, yt will be a fat matrix with points stored as columns. All directions
            % will be frozen to their average value, but the two of the local cut
            % (i.e., all rows but two will be constant).
            % We begin by making it all constant rows and then changing the rows we need
            %
            % PTS = [ avg_dir1; avg_dir2; avg_dir3 ...]
            PTS = diag(avg_vec)*ones(N,nb_pts);
            
            % replace lines of non-constant directions
            PTS(v1,:)= XP(:)';
            PTS(v2,:)= YP(:)';
            
            % interpolate on sparse grid
            f_interp_eval = f_interp(PTS);
            
            % reshape to use surf
            FIP = reshape(f_interp_eval,size(XP));
            
            h(ii)=figure;
            surf(XP,YP,FIP);
            if with_f_values
                hold on
                plot3(Sr.knots(v1,:),Sr.knots(v2,:),f_values,'ok','MarkerSize',16,'MarkerFaceColor','r')
            end
            title(['cut ',num2str(ii),' of ',num2str(CUTS),' over directions ',num2str(v1),' and ',num2str(v2)])
            xlabel(['y_',num2str(v1)])
            ylabel(['y_',num2str(v2)])
            
        end
        
end

end


function v =value_of(string,cell)

% if string is found in cell, return the value in the next cell

% logical array, 1 if string is found in cell
found = cellfun(@(in) strcmp(in,string), cell);

if any(found)
    % find the location of 1 in found
    pos = find(found);
    % the next input is our guy
    v = cell{pos+1};
else
    v= [];
end

end