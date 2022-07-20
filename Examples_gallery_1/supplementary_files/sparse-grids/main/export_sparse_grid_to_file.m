function export_sparse_grid_to_file(Sr,filename,with_weights)

% EXPORT_SPARSE_GRID_TO_FILE saves knots of a reduced sparse grid to an ASCII file.
% The first line of the file shows the number of points in the grid and their dimension;
% then, points are stored as lines.
%
% Therefore, a sparse grid with 50 points in 2 dimensions will be stored in a file that
% looks as follows
% 
% 50 2
% coord1(P1) coord2(P1)
% coord1(P2) coord2(P2)
% coord1(P3) coord2(P3)
% ...
%
% EXPORT_SPARSE_GRID_TO_FILE(SR) saves the points in Sr in a file called POINTS.DAT
%
% EXPORT_SPARSE_GRID_TO_FILE(SR,FILENAME) saves the points in Sr in a file called FILENAME 
%       (extension should be provided in FILENAME)
%
% EXPORT_SPARSE_GRID_TO_FILE(SR,FILENAME,'with_weights') saves the points in Sr in a file called FILENAME 
%       (extension should be provided in FILENAME) and adds also the corresponding weight as
%       last entry of the row,  ie.
% 
% 50 2
% coord1(P1) coord2(P1) weight(P1)
% coord1(P2) coord2(P2) weight(P2)
% coord1(P3) coord2(P3) weight(P3)
% ...

%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------




switch nargin
    case 1
        filename = 'points.dat';
        with_weights = false;
    case 2
        with_weights = false;
    case 3
        if strcmp(with_weights,'with_weights')
            with_weights=true;
        else
            error('unknown second parameter')
        end
end

% first row, number of points and dimension
fid = fopen(filename,'w');
fprintf(fid,'%d %d\n',fliplr(size(Sr.knots)));
fclose(fid);

% then add all points, one per row,  and possibly weights if requested
if with_weights
    out = [Sr.knots' Sr.weights']; %#ok<*NASGU>
    save(filename,'out','-ascii', '-double','-append')    
else 
    out = Sr.knots';
    save(filename,'out','-ascii', '-double','-append')
end