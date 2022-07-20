%----------------------------------------------------------------------------------
% Sparse Grids Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------------------------------------



%---------------------------
% Contributors
%---------------------------
%
% 1) Lorenzo Tamellini
% 2) Diane Guignard
% 3) Fabio Nobile
% 4) Giovanni Porta
% 5) Bjorn Sprungk
% 6) Francesco Tesei


%---------------------------
% 1) QUICK START
%
% Take a look at:
%
% -- RELEASE_NOTE.m for a detailed list of diffrences with the previous version
%
% -- SPARSE_GRIDS_TUTORIAL.m in the folder docs-examples for a gradual introduction to the code
%
% -- the slides SPARSE_GRIDS_MATLAB_KIT_TALK.pdf in the folder docs-examples also discuss basic and some advanced features of the code



%---------------------------
% 2) HOW TO INSTALL THE TOOLKIT
% 
% add to path this folder (with subfolders). Alternatively, run 
%
% addpath(genpath(pwd))
%
% at the beginning of the Matlab session. Consider launching Matlab as  
%
% matlab -r addpath(genpath(pwd))
%
% under Unix, it will run add_to_path during the start-up




%---------------------------
% 3) HOW TO USE THE TOOLKIT
%
% To create a sparse grids, use any of these two functions:
%
% -> smolyak_grid
% -> smolyak_grid_multiidx_set
%
% A number of functionalities to operate on sparse grids are provided: 
%
% -> adapt_sparse_grid
% -> convert_to_modal
% -> compute_sobol_indices_from_sparse_grid
% -> derive_sparse_grid
% -> export_sparse_grid_to_file
% -> evaluate_on_sparse_grid
% -> interpolate_on_sparse_grid
% -> plot_sparse_grid
% -> plot_sparse_grid_interpolant
% -> quadrature_on_sparse_grid
% -> reduce_sparse_grid
% 
% type HELP FUNCTIONAME to access help for each function.
% Documentation and examples can be found in docs-examples folder. 




%---------------------------
% 4) REPORT BUGS to tamellini@imati.cnr.it
% 
% also, let us know your email address if you want us to warn you whenever we release a new version of the toolbox




%---------------------------
% 5) PLEASE CITE US

% Please cite our toolbox by mentioning the webpage containing the package (http://csqi.epfl.ch) 
% and adding the following reference to your work:
%
% @InCollection{back.nobile.eal:comparison,
%  author    = {B\"ack, J. and Nobile, F. and Tamellini, L. and Tempone, R.},
%  title        = {Stochastic spectral {G}alerkin and collocation methods for {PDE}s with random coefficients: a numerical comparison},
%  booktitle    = {Spectral and High Order Methods for Partial Differential Equations},
%  pages        = {43--62},
%  publisher    = {Springer},
%  year        = 2011,
%  volume    = 76,
%  series    = {Lecture Notes in Computational Science and Engineering}, 
%  editor    = {Hesthaven, J.S. and Ronquist, E.M.},
%  note        = {Selected papers from the ICOSAHOM '09 conference, June 22-26, Trondheim, Norway}
%}

