function [res,evals] = quadrature_on_sparse_grid(f,S,Sr,evals_old,S_old,Sr_old,paral,tol)

% QUADRATURE_ON_SPARSE_GRID uses a sparse grid to compute the integral of a function.
% It provides evaluation recycling and parallel toolbox support. This function behaves as
% EVALUATE_ON_SPARSE_GRID, except that it return the value of the approximated integral 
% of the function. See EVALUATE_ON_SPARSE_GRID for more information on inputs. Possible calls:
%
% res = QUADRATURE_ON_SPARSE_GRID(F,SR)
%
%
% res = QUADRATURE_ON_SPARSE_GRID(F,S,SR,EVALS_OLD,S_OLD,SR_OLD) 
%
% res = QUADRATURE_ON_SPARSE_GRID(F,S,SR,[],[],[]) 
%
% res = QUADRATURE_ON_SPARSE_GRID(F,S,SR,EVALS_OLD,[],SR_OLD) 
%
%
% res = QUADRATURE_ON_SPARSE_GRID(F,S,SR,EVALS_OLD,S_OLD,SR_OLD,PARAL) 
%
%
% res = QUADRATURE_ON_SPARSE_GRID(F,S,SR,EVALS_OLD,S_OLD,SR_OLD,PARAL,TOL)
%
%
% [res,evals] = QUADRATURE_ON_SPARSE_GRID(...) returns the evaluations of the function F 
%       over the points of the sparse grid S


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

% declare a global variable controlling verbosity
global MATLAB_SPARSE_KIT_VERBOSE
if isempty(MATLAB_SPARSE_KIT_VERBOSE)
    MATLAB_SPARSE_KIT_VERBOSE=1;
end


switch nargin
    
    case 1
        error('SparseGKit:WrongInput','not enough input arguments')
        
    case 2
        % res = QUADRATURE_ON_SPARSE_GRID(f,S), S being a reduced sparse grid. 
        if ~isreduced(S)
            error('SparseGKit:WrongInput','when quadrature_on_sparse_grid is called with two inputs, the second one must be a reduced sparse grid')
        end
        evals = evaluate_on_sparse_grid(f,S);
        res = evals*S.weights';
        return
        
    case {3,4}
        errmsg = ['QUADRATURE_ON_SPARSE_GRID does not accept ',num2str(nargin),' inputs. ' ... 
            'Observe that QUADRATURE_ON_SPARSE_GRID has been changed after release 15.8 and ' ...
            'now takes as input also the non-reduced versions of the sparse grids on which ' ...
            'one wants to evaluate the function. This allows for significant savings in computational '...
            'time, especially for N large. See QUADRATURE_ON_SPARSE_GRID for more information and '...
            'QUADRATURE_ON_SPARSE_GRID_LEGACY if you are still using a version of the Sparse Grids Matlab kit older than 14.4'];
        error('SparseGKit:WrongInput',errmsg)
        
    case 5
        errmsg = ['QUADRATURE_ON_SPARSE_GRID does not accept ',num2str(nargin),' inputs. ' ... 
            'Observe that QUADRATURE_ON_SPARSE_GRID has been changed after release 15.8 and ' ...
            'now takes as input also the non-reduced versions of the sparse grids on which ' ...
            'one wants to evaluate the function. This allows for significant savings in computational '...
            'time, especially for N large. See QUADRATURE_ON_SPARSE_GRID for more information. '...
            'As a quick fix, you can use QUADRATURE_ON_SPARSE_GRID_LEGACY '...
            'which is the old version of QUADRATURE_ON_SPARSE_GRID, see help QUADRATURE_ON_SPARSE_GRID_LEGACY. '...
            'This function is however deprecated and will disappear from future relesases of the Sparse Grid Matlab Kit. '];
        error('SparseGKit:WrongInput',errmsg)
        
    case 6
        % in the previous versions, this was evals = quadrature_on_sparse_grid(f,Sr,evals_old,Sr_old,paral,tol);
        %
        % while now it is:
        %
        % quadrature_on_sparse_grid(f,S,Sr,evals_old,S_old,Sr_old)
        % or
        % quadrature_on_sparse_grid(f,S,Sr,[],[],[])
        % or  
        % quadrature_on_sparse_grid(f,S,Sr,evals_old,[],Sr_old)
        
        if ~issmolyak(S)
            errmsg = ['The second input of QUADRATURE_ON_SPARSE_GRID must be a non-reduced sparse grid if the function is called with 6 inputs. ' ...
                'Observe that QUADRATURE_ON_SPARSE_GRID has been changed after release 15.8 and ' ...
                'now takes as input also the non-reduced versions of the sparse grids on which ' ...
                'one wants to evaluate the function. This allows for significant savings in computational '...
                'time, especially for N large. See QUADRATURE_ON_SPARSE_GRID for more information. '...
                'As a quick fix, you can use QUADRATURE_ON_SPARSE_GRID_LEGACY '...
                'which is the old version of QUADRATURE_ON_SPARSE_GRID, see help QUADRATURE_ON_SPARSE_GRID_LEGACY. '...
                'This function is however deprecated and will disappear from future relesases of the Sparse Grid Matlab Kit.'];
            error('SparseGKit:WrongInput',errmsg)
        end
        evals = evaluate_on_sparse_grid(f,S,Sr,evals_old,S_old,Sr_old);
    case 7
        evals = evaluate_on_sparse_grid(f,S,Sr,evals_old,S_old,Sr_old,paral);
    case 8
        evals = evaluate_on_sparse_grid(f,S,Sr,evals_old,S_old,Sr_old,paral,tol);
end

res = evals*Sr.weights';