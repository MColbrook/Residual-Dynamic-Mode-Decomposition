function [res,evals] = quadrature_on_sparse_grid_legacy(f,S,evals_old,S_old,paral,tol)

% QUADRATURE_ON_SPARSE_GRID_LEGACY uses a sparse grid to compute the integral of a function.
% It provides evaluation recycling and parallel toolbox support, see also EVALUATE_ON_SPARSE_GRID.
% This function is deprecated and will disappear from future relesases of the Sparse Grid Matlab Kit.
% Use QUADRATURE_ON_SPARSE_GRID instead.
%
% res = QUADRATURE_ON_SPARSE_GRID_LEGACY(F,S) computes the integral of F using the sparse grid S.
%       S is a reduced sparse grid; F is a function that takes as input a column vector 
%       (i.e. a sparse grid point) and returns either a scalar or a
%       column vector. F will be evaluated one point at a time so there's no need for F to accept as input 
%       matrices as well.
%
% res = QUADRATURE_ON_SPARSE_GRID_LEGACY(F,S,EVALS_OLD,S_OLD) recycles available evaluations 
%       of F on a different sparse grid (S_OLD). EVALS_OLD is a matrix storing the evaluations 
%       of F on S_OLD, where each evaluation stored as a column vector; S_OLD is a reduced sparse_grid 
%
% res = QUADRATURE_ON_SPARSE_GRID_LEGACY(F,S,[],[]) is equivalent to res = QUADRATURE_ON_SPARSE_GRID(F,S)
%
% res = QUADRATURE_ON_SPARSE_GRID_LEGACY(F,S,EVALS_OLD,S_OLD,PARAL) specifies the parameter
%       controlling the use of the parallel computing toolbox, see EVALUATE_ON_SPARSE_GRID. 
%
% res = QUADRATURE_ON_SPARSE_GRID_LEGACY(F,S,EVALS_OLD,S_OLD,PARAL,TOL), specifies the tolerance
%       used when checking equality between points of S and S_OLD. 
%
% [res,evals] = QUADRATURE_ON_SPARSE_GRID_LEGACY(...) returns the evaluations of the function F 
%       over the points of the sparse grid S


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

warning('SparseGKit:deprecated','quadrature_on_sparse_grid_legacy is a deprecated function')


switch nargin
    case 1
        error('SparseGKit:WrongInput','not enough input arguments')
    case 2
        % res = QUADRATURE_ON_SPARSE_GRID(f,S)
        evals = evaluate_on_sparse_grid_legacy(f,S);
    case 3
        % this was once this call, no longer valid
        % res = QUADRATURE_ON_SPARSE_GRID(f,S,map)
        errmsg=[' Starting from release 14.4 QUADRATURE_ON_SPARSE_GRID does not accept 3 input arguments. '...
                'In particular, QUADRATURE_ON_SPARSE_GRID does no longer accept unput argument MAP. '...
                'To fix this, either regenerate your sparse grid using MAP as input (help SMOLYAK_GRID) '...
                'or modify the field KNOTS of your reduced sparse grid applying MAP '...
                '(e.g. Sr.knots = map(Sr.knots) ). '...
                'This message will disappear in future releases of SPARSE_GRID_MATLAB_KIT'];
        error('SparseGKit:WrongInput',errmsg);
    case 4
        % this is either
        % % res = QUADRATURE_ON_SPARSE_GRID(f,S,map,weights_fact) or QUADRATURE_ON_SPARSE_GRID(f,S,[],weights_fact)
        % which are no longer valid, or
        % res = QUADRATURE_ON_SPARSE_GRID(f,S,evals_old,S_old)
        % or
        % res = QUADRATURE_ON_SPARSE_GRID(f,S,[],[])
        if ( isa(evals_old,'function_handle') || isempty(evals_old)) && ~isempty(S_old) 
            errmsg=['Unrecognized inputs. Note that starting from release 14.4 '...
                'QUADRATURE_ON_SPARSE_GRID does no longer accept input arguments MAP and WEIGHTS_FACT. '...
                'To fix this, either regenerate your sparse grid using MAP and WEIGHTS_FACT as input (help SMOLYAK_GRID) '...
                'or modify the fields KNOTS and WEIGHTS of your reduced sparse grid applying MAP and WEIGHTS_FACT'...
                '(e.g. Sr.knots = map(Sr.knots),  Sr.weights=weights_fact*Sr.weights). '...
                'This message will disappear in future releases of SPARSE_GRID_MATLAB_KIT'];
            error('SparseGKit:WrongInput',errmsg);
        end
        
        evals = evaluate_on_sparse_grid_legacy(f,S,evals_old,S_old);
        
    case 5
        evals = evaluate_on_sparse_grid_legacy(f,S,evals_old,S_old,paral);
    case 6
        evals = evaluate_on_sparse_grid_legacy(f,S,evals_old,S_old,paral,tol);

end

res = evals*S.weights';