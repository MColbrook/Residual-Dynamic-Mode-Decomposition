function adapted = adapt_sparse_grid(f,N_full,knots,lev2knots,prev_adapt,controls)

% ADAPT_SPARSE_GRID computes an adaptive sparse grid approximation of a function f:
%
% ADAPTED = ADAPT_SPARSE_GRID(F,N_FULL,KNOTS,LEV2KNOTS,PREV_ADAPT,CONTROLS)
%
% computes a sparse grid adapted to the approximation of a function F: \Gamma subset R^N_full \to R^V 
% according to some profit indicator. The implementation closely resembles the algorithm proposed in
%
% T. Gerstner, M. Griebel, Dimension–Adaptive Tensor–Product Quadrature, Computing 71, 65–87 (2003)
%
% and the various definitions of profits are discussed in 
%
% F. Nobile, L. Tamellini, F. Tesei, and R. Tempone, An Adaptive Sparse Grid Algorithm for Elliptic PDEs 
% with Lognormal Diffusion Coefficient, Sparse Grids and Applications – Stuttgart 2014,
% Lecture Notes in Computational Science and Engineering 109, p. 191-220, 
% Springer International Publishing Switzerland 2016
%
% Roughly speaking, at each iteration the algorithm considers the multiidx with the highest profit
% among those already explored, adds it to the sparse grid, computes the profits of all its neighbouring 
% indices (see controls.prof) and adds them to the list of explored indices.
% 
%
% --> F is a function-handle, taking as input a column vector (N_full components) and returning a column-vector value (V components). 
% 
% --> N_full is a scalar, denotes the dimension of the space \Gamma over which F and the sparse grid are defined.
%       Note that this is the full dimension, yet the algorithm might explore a smaller subset of variables to
%       start with, N_curr, as specified in CONTROLS.var_buffer_size (see below)
%
% --> KNOTS, LEV2KNOTS are the quantities defining the knots to be used in each direction (same family of
%       knots assumed in each direction of the space Gamma)
%
% --> PREV_ADAPT: by setting PREV_ADAPT as the output ADAPTED of a previous computation, the new computation
%       will resume from where the previous one left. Set PREV_ADAPT = [] for a fresh start
%
% --> CONTROLS is a struct defining several parameters that control the algorithm flow. It is possible to
%       define only some (or none) of these parameters, the others will be set to default value. Only
%       the parameter 'nested' is mandatory. The parameter 'pdf' is mandatory only if controls.prof is
%       set to any of the 'weighted*' choices (see below). For safety, 'pdf' is not set to default and
%       an error will be raised if the field is not set
%
%       controls.nested   : (mandatory) true if nested points will be used, false otherwise
%
%       controls.max_pts  : is the maximum number of points the adapted sparse grid can have (default 1000).
%                           Observe that if the max nb of points is reached while checking the neighbours 
%                           of the current idx, the algorithm will not be interrupted.
%
%       controls.prof_tol : the algorithm will stop as soon as profit of the considered idx is smaller
%                           then prof_tol (default 1e-14)
%
%       controls.paral    : parallel computation is used as soon as at least <paral> function evaluation are
%                           needed
%
%       controls.recycling: can be set to 'priority_to_evaluation' (default) or 'priority_to_recycling'. This affects
%                           only when non-nested points are used. Because in this case points enter and exit
%                           from the sparse grid, we should keep a record of all points that have been
%                           evaluated and recycle from it whenever possible. However, doing this exactly is **very
%                           expensive** in terms of computational cost with the current data structure if N is
%                           large. So, unless your function evaluation is **very** expensive, we recommend to
%                           leave the default setting, which does a "partial search", i.e. searches for
%                           recycling only in the previous grid (instead than on the entire sparse grid
%                           history) and therefore implies that the same point might be evaluated multiple times.
%
%       controls.prof     : chooses how to compute the profit of an idx. In general, this entails comparing
%                           the sparse grids approximation of f *before* adding idx to the sparse grids (call 
%                           this approximation S) and *after* having added it (call this approximation T).
%                           Note that S, T are actually functions of y in \Gamma and return a vector with V components
%                           (same as f, the function to be approximated). In general, profit takes either of the two forms
%                                 
%                                 P(idx) = Op_y ( Op_vect(S,T) ) ,  P(idx) = Op_vect ( Op_y(S,T) )
%             
%                           where Op_vect acts as a "vector norm" over the V components of S,T (e.g., the euclidean norm for
%                           fixed y), and Op_y acts as a norm over Gamma (e.g. the maximum over a discrete set of values of y).
%                           Currently, the choices below are available. In all of these choices, Op_vect can
%                           be changed by setting the field controls.op_vect (see below)
%
%                               'Linf'                      :   P(idx) = Op_y( Op_vect(S(y),T(y) ) where
%                                                               Op_vect(S(y),T(y)) = euclidean norm of S(y) - T(y) for a fixed y; 
%                                                               Op_y = max over a discrete set of points of Gamma.
%                                                               For nested points, the algorithm considers the points that would
%                                                               be added to the grid by the considered idx,
%                                                               evalutes the function and the previous sparse grid on such
%                                                               points, and takes the max of such difference
%                                                               as profit estimate. For non-nested points, the sparse grid approx
%                                                               is not interpolant, hence we consider the difference between both the
%                                                               *new* and the previous sparse grids on the new points.
%
%                               'Linf/new_points' (default) :   the 'Linf' profit estimate above is further divided
%                                                               by the number of new points considered 
%                                                               (i.e. gain-cost ratio). For non-nested points, this is actually
%                                                               the cost of the tensor grid associated to the current midx (see
%                                                               Nobile Tamellini Tempone, ``Convergence of quasi optimal sparse
%                                                               grid approximation ...''
%
%                               'weighted Linf'             :   same as 'Linf', but the difference between
%                                                               previous sparse grid and function evaluation (or between new and old
%                                                               sparse grids, depending whether nested or non-nested points are used) is 
%                                                               multiplied by the probability density function (see below for details). 
%                                                               Useful when considering sparse grids on unbounded domains.
%
%                               'weighted Linf/new_points'  :   the 'weighted Linf' profit estimate is divided by the number of new points
%
%                               'deltaint'                  :   P(idx) = Op_vect ( Op_y(S,T) ) = euclidean_norm( expct(S) - expct(T))
%                                                               i.e. the profit estimate is the euclidean norm of the difference of the
%                                                               quadratures using the sparse grids with and without the given idx
%
%                               'deltaint/new_points'       :   the 'deltaint' profit estimate is divided by the number of new points            
%
%
%       controls.op_vect   : changes the operator op_vect above. It is defined as a handle function
%                            @(A,B) some_function_of_(A,B)   
%                            where A B are two matrices with V rows, containing the evaluation of the two operators S T discussed above 
%                            as columns, e.g. A=[S(y1) S(y2) S(y3) ...], B=[T(y1) T(y2) T(y3) ...]
%
%       controls.pts_tol   : is the tolerance used by evaluate_on_sparse_grid to check for equal points (default 1e-14)
%
%		controls.pdf       : pdf over Gamma \subset R^N_full, the set over which the sparse 
%                            of the random variables, to be used in weighted* profits. It must be provided
%                            as a function handle that takes as input a matrix with N rows (for generic N), i.e. where each column is a point
%                            where the pdf must be evaluated. For instance, if the weight is the standard gaussian
%                            pdf on each direction, controls.pdf = @(Y) prod(normpdf(Y,0,1),1), with Y matrix with N
%                            rows. For safety, 'pdf' is not set to default and an error will be raised if the field is not set
%
%       controls.var_buffer_size : the algorithm starts exploring N_curr = <var_buffer_size> dimensions. As soon as
%                           points are placed in one dimension, a new dimension is added to the set of
%                           explored variables, i.e., N_curr = N_curr+1. In this way we ensure that there are always <var_buffer_size>
%                           explored but "non-activated" variables, i.e., along which no point is placed (default min(N_full,5) )                         
%
%
%       controls.plot      : plot multiidx set and pause (default false)  
%
%       
%
%
% --> ADAPTED is a struct containing the results of the algorithms. The same structure can be passed
%     as input (PREV_ADAPT) to resume the computation. Its fields are 
%
%       adapted.S       : the adapted sparse grid (built over both all indices whose profit has been computed,
%                           even those whose profit hasn't yet been chosen as the best one, i.e. those whose
%                           neighbours haven't been explored yet).
%
%       adapted.Sr      : its reduced version;
%
%       adapted.f_on_Sr : the evaluations of f over the list of points contained in Sr;
%
%       adapted.nb_pts  : the number of points in Sr;
%
%       adapted.nb_pts_visited: the number of points visited while building the sparse grid. For non-nested
%                           points, this will be larger than nb_pts, because some points enter and then exit
%                           the grid when the corresponding idx exits from the combination technique.
%
%       adapted.num_evals: the number of function evaluations done to compute the sparse grid. This is not
%                           necessarily equal to the previous one, because for speed reasons (looking for
%                           points in expensive for N large) sometimes one point might be recomputed twice
%
%       adapted.N       : the current number of dimensions considered for exploration
%
%       adapted.private : a structure contained more detailed information on the status of the adaptive algorithm, that is needed
%                         to resume the computation.  In particular, the needed data structure consists of:
%
%                           private.G       :   the multi-idx set used to build the sparse grid. This is called
%                                               \mathcal{I} in Gerstner-Griebel paper.
%
%                           private.I       :   the set of explored indices; G = I plus the neighbours of I.
%                                               This is called \mathcal{O} in Gerstner-Griebel paper. I is sorted
%                                               lexicographically
%
%                           private.I_log   :   same as private.I, but sorted in the order with which indices
%                                               are added to I insted of lexicographic
%
%                           private.N_log   :   for each iteration, the value of N_curr
%
%                           private.idx=idx :   the idx whose neighbour is the next to be explored 
%
%                           private.maxprof :   the corresponding profit
%
%                           private.profits :   the set of idx whose profit has been computed. They have been added to the grid 
%                                               but their neighbour is yet to be explored.
%
%                           private.idx_bin :   the corresponding set of profits. This is called A in
%                                               Gerstner-Griebel paper
%                           private.var_with_pts: vector of variables in which we have actually put points. 
%                                               length(var_with_pts) + controls.var_buffer_size = N_curr
%                       
%                           private.nb_pts_log : for each iteration, the current nb_pts
%
%                           private.Hr      :   for non-nested points, the list of points visited by the
%                                               algorithm. Empty for nested-points
%
%                           private.f_on_Hr :   for non-nested points, the evaluations of f on Hr. Empty for nested-points

%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile, D. Guignard, F. Tesei, B. Sprungk
% See LICENSE.txt for license
%----------------------------------------------------


% declare a global variable controlling verbosity
global MATLAB_SPARSE_KIT_VERBOSE
if isempty(MATLAB_SPARSE_KIT_VERBOSE)
    MATLAB_SPARSE_KIT_VERBOSE=1;
end



% set control fields
controls = default_controls(controls,N_full);

% init / resume adaptive algorithm.
%
% Here's the data structure we need
%
% --> idx       : is the idx with the highest profit, whose neighbour we will next explore
% --> maxprof   : the corresponding profit
% --> I         : is the set of explored indices (the grid is actually larger, it includes as well their
%                   neighbours), corresponds to O in Gerstner-Griebel
% --> I_log     : I must be sorted lexicographically for software reasons. I_log is the same set of indices, but
%                   sorted in the order in which they are chosen by the algorithm
% --> idx_bin   : is the set of idx whose profit has been computed. They have been added to the grid but their 
%                 neighbour is yet to be explored. Corresponds to A in Gerstner-Griebel
% --> profits   : is the corresponding set of profits
% --> G         : is the set of the grid. Corresponds to I in Gerstner-Griebel
% --> nb_pts    : the number of points in the grid
% --> nb_pts_visited : the number of points visited during the sparse grid construction
% --> num_evals : the number of function evaluations performed
% --> nb_pts_log: for each iteration, the current nb_pts
% --> S         : smolyak_grid_multiidx_set(G,knots,lev2knots);
% --> Sr        : reduce_sparse_grid(S);
% --> f_on_Sr   : evaluate_on_sparse_grid(f,Sr);
% --> var_with_pts: vector of variables in which we have actually put points. 
%                   length(var_with_pts) + controls.var_buffer_size = N_curr
% --> N_log     : for every iteration, the value of N_curr

% we need to distinguish the full dimension of the parameter space (N_tot) and the dimensional of the subspace
% currently explored (N_curr = length(expl_var)). For consistency with previous code, we actually use N instead of N_curr

% this is the current number of dimensions activated. Its value might change at the next line if we are
% actually resuming a previous computation. Also, we make sure that this is actually no larger than N_full
% already!

N = controls.var_buffer_size;


[N,N_log,var_with_pts,S,Sr,f_on_Sr,I,I_log,idx,maxprof,idx_bin,profits,G,G_log,Hr,f_on_Hr,...
    nb_pts,nb_pts_log,num_evals,intf] = start_adapt(f,N,knots,lev2knots,prev_adapt); 


% here's the adapt algo


while nb_pts < controls.max_pts   %while nb_pts_wrong_count < controls.max_pts

    if maxprof < controls.prof_tol
        break
    end
    
    %-------------------------------------------
    % compute neighbours of the idx with highest profit
    %-------------------------------------------
    
    
    Ng=ones(N,1)*idx + eye(N);  
    
    %-------------------------------------------
    % remove those unadmissible
    %-------------------------------------------

    m=1;
    while m<=size(Ng,1)
        if check_index_admissibility(Ng(m,:),I)
            m=m+1;
        else
            Ng(m,:)=[];
        end
    end
    
    %-------------------------------------------
    % for those safe, add them to the grid & compute profit
    %-------------------------------------------
    
    M=size(Ng,1);                                               
    
    Prof_temp=zeros(1,M);
    
    
    for m=1:M
        
        G = sortrows([G; Ng(m,:)]);
        G_log = [G_log; Ng(m,:)]; %#ok<AGROW>
        %T = smolyak_grid_multiidx_set(G,knots,lev2knots);
        T = smolyak_grid_multiidx_set(G,knots,lev2knots,S); % recycle tensor grids from previous sparse
        Tr = reduce_sparse_grid(T);

        [nb_pts,num_evals,nb_pts_log,Prof_temp(m),f_on_Tr,Hr,f_on_Hr,intnew] = ...
            compute_profit_idx(Ng(m,:),f,S,T,Tr,Sr,Hr,f_on_Sr,f_on_Hr,intf,nb_pts,num_evals,nb_pts_log,knots,lev2knots,controls);
                
        S=T;
        Sr = Tr;
        f_on_Sr = f_on_Tr;
        intf = intnew;
        
    end
    
    %-------------------------------------------
    % update the list of profits
    %-------------------------------------------

    profits= [profits, Prof_temp]; %#ok<AGROW>
    
    % observe that any of the indices in Ng can be already in idx_bin, i.e. there are no duplicates. Indeed, if ii is in Ng, this
    % means that ii is admissible wrt I, therefore all previous indices jj<=ii are already included in I, therefore they have already
    % left idx_bin therefore no jj will be considered and ii will no longer be generated in Ng. 
    % Anyway, let's leave a check for the moment
    
    idx_bin = [idx_bin; Ng];  %#ok<AGROW>
    if size(unique(idx_bin,'rows'),1)~=size(idx_bin,1)
        error('an index is appearing twice in idx_bin')
    end
    
    
    %-------------------------------------------
    % take the next idx with highest profit. I'd like to remove it but first I need to add (possibly) the new dimension.
    % Note that idx it's already in G
    %-------------------------------------------
            
    [maxprof,k]=max(profits);
    idx=idx_bin(k,:);

    %-------------------------------------------
    % now we take care of possibly changing the variables buffer, i.e., adding a new variable to be explored
    %-------------------------------------------
    
    % find the list of variables in which the current choice of idx wants to add points    
    to_be_explored=find(idx>1);
    
    new_var = setdiff(to_be_explored,var_with_pts); % i.e., the variables activated by the new profit in which 
                                                % no points have still been placed. Note that obviously new_var <= N_full
    
    switch length(new_var)
        case 0
            % do nothing here, we keep working in the same subspace
            if MATLAB_SPARSE_KIT_VERBOSE,
                disp('keep number of dimensions as is')
            end
        case 1
            % in this case, we add the new proposed variable to the ones where we have added points and
            % we also add one variable to the explored one, to maintain the balance
            % length(var_with_pts) + controls.var_buffer_size = N_curr

            if MATLAB_SPARSE_KIT_VERBOSE,
                disp('adding points a new variable')
            end
            var_with_pts = union(var_with_pts,new_var);
            
            % the new variable to be explored is necessarily the N_curr+1, so adding it is just a matter
            % of adding a column to the proper containers, unless the hard limit of N_full, 
            % i.e. the total number of variables of f has been reached; then, there are no more variables to explore.
            % Observe though that even if all variables are explored, not all variables necessarily have points,
            % so the previous line is not to be put inside the following if
            
            if N < N_full            

                if MATLAB_SPARSE_KIT_VERBOSE,
                    disp('adding a new variable to the explored dimensions')
                end

                % let's add one variable and increase the containers. We first add one dimension to the index containers
                N = N+1;               
                I=[I,ones(size(I,1),1)]; %#ok<AGROW>
                I_log=[I_log,ones(size(I_log,1),1)]; %#ok<AGROW>
                G=[G,ones(size(G,1),1)]; %#ok<AGROW>
                G_log=[G_log,ones(size(G_log,1),1)]; %#ok<AGROW>
                idx_bin=[idx_bin,ones(size(idx_bin,1),1)]; %#ok<AGROW>
                
                % then we add one coordinate to the sparse grid points that we have already generated. The new
                % coordinate is the midpoint of the parameter space along the new direction. Here it's crucial that
                % f(x)==f([x mp]), otherwise I have to reevaluate everything, which I do not want to do. I'll just raise a 
                % warning for the time being
                if MATLAB_SPARSE_KIT_VERBOSE,
                    disp('adding a new variable, hence a new coordinate to points. Does this change f evaluations? If so, the code does not work because it does not recompute function evaluations')                     
                end
                
                mp = knots(1);
                for j=1:size(S,2)
                    S(j).knots=[S(j).knots;mp*ones(1,S(j).size)];
                    S(j).knots_per_dim{end+1}=mp;
                    S(j).m(end+1)=1;
                    S(j).idx = [S(j).idx 1];
                end
%                 if isfield(S(1),'knots_per_dim')
%                     for j=1:size(S,2)
%                         S(j).knots_per_dim{end+1}=mp;
%                         S(j).m(end+1)=1;
%                     end
%                 end
                Sr.knots=[Sr.knots;mp*ones(1,size(Sr.knots,2))];
                if ~controls.nested
                    Hr.knots=[Hr.knots;mp*ones(1,size(Hr.knots,2))];
                end
                
                % because we have added the new variable, we need to add right away the profit of the first index in 
                % its direction, as if we had started with Ng=ones(N+1,1)*idx + eye(N+1) in the first place. In this way, 
                % we trigger exploration in that direction too

                Ng=ones(1,N); Ng(end)=2;
                
                G = sortrows([G; Ng]);
                G_log = [G_log; Ng];  %#ok<AGROW>
                %T = smolyak_grid_multiidx_set(G,knots,lev2knots);
                T = smolyak_grid_multiidx_set(G,knots,lev2knots,S); % recycle tensor grids from previos sparse
                Tr = reduce_sparse_grid(T);

                % [nb_pts,nb_pts_log,nb_pts_wrong_count,Prof_temp,f_on_Tr,Hr,f_on_Hr,intnew] = ...
                %   compute_profit_idx(Ng,f,S,T,Tr,Sr,Hr,f_on_Sr,f_on_Hr,intf,nb_pts,nb_pts_log,knots,lev2knots,controls);
                [nb_pts,num_evals,nb_pts_log,Prof_temp,f_on_Tr,Hr,f_on_Hr,intnew] = ...
                    compute_profit_idx(Ng,f,S,T,Tr,Sr,Hr,f_on_Sr,f_on_Hr,intf,nb_pts,num_evals,nb_pts_log,knots,lev2knots,controls);
                    
                S=T;
                Sr = Tr;
                f_on_Sr = f_on_Tr;
                intf = intnew;
                
                profits= [profits, Prof_temp]; %#ok<AGROW>                
                
                idx_bin = [idx_bin; Ng];  %#ok<AGROW>
                if size(unique(idx_bin,'rows'),1)~=size(idx_bin,1)
                    error('an index is appearing twice in idx_bin')
                end

                % now, this new profit might already be the best one, so I need to reconsider my previous decision. 
                % I recompute the best profit and now I can also remove it from idx_bin (a few lines below)
                
                [maxprof,k]=max(profits);
                    
                idx=idx_bin(k,:);
                
            else
                if MATLAB_SPARSE_KIT_VERBOSE,
                    disp('maximum number of variables to be explored reached, continuing as is') 
                end
            end
            
        otherwise
            error('still don''t know how to deal with the case where 2 or more new variables need to be explored')
    
    end
    %-------------------------------------------
    % end of the section where we change the variables buffer
    %-------------------------------------------

    I_log = [I_log; idx]; %#ok<AGROW>
    I=sortrows([I; idx]); % note that I must be lexicog sorted for check_index_admissibility(Ng(m,:),I) to work

    idx_bin(k,:)=[];
    
    profits(k)=[];

    N_log(end+1)=N; %#ok<AGROW>
    
    if controls.plot
        plot_idx_status(G,I,idx_bin,idx)
        pause
    end
end

% done with the loop on idx, we can reduce Hr to squeeze out multiple occurrencies of same point (if any) and
% fix nb_pts_visited
if ~controls.nested && strcmp(controls.recycling,'priority_to_evaluation')
    % I need to make Hr look like a sequence of tensor grids. I actually only need to add to it a fake weights field
    Hr.weights = zeros(1,size(Hr.knots,2));
    Hr= reduce_sparse_grid(Hr);
    f_on_Hr = f_on_Hr(:,Hr.m); % here I remove duplicates in f_on_Hr too
end


adapted.N=N;
adapted.S=S;
adapted.Sr=Sr;
adapted.f_on_Sr=f_on_Sr;
adapted.nb_pts=nb_pts;
if controls.nested
    adapted.nb_pts_visited = nb_pts;
else
    adapted.nb_pts_visited = size(Hr.knots,2);
end
adapted.num_evals = num_evals;
if num_evals>adapted.nb_pts_visited, 
    disp([  'Some points have been evaluated more than once. Total: ',num2str(num_evals-adapted.nb_pts_visited),...
            ' extra evaluations over ',num2str(adapted.nb_pts_visited),' function evaluations']); 
end
adapted.intf=intf;


private.G=G;
private.G_log=G_log;
private.I=I;
private.I_log = I_log;
private.maxprof=maxprof;
private.idx=idx;
private.profits=profits;
private.idx_bin=idx_bin;
private.Hr=Hr;
private.f_on_Hr=f_on_Hr;
private.var_with_pts=var_with_pts;
private.N_log = N_log;
private.nb_pts_log=nb_pts_log;  

adapted.private = private;


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   --------------------------  end of main function -----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%
%
%---------------------------------------------------------------------------------
% auxiliary function, set default controls
%---------------------------------------------------------------------------------


function controls = default_controls(controls,N_full)


if ~isfield(controls,'pts_tol')
    controls.pts_tol=1e-14;
end
if ~isfield(controls,'max_pts')
    controls.max_pts=1000;
end
if ~isfield(controls,'prof_tol')
    controls.prof_tol=1e-14;
end
if ~isfield(controls,'paral')
    controls.paral=NaN;
end
if ~isfield(controls,'prof')
    controls.prof='Linf/new_points';
end
if ~isfield(controls,'plot')
    controls.plot=false;
end
if ~isfield(controls,'nested')
    error('controls must specify the value of ''nested'' field')
end
if ~isfield(controls,'op_vect')
    controls.op_vect = @(A,B) sqrt(sum((A - B).^2,1));
end
if strcmp(controls.prof,'weighted Linf/new_points') || strcmp(controls.prof,'weighted Linf')
    if ~isfield(controls,'pdf')
        error('you need to set the field ''pdf'' to use ''weighted Linf'' and ''weighted Linf/new_points'' profits')
    end
end
if ~isfield(controls,'var_buffer_size')
    controls.var_buffer_size = min(N_full,5);
elseif isfield(controls,'var_buffer_size') && controls.var_buffer_size > N_full
    controls.var_buffer_size = N_full;
    warning('SparseGKit:BuffGTNfull','controls.var_buffer_size cannot be greater than N_full. The code will proceed with controls.var_buffer_size = N_full;') 
    pause
end

if ~isfield(controls,'recycling')
    controls.recycling = 'priority_to_evaluation';
end
switch controls.recycling
    case {'priority_to_evaluation','priority_to_recycling'}
    otherwise
        error('unknown value of field controls.recycling')
end
end

%
%
%
%
%
%
%
%---------------------------------------------------------------------------------
% auxiliary function, init/resume algorithm
%---------------------------------------------------------------------------------


function [N,N_log,var_with_pts,S,Sr,f_on_Sr,I,I_log,idx,maxprof,idx_bin,profits,G,G_log,Hr,f_on_Hr,...
    nb_pts,nb_pts_log,num_evals,intf] = start_adapt(f,N,knots,lev2knots,prev_adapt)


% --> I         : is the set of explored indices (the grid is actually larger, it includes as well their neighbours)   
% --> I_log     : I must be sorted lexicographically for software reasons. I_log is the same set of indices, but
%                   sorted in the order in which they are chosen by the algorithm
% --> idx       : is the idx with the highest profit, whose neighbour we will next explore
% --> maxprof   : the corresponding profit
% --> idx_bin   : is the set of idx whose profit has been computed. They have been added to the grid but their neighbour is yet to be explored
% --> profits   : is the corresponding set of profits
% --> G         : is the set of the grid.
% --> nb_pts    : the number of points in the grid
% --> num_evals : the number of function evaluations
% --> nb_pts_log: for each iteration, the current nb_pts
% --> S         : smolyak_grid_multiidx_set(G,knots,lev2knots);
% --> Sr        : reduce_sparse_grid(S);
% --> f_on_Sr   : evaluate_on_sparse_grid(f,Sr);
% --> intf      : approx of integral of f using Sr
% --> Hr        : all the points visited by the algo, stored as a reduced grid to be able to use ; only useful for non-nested grids, where it differs from Sr.
% --> var_with_pts: vector of variables in which we have actually put points. 
%                   length(var_with_pts) + controls.var_buffer_size = N_curr
% --> N_log     : for every iteration, the value of N_curr

% declare a global variable controlling verbosity
global MATLAB_SPARSE_KIT_VERBOSE

if isempty(MATLAB_SPARSE_KIT_VERBOSE)
    MATLAB_SPARSE_KIT_VERBOSE=1;
end


if isempty(prev_adapt)

    % it's a fresh start
    %--------------------------------------------

    var_with_pts = []; % we have put no points in no variables for now
    N_log = N;
    I=ones(1,N);
    I_log = ones(1,N);
    idx = ones(1,N);
    maxprof = Inf;

    idx_bin=[];
    profits=[];
    
    G = I;    
    G_log = G;
    S  = smolyak_grid_multiidx_set(G,knots,lev2knots);
    Sr = reduce_sparse_grid(S);
    f_on_Sr = evaluate_on_sparse_grid(f,Sr);

    Hr = Sr;
    f_on_Hr = f_on_Sr; % it is a matrix of size VxM where M is the number of points and f:R^N->R^V

    intf=f_on_Sr*Sr.weights';
    nb_pts=size(f_on_Sr,2); 
    nb_pts_log = nb_pts;
    num_evals = nb_pts;
else

    % we are resuming from a previous run
    %--------------------------------------------

    if MATLAB_SPARSE_KIT_VERBOSE,
        disp('adapt--recycling')
    end

    N = prev_adapt.N;
    N_log = prev_adapt.private.N_log;
    var_with_pts = prev_adapt.private.var_with_pts;
    I=prev_adapt.private.I;
    I_log = prev_adapt.private.I_log;
    idx = prev_adapt.private.idx;
    maxprof = prev_adapt.private.maxprof;

    idx_bin=prev_adapt.private.idx_bin;
    profits=prev_adapt.private.profits;
    
    G = prev_adapt.private.G;    
    G_log = prev_adapt.private.G_log;    
    S = prev_adapt.S;
    Sr= prev_adapt.Sr;
    f_on_Sr = prev_adapt.f_on_Sr;
    
    Hr=prev_adapt.private.Hr;
    f_on_Hr = prev_adapt.private.f_on_Hr;
    
    intf = prev_adapt.intf;
    nb_pts=prev_adapt.nb_pts;  
    nb_pts_log=prev_adapt.private.nb_pts_log;
    num_evals = prev_adapt.num_evals;
    
end

end


%
%
%
%
%
%
%
%---------------------------------------------------------------------------------
% auxiliary function, compute profits
%---------------------------------------------------------------------------------

function [nb_pts,num_evals,nb_pts_log,Prof_temp,f_on_Tr,Hr,f_on_Hr,intnew] = ...
      compute_profit_idx(ng_idx,f,S,T,Tr,Sr,Hr,f_on_Sr,f_on_Hr,intf,nb_pts,num_evals,nb_pts_log,knots,lev2knots,controls)


N = size(ng_idx,2);

if controls.nested

    % here we we evaluate on new points only. Note that finding which points have been evaluated already
    % relies on multiindex info almost exclusively (because the points are nested) so this is quite efficient
    
    [f_on_Tr,new_points,idx_newp] = evaluate_on_sparse_grid(f,T,Tr,f_on_Sr,S,Sr,controls.paral,controls.pts_tol);
    intnew = f_on_Tr*Tr.weights';
    
    newp = length(idx_newp);
    nb_pts=nb_pts + newp;
    nb_pts_log(end+1)=nb_pts; 
    num_evals = nb_pts;
    Hr=[];
    f_on_Hr=[];
    
else
    
    % in this case, we need to keep track of all the points explored, even those that have been discarded in
    % previous iterations
    
    if strcmp(controls.recycling,'priority_to_evaluation') 
        %here we allow for multiple evaluations of the same point because we recycle from the previous grid only.
        % if the function evaluation is "cheap" this is much faster, because the search for common points relies
        % on multiindices and not on comparison of coordinates
        [f_on_Tr,~,idx_newp] = evaluate_on_sparse_grid(f,T,Tr,f_on_Sr,S,Sr,controls.paral,controls.pts_tol);
        intnew = f_on_Tr*Tr.weights';
        
        Hr.knots = [Hr.knots, Tr.knots(:,idx_newp)];
        f_on_Hr = [f_on_Hr, f_on_Tr(:,idx_newp)];

        newp = prod(lev2knots(ng_idx));
        nb_pts = size(Tr.knots,2); 
        num_evals = num_evals + length(idx_newp);
        nb_pts_log(end+1)=nb_pts;
        
    else
        % here we want to make sure no multiple evaluations of the same point occur. Thus we look in 
        % all points ever visited, but this is expensive because we rely on point coordinates only!

        [f_on_Tr,~,idx_newp] = evaluate_on_sparse_grid(f,T,Tr,f_on_Hr,[],Hr.knots,controls.paral,controls.pts_tol); 
        intnew = f_on_Tr*Tr.weights';
        
        Hr.knots = [Hr.knots, Tr.knots(:,idx_newp)];
        f_on_Hr = [f_on_Hr, f_on_Tr(:,idx_newp)];
        
        newp = prod(lev2knots(ng_idx));
        nb_pts = size(f_on_Tr,2); 
        num_evals = size(Hr.knots,2);
        nb_pts_log(end+1)=nb_pts;
    end
    % moreover, if profit is of type Linf, we need to evaluate the new grid on the ``nominally new points'',
    
    switch controls.prof
        
        case {'Linf/new_points','Linf','weighted Linf/new_points','weighted Linf'}
            
            Tx = tensor_grid(N,lev2knots(ng_idx),knots);
            new_points= Tx.knots;
            Tr_on_new_pts = interpolate_on_sparse_grid(T,Tr,f_on_Tr,new_points);
            
            
        case {'deltaint/new_points','deltaint'}
            
            % no need of new points
            
        otherwise
            
            error('do we need new points in this case? fix code here')
    end
    
end % closes if controls.nested


switch controls.prof
    
    case 'Linf/new_points'
        Sr_on_new_pts = interpolate_on_sparse_grid(S,Sr,f_on_Sr,new_points);
        if controls.nested
            Prof_temp = max(controls.op_vect(f_on_Tr(:,idx_newp),Sr_on_new_pts))/newp;
        else
            Prof_temp = max(controls.op_vect(Tr_on_new_pts,Sr_on_new_pts))/newp;
        end
        
    case 'Linf'
        Sr_on_new_pts = interpolate_on_sparse_grid(S,Sr,f_on_Sr,new_points);
        if controls.nested
            Prof_temp = max(controls.op_vect(f_on_Tr(:,idx_newp),Sr_on_new_pts));
        else
            Prof_temp = max(controls.op_vect(Tr_on_new_pts,Sr_on_new_pts));
        end
        
    case 'deltaint/new_points'
        delta_int  = controls.op_vect(intnew,intf);
        Prof_temp = delta_int / newp;
        
    case 'deltaint'
        delta_int  = controls.op_vect(intnew,intf);
        Prof_temp = delta_int;
        
    case 'weighted Linf/new_points'
        Sr_on_new_pts = interpolate_on_sparse_grid(S,Sr,f_on_Sr,new_points);
        if controls.nested
            Prof_temp = max( controls.op_vect(f_on_Tr(:,idx_newp),Sr_on_new_pts).*controls.pdf(new_points) )/newp;
        else
            Prof_temp = max( controls.op_vect(Tr_on_new_pts,Sr_on_new_pts).*controls.pdf(new_points) )/newp;
        end
        
    case 'weighted Linf'
        Sr_on_new_pts = interpolate_on_sparse_grid(S,Sr,f_on_Sr,new_points);
        if controls.nested
            Prof_temp = max( controls.op_vect(f_on_Tr(:,idx_newp),Sr_on_new_pts).*controls.pdf(new_points) );
        else
            Prof_temp = max( controls.op_vect(Tr_on_new_pts,Sr_on_new_pts).*controls.pdf(new_points) );
        end
    otherwise
        error('unknown profit indicator. Check spelling')
end

end