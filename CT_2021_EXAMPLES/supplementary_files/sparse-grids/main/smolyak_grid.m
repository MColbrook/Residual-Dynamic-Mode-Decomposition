function [S,C] = smolyak_grid(N,w,knots,lev2knots,idxset,arg6,weights_coeff)


%  SMOLYAK_GRID generates a Smolyak sparse grid (and corresponding quadrature weights)
%  as a linear combination of full tensor grids, employing formula (2.9)
%  of [Nobile-Tempone-Webster, SINUM 46/5, pages 2309-2345] 
%
%  [S,C] = SMOLYAK_GRID(N,W,KNOTS,LEV2KNOTS,IDXSET) creates a sparse grid in N dimensions
%       using the multiindex set defined by
%
%       IDXSET(I) <= W,   
%
%       where  I is N-dimensional multiindex, W is an integer non-negative value and IDXSET is a function.
%       IDXSET is an optional argument, the default being IDXSET = @(i) sum(i-1).
%
%       KNOTS is either a cell array containing the functions to be used to generate the knots 
%       in each direction, i.e. 
%   
%       KNOTS={@knots_function1, @knots_function2, ... }
%
%       or a single function to be used in every direction, i.e.  KNOTS=@knots_function1.
%       In both cases, the header of knots_function is [x,w]=knots_function(m)
%
%       LEV2KNOTS is either a cell array containing the functions defining the relation between level 
%       and number of knots to be used in each direction, i.e. 
%
%       LEV2KNOTS={@m_function1, @m_function2, ... }
%
%       or a single function to be used in every direction, i.e. LEV2KNOTS=@m_function1
%       In both cases, the header of m_function is m=m_function(i)
%
%       The sparse grid information is stored as a vector of "tensor grids", 
%       each "tensor grid" S(j) is a four fields structure:
%           S(j).knots: vector containing the tensor grid knots
%           S(j).weights: vector containing the corresponding weights
%           S(j).size: size of the tensor grid, S(j).size = prod(m)
%           S(j).knots_per_dim: cell array (N components), each component is the set of 1D knots used 
%               to build the tensor grid    
%           S(j).m: the input vector m, m == lev2knots(idx), m(i)==length(S(j).knots_per_dim(i))
%           S(j).coeff: how many times the tensor grid appears in the sparse grid (with sign)
%               the index j runs over all points in the level set Y(W,N) 
%           S(j).idx: the multiidx vector corresponding to the current grid, whose number of points
%               is defined by lev2knots(i)
%
%       The outputs of SMOLYAK_GRID are 
%       S: structure containing the information on the sparse grid (vector of tensor grids; see above)
%       C: multi-index set used to generate the sparse grid  
%
%
%
%
% [S,C] = SMOLYAK_GRID(N,W,KNOTS,LEV2KNOTS,IDXSET,S2), where S2 is another Smolyak grid, 
%       tries to recycle tensor grids from S2 to build those of S instead of recomputing them.
%       This can be helpful whenever sequences of Smolyak grid are generates. Note that *NO* check
%       will performed whether S2 was generated with the same lev2knots as the one given as input.
%       S2 can also be empty, S2=[]
%
%
%
%  [S,C] = SMOLYAK_GRID(N,W,KNOTS,LEV2KNOTS,IDXSET,MAP,WEIGHTS_COEFF) can be used as an alternative
%       to generate a sparse grid on a hyper-rectangle. 
%       Instead of typing out one KNOTS function and one LEV2KNOTS for each dimension, like in
%
%       [S,C] = SMOLYAK_GRID(N,W,{@knots1, @knots2, ...},{@m1, @m2 ...}),
%
%       one can use 
%
%       [S,C] = SMOLYAK_GRID(N,W,@KNOTS,@M,MAP,WEIGHTS_COEFF),
%
%       which generates the sparse grid on the hypercube corresponding to @KNOTS and and then shifts it 
%       according to the mapping defined by MAP, e.g. from (-1,1)^N to (a1,b1)x(a2,b2)x...x(a_N,b_N). 
%       See also GET_INTERVAL_MAP. Specifying the scalar value WEIGHTS_COEFF will also multiply the 
%       quadrature weights by WEIGHTS_COEFF. Use IDXSET=[] to use the default value, IDXSET = @(i) sum(i-1).    



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


%----------------------------------------------------------
% input handling

if nargin==4 || (nargin==6 && isempty(idxset)) || (nargin==7 && isempty(idxset)) 
    idxset=@(i) sum(i-1);
end

% discriminate between calls:
if exist('arg6','var') 
    if isa(arg6,'function_handle') % we are in the case smolyak_grid(N,w,knots,lev2knots,idxset,map,weights_coeff)
        map = arg6;
    elseif issmolyak(arg6) || isempty(arg6) % we are in the case smolyak_grid(N,w,knots,lev2knots,idxset,S2)
        S2 = arg6;
    else
        error('SparseGKit:WrongInput','unknown type for 6th input')
    end
    clear arg6
end



% if knots and  lev2knots are simple function, we replicate them in a cell
if isa(knots,'function_handle')
    fknots=knots;
    knots=cell(1,N);
    for i=1:N
        knots{i}=fknots;
    end
end
if isa(lev2knots,'function_handle')
    f_lev2knots=lev2knots;
    lev2knots=cell(1,N);
    for i=1:N
        lev2knots{i}=f_lev2knots;
    end
end



if w==0

    %----------------------------------------------------------
    % the trivial case

    
    i = ones(1,N);
    m = apply_lev2knots(i,lev2knots,N);
    S(1) = tensor_grid(N,m,knots);
    S(1).coeff=1;
    S(1).idx = i;
    %disp('using 1 multiindex')
    C=i;

else

    
    %----------------------------------------------------------
    % let's go with the sparse construction
    
    
    % build the list of multiindices in the set: idxset(i)<=w.
    C=multiidx_gen(N,idxset,w,1);

    
    % we also build the set of S2 if provided
    if exist('S2','var') && ~isempty(S2)
        % get index set of S2. Note that if S2 has empty fields, C2==[]
        nb_idx_C2 = length(S2);
        C2 = reshape([S2.idx],N,nb_idx_C2)'; % [S2.idx] puts all indices on a row. I reshape them and each N elements I get one column.
        % Then I need to traspose the matrix
    end
    
    
    
    % now compute the tensor grids out of the delta operators listed in C.
    % Exploit partial ordering of the sequence of multiindices
    %-----------------------------------------------
    %
    % each multiindex of C is a delta grid. Given say [i1 i2 i3] i will get
    % (i1 - i1-1) x (i2 - i2-1) x (i3 - i3-1)
    % that is the grids associated to
    % [i1 i2 i3],[i1-1 i2 13],[i1 i2-1 i3],[i1-1 i2-1 i3] and so on
    %
    % note that all of these multiindeces are also included in the initial set C
    % (if [i1 i2 i3] ratisfies the rule, all other will, because their indeces are equal or lower)
    %
    % So the set of all possible grids (non delta grids) is the same set C, but some of them will cancel out
    %
    % C is partially ordered (lexicographis): [x x x x i], are listed increasing with i,
    % [x x x j i] are listed increasing first with j and then with i ...
    % Now take a row of C, c. Because of the ordering, if you take c as a grid index
    % the same grid can appear again only from delta grids coming from rows following.
    %
    % Now we scroll these rows following (say c2) and compute how many times c will be generated, and with
    % what sign. It has to be that d=c2-c has only 0 or 1
    %
    % if this is the case, the sign of this new appearence of c could be both + or - .
    % To determine the sign, start with + and switch it every time a 1 appears in d
    %
    % d=[0 1 0 0 1] => sign= +
    %
    % the formula is (-1)^sum(d) ( d with even appearences of 1 gives a + , odd gives a -)
    %
    % and just sum all the coefficients to see if c will survive or cancel out
    
    nn=size(C,1);
    coeff=ones(1,nn); % initialize coefficients to 1: all c survive
    
    % I can at least restrict the search to multiindices whose first component is c(i) + 2, so I define
    [~,bookmarks]=unique(C(:,1),'first');
    bk = [bookmarks(3:end)'-1 nn nn];
    % i.e. those who begin with 1 end at bookmark(3)-1, those who begin with 2-1 end at bookmark(4) and so on, 
    % until there's no multiindex with c(i)+2

    for i=1:nn % scroll c
        cc = C(i,:);
        % recover the range in which we have to look. Observe that the first column of C contains necessarily 1,2,3 ...
        % so we can use them to access bk
        range=bk(cc(1));    
        for j=i+1:range
            % scroll c2, the following rows
            d=C(j,:)-cc;
            if max(d)<=1 && min(d)>=0  % much faster to check then if only 0 and 1 appears. Also, && is short-circuited,
                % so if max(d)<=1 is false the other condition is not even checked
                coeff(i)=coeff(i) + (-1)^sum(d);
            end
        end 
    end

        
    % now we can store only those grids who survived, i.e. coeff~=0
    %------------------------------------------------------
    
    nb_grids=sum(coeff~=0);
    empty_cells=cell(1,nb_grids);
    S=struct('knots',empty_cells,'weights',empty_cells,'size',empty_cells,'knots_per_dim',empty_cells,'m',empty_cells);
    fieldnms = fieldnames(S)'; % I'll later need to loop over these fields - note the transpose, it has to be a row cell
    coeff_condensed=zeros(1,nb_grids);
    ss=1;

    % for each nonzero coeff, generate the tensor grid and store it. If possible, recycle from S2.
    if exist('C2','var') && ~isempty(C2)
        if MATLAB_SPARSE_KIT_VERBOSE
            disp('build smolyak grid with tensor grid recycling')
        end
        for j=1:nn
            if coeff(j)~=0
                i = C(j,:);
                [found,pos] = find_lexicographic(i,C2,'nocheck'); % this exploits that C2 is lexicographic, so it's efficient (cost ~ log(nb_rows_C2))
                if found
                    %disp('found')
                    % Note that at this point elements of S are tensor grids while S2 is a sparse grid therefore it has additional fields
                    % (coeff, idx). We thus need to copy field by field otherwise we'll have "assignment between dissimilar
                    % structures" error. We use dynamic filed names to this end
                    for fn = fieldnms
                        S(ss).(fn{1}) = S2(pos).(fn{1}); % note that each fn is a 1-element cell, so to access its value I need the notation fn{1}
                    end
                    % however we need to fix the weights. Indeed, they are stored in S2 as weights*coeff, so we need to reverse
                    % that multiplication
                    S(ss).weights = S(ss).weights/S2(pos).coeff;
                else
                    m =apply_lev2knots(i,lev2knots,N); % n. of points in each direction
                    S(ss) = tensor_grid(N,m,knots);
                end
                S(ss).weights=S(ss).weights*coeff(j);
                coeff_condensed(ss)=coeff(j);
                ss=ss+1;
            end
        end
    else
        for j=1:nn
            if coeff(j)~=0
                i = C(j,:);
                m =apply_lev2knots(i,lev2knots,N); % n. of points in each direction
                S(ss) = tensor_grid(N,m,knots);
                S(ss).weights=S(ss).weights*coeff(j);
                coeff_condensed(ss)=coeff(j);
                ss=ss+1;
            end
        end
    end
        
    % finally, shift the points according to map if needed
    if exist('map','var') && ~isempty(map)
        for ss=1:nb_grids
            S(ss).knots = map(S(ss).knots);
        end
    end
    
    % and possibly fix weights
    if exist('weights_coeff','var') && ~isempty(weights_coeff)
        for ss=1:nb_grids
            S(ss).weights = S(ss).weights*weights_coeff;
        end
    end
    
    % now store the coeff value. It has to be stored after the first loop, becuase tensor_grid returns a grid
    % WITHOUT coeff field, and Matlab would throw an error (Subscripted assignment between dissimilar structures)
    
    for ss=1:nb_grids
        S(ss).coeff=coeff_condensed(ss);
    end
     
    % similarly for the multiidx generating each tensor grid
    ss=1;
    for j=1:nn
        if coeff(j)~=0
            i = C(j,:);
            S(ss).idx = i; 
            ss=ss+1;
        end
    end

    
    

end % this end closes if w == 0


end % this end closes the function









function m = apply_lev2knots(i,lev2knots,N)
    
% m = apply_lev2knots(i,lev2knots,N)
%
% return a vector m s.t. m(n) = m(i(n)). m has to be single (or double) precision or following functions won't work 
% (more specifically, prod is not defined for intXY classes)


% N could be deduced by N but it's better passed as an input, to speed computation
% init m to zero vector
m=0*i;

% next, iterate on each direction 1,...,N. 

for n=1:N
    m(n) = lev2knots{n}(i(n)); 
end

end

