%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


%%
clear

% dimension 
N=6;

% vector-valued function to be interpolated. input column points, out
% 2-row matrix, one per output component 

f1=@(x) 1./(1+0.5/N*sum(x)); 
f2=@(x) sum(abs(x).^3); 

f=@(x) [f1(x); f2(x)];

% define sparse grid
[lev2knots,idxset]=define_functions_for_rule('TD',N);
knots=@(n) knots_uniform(n,-1,1,'nonprob');
a=-1; b=1;



% for loop

w_max=6;
interp_error=zeros(2,w_max+1);
work=zeros(1,w_max+1);


non_grid_points=rand(N,2000)*(b-a)+a;
%non_grid_points=[0.5*ones(N,1), zeros(N,1)];
S_old = [];  % we also recycle previous grids

tic
for w=0:w_max

    disp(w)

    % create grid
    [S,C]=smolyak_grid(N,w,knots,lev2knots,idxset,S_old);
    S_old = S;
    
    Sr=reduce_sparse_grid(S);

    
    % move points to actual interval here point are columns
    % Sr.knots=interval_map(Sr.knots);

    % compute work
    work(w+1)=Sr.size;
    
    % compute estimate of polynomial size
    pol_size(w+1)=size(C,1);
    
    % compute the nodal values to be used to interpolate. It has to be
    % row vector (more rows for vector-valued output functions)
    function_on_grid=f(Sr.knots);   
    
    % compute interpolated values. Here f_values is row
    f_values = interpolate_on_sparse_grid(S,Sr,function_on_grid,non_grid_points);

    % compute error
    interp_error(:,w+1)=max( abs((f(non_grid_points) - f_values)), [],2 ) ;    

end
toc

%%
figure
semilogy(0:w,interp_error(1,:),'-or','DisplayName','Numerical Error component 1, w.r.t. grid level');
hold on
semilogy(0:w,interp_error(2,:),'-ok','DisplayName','Numerical Error component 2, w.r.t. grid level');
semilogy(0:w_max,1./exp(0:w_max),'--o','DisplayName','exp(-level)')
legend show

%%
figure
loglog(work,interp_error(1,:),'-or','DisplayName','Numerical Error component 1, w.r.t. #points');
hold on
loglog(work,interp_error(2,:),'-ok','DisplayName','Numerical Error component 2, w.r.t. #points');
legend show

