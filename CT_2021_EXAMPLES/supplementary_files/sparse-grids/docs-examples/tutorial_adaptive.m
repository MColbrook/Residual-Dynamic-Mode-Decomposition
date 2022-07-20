
%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


% an example 2D

clear

f=@(x) 1./(x(1)^2+x(2)^2 + 0.3);
N=2;
a=-1;
b=1;
knots=@(n) knots_CC(n,a,b);
lev2knots=@lev2knots_doubling;

controls.paral=NaN; %no parallel evaluation of f over grids
controls.max_pts=200;
controls.prof_toll = 1e-10;
prev_adapt = [];
controls.nested=true;
controls.plot=false;
adapt1 = adapt_sparse_grid(f,N,knots,lev2knots,prev_adapt,controls);

% plot the indices used 
plot_idx_status(adapt1.private.G,adapt1.private.I,adapt1.private.idx_bin,adapt1.private.idx)
axis([0 12 0 12])
set(gca,'FontSize',14)

%pdfsaving('idx-status')

%% --------------------------------------------------------------
% set a profit indicator other than the default one


clear

f=@(x) 1./(x(1)^2+x(2)^2 + 0.3) ;
N=2;
a=-1;
b=1;
knots=@(n) knots_CC(n,a,b);
lev2knots=@lev2knots_doubling;

controls.paral=NaN;
controls.max_pts=200;
controls.prof_toll = 1e-10;
controls.prof='deltaint/new_points';
controls.nested=true;

adapt_prev = [];
adapt2 = adapt_sparse_grid(f,N,knots,lev2knots,adapt_prev,controls);


% plot the indices used 
plot_idx_status(adapt2.private.G,adapt2.private.I,adapt2.private.idx_bin,adapt2.private.idx)


%% ----------------------------------------------------------------
% increase the number of points: can recycle previous run

controls.max_pts=1500;

adapt3 = adapt_sparse_grid(f,N,knots,lev2knots,adapt2,controls);


% plot the indices used 
plot_idx_status(adapt3.private.G,adapt3.private.I,adapt3.private.idx_bin,adapt3.private.idx)




%% ---------------------------------------------------
% another example 2D, on an unbounded interval, with both nested and non-nested

clear

f=@(x) 1./(2+exp(x(1)) + exp(x(2))); 
N=2;

knots=@(n) knots_kpn(n);
lev2knots=@lev2knots_kpn;
controls.paral=NaN;
controls.max_pts=150;
controls.prof_toll = 1e-10;
controls.prof='weighted Linf/new_points';
%controls.prof='deltaint';
prev_adapt = [];
controls.nested=true;
controls.pdf = @(Y) prod(normpdf(Y,0,1),1); % note that we need to define a pdf here
adapt1 = adapt_sparse_grid(f,N,knots,lev2knots,prev_adapt,controls);
plot_idx_status(adapt1.private.G,adapt1.private.I,adapt1.private.idx_bin,adapt1.private.idx)


knots=@(n) knots_gaussian(n,0,1);
lev2knots=@lev2knots_lin;
controls.nested=false; % changing to nested false for gaussian
% here's the adapt non-nested. You will see some message like:
% "Some points have been evaluated more than once. Total: 191 extra evaluations over 295 function evaluations"
% you can change this behaviour by changing the default value of controls.recycling from 
%
% controls.recycling = 'priority_to_evaluation'
%
% to
%
% controls.recycling = 'priority_to_recycling'
%
% however, this is **not** recommended if N is large and evaluating f is cheap.
% see help ADAPT_SPARSE_GRID > CONTROLS.RECYCLING for more information

adapt2 = adapt_sparse_grid(f,N,knots,lev2knots,prev_adapt,controls);

plot_idx_status(adapt2.private.G,adapt2.private.I,adapt2.private.idx_bin,adapt2.private.idx)

adapt1.intf
adapt2.intf

controls.recycling = 'priority_to_recycling';
adapt3 = adapt_sparse_grid(f,N,knots,lev2knots,prev_adapt,controls);

% this will be false, because num_evals and nb_pts_visited are now identical
isequal(adapt2,adapt3)

%% observe that there's a maximum number of tabulated points with KPN that one will hit sooner or later, so asking too many points
% will result in an error


clear

f=@(x) 1./(2+exp(x(1)) + exp(x(2))) 
N=2;
knots=@(n) knots_kpn(n);
lev2knots=@lev2knots_kpn;

controls.paral=NaN; 
controls.max_pts=1500;
controls.prof_toll = 1e-10;
controls.pdf = @(Y) exp(-0.5*sum(Y.^2,1)); % define the weight for the profit. Not doing this will raise an error
controls.prof='weighted Linf/new_points';
prev_adapt = [];
controls.nested=true;

adapt1 = adapt_sparse_grid(f,N,knots,lev2knots,prev_adapt,controls);


plot_idx_status(adapt1.private.G,adapt1.private.I,adapt1.private.idx_bin,adapt1.private.idx)





%% also, using a non-weighted profit estimate will cause troubles, the interpolation is not converging in Linf sense on the 
% whole real axis so profit estimates are unreliable and they soon lead to hit the maximum number of points allowed

clear

f=@(x) 1./(2+exp(x(1)) + exp(x(2))) 
N=2;
knots=@(n) knots_kpn(n);
lev2knots=@lev2knots_kpn;

controls.paral=NaN;
controls.max_pts=1500;
controls.prof_toll = 1e-10;
controls.prof='Linf/new_points';
prev_adapt = [];
controls.nested=true;

adapt1 = adapt_sparse_grid(f,N,knots,lev2knots,prev_adapt,controls);


plot_idx_status(adapt1.private.G,adapt1.private.I,adapt1.private.idx_bin,adapt1.private.idx)




%% an example 2D with a vector-valued function

clear

f=@(x) [1./(x(1)^2+x(2)^2 + 0.3); 1./(x(1)^2+0.1*x(2)^2 + 2)];

N=2;
a=-1;
b=1;
knots=@(n) knots_CC(n,a,b);
lev2knots=@lev2knots_doubling;

controls.paral=NaN; %no parallel evaluation of f over grids
controls.max_pts=200;
controls.prof_toll = 1e-10;
prev_adapt = [];
controls.nested=true;
controls.plot=false;
adapt1 = adapt_sparse_grid(f,N,knots,lev2knots,prev_adapt,controls);

% adaptivity for vector-valued quantites is done by computing by default the euclidean norm of the output.
% by changing the way in which profit of vector valued quantities are computed, I can recover exactly the same
% behaviour as if I was considering only one of the two components (see help for more details)

% use e.g. this one to recover the scalar result for the first function only
controls.op_vect = @(A,B) sqrt(sum((A(1,:) - B(1,:)).^2,1));
adapt_f1 = adapt_sparse_grid(f,N,knots,lev2knots,prev_adapt,controls);

f1=@(x) 1./(x(1)^2+x(2)^2 + 0.3);
controls = rmfield(controls,'op_vect'); %restore default
adapt_f1_check = adapt_sparse_grid(f1,N,knots,lev2knots,prev_adapt,controls);


% use e.g. this one to recover the scalar result for the second function only
controls.op_vect = @(A,B) sqrt(sum((A(2,:) - B(2,:)).^2,1));
adapt_f2 = adapt_sparse_grid(f,N,knots,lev2knots,prev_adapt,controls);

f2=@(x) 1./(x(1)^2+0.1*x(2)^2 + 2);
controls = rmfield(controls,'op_vect'); %restore default
adapt_f2_check = adapt_sparse_grid(f2,N,knots,lev2knots,prev_adapt,controls);




%% ---------------------------------------------------
% an example 4D

clear

f=@(x) 1./(x(1)^2+x(2)^2 + 0.3 + 0.1*sin(x(3)).*exp(0.4*x(4)))
N=4;
a=-1;
b=1;
knots=@(n) knots_CC(n,a,b);
lev2knots=@lev2knots_doubling;

controls.paral=NaN; 
controls.max_pts=400;
controls.prof_toll = 1e-10;
prev_adapt = [];
controls.nested=true;
tic
adapt1 = adapt_sparse_grid(f,N,knots,lev2knots,prev_adapt,controls);
toc
% plot indices. PS those idx_bin in the corner next to the origin are 2D projection of 4D indices not in idx_bin, like
% [1 1 4 1], so that's ok
plot_idx_status(adapt1.private.G(:,1:2),adapt1.private.I(:,1:2),adapt1.private.idx_bin(:,1:2),adapt1.private.idx(:,1:2))
plot_idx_status(adapt1.private.G(:,3:4),adapt1.private.I(:,3:4),adapt1.private.idx_bin(:,3:4),adapt1.private.idx(:,3:4))



%% now we verify that using partial exploration we gain in computational work

clear


% the target function is a 25-dim function, only the first 3 are meaningful:
% f=@(y) 1./(4 + y(1) + 0.2*y(2) + 0.04*y(3) + 0*y(4:25) );
%
% so with a buffer of 2 we expect that dimadapt will explore only the first 5, with a gain
% of 20x2= 40 points (20 unexplored dim x 2 points to explore each of them to realize that they 
% are meaningless)
%
% Because f must be ablo to accept multiple dimensions of y, so we hack it as

aaa=[1 0.2 0.04 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
f=@(y) 1./(4 + dot(aaa(1:length(y)),y));

N=length(aaa);
a=-1;
b=1;
knots=@(n) knots_CC(n,a,b);
lev2knots=@lev2knots_doubling;

controls.paral=NaN; 
controls.prof_toll = 1e-10;
prev_adapt = [];
controls.nested=true;

% we do a convergence study for standard adapt, with 30 values of work between 0 and 300. The reference is set
% at 500 points
nb_pts = 30; 
max_pts = [ceil(logspace(1,log10(300),nb_pts)) 500];
intf_vals = zeros(1,nb_pts);
true_pts = zeros(1,nb_pts);

k=1;
for p = max_pts
    controls.max_pts=p;
    adapt1 = adapt_sparse_grid(f,N,knots,lev2knots,prev_adapt,controls);
    intf_vals(k)=adapt1.intf;
    true_pts(k) = adapt1.nb_pts;
    %figure; spy(adapt1.private.I_log-1)
    if p<max_pts(end)
        prev_adapt=adapt1;
        k=k+1;
    end
end

% move the last solve to reference and delete it from the convergence study
adapt1 = prev_adapt;
intf_ex = intf_vals(end);
intf_vals(end)=[];
max_pts(end)=[];
true_pts(end)=[];

% now repeat for dimension adapt

controls.var_buffer_size = 2;
N_full = N;
dimad_intf_vals = zeros(1,nb_pts);
dimad_true_pts = zeros(1,nb_pts);
dimad_N = zeros(1,nb_pts);

prev_adapt=[];
k=1;
for p = max_pts
    controls.max_pts=p;
    adapt2 = adapt_sparse_grid(f,N_full,knots,lev2knots,prev_adapt,controls);
    dimad_intf_vals(k)=adapt2.intf;
    dimad_true_pts(k) = adapt2.nb_pts;
    dimad_N(k)=adapt2.N;
    %figure; spy(adapt2.private.I_log-1)
    prev_adapt=adapt2;
    k=k+1;
end


% we now plot the convergence. The 40 points gain is confirmed
figure
loglog(true_pts,abs(intf_vals-intf_ex),'-ob','LineWidth',2,'MarkerFaceColor','b','DisplayName','standard adapt')
hold on
loglog(dimad_true_pts,abs(dimad_intf_vals-intf_ex),'-or','LineWidth',2,'MarkerFaceColor','r','DisplayName','dim-adapt')
legend show
grid on

% we also plot the sequence with which indices are added to the grid
figure
subplot(1,2,1)
spy(adapt1.private.G_log-1)
subplot(1,2,2)
spy(adapt2.private.G_log-1)
 
% observe that after the initial part in which dimadapt gains its advantage, they then continue in the same order. 
figure
subplot(1,2,1)
spy(adapt1.private.G_log(30:end,1:3)-1)
subplot(1,2,2)
spy(adapt2.private.G_log(10:end,1:3)-1)

% the same can be seen from the pts count, which is parallel, 40 pts diff
figure
plot(adapt1.private.nb_pts_log(30:end),'b','DisplayName','adapt')
hold on
plot(adapt2.private.nb_pts_log(10:end),'r','DisplayName','dim adapt')
legend show


%% something similar happens also for non-nested points. Dspite the linear growth, using non-nested points requires 2 new points to assess that
% a dimension in useless (2 pts at level 2, different from the point at level 1) so the gain is still 20x2 =
% 40 points. Things are less evident though, because nb_pts_log reports the count of Tr and not of unique(Hr),
% even if I set priority to recycling the lines won't be parallel

clear

aaa=[1 0.2 0.04 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
f=@(y) 1./(4 + dot(aaa(1:length(y)),y));

N=length(aaa);
a=-1;
b=1;
knots=@(n) knots_uniform(n,a,b);
lev2knots=@lev2knots_lin;

controls.paral=NaN; 
controls.prof_toll = 1e-10;
prev_adapt = [];
controls.nested=false;

% we do a convergence study for standard adapt, with 30 values of work between 0 and 300. The reference is set
% at 500 points
nb_pts = 30; 
max_pts = [ceil(logspace(1,log10(300),nb_pts)) 500];
intf_vals = zeros(1,nb_pts);
true_pts = zeros(1,nb_pts);

k=1;
for p = max_pts
    controls.max_pts=p;
    adapt1 = adapt_sparse_grid(f,N,knots,lev2knots,prev_adapt,controls);
    intf_vals(k)=adapt1.intf;
    true_pts(k) = adapt1.nb_pts;
    %figure; spy(adapt1.private.I_log-1)
    if p<max_pts(end)
        prev_adapt=adapt1;
        k=k+1;
    end
end

% move the last solve to reference and delete it from the convergence study
adapt1 = prev_adapt;
intf_ex = intf_vals(end);
intf_vals(end)=[];
max_pts(end)=[];
true_pts(end)=[];

% now repeat for dimension adapt
controls.var_buffer_size = 2;
%controls.recycling = 'priority_to_recycling';
N_full = N;
dimad_intf_vals = zeros(1,nb_pts);
dimad_true_pts = zeros(1,nb_pts);
dimad_N = zeros(1,nb_pts);

prev_adapt=[];
k=1;
for p = max_pts
    controls.max_pts=p;
    adapt2 = adapt_sparse_grid(f,N_full,knots,lev2knots,prev_adapt,controls);
    dimad_intf_vals(k)=adapt2.intf;
    dimad_true_pts(k) = adapt2.nb_pts;
    dimad_N(k)=adapt2.N;
    %figure; spy(adapt2.private.I_log-1)
    prev_adapt=adapt2;
    k=k+1;
end

%% we now plot the convergence. The 40 points gain is confirmed
figure
loglog(true_pts,abs(intf_vals-intf_ex),'-ob','LineWidth',2,'MarkerFaceColor','b','DisplayName','standard adapt')
hold on
loglog(dimad_true_pts,abs(dimad_intf_vals-intf_ex),'-or','LineWidth',2,'MarkerFaceColor','r','DisplayName','dim-adapt')
legend show
grid on

% we also plot the sequence with which indices are added to the grid
figure
subplot(1,2,1)
spy(adapt1.private.G_log-1)
subplot(1,2,2)
spy(adapt2.private.G_log-1)
 
% observe that after the initial part in which dimadapt gains its advantage, they then continue in the same order. 

figure
subplot(1,2,1)
spy(adapt1.private.G_log(31:end,1:3)-1)
subplot(1,2,2)
spy(adapt2.private.G_log(11:end,1:3)-1)

% the same can be seen from the pts count, which is parallel, 40 pts diff
figure
plot(adapt1.private.nb_pts_log(31:end),'b','DisplayName','adapt')
hold on
plot(adapt2.private.nb_pts_log(11:end),'r','DisplayName','dim adapt')
legend show


%% a convergence study in L-inf norm (max in interpolation error)


clear

% function to be interpolated
f=@(x) 1./(x(1,:).^2+x(2,:).^2 + 0.3);

% domain is [a,b]^N
N=2;
a=-1;
b=1;

% settings for sparse grids
knots=@(n) knots_CC(n,a,b);
lev2knots=@lev2knots_doubling;
controls.paral=NaN; 
controls.prof_toll = 1e-10;
prev_adapt = [];
controls.nested=true;

% evaluate error as max error over 100 random points in [a,b]^2. Note that here we have hard-coded that a=-1,  b=1
nb_rand_pts = 100;
Rand_pts = 2*rand(2,nb_rand_pts)-1;
truef_evals = f(Rand_pts);

% generate a sequence of sparse grids with these many points (approximately), for each save values of interest 
% (exact nb pts, error,  approximation of integral of f)
max_pts = [5 7 13 21 29 50 80 200 400 600 1000];
PP = length(max_pts);
quadf_vals = zeros(1,PP);
sg_pts = zeros(1,PP);
sg_err =zeros(1,PP);


% the loop over the sparse grids
k=1;
for p = max_pts
    controls.max_pts=p;
    adapt1 = adapt_sparse_grid(f,N,knots,lev2knots,prev_adapt,controls);
    sg_eval = interpolate_on_sparse_grid(adapt1.S,adapt1.Sr,adapt1.f_on_Sr,Rand_pts);
    sg_err(k) = max(abs(sg_eval - truef_evals));
    quadf_vals(k)=adapt1.intf;
    sg_pts(k) = adapt1.nb_pts;
    if p<max_pts(end)
        prev_adapt=adapt1;
        k=k+1;
    end
end

% take the last computed integral as reference integral
quadf_ref = quadf_vals(end);
quadf_vals(end)=[];
max_pts(end)=[];

% error plots
figure
loglog(sg_pts(1:end-1),abs(quadf_vals-quadf_ref),'-ob','LineWidth',2,'MarkerFaceColor','b','DisplayName','adaptive sg')
title('quadrature error')
grid on

figure
loglog(sg_pts,sg_err,'-ob','LineWidth',2,'MarkerFaceColor','b','DisplayName','adaptive sg')
title('interp error')
grid on