%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

%% a more complex example of recycle

clc
clear

fs=@(x) sum(x);

N=2; 


previous_evals=[];
S_old=[];
Sr_old=[];


for w=0:4
    
    S=smolyak_grid(N,w,@(n) knots_uniform(n,-1,1),@lev2knots_lin);
    Sr= reduce_sparse_grid(S);

    evals_rec = evaluate_on_sparse_grid(fs,S,Sr,previous_evals,S_old,Sr_old);    
    evals_non_rec = evaluate_on_sparse_grid(fs,Sr);
    
    previous_evals=evals_rec;
    S_old = S;
    Sr_old=Sr;

    max(abs(evals_non_rec(:)-evals_rec(:))) 

end



%% it is also possible to extract information about what points are new and what points are in common between to grids, 
% and what points are in the old sparse grid only

clc
clear

f=@(x) sum(x);

N=2; w=4;
S=smolyak_grid(N,w,@(n) knots_uniform(n,-1,1),@lev2knots_lin);
Sr= reduce_sparse_grid(S);

w=5;
T=smolyak_grid(N,w,@(n) knots_uniform(n,-1,1),@lev2knots_lin);
Tr= reduce_sparse_grid(T);


tol=1e-16;
evals_old=evaluate_on_sparse_grid(f,Sr);
[evals_par,new_points,tocomp_list,discard_points,discard_list] =evaluate_on_sparse_grid(f,T,Tr,evals_old,S,Sr);




%% an extreme test: we recycle from a larger grid to a smaller

clc
clear

f=@(x) sum(x);

N=2; w=3;
S=smolyak_grid(N,w,@(n) knots_uniform(n,-1,1),@lev2knots_lin);
Sr= reduce_sparse_grid(S);

w=4;
T=smolyak_grid(N,w,@(n) knots_uniform(n,-1,1),@lev2knots_lin);
Tr= reduce_sparse_grid(T);

evals_nr=evaluate_on_sparse_grid(f,Sr);
evals_r=evaluate_on_sparse_grid(f,S,Sr,evaluate_on_sparse_grid(f,Tr),T,Tr);

%[i,j]=max(abs(evals_nr(:)-evals_r(:))) 
max(abs(evals_nr(:)-evals_r(:))) 