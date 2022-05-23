close all; clear all; clc;
dim = [2000, 1000];      
m = dim(1);
n = dim(2); 
r = 5; %rank of the optimal matrix
p = 0.3; %uniform observation rate on the submatrices
delta = 0.2;%percentage of sampled columns or rows
max_ite = 500;
TOL = 1e-4;
eta = [1/p, 1/p, 1/(2*p)];
params.eta = eta;
params.TOL = TOL;
params.max_ite = max_ite;

%%
%Generate the underlying matrix with rank = r
A_generater = randn(m,r);
B_generater = randn(r,n);
X = A_generater * B_generater;
%Generate observed data under CCS with give p and delta
[X_Omega_UR, Ind_I, Ind_J] = CCS(X, p, delta);

fprintf('Running ICURC with m=%d, n=%d, r=%f...\n',m,n,r);
%Run ICURC 
[C,U_r,R, fct_all_time, ite] = ICURC(X_Omega_UR, Ind_I, Ind_J, r, params);
Mount_CURf = C*U_r*R; 
Error = norm(Mount_CURf - X,'fro') / norm(X,'fro');

fprintf('ICURC finished with relative error in frobenius norm=%f with %f iterations in time t=%f \n',Error, ite,fct_all_time);
