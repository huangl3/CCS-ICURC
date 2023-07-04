close all; clear all; clc;
dim = [3000, 2000];      
m = dim(1);
n = dim(2); 
r = 5; %rank of the optimal matrix

params_CCS.p = 0.3; %uniform observation rate on the submatrices
params_CCS.delta = 0.2;%percentage of sampled columns or rows

params_ICURC.eta = [1/params_CCS.p, 1/params_CCS.p, 1/(2*params_CCS.p)]; %step sizes for updating C, R, and U
params_ICURC.TOL = 1e-4;
params_ICURC.max_ite = 500;

%%
disp('Generating a CCS matrix completion problem.')

% Generate the underlying matrix with rank = r
A_generater = randn(m,r);
B_generater = randn(r,n);
X = A_generater * B_generater;

% Generate observed data under CCS with give p and delta
[X_Omega_ccs, I_ccs, J_ccs] = CCS(X, params_CCS); 

fprintf('Running ICURC with m=%d, n=%d, r=%d \n',m,n,r);

%% Run ICURC 
[C, U_pinv, R, ICURC_time] = ICURC(X_Omega_ccs, I_ccs, J_ccs, r, params_ICURC);
Mout_CURf = C* U_pinv *R; 
Error = norm(Mout_CURf - X,'fro') / norm(X,'fro');

fprintf('ICURC recovers the ground truth with relative error %f \n', Error);

