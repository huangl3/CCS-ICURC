# ICURC for Cross-Concentrated Sampling Model 
This is Matlab repo for an efficient non-convex algorithm, termed Iterative CUR Completion (ICURC), for Cross-Concentrated Sampling based matrix completion. For the details of Cross-Concentrated Sampling and ICURC, read and cite our latest TPAMI paper:

[1] HanQin Cai, Longxiu Huang, Pengyu Li, and Deanna Needell. <a href=https://doi.org/10.1109/TPAMI.2023.3261185>Matrix Completion with Cross-Concentrated Sampling: Bridging Uniform Sampling and CUR Sampling</a>. *IEEE Transactions on Pattern Analysis and Machine Intelligence*, 2023. doi: 10.1109/TPAMI.2023.3261185.

###### To display math symbols properly, one may have to install a MathJax plugin. For example, [MathJax Plugin for Github](https://chrome.google.com/webstore/detail/mathjax-plugin-for-github/ioemnmodlmafdkllaclgeombjnmnbima?hl=en).


## Introduction
In this work, we propose a novel and easy-to-implement sampling strategy, coined Cross-Concentrated Sampling (CCS). Additionally, we propose a highly efficient non-convex algorithm, named Iterative CUR Completion (ICURC), for the proposed CCS model. More details can be found in our paper [1].  


## Syntex
Using all default parameters:
```
[X_Omega_css, I_css, J_css] = CCS(X, '');

[C,U_pinv,R, ICURC_time, ICURC_ite] = ICURC(X_Omega_css, I_css, J_css, r,'');
```

Using custom parameters:
```
params_ccs.delta = 0.2;
params_ccs.p = 0.3;
[X_Omega_css, I_css, J_css] = CCS(X, params_ccs);

p = params_ccs.p;
params_icurc.eta = [1/p, 1/p, 1/(2*p)];
params_icurc.TOL = 1e-4;
params_icurc.max_ite = 500;
[C,U_pinv,R, ICURC_time] = ICURC(X_Omega_css, I_css, J_css, r,params_icurc);
```

## Input Description for CCS
1. X : low rank matrix.  
2. params_ccs : parameters for the algorithm CCS
   * .delta : rate of sampled columns or rows. (default 0.2)
   * .p : observation rate on the selected submatrices. (default 0.2)


## Output Description for CCS
1. X_Omega_css : observed data matrix based on CSS sampling model
2. I_css : row indices of the selected row submatrix
3. J_css : column indices of the selected column submatrix

## Input Description for ICURC
1. X_Omega_css : observed data matrix based on CSS sampling model
2. I_css : row indices of the selected row submatrix
3. J_css : column indices of the selected column submatrix
4. r : rank of X.
5. params_icurc : parameters for the algorithm ICURC
   * .max_iter : Maximum number of iterations. (default 500)
   * .TOL : Desired Frobenius norm error. (default 1e-4)
   * .eta :  eta(1), eta(2), and eta(3) the step sizes for updating C, R, and U. (default [1 1 1])

## Output Description for ICURC
1. C，U_pinv，R : CUR decomposition of $X = C U^\dagger R$, where $U^\dagger$ is the pseudo-inverse of $U$.
2. ICURC_time : runtime for ICURC.

