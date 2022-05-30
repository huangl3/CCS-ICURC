# ICURC for Cross-Concentrated Sampling Model 
This is Matlab repo for an efficient non-convex matrix completion algorithm, termed Iterative CUR Completion (ICURC):

[1] HanQin Cai, Longxiu Huang, Pengyu Li, and Deanna Needell. Matrix Completion with Cross-Concentrated Sampling: Bridging Uniform Sampling and CUR Sampling

###### To display math symbols properly, one may have to install a MathJax plugin. For example, [MathJax Plugin for Github](https://chrome.google.com/webstore/detail/mathjax-plugin-for-github/ioemnmodlmafdkllaclgeombjnmnbima?hl=en).


## Introduction
In this work, we propose a novel and easy-to-implement sampling strategy, coined Cross-Concentrated Sampling (CCS). Additionally, we propose a highly efficient non-convex algorithm, named Iterative CUR Completion (ICURC), for the proposed CCS model. More details can be found in our paper [1].  


## Syntex
Using all default parameters:
```
[X_Omega_css, I_css, J_css] = CCS(X, p, delta);

[C,U_r,R, fct_all_time, ite] = ICURC(X_Omega_css, I_css, J_css, r,'');
```

Using custom parameters:
```
params.eta = [1/p, 1/p, 1/(2*p)];
params.TOL = 1e-4;
params.max_ite = 500;
[X_Omega_css, I_css, J_css] = CCS(X, p, delta);

[C,U_pinv,R, ICURC_time, ICURC_ite] = ICURC(X_Omega_css, I_css, J_css, r,params);
```

