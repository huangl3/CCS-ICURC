function [X_Omega_css, I_css, J_css] = CCS(X, params_CCS)
%Return the observed data matrix based on CSS sampling model along with 
%row and column indices of the selected row submatrix 
%
%Inputs:
%   X - low rank matrix.
%   params_CCS - parameter structure containing the following fields:
%       delta - rate of sampled columns or rows
%       p - observation rate on the selected submatrices

%Outputs:
%   X_Omega_css - observed data matrix based on CSS sampling model
%   I_css - row indices of the selected row submatrix
%   J_css - column indices of the selected column submatrix

    if(~exist('params_CCS','var'))
        params_CCS=struct();
    end
    params_CCS = SetDefaultParams_CCS(params_CCS);
    p=params_CCS.p;
    delta=params_CCS.delta;
    
    [m, n] = size(X);
    % Generate I_css and J_css based on delta

    num_r = round(m*delta);
    num_c = round(n*delta);
    I_css = randsample(m,num_r,false);
    J_css = randsample(n,num_c,false);     
    
    C = X(:,J_css);
    R = X(I_css,:);
    
    ubc = min(num_c * m, ceil(p*num_c * m));
    ubr = min(num_r * n, ceil(p*num_r * n));
    
    % Randomly sample entrice for column and row respectively.
    C_obs_ind = randsample(num_c * m, ubc); 
    R_obs_ind = randsample(num_r * n, ubr); 
    
    C_Obs = zeros(m,num_c);
    R_Obs = zeros(num_r,n);
    
    
    C_Obs(C_obs_ind) = C(C_obs_ind);
    R_Obs(R_obs_ind) = R(R_obs_ind);

    % Generate observed data
    X_Omega_css = zeros(m, n);
    X_Omega_css(I_css, :) = R_Obs;
    X_Omega_css(:, J_css) = C_Obs;
end

