function [X_Omega_UR, Ind_I, Ind_J] = CCS(X, p, delta)
  
    [m, n] = size(X);
    % Generate Ind_I and Ind_J based on delta

    num_r = round(m*delta);
    num_c = round(n*delta);
    Ind_I = randsample(m,num_r,false);
    Ind_J = randsample(n,num_c,false);     
    
    C = X(:,Ind_J);
    R = X(Ind_I,:);
    
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
    X_Omega_UR = zeros(m, n);
    X_Omega_UR(Ind_I, :) = R_Obs;
    X_Omega_UR(:, Ind_J) = C_Obs;
end

