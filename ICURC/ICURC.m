function [C,U_r,R, fct_all_time, ite] = ICURC(X_Omega_UR, Ind_I, Ind_J, r, params)
%Return the CUR components of result of solved Matrix Completion problem
%under Cross-Concentrated Sampling (CCS). 
%
%Inputs:
%   X_Omega_UR - observed data
%   Ind_I - row indices that define R
%   Ind_J - column indices that define C
%   r - target rank
%   params - parameter structure containing the following fields:
%       TOL,max_ite - stopping criteria
%           defaults = 1e-4 (TOL), 500 (mxitr)
%       eta - step size
%           default = [1, 1, 1] for [eta_C, eat_R, era_U]
%       steps_are1 - if all step sizes are 1
%           default = true


%Outputs:
%   C, U_r, R - The CUR components of recovered X
%   fct_all_time -Function execution time
%   ite - Number of iterations required

    if(~exist('params','var'))
        params=struct();
    end
    params = SetDefaultParams(params);
    eta=params.eta;
    TOL=params.TOL;
    max_ite=params.max_ite;
    steps_are1=params.steps_are1;
    
    %This step is to extract observed C, U, and R 
    Obs_U = X_Omega_UR(Ind_I, Ind_J);  
    Obs_C = X_Omega_UR(:, Ind_J);
    Obs_R = X_Omega_UR(Ind_I, :);
    C_size = size(Obs_C);
    R_size = size(Obs_R);

    all_row_ind = 1:C_size(1); 
    all_col_ind = 1:R_size(2); 

    Ind_I_comp = setdiff(all_row_ind, Ind_I); 
    Ind_J_comp = setdiff(all_col_ind, Ind_J); 
    
    C = Obs_C(Ind_I_comp, :); 
    R = Obs_R(:, Ind_J_comp); 
    %C and R are C_obs\U, R_obs\U.
    Smp_C = C ~= 0;
    Smp_R = R ~= 0; 
    Smp_U = Obs_U ~= 0; 
    L_obs_only_row = R;
    L_obs_only_col = C;  
    L_obs_only_U = Obs_U; 
    
    Omega_row = find(Smp_R);    
    Omega_col = find(Smp_C);  
    Omega_U = find(Smp_U);  
    
    L_obs_row_vec = L_obs_only_row(Omega_row);
    L_obs_col_vec = L_obs_only_col(Omega_col);   
    L_obs_U_vec = L_obs_only_U(Omega_U); 
    
    normC_obs = norm(L_obs_col_vec,'fro');
    normU_obs = norm(L_obs_U_vec,'fro');
    normR_obs = norm(L_obs_row_vec,'fro'); 

    col_row_norm_sum = normC_obs + normU_obs + normR_obs;    
    
    %Initializing U
    U_i = Obs_U;
    [u,s,v] = svd(U_i);  
    u = u(:, 1:r);
    v = v(:, 1:r);
    s = s(1:r, 1:r);  
    U_i = u*s*v'; 
    
    %Calculating error
    New_Error = (norm(R(Omega_row) - L_obs_row_vec,'fro') + ... 
                            norm(C(Omega_col)- L_obs_col_vec, 'fro') + ... 
                            norm(U_i(Omega_U)- L_obs_U_vec))/col_row_norm_sum;
    %If all step sizes are 1
    if steps_are1 
        fct_time = tic;
        for ite = 1:max_ite

            R = u*u'*R;
            C = C*v*v';        
            Old_error = New_Error;

            New_Error = (norm(R(Omega_row) - L_obs_row_vec,'fro') + ... 
                            norm(C(Omega_col)- L_obs_col_vec, 'fro') + ... 
                            norm(U_i(Omega_U)- L_obs_U_vec))/col_row_norm_sum;  
            if New_Error < TOL || ite == max_ite 
                fct_all_time = toc(fct_time); 
                R(Omega_row) = L_obs_row_vec;
                C(Omega_col) = L_obs_col_vec;
                U_i(Omega_U) = L_obs_U_vec; 

                Final_C = zeros(C_size);
                Final_R = zeros(R_size); 
                Final_C(Ind_I_comp, :) = C;
                Final_R(:, Ind_J_comp) = R;
                Final_C(Ind_I, :) = U_i; 
                Final_R(:, Ind_J) = U_i; 

                C = Final_C;
                R = Final_R; 
                
                [u,s,v] = svd(U_i);  
                u = u(:, 1:r);
                v = v(:, 1:r);
                s = s(1:r, 1:r);  
                U_r = v*pinv(s)*u';
                return
            end
            %Updating R, C, and U
            R(Omega_row) = L_obs_row_vec;
            C(Omega_col) = L_obs_col_vec;
            U_i(Omega_U) =  L_obs_U_vec; 

            [u,s,v] = svd(U_i);  
            u = u(:, 1:r);
            v = v(:, 1:r);
            s = s(1:r, 1:r);  
            U_i = u*s*v';
        end
        
    %If all step sizes are not 1    
    else
        step_c = eta(1);
        step_r = eta(2);
        step_u = eta(3);
        
        fct_time = tic;
        for ite = 1:max_ite

            R = u*u'*R;
            C = C*v*v';        
            Old_error = New_Error;
            New_Error = (norm(R(Omega_row) - L_obs_row_vec,'fro') + ... 
                            norm(C(Omega_col)- L_obs_col_vec, 'fro') + ... 
                            norm(U_i(Omega_U)- L_obs_U_vec))/col_row_norm_sum;   
                        
            if New_Error < TOL || ite == max_ite
                fct_all_time = toc(fct_time); 

                R(Omega_row) = L_obs_row_vec;
                C(Omega_col) = L_obs_col_vec;
                U_i(Omega_U) = L_obs_U_vec; 

                Final_C = zeros(C_size);
                Final_R = zeros(R_size); 
                Final_C(Ind_I_comp, :) = C;
                Final_R(:, Ind_J_comp) = R;
                Final_C(Ind_I, :) = U_i; 
                Final_R(:, Ind_J) = U_i; 

                C = Final_C;
                R = Final_R; 

                [u,s,v] = svd(U_i);  
                u = u(:, 1:r);
                v = v(:, 1:r);
                s = s(1:r, 1:r);  
                U_r = v*pinv(s)*u';
                return
            end  
            %Updating R, C, and U
            R(Omega_row) = (1 - step_r)*R(Omega_row) + step_r*(L_obs_row_vec);     
            C(Omega_col) = (1 - step_c)*C(Omega_col) + step_c*(L_obs_col_vec);
            U_i(Omega_U) =  (1 - step_u)*U_i(Omega_U) + step_u*(L_obs_U_vec); 

            [u,s,v] = svd(U_i);  
            u = u(:, 1:r);
            v = v(:, 1:r);
            s = s(1:r, 1:r);  
            U_i = u*s*v';
         end
     end
end