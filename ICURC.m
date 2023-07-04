function [C,U_pinv,R, ICURC_time] = ICURC(X_Omega_UR, I_ccs, J_ccs, r, params_ICURC)
%Return the CUR components of result of solved Matrix Completion problem
%under Cross-Concentrated Sampling (CCS). 
%
%Inputs:
%   X_Omega_ccs - observed data matrix based on CCS sampling model
%   I_ccs - row indices of the selected row submatrix
%   J_ccs - column indices of the selected column submatrix
%   r : rank of targt X
%   params_ICURC - parameter structure containing the following fields:
%       TOL,max_ite - stopping crICURC_iteria
%           defaults = 1e-4 (TOL), 500 (mxitr)
%       eta - step size
%           default = [1, 1, 1] for [eta_C, eat_R, era_U], step sizes for
%           updating C, R, and U respectively 
%       steps_are1 - if all step sizes are 1
%           default = true

%Outputs:
%   C, U_pinv, R - The CUR components of recovered X
%   ICURC_time -Function execution time
%   ICURC_ite - Number of ICURC_iterations required

    if(~exist('params_ICURC','var'))
        params_ICURC=struct();
    end
    params_ICURC = SetDefaultParams_ICURC(params_ICURC);
    
    eta=params_ICURC.eta; 
    fprintf('using stepsize eta_C = %f.\n', eta(1));
    fprintf('using stepsize eta_R = %f.\n', eta(2));
    fprintf('using stepsize eta_U = %f.\n', eta(3));
    
    TOL=params_ICURC.TOL;
    max_ite=params_ICURC.max_ite;
    steps_are1=params_ICURC.steps_are1;
    
    %This step is to extract observed C, U, and R 
    Obs_U = X_Omega_UR(I_ccs, J_ccs);  
    Obs_C = X_Omega_UR(:, J_ccs);
    Obs_R = X_Omega_UR(I_ccs, :);
    C_size = size(Obs_C);
    R_size = size(Obs_R);

    all_row_ind = 1:C_size(1); 
    all_col_ind = 1:R_size(2); 

    I_ccs_comp = setdiff(all_row_ind, I_ccs); 
    J_ccs_comp = setdiff(all_col_ind, J_ccs); 
    
    C = Obs_C(I_ccs_comp, :); 
    R = Obs_R(:, J_ccs_comp); 
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
        fprintf('running ICURC with stepsize all equal to 1\n');
        fct_time = tic;
        for ICURC_ite = 1:max_ite
            ite_itme = tic;
            R = u*(u'*R);
            C = C*v*v';        
            Old_error = New_Error;
            
            New_Error = (norm(R(Omega_row) - L_obs_row_vec,'fro') + ... 
                            norm(C(Omega_col)- L_obs_col_vec, 'fro') + ... 
                            norm(U_i(Omega_U)- L_obs_U_vec))/col_row_norm_sum;  
            if New_Error < TOL || ICURC_ite == max_ite 
                ICURC_time = toc(fct_time); 
                R(Omega_row) = L_obs_row_vec;
                C(Omega_col) = L_obs_col_vec;
                U_i(Omega_U) = L_obs_U_vec; 

                Final_C = zeros(C_size);
                Final_R = zeros(R_size); 
                Final_C(I_ccs_comp, :) = C;
                Final_R(:, J_ccs_comp) = R;
                Final_C(I_ccs, :) = U_i; 
                Final_R(:, J_ccs) = U_i; 

                C = Final_C;
                R = Final_R; 
                
                [u,s,v] = svd(U_i);  
                u = u(:, 1:r);
                v = v(:, 1:r);
                s = s(1:r, 1:r);  
                U_pinv = v*pinv(s)*u';
                fprintf('ICURC finished at  %d-th iteration \n', ICURC_ite);
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
            
            fprintf('Iteration %d: error: %e, timer: %f \n', ICURC_ite, New_Error, toc(ite_itme));
        end
        
    %If all step sizes are not 1    
    else
        fprintf('running ICURC with given stepsizes \n');
        step_c = eta(1);
        step_r = eta(2);
        step_u = eta(3);
        
        fct_time = tic;
        for ICURC_ite = 1:max_ite
            ite_itme = tic;
            R = u*(u'*R);
            C = C*v*v';        
            Old_error = New_Error;
            New_Error = (norm(R(Omega_row) - L_obs_row_vec,'fro') + ... 
                            norm(C(Omega_col)- L_obs_col_vec, 'fro') + ... 
                            norm(U_i(Omega_U)- L_obs_U_vec))/col_row_norm_sum;   
                        
            if New_Error < TOL || ICURC_ite == max_ite
                ICURC_time = toc(fct_time); 

                R(Omega_row) = L_obs_row_vec;
                C(Omega_col) = L_obs_col_vec;
                U_i(Omega_U) = L_obs_U_vec; 

                Final_C = zeros(C_size);
                Final_R = zeros(R_size); 
                Final_C(I_ccs_comp, :) = C;
                Final_R(:, J_ccs_comp) = R;
                Final_C(I_ccs, :) = U_i; 
                Final_R(:, J_ccs) = U_i; 

                C = Final_C;
                R = Final_R; 

                [u,s,v] = svd(U_i);  
                u = u(:, 1:r);
                v = v(:, 1:r);
                s = s(1:r, 1:r);  
                U_pinv = v*pinv(s)*u';
                fprintf('ICURC finished in %d iterations, final error: %e, total runtime: %f \n', ICURC_ite, New_Error, ICURC_time);
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
            fprintf('Iteration %d: error: %e, timer: %f \n', ICURC_ite, New_Error, toc(ite_itme));

         end
     end
end
