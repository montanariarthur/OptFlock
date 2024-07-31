function mtle = dde_rightmost_eig(J1,J2,tau)
%   Computes the maximum eigenvalue of the generalized eigenvalue problem.
%   J1, J2                    -   Jacobians in the generalized eigenvalue problem 
%   tau                       -   Time delay
%   mtle                      -   Maximum eigenvalue

if tau~=0
    flag_newhheur = 1;
    method = df_mthod('stst',flag_newhheur);
    method.stability.minimal_real_part = -10;
    
    stability = stst_stabil_mod(J1,J2,tau,method);
    
    mtle = unique(max(real(stability.l1)));
    if abs(mtle)<1e-6
        lamIndex = find(abs(stability.l1)>10^-6 & imag(stability.l1)<=0);
        mtle = unique(max(real(stability.l1(lamIndex))));
    end
    if any(size(mtle) == 0)
        mtle = unique(max(real(stability.l1)));
        if any(size(mtle) == 0)
            disp("problem")
            mtle = 0;
        end
    end
end

if tau == 0 
    e = eig(J1+J2);
    e_max = max(real(e));
    e_max_index = find(real(e)==e_max);
    if abs(e_max)<1e-4
        e(e_max_index) = [];
        e_max = max(real(e));
        e_max_index = find(real(e)==e_max & imag(e)>=0);
    end
    e_lambda = e(e_max_index);
    mtle = unique(real(e_lambda));
end
