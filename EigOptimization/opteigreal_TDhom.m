function lambdamax = opteigreal_TDhom(b,L,tau,n)
%   Computes the Lyapunov exponent of the time-delay consensus model.
%   b                         -   optimization parameter
%   L, tau, n                 -   system parameters
%   J1, J2                    -   Jacobians in the generalized eigenvalue problem 
%   lambdamax                 -   maximum eigenvalue

if length(b) == 1
    beta = b;
    J1 = [zeros(n,n), eye(n); zeros(n,n), zeros(n,n)];
    J2 = [zeros(n,n), zeros(n,n); -beta*L , -beta*L];
elseif length(b) == 2
    alpha = b(1:n);
    beta = b(n+1:2*n);
    J1 = [zeros(n,n), eye(n); zeros(n,n), zeros(n,n)];
    J2 = [zeros(n,n), zeros(n,n);-alpha*L , -beta*L];
end

% Eigenvalue
lambdamax = dde_rightmost_eig(J1,J2,tau);  % sort eigenvalues in descending order

end

