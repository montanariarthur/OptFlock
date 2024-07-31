function lambdamax = opteigreal_CSDhet(b,L,tau,n)
%   Computes the Lyapunov exponent of the time-delay consensus model.
%   b                         -   nx1 (heterogeneous) optimization parameter
%   L, tau, n                 -   system parameters
%   J1, J2                    -   Jacobians in the generalized eigenvalue problem 
%   lambdamax                 -   maximum eigenvalue

if length(b) == n
    beta = b;
    J1 = [zeros(n,n), eye(n); zeros(n,n), zeros(n,n)];
    J2 = [zeros(n,n), zeros(n,n); -L .* beta, -L .* beta];
elseif length(b) == 2*n
    alpha = b(1:n);
    beta = b(n+1:2*n);
    J1 = [zeros(n,n), eye(n); zeros(n,n), zeros(n,n)];
    J2 = [zeros(n,n), zeros(n,n);-L .* alpha, -L .* beta];
end

% Eigenvalue
lambdamax = dde_rightmost_eig(J1,J2,tau);  % sort eigenvalues in descending order

end
