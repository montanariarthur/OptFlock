function lambdamax = opteigreal_OShom(b,N,dim,Lap_alpha,Lap_beta,param)
%   Computes the Lyapunov exponent of the OS flocking model with
%   homogeneous agents.
%
%   b                     -  optimization parameter (2-dimensional)
%   N, dim, param         -  system parameters
%   Lap_alpha, Lap_beta   -  network structures (alpha-alpha, alpha-beta)

% Parameters
n = dim*N;

% Jacobian matrix
Jac11 = zeros(n,n);
Jac12 = eye(n);

Jac21_11 = - b(1)*eye(N);
Jac21 = kron(eye(dim), Jac21_11);

% Checks if there are obstacles in the system
if param.Nbeta == 0
    Jac22_11 = - (b(2)*eye(N) - param.C2alpha*Lap_alpha);
    Jac22 = kron(eye(dim), Jac22_11);
else
    Jac22_11 = - (b(2)*eye(N) - param.C2alpha*Lap_alpha);
    Jac22 = kron(eye(dim), Jac22_11) - param.C2beta*Lap_beta;
end


Jac = [Jac11, Jac12; Jac21, Jac22];

% Eigenvalue
lambdamax = max(real(eig(Jac)));  % sort eigenvalues in descending order

end