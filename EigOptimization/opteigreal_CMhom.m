function lambdamax = opteigreal_CMhom(b,Lap,gamma,N,dim,mass,optB,b1,b2)
%   Computes the Lyapunov exponent of the CM flocking model with
%   homogeneous agents.
%
%   b                         -  optimization parameter (homogeneous)
%   Lap, gamma, N, dim, mass  -  system parameters
%
%   optB                      -  defines the optimization parameter
%        'argmax = B1' -> position feedback gain, b is scalar
%        'argmax = B2' -> velocity feedback gain, b is scalar
%        'argmax = B'  -> position and velocity gain 
%                         (constrained to be equal, b is scalar)
%    'argmax = B1, B2' -> position and velocity gains 
%                         (independent, b is 2-dimensional)
%
%   b1, b2                    -  constant feedback gain (scalar)
%   (e.g., for cases in which only b1 is optimized we need to define b2)

% Parameters
n = dim*N;

% Jacobian
Jac11 = zeros(n,n);
Jac12 = eye(n);

switch optB
    case 'argmax = B1'
        Jac21_11 = -(b*eye(N) - Lap);
        Jac22_11 = -gamma*(b2*eye(N) - Lap);

    case 'argmax = B2'
        Jac21_11 = -(b1*eye(N) - Lap);
        Jac22_11 = -gamma*(b*eye(N) - Lap);

    case 'argmax = B'
        Jac21_11 = -(b*eye(N) - Lap);
        Jac22_11 = -gamma*(b*eye(N) - Lap);

    case 'argmax = B1, B2'
        Jac21_11 = -(b(1)*eye(N) - Lap);
        Jac22_11 = -gamma*(b(2)*eye(N) - Lap);
end

Jac21 = kron(eye(dim), Jac21_11);
Jac22 = kron(eye(dim), Jac22_11);

Jac = [Jac11, Jac12; (1/mass)*Jac21, (1/mass)*Jac22];

% Eigenvalue
lambdamax = max(real(eig(Jac)));  % sort eigenvalues in descending order

end