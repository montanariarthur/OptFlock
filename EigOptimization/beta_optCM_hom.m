function [b_opt,lyapexp] = beta_optCM_hom(Adj,param,optB,b1,b2,b0,blim)
%   Finds the optimal (homogeneous) parameter 'b_opt' that minimizes the
%   Lyapunov exponent 'lyapexp' of the CM flocking model.
%
%   Adj         -  adjacency matrix
%   param       -  system parameters 
%
%   optB        -  defines the optimization parameter b
%        'argmax = B1' -> position feedback gain, b is scalar
%        'argmax = B2' -> velocity feedback gain, b is scalar
%        'argmax = B'  -> position and velocity gain 
%                         (constrained to be equal, b is scalar)
%    'argmax = B1, B2' -> position and velocity gains 
%                         (independent, b is 2-dimensional)
%
%   b1, b2      -  constant feedback gain (scalar)
%   (e.g., for cases in which only b1 is optimized we need to define b2)
%
%   b0          -  initial guess for the optimization (same size as b)
%   blim        -  upper bound for b
%   Nopt        -  number of realizations for the optimization procedure

% Parameters
N = param.N;
gamma = param.gamma;
mass = param.mass;
dim = param.dim;

% Solves optimization problem 
switch optB
    
    % For this special case, the optimization is solved analytically
    case 'argmax = B'
        D = diag(sum(Adj));
        L = D - Adj;
        lambda = sort(eig(L));
        b_opt = (2*mass)/gamma^2 - lambda(N) + sqrt(lambda(N)^2 + (4*mass^2)/gamma^4);
        lyapexp = -gamma*b_opt/(2*mass);

    % Optimization routine using fmincon (constrained optimization)
    otherwise
        D = diag(sum(Adj));
        Lap = Adj - D;
        b_opt = fmincon(@(b)opteigreal_CMhom(b,Lap,gamma,N,dim,mass,optB,b1,b2),b0,[],[],[],[],0*ones(size(b0)),blim*ones(size(b0)));
        lyapexp = opteigreal_CMhom(b_opt,Lap,gamma,N,dim,mass,optB,b1,b2);
end

end