function [b_opt,lyapexp] = beta_optCM_het(Adj,param,optB,b1,b2,b0,blim,Nopt)
%   Finds the optimal (heterogeneous) parameter 'b_opt' that minimizes the
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

% Adjacency matrix
D = diag(sum(Adj));
Lap = Adj - D;

% Solves optimization problem 
switch optB

    case 'argmax = B1, B2'
        B_het_mc = zeros(2*N,Nopt);
        lambda_het_mc = zeros(Nopt,1);
        for MCopt = 1:Nopt      % for over independent realizations
            B_het_mc(:,MCopt) = fmincon(@(b)opteigreal_CMhet(b,Lap,gamma,N,dim,mass,optB,b1,b2),[b0(1)+0.1*randn(N,1); b0(2)+0.1*randn(N,1)],[],[],[],[],0*ones(2*N,1),blim*ones(2*N,1));
            lambda_het_mc(MCopt) = opteigreal_CMhet(B_het_mc(:,MCopt),Lap,gamma,N,dim,mass,optB,b1,b2);
        end

    otherwise
        B_het_mc = zeros(N,Nopt);
        lambda_het_mc = zeros(Nopt,1);
        for MCopt = 1:Nopt      % for over independent realizations
            B_het_mc(:,MCopt) = fmincon(@(b)opteigreal_CMhet(b,Lap,gamma,N,dim,mass,optB,b1,b2),b0+0.1*randn(N,1),[],[],[],[],0*ones(N,1),blim*ones(N,1));
            lambda_het_mc(MCopt) = opteigreal_CMhet(B_het_mc(:,MCopt),Lap,gamma,N,dim,mass,optB,b1,b2);
        end     
end

% Picks the best heterogeneous solution
[~,minindex] = min(lambda_het_mc);
lyapexp = lambda_het_mc(minindex);
b_opt = B_het_mc(:,minindex);

end