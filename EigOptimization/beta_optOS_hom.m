function [b_opt,lyapexp] = beta_optOS_hom(Lap_alpha,Lap_beta,param,b0,blim)
%   Finds the optimal (homogeneous) parameter 'b_opt' that minimizes the
%   Lyapunov exponent 'lyapexp' of the OS flocking model.
%
%   Lap_alpha, Lap_beta   -  network structures (alpha-alpha, alpha-beta)
%   param                 -  system parameters 
%
%   b0          -  initial guess for the optimization (same size as b)
%   blim        -  upper bound for b
%   Nopt        -  number of realizations for the optimization procedure

% Parameters
N = param.N;
dim = param.dim;

% Optimization routine
b_opt = fmincon(@(b)opteigreal_OShom(b,N,dim,Lap_alpha,Lap_beta,param),b0,[],[],[],[],0*ones(2,1),[blim(1)*ones(N,1); blim(2)*ones(N,1)]);
lyapexp = opteigreal_OShom(b_opt,N,dim,Lap_alpha,Lap_beta,param);

end