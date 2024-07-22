function [b_opt,lyapexp] = beta_optOS_het(Lap_alpha,Lap_beta,param,b0,blim,Nopt)
%   Finds the optimal (heterogeneous) parameter 'b_opt' that minimizes the
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

% Optimization routine over independent realizations
for MCopt = 1:Nopt
    B_het_mc(:,MCopt) = fmincon(@(b)opteigreal_OShet(b,N,dim,Lap_alpha,Lap_beta,param),[b0(1)+0.1*randn(N,1); b0(2)+0.1*randn(N,1)],[],[],[],[],0*ones(2*N,1),[blim(1)*ones(N,1); blim(2)*ones(N,1)]);
    lambda_het_mc(MCopt) = opteigreal_OShet(B_het_mc(:,MCopt),N,dim,Lap_alpha,Lap_beta,param);
end

% Picks the best heterogeneous solution
[~,minindex] = min(lambda_het_mc);
lyapexp = lambda_het_mc(minindex);
b_opt = B_het_mc(:,minindex);

end