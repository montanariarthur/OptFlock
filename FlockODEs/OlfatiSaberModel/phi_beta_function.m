function phi_beta = phi_beta_function(z,param)
% Equation (15) of Reference [2].

dist_beta = param.dist_beta;

d_beta = sigma_norm(dist_beta,param);
arg = z - d_beta;
sigma1 = arg/sqrt(1 + arg^2);
phi_beta = bump_function(z/d_beta,param) * (sigma1 - 1);

end