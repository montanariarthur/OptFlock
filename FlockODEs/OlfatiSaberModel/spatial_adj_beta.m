function Bij = spatial_adj_beta(qi,qj,param)
% Equation (55) of Reference [2].

dist_beta = param.dist_beta;
d_beta = sigma_norm(dist_beta,param);

Bij = bump_function( sigma_norm(qj - qi,param) / d_beta , param);

end