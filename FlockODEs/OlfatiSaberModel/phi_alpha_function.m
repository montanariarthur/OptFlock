function phi_alpha = phi_alpha_function(z,param)
% Equation (15) of Reference [2].

dist = param.dist;
range = param.range;

d_alpha = sigma_norm(dist,param);
r_alpha = sigma_norm(range,param);
phi_alpha = bump_function(z/r_alpha,param) * phi_function(z - d_alpha,param);

end