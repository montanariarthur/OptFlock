function Aij = spatial_adj_entry(qi,qj,param)
% Equation (11) of Reference [2].

range = param.range;
r_alpha = sigma_norm(range,param);

Aij = bump_function( sigma_norm(qj - qi,param) / r_alpha , param);

end