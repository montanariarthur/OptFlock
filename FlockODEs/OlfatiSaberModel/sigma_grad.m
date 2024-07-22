function sigma_eps = sigma_grad(z,param)
% Equation (9) of Reference [2].

epsilon = param.epsilon;

sigma_eps = z / (1 + epsilon*sigma_norm(z,param));

end