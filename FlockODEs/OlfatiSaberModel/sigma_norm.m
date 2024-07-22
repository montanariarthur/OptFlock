function sigma = sigma_norm(z,param)
% Equation (8) of Reference [2].

epsilon = param.epsilon;

sigma = (1/epsilon) * ( sqrt(1 + epsilon*norm(z)^2) - 1 );

end