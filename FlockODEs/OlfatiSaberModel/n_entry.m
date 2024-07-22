function nij = n_entry(qi,qj,param)
% Equation (11) of Reference [2].

nij = sigma_grad(qj - qi,param);

end