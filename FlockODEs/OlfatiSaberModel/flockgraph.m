function [G,Q] = flockgraph(q,param)
% Equation (4) of Reference [2].

N = param.N;
dim = param.dim;
range = param.range;

% Matrix Q(i,j) = qj_k - qi_k distance along each dimension k
for k = 1:dim
    qk = q(N*(k-1) + 1 : N*k , 1);
    Qk{k} = ones(N,1)*qk' - qk*ones(1,N);
end

% Distance matrix Q(i,j) = norm(qj - qi)
Q = zeros(N,N);
for k = 1:dim
    Q = Q + Qk{k}.^2;
end
Q = sqrt(Q);

% Undirected, binary graph
G = (Q < range);
G = G - eye(N);

end