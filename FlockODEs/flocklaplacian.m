function [Lap,Adj] = flocklaplacian(p,sigma,beta,K,N)
% Laplacian matrix 'Lap' and adjacency matrix 'Adj' of the CM flocking.
%
% p               -  agents' positions
% N               -  number of agents
% beta, K, sigma  -  parameters

for i = 1:N
    for j = 1:N
        if i == j
            Adj(i,j) = 0;
        else
            pdiff = p(i:N:end) - p(j:N:end);
            Adj(i,j) = K / ( (sigma^2 + pdiff'*pdiff)^beta );
        end
    end
end
Lap = Adj - diag(sum(Adj,2));