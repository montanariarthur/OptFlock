function  [S,RMSE,RMSE_pos] = agentcoord_CM(x,dim,N,pg)
% Computes, for each instant of time k, the corresponding coordinates of
% the agents 'S' and the RMSE of position and velocity 'RMSE'.
%
%  x      -  state vector
%  dim    -  Euclidean space dimension
%  N      -  number of agents
%  pg     -  relative desired position in the formation

% Determines x,y,z coordinates of agents (position and velocity)
if dim == 2
    % Position
    S.px = x(1:N,:);
    S.py = x(N+1:2*N,:);
    S.p = [S.px; S.py];
    % Velocity
    S.qx = x(2*N+1:3*N,:);
    S.qy = x(3*N+1:4*N,:);
    S.q = [S.qx; S.qy];

    % Reference position
    S.pref_x = x(2*N*dim+1,:);
    S.pref_y = x(2*N*dim+2,:);
    S.pref = [S.pref_x; S.pref_y];
    % Reference velocity
    S.qref_x = x(2*N*dim+3,:);
    S.qref_y = x(2*N*dim+4,:);
    S.qref = [S.qref_x; S.qref_y];
    
    % Position difference
    Gcomplete = ones(N);    % complete graph
    Lcomplete = Gcomplete - diag(sum(Gcomplete,2));
    S.pdiffx = Lcomplete*S.px;
end

% Root mean square error of position and velocity
for k = 1:length(x)
    aux = [kron(S.pref(:,k),ones(N,1)); kron(S.qref(:,k),ones(N,1))];
    RMSE(k,1) = norm(x(1:2*dim*N,k) - aux - [pg; zeros(dim*N,1)] );
    RMSE_pos(k,1) = norm(x(1:dim*N,k) - aux(1:dim*N,1) - [pg] );
end

end