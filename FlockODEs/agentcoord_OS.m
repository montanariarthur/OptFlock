function  [S,Edev,Conn,Vmis,Terror] = agentcoord_OS(x,dim,N,param,knob)
% Computes, for each instant of time k, the corresponding coordinates of
% the agents 'S', the lattice deviation energy 'Edev', the relative network
% connectivity 'Conn, the velocity mismatch 'Vmis', and the centering error
% with respect to the target 'Terror'
%
%  x      -  state vector
%  dim    -  Euclidean space dimension
%  N      -  number of agents
%  param  -  system parameters
%  knob   -  knob controlling internal routines
%     knob.protocol - implements 'Algorithm 2' or 'Algorithm 3' of Ref. [2]

% Determines x,y,z coordinates of agents (position and velocity)
q(:,1) = x(1:N*dim);
p(:,1) = x(N*dim+1:2*N*dim);

% Agent positions
switch dim
    case 2
        % Position
        S.qx = x(1:N,:);
        S.qy = x(N+1:2*N,:);
        S.q = [S.qx; S.qy];
        % Velocity
        S.px = x(2*N+1:3*N,:);
        S.py = x(3*N+1:4*N,:);
        S.p = [S.px; S.py];

        switch knob.protocol
            case {'Algorithm 2','Algorithm 3'}
                % Reference position
                S.qref_x = x(2*N*dim+1,:);
                S.qref_y = x(2*N*dim+2,:);
                S.qref = [S.qref_x; S.qref_y];
                % Reference velocity
                S.pref_x = x(2*N*dim+3,:);
                S.pref_y = x(2*N*dim+4,:);
                S.pref = [S.pref_x; S.pref_y];
        end


    case 3
        % Position
        S.qx = x(1:N,:);
        S.qy = x(N+1:2*N,:);
        S.qz = x(2*N+1:3*N,:);
        S.q = [S.qx; S.qy; S.qz];
        % Velocity
        S.px = x(3*N+1:4*N,:);
        S.py = x(4*N+1:5*N,:);
        S.pz = x(5*N+1:6*N,:)
        S.p = [S.px; S.py; S.pz];

        switch knob.protocol
            case {'Algorithm 2','Algorithm 3'}
                % Reference position
                S.qref_x = x(2*N*dim+1,:);
                S.qref_y = x(2*N*dim+2,:);
                S.qref_z = x(2*N*dim+3,:);
                S.qref = [S.qref_x; S.qref_y; S.qref_z];
                % Reference velocity
                S.pref_x = x(2*N*dim+4,:);
                S.pref_y = x(2*N*dim+5,:);
                S.pref_z = x(2*N*dim+6,:);
                S.pref = [S.pref_x; S.pref_y; S.pref_z];
        end
end

% Measures
param.h = param.h_alpha;
for k = 1:length(x)

    % Deviation energy (Eq. 7 of Ref [2])
    E = 0;
    Nedges = 0;
    for i = 1:N
        qi = S.q(i:N:i+N*(dim-1),k);
        for j = 1:N
            qj = S.q(j:N:j+N*(dim-1),k);

            distance = norm(qi - qj);
            if distance < param.range && i ~= j
                E = E + (norm(qj - qi) - param.dist)^2;
                Nedges = Nedges + 1;
                Adj_k(i,j) = spatial_adj_entry(qi,qj,param);
            else
                Adj_k(i,j) = 0;
            end
        end
    end
    Edev(k,1) = (E/(Nedges + 1))/(param.dist^2);
    
    % Relative connectivity (page 418 of Ref. [2])
    Lap_k = Adj_k - diag(sum(Adj_k,2));
    Conn(k) = (1/(N-1))*rank(Lap_k);
end

% Velocity mismatch
Vmis = vecnorm(S.p - kron(S.pref,ones(N,1)),2);

% Centering error with respect to the target position
Terror = vecnorm([mean(S.qx,1); mean(S.qy,1)] - S.qref,2);