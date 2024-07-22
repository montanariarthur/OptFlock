function xdot = OSflock(t,x,C1alpha,C2alpha,C1gamma,C2gamma,C1beta,C2beta,qbeta,param,knob)
% ODE of Olfati-Saber model of flocks. By convention, flock agents are
% called alpha-agents, the virtual target is called a gamma-agent, and
% obstacles are called beta-agents.
%
%  x                          -  state vector
%  C1alpha, C1gamma, C1beta   -  gain of the gradient term for alpha-alpha,
%                                alpha-gamma, & alpha-beta interaction
%  C2alpha, C2gamma, C2beta   -  gain of the velocity consensus
%  qbeta                      -  obstacle positions
%  param                      -  system parameters
%  knob                       -  knob controlling internal routines
%     knob.protocol - implements 'Algorithm 2' or 'Algorithm 3' of Ref. [2]
%     knob.ref_traj - target trajectory ('stationary' or 'constant')

%% Parameters
N = param.N;
dim = param.dim;
Nbeta = param.Nbeta;

% State vector (alpha-agent)
q(:,1) = x(1:N*dim);
p(:,1) = x(N*dim+1:2*N*dim);

% Reference trajectory (gamma-agent)
switch knob.protocol
    case {'Algorithm 2','Algorithm 3'}
        % State vector (gamma-agent)
        qref(:,1) = x(2*N*dim+1:2*N*dim+dim);
        pref(:,1) = x(2*N*dim+dim+1:2*N*dim+2*dim);
        switch knob.ref_traj
            case 'stationary'
                qrefdot = zeros(dim,1);
                prefdot = zeros(dim,1);
            case 'constant'
                qrefdot = pref;
                prefdot = zeros(dim,1);
        end
end

%% For each alpha-agent in the network

u_alpha = zeros(N*dim,1);
u_beta = zeros(N*dim,1);
for i = 1:N
    qi = q(i:N:i+N*(dim-1));
    pi = p(i:N:i+N*(dim-1));
    param.h = param.h_alpha;

    % Interaction term through each neighbor (alpha-alpha)
    for j = 1:N
        
        % except the node itself
        if i~=j                  
            qj = q(j:N:j+N*(dim-1));
            pj = p(j:N:j+N*(dim-1));

            % if it is within the interaction range
            distance = norm(qi - qj);
            if distance < param.range
                u_alpha(i:N:i+N*(dim-1),1) = u_alpha(i:N:i+N*(dim-1),1) ...
                    + C1alpha*phi_alpha_function(sigma_norm(qj-qi,param),param)*n_entry(qi,qj,param) ...
                    + C2alpha*spatial_adj_entry(qi,qj,param)*(pj - pi);
            end
        end
    end

    % Obstacle avoidance (alpha-beta)
    if Nbeta > 0
        param.h = param.h_beta;

        for k = 1:Nbeta
            % obstacle center and radius
            yk = qbeta(k:Nbeta:k+Nbeta*(dim-1),1);
            Rk = param.radius_beta(k);

            % calculates position and velocity according to Sec. VII-D
            ak = (qi - yk)/norm(qi - yk);
            mu = Rk / norm(qi - yk);
            P = eye(dim) - ak*ak';
            qk = mu*qi + (1 - mu)*yk;

            % if it is within the interaction range
            distance_beta = norm(qi - qk);
            if distance_beta < param.range_beta
                pk = mu*P*pi;
                u_beta(i:N:i+N*(dim-1),1) = u_beta(i:N:i+N*(dim-1),1) ...
                    + C1beta*phi_beta_function(sigma_norm(qk-qi,param),param)*n_entry(qi,qk,param) ... 
                    + C2beta*spatial_adj_beta(qi,qk,param)*(pk - pi);
            end
        end
    end
end

% Target tracking (alpha-gamma)
switch knob.protocol
    case {'Algorithm 2','Algorithm 3'}
        u_gamma = - kron(eye(dim),C1gamma)*(q - kron(qref,ones(N,1))) - kron(eye(dim),C2gamma)*(p - kron(pref,ones(N,1)));
end


%% Output
switch knob.protocol
    case 'Algorithm 1'
        qdot = p;
        pdot = u_alpha;
        xdot = [qdot; pdot];
    case {'Algorithm 2','Algorithm 3'}
        qdot = p;
        pdot = u_alpha + u_beta + u_gamma;
        xdot = [qdot; pdot; qrefdot; prefdot];
end


end