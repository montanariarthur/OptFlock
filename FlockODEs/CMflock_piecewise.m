function xdot = CMflock_piecewise(t,x,N,dim,B1,B2,gamma,pg,mass,noise,Lap,reference_trajectory)
% ODE of the "centroid model" of flocks.
% In this model, agents are tasked to main a pre-specified formation while
% tracking a target in space.
%
%  x      -  state vector
%  N      -  number of agents
%  dim    -  Euclidean space dimension
%  B1     -  position feedback matrix
%  B2     -  velocity feedback matrix
%  gamma  -  tuning parameter
%  pg     -  intended relative positions for each agent in the formation
%  mass   -  agent mass
%  noise  -  standard deviation of the noise
%  Lap    -  Laplacian matrix
%  reference_trajectory - five trajectories are implemented for the target 
%                         (stationary, constant speed, sinusoidal/zigzag, 
%                         circular, fuzzy/random walk)

% State vector
p(:,1) = x(1:N*dim);
q(:,1) = x(N*dim+1:2*N*dim);
pref(:,1) = x(2*N*dim+1:2*N*dim+dim);
qref(:,1) = x(2*N*dim+dim+1:2*N*dim+2*dim);

% Reference trajectory
switch reference_trajectory
    case 'stationary'
        prefdot = zeros(dim,1);
        qrefdot = zeros(dim,1);
    case 'constant'
        prefdot = qref;
        qrefdot = zeros(dim,1);
    case 'zigzag'
        omega = 1;                           
        radius = 10;
        prefdot = [qref(1); radius*sin(omega*t)]; 
        qrefdot = [0;       radius*cos(omega*t)*omega];
    case 'circular' 
        omega  = 2;
        lambda = 1;
        radius = 10;
        prefdot = [- omega*pref(2) - lambda*(pref(1)^2 + pref(2)^2 - radius^2)*pref(1);
                   + omega*pref(1) - lambda*(pref(1)^2 + pref(2)^2 - radius^2)*pref(2)];
        qrefdot = [lambda*qref(1)*(radius^2 - 3*pref(1)^2 - pref(2)^2) - qref(2)*(+ 2*lambda*pref(1)*pref(2) + omega);
                   lambda*qref(2)*(radius^2 - 3*pref(2)^2 - pref(1)^2) + qref(1)*(- 2*lambda*pref(1)*pref(2) + omega)];
    case 'fuzzy'
        sigma = 30;
        gamma = 0.1;
        m = 0.1;
        prefdot = qref;
        qrefdot = (1/m)*(- gamma*qref + sigma*randn(2,1));
end

% ODEs of agents
pdot = q;
qdot = (1./mass).*(- kron(eye(dim),B1)*(p - kron(pref,ones(N,1)) - pg) - gamma*kron(eye(dim),B2)*(q - kron(qref,ones(N,1))) ...  % 1st term (position consensus)
       + kron(eye(dim),Lap)*((p - pg) + gamma*q) + kron(qrefdot,ones(N,1))  ...                                                  % 2nd term (velocity consensus)
       + 0*kron(qrefdot,ones(N,1))) + noise*randn(N*dim,1);                                                                      % target acceleration + noise

% Output
xdot = [pdot; qdot; prefdot; qrefdot];

end