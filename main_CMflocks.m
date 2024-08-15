%% Real-time optimization of the flocking centroid model (CM)
% This code compares the performance of a flock of agents optimized in real
% time (over regular optimization windows) with:
%  - constrained parameters (homogeneous flock)
%  - unconstrained parameters (heterogeneous flock)

clear all; close all; clc;
addpath([pwd,'/EigOptimization/'])
addpath([pwd,'/FlockODEs/'])

%% Setup simulation

% Parameters to be optimized
knob.optB = 'argmax = B1, B2';
    % argmax = B1     - only position feedback gains are optimized
    % argmax = B2     - only velocity feedback gains are optimized
    % argmax = B      - both feedback gains are constrained to be equal
    % argmax = B1, B2 - both feedback gains are optimized independently

% Desired formation shape
knob.formation = 'sun';
    % 'random' -  uniformly random in a box 
    % 'circle' -  single circle with fixed radius
    % 'sun'    -  circular patterns with increasing radius

% Update rule for the adjacency matrix ('continuous','piecewise')
knob.adjupdate = 'piecewise';

% Target trajectory
knob.ref_traj = 'constant';
    % 'stationary' - does not move
    % 'constant'   - moves with constant speed along x-axis
    % 'zigzag'     - sinusoidal pattern along x-axis
    % 'circular'   - moves in circles
    % 'fuzzy'      - random walk

disp(knob)

%% Parameters

% System parameters
N = 30;        param.N = N;         % number of agents
dim = 2;       param.dim = dim;     % dimension of the Euclidean space
beta = 0.8;         % exponent of the distance function in Adj
sigma = 0.1;        % repulsion force
K = 2;              % coupling strength
gamma = 1;     param.gamma = gamma;      
mass = 1;      param.mass = mass;
noise = 0.1;

% Damping coefficients
b1 = 5;              % position-integrator feedback gain
b2 = 5;              % velocity-integrator feedback gain
blim = Inf;          % upper bound for feedback gain

% Initial conditions
p0 = unifrnd(-2,2,N*dim,1);        % initial position
q0 = unifrnd(0,0,N*dim,1);         % initial velocity
switch knob.ref_traj               % initial reference pos/velocity
    case 'stationary'
        pref0 = 100*[1; 1];
        qref0 = 0*[1; 1];
    case 'constant'
        pref0 = 5000*[1; 1];
        qref0 = 100*[1; 0];
    case 'zigzag'
        pref0 = 100*[1; 1];
        qref0 = 10*[1; 1];
    case 'circular'
        pref0 = 10*[1; 0];
        qref0 = 10*[0; 1];
    case 'fuzzy'
        pref0 = 0*[1; 1];
        qref0 = 0*[1; 1];
end
x0 = [p0;q0;pref0;qref0];          % state vector

% Relative positions
switch knob.formation
    case 'random'
    % uniform random distribution
        pg_x = unifrnd(-5,5,N,1);
        pg_y = unifrnd(-5,5,N,1);
    case 'circle'
        % (circle of radius R around reference point, evenly spaced agents)
        radius = 5;
        angle = linspace(0,2*pi-2*pi/N,N)';
        pg_x = [radius*cos(angle)];
        pg_y = [radius*sin(angle)];
    case 'sun'
        radius = 1; angle = linspace(0,2*pi-2*pi/5,5)';
        pg_x1 = [radius*cos(angle)];
        pg_y1 = [radius*sin(angle)];

        radius = 3; angle = linspace(0,2*pi-2*pi/10,10)';
        pg_x2 = [radius*cos(angle)];
        pg_y2 = [radius*sin(angle)];

        radius = 5; angle = linspace(0,2*pi-2*pi/15,15)';
        pg_x3 = [radius*cos(angle)];
        pg_y3 = [radius*sin(angle)];

        pg_x = [pg_x1; pg_x2; pg_x3];
        pg_y = [pg_y1; pg_y2; pg_y3];
end
if dim == 2
    pg = [pg_x; pg_y];
else
    pg = unifrnd(-5,5,dim*N,1);
end
% figure(1); scatter(pg_x,pg_y)

% Simulation time
dt = 1e-3;
tf = 30;

% Real-time optimization
tstep = 1;
Nsteps = tf/tstep;
Nopt = 10;                      % MC iterations for optimization

%% Online simulation (homogeneous optimization)
t_hom = 0; x_hom = x0;
B1 = b1*eye(N); B2 = b2*eye(N);
for k = 1:Nsteps
    % Adjacency matrix computed at time k
    p_k = x_hom(1:N,end);
    [Lap_k,Adj_k] = flocklaplacian(p_k,sigma,beta,K,N);

    % Eigenvalue optimization
    switch knob.optB
        case 'argmax = B1'
            [B_hom(:,k),lambda_hom(k)] = beta_optCM_hom(Adj_k,param,knob.optB,b1,b2,B1(1,1),blim);
            B1 = B_hom(k)*eye(N);
            B2 = b2*eye(N);
        case 'argmax = B2'
            [B_hom(:,k),lambda_hom(k)] = beta_optCM_hom(Adj_k,param,knob.optB,b1,b2,B2(1,1),blim);
            B1 = b1*eye(N);
            B2 = B_hom(k)*eye(N);
        case 'argmax = B'
            [B_hom(:,k),lambda_hom(k)] = beta_optCM_hom(Adj_k,param,knob.optB,b1,b2,B1(1,1),blim);
            B1 = B_hom(k)*eye(N);
            B2 = B_hom(k)*eye(N);
        case 'argmax = B1, B2'
            [B_hom(:,k),lambda_hom(k)] = beta_optCM_hom(Adj_k,param,knob.optB,b1,b2,[B1(1,1); B2(1,1)],blim);
            B1 = B_hom(1,k)*eye(N);
            B2 = B_hom(2,k)*eye(N);
    end

    % ODE integration
    switch knob.adjupdate
        case 'continuous'
            [ttemp,xtemp] = odeRK(@(t,x)CMflock_continuous(t,x,N,dim,beta,sigma,K,B1,B2,gamma,pg,mass,noise,knob.ref_traj),[t_hom(end) dt t_hom(end)+tstep],x_hom(:,end)');
        case 'piecewise'
            [ttemp,xtemp] = odeRK(@(t,x)CMflock_piecewise(t,x,N,dim,B1,B2,gamma,pg,mass,noise,Lap_k,knob.ref_traj),[t_hom(end) dt t_hom(end)+tstep],x_hom(:,end)');
    end
    xtemp = xtemp';
    t_hom = [t_hom; ttemp(2:end)'];
    x_hom = [x_hom  xtemp(:,2:end)];
    xtemp = []; ttemp = [];
end

%% Online simulation (heterogeneous optimization, BFGS method)
t_het = 0; x_het = x0;
B1 = b1*eye(N); B2 = b2*eye(N);

for k = 1:Nsteps
    k
    % Adjacency matrix computed at time k
    p_k = x_het(1:N,end);
    [Lap_k,Adj_k] = flocklaplacian(p_k,sigma,beta,K,N);

    % Eigenvalue optimization
    switch knob.optB
        case 'argmax = B1'
            [B_hom_het(:,k),lambda_hom_het(k)] = beta_optCM_hom(Adj_k,param,knob.optB,b1,b2,B1(1,1),blim);
            [B_het(:,k),lambda_het(k)] = beta_optCM_het(Adj_k,param,knob.optB,b1,b2,B_hom_het(:,k),blim,Nopt);
            B1 = diag(B_het(:,k));
            B2 = b2*eye(N);
        case 'argmax = B2'
            [B_hom_het(:,k),lambda_hom_het(k)] = beta_optCM_hom(Adj_k,param,knob.optB,b1,b2,B2(1,1),blim);
            [B_het(:,k),lambda_het(k)] = beta_optCM_het(Adj_k,param,knob.optB,b1,b2,B_hom_het(:,k),blim,Nopt);
            B1 = b1*eye(N);
            B2 = diag(B_het(:,k));
        case 'argmax = B'
            [B_hom_het(:,k),lambda_hom_het(k)] = beta_optCM_hom(Adj_k,param,knob.optB,b1,b2,B1(1,1),blim);
            [B_het(:,k),lambda_het(k)] = beta_optCM_het(Adj_k,param,knob.optB,b1,b2,B_hom_het(:,k),blim,Nopt);
            B1 = diag(B_het(:,k));
            B2 = diag(B_het(:,k));
        case 'argmax = B1, B2'
            [B_hom_het(:,k),lambda_hom_het(k)] = beta_optCM_hom(Adj_k,param,knob.optB,b1,b2,[B1(1,1); B2(1,1)],blim);
            [B_het(:,k),lambda_het(k)] = beta_optCM_het(Adj_k,param,knob.optB,b1,b2,B_hom_het(:,k),blim,Nopt);
            B1 = diag(B_het(1:N,k));
            B2 = diag(B_het(N+1:2*N,k));
    end

    % ODE integration
    switch knob.adjupdate
        case 'continuous'
            [ttemp,xtemp] = odeRK(@(t,x)CMflock_continuous(t,x,N,dim,beta,sigma,K,B1,B2,gamma,pg,mass,noise,knob.ref_traj),[t_het(end) dt t_het(end)+tstep],x_het(:,end)');
        case 'piecewise'
            [ttemp,xtemp] = odeRK(@(t,x)CMflock_piecewise(t,x,N,dim,B1,B2,gamma,pg,mass,noise,Lap_k,knob.ref_traj),[t_het(end) dt t_het(end)+tstep],x_het(:,end)');
    end
    xtemp = xtemp';
    t_het = [t_het; ttemp(2:end)'];
    x_het = [x_het  xtemp(:,2:end)];
    xtemp = []; ttemp = [];
end


%% Performance analysis

% Homogeneous optimization
[Xhom,RMSE_hom] = agentcoord_CM(x_hom,dim,N,pg);
figure(1)
subplot(231); plot(t_hom,Xhom.px - Xhom.pref_x - pg_x)
subplot(232); plot(t_hom,Xhom.qx - Xhom.qref_x)
subplot(233); plot(t_hom,RMSE_hom )

% Heterogeneous optimization
[Xhet,RMSE_het] = agentcoord_CM(x_het,dim,N,pg);
subplot(234); plot(t_het,Xhet.px - Xhet.pref_x - pg_x)
subplot(235); plot(t_het,Xhet.qx - Xhet.qref_x)
subplot(236); plot(t_het,RMSE_het )

% RMSE comparison
figure(2)
subplot(131); semilogy(t_hom,RMSE_hom,t_het,RMSE_het)
xlabel('time'); ylabel('RMSE')
legend('homogeneous','heterogeneous')
title('Convergence time')

subplot(132); plot(1:Nsteps,lambda_hom,1:Nsteps,lambda_het)
xlabel('time'); ylabel('\lambda_m_a_x')
legend('homogeneous','heterogeneous')
title('Stability hom vs het')

subplot(133); plot(1:Nsteps,lambda_hom_het,1:Nsteps,lambda_het)
xlabel('time'); ylabel('\lambda_m_a_x')
legend('homogeneous','heterogeneous')
title('Stability relative improvement')

