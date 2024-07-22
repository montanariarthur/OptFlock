%% Real-time optimization of the Olfati-Saber flocking model (OS)
% This code compares the performance of a flock of agents optimized in real
% time (over regular optimization windows) with:
%  - constrained parameters (homogeneous flock)
%  - unconstrained parameters (heterogeneous flock)

clear all; close all; clc;

% Eigenvalue optimization functions
addpath([pwd,'/EigJacobian/'])
addpath([pwd,'/EigOptimization/'])

% Differential equations and data processing
addpath([pwd,'/FlockODEs/'])
addpath([pwd,'/FlockODEs/OlfatiSaberModel'])

% Plots
addpath([pwd,'/ShadedPlots/'])

%% Setup simulation

% Target trajectory
knob.ref_traj = 'constant';
    % 'stationary' - does not move
    % 'constant'   - moves with constant speed along x-axis

% Protocol ('Algorithm 2' or 'Algorithm 3' according to Ref. 2)
knob.protocol = 'Algorithm 3';

% Distribution of obstacles in space
knob.obstacles = 'sobolseq';          
    % 'random'    - uniform distribution
    % 'specified' - specified manually
    % 'sobolseq'  - Sobol sequence

disp(knob)

%% Parameters

% Network interaction parameters
param.epsilon = 0.1;            % for sigma norm
param.h_alpha = 0.2;            % for bump function of alpha-agents
param.h_beta = 0.9;             % for bump function of beta-agents
param.a = 5; param.b = 5;       % for phi function
param.dist = 7;                 % lattice constrained distance
param.range = 1.2*param.dist;   % interaction range
param.dist_beta = 0.6*param.dist;           % lattice constrained distance 
param.range_beta = 1.2*param.dist_beta;     % to beta agent

% System parameters
N = 30;        param.N = N;         % number of agents
dim = 2;       param.dim = dim;     % dimension of the Euclidean space

%% Gain coefficients 

% (alpha-gamma)
b1 = 1;   %1           % position-integrator feedback gain 
b2 = 1;   %2*sqrt(b1);     % velocity-integrator feedback gain 
blim = [5; 5];            % upper bound for feedback gain

% (alpha-alpha)
C1alpha = 30;
C2alpha = 2*sqrt(C1alpha);

% (alpha-beta)
C1beta = 300;
C2beta = 2*sqrt(C1beta);

param.C1alpha = C1alpha;
param.C2alpha = C2alpha;
param.C1beta = C1beta;
param.C2beta = C2beta;

%% Initial conditions
q0 = unifrnd(-50,50,N*dim,1);    % initial position
p0 = unifrnd(0,0,N*dim,1);         % initial velocity

switch knob.protocol
    case {'Algorithm 1'}
        x0 = [q0;p0];
    case {'Algorithm 2','Algorithm 3'}
        switch knob.ref_traj
            case 'stationary'
                qref0 = 1000*[1; 0];
                pref0 = 0*[1; 1];
            case 'constant'
                qref0 = [0; 0];
                pref0 = [20; 0];
        end
        x0 = [q0;p0;qref0;pref0];
end

% Simulation time
dt = 1e-2;
tf = 30;

% Real-time optimization
tstep = 0.5;
Nsteps = tf/tstep;
Nopt = 10;                      % MC iterations for optimization

%% Define obstacles
switch knob.protocol
    case 'Algorithm 3'
        switch knob.obstacles
            case 'random'
                param.Nbeta = 30;
                switch dim
                    case 2
                        qbeta_x = unifrnd(150,600,param.Nbeta,1);
                        qbeta_y = unifrnd(-20,20,param.Nbeta,1);
                        qbeta = [qbeta_x; qbeta_y];
                end
                pbeta = unifrnd(0,0,param.Nbeta*dim,1);
                %x0 = [x0; qbeta; pbeta];
                param.radius_beta = unifrnd(1,4,param.Nbeta,1);
            case 'specified'
                param.Nbeta = 0;
                qbeta_x = [100 110 120 130 150 160]';
                qbeta_y = [20 60 40 -20 40 0]';
                qbeta = [qbeta_x; qbeta_y];
                pbeta = unifrnd(0,0,param.Nbeta*dim,1);
                param.radius_beta = [10 4 2 5 5 3]';
            case 'sobolseq'
                param.Nbeta = 30/2;
                % Define the bounds of the box [xmin xmax; ymin ymax; ...]
                bounds = [150 600; -20 20]; % example bounds for a 2D box
                % Create a Sobol sequence
                p = sobolset(dim);
                % Generate the points
                points = net(p, param.Nbeta);
                % Scale the points to fit within the specified bounds
                for i = 1:dim
                    qbeta(param.Nbeta*(i-1)+1:param.Nbeta*i,1) = bounds(i, 1) + points(:, i) * (bounds(i, 2) - bounds(i, 1));
                end
                pbeta = unifrnd(0,0,param.Nbeta*dim,1);
                param.radius_beta = unifrnd(1,4,param.Nbeta,1);
        end
    case {'Algorithm 1','Algorithm 2'}
        param.Nbeta = 0;
        qbeta = [];
end

%% Online simulation (homogeneous optimization)
t_hom = 0; x_hom = x0;
B1gamma = b1*eye(N); B2gamma = b2*eye(N);
for k = 1:Nsteps

    % Network structure: alpha-alpha
    q = x_hom(1:N*dim,end);
    param_hom = param; param_hom.h = param_hom.h_alpha;
    for i = 1:N
        qi = q(i:N:i+N*(dim-1));
        for j = 1:N
            qj = q(j:N:j+N*(dim-1));
            distance = norm(qi - qj);
            if distance < param.range && i ~= j
                Adj_k(i,j) = spatial_adj_entry(qi,qj,param_hom);
            else
                Adj_k(i,j) = 0;
            end
        end
    end
    Lap_alpha_k = Adj_k - diag(sum(Adj_k,2));

    % Network structure: alpha-beta
    Lap_beta_k = zeros(N*dim,N*dim);
    for i = 1:N
        qi = q(i:N:i+N*(dim-1));
        Adj_beta_k = zeros(dim,dim);
        for kk = 1:param.Nbeta
            % obstacle center and radius
            yk = qbeta(kk:param.Nbeta:kk+param.Nbeta*(dim-1),1);
            Rk = param.radius_beta(kk);

            % calculates position and velocity according to Sec. VII-D
            ak = (qi - yk)/norm(qi - yk);
            mu = Rk / norm(qi - yk);
            P = eye(dim) - ak*ak';
            qk = mu*qi + (1 - mu)*yk;

            % if it is within the interaction range
            distance_beta = norm(qi - qk);
            if distance_beta < param.range_beta
                pk = mu*P*pi;
                Adj_beta_k = Adj_beta_k + spatial_adj_beta(qi,qk,param_hom)*(mu*P - eye(dim));
            end
        end
        Lap_beta_k(i:N:N*dim,i:N:N*dim) = Adj_beta_k;
    end


    % Eigenvalue optimization
    [B_hom(:,k),lambda_hom(k)] = beta_optOS_hom(Lap_alpha_k,Lap_beta_k,param,[B1gamma(1,1); B2gamma(1,1)],blim);
    B1gamma = B_hom(1,k)*eye(N);
    B2gamma = B_hom(2,k)*eye(N);

    % ODE integration
    [ttemp,xtemp] = odeRK(@(t,x)OSflock(t,x,C1alpha,C2alpha,B1gamma,B2gamma,C1beta,C2beta,qbeta,param,knob),[t_hom(end) dt t_hom(end)+tstep],x_hom(:,end)');
    xtemp = xtemp';
    t_hom = [t_hom; ttemp(2:end)'];
    x_hom = [x_hom  xtemp(:,2:end)];
    xtemp = []; ttemp = [];
end


%% Online simulation (heterogeneous optimization)
t_het = 0; x_het = x0;
for k = 1:Nsteps
    % Network structure: alpha-alpha
    q = x_het(1:N*dim,end);
    param_het = param; param_het.h = param_het.h_alpha;
    for i = 1:N
        qi = q(i:N:i+N*(dim-1));
        for j = 1:N
            qj = q(j:N:j+N*(dim-1));
            distance = norm(qi - qj);
            if distance < param.range && i ~= j
                Adj_k(i,j) = spatial_adj_entry(qi,qj,param_het);
            else
                Adj_k(i,j) = 0;
            end
        end
    end
    Lap_alpha_k = Adj_k - diag(sum(Adj_k,2));

    % Network structure: alpha-beta
    Lap_beta_k = zeros(N*dim,N*dim);
    for i = 1:N
        qi = q(i:N:i+N*(dim-1));
        Adj_beta_k = zeros(dim,dim);
        for kk = 1:param.Nbeta
            % obstacle center and radius
            yk = qbeta(kk:param.Nbeta:kk+param.Nbeta*(dim-1),1);
            Rk = param.radius_beta(kk);

            % calculates position and velocity according to Sec. VII-D
            ak = (qi - yk)/norm(qi - yk);
            mu = Rk / norm(qi - yk);
            P = eye(dim) - ak*ak';
            qk = mu*qi + (1 - mu)*yk;

            % if it is within the interaction range
            distance_beta = norm(qi - qk);
            if distance_beta < param.range_beta
                pk = mu*P*pi;
                Adj_beta_k = Adj_beta_k + spatial_adj_beta(qi,qk,param_het)*(mu*P - eye(dim));
            end
        end
        Lap_beta_k(i:N:N*dim,i:N:N*dim) = Adj_beta_k;
    end

    % Eigenvalue optimization
    [B_hom_bfgs(:,k),lambda_hom_bfgs(k)] = beta_optOS_hom(Lap_alpha_k,Lap_beta_k,param,[B1gamma(1,1); B2gamma(1,1)],blim);
    [B_het(:,k),lambda_bfgs(k)] = beta_optOS_het(Lap_alpha_k,Lap_beta_k,param,B_hom_bfgs(:,k),blim,Nopt);
    B1gamma = diag(B_het(1:N,k));
    B2gamma = diag(B_het(N+1:2*N,k));

    % ODE integration
    [ttemp,xtemp] = odeRK(@(t,x)OSflock(t,x,C1alpha,C2alpha,B1gamma,B2gamma,C1beta,C2beta,qbeta,param,knob),[t_het(end) dt t_het(end)+tstep],x_het(:,end)');
    xtemp = xtemp';
    t_het = [t_het; ttemp(2:end)'];
    x_het = [x_het  xtemp(:,2:end)];
    xtemp = []; ttemp = [];
end

%% Time-series data
[Xhom,Ehom,Chom,Vhom,Thom] = agentcoord_OS(x_hom,dim,N,param,knob);
[Xhet,Ehet,Chet,Vhet,Thet] = agentcoord_OS(x_het,dim,N,param,knob);

%% Time-series plot
figure(1);
tinst = tf/dt;
subplot(231); plot(t_hom,Xhom.qx); title('Agent position')
subplot(232); plot(t_hom,Xhom.px); title('Agent velocity')
subplot(233); scatter(Xhom.qx(:,tinst),Xhom.qy(:,tinst)); title('Final formation')
subplot(234); plot(t_het,Xhet.qx);
subplot(235); plot(t_het,Xhet.px);
subplot(236); scatter(Xhet.qx(:,tinst),Xhet.qy(:,tinst));

% Damping
figure(2)
subplot(221); stairs(B_hom(1,:)); title('Position feedback gain')
subplot(222); stairs(B_hom(2,:)); title('Velocity feedback gain')
subplot(223); stairs(B_het(1:N,:));
subplot(224); stairs(B_het(N+1:2*N,:)); 

% Performance plot
figure(3)
subplot(141); plot(t_hom,Ehom,t_het,Ehet)
legend('hom','het')
subplot(142); plot(t_hom,Chom,t_het,Chet)
subplot(143); plot(t_hom,Vhom,t_het,Vhet)
% set(gca,'YScale','log')
subplot(144); plot(t_hom,Thom,t_het,Thet)
% set(gca,'YScale','log')

