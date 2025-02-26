%% Distributed optimization of heterogeneous agents
%  using the Ren (centroid) model of flocking.
clear all; clc; close all;

% Eigenvalue optimization functions
addpath([pwd,'/EigJacobian/'])
addpath([pwd,'/EigOptimization/'])

% Differential equations and data processing
addpath([pwd,'/FlockODEs/'])

% Parallel computation
poolobj = gcp('nocreate');
delete(poolobj);
ncores = input('Number of cores: ')
parpool(ncores);

%% Setup simulation

% Parameter optimization
knob.optB = 'argmax = B1, B2'

% Desired formation
knob.formation = 'random';
    % random; lattice

% Update rule for the adjacency matrix
knob.adjupdate = 'piecewise';
    % continuous; piecewise

knob.ref_traj = 'constant';
    % stationary; constant; zigzag; circular; fuzzy

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

% Initial conditions
p0 = unifrnd(-2,2,N*dim,1);        % initial position
q0 = unifrnd(0,0,N*dim,1);         % initial velocity
switch knob.ref_traj               % initial reference pos/velocity
    case 'stationary'
        pref0 = 100*[1; 1];
        qref0 = 0*[1; 1];
    case 'constant'
        pref0 = 100*[1; 1];
        qref0 = 100*[1; 0];
end
x0 = [p0;q0;pref0;qref0];          % state vector

% Relative positions
switch knob.formation
    case 'random'
    % uniform random distribution
        pg_x = unifrnd(-5,5,N,1);
        pg_y = unifrnd(-5,5,N,1);
    case 'lattice'
        
end
if dim == 2
    pg = [pg_x; pg_y];
else
    pg = unifrnd(-5,5,dim*N,1);
end
% figure(1); scatter(pg_x,pg_y)

%% Damping coefficients
b1 = 10;              % position-integrator feedback gain
b2 = 10;              % velocity-integrator feedback gain
blim = 30;          % upper bound for feedback gain

% Simulation time
dt = 1e-3;
tf = 30;

% Real-time optimization
tstep = tf;
Nsteps = tf/tstep;
Nopt = 10;                      % MC iterations for optimization

%% Eigenvalue optimization
Nmc = 100;             % number of realizations
Nradius = [0.1 0.25 0.5:0.25:10];
radius_length = length(Nradius);

% Initialization
lambda_hom_global = zeros(Nmc,1);
lambda_het_global = zeros(Nmc,1);
lambda_het = zeros(Nmc,radius_length);
N_neighbors_mean = zeros(Nmc,radius_length);
N_neighbors_std = zeros(Nmc,radius_length);

% Parallel computation
parfor i = 1:Nmc

    % Initial conditions
    p0 = unifrnd(-2,2,N*dim,1);        % initial position
    q0 = unifrnd(0,0,N*dim,1);         % initial velocity
    switch knob.ref_traj               % initial reference pos/velocity
        case 'stationary'
            pref0 = 100*[1; 1];
            qref0 = 0*[1; 1];
        case 'constant'
            pref0 = 100*[1; 1];
            qref0 = 100*[1; 0];
    end
    x0 = [p0;q0;pref0;qref0];          % state vector

    % Relative positions
    switch knob.formation
        case 'random'
            % uniform random distribution
            pg_x = unifrnd(-5,5,N,1);
            pg_y = unifrnd(-5,5,N,1);
        case 'lattice'

    end
    pg = [pg_x; pg_y];
    
    % State of the system at final formation
    x_bfgs = [pg; kron(qref0,ones(30,1)); pref0; qref0];

    % Maximum radius for optimization
    for j = 1:radius_length
        radius = Nradius(j);

        % Feedback coefficients
        B1_het = zeros(N,N);
        B2_het = zeros(N,N);

        % Adjacency matrix computed for each agent
        p_k = x_bfgs(1:2*N,end);
        [Lap_k,Adj_k] = flocklaplacian(p_k,sigma,beta,K,N);
        N_neighbors = zeros(N,1);
        for ii = 1:N
            pii = p_k([ii ii+N],1);
            Graph_bin = zeros(N,N);

            % Neighborhood radius
            neighborhood_ii = [];
            for jj = 1:N
                pjj = p_k([jj jj+N],1);
                if vecnorm(pii - pjj) < radius
                    neighborhood_ii = [neighborhood_ii jj];
                    if ii == jj
                        neighborhood_index = length(neighborhood_ii);
                    end
                end
            end

            % Local adjacency matrix
            Adj_ii = Adj_k(neighborhood_ii,neighborhood_ii);
            N_neighbors(ii) = size(Adj_ii,2);
            param_ii = param; param_ii.N = N_neighbors(ii);

            % Eigenvalue optimization
            [B_hom_bfgs,lambda_hom_bfgs] = beta_optCM_hom(Adj_ii,param_ii,knob.optB,10,10,[10; 10],blim);
            [B_bfgs,lambda_bfgs] = beta_optCM_bfgs(Adj_ii,param_ii,knob.optB,10,10,B_hom_bfgs,blim,Nopt);
            B1_het(ii,ii) = B_bfgs(neighborhood_index);
            B2_het(ii,ii) = B_bfgs(neighborhood_index+N_neighbors(ii));
        end

        % Distributed optimization
        b_het = [diag(B1_het); diag(B2_het)];
        lambda_het(i,j) = opteigreal_CMhet(b_het,Lap_k,gamma,N,dim,mass,knob.optB,[],[]);
        N_neighbors_mean(i,j) = mean(N_neighbors);
        N_neighbors_std(i,j) = std(N_neighbors);
    end

    % Global optimization
    [B_hom_global,lambda_hom_global(i)] = beta_optCM_hom(Adj_k,param,knob.optB,10,10,[10; 10],blim);
    [B_het_global,lambda_het_global(i)] = beta_optCM_bfgs(Adj_k,param,knob.optB,b1,b2,B_hom_global,blim,3*Nopt);
end

%% Plot
figure(1);
subplot(121); 
plot(Nradius,median(lambda_hom_global,1)*ones(radius_length,1)); hold on;
plot(Nradius,median(lambda_het_global,1)*ones(radius_length,1));
plot(Nradius,median(lambda_het,1));
legend('homogeneous, centralized', 'heterogeneous, centralized', 'heterogeneous, distributed')
xlabel('sensing radius R')
ylabel('Lyapunov exponent \lambda_m_a_x')

subplot(122); plot(Nradius,median(N_neighbors_mean,1));
xlabel('sensing radius R')
ylabel('neighborhood size')

% save data_distributedeig.mat
