%% Stability optimization of the time-delay (TD) consensus model.
% This code computes the optimal choice of homogeneous and heterogeneous
% parameter 'b' for a given time-delay 'tau'

% Declare parameters
n = 4; % Number of agents
Nopt = 100; % Number of optimization trials
tau = 0.6; % Time delay
L =  [1 0 -1 0 ;-1 1 0 0; 0 -1 1 0; -1 0 0 1]; % Laplacian matrix

% Find best homogeneous solution
b0_hom = 1;
blim = 1;
B_bfgs_mc_hom = zeros(1,Nopt);
lambda_bfgs_mc_hom = zeros(Nopt,1);
for MCopt = 1:Nopt
    B_bfgs_mc_hom(1,MCopt) = fmincon(@(b)opteigreal_CSDhom(b,L,tau,n),[b0_hom(1)+1*randn(1,1)],[],[],[],[],0.1,blim);
    lambda_bfgs_mc_hom(MCopt) = opteigreal_CSDhom(B_bfgs_mc_hom(:,MCopt),L,tau,n);
end

% Picks the best homogeneous solution
[~,minindex_hom] = min(lambda_bfgs_mc_hom);
lyapexp_hom = lambda_bfgs_mc_hom(minindex_hom);
b_opt_hom = B_bfgs_mc_hom(:,minindex_hom);


% Find best heterogeneous solution
b0_het = b_opt_hom;
B_bfgs_mc_het = zeros(n,Nopt);
lambda_bfgs_mc_het = zeros(Nopt,1);
for MCopt = 1:Nopt
    B_bfgs_mc_het(:,MCopt) = fmincon(@(b)opteigreal_CSDhet(b,L,tau,n),[b0_het(1)+0.1*randn(n,1)],[],[],[],[],0.1*ones(n,1),blim*ones(n,1));
    lambda_bfgs_mc_het(MCopt) = opteigreal_CSDhet(B_bfgs_mc_het(:,MCopt),L,tau,n);
end

% Picks the best heterogeneous solution
[~,minindex_het] = min(lambda_bfgs_mc_het);
lyapexp_het = lambda_bfgs_mc_het(minindex_het);
b_opt_het = B_bfgs_mc_het(:,minindex_het);

%% TimeDelay_dynamics
% Calculates the time series of the dynamics of the time delay consensus 
% model for the optimal choices in the homogeneous and heterogeneous

beta_hom = b_opt_hom*ones(n,1);  % Homogeneous optimal 
beta_het = b_opt_het; % Heterogeneous optimal 

%%%%--- Parameters vector ---%%%%%
param.n = n;       
param.L = L;              
param.tau = tau;  

% simulating the system for T time units
T = 80;
time = [0 T];
opts = ddeset('RelTol',1e-3,'AbsTol',1e-6, 'MaxStep', 1e-2);

x0 = unifrnd(-2,2,n,1);
y0 = unifrnd(-2,2,n,1);
vx0 = unifrnd(0,2,n,1);
vy0 = unifrnd(0,2,n,1);
state0 = [x0;y0;vx0;vy0];

% Solving for heterogeneous beta
rhs_het = @(t,x,Z)rhs_aux(t,x,Z,n,beta_het,L);
rhshist = @(t)rhshist_aux(t,x0,y0,vx0,vy0);

% solving the DDE
sol_het = dde23(rhs_het,tau,rhshist,time,opts);

posx_het = sol_het.y(1:n,:);
posy_het = sol_het.y(n+1:2*n,:);
velx_het = sol_het.y(2*n+1:3*n,:);
vely_het = sol_het.y(3*n+1:4*n,:);

% Solving for homogeneous beta
rhs_hom = @(t,x,Z)rhs_aux(t,x,Z,n,beta_hom,L);

% solving the DDE
sol_hom = dde23(rhs_hom,tau,rhshist,time,opts);

posx_hom = sol_hom.y(1:n,:);
posy_hom = sol_hom.y(n+1:2*n,:);
velx_hom = sol_hom.y(2*n+1:3*n,:);
vely_hom = sol_hom.y(3*n+1:4*n,:);

%% Ploting the figures

% Homogeneous position
figure();
plot(sol_hom.x,posx_hom(1,:))
hold on
for iim= 2:n
plot(sol_hom.x,posx_hom(iim,:))
end

figure();
plot(sol_hom.x,posy_hom(1,:))
hold on
for iim= 2:n
plot(sol_hom.x,posy_hom(iim,:))
end

% Homogeneous velocity
figure();
plot(sol_hom.x,velx_hom(1,:))
hold on
for iim= 2:n
plot(sol_hom.x,velx_hom(iim,:))
end

figure();
plot(sol_hom.x,vely_hom(1,:))
hold on
for iim= 2:n
plot(sol_hom.x,vely_hom(iim,:))
end

% Heterogeneous position
figure();
plot(sol_het.x,posx_het(1,:))
hold on
for iim= 2:n
plot(sol_het.x,posx_het(iim,:))
end

figure();
plot(sol_het.x,posy_het(1,:))
hold on
for iim= 2:n
plot(sol_het.x,posy_het(iim,:))
end

% Heterogeneous velocity
figure();
plot(sol_het.x,velx_het(1,:))
hold on
for iim= 2:n
plot(sol_het.x,velx_het(iim,:))
end

figure();
plot(sol_het.x,vely_het(1,:))
hold on
for iim= 2:n
plot(sol_het.x,vely_het(iim,:))
end


%%%%%--- Defining the DDEs ---%%%%%
function dxdt = rhs_aux(t,state,Z,n,beta_v,L) 

    dxdt = [state(2*n+1:4*n);-diag(beta_v)*L*(Z(1:n) + Z(2*n+1:3*n));-diag(beta_v)*L*(Z(n+1:2*n) + Z(3*n+1:4*n))];
end

%%%%%--- Defining the initial conditions ---%%%%%
function S = rhshist_aux(t,x0,y0,vx0,vy0)	
    %%%%--- Parameters ---%%%%%
	S = [x0;y0;vx0;vy0];
end


