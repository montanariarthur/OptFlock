function stability=stst_stabil_mod(J1,J2,tau,method) 
% This is a modified version of original stst_stabil function
% that directly returns the arguments J1,J2,tau instead of the
% variable stst (steady state point).
% function stability=stst_stabil(stst,method)
% INPUT:
%   J1, J2                    -   Jacobians in the generalized eigenvalue problem
%   tau                       -   Time delay
%	method                    -   Method parameters 
% OUTPUT:
%	stability                 -   Stability information
% COMMENT:
%       Assumes (imag(method.lms_parameter_rho)~=0)
%       This condition is tested in the code.

% Original code:  
% (c) DDE-BIFTOOL v. 2.03, 05/03/2007
% Added on 05/03/2007 
% Modified by Ana Barioni, 07/30/2024

if imag(method.stability.lms_parameter_rho)==0,
  error(['STST_STABIL: imag(method.lms_parameter_rho)==0 :' ...
	 ' hence, method contains inapproriate paramters for this function. ' ...
	 'Use method=df_mthod(''stst'',1);stst_stabil(stst,' ...
	 ' method.stability), to obtain an appropriate method.stability' ...
	 ' structure for this function.']);
end;

nb_nu = 20;
delta_real_min = 0.1;

% MODIFICATION: Does not need if statement or for loop since our system only has one time delay, given in the arguments

taumin=min(tau);
taumax=max(tau);

% MODIFICATION: Takes the Jacobians from the argument instead of calculating the derivatives of the dynamical equations
AA{1} = J1;
AA{2} = J2;


n=size(AA{1},2);

alpha=method.stability.lms_parameter_alpha;
beta=method.stability.lms_parameter_beta;
k_lms = length(alpha) - 1;

interp_order=method.stability.interpolation_order;
order_method=interp_order;
s_min = floor((interp_order-1)/2);
s_plus = (interp_order-1) - s_min;
real_min=method.stability.minimal_real_part;
if isempty(real_min) | (real_min==-Inf) 
  if taumax>0 & taumin>0.1
    real_min=-1/taumin;
  else
    real_min=-1;
  end
end

a_ellipse=real(method.stability.lms_parameter_rho);
b_ellipse=imag(method.stability.lms_parameter_rho);
  
% 
eigA0 = eig(AA{1}).';
pts=get_pts_h_new(AA,tau,real_min,delta_real_min,nb_nu,eigA0);
hh=0.9/sqrt(max((real(pts)/a_ellipse).^2+(imag(pts)/b_ellipse).^2));

hh=min([hh,min(tau(tau>=hh*interp_order*1e-2))/s_plus,taumax*method.stability.maximal_time_step]);
hh=max([hh,taumax*method.stability.minimal_time_step]);
nn=n*(k_lms + ceil(taumax/hh) + s_min);

[mu,nL]=help_stst_stabil(AA,tau,hh,alpha,beta,interp_order);   

% Note: nn==nL,
% except if some interpolation can be avoided ...

% Throw away mu too close to the origin
ss=exp(real_min*hh);
mu=mu(ss<=abs(mu));

% % Throw away zeros
% mu=find(abs(mu));
% % Zeros were already discarded in the previous step

lambda=mu_to_lambda(mu,hh);

% "A posteriori safeguard", 
% but here we approximate (0.9/h)*LMS^{-1}(T_delta) 
% by (0.9/h)*[ellipse(a_ellipse,b_ellipse)] ...:
lambda=lambda((real(lambda)/a_ellipse).^2+(imag(lambda)/b_ellipse).^2<=(0.9/hh)^2);

[dummy,idx]=sort(real(lambda));
lambda=lambda(idx(end:-1:1));

if length(lambda)>method.stability.max_number_of_eigenvalues
  lambda=lambda(1:method.stability.max_number_of_eigenvalues);
end;

stability=struct('h',hh,'l0',lambda,'l1',[],'n1',[]);

stability=stst_stabil_nwt_corr(stability,AA,tau,method);

return;
