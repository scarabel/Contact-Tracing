% Copyright (c) 2021 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE for 
% licensing information. 
% 
% If you use this code, please cite 
% Scarabel, Pellis, Ogden, Wu, 'A renewal equation model to assess roles and
% limitations of contact tracing for disease outbreak control',
% Royal Society Open Science, 2021.
% 
%% Script R0_ec.m
% Computes the contour lines of R_{d,c}=1 in the plane (R0,epsilon_c).
% Uses the function linear_contact_tracing.m

clear
close all

step = 0.05; % stepsize for numerical solution

% Epidemiological parameters

% Distribution of incubation time: Gamma distribution (Overton et al, 2020)
mean_incubation = 4.84;
std_incubation = 2.79;

shape_incubation = (mean_incubation/std_incubation)^2;
scale_incubation = std_incubation^2/mean_incubation;

density_incubation_f = @(s) gampdf(s,shape_incubation,scale_incubation);
surv_symptoms_f = @(s) 1-integral(@(y) density_incubation_f(y),0,s);

% infectiousness profile: Gamma distribution (Ferretti et al, 2020)
bmax = 20; % maximal bound to infectiousness period

mean_beta = 5;
std_beta = 1.9;

shape_beta = (mean_beta/std_beta)^2;
scale_beta = std_beta^2/mean_beta;

% percentage symptomatic from He et al, 2020, Systematic review: 85%
epsilon_s = 0.85;

% diagnosis parameters
dmax = 20;
delay_diagnosis = 2;
epsilon_d = 1; % fraction of symptomatic individuals diagnosed

density_diagnosis = @(x) epsilon_d*epsilon_s*(x<=dmax).*gampdf(x-delay_diagnosis,shape_incubation,scale_incubation);
surv_diagnosis_f = @(x) 1-integral(@(y) density_diagnosis(y),0,x);

% Tracing window
cmax = 5;

%% Iterations to explore different combinations of parameters

legendname = {};
epsilon_c_vector = 0:0.1:1;
R0_vector = 0.2:0.2:3;

for index_x = 1:length(R0_vector)
for index_y = 1:length(epsilon_c_vector)

    epsilon_c = epsilon_c_vector(index_y);
    R0 = R0_vector(index_x);

    beta_transm = @(x) R0*(x<=bmax).*gampdf(x,shape_beta,scale_beta);

    % Solution of the linear system
    nd = dmax/step;
    nc = cmax/step;
    nb = bmax/step;

    N = max([nb,nd,nc,nc+nb]);
    dgrid = step*(1:nd);
    bgrid = step*(1:nb);
    Ngrid = step*(1:N);

    % Initialization of known parameters (discretization of functions)
    beta_mat = zeros(N,1);
    h_d = zeros(N,1);
    surv_d = (1-epsilon_d*epsilon_s)*ones(N,1); % survival diagnosis
    dens_d = zeros(N,1);

    for itau = 1:N
        tau = itau*step;
        beta_mat(itau) = beta_transm(tau);
    end

    dens_d(1) = density_diagnosis(step);
    surv_d(1) = surv_diagnosis_f(step);
    h_d(1) = -log(surv_d(1))/step;
    for itau = 2:(dmax/step)
        dens_d(itau) = density_diagnosis(step*itau);
        surv_d(itau) = surv_diagnosis_f(itau*step);
        h_d(itau) = - (log(surv_d(itau))-log(surv_d(itau-1)))/step;
    end
    surv_d(nd+1:end)=surv_d(nd);

    % calculation of reproduction number via quadrature formulas
    R0 = step*trapz(beta_mat);
    Rd = step*trapz(beta_mat.*surv_d);
    
    r0 = fsolve(@(x) 1- step*trapz(beta_mat.*exp(-x*step*(1:N)')), 0.1);
    rd = fsolve(@(x) 1- step*trapz(beta_mat.*surv_d.*exp(-x*step*(1:N)')), r0);
    
    % initialize probability of contact tracing
    x0 = zeros(N+1,1); % the last entry will represent the exponential growth rate
    x0(1:nc)=zeros(1,nc);
    x0(end)=rd;
    % linear_contact_tracing(x0(1:N),x0(N+1),step,nc,nd,epsilon_c,beta_mat,h_d,surv_d)

    options = optimoptions('fsolve','Display','none','MaxIter',100000);
    Sol = fsolve(@(x) [x(1:N);1] - linear_contact_tracing(x(1:N),x(N+1),step,nc,nd,epsilon_c,beta_mat,h_d,surv_d), x0);
    h_ct = Sol(1:N);
    rct = Sol(N+1);

    % calculation of Rct, reproduction number with contact tracing
    Rct = 0;
    surv_ct = zeros(N,1);
    for itau = 1:N
        surv_ct(itau) = exp(-step*sum(h_ct(1:itau)));
        Rct = Rct + step*beta_mat(itau)*surv_d(itau)*surv_ct(itau);
    end
 
    Rct_matrix(index_y,index_x) = Rct;
    rct_matrix(index_y,index_x) = rct;
    
end

end

%% Contour plot

[XX,YY]=meshgrid(R0_vector,epsilon_c_vector);

figure
contour(XX,YY,Rct_matrix,0.2:0.2:3,'ShowText','on','LineWidth',2); hold on
xlabel('$R_0$','Interpreter','latex');
ylabel('fraction of contacts traced','Interpreter','latex');
title('Effective reproduction number $R_{d,c}$','Interpreter','latex');
set(gca,'fontsize',14)
