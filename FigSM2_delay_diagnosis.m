% Copyright (c) 2021 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE for 
% licensing information. 
% 
% If you use this code, please cite 
% Scarabel, Pellis, Ogden, Wu, 'A renewal equation model to assess roles and
% limitations of contact tracing for disease outbreak control',
% Royal Society Open Science, 2021.
% 
%% script delay_diagnosis.m
% Computes R_{d,c} and related quantities varying epsilon_c and the diagnosis delay.
% Uses the function linear_contact_tracing.m
% Figure SM2 in the Supplementary Material is obtained with
% R0 = 2.5 and R0 = 1.5, respectively

% Iterations to explore different combinations of parameters

clearvars

step = 0.05; % stepsize for numerical solution

% Epidemiological parameters

% Basic reproduction number
R0 = 2.5; % 1.5; % 2.5

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

beta_transm = @(x) R0*(x<=bmax).*gampdf(x,shape_beta,scale_beta);

% percentage symptomatic from He et al, 2020, Systematic review: 85%
epsilon_s = 0.85; % percentage of symptomatic individuals

% diagnosis parameters
dmax = 20;
epsilon_d = 1; % fraction of symptomatic individuals diagnosed

cmax = 5;
legendname = {};
epsilon_c_vector = 0:0.1:1;
delay_vector = 0:2;

figure(20); clf; hold on
plot([0 1],[1 1],'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',1,'HandleVisibility','off'); % reference line R=1

colorcode = lines(length(delay_vector));
opt_fsolve = optimoptions('fsolve','Display','none','MaxIter',100000);

Lx = length(delay_vector);
Ly = length(epsilon_c_vector);

Rct_matrix = zeros(Ly,Lx);
rct_matrix = Rct_matrix;

Theta_d = Rct_matrix;
Theta_ct = Rct_matrix;

for index_x = 1:Lx
for index_y = 1:Ly

    delay_diagnosis = delay_vector(index_x);
    epsilon_c = epsilon_c_vector(index_y);
    
    density_diagnosis = @(x) epsilon_d*epsilon_s*(x<=dmax).*gampdf(x-delay_diagnosis,shape_incubation,scale_incubation);
    surv_diagnosis_f = @(x) 1-integral(@(y) density_diagnosis(y),0,x);

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
    surv_d = (1-epsilon_d*epsilon_s)*ones(N,1); % probability of not being diagnosed
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
    
    r0 = fsolve(@(x) 1- step*trapz(beta_mat.*exp(-x*step*(1:N)')), 0.1, opt_fsolve);
    rd = fsolve(@(x) 1- step*trapz(beta_mat.*surv_d.*exp(-x*step*(1:N)')), r0, opt_fsolve);
    
    % initialize probability of contact tracing
    x0 = zeros(N+1,1); % the last entry will represent the exponential growth rate
    x0(1:nc)=ones(1,nc);
    x0(end)=rd;

    % Solve system (3.2)-(3.3) in the main text
    Sol = fsolve(@(x) [x(1:N);1] - linear_contact_tracing(x(1:N),x(N+1),step,nc,nd,epsilon_c,beta_mat,h_d,surv_d), x0, opt_fsolve);
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

    % Prevented transmission
    Theta_d(index_y,index_x) = Rd/R0;
    Theta_ct(index_y,index_x) = Rct/R0;
    
end

    % Plot of effective reproduction number as a function of epsilon_c
    figure(20); hold on
    plot(epsilon_c_vector,Rct_matrix(:,index_x),'LineWidth',2,'Color',colorcode(index_x,:))
    legendname{index_x} = [num2str(delay_diagnosis),'d diagnosis delay'];
    
end

legendname{1} = 'no diagnosis delay';
figure(20)
legend(legendname)
xlabel('tracing coverage','Interpreter','latex');
ylabel('$R_{d,c}$','Interpreter','latex');
title('Effect of the tracing coverage on $R_{d,c}$','Interpreter','latex');
set(gca,'fontsize',14)
