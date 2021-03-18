% Copyright (c) 2021 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE for 
% licensing information. 
% 
% If you use this code, please cite 
% Scarabel, Pellis, Ogden, Wu, 'A renewal equation model to assess roles and
% limitations of contact tracing for disease outbreak control',
% Royal Society Open Science, 2021.
% 
%% Script controllability.m 
% Computes the lines R_d=1 and R_{d,c}=1 in the plane (epsilon_d,delay),
% for given value of R0.
% Uses the function linear_contact_tracing.m
% The two panels in the figure are obtained for R0 = 1.5 and R0 = 2.5

clear
close all

step = 0.05; % stepsize for numerical solution

% Epidemiological parameters

% Basic reproduction number
R0 = 2.5; % 2.5; % 1.5;

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

% Plot of disease-related parameters
bgrid = 0:0.1:bmax;
for jj=1:length(bgrid)
    B(jj) = beta_transm(bgrid(jj))/R0;
end

figure(10); hold on
plot(bgrid,density_incubation_f(bgrid),bgrid,B,'LineWidth',2)
legend('incubation time','generation time')
xlabel('age since infection','Interpreter','latex');
title('Disease-related parameters','Interpreter','latex')
set(gca,'fontsize',14)

% for controllability: 
epsilon_s = 1; % percentage of symptomatic individuals

% max diagnosis time
dmax = 20;

% max tracing window
cmax = 5;

%% Iterations

epsilon_d_vector = 0.1:0.3:1; 
delay_diagnosis_vector = 0:1:5;
epsilon_c_vector = 0.2:0.2:1;

xgrid = epsilon_d_vector;
ygrid = delay_diagnosis_vector;
zgrid = epsilon_c_vector;

Lx = length(xgrid);
Ly = length(ygrid);
Lz = length(zgrid);

Rd_matrix = zeros(Ly,Lx);
Rct_matrix = zeros(Ly,Lx,Lz);

for index_z = 1:Lz
for index_x = 1:Lx
for index_y = 1:Ly
   
    epsilon_d = epsilon_d_vector(index_x);
    delay_diagnosis = delay_diagnosis_vector(index_y);
    epsilon_c = epsilon_c_vector(index_z);

    % Diagnosis process
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
    x0(1:nc)=ones(1,nc);
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

    Rd_matrix(index_y,index_x) = Rd;
    Rct_matrix(index_y,index_x,index_z) = Rct;

end
end
end

%% Plot

[XX,YY]=meshgrid(xgrid,ygrid); 
colorscode = jet(Lz);

figure
contour(XX,YY,Rd_matrix,[1 1],'ShowText','off','LineWidth',2,'LineColor','b'); hold on %,'DisplayName','R_d'
for index_z = 1:Lz
    contour(XX,YY,Rct_matrix(:,:,index_z),[1 1],'ShowText','off','LineWidth',2,'LineColor',colorscode(index_z,:)); hold on
end

xlabel('fraction effectively diagnosed','Interpreter','latex');
ylabel('delay from symptoms to diagnosis','Interpreter','latex');
title(['$R_0=$',num2str(R0)],'Interpreter','latex');
legend('$R_d=1$','$R_{d,c}=1$, 20\%', '$R_{d,c}=1$, 40\%', '$R_{d,c}=1$, 60\%','$R_{d,c}=1$, 80\%','$R_{d,c}=1$, 100\%', 'Location', 'NorthWest','Interpreter','latex')
set(gca,'fontsize',14)