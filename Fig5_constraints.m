% Copyright (c) 2021 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE for 
% licensing information. 
% 
% If you use this code, please cite 
% Scarabel, Pellis, Ogden, Wu, 'A renewal equation model to assess roles and
% limitations of contact tracing for disease outbreak control',
% Royal Society Open Science, 2021.
% 
%% Script constraints.m
% Simulations with maximal capacity on tracing or diagnosis
% The two panels in Fig 5 in the main text are obtained setting 
% max_capacity = 'CT' or 'Diag', respectively

clearvars

step = 0.05; % stepsize for numerical solution

% Epidemiological parameters

% Basic reproduction number
R0 = 1.5; 

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
epsilon_s = 0.85;

% Max diagnosis
dmax = 20;
delay_diagnosis = 2;

% Contact tracing process
cmax = 5;

%% Simulation with constraints

colorscode = lines(5);

% switch max_capacity for constraints either on diagnosis or tracing:
% 'CT' - constraints on tracing
% 'Diag' - constraints on diagnosis
max_capacity = 'CT'; %'Diag', 'CT'

Max_Test = 1e-7;
I0 = 1e-7; % multiplicative factor for initial condition

Tf = 100; % final time
T1 = 30; % time of interruption of contact tracing

epsilon_d_vector = [0.5;0.6;0.6];
duration_vector = [30;30;50]; % 50, 30  % duration of interruption of contact tracing

epsilon_c0 = 1; %0; %0.5;
epsilon_c_reduced = 0; % epsilon at interruption

% discretization
nd = dmax/step;
nc = cmax/step;
nb = bmax/step;

nt = Tf/step;
na = nt+nb;

N = max([nb,nd,nc,nc+nb]);

tgrid = step*(1:nt);
dgrid = step*(1:nd);
bgrid = step*(1:nb);
Ngrid = step*(1:N);

for ind_test = 1:3
    ind_color=1;
    
    duration = duration_vector(ind_test);
    epsilon_d = epsilon_d_vector(ind_test);
    
    % Diagnosis process
    density_diagnosis = @(x) epsilon_d*epsilon_s*(x<=dmax).*gampdf(x-delay_diagnosis,shape_incubation,scale_incubation);
    surv_diagnosis_f = @(x) 1-integral(@(y) density_diagnosis(y),0,x);

    % Initialization of known parameters
    beta_mat = zeros(N,1);
    h_d = zeros(N,1);
    surv_d = (1-epsilon_d)*ones(N,1); % survival diagnosis
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

    r0 = fzero(@(x) 1- step*trapz(beta_mat.*exp(-x*step*(1:N)')), 0.1);
    rd = fzero(@(x) 1- step*trapz(beta_mat.*surv_d.*exp(-x*step*(1:N)')), r0);

    %% Simulation of the nonlinear model (2.9) and (2.11) in the main text

    % definition of the initial history function
    incidence_f = @(x) I0*exp(rd*x);

    inc_history = incidence_f(-step*N:step:0);
    incidence_det = incidence_f(tgrid)';

    h_ct = zeros(nt,N);
    surv_ct = ones(nt,N);
    epsilon_c = epsilon_c0;
    Epsilon = zeros(nt,1);

    h_d_effective = repmat(h_d',nt,1);
    surv_d_effective = repmat(surv_d',nt,1);

    susceptibles_det = ones(nt,1);
    prevalence_det = zeros(nt,1);
    num_traced_det = zeros(nt,1);
    num_diagnosed_det = zeros(nt,1);

    % Loops
    for it = 2:nt % loop over time

        if it*step>T1 && it*step<=(T1+duration) % Interruption contact tracing
            epsilon_c_baseline = epsilon_c_reduced;
            epsilon_c = epsilon_c_baseline;
        else
            epsilon_c_baseline = epsilon_c0;
        end

        if it>1/step

            switch max_capacity
                % "ratio" is the proportion by which we divide
                % the proportion of contacts diagnosed/traced

                case 'CT' % max only on tracing
                    % ratio = min([1; Max_Test/(step*sum(num_traced_det(it-ceil(1/step):it-1)))]);
                    ratio = min([1; Max_Test/(0.5*sum(num_traced_det(it-2:it-1)))]);
                    epsilon_c = epsilon_c_baseline*ratio;   

                case 'Diag' % max only on diagnosis
                    % ratio = min([1; Max_Test/(step*sum(num_diagnosed_det(it-ceil(1/step):it-1)))]);
                    ratio = min([1; Max_Test/(0.5*sum(num_diagnosed_det(it-2:it-1)))]);

                    h_d_effective(it,:) = h_d'*ratio;
                    surv_d_effective(it,2:end) = surv_d_effective(it-1,1:end-1).*exp(-step*h_d_effective(it-1,1:end-1));
                    epsilon_c = epsilon_c_baseline;

                otherwise
                    error('unknown')                
            end

        end

        % update survival to contact tracing
        surv_ct(it,2:end)=surv_ct(it-1,1:end-1).*exp(-step*h_ct(it-1,1:end-1));

        for itau = nc:-1:1 % loop over age (decreasing) - we need to solve it only up to nc as the hazard of being traced is zero afterwards

            if it>itau
                FOI_delay = incidence_det(it-itau)/susceptibles_det(it-itau);
            else
                FOI_delay = inc_history(N+1+(it-itau));
            end

            integral_hc = 0; % initialization of integral
            for isigma = itau+1:N
                if it>isigma
                    inc_delay = incidence_det(it-isigma);
                else
                    inc_delay = inc_history(N+1+(it-isigma));
                end

                % rectangles quadrature formula
                integral_hc = integral_hc + step*(beta_mat(isigma-itau)*inc_delay...
                    *surv_d_effective(it,isigma)*surv_ct(it,isigma)*(h_d_effective(it,isigma)+h_ct(it,isigma)) );

            end

            % compute hazard of contact tracing
            % Assuming beta(0)=0
            if FOI_delay>0
                h_ct(it,itau) = epsilon_c*integral_hc/FOI_delay;
            else
                disp('zero FOI');
            end

        end

        % solution of equation for the incidence at time t
        int_inc = 0;
        int_prev = 0;
        int_traced = 0;
        int_diagnosed = 0;
        for is = 1:nb 
            idiff = it-is;
            if it>is 
                inc_delay = incidence_det(idiff);
            else
                inc_delay = inc_history(N+1+idiff);
            end

            int_inc = int_inc + step*...
                beta_mat(is)*inc_delay*surv_d_effective(it,is)*surv_ct(it,is); %*surv_d(is)*surv_ct(it,is);

            int_prev = int_prev + step*inc_delay; % not necessary for solution
            int_diagnosed = int_diagnosed + step*inc_delay*surv_d_effective(it,is)*surv_ct(it,is)*h_d_effective(it,is);
            int_traced = int_traced + step*inc_delay*surv_d_effective(it,is)*surv_ct(it,is)*h_ct(it,is);

        end

        prevalence_det(it) = int_prev;

        % additional piece for number of traced at each time 
        for is = nb+1:N
            idiff = it-is;
            if idiff>0
                inc_delay = incidence_det(idiff);
            else
                inc_delay = inc_history(N+1+idiff);
            end
            int_diagnosed = int_diagnosed + step*inc_delay*surv_d_effective(it,is)*surv_ct(it,is)*h_d_effective(it,is);
            int_traced = int_traced + step*inc_delay*surv_d_effective(it,is)*surv_ct(it,is)*h_ct(it,is);
        end

        num_traced_det(it) = int_traced;    
        num_diagnosed_det(it) = int_diagnosed;

        % update incidence and compute susceptibles
        if it==1
            incidence_det(it) = int_inc;
            susceptibles_det(it) = 1 - step*incidence_det(it);
        else
            incidence_det(it) = susceptibles_det(it-1)*int_inc;
            susceptibles_det(it) = susceptibles_det(it-1) - step*incidence_det(it);
        end
        
        % Plot profile h_ct and distribution
        if it==(10/step) || it==(T1+duration+10)/step
            
        cdf_ct = 1-exp(-step*cumsum(h_ct(it,:)))'; % Cumulative distribution of contact trac
        pdf_ct = diff([0;cdf_ct])/step;
        
        figure(50); subplot(3,1,ind_test);
        hold on 
        plot(step*(1:nc),h_ct(it,1:nc),'LineWidth',2,'LineStyle','--','Color',colorscode(ind_color,:),'HandleVisibility','off');
        plot(step*[nc,nc],[h_ct(it,nc),0],'LineWidth',2,'LineStyle','--','Color',colorscode(ind_color,:),'HandleVisibility','off');
        plot(step*(1:nc),cdf_ct(1:nc),'LineWidth',2,'Color',colorscode(ind_color,:));
        plot(step*[nc,N],[cdf_ct(nc),cdf_ct(N)],'LineWidth',2,'Color',colorscode(ind_color,:),'HandleVisibility','off');
        axis([0 cmax 0 0.3])
        ind_color = ind_color+1;
        end
        
    end
        
    figure(50)
    subplot(3,1,ind_test)
    legend('time $=$10', ['time $=$',num2str(T1+duration+10)],'Interpreter','latex','Location','NorthWest');
    set(gca,'fontsize',14)

    figure(5); hold on
    subplot(3,2,2*(ind_test-1)+1)
    plot(tgrid, incidence_det,'LineWidth',2); hold on
    plot(tgrid, num_diagnosed_det,'LineWidth',2); hold on
    plot(tgrid, num_traced_det,'LineWidth',2); hold on
    legend('Incidence','Diagnosed','Traced','Location','NorthWest');
    xlabel('time','Interpreter','latex')
    plot([0 Tf],[Max_Test Max_Test],'-.','HandleVisibility','off')

    subplot(3,2,2*ind_test)
    plot(tgrid(2:end), log(incidence_det(2:end)./incidence_det(1:end-1))/step,'LineWidth',2); hold on
    plot([0 Tf], [0 0],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle','--','HandleVisibility','off'); hold on
    xlabel('time','Interpreter','latex')
    ylabel('growth rate','Interpreter','latex')

end
        
switch max_capacity
    case 'CT'
        figure(5)
        sgtitle('Limit on tracing','Interpreter','latex');
        figure(50)
        sgtitle('Limit on tracing','Interpreter','latex');
    case 'Diag'
        figure(5)
        sgtitle('Limit on diagnosis','Interpreter','latex');
        figure(50)
        sgtitle('Limit on diagnosis','Interpreter','latex');
    otherwise
        display('unknown');
end

figure(50)
xlabel('age since infection','Interpreter','latex');