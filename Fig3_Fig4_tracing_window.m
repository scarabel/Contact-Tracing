% Copyright (c) 2021 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE for 
% licensing information. 
% 
% If you use this code, please cite 
% Scarabel, Pellis, Ogden, Wu, 'A renewal equation model to assess roles and
% limitations of contact tracing for disease outbreak control',
% Royal Society Open Science, 2021.
% 
%% script tracing_window.m
% Computes R_{d,c} and related quantities varying epsilon_c and cmax.
% Uses the function linear_contact_tracing.m
% Fig 3 in the main text is obtained with R0 = 2.5 and delay_diagnosis = 0

clear
close all

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
delay_diagnosis = 0; % 2; % 0;
epsilon_d = 1; % fraction of symptomatic individuals diagnosed

density_diagnosis = @(x) epsilon_d*epsilon_s*(x<=dmax).*gampdf(x-delay_diagnosis,shape_incubation,scale_incubation);
surv_diagnosis_f = @(x) 1-integral(@(y) density_diagnosis(y),0,x);

%% Iterations to explore different combinations of parameters

legendname = {};
epsilon_c_vector = 0:0.1:1;
cmax_vector = 2:5;

Rct_matrix = zeros(length(epsilon_c_vector),length(cmax_vector));
rct_matrix = Rct_matrix;
Theta_d_matrix = Rct_matrix;
Thera_ct_matrix = Rct_matrix;

figure(5); hold on
plot([0 1],[1 1],'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',1,'HandleVisibility','off'); % reference line R=1

colorscode = lines(length(cmax_vector));

for index_x = 1:length(cmax_vector)
for index_y = 1:length(epsilon_c_vector)

    cmax = cmax_vector(index_x);
    epsilon_c = epsilon_c_vector(index_y);

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

    % Prevented transmission
    Theta_d(index_y,index_x) = Rd/R0;
    Theta_ct(index_y,index_x) = Rct/R0;
    
end

    % Plot of effective reproduction number as a function of epsilon_c
    figure(5); hold on
    plot(epsilon_c_vector,Rct_matrix(:,index_x),'LineWidth',2,'Color',colorscode(index_x,:))
    legendname{index_x} = [num2str(cmax),'d tracing'];
    
    % Plot of profile of h_ct for epsilon_c=1
    cdf_ct = 1-exp(-step*cumsum(h_ct)); % Cumulative distribution of contact trac
    pdf_ct = diff([0;cdf_ct])/step;
    cdf_diag = 1-exp(-step*cumsum(h_d)); % Cumulative distribution of contact trac
    cdf_diag_ct = 1-exp(-step*cumsum(h_d+h_ct)); % Cumulative distribution of contact trac
    pdf_diag_ct = diff([0;cdf_diag_ct])/step;
    
    figure(7); 
    subplot(2,2,1)
    hold on
    plot(step*(1:nc),h_ct(1:nc),'LineWidth',2,'Color',colorscode(index_x,:));
    plot(step*[nc,nc],[h_ct(nc),0],'LineWidth',2,'Color',colorscode(index_x,:),'HandleVisibility','off');
    axis([0 cmax 0 0.2])
    title('Tracing rate','Interpreter','latex')

    subplot(2,2,2)
    hold on
    plot(step*(1:nd),h_d(1:nd),'LineWidth',2,'Color',[0.5 0.5 0.5]);
    plot(step*[nd,nd],[h_d(nd),0],'LineWidth',2,'Color',[0.5 0.5 0.5],'HandleVisibility','off');
    plot(step*(1:nd),h_d(1:nd)+h_ct(1:nd),'LineWidth',2,'Color',colorscode(index_x,:));
    plot(step*[nd,nd],[h_d(nd)+h_ct(nd),0],'LineWidth',2,'Color',colorscode(index_x,:),'HandleVisibility','off');
    title('Detection rate','Interpreter','latex')

    subplot(2,2,3)
    hold on
    plot(step*(1:nc),cdf_ct(1:nc),'LineWidth',2,'Color',colorscode(index_x,:));
    plot(step*[nc,N],[cdf_ct(nc),cdf_ct(N)],'LineWidth',2,'Color',colorscode(index_x,:),'HandleVisibility','off');
    axis([0 cmax 0 0.5])
    title('Tracing cdf','Interpreter','latex')
    xlabel('age since infection','Interpreter','latex');

    subplot(2,2,4); 
    hold on
    if index_x==1
        plot(step*(1:nd),cdf_diag(1:nd),'LineWidth',2,'Color',[0.5 0.5 0.5]);
    end
    plot(step*(1:nd),cdf_diag_ct(1:nd),'LineWidth',2,'Color',colorscode(index_x,:));
    title('Detection cdf','Interpreter','latex')
    xlabel('age since infection','Interpreter','latex');

end

%%
figure(5)
legend(legendname)
xlabel('tracing coverage','Interpreter','latex');
ylabel('$R_{d,c}$','Interpreter','latex');
title('Effect of the tracing coverage on $R_{d,c}$','Interpreter','latex');
set(gca,'fontsize',14)

figure(7)
legenddiag{1} = 'no tracing';
for jj=1:4
    legenddiag{jj+1} = legendname{jj};
end
subplot(2,2,4)
legend(legenddiag,'Location','SouthEast')
set(gca,'fontsize',14)


%% Contour plots

[XX,YY]=meshgrid(cmax_vector,epsilon_c_vector);

% colorscode = jet(Lz);
figure
contour(XX,YY,Rct_matrix,0.2:0.1:2,'ShowText','on','LineWidth',2); hold on
xlabel('tracing window','Interpreter','latex');
ylabel('tracing coverage','Interpreter','latex');
title('Effective reproduction number','Interpreter','latex');
set(gca,'fontsize',14)

figure
contour(XX,YY,rct_matrix,-0.5:0.02:0.5,'ShowText','on','LineWidth',2); hold on
xlabel('tracing window','Interpreter','latex');
ylabel('tracing coverage','Interpreter','latex');
title('Growth rate','Interpreter','latex');
set(gca,'fontsize',14)

figure
[C,h] = contour(XX,YY,rct_matrix,[-0.5:0.02:0.5],'LineWidth',2); hold on
% % Uncomment the following for manual labelling:
% hcl = clabel(C,h,'manual','FontSize',12,'Color','k');
% for i = 1:length(hcl)
%     oldLabelText = get(hcl(i), 'String');
%     doubl_time = round(log(2)./str2double(oldLabelText));
%     newLabelText = num2str(doubl_time);
%     set(hcl(i), 'String', newLabelText);
% end
xlabel('tracing window','Interpreter','latex');
ylabel('tracing coverage','Interpreter','latex');
title('Doubling time','Interpreter','latex');
set(gca,'fontsize',14)

figure
[C,h] = contour(XX,YY,Theta_d-Theta_ct,'LineWidth',2); hold on
% % Uncomment the following for manual labelling:
% hcl = clabel(C,h,'manual','FontSize',12,'Color','k');
% for i = 1:length(hcl)
%     oldLabelText = get(hcl(i), 'String');
%     percentage = str2double(oldLabelText)*100;
%     newLabelText = [num2str(percentage) ' %'];
%     set(hcl(i), 'String', newLabelText);
% end
xlabel('tracing window','Interpreter','latex');
ylabel('tracing coverage','Interpreter','latex');
title('Transmission prevented by contact tracing only','Interpreter','latex');
set(gca,'fontsize',14)
