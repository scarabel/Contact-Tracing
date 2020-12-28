function output = linear_contact_tracing(h_c,r,step,nc,nd,epsilon_c,beta_mat,h_d,surv_d)
% Copyright (c) 2020 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE.txt for 
% licensing information. 
% 
% If using this code, please cite 
% Scarabel, Pellis, Ogden, Wu, 'A renewal equation model to assess roles and
% limitations of contact tracing for disease outbreak control'
% 
%% function linear_contact_tracing.m
% computes the right-hand side of the system in the linear approximation
% h_c vector of the contact tracing hazard rate (unknown in the system)
% r - rate of exponential growth (unknown in the system)
% step - stepsize
% nc - index corresponding to maximal contact tracing age
% nd - index corresponding to maximal diagnosis (nd>nc)
% epsilon_c - fraction of traced contacts
% beta_mat - vector with transmissibility values
% h_d - vector with diagnosis hazard rate
% surv_d - vector with survival probabilities to diagnosis


    N=length(h_c);
    output_hc = zeros(N,1);
    surv_ct = zeros(N,1);

    for itau = 1:N
        surv_ct(itau) = exp(-step*sum(h_c(1:itau)));
    end

    for itau = nc:-1:1 % loop over age (decreasing)

        int_ct = 0;
        for isigma = itau+1:nc % integral up to nc*step (with nd>nc)

            idiff = isigma-itau;
            % rectangles quadrature formula
            int_ct = int_ct + step*(beta_mat(idiff)*exp(-r*(step*idiff))...
                *surv_d(isigma)*surv_ct(isigma)*(h_d(isigma)+h_c(isigma)) );

        end
        
        for isigma = nc+1:nd % integral up to nd*step (with nd>nc) and hc = 0

            idiff = isigma-itau;
            % rectangles quadrature formula
            int_ct = int_ct + step*(beta_mat(idiff)*exp(-r*(step*idiff))...
                *surv_d(isigma)*surv_ct(isigma)*h_d(isigma));
            
        end
        output_hc(itau) = epsilon_c*int_ct;
                        
    end
    
    % integral for Lotka-Euler equation
    int_inc = 0;
    for is=1:N
        int_inc = int_inc + step*...
            (beta_mat(is)*exp(-r*(step*is))*surv_d(is)*surv_ct(is));
    end
    
    output = [output_hc; int_inc];

end