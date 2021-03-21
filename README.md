# Contact-Tracing
Matlab codes for the renewal equation model of contact tracing

The codes in this repository are used to reproduce the figures in the paper:
Scarabel F, Pellis L, Ogden NH, Wu J. A renewal equation model to assess roles and limitations of contact tracing for disease outbreak control. Royal Society Open Science, 2021.
Please cite this reference if you use the codes.
The codes are distributed under the MIT license, see LICENSE for licensing information. 

The following scripts are available:
- 'Fig2_controllability.m': reproduces Fig 2;
- 'Fig3_Fig4_tracing_window.m': reproduces Fig 3 and Fig 4;
- 'Fig5_constraints.m': reproduces Fig 5;
- 'Fig6_R0_ec.m': reproduces Fig 6;
- 'FigSM2_delay_diagnosis.m': reproduces Supplementary Material, Fig 2.

The script use the function 'linear_contact_tracing.m', which solves the linearized equation.
The folder 'Figures' contains the figures included in the paper.

The codes are tested on MATLAB R2020b and use the function fsolve.m from Matlab Optimization Toolbox for the solution of a nonlinear system.

## Method
System (2.9) with (2.11) and system (3.2)-(3.3) are solved numerically by discretizing the domain of the tau variable with equidistant nodes in the interval [0,B] containing the (compact) support of the function beta. The integrals are then approximated using rectangles quadrature formulas.

The linearised system (3.2)-(3.3) is solved by using a numerical solver (fsolve for Matlab) to find a solution of the nonlinear algebraic system of equations for the vector approximating h_c and for the approximation of the variable r. 

The nonlinear system (2.9) with (2.11) is solved iteratively for increasing time, using equidistant nodes in the time interval of the same size of the tau-stepsize. For each time-iteration, the tau-iteration is solved for decreasing indices, from B to 0.
