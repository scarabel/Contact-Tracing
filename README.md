# Contact-Tracing
Codes for the renewal equation model of contact tracing

The codes in this repository are used to reproduce the figures in the paper:
Scarabel F, Pellis L, Ogden NH, Wu J. A renewal equation model to assess roles and limitations of contact tracing for disease outbreak control. Royal Society Open Science, 2021.
Please cite this reference if you use the codes.
The codes are distributed under the MIT license, see LICENSE for licensing information. 

The following scripts are available:
- 'Fig2_controllability.m': reproduces Fig 2;
- 'Fig3_Fig4_tracing_window.m': reproduces Fig 3 and Fig 4;
- 'Fig5_constraints.m': reproduces Fig 5;
- 'Fig6_R0_ec.m': reproduces Fig 6;
- 'FigSM2_delay_diagnosis.m': reproduces Supplementary Material, Fig 2;
The script use the function 'linear_contact_tracing.m', which solves the linearized equation.
The folder 'Figures' contains the figures included in the paper.

The codes are tested on MATLAB R2020b.