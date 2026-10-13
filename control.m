%% ========================================================================
%  1. NMPC DIMENSIONS
%  ========================================================================

nx = 6;       % states
ny = 2;       % controlled outputs
nu = 2;       % manipulated variables: [beta_l; beta_f]

nlobj = nlmpc(nx,ny,nu);


%% ========================================================================
%  2. SAMPLE TIME AND HORIZONS
%  ========================================================================

% Values retained from the tested orbit NMPC setup.
nlobj.Ts = 40;                    % controller sample time [s]
nlobj.PredictionHorizon = 60;
nlobj.ControlHorizon = 3;


%% ========================================================================
%  3. NONLINEAR PREDICTION MODEL
%  ========================================================================

nlobj.Model.StateFcn = 'nmpc_orbit_state';
nlobj.Model.OutputFcn = 'nmpc_orbit_output';
nlobj.States(1).Max = inf;
nlobj.States(1).Min = -inf;
% nmpc_orbit_state returns dx/dt, therefore the prediction model is
% continuous-time and nlmpc performs the prediction integration internally.
nlobj.Model.IsContinuousTime = true;

% One optional model parameter is used. It is a 5x1 VECTOR:
%
%   p(k) = [a_f(k); e_f(k); i_f(k); omega_f(k); rho_f(k)]
%
% NumberOfParameters = 1 because this entire vector is ONE parameter.
% Its numerical value is supplied online through the Simulink params port.
nlobj.Model.NumberOfParameters = 1;


%% ========================================================================
%  4. MANIPULATED-VARIABLE CONSTRAINTS
%  ========================================================================

% u(1) = beta_l
nlobj.MV(1).Min = B_l_min;
nlobj.MV(1).Max = B_l_max;

% u(2) = beta_f
nlobj.MV(2).Min = B_f_min;
nlobj.MV(2).Max = B_f_max;

% Sanity checks
assert(B_l_min < B_l_max, 'Leader beta bounds are invalid.');
assert(B_f_min < B_f_max, 'Follower beta bounds are invalid.');


%% ========================================================================
%  5. NMPC WEIGHTS
%  ========================================================================

% Output weights:
% y(1) = delta_a
% y(2) = delta_lambda
%
nlobj.Weights.OutputVariables = [2.8e3 1e3];

% Manipulated-variable magnitude weights.
nlobj.Weights.ManipulatedVariables = [0.009 0.009];

% Manipulated-variable RATE weights
nlobj.Weights.ManipulatedVariablesRate = [9.4e-3 9.4e-3];




% Initial state and manipulated variables
x0 = [-5.423641981982115e-09;
1.746083681375150e-05; 
1.834066508843608e-08; 
6.129582174729115e-09; 
-7.527684364738718e-09; 
2.745498157087124e-09];

u0 = [0.04; 0.04];

% NMPC model parameters
a_f_nmpc     = a_follower;
e_f_nmpc     = e_follower;
i_f_nmpc     = deg2rad(i_follower);
omega_f_nmpc = deg2rad(omega_follower);
rho_f        = rho_follower;

% Single parameter vector used by nmpc_orbit_state and the Simulink bus
P = [a_f_nmpc;
     e_f_nmpc;
     i_f_nmpc;
     omega_f_nmpc;
     rho_f];

% Test the state function with the same parameter structure used by NMPC
dx0 = nmpc_orbit_state(x0,u0,P);

% Since the model has ONE parameter argument (P), pass one cell containing P
validateFcns(nlobj, x0, u0, [], {P});

% Create the parameter bus for the Nonlinear MPC Controller block
open('system_model_NMPC.slx');
mdl = bdroot;
createParameterBus(nlobj, [mdl '/Model/control_nmpc/Nonlinear MPC Controller'], 'BusObject', {P});


