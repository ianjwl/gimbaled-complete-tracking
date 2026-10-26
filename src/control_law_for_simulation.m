function [u,controller_state_tpdt,outputs] = control_law_for_simulation(t,dt,continuous_system_state,discrete_system_state,parameters)
%CONTROL_LAW_FOR_SIMULATION  Complete control law, in the sense of the theory.
%   [u,xd_tpdt,outputs] = CONTROL_LAW_FOR_SIMULATION(t,dt,xc,xd,parameters)
%   runs the control law against a user command and returns the net thrust
%   and torque..
%
%   This is the version of the control law used in the simulation.  It
%   consists of three parts:
%       1. user_input_function   builds the user command: the reference
%                                position history and the desired pointing
%                                vector nu
%       2. control_law           runs the controller cascade and returns
%                                the individual actuator forces f
%       3. mixing_matrix         recombines f into the net thrust and
%                                torque u, the variables in which the
%                                theory states the control law
%
%   The first and third parts are what the theory treats as given by
%   context rather than designed, and they depend on the simulated example
%   (the reference history and the actuator geometry).  The controller
%   cascade itself is in CONTROL_LAW, which is the portion that could be
%   implemented on a vehicle.  Neither function performs controller state
%   switching; the MRP switching rules live in the hybrid jump guard and
%   jump map.
%
%   u is ordered as [u1; u2]: u(1) is the collective thrust and u(2:4) is
%   the torque, so the output is directly compatible with the vehicle
%   dynamics in MULTIROTOR_FLOW_MAP.
%
%   Inputs:
%       t, dt              current time and flow step length
%       xc                 continuous vehicle state, layout as in
%                          MULTIROTOR_FLOW_MAP
%       xd                 discrete controller state, layout as in
%                          CONTROL_LAW
%       parameters         struct of system and controller parameters
%
%   Outputs:
%       u                  4x1 net thrust and torque [u1; u2]
%       xd_tpdt            controller state at t + dt
%       outputs            diagnostics struct; outputs.initialized is set
%                          true on entry and outputs.user_input on exit
%
%   Requires on the MATLAB path:
%       user_input_function, control_law.
%
%   Called by MULTIROTOR_FLOW_MAP.

% Load parameter values...
mixing_matrix = parameters.mixing_matrix;

% Separate system states...
vehicle_state = continuous_system_state;
controller_state = discrete_system_state;

% Initialize outputs structure...
outputs.initialized = true;

% Call user input function...
user_input = user_input_function(t,vehicle_state,parameters);

% Call control law...
[f,controller_state_tpdt,outputs] = ...
    control_law(t,dt,user_input,vehicle_state,controller_state,parameters,outputs);

% Calculate control action...
u = mixing_matrix*f;

% Manage outputs...
outputs.user_input = user_input;

end