function [xcdottj,xdtpdtj,outputs] = multirotor_flow_map(tj,dt,xctj,xdtj,parameters)
%MULTIROTOR_FLOW_MAP  Flow map for the multirotor hybrid system.
%   [xcdot,xd_tpdt,outputs] = MULTIROTOR_FLOW_MAP(tj,dt,xc,xd,parameters)
%   returns the time derivative of the continuous state at the start of a 
%   flow step and the discrete state at the end of a flow step of length 
%   dt, for a multirotor under the appended practical complete tracking 
%   control law. 
%
%   Continuous state layout (NED):
%       xc(1:3)    x_I       inertial position                     [m]
%       xc(4:6)    v_I       inertial velocity                     [m/s]
%       xc(7:10)   q_IB      unit quaternion, scalar first
%       xc(11:13)  omega_B   body angular velocity                 [rad/s]
%
%   Inputs:
%       tj         1x2 hybrid time [t j]
%       dt         length of the trial flow step
%       xc         13x1 continuous state, layout above
%       xd         discrete controller state, layout as in CONTROL_LAW
%       parameters struct of system and controller parameters; uses
%                  parameters.m, parameters.In, parameters.g
%
%   Outputs:
%       xcdot      13x1 continuous state derivative
%       xd_tpdt    discrete controller state at the end of the step
%       outputs    struct of controller diagnostics, accumulated by
%                  CONTROL_LAW_FOR_SIMULATION at this state
%
%   Requires on the MATLAB path:
%       control_law_for_simulation, attitude_quaternion_to_matrix,
%       attitude_quaternion_kinematic_matrix, unit_vector.

% Parse inputs...
t = tj(1);
j = tj(2);
continuous_system_state = xctj;
discrete_system_state = xdtj;

% Load parameter values...
m = parameters.m;
In = parameters.In;
g = parameters.g;

% Separate vehicle states...
vehicle_state = continuous_system_state;
x_I = vehicle_state(1:3); % Inertial frame cartesian position
v_I = vehicle_state(4:6); % Velocity of the center of mass relative to the inertial frame resolved in the inertial frame
q_IB = vehicle_state(7:10); q_IB = q_IB/norm(q_IB); % Unit quaternion corresponding to R_IB
omega_B = vehicle_state(11:13); % Angular velocity of the body frame relative to the inertial frame resolved in the body frame

% Call simulation version of control law...
[u,controller_state_tpdt,outputs] = ...
    control_law_for_simulation(t,dt,continuous_system_state,discrete_system_state,parameters);

% Vehicle state derivatives...
ddt_x_I = v_I;
ddt_v_I = [0;0;g] + attitude_quaternion_to_matrix(q_IB)*(-unit_vector(3,3))*u(1,1)/m;
ddt_q_IB = attitude_quaternion_kinematic_matrix(q_IB)*omega_B;
ddt_omega_B = (In\(cross(In*omega_B,omega_B) + u(2:4,1)));
ddt_vehicle_state = [
    ddt_x_I
    ddt_v_I
    ddt_q_IB
    ddt_omega_B
];

% Assemble continuous system state derivative...
ddt_continuous_system_state = ddt_vehicle_state;
xcdottj = ddt_continuous_system_state;

% Assemble discrete system state increment...
discrete_system_state_tpdt = controller_state_tpdt;
xdtpdtj = discrete_system_state_tpdt;

end