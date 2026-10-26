function [f,controller_state_tpdt,outputs] = control_law(t,dt,user_input,vehicle_state,controller_state,parameters,outputs)
%CONTROL_LAW  The control law proper, as it would be implemented on a vehicle.
%   [f,xd_tpdt,outputs] = CONTROL_LAW(t,dt,user_input,vehicle_state, ...
%       controller_state,parameters,outputs)
%   runs the full controller cascade once with a step of length dt and
%   returns the resulting individual actuator forces f.
%
%   The cascade mirrors Sections III and IV:
%       reference_model          generates the reference outputs x* and the
%                                pointing vector nu from user_input
%       outer_loop_controller    computes the target actuation a* = kappa_x*
%                                and the matrices that supply ddt_a* and
%                                ddt2_a*, needed by the inner loop
%       inner_loop_controller    computes the thrust and angular
%                                acceleration prescriptions u1 and u2
%       control_allocation       distributes [u1; u2] across the individual
%                                actuators to give f
%
%   This function holds the portion of the control system that would be
%   implemented on the vehicle.  Generating the user command from a
%   mission-level input, and mapping f back to the net thrust and torque
%   that the vehicle model in this paper assumes, are handled by
%   CONTROL_LAW_FOR_SIMULATION, which is the version used in simulation.
%   Neither function performs controller state switching; the MRP
%   switching rules live in the hybrid jump guard and jump map.
%
%   Controller state layout (fixed order; the returned increment is
%   assembled in the same order as it is parsed):
%       xd(1:7)      inner-loop state:     [mu_star; qeta_lp]
%       xd(8:11)     control allocation state
%       xd(12:29)    reference-model state
%       xd(30:end)   outer-loop state
%
%   Inputs:
%       t, dt              current time and flow step length
%       user_input         6x1 command: reference position, then the
%                          desired pointing vector nu
%       vehicle_state      outer-loop and attitude states, layout as in
%                          MULTIROTOR_FLOW_MAP
%       xd                 controller state, layout above
%       parameters         struct of system and controller parameters
%       outputs            struct of diagnostics, threaded through so that
%                          one struct accumulates every stage's outputs
%
%   Outputs:
%       f                  4x1 individual actuator forces
%       xd_tpdt            controller state at t + dt
%       outputs            diagnostics struct, extended by every stage
%
%   Requires on the MATLAB path:
%       reference_model, outer_loop_controller, inner_loop_controller,
%       control_allocation.
%
%   Called by CONTROL_LAW_FOR_SIMULATION.

% Parse controller states...
IL_controller_state = controller_state(1:7,1);
CA_controller_state = controller_state(8:11,1);
reference_state = controller_state(12:29,1);
OL_controller_state = controller_state(30:end,1);

% Run reference model...
[discrete_reference_state_tpdt,reference_output,outputs] = ...
    reference_model(t,dt,reference_state,user_input,parameters,outputs);

% Run outer-loop controller...
[a_star,ddt_a_star_MATRIX,ddt2_a_star_MATRIX,OL_controller_state_tpdt,outputs] = ...
    outer_loop_controller(t,dt,vehicle_state,reference_output,OL_controller_state,parameters,outputs);

% Run inner-loop controller...
[u_baseline,discrete_IL_controller_state_tpdt,outputs] = ...
    inner_loop_controller(t,dt,vehicle_state,reference_output,IL_controller_state,a_star,ddt_a_star_MATRIX,ddt2_a_star_MATRIX,parameters,outputs);

% Run control-allocation function...
[f,discrete_CA_controller_state_tpdt,outputs] = ...
    control_allocation(u_baseline,CA_controller_state,parameters,outputs);

% Assemble controller state increment...
controller_state_tpdt = [discrete_IL_controller_state_tpdt; discrete_CA_controller_state_tpdt; discrete_reference_state_tpdt; OL_controller_state_tpdt];

end