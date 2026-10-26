function [a_star,ddt_a_star_MATRIX,ddt2_a_star_MATRIX,OL_controller_state_tpdt,outputs] = outer_loop_controller(t,dt,vehicle_state,reference_output,OL_controller_state,parameters,outputs)
%OUTER_LOOP_CONTROLLER  Outer-loop control law and exact derivative data for a*.
%   [a_star,ddt_MATRIX,ddt2_MATRIX,xd_OL_tpdt,outputs] = ...
%       OUTER_LOOP_CONTROLLER(t,dt,vehicle_state,reference_output, ...
%       OL_controller_state,parameters,outputs)
%   computes the target actuation a* = kappa_x* that the inner loop is asked
%   to realize, along with the coefficient matrices that give its first two
%   time derivatives.
%
%   In state-feedback mode ("P" in RUN_ALL) the law is exactly the one given 
%   in the paper,
%       a* = m ddot_x* - m g e_3 - [k_P I  k_D I] e,   e = x_OL - x*,
%   with the mass scaling supplied by invCAB_OL rather than by an explicit
%   1/m. In the integral and adaptive modes the controller state
%   alpha_hat_OL adds the remaining terms, and Gamma, gamma_nu_OL, and
%   gamma_sigma_OL are nonzero.
%
%   Why matrices rather than values: The inner loop needs ddt_a* and
%   ddt2_a* to evaluate omega* and ddt_omega*. Both are affine in the
%   actuation actually applied, a = lambda u1, and in its derivative,
%       ddt_a*  = ddt_a_star_MATRIX * [1; a]            (3x4 acting on 4x1)
%       ddt2_a* = ddt2_a_star_MATRIX * [1; a; ddt_a]    (3x7 acting on 7x1)
%   so this function builds the two coefficient matrices, and the inner
%   loop supplies a and ddt_a at the point where they are known. This keeps
%   the derivatives exact, with no numerical differentiation and no
%   symbolic algebra at run time, at the cost of the inner loop having to
%   evaluate the products itself.
%
%   Outer-loop state and jets:
%       x_OL = [x_I; v_I]              outer-loop state, 6x1
%       y_jet_01 = ref_out(1:6)        [x*; dx*]
%       y_jet_02 = ref_out(1:9)        [x*; dx*; ddx*]
%       y_jet_12 = ref_out(4:9)        [dx*; ddx*]
%       y_jet_13 = ref_out(4:12)       [dx*; ddx*; dddx*]
%       y_jet_24 = ref_out(7:15)       [ddx*; dddx*; ddddx*]
%   The jet notation is the same as in REFERENCE_MODEL: the subscript pair
%   in y_jet_ij names the derivative orders spanned.
%
%   Controller state:
%       alpha_hat_OL   integral or adaptive state, 3x1. Identically zero in
%                      state-feedback mode, where ddt_alpha_hat_OL = 0.
%
%   The returned state is advanced with a single forward-Euler increment,
%       alpha_hat_OL_tpdt = alpha_hat_OL + ddt_alpha_hat_OL*dt,
%   which is exact in state-feedback mode (the rate is zero) and a
%   first-order step otherwise. This differs from REFERENCE_MODEL, which
%   uses an exact zero-order-hold discretization.
%
%   Inputs:
%       t, dt              current time and flow step length; t is unused
%       vehicle_state      layout as in MULTIROTOR_FLOW_MAP; uses
%                          vehicle_state(1:6)
%       reference_output   21x1 output of REFERENCE_MODEL; uses 1:15
%       OL_controller_state 3x1 alpha_hat_OL
%       parameters         struct of system and controller parameters; uses
%                          parameters.m, parameters.g, and the outer-loop
%                          matrices and gains A_OL, B_OL, C_OL, invCAB_OL,
%                          K_P_OL, K_D_OL, Gamma_OL, gamma_nu_OL,
%                          gamma_sigma_OL
%       outputs            diagnostics struct, threaded through
%
%   Outputs:
%       a_star             3x1 target actuation, the signal a* of the paper
%       ddt_MATRIX         3x4 coefficient matrix for ddt_a*
%       ddt2_MATRIX        3x7 coefficient matrix for ddt2_a*
%       xd_OL_tpdt         3x1 alpha_hat_OL at t + dt
%       outputs            diagnostics struct, extended with alpha_hat_OL,
%                          OL_controller_state_tpdt, target_thrust, a_star,
%                          and the five reference-jet entries
%
%   Requires on the MATLAB path: nothing beyond base MATLAB.
%
%   Called by CONTROL_LAW.
%
%   Note: MULTIROTOR_JUMP_GUARD and MULTIROTOR_LYAPUNOV_FUNCTION both
%   reimplement this law inline rather than calling it, so that they can
%   evaluate V and the guard without invoking the rest of the cascade. Any
%   change to the law here must be mirrored in both of those files, or the
%   guard and the Lyapunov function will silently disagree with the control
%   actually applied. 
    
% Load parameter values...
m = parameters.m;
g = parameters.g;
d_OL = [0;0;m*g];
A_OL = parameters.A_OL;
B_OL = parameters.B_OL;
C_OL = parameters.C_OL;
invCAB_OL = parameters.invCAB_OL;
K_P = parameters.K_P_OL;
K_D = parameters.K_D_OL;
Gamma = parameters.Gamma_OL;
gamma_nu_OL = parameters.gamma_nu_OL;
gamma_sigma_OL = parameters.gamma_sigma_OL;
C_jet_2 = [C_OL; C_OL*A_OL];
C_jet_3 = [C_OL; C_OL*A_OL; C_OL*A_OL^2];
KPD1 = [K_P K_D eye(3)];


% Separate vehicle states...
x_I = vehicle_state(1:3); % Inertial frame cartesian position
v_I = vehicle_state(4:6); % Velocity of the center of mass relative to the inertial frame resolved in the inertial frame

% Separate reference states...
y_ref = reference_output(1:3);
ddt_y_ref = reference_output(4:6);
ddt2_y_ref = reference_output(7:9);
ddt3_y_ref = reference_output(10:12);
ddt4_y_ref = reference_output(13:15);
y_jet_01 = reference_output(1:6);
y_jet_02 = reference_output(1:9);
y_jet_12 = reference_output(4:9);
y_jet_13 = reference_output(4:12);
y_jet_24 = reference_output(7:15);

% Parse controller states...
alpha_hat_OL = OL_controller_state;

%% Outer-loop Controller

% Load outer-loop control states...
x_OL = [x_I; v_I];

% Compute outer-loop control signal (target control)...
a_star = -d_OL + invCAB_OL*(alpha_hat_OL - KPD1*C_jet_3*x_OL + KPD1*y_jet_02);

% Prepare to compute ddt_x_OL...
ddt_x_OL_MATRIX = [A_OL*x_OL + B_OL*d_OL   B_OL];
    % ddt_x_OL = ddt_x_OL_MATRIX*[1; a];
    
% Prepare to compute ddt2_x_OL...
ddt2_x_OL_MATRIX = [A_OL^2*x_OL + A_OL*B_OL*d_OL   A_OL*B_OL   B_OL];
    % ddt2_x_OL = ddt2_x_OL_MATRIX*[1; a; ddt_a];

% Prepare to compute ddt_a_star...
ddt_a_star_MATRIX = invCAB_OL*(...
    [Gamma*(y_jet_01 - C_jet_2*x_OL) + KPD1*y_jet_13 zeros(3,3)] ...
    -KPD1*C_jet_3*ddt_x_OL_MATRIX ...
    );
    % ddt_a_star = ddt_a_star_MATRIX*[1; a];

% Prepare to compute ddt2_a_star...
ddt2_a_star_MATRIX = invCAB_OL*(...
    [Gamma*y_jet_12 + KPD1*y_jet_24 zeros(3,6)] ...
    - Gamma*C_jet_2*[ddt_x_OL_MATRIX zeros(6,3)] ...
    - KPD1*C_jet_3*ddt2_x_OL_MATRIX ...
    );
    % ddt2_a_star = ddt2_a_star_MATRIX*[1; a; ddt_a];

% Assemble outer-loop controller state rate...
ddt_alpha_hat_OL = -Gamma*(x_OL - y_jet_01) - (gamma_sigma_OL + gamma_nu_OL*norm(x_OL - y_jet_01))*alpha_hat_OL; 

% Calculate outer-loop controller state increment...
OL_controller_state_tpdt = alpha_hat_OL + ddt_alpha_hat_OL*dt;

% Manage outputs...
outputs.alpha_hat_OL = alpha_hat_OL;
outputs.OL_controller_state_tpdt = OL_controller_state_tpdt;
outputs.target_thrust = a_star;
outputs.a_star = a_star;
outputs.x_I_ref = y_ref;
outputs.ddt_x_I_ref = ddt_y_ref;
outputs.ddt2_x_I_ref = ddt2_y_ref;
outputs.ddt3_x_I_ref = ddt3_y_ref;
outputs.ddt4_x_I_ref = ddt4_y_ref;

end