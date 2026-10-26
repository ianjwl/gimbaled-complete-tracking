function sim = run_all(save_figures,OLcontrolMode)
%RUN_ALL  Single entry point: simulate the closed-loop system and build every figure.
%   sim = RUN_ALL(save_figures,OLcontrolMode)
%   runs the multirotor simulation of the accompanying paper, postprocesses
%   the diagnostics, generates the paper figures and five diagnostic
%   figures, and returns the postprocessed data as a struct.
%
%   RUN_ALL sets the path (src/, scripts/, third_party/), so it must be
%   called from the repository root or with the repository root on the
%   path. Nothing else in the repository adds paths.
%
%   Sections, in order:
%       0   User configuration: selects the outer-loop control mode
%       1   Initial conditions: continuous state and controller state
%       2   Parameters: vehicle, reference model, controller gains, and
%           the mixing matrix; everything downstream reads from the
%           `parameters` struct built here
%       3   Simulation: HYBRID_RK with RK4 flows and the Lyapunov
%           step-size control active
%       4   Postprocessing: re-evaluates the flow map at every recorded
%           state to accumulate the diagnostics, then extracts the
%           histories the figures need and bundles them into `sim`
%       5   Plotting: five diagnostic figures, then the six paper figures
%
%   SIMULATION SETTINGS: Base step dT = 0.01 s on [0, 25] s, classical RK4
%   flows, event tolerance 1e-12, jump index capped at 1e5, and Lyapunov
%   step-size control with tol = 0.1*min(k_a,k_omega) and a relative
%   bound. The Lyapunov control shortens steps where the computed solution
%   violates Vdot <= -kV; it enforces accuracy by enforcing the certificate.
%
%   PARAMETERS: Every numeric value in Section 2 is a control parameter of
%   the paper. The vehicle is m = 2.09 kg with J = diag(0.0297, 0.0320,
%   0.0488) kg m^2. The controller uses r = 5 N and U1 = 0.55*r, so
%   rho = U1/r = 0.55. The reference model is a cascade of two
%   second-order filters per axis for the position channel plus a
%   second-order filter for the pointing channel, assembled with tf and
%   tf2ss and then replicated across the three axes with kron. The
%   outer-loop gains come from lqr. NOTE: the schematics in figure_2,
%   figure_3, and figure_4 use rho = 0.6 for illustration, which differs
%   from the simulation value above.
%
%   THE `sim` STRUCT: Section 4 bundles every variable in the workspace
%   into `sim` with a `who`/`eval` sweep, so `sim` contains far more than
%   the figures need, including the `parameters` struct, the four function
%   handles passed to HYBRID_RK, and the loop counters. The fields the
%   figure functions actually consume are:
%       sim.t                       time history                      [Kx1]
%       sim.positions               inertial position, NED            [Kx3]
%       sim.reference_positions     reference position, NED           [Kx3]
%       sim.attitudes               attitude quaternions              [Kx4]
%       sim.ENU_positions           position in East-North-Up         [Kx3]
%       sim.ENU_reference_positions reference in East-North-Up        [Kx3]
%       sim.ENU_Thrust              Tbar(t) = lambda(t) u1(t)/m       [Kx3]
%       sim.ENU_Thrust_star         reference Tbar(t)                 [Kx3]
%       sim.theta, sim.theta_star   angular offset and its bound      [Kx1]
%       sim.uprightness             (-e_3)'*lambda(t)                 [Kx1]
%   The remaining simulation outputs are available for inspection but are
%   not read by any figure.
%
%   FRAME CONVENTION: The plant state is NED, but the reference returned by
%   ALTERNATE_REFERENCE_FUNCTION is already in ENU as [Y; X; -Z], so it is
%   translated here by hand for plotting. The component swap [2 1 3] with a
%   sign flip on the third column appears throughout the extraction block
%   and is what puts both signals in the same frame.
%
%   OUTPUT FIGURES: Five diagnostic figures are drawn from the simulation
%   histories and exported as figures/diagnostic_figure_1.pdf through
%   diagnostic_figure_5.pdf: the timing history, the MRP error against its
%   bound, the alignment of the pointing vector, the Lyapunov value and its
%   attitude component, and the step-size history. These are not paper
%   figures; they are kept as a record of the solver's behavior, and the
%   step-size and Lyapunov traces are the two that matter if you are
%   auditing the integration.
%
%   Inputs:
%       save_figures   logical. When true, every figure is exported to PDF
%                      under figures/. When false, the figures are drawn
%                      and left open.
%       OLcontrolMode  string. Selects the outer-loop law. The valid 
%                      options are:
%                       "P"    state feedback. This is the mode used for 
%                              the results in the paper, and it is the 
%                              DEFAULT branch.
%                       "PI"   proportional-integral.
%                       "A"    adaptive.
%                      All three share the same inner loop, reference 
%                      model, and vehicle model; only the outer-loop gains 
%                      and the compensation state alpha_hat_OL differ. In 
%                      "P" mode the compensation rate is identically zero, 
%                      so the returned state is exact under the forward-
%                      Euler increment used in OUTER_LOOP_CONTROLLER.
%
%   Outputs:
%       sim            struct of postprocessed simulation data; see THE
%                      `sim` STRUCT above for the fields the figures use.
%                      Documentation on the contents of `sim` is available
%                      below under `Contents of sim.`  
%
%   Requires on the MATLAB path:
%       Everything in src/ and scripts/ (added here), plus
%       third_party/FlightTest-Tools-main for plotObjectTrajectory and the
%       quad.stl mesh used by figure_5 and figure_10. Also requires the
%       Control System Toolbox (tf, tf2ss, lqr) and, as currently written,
%       the Optimization Toolbox, because optimoptions is called for the
%       quadratic-program allocators even though apply_QP is false and none
%       of them is used. See KNOWN_ISSUES.md.
%
%   See also MULTIROTOR_FLOW_MAP, CONTROL_LAW_FOR_SIMULATION, CONTROL_LAW.
%
%
%% Contents of sim
%
%   RUN_ALL returns every variable in its workspace as a field of sim. The
%   fields fall into seven groups: the vehicle, the controller parameters,
%   the outer-loop design, the reference-model design, the control
%   allocation settings, the simulation settings, and the simulation output.
%   The output group is much the largest; it is split below into what the
%   integrator returns, what the diagnostic pass through the flow map adds,
%   and the coordinate histories that the figure functions read.
%
%   K = size(sim.tj_history,1) is the number of recorded states: one row per
%   accepted flow step, plus one row per jump. Every history below has K
%   rows. Time-varying quantities are Kx3 or Kx1 sampled at those rows.
%   A quantitity that is 1x1 is a scalar setting or parameter, constant for
%   the whole run.
%
%   =====================================================================
%   VEHICLE
%   =====================================================================
%     sim.m            scalar   mass                               [kg]
%     sim.In           3x3      inertia matrix, body frame         [kg m^2]
%     sim.g            scalar   gravitational acceleration         [m/s^2]
%
%   =====================================================================
%   CONTROLLER PARAMETERS (inner loop)
%   =====================================================================
%   -- Set-based objective
%     sim.r            scalar   radius r of the target ball A(t)        [N]
%     sim.U1min        scalar   lower bound U_1 on the thrust           [N]
%     sim.rho          scalar   U_1/r, the normalized lower bound       [-]
%
%   -- Angle geometry: breakpoints of theta*
%     sim.base         scalar   (1-rho)/2, the common scale of the four
%                               breakpoints                             [-]
%     sim.alpha_2pi    scalar   (1-rho)/8, top of the 2*pi plateau of
%                               theta*, below which the thrust direction
%                               is unconstrained                        [-]
%
%   -- Thrust law
%     sim.delta0       scalar   lower edge of the alignment transition
%                               in kappa                              [rad]
%     sim.delta1       scalar   upper edge of the same transition     [rad]
%     sim.a0           scalar   r*alpha_2pi/5, lower edge of the k_eta
%                               gate                                  [N]
%     sim.a1           scalar   r/5, upper edge of that gate          [N]
%
%   -- Backstepping and Lyapunov function
%     sim.varrho       scalar   containment gate fraction: the terms of
%                               omega* proportional to ddt_a*/||a*|| are
%                               active only where ||eta||/n* > varrho [-]
%     sim.Delta_a      scalar   width of the ramp in V_a, and of the
%                               vanishing step below it               [N]
%     sim.Delta_ddt_n_star  scalar   width of the ramp on ddt_n_star  [1/s]
%     sim.k_a          scalar   dissipation gain in Omega_a           [1/s]
%     sim.k_eta        scalar   attitude trim gain for keeping eta
%                               small once it reaches N*              [rad/s]
%     sim.k_mu         scalar   auxiliary-dir gain in Omega_mu        [rad/m]
%     sim.p_omega      scalar   weight on the angular rate-error term in 
%                               the Lyapunov function V               [N s^2/rad^2]
%     sim.k_omega      scalar   rate feedback gain in the u_2 law     [1/s]
%
%   -- MRP switching
%     sim.epsilon      scalar   required drop in V for the anti-unwinding
%                               jump to fire                          [N]
%
%   =====================================================================
%   OUTER-LOOP DESIGN
%   =====================================================================
%   -- Selected law
%     sim.OLcontrolMode  string  "P", "PI", or "A". "P" is the state
%                                feedback law used in the paper.
%
%   -- Plant matrices used by the outer loop, 6x6 / 6x3 / 3x6
%     sim.A_OL, sim.B_OL, sim.C_OL      the outer-loop model
%     sim.invCAB_OL     3x3      (C_OL*A_OL*B_OL)^-1, which equals m*I.
%                                Supplies the mass scaling in a*, so a*
%                                comes out as a force with no explicit 1/m.
%
%   -- Gains
%     sim.K_P_OL, sim.K_D_OL   3x3   k_P*I and k_D*I
%     sim.Gamma_OL      3x6     integral gain; zeros in "P" mode
%     sim.gamma_nu_OL, sim.gamma_sigma_OL   scalar   adaptive terms;
%                                 zero outside "A" mode
%
%   -- Design intermediates, kept because they record how the gains were
%      computed. All are reconstructible from the lines above and are not
%      read anywhere.
%     sim.A_OL_scalar, sim.B_OL_scalar, sim.C_OL_scalar
%     sim.A_K_scalar, sim.B_K_scalar, sim.Q_K_scalar, sim.R_K_scalar
%     sim.K_OL_scalar   1x2     [k_P/m, k_D/m] from the scalar LQR solve
%
%   =====================================================================
%   REFERENCE MODEL DESIGN
%   =====================================================================
%   The reference model is a cascade of two second-order low-pass filters
%   per axis for the position channel, plus one second-order filter for the
%   pointing channel. Both were designed as transfer functions and converted
%   to state space, then replicated across the three spatial axes.
%
%   -- Filter parameters
%     sim.wx1, sim.wx2       scalar   break frequencies of the two
%                                     position filters               [rad/s]
%     sim.zx1, sim.zx2       scalar   their damping ratios           [-]
%     sim.wnu, sim.znu       scalar   frequency and damping of the
%                                     pointing-vector filter         [rad/s]
%
%   -- Transfer-function objects (design only, not simulation)
%     sim.hx1, sim.hx2, sim.hx   1x1 tf   the two position filters and
%                                         their cascade; hx is fourth order
%     sim.hnu                    1x1 tf   the pointing filter, second order
%
%   -- State-space matrices, position channel (4 states per axis)
%     sim.Ax0, sim.Bx0, sim.Cx0, sim.Dx0   from tf2ss
%     sim.Ax, sim.Bx, sim.Cx, sim.Dx       after the similarity transform
%     sim.Cx1, sim.Dx1                     4x4 and 4x1, used in that
%                                          transform; Dx1 is unused
%
%   -- State-space matrices, pointing channel (2 states per axis)
%     sim.Anu0, sim.Bnu0, sim.Cnu0, sim.Dnu0  from tf2ss
%     sim.Anu, sim.Bnu, sim.Cnu, sim.Dnu      after the transform
%     sim.Cnu1, sim.Dnu1                       used in the transform;
%                                              Dnu1 is unused
%
%   -- Assembled single-axis and three-axis models
%     sim.A, sim.B, sim.C, sim.D          single axis, 6 states in,
%                                         7 outputs out
%     sim.A_ref 18x18, sim.B_ref 18x6     replicated with kron(A,eye(3)), etc.
%     sim.C_ref 21x18, sim.D_ref 21x6
%
%   The similarity transform on the position channel is what makes the first
%   twelve reference states equal the position and its first three time
%   derivatives. The derivative rows of C_ref index that ordering, so the
%   two are not independent choices.
%
%   =====================================================================
%   CONTROL ALLOCATION
%   =====================================================================
%     sim.mixing_matrix  4x4     M, mapping motor forces to the net wrench
%                                [thrust; roll; pitch; yaw]
%     sim.ell            scalar  arm length used in M                   [m]
%     sim.kt             scalar  yaw moment coefficient used in M       [-]
%     sim.arm_length, sim.yaw_moment_coefficient
%                        scalar  the two values ell and kt were built from
%     sim.max_individual_motor_thrust  scalar  per-motor bound          [N]
%     sim.apply_QP      logical  false. enables control allocation if true.
%     sim.saturate_controls  
%                       logical  false. enables control saturation if true.
%     sim.qp_version     scalar  1 through 7; selects a QP allocator when
%                                apply_QP is true. Inactive as shipped.
%     sim.qp_slack, sim.lambda_qp   QP parameters; inactive as shipped
%     sim.qp_options     1x1 optim.options.Quadprog   QP options object;
%                                inactive as shipped, but its construction
%                                is what makes the Optimization Toolbox a
%                                nominal dependency. 
%
%   =====================================================================
%   USER COMMAND SETTINGS
%   =====================================================================
%     sim.override_reference  logical  true. The position channel of the
%                                reference output is replaced by the
%                                analytic history of
%                                ALTERNATE_REFERENCE_FUNCTION rather than
%                                filtered; the pointing channel still runs
%                                through its filter.
%     sim.user_evaluates_error logical  false. Enables a return-to-user.
%     sim.user_position_error_threshold scalar.  Return threshold.      [m]
%
%   =====================================================================
%   SIMULATION SETTINGS
%   =====================================================================
%   -- Horizon and step
%     sim.initial_time, sim.final_time     scalar   t0 and t1          [s]
%     sim.fundamental_time_increment       scalar   the base step dT   [s]
%     sim.t0, sim.t1, sim.dT               scalar   the same three values
%                                          again, under the names
%                                          HYBRID_RK is called with
%
%   -- Initial conditions, also stored twice
%     sim.initial_continuous_system_state  13x1   x(0), layout as in
%                                                 MULTIROTOR_FLOW_MAP
%     sim.initial_discrete_system_state    32x1   xd(0), layout as in
%                                                 CONTROL_LAW
%     sim.xc0, sim.xd0                      the same two vectors again
%
%   -- Integrator controls
%     sim.j0             scalar   initial jump index
%     sim.Jmax           scalar   jump index cap; the run never approaches it
%     sim.event_tol      scalar   event-location tolerance in t         [s]
%     sim.method         string   "rk4"
%     sim.display_flag   logical  true
%     sim.lyap           1x1 struct   Lyapunov step-size control settings;
%                                     fields tol, V, k, relative, V_floor,
%                                     min_step
%
%   -- Function handles passed to HYBRID_RK
%     sim.flow_map, sim.jump_map, sim.jump_guard, sim.manifold_map
%     sim.lyapunov_function
%
%   -- Output control
%     sim.save_figures   logical  the argument RUN_ALL was called with
%
%   =====================================================================
%   SIMULATION OUTPUT
%   =====================================================================
%   -- What HYBRID_RK returns
%     sim.tj_history    Kx2     hybrid time [t j] of each recorded state
%     sim.xc_history    Kx13    continuous state at each recorded state
%     sim.xd_history    Kx32    discrete controller state at each recorded
%                               state
%     sim.return_condition  scalar  1 through 4; 4 means the flow reached
%                               t = tf with the jump index below Jmax,
%                               which is the normal termination
%     sim.stats         1x1 struct   work counters and diagnostics from the
%                               integrator: accepted_steps, jumps,
%                               lyap_rejections, defect_integral,
%                               min_step, and others
%     sim.K             scalar  number of recorded states, size of every
%                               history below
%
%   -- What the diagnostic pass through the flow map adds
%      RUN_ALL re-evaluates MULTIROTOR_FLOW_MAP at every recorded state and
%      harvests the `outputs` struct each stage of the control law writes
%      into. The result is one row per recorded state for each diagnostic.
%
%     Lyapunov function
%       sim.lyapunov_value   Kx1   V(t, eta, omega), the quantity the
%                                  step-size control monitors
%       sim.V                Kx1   the same values again, under the name
%                                  the figures use; identical to
%                                  sim.lyapunov_value
%       sim.V_a              Kx1   the actuation-error part, computed as
%                                  Z(||a*||(||eta|| - n*); 0, Delta_a)
%       sim.ddt_V            Kx1   the certificate bound on dV/dt
%       sim.ddt_V_a          Kx1   the actuation part of ddt_V
%
%     Attitude error
%       sim.eta              Kx3   MRP error vector
%       sim.ddt_eta          Kx3   its time derivative
%       sim.n_star           Kx1   bound on ||eta||, tan(theta*/4); inf
%                                  on the 2*pi plateau
%       sim.ddt_n_star       Kx1   its time derivative
%       sim.ddt2_n_star      Kx1   its second time derivative
%       sim.R_tilde_col      Kx9   R_tilde = R*'R, stored column-major
%       sim.ddt_R_tilde_col  Kx9   its time derivative, same packing
%       sim.switching_variable  Kx1  always false; the switching logic
%                                  lives in MULTIROTOR_JUMP_GUARD
%
%     Target angular velocity omega*
%       sim.omega_star             Kx3   the five-term target rate
%       sim.omega_star_B           Kx3   the same, in the body frame
%       sim.ddt_omega_star         Kx3   its time derivative
%       sim.ddt_omega_R_star       Kx3   derivative of the target-frame
%                                        angular velocity omega_R*
%       sim.omega_star_term_1      Kx3   dissipation, through Omega_a
%       sim.omega_star_term_2      Kx3   Omega_mu feedforward
%       sim.omega_star_term_3      Kx3   ddt_a* feedforward
%       sim.omega_star_term_4      Kx3   ddt_n* containment
%       sim.omega_star_term_5      Kx3   k_eta trim
%       sim.ddt_omega_star_term_1  Kx3   the five derivatives of the above
%       sim.ddt_omega_star_term_2  Kx3
%       sim.ddt_omega_star_term_3  Kx3
%       sim.ddt_omega_star_term_4  Kx3
%       sim.ddt_omega_star_term_5  Kx3
%
%       Terms 1 through 4 are gated by zeta(||eta||/n*; varrho, 1), which
%       is zero wherever ||eta||/n* <= varrho; term 5 is gated separately
%       by zeta(||a*||; a0, a1). The gating is why omega* stays continuous
%       through the free-fall singularity.
%
%     Dissipation
%       sim.Omega_a          Kx1   the dissipation term in omega*; zero
%                                  when V_a is zero or when the ratio is
%                                  not finite
%       sim.ddt_Omega_a      Kx1   its time derivative
%
%     Angular offset
%       sim.theta            Kx1   angle between lambda and lambda*
%       sim.theta_star       Kx1   the smooth bound on it, theta*
%       sim.ddt_theta        Kx1   time derivative of theta
%       sim.ddt_theta_star   Kx1   time derivative of theta*
%       sim.ddt2_theta_star  Kx1   second time derivative of theta*
%
%     Thrust law and attitude directions
%       sim.a_star           Kx3   target actuation a* = kappa_x*
%       sim.ddt_a_star       Kx3   its time derivative
%       sim.ddt2_a_star      Kx3   its second time derivative
%       sim.a_               Kx1   ||a*||, the norm of the target actuation
%       sim.ddt_a_           Kx1   its time derivative
%       sim.ddt2_a_          Kx1   its second time derivative
%       sim.lambda           Kx3   actual thrust direction, third column
%                                  of R_IC
%       sim.lambda_star      Kx3   target direction a*/||a*||
%       sim.ddt_lambda       Kx3   time derivative of lambda
%       sim.ddt_lambda_star  Kx3   time derivative of lambda_star
%       sim.saturation_ratio_term      Kx1  zeta(||eta||/n*; varrho, 1)
%       sim.ddt_saturation_ratio_term  Kx1  its time derivative
%       sim.cross_term       Kx3   the eta-proportional term added to
%                                  u_2 when ||eta|| > n*
%
%     Auxiliary direction
%       sim.mu_star          Kx3   the auxiliary direction
%       sim.ddt_mu_star      Kx3   its time derivative
%       sim.mu_star_tpdt     Kx3   its value at the end of the step
%
%     Controls
%       sim.u_baseline       Kx4   [u1; u2], thrust then body torque
%       sim.f                Kx4   individual motor forces
%       sim.f_baseline       Kx4   the inverse allocation of u_baseline
%       sim.u                Kx4   the net wrench actually realized
%       sim.ddt_omega_B      Kx3   body angular acceleration implied by u2
%       sim.initialized      Kx1   always 1; a flag set by the control law
%       sim.target_thrust    Kx3   a copy of a* written by the outer loop
%
%     Reference and outer-loop state
%       sim.reference_state       Kx18   reference state at each step
%       sim.reference_output      Kx21   reference output at each step
%       sim.user_input            Kx6    user command at each step
%       sim.x_I_ref               Kx3    reference position
%       sim.ddt_x_I_ref           Kx3    its first derivative
%       sim.ddt2_x_I_ref          Kx3    second derivative
%       sim.ddt3_x_I_ref          Kx3    third derivative
%       sim.ddt4_x_I_ref          Kx3    fourth derivative
%       sim.alpha_hat_OL          Kx3    outer-loop compensation state
%       sim.OL_controller_state_tpdt  Kx3  its value at the end of the step
%
%   -- Coordinate histories assembled for the figures
%      These are derived in the extraction block at the end of the
%      postprocessing section, from the quantities above.
%
%     sim.t                Kx1   time at each recorded state            [s]
%     sim.j                Kx1   jump index at each recorded state
%     sim.states           Kx13  continuous state, same as xc_history
%     sim.system_state     Kx45  [xc_history xd_history]
%     sim.mu_star_state    Kx3   columns 14:16 of system_state
%     sim.controls         Kx4   same as sim.u
%     sim.reference_states Kx18  same as sim.reference_state
%     sim.user_input_history  Kx6  same as sim.user_input
%     sim.target_controls  Kx3   same as a*, under a third name
%     sim.u1               Kx1   collective thrust, column 1 of control [N]
%     sim.u2               Kx3   torque, columns 2:4 of controls        [Nm]
%
%     Position and velocity, NED and ENU
%       sim.positions                 Kx3  actual position, NED         [m]
%       sim.reference_positions       Kx3  reference position, NED      [m]
%       sim.velocities                Kx3  actual velocity, NED         [m/s]
%       sim.reference_velocities      Kx3  reference velocity, NED      [m/s]
%       sim.ENU_positions             Kx3  actual, reordered to East-North-Up
%       sim.ENU_reference_positions   Kx3  reference, same reordering
%       sim.ENU_velocities            Kx3  actual velocity, same reordering
%       sim.ENU_reference_velocities  Kx3  reference velocity, same reordering
%       sim.user_positions            Kx3  reference position from the user
%                                          command, NED
%       sim.ENU_user_positions        Kx3  the same, reordered
%
%       The reordering is [N, E, -D] -> [E, N, U]: columns 1 and 2 swap and
%       column 3 changes sign. Because override_reference is true,
%       user_positions and reference_positions hold the same values.
%
%     Attitude and rotation
%       sim.attitudes            Kx4  attitude quaternions, scalar first
%       sim.angular_velocities   Kx3  body angular velocity             [rad/s]
%
%     Thrust vectors and the uprightness trace
%       sim.lambda      Kx3  actual thrust direction, direction only
%       sim.mu          Kx3  the vehicle's forward direction
%       sim.Thrust      Kx3  lambda(t)*u1(t), the actuation vector, NED [N]
%       sim.ENU_Thrust  Kx3  the same, reordered
%       sim.Thrust_star Kx3  the reference actuation, NED               [N]
%       sim.ENU_Thrust_star  Kx3  the same, reordered
%       sim.uprightness Kx1  (-e_3)'*lambda(t); +1 upright, -1 inverted
%
%     Auxiliary-direction alignment
%       sim.straightness              Kx1  alignment of mu with the
%                                          ground-track direction of motion
%       sim.straightness_star         Kx1  alignment of mu* with the same
%       sim.ip_of_mu_and_mu_star      Kx1  mu'*mu*
%
%       straightness and straightness_star are set to NaN wherever the
%       horizontal speed falls below 1e-3 m/s, so the ratio is undefined.
%
%     Reference derivatives
%       sim.ddt_reference_states  Kx18  time derivative of the reference
%                                      state, from A_ref and B_ref
%       sim.reference_outputs     Kx21  reference output recomputed from
%                                      C_ref and D_ref; in override mode
%                                      the first three columns repeat
%                                      user_input_history(:,1:3)
%
%
%% Section 0: Initialization
root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root,'src')));
addpath(genpath(fullfile(root,'scripts')));
addpath(genpath(fullfile(root,'third_party')));
if nargin < 2
    OLcontrolMode = "P"; 
end
if nargin < 1
    save_figures = false;
end

%% Section 1: Initial Conditions
initial_time = 0;
final_time = 25; % Total simulation time
fundamental_time_increment = 0.01;%0.005;%0.0005;%
initial_continuous_system_state = [
    [0;0;0] % Position [m]
    [0;0;0] % Velocity [m/s]
    [1;0;0;0] % Quaternion attitude
    [0;0;0] % Angular velocity [rad/s]
];
initial_discrete_system_state = [
    [1;0;0] % mu_star
    [1;0;0;0] % qeta_lp
    [1;1;1;1] % s_lp
    [0;0;0] % alpha_hat_OL
    zeros(12,1) % Position reference state
    zeros(6,1) % Attitude reference state
];

%% Section 2: Parameters
% Vehicle parameters are based on AirSim's stock quadcopter model.
m = 2.09;
In = diag([0.0297 0.0320 0.0488]);
g = 9.81;
s = tf('s');
wx1 = 7;
wx2 = 8;
zx1 = 1.01;
zx2 = 0.99;
hx1 = 1/(s^2/wx1^2 + 2*zx1*s/wx1 + 1);
hx2 = 1/(s^2/wx2^2 + 2*zx2*s/wx2 + 1);
hx = hx1*hx2;
[Ax0,Bx0,Cx0,Dx0] = tf2ss(hx.Numerator{1},hx.Denominator{1});
Dx1 = [Dx0; Cx0*Bx0; Cx0*Ax0*Bx0; Cx0*Ax0^2*Bx0]; 
Cx1 = [Cx0; Cx0*Ax0; Cx0*Ax0^2; Cx0*Ax0^3];
Ax = (Cx1*(Ax0/Cx1)); 
Ax = [zeros(3,1) eye(3); Ax(4,:)];
Bx = Cx1*Bx0;
Bx = [zeros(3,1); Bx(4)];
Cx = [eye(4); Ax(4,:)];
Dx = [zeros(4,1); Bx(4)];
wnu = 7;
znu = 1.1;
hnu = 1/(s^2/wnu^2 + 2*znu*s/wnu + 1);
[Anu0,Bnu0,Cnu0,Dnu0] = tf2ss(hnu.Numerator{1},hnu.Denominator{1});
Dnu1 = zeros(2,1);
Cnu1 = [Cnu0; Cnu0*Anu0];
Anu = (Cnu1*(Anu0/Cnu1));
Anu = [0 1; Anu(2,:)];
Bnu = Cnu1*Bnu0;
Bnu = [0; Bnu(2)];
Cnu = eye(2);
Dnu = zeros(2,1);
A = [Ax zeros(4,2); zeros(2,4) Anu];
B = [Bx zeros(4,1); zeros(2,1) Bnu];
C = [Cx zeros(5,2); zeros(2,4) Cnu];
D = [Dx zeros(5,1); zeros(2,1) Dnu];
A_ref = kron(A,eye(3));
B_ref = kron(B,eye(3));
C_ref = kron(C,eye(3));
D_ref = kron(D,eye(3));
switch OLcontrolMode
    case "PI"
    % PI control
    A_OL_scalar = [0  1; 0 0]; 
    B_OL_scalar = [0; 1/m];
    C_OL_scalar = [1  0];
    A_OL = kron(A_OL_scalar,eye(3));
    B_OL = kron(B_OL_scalar,eye(3));
    C_OL = kron(C_OL_scalar,eye(3));
    invCAB_OL = inv(C_OL*A_OL*B_OL);
    A_K_scalar = [0 1 0; 0 0 1; 0 0 0]; 
    B_K_scalar = [0; 0; 1];
    Q_K_scalar = diag([0.2 2 2]);
    R_K_scalar = 1;
    [K_OL_scalar,~] = lqr(A_K_scalar,B_K_scalar,Q_K_scalar,R_K_scalar);
    K_I_OL = K_OL_scalar(1)*eye(3);
    K_P_OL = K_OL_scalar(2)*eye(3);
    K_D_OL = K_OL_scalar(3)*eye(3);
    Gamma_OL = [K_I_OL zeros(3,3)];
    gamma_nu_OL = 0;
    gamma_sigma_OL = 0;
    case "A"
    % Adaptive control
    A_OL_scalar = [0  1; 0 0]; 
    B_OL_scalar = [0; 1/m];
    C_OL_scalar = [1  0];
    A_OL = kron(A_OL_scalar,eye(3));
    B_OL = kron(B_OL_scalar,eye(3));
    C_OL = kron(C_OL_scalar,eye(3));
    invCAB_OL = inv(C_OL*A_OL*B_OL);
    A_K_scalar = [0 1; 0 0]; 
    B_K_scalar = [0; 1];
    Q_K_scalar = diag([2 2]);
    R_K_scalar = 1;
    [K_OL_scalar,P_OL_scalar] = lqr(A_K_scalar,B_K_scalar,Q_K_scalar,R_K_scalar);
    K_P_OL = K_OL_scalar(1)*eye(3);
    K_D_OL = K_OL_scalar(2)*eye(3);
    R_OL_scalar = P_OL_scalar\(Q_K_scalar/P_OL_scalar)/2 - B_K_scalar*(R_K_scalar\B_K_scalar')/2;
    [~,R_is_positive_definite] = chol(R_OL_scalar);
    if ~R_is_positive_definite; error('main: R was not positive definite!'); end
    gamma_OL = 10;
    Gamma_OL = gamma_OL*[P_OL_scalar(1,2)*eye(3) P_OL_scalar(2,2)*eye(3)];
    nu_OL = 10;
    sigma_OL = 0;
    gamma_nu_OL = gamma_OL*nu_OL;
    gamma_sigma_OL = gamma_OL*sigma_OL;
    otherwise % Default to "P"
    % P control
    A_OL_scalar = [0  1; 0 0]; 
    B_OL_scalar = [0; 1/m];
    C_OL_scalar = [1  0];
    A_OL = kron(A_OL_scalar,eye(3));
    B_OL = kron(B_OL_scalar,eye(3));
    C_OL = kron(C_OL_scalar,eye(3));
    invCAB_OL = inv(C_OL*A_OL*B_OL);
    A_K_scalar = [0 1; 0 0]; 
    B_K_scalar = [0; 1];
    Q_K_scalar = diag([2 2]);
    R_K_scalar = 1;
    K_OL_scalar = lqr(A_K_scalar,B_K_scalar,Q_K_scalar,R_K_scalar);
    K_P_OL = K_OL_scalar(1)*eye(3);
    K_D_OL = K_OL_scalar(2)*eye(3);
    Gamma_OL = zeros(3,6);
    gamma_nu_OL = 0;
    gamma_sigma_OL = 0;
end
r = 5;
U1min = 0.55*r;%0;%
rho = U1min/r;
base = (1-rho)/2;
alpha_2pi = base/4;
delta0 = 2*pi/24;
delta1 = 2*delta0;
varrho = 0.2;
Delta_a = 2;
Delta_ddt_n_star = 2;
epsilon = 10;
k_a = 2;
k_eta = 10;
a0 = alpha_2pi*r/5;
a1 = r/5;
k_mu = 2;
p_omega = 1;
k_omega = 40;
arm_length = 0.255;
yaw_moment_coefficient = 0.1;
ell = arm_length;
kt = yaw_moment_coefficient;
mixing_matrix = [
    [  1    1    1    1 ]
    [-ell  ell  ell -ell]/sqrt(2)
    [-ell -ell  ell  ell]/sqrt(2)
    [ kt  -kt   kt  -kt ]
];
max_individual_motor_thrust = 12;
qp_slack = 1e-3;
lambda_qp = [0.01; 0.01];
qp_options = optimoptions('quadprog','Display','off','Algorithm','active-set');
qp_version = 7;
apply_QP = false;
saturate_controls = false;
user_evaluates_error = false;
user_position_error_threshold = 5;
override_reference = true;

parameters.initial_time = initial_time;
parameters.final_time = final_time;
parameters.sampling_time = fundamental_time_increment;
parameters.m = m;
parameters.In = In;
parameters.g = g;
parameters.A_ref = A_ref;
parameters.B_ref = B_ref;
parameters.C_ref = C_ref;
parameters.D_ref = D_ref;
parameters.A_OL = A_OL;
parameters.B_OL = B_OL;
parameters.C_OL = C_OL;
parameters.invCAB_OL = invCAB_OL;
parameters.K_P_OL = K_P_OL;
parameters.K_D_OL = K_D_OL;
parameters.Gamma_OL = Gamma_OL;
parameters.gamma_nu_OL = gamma_nu_OL;
parameters.gamma_sigma_OL = gamma_sigma_OL;
parameters.U1min = U1min;
parameters.r = r;
parameters.rho = rho;
parameters.delta0 = delta0;
parameters.delta1 = delta1;
parameters.varrho = varrho;
parameters.Delta_a = Delta_a;
parameters.Delta_ddt_n_star = Delta_ddt_n_star;
parameters.epsilon = epsilon;
parameters.k_a = k_a;
parameters.k_eta = k_eta;
parameters.a0 = a0;
parameters.a1 = a1;
parameters.k_mu = k_mu;
parameters.k_omega = k_omega;
parameters.p_omega = p_omega;
parameters.mixing_matrix = mixing_matrix;
parameters.max_individual_motor_thrust = max_individual_motor_thrust;
parameters.apply_QP = apply_QP;
parameters.qp_slack = qp_slack;
parameters.lambda_qp = lambda_qp;
parameters.qp_options = qp_options;
parameters.saturate_controls = saturate_controls;
parameters.qp_version = qp_version;
parameters.user_evaluates_error = user_evaluates_error;
parameters.user_position_error_threshold = user_position_error_threshold;
parameters.override_reference = override_reference;

%% Section 3: System Simulation with Continuous Control

% Prepare for integration...
flow_map = @multirotor_flow_map;
jump_map = @multirotor_jump_map;
jump_guard = @multirotor_jump_guard;
manifold_map = @multirotor_manifold_map;
t0 = initial_time;
t1 = final_time;
dT = fundamental_time_increment;
j0 = 0;
Jmax = 1e5;
event_tol = 1e-12;
xc0 = initial_continuous_system_state;
xd0 = initial_discrete_system_state;
display_flag = true;
method = "rk4";
lyapunov_function = @multirotor_lyapunov_function;
k = min(k_a,k_omega);
lyap = struct('tol', 0.1*k, 'V', lyapunov_function, 'k', k, 'relative', true, 'V_floor', 1e-6, 'min_step', 1e-8*dT);

% Integrate simulation...
[tj_history,xc_history,xd_history,return_condition,stats] = ...
    hybrid_RK(flow_map,jump_map,jump_guard,manifold_map,...
    t0,t1,dT,j0,Jmax,event_tol,xc0,xd0,parameters,display_flag,method,lyap);

%% Section 4: Postprocess Data

% Extract output structures...
K = size(tj_history,1);
dK = 1e3;
outc = cell(K,1);
for k = 1:K
    if k < K
        dtk = tj_history(k+1,1)-tj_history(k,1);
    else
        dtk = dT;
    end
    [~,~,outk] = multirotor_flow_map(tj_history(k,:),dtk,xc_history(k,:)',xd_history(k,:)',parameters);
    outc{k,1} = outk;
    if (k == 1 || mod(k,dK) == 0 || k == K); disp("Postprocessing Progress: loop 1/4, k = "+k+"/"+K); end
end
fields = fieldnames(outk);
clear k dtk dK
outputs = struct();
for jj = 1:length(fields)
    outputs.(fields{jj}) = nan(K,numel(outk.(fields{jj})));
    disp("Postprocessing Progress: loop 2/4, jj1 = "+jj+"/"+length(fields))
end
clear outk
for jj = 1:length(fields)
    for ii = 1:K
        outputs.(fields{jj})(ii,:) = reshape(outc{ii,1}.(fields{jj}),1,[]);
    end
    disp("Postprocessing Progress: loop 3/4, jj2 = "+jj+"/"+length(fields))
end
clear jj outc
fields = fieldnames(outputs);
for ii = 1:numel(fields)
    eval([fields{ii} ' = outputs.(fields{ii});'])
    disp("Postprocessing Progress: loop 4/4, ii = "+ii+"/"+numel(fields))
end
clear ii fields

% Extract output histories...
t = tj_history(:,1);
j = tj_history(:,2);
system_state = [xc_history xd_history];
states = system_state(:,1:13);
mu_star_state = system_state(:,14:16);
controls = u;
reference_states = reference_state;
user_input_history = user_input;
target_controls = a_star;
V = lyapunov_value;
positions = states(:,1:3);
ENU_positions = [positions(:,2) positions(:,1) -positions(:,3)];
velocities = states(:,4:6);
ENU_velocities = [velocities(:,2) velocities(:,1) -velocities(:,3)];
attitudes = states(:,7:10);
angular_velocities = states(:,11:13);
ddt_reference_states = reference_states*A_ref' + user_input_history*B_ref';
reference_outputs = reference_states*C_ref' + user_input_history*D_ref';
reference_positions = reference_outputs(:,1:3);
ENU_reference_positions = [reference_positions(:,2) reference_positions(:,1) -reference_positions(:,3)];
user_positions = user_input_history(:,1:3);
ENU_user_positions = [user_positions(:,2) user_positions(:,1) -user_positions(:,3)];
reference_velocities = reference_outputs(:,4:6);
ENU_reference_velocities = [reference_velocities(:,2) reference_velocities(:,1) -reference_velocities(:,3)];
u1 = controls(:,1);
u2 = controls(:,2:4);
uprightness = nan(K,1);
straightness = nan(K,1);
straightness_star = nan(K,1);
ip_of_mu_and_mu_star = nan(K,1);
lambda = nan(K,3);
mu = nan(K,3);
Thrust = nan(K,3);
V_a = nan(K,1);
for kk = 1:K
    lambda(kk,:) = (attitude_quaternion_to_matrix(states(kk,7:10)'/norm(states(kk,7:10)))*(-unit_vector(3,3)))';
    mu(kk,:) = (attitude_quaternion_to_matrix(states(kk,7:10)'/norm(states(kk,7:10)))*(unit_vector(1,3)))';
    Thrust(kk,:) = u1(kk)*lambda(kk,:);
    uprightness(kk) = lambda(kk,:)*(-unit_vector(3,3));
    temp = norm(user_input_history(kk,4:6));
    temp(temp<0.001) = nan;
    straightness(kk) = user_input_history(kk,4:6)*(mu(kk,:))'/temp;
    straightness_star(kk) = user_input_history(kk,4:6)*(mu_star(kk,:))'/temp;
    ip_of_mu_and_mu_star(kk) = mu_star(kk,:)*(mu(kk,:)');
    V_a(kk) = intssat(norm(a_star(kk,:))*(norm(eta(kk,:))-n_star(kk)),0,Delta_a);
end
clear kk temp
ENU_Thrust = [Thrust(:,2) Thrust(:,1) -Thrust(:,3)];
Thrust_star = target_controls;
tau = vecnorm(Thrust_star')';
lambda_star = Thrust_star./vecnorm(Thrust_star')';
ENU_Thrust_star = [Thrust_star(:,2) Thrust_star(:,1) -Thrust_star(:,3)];

clear outputs initialized switching_variable s root
list = who;
sim = struct();
for ii = 1:numel(list)
    eval(['sim.' list{ii} ' = ' list{ii} ';'])
end

%% Section 5: Plot Results

% Timing history
fig = figure(Theme='light');
yyaxis left
plot(t)
xlabel('Index, $k$','Interpreter','latex')
ylabel('Time, $t$','Interpreter','latex')
yyaxis right
plot(j)
ylabel('Jump index, $j$','Interpreter','latex')
title('Timing History')
xlim([1 K])
fig.Position = [2 50 718 398];
if save_figures; exportgraphics(fig,fullfile('figures','diagnostic_figure_1.pdf'), 'ContentType','vector'); end

% MRP history
fig = figure(Theme='light');
plot(t,eta,'.','DisplayName','$\eta(t)$','MarkerSize',3)
hold on
plot(t,vecnorm(eta')','DisplayName','$\|\eta(t)\|$')
yL = ylim;
plot(t,n_star,'DisplayName','$n^{\star}(t)$')
hold off
ylim(yL)
xlabel('Time [sec]','Interpreter','latex')
ylabel('Error MRP','Interpreter','latex')
grid on
legend('Location','best','Interpreter','latex')
title('MRP History')
fig.Position = [2 50 718 398];
if save_figures; exportgraphics(fig,fullfile('figures','diagnostic_figure_2.pdf'), 'ContentType','vector'); end

% Straightness history
fig = figure(Theme='light');
plot(t,straightness,'DisplayName','$((I - e_3 e_3^T){\rm v}(t))^T\mu(t)/\|{\rm v}(t)\|$')
hold on
plot(t,straightness_star,'DisplayName','$((I - e_3 e_3^T){\rm v}(t))^T\mu^{\star}(t)/\|{\rm v}(t)\|$')
plot(t,ip_of_mu_and_mu_star,'DisplayName','$\mu(t)^T\mu^{\star}(t)$')
xlabel('Time [sec]','Interpreter','latex')
ylabel('Cosine','Interpreter','latex')
grid on
legend('Location','southwest','Interpreter','latex')
ylim([-1.2 1.2])
title('Alignment of FPV camera direction')
fig.Position = [2 50 718 398];
if save_figures; exportgraphics(fig,fullfile('figures','diagnostic_figure_3.pdf'), 'ContentType','vector'); end

% Lyapunov function history
fig = figure(Theme='light');
semilogy(t,V,'DisplayName','$V(t)$')
hold on
semilogy(t,V_a,'DisplayName','$V_{\rm a}(t)$')
hold off
xlabel('Time, $t$','Interpreter','latex')
ylabel('Lyapunov value','Interpreter','latex')
grid on
legend('Location','best','Interpreter','latex')
title('Lyapunov Function History')
fig.Position = [2 50 718 398];
if save_figures; exportgraphics(fig,fullfile('figures','diagnostic_figure_4.pdf'), 'ContentType','vector'); end

% Time step history
fig = figure(Theme='light');
semilogy(t(1:end-1),diff(t))
xlabel('Time, $t$','Interpreter','latex')
ylabel('Time Step, $dt$','Interpreter','latex');
grid on;
legend('Location','best','Interpreter','latex');
title('Simulation Step Size History')
fig.Position = [2 50 718 398];
if save_figures; exportgraphics(fig,fullfile('figures','diagnostic_figure_5.pdf'), 'ContentType','vector'); end

% Figures included in the paper
figure_2(save_figures)
figure_3(save_figures)
figure_3_alt(save_figures)
figure_4(save_figures)
figure_5(sim,save_figures)
figure_6_through_9(sim,save_figures)
figure_10(sim,save_figures)

end