% PCT for GIMBALED UNIDIRECTIONAL ACTUATION
%   Control and simulation code accompanying
%       "Practical Complete Tracking Through Gimbaled Unidirectional
%        Actuation with Application to Multirotor UAVs"
%       I. J. Willebeek-LeMair and C. A. Woolsey
%
%   The single entry point is RUN_ALL, one level up from this folder. It
%   sets the parameters and initial conditions, integrates the closed-loop
%   system with HYBRID_RK, postprocesses the diagnostics into a struct, and
%   calls each figure function in ../scripts. The functions here implement
%   recently developed theory.
%
%   Integration
%     hybrid_RK                       - Hybrid system integrator: Runge-Kutta
%                                       flows, event location, and Lyapunov
%                                       step-size control
%
%   Vehicle model and hybrid-system data
%     multirotor_flow_map             - Continuous-state derivative and discrete-
%                                       state update for one flow step
%     multirotor_jump_guard           - Jump set: evaluates V and the switching
%                                       test (rebuilds the inner-loop law inline)
%     multirotor_jump_map             - Jump map: flips the MRP-hemisphere memory
%     multirotor_lyapunov_function    - Lyapunov value V for the inner loop
%                                       (rebuilds the inner-loop law inline)
%     multirotor_manifold_map         - Projects the state onto the system
%                                       manifold
%
%   Controller cascade
%     control_law_for_simulation      - Complete control law: user command,
%                                       controller cascade, and mixing
%     control_law                     - The cascade itself, as it would be
%                                       implemented on a vehicle
%     reference_model                 - Linear reference model, exact ZOH step
%     outer_loop_controller           - Target actuation a* and the coefficient
%                                       matrices for its time derivatives
%     inner_loop_controller           - Gimbal control law: thrust u1 and
%                                       torque u2
%     control_allocation              - Distributes thrust and torque across
%                                       the individual rotors
%
%   Reference history and user command
%     alternate_reference_function    - Obstacle-course reference position and
%                                       its exact derivatives through 4th order
%     user_input_function             - Builds the user command at time t
%
%   Angular geometry
%     theta_bar                       - Exact admissible half-angle (used by the
%                                       schematic figures)
%     theta_star_function             - Smooth underestimate of theta_bar
%     dtheta_star                     - First derivative of theta_star
%     ddtheta_star                    - Second derivative of theta_star
%
%   Smooth step machinery (Appendix A)
%     ssat                            - Smooth step function zeta
%     dssat                           - First derivative of zeta
%     ddssat                          - Second derivative of zeta
%     dddssat                         - Third derivative of zeta
%     ddddssat                        - Fourth derivative of zeta
%     intssat                         - Integral of zeta: the ramp function Z
%
%   Attitude utilities
%     attitude_matrix_to_quaternion   - Rotation matrix to quaternion, with
%                                       hemisphere continuity
%     attitude_quaternion_kinematic_matrix
%                                     - The matrix E(q) in qdot = E(q)*omega_B
%     attitude_quaternion_to_matrix   - Quaternion to rotation matrix
%     hat_map                         - Cross-product equivalent matrix
%     mrp_kinematic_matrix            - The matrix T(eta) in etadot = T(eta)*omega
%     quat2mrp                        - Quaternion to modified Rodrigues parameter
%     unit_vector                     - The unit vector e_i in R^n
%
%   Figure styling (used by the schematic figures)
%     tacStyle                        - Shared figure style constants
%     tacAxes                         - Applies the shared styling to an axes
%
%   See also RUN_ALL, figure_2, figure_3, figure_4, figure_5,
%   figure_6_through_9, figure_10.