function [reference_state_tpdt,reference_output,outputs] = reference_model(t,dt,reference_state,user_input,parameters,outputs)
%REFERENCE_MODEL  Linear reference model, stepped with an exact ZOH discretization.
%   [xd_ref_tpdt,ref_out,outputs] = REFERENCE_MODEL(t,dt,xd_ref,user_input, ...
%       parameters,outputs)
%   advances the reference state by one flow step of length dt and returns
%   the reference outputs evaluated at the current time.
%
%   The model is the linear system (A_ref, B_ref, C_ref, D_ref), built in
%   RUN_ALL from cascaded second-order filters: a fourth-order low-pass per
%   component for the position reference, and the second-order filter of
%   eq. (outer_claw_2) for the pointing vector nu. The step is taken with
%   the exact zero-order-hold discretization
%       [Phi Gamma; 0 I] = expm([A_ref B_ref; 0 0]*dt)
%   so the propagation is exact for a user input that is constant over the 
%   step, with no truncation error, and no accuracy penalty for the large 
%   dT used in the simulation. 
%
%   Time convention: the reference outputs are evaluated at the current
%   time t from reference_state, while the returned state has advanced to
%   t + dt. The two are therefore offset by one step by design.  
%
%   Reference output layout (21x1):
%       ref_out(1:15)     position and its first four time derivatives, as
%                         three 3x1 blocks of [x*; dx*; ddx*; dddx*; ddddx*]
%                         in the component ordering returned by
%                         ALTERNATE_REFERENCE_FUNCTION (see below)
%       ref_out(16:18)    pointing vector nu
%       ref_out(19:21)    ddt_nu, the time derivative of nu
%
%   Reference state layout (18x1):
%       xd_ref(1:12)      position reference: four states per component
%       xd_ref(13:18)     nu filter: two states per component
%
%   When parameters.override_reference is true, the position blocks are
%   replaced by the analytic history of ALTERNATE_REFERENCE_FUNCTION: the
%   first fifteen entries of ref_out are overwritten with x_user and its
%   first four derivatives at t, and the first twelve entries of the
%   returned state are overwritten with x_user and its first three
%   derivatives at t + dt. The nu filter states are NOT overridden, so nu
%   remains a filtered version of the direction channel of user_input even
%   in override mode. The component ordering of the position block is
%   [Y; X; -Z], which is a component swap relative to the ordering of the
%   vehicle state; see ALTERNATE_REFERENCE_FUNCTION for details.
%
%   RUN_ALL sets override_reference to true, which matches the paper's
%   statement that the target outer-loop trajectory is not processed
%   through a reference model. With override_reference false, the model is
%   used as designed: a linear reference model whose initial conditions
%   agree with the vehicle state, as discussed in the simulation section.
%
%   Inputs:
%       t, dt              current time and flow step length
%       xd_ref             18x1 reference state, layout above
%       user_input         6x1 user command from USER_INPUT_FUNCTION:
%                          reference position, then the desired pointing
%                          vector; only the direction channel drives the nu
%                          filter
%       parameters         struct of system and controller parameters; uses
%                          parameters.A_ref, parameters.B_ref,
%                          parameters.C_ref, parameters.D_ref, and
%                          parameters.override_reference
%       outputs            diagnostics struct, threaded through
%
%   Outputs:
%       xd_ref_tpdt        18x1 reference state at t + dt
%       ref_out            21x1 reference output at t, layout above
%       outputs            diagnostics struct, extended with
%                          outputs.reference_state and
%                          outputs.reference_output
%
%   Requires on the MATLAB path:
%       alternate_reference_function (only when override_reference is true).
%
%   Called by CONTROL_LAW.

% Load parameter values...
A_ref = parameters.A_ref;
B_ref = parameters.B_ref;
C_ref = parameters.C_ref;
D_ref = parameters.D_ref;

% Exact zero-order-hold discretization over [t, t+dt]
n = size(A_ref,1);
m = size(B_ref,2);
PhiGamma = expm([A_ref, B_ref; zeros(m,n), zeros(m,m)]*dt);
Phi   = PhiGamma(1:n, 1:n);
Gamma = PhiGamma(1:n, n+1:n+m);

% Propagate reference state to t+dt, then evaluate the output there
reference_state_tpdt = Phi*reference_state + Gamma*user_input;
reference_output     = C_ref*reference_state + D_ref*user_input;

% Handle the case when reference is overridden...
if parameters.override_reference
    [x_user,dx_user,ddx_user,dddx_user,ddddx_user,~] = alternate_reference_function(t);
    reference_output(1:15,1) = [x_user;dx_user;ddx_user;dddx_user;ddddx_user];
    [x_user,dx_user,ddx_user,dddx_user,~,~] = alternate_reference_function(t+dt);
    reference_state_tpdt(1:12,1) = [x_user;dx_user;ddx_user;dddx_user];
end
    
% Manage outputs...
outputs.reference_state = reference_state; 
outputs.reference_output = reference_output; 

end