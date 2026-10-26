function [tj_history,xc_history,xd_history,return_condition,stats] = hybrid_RK(flow_map,jump_map,jump_guard,manifold_map,t0,tf,dT,j0,Jmax,event_tol,xc0,xd0,parameters,display_flag,method,lyap)
%HYBRID_RK Simulate a hybrid system using explicit Runge-Kutta flows.
%
%   [tj_history,xc_history,xd_history,return_condition] = ...
%       hybrid_RK(flow_map,jump_map,jump_guard, ...
%       manifold_map,t0,tf,dT,j0,Jmax,event_tol,xc0,xd0,parameters)
%
%   [...] = hybrid_RK(...,parameters,display_flag)
%   also prints progress messages to the Command Window when display_flag
%   is true.
%
%   [...] = hybrid_RK(...,parameters,display_flag,method)
%   selects the Runge-Kutta method used during flows. The default is the
%   classical fourth-order method, "rk4".
%
%   [...] = hybrid_RK(...,display_flag,method,lyap)
%   shortens flow steps where the computed solution violates a Lyapunov
%   certificate Vdot <= -k*V of the control law (see LYAPUNOV STEP-SIZE
%   CONTROL).
%
%   [...,return_condition,stats] = hybrid_RK(...) also
%   returns a struct of work counters and diagnostics (see OUTPUTS).
%
%   Computes an approximate solution of a hybrid system with continuous
%   state xc and discrete state xd on a hybrid time domain indexed by
%   (t,j). The jump set is the sublevel set
%       D = {(t,j,xc,xd) : g(t,j,xc,xd) <= 0}
%   of a scalar guard function g supplied by jump_guard. During flows, xc
%   is propagated with an explicit Runge-Kutta method, and flow_map
%   supplies the value of xd at the end of each flow step, given the step
%   length dt. This lets xd carry discrete-time states, such as
%   zero-order-hold memory or quantities updated over each step, that
%   change once per accepted flow step. When a flow step ends in D, the
%   time at which the flow enters D is located by root-finding on g, the
%   flow is stopped there, and the jump map is applied: the jump index j
%   increments and t is held fixed. After every flow step and every jump,
%   the state is passed through manifold_map to map it back onto the
%   system manifold.
%
%   HYBRID CONVENTIONS
%   - Jumps have priority. The guard is checked first, and the system
%     flows only from states with g > 0. The flow set is therefore the
%     complement of the jump set, so together they cover the entire state
%     space and a solution never terminates by leaving their union.
%   - D is closed: a state with g = 0 is in D. The jump map should
%     therefore map into g > 0 (or to a state that jumps again on purpose);
%     a post-jump state left on g = 0 jumps again immediately.
%   - Consecutive jumps at the same t are allowed. Each jump sequence runs
%     until the state leaves D or j reaches Jmax. This includes jumps at
%     t = tf.
%   - Several guards g_1,...,g_m combine into one: use min(g_1,...,g_m)
%     for the union of their sublevel sets (jump when any guard is
%     reached) and max(g_1,...,g_m) for the intersection (jump only when
%     all are). For example, a bouncing ball with height h and vertical
%     velocity v uses g = max(h,v), so D = {h <= 0 and v <= 0} and the
%     ball leaves D at every bounce.
%
%   RUNGE-KUTTA FLOW STEP
%   A method with s stages and Butcher tableau (A,b,c) takes a trial step
%   of length dt from the state (xc,xd) at hybrid time [t j] as follows.
%   For each stage i = 1,...,s,
%       xc_i = xc + dt*sum_{m<i} A(i,m)*K_m,
%       K_i  = xcdot returned by flow_map([t+c(i)*dt j],dt,xc_i,xd,...),
%   and the trial state is
%       xc_next = xc + dt*sum_i b(i)*K_i,
%       xd_next = discrete-state output of the stage-1 call,
%   after which (xc_next,xd_next) is passed through manifold_map. So xd is
%   held at its value at the start of the step for every stage
%   (zero-order hold within the step), and only the stage-1 call, which is
%   made at the start state, defines the discrete update. The
%   discrete-state outputs of stages 2,...,s are discarded. Stage states
%   are not passed through manifold_map; only the end-of-step state is.
%   This is the standard projection method, and it keeps the order of the
%   underlying Runge-Kutta method when flow_map is smooth in a
%   neighborhood of the manifold.
%
%   STEP-SIZE AND EVENT STRATEGY
%   - Flow steps have length dT unless Lyapunov step-size control
%     shortens them (see LYAPUNOV STEP-SIZE CONTROL), and they are always
%     clamped so t never passes tf. A step clamped by tf ends exactly at
%     t = tf, and a step that would end within 1e-10*dT of tf is
%     lengthened to end there, so roundoff accumulated in t (e.g., ten
%     steps of 0.1 summing to 1 - 1.1e-16) cannot leave a sliver step just
%     before tf.
%   - If a trial step of length dt ends with g > 0, it is accepted. If it
%     ends with g <= 0, then g changes sign over the step, because the
%     step started with g > 0. Consider the trial-step map
%         tau -> g at the end of a trial step of length tau, 0 < tau <= dt,
%     taken from the same start state with the same method. The event
%     search shrinks the bracket [0,dt], on which this map goes from
%     positive to nonpositive, until its width is at most event_tol, using
%     the Illinois variant of regula falsi with a bisection safeguard.
%     Every trial point is a complete Runge-Kutta step (s flow_map calls,
%     one manifold_map call, and one jump_guard call), so the accepted
%     step is exactly a step of the method, with flow_map's discrete
%     update computed for its true length.
%   - The step accepted at the end of the search is the right end of the
%     final bracket. It ends in D (g <= 0), at most event_tol in t after
%     the point where the computed flow enters D. The jump sequence is
%     then applied, and steps of length dT resume from the post-jump time.
%   - The search converges superlinearly when the trial-step map is
%     smooth near its root, typically in a handful of trial steps. The
%     bisection safeguard halves the bracket at least every second trial,
%     so a search never takes more than about 2*log2(dT/event_tol) trials.
%   - If the discrete update xd_next does not approach xd as the step
%     length shrinks, and that update alone places the end of the step in
%     D, the trial-step map is discontinuous at tau = 0. The search then
%     converges to tau = 0, and a step no longer than event_tol is
%     accepted before the jump, at the cost of about 2*log2(dT/event_tol)
%     trial steps.
%   - Jump detection checks the guard only at the ends of steps. A
%     trajectory that enters and leaves D within a single step (g dips
%     below zero and recovers) is not detected, so D should have thickness
%     along the flow. A convergence study in dT is the intended check that
%     no jumps are being missed.
%
%   LYAPUNOV STEP-SIZE CONTROL
%   For a control law certified by a Lyapunov function V with
%   Vdot <= -k*V (k >= 0), the solver can shorten flow steps where the
%   computed solution violates the certificate. The certificate is
%   supplied as the function lyap.V and the rate lyap.k. Because
%       d/dt( exp(k*(t - t_k))*V ) = exp(k*(t - t_k))*(Vdot + k*V),
%   the integral of Vdot + k*V over a step from t_k to t_{k+1} = t_k + dt,
%   weighted by exp(k*(t - t_k)) and taken along the computed solution,
%   is exactly the certificate defect
%       delta = exp(k*dt)*V(t_{k+1}) - V(t_k),
%   with no quadrature error. The defect is computed from values of V,
%   not from the analytical Vdot: wherever the certificate holds as a
%   pointwise inequality, Vdot + k*V <= 0 at every state, accurate or not,
%   so evaluating Vdot + k*V can never reveal integration error. delta > 0
%   means the computed step increased V faster than the certificate
%   allows, whatever the cause: integration error (including the
%   instability of an explicit method on a stiff closed loop), stale
%   discrete states, or roundoff in the control law.
%
%   lyap.V is evaluated once per checked trial step, at the end of the
%   step after manifold_map; V(t_k) is carried over from the end of the
%   previous step (or evaluated at the start of each flow interval, after
%   the initial state and after every jump sequence). flow_map is not
%   called any additional times.

%   A step's violation is max(delta,0). With lyap.relative = true, it is
%   then divided by max(V(t_k),lyap.V_floor), so it measures growth of V
%   relative to its current size (see SCALE OF THE BOUND). A trial step is
%   rejected when its violation exceeds
%       c = lyap.tol*dT.
%   The step length is then multiplied by 0.9*c/violation, limited to the
%   range [0.1,0.5], and the step is retried, but never with a length
%   below lyap.min_step. A step that still violates the bound at
%   lyap.min_step is accepted anyway: a warning is issued the first time
%   this happens in a call, and every such step is counted in
%   stats.lyap_floor_steps. After an accepted step of length dt with
%   violation v > 0, the next step is at most 0.9*c*dt/v; otherwise steps
%   double back toward dT. After a jump, the step resets to dT. Trial
%   steps of an event search are shorter than the step that triggered the
%   search and are not checked.
%
%   The bound c scales with dT so that convergence studies stay clean.
%   When the violation accrues at a rate set by the state (delta about
%   dt*rbar, with rbar the weighted mean of Vdot + k*V along the computed
%   step), the controlled step is about dT*lyap.tol/rbar, which is
%   proportional to dT, so halving dT halves every step, shortened ones
%   included. When the violation comes from integration error (delta about
%   C*dt^(p+1) for a method of order p), the controlled step is about
%   (lyap.tol*dT/C)^(1/(p+1)); it still shrinks as dT is halved, but more
%   slowly, so a study in dT refines those stretches less than the rest of
%   the solution. lyap.tol is the violation rate that a step of length dT
%   may carry before it is shortened.
%
%   SCALE OF THE BOUND
%   By default (lyap.relative = false), c is an absolute amount of V, and
%   lyap.tol has units of V per unit time. Because delta is proportional
%   to V, an absolute bound is strict far from the equilibrium and
%   permissive near it: once V << c, even a step that multiplies V several
%   times over passes. With lyap.relative = true, the violation is divided
%   by max(V(t_k),lyap.V_floor), lyap.tol has units of 1/time like k, and
%   a step passes when, roughly, V grows over the step by at most a factor
%   exp(-k*dt)*(1 + lyap.tol*dT). For example, lyap.tol = 0.1*k lets each
%   step of length dT decay at most 10% slower than certified. Below
%   lyap.V_floor the bound becomes absolute again, which keeps roundoff in
%   V near the equilibrium from forcing tiny steps.
%
%   The control promotes greater accuracy in regions where lower temporal
%   resolution might wrongly suggest the certificate does not hold. In some
%   circumstances, the certificate may not hold exactly for any step size 
%   since numerical inacuracies constitute a disturbance to the control 
%   system (and the certificate Vdot <= -k*V only holds exactly in the 
%   absence of disturbances). Reducing the step size whenever the 
%   certificate is violated reduces the extent of the violation if the 
%   source of the violation is a numerical inaccuracy.  If the source of a 
%   violation were within the dynamical system (i.e., not due to numerical 
%   integration), then reducing the step size would not reduce the extent 
%   of the violation.  Therefore, the Lyapunov step-size control promotes 
%   accuracy locally (in time) by enforcing the certificate and does not 
%   hide any performance that contradicts the claimed certificate. 
%   Globally, accuracy is set by the base step size `dT`.  Global accuracy 
%   remains the job of the convergence study in dT. 
%
%   stats.defect_integral sums the violations of accepted steps (in the
%   units of the bound). If it shrinks as dT is halved, the violation
%   depends on the step length, and the control is removing it. If it
%   stays roughly constant, the violation is a property of the computed
%   states (e.g., roundoff in the control law at stiff points, or a
%   certificate that does not hold there), and shorter steps only spread
%   it over more steps.
%
%   ACCURACY ACROSS JUMPS
%   Between jumps, a method of order p has global error O(dT^p). Each jump
%   time is located to within event_tol of the jump time of the computed
%   solution, which itself differs from the exact jump time by O(dT^p)
%   when the flow crosses g = 0 transversally. The error after a jump is
%   therefore O(dT^p) + O(event_tol), so order p is kept through jumps as
%   long as event_tol is well below the integration error. Because the
%   search converges superlinearly, a tight tolerance (e.g., 1e-12 for t
%   of order 1 to 100) costs only a few extra trial steps per jump.
%
%   INPUTS
%   flow_map     - Function handle,
%                  [xcdot,xd_next] = flow_map(tj,dt,xc,xd,parameters).
%                  Given a trial step of length dt from the state (xc,xd)
%                  at hybrid time tj, returns the time derivative of xc as
%                  a column vector the same size as xc, and the discrete
%                  state at the end of the step as a column vector the
%                  same size as xd. Return xd unchanged for any discrete
%                  state that is constant during flows. flow_map is called
%                  once per stage (s times) for every trial step,
%                  including those of an event search and those rejected
%                  by the Lyapunov control, and always with two outputs
%                  requested, so it must not have side effects that assume
%                  each call is an accepted step. In stage calls, tj(1) is
%                  the stage time, dt is the full trial step length, and xd
%                  is the start-of-step discrete state.
%   jump_map     - Function handle,
%                  [xc_plus,xd_plus] = jump_map(tj,xc,xd,parameters).
%                  Returns the post-jump states as column vectors. It is
%                  evaluated at the pre-jump hybrid time.
%   jump_guard   - Function handle, g = jump_guard(tj,xc,xd,parameters).
%                  Returns the guard value g as a real double scalar; the
%                  state is in the jump set when g <= 0. g must not be NaN,
%                  and it must not be a logical (a former jump_set
%                  function returning g <= 0 must instead return g). For
%                  accurate event location, g should be continuous along
%                  the flow, ideally smooth, and cross zero transversally.
%   manifold_map - Function handle,
%                  [xc,xd] = manifold_map(tj,xc,xd,parameters).
%                  Maps the state onto the system manifold (e.g., keep
%                  quaternions on $\mathbb{S}^3$), overwriting both states.
%                  Use an identity map if no projection is needed.
%   t0           - Initial time. Real, finite scalar.
%   tf           - Final time. Real, finite scalar with tf >= t0.
%   dT           - Baseline (maximum) flow step size. Real, finite,
%                  positive scalar.
%   j0           - Initial jump index. Nonnegative integer.
%   Jmax         - Maximum jump index; the simulation stops when j reaches
%                  Jmax. Finite nonnegative integer (it is used for
%                  preallocation).
%   event_tol    - Event-location tolerance in t. Real, finite, positive
%                  scalar. Each event search ends once the bracket on the
%                  step length is no wider than event_tol. A value of
%                  event_tol >= dT disables the search: a step that ends
%                  in D is accepted as is, and the jump is applied at its
%                  end.
%   xc0          - Initial continuous state, nc-by-1 real, finite vector.
%                  A row vector is converted to a column.
%   xd0          - Initial discrete state, nd-by-1 real, finite vector.
%                  A row vector is converted to a column. Pass zeros(0,1)
%                  if there is no discrete state.
%   parameters   - Passed unchanged to every function handle. Any type.
%   display_flag - (Optional, default false) Logical scalar. When true,
%                  prints the current hybrid time [t j] and the row count
%                  k at the start, every 1000 recorded rows, and at the
%                  end, followed by "Simulation Complete!". Progress is
%                  counted in recorded rows, not in t, so the messages are
%                  not evenly spaced in time.
%   method       - (Optional, default "rk4") Runge-Kutta method for flows.
%                  Either one of the names below, or a struct with fields
%                  A (s-by-s, strictly lower triangular), b (length s,
%                  summing to 1), and optionally c (length s, equal to the
%                  row sums of A; computed from A if omitted) and order
%                  (reported in stats; NaN if omitted).
%                    "euler"    - forward Euler,              s = 1, order 1
%                    "midpoint" - explicit midpoint,          s = 2, order 2
%                    "heun"     - Heun (explicit trapezoid),  s = 2, order 2
%                    "ralston"  - Ralston,                    s = 2, order 2
%                    "ssprk3"   - Shu-Osher SSP RK3,          s = 3, order 3
%                    "rk4"      - classical Runge-Kutta,      s = 4, order 4
%                    "rk38"     - Kutta's 3/8 rule,           s = 4, order 4
%                  Names are case-insensitive.
%   lyap         - (Optional, default struct(), which disables the
%                  control) Struct of Lyapunov step-size control settings.
%                  Fields:
%       tol          - Required. Real scalar tol >= 0, in units of V per
%                      unit time (1/time when relative is true); each
%                      step's violation is bounded by c = tol*dT. tol = Inf
%                      computes the diagnostics in stats without
%                      shortening any step.
%       V            - Required. Function handle,
%                      v = V(tj,xc,xd,parameters), returning the Lyapunov
%                      function value as a real, finite scalar.
%       k            - Required. Certificate decay rate, real, finite
%                      scalar k >= 0, from Vdot <= -k*V.
%       min_step     - (Default 1e-4*dT) Shortest step length the control
%                      may impose, 0 < min_step <= dT. If you hit it in a
%                      convergence study, scale it with dT as well.
%       relative     - (Default false) Logical scalar. When true, the
%                      violation is divided by max(V(t_k),V_floor), with
%                      V(t_k) the value at the start of the step (see
%                      SCALE OF THE BOUND).
%       V_floor      - (Default 1e-10 times V at the start of the first
%                      flow step, or realmin if that is 0) Real scalar
%                      V_floor > 0. Level of V below which the relative
%                      bound becomes absolute. Set it to the size of V
%                      that roundoff makes meaningless, and set it
%                      explicitly if the flow starts at the equilibrium.
%
%   flow_map must be a scalar function handle that accepts the five inputs
%   (tj,dt,xc,xd,parameters). Each of jump_map, jump_guard, and
%   manifold_map (and lyap.V) must be a scalar function handle that
%   accepts the four inputs (tj,xc,xd,parameters). In every call, tj = [t j] is a 1-by-2
%   row vector holding the hybrid time of the state being evaluated.
%
%   OUTPUTS
%   tj_history       - K-by-2 array. Row k is the hybrid time [t j] of the
%                      k-th recorded state.
%   xc_history       - K-by-nc array. Row k is the continuous state at
%                      tj_history(k,:), stored as a row.
%   xd_history       - K-by-nd array. Row k is the discrete state at
%                      tj_history(k,:), stored as a row.
%   return_condition - Integer code giving the reason the simulation
%                      stopped:
%       1 - The initial hybrid time already satisfies t0 = tf or
%           j0 >= Jmax. Only the mapped initial state is returned, and no
%           jumps are applied, even if the initial state is in the jump
%           set.
%       2 - j reached Jmax during the initial jump sequence at t0, before
%           any flow.
%       3 - j reached Jmax during a jump sequence that followed a flow.
%       4 - The flow reached t = tf with j below Jmax.
%   stats            - Struct of work counters and diagnostics:
%       method            - Method name ("custom" for a struct tableau).
%       stages            - Number of stages s.
%       order             - Classical order of the method.
%       accepted_steps    - Accepted flow steps.
%       event_searches    - Event searches performed (flow steps that
%                           ended in D).
%       search_trials     - Trial steps taken inside event searches, in
%                           addition to the step that triggered each
%                           search.
%       lyap_rejections   - Trial steps rejected by the Lyapunov control.
%       lyap_floor_steps  - Steps accepted at lyap.min_step while still
%                           violating the bound.
%       defect_integral   - Sum of the violations of accepted steps, in
%                           the units of the bound (0 when lyap is not
%                           given; steps that end an event search are
%                           not checked and add 0).
%       jumps             - Jumps applied.
%       flow_map_calls    - Total flow_map evaluations (s per trial step).
%       V_calls           - Total lyap.V evaluations.
%       guard_calls       - Total jump_guard evaluations.
%       min_step          - Shortest accepted flow step (Inf if none).
%
%   One row is recorded per accepted flow step and one per jump, so the
%   pre-jump and post-jump states of each jump share the same t and
%   differ by one in j. Trial steps that are not accepted are not
%   recorded.
%
%   NOTES
%   - The requirements listed under INPUTS and OUTPUTS are enforced by the
%     arguments blocks below, except the requirements on the values
%     returned by the function handles, which are checked at every call,
%     and the fields of lyap, which are checked when the function starts.
%     The Output block requires MATLAB R2022b or later; the Input block
%     requires R2019b or later.
%   - If more than the preallocated number of rows is needed, MATLAB
%     grows the arrays automatically. The results are unaffected; only
%     performance is.

%% Validate inputs and outputs
arguments (Input)
    % Hybrid system data. flow_map is called as f(tj,dt,xc,xd,parameters);
    % the other three maps are called as f(tj,xc,xd,parameters).
    flow_map     (1,1) function_handle {mustAcceptInputs(flow_map,5)}
    jump_map     (1,1) function_handle {mustAcceptInputs(jump_map,4)}
    jump_guard   (1,1) function_handle {mustAcceptInputs(jump_guard,4)}
    manifold_map (1,1) function_handle {mustAcceptInputs(manifold_map,4)}
    % Time horizon and baseline step size.
    t0           (1,1) double {mustBeReal, mustBeFinite}
    tf           (1,1) double {mustBeReal, mustBeFinite, mustBeGreaterThanOrEqual(tf,t0)}
    dT           (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
    % Jump indices and event-location tolerance.
    j0           (1,1) double {mustBeInteger, mustBeNonnegative}
    Jmax         (1,1) double {mustBeFinite, mustBeInteger, mustBeNonnegative}
    event_tol    (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
    % Initial states (row vectors are converted to columns).
    xc0          (:,1) double {mustBeReal, mustBeFinite}
    xd0          (:,1) double {mustBeReal, mustBeFinite}
    % Problem parameters, passed through unchanged.
    parameters
    % Progress display switch; may be omitted.
    display_flag (1,1) logical = false
    % Runge-Kutta method: a name or a Butcher-tableau struct; may be
    % omitted.
    method                     {mustBeRKMethod(method)} = "rk4"
    % Lyapunov step-size control settings; may be omitted. The fields
    % are checked by lyap_settings below.
    lyap         (1,1) struct = struct()
end
arguments (Output)
    % One row per recorded state.
    tj_history       (:,2) double
    xc_history       (:,:) double
    xd_history       (:,:) double
    % Termination code, see the header.
    return_condition (1,1) double {mustBeMember(return_condition,[1 2 3 4])}
    % Work counters and diagnostics, see the header.
    stats            (1,1) struct
end

%% Set up the Runge-Kutta method and the step-size control
% Butcher tableau: stage matrix A (strictly lower triangular), weights b
% (column), and nodes c (column). s is the number of stages.
[A,b,c,method_name,method_order] = butcher_tableau(method);
s = numel(b);

% Lyapunov step-size control settings, with defaults filled in. L.c is
% the per-step violation bound (Inf when the control is off).
L = lyap_settings(lyap,dT);

% Work counters and diagnostics returned in stats.
stats = struct( ...
    'method',            method_name, ...
    'stages',            s, ...
    'order',             method_order, ...
    'accepted_steps',    0, ...
    'event_searches',    0, ...
    'search_trials',     0, ...
    'lyap_rejections',   0, ...
    'lyap_floor_steps',  0, ...
    'defect_integral',   0, ...
    'jumps',             0, ...
    'flow_map_calls',    0, ...
    'V_calls',           0, ...
    'guard_calls',       0, ...
    'min_step',          Inf);
% Set once the min_step warning has been issued in this call.
floor_warned = false;

%% Preallocate the output histories
% Estimated row count: one row per nominal flow step over [t0,tf], plus,
% for each allowed jump, the jump itself and the shortened step that ends
% at the event. Steps shortened by the Lyapunov control are not counted.
K0 = ceil((tf-t0)/dT) + 2*Jmax;
% Display interval in recorded rows. Only used when progress is
% displayed; every use below is guarded by display_flag && ..., which
% short-circuits.
dK = 1000;
k = 0;                                  % number of rows recorded so far

tj_history = nan(K0,2);
xc_history = nan(K0,size(xc0,1));
xd_history = nan(K0,size(xd0,1));

%% Record the initial state
tj0 = [t0 j0];
% Report the starting hybrid time.
if display_flag
    report_progress(tj0,k)
end
% Map the initial condition onto the manifold before recording it.
[xc0,xd0] = manifold_map(tj0,xc0,xd0,parameters);

k = k + 1;
tj_history(k,:) = tj0;
xc_history(k,:) = xc0';
xd_history(k,:) = xd0';

% Current hybrid time and state. States are stored as rows in the
% histories but handled as column vectors in the computations.
tj = tj_history(k,:);
xctj = xc_history(k,:)';
xdtj = xd_history(k,:)';

% Stop immediately if the initial hybrid time is already at or beyond the
% horizon in either t or j.
if any(tj >= [tf Jmax])
    % Remove the unused preallocated rows.
    tj_history(k+1:end,:) = [];
    xc_history(k+1:end,:) = [];
    xd_history(k+1:end,:) = [];
    return_condition = 1;
    % Report the final hybrid time.
    if display_flag
        report_progress(tj,k)
        disp("Simulation Complete!")
    end
    return
end

%% Initial jump sequence
% Guard value at the current state. After every jump sequence, it is
% positive, which is what the event search relies on.
gtj = eval_guard(jump_guard,tj,xctj,xdtj,parameters);
stats.guard_calls = stats.guard_calls + 1;
% Apply jumps at t0 until the state leaves the jump set or j reaches Jmax.
while gtj <= 0
    k = k + 1;
    stats.jumps = stats.jumps + 1;
    % A jump holds t fixed and increments j.
    tjp1 = tj + [0 1];
    tj_history(k,:) = tjp1;
    % Evaluate the jump map at the pre-jump hybrid time, then map the
    % post-jump state onto the manifold at the post-jump hybrid time.
    [xctjp1,xdtjp1] = jump_map(tj,xctj,xdtj,parameters);
    [xctjp1,xdtjp1] = manifold_map(tjp1,xctjp1,xdtjp1,parameters);
    xc_history(k,:) = xctjp1';
    xd_history(k,:) = xdtjp1';

    % Advance the current hybrid time and state to the post-jump values.
    tj = tjp1;
    xctj = xctjp1;
    xdtj = xdtjp1;

    % Report progress every dK recorded rows.
    if display_flag && mod(k, dK) == 0
        report_progress(tj,k)
    end

    % Stop if the jump limit has been reached.
    if tj(2) >= Jmax
        % Remove the unused preallocated rows.
        tj_history(k+1:end,:) = [];
        xc_history(k+1:end,:) = [];
        xd_history(k+1:end,:) = [];
        return_condition = 2;
        % Report the final hybrid time.
        if display_flag
            report_progress(tj,k)
            disp("Simulation Complete!")
        end
        return
    end

    % Guard value at the post-jump state.
    gtj = eval_guard(jump_guard,tj,xctj,xdtj,parameters);
    stats.guard_calls = stats.guard_calls + 1;
end

%% Flow with event location and Lyapunov step-size control
% Lyapunov function value at the current state, used by the certificate
% defect of the next step (NaN when the control is off).
[Vtj,stats] = eval_V(L,tj,xctj,xdtj,parameters,stats);
% Default floor of the relative bound: a small fraction of V at the start
% of the first flow step.
if L.relative && isnan(L.V_floor)
    L.V_floor = max(1e-10*Vtj,realmin);
end
% Step length planned for the next step, before the tf clamp. It is dT
% except while the Lyapunov control is shortening steps.
dt_next = dT;
while tj(1) < tf

    % Trial step of the planned length, clamped so it does not pass tf.
    % The Lyapunov control shortens and retries it until its violation is
    % within the bound or its length reaches the floor.
    dt = clamp_to_tf(dt_next,tj(1),tf,dT);
    while true
        T = rk_trial(flow_map,manifold_map,jump_guard,parameters, ...
            A,b,c,L,tf,tj,dt,xctj,xdtj,Vtj,true);
        stats.flow_map_calls = stats.flow_map_calls + s;
        stats.guard_calls = stats.guard_calls + 1;
        stats.V_calls = stats.V_calls + L.enabled;
        if T.violation <= L.c
            % Within the bound (always, when the control is off).
            break
        end
        if dt > L.min_step
            % Reject the trial step and retry with a shorter one.
            stats.lyap_rejections = stats.lyap_rejections + 1;
            dt = max(L.min_step, dt*max(0.1,min(0.5,0.9*L.c/T.violation)));
            continue
        end
        % The step is already at the floor: accept it anyway and warn
        % the first time this happens.
        stats.lyap_floor_steps = stats.lyap_floor_steps + 1;
        if ~floor_warned
            warning('hybrid_RK:lyapMinStep', ...
                ['At t = %.10g, a step of lyap.min_step = %g still ' ...
                 'violates the Lyapunov bound (violation %g > c = %g). ' ...
                 'Continuing at min_step while the violation persists, ' ...
                 'which can be slow; later occurrences are counted in ' ...
                 'stats.lyap_floor_steps. If this persists, check that ' ...
                 'the certificate holds along the true flow (e.g., that ' ...
                 'k is not too large).'], ...
                tj(1),L.min_step,T.violation,L.c);
            floor_warned = true;
        end
        break
    end

    % If the step ends in the jump set, locate where the flow enters it.
    % The step started with gtj > 0 and ended with T.g <= 0, so [0,dt]
    % brackets the event. The search returns the shortest step it found
    % that still ends in D, at most event_tol longer than the event.
    if T.g <= 0
        stats.event_searches = stats.event_searches + 1;
        % Search trials are not checked by the Lyapunov control.
        trial_tau = @(tau) rk_trial(flow_map,manifold_map,jump_guard, ...
            parameters,A,b,c,L,tf,tj,tau,xctj,xdtj,Vtj,false);
        [T,n_trials] = locate_event(trial_tau,gtj,dt,T,event_tol);
        stats.search_trials = stats.search_trials + n_trials;
        stats.flow_map_calls = stats.flow_map_calls + s*n_trials;
        stats.guard_calls = stats.guard_calls + n_trials;
    end

    % Accept the step.
    dt_accepted = T.tj(1) - tj(1);
    k = k + 1;
    stats.accepted_steps = stats.accepted_steps + 1;
    stats.min_step = min(stats.min_step,dt_accepted);
    stats.defect_integral = stats.defect_integral + T.violation;
    tj_history(k,:) = T.tj;
    xc_history(k,:) = T.xc';
    xd_history(k,:) = T.xd';

    % Advance the current hybrid time and state to the end of the step.
    tj = T.tj;
    xctj = T.xc;
    xdtj = T.xd;
    gtj = T.g;
    Vtj = T.V;

    % Report progress every dK recorded rows.
    if display_flag && mod(k, dK) == 0
        report_progress(tj,k)
    end

    % Plan the next step from the trial length dt.
    dt_next = plan_next_step(dt,T.violation,L,dT);

    % Apply jumps until the state leaves the jump set or j reaches Jmax.
    % This runs only after a step that ended in D. (Rows added here are
    % counted in k but do not trigger progress messages.)
    jumped = false;
    while gtj <= 0
        jumped = true;
        k = k + 1;
        stats.jumps = stats.jumps + 1;
        % A jump holds t fixed and increments j.
        tjp1 = tj + [0 1];
        tj_history(k,:) = tjp1;
        % Evaluate the jump map at the pre-jump hybrid time, then map the
        % post-jump state onto the manifold at the post-jump hybrid time.
        [xctjp1,xdtjp1] = jump_map(tj,xctj,xdtj,parameters);
        [xctjp1,xdtjp1] = manifold_map(tjp1,xctjp1,xdtjp1,parameters);
        xc_history(k,:) = xctjp1';
        xd_history(k,:) = xdtjp1';

        % Advance the current hybrid time and state to the post-jump
        % values.
        tj = tjp1;
        xctj = xctjp1;
        xdtj = xdtjp1;

        % Stop if the jump limit has been reached.
        if tj(2) >= Jmax
            % Remove the unused preallocated rows.
            tj_history(k+1:end,:) = [];
            xc_history(k+1:end,:) = [];
            xd_history(k+1:end,:) = [];
            return_condition = 3;
            % Report the final hybrid time.
            if display_flag
                report_progress(tj,k)
                disp("Simulation Complete!")
            end
            return
        end

        % Guard value at the post-jump state.
        gtj = eval_guard(jump_guard,tj,xctj,xdtj,parameters);
        stats.guard_calls = stats.guard_calls + 1;
    end
    if jumped
        % The state changed discontinuously, so restart from a full step
        % and evaluate V at the post-jump state.
        dt_next = dT;
        [Vtj,stats] = eval_V(L,tj,xctj,xdtj,parameters,stats);
    end

end

%% Final time reached by flow
% Remove the unused preallocated rows.
tj_history(k+1:end,:) = [];
xc_history(k+1:end,:) = [];
xd_history(k+1:end,:) = [];
return_condition = 4;
% Report the final hybrid time.
if display_flag
    report_progress(tj,k)
    disp("Simulation Complete!")
end

end

%% Local functions

function T = rk_trial(flow_map,manifold_map,jump_guard,parameters,A,b,c,L,tf,tj,dt,xctj,xdtj,Vtj,check)
%RK_TRIAL One Runge-Kutta trial step, mapped to the manifold and checked.
%   Takes a step of length dt from (xctj,xdtj) at hybrid time tj with the
%   tableau (A,b,c), passes the end state through manifold_map, and
%   evaluates the guard there. When the Lyapunov control is on and check
%   is true, it also evaluates lyap.V at the end of the step and the
%   certificate defect against Vtj, the value at the start. See
%   RUNGE-KUTTA FLOW STEP and LYAPUNOV STEP-SIZE CONTROL in the header.
%
%   Returns a struct with fields
%       tj, xc, xd - hybrid time and states at the end of the step,
%       g          - guard value at the end of the step,
%       V          - V at the end of the step (NaN when not evaluated),
%       violation  - max(delta,0), scaled for the relative bound (0 when
%                    not checked).
s = numel(b);
Kc = zeros(numel(xctj),s);
% Stage 1 is evaluated at the start state, and it alone supplies the
% discrete state at the end of the step.
[Kc(:,1),xdtpdtj] = flow_map(tj,dt,xctj,xdtj,parameters);
for i = 2:s
    % Stage state: explicit combination of the earlier stages. The
    % discrete state is held at its start-of-step value.
    xci = xctj + dt*(Kc(:,1:i-1)*A(i,1:i-1)');
    % Both outputs are requested so that flow_map is always called the
    % same way; the stage's discrete-state output is discarded.
    [Kc(:,i),~] = flow_map(tj + [c(i)*dt 0],dt,xci,xdtj,parameters);
end
% End of the step. A step clamped by tf ends exactly at tf, even if
% tj(1) + dt rounds to a neighbor of tf.
tpdtj = tj + [dt 0];
if dt >= tf - tj(1)
    tpdtj(1) = tf;
end
xctpdtj = xctj + dt*(Kc*b);
[xctpdtj,xdtpdtj] = manifold_map(tpdtj,xctpdtj,xdtpdtj,parameters);
gtpdtj = eval_guard(jump_guard,tpdtj,xctpdtj,xdtpdtj,parameters);

% Certificate defect, from V at both ends of the step.
Vtpdtj = NaN;
violation = 0;
if L.enabled && check
    Vtpdtj = check_V(L.V(tpdtj,xctpdtj,xdtpdtj,parameters));
    % delta = exp(k*dt)*V(t_{k+1}) - V(t_k), the exp(k*(t - t_k))-weighted
    % integral of Vdot + k*V along the step.
    violation = max(exp(L.k*(tpdtj(1) - tj(1)))*Vtpdtj - Vtj, 0);
    if L.relative
        % Measure against the current size of V.
        violation = violation/max(max(Vtj,L.V_floor),realmin);
    end
end

T = struct('tj',tpdtj,'xc',xctpdtj,'xd',xdtpdtj,'g',gtpdtj, ...
    'V',Vtpdtj,'violation',violation);
end

function dt_next = plan_next_step(dt,violation,L,dT)
%PLAN_NEXT_STEP Length of the next step after an accepted step of length dt.
%   A violating step predicts the length at which the next step's
%   violation would be 0.9*c, assuming the violation scales with the step
%   length; otherwise steps double back toward dT. The result is never
%   below lyap.min_step (0 when the control is off). It is not yet clamped
%   by tf.
dt_next = min(2*dt,dT);
if violation > 0 && isfinite(L.c)
    dt_next = min(dt_next,0.9*L.c*dt/violation);
end
dt_next = max(dt_next,min(L.min_step,dT));
end

function [Tb,n_trials] = locate_event(trial_tau,g0,dt,Tb,tol)
%LOCATE_EVENT Locate where the flow enters the jump set along one step.
%   trial_tau(tau) returns the trial-step struct of a step of length tau
%   from the current state (see rk_trial). Shrinks a bracket [tau_a,tau_b]
%   on the step length tau, starting from [0,dt], over which the guard at
%   the end of a trial step goes from positive (g0 > 0 at tau = 0, the
%   current state) to nonpositive (Tb.g <= 0 at tau = dt). Stops when
%   tau_b - tau_a <= tol and returns the trial step of length tau_b,
%   which ends in the jump set, along with the number of trial steps
%   taken.
%
%   Each new point comes from the Illinois variant of regula falsi: the
%   secant through the two bracket ends, with the value at an end that has
%   been kept for two iterations in a row halved so that the method keeps
%   converging superlinearly instead of stalling at one end. Two
%   safeguards make the search reliable for any continuous guard:
%   - Each point is kept at least tol/2 inside the bracket. Once the
%     secant estimate is accurate, the next point therefore lands within
%     tol of the root on one side or the other, and the search ends.
%   - If two iterations together have not halved the bracket, the next
%     point is the midpoint. The bracket therefore halves at least every
%     second iteration, which bounds the search at about
%     2*log2(dt/tol) trial steps.
tau_a = 0;
tau_b = dt;
fa = g0;                % interpolation values (Illinois-scaled copies)
fb = Tb.g;
kept = 0;               % +1 if tau_b was kept last time, -1 if tau_a was
width_hist = [Inf Inf]; % bracket widths before the last two iterations
n_trials = 0;
while tau_b - tau_a > tol
    width = tau_b - tau_a;
    if width > 0.5*width_hist(1)
        % Two iterations have not halved the bracket, so bisect.
        tau_m = tau_a + 0.5*width;
    else
        % Secant (regula falsi) point. fb <= 0 < fa, so the denominator
        % is negative and the point lies in [tau_a,tau_b].
        tau_m = tau_b - fb*(tau_b - tau_a)/(fb - fa);
        tau_m = min(max(tau_m,tau_a + 0.5*tol),tau_b - 0.5*tol);
    end
    width_hist = [width_hist(2) width];

    Tm = trial_tau(tau_m);
    n_trials = n_trials + 1;

    if Tm.g <= 0
        % The new point is in the jump set: it becomes the right end.
        tau_b = tau_m;
        Tb = Tm;
        fb = Tm.g;
        if kept == -1
            % tau_a has now been kept twice in a row, so halve its value.
            fa = 0.5*fa;
        end
        kept = -1;
    else
        % The new point is outside the jump set: it becomes the left end.
        tau_a = tau_m;
        fa = Tm.g;
        if kept == +1
            % tau_b has now been kept twice in a row, so halve its value.
            fb = 0.5*fb;
        end
        kept = +1;
    end
end
end

function L = lyap_settings(lyap,dT)
%LYAP_SETTINGS Check the lyap struct and fill in defaults.
%   Returns a struct with fields enabled, c (per-step violation bound,
%   Inf when disabled), min_step, V, k, relative, and V_floor (NaN until
%   its default is set at the first flow step). The checks run here rather
%   than in the arguments block so that they also run where arguments
%   blocks are not supported.
allowed = {'tol','V','k','min_step','relative','V_floor'};
names = fieldnames(lyap);
unknown = setdiff(names,allowed);
if ~isempty(unknown)
    error('hybrid_RK:badLyap', ...
        'Unknown lyap field "%s". Valid fields are: %s.', ...
        unknown{1},strjoin(allowed,', '));
end

L = struct('enabled',false,'c',Inf,'min_step',0,'V',[],'k',0, ...
    'relative',false,'V_floor',NaN);
if isempty(names)
    return
end

if ~all(isfield(lyap,{'tol','V','k'}))
    error('hybrid_RK:badLyap', ...
        'lyap.tol, lyap.V, and lyap.k are required to enable the Lyapunov control.');
end
if ~(isa(lyap.V,'function_handle') && isscalar(lyap.V))
    error('hybrid_RK:badLyap', ...
        'lyap.V must be a function handle V(tj,xc,xd,parameters).');
end
n = nargin(lyap.V);
if ~(n == 4 || (n < 0 && abs(n) - 1 <= 4))
    error('hybrid_RK:badLyap', ...
        'lyap.V must accept the four inputs (tj,xc,xd,parameters).');
end
if ~(isnumeric(lyap.tol) && isscalar(lyap.tol) && isreal(lyap.tol) ...
        && lyap.tol >= 0)
    error('hybrid_RK:badLyap', ...
        'lyap.tol must be a real scalar >= 0 (Inf is allowed).');
end
kk = lyap.k;
if ~(isnumeric(kk) && isscalar(kk) && isreal(kk) && isfinite(kk) && kk >= 0)
    error('hybrid_RK:badLyap', ...
        'lyap.k must be a real, finite scalar >= 0.');
end
L.enabled = true;
L.c = double(lyap.tol)*dT;
L.V = lyap.V;
L.k = double(kk);

L.min_step = 1e-4*dT;
if isfield(lyap,'min_step')
    m = lyap.min_step;
    if ~(isnumeric(m) && isscalar(m) && isreal(m) && m > 0 && m <= dT)
        error('hybrid_RK:badLyap', ...
            'lyap.min_step must be a real scalar with 0 < min_step <= dT.');
    end
    L.min_step = double(m);
end

if isfield(lyap,'relative')
    u = lyap.relative;
    if ~((islogical(u) || isnumeric(u)) && isscalar(u) && (u == 0 || u == 1))
        error('hybrid_RK:badLyap', ...
            'lyap.relative must be a logical scalar.');
    end
    L.relative = logical(u);
end
if isfield(lyap,'V_floor')
    v = lyap.V_floor;
    if ~(isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v) && v > 0)
        error('hybrid_RK:badLyap', ...
            'lyap.V_floor must be a real, finite scalar > 0.');
    end
    L.V_floor = double(v);
end
end

function [V,stats] = eval_V(L,tj,xc,xd,parameters,stats)
%EVAL_V Lyapunov function value at a state, or NaN when the control is off.
if ~L.enabled
    V = NaN;
    return
end
V = check_V(L.V(tj,xc,xd,parameters));
stats.V_calls = stats.V_calls + 1;
end

function V = check_V(V)
%CHECK_V Check a value returned by lyap.V.
if ~(isnumeric(V) && isscalar(V) && isreal(V) && isfinite(V))
    error('hybrid_RK:badLyapunovValue', ...
        'lyap.V must return a real, finite numeric scalar.');
end
V = double(V);
end

function g = eval_guard(jump_guard,tj,xc,xd,parameters)
%EVAL_GUARD Evaluate the guard and check that its value is usable.
g = jump_guard(tj,xc,xd,parameters);
if islogical(g)
    error('hybrid_RK:logicalGuard', ...
        ['jump_guard returned a logical. It must return the guard value ' ...
         'g itself; the jump set is {g <= 0}.']);
end
if ~(isnumeric(g) && isscalar(g) && isreal(g) && ~isnan(g))
    error('hybrid_RK:badGuard', ...
        'jump_guard must return a real, non-NaN numeric scalar.');
end
g = double(g);
end

function report_progress(tj,k)
%REPORT_PROGRESS Print the current hybrid time and recorded row count.
disp("Simulation Progress: [t j k] = ["+tj(1)+" "+tj(2)+" "+k+"]")
end

function dt = clamp_to_tf(dt,t,tf,dT)
%CLAMP_TO_TF Clamp a step so it does not pass tf.
%   Returns tf - t if dt would reach, pass, or stop within 1e-10*dT of tf;
%   otherwise returns dt unchanged.
remaining = tf - t;
if remaining - dt <= 1e-10*dT
    dt = remaining;
end
end

function [A,b,c,name,order] = butcher_tableau(method)
%BUTCHER_TABLEAU Return the Butcher tableau of an explicit RK method.
%   method is a method name or a struct with fields A, b, and optionally
%   c and order (already checked by mustBeRKMethod). b and c are returned
%   as columns.
if isstruct(method)
    A = method.A;
    b = method.b(:);
    if isfield(method,'c')
        c = method.c(:);
    else
        c = sum(A,2);
    end
    if isfield(method,'order')
        order = method.order;
    else
        order = NaN;
    end
    name = "custom";
    return
end

name = lower(string(method));
switch name
    case "euler"
        A = 0;
        b = 1;
        order = 1;
    case "midpoint"
        A = [0   0
             1/2 0];
        b = [0; 1];
        order = 2;
    case "heun"
        A = [0 0
             1 0];
        b = [1/2; 1/2];
        order = 2;
    case "ralston"
        A = [0   0
             2/3 0];
        b = [1/4; 3/4];
        order = 2;
    case "ssprk3"
        A = [0   0   0
             1   0   0
             1/4 1/4 0];
        b = [1/6; 1/6; 2/3];
        order = 3;
    case "rk4"
        A = [0   0   0 0
             1/2 0   0 0
             0   1/2 0 0
             0   0   1 0];
        b = [1/6; 1/3; 1/3; 1/6];
        order = 4;
    case "rk38"
        A = [ 0   0  0 0
              1/3 0  0 0
             -1/3 1  0 0
              1  -1  1 0];
        b = [1/8; 3/8; 3/8; 1/8];
        order = 4;
end
c = sum(A,2);
end

function mustBeRKMethod(method)
%MUSTBERKMETHOD Validate a Runge-Kutta method name or Butcher tableau.
%   A name must be one of the built-in methods. A struct must define an
%   explicit, consistent tableau: A square and strictly lower triangular,
%   b with one weight per stage summing to 1, and, if given, c equal to
%   the row sums of A.
names = ["euler" "midpoint" "heun" "ralston" "ssprk3" "rk4" "rk38"];
if (ischar(method) && isrow(method)) || (isstring(method) && isscalar(method))
    if ~any(lower(string(method)) == names)
        error('hybrid_RK:badMethod', ...
            'Unknown method "%s". Valid names are: %s.', ...
            method, strjoin(names,', '));
    end
    return
end
if ~(isstruct(method) && isscalar(method))
    error('hybrid_RK:badMethod', ...
        'method must be a method name or a struct with fields A and b.');
end
if ~all(isfield(method,{'A','b'}))
    error('hybrid_RK:badMethod', ...
        'A method struct must have fields A and b.');
end
A = method.A;
b = method.b;
s = size(A,1);
if ~(isnumeric(A) && isreal(A) && ismatrix(A) && size(A,2) == s && s >= 1 ...
        && all(isfinite(A),'all'))
    error('hybrid_RK:badMethod', ...
        'method.A must be a real, finite, square matrix.');
end
if any(triu(A) ~= 0,'all')
    error('hybrid_RK:badMethod', ...
        'method.A must be strictly lower triangular (explicit method).');
end
if ~(isnumeric(b) && isreal(b) && isvector(b) && numel(b) == s ...
        && all(isfinite(b)))
    error('hybrid_RK:badMethod', ...
        'method.b must be a real, finite vector with one weight per stage.');
end
if abs(sum(b) - 1) > 1e-12
    error('hybrid_RK:badMethod', ...
        'method.b must sum to 1 (consistency).');
end
if isfield(method,'c')
    c = method.c;
    if ~(isnumeric(c) && isreal(c) && isvector(c) && numel(c) == s)
        error('hybrid_RK:badMethod', ...
            'method.c must be a real vector with one node per stage.');
    end
    if any(abs(c(:) - sum(A,2)) > 1e-12)
        error('hybrid_RK:badMethod', ...
            'method.c must equal the row sums of method.A.');
    end
end
if isfield(method,'order') && ~(isnumeric(method.order) && isscalar(method.order))
    error('hybrid_RK:badMethod', ...
        'method.order must be a numeric scalar.');
end
end

function mustAcceptInputs(f,n_required)
%MUSTACCEPTINPUTS Validate that f can be called with n_required inputs.
%   nargin(f) is the number of declared inputs, or -(n+1) when f has n
%   fixed inputs followed by varargin. Either form must admit n_required
%   inputs.
n = nargin(f);
if ~(n == n_required || (n < 0 && abs(n) - 1 <= n_required))
    error('hybrid_RK:badHandle', ...
        'Function handle %s must accept %d inputs.', ...
        func2str(f), n_required);
end
end