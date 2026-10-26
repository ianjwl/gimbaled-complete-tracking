# Known issues and implementation notes

Limitations, deliberate design choices, and fragilities in the code. Items
marked **caution** can silently produce wrong results if the code is edited
carelessly.

---

## 1. The inner-loop control law is written out in three places

**Caution.** The control law — the outer-loop expression for `a*`, the target
angular velocity `omega*` and its five terms, and the Lyapunov function `V` —
is implemented independently in:

| File | Purpose |
|---|---|
| `src/inner_loop_controller.m` | the law that is applied |
| `src/multirotor_jump_guard.m` | evaluates the switching test |
| `src/multirotor_lyapunov_function.m` | evaluates `V` for the step-size control |

The latter two rebuild the law inline rather than calling
`inner_loop_controller`, so that they can run without executing the full
controller cascade.

**Consequence.** Editing the law in one file and not the others leaves the
switching test and the Lyapunov step-size control evaluating a law that is no
longer being applied. The symptom is subtle: the guard fires at slightly
wrong times, or the step-size control shortens steps against a certificate
that no longer matches the dynamics.

**If you change the control law, change it in all three files.** A quick
sanity check after an edit is to evaluate `V` through both paths on a fixed
state and confirm that the two agree to machine precision.

---

## 2. Optimization Toolbox is a nominal dependency

`run_all.m` constructs `parameters.qp_options` with
`optimoptions('quadprog', ...)` unconditionally, so the Optimization Toolbox
must be installed to run anything — even though the configuration used for
the results sets `apply_QP = false` and exercises none of the seven
quadratic-program allocators in `src/control_allocation.m`.

The dependency is included to enable users to easily explore how the 
control law behaves when inputs are constrained. 

To remove the dependency:

1. Delete the seven `qp_version` branches in `src/control_allocation.m`, and
2. Remove `qp_options`, `qp_slack`, `lambda_qp`, and `qp_version` from the
   parameter block and from `parameters`.

The seven allocators trade thrust accuracy against torque accuracy under
per-motor bounds; none is used in the accompanying paper. 

---

## 3. Parameters hard-coded in the figure scripts

| Constant | Where | Value |
|---|---|---|
| Origin disk radius `rho` in the schematics | `figure_2.m`, `figure_3.m`, `figure_4.m` | 0.6 |
| `rho` in the simulation | `run_all.m` | 0.55 |
| Commanded free-fall window | `figure_6_through_9.m`, `figure_10.m` | [14.5, 17] s |
| Flip window | `figure_10.m` | [7.85, 8.85] s |

The schematic figures deliberately use `rho = 0.6` for illustration, which
differs from the simulation value `rho = U1/r = 2.75/5 = 0.55`. 

---

## 4. Deliberate numerical choices

**Forward-Euler increments on controller states.** The auxiliary direction
`mu_star` and the outer-loop controller state `alpha_hat_OL` are advanced with
a single forward-Euler increment, while the reference model uses an exact
zero-order-hold step. This is intentional: it mirrors a fixed-sample-rate
update on real hardware. The outer-loop increment is exact in the `"P"` mode
used for the results, because that mode's state rate is identically zero. In
`"PI"` and `"A"` modes the rate is nonzero and the step is first order. The 
forward euler increments could be replaced with high precision integration by 
moving the relevant controller states from the discrete time part of the 
system state to the continuous time part of the system state and modifying 
the control law so that it returns the derivatives of the states, rather than
their increments.

**Degenerate transition intervals.** The smooth step `zeta(x;a,b)` is
undefined when `a = b`, since the normalization divides by `b - a`. The only 
case where this can happen in the code is `a = b = 0`. When `a = b`, the code
uses `zeta(x;a,a) = sign(x - a)`, with the derivatives of `zeta` vanishing, 
while the integral `Z` uses `Z(x;a,a) = max(0, x - a)`. These two
agree on the whole domain of use, because the schedule argument
`alpha = ||a*||/r` is nonnegative (i.e., greater than or equal to `a = b = 0`) 
— but they are not a literal derivative/integral pair at that single boundary 
value. This case is reachable at `rho = 0`, i.e. `U1 = 0`.

**Attitude error extraction.** `attitude_matrix_to_quaternion.m` carries the
previous quaternion and selects the antipode that keeps the extraction
continuous across steps. If that previous value goes stale, the MRP error can
jump between coverings of `SO(3)`. The memory is written unconditionally on
every call, so it should not go stale, but any refactor that makes the write
conditional would reintroduce the risk.

---

## 5. Branches the reference trajectory never reaches

- The `a* = 0` branch of `inner_loop_controller.m` is not entered, because the
  outer loop never drives the commanded actuation to exactly zero. The branch
  exists for correctness and is entered only if you construct a trajectory
  that reaches it.
- The saturation level `s_lp` is reset on every call to
  `control_allocation.m` in the active configuration, so it carries no
  information between steps. The state slots are reserved so the controller
  state layout is unchanged if the allocators are enabled.

---

## 6. What the Lyapunov step-size control does and does not do

`hybrid_RK.m` can shorten flow steps where the computed solution violates
`Vdot <= -k V`. This promotes greater accuracy in regions where lower 
precision might wrongly suggest the certificate does not hold.  The 
certificate may not hold exactly for any step size since numerical 
inacuracies constitute a disturbance to the control system (and 
`Vdot <= -k V` only holds exactly in the absence of disturbances).  Reducing 
the step size whenever numerical inaccuracies cause significant violations 
to `Vdot <= -k V`, however, effectively reduces the size of the disturbance 
on that time period.  Therefore, the Lyapunov step-size control promotes 
accuracy locally (in time) by enforcing the certificate. Globally, accuracy 
is set by the base step size `dT`.  For more details, see the header of 
`hybrid_RK.m`.

If the certificate does not hold along the true flow — for instance if `k` is
set larger than the design supports — every step violates the bound and the
solver runs at its minimum step, which is slow and issues a warning.
`stats.defect_integral` accumulates the violations of accepted steps: if it
shrinks when `dT` is halved, the violation is step-size dependent and the
control is removing it.