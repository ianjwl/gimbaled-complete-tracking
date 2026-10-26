# gimbaled-complete-tracking
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23669210.svg)](https://doi.org/10.5281/zenodo.23669210)

Simulation and figure-generation code for the accompanying paper:

> I. J. Willebeek-LeMair and C. A. Woolsey, "Practical Complete Tracking
> Through Gimbaled Unidirectional Actuation with Application to Multirotor
> UAVs," [arXiv:2610.31962](https://doi.org/10.48550/arXiv.2610.31962), 2026.

The repository contains a MATLAB implementation of the control design, the
reference trajectory used in the simulation study, and the scripts that
generate every data and geometry figure.

---

## What this code does

This code simulates a multirotor unmanned aerial vehicle (UAV) flying with a 
recently developed control law.  The control law drives the actuation (the 
thrust) of a gimbaled unidirectional actuator (the UAV) to the
ball `A(t) = a*(t) + r B^3` centered on the target actuation, rather than to
an exact target value. This set-based objective removes the free-fall
singularity that leaves the target attitude undefined when the commanded
thrust passes through zero.

The implementation covers:

- the inner-loop backstepping design, including the target angular velocity
  and its five gated terms
- the smooth thrust schedule, in its radially symmetric and radially
  asymmetric variants
- the modified-Rodrigues-parameter attitude error, its Lyapunov function,
  and the switching condition that prevents unwinding
- the multirotor UAV specialization and the outer-loop state-feedback law
- a hybrid integrator with Runge-Kutta flows, event location, and Lyapunov
  step-size control

The simulation is a multirotor flying an obstacle course that includes a
vertical maneuver with two commanded free falls and a ballistic arc with a
commanded free fall.

---

## Attribution

The original version of every file in this repository was written entirely 
by Ian Willebeek-LeMair (excluding the `third_party/` files written by Dr. 
Jeremy Hopwood). Artificial intelligence tools (Claude Opus 5.5, 
Claude Fable 5.1, and DeepSeek-V4.1-Flash [hosted by the Virginia Tech ARC 
cluster]) were used to revise several files and to improve/expand/add 
headers and comments. Some functions (e.g., the `ssat.m` family and 
`hybrid_RK.m` utilities) were thoroughly modified while others (e.g., 
`inner_loop_controller.m`) received only minor changes. For the most part, 
AI-supported changes were made to improve numerical accuracy. Excluding 
`third_party/` files, the entire repository was checked for typos manually 
by the author and autonomously using the DeepSeek-V4.1-Flash model hosted 
by the Virginia Tech ARC cluster.  

---

## Requirements

| Requirement | Notes |
|---|---|
| MATLAB R2022b or later | `hybrid_RK.m` uses `arguments` blocks whose `Output` block validation requires R2022b |
| Control System Toolbox | Builds the reference model (`tf`, `tf2ss`) and the outer-loop gains (`lqr`) |
| Optimization Toolbox | Required as shipped, but not exercised. See `KNOWN_ISSUES.md` |

No other toolboxes are needed. 

---

## Quick start

Run everything from the repository root:

    % Simulate in the state-feedback mode used in the paper and draw every
    % figure, without writing any files:
    run_all

    % Write every figure to figures/ as PDF:
    run_all(true)

    % Same, but with the proportional-integral outer loop:
    run_all(true, "PI")

    % Same, but with the adaptive outer loop:
    run_all(true, "A")

    % Simulate only, returning the data with no figures drawn:
    sim = run_all();

`run_all` adds `src/`, `scripts/`, and `third_party/` to the MATLAB path,
builds the parameters, integrates the closed-loop system, postprocesses the
diagnostics, and then calls each figure function in turn. It returns the
postprocessed simulation data as a struct.

The `figures/` folder ships with the repository. `exportgraphics` does not
create it, so if you delete it, recreate it before running with
`save_figures = true`.

### `run_all` inputs

| Argument | Default | Meaning |
|---|---|---|
| `save_figures` | `false` | When `true`, every figure is exported to PDF under `figures/` |
| `OLcontrolMode` | `"P"` | Outer-loop law: `"P"` state feedback (used in the paper), `"PI"` proportional-integral, `"A"` adaptive |

All three outer-loop modes share the same inner loop, reference model, and
vehicle model. The `"P"` mode reproduces the gains reported in the paper. 
More information is available in the header for `run_all`.

### Working with the returned data

    sim = run_all();                       % simulate only, no figures drawn
    err = vecnorm(sim.ENU_positions - sim.ENU_reference_positions, 2, 2);
    plot(sim.t, err)                       % position-error norm history

`run_all.m` documents the contents of `sim` in full, grouped by role. The
fields the figure functions read are `t`, `positions`, `reference_positions`,
`attitudes`, `ENU_positions`, `ENU_reference_positions`, `ENU_Thrust`,
`ENU_Thrust_star`, `theta`, `theta_star`, and `uprightness`.

---

## Simulation settings

| Setting | Value |
|---|---|
| Integrator | `hybrid_RK` with classical RK4 flows |
| Base step `dT` | 0.01 s |
| Time horizon | 0 to 25 s |
| Event tolerance | 1e-12 |
| Jump index cap | 1e5 |
| Lyapunov step control | active, `tol = 0.1*min(k_a,k_omega)`, relative bound |
| Termination | `return_condition = 4` |

The vehicle parameters are based on AirSim's stock quadcopter model:
`m = 2.09` kg and `J = diag(0.0297, 0.0320, 0.0488)` kg m^2.

The Lyapunov step-size control shortens flow steps where the computed
solution violates the certificate `Vdot <= -k V`. It promotes accuracy by
enforcing the certificate.

---

## Repository layout

    .
    |-- run_all.m              Single entry point: data, diagnostics, figures
    |-- README.md
    |-- LICENSE
    |-- CITATION.cff
    |-- KNOWN_ISSUES.md
    |-- .gitignore
    |-- src/                   Controller, integrator, vehicle model, utilities
    |   `-- Contents.m         Index of every library file, grouped by role
    |-- scripts/               One function per figure
    |-- figures/               Rendered PDFs, committed for reference
    `-- third_party/           Vendored mesh and 3-D plotting utility

From the repository root, `help src` prints the `Contents.m` index.

---

## Figure map

| Script | Files written | Content |
|---|---|---|
| `scripts/figure_2.m` | `figure_2_a.pdf` ... `figure_2_d.pdf` | Geometry of `D(alpha)` and `rho B^2`, four configurations |
| `scripts/figure_3.m` | `figure_3.pdf` | Graphs of `theta_bar(alpha)` and `theta*(alpha)` |
| `scripts/figure_4.m` | `figure_4_a.pdf` ... `figure_4_d.pdf` | Annulus sections `C(alpha)`, four configurations |
| `scripts/figure_5.m` | `figure_5.pdf` | Simulated trajectory with rendered vehicle poses |
| `scripts/figure_6_through_9.m` | `figure_6.pdf` ... `figure_9.pdf` | Position, position-error norm, thrust, angular offset |
| `scripts/figure_10.m` | `figure_10.pdf` | Flip sequence and uprightness history |

`run_all(true)` additionally writes five diagnostic figures from the
simulation histories, `figures/diagnostic_figure_1.pdf` through
`diagnostic_figure_5.pdf`: the timing history, the MRP error against its
bound, the auxiliary-direction alignment, the Lyapunov value and its
actuation component, and the step-size history. These are not paper figures.
They are a record of the solver's behavior, and the step-size and Lyapunov
traces are the two to look at when auditing the integration.

---

## Reproducing the results

1. Clone the repository and start MATLAB in the repository root.
2. Run `run_all(true)`.
3. Compare the PDFs written to `figures/` against the committed copies.

The pipeline is deterministic: there is no random number generation anywhere,
so a fixed MATLAB release reproduces the same trajectories exactly. Record
the release you used if you report numbers derived from it.

---

## Third-party code

`third_party/FlightTest-Tools-main/` is not part of this work. It is
distributed under the MIT License, Copyright (c) 2023 Jeremy W. Hopwood, and
its license file is retained unmodified in that folder. Two components are
used:

- `plotObjectTrajectory.m`, for rendering vehicle poses in the trajectory and
  flip figures
- `quad.stl`, the multirotor mesh

Nothing in `third_party/` is on the control or simulation path. It affects
figures only.

---

## Citation

If you use this code, please cite the paper:

> I. J. Willebeek-LeMair and C. A. Woolsey, "Practical Complete Tracking Through 
> Gimbaled Unidirectional Actuation with Application to Multirotor UAVs," 
> Oct. 22, 2026, arXiv:2610.31962. doi: 10.48550/arXiv.2610.31962.

```bibtex
@misc{willebeek_lemair_practical_2026,
    title = {Practical {Complete} {Tracking} {Through} {Gimbaled} {Unidirectional} {Actuation} with {Application} to {Multirotor} {UAVs}},
    author = {Willebeek-LeMair, Ian J. and Woolsey, Craig A.},
    year = {2026},
    eprint = {2610.31962},
    archivePrefix = {arXiv},
    primaryClass = {eess.SY},
    url = {https://arxiv.org/abs/2610.31962},
    doi = {10.48550/arXiv.2610.31962},
}
```

The citation above will be updated once the journal version of the paper is released.

The code itself is archived on Zenodo. To cite it directly:

> I. J. Willebeek-LeMair, "Simulation code for Practical Complete Tracking
> Through Gimbaled Unidirectional Actuation with Application to Multirotor
> UAVs," Zenodo, 2026. doi: 10.5281/zenodo.23669210

```bibtex
@software{willebeek_lemair_code_2026,
    title     = {Simulation code for Practical Complete Tracking Through Gimbaled Unidirectional Actuation with Application to Multirotor {UAVs}},
    author    = {Willebeek-LeMair, Ian J.},
    year      = {2026},
    publisher = {Zenodo},
    doi       = {10.5281/zenodo.23669210},
    url       = {https://doi.org/10.5281/zenodo.23669210}
}
```

`CITATION.cff` supplies the metadata for GitHub's "Cite this repository" panel.

---

## License

This project is licensed under the MIT License; see `LICENSE`. The vendored
code under `third_party/` carries its own MIT license and its own copyright
holder.

---

## Contact

- Ian J. Willebeek-LeMair, Virginia Tech — ianwl@vt.edu — ianjwlemair@gmail.com

All reports of possible issues with the code are greatly appreciated.
