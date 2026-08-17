# Example Gallery

Runnable examples live under each crate's `examples/` directory. Build/run any
of them with `cargo run --example <name> -p <crate>` (add `--release` for speed).

| Example | Crate | Domain | What it shows |
| --- | --- | --- | --- |
| `beam` | `tpt-physics-fea` | FEA | Cantilever Euler–Bernoulli beam tip deflection vs. analytic. |
| `cavity` | `tpt-physics-cfd` | CFD | Lid-driven cavity flow (experimental vortex). |
| `granular_pile` | `tpt-physics-dem` | DEM | Spheres settling under gravity to a quiescent pile. |
| `rl_pendulum` | `tpt-physics-ai` | AI | Differentiable pendulum as a Gym env + autodiff Jacobians. |
| `pile_cage_spacer` | `tpt-physics-fea` | FEA | End-to-end milestone: CAD→mesh→elasticity. |
| `spacer_benchmark` | `tpt-physics-fea` | FEA | Full-stack spacer solve timed end-to-end. |

Run the whole gallery with `./scripts/run_gallery.sh` (or
`scripts/run_gallery.ps1` on Windows). The runner executes every example and
prints its stdout so you can eyeball the physics.

## Validation tests (the real "gallery" of correctness)

The `tests/` folders are the authoritative demonstrations — each asserts
*qualitative physics*, not just "it compiles":

- **FEA** — `tet10` rigid-body / symmetry / PD; `beam3d` cantilever (L≠1);
  `shell4` rigid-body + simply-supported plate; nonlinear Total-Lagrangian
  force/tangent; J2 plasticity Cook's membrane.
- **DEM** — Hertz contact, spatial-hash pairs, SIMD narrow-phase, settling,
  pile-cage flow, SSI spacer, hopper discharge, random close packing,
  cohesive bonds + inter-particle heat, bincode checkpoint resume.
- **CFD** — Poiseuille (analytic), lid-driven cavity (experimental),
  flow-past-cylinder vortex shedding.
- **Solver** — CG / preconditioned CG, GMRES / preconditioned GMRES, Newmark-β,
  RK4, hardware dispatch.
- **AI** — differentiable harmonic oscillator + pendulum matching analytic
  Jacobians.
