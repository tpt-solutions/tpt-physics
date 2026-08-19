# Example Gallery

Runnable examples live under each crate's `examples/` directory. Build/run any
of them with `cargo run --example <name> -p <crate>` (add `--release` for speed).

| Example | Crate | Domain | What it shows |
| --- | --- | --- | --- |
| `cavity` | `tpt-phys-cfd` | CFD | Lid-driven cavity flow (experimental vortex). |
| `granular_pile` | `tpt-phys-dem` | DEM | Spheres settling under gravity to a quiescent pile. |
| `rl_pendulum` | `tpt-phys-orchestrator` | AI | Differentiable pendulum as a Gym env + autodiff Jacobians. |
| `uq_cantilever` | `tpt-phys-core` | UQ | Monte-Carlo scatter on a cantilever tip deflection (`--features uq`). |

Run the whole gallery with `./scripts/run_gallery.sh` (or
`scripts/run_gallery.ps1` on Windows). The runner executes every example and
prints its stdout so you can eyeball the physics.

## Validation tests (the real "gallery" of correctness)

The `tests/` folders are the authoritative demonstrations — each asserts
*qualitative physics*, not just "it compiles":

- **DEM** — Hertz contact, spatial-hash pairs, SIMD narrow-phase, settling,
  pile-cage flow, SSI spacer, hopper discharge, random close packing,
  cohesive bonds + inter-particle heat, bincode checkpoint resume.
- **CFD** — Poiseuille (analytic), lid-driven cavity (experimental),
  flow-past-cylinder vortex shedding, SPH dam break (free-surface).
- **FSI** — partitioned coupling driver: structure displaces under steady flow
  and relaxes when the flow stops.
- **Electro-thermal** — Joule heating rod: heats under voltage, self-limits.
- **Thermal-struct** — `thermal_load_vector` integration test.
- **Orchestrator** — `SubModel` adapters drive a coupled `Simulation`.
