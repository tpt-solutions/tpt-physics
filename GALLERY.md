# tpt-physics Example Gallery

A one-stop tour of every physics domain in the workspace, runnable as a single
binary:

```text
cargo run -p tpt-phys-gallery
```

Each line printed is a self-contained "hello world" for one crate. The gallery
binary lives in `crates/tpt-phys-gallery/` and calls into every workspace
crate end-to-end.

## Per-domain hello-world examples

| Domain | Crate | Example | What it shows |
| --- | --- | --- | --- |
| Materials | `tpt-phys-core` | `MaterialRegistry::with_defaults()` | Typed material database + serde |
| DEM | `tpt-phys-dem` | `cargo run --example granular_pile` | Granular settling under gravity |
| CFD | `tpt-phys-cfd` | `cargo run --example cavity` | Lid-driven cavity / channel flow |
| FSI | `tpt-phys-fsi` | `crate::coupling` | Partitioned fluid–structure coupling driver |
| Thermal-struct | `tpt-phys-thermal-struct` | `thermal_load_vector` | Thermal-strain load on a tet mesh |
| Electro-thermal | `tpt-phys-electro-thermal` | `ElectroThermalRod` | Joule heating of a resistive rod |
| Orchestrator | `tpt-phys-orchestrator` | `build_demo_simulation` | Multiphysics co-simulation (`SubModel`s) |
| AI | `tpt-phys-orchestrator` | `cargo run --example rl_pendulum` | Differentiable pendulum + RL wrapper |

## What the gallery demonstrates

- **core** — looks up *Structural Steel* and prints its Young's and shear moduli.
- **dem** — drops two particles under gravity and reports the settled kinetic
  energy (Hertz–Mindlin contacts).
- **cfd** — drives a 2-D channel with the LBM solver (non-reflective open
  boundary) and prints the centre-line velocity.
- **fsi** — advances a partitioned fluid–structure loop (fluid traction →
  lumped-structure displacement → moving-wall feedback).
- **thermal-struct** — assembles the thermal-strain load vector for a tet.
- **electro-thermal** — heats a resistive rod under voltage (Joule heating).
- **orchestrator** — drives a `Simulation` coupling electro-thermal →
  thermal-structural models.
- **ai** — evaluates the forward-mode-autodiff Jacobian of the harmonic
  oscillator and checks it against the analytic value.

## Benchmarks

Long-running performance tracking lives under each crate's `benches/`
directory and uses `criterion`:

```text
cargo bench -p tpt-phys-dem      # >100k-particle step_par
cargo bench -p tpt-phys-cfd      # LBM lattice sweep + SPH step
cargo bench -p tpt-phys-fsi      # FSI coupling-iteration driver
```
