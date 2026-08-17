# tpt-physics Example Gallery

A one-stop tour of every physics domain in the workspace, runnable as a single
binary:

```text
cargo run -p tpt-physics-gallery
```

Each line printed is a self-contained "hello world" for one crate. The gallery
binary lives in `crates/tpt-physics-gallery/` and calls into every workspace
crate end-to-end.

## Per-domain hello-world examples

| Domain | Crate | Example | What it shows |
| --- | --- | --- | --- |
| Materials | `tpt-physics-core` | `MaterialRegistry::with_defaults()` | Typed material database + serde |
| Solvers | `tpt-physics-solver` | `examples` / `cargo bench -p tpt-physics-solver` | CG / preconditioned-CG / GMRES |
| FEA | `tpt-physics-fea` | `cargo run --example beam` | Cantilever tip deflection (3-D beam) |
| DEM | `tpt-physics-dem` | `cargo run --example granular_pile` | Granular settling under gravity |
| CFD | `tpt-physics-cfd` | `cargo run --example cavity` | Lid-driven cavity / channel flow |
| AI | `tpt-physics-ai` | `cargo run --example rl_pendulum` | Differentiable pendulum + RL wrapper |

## What the gallery demonstrates

- **core** — looks up *Structural Steel* and prints its Young's and shear moduli.
- **solver** — solves a 2×2 SPD system with Conjugate Gradient.
- **fea** — builds a `ProblemSpec` (JSON-describable) box, solves static
  linear elasticity under self-weight, prints the top settlement.
- **dem** — drops two particles under gravity and reports the settled kinetic
  energy (Hertz–Mindlin contacts).
- **cfd** — drives a 2-D channel with the LBM solver (non-reflective open
  boundary) and prints the centre-line velocity.
- **ai** — evaluates the forward-mode-autodiff Jacobian of the harmonic
  oscillator and checks it against the analytic value.

## Benchmarks

Long-running performance tracking lives under each crate's `benches/`
directory and uses `criterion`:

```text
cargo bench -p tpt-physics-solver   # CG / preconditioned-CG / GMRES
cargo bench -p tpt-physics-dem      # >100k-particle step_par
cargo bench -p tpt-physics-cfd      # LBM lattice sweep
cargo bench -p tpt-physics-fea      # spacer benchmark
```
