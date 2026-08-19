# tpt-physics Example Gallery

A one-stop tour of every physics domain in the workspace, runnable as a single
binary:

```text
cargo run -p tpt-phys-gallery
```

Each line printed is a self-contained "hello world" for one crate. The gallery
binary (`crates/tpt-phys-gallery/src/main.rs`) calls into every workspace crate
end-to-end. Expanded, interactive versions live in each crate's `examples/`.

## Gallery runner (per-domain)

| Domain | Crate | What it shows |
| --- | --- | --- |
| Materials | `tpt-phys-core` | `MaterialRegistry::with_defaults()` — typed material DB + serde |
| DEM | `tpt-phys-dem` | two particles settling under gravity (Hertz–Mindlin) |
| CFD | `tpt-phys-cfd` | driven channel (LBM, non-reflective open boundary) |
| FSI | `tpt-phys-fsi` | partitioned fluid→structure loop (`LumpedStructure`) |
| Thermal-struct | `tpt-phys-thermal-struct` | `thermal_load_vector` on a tet mesh |
| Electro-thermal | `tpt-phys-electro-thermal` | `ElectroThermalRod` Joule heating |
| Orchestrator | `tpt-phys-orchestrator` | `build_demo_simulation` co-simulation |

## Per-crate examples

Run any with `cargo run --example <name> -p <crate>` (add `--release` for speed).

| Example | Crate | Domain | What it shows |
| --- | --- | --- | --- |
| `cavity` | `tpt-phys-cfd` | CFD | Lid-driven cavity (⚠️ experimental vortex) |
| `granular_pile` | `tpt-phys-dem` | DEM | Spheres settling under gravity |
| `rl_pendulum` | `tpt-phys-orchestrator` | AI | Differentiable pendulum + RL wrapper + AD Jacobians |
| `uq_coupled` | `tpt-phys-orchestrator` | UQ + co-sim | Monte-Carlo material scatter through the coupled `Simulation` |
| `uq_cantilever` | `tpt-phys-core` | UQ | Monte-Carlo scatter on a cantilever tip (`--features uq`) |

The full runner is `cargo run -p tpt-phys-gallery`.

## Validation tests (the real "gallery" of correctness)

`tests/` folders are the authoritative demonstrations — each asserts *qualitative
physics*, not just "it compiles":

- **DEM** — Hertz contact, spatial-hash pairs, SIMD narrow-phase, settling,
  pile-cage flow, SSI spacer, hopper discharge, random close packing, cohesive
  bonds + inter-particle heat, bincode checkpoint resume.
- **CFD** — Poiseuille (analytic), lid-driven cavity (⚠️ experimental),
  flow-past-cylinder vortex shedding, SPH dam break (free-surface).
- **FSI** — partitioned coupling driver: structure displaces under steady flow
  and relaxes when the flow stops (lumped model).
- **Electro-thermal** — Joule heating rod: heats under voltage, self-limits.
- **Thermal-struct** — `thermal_load_vector` integration test.
- **Orchestrator** — `SubModel` adapters drive a coupled `Simulation`.

## Benchmarks

Long-running performance tracking lives under each crate's `benches/`
directory and uses `criterion`:

```text
cargo bench -p tpt-phys-dem      # >100k-particle step_par
cargo bench -p tpt-phys-cfd      # LBM lattice sweep + SPH step
cargo bench -p tpt-phys-fsi      # FSI coupling-iteration driver
```
