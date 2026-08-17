# tpt-physics

**A pure-Rust, AI-native engineering & physics simulation framework.**

TPT Solutions · Dual-licensed [MIT](LICENSE-MIT) / [Apache-2.0](LICENSE-APACHE-2.0)

`tpt-physics` is a ground-up, 100% pure-Rust rewrite of core physics
simulation paradigms — Finite Element Analysis (FEA), Discrete Element Method
(DEM), and Computational Fluid Dynamics (CFD) — with no C/C++ FFI, no
`bindgen`, and no system-level dependency hell. It is data-oriented,
unit-safe, and designed from day one to expose differentiable simulation state
to reinforcement-learning agents.

## Why

The engineering simulation industry is dominated by monolithic, decades-old
C++/Fortran codebases (OpenFOAM, CalculiX, ANSYS) that suffer from
build-system hell, opaque memory management, restrictive (GPL) licenses, and
poor integration with modern AI/ML workflows. `tpt-physics` replaces that with
a modular, composable, permissively-licensed Rust workspace.

## Design principles

- **100% Pure Rust.** Zero C/C++ FFI. Cross-compiles trivially to WASM,
  Linux, Windows, and macOS.
- **Strict MIT / Apache-2.0 licensing.** All dependencies are audited via
  [`cargo-deny`](deny.toml); any copyleft (GPL/LGPL/AGPL/SSPL) crate is
  rejected.
- **Compile-time unit safety.** SI units are typed through
  [`tpt-math-units`](https://github.com/tpt-solutions/tpt-math) (a wrap over
  `uom`), preventing catastrophic unit mismatches at compile time.
- **Data-oriented design.** Structure-of-arrays layouts and `rayon` parallelism
  for fearless concurrency on large problems.
- **AI-native & differentiable.** Simulation state is exposed as Gymnasium-like
  environments backed by [`tpt-math-autodiff`].

## Architecture

```
tpt-physics/
├── tpt-physics-core      # Material database + CAD ingestion adapter
├── tpt-physics-solver    # Iterative solvers (CG, GMRES*), time integration, HW dispatch*
├── tpt-physics-fea       # FEA elements (Tet10, 3D beam*, shell*), nonlinear (Tet4 only), thermal coupling
├── tpt-physics-dem       # Granular physics: Hertz-Mindlin, spatial hashing, SIMD
├── tpt-physics-cfd       # Lattice Boltzmann (D2Q9) incompressible flow
└── tpt-physics-ai        # Differentiable physics wrappers for RL agents
```
(\* = experimental: preconditioned GMRES, GPU dispatch backend, 3-D beam and
shell elements. Validated: Tet10, CG(+precond), nonlinear Tet4, J2 plasticity,
DEM, LBM cavity/cylinder, differentiable plants.)

## Crate-reuse map

Per the 2026-08-15 re-scope, `tpt-physics` reuses the sibling
[`tpt-math`](https://github.com/tpt-solutions/tpt-math) and
[`tpt-fem`](https://github.com/tpt-solutions/tpt-fem) workspaces **directly**
(no wrapper crates, no re-export shims). Each `tpt-physics-*` crate depends on
exactly the sibling crates it needs:

| Capability | Directly-used sibling crate |
| --- | --- |
| SI units | `tpt-math-units`, `tpt-math-units-dyn` |
| Dense linear algebra | `tpt-math-linalg-dense`, `tpt-math-linalg-fixed` |
| Geometry / autodiff | `tpt-math-geometry`, `tpt-math-autodiff` |
| Mesh + generation | `tpt-fem-mesh`, `tpt-fem-mesh-gen` |
| Reference elements / quadrature | `tpt-fem-element`, `tpt-fem-quadrature` |
| Sparse assembly + solve | `tpt-fem-sparse` |
| BCs / elasticity / thermal / eigen / solve | `tpt-fem-assembly`, `tpt-fem-elasticity`, `tpt-fem-thermal`, `tpt-fem-eigen`, `tpt-fem-solve` |

Only the genuinely net-new code lives in `tpt-physics-*`: the material
database, the CAD→`tpt-fem-mesh` ingestion adapter, the iterative solvers
(CG/GMRES) and time-integration schemes (no iterative/dynamic solver exists in
the siblings), the quadratic/3D-beam/shell elements, the geometric-nonlinear
framework, thermal-to-structural coupling, the DEM contact model, the LBM CFD
solver, and the RL environment wrappers.

## Building

```sh
cargo build --workspace
cargo test  --workspace
```

Requires a Rust toolchain ≥ 1.84. The sibling `tpt-math` and `t-fem` repos
must be checked out as siblings of this directory (the `[workspace.dependencies]`
in `Cargo.toml` reference them by relative path).

## Quickstart

```sh
# 1. Verify the sibling dependency workspaces are present.
./scripts/bootstrap.sh        # or: just setup

# 2. Build, test, and run the example gallery.
just test                     # cargo test --workspace
just examples                 # build every example
cargo run --release --example beam -p tpt-physics-fea
```

A minimal end-to-end FEA solve via the declarative JSON problem spec
(`tpt-physics-fea::spec`):

```rust
use tpt_physics_fea::spec::{DomainSpec, LoadSpec, ProblemSpec, SolverSpec};

let spec = ProblemSpec {
    materials: None,
    material: tpt_physics_fea::spec::MaterialRef::Inline(
        tpt_physics_core::Material::new("Steel", 200e9, 0.3, 7850.0, 12e-6),
    ),
    domain: DomainSpec::Box { min: [0.0,0.0,0.0], max: [1.0,1.0,1.0], n: [4,4,4] },
    boundary_conditions: tpt_physics_fea::spec::BcSpec {
        fixed_planes: vec!["y_min".to_string()],
    },
    loads: LoadSpec { self_weight: true, gravity: 9.81 },
    solver: SolverSpec::StaticLinear,
};
let solved = spec.solve(&tpt_physics_core::MaterialRegistry::new()).unwrap();
println!("top settlement = {:.3e} m", solved.free_top_settlement_y);
```

The spec round-trips through JSON (`ProblemSpec::from_json` /
`to_json`), so problems can be authored as data files rather than code.

## Validated vs. experimental

Capabilities are split by verification maturity. Anything *experimental* is
useful and tested for the cases noted, but has **not** been validated against a
full benchmark suite and may need more work before production use.

| Capability | Status | Notes |
| --- | --- | --- |
| Tet10 element, CG (+precond), GMRES | ✅ Validated | unit + integration tests |
| Geometric-nonlinear Tet4, J2 plasticity | ✅ Validated | Total-Lagrangian; Tet4 continuum only |
| DEM (Hertz–Mindlin, SIMD, >100k) | ✅ Validated | multiple physics tests |
| LBM (Poiseuille, flow-past-cylinder) | ✅ Validated | analytic + shedding benchmarks |
| Differentiable plants (oscillator, pendulum) | ✅ Validated | analytic Jacobian match |
| Declarative JSON problem spec | ✅ Validated | self-weight box solves |
| Cohesive bonds + inter-particle heat | ✅ Validated | new in this pass |
| 3-D beam (Euler–Bernoulli) element | ⚠️ Experimental | cantilever verified; no shear/flexible-beam benchmark |
| Mindlin–Reissner shell4 | ⚠️ Experimental | rigid-body + simply-supported plate; no Scordelis–Lo |
| Preconditioned GMRES | ⚠️ Experimental | Jacobi preconditioner; no AMG yet |
| GPU hardware-dispatch backend | ⚠️ Experimental | real `wgpu` WGSL `matvec` kernel behind `--features gpu`; falls back to `BackendUnavailable` when no adapter |
| Lid-driven cavity | ⚠️ Experimental | `#[ignore]`d — primary-vortex convergence issue |

## Roadmap

- **Phase 1 — Foundation & FEA MVP:** workspace, solvers, linear/nonlinear
  FEA. Milestone: simulate the 3D-printed pile cage spacer.
- **Phase 2 — Granular & performance:** DEM (concrete aggregate flow),
  GPU acceleration for >100k particles, nonlinear FEA extensions.
- **Phase 3 — Fluids & AI:** LBM CFD, differentiable gym environments,
  documentation, benchmarks, and the "Spacer Benchmark" case study.

See [`todo.md`](todo.md) for the full checklist.

## Validation status

Every non-trivial capability is backed by an integration test (or unit test)
that asserts the *qualitative physics*, not just "it compiles". Highlights:

| Area | Validation | What it proves |
| --- | --- | --- |
| DEM | `tests/granular_settling.rs` | a poured pile settles without blow-up or interpenetration |
| DEM | `tests/ssi_spacer.rs` | soil settles around an embedded cylindrical spacer, no penetration |
| DEM | `tests/pile_cage_flow.rs` | fluidized aggregate flows past a pile cage and comes to rest |
| DEM | `tests/hopper_discharge.rs` | discharge rate follows the Beverloo trend (arching when `D < d`) |
| DEM | `tests/random_close_packing.rs` | mono-disperse bed packs to the RCP fraction (~0.64) |
| DEM | `tests/large_scale.rs` | the `rayon` stepper advances **>100k** particles stably |
| FEA | `tests/cooks_membrane.rs` | J2 plasticity: more compliant than elastic, hardening stiffens |
| CFD | `tests/lid_driven_cavity.rs` | lid-driven cavity: **known-failing / `#[ignore]`d** (primary-vortex convergence issue; run `cargo test --release -- --ignored`) |
| CFD | `tests/flow_past_cylinder.rs` | steady symmetric wake (low Re) and von Kármán shedding (moderate Re) |
| AI | `lib.rs` | second differentiable plant (pendulum) matches analytic Jacobians |
| FEA | `examples/pile_cage_spacer.rs` | end-to-end spacer milestone (CAD→mesh→elasticity) |

## Benchmarks & examples

Long-term performance tracking uses a [`criterion`](https://github.com/japaric/criterion.rs)
harness (replacing the old `eprintln!`-timing examples):

```sh
cargo bench --workspace            # all crates
cargo bench -p tpt-physics-dem    # DEM parallel stepper (>100k particles)
cargo bench -p tpt-physics-solver # CG / preconditioned-CG / GMRES
cargo bench -p tpt-physics-fea    # declarative-spec box solve
cargo bench -p tpt-physics-cfd     # LBM lattice step
```

Runnable examples (see [`docs/GALLERY.md`](docs/GALLERY.md) for the full
index, or run `./scripts/run_gallery.sh`):

```sh
cargo run --release --example beam          -p tpt-physics-fea   # cantilever tip deflection
cargo run --release --example cavity        -p tpt-physics-cfd   # lid-driven cavity
cargo run --release --example granular_pile -p tpt-physics-dem   # settling pile
cargo run --release --example rl_pendulum  -p tpt-physics-ai    # differentiable env + Jacobians
cargo run --release --example gradient_opt -p tpt-physics-ai    # gradient-based control opt
```

The DEM `rayon` stepper ([`World::step_par`]) is the CPU-acceleration path for
large counts; the hardware-dispatch API in `tpt-physics-solver` selects the GPU
target for problems above its size threshold. With `--features gpu` a real WGSL
`matvec` compute kernel (`tpt-physics-solver::gpu`) runs on the first available
adapter; without the feature, or with no GPU present, it returns
`BackendUnavailable` and callers fall back to the CPU path.

## Troubleshooting

- **`failed to load source for dependency tpt-math-...`** — the sibling repos
  are missing. Run `./scripts/bootstrap.sh` (or `scripts/bootstrap.ps1`) and
  clone `tpt-math` / `tpt-fem` as siblings of this directory.
- **`lid_driven_cavity` test "fails"** — it is intentionally `#[ignore]`d
  (experimental). Run it explicitly with
  `cargo test --release -- --ignored` if you want to inspect it.

## License

Dual-licensed under [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE-2.0) at
your option. Copyright TPT Solutions.
