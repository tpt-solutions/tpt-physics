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
├── tpt-physics-solver    # Iterative solvers (CG, GMRES), time integration, HW dispatch
├── tpt-physics-fea       # FEA elements (Tet10, 3D beam, shell), nonlinear, thermal coupling
├── tpt-physics-dem       # Granular physics: Hertz-Mindlin, spatial hashing, SIMD
├── tpt-physics-cfd       # Lattice Boltzmann (D2Q9) incompressible flow
└── tpt-physics-ai        # Differentiable physics wrappers for RL agents
```

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

Requires a Rust toolchain ≥ 1.84. The sibling `tpt-math` and `tpt-fem` repos
must be checked out as siblings of this directory (the `[workspace.dependencies]`
in `Cargo.toml` reference them by relative path).

## Roadmap

- **Phase 1 — Foundation & FEA MVP:** workspace, solvers, linear/nonlinear
  FEA. Milestone: simulate the 3D-printed pile cage spacer.
- **Phase 2 — Granular & performance:** DEM (concrete aggregate flow),
  GPU acceleration for >100k particles, nonlinear FEA extensions.
- **Phase 3 — Fluids & AI:** LBM CFD, differentiable gym environments,
  documentation, benchmarks, and the "Spacer Benchmark" case study.

See [`todo.md`](todo.md) for the full checklist.
