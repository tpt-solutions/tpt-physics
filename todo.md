# tpt-physics — Project Checklist

TPT Solutions | Dual-licensed MIT / Apache-2.0

## Phase 0: Project Setup & Tooling
- [x] Initialize git repository
- [x] Add .gitignore (Rust workspace)
- [ ] Add dual license files (LICENSE-MIT, LICENSE-APACHE-2.0; copyright TPT Solutions)
- [ ] Set `license = "MIT OR Apache-2.0"` across all crate manifests
- [ ] Scaffold Cargo workspace with member crates: tpt-physics-core, tpt-physics-solver,
      tpt-physics-fea, tpt-physics-dem, tpt-physics-cfd, tpt-physics-ai (optional)
- [ ] Set up cargo-deny pipeline rejecting copyleft deps (GPL/LGPL/AGPL/SSPL)
- [ ] Add approved dependency stack (faer, nalgebra, uom, rayon, wgpu, serde, bincode,
      tracing, proptest), vetted via cargo-deny
- [ ] Add root README (executive summary, architecture overview)

## Phase 1: Foundation & FEA MVP (Months 1-3)

### tpt-physics-core
- [ ] Unit system: wrapper around `uom` enforcing strict SI unit typing
- [ ] Mesh abstraction: topology-aware mesh (Nodes, Elements, Faces) ingesting
      tpt-cad/biocad output without heavy STEP/IGES parsers
- [ ] Material database: type-safe registry (Young's Modulus, Poisson's ratio,
      density) with JSON/serde serialization

### tpt-physics-solver
- [ ] Sparse linear algebra on `faer`/`sprs`
- [ ] Conjugate Gradient (CG) solver
- [ ] GMRES solver
- [ ] Time integration: Newmark-beta scheme
- [ ] Time integration: Runge-Kutta scheme
- [ ] Hardware dispatch API: route matrix ops to CPU (rayon) or GPU (wgpu/spark)
      by problem size, with zero code changes for the end user

### tpt-physics-fea
- [ ] Element types: linear and quadratic tetrahedrons
- [ ] Element types: hexahedrons
- [ ] Element types: beam/shell elements
- [ ] Static linear stress/strain analysis
- [ ] Modal analysis (natural frequencies)
- [ ] Basic non-linear large-deformation analysis
- [ ] Boundary conditions: fixed supports (builder API)
- [ ] Boundary conditions: point loads
- [ ] Boundary conditions: pressure loads
- [ ] Boundary conditions: thermal gradients

### Milestone
- [ ] Successfully simulate the 3D-printed pile cage spacer

## Phase 2: Granular Physics & Performance (Months 4-6)

### tpt-physics-dem
- [ ] Hertz-Mindlin contact model
- [ ] Friction and damping
- [ ] Spatial hashing (broad-phase collision)
- [ ] SIMD-accelerated narrow-phase contact resolution
- [ ] Validate: wet concrete aggregate flow through a pile cage
- [ ] Validate: soil-structure interaction around a 3D-printed spacer

### Performance
- [ ] GPU acceleration via spark/wgpu for large DEM particle counts (>100k particles)

### tpt-physics-fea (extensions)
- [ ] Non-linear FEA: large deformations
- [ ] Non-linear FEA: plasticity

## Phase 3: Fluid Dynamics & AI Integration (Months 7-12)

### tpt-physics-cfd
- [ ] Lattice Boltzmann Method (LBM) solver for incompressible flow

### tpt-physics-ai
- [ ] Differentiable physics wrappers exposing simulation state as a
      Gymnasium-like environment for tpt-anima reinforcement-learning agents

### Publishing
- [ ] Comprehensive documentation
- [ ] Benchmarks
- [ ] "Spacer Benchmark" case study
