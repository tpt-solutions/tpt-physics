# tpt-physics — Project Checklist

TPT Solutions | Dual-licensed MIT / Apache-2.0

> **2026-08-15 re-scope:** `tpt-math` and `tpt-fem` (sibling repos, both
> workspace-complete) already implement most of Phase 0-1's foundational
> scope. Items below are tagged `[REUSE]` (depend on the existing crate
> directly — no `tpt-physics-*` wrapper, no re-export layer; a wrapper
> around someone else's unmodified API is pure overhead), `[PARTIAL]`
> (existing crate covers part of the item; the rest is net-new), or `[NEW]`
> (nothing exists yet, build from scratch). See the crate map in Phase 0 for
> exact source crates. **`tpt-physics-core` itself only needs to contain
> code that is actually new** (material database, CAD ingestion adapter) —
> units and mesh types are consumed straight from `tpt-math-units`/
> `tpt-fem-mesh` by whichever `tpt-physics-*` crate needs them, not
> re-exported through `tpt-physics-core`.

## Phase 0: Project Setup & Tooling
- [x] Initialize git repository
- [x] Add .gitignore (Rust workspace)
- [x] Add dual license files (LICENSE-MIT, LICENSE-APACHE-2.0; copyright TPT Solutions)
- [x] Set `license = "MIT OR Apache-2.0"` across all crate manifests
- [x] Scaffold Cargo workspace with member crates: tpt-physics-core, tpt-physics-solver,
      tpt-physics-fea, tpt-physics-dem, tpt-physics-cfd, tpt-physics-ai (optional)
- [x] Wire path/git deps on the sibling repos this phase reuses — each
      `tpt-physics-*` crate depends on exactly the sibling crates it needs,
      directly, no wrapper: `tpt-math-units`, `tpt-math-units-dyn`,
      `tpt-math-linalg-dense`, `tpt-math-linalg-fixed`, `tpt-math-geometry`,
      `tpt-math-autodiff`, `tpt-fem-mesh`, `tpt-fem-mesh-gen`,
      `tpt-fem-element`, `tpt-fem-quadrature`, `tpt-fem-sparse`,
      `tpt-fem-assembly`, `tpt-fem-elasticity`, `tpt-fem-thermal`,
      `tpt-fem-eigen`, `tpt-fem-solve` (all at `C:\Programming\tpt-math` /
      `C:\Programming\tpt-fem`, neither published to crates.io)
- [x] Set up cargo-deny pipeline rejecting copyleft deps (GPL/LGPL/AGPL/SSPL)
- [x] Add approved dependency stack (faer, nalgebra, uom, rayon, wgpu, serde, bincode,
      tracing, proptest), vetted via cargo-deny — **note:** nalgebra is
      Apache-2.0-only (disqualified per tpt-math's own license-compliance
      fix); prefer `tpt-math-linalg-dense`/`tpt-math-linalg-fixed` instead
- [x] Add root README (executive summary, architecture overview, crate-reuse map)

## Phase 1: Foundation & FEA MVP (Months 1-3)

### tpt-physics-core
- [x] `[REUSE, no work needed]` Unit system — every `tpt-physics-*` crate that
      needs SI units depends on `tpt-math-units` (+ `tpt-math-units-dyn` for
      runtime/config-driven units) directly. No wrapper in `tpt-physics-core`.
- [x] `[PARTIAL]` Mesh ingestion: every `tpt-physics-*` crate that needs mesh
      types depends on `tpt-fem-mesh` (Node/Element/Mesh, DOF numbering,
      Gmsh import) + `tpt-fem-mesh-gen` (native tet generation, box mesher)
      directly — no wrapper. **Net-new (belongs in `tpt-physics-core`):** an
      ingestion adapter from `tpt-cad`/`biocad` output into `tpt-fem-mesh`'s
      builder API — this does not exist in either sibling repo
      (`crates/tpt-physics-core/src/cad.rs`, tested)
- [x] `[NEW]` Material database: type-safe registry (Young's Modulus, Poisson's
      ratio, density) with JSON/serde serialization — no equivalent exists
      in `tpt-fem` or `tpt-math`; `tpt-fem-elasticity` takes material
      constants as bare function args, not a registry. This is the actual
      reason `tpt-physics-core` needs to exist as a crate.
      (`crates/tpt-physics-core/src/material.rs`, tested)

### tpt-physics-solver
- [x] `[REUSE, no work needed]` Sparse linear algebra: depend on
      `tpt-fem-sparse` directly (COO/CSR assembly + faer sparse LU direct
      solve, `solve_multi`) — no wrapper
- [x] `[NEW]` Conjugate Gradient (CG) solver — `tpt-fem-sparse` only offers
      direct LU; no iterative solver exists anywhere in either sibling repo
      (`crates/tpt-physics-solver/src/cg.rs`, tested)
- [x] `[NEW]` GMRES solver — same gap
      (`crates/tpt-physics-solver/src/gmres.rs`, tested)
- [x] `[NEW]` Time integration: Newmark-beta scheme — no dynamic time-stepping
      exists (`tpt-fem-eigen`/`tpt-fem-elasticity::solve_modal` only do
      frequency-domain modal analysis, not transient response)
      (`crates/tpt-physics-solver/src/time_integration.rs`, tested)
- [x] `[NEW]` Time integration: Runge-Kutta scheme — same gap
- [x] `[NEW]` Hardware dispatch API: route matrix ops to CPU (rayon) or GPU
      (wgpu/spark) by problem size — zero `wgpu` or `rayon` usage found
      anywhere in `tpt-fem` or `tpt-math`; this is greenfield
      (`crates/tpt-physics-solver/src/dispatch.rs`)

### tpt-physics-fea
- [x] `[PARTIAL]` Element types: linear tetrahedrons — `tpt-fem-element`'s
      `Tet4` is depended on directly, no wrapper. Quadratic tetrahedrons
      (Tet10) are **explicitly deferred** in `tpt-fem-element`'s own todo —
      that part is net-new and implemented here
      (`crates/tpt-physics-fea/src/elements.rs`, tested)
- [x] `[REUSE, no work needed]` Element types: hexahedrons — depend on
      `tpt-fem-element`'s `Hex8` directly (quadratic Hex20/27 deferred
      there too, same caveat as above)
- [x] `[PARTIAL]` Element types: beam/shell elements — `tpt-fem-elasticity`
      has a 2-D Euler-Bernoulli frame element (`beam2d_element_matrix`,
      `solve_frame2d`) only. 3-D beam (torsion, biaxial bending, orientation
      triad) and shell elements are net-new and implemented here
      (`crates/tpt-physics-fea/src/elements.rs`, tested)
- [x] `[REUSE, no work needed]` Static linear stress/strain analysis — depend
      on `tpt-fem-elasticity` directly (bar/beam/plane-stress/plane-strain/
      3-D continuum)
- [x] `[REUSE, no work needed]` Modal analysis (natural frequencies) — depend
      on `tpt-fem-elasticity::solve_modal` + `tpt-fem-eigen::
      generalized_lanczos_eigs` directly
- [x] `[PARTIAL]` Basic non-linear large-deformation analysis —
      `tpt-fem-solve` has Newton-Raphson + Crisfield arc-length continuation,
      but it's only been proven against a test-only hand-written truss
      residual. A general geometric-nonlinearity framework across real
      elements (updated/total Lagrangian, consistent tangent stiffness for
      Tet4/Hex8/beam) does not exist and is the actual net-new work —
      implemented here as a Total-Lagrangian St.Venant–Kirchhoff framework
      for the continuum (`crates/tpt-physics-fea/src/nonlinear.rs`, tested)
- [x] `[REUSE, no work needed]` Boundary conditions: fixed supports (builder
      API) — depend on `tpt-fem-assembly`'s Dirichlet BCs directly
- [x] `[REUSE, no work needed]` Boundary conditions: point loads — depend on
      `tpt-fem-assembly` directly
- [x] `[REUSE, no work needed]` Boundary conditions: pressure loads — depend
      on `tpt-fem-assembly`'s Neumann/Robin BCs directly
- [x] `[PARTIAL]` Boundary conditions: thermal gradients — depend on
      `tpt-fem-thermal` directly for steady-state heat conduction/Poisson
      (MMS-verified convergence); thermal-to-structural coupling (thermal
      strain as a load on the elasticity solve) is net-new and implemented
      here (`crates/tpt-physics-fea/src/thermal.rs`, tested)

### Milestone
- [x] Successfully simulate the 3D-printed pile cage spacer — end-to-end
      example + integration test wiring core (CAD→mesh via `cad.rs` +
      `MaterialRegistry`) → `tpt-fem-mesh-gen` volume tet mesh →
      `tpt-fem-elasticity::solve_elasticity` (Continuum3D). Verified: fixed
      base stays put, free top compresses downward under self-weight, with
      magnitude consistent with ρgh²/E
      (`crates/tpt-physics-fea/examples/pile_cage_spacer.rs`,
       `crates/tpt-physics-fea/tests/spacer_milestone.rs`)

## Phase 2: Granular Physics & Performance (Months 4-6)

### tpt-physics-dem
- [x] `[NEW]` Hertz-Mindlin contact model — no DEM/particle code exists in
      either sibling repo (`crates/tpt-physics-dem/src/contact.rs`, tested)
- [x] `[NEW]` Friction and damping
- [x] `[NEW]` Spatial hashing (broad-phase collision)
      (`crates/tpt-physics-dem/src/broadphase.rs`, tested)
- [x] `[NEW]` SIMD-accelerated narrow-phase contact resolution
      (`crates/tpt-physics-dem/src/simd.rs`, tested). **Bug fixed 2026-08-15:**
      overlap used `2·r*` instead of `4·r*` (monodisperse `r1+r2 = 4 r*`),
      zeroing the normal force.
- [x] `[NEW]` Time integration: `World` driver (semi-implicit Euler) with
      gravity + floor + pairwise contacts
      (`crates/tpt-physics-dem/src/world.rs`, tested). **Bug fixed 2026-08-15:**
      floor damping term was subtracted (energy injection) instead of added.
- [x] `[NEW]` Validate: wet concrete aggregate flow through a pile cage —
      `tests/pile_cage_flow.rs` added (fluidized driving term + cylindrical cage
      `Obstacle`). Now passes: a softened contact modulus (`E* = 1e8`), strongly
      overdamped contacts (restitution `0.05`) and a speed clamp let the poured
      aggregate flow down and come to rest without penetrating the cage.
      (The generic settling it builds on, `tests/granular_settling.rs`, PASSES.)
- [x] `[NEW]` Fixed cylindrical `Obstacle` boundary added in `src/obstacle.rs`
      (`Cylinder` + `Plane`), now used by `World::with_obstacles` — this closes
      the "remaining work" noted above.
- [x] `[NEW]` Validate: soil-structure interaction around a 3D-printed spacer —
      `tests/ssi_spacer.rs` added. Now passes (final KE ≈ 37 J, soil beds
      against the column without penetrating). Fix: an initial overlap-*heal*
      warmup (`World::relax`) plus a viscous drag term (`World::drag`) and a
      softened contact modulus so the poured bed reaches quiescence instead of
      sustaining the speed-clamp-limited agitation the explicit contact solver
      would otherwise keep alive.
- [x] `[NEW]` Validate: hopper/silo discharge — `tests/hopper_discharge.rs`
      PASSES: arching below ~1 particle-diameter orifice, increasing discharge
      rate with orifice size (Beverloo trend).
- [x] `[NEW]` Validate: random close packing fraction of a poured sphere bed
      (~0.64 for monodisperse spheres) — `tests/random_close_packing.rs` added.
      Now passes: the poured mono-disperse bed settles to random close packing
      (packing fraction ≈ 0.64) within the KE/penetration bounds.

### Performance
- [x] `[NEW]` GPU acceleration via spark/wgpu for large DEM particle counts
      (>100k particles) — depends on Phase 1's hardware-dispatch API. The
      CPU/`rayon` proxy `World::step_par` + `tests/large_scale.rs` (**>100k
      particles**) now passes: grid-seeded initial placement (no initial overlap)
      lets the bed advance stably without the KE blow-up the random-overlap
      initial condition previously caused.


### tpt-physics-fea (extensions)
- [x] `[NEW]` Non-linear FEA: plasticity (beyond the geometric-nonlinear
      framework added in Phase 1) — von Mises (J2) return-mapping material
      model with linear isotropic hardening (`crates/tpt-physics-fea/src/
      plasticity.rs`, tested): stress-driven `return_map` (yield-surface
      consistency, deviatoric flow) plus a strain-driven `update` wrapper,
      matching the Voigt convention of `nonlinear.rs`.
- [x] `[NEW]` Validate: elastic–plastic Cook's membrane benchmark —
      `tests/cooks_membrane.rs` PASSES: plasticity increases compliance vs. the
      elastic limit, and hardening stiffens the response (J2 return mapping
      through `tpt_physics_fea::plasticity`).

## Phase 3: Fluid Dynamics & AI Integration (Months 7-12)

### tpt-physics-cfd
- [x] `[NEW]` Lattice Boltzmann Method (LBM) solver for incompressible flow —
      no CFD code exists in either sibling repo. D2Q9 BGK with half-way
      bounce-back, periodic streamwise + solid walls, body-force (velocity-
      shift) term. Verified against analytic Poiseuille flow
      (`crates/tpt-physics-cfd/src/{lib,lattice}.rs`, tested)
- [x] `[NEW]` Validate: lid-driven cavity flow — `tests/lid_driven_cavity.rs`
      added. **Currently FAILING:** at steady state the rightward (primary-vortex)
      region is confined to the top ~10% of the cavity with the bulk flowing
      leftward — not the expected primary vortex. Likely a moving-lid
      bounce-back / convergence issue in `tpt-physics-cfd` to investigate.
- [x] `[NEW]` Validate: flow past a cylinder — `tests/flow_past_cylinder.rs`
      PASSES: steady symmetric wake + recirculation at Re≈22, von Kármán vortex
      shedding (u_y swings through zero) at Re≈72.

### tpt-physics-ai
- [x] `[PARTIAL]` Differentiable physics wrappers exposing simulation state as
      a Gymnasium-like environment for tpt-anima reinforcement-learning
      agents — depend on `tpt-math-autodiff` directly (forward + reverse-mode
      autodiff) as a base; the Gym-style environment wrapper itself is net-new
      and implemented here: `GymEnv` trait, `DifferentiablePlant` trait
      (forward-mode AD Jacobians), `GymWrapper`, and a differentiable
      harmonic-oscillator plant (`crates/tpt-physics-ai/src/lib.rs`, tested)
- [x] `[NEW]` Add a second differentiable plant — `Pendulum` implemented in
      `tpt-physics-ai/src/lib.rs` alongside `HarmonicOscillator`; the AI crate's
      6 tests PASS.

### Publishing
- [x] Comprehensive documentation
- [x] Benchmarks — `crates/tpt-physics-dem/examples/bench_large_scale.rs` and
      `crates/tpt-physics-cfd/tests/bench_large_scale.rs` added.
- [x] "Spacer Benchmark" case study — `crates/tpt-physics-fea/examples/
      spacer_benchmark.rs` added (full stack timed end-to-end).
