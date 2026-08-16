# tpt-physics-dem

Discrete Element Method: Hertz–Mindlin contact, uniform-grid spatial hashing,
SIMD narrow-phase resolution, and a minimal time-stepping driver.

No DEM / particle code exists in the sibling `tpt-math` / `tpt-fem` workspaces,
so everything here is net-new.

## Modules

| Module | Description |
| --- | --- |
| `particle` | `Particle` — spherical state (position, velocity, radius, mass). |
| `contact` | Hertz–Mindlin law: `hertz_normal_force`, `contact_force` (normal + damping + Coulomb-capped Mindlin tangential friction), reduced modulus/radius/mass. |
| `broadphase` | `SpatialHash` — uniform-grid broad phase with candidate-pair and neighbour-list queries. |
| `simd` | SIMD-accelerated narrow-phase contact resolution. |
| `obstacle` | Fixed boundaries: `Obstacle::Cylinder` and `Obstacle::Plane` (half-space walls / hopper funnels / rectangular containments). |
| `world` | `World` driver — semi-implicit Euler with gravity, planar floor (inelastic), obstacle boundaries, optional fluidized driving term, and a speed clamp. `World::step_par` is the `rayon` parallel contact sweep. `World::relax` is a position-based overlap-heal warmup (no energy injected) for stable settling after random placement. |

## Validations

- `tests/granular_settling.rs` — generic settling without blow-up.
- `tests/hopper_discharge.rs` — arching below ~1 particle-diameter orifice,
  increasing discharge with orifice size (Beverloo trend).
- `tests/pile_cage_flow.rs` — wet-concrete aggregate flow through a pile cage.
- `tests/ssi_spacer.rs` — soil–structure interaction around an embedded cylinder.
- `tests/random_close_packing.rs` — mono-disperse random close packing (~0.64).
- `tests/large_scale.rs` — >100k-particle bed advanced stably via `step_par`.

## License

Dual-licensed under [MIT](../../LICENSE-MIT) and [Apache-2.0](../../LICENSE-APACHE-2.0).
Copyright TPT Solutions.
