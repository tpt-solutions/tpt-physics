# tpt-phys-cfd

Computational Fluid Dynamics for `tpt-physics`, via a pure-Rust, data-oriented
Lattice Boltzmann Method (LBM).

No CFD code exists in the sibling `tpt-math` / `tpt-fem` workspaces, so this is
net-new.

## Solver

`Lbm2D` is a D2Q9 (9-velocity, 2-D) BGK lattice with:

- half-way bounce-back solid boundaries (stationary or moving lids);
- periodic, or velocity-inlet / zero-gradient-outlet `x` boundaries;
- a body-force (Guo-style velocity-shift) term for pressure-gradient /
  gravity-driven flow;
- a circular obstacle primitive for bluff-body flows;
- `rayon`-friendly bulk loops (`step_par`).

## Modules

| Module | Description |
| --- | --- |
| `lattice` | `D2Q9` definition: velocity set, weights, equilibrium distribution. |
| `lib` | `Lbm2D` solver, streaming/bounce-back, macroscopic field recovery. |

## Validations

- `tests/lid_driven_cavity.rs` — lid-driven cavity: primary recirculating
  vortex (lid drags fluid `+x` near the top, return flow near the floor).
- `tests/flow_past_cylinder.rs` — steady symmetric wake + recirculation at
  Re≈22, von Kármán vortex shedding at Re≈72.

## License

Dual-licensed under [MIT](../../LICENSE-MIT) and [Apache-2.0](../../LICENSE-APACHE-2.0).
Copyright TPT Solutions.
