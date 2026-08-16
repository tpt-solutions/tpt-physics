# tpt-physics-fea

Finite Element Analysis for `tpt-physics`.

Linear stress/strain, boundary conditions, steady-state heat conduction, modal
analysis, and the Newton–Raphson / arc-length nonlinear solver are **reused
directly** from the sibling `tpt-fem` crates (no wrapper). What is net-new and
lives only here:

## Net-new modules

| Module | Description |
| --- | --- |
| `elements` | Quadratic tetrahedron `tet10_stiffness` (deferred in `tpt-fem-element`), 3-D beam/frame `beam3d_global_stiffness` (torsion, biaxial bending, orientation triad), and Mindlin–Reissner plate/shell `shell4_stiffness`. |
| `nonlinear` | Total-Lagrangian St. Venant–Kirchhoff geometric-nonlinearity framework with consistent internal-force (`tet4_internal_force`) and tangent (`tet4_tangent`) operators for the continuum. |
| `thermal` | Thermal-to-structural coupling: a thermal-strain load vector on top of `tpt-fem-thermal`'s temperature field. |
| `plasticity` | von Mises (J2) return-mapping with linear isotropic hardening — `return_map` (stress-driven) and `update` (strain-driven), matching the Voigt convention of `nonlinear`. |
| `assemble` | Element-to-global scatter helper. |

The directly-reused solvers are reachable as `tpt_physics_fea::elasticity::*`,
`::assembly::*`, `::eigen::*`, `::solve::*`, `::thermal_solver::*`.

## Validations

- `tests/cooks_membrane.rs` — elastic–plastic Cook's membrane: plasticity
  increases compliance vs. the elastic limit; hardening stiffens the response.
- `tests/spacer_milestone.rs` — end-to-end 3D-printed pile-cage spacer:
  CAD→mesh→`tpt-fem-elasticity::solve_elasticity` (Continuum3D), fixed base
  stays put, free top compresses under self-weight.

## License

Dual-licensed under [MIT](../../LICENSE-MIT) and [Apache-2.0](../../LICENSE-APACHE-2.0).
Copyright TPT Solutions.
