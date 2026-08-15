//! Finite Element Analysis for `tpt-physics`.
//!
//! This crate reuses the sibling `tpt-fem` crates **directly** for the
//! capabilities they already provide: linear stress/strain
//! ([`elasticity`]), boundary conditions ([`assembly`]), steady-state
//! heat conduction ([`thermal`]), modal/eigen analysis ([`eigen`]), and the
//! Newton–Raphson / arc-length nonlinear solver ([`solve`]).
//!
//! What is net-new here (and lives only in this crate):
//!
//! * [`elements`] — quadratic tetrahedron [`elements::tet10_stiffness`],
//!   3-D beam/frame [`elements::beam3d_global_stiffness`], and a Mindlin–
//!   Reissner plate/shell [`elements::shell4_stiffness`];
//! * [`nonlinear`] — a Total-Lagrangian, St. Venant–Kirchhoff geometric-
//!   nonlinearity framework with consistent internal-force and tangent
//!   operators for the continuum [`nonlinear::tet4_internal_force`] /
//!   [`nonlinear::tet4_tangent`];
//! * [`thermal`] — thermal-to-structural coupling (thermal-strain load
//!   vector) on top of `tpt-fem-thermal`'s temperature field;
//! * [`assemble`] — element-to-global scatter helper.

pub mod assemble;
pub mod elements;
pub mod nonlinear;
pub mod plasticity;
pub mod thermal;

// Directly-reused sibling crates (no wrapper, no re-export shim beyond this
// namespace). Consumers reach `tpt-fem-elasticity`'s solver through
// `tpt_physics_fea::elasticity::*`.
pub use tpt_fem_assembly as assembly;
pub use tpt_fem_eigen as eigen;
pub use tpt_fem_elasticity as elasticity;
pub use tpt_fem_solve as solve;
pub use tpt_fem_thermal as thermal_solver;
