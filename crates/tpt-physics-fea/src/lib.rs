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
//!   Reissner plate/shell [`elements::shell4_stiffness`]. **Note:** the
//!   `shell4` element is *experimental* — it passes rigid-body-exactness and a
//!   simply-supported-plate sanity checks, but has not been validated against
//!   a curved-shell benchmark (Scordelis–Lo). The 3-D beam is Euler–Bernoulli
//!   (no shear deformation) and likewise experimental.
//! * [`nonlinear`] — a Total-Lagrangian, St. Venant–Kirchhoff geometric-
//!   nonlinearity framework with consistent internal-force and tangent
//!   operators for the **linear-tetrahedral continuum**
//!   [`nonlinear::tet4_internal_force`] / [`nonlinear::tet4_tangent`]. This is
//!   the only continuum nonlinear element implemented; there is **no** Hex8 or
//!   beam geometric-nonlinear framework yet (those would reuse this same
//!   Total-Lagrangian pattern with their own reference gradients/quadrature).
//! * [`plasticity`] — J2 (von Mises) return-mapping with linear isotropic
//!   hardening, matching the nonlinear Voigt convention. This is the only
//!   plasticity model shipped.
//! * [`thermal`] — thermal-to-structural coupling (thermal-strain load
//!   vector) on top of `tpt-fem-thermal`'s temperature field;
//! * [`assemble`] — element-to-global scatter helper;
//! * [`spec`] — a declarative JSON problem specification (material + domain +
//!   boundary conditions + loads + solver) solved end-to-end through the same
//!   pipeline, reusing [`tpt_physics_core::MaterialRegistry`].
//!
//! **Status:** the JSON spec is a thin façade over the validated linear-elastic
//! solve — it is *experimental* in the sense that only the `static_linear`
//! (continuum 3-D) path is wired; modal / nonlinear problems must still be
//! driven by the direct APIs.

pub mod assemble;
pub mod elements;
pub mod nonlinear;
pub mod plasticity;
pub mod spec;
pub mod thermal;

// Directly-reused sibling crates (no wrapper, no re-export shim beyond this
// namespace). Consumers reach `tpt-fem-elasticity`'s solver through
// `tpt_physics_fea::elasticity::*`.
pub use tpt_fem_assembly as assembly;
pub use tpt_fem_eigen as eigen;
pub use tpt_fem_elasticity as elasticity;
pub use tpt_fem_solve as solve;
pub use tpt_fem_thermal as thermal_solver;
