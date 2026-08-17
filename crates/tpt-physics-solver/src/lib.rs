//! Iterative linear solvers, time-integration schemes, and CPU/GPU hardware
//! dispatch for `tpt-physics`.
//!
//! The sibling `tpt-fem-sparse` crate provides only a direct LU solver. This
//! crate fills the iterative-solver and transient-simulation gap that the
//! re-scope identified as net-new:
//!
//! * [`cg`] / [`cg_pc`] — (preconditioned) Conjugate Gradient
//! * [`gmres`] — restarted GMRES for non-symmetric systems (**preconditioned
//!   GMRES is experimental**; only the unpreconditioned variant is currently
//!   wired, see the crate roadmap)
//! * [`rk4`] / [`NewmarkBeta`] — Runge–Kutta and Newmark-beta time integration
//! * [`HardwareDispatch`] — route kernels to CPU (`rayon`) or GPU (`wgpu`/`spark`).
//!   **The GPU backend is experimental** — it selects the target and runs the
//!   CPU path, but the actual `wgpu`/`spark` kernel returns
//!   `BackendUnavailable` until linked.
//!
//! All solvers operate on the [`LinearOperator`] trait, so they work on
//! [`tpt_fem_sparse::Csr`] (reused from `tpt-fem`), dense matrices, or
//! matrix-free operators alike.

pub mod cg;
pub mod dispatch;
pub mod error;
pub mod gmres;
pub mod linalg;
pub mod time_integration;

pub use cg::{cg, cg_pc};
pub use dispatch::{parallel_matvec_csr, ComputeTarget, HardwareDispatch};
pub use error::{SolveReport, SolverError};
pub use gmres::gmres;
pub use linalg::{axpy, axpy_inplace, csr_from_dense, dot, norm2, DenseMat, LinearOperator};
pub use time_integration::{rk4, NewmarkBeta};
