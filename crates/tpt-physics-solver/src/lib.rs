//! Iterative linear solvers, time-integration schemes, and CPU/GPU hardware
//! dispatch for `tpt-physics`.
//!
//! The sibling `tpt-fem-sparse` crate provides only a direct LU solver. This
//! crate fills the iterative-solver and transient-simulation gap that the
//! re-scope identified as net-new:
//!
//! * [`cg`] / [`cg_pc`] — (preconditioned) Conjugate Gradient
//! * [`gmres`] / [`gmres_pc`] — restarted GMRES for non-symmetric systems.
//!   Preconditioned GMRES ([`gmres_pc`]) is **experimental** — it is correct for
//!   `M = I` and for a Jacobi (diagonal) preconditioner, but a multigrid/AMG
//!   cycle has not been benchmarked yet (see [`multigrid`]).
//! * [`rk4`] / [`NewmarkBeta`] — Runge–Kutta and Newmark-beta time integration
//! * [`HardwareDispatch`] — route kernels to CPU (`rayon`) or GPU (`wgpu`).
//!   **The GPU backend is experimental** — it is enabled by the `gpu` feature
//!   and runs a real WGSL `matvec` compute kernel on the first available
//!   adapter (see [`crate::gpu`]). Without the feature, or with no GPU adapter
//!   present, the GPU target returns `BackendUnavailable` and callers can fall
//!   back to the CPU path.
//!
//! All solvers operate on the [`LinearOperator`] trait, so they work on
//! [`tpt_fem_sparse::Csr`] (reused from `tpt-fem`), dense matrices, or
//! matrix-free operators alike.

pub mod cg;
pub mod dispatch;
pub mod error;
pub mod gmres;
pub mod linalg;
pub mod multigrid;
pub mod time_integration;

#[cfg(feature = "gpu")]
pub mod gpu;

pub use cg::{cg, cg_pc};
pub use dispatch::{parallel_matvec_csr, ComputeTarget, HardwareDispatch};
pub use error::{SolveReport, SolverError};
pub use gmres::gmres;
pub use linalg::{axpy, axpy_inplace, csr_from_dense, dot, norm2, DenseMat, LinearOperator};
pub use multigrid::{Grid2D, Multigrid};
pub use time_integration::{rk4, NewmarkBeta};
