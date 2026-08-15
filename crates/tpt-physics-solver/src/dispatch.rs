//! Hardware dispatch: route matrix operations to multi-threaded CPU (`rayon`)
//! or GPU (`wgpu`/`spark`) by problem size, with a single API for callers.
//!
//! This is greenfield — neither `tpt-fem` nor `tpt-math` uses `rayon`/`wgpu`
//! for matrix work. The CPU path uses `rayon` parallel row-loop matvecs. The
//! GPU path is declared and selected, but the actual `wgpu`/`spark` kernel
//! backend is a tracked follow-up; until it is linked this returns
//! [`SolverError::BackendUnavailable`] rather than silently falling back.

use crate::error::SolverError;
use rayon::prelude::*;
use tpt_fem_sparse::Csr;

/// Where a compute kernel should run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ComputeTarget {
    /// Multi-threaded CPU via `rayon`.
    Cpu,
    /// GPU via `wgpu`/`spark` (not yet linked).
    Gpu,
}

impl std::fmt::Display for ComputeTarget {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ComputeTarget::Cpu => write!(f, "cpu"),
            ComputeTarget::Gpu => write!(f, "gpu"),
        }
    }
}

/// Decides and performs the routing of linear-algebra kernels.
#[derive(Debug, Clone)]
pub struct HardwareDispatch {
    /// Preferred target when the problem is below `gpu_size_threshold`.
    pub preference: ComputeTarget,
    /// Problems with `n ≥ gpu_size_threshold` are routed to the GPU.
    pub gpu_size_threshold: usize,
}

impl Default for HardwareDispatch {
    fn default() -> Self {
        HardwareDispatch {
            preference: ComputeTarget::Cpu,
            gpu_size_threshold: 50_000,
        }
    }
}

impl HardwareDispatch {
    /// A CPU-first dispatcher that promotes large problems (≥ `threshold`
    /// degrees of freedom) to the GPU.
    pub fn new(gpu_size_threshold: usize) -> Self {
        HardwareDispatch {
            preference: ComputeTarget::Cpu,
            gpu_size_threshold,
        }
    }

    /// Choose the target for a problem of `problem_size` degrees of freedom.
    pub fn select(&self, problem_size: usize) -> ComputeTarget {
        if problem_size >= self.gpu_size_threshold {
            ComputeTarget::Gpu
        } else {
            self.preference
        }
    }

    /// Sparse matrix-vector product `y = A x` routed by problem size.
    ///
    /// On the CPU this uses a `rayon` parallel row loop. On the GPU it returns
    /// [`SolverError::BackendUnavailable`] until the `wgpu`/`spark` kernel
    /// backend is integrated.
    pub fn matvec(&self, a: &Csr, x: &[f64], y: &mut [f64]) -> Result<(), SolverError> {
        match self.select(a.nrows) {
            ComputeTarget::Cpu => {
                parallel_matvec_csr(a, x, y);
                Ok(())
            }
            ComputeTarget::Gpu => Err(SolverError::BackendUnavailable(
                "GPU (wgpu/spark) kernel backend not linked in this build".into(),
            )),
        }
    }
}

/// Parallel sparse matrix-vector product using `rayon`.
///
/// Each row's dot product is independent, so the row loop is distributed across
/// the global thread pool. This is the CPU work-horse that the hardware-dispatch
/// API routes to.
pub fn parallel_matvec_csr(a: &Csr, x: &[f64], y: &mut [f64]) {
    debug_assert_eq!(x.len(), a.ncols);
    debug_assert_eq!(y.len(), a.nrows);
    y.par_iter_mut().enumerate().for_each(|(r, yr)| {
        let mut s = 0.0;
        for idx in a.row_ptrs[r]..a.row_ptrs[r + 1] {
            s += a.values[idx] * x[a.col_ind[idx]];
        }
        *yr = s;
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::csr_from_dense;
    use crate::linalg::LinearOperator;

    #[test]
    fn dispatch_picks_gpu_for_large() {
        let d = HardwareDispatch::new(100);
        assert_eq!(d.select(10), ComputeTarget::Cpu);
        assert_eq!(d.select(100), ComputeTarget::Gpu);
        assert_eq!(d.select(1_000_000), ComputeTarget::Gpu);
    }

    #[test]
    fn parallel_matvec_matches_sequential() {
        let data = vec![
            4.0, -1.0, 0.0, -1.0, -1.0, 0.0, //
            0.0, 0.0, 4.0, -1.0, -1.0, 0.0, //
            -1.0, 0.0, -1.0, 4.0, 0.0, -1.0, //
            0.0, -1.0, 0.0, -1.0, 4.0, -1.0, //
            -1.0, 0.0, 0.0, 0.0, -1.0, 4.0, //
            0.0, 0.0, -1.0, 0.0, -1.0, 4.0, //
        ];
        let csr = csr_from_dense(6, 6, &data);
        let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let mut y_seq = vec![0.0; 6];
        csr.apply(&x, &mut y_seq);
        let mut y_par = vec![0.0; 6];
        parallel_matvec_csr(&csr, &x, &mut y_par);
        for i in 0..6 {
            assert!((y_seq[i] - y_par[i]).abs() < 1e-12, "row {i}");
        }
    }

    #[test]
    fn gpu_backend_reports_unavailable() {
        let d = HardwareDispatch::new(2);
        let csr = csr_from_dense(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let mut y = [0.0; 2];
        assert!(matches!(
            d.matvec(&csr, &[1.0, 1.0], &mut y),
            Err(SolverError::BackendUnavailable(_))
        ));
    }
}
