//! Errors produced by the iterative and time-integration solvers.

use std::fmt;

/// Solver failures.
#[derive(Debug)]
pub enum SolverError {
    /// The iteration did not reach the requested tolerance within the
    /// iteration budget.
    NotConverged {
        /// Iterations performed.
        iterations: usize,
        /// Final relative residual.
        residual: f64,
    },
    /// The problem was not square (`nrows != ncols`).
    NotSquare {
        /// Rows.
        nrows: usize,
        /// Columns.
        ncols: usize,
    },
    /// A zero or non-finite value was encountered where a positive quantity
    /// was required (e.g. a preconditioner solve).
    Numerical(String),
    /// A requested backend (e.g. GPU) is not available in this build.
    BackendUnavailable(String),
}

impl fmt::Display for SolverError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SolverError::NotConverged {
                iterations,
                residual,
            } => write!(
                f,
                "solver did not converge in {iterations} iterations (residual {residual:.3e})"
            ),
            SolverError::NotSquare { nrows, ncols } => {
                write!(f, "operator is not square: {nrows} x {ncols}")
            }
            SolverError::Numerical(m) => write!(f, "numerical failure: {m}"),
            SolverError::BackendUnavailable(m) => write!(f, "backend unavailable: {m}"),
        }
    }
}

impl std::error::Error for SolverError {}

/// Result of a converged (or stopped) iterative solve.
#[derive(Debug, Clone)]
pub struct SolveReport {
    /// Iterations performed.
    pub iterations: usize,
    /// Final relative residual `||r|| / ||b||`.
    pub residual: f64,
    /// Whether the requested tolerance was reached.
    pub converged: bool,
}
