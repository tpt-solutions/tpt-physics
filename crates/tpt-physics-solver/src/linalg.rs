//! Linear-operator abstraction and dense vector helpers shared by the
//! iterative solvers.
//!
//! The solvers operate on a [`LinearOperator`] trait rather than a concrete
//! matrix type, so they work equally well on a [`tpt_fem_sparse::Csr`] (the
//! sparse assembly format reused from `tpt-fem`), a dense matrix, or a
//! matrix-free operator (e.g. a fast-multiple or FFT-based action).

use tpt_fem_sparse::Csr;

/// A linear operator `y = A x` acting on `f64` vectors.
pub trait LinearOperator {
    /// Number of rows of `A`.
    fn nrows(&self) -> usize;
    /// Number of columns of `A`.
    fn ncols(&self) -> usize;
    /// Compute `y = A x` (overwrites `y`).
    fn apply(&self, x: &[f64], y: &mut [f64]);
}

impl LinearOperator for Csr {
    fn nrows(&self) -> usize {
        self.nrows
    }
    fn ncols(&self) -> usize {
        self.ncols
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        debug_assert_eq!(x.len(), self.ncols);
        debug_assert_eq!(y.len(), self.nrows);
        for (r, y_r) in y.iter_mut().enumerate() {
            let mut s = 0.0;
            for idx in self.row_ptrs[r]..self.row_ptrs[r + 1] {
                s += self.values[idx] * x[self.col_ind[idx]];
            }
            *y_r = s;
        }
    }
}

/// A dense column-major-free dense operator backed by a row-major `n×m` buffer.
///
/// Used mainly for small test/teaching systems and for the effective matrices
/// of the time integrators.
#[derive(Debug, Clone)]
pub struct DenseMat {
    /// Row count.
    pub nrows: usize,
    /// Column count.
    pub ncols: usize,
    /// Row-major entries.
    pub data: Vec<f64>,
}

impl DenseMat {
    /// Build a dense matrix from a row-major buffer of length `nrows*ncols`.
    pub fn from_row_major(nrows: usize, ncols: usize, data: Vec<f64>) -> Self {
        assert_eq!(data.len(), nrows * ncols);
        DenseMat { nrows, ncols, data }
    }

    /// `A[i][j]`.
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.ncols + j]
    }

    /// Mutable `A[i][j]`.
    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        self.data[i * self.ncols + j] = v;
    }
}

impl LinearOperator for DenseMat {
    fn nrows(&self) -> usize {
        self.nrows
    }
    fn ncols(&self) -> usize {
        self.ncols
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        for (i, y_i) in y.iter_mut().enumerate() {
            let base = i * self.ncols;
            let s: f64 = self.data[base..base + self.ncols]
                .iter()
                .zip(x)
                .map(|(a, b)| a * b)
                .sum();
            *y_i = s;
        }
    }
}

/// `a = alpha * x + y` (vector saxpy into `a`).
pub fn axpy(alpha: f64, x: &[f64], y: &[f64], a: &mut [f64]) {
    for i in 0..a.len() {
        a[i] = alpha * x[i] + y[i];
    }
}

/// `y = alpha * x + y` (in place).
pub fn axpy_inplace(alpha: f64, x: &[f64], y: &mut [f64]) {
    for i in 0..y.len() {
        y[i] += alpha * x[i];
    }
}

/// Dot product `Σ x_i y_i`.
pub fn dot(x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(y).map(|(a, b)| a * b).sum()
}

/// Euclidean (`L²`) norm.
pub fn norm2(x: &[f64]) -> f64 {
    x.iter().map(|v| v * v).sum::<f64>().sqrt()
}

/// Build a `Csr` from a dense `(rows, cols, data)` triple, dropping zeros.
pub fn csr_from_dense(nrows: usize, ncols: usize, data: &[f64]) -> Csr {
    let mut coo_rows = Vec::new();
    let mut coo_cols = Vec::new();
    let mut coo_vals = Vec::new();
    for i in 0..nrows {
        for j in 0..ncols {
            let v = data[i * ncols + j];
            if v != 0.0 {
                coo_rows.push(i);
                coo_cols.push(j);
                coo_vals.push(v);
            }
        }
    }
    let mut coo = tpt_fem_sparse::Coo::new();
    for k in 0..coo_rows.len() {
        coo.push(coo_rows[k], coo_cols[k], coo_vals[k]);
    }
    coo.to_csr()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn csr_matvec_matches_dense() {
        // [[4,-1,0],[-1,4,-1],[0,-1,4]]
        let csr = csr_from_dense(3, 3, &[4.0, -1.0, 0.0, -1.0, 4.0, -1.0, 0.0, -1.0, 4.0]);
        let x = [1.0, 2.0, 3.0];
        let mut y = [0.0; 3];
        csr.apply(&x, &mut y);
        assert_eq!(y, [2.0, 4.0, 10.0]);

        let dense =
            DenseMat::from_row_major(3, 3, vec![4.0, -1.0, 0.0, -1.0, 4.0, -1.0, 0.0, -1.0, 4.0]);
        let mut yd = [0.0; 3];
        dense.apply(&x, &mut yd);
        assert_eq!(y, yd);
    }

    #[test]
    fn vector_helpers() {
        let x = [1.0, 2.0, 3.0];
        let y = [4.0, -1.0, 0.5];
        assert!((dot(&x, &y) - (4.0 - 2.0 + 1.5)).abs() < 1e-12);
        assert!((norm2(&x) - 14.0_f64.sqrt()).abs() < 1e-12);
        let mut a = [0.0; 3];
        axpy(2.0, &x, &y, &mut a);
        assert_eq!(a, [6.0, 3.0, 6.5]);
    }
}
