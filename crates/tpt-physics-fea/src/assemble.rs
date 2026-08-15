//! Scatter helpers for assembling global sparse matrices from element matrices.

use tpt_fem_sparse::Coo;

/// Add an `n×n` element stiffness `ke` (row-major) into `coo` at the global DOF
/// indices given by `dofs` (length `n`).
pub fn scatter(coo: &mut Coo, dofs: &[usize], ke: &[f64], n: usize) {
    debug_assert_eq!(ke.len(), n * n);
    for i in 0..n {
        for j in 0..n {
            coo.push(dofs[i], dofs[j], ke[i * n + j]);
        }
    }
}
