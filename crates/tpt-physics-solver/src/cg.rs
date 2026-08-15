//! Conjugate Gradient (CG) and preconditioned CG for symmetric positive
//! definite systems.
//!
//! `tpt-fem-sparse` (reused from `tpt-fem`) only offers a direct LU factorise-
//! and-solve path; no iterative solver exists in either sibling repo. CG is
//! the workhorse for the large, sparse, SPD stiffness matrices that arise in
//! FEA, so it is implemented here.

use crate::error::{SolveReport, SolverError};
use crate::linalg::{dot, norm2, LinearOperator};

/// Plain (unpreconditioned) Conjugate Gradient.
///
/// Solves `A x = b` for SPD `A`. `x0`, if `None`, starts from the zero vector.
pub fn cg<A: LinearOperator + ?Sized>(
    a: &A,
    b: &[f64],
    x0: Option<&[f64]>,
    tol: f64,
    max_iter: usize,
) -> Result<(Vec<f64>, SolveReport), SolverError> {
    cg_pc(a, b, x0, tol, max_iter, None)
}

/// Preconditioned Conjugate Gradient.
///
/// `apply_p` (if provided) applies a preconditioner `M⁻¹` to a vector
/// in-place: `z = M⁻¹ r`. A common choice is Jacobi (diagonal) scaling.
pub fn cg_pc<A: LinearOperator + ?Sized>(
    a: &A,
    b: &[f64],
    x0: Option<&[f64]>,
    tol: f64,
    max_iter: usize,
    apply_p: Option<&dyn Fn(&[f64], &mut [f64])>,
) -> Result<(Vec<f64>, SolveReport), SolverError> {
    let n = b.len();
    if a.nrows() != n || a.ncols() != n {
        return Err(SolverError::NotSquare {
            nrows: a.nrows(),
            ncols: a.ncols(),
        });
    }
    let bnorm = norm2(b);
    if bnorm == 0.0 {
        return Ok((vec![0.0; n], SolveReport { iterations: 0, residual: 0.0, converged: true }));
    }

    let mut x = match x0 {
        Some(v) => v.to_vec(),
        None => vec![0.0; n],
    };
    let mut r = vec![0.0; n];
    a.apply(&x, &mut r);
    for i in 0..n {
        r[i] = b[i] - r[i];
    }

    let mut z = vec![0.0; n];
    match apply_p {
        Some(p) => p(&r, &mut z),
        None => z.copy_from_slice(&r),
    }

    let mut p = z.clone();
    let mut rz_old = dot(&r, &z);

    let mut residual = norm2(&r) / bnorm;
    if residual < tol {
        return Ok((x, SolveReport { iterations: 0, residual, converged: true }));
    }

    let mut ap = vec![0.0; n];
    let mut iterations = 0;
    for k in 0..max_iter {
        iterations = k + 1;
        a.apply(&p, &mut ap);
        let pap = dot(&p, &ap);
        if pap.abs() < 1e-300 {
            return Err(SolverError::Numerical("pᵀAp vanished".into()));
        }
        let alpha = rz_old / pap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        residual = norm2(&r) / bnorm;
        if residual < tol {
            return Ok((x, SolveReport { iterations, residual, converged: true }));
        }
        let mut z_new = vec![0.0; n];
        match apply_p {
            Some(pfn) => pfn(&r, &mut z_new),
            None => z_new.copy_from_slice(&r),
        }
        let rz_new = dot(&r, &z_new);
        let beta = rz_new / rz_old;
        for i in 0..n {
            p[i] = z_new[i] + beta * p[i];
        }
        rz_old = rz_new;
    }

    Err(SolverError::NotConverged { iterations, residual })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::csr_from_dense;
    use tpt_fem_sparse::Csr;

    fn spd() -> (CsrMock, Vec<f64>) {
        // 2D Poisson 3x3 (Laplacian) is SPD. Use [[4,-1,0,...],...].
        let n = 3;
        let mut data = vec![0.0; n * n];
        for i in 0..n {
            data[i * n + i] = 4.0;
            if i > 0 {
                data[i * n + (i - 1)] = -1.0;
                data[(i - 1) * n + i] = -1.0;
            }
        }
        let csr = csr_from_dense(n, n, &data);
        (CsrMock(csr), vec![3.0, 2.0, 3.0])
    }

    struct CsrMock(Csr);
    impl LinearOperator for CsrMock {
        fn nrows(&self) -> usize { self.0.nrows }
        fn ncols(&self) -> usize { self.0.ncols }
        fn apply(&self, x: &[f64], y: &mut [f64]) { self.0.apply(x, y) }
    }

    #[test]
    fn cg_solves_spd() {
        let (a, b) = spd();
        let (x, rep) = cg(&a, &b, None, 1e-10, 100).expect("cg");
        // Hand solution is x = [1,1,1].
        for v in &x {
            assert!((v - 1.0).abs() < 1e-8, "got {v}");
        }
        assert!(rep.converged);
    }

    #[test]
    fn cg_zero_rhs() {
        let (a, _) = spd();
        let b = vec![0.0; 3];
        let (x, rep) = cg(&a, &b, None, 1e-10, 10).unwrap();
        assert!(rep.converged);
        assert!(x.iter().all(|v| *v == 0.0));
    }

    #[test]
    fn cg_rejects_rectangular() {
        // 2x3 operator that actually spans three columns.
        let a = CsrMock(csr_from_dense(2, 3, &[1.0, 0.0, 1.0, 0.0, 1.0, 0.0]));
        let b = vec![1.0, 1.0];
        assert!(matches!(cg(&a, &b, None, 1e-9, 5), Err(SolverError::NotSquare { .. })));
    }
}
