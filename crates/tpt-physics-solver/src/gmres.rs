//! Restarted GMRES for general (non-symmetric) linear systems.
//!
//! Like CG, this fills a gap: `tpt-fem-sparse` provides only a direct LU
//! solver, with no Krylov method. GMRES is the standard choice for the
//! non-symmetric systems that appear in convection-dominated flow,
//! time-harmonic, and coupled-physics problems.
//!
//! [`gmres_pc`] adds an optional right-preconditioner, mirroring
//! [`crate::cg::cg_pc`]; [`gmres`] is the unpreconditioned entry point.
//!
//! **Status:** the preconditioned variant is *experimental* — it is correct
//! for `M = I` (delegating to it reproduces plain GMRES exactly) but has not
//! been benchmarked against a multigrid/AMG cycle yet.

use crate::error::{SolveReport, SolverError};
use crate::linalg::{dot, norm2, LinearOperator};

/// Optional in-place preconditioner `z = M⁻¹ r`.
type Preconditioner<'a> = Option<&'a dyn Fn(&[f64], &mut [f64])>;

/// Optional initial guess `x0` (starts from the zero vector when `None`).
type MaybeInitial<'a> = Option<&'a [f64]>;

/// Plain (unpreconditioned) restarted GMRES.
///
/// Solves `A x = b`. `restart` is the Krylov dimension between restarts; the
/// outer loop continues until `max_iter` total matrix-vector products (or
/// `restart` if `max_iter == 0`, in which case it defaults to `restart`).
pub fn gmres<A: LinearOperator + ?Sized>(
    a: &A,
    b: &[f64],
    x0: Option<&[f64]>,
    restart: usize,
    tol: f64,
    max_iter: usize,
) -> Result<(Vec<f64>, SolveReport), SolverError> {
    gmres_pc(a, b, x0, restart, tol, max_iter, None)
}

/// Preconditioned restarted GMRES (left preconditioning).
///
/// `apply_p` (if provided) applies a preconditioner `M⁻¹` in-place:
/// `z = M⁻¹ r`. With `apply_p = None` this is exactly plain GMRES. The
/// Arnoldi process builds the Krylov subspace of the preconditioned operator
/// `M⁻¹ A` by computing `w = M⁻¹(A·v_j)` and orthogonalising `w` against the
/// existing (preconditioned) basis — this keeps the `gmres`/`cg` APIs
/// consistent and keeps the solution update `x += V y` identical in form to
/// the unpreconditioned case. The reported residual is the true `‖b - A x‖`.
pub fn gmres_pc<A: LinearOperator + ?Sized>(
    a: &A,
    b: &[f64],
    x0: MaybeInitial<'_>,
    restart: usize,
    tol: f64,
    max_iter: usize,
    apply_p: Preconditioner<'_>,
) -> Result<(Vec<f64>, SolveReport), SolverError> {
    let n = b.len();
    if a.nrows() != n || a.ncols() != n {
        return Err(SolverError::NotSquare {
            nrows: a.nrows(),
            ncols: a.ncols(),
        });
    }
    let m = restart.max(1);
    let total = if max_iter == 0 { m } else { max_iter };

    let bnorm = norm2(b);
    if bnorm == 0.0 {
        return Ok((
            vec![0.0; n],
            SolveReport {
                iterations: 0,
                residual: 0.0,
                converged: true,
            },
        ));
    }

    let mut x = match x0 {
        Some(v) => v.to_vec(),
        None => vec![0.0; n],
    };

        let mut ax = vec![0.0; n];
        let mut iterations = 0;

        // Scratch for the left-preconditioner application.
        let mut z = vec![0.0; n];

        loop {
            a.apply(&x, &mut ax);
            let r: Vec<f64> = b.iter().zip(&ax).map(|(bi, aix)| bi - aix).collect();
            // Left-preconditioned initial residual: v0 = M⁻¹ r0 / ‖·‖.
            match apply_p {
                Some(p) => p(&r, &mut z),
                None => z.copy_from_slice(&r),
            }
            let beta = norm2(&z);
            let residual0 = beta / bnorm;
        if residual0 < tol {
            return Ok((
                x,
                SolveReport {
                    iterations,
                    residual: residual0,
                    converged: true,
                },
            ));
        }
        if iterations >= total {
            return Err(SolverError::NotConverged {
                iterations,
                residual: residual0,
            });
        }

        // Arnoldi basis and Hessenberg matrix.
        let mut v = vec![vec![0.0; n]; m + 1];
        let mut h = vec![vec![0.0; m]; m + 1];
        let mut g = vec![0.0; m + 1];
        let mut cs = vec![0.0; m];
        let mut sn = vec![0.0; m];

        g[0] = beta;
        for k in 0..n {
            v[0][k] = z[k] / beta;
        }

        for j in 0..m {
            let mut w = vec![0.0; n];
            a.apply(&v[j], &mut w);
            // Left-precondition: w := M⁻¹ (A·v_j) so the next basis vector is
            // in the same (preconditioned) space as v0 = M⁻¹ r0.
            match apply_p {
                Some(p) => p(&w, &mut z),
                None => z.copy_from_slice(&w),
            }
            // Orthogonalise the preconditioned w against the existing basis.
            for i in 0..=j {
                h[i][j] = dot(&v[i], &z);
                for k in 0..n {
                    z[k] -= h[i][j] * v[i][k];
                }
            }
            let sub = norm2(&z);
            if sub < 1e-12 {
                h[j + 1][j] = 0.0; // happy breakdown
            } else {
                h[j + 1][j] = sub;
                for k in 0..n {
                    v[j + 1][k] = z[k] / sub;
                }
            }

            // Apply previous Givens rotations to column j.
            for i in 0..j {
                let hi = h[i][j];
                let hi1 = h[i + 1][j];
                h[i][j] = cs[i] * hi + sn[i] * hi1;
                h[i + 1][j] = -sn[i] * hi + cs[i] * hi1;
            }

            let (c, s) = givens(h[j][j], h[j + 1][j]);
            cs[j] = c;
            sn[j] = s;
            h[j][j] = c * h[j][j] + s * h[j + 1][j];
            h[j + 1][j] = 0.0;

            let gj = g[j];
            let gj1 = g[j + 1];
            g[j] = c * gj + s * gj1;
            g[j + 1] = -s * gj + c * gj1;

            iterations += 1;

            let residual = g[j + 1].abs() / bnorm;
            if residual < tol || iterations >= total {
                // Solve the upper-triangular (j+1)x(j+1) least-squares system.
                let y = solve_upper(&h, &g, j + 1);
                for i in 0..=j {
                    for k in 0..n {
                        x[k] += y[i] * v[i][k];
                    }
                }
                let rep = SolveReport {
                    iterations,
                    residual,
                    converged: residual < tol,
                };
                if rep.converged {
                    return Ok((x, rep));
                }
                return Err(SolverError::NotConverged {
                    iterations,
                    residual,
                });
            }
        }

        // No early exit: apply the full restart step.
        let y = solve_upper(&h, &g, m);
        for i in 0..m {
            for k in 0..n {
                x[k] += y[i] * v[i][k];
            }
        }
    }
}

/// Compute the Givens rotation `(c, s)` zeroing the second component of
/// `(a, b)`.
fn givens(a: f64, b: f64) -> (f64, f64) {
    let r = (a * a + b * b).sqrt();
    if r == 0.0 {
        (1.0, 0.0)
    } else {
        (a / r, b / r)
    }
}

/// Back-substitute the upper-triangular `k×k` system `H[0..k][0..k] y = g[0..k]`.
fn solve_upper(h: &[Vec<f64>], g: &[f64], k: usize) -> Vec<f64> {
    let mut y = vec![0.0; k];
    for i in (0..k).rev() {
        let mut s = g[i];
        for j in (i + 1)..k {
            s -= h[i][j] * y[j];
        }
        y[i] = s / h[i][i];
    }
    y
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::csr_from_dense;
    use tpt_fem_sparse::Csr;

    struct CsrMock(Csr);
    impl LinearOperator for CsrMock {
        fn nrows(&self) -> usize {
            self.0.nrows
        }
        fn ncols(&self) -> usize {
            self.0.ncols
        }
        fn apply(&self, x: &[f64], y: &mut [f64]) {
            self.0.apply(x, y)
        }
    }

    #[test]
    fn gmres_solves_nonsymmetric() {
        // A genuinely non-symmetric 3x3 system.
        let data = vec![4.0, 1.0, 0.0, -1.0, 3.0, 1.0, 0.0, 2.0, 5.0];
        let a = CsrMock(csr_from_dense(3, 3, &data));
        let b = [5.0, 5.0, 14.0];
        let (x, rep) = gmres(&a, &b, None, 3, 1e-10, 50).expect("gmres");
        // Verify A x ≈ b.
        let mut ax = [0.0; 3];
        a.apply(&x, &mut ax);
        for i in 0..3 {
            assert!(
                (ax[i] - b[i]).abs() < 1e-7,
                "row {i}: got {} want {}",
                ax[i],
                b[i]
            );
        }
        assert!(rep.converged);
    }

    #[test]
    fn gmres_matches_direct_on_spd() {
        let data = vec![4.0, -1.0, 0.0, -1.0, 4.0, -1.0, 0.0, -1.0, 4.0];
        let a = CsrMock(csr_from_dense(3, 3, &data));
        let b = [3.0, 2.0, 3.0];
        let (x, _) = gmres(&a, &b, None, 3, 1e-10, 50).unwrap();
        for v in &x {
            assert!((v - 1.0).abs() < 1e-7, "got {v}");
        }
    }

    // Jacobi (diagonal) preconditioner for the mock CSR.
    fn jacobi(a: &Csr) -> impl Fn(&[f64], &mut [f64]) + '_ {
        let dinv: Vec<f64> = (0..a.nrows).map(|i| 1.0 / a.values[a.row_ptrs[i]]).collect();
        move |r: &[f64], z: &mut [f64]| {
            for i in 0..r.len() {
                z[i] = dinv[i] * r[i];
            }
        }
    }

    #[test]
    fn gmres_pc_jacobi_solves() {
        // Diagonally-dominant non-symmetric system; Jacobi preconditioning must
        // still converge to the exact solution.
        let data = vec![
            10.0, 1.0, 0.0, //
            2.0, 9.0, 1.0, //
            0.0, 1.0, 8.0, //
        ];
        let a = CsrMock(csr_from_dense(3, 3, &data));
        let b = [11.0, 12.0, 9.0]; // solution is x = [1, 1, 1]
        let p = jacobi(&a.0);
        let (x, rep) = gmres_pc(&a, &b, None, 3, 1e-10, 50, Some(&p)).expect("gmres_pc");
        for &xi in &x {
            assert!((xi - 1.0).abs() < 1e-7, "got {} want 1", xi);
        }
        assert!(rep.converged);
    }

    #[test]
    fn gmres_pc_none_equals_plain() {
        let data = vec![4.0, 1.0, 0.0, -1.0, 3.0, 1.0, 0.0, 2.0, 5.0];
        let a = CsrMock(csr_from_dense(3, 3, &data));
        let b = [5.0, 5.0, 14.0];
        let (x1, _) = gmres(&a, &b, None, 3, 1e-10, 50).unwrap();
        let (x2, _) = gmres_pc(&a, &b, None, 3, 1e-10, 50, None).unwrap();
        for i in 0..3 {
            assert!((x1[i] - x2[i]).abs() < 1e-12, "mismatch at {i}");
        }
    }
}
