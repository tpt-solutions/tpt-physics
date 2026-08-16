//! Restarted GMRES for general (non-symmetric) linear systems.
//!
//! Like CG, this fills a gap: `tpt-fem-sparse` provides only a direct LU
//! solver, with no Krylov method. GMRES is the standard choice for the
//! non-symmetric systems that appear in convection-dominated flow,
//! time-harmonic, and coupled-physics problems.

use crate::error::{SolveReport, SolverError};
use crate::linalg::{dot, norm2, LinearOperator};

/// Restarted GMRES.
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

    loop {
        a.apply(&x, &mut ax);
        let r: Vec<f64> = b.iter().zip(&ax).map(|(bi, aix)| bi - aix).collect();
        let beta = norm2(&r);
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
            v[0][k] = r[k] / beta;
        }

        for j in 0..m {
            let mut w = vec![0.0; n];
            a.apply(&v[j], &mut w);
            for i in 0..=j {
                h[i][j] = dot(&v[i], &w);
                for k in 0..n {
                    w[k] -= h[i][j] * v[i][k];
                }
            }
            let sub = norm2(&w);
            if sub < 1e-12 {
                h[j + 1][j] = 0.0; // happy breakdown
            } else {
                h[j + 1][j] = sub;
                for k in 0..n {
                    v[j + 1][k] = w[k] / sub;
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
}
