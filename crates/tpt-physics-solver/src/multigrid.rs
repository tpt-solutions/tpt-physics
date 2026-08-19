//! Lightweight two-level geometric multigrid preconditioner for structured
//! 2-D Poisson (5-point Laplacian) problems.
//!
//! This is the "lightweight AMG/multigrid" follow-up to the preconditioned
//! CG/GMRES solvers. It is **experimental**: it assumes the operator is a grid-
//! structured 2-D Laplacian (the common SPD diffusion system from FEA/CFD
//! thermal/viscous solves), builds a half-resolution coarse grid by rediscreti-
//! zation, and applies a V-cycle (weighted-Jacobi smoothing + direct coarse
//! solve) as the preconditioner `M⁻¹`.
//!
//! It is *not* a general black-box algebraic multigrid yet — for an arbitrary
//! sparse matrix use the Jacobi preconditioner in [`crate::cg_pc`] /
//! [`crate::gmres_pc`].

use crate::linalg::LinearOperator;
use tpt_fem_sparse::{Coo, Csr};

/// A structured 2-D grid with `nx × ny` nodes, indexed `j*nx + i`.
#[derive(Debug, Clone, Copy)]
pub struct Grid2D {
    /// Number of nodes in `x`.
    pub nx: usize,
    /// Number of nodes in `y`.
    pub ny: usize,
}

impl Grid2D {
    /// Build a grid.
    pub fn new(nx: usize, ny: usize) -> Self {
        Grid2D { nx, ny }
    }

    /// Linear index of node `(i, j)`.
    #[inline]
    pub fn idx(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }

    /// Total node count.
    pub fn len(&self) -> usize {
        self.nx * self.ny
    }

    /// `true` when the grid has no nodes.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Assemble the standard 5-point Laplacian `A` for the operator `-Δ` with
    /// unit spacing (the interior stencil `4·center − 4·neighbours`). Boundary
    /// rows are left as the stencil sees them; callers pin Dirichlet rows as
    /// needed.
    pub fn laplacian_csr(&self) -> Csr {
        let mut coo = Coo::with_capacity(self.len() * 5);
        for j in 0..self.ny {
            for i in 0..self.nx {
                let r = self.idx(i, j);
                coo.push(r, r, 4.0);
                if i > 0 {
                    coo.push(r, self.idx(i - 1, j), -1.0);
                }
                if i + 1 < self.nx {
                    coo.push(r, self.idx(i + 1, j), -1.0);
                }
                if j > 0 {
                    coo.push(r, self.idx(i, j - 1), -1.0);
                }
                if j + 1 < self.ny {
                    coo.push(r, self.idx(i, j + 1), -1.0);
                }
            }
        }
        coo.to_csr()
    }
}

/// A two-level geometric multigrid V-cycle preconditioner for a 2-D Laplacian.
pub struct Multigrid {
    fine: Csr,
    coarse: Csr,
    coarse_nx: usize,
    coarse_ny: usize,
    fine_nx: usize,
    fine_ny: usize,
    /// Pre-smoothing sweeps.
    nu_pre: usize,
    /// Post-smoothing sweeps.
    nu_post: usize,
}

impl Multigrid {
    /// Build a V-cycle preconditioner from a fine grid. The coarse grid is the
    /// half-resolution Laplacian. Both dimensions must be even.
    pub fn new(fine: Grid2D, nu_pre: usize, nu_post: usize) -> Self {
        assert!(
            fine.nx % 2 == 0 && fine.ny % 2 == 0,
            "grid must be even-sized"
        );
        let coarse = Grid2D::new(fine.nx / 2, fine.ny / 2);
        Multigrid {
            fine: fine.laplacian_csr(),
            coarse: coarse.laplacian_csr(),
            coarse_nx: coarse.nx,
            coarse_ny: coarse.ny,
            fine_nx: fine.nx,
            fine_ny: fine.ny,
            nu_pre,
            nu_post,
        }
    }

    /// Apply the V-cycle: `z = M⁻¹ b` (zero initial guess).
    pub fn apply(&self, z: &mut [f64], b: &[f64]) {
        let n = self.fine.nrows;
        let mut x = vec![0.0; n];
        for _ in 0..self.nu_pre {
            jacobi(&self.fine, b, &mut x);
        }
        // Residual r = b - A x.
        let mut ax = vec![0.0; n];
        self.fine.apply(&x, &mut ax);
        let r: Vec<f64> = (0..n).map(|i| b[i] - ax[i]).collect();

        // Restrict to coarse grid and solve the coarse system directly.
        let rc = self.restrict(&r);
        let xc = solve_dense_csr(&self.coarse, &rc);

        // Prolongate (interpolate) the coarse correction into `x`.
        self.prolongate_add(&xc, &mut x);

        for _ in 0..self.nu_post {
            jacobi(&self.fine, b, &mut x);
        }
        z.copy_from_slice(&x);
    }

    /// Full-weighting restriction of a fine-grid vector to the coarse grid.
    fn restrict(&self, fine: &[f64]) -> Vec<f64> {
        let nx = self.fine_nx;
        let ny = self.fine_ny;
        let cnx = self.coarse_nx;
        let cny = self.coarse_ny;
        let mut out = vec![0.0; cnx * cny];
        for (jj, j) in (0..cny).enumerate() {
            for (ii, i) in (0..cnx).enumerate() {
                let fi = 2 * ii;
                let fj = 2 * jj;
                let mut acc = 0.0;
                let mut wsum = 0.0;
                for (dj, dw) in [(0i64, 4.0), (-1, 2.0), (1, 2.0)] {
                    for (di, dw2) in [(0i64, 4.0), (-1, 2.0), (1, 2.0)] {
                        let wi = fi as i64 + di;
                        let wj = fj as i64 + dj;
                        if wi < 0 || wi >= nx as i64 || wj < 0 || wj >= ny as i64 {
                            continue;
                        }
                        let w = (dw * dw2) / 16.0;
                        acc += w * fine[(wj as usize) * nx + (wi as usize)];
                        wsum += w;
                    }
                }
                if wsum > 0.0 {
                    out[j * cnx + i] = acc / wsum;
                }
            }
        }
        out
    }

    /// Bilinear interpolation of a coarse-grid correction into the fine grid
    /// (added to `fine`).
    fn prolongate_add(&self, coarse: &[f64], fine: &mut [f64]) {
        let nx = self.fine_nx;
        let ny = self.fine_ny;
        let cnx = self.coarse_nx;
        let cny = self.coarse_ny;
        for j in 0..ny {
            for i in 0..nx {
                let ci = (i / 2).min(cnx - 1);
                let cj = (j / 2).min(cny - 1);
                let fx = (i % 2) as f64;
                let fy = (j % 2) as f64;
                let w00 = (1.0 - fx) * (1.0 - fy);
                let w10 = fx * (1.0 - fy);
                let w01 = (1.0 - fx) * fy;
                let w11 = fx * fy;
                let mut v = 0.0;
                v += w00 * coarse[cj * cnx + ci];
                if ci + 1 < cnx {
                    v += w10 * coarse[cj * cnx + (ci + 1)];
                }
                if cj + 1 < cny {
                    v += w01 * coarse[(cj + 1) * cnx + ci];
                }
                if ci + 1 < cnx && cj + 1 < cny {
                    v += w11 * coarse[(cj + 1) * cnx + (ci + 1)];
                }
                fine[j * nx + i] += v;
            }
        }
    }
}

/// Weighted-Jacobi smoothing in place: `x ← x + ω D⁻¹ (b − A x)`.
fn jacobi(a: &Csr, b: &[f64], x: &mut [f64]) {
    let w = 2.0 / 3.0;
    let n = a.nrows;
    let mut off = vec![0.0; n];
    let mut diag = vec![1.0; n];
    for r in 0..n {
        let mut s = 0.0;
        for idx in a.row_ptrs[r]..a.row_ptrs[r + 1] {
            let c = a.col_ind[idx];
            if c == r {
                diag[r] = a.values[idx];
            } else {
                s += a.values[idx] * x[c];
            }
        }
        off[r] = s;
    }
    for r in 0..n {
        x[r] = (1.0 - w) * x[r] + w * (b[r] - off[r]) / diag[r];
    }
}

/// Direct dense solve of a small SPD `Csr` system via partial-pivot Gaussian
/// elimination (used for the coarse grid).
fn solve_dense_csr(a: &Csr, b: &[f64]) -> Vec<f64> {
    let n = a.nrows;
    let mut m = vec![vec![0.0; n]; n];
    for (r, m_r) in m.iter_mut().enumerate() {
        for idx in a.row_ptrs[r]..a.row_ptrs[r + 1] {
            m_r[a.col_ind[idx]] = a.values[idx];
        }
    }
    // Augment with RHS and eliminate.
    let mut aug = vec![vec![0.0; n + 1]; n];
    for i in 0..n {
        for j in 0..n {
            aug[i][j] = m[i][j];
        }
        aug[i][n] = b[i];
    }
    for col in 0..n {
        let mut piv = col;
        let mut best = aug[col][col].abs();
        for (r, aug_r) in aug.iter().enumerate().skip(col + 1) {
            if aug_r[col].abs() > best {
                best = aug_r[col].abs();
                piv = r;
            }
        }
        aug.swap(col, piv);
        let d = aug[col][col];
        for aug_col in aug[col].iter_mut().skip(col).take(n + 1 - col) {
            *aug_col /= d;
        }
        for r in 0..n {
            if r != col {
                let f = aug[r][col];
                let col_row: Vec<f64> = aug[col][col..=n].to_vec();
                for (j, aug_rj) in aug[r].iter_mut().enumerate().skip(col).take(n + 1 - col) {
                    *aug_rj -= f * col_row[j - col];
                }
            }
        }
    }
    (0..n).map(|i| aug[i][n]).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gmres::gmres;
    use crate::gmres::gmres_pc;
    use crate::linalg::csr_from_dense;

    // A small 2×2 SPD system used to check `solve_dense_csr`.
    #[test]
    fn dense_csr_solve_matches_reference() {
        let a = csr_from_dense(2, 2, &[4.0, -1.0, -1.0, 3.0]);
        let x = solve_dense_csr(&a, &[3.0, 2.0]);
        // solution is x = [1, 1]
        assert!((x[0] - 1.0).abs() < 1e-9 && (x[1] - 1.0).abs() < 1e-9);
    }

    // 2-D Poisson: GMRES preconditioned with the multigrid V-cycle should
    // converge in no more iterations than a plain Jacobi preconditioner,
    // demonstrating that the V-cycle is an effective coarse correction.
    #[test]
    fn mg_pc_beats_jacobi_on_poisson() {
        let grid = Grid2D::new(32, 32);
        let a = grid.laplacian_csr();
        // RHS: a smooth source so the solution is smooth (ideal for MG).
        let n = grid.len();
        let b: Vec<f64> = (0..n)
            .map(|idx| {
                let i = idx % grid.nx;
                let j = idx / grid.nx;
                let x = i as f64 / (grid.nx as f64 - 1.0);
                let y = j as f64 / (grid.ny as f64 - 1.0);
                (x - 0.5) * (y - 0.5)
            })
            .collect();

        // Reference: unpreconditioned GMRES iteration count.
        let (_, rep_none) = gmres(&a, &b, None, 50, 1e-8, 200).unwrap();

        // Jacobi preconditioner.
        let dinv: Vec<f64> = (0..a.nrows)
            .map(|i| 1.0 / a.values[a.row_ptrs[i]])
            .collect();
        let jacobi = |r: &[f64], z: &mut [f64]| {
            for i in 0..r.len() {
                z[i] = dinv[i] * r[i];
            }
        };
        let (_, rep_jac) = gmres_pc(&a, &b, None, 50, 1e-8, 200, Some(&jacobi)).unwrap();

        // Multigrid V-cycle preconditioner.
        let mg = Multigrid::new(grid, 2, 2);
        let mg_closure = |r: &[f64], z: &mut [f64]| mg.apply(z, r);
        let (x_mg, rep_mg) =
            gmres_pc(&a, &b, None, 50, 1e-8, 200, Some(&mg_closure)).unwrap();

        // Verify the MG solution is correct (A x ≈ b).
        let mut ax = vec![0.0; n];
        a.apply(&x_mg, &mut ax);
        let mut res = 0.0;
        for i in 0..n {
            res += (ax[i] - b[i]).powi(2);
        }
        assert!(res.sqrt() / (n as f64).sqrt() < 1e-6, "mg residual {res}");

        // The multigrid V-cycle must be at least as effective as Jacobi (and on
        // a smooth Laplacian it strictly reduces the iteration count).
        assert!(
            rep_mg.iterations <= rep_jac.iterations,
            "mg {} iters vs jacobi {} iters vs none {}",
            rep_mg.iterations,
            rep_jac.iterations,
            rep_none.iterations
        );
        assert!(rep_mg.converged);
    }
}
