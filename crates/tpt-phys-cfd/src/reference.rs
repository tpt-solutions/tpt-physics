//! Test-only scalar reference implementation of the BGK step.
//!
//! This is the original Array-of-Structures collide / stream / moment scheme
//! the optimised SoA + SIMD + Rayon kernels in [`crate::kernel`] replaced. It
//! is kept (verbatim in structure) as an independent oracle: the equivalence
//! tests assert the production kernels reproduce it.

use crate::lattice::D2Q9;
use crate::{clamp_wall, wrap, Lbm2D, XBoundary};

impl Lbm2D {
    /// One BGK step using the scalar reference scheme.
    pub(crate) fn step_reference(&mut self, force: [f64; 2]) {
        let (nx, ny) = (self.nx, self.ny);
        let n = nx * ny;
        let tau = self.tau;
        let mut f: Vec<[f64; 9]> = (0..n).map(|i| self.node_f(i)).collect();

        // Collide.
        let mut boundary_force = [0.0_f64; 2];
        for iy in 0..ny {
            for ix in 0..nx {
                let i = iy * nx + ix;
                if self.solid[i] {
                    continue;
                }
                let fi = f[i];
                let rho = fi.iter().copied().sum::<f64>();
                let ux = (0..9).map(|k| fi[k] * D2Q9::E[k].0 as f64).sum::<f64>() / rho;
                let uy = (0..9).map(|k| fi[k] * D2Q9::E[k].1 as f64).sum::<f64>() / rho;
                let feq = D2Q9::equilibrium(rho, [ux + 0.5 * force[0], uy + 0.5 * force[1]]);
                for k in 0..9 {
                    let f_post = fi[k] - (1.0 / tau) * (fi[k] - feq[k]);
                    f[i][k] = f_post;
                    let (ex, ey) = D2Q9::E[k];
                    let sx = ix as isize + ex as isize;
                    let sy = iy as isize + ey as isize;
                    if sx >= 0
                        && sx < nx as isize
                        && sy >= 0
                        && sy < ny as isize
                        && self.solid[sy as usize * nx + sx as usize]
                    {
                        boundary_force[0] += 2.0 * f_post * ex as f64;
                        boundary_force[1] += 2.0 * f_post * ey as f64;
                    }
                }
            }
        }
        self.boundary_force = boundary_force;

        // Stream.
        let prev = f.clone();
        let bounce = |i: usize, k: usize, sidx: usize| -> f64 {
            let wv = self.wall_vel[sidx];
            let base = prev[i][D2Q9::OPP[k]];
            if wv[0] == 0.0 && wv[1] == 0.0 {
                base
            } else {
                let eu = D2Q9::E[k].0 as f64 * wv[0] + D2Q9::E[k].1 as f64 * wv[1];
                base + 2.0 * D2Q9::W[k] * self.rho[i] * eu / D2Q9::CS2
            }
        };
        for iy in 0..ny {
            for ix in 0..nx {
                let i = iy * nx + ix;
                if self.solid[i] {
                    continue;
                }
                let zh = if ix == nx - 1 && matches!(self.x_boundary, XBoundary::Inlet(_)) {
                    let rs = iy.saturating_sub(1);
                    let rn = (iy + 1).min(ny - 1);
                    let (f0, f1, f2, f4, f5, f8) = (
                        prev[i][0],
                        prev[iy * nx + ix - 1][1],
                        prev[rs * nx + ix][2],
                        prev[rn * nx + ix][4],
                        prev[rs * nx + ix - 1][5],
                        prev[rn * nx + ix - 1][8],
                    );
                    let ux = ((2.0 * (f1 + f5 + f8) + f0 + f2 + f4 - 1.0) / 1.0).clamp(-0.3, 0.3);
                    let feq = D2Q9::equilibrium(1.0, [ux, 0.0]);
                    let f6 = feq[6] - (f5 - feq[5]);
                    let f7 = feq[7] - (f8 - feq[8]);
                    let f3 = 1.0 - (f0 + f1 + f2 + f4 + f5 + f8) - f6 - f7;
                    Some([f3, f6, f7])
                } else {
                    None
                };
                for k in 0..9 {
                    let (ex, ey) = D2Q9::E[k];
                    let sx = ix as isize - ex as isize;
                    let sy = iy as isize - ey as isize;
                    let val = if sx < 0 || sx >= nx as isize {
                        match self.x_boundary {
                            XBoundary::Periodic => {
                                prev[clamp_wall(sy, ny as isize) * nx + wrap(sx, nx as isize)][k]
                            }
                            XBoundary::Inlet(u) => {
                                if sx < 0 {
                                    D2Q9::equilibrium(1.0, [u, 0.0])[k]
                                } else {
                                    match k {
                                        3 => zh.unwrap()[0],
                                        6 => zh.unwrap()[1],
                                        7 => zh.unwrap()[2],
                                        _ => prev[iy * nx + (nx - 2)][k],
                                    }
                                }
                            }
                            XBoundary::Open => {
                                if sx < 0 {
                                    prev[clamp_wall(sy, ny as isize) * nx + wrap(sx, nx as isize)]
                                        [k]
                                } else {
                                    prev[iy * nx + (nx - 2)][k]
                                }
                            }
                        }
                    } else if sy < 0 || sy >= ny as isize {
                        let sidx = clamp_wall(sy, ny as isize) * nx + wrap(sx, nx as isize);
                        if self.solid[sidx] {
                            bounce(i, k, sidx)
                        } else {
                            prev[sidx][k]
                        }
                    } else {
                        let sidx = sy as usize * nx + sx as usize;
                        if self.solid[sidx] {
                            bounce(i, k, sidx)
                        } else {
                            prev[sidx][k]
                        }
                    };
                    f[i][k] = val;
                }
            }
        }

        // Moments.
        for (i, &fi) in f.iter().enumerate() {
            let (rho, ux, uy) = if self.solid[i] {
                (0.0, 0.0, 0.0)
            } else {
                let rho = fi.iter().copied().sum::<f64>();
                (
                    rho,
                    (0..9).map(|k| fi[k] * D2Q9::E[k].0 as f64).sum::<f64>() / rho,
                    (0..9).map(|k| fi[k] * D2Q9::E[k].1 as f64).sum::<f64>() / rho,
                )
            };
            self.rho[i] = rho;
            self.ux[i] = ux;
            self.uy[i] = uy;
        }
        for (i, fi) in f.into_iter().enumerate() {
            self.set_node_f(i, fi);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn build(
        nx: usize,
        ny: usize,
        tau: f64,
        xb: u8,
        walls: u8,
        blocks: &[(usize, usize)],
        lid: f64,
    ) -> Lbm2D {
        let mut s = Lbm2D::new(nx, ny, tau);
        s.set_x_boundary(match xb {
            0 => XBoundary::Periodic,
            1 => XBoundary::Inlet(0.04),
            _ => XBoundary::Open,
        });
        match walls {
            1 => s.set_horizontal_walls(),
            2 => {
                s.set_box_walls();
                s.set_moving_lid(ny - 1, lid);
            }
            _ => {}
        }
        for &(bx, by) in blocks {
            let (x, y) = (bx % nx, by % ny);
            s.add_rect(x, y, x + 1, y + 1);
        }
        s.initialise(1.0, [0.03, 0.01]);
        s
    }

    fn assert_same(a: &Lbm2D, b: &Lbm2D, tol: f64) {
        for i in 0..a.nx * a.ny {
            for k in 0..9 {
                let (x, y) = (a.f[k][i], b.f[k][i]);
                assert!((x - y).abs() <= tol, "f[{k}][{i}]: {x} vs {y}");
            }
            for (p, q) in [(a.rho[i], b.rho[i]), (a.ux[i], b.ux[i]), (a.uy[i], b.uy[i])] {
                assert!((p - q).abs() <= tol, "moment at {i}: {p} vs {q}");
            }
        }
        for c in 0..2 {
            assert!((a.boundary_force[c] - b.boundary_force[c]).abs() <= 1e-10);
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(24))]

        /// The SoA/SIMD kernels (serial and Rayon) reproduce the scalar
        /// reference across lattice sizes (incl. non-multiples of the SIMD
        /// width), relaxation times, forcing, obstacles, moving walls and all
        /// three `x` boundary treatments.
        #[test]
        fn kernels_match_scalar_reference(
            nx in 4usize..23,
            ny in 3usize..14,
            tau in 0.52f64..1.4,
            xb in 0u8..3,
            walls in 0u8..3,
            blocks in proptest::collection::vec((0usize..64, 0usize..64), 0..4),
            fx in -1e-4f64..1e-4,
            fy in -1e-4f64..1e-4,
            lid in 0.0f64..0.08,
            steps in 1usize..12,
        ) {
            let mut r = build(nx, ny, tau, xb, walls, &blocks, lid);
            let mut a = r.clone();
            let mut b = r.clone();
            for _ in 0..steps {
                r.step_reference([fx, fy]);
                a.step([fx, fy]);
                b.step_par([fx, fy]);
            }
            assert_same(&r, &a, 1e-12);
            assert_same(&a, &b, 0.0);
        }
    }

    #[test]
    fn trt_with_bgk_magic_matches_bgk() {
        let tau = 0.7;
        let mut a = build(13, 9, tau, 0, 1, &[(5, 4)], 0.0);
        let mut b = a.clone();
        b.set_collision(Collision::Trt {
            magic: (tau - 0.5) * (tau - 0.5),
        });
        for _ in 0..20 {
            a.step([1e-5, 0.0]);
            b.step([1e-5, 0.0]);
        }
        assert_same(&a, &b, 1e-10);
    }

    use crate::Collision;

    #[test]
    fn mrt_with_uniform_rates_is_bgk() {
        let tau = 0.8;
        let w = 1.0 / tau;
        let mut a = build(13, 9, tau, 0, 2, &[(5, 4)], 0.05);
        let mut b = a.clone();
        b.set_collision(Collision::Mrt {
            s_e: w,
            s_eps: w,
            s_q: w,
        });
        for _ in 0..30 {
            a.step([1e-5, 0.0]);
            b.step([1e-5, 0.0]);
        }
        assert_same(&a, &b, 1e-10);
    }

    #[test]
    fn mrt_simd_matches_serial_and_conserves_mass() {
        let mut a = build(23, 11, 0.6, 0, 2, &[], 0.08);
        a.set_collision(Collision::MRT_DEFAULT);
        let mut b = a.clone();
        for _ in 0..50 {
            a.step([0.0, 0.0]);
            b.step_par([0.0, 0.0]);
        }
        assert_same(&a, &b, 0.0);
        let mass: f64 = (0..a.nx * a.ny)
            .filter(|&i| !a.solid[i])
            .map(|i| a.rho[i])
            .sum();
        let fluid = a.solid.iter().filter(|&&s| !s).count() as f64;
        assert!(mass.is_finite() && (mass - fluid).abs() / fluid < 1e-3);
    }
}
