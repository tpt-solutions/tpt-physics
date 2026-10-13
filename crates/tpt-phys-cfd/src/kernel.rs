//! Collide / stream / moment kernels for [`Lbm2D`].
//!
//! The distributions are stored Structure-of-Arrays (`f[k][node]`), so four
//! neighbouring nodes of one velocity are contiguous and a `wide::f64x4`
//! lane processes four nodes at once. Each step is three row-parallel passes
//! over two buffers (no per-step copy):
//!
//! 1. **collide** `f → g`: relax toward equilibrium (SIMD on nodes whose whole
//!    3×3 neighbourhood is fluid, scalar next to solids / lattice edges),
//!    accumulating the boundary momentum exchange;
//! 2. **stream** `g → f`: pull each distribution from its upwind neighbour
//!    (plain slice copies in the bulk; scalar bounce-back / inlet / outlet
//!    logic at boundaries);
//! 3. **moments**: refresh `ρ`, `u` from `f`.
//!
//! All arithmetic is written once, generic over [`Lane`] (`f64` or `f64x4`),
//! in the exact operation order of the original scalar scheme, so the SIMD and
//! scalar paths agree bit-for-bit and `step` equals `step_par`.

// Index loops over the nine velocity slots mirror the lattice tables
// (`E`, `W`, `OPP`) and are clearer than iterator chains here.
#![allow(clippy::needless_range_loop)]

use std::ops::{Add, Div, Mul, Sub};

use rayon::prelude::*;
use wide::f64x4;

use crate::collision::Collision;
use crate::lattice::D2Q9;
use crate::{clamp_wall, wrap, Lbm2D, XBoundary};

/// SIMD width in nodes.
const W: usize = 4;

/// A scalar `f64` or a `f64x4` of four nodes: the arithmetic the kernels need.
pub(crate) trait Lane:
    Copy + Add<Output = Self> + Sub<Output = Self> + Mul<Output = Self> + Div<Output = Self>
{
    fn splat(x: f64) -> Self;
}

impl Lane for f64 {
    #[inline(always)]
    fn splat(x: f64) -> Self {
        x
    }
}

impl Lane for f64x4 {
    #[inline(always)]
    fn splat(x: f64) -> Self {
        f64x4::splat(x)
    }
}

/// Density and velocity of one node (or four) from its nine distributions.
#[inline(always)]
fn moments<L: Lane>(f: &[L; 9]) -> (L, L, L) {
    let rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
    let ux = (f[1] - f[3] + f[5] - f[6] - f[7] + f[8]) / rho;
    let uy = (f[2] - f[4] + f[5] + f[6] - f[7] - f[8]) / rho;
    (rho, ux, uy)
}

/// Second-order D2Q9 equilibrium.
#[inline(always)]
fn equilibrium<L: Lane>(rho: L, ux: L, uy: L) -> [L; 9] {
    let zero = L::splat(0.0);
    let u2 = ux * ux + uy * uy;
    let eu = [
        zero,
        ux,
        uy,
        zero - ux,
        zero - uy,
        ux + uy,
        uy - ux,
        zero - ux - uy,
        ux - uy,
    ];
    let one = L::splat(1.0);
    let mut feq = [zero; 9];
    for k in 0..9 {
        let e = eu[k];
        feq[k] = L::splat(D2Q9::W[k])
            * rho
            * (one + L::splat(3.0) * e + L::splat(4.5) * e * e - L::splat(1.5) * u2);
    }
    feq
}

/// Relaxation rates and operator for one step.
#[derive(Clone, Copy)]
pub(crate) struct Relax {
    pub bgk: bool,
    pub wp: f64,
    pub wm: f64,
}

impl Relax {
    pub fn new(collision: Collision, tau: f64) -> Self {
        let (wp, wm) = collision.rates(tau);
        Relax {
            bgk: matches!(collision, Collision::Bgk),
            wp,
            wm,
        }
    }
}

/// Collide one node (or four): body-force velocity shift, equilibrium, then
/// BGK or TRT relaxation.
#[inline(always)]
fn collide<L: Lane>(f: [L; 9], force: [f64; 2], r: Relax) -> [L; 9] {
    let (rho, ux, uy) = moments(&f);
    // Velocity used in the equilibrium is shifted by half the force
    // (explicit Guo-style forcing): u* = u + ½ F.
    let uxs = ux + L::splat(0.5) * L::splat(force[0]);
    let uys = uy + L::splat(0.5) * L::splat(force[1]);
    let feq = equilibrium(rho, uxs, uys);
    let wp = L::splat(r.wp);
    let mut out = [L::splat(0.0); 9];
    if r.bgk {
        for k in 0..9 {
            out[k] = f[k] - wp * (f[k] - feq[k]);
        }
    } else {
        let wm = L::splat(r.wm);
        let half = L::splat(0.5);
        out[0] = f[0] - wp * (f[0] - feq[0]);
        // Opposite pairs: (1,3) (2,4) (5,7) (6,8).
        for (a, b) in [(1usize, 3usize), (2, 4), (5, 7), (6, 8)] {
            let fs = (f[a] + f[b]) * half;
            let fa = (f[a] - f[b]) * half;
            let es = (feq[a] + feq[b]) * half;
            let ea = (feq[a] - feq[b]) * half;
            let ds = wp * (fs - es);
            let da = wm * (fa - ea);
            out[a] = f[a] - ds - da;
            out[b] = f[b] - ds + da;
        }
    }
    out
}

/// Row-parallel (or serial) map over per-row work items, order-preserving.
fn map_rows<T: Send, R: Send>(
    rows: Vec<T>,
    par: bool,
    f: impl Fn(usize, T) -> R + Sync + Send,
) -> Vec<R> {
    if par {
        rows.into_par_iter()
            .enumerate()
            .map(|(i, t)| f(i, t))
            .collect()
    } else {
        rows.into_iter().enumerate().map(|(i, t)| f(i, t)).collect()
    }
}

/// Split a nine-lane SoA buffer into per-row mutable slices.
fn split_rows(buf: &mut [Vec<f64>; 9], nx: usize) -> Vec<[&mut [f64]; 9]> {
    let ny = buf[0].len() / nx;
    let mut its = buf.each_mut().map(|v| v.chunks_mut(nx));
    (0..ny)
        .map(|_| its.each_mut().map(|it| it.next().expect("row")))
        .collect()
}

fn as_refs(buf: &[Vec<f64>; 9]) -> [&[f64]; 9] {
    buf.each_ref().map(|v| v.as_slice())
}

/// Read-only data shared by the stream kernel.
struct StreamCtx<'a> {
    nx: usize,
    ny: usize,
    solid: &'a [bool],
    clear: &'a [bool],
    rho: &'a [f64],
    wall_vel: &'a [[f64; 2]],
    xb: XBoundary,
    /// Post-collision distributions `g[k][node]`.
    g: [&'a [f64]; 9],
}

impl Lbm2D {
    /// Recompute the "clear" mask (fluid node whose in-lattice 3×3
    /// neighbourhood is all fluid) when the public `solid` array changed.
    fn refresh_clear(&mut self) {
        if self.clear.len() == self.solid.len() && self.solid_key == self.solid {
            return;
        }
        let (nx, ny) = (self.nx, self.ny);
        self.solid_key.clone_from(&self.solid);
        self.clear.clear();
        self.clear.resize(nx * ny, false);
        for iy in 0..ny {
            for ix in 0..nx {
                let mut ok = true;
                'n: for dy in -1_isize..=1 {
                    for dx in -1_isize..=1 {
                        let (sx, sy) = (ix as isize + dx, iy as isize + dy);
                        if sx >= 0
                            && sx < nx as isize
                            && sy >= 0
                            && sy < ny as isize
                            && self.solid[sy as usize * nx + sx as usize]
                        {
                            ok = false;
                            break 'n;
                        }
                    }
                }
                self.clear[iy * nx + ix] = ok;
            }
        }
    }

    /// One full step (collide, stream, moments); `par` selects Rayon rows.
    pub(crate) fn step_impl(&mut self, force: [f64; 2], par: bool) {
        self.refresh_clear();
        let (nx, ny) = (self.nx, self.ny);
        let relax = Relax::new(self.collision, self.tau);

        // 1. collide: f -> g.
        {
            let f = as_refs(&self.f);
            let solid = &self.solid[..];
            let clear = &self.clear[..];
            let rows = split_rows(&mut self.g, nx);
            let parts = map_rows(rows, par, |iy, dst| {
                collide_row(nx, ny, iy, solid, clear, &f, dst, force, relax)
            });
            // Fixed (row-major) summation order: deterministic for any thread
            // count, and identical for `step` and `step_par`.
            self.boundary_force = parts
                .iter()
                .fold([0.0; 2], |a, p| [a[0] + p[0], a[1] + p[1]]);
        }

        // 2. stream: g -> f.
        {
            let ctx = StreamCtx {
                nx,
                ny,
                solid: &self.solid,
                clear: &self.clear,
                rho: &self.rho,
                wall_vel: &self.wall_vel,
                xb: self.x_boundary,
                g: as_refs(&self.g),
            };
            let rows = split_rows(&mut self.f, nx);
            map_rows(rows, par, |iy, dst| stream_row(&ctx, iy, dst));
        }

        // 3. moments.
        self.moments_impl(par);
    }

    /// Refresh `ρ, u` from `f` (solid nodes are zeroed).
    pub(crate) fn moments_impl(&mut self, par: bool) {
        let nx = self.nx;
        let f = as_refs(&self.f);
        let solid = &self.solid[..];
        let rows: Vec<_> = self
            .rho
            .chunks_mut(nx)
            .zip(self.ux.chunks_mut(nx))
            .zip(self.uy.chunks_mut(nx))
            .map(|((r, u), v)| (r, u, v))
            .collect();
        map_rows(rows, par, |iy, (rho, ux, uy)| {
            let base = iy * nx;
            let mut ix = 0;
            while ix < nx {
                if ix + W <= nx && !solid[base + ix..base + ix + W].iter().any(|&s| s) {
                    let fl: [f64x4; 9] = f.map(|k| load4(&k[base + ix..]));
                    let (r, u, v) = moments(&fl);
                    rho[ix..ix + W].copy_from_slice(&r.to_array());
                    ux[ix..ix + W].copy_from_slice(&u.to_array());
                    uy[ix..ix + W].copy_from_slice(&v.to_array());
                    ix += W;
                } else {
                    let i = base + ix;
                    if solid[i] {
                        (rho[ix], ux[ix], uy[ix]) = (0.0, 0.0, 0.0);
                    } else {
                        let fl: [f64; 9] = f.map(|k| k[i]);
                        (rho[ix], ux[ix], uy[ix]) = moments(&fl);
                    }
                    ix += 1;
                }
            }
        });
    }
}

#[inline(always)]
fn load4(s: &[f64]) -> f64x4 {
    f64x4::new([s[0], s[1], s[2], s[3]])
}

/// Collide row `iy`; returns the boundary momentum exchange of its nodes.
#[allow(clippy::too_many_arguments)]
fn collide_row(
    nx: usize,
    ny: usize,
    iy: usize,
    solid: &[bool],
    clear: &[bool],
    f: &[&[f64]; 9],
    dst: [&mut [f64]; 9],
    force: [f64; 2],
    relax: Relax,
) -> [f64; 2] {
    let base = iy * nx;
    let mut bf = [0.0_f64; 2];
    let mut ix = 0;
    while ix < nx {
        if ix + W <= nx && clear[base + ix..base + ix + W].iter().all(|&c| c) {
            // Whole 3×3 neighbourhoods are fluid: no boundary exchange.
            let fl: [f64x4; 9] = f.map(|k| load4(&k[base + ix..]));
            let post = collide(fl, force, relax);
            for k in 0..9 {
                dst[k][ix..ix + W].copy_from_slice(&post[k].to_array());
            }
            ix += W;
            continue;
        }
        let i = base + ix;
        if !solid[i] {
            let fl: [f64; 9] = f.map(|k| k[i]);
            let post = collide(fl, force, relax);
            for k in 0..9 {
                dst[k][ix] = post[k];
                // Momentum exchange: the distribution heading into a solid
                // neighbour bounces straight back, handing the body
                // `2 f_post e_k` of momentum this step (half-way bounce-back,
                // Ladd). Neighbours outside the lattice are never solid.
                let (ex, ey) = D2Q9::E[k];
                let sx = ix as isize + ex as isize;
                let sy = iy as isize + ey as isize;
                if sx >= 0
                    && sx < nx as isize
                    && sy >= 0
                    && sy < ny as isize
                    && solid[sy as usize * nx + sx as usize]
                {
                    bf[0] += 2.0 * post[k] * ex as f64;
                    bf[1] += 2.0 * post[k] * ey as f64;
                }
            }
        }
        ix += 1;
    }
    bf
}

/// Stream into row `iy`: a pull of every distribution from its upwind node.
fn stream_row(c: &StreamCtx, iy: usize, mut out: [&mut [f64]; 9]) {
    let (nx, ny) = (c.nx, c.ny);
    let base = iy * nx;
    let interior_row = iy >= 1 && iy + 1 < ny;
    let fast = |ix: usize| interior_row && ix >= 1 && ix + 1 < nx && c.clear[base + ix];
    let mut ix = 0;
    while ix < nx {
        if fast(ix) {
            // Bulk run: every source is an in-lattice fluid node, so
            // streaming is a shifted slice copy per velocity.
            let mut end = ix + 1;
            while end < nx && fast(end) {
                end += 1;
            }
            for k in 0..9 {
                let (ex, ey) = D2Q9::E[k];
                let sy = (iy as isize - ey as isize) as usize;
                let s0 = sy * nx + (ix as isize - ex as isize) as usize;
                out[k][ix..end].copy_from_slice(&c.g[k][s0..s0 + (end - ix)]);
            }
            ix = end;
        } else {
            stream_node(c, ix, iy, &mut out);
            ix += 1;
        }
    }
}

/// Scalar streaming of one node with the full boundary logic: half-way
/// bounce-back (with the Ladd moving-wall correction), periodic / inlet /
/// open `x` boundaries.
fn stream_node(c: &StreamCtx, ix: usize, iy: usize, out: &mut [&mut [f64]; 9]) {
    let (nx, ny) = (c.nx, c.ny);
    let i = iy * nx + ix;
    if c.solid[i] {
        return;
    }
    let g = |idx: usize, k: usize| c.g[k][idx];

    // Zou–Hou reconstruction for the east outflow column of an `Inlet` run:
    // the three west-pointing distributions (k = 3, 6, 7) are missing (their
    // source lies outside the lattice), so they are rebuilt from the known
    // post-stream set with the static density pinned to the reference — a
    // mass-consistent pressure outlet. (A plain zero-gradient copy here let
    // the domain fill up over time, decelerating the flow and decaying the
    // body force.)
    let zh = if ix == nx - 1 && matches!(c.xb, XBoundary::Inlet(_)) {
        let row_s = iy.saturating_sub(1);
        let row_n = (iy + 1).min(ny - 1);
        let known = (
            g(i, 0),
            g(iy * nx + ix - 1, 1),
            g(row_s * nx + ix, 2),
            g(row_n * nx + ix, 4),
            g(row_s * nx + ix - 1, 5),
            g(row_n * nx + ix - 1, 8),
        );
        Some(zou_hou_east(known))
    } else {
        None
    };

    for k in 0..9 {
        let (ex, ey) = D2Q9::E[k];
        let sx = ix as isize - ex as isize;
        let sy = iy as isize - ey as isize;
        let val = if sx < 0 || sx >= nx as isize {
            match c.xb {
                XBoundary::Periodic => {
                    let wx = wrap(sx, nx as isize);
                    let wy = clamp_wall(sy, ny as isize);
                    g(wy * nx + wx, k)
                }
                XBoundary::Inlet(u) => {
                    if sx < 0 {
                        // Entering from the west inlet.
                        D2Q9::equilibrium(1.0, [u, 0.0])[k]
                    } else {
                        // East open (outflow) boundary: a density-pinned
                        // Zou–Hou pressure outlet for the west-pointing
                        // distributions, zero-gradient copy for the rest.
                        match k {
                            3 => zh.unwrap()[0],
                            6 => zh.unwrap()[1],
                            7 => zh.unwrap()[2],
                            _ => g(iy * nx + (nx - 2), k),
                        }
                    }
                }
                XBoundary::Open => {
                    // Downstream outflow: the west edge wraps (periodic,
                    // well-posed inflow); the east edge uses first-order
                    // upwind / zero-gradient, copying its neighbouring
                    // interior column so that mass and waves leave the domain
                    // instead of reflecting.
                    if sx < 0 {
                        let wx = wrap(sx, nx as isize);
                        let wy = clamp_wall(sy, ny as isize);
                        g(wy * nx + wx, k)
                    } else {
                        g(iy * nx + (nx - 2), k)
                    }
                }
            }
        } else if sy < 0 || sy >= ny as isize {
            let wx = wrap(sx, nx as isize);
            let wy = clamp_wall(sy, ny as isize);
            let sidx = wy * nx + wx;
            if c.solid[sidx] {
                bounce(c, i, k, sidx)
            } else {
                g(sidx, k)
            }
        } else {
            let sidx = sy as usize * nx + sx as usize;
            if c.solid[sidx] {
                bounce(c, i, k, sidx)
            } else {
                g(sidx, k)
            }
        };
        out[k][ix] = val;
    }
}

/// Zou–Hou (1997) density-specified outflow reconstruction for the east
/// column: given the known post-stream distributions at the boundary node
/// `(f0, f1, f2, f4, f5, f8)` and the pinned outlet density (the lattice
/// reference ρ = 1, with `u_y = 0`), returns the three unknown west-pointing
/// distributions `[f3, f6, f7]`. `u_x` follows from combining mass
/// conservation with x-momentum; the diagonal unknowns use the
/// non-equilibrium bounce-back closure and the normal unknown is set by exact
/// mass conservation (which is what pins the domain's density).
#[inline]
fn zou_hou_east(known: (f64, f64, f64, f64, f64, f64)) -> [f64; 3] {
    let (f0, f1, f2, f4, f5, f8) = known;
    let rho = 1.0;
    let ux = (2.0 * (f1 + f5 + f8) + f0 + f2 + f4 - rho) / rho;
    let ux = ux.clamp(-0.3, 0.3);
    let feq = D2Q9::equilibrium(rho, [ux, 0.0]);
    let f6 = feq[6] - (f5 - feq[5]);
    let f7 = feq[7] - (f8 - feq[8]);
    let f3 = rho - (f0 + f1 + f2 + f4 + f5 + f8) - f6 - f7;
    [f3, f6, f7]
}

/// Half-way bounce-back value for fluid node `i`, direction `k`, whose source
/// is the solid node `sidx`. Applies the Ladd moving-wall correction when the
/// solid node carries a non-zero wall velocity.
#[inline]
fn bounce(c: &StreamCtx, i: usize, k: usize, sidx: usize) -> f64 {
    let wv = c.wall_vel[sidx];
    let base = c.g[D2Q9::OPP[k]][i];
    if wv[0] == 0.0 && wv[1] == 0.0 {
        base
    } else {
        let rho_i = c.rho[i];
        let eu = D2Q9::E[k].0 as f64 * wv[0] + D2Q9::E[k].1 as f64 * wv[1];
        base + 2.0 * D2Q9::W[k] * rho_i * eu / D2Q9::CS2
    }
}
