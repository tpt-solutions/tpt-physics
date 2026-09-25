//! Computational Fluid Dynamics for `tpt-physics`.
//!
//! This crate provides a pure-Rust, data-oriented Lattice Boltzmann Method
//! (LBM) solver for incompressible flow. No equivalent CFD code exists in the
//! sibling `tpt-math` / `tpt-fem` workspaces, so this is genuinely net-new.
//!
//! The current implementation is the canonical **D2Q9** (9-velocity,
//! two-dimensional) BGK-collision lattice with:
//!
//! * half-way bounce-back solid boundaries (stationary or moving),
//! * periodic, or velocity-inlet / zero-gradient-outlet `x` boundaries,
//! * a body-force (Guo-style velocity-shift) term so that pressure-gradient /
//!   gravity-driven flows can be simulated,
//! * a circular obstacle primitive for bluff-body flows,
//! * `rayon`-friendly bulk loops.

pub mod lattice;
pub mod sph;

use lattice::D2Q9;
use rayon::prelude::*;

/// Streamwise (`x`) boundary treatment.
#[derive(Debug, Clone, Copy)]
pub enum XBoundary {
    /// Wrap around in `x` (periodic channel / recirculating wake).
    Periodic,
    /// Velocity inlet at the west edge (`u = (U, 0)`) and a zero-gradient /
    /// convective outlet at the east edge — used for external
    /// (flow-past-body) simulations.
    Inlet(f64),
    /// Open (non-reflective) *downstream* boundary: the upstream (west) edge
    /// wraps (periodic) while the downstream (east) edge adopts the
    /// distributions of its neighbouring interior column (first-order upwind /
    /// zero-gradient). Waves and mass therefore leave the domain at the east
    /// instead of reflecting, while the west keeps a well-posed (periodic)
    /// inflow. Use this for free-outflow / recirculating problems (plumes,
    /// wakes) where a fully periodic domain is undesirable.
    Open,
}

/// A two-dimensional D2Q9 Lattice Boltzmann solver.
///
/// The lattice holds `nx * ny` nodes. A node is either fluid or solid
/// (boundary). Distributions are stored Structure-of-Arrays as a flat `Vec` of
/// length `nx * ny` of `[f64; 9]` (one lane per lattice velocity).
pub struct Lbm2D {
    /// Number of nodes in the `x` (streamwise) direction.
    pub nx: usize,
    /// Number of nodes in the `y` (wall-normal) direction.
    pub ny: usize,
    /// BGK relaxation time `τ` (>`0.5` for stability).
    pub tau: f64,
    /// Post-collision distributions `f_i(x)` for every node.
    f: Vec<[f64; 9]>,
    /// Macroscopic density `ρ`.
    pub rho: Vec<f64>,
    /// Macroscopic x-velocity `u`.
    pub ux: Vec<f64>,
    /// Macroscopic y-velocity `v`.
    pub uy: Vec<f64>,
    /// `true` for solid (no-slip) nodes.
    pub solid: Vec<bool>,
    /// Wall velocities `[u, v]` for moving (lid) boundaries; zero = stationary.
    wall_vel: Vec<[f64; 2]>,
    /// Streamwise (`x`) boundary treatment.
    x_boundary: XBoundary,
    /// Fluid force on the solid boundaries during the most recent step,
    /// `[Fx, Fy]`, accumulated by boundary momentum exchange in
    /// [`Lbm2D::collide`].
    boundary_force: [f64; 2],
    /// Reused streaming source buffer (a copy of the previous `f`): keeps
    /// `stream` allocation-free — a fresh `Vec` per step showed up as real
    /// cost in the throughput budget.
    f_swap: Vec<[f64; 9]>,
}

impl Lbm2D {
    /// Construct an `nx × ny` lattice with relaxation time `tau`.
    ///
    /// ```
    /// use tpt_phys_cfd::Lbm2D;
    /// let mut lbm = Lbm2D::new(24, 24, 0.53);
    /// lbm.step([0.0, 0.0]); // one collision+stream step
    /// ```
    pub fn new(nx: usize, ny: usize, tau: f64) -> Self {
        assert!(nx > 2 && ny > 1, "lattice must be at least 3x2");
        assert!(tau > 0.5, "tau must exceed 0.5 for stability");
        let n = nx * ny;
        Lbm2D {
            nx,
            ny,
            tau,
            f: vec![[0.0; 9]; n],
            rho: vec![1.0; n],
            ux: vec![0.0; n],
            uy: vec![0.0; n],
            solid: vec![false; n],
            wall_vel: vec![[0.0; 2]; n],
            x_boundary: XBoundary::Periodic,
            boundary_force: [0.0; 2],
            f_swap: Vec::new(),
        }
    }

    /// Linear index of node `(ix, iy)`.
    #[inline]
    pub fn idx(&self, ix: usize, iy: usize) -> usize {
        iy * self.nx + ix
    }

    /// Mark the first and last rows (`iy = 0` and `iy = ny-1`) as solid walls.
    /// Together with periodic x this forms a 2-D channel.
    pub fn set_horizontal_walls(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        for ix in 0..nx {
            self.solid[ix] = true;
            self.solid[(ny - 1) * nx + ix] = true;
        }
    }

    /// Mark all four borders as solid walls (fully enclosed cavity).
    pub fn set_box_walls(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        for ix in 0..nx {
            self.solid[ix] = true;
            self.solid[(ny - 1) * nx + ix] = true;
        }
        for iy in 0..ny {
            self.solid[iy * nx] = true;
            self.solid[iy * nx + (nx - 1)] = true;
        }
    }

    /// Set the streamwise (`x`) boundary treatment.
    pub fn set_x_boundary(&mut self, b: XBoundary) {
        self.x_boundary = b;
    }

    /// Make a horizontal strip of wall nodes move with velocity `v` in `x`
    /// (used as the lid of a lid-driven cavity). `iy` is the row (e.g. `ny-1`).
    pub fn set_moving_lid(&mut self, iy: usize, v: f64) {
        let nx = self.nx;
        for ix in 0..nx {
            let i = self.idx(ix, iy);
            self.solid[i] = true;
            self.wall_vel[i] = [v, 0.0];
        }
    }

    /// Set a single node `(ix, iy)` as a moving wall with lattice velocity
    /// `v = [u, w]`. The node is marked solid so the Ladd bounce-back
    /// correction applies; used by FSI back-coupling to drag the fluid with a
    /// moving structural boundary.
    pub fn set_wall_velocity(&mut self, ix: usize, iy: usize, v: [f64; 2]) {
        let i = self.idx(ix, iy);
        self.solid[i] = true;
        self.wall_vel[i] = v;
    }

    /// Add a solid circular obstacle centred at `(cx, cy)` with lattice radius
    /// `r` (no-slip bounce-back on every node inside).
    pub fn add_circle(&mut self, cx: f64, cy: f64, r: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let r2 = r * r;
        for iy in 0..ny {
            for ix in 0..nx {
                let dx = ix as f64 - cx;
                let dy = iy as f64 - cy;
                if dx * dx + dy * dy <= r2 {
                    self.solid[iy * nx + ix] = true;
                }
            }
        }
    }

    /// Add an axis-aligned solid rectangle `[x0..x1] × [y0..y1]` (inclusive).
    pub fn add_rect(&mut self, x0: usize, y0: usize, x1: usize, y1: usize) {
        let nx = self.nx;
        for iy in y0..=y1.min(self.ny - 1) {
            for ix in x0..=x1.min(self.nx - 1) {
                self.solid[iy * nx + ix] = true;
            }
        }
    }

    /// Initialise the whole lattice to equilibrium for density `rho0` and
    /// velocity `u0`.
    pub fn initialise(&mut self, rho0: f64, u0: [f64; 2]) {
        for i in 0..self.rho.len() {
            if self.solid[i] {
                continue;
            }
            self.rho[i] = rho0;
            self.ux[i] = u0[0];
            self.uy[i] = u0[1];
            self.f[i] = D2Q9::equilibrium(rho0, u0);
        }
    }

    /// Add a localised transverse velocity "bump" — a cos²-enveloped
    /// perturbation centred at `(cx, cy)` with lattice radius `r`. Used to
    /// seed the Kármán instability behind bluff bodies: a perfectly
    /// symmetric initial field can hold the (unstable) steady wake branch
    /// for tens of thousands of steps before round-off asymmetry grows, so
    /// unsteady studies (and honest time-averaged drag) start from a
    /// slightly broken symmetry instead.
    pub fn add_velocity_bump(&mut self, cx: f64, cy: f64, r: f64, dv: [f64; 2]) {
        let nx = self.nx;
        let ny = self.ny;
        let r2 = r * r;
        for iy in 0..ny {
            for ix in 0..nx {
                let i = iy * nx + ix;
                if self.solid[i] {
                    continue;
                }
                let dx = ix as f64 - cx;
                let dy = iy as f64 - cy;
                let d2 = dx * dx + dy * dy;
                if d2 >= r2 {
                    continue;
                }
                let w = (1.0 - d2 / r2).powi(2); // cos²-style envelope
                let u = [self.ux[i] + dv[0] * w, self.uy[i] + dv[1] * w];
                self.f[i] = D2Q9::equilibrium(self.rho[i], u);
            }
        }
    }

    /// Kinematic viscosity `ν = cs² (τ − ½)`, with `cs² = 1/3` for D2Q9.
    pub fn viscosity(&self) -> f64 {
        D2Q9::CS2 * (self.tau - 0.5)
    }

    /// Advance the simulation by one collision+streaming step under body-force
    /// density `force = (Fx, Fy)` (acceleration per unit mass, lattice units).
    pub fn step(&mut self, force: [f64; 2]) {
        self.collide(force);
        self.stream();
        self.macro_fields();
    }

    /// BGK collision with the velocity-shift body-force scheme. While
    /// post-collision distributions are written, every link whose neighbour
    /// is solid contributes its momentum exchange `2 f_post e_k` to the
    /// per-step [`Lbm2D::boundary_force`] (half-way bounce-back, Ladd).
    fn collide(&mut self, force: [f64; 2]) {
        let tau = self.tau;
        let nx = self.nx;
        let ny = self.ny;
        let mut boundary_force = [0.0_f64; 2];
        for iy in 0..ny {
            for ix in 0..nx {
                let i = iy * nx + ix;
                if self.solid[i] {
                    continue;
                }
                let f = self.f[i];
                let rho = f.iter().copied().sum::<f64>();
                let ux = (0..9).map(|k| f[k] * D2Q9::E[k].0 as f64).sum::<f64>() / rho;
                let uy = (0..9).map(|k| f[k] * D2Q9::E[k].1 as f64).sum::<f64>() / rho;

                // Velocity used in the equilibrium is shifted by half the force
                // (Guo forcing, explicit form): u* = u + ½ F.
                let uxs = ux + 0.5 * force[0];
                let uys = uy + 0.5 * force[1];

                let feq = D2Q9::equilibrium(rho, [uxs, uys]);
                for k in 0..9 {
                    let f_post = f[k] - (1.0 / tau) * (f[k] - feq[k]);
                    self.f[i][k] = f_post;
                    // Momentum exchange: the distribution heading into a
                    // solid neighbour bounces straight back, handing the
                    // body `2 f_post e_k` of momentum this step. Neighbours
                    // outside the lattice are never solid, so no wrap
                    // lookup is needed here. (The moving-wall correction is
                    // ignored — force tracking targets stationary bodies.)
                    let (ex, ey) = D2Q9::E[k];
                    let sx = ix as isize + ex as isize;
                    let sy = iy as isize + ey as isize;
                    if sx >= 0 && sx < nx as isize && sy >= 0 && sy < ny as isize {
                        if self.solid[sy as usize * nx + sx as usize] {
                            boundary_force[0] += 2.0 * f_post * ex as f64;
                            boundary_force[1] += 2.0 * f_post * ey as f64;
                        }
                    }
                }
            }
        }
        self.boundary_force = boundary_force;
    }

    /// Streaming with half-way bounce-back at solid nodes, moving-wall
    /// corrections, wrap/periodicity or inlet–outlet in `x`.
    fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        self.f_swap.clear();
        self.f_swap.extend_from_slice(&self.f);
        let f_prev = &self.f_swap;
        for iy in 0..ny {
            for ix in 0..nx {
                let i = self.idx(ix, iy);
                if self.solid[i] {
                    continue;
                }
                // Zou–Hou reconstruction for the east outflow column of an
                // `Inlet` run: the three west-pointing distributions (k = 3,
                // 6, 7) are missing (their source lies outside the lattice),
                // so they are rebuilt from the known post-stream set with the
                // static density pinned to the reference — a mass-consistent
                // pressure outlet. (A plain zero-gradient copy here let the
                // domain fill up over time, decelerating the flow and
                // decaying the body force.)
                let zh = if ix == nx - 1 && matches!(self.x_boundary, XBoundary::Inlet(_)) {
                    let row_s = iy.saturating_sub(1);
                    let row_n = (iy + 1).min(ny - 1);
                    let known = (
                        f_prev[i][0],
                        f_prev[iy * nx + ix - 1][1],
                        f_prev[row_s * nx + ix][2],
                        f_prev[row_n * nx + ix][4],
                        f_prev[row_s * nx + ix - 1][5],
                        f_prev[row_n * nx + ix - 1][8],
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
                        match self.x_boundary {
                            XBoundary::Periodic => {
                                let wx = wrap(sx, nx as isize);
                                let wy = clamp_wall(sy, ny as isize);
                                f_prev[wy * nx + wx][k]
                            }
                            XBoundary::Inlet(u) => {
                                if sx < 0 {
                                    // Entering from the west inlet.
                                    D2Q9::equilibrium(1.0, [u, 0.0])[k]
                                } else {
                                    // East open (outflow) boundary: a
                                    // density-pinned Zou–Hou pressure outlet
                                    // for the west-pointing distributions,
                                    // zero-gradient copy for the rest.
                                    match k {
                                        3 => zh.unwrap()[0],
                                        6 => zh.unwrap()[1],
                                        7 => zh.unwrap()[2],
                                        _ => f_prev[iy * nx + (nx - 2)][k],
                                    }
                                }
                            }
                            XBoundary::Open => {
                                // Downstream outflow: the west edge wraps
                                // (periodic, well-posed inflow); the east edge
                                // uses first-order upwind / zero-gradient,
                                // copying its neighbouring interior column so
                                // that mass and waves leave the domain instead
                                // of reflecting.
                                if sx < 0 {
                                    let wx = wrap(sx, nx as isize);
                                    let wy = clamp_wall(sy, ny as isize);
                                    f_prev[wy * nx + wx][k]
                                } else {
                                    let c = nx - 2;
                                    let ci = iy * nx + c;
                                    f_prev[ci][k]
                                }
                            }
                        }
                    } else if sy < 0 || sy >= ny as isize {
                        let wx = wrap(sx, nx as isize);
                        let wy = clamp_wall(sy, ny as isize);
                        let sidx = wy * nx + wx;
                        if self.solid[sidx] {
                            bounce(f_prev, self, i, k, sidx)
                        } else {
                            f_prev[sidx][k]
                        }
                    } else {
                        let sidx = sy as usize * nx + sx as usize;
                        if self.solid[sidx] {
                            bounce(f_prev, self, i, k, sidx)
                        } else {
                            f_prev[sidx][k]
                        }
                    };
                    self.f[i][k] = val;
                }
            }
        }
    }

    /// Recompute macroscopic density and velocity (used as a post-step check and
    /// externally after stepping).
    pub fn macro_fields(&mut self) -> (&Vec<f64>, &Vec<f64>, &Vec<f64>) {
        for i in 0..self.rho.len() {
            if self.solid[i] {
                self.rho[i] = 0.0;
                self.ux[i] = 0.0;
                self.uy[i] = 0.0;
                continue;
            }
            let f = self.f[i];
            let rho = f.iter().copied().sum::<f64>();
            let ux = (0..9).map(|k| f[k] * D2Q9::E[k].0 as f64).sum::<f64>() / rho;
            let uy = (0..9).map(|k| f[k] * D2Q9::E[k].1 as f64).sum::<f64>() / rho;
            self.rho[i] = rho;
            self.ux[i] = ux;
            self.uy[i] = uy;
        }
        (&self.rho, &self.ux, &self.uy)
    }
}

/// Zou–Hou (1997) density-specified outflow reconstruction for the east
/// column: given the known post-stream distributions at the boundary node
/// `(f0, f1, f2, f4, f5, f8)` and the pinned outlet density (the lattice
/// reference ρ = 1, with `u_y = 0`), returns the three unknown
/// west-pointing distributions `[f3, f6, f7]`. `u_x` follows from combining
/// mass conservation with x-momentum; the diagonal unknowns use the
/// non-equilibrium bounce-back closure and the normal unknown is set by
/// exact mass conservation (which is what pins the domain's density).
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
fn bounce(f_prev: &[[f64; 9]], sim: &Lbm2D, i: usize, k: usize, sidx: usize) -> f64 {
    let wv = sim.wall_vel[sidx];
    let base = f_prev[i][D2Q9::OPP[k]];
    if wv[0] == 0.0 && wv[1] == 0.0 {
        base
    } else {
        let rho_i = sim.rho[i];
        let eu = D2Q9::E[k].0 as f64 * wv[0] + D2Q9::E[k].1 as f64 * wv[1];
        base + 2.0 * D2Q9::W[k] * rho_i * eu / D2Q9::CS2
    }
}

#[inline]
fn wrap(x: isize, n: isize) -> usize {
    let mut r = x % n;
    if r < 0 {
        r += n;
    }
    r as usize
}

#[inline]
fn clamp_wall(y: isize, n: isize) -> usize {
    if y < 0 {
        0
    } else if y >= n {
        (n - 1) as usize
    } else {
        y as usize
    }
}

impl Lbm2D {
    /// Fluid force on the solid boundaries during the most recent
    /// [`Lbm2D::step`] / [`Lbm2D::step_par`], `[Fx, Fy]`, in lattice units
    /// per unit depth (δx = δt = 1), by boundary momentum exchange
    /// (half-way bounce-back). For an `XBoundary::Inlet(u)` external flow
    /// along +x this is the hydrodynamic force on the body: drag is `Fx`,
    /// lift `Fy`. A drag coefficient follows as
    /// `Fx / (½ ρ U² H)` with ρ ≈ 1 and `H` the body's lattice frontal
    /// height (number of solid-node rows).
    pub fn boundary_force(&self) -> [f64; 2] {
        self.boundary_force
    }

    /// Average x-velocity profile along `y` (averaged over the `x` interior).
    pub fn x_velocity_profile(&self) -> Vec<f64> {        (0..self.ny)
            .map(|iy| {
                let mut s = 0.0;
                for ix in 0..self.nx {
                    s += self.ux[self.idx(ix, iy)];
                }
                s / self.nx as f64
            })
            .collect()
    }

    /// Mean x-velocity along the vertical line `ix` (for wake / cavity studies).
    pub fn x_velocity_at_column(&self, ix: usize) -> Vec<f64> {
        (0..self.ny).map(|iy| self.ux[self.idx(ix, iy)]).collect()
    }
}

/// Convenience helpers that parallelise the collision sweep (used to exercise
/// the `rayon` dependency path without changing numerics).
impl Lbm2D {
    /// Like [`Lbm2D::step`] but runs the collision sweep in parallel.
    pub fn step_par(&mut self, force: [f64; 2]) {
        let tau = self.tau;
        let nx = self.nx;
        let ny = self.ny;
        let solid = &self.solid;
        let (f_new, boundary_force): (Vec<[f64; 9]>, Vec<[f64; 2]>) = (0..self.rho.len())
            .into_par_iter()
            .map(|i| {
                if solid[i] {
                    return ([0.0; 9], [0.0; 2]);
                }
                let iy = i / nx;
                let ix = i % nx;
                let f = self.f[i];
                let rho = f.iter().copied().sum::<f64>();
                let ux = (0..9).map(|k| f[k] * D2Q9::E[k].0 as f64).sum::<f64>() / rho;
                let uy = (0..9).map(|k| f[k] * D2Q9::E[k].1 as f64).sum::<f64>() / rho;
                let uxs = ux + 0.5 * force[0];
                let uys = uy + 0.5 * force[1];
                let feq = D2Q9::equilibrium(rho, [uxs, uys]);
                let mut out = [0.0; 9];
                let mut link_force = [0.0; 2];
                for k in 0..9 {
                    let f_post = f[k] - (1.0 / tau) * (f[k] - feq[k]);
                    out[k] = f_post;
                    // Same boundary momentum exchange as `collide`.
                    let (ex, ey) = D2Q9::E[k];
                    let sx = ix as isize + ex as isize;
                    let sy = iy as isize + ey as isize;
                    if sx >= 0 && sx < nx as isize && sy >= 0 && sy < ny as isize {
                        if solid[sy as usize * nx + sx as usize] {
                            link_force[0] += 2.0 * f_post * ex as f64;
                            link_force[1] += 2.0 * f_post * ey as f64;
                        }
                    }
                }
                (out, link_force)
            })
            .unzip();
        self.f = f_new;
        self.boundary_force = boundary_force
            .into_iter()
            .fold([0.0; 2], |mut acc, f| {
                acc[0] += f[0];
                acc[1] += f[1];
                acc
            });
        self.stream();
        self.macro_fields();
    }
}

/// Re-export so callers can build equilibria / inspect the lattice constants.
pub use lattice::D2Q9 as Lattice;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equilibrium_conserves_mass_and_momentum() {
        let f = D2Q9::equilibrium(1.3, [0.1, -0.05]);
        let rho: f64 = f.iter().sum();
        assert!((rho - 1.3).abs() < 1e-12, "rho = {rho}");
        let ux = f
            .iter()
            .enumerate()
            .map(|(k, fk)| fk * D2Q9::E[k].0 as f64)
            .sum::<f64>()
            / rho;
        let uy = f
            .iter()
            .enumerate()
            .map(|(k, fk)| fk * D2Q9::E[k].1 as f64)
            .sum::<f64>()
            / rho;
        assert!((ux - 0.1).abs() < 1e-12, "ux = {ux}");
        assert!((uy + 0.05).abs() < 1e-12, "uy = {uy}");
    }

    #[test]
    fn poiseuille_is_parabolic() {
        let nx = 48;
        let ny = 24;
        let tau = 0.6;
        let mut sim = Lbm2D::new(nx, ny, tau);
        sim.set_horizontal_walls();
        sim.initialise(1.0, [0.0, 0.0]);
        let fx = 1e-5;
        for _ in 0..20000 {
            sim.step([fx, 0.0]);
        }

        let prof = sim.x_velocity_profile();
        // No-slip at the walls.
        assert!(prof[0].abs() < 1e-9, "wall u = {}", prof[0]);
        assert!(prof[ny - 1].abs() < 1e-9, "wall u = {}", prof[ny - 1]);
        // Positive flow through the interior.
        for (y, &u) in prof.iter().enumerate().take(ny - 1).skip(1) {
            assert!(u > 0.0, "u[{y}] = {u}");
        }
        // Symmetric about the centre row.
        for y in 1..ny / 2 {
            assert!((prof[y] - prof[ny - 1 - y]).abs() < 1e-6, "asym at {y}");
        }
        // Poiseuille: d²u/dy² = -Fx/ν (independent of boundary convention).
        // The velocity-shift forcing scheme is first-order accurate, so the
        // measured curvature is ~5/6 of the ideal; we allow 20%.
        let nu = sim.viscosity();
        let yc = ny / 2;
        let d2 = prof[yc - 1] - 2.0 * prof[yc] + prof[yc + 1];
        let expected = -fx / nu;
        assert!(
            (d2 - expected).abs() / expected.abs() < 0.20,
            "d2 = {d2}, expected {expected}"
        );
    }

    #[test]
    fn fields_stay_finite_under_forcing() {
        let mut sim = Lbm2D::new(20, 16, 0.55);
        sim.set_horizontal_walls();
        sim.initialise(1.0, [0.0, 0.0]);
        for _ in 0..500 {
            sim.step([1e-6, 1e-7]);
        }
        let (rho, ux, uy) = sim.macro_fields();
        assert!(rho.iter().chain(ux).chain(uy).all(|v| v.is_finite()));
    }

    #[test]
    fn open_boundary_outflow_is_nonreflecting() {
        // A rightward-moving pulse in a periodic-upstream / open-downstream
        // domain must travel east and *leave* through the outflow instead of
        // reflecting back. After enough steps the upstream (west) region should
        // have relaxed toward rest while the fields stay finite and bounded.
        let nx = 80;
        let ny = 16;
        let tau = 0.6;
        let mut sim = Lbm2D::new(nx, ny, tau);
        sim.set_x_boundary(XBoundary::Open);
        // Rest everywhere...
        sim.initialise(1.0, [0.0, 0.0]);
        // ...with a localized rightward-moving blob in the left-centre.
        for iy in 0..ny {
            for ix in 30..50 {
                let i = sim.idx(ix, iy);
                sim.f[i] = D2Q9::equilibrium(1.0, [0.1, 0.0]);
            }
        }
        for _ in 0..500 {
            sim.step([0.0, 0.0]);
        }
        let (rho, ux, uy) = sim.macro_fields();
        assert!(rho.iter().chain(ux).chain(uy).all(|v| v.is_finite()));
        // The pulse has advected ~0.1·steps ≈ 50 cells to the right and exited
        // through the east outflow, so the upstream (west) quarter is at rest.
        let west_u: f64 = (0..nx / 4)
            .map(|ix| (0..ny).map(|iy| ux[iy * nx + ix].abs()).sum::<f64>())
            .sum();
        assert!(
            west_u < (nx / 4) as f64 * ny as f64 * 0.05,
            "upstream not relaxed: {west_u}"
        );
    }

    #[test]
    fn boundary_force_vanishes_without_solids() {
        // Uniform streaming with no solid nodes: nothing exchanges momentum
        // with a boundary, so the per-step force must stay at zero.
        let mut sim = Lbm2D::new(48, 24, 0.6);
        sim.initialise(1.0, [0.1, 0.0]);
        for _ in 0..500 {
            sim.step([0.0, 0.0]);
        }
        let [fx, fy] = sim.boundary_force();
        assert!(fx.abs() < 1e-12 && fy.abs() < 1e-12, "force = [{fx}, {fy}]");
    }

    #[test]
    fn cylinder_drag_is_positive_and_symmetric() {
        // Flow past a circle at Re = 20 (steady regime) through a velocity
        // inlet: the momentum-exchange force must point downstream, sit
        // symmetrically about the centreline (|lift| << drag), and give a
        // drag coefficient in the right ballpark (experimental Cd(Re=20) ≈
        // 2.0–2.2). Loose bounds here — this asserts the mechanism; the CFD
        // app's engine tests hold the quantitative Re sweep.
        let nx = 200;
        let ny = 100;
        let u0 = 0.05;
        let h = 24.0; // lattice frontal height of the disc below
        let re = 20.0;
        let nu = u0 * h / re;
        let mut sim = Lbm2D::new(nx, ny, 3.0 * nu + 0.5);
        sim.set_x_boundary(XBoundary::Inlet(u0));
        sim.add_circle(50.0, (ny - 1) as f64 / 2.0, h / 2.0);
        sim.initialise(1.0, [u0, 0.0]);
        for _ in 0..4000 {
            sim.step([0.0, 0.0]);
        }
        let [fx, fy] = sim.boundary_force();
        assert!(fx > 0.0, "drag must point downstream, fx = {fx}");
        assert!(
            fy.abs() < 0.05 * fx,
            "symmetric body cannot carry lift: fx = {fx}, fy = {fy}"
        );
        let cd = fx / (0.5 * u0 * u0 * h);
        assert!(
            (1.2..3.2).contains(&cd),
            "Cd(Re=20) = {cd} outside the sanity band"
        );
    }

    #[test]
    fn velocity_bump_seeds_shedding_promptly() {
        // Same Re = 100 cylinder, but with the Kármán instability seeded by
        // a small transverse bump right after initialisation: the lift must
        // be visibly oscillating within a few thousand steps instead of
        // after ~100k of metastable symmetric hold.
        let nx = 160;
        let ny = 80;
        let u0 = 0.075;
        let h = 16.0;
        let nu = u0 * h / 100.0;
        let mut sim = Lbm2D::new(nx, ny, 3.0 * nu + 0.5);
        sim.set_x_boundary(XBoundary::Inlet(u0));
        sim.add_circle(40.0, (ny - 1) as f64 / 2.0, h / 2.0);
        sim.initialise(1.0, [u0, 0.0]);
        sim.add_velocity_bump(48.0, (ny - 1) as f64 / 2.0, 6.0, [0.0, 0.1 * u0]);

        let mut cl_min = f64::INFINITY;
        let mut cl_max = f64::NEG_INFINITY;
        for step in 0..12_000 {
            sim.step([0.0, 0.0]);
            if step % 4 == 0 {
                let [_, fy] = sim.boundary_force();
                let cl = fy / (0.5 * u0 * u0 * h);
                cl_min = cl_min.min(cl);
                cl_max = cl_max.max(cl);
            }
        }
        assert!(
            cl_max - cl_min > 0.1,
            "seeded wake not shedding: Cl swung only {:.4}",
            cl_max - cl_min
        );
    }

    #[test]
    fn inlet_outflow_keeps_mass_and_drag_stationary() {
        // The Zou–Hou pressure outlet must not let the domain fill up or
        // drain: over a long Re = 20 run the mean fluid density stays put,
        // and the drag coefficient settles instead of drifting (a copy-based
        // outlet decays Cd from ~2.1 toward ~1.2 as back-pressure strangles
        // the through-flow). Doubles as the quantitative cylinder check:
        // settled Cd(Re=20) ≈ 2.0–2.2 (Tritton 1959, Braza et al. 1986).
        let nx = 160;
        let ny = 80;
        let u0 = 0.05;
        let h = 16.0;
        let re = 20.0;
        let nu = u0 * h / re;
        let mut sim = Lbm2D::new(nx, ny, 3.0 * nu + 0.5);
        sim.set_x_boundary(XBoundary::Inlet(u0));
        sim.add_circle(40.0, (ny - 1) as f64 / 2.0, h / 2.0);
        sim.initialise(1.0, [u0, 0.0]);

        let mean_rho = |sim: &Lbm2D| -> f64 {
            let n = sim.solid.iter().filter(|&&s| !s).count() as f64;
            sim.rho
                .iter()
                .zip(sim.solid.iter())
                .filter(|(_, &s)| !s)
                .map(|(&r, _)| r)
                .sum::<f64>()
                / n
        };
        let rho_start = mean_rho(&sim);
        let cd_at = |sim: &Lbm2D| -> f64 {
            let [fx, _] = sim.boundary_force();
            fx / (0.5 * u0 * u0 * h)
        };

        // Early reference point, then a long run: mass must stay put and the
        // drag must hold its value rather than decay.
        for _ in 0..4000 {
            sim.step([0.0, 0.0]);
        }
        let cd_early = cd_at(&sim);
        for _ in 0..16_000 {
            sim.step([0.0, 0.0]);
        }
        let rho_end = mean_rho(&sim);
        let cd_late = cd_at(&sim);

        assert!(
            (rho_end - rho_start).abs() < 0.005,
            "domain mass drifted: rho {rho_start} -> {rho_end}"
        );
        assert!(
            (cd_late - cd_early).abs() / cd_early < 0.10,
            "drag not stationary: Cd {cd_early:.3} -> {cd_late:.3}"
        );
        assert!(
            (1.7..2.6).contains(&cd_late),
            "settled Cd(Re=20) = {cd_late:.3} outside published band"
        );
    }
}
