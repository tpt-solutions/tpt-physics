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
//! * a Structure-of-Arrays, `wide::f64x4`-vectorised and Rayon row-parallel
//!   stepping kernel ([`Lbm2D::step`] serial, [`Lbm2D::step_par`] parallel —
//!   bit-identical results),
//! * selectable collision operators ([`Collision`]: BGK, TRT),
//! * lid-driven-cavity reference data and error metrics ([`ghia`]).

pub mod collision;
pub mod ghia;
mod kernel;
pub mod lattice;
pub mod sph;

#[cfg(test)]
mod reference;

pub use collision::Collision;
use lattice::D2Q9;

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
/// (boundary). Distributions are stored Structure-of-Arrays: nine flat
/// `Vec<f64>` of length `nx * ny` (one per lattice velocity), so neighbouring
/// nodes are contiguous for SIMD.
#[derive(Clone)]
pub struct Lbm2D {
    /// Number of nodes in the `x` (streamwise) direction.
    pub nx: usize,
    /// Number of nodes in the `y` (wall-normal) direction.
    pub ny: usize,
    /// BGK relaxation time `τ` (>`0.5` for stability).
    pub tau: f64,
    /// Post-stream distributions `f_k(x)`, indexed `f[k][node]`.
    f: [Vec<f64>; 9],
    /// Post-collision distributions (the streaming source). Solid nodes are
    /// never written and stay zero.
    g: [Vec<f64>; 9],
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
    /// Collision operator.
    collision: Collision,
    /// Fluid force on the solid boundaries during the most recent step,
    /// `[Fx, Fy]`, accumulated by boundary momentum exchange in the collision
    /// pass.
    boundary_force: [f64; 2],
    /// Fluid nodes whose whole in-lattice 3×3 neighbourhood is fluid (SIMD /
    /// bulk-copy fast path). Rebuilt when `solid` no longer equals
    /// `solid_key`.
    clear: Vec<bool>,
    /// Snapshot of `solid` that `clear` was computed from.
    solid_key: Vec<bool>,
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
            f: std::array::from_fn(|_| vec![0.0; n]),
            g: std::array::from_fn(|_| vec![0.0; n]),
            rho: vec![1.0; n],
            ux: vec![0.0; n],
            uy: vec![0.0; n],
            solid: vec![false; n],
            wall_vel: vec![[0.0; 2]; n],
            x_boundary: XBoundary::Periodic,
            collision: Collision::Bgk,
            boundary_force: [0.0; 2],
            clear: Vec::new(),
            solid_key: Vec::new(),
        }
    }

    /// Select the collision operator (default [`Collision::Bgk`]).
    pub fn set_collision(&mut self, collision: Collision) {
        if let Collision::Trt { magic } = collision {
            assert!(magic > 0.0, "TRT magic parameter must be positive");
        }
        self.collision = collision;
    }

    /// The active collision operator.
    pub fn collision(&self) -> Collision {
        self.collision
    }

    /// Distributions of node `i` as `[f_0 … f_8]`.
    #[cfg(test)]
    #[inline]
    fn node_f(&self, i: usize) -> [f64; 9] {
        std::array::from_fn(|k| self.f[k][i])
    }

    /// Overwrite the distributions of node `i`.
    #[inline]
    fn set_node_f(&mut self, i: usize, v: [f64; 9]) {
        for (k, x) in v.into_iter().enumerate() {
            self.f[k][i] = x;
        }
    }

    /// A ready-to-run lid-driven cavity: an `n × n` lattice enclosed by solid
    /// walls with the top wall moving in `+x` at `u_lid`, initialised at rest.
    ///
    /// The fluid occupies `L = n − 2` cells (the border nodes are the walls,
    /// whose half-way bounce-back surfaces sit half a cell inside), so the
    /// Reynolds number is `Re = u_lid · L / ν` and the relaxation time is
    /// `τ = ½ + 3 ν` with `ν = u_lid · L / Re`. Keep `u_lid ≲ 0.1` (Mach ≪ 1)
    /// and `τ` comfortably above `½`; the constructor panics if `Re` is so high
    /// for this `n` that `τ ≤ ½`. Compare against [`ghia`] reference data.
    pub fn cavity(n: usize, re: f64, u_lid: f64) -> Self {
        assert!(re > 0.0 && u_lid > 0.0, "re and u_lid must be positive");
        let nu = u_lid * (n as f64 - 2.0) / re;
        let mut sim = Lbm2D::new(n, n, 0.5 + 3.0 * nu);
        sim.set_box_walls();
        sim.set_moving_lid(n - 1, u_lid);
        sim.initialise(1.0, [0.0, 0.0]);
        sim
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
            self.set_node_f(i, D2Q9::equilibrium(rho0, u0));
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
                self.set_node_f(i, D2Q9::equilibrium(self.rho[i], u));
            }
        }
    }

    /// Kinematic viscosity `ν = cs² (τ − ½)`, with `cs² = 1/3` for D2Q9.
    pub fn viscosity(&self) -> f64 {
        D2Q9::CS2 * (self.tau - 0.5)
    }

    /// Advance the simulation by one collision+streaming step under body-force
    /// density `force = (Fx, Fy)` (acceleration per unit mass, lattice units).
    ///
    /// Single-threaded; see [`Lbm2D::step_par`] for the Rayon version (same
    /// result, bit for bit).
    pub fn step(&mut self, force: [f64; 2]) {
        self.step_impl(force, false);
    }

    /// Like [`Lbm2D::step`] but runs every pass row-parallel on the Rayon
    /// pool. Results are bit-identical to the serial step.
    pub fn step_par(&mut self, force: [f64; 2]) {
        self.step_impl(force, true);
    }

    /// Recompute macroscopic density and velocity (used as a post-step check and
    /// externally after stepping).
    pub fn macro_fields(&mut self) -> (&Vec<f64>, &Vec<f64>, &Vec<f64>) {
        self.moments_impl(false);
        (&self.rho, &self.ux, &self.uy)
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
    pub fn x_velocity_profile(&self) -> Vec<f64> {
        (0..self.ny)
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
                sim.set_node_f(i, D2Q9::equilibrium(1.0, [0.1, 0.0]));
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
