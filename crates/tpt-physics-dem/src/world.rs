//! Minimal time-stepping DEM driver.

use crate::broadphase::SpatialHash;
use crate::contact::{contact_force, hertz_normal_force, reduced_modulus};
use crate::obstacle::Obstacle;
use crate::particle::Particle;
use rayon::prelude::*;

/// A granular system of spherical particles under gravity, with a planar floor
/// at `y = floor_y`, fixed [`Obstacle`] boundaries, and pairwise
/// Hertz–Mindlin contacts.
pub struct World {
    /// The particles.
    pub particles: Vec<Particle>,
    /// Time step.
    pub dt: f64,
    /// Gravitational acceleration `(gx, gy, gz)`.
    pub gravity: [f64; 3],
    /// Reduced (harmonic) contact modulus `E*`.
    pub e_star: f64,
    /// Coefficient of friction `μ`.
    pub friction: f64,
    /// Coefficient of restitution `e` (sets contact damping).
    pub restitution: f64,
    /// Height of the floor plane.
    pub floor_y: f64,
    /// Fixed obstacle boundaries (cylinders, walls, containers).
    pub obstacles: Vec<Obstacle>,
    /// Extra downward body acceleration (modelling a fluidized/driven flow such
    /// as wet concrete being poured through a cage). Zero by default.
    pub fluidization: f64,
    /// Maximum speed magnitude (m/s). A positive value clamps particle speed
    /// after each step as a stability guard against rare wedged-contact
    /// explosions; `0.0` disables clamping.
    pub max_speed: f64,
    /// Viscous drag coefficient (1/s). A positive value applies a velocity-
    /// proportional damping `v *= 1/(1 + drag·dt)` each step, modelling
    /// particle–fluid drag and guaranteeing asymptotic settling of poured beds
    /// (it removes residual agitation that an explicit contact solver with a
    /// speed clamp would otherwise sustain). `0.0` disables drag.
    pub drag: f64,
}

impl World {
    /// Build a world with steel-like contact stiffness (`E = 200 GPa`, `ν = 0.3`)
    /// and a downward gravity of `9.81`.
    pub fn new(particles: Vec<Particle>, dt: f64) -> Self {
        let e_star = reduced_modulus(200e9, 0.3, 200e9, 0.3);
        World {
            particles,
            dt,
            gravity: [0.0, -9.81, 0.0],
            e_star,
            friction: 0.5,
            restitution: 0.2,
            floor_y: 0.0,
            obstacles: Vec::new(),
            fluidization: 0.0,
            max_speed: 0.0,
            drag: 0.0,
        }
    }

    /// Build a world and attach `obstacles`.
    pub fn with_obstacles(particles: Vec<Particle>, dt: f64, obstacles: Vec<Obstacle>) -> Self {
        World {
            obstacles,
            ..Self::new(particles, dt)
        }
    }

    /// Advance the simulation by one time step (semi-implicit Euler).
    pub fn step(&mut self) {
        let n = self.particles.len();
        let mut force = vec![[0.0_f64; 3]; n];

        for i in 0..n {
            for k in 0..3 {
                force[i][k] += self.gravity[k] * self.particles[i].mass;
            }
            // Optional fluidized-driving term (extra downward body force).
            force[i][1] += self.fluidization * self.particles[i].mass;
        }

        let max_r = self
            .particles
            .iter()
            .map(|p| p.radius)
            .fold(0.0_f64, f64::max)
            .max(1e-6);
        let hash = SpatialHash::build(&self.particles, 2.0 * max_r);
        for (i, j) in hash.candidate_pairs() {
            let f = contact_force(
                &self.particles[i],
                &self.particles[j],
                self.e_star,
                self.friction,
                self.restitution,
            );
            for k in 0..3 {
                force[i][k] += f[k];
                force[j][k] -= f[k];
            }
        }

        for i in 0..n {
            let pen = self.floor_y - (self.particles[i].position[1] - self.particles[i].radius);
            if pen > 0.0 {
                let f_n = hertz_normal_force(self.e_star, self.particles[i].radius, pen);
                let vn = self.particles[i].velocity[1];
                let damp = -2.0 * (self.e_star * self.particles[i].radius * self.particles[i].mass).sqrt()
                    * self.restitution
                    * vn.min(0.0);
                // `damp` is already the upward damping force magnitude (it is
                // positive when the particle moves downward into the floor), so
                // it must be *added* to the Hertz normal force to dissipate.
                force[i][1] += f_n + damp;
                // Inelastic floor: remove any remaining downward (into-floor)
                // velocity so particles come to rest on the floor instead of
                // bouncing — the damping term alone is sub-critical for an open
                // pile and leaves a perpetually-agitated layer.
                if self.particles[i].velocity[1] < 0.0 {
                    self.particles[i].velocity[1] = 0.0;
                }
            }
        }

        // Fixed obstacle boundaries (cylinders, walls, containers).
        for i in 0..n {
            for obs in &self.obstacles {
                if let Some((f, corr)) = obs.resolve(
                    &self.particles[i],
                    self.e_star,
                    self.friction,
                    self.restitution,
                ) {
                    for k in 0..3 {
                        force[i][k] += f[k];
                        self.particles[i].position[k] += corr[k];
                    }
                    // A fixed obstacle is inelastic: after the de-penetration
                    // snap, remove any remaining inward normal velocity so the
                    // particle settles against it instead of limit-cycling
                    // (position-only correction leaves a residual inward
                    // velocity that gravity re-injects every step).
                    let cl = corr[0] * corr[0] + corr[1] * corr[1] + corr[2] * corr[2];
                    if cl > 1e-12 {
                        let inv = 1.0 / cl.sqrt();
                        let (nx, ny, nz) = (corr[0] * inv, corr[1] * inv, corr[2] * inv);
                        let vn = self.particles[i].velocity[0] * nx
                            + self.particles[i].velocity[1] * ny
                            + self.particles[i].velocity[2] * nz;
                        if vn < 0.0 {
                            self.particles[i].velocity[0] -= vn * nx;
                            self.particles[i].velocity[1] -= vn * ny;
                            self.particles[i].velocity[2] -= vn * nz;
                        }
                    }
                }
            }
        }

        let dt = self.dt;
        let vmax = self.max_speed;
        let drag_f = 1.0 / (1.0 + self.drag * dt);
        for i in 0..n {
            let inv = self.particles[i].inv_mass();
            for k in 0..3 {
                let v = (self.particles[i].velocity[k] + force[i][k] * inv * dt) * drag_f;
                self.particles[i].velocity[k] = v;
                self.particles[i].position[k] += v * dt;
            }
            if vmax > 0.0 {
                let v = &self.particles[i].velocity;
                let speed = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                if speed > vmax {
                    let s = vmax / speed;
                    for k in 0..3 {
                        self.particles[i].velocity[k] *= s;
                    }
                }
            }
        }
    }

    /// Advance the simulation by one time step using a `rayon`-parallel
    /// contact sweep, suitable for very large particle counts (≫ 100k).
    ///
    /// Each particle's force is computed independently by summing its
    /// neighbour contacts, so there are no cross-particle write races and the
    /// per-particle loop is distributed across the global thread pool. Floor
    /// and gravity handling are included; fixed obstacles contribute their
    /// contact force but *not* the positional de-penetration (which would be a
    /// race) — for obstacle-laden scenes use the sequential [`World::step`].
    pub fn step_par(&mut self) {
        let n = self.particles.len();
        let max_r = self
            .particles
            .iter()
            .map(|p| p.radius)
            .fold(0.0_f64, f64::max)
            .max(1e-6);
        let hash = SpatialHash::build(&self.particles, 2.0 * max_r);
        let neighbours = hash.neighbour_lists();

        let e_star = self.e_star;
        let mu = self.friction;
        let rest = self.restitution;
        let grav = self.gravity;
        let floor_y = self.floor_y;
        let fluid = self.fluidization;
        let obstacles = &self.obstacles;

        let forces: Vec<[f64; 3]> = (0..n)
            .into_par_iter()
            .map(|i| {
                let pi = &self.particles[i];
                let mut f = [0.0_f64; 3];
                for k in 0..3 {
                    f[k] += grav[k] * pi.mass;
                }
                f[1] += fluid * pi.mass;
                // Floor.
                let pen = floor_y - (pi.position[1] - pi.radius);
                if pen > 0.0 {
                    let f_n = hertz_normal_force(e_star, pi.radius, pen);
                    let vn = pi.velocity[1];
                    let damp = -2.0
                        * (e_star * pi.radius * pi.mass).sqrt()
                        * rest
                        * vn.min(0.0);
                    f[1] += f_n + damp;
                }
                // Pairwise contacts (force on i from each neighbour j).
                for &j in &neighbours[i] {
                    if j == i {
                        continue;
                    }
                    let pj = &self.particles[j];
                    let cf = contact_force(pi, pj, e_star, mu, rest);
                    for k in 0..3 {
                        f[k] += cf[k];
                    }
                }
                // Fixed obstacles (force only; no de-penetration).
                for obs in obstacles {
                    if let Some((of, _corr)) = obs.resolve(pi, e_star, mu, rest) {
                        for k in 0..3 {
                            f[k] += of[k];
                        }
                    }
                }
                f
            })
            .collect();

        let dt = self.dt;
        let vmax = self.max_speed;
        let drag_f = 1.0 / (1.0 + self.drag * dt);
        for i in 0..n {
            let inv = self.particles[i].inv_mass();
            for k in 0..3 {
                let v = (self.particles[i].velocity[k] + forces[i][k] * inv * dt) * drag_f;
                self.particles[i].velocity[k] = v;
                self.particles[i].position[k] += v * dt;
            }
            if vmax > 0.0 {
                let v = &self.particles[i].velocity;
                let speed = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                if speed > vmax {
                    let s = vmax / speed;
                    for k in 0..3 {
                        self.particles[i].velocity[k] *= s;
                    }
                }
            }
        }
    }

    /// Resolve initial overlaps via position-based dynamics (no velocity change),
    /// so a randomly-placed cloud starts from a feasible near-contact state and
    /// the dynamic phase settles instead of rattling at the speed clamp.
    ///
    /// Each iteration pushes every overlapping pair apart by half the penetration
    /// and lifts every floor-/obstacle-penetrating particle out of the boundary.
    /// Because positions only move and velocities are untouched, no kinetic energy
    /// is injected — this is an initial-condition *heal*, not a dynamic step.
    pub fn relax(&mut self, iterations: usize) {
        let max_r = self
            .particles
            .iter()
            .map(|p| p.radius)
            .fold(0.0_f64, f64::max)
            .max(1e-6);
        let e_star = self.e_star;
        let mu = self.friction;
        let rest = self.restitution;
        for _ in 0..iterations {
            let hash = SpatialHash::build(&self.particles, 2.0 * max_r);
            for (i, j) in hash.candidate_pairs() {
                let dx = self.particles[i].position[0] - self.particles[j].position[0];
                let dy = self.particles[i].position[1] - self.particles[j].position[1];
                let dz = self.particles[i].position[2] - self.particles[j].position[2];
                let d = (dx * dx + dy * dy + dz * dz).sqrt();
                let rsum = self.particles[i].radius + self.particles[j].radius;
                if d >= rsum {
                    continue;
                }
                let dd = d.max(1e-9);
                let pen = rsum - d;
                let corr = 0.5 * pen;
                let nx = dx / dd;
                let ny = dy / dd;
                let nz = dz / dd;
                self.particles[i].position[0] -= corr * nx;
                self.particles[i].position[1] -= corr * ny;
                self.particles[i].position[2] -= corr * nz;
                self.particles[j].position[0] += corr * nx;
                self.particles[j].position[1] += corr * ny;
                self.particles[j].position[2] += corr * nz;
            }
            for p in &mut self.particles {
                let pen = self.floor_y - (p.position[1] - p.radius);
                if pen > 0.0 {
                    p.position[1] += pen;
                }
            }
            for i in 0..self.particles.len() {
                for obs in &self.obstacles {
                    if let Some((_f, corr)) =
                        obs.resolve(&self.particles[i], e_star, mu, rest)
                    {
                        self.particles[i].position[0] += corr[0];
                        self.particles[i].position[1] += corr[1];
                        self.particles[i].position[2] += corr[2];
                    }
                }
            }
        }
    }

    /// Total kinetic energy of the system.
    pub fn kinetic_energy(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| {
                0.5 * p.mass
                    * (p.velocity[0].powi(2) + p.velocity[1].powi(2) + p.velocity[2].powi(2))
            })
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn particles_settle_without_blowup() {
        let ps = vec![
            Particle::new([0.0, 1.0, 0.0], [0.0; 3], 0.5, 1000.0),
            Particle::new([0.2, 2.0, 0.0], [0.0; 3], 0.5, 1000.0),
        ];
        let mut w = World::new(ps, 2e-4);
        for _ in 0..3000 {
            w.step();
        }
        // No NaN, particles rest on/above the floor, energy bounded.
        for p in &w.particles {
            assert!(p.position.iter().all(|v| v.is_finite()));
            assert!(p.position[1] >= w.floor_y - 1e-6);
        }
        assert!(w.kinetic_energy() < 10.0, "KE = {}", w.kinetic_energy());
    }

    #[test]
    fn overlapping_particles_repel() {
        // Two particles overlapping should push apart over a few steps.
        let ps = vec![
            Particle::new([0.0, 0.6, 0.0], [0.0; 3], 0.5, 1000.0),
            Particle::new([0.6, 0.6, 0.0], [0.0; 3], 0.5, 1000.0),
        ];
        let mut w = World::new(ps, 1e-4);
        for _ in 0..50 {
            w.step();
        }
        let dist = (w.particles[0].position[0] - w.particles[1].position[0]).abs();
        assert!(dist > 0.2, "particles separated to {dist}");
    }
}
