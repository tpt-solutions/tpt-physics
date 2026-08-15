//! Minimal time-stepping DEM driver.

use crate::broadphase::SpatialHash;
use crate::contact::{contact_force, hertz_normal_force, reduced_modulus};
use crate::particle::Particle;

/// A granular system of spherical particles under gravity, with a planar floor
/// at `y = floor_y` and pairwise Hertz–Mindlin contacts.
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
            }
        }

        let dt = self.dt;
        for i in 0..n {
            let inv = self.particles[i].inv_mass();
            for k in 0..3 {
                self.particles[i].velocity[k] += force[i][k] * inv * dt;
                self.particles[i].position[k] += self.particles[i].velocity[k] * dt;
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
