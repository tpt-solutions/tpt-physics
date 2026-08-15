//! Rigid-sphere particle state.

/// A single spherical DEM particle.
#[derive(Debug, Clone)]
pub struct Particle {
    /// Position `(x, y, z)`.
    pub position: [f64; 3],
    /// Velocity `(vx, vy, vz)`.
    pub velocity: [f64; 3],
    /// Radius.
    pub radius: f64,
    /// Mass.
    pub mass: f64,
}

impl Particle {
    /// Construct a particle. Mass is derived from radius and density if not
    /// supplied explicitly.
    pub fn new(position: [f64; 3], velocity: [f64; 3], radius: f64, density: f64) -> Self {
        let mass = density * (4.0 / 3.0) * std::f64::consts::PI * radius * radius * radius;
        Particle {
            position,
            velocity,
            radius,
            mass,
        }
    }

    /// Inverse mass (zero for an immovable particle).
    pub fn inv_mass(&self) -> f64 {
        if self.mass > 0.0 {
            1.0 / self.mass
        } else {
            0.0
        }
    }
}
