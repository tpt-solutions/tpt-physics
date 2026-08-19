//! Python bindings for `tpt-physics` (thin PyO3 layer over the DEM [`World`]
//! and the material registry).
//!
//! Built with `maturin develop` / `maturin build` from this directory. The
//! crate is intentionally *outside* the Rust workspace (`exclude`d in the root
//! `Cargo.toml`) so the default `cargo build --workspace` does not require
//! Python development headers.
//!
//! The FEA declarative `ProblemSpec` binding was removed when `tpt-physics-fea`
//! was retired (FEM is now delegated to `tpt-fem`); the DEM and material APIs
//! remain.

use pyo3::prelude::*;
use tpt_phys_dem::particle::Particle;
use tpt_phys_dem::world::World;

/// A Python-facing granular world. Wraps [`tpt_physics_dem::world::World`].
#[pyclass]
struct DemWorld {
    inner: World,
}

#[pymethods]
impl DemWorld {
    /// Create an empty world with the given time step (s).
    #[new]
    fn new(dt: f64) -> Self {
        DemWorld {
            inner: World::new(Vec::new(), dt),
        }
    }

    /// Add a spherical particle at `(x, y, z)` with radius `r` and `density`.
    fn add_particle(&mut self, x: f64, y: f64, z: f64, r: f64, density: f64) {
        self.inner
            .particles
            .push(Particle::new([x, y, z], [0.0; 3], r, density));
    }

    /// Advance the simulation by one time step.
    fn step(&mut self) {
        self.inner.step();
    }

    /// Total kinetic energy of the system (J).
    fn kinetic_energy(&self) -> f64 {
        self.inner.kinetic_energy()
    }

    /// Flattened `[x0, y0, z0, x1, y1, z1, ...]` particle positions.
    fn positions(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.inner.particles.len() * 3);
        for p in &self.inner.particles {
            out.extend_from_slice(&p.position);
        }
        out
    }

    /// Number of particles.
    fn n_particles(&self) -> usize {
        self.inner.particles.len()
    }
}

/// The `tpt_physics` Python module.
#[pymodule]
fn tpt_physics_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<DemWorld>()?;
    Ok(())
}
