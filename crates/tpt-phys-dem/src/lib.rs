//! Discrete Element Method for `tpt-physics`.
//!
//! Granular / particle physics with no equivalent in the sibling repos. This
//! crate provides:
//!
//! * [`Particle`] — a rigid sphere with position, velocity, mass and radius;
//! * [`hertz_normal_force`] / [`contact_force`] — the Hertz–Mindlin normal
//!   contact law with Coulomb-capped tangential friction and damping;
//! * [`SpatialHash`] — uniform-grid broad-phase collision detection;
//! * SIMD-accelerated narrow-phase contact resolution (see [`simd`]);
//! * [`World`] — a small time-stepping DEM driver.

pub mod broadphase;
pub mod contact;
pub mod obstacle;
pub mod particle;
pub mod simd;
pub mod world;

pub use obstacle::Obstacle;
