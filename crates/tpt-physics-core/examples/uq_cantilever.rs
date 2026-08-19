//! Monte-Carlo uncertainty quantification of a cantilever beam.
//!
//! Run with:
//!
//! ```text
//! cargo run -p tpt-physics-core --features uq --example uq_cantilever
//! ```
//!
//! We take a nominal "Structural Steel", scatter its Young's modulus by ±10%
//! (a typical datasheet spread between heats) and its density by ±5%, then
//! Monte-Carlo-sweep two forward models:
//!
//! * tip deflection `δ = P·L³ / (3·E·I)` — depends only on `E`;
//! * fundamental natural frequency `f₁` — depends on *both* `E` and `ρ`.
//!
//! The sweep prints the resulting response statistics, demonstrating how
//! material scatter propagates into engineering quantities of interest.

use tpt_physics_core::material::Material;
use tpt_physics_core::uq::{
    cantilever_natural_frequency, cantilever_tip_deflection, monte_carlo, tol_band, Statistics,
};

fn main() {
    let steel = Material::new("Structural Steel", 200e9, 0.30, 7850.0, 12e-6);

    // Geometry of the (known) cantilever cross-section and span.
    let load = 1000.0; // N, end load
    let length = 2.0; // m
    let inertia = 1e-6; // m⁴, second moment of area
    let area = 1e-3; // m², cross-sectional area

    // ±10% on E, ±5% on ρ, no scatter on ν / α.
    let strategy = tol_band(&steel, 0.10, 0.0, 0.05, 0.0);

    let n = 20_000;

    let deflection = monte_carlo(&strategy, |m| {
        cantilever_tip_deflection(m, load, length, inertia)
    }, n);
    let frequency = monte_carlo(&strategy, |m| {
        cantilever_natural_frequency(m, length, inertia, area)
    }, n);

    println!("Monte-Carlo UQ — Structural Steel cantilever (n = {n})");
    println!("Scatter: ±10% E, ±5% ρ\n");
    print_stats("Tip deflection δ (m)", deflection);
    print_stats("Natural frequency f₁ (Hz)", frequency);

    // Sanity: a narrower band on E should give a lower deflection COV.
    let tight = tol_band(&steel, 0.02, 0.0, 0.0, 0.0);
    let tight_def = monte_carlo(&tight, |m| {
        cantilever_tip_deflection(m, load, length, inertia)
    }, n);
    println!(
        "Deflection COV @ ±2% E = {:.3}%  (vs ±10% E = {:.3}%)",
        tight_def.cov * 100.0,
        deflection.cov * 100.0
    );
    assert!(tight_def.cov < deflection.cov, "tighter band must reduce COV");
}

fn print_stats(label: &str, s: Statistics) {
    println!("{label}");
    println!("  n    = {}", s.n);
    println!("  mean = {:.6e}", s.mean);
    println!("  std  = {:.6e}", s.std);
    println!("  cov  = {:.3}%", s.cov * 100.0);
    println!("  min  = {:.6e}", s.min);
    println!("  p05  = {:.6e}", s.p05);
    println!("  p50  = {:.6e}", s.p50);
    println!("  p95  = {:.6e}", s.p95);
    println!("  max  = {:.6e}\n", s.max);
}
