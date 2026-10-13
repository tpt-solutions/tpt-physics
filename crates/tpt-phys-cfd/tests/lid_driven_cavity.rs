//! CFD validation: lid-driven cavity flow against Ghia, Ghia & Shin (1982).
//!
//! A square cavity with three stationary no-slip walls and a top lid moving in
//! `+x` develops a primary recirculating vortex. [`Lbm2D::cavity`] builds the
//! benchmark at a given Reynolds number; [`tpt_phys_cfd::ghia`] samples the
//! lattice at the published centreline coordinates and reports the error.
//!
//! Both tests here run by default (a 34×34 lattice is enough to land within
//! ~1% of the published profiles); the finer Re = 100 and Re = 400 runs live in
//! `refined_benchmarks.rs` and are `#[ignore]`d for speed.

use tpt_phys_cfd::ghia::{self, u_at, v_at};
use tpt_phys_cfd::Lbm2D;

const U_LID: f64 = 0.1;

fn run_to_steady(sim: &mut Lbm2D, max_steps: usize) {
    // Converge on the centreline velocity rather than a fixed step count.
    let mut last = f64::INFINITY;
    for chunk in 0..max_steps / 1000 {
        for _ in 0..1000 {
            sim.step_par([0.0, 0.0]);
        }
        let u = u_at(sim, U_LID, 0.5, 0.25);
        if chunk > 4 && (u - last).abs() < 1e-8 {
            return;
        }
        last = u;
    }
}

#[test]
fn cavity_smoke_primary_vortex_forms() {
    let mut sim = Lbm2D::cavity(34, 100.0, U_LID);
    for _ in 0..8000 {
        sim.step_par([0.0, 0.0]);
    }
    let (rho, ux, uy) = sim.macro_fields();
    assert!(rho.iter().chain(ux).chain(uy).all(|v| v.is_finite()));
    // The lid drags fluid in +x right beneath it ...
    assert!(u_at(&sim, U_LID, 0.5, 0.95) > 0.4);
    // ... and a return flow develops in the lower half (primary vortex).
    assert!(u_at(&sim, U_LID, 0.5, 0.45) < -0.1);
    // Closed box: total mass is conserved (initial fluid mass = fluid nodes).
    let fluid = sim.solid.iter().filter(|&&s| !s).count() as f64;
    let mass: f64 = sim
        .rho
        .iter()
        .zip(&sim.solid)
        .filter(|(_, &s)| !s)
        .map(|(r, _)| r)
        .sum();
    assert!(
        (mass - fluid).abs() / fluid < 1e-9,
        "mass {mass} vs {fluid}"
    );
}

#[test]
fn lid_driven_cavity_re100_matches_ghia() {
    let mut sim = Lbm2D::cavity(34, 100.0, U_LID);
    run_to_steady(&mut sim, 30_000);
    let eu = ghia::u_centerline_error(&sim, U_LID, ghia::U_RE100);
    let ev = ghia::v_centerline_error(&sim, U_LID, ghia::V_RE100);
    assert!(eu.rms < 0.012 && eu.max < 0.02, "u(y) vs Ghia: {eu:?}");
    assert!(ev.rms < 0.012 && ev.max < 0.02, "v(x) vs Ghia: {ev:?}");
    // Primary-vortex signature (Ghia: u_min = -0.2109 @ y = 0.4531).
    assert!((u_at(&sim, U_LID, 0.5, 0.4531) + 0.2109).abs() < 0.015);
    assert!((v_at(&sim, U_LID, 0.8047, 0.5) + 0.2453).abs() < 0.015);
}
