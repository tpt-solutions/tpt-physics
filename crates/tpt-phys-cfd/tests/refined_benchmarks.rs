//! Refined lid-driven-cavity benchmarks against Ghia, Ghia & Shin (1982).
//!
//! Finer lattices and the Re = 400 case need many steps to reach steady state,
//! so they are `#[ignore]`d; run them with
//! `cargo test --release -p tpt-phys-cfd --test refined_benchmarks -- --ignored`.
//! The coarse Re = 100 check that runs by default lives in `lid_driven_cavity.rs`.

use tpt_phys_cfd::ghia::{self, u_at, v_at};
use tpt_phys_cfd::Lbm2D;

const U_LID: f64 = 0.1;

fn run_to_steady(sim: &mut Lbm2D, max_steps: usize) {
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
#[ignore = "slow steady-state LBM; run with `cargo test --release -- --ignored`"]
fn re100_66x66_matches_ghia() {
    let mut sim = Lbm2D::cavity(66, 100.0, U_LID);
    run_to_steady(&mut sim, 200_000);
    let eu = ghia::u_centerline_error(&sim, U_LID, ghia::U_RE100);
    let ev = ghia::v_centerline_error(&sim, U_LID, ghia::V_RE100);
    eprintln!(
        "Re=100 66²: u rms {:.4} max {:.4} | v rms {:.4} max {:.4}",
        eu.rms, eu.max, ev.rms, ev.max
    );
    assert!(eu.rms < 0.01 && eu.max < 0.02, "u(y) vs Ghia: {eu:?}");
    assert!(ev.rms < 0.01 && ev.max < 0.02, "v(x) vs Ghia: {ev:?}");
    // Extrema of the centreline profiles (Ghia: u_min = -0.2109 @ y = 0.4531).
    assert!((u_at(&sim, U_LID, 0.5, 0.4531) + 0.2109).abs() < 0.01);
    assert!((v_at(&sim, U_LID, 0.8047, 0.5) + 0.2453).abs() < 0.01);
}

#[test]
#[ignore = "slow steady-state LBM; run with `cargo test --release -- --ignored`"]
fn re400_130x130_matches_ghia() {
    let mut sim = Lbm2D::cavity(130, 400.0, U_LID);
    run_to_steady(&mut sim, 200_000);
    let eu = ghia::u_centerline_error(&sim, U_LID, ghia::U_RE400);
    let ev = ghia::v_centerline_error(&sim, U_LID, ghia::V_RE400);
    eprintln!(
        "Re=400 130²: u rms {:.4} max {:.4} | v rms {:.4} max {:.4}",
        eu.rms, eu.max, ev.rms, ev.max
    );
    assert!(eu.rms < 0.01 && eu.max < 0.02, "u(y) vs Ghia: {eu:?}");
    assert!(ev.rms < 0.02 && ev.max < 0.05, "v(x) vs Ghia: {ev:?}");
}
