//! Hello-world for `tpt-phys-cfd`: the lid-driven cavity benchmark.
//!
//! A fully-enclosed square cavity whose top lid moves in `+x` drives a
//! primary recirculating vortex. This runs the D2Q9 LBM solver to steady
//! state and compares the centreline velocity profiles with Ghia, Ghia &
//! Shin (1982).
//!
//! Run with: `cargo run --example cavity --release [n] [re]`

use std::time::Instant;
use tpt_phys_cfd::{ghia, Lbm2D};

fn main() {
    let mut args = std::env::args().skip(1);
    let n: usize = args.next().and_then(|a| a.parse().ok()).unwrap_or(66);
    let re: f64 = args.next().and_then(|a| a.parse().ok()).unwrap_or(100.0);
    let u_lid = 0.1;

    let mut lbm = Lbm2D::cavity(n, re, u_lid);
    let steps = 40 * n * n; // several diffusion times
    let t0 = Instant::now();
    for _ in 0..steps {
        lbm.step_par([0.0, 0.0]);
    }
    let dt = t0.elapsed().as_secs_f64();
    let mlups = (steps * n * n) as f64 / dt / 1e6;

    let (ru, rv) = if re < 250.0 {
        (ghia::U_RE100, ghia::V_RE100)
    } else {
        (ghia::U_RE400, ghia::V_RE400)
    };
    let eu = ghia::u_centerline_error(&lbm, u_lid, ru);
    let ev = ghia::v_centerline_error(&lbm, u_lid, rv);

    println!(
        "Lid-driven cavity ({n}x{n}, Re={re}, tau={:.4}, u_lid={u_lid})",
        lbm.tau
    );
    println!("  steps / speed        : {steps} steps, {mlups:.1} MLUPS");
    println!(
        "  u(y) at x=0.5 vs Ghia: rms {:.4}, max {:.4}",
        eu.rms, eu.max
    );
    println!(
        "  v(x) at y=0.5 vs Ghia: rms {:.4}, max {:.4}",
        ev.rms, ev.max
    );
    println!("  y        u_lbm    u_ghia");
    for &(y, r) in ru {
        println!("  {y:.4}  {:8.4}  {r:8.4}", ghia::u_at(&lbm, u_lid, 0.5, y));
    }
}
