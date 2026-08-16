//! CFD validation: lid-driven cavity flow (Ghia et al. reference).
//!
//! A square cavity with three stationary no-slip walls and a top lid moving in
//! `+x` at speed `U` develops a primary recirculating vortex: the lid drags
//! fluid in `+x` just below it, the fastest streamwise flow sits in the upper
//! half of the cavity, and a return flow (`u_x < 0`) appears near the floor. On
//! the coarse lattice used here the lid-dominated region spans a substantial
//! portion (≈40%) of the centre-line height; the test asserts the qualitative
//! vortex structure rather than an exact Ghia centre-line match.

use tpt_physics_cfd::Lbm2D;

#[test]
fn lid_driven_cavity_primary_vortex() {
    let nx = 48;
    let ny = 48;
    // Relaxation time `τ` sets the kinematic viscosity `ν = cs²(τ−½)`; with
    // `u_lid = 0.1` and `L = 48` this gives `Re = u·L/ν ≈ 720`. On this coarse
    // lattice the lid-driven primary recirculating vortex has its centre in the
    // upper portion of the cavity, so `u_x > 0` over a substantial (≈40%) span
    // of the vertical centre-line just beneath the lid, with a return flow
    // (`u_x < 0`) near the floor.
    let tau = 0.52;
    let u_lid = 0.1;

    let mut sim = Lbm2D::new(nx, ny, tau);
    sim.set_box_walls();
    sim.set_moving_lid(ny - 1, u_lid);
    sim.initialise(1.0, [0.0, 0.0]);

    for _ in 0..250000 {
        sim.step([0.0, 0.0]);
    }

    let col = sim.x_velocity_at_column(nx / 2);
    let interior: Vec<f64> = col[1..ny - 1].to_vec();
    for v in &interior {
        assert!(v.is_finite());
    }

    let beneath_lid = interior[interior.len() - 1];
    let max_u = interior.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let argmax = interior
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .unwrap()
        .0;
    let min_u = interior.iter().cloned().fold(f64::INFINITY, f64::min);
    let positive = interior.iter().filter(|&&v| v > 0.0).count();
    let pos_frac = positive as f64 / interior.len() as f64;
    eprintln!(
        "DBG cavity u_lid={:.3}: beneath={:.4} max={:.4}@{} min={:.4} pos%={:.1} prof={:?}",
        u_lid,
        beneath_lid,
        max_u,
        argmax,
        min_u,
        100.0 * pos_frac,
        interior.iter().rev().step_by(4).collect::<Vec<_>>()
    );

    // Primary-vortex validation (lid-driven cavity, Ghia et al. reference
    // structure): the moving lid drags fluid in `+x` just beneath it, the
    // fastest streamwise flow sits in the upper half of the cavity, and a
    // recirculating vortex forms with `u_x > 0` over a substantial span of the
    // centre-line (the lid-dominated upper region) and a return flow
    // (`u_x < 0`) near the floor.
    assert!(beneath_lid > 0.01, "lid did not drag fluid: {beneath_lid}");
    assert!(
        argmax >= interior.len() / 2,
        "peak velocity not in upper half: @{}",
        argmax
    );
    assert!(
        pos_frac > 0.30,
        "primary vortex too small: u_x>0 over {:.0}% of height",
        100.0 * pos_frac
    );
    assert!(min_u < 0.0, "no return flow at the floor: min_u = {min_u}");
}
