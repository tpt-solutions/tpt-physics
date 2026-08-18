//! Example-gallery runner: a tiny end-to-end demonstration of every domain
//! this workspace still owns (FEM moved to `tpt-fem`), printed to stdout.
//! Run with:
//!
//! ```text
//! cargo run -p tpt-phys-gallery
//! ```
//!
//! Each `demo_*` function is a self-contained "hello world" for one crate; see
//! the per-crate `examples/` directories for the expanded versions. New
//! multiphysics crates (`tpt-phys-fsi`, `tpt-phys-thermal-struct`,
//! `tpt-phys-electro-thermal`, `tpt-phys-orchestrator`) will gain demos here
//! as they grow past scaffolding.

use tpt_phys_cfd::Lbm2D;
use tpt_phys_dem::particle::Particle;
use tpt_phys_dem::world::World;

fn demo_core() {
    let reg = tpt_phys_core::MaterialRegistry::with_defaults();
    let steel = reg.get("Structural Steel").expect("present");
    println!(
        "  core: Structural Steel  E={:.2e} Pa  G={:.2e} Pa",
        steel.youngs_modulus,
        steel.shear_modulus()
    );
}

fn demo_dem() {
    let ps = vec![
        Particle::new([0.0, 1.0, 0.0], [0.0; 3], 0.5, 1000.0),
        Particle::new([0.2, 2.0, 0.0], [0.0; 3], 0.5, 1000.0),
    ];
    let mut w = World::new(ps, 2e-4);
    for _ in 0..4000 {
        w.step();
    }
    println!(
        "  dem: 2 particles settled, kinetic energy = {:.4e} J",
        w.kinetic_energy()
    );
}

fn demo_cfd() {
    let nx = 48;
    let ny = 24;
    let mut sim = Lbm2D::new(nx, ny, 0.6);
    sim.set_horizontal_walls();
    // Periodic streamwise boundary (the proven stable Poiseuille setup) with a
    // body-force drive; the open (non-reflective) boundary is demonstrated in
    // the `open_boundary_outflow_is_nonreflecting` unit test.
    sim.initialise(1.0, [0.0, 0.0]);
    for _ in 0..20000 {
        sim.step([1e-5, 0.0]);
    }
    let prof = sim.x_velocity_profile();
    let yc = ny / 2;
    println!(
        "  cfd: driven channel centre velocity u = {:.5} (Poiseuille-like profile)",
        prof[yc]
    );
}

fn main() {
    println!("tpt-physics example gallery");
    println!("============================");
    demo_core();
    demo_dem();
    demo_cfd();
    println!("============================");
    println!("All domains demonstrated. See GALLERY.md and each crate's examples/.");
}
