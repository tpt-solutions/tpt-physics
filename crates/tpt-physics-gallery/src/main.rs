//! Example-gallery runner: a tiny end-to-end demonstration of every physics
//! domain in the workspace, printed to stdout. Run with:
//!
//! ```text
//! cargo run -p tpt-physics-gallery
//! ```
//!
//! Each `demo_*` function is a self-contained "hello world" for one crate; see
//! the per-crate `examples/` directories for the expanded versions.

use tpt_physics_core::Material;
use tpt_physics_dem::particle::Particle;
use tpt_physics_dem::world::World;
use tpt_physics_fea::spec::{DomainSpec, LoadSpec, ProblemSpec, SolverSpec};
use tpt_physics_cfd::Lbm2D;

fn demo_core() {
    let reg = tpt_physics_core::MaterialRegistry::with_defaults();
    let steel = reg.get("Structural Steel").expect("present");
    println!(
        "  core: Structural Steel  E={:.2e} Pa  G={:.2e} Pa",
        steel.youngs_modulus,
        steel.shear_modulus()
    );
}

fn demo_solver() {
    use tpt_physics_solver::cg::cg;
    use tpt_physics_solver::linalg::csr_from_dense;
    let a = csr_from_dense(2, 2, &[4.0, -1.0, -1.0, 4.0]);
    let (x, rep) = cg(&a, &[3.0, 3.0], None, 1e-9, 100).unwrap();
    println!(
        "  solver: CG solved 2x2 SPD in {} iters, x = [{:.3}, {:.3}]",
        rep.iterations, x[0], x[1]
    );
}

fn demo_fea() {
    let spec = ProblemSpec {
        materials: None,
        material: tpt_physics_fea::spec::MaterialRef::Inline(Material::new(
            "PLA", 3.5e9, 0.36, 1240.0, 68e-6,
        )),
        domain: DomainSpec::Box {
            min: [0.0, 0.0, 0.0],
            max: [0.04, 0.05, 0.04],
            n: [4, 5, 4],
        },
        boundary_conditions: tpt_physics_fea::spec::BcSpec {
            fixed_planes: vec!["y_min".to_string()],
        },
        loads: LoadSpec {
            self_weight: true,
            gravity: 9.81,
        },
        solver: SolverSpec::StaticLinear,
    };
    let solved = spec.solve(&tpt_physics_core::MaterialRegistry::new()).unwrap();
    println!(
        "  fea: box self-weight → top settles {:.3e} m ({} nodes)",
        solved.free_top_settlement_y, solved.n_nodes
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

fn demo_ai() {
    use tpt_physics_ai::{DifferentiablePlant, HarmonicOscillator};
    let plant = HarmonicOscillator::new();
    let (ds, _da) = plant.jacobians(&[0.5, 0.25], &[0.1]);
    // Forward-mode AD Jacobian: dx'/dx = 1 - dt^2 k/m.
    let expected = 1.0 - plant.dt * plant.dt * plant.k / plant.m;
    println!(
        "  ai: HarmonicOscillator ∂x'/∂x = {:.4} (autodiff, expect {:.4})",
        ds[0][0], expected
    );
}

fn main() {
    println!("tpt-physics example gallery");
    println!("============================");
    demo_core();
    demo_solver();
    demo_fea();
    demo_dem();
    demo_cfd();
    demo_ai();
    println!("============================");
    println!("All domains demonstrated. See GALLERY.md and each crate's examples/.");
}
