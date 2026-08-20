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
//! `tpt-phys-electro-thermal`, `tpt-phys-orchestrator`) each gain a demo here.

use tpt_phys_cfd::Lbm2D;
use tpt_phys_dem::particle::Particle;
use tpt_phys_dem::world::World;
use tpt_phys_fsi::StructuralModel;

fn demo_core() {
    let reg = tpt_phys_core::MaterialRegistry::with_defaults();
    let steel = reg
        .get("Structural Steel")
        .expect("MaterialRegistry::with_defaults() must seed \"Structural Steel\"");
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

fn demo_fsi() {
    // Partitioned FSI: a driven channel whose right wall is the fluid–structure
    // interface, coupled to a lumped structural oscillator.
    let nx = 32;
    let ny = 16;
    let mut sim = Lbm2D::new(nx, ny, 0.6);
    sim.set_horizontal_walls();
    sim.add_rect(nx - 2, 1, nx - 2, ny - 2);
    sim.initialise(1.0, [0.1, 0.0]); // steady rightward flow

    let mut b = tpt_fem_mesh::MeshBuilder::new();
    b.add_node(vec![0.0, 0.0, 0.0]);
    let mesh = b.build();
    let mut structure = tpt_phys_fsi::LumpedStructure::new(1, 1.0, 10.0, 2.0);
    let mut driver = tpt_phys_fsi::FsiDriver::new(&sim, &structure, &mesh);
    for _ in 0..200 {
        driver.step(&mut sim, &mut structure, 1.0);
    }
    let disp = structure.displacement(0)[0];
    println!(
        "  fsi: structure pushed downstream by flow, displacement x = {:.3e} m, energy = {:.3e} J (lumped model)",
        disp,
        structure.energy()
    );
}

fn demo_thermal_struct() {
    // Assemble the thermal-strain load vector for a single tet with a uniform
    // temperature rise — the coupling this crate ports from the old FEA crate.
    let mut b = tpt_fem_mesh::MeshBuilder::new();
    let n0 = b.add_node(vec![0.0, 0.0, 0.0]);
    let n1 = b.add_node(vec![1.0, 0.0, 0.0]);
    let n2 = b.add_node(vec![0.0, 1.0, 0.0]);
    let n3 = b.add_node(vec![0.0, 0.0, 1.0]);
    b.add_element(tpt_fem_mesh::CellType::Tet, vec![n0, n1, n2, n3]);
    let mesh = b.build();
    let steel = tpt_phys_core::Material::new("Steel", 200e9, 0.3, 7850.0, 12e-6);
    let temps = vec![50.0, 50.0, 50.0, 50.0];
    let load = tpt_phys_thermal_struct::thermal_load_vector(&mesh, 3, &steel, &temps, 20.0);
    let norm: f64 = load.iter().map(|v| v * v).sum::<f64>().sqrt();
    println!(
        "  thermal-struct: uniform-temp load-vector norm = {:.4e} (zero by self-equilibrium)",
        norm
    );
}

fn demo_electro_thermal() {
    // A resistive rod heating under an applied voltage (Joule heating).
    let mut rod = tpt_phys_electro_thermal::ElectroThermalRod::new(21, 300.0);
    rod.dx = 0.01;
    rod.set_voltage(10.0);
    rod.convection = 50.0;
    let t0 = rod.temperatures()[10];
    for _ in 0..2000 {
        rod.step(1e-4);
    }
    let t1 = rod.temperatures()[10];
    println!(
        "  electro-thermal: rod centre heats under load {:.1} K -> {:.1} K (Joule + self-limiting σ(T))",
        t0, t1
    );
}

fn demo_orchestrator() {
    // Multiphysics co-simulation: electro-thermal -> thermal-structural, plus
    // the FSI sub-model, driven by tpt-sci-sim-core's Simulation engine.
    let mut sim = tpt_phys_orchestrator::build_demo_simulation();
    for _ in 0..200 {
        sim.step(1e-4);
    }
    let mut buf = vec![0.0; sim.submodel(0).state_dim()];
    sim.submodel(0).gather_state(&mut buf);
    println!(
        "  orchestrator: demo Simulation stepped, electro-thermal node[0] T = {:.1} K",
        buf[0]
    );
}

fn main() {
    println!("tpt-physics example gallery");
    println!("============================");
    demo_core();
    demo_dem();
    demo_cfd();
    demo_fsi();
    demo_thermal_struct();
    demo_electro_thermal();
    demo_orchestrator();
    println!("============================");
    println!("All domains demonstrated. See GALLERY.md and each crate's examples/.");
}
