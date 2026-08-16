//! "Spacer Benchmark" case study: end-to-end simulation of the 3D-printed pile
//! cage spacer, timed.
//!
//! This is the Phase-1 milestone wrapped as a benchmark: it wires the full
//! `tpt-physics` stack (material database → CAD ingestion adapter → volume tet
//! mesh → linear elasticity with self-weight) and reports the solve time and the
//! physics result (fixed base stays put, free top compresses under its own
//! weight, magnitude consistent with `ρ g h² / E`).

use std::time::Instant;

use tpt_fem_elasticity::{solve_elasticity, ElasticModel};
use tpt_fem_mesh_gen::box_mesh;
use tpt_physics_core::cad::{CadIngestor, CadSolid, CadVertex};
use tpt_physics_core::MaterialRegistry;

fn main() {
    // 1. Material — PLA from the net-new material database.
    let registry = MaterialRegistry::with_defaults();
    let material = registry
        .get("PLA (3D print, ~50% infill)")
        .expect("PLA present in default registry");
    let (young, poisson, density) = (
        material.youngs_modulus,
        material.poissons_ratio,
        material.density,
    );
    println!(
        "Material: {}  E={:.3e} Pa  nu={:.3}  rho={:.1} kg/m^3",
        material.name, young, poisson, density
    );

    // 2. CAD ingestion adapter: spacer footprint as a tpt-cad solid lowered to a
    //    tpt-fem-mesh surface mesh.
    let mut cad = CadSolid::new();
    for corner in [
        [0.0, 0.0, 0.0],
        [0.04, 0.0, 0.0],
        [0.04, 0.05, 0.0],
        [0.0, 0.05, 0.0],
        [0.0, 0.0, 0.04],
        [0.04, 0.0, 0.04],
        [0.04, 0.05, 0.04],
        [0.0, 0.05, 0.04],
    ] {
        cad.add_vertex(CadVertex::new(corner[0], corner[1], corner[2]));
    }
    cad.add_face(vec![0, 1, 2, 3], Some(1));
    cad.add_face(vec![4, 5, 6, 7], Some(2));
    cad.add_face(vec![0, 3, 7, 4], Some(3));
    cad.add_face(vec![1, 2, 6, 5], Some(4));
    cad.add_face(vec![0, 1, 5, 4], Some(5));
    cad.add_face(vec![3, 2, 6, 7], Some(6));
    let surface = cad.ingest().expect("ingest CAD solid");
    println!(
        "CAD adapter → surface mesh: {} nodes, {} tris",
        surface.node_count(),
        surface.element_count()
    );

    // 3. Volume tet mesh (40×50×40 mm).
    let (w, h, d) = (0.04_f64, 0.05_f64, 0.04_f64);
    let mesh = box_mesh([0.0, 0.0, 0.0], [w, h, d], [4, 5, 4]);
    let n_nodes = mesh.node_count();
    let n_elems = mesh.element_count();
    println!("Volume tet mesh: {n_nodes} nodes, {n_elems} elements");

    // 4. Boundary conditions: fix the bottom face (y ≈ 0).
    let mut dirichlet = Vec::new();
    for n in 0..n_nodes {
        let c = mesh.node_coords(n);
        if c[1].abs() < 1e-9 {
            for comp in 0..3 {
                dirichlet.push((n * 3 + comp, 0.0));
            }
        }
    }
    println!("Fixed base nodes: {}", dirichlet.len() / 3);

    // Self-weight body force.
    let g = 9.81;
    let body = |_x: &[f64]| vec![0.0, -density * g, 0.0];

    // 5. Solve, timed.
    let t0 = Instant::now();
    let u = solve_elasticity(
        &mesh,
        ElasticModel::Continuum3D,
        young,
        poisson,
        2,
        body,
        &dirichlet,
    )
    .expect("elasticity solve");
    let solve_time = t0.elapsed();

    // 6. Post-process.
    let mut bottom_disp = 0.0_f64;
    let mut top_node = 0usize;
    let mut top_y = -1e9_f64;
    let mut max_disp = 0.0_f64;
    for n in 0..n_nodes {
        let c = mesh.node_coords(n);
        let mag = (u[n * 3].powi(2) + u[n * 3 + 1].powi(2) + u[n * 3 + 2].powi(2)).sqrt();
        max_disp = max_disp.max(mag);
        if c[1].abs() < 1e-9 {
            bottom_disp = bottom_disp.max(mag);
        }
        if c[1] > top_y {
            top_y = c[1];
            top_node = n;
        }
    }
    let top_disp_y = u[top_node * 3 + 1];
    let top_disp_mag =
        (u[top_node * 3].powi(2) + u[top_node * 3 + 1].powi(2) + u[top_node * 3 + 2].powi(2))
            .sqrt();
    let strain_est = density * g * h / young;
    let settle_est = strain_est * h;

    println!("--- Results ---");
    println!("Solve time            : {solve_time:?}");
    println!("Max |displacement|     : {max_disp:.3e} m");
    println!("Bottom (fixed) max |u| : {bottom_disp:.3e} m");
    println!("Top-centre disp (y)    : {top_disp_y:.3e} m  (mag {top_disp_mag:.3e})");
    println!("Settlement estimate ρgh²/E : {settle_est:.3e} m");

    assert!(bottom_disp < 1e-12, "fixed base must not move");
    assert!(top_disp_y < 0.0, "free top must compress downward");
    assert!(top_disp_y.abs() < settle_est * 10.0 && top_disp_y.abs() > settle_est * 0.01);
    println!("\nSpacer Benchmark PASSED.");
}
