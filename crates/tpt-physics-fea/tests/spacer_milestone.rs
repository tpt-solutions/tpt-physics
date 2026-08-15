//! Milestone integration test: simulate the 3D-printed pile cage spacer.
//!
//! Exercises the full `tpt-physics` stack — material database (core), CAD
//! ingestion adapter (core), volume mesh generation (`tpt-fem-mesh-gen`), and
//! static linear FEA (`tpt-fem-elasticity`) — and checks the solution is
//! physical: the fixed base stays put while the free top compresses downward
//! under self-weight.

use tpt_fem_elasticity::{solve_elasticity, ElasticModel};
use tpt_fem_mesh_gen::box_mesh;
use tpt_physics_core::cad::{CadIngestor, CadSolid, CadVertex};
use tpt_physics_core::MaterialRegistry;

#[test]
fn pile_cage_spacer_static_solve_is_physical() {
    // Material from the registry.
    let registry = MaterialRegistry::with_defaults();
    let material = registry
        .get("PLA (3D print, ~50% infill)")
        .expect("PLA present");
    let (young, poisson, density) =
        (material.youngs_modulus, material.poissons_ratio, material.density);

    // CAD ingestion adapter produces a valid surface mesh of the footprint.
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
    cad.add_face(vec![0, 1, 2, 3], None);
    cad.add_face(vec![4, 5, 6, 7], None);
    cad.add_face(vec![0, 3, 7, 4], None);
    cad.add_face(vec![1, 2, 6, 5], None);
    cad.add_face(vec![0, 1, 5, 4], None);
    cad.add_face(vec![3, 2, 6, 7], None);
    let surface = cad.ingest().expect("ingest CAD solid");
    assert_eq!(surface.node_count(), 8);
    assert_eq!(surface.element_count(), 12);

    // Volume tet mesh (vertical axis = y, 50 mm).
    let (w, h, d) = (0.04_f64, 0.05_f64, 0.04_f64);
    let mesh = box_mesh([0.0, 0.0, 0.0], [w, h, d], [4, 5, 4]);
    let n_nodes = mesh.node_count();
    assert!(n_nodes > 0 && mesh.element_count() > 0);

    // Fix the bottom face (y ≈ 0) in all 3 DOFs.
    let mut dirichlet = Vec::new();
    for n in 0..n_nodes {
        let c = mesh.node_coords(n);
        if c[1].abs() < 1e-9 {
            for comp in 0..3 {
                dirichlet.push((n * 3 + comp, 0.0));
            }
        }
    }

    // Self-weight (gravity) body force.
    let g = 9.81;
    let body = |_x: &[f64]| vec![0.0, -density * g, 0.0];

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

    // Post-process.
    let mut bottom_max = 0.0_f64;
    let mut top_node = 0usize;
    let mut top_y = -1e9_f64;
    let mut max_disp = 0.0_f64;
    for n in 0..n_nodes {
        let c = mesh.node_coords(n);
        let mag = (u[n * 3].powi(2) + u[n * 3 + 1].powi(2) + u[n * 3 + 2].powi(2)).sqrt();
        assert!(mag.is_finite());
        max_disp = max_disp.max(mag);
        if c[1].abs() < 1e-9 {
            bottom_max = bottom_max.max(mag);
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

    // Physical sanity.
    assert!(bottom_max < 1e-12, "fixed base must not move");
    assert!(top_disp_y < 0.0, "free top must compress downward");
    assert!(top_disp_mag > 0.0, "load must produce deformation");
    let settle_est = (density * g * h / young) * h;
    assert!(
        top_disp_y.abs() < settle_est * 10.0 && top_disp_y.abs() > settle_est * 0.01,
        "settlement {top_disp_y:e} far from ρgh²/E estimate {settle_est:e}"
    );
    let _ = max_disp;
}
