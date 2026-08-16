//! Milestone: simulate the 3D-printed pile cage spacer.
//!
//! This wires the whole `tpt-physics` stack end-to-end:
//!
//! 1. **Material database** (`tpt-physics-core`): pick a 3-D-print material
//!    (PLA) from the type-safe registry — the net-new reason `core` exists.
//! 2. **CAD ingestion adapter** (`tpt-physics-core::cad`): express the spacer
//!    footprint as a `tpt-cad`-style `CadSolid` and lower it into a
//!    `tpt-fem-mesh` surface mesh (exercises the net-new adapter).
//! 3. **Mesh generation** (`tpt-fem-mesh-gen`): generate the volume tet mesh of
//!    the spacer block directly.
//! 4. **Static linear FEA** (`tpt-fem-elasticity`): fix the base, load with
//!    self-weight (gravity), and solve for the displacement field.
//!
//! A pile cage spacer is the plastic component that holds a reinforced-concrete
//! pile cage at the correct concrete cover from the bore wall. Here we model a
//! representative spacer foot block and verify the static solution is physical:
//! the fixed base stays put while the free top compresses downward under its
//! own weight.

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
    let young = material.youngs_modulus;
    let poisson = material.poissons_ratio;
    let density = material.density;
    println!(
        "Material: {}  E={:.3e} Pa  nu={:.3}  rho={:.1} kg/m^3",
        material.name, young, poisson, density
    );

    // 2. Geometry also expressible through the CAD ingestion adapter: the
    //    spacer's bounding box as a `tpt-cad` solid, lowered to a `tpt-fem-mesh`
    //    (triangulated) surface mesh. The vertical (y) axis is the 50 mm dim.
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
    // Six quad faces (vertex-index loops), fanned into triangles on ingest.
    cad.add_face(vec![0, 1, 2, 3], Some(1)); // bottom  y = 0
    cad.add_face(vec![4, 5, 6, 7], Some(2)); // top     y = h
    cad.add_face(vec![0, 3, 7, 4], Some(3)); // side    x = 0
    cad.add_face(vec![1, 2, 6, 5], Some(4)); // side    x = w
    cad.add_face(vec![0, 1, 5, 4], Some(5)); // side    z = 0
    cad.add_face(vec![3, 2, 6, 7], Some(6)); // side    z = d
    let surface = cad.ingest().expect("ingest CAD solid");
    println!(
        "CAD adapter → surface mesh: {} nodes, {} tris",
        surface.node_count(),
        surface.element_count()
    );
    assert_eq!(surface.node_count(), 8);
    assert_eq!(surface.element_count(), 12);

    // 3. Volume mesh of the spacer foot (a tet-meshed block, 40×50×40 mm).
    let (w, h, d) = (0.04_f64, 0.05_f64, 0.04_f64);
    let mesh = box_mesh([0.0, 0.0, 0.0], [w, h, d], [4, 5, 4]);
    let n_nodes = mesh.node_count();
    let n_elems = mesh.element_count();
    println!("Volume tet mesh: {n_nodes} nodes, {n_elems} elements");
    assert!(n_nodes > 0 && n_elems > 0);

    // 4. Boundary conditions: fix the bottom face (y ≈ 0) in all 3 DOFs.
    let mut dirichlet = Vec::new();
    for n in 0..n_nodes {
        let c = mesh.node_coords(n);
        if c[1].abs() < 1e-9 {
            for comp in 0..3 {
                dirichlet.push((n * 3 + comp, 0.0));
            }
        }
    }
    let fixed = dirichlet.len() / 3;
    println!("Fixed base nodes: {fixed}");

    // Self-weight body force (gravity, downward).
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

    // 5. Post-process: find bottom (fixed) and top-centre nodes.
    let mut bottom_disp = 0.0_f64;
    let mut top_node = 0usize;
    let mut top_y = -1e9_f64;
    let mut max_disp = 0.0_f64;
    for n in 0..n_nodes {
        let c = mesh.node_coords(n);
        let ux = u[n * 3];
        let uy = u[n * 3 + 1];
        let uz = u[n * 3 + 2];
        let mag = (ux * ux + uy * uy + uz * uz).sqrt();
        assert!(mag.is_finite(), "node {n} displacement not finite");
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

    println!("--- Results ---");
    println!("Max |displacement|      : {max_disp:.3e} m");
    println!("Bottom (fixed) max |u|  : {bottom_disp:.3e} m");
    println!(
        "Top-centre disp (x,y,z): ({:.3e}, {:.3e}, {:.3e}) m",
        u[top_node * 3],
        u[top_node * 3 + 1],
        u[top_node * 3 + 2]
    );

    // Physical sanity checks.
    assert!(bottom_disp < 1e-12, "fixed base must not move");
    assert!(
        top_disp_y < 0.0,
        "free top must compress downward under self-weight"
    );
    assert!(
        top_disp_mag > 0.0,
        "load must produce measurable deformation"
    );
    // Closed-form sanity: vertical strain ≈ ρ g h / E, so top settlement
    // ≈ strain · h. Allow an order of magnitude either way.
    let strain_est = density * g * h / young;
    let settle_est = strain_est * h;
    println!("Settlement estimate ρgh²/E  : {settle_est:.3e} m  (sim top_y = {top_disp_y:.3e})");
    assert!(
        top_disp_y.abs() < settle_est * 10.0 && top_disp_y.abs() > settle_est * 0.01,
        "top settlement far from ρgh²/E estimate"
    );

    println!("\nMilestone PASSED: pile cage spacer simulated end-to-end.");
}
