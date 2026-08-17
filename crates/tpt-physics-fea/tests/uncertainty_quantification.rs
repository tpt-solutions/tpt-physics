//! Uncertainty quantification: Monte-Carlo sampling of material properties.
//!
//! Uses `proptest` to draw Young's modulus and Poisson's ratio from realistic
//! ranges and verifies the declarative FEA solve remains physical (finite,
//! bounded settlement) across the sampled population — a lightweight UQ smoke
//! test over the material parameters.

use proptest::prelude::*;
use tpt_physics_core::Material;
use tpt_physics_fea::spec::{DomainSpec, LoadSpec, ProblemSpec, SolverSpec};

fn strategy() -> impl Strategy<Value = (f64, f64)> {
    // E ∈ [1 GPa, 300 GPa], ν ∈ [0.1, 0.45] — a broad engineering range.
    (1e9..3e11f64, 0.1f64..0.45f64)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(40))]

    #[test]
    fn settlement_stays_physical_under_material_uncertainty((e, nu) in strategy()) {
        let spec = ProblemSpec {
            materials: None,
            material: tpt_physics_fea::spec::MaterialRef::Inline(Material::new(
                "UQ", e, nu, 1240.0, 68e-6,
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

        // Settlement must be finite and downward.
        prop_assert!(solved.free_top_settlement_y.is_finite());
        prop_assert!(solved.free_top_settlement_y < 0.0);

        // Bounded by the analytic ρgh²/E estimate to within two orders of
        // magnitude on either side (a wide band that should hold for any
        // realistic (E, ν) in range).
        let dens = 1240.0;
        let h = 0.05;
        let est = dens * 9.81 * h * h / e;
        let s = solved.free_top_settlement_y.abs();
        prop_assert!(
            s > est * 0.01 && s < est * 100.0,
            "settlement {s:.3e} outside band around ρgh²/E = {est:.3e}"
        );
    }
}
