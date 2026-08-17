//! Criterion benchmark harness for the FEA stack (declarative problem spec).
//!
//! Run `cargo bench -p tpt-physics-fea` for long-term tracking.

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use tpt_physics_core::MaterialRegistry;
use tpt_physics_fea::spec::{DomainSpec, LoadSpec, ProblemSpec, SolverSpec};

fn spacer_spec(n: usize) -> ProblemSpec {
    ProblemSpec {
        materials: None,
        material: tpt_physics_fea::spec::MaterialRef::Inline(
            tpt_physics_core::Material::new("PLA", 3.5e9, 0.36, 1240.0, 68e-6),
        ),
        domain: DomainSpec::Box {
            min: [0.0, 0.0, 0.0],
            max: [0.04, 0.05, 0.04],
            n: [n, n + 1, n],
        },
        boundary_conditions: tpt_physics_fea::spec::BcSpec {
            fixed_planes: vec!["y_min".to_string()],
        },
        loads: LoadSpec {
            self_weight: true,
            gravity: 9.81,
        },
        solver: SolverSpec::StaticLinear,
    }
}

fn bench_spec(c: &mut Criterion) {
    let mut group = c.benchmark_group("fea_spec_solve");
    let reg = MaterialRegistry::new();
    for &n in &[4usize, 8usize] {
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, &n| {
            let spec = spacer_spec(n);
            b.iter(|| {
                let solved = spec.solve(&reg).unwrap();
                solved.free_top_settlement_y
            });
        });
    }
    group.finish();
}

criterion_group!(benches, bench_spec);
criterion_main!(benches);
