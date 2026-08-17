//! Criterion benchmark harness for the iterative solvers.
//!
//! Builds a 2-D Poisson (5-point Laplacian) system of adjustable size and
//! benchmarks the (preconditioned) CG and GMRES solves. Run
//! `cargo bench -p tpt-physics-solver` for long-term tracking.

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use tpt_physics_solver::cg::cg;
use tpt_physics_solver::cg::cg_pc;
use tpt_physics_solver::gmres::gmres;
use tpt_physics_solver::linalg::csr_from_dense;

/// Assemble a dense `side×side` 2-D Poisson matrix (Dirichlet) and return it
/// plus a unit RHS.
fn poisson(side: usize) -> (tpt_fem_sparse::Csr, Vec<f64>) {
    let n = side * side;
    let mut dense = vec![0.0; n * n];
    let idx = |i: usize, j: usize| i * side + j;
    for i in 0..side {
        for j in 0..side {
            let k = idx(i, j);
            dense[k * n + k] = 4.0;
            if i > 0 {
                dense[k * n + idx(i - 1, j)] = -1.0;
            }
            if i + 1 < side {
                dense[k * n + idx(i + 1, j)] = -1.0;
            }
            if j > 0 {
                dense[k * n + idx(i, j - 1)] = -1.0;
            }
            if j + 1 < side {
                dense[k * n + idx(i, j + 1)] = -1.0;
            }
        }
    }
    let b = vec![1.0; n];
    (csr_from_dense(n, n, &dense), b)
}

fn bench_solvers(c: &mut Criterion) {
    let mut group = c.benchmark_group("solvers");
    for &side in &[40usize, 80usize] {
        let (a, b) = poisson(side);
        let n = side * side;
        group.bench_with_input(BenchmarkId::new("cg", n), &n, |bc, _| {
            bc.iter(|| {
                let (x, _) = cg(&a, &b, None, 1e-8, 2000).unwrap();
                x
            });
        });
        group.bench_with_input(BenchmarkId::new("cg_jacobi", n), &n, |bc, _| {
            let dinv: Vec<f64> = (0..a.nrows)
                .map(|i| 1.0 / a.values[a.row_ptrs[i]])
                .collect();
            bc.iter(|| {
                let (x, _) = cg_pc(&a, &b, None, 1e-8, 2000, Some(&|r: &[f64], z: &mut [f64]| {
                    for i in 0..r.len() {
                        z[i] = dinv[i] * r[i];
                    }
                }))
                .unwrap();
                x
            });
        });
        group.bench_with_input(BenchmarkId::new("gmres", n), &n, |bc, _| {
            bc.iter(|| {
                let (x, _) = gmres(&a, &b, None, 30, 1e-8, 2000).unwrap();
                x
            });
        });
    }
    group.finish();
}

criterion_group!(benches, bench_solvers);
criterion_main!(benches);
