//! Criterion benchmark harness for the LBM solver.
//!
//! Run `cargo bench -p tpt-phys-cfd` for long-term tracking. Throughput is
//! reported in lattice-node updates per second (divide by 1e6 for MLUPS).

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use tpt_phys_cfd::{Collision, Lbm2D};

fn bench_lbm(c: &mut Criterion) {
    let mut group = c.benchmark_group("lbm_step");
    for &(nx, ny) in &[(128usize, 128usize), (256usize, 128usize)] {
        group.throughput(Throughput::Elements((nx * ny) as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{nx}x{ny}")),
            &(nx, ny),
            |b, &(nx, ny)| {
                let mut lat = Lbm2D::new(nx, ny, 0.53);
                b.iter(|| lat.step([0.0, 0.0]));
            },
        );
    }
    group.finish();
}

fn bench_cavity(c: &mut Criterion) {
    let mut group = c.benchmark_group("cavity_mlups");
    for &n in &[128usize, 256, 512] {
        group.throughput(Throughput::Elements((n * n) as u64));
        for (name, par, col) in [
            ("serial_bgk", false, Collision::Bgk),
            ("parallel_bgk", true, Collision::Bgk),
            ("parallel_trt", true, Collision::TRT_WALL_EXACT),
        ] {
            group.bench_with_input(BenchmarkId::new(name, n), &n, |b, &n| {
                let mut lat = Lbm2D::cavity(n, 1000.0, 0.1);
                lat.set_collision(col);
                b.iter(|| {
                    if par {
                        lat.step_par([0.0, 0.0])
                    } else {
                        lat.step([0.0, 0.0])
                    }
                });
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_lbm, bench_cavity);
criterion_main!(benches);
