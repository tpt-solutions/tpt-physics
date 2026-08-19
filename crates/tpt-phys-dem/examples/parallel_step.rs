//! Parallel contact sweep: `World::step` vs `World::step_par`.
//!
//! [`World::step`] is a straightforward sequential semi-implicit Euler step.
//! [`World::step_par`] instead computes each particle's force independently
//! (summing its own neighbour contacts from the spatial-hash neighbour lists),
//! so the sweep distributes across the `rayon` thread pool with no
//! cross-particle write races. This is the CPU acceleration path for very large
//! particle counts (the crate's `large_scale` test drives >100k particles
//! through it).
//!
//! Because the pairwise Hertz–Mindlin contact force is exactly antisymmetric,
//! the two paths agree numerically on pairwise contacts — only the floating-point
//! summation order differs. The two steppers do differ in *boundary* handling;
//! see the note printed at the end. This example therefore uses a floor-free,
//! obstacle-free, bond-free cloud so the comparison is apples to apples.
//!
//! Run with (release is essential for a meaningful timing):
//!
//! ```text
//! cargo run --release --example parallel_step -p tpt-phys-dem
//! ```

use std::time::Instant;

use tpt_phys_dem::particle::Particle;
use tpt_phys_dem::world::World;

const NX: usize = 40;
const NY: usize = 12;
const NZ: usize = 40;
const R: f64 = 0.05;
/// Slight initial overlap so every particle carries real contact work.
const OVERLAP: f64 = 1.0e-3;
const STEPS: usize = 100;

fn build_world() -> World {
    let spacing = 2.0 * R - OVERLAP;
    let mut particles = Vec::with_capacity(NX * NY * NZ);
    for i in 0..NX {
        for j in 0..NY {
            for k in 0..NZ {
                particles.push(Particle::new(
                    [i as f64 * spacing, j as f64 * spacing, k as f64 * spacing],
                    [0.0; 3],
                    R,
                    1000.0,
                ));
            }
        }
    }
    let mut world = World::new(particles, 1e-5);
    world.e_star = 1e7; // soft contacts so the explicit step stays stable
    world.drag = 5.0; // damp the expanding cloud
    world.max_speed = 2.0;
    // No floor, no obstacles, no bonds ⇒ the two steppers are directly
    // comparable (see the note at the end of the output).
    world.floor_y = -1.0e3;
    world
}

fn main() {
    let mut sequential = build_world();
    let mut parallel = build_world();
    let n = sequential.particles.len();

    println!("DEM parallel contact sweep");
    println!("  particles         : {n}");
    println!("  steps             : {STEPS}");
    println!(
        "  threads available : {}",
        std::thread::available_parallelism()
            .map(|p| p.get())
            .unwrap_or(1)
    );

    let t0 = Instant::now();
    for _ in 0..STEPS {
        sequential.step();
    }
    let seq_elapsed = t0.elapsed();

    let t1 = Instant::now();
    for _ in 0..STEPS {
        parallel.step_par();
    }
    let par_elapsed = t1.elapsed();

    let total = (n * STEPS) as f64;
    let seq_rate = total / seq_elapsed.as_secs_f64();
    let par_rate = total / par_elapsed.as_secs_f64();

    println!();
    println!(
        "  {:<12} {:>12} {:>20}",
        "stepper", "wall time", "particle-steps/s"
    );
    println!("  {}", "-".repeat(46));
    println!("  {:<12} {:>12.3?} {:>20.3e}", "step()", seq_elapsed, seq_rate);
    println!(
        "  {:<12} {:>12.3?} {:>20.3e}",
        "step_par()", par_elapsed, par_rate
    );
    println!("  step_par speedup : {:.2}x", par_rate / seq_rate);

    // Numerical agreement: identical pairwise physics, only summation order
    // differs, so the energies must agree to round-off.
    let ke_seq = sequential.kinetic_energy();
    let ke_par = parallel.kinetic_energy();
    let rel = (ke_seq - ke_par).abs() / ke_seq.max(1e-30);
    println!();
    println!("  kinetic energy   : step() {ke_seq:.6e} J vs step_par() {ke_par:.6e} J");
    println!("  relative delta   : {rel:.2e} (summation-order round-off only)");

    assert!(
        ke_seq.is_finite() && ke_par.is_finite(),
        "simulation diverged"
    );
    assert!(
        rel < 1e-6,
        "parallel sweep disagreed with the sequential one by {rel}"
    );

    println!();
    println!("Note on boundaries: `step` additionally applies an *inelastic floor* —");
    println!("it zeroes any remaining downward velocity of a floor-contacting particle —");
    println!("which `step_par` does not, so a scene with an active floor settles");
    println!("slightly differently between the two steppers. Fixed-obstacle");
    println!("de-penetration *is* applied by both. Keep this in mind when switching an");
    println!("existing floor-bounded scene over to `step_par`.");
}
