//! Benchmark: DEM parallel stepper throughput.
//!
//! Times [`World::step_par`] (the `rayon` CPU-acceleration path) on 10k and
//! 100k-particle beds and reports particles processed per second, demonstrating
//! the scale at which the hardware-dispatch API would route to a GPU.

use tpt_physics_dem::particle::Particle;
use tpt_physics_dem::world::World;

fn make(n: usize) -> World {
    let r = 0.25;
    let l = (n as f64).sqrt().ceil() * 1.0 + 4.0;
    let mut rng = Lcg::new(0xBEEF);
    let mut particles = Vec::with_capacity(n);
    for _ in 0..n {
        let x = r + rng.next_f64() * (l - 2.0 * r);
        let y = r + rng.next_f64() * (l - 2.0 * r);
        let z = r + rng.next_f64() * (l - 2.0 * r);
        particles.push(Particle::new([x, y, z], [0.0; 3], r, 1000.0));
    }
    let mut w = World::new(particles, 1e-4);
    w.e_star = 5e7;
    w.max_speed = 15.0;
    w
}

fn main() {
    for &n in &[10_000usize, 100_000usize] {
        let mut world = make(n);
        let steps = 20;
        let t0 = std::time::Instant::now();
        for _ in 0..steps {
            world.step_par();
        }
        let elapsed = t0.elapsed();
        let per_step = elapsed / steps as u32;
        let rate = n as f64 / per_step.as_secs_f64();
        println!(
            "{:>7} particles: {:>8.2?}/step  =>  {:.0} particles/s",
            n, per_step, rate
        );
    }
}

struct Lcg(u64);
impl Lcg {
    fn new(seed: u64) -> Self {
        Lcg(seed | 1)
    }
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0
    }
    fn next_f64(&mut self) -> f64 {
        ((self.next_u64() >> 11) as f64) / (1u64 << 53) as f64
    }
}
