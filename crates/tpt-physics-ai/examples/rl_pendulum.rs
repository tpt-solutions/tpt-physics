//! Hello-world for `tpt-physics-ai`: a differentiable pendulum as an RL env.
//!
//! Wraps the net-new [`Pendulum`] differentiable plant in a Gymnasium-like
//! [`GymWrapper`], steps it with a (zero) action, and prints the
//! forward-mode-autodiff Jacobians `∂f/∂s` and `∂f/∂a`. Run with:
//! `cargo run --example rl_pendulum`

use tpt_math_linalg_dense::DVector;
use tpt_physics_ai::{GymEnv, GymWrapper, Pendulum};

fn main() {
    let plant = Pendulum::new(); // g = L = m = 1, light damping, dt = 0.01
    let mut env = GymWrapper::new(plant, vec![0.5, 0.0], 100); // start at θ = 0.5 rad

    // Step the environment a few times with a zero control action.
    let action = DVector::from_fn(1, |_| 0.0);
    for step in 0..5 {
        let obs = env.state().to_vec();
        let out = env.step(&action);
        println!(
            "step {step}: θ = {:.4}, ω = {:.4}  (reward = {:.4}, terminated = {})",
            obs[0], obs[1], out.reward, out.terminated
        );
    }

    // The differentiable plant: Jacobians at the current state (zero action).
    let (dfs, dfa) = env.gradients();
    println!("\n∂f/∂s (state transition w.r.t. state):");
    for row in &dfs {
        println!("  [{:.4} {:.4}]", row[0], row[1]);
    }
    println!("∂f/∂a (state transition w.r.t. action):");
    for row in &dfa {
        println!("  [{:.4}]", row[0]);
    }
}
