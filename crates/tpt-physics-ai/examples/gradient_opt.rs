//! Gradient-based design/control optimization on the differentiable path.
//!
//! Demonstrates using the forward-mode-autodiff Jacobians from
//! [`DifferentiablePlant::jacobians`] to do gradient descent: here we find a
//! *constant* control `u` that drives the [`Pendulum`] to rest (minimises the
//! terminal cost `θ² + ω²`). The same chain-rule machinery underpins the
//! topology / shape-optimization workflows the differentiable physics enables
//! (sensitivity analysis through the simulated plant).
//!
//! Run with: `cargo run --example gradient_opt`

use tpt_physics_ai::{DifferentiablePlant, Pendulum};

fn main() {
    let plant = Pendulum::new();
    let n = Pendulum::S; // state dim = 2
    let s0 = vec![0.8, 0.0]; // start displaced from upright
    let n_steps = 12;
    let lr = 0.3;

    let mut u = 0.0_f64;
    println!("Gradient-based control optimization (minimise terminal θ²+ω²):");
    for iter in 0..60 {
        // Forward pass + backward pass (chain rule) for constant control `u`.
        let mut s = s0.clone();
        let mut ds_du = vec![0.0; n]; // ∂s_N/∂u
        for _ in 0..n_steps {
            let (dfs, dfa) = plant.jacobians(&s, &[u]);
            let mut next = vec![0.0; n];
            for i in 0..n {
                let mut v = dfa[i][0]; // ∂f_i/∂u  (single action)
                for k in 0..n {
                    v += dfs[i][k] * ds_du[k];
                }
                next[i] = v;
            }
            ds_du = next;
            s = plant.step_prim(&s, &[u]);
        }
        let cost = s[0] * s[0] + s[1] * s[1];
        let dc_du = 2.0 * s[0] * ds_du[0] + 2.0 * s[1] * ds_du[1];
        if iter % 10 == 0 || iter == 59 {
            println!("  iter {iter:>3}: u = {u:+.4}  terminal_cost = {cost:.5}");
        }
        u -= lr * dc_du;
        u = u.clamp(-5.0, 5.0);
    }
}
