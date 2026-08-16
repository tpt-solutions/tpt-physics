# tpt-physics-ai

Differentiable physics wrappers exposing simulation state as Gymnasium-like
RL environments for `tpt-anima` reinforcement-learning agents.

Numerical differentiability is provided by `tpt-math-autodiff` (forward mode);
the Gym-style wrapper itself is net-new.

## Core traits

- `GymEnv` — minimal Gymnasium-like environment (`reset` → observation,
  `step(action)` → `(observation, reward, terminated)`).
- `DifferentiablePlant` — a plant whose transition `s' = f(s, a)` is
  differentiable; `jacobians` returns `∂f/∂s` and `∂f/∂a` via a single
  forward-mode AD pass.
- `GymWrapper` — adapts any `DifferentiablePlant` into a `GymEnv`, so an RL
  agent can both act on and differentiate through the physics.

## Plants

- `HarmonicOscillator` — differentiable spring–mass–damper.
- `Pendulum` — differentiable pendulum dynamics.

## License

Dual-licensed under [MIT](../../LICENSE-MIT) and [Apache-2.0](../../LICENSE-APACHE-2.0).
Copyright TPT Solutions.
