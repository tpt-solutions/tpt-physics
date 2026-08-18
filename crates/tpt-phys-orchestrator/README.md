# tpt-phys-orchestrator

AI-native coupling layer: adaptive time-stepping, solver coordination,
co-simulation API, and reinforcement-learning hooks.

Built on top of `tpt-sci-sim-core` (the sibling `tpt-science` repo's
multi-scale simulation engine — heterogeneous sub-model time-stepping,
cross-scale coupling, checkpointing) rather than reimplementing it. This
crate adds the physics-domain layer on top: RL/differentiable hooks (the
`rl` module, ported from the removed `tpt-physics-ai` crate) and, as
`tpt-phys-fsi`/`tpt-phys-thermal-struct`/`tpt-phys-electro-thermal` mature,
`SubModel` adapters wiring them into a `Simulation`.
