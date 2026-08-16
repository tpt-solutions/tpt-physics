# tpt-physics-solver

Iterative linear solvers, transient time-integration schemes, and CPU/GPU
hardware dispatch for `tpt-physics`.

The sibling `tpt-fem-sparse` crate provides only a **direct** LU solver. This
crate fills the net-new iterative-solver and transient-simulation gap.

## Modules

| Module | Description |
| --- | --- |
| `cg` | (Preconditioned) Conjugate Gradient — `cg`, `cg_pc`. The only iterative solver in the workspace; `tpt-fem-sparse` offers direct LU only. |
| `gmres` | Restarted GMRES for non-symmetric systems — `gmres`. |
| `time_integration` | `rk4` (Runge–Kutta) and `NewmarkBeta` (transient structural dynamics). `tpt-fem-eigen` only does frequency-domain modal analysis. |
| `dispatch` | `HardwareDispatch` / `ComputeTarget` route matrix kernels to CPU (`rayon`) or GPU (`wgpu` / `spark`) by problem size. `parallel_matvec_csr` is the CPU proxy. |
| `linalg` | `LinearOperator` trait (works on `tpt_fem_sparse::Csr`, dense matrices, or matrix-free operators), `DenseMat`, BLAS-like helpers. |
| `error` | `SolverError`, `SolveReport`. |

## License

Dual-licensed under [MIT](../../LICENSE-MIT) and [Apache-2.0](../../LICENSE-APACHE-2.0).
Copyright TPT Solutions.
