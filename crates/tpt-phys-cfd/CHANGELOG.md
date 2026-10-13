# Changelog

All notable changes to this crate are documented here. This project
adheres to [Semantic Versioning](https://semver.org).

## [0.1.0] - Unreleased

### Added
- `Lbm2D::cavity(n, re, u_lid)` benchmark constructor, `ghia` module (Ghia, Ghia
  & Shin 1982 Re 100/400 reference profiles + centreline error metrics),
  `Collision::{Bgk, Trt, Mrt}` operators (`Lbm2D::set_collision`), and a
  `LatticeModel` trait implemented for `D2Q9`.
- Cavity validation: Re = 100 centrelines within ~1% of Ghia et al. (default
  test) plus `#[ignore]`d refined Re = 100 / 400 runs.
- `D3Q19` lattice data behind `LatticeModel` (isotropy-tested; no 3-D solver yet).
- `cavity_mlups` criterion benchmark (serial / parallel, BGK / TRT).

### Fixed
- `cylinder_vortex_shedding` diverged to NaN with BGK at τ = 0.55; it now runs
  at τ = 0.6 with the cylinder scaled to Re ≈ 108 (supercritical at this
  blockage). Note: `Collision::Trt` is *less* stable than BGK for open
  inlet/outlet flows at low τ (under-relaxed odd modes), so keep BGK there.

### Changed
- Solver internals rewritten as Structure-of-Arrays with a fused,
  `wide::f64x4`-vectorised collide/stream kernel and row-parallel `step_par`
  (bit-identical to `step`, and to the previous scalar scheme to 1e-12 —
  checked by a proptest against a retained scalar reference). Measured on the
  driven-cavity bench: `step` 1.1–1.4× faster, `step_par` 1.8–3.2× faster
  (e.g. 128²: 10.9 → 14.3 / 12.8 → 22.9 MLUPS; 512²: 7.0 → 7.7 / 10.3 → 32.1).
- The public API (`Lbm2D`, `XBoundary`, public fields) is unchanged; `Lbm2D`
  is now `Clone`.
- The cavity no longer fails: the earlier "primary-vortex bug" was in the test
  (off-centre sampling, loose thresholds), not the solver — see the README.

### Added (initial)
- `Lbm2D` D2Q9 Lattice Boltzmann solver: BGK collision, half-way bounce-back
  boundaries (stationary and moving lids), periodic / inlet–outlet / open `x`
  boundaries, Guo body-force term, `add_circle` / `add_rect` obstacles, and a
  `rayon` `step_par` variant.
- `Sph2D` weakly-compressible SPH solver (Müller et al. 2003): Poly6 / spiky /
  viscous kernels, Tait EOS, uniform-grid neighbour search, density auto-
  calibration.
- Examples: `cavity`, `poiseuille`, `flow_past_cylinder`, `dam_break_sph`.
- Validation tests: lid-driven cavity, flow-past-cylinder (steady + shedding),
  SPH dam-break finiteness and settling.

[0.1.0]: https://github.com/tpt-physics/tpt-physics/releases/tag/tpt-phys-cfd-0.1.0
