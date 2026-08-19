# tpt-physics-wasm

**WebAssembly bindings for the `tpt-physics` solvers** — a browser playground
that runs the pure-Rust DEM and CFD engines directly in WebGL, with no server.

This crate is the repo's WebAssembly front door. It exposes two solvers to
JavaScript via [`wasm-bindgen`](https://crates.io/crates/wasm-bindgen):

* **DEM** ([`tpt_phys_dem`]) — a granular `World` you drive step-by-step.
* **CFD** ([`tpt_phys_cfd`]) — the D2Q9 lattice-Boltzmann (`Lbm2D`) solver.

Load a scene (particles/obstacles or a lattice setup) as JSON, advance it
step-by-step, and pull state back out as flat `Float32Array`s ready to upload
to a WebGL buffer.

> **Note (2026-08-19):** the playground currently exposes **DEM**, **CFD**, and
> **electro-thermal**. FSI / orchestrator bindings are planned but not yet wired
> into the WebGL frontend — see `src/lib.rs` for the exact surface.

## Build

Requires `wasm-pack` (or `cargo` + `wasm-bindgen-cli`) on `PATH`.

```sh
# from the workspace root
just wasm          # runs scripts/build_wasm.ps1 (or .sh)
just serve-wasm    # serves crates/tpt-physics-wasm/www on :8080
```

`just wasm` produces `crates/tpt-physics-wasm/www/pkg`; open
`crates/tpt-physics-wasm/www/index.html` (served by `just serve-wasm`) to use
the playground.

## Public surface

* `DemSimulation` — `new(json)`, `step()`, `count()`, `positions()`,
  `velocities()`, `temperatures()`, `kinetic_energy()`.
* `CfdSimulation` — `new(json)`, `step()`, `nx()`, `ny()`, `velocity()`,
  `solid()`.

See the doc comments on each constructor for the exact JSON scene schema
(particles + obstacles for DEM; lattice size, walls, obstacles, moving lid for
CFD).
