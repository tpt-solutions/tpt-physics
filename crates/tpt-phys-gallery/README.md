# tpt-phys-gallery

End-to-end example-gallery runner for the `tpt-physics` workspace.

This is a **binary crate** (not a library): `cargo run -p tpt-phys-gallery`
prints a one-line "hello world" for every physics domain the workspace owns,
exercising the real solvers end to end. It is the smoke test that proves the
whole stack links and runs together.

> FEM itself lives in the sibling `tpt-fem` workspace; this crate only couples
> the physics crates (`tpt-phys-*`) that remain here.

## Run

```sh
cargo run -p tpt-phys-gallery
```

## Demos

Each `demo_*` function is a self-contained showcase of one crate:

| Demo | Crate | Shows |
| --- | --- | --- |
| `demo_core` | `tpt-phys-core` | Material registry lookup (Young's / shear modulus). |
| `demo_dem` | `tpt-phys-dem` | Two particles settling under gravity (DEM `World`). |
| `demo_cfd` | `tpt-phys-cfd` | Driven channel developing a Poiseuille-like profile (D2Q9 LBM). |
| `demo_fsi` | `tpt-phys-fsi` | Partitioned FSI: a wall pushed downstream by channel flow. |
| `demo_thermal_struct` | `tpt-phys-thermal-struct` | Self-equilibrated thermal-strain load vector (uniform temp ⇒ 0). |
| `demo_electro_thermal` | `tpt-phys-electro-thermal` | Resistive rod Joule-heating under voltage. |
| `demo_orchestrator` | `tpt-phys-orchestrator` | Multi-crate `Simulation` co-simulation stepping. |

For the expanded, standalone versions of each demo, see the `examples/`
directories of the individual crates and [`GALLERY.md`](./GALLERY.md).

## License

Dual-licensed under [MIT](../../LICENSE-MIT) and [Apache-2.0](../../LICENSE-APACHE-2.0).
Copyright TPT Solutions.
