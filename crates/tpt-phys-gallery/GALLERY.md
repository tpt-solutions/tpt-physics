# Gallery

A tour of every physics domain in the `tpt-physics` workspace, with the
standalone example to read for each. Run the whole tour with:

```sh
cargo run -p tpt-phys-gallery
```

| Domain | Crate | Standalone examples | Gallery demo |
| --- | --- | --- | --- |
| Materials & CAD | `tpt-phys-core` | `material_registry`, `cad_ingest`, `uq_natural_frequency` (feature `uq`) | `demo_core` |
| Discrete elements | `tpt-phys-dem` | `hopper_discharge`, `cohesive_bonds`, `obstacles_ssi`, `heat_conduction`, `checkpoint`, `parallel_step` | `demo_dem` |
| Fluid (LBM + SPH) | `tpt-phys-cfd` | `cavity`, `poiseuille`, `flow_past_cylinder`, `dam_break_sph` | `demo_cfd` |
| Fluid–structure | `tpt-phys-fsi` | `compliant_wall`, `mesh_mapping` | `demo_fsi` |
| Thermal–structural | `tpt-phys-thermal-struct` | `uniform_expansion`, `layered_strip` | `demo_thermal_struct` |
| Electro-thermal | `tpt-phys-electro-thermal` | `heated_rod`, `self_limiting` | `demo_electro_thermal` |
| Orchestration | `tpt-phys-orchestrator` | `rl_pendulum`, `uq_coupled`, `coupled_simulation` | `demo_orchestrator` |

## Running a standalone example

```sh
cargo run --release --example poiseuille -p tpt-phys-cfd
```

Each example prints diagnostic output and asserts the expected physical
behaviour (e.g. a parabolic profile, vortex shedding, self-equilibrated thermal
load, finite coupled stepping).
