# tpt-physics-core

Foundational, net-new data structures for the `tpt-physics` workspace.

Per the 2026-08-15 re-scope, `tpt-physics-core` contains **only code that is
genuinely new to the physics workspace**. SI units and mesh containers are
consumed *directly* from the sibling `tpt-math` / `tpt-fem` crates — there is no
wrapper or re-export shim.

## What lives here

| Module | Status | Description |
| --- | --- | --- |
| `material` | **[NEW]** | Type-safe material registry (`Material`, `MaterialRegistry`) with `serde` JSON serialization — Young's modulus, Poisson's ratio, density. No equivalent exists in `tpt-fem`/`tpt-math` (which pass material constants as bare function args). |
| `cad` | **[NEW]** | CAD→mesh ingestion adapter (`CadIngestor`) that streams `tpt-cad` / `biocad` solids into `tpt-fem-mesh`'s builder API. No equivalent exists in either sibling repo. |

The sibling units crate is re-exported as `tpt_physics_core::units` for
ergonomic single-import consumption.

## Example

```rust
use tpt_physics_core::MaterialRegistry;

let reg = MaterialRegistry::with_defaults();
let steel = reg.get("Structural Steel").unwrap();
assert!(steel.youngs_modulus > 100e9);
```

## License

Dual-licensed under [MIT](../../LICENSE-MIT) and [Apache-2.0](../../LICENSE-APACHE-2.0).
Copyright TPT Solutions.
