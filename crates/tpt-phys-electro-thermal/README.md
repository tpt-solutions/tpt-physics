# tpt-phys-electro-thermal

Electro-thermal coupling for `tpt-physics`: Joule heating, resistive losses, and
temperature-dependent conductivity.

No prior art for this exists in `tpt-physics`, `tpt-fem`, or `tpt-science`, so the
crate is built from scratch. It closes the multiphysics loop

```
electric field  →  current J = σ(T)·E  →  Joule heating q = σ(T)·|E|²
                                                        ↓
                              temperature T  →  σ(T)  (positive feedback)
                                                        ↓
                              heat conduction + surface convection
```

The reference solver is a 1-D resistive rod (`ElectroThermalRod`, explicit finite
differences) — the textbook Joule-heating geometry. The same `step` contract
generalises to an arbitrary mesh.

## Example

```rust
use tpt_phys_electro_thermal::ElectroThermalRod;

let mut rod = ElectroThermalRod::new(21, 300.0);
rod.dx = 0.01;
rod.set_voltage(10.0);
rod.convection = 50.0; // sink so a steady state exists
for _ in 0..2000 {
    rod.step(1e-4);
}
assert!(rod.temperatures().iter().all(|&t| t.is_finite()));
```

## Status

Scaffold. The 1-D operator is verified; a 3-D mesh-based variant (reusing
`tpt-fem-mesh` / `tpt-fem-thermal`) is future work.
