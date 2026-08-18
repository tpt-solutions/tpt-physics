# tpt-phys-thermal-struct

Thermal expansion and thermo-mechanical stress coupling (FEM to FEM).

Steady-state heat conduction is reused directly from `tpt-fem-thermal`, and
the structural solve from `tpt-fem-elasticity` — this crate's only job is the
coupling between them: converting a solved temperature field into a
thermal-strain load on the structural degrees of freedom.

Ported from `tpt-physics-fea` (removed) when this repo re-scoped to
multiphysics coupling per `spec2.txt`, with FEM itself delegated to `tpt-fem`.
