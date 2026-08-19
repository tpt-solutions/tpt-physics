//! `SubModel` adapters wiring the `tpt-phys` domain coupling crates into a
//! [`tpt_sci_sim_core::Simulation`].
//!
//! Each adapter implements [`tpt_sci_sim_core::SubModel`] over a concrete
//! solver so the generic co-simulation engine can drive it:
//!
//! * [`FsiSubModel`] — partitioned fluid–structure interaction
//!   (`tpt-phys-fsi` + `tpt-phys-cfd`),
//! * [`ElectroThermalSubModel`] — Joule heating (`tpt-phys-electro-thermal`),
//! * [`ThermalStructSubModel`] — thermal-to-structural coupling
//!   (`tpt-phys-thermal-struct` + `tpt-phys-core`).
//!
//! [`build_demo_simulation`] registers one of each and couples the
//! electro-thermal temperature field into the thermal-structural model,
//! demonstrating the multi-crate orchestration end to end.

use std::fmt::Debug;

use tpt_fem_mesh::{CellType, Mesh, MeshBuilder};
use tpt_phys_cfd::Lbm2D;
use tpt_phys_core::Material;
use tpt_phys_electro_thermal::ElectroThermalRod;
use tpt_phys_fsi::{FsiDriver, LumpedStructure, StructuralModel};
use tpt_phys_thermal_struct::thermal_load_vector;
use tpt_sci_sim_core::{CouplingFn, Simulation, SubModel};

/// Partitioned FSI sub-model for the co-simulation engine.
pub struct FsiSubModel {
    name: String,
    fluid: Lbm2D,
    structure: LumpedStructure,
    driver: FsiDriver,
}

impl std::fmt::Debug for FsiSubModel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FsiSubModel")
            .field("name", &self.name)
            .field("structure", &self.structure)
            .field("driver", &self.driver)
            .finish_non_exhaustive()
    }
}

impl FsiSubModel {
    /// Build from a fluid lattice, a structural model, and its interface mesh.
    pub fn new(fluid: Lbm2D, structure: LumpedStructure, structural_mesh: &Mesh) -> Self {
        let driver = FsiDriver::new(&fluid, &structure, structural_mesh);
        FsiSubModel {
            name: "fsi".to_string(),
            fluid,
            structure,
            driver,
        }
    }

    /// Immutable access to the fluid lattice (for inspection).
    pub fn fluid(&self) -> &Lbm2D {
        &self.fluid
    }
}

impl SubModel for FsiSubModel {
    fn name(&self) -> &str {
        &self.name
    }
    fn step(&mut self, dt: f64) {
        self.driver.step(&mut self.fluid, &mut self.structure, dt);
    }
    fn state_dim(&self) -> usize {
        self.structure.n_interface_nodes() * 3
    }
    fn gather_state(&self, out: &mut [f64]) {
        for k in 0..self.structure.n_interface_nodes() {
            let d = self.structure.displacement(k);
            out[3 * k] = d[0];
            out[3 * k + 1] = d[1];
            out[3 * k + 2] = d[2];
        }
    }
    fn apply_input(&mut self, _input: &[f64]) {
        // The partitioned loop already feeds structure motion back to the fluid
        // internally; external interface data is not consumed here.
    }
}

/// Electro-thermal sub-model (Joule heating) for the co-simulation engine.
#[derive(Debug)]
pub struct ElectroThermalSubModel {
    name: String,
    rod: ElectroThermalRod,
}

impl ElectroThermalSubModel {
    /// Wrap an [`ElectroThermalRod`].
    pub fn new(rod: ElectroThermalRod) -> Self {
        ElectroThermalSubModel {
            name: "electro-thermal".to_string(),
            rod,
        }
    }

    /// Immutable access to the rod (for inspection).
    pub fn rod(&self) -> &ElectroThermalRod {
        &self.rod
    }
}

impl SubModel for ElectroThermalSubModel {
    fn name(&self) -> &str {
        &self.name
    }
    fn step(&mut self, dt: f64) {
        self.rod.step(dt);
    }
    fn state_dim(&self) -> usize {
        self.rod.len()
    }
    fn gather_state(&self, out: &mut [f64]) {
        out.copy_from_slice(self.rod.temperatures());
    }
    fn apply_input(&mut self, _input: &[f64]) {
        // Heat source is internal (applied voltage); external coupling feeds in
        // via the rod API directly if needed.
    }
}

/// Thermal-to-structural sub-model: a lumped structural response driven by the
/// thermal-strain load [`thermal_load_vector`] computes for a tetrahedral mesh.
#[derive(Debug)]
pub struct ThermalStructSubModel {
    name: String,
    mesh: Mesh,
    material: Material,
    t_ref: f64,
    mass: f64,
    stiffness: f64,
    damping: f64,
    pos: Vec<[f64; 3]>,
    vel: Vec<[f64; 3]>,
    force: Vec<[f64; 3]>,
    /// Per-node temperature (K), updated by the electro-thermal coupling.
    temp: Vec<f64>,
}

impl ThermalStructSubModel {
    /// Build from a structural tet mesh, material, and lumped oscillator params.
    pub fn new(
        mesh: Mesh,
        material: Material,
        t_ref: f64,
        mass: f64,
        stiffness: f64,
        damping: f64,
    ) -> Self {
        let n = mesh.node_count();
        ThermalStructSubModel {
            name: "thermal-struct".to_string(),
            mesh,
            material,
            t_ref,
            mass,
            stiffness,
            damping,
            pos: vec![[0.0; 3]; n],
            vel: vec![[0.0; 3]; n],
            force: vec![[0.0; 3]; n],
            temp: vec![t_ref; n],
        }
    }

    /// Immutable temperatures (K).
    pub fn temperatures(&self) -> &[f64] {
        &self.temp
    }
}

impl SubModel for ThermalStructSubModel {
    fn name(&self) -> &str {
        &self.name
    }
    fn step(&mut self, dt: f64) {
        // Thermal strain from the current temperature field becomes a nodal
        // load through the ported coupling primitive.
        let load = thermal_load_vector(&self.mesh, 3, &self.material, &self.temp, self.t_ref);
        for (k, f) in self.force.iter_mut().enumerate() {
            f[0] += load[k * 3];
            f[1] += load[k * 3 + 1];
            f[2] += load[k * 3 + 2];
        }
        for k in 0..self.pos.len() {
            for i in 0..3 {
                let a = (self.force[k][i]
                    - self.stiffness * self.pos[k][i]
                    - self.damping * self.vel[k][i])
                    / self.mass;
                self.vel[k][i] += a * dt;
                self.pos[k][i] += self.vel[k][i] * dt;
            }
            self.force[k] = [0.0; 3];
        }
    }
    fn state_dim(&self) -> usize {
        self.pos.len() * 3
    }
    fn gather_state(&self, out: &mut [f64]) {
        for k in 0..self.pos.len() {
            out[3 * k] = self.pos[k][0];
            out[3 * k + 1] = self.pos[k][1];
            out[3 * k + 2] = self.pos[k][2];
        }
    }
    fn apply_input(&mut self, input: &[f64]) {
        // Electro-thermal coupling delivers a temperature field of equal length.
        if input.len() == self.temp.len() {
            for (t, v) in self.temp.iter_mut().zip(input) {
                *t = *v;
            }
        }
    }
}

/// A single tetrahedral element spanning the unit cube, used as the structural
/// domain for the demo [`ThermalStructSubModel`].
fn unit_tet_mesh() -> Mesh {
    let mut b = MeshBuilder::new();
    let n0 = b.add_node(vec![0.0, 0.0, 0.0]);
    let n1 = b.add_node(vec![1.0, 0.0, 0.0]);
    let n2 = b.add_node(vec![0.0, 1.0, 0.0]);
    let n3 = b.add_node(vec![0.0, 0.0, 1.0]);
    b.add_element(CellType::Tet, vec![n0, n1, n2, n3]);
    b.build()
}

/// Build a coupled demo simulation wiring electro-thermal heating into the
/// thermal-structural model (temperature → thermal-strain load).
///
/// `material` is the structural material used by the thermal-structural model and
/// `voltage` is the applied electro-thermal rod voltage — both are parameters so
/// callers (e.g. a Monte-Carlo UQ sweep) can vary them.
pub fn build_demo_simulation_for(material: &Material, voltage: f64) -> Simulation {
    // FSI: a driven channel whose right wall is the fluid–structure interface.
    let mut fluid = Lbm2D::new(32, 16, 0.6);
    fluid.set_horizontal_walls();
    fluid.add_rect(30, 1, 30, 14);
    fluid.initialise(1.0, [0.1, 0.0]);
    let mut fsi_b = MeshBuilder::new();
    fsi_b.add_node(vec![0.0, 0.0, 0.0]);
    let fsi_mesh = fsi_b.build();
    let fsi = FsiSubModel::new(fluid, LumpedStructure::new(1, 1.0, 10.0, 2.0), &fsi_mesh);

    // Electro-thermal: a rod under voltage (steady Joule heating).
    let mut rod = ElectroThermalRod::new(11, 300.0);
    rod.dx = 0.01;
    rod.set_voltage(voltage);
    rod.convection = 50.0;
    let et = ElectroThermalSubModel::new(rod);

    // Thermal-structural: the unit tet driven by thermal strain.
    let ts = ThermalStructSubModel::new(unit_tet_mesh(), material.clone(), 300.0, 1.0, 10.0, 2.0);

    let mut sim = Simulation::new();
    let i_et = sim.add_submodel(Box::new(et));
    let i_ts = sim.add_submodel(Box::new(ts));
    let _ = sim.add_submodel(Box::new(fsi));
    // Electro-thermal temperature field → thermal-structural temperature.
    let coupling: CouplingFn = Box::new(|src: &[f64], dst: &mut [f64]| {
        let n = src.len().min(dst.len());
        dst[..n].copy_from_slice(&src[..n]);
    });
    sim.add_coupling(i_et, i_ts, coupling);
    sim
}

/// Build the default coupled demo simulation (nominal steel, 10 V rod).
pub fn build_demo_simulation() -> Simulation {
    let mat = Material::new("Demo", 200e9, 0.3, 7850.0, 12e-6);
    build_demo_simulation_for(&mat, 10.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn demo_simulation_steps_and_stays_finite() {
        let mut sim = build_demo_simulation();
        for _ in 0..200 {
            sim.step(1e-4);
        }
        let mut buf = vec![0.0; sim.submodel(0).state_dim()];
        sim.submodel(0).gather_state(&mut buf);
        assert!(buf.iter().all(|v| v.is_finite()), "electro-thermal NaN");
        let mut buf2 = vec![0.0; sim.submodel(1).state_dim()];
        sim.submodel(1).gather_state(&mut buf2);
        assert!(buf2.iter().all(|v| v.is_finite()), "thermal-struct NaN");
    }
}
