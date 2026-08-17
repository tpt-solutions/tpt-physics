//! Declarative problem specification: describe a static linear-elasticity
//! problem (material + geometry + boundary conditions + loads + solver) as a
//! single JSON document and solve it end-to-end.
//!
//! The spec reuses [`tpt_physics_core::MaterialRegistry`] (so a registry can be
//! embedded in the document via `MaterialRegistry::from_json`) for materials
//! and `tpt-fem-mesh-gen` / `tpt-fem-elasticity` for meshing and the actual
//! solve. It is intentionally a thin, serialization-first façade over the same
//! pipeline the hand-written examples use — no new numerics.
//!
//! ```json
//! {
//!   "materials": { "materials": [ { "name": "Steel", "youngs_modulus": 200e9,
//!                                  "poissons_ratio": 0.3, "density": 7850.0 } ] },
//!   "material": { "name": "Steel", "youngs_modulus": 200e9,
//!                 "poissons_ratio": 0.3, "density": 7850.0 },
//!   "domain":   { "type": "box", "min": [0,0,0], "max": [0.04,0.05,0.04],
//!                 "n": [4,5,4] },
//!   "boundary_conditions": { "fixed_planes": ["y_min"] },
//!   "loads":    { "self_weight": true, "gravity": 9.81 },
//!   "solver":   { "type": "static_linear" }
//! }
//! ```

use serde::{Deserialize, Serialize};
use tpt_fem_elasticity::{solve_elasticity, ElasticModel};
use tpt_fem_mesh_gen::box_mesh;
use tpt_physics_core::{Material, MaterialRegistry};

/// A fully-described static-linear-elasticity problem.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ProblemSpec {
    /// Optional embedded material registry (reuses `MaterialRegistry::from_json`).
    /// Materials referenced by name via [`ProblemSpec::material`] are resolved
    /// here first, then in the registry passed to [`ProblemSpec::solve`].
    #[serde(default)]
    pub materials: Option<MaterialRegistry>,
    /// The material to use: an inline definition or a name to resolve.
    pub material: MaterialRef,
    /// Geometry to mesh.
    pub domain: DomainSpec,
    /// Boundary conditions (essential / Dirichlet).
    #[serde(default)]
    pub boundary_conditions: BcSpec,
    /// Loads (body forces, etc.).
    #[serde(default)]
    pub loads: LoadSpec,
    /// Solver selection.
    #[serde(default)]
    pub solver: SolverSpec,
}

/// How the active material is supplied.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(untagged)]
pub enum MaterialRef {
    /// A complete inline [`Material`] definition.
    Inline(Material),
    /// A name resolved from an embedded / passed [`MaterialRegistry`].
    Ref(String),
}

/// Geometry to be meshed into a volume tet mesh.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum DomainSpec {
    /// An axis-aligned box meshed with `tpt-fem-mesh-gen::box_mesh`.
    Box {
        /// Minimum corner `[x, y, z]`.
        min: [f64; 3],
        /// Maximum corner `[x, y, z]`.
        max: [f64; 3],
        /// Element counts per axis `[nx, ny, nz]`.
        n: [usize; 3],
    },
}

/// Boundary conditions applied as essential (Dirichlet) constraints.
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct BcSpec {
    /// Planes whose nodes are fully fixed. Values are one of `x_min`, `x_max`,
    /// `y_min`, `y_max`, `z_min`, `z_max`.
    #[serde(default)]
    pub fixed_planes: Vec<String>,
}

/// Body loads.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct LoadSpec {
    /// Apply self-weight (gravity acting on the material density).
    #[serde(default)]
    pub self_weight: bool,
    /// Gravitational acceleration magnitude (m/s²). Default `9.81`.
    #[serde(default = "default_gravity")]
    pub gravity: f64,
}

fn default_gravity() -> f64 {
    9.81
}

impl Default for LoadSpec {
    fn default() -> Self {
        LoadSpec {
            self_weight: false,
            gravity: default_gravity(),
        }
    }
}

/// Solver selection. Only `static_linear` (continuum 3-D) is wired today; other
/// variants are accepted for forward-compatibility but resolve to the same
/// linear solve.
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "snake_case")]
pub enum SolverSpec {
    /// Linear static continuum elasticity (the only implemented path).
    #[default]
    StaticLinear,
}

/// Errors produced while resolving or solving a [`ProblemSpec`].
#[derive(Debug, Clone)]
pub enum SpecError {
    /// The named material could not be found in any supplied registry.
    MaterialNotFound(String),
    /// The underlying elasticity solve failed.
    Solve(String),
}

impl std::fmt::Display for SpecError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SpecError::MaterialNotFound(n) => write!(f, "material '{n}' not found"),
            SpecError::Solve(e) => write!(f, "elasticity solve failed: {e}"),
        }
    }
}

impl std::error::Error for SpecError {}

impl ProblemSpec {
    /// Serialize to a pretty JSON string.
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }

    /// Deserialize from a JSON string (also accepts embedded `materials` via
    /// [`MaterialRegistry::from_json`]'s schema).
    ///
    /// ```
    /// use tpt_physics_fea::spec::ProblemSpec;
    /// let json = r#"{
    ///   "material": { "name": "PLA", "youngs_modulus": 3.5e9, "poissons_ratio": 0.36, "density": 1240.0 },
    ///   "domain": { "type": "box", "min": [0,0,0], "max": [0.04,0.05,0.04], "n": [4,5,4] },
    ///   "boundary_conditions": { "fixed_planes": ["y_min"] },
    ///   "loads": { "self_weight": true }
    /// }"#;
    /// let spec = ProblemSpec::from_json(json).unwrap();
    /// let solved = spec.solve(&tpt_physics_core::MaterialRegistry::new()).unwrap();
    /// assert!(solved.free_top_settlement_y < 0.0); // compresses under self-weight
    /// ```
    pub fn from_json(s: &str) -> Result<Self, serde_json::Error> {
        serde_json::from_str(s)
    }

    /// Resolve the active [`Material`], searching the embedded registry (if
    /// any) then `external`.
    fn resolve_material(&self, external: &MaterialRegistry) -> Result<Material, SpecError> {
        match &self.material {
            MaterialRef::Inline(m) => Ok(m.clone()),
            MaterialRef::Ref(name) => external
                .get(name)
                .or_else(|| self.materials.as_ref().and_then(|r| r.get(name)))
                .cloned()
                .ok_or_else(|| SpecError::MaterialNotFound(name.clone())),
        }
    }

    /// Solve the described problem.
    ///
    /// `external` is an optional registry used to resolve materials given by
    /// name (the embedded `materials` map takes precedence). Returns the full
    /// displacement vector (`node_count × 3`) plus a small physical summary.
    pub fn solve(&self, external: &MaterialRegistry) -> Result<SolvedProblem, SpecError> {
        let material = self.resolve_material(external)?;

        let (min, max, n) = match &self.domain {
            DomainSpec::Box { min, max, n } => (*min, *max, *n),
        };
        let mesh = box_mesh(min, max, n);
        let n_nodes = mesh.node_count();

        // Bounding box → fixed-plane detection tolerances.
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for nd in 0..n_nodes {
            let c = mesh.node_coords(nd);
            for a in 0..3 {
                lo[a] = lo[a].min(c[a]);
                hi[a] = hi[a].max(c[a]);
            }
        }
        let tol = |a: usize| (hi[a] - lo[a]).abs() * 1e-9 + 1e-12;

        let plane_axis = |tag: &str| -> Option<(usize, bool)> {
            // Returns (axis, is_max).
            match tag {
                "x_min" => Some((0, false)),
                "x_max" => Some((0, true)),
                "y_min" => Some((1, false)),
                "y_max" => Some((1, true)),
                "z_min" => Some((2, false)),
                "z_max" => Some((2, true)),
                _ => None,
            }
        };

        let mut dirichlet = Vec::new();
        for nd in 0..n_nodes {
            let c = mesh.node_coords(nd);
            for tag in &self.boundary_conditions.fixed_planes {
                if let Some((axis, is_max)) = plane_axis(tag) {
                    let on_plane = if is_max {
                        (c[axis] - hi[axis]).abs() < tol(axis)
                    } else {
                        (c[axis] - lo[axis]).abs() < tol(axis)
                    };
                    if on_plane {
                        for comp in 0..3 {
                            dirichlet.push((nd * 3 + comp, 0.0));
                        }
                    }
                }
            }
        }

        let density = material.density;
        let g = self.loads.gravity;
        let body = move |_x: &[f64]| {
            if self.loads.self_weight {
                vec![0.0, -density * g, 0.0]
            } else {
                vec![0.0, 0.0, 0.0]
            }
        };

        let u = solve_elasticity(
            &mesh,
            ElasticModel::Continuum3D,
            material.youngs_modulus,
            material.poissons_ratio,
            2,
            body,
            &dirichlet,
        )
        .map_err(|e| SpecError::Solve(format!("{e:?}")))?;

        // Physical summary: fixed-base residual and free-top settlement.
        let mut bottom_max = 0.0_f64;
        let mut top_y = -1e9_f64;
        let mut top_node = 0usize;
        let mut max_abs = 0.0_f64;
        for nd in 0..n_nodes {
            let c = mesh.node_coords(nd);
            let mag = (u[nd * 3].powi(2) + u[nd * 3 + 1].powi(2) + u[nd * 3 + 2].powi(2)).sqrt();
            max_abs = max_abs.max(mag);
            if c[1] <= lo[1] + tol(1) {
                bottom_max = bottom_max.max(mag);
            }
            if c[1] > top_y {
                top_y = c[1];
                top_node = nd;
            }
        }
        let top_settlement_y = u[top_node * 3 + 1];

        Ok(SolvedProblem {
            displacements: u,
            n_nodes,
            n_elements: mesh.element_count(),
            max_abs_displacement: max_abs,
            fixed_base_max_displacement: bottom_max,
            free_top_settlement_y: top_settlement_y,
        })
    }
}

/// The result of solving a [`ProblemSpec`].
#[derive(Debug, Clone)]
pub struct SolvedProblem {
    /// Per-DOF displacement vector (`node_count × 3`).
    pub displacements: Vec<f64>,
    /// Number of mesh nodes.
    pub n_nodes: usize,
    /// Number of mesh elements.
    pub n_elements: usize,
    /// Maximum absolute nodal displacement magnitude (m).
    pub max_abs_displacement: f64,
    /// Maximum absolute displacement on the fixed base (should be ~0).
    pub fixed_base_max_displacement: f64,
    /// Vertical displacement of the highest node (negative ⇒ compression).
    pub free_top_settlement_y: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn json_roundtrip_preserves_spec() {
        let json = r#"{
            "material": { "name": "PLA", "youngs_modulus": 3.5e9, "poissons_ratio": 0.36, "density": 1240.0 },
            "domain": { "type": "box", "min": [0,0,0], "max": [0.04,0.05,0.04], "n": [4,5,4] },
            "boundary_conditions": { "fixed_planes": ["y_min"] },
            "loads": { "self_weight": true },
            "solver": "static_linear"
        }"#;
        let spec = ProblemSpec::from_json(json).expect("parse");
        let back = ProblemSpec::from_json(&spec.to_json().unwrap()).expect("re-parse");
        assert_eq!(spec.domain, back.domain);
        assert_eq!(spec.boundary_conditions, back.boundary_conditions);
    }

    #[test]
    fn solve_box_self_weight_is_physical() {
        let json = r#"{
            "material": { "name": "PLA", "youngs_modulus": 3.5e9, "poissons_ratio": 0.36, "density": 1240.0 },
            "domain": { "type": "box", "min": [0,0,0], "max": [0.04,0.05,0.04], "n": [4,5,4] },
            "boundary_conditions": { "fixed_planes": ["y_min"] },
            "loads": { "self_weight": true },
            "solver": "static_linear"
        }"#;
        let spec = ProblemSpec::from_json(json).expect("parse");
        let solved = spec.solve(&MaterialRegistry::new()).expect("solve");

        assert!(solved.n_nodes > 0 && solved.n_elements > 0);
        // Fixed base must not move.
        assert!(
            solved.fixed_base_max_displacement < 1e-9,
            "fixed base moved: {}",
            solved.fixed_base_max_displacement
        );
        // Free top must compress downward under self-weight.
        assert!(
            solved.free_top_settlement_y < 0.0,
            "top did not compress: {}",
            solved.free_top_settlement_y
        );
        // And be within an order of magnitude of ρgh²/E.
        let dens = 1240.0;
        let e = 3.5e9;
        let h = 0.05;
        let g = 9.81;
        let est = dens * g * h * h / e;
        assert!(
            solved.free_top_settlement_y.abs() < est * 10.0
                && solved.free_top_settlement_y.abs() > est * 0.01,
            "settlement {:.3e} far from estimate {:.3e}",
            solved.free_top_settlement_y,
            est
        );
    }

    #[test]
    fn named_material_resolves_from_embedded_registry() {
        let json = r#"{
            "materials": { "materials": [
                { "name": "Steel", "youngs_modulus": 200e9, "poissons_ratio": 0.3, "density": 7850.0 }
            ] },
            "material": "Steel",
            "domain": { "type": "box", "min": [0,0,0], "max": [1,1,1], "n": [2,2,2] },
            "boundary_conditions": { "fixed_planes": ["y_min", "y_max"] }
        }"#;
        let spec = ProblemSpec::from_json(json).expect("parse");
        let solved = spec.solve(&MaterialRegistry::new()).expect("solve");
        // Fully clamped box → essentially no movement.
        assert!(solved.max_abs_displacement < 1e-6, "got {}", solved.max_abs_displacement);
    }

    #[test]
    fn missing_named_material_errors() {
        let json = r#"{
            "material": "Ghost",
            "domain": { "type": "box", "min": [0,0,0], "max": [1,1,1], "n": [2,2,2] }
        }"#;
        let spec = ProblemSpec::from_json(json).expect("parse");
        assert!(matches!(
            spec.solve(&MaterialRegistry::new()),
            Err(SpecError::MaterialNotFound(_))
        ));
    }
}
