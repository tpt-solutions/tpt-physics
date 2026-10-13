//! Collision operators for the D2Q9 lattice.
//!
//! The collision step relaxes the distributions toward their local
//! equilibrium. Every operator here recovers the same Navier–Stokes
//! viscosity `ν = cs² (τ − ½)`; they differ in how the *non-hydrodynamic*
//! ("ghost") modes are damped, which is what decides stability at low
//! viscosity and how accurate bounce-back walls are.

/// Collision operator selection for [`crate::Lbm2D`].
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum Collision {
    /// Single-relaxation-time Bhatnagar–Gross–Krook: every moment relaxes at
    /// `1/τ`. The default; simplest and cheapest.
    #[default]
    Bgk,
    /// Two-relaxation-time (Ginzburg). The symmetric part relaxes at `1/τ`
    /// (fixing the viscosity) and the antisymmetric part at `1/τ⁻`, where
    /// `τ⁻` is chosen so that the "magic" parameter
    /// `Λ = (τ − ½)(τ⁻ − ½)` takes the given value. `Λ = 3/16` places
    /// half-way bounce-back walls exactly for Poiseuille flow, `Λ = 1/4`
    /// maximises stability. `Λ = (τ − ½)²` reproduces [`Collision::Bgk`].
    Trt {
        /// The magic parameter `Λ > 0`.
        magic: f64,
    },
    /// Multiple-relaxation-time (Lallemand & Luo 2000) in the D2Q9 moment
    /// basis `(ρ, e, ε, jx, qx, jy, qy, pxx, pxy)`. The shear moments
    /// (`pxx`, `pxy`) relax at `1/τ` (fixing the viscosity) and the momentum
    /// moments at `1/τ` too (so body forcing behaves exactly as in BGK); the
    /// energy, energy-square and heat-flux ghost modes relax at the given
    /// free rates, each in `(0, 2)`. With every rate `1/τ` this is BGK.
    Mrt {
        /// Energy `e` relaxation rate (bulk viscosity).
        s_e: f64,
        /// Energy-square `ε` relaxation rate.
        s_eps: f64,
        /// Heat-flux `qx`, `qy` relaxation rate.
        s_q: f64,
    },
}

impl Collision {
    /// TRT with the wall-exact magic parameter `Λ = 3/16`.
    pub const TRT_WALL_EXACT: Collision = Collision::Trt { magic: 3.0 / 16.0 };

    /// MRT with the commonly used ghost-mode rates (Lallemand & Luo 2000).
    pub const MRT_DEFAULT: Collision = Collision::Mrt {
        s_e: 1.64,
        s_eps: 1.54,
        s_q: 1.9,
    };

    /// Diagonal moment-space rates `(ρ, e, ε, jx, qx, jy, qy, pxx, pxy)` for
    /// [`Collision::Mrt`]; `None` for the other operators.
    pub(crate) fn mrt_rates(self, tau: f64) -> Option<[f64; 9]> {
        match self {
            Collision::Mrt { s_e, s_eps, s_q } => {
                for r in [s_e, s_eps, s_q] {
                    assert!(r > 0.0 && r < 2.0, "MRT rates must lie in (0, 2)");
                }
                let w = 1.0 / tau;
                Some([0.0, s_e, s_eps, w, s_q, w, s_q, w, w])
            }
            _ => None,
        }
    }

    /// Relaxation rates `(ω⁺, ω⁻)` for the symmetric and antisymmetric parts
    /// at the viscosity-setting relaxation time `tau`.
    pub(crate) fn rates(self, tau: f64) -> (f64, f64) {
        match self {
            Collision::Bgk => (1.0 / tau, 1.0 / tau),
            Collision::Trt { magic } => {
                assert!(magic > 0.0, "TRT magic parameter must be positive");
                let tau_minus = 0.5 + magic / (tau - 0.5);
                (1.0 / tau, 1.0 / tau_minus)
            }
            Collision::Mrt { .. } => (1.0 / tau, 1.0 / tau),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bgk_rates_are_equal() {
        let (wp, wm) = Collision::Bgk.rates(0.8);
        assert_eq!(wp, wm);
        assert!((wp - 1.25).abs() < 1e-15);
    }

    #[test]
    fn trt_with_squared_offset_magic_is_bgk() {
        let tau: f64 = 0.7;
        let (wp, wm) = Collision::Trt {
            magic: (tau - 0.5) * (tau - 0.5),
        }
        .rates(tau);
        assert!((wp - wm).abs() < 1e-12, "{wp} vs {wm}");
    }
}
