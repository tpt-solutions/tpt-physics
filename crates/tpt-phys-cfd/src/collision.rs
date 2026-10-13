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
}

impl Collision {
    /// TRT with the wall-exact magic parameter `Λ = 3/16`.
    pub const TRT_WALL_EXACT: Collision = Collision::Trt { magic: 3.0 / 16.0 };

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
