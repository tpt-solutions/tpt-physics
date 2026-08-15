//! Material-point plasticity for the FEA crate.
//!
//! This is the net-new *material* nonlinearity complementing the geometric
//! (Total-Lagrangian) framework in [`crate::nonlinear`]. It implements an
//! associative **von Mises (J2) plasticity** model with **linear isotropic
//! hardening**, integrated by the standard elastic-predictor / plastic-corrector
//! (return mapping) algorithm.
//!
//! The model operates on a single integration point in Voigt notation
//! `[σxx, σyy, σzz, σxy, σyz, σxz]` with *engineering* shear strains/rates,
//! matching the convention used by [`crate::nonlinear`]. It is stress-driven
//! ([`PlasticMaterial::return_map`]) so it can be coupled into any element's
//! consistent-tangent driver, and also exposes a strain-driven convenience
//! wrapper ([`PlasticMaterial::update`]).
//!
//! ```
//! use tpt_physics_fea::plasticity::{PlasticMaterial, PlasticState};
//! let mat = PlasticMaterial::new(200e9, 0.3, 250e6, 1e9);
//! // Uniaxial tensile trial stress well past first yield.
//! let trial = [400e6, 0.0, 0.0, 0.0, 0.0, 0.0];
//! let (sigma, state) = mat.return_map(&trial, &PlasticState::default());
//! // Stress is capped at the (hardened) yield surface and plastic strain grew.
//! assert!(state.eq_plastic_strain > 0.0);
//! assert!((mat.von_mises(&sigma) - (mat.sigma_y0 + mat.hard * state.eq_plastic_strain)).abs() < 1.0);
//! ```

/// Isotropic von Mises material with linear isotropic hardening.
#[derive(Debug, Clone, Copy)]
pub struct PlasticMaterial {
    /// Young's modulus `E`.
    pub young: f64,
    /// Poisson's ratio `ν`.
    pub poisson: f64,
    /// Initial uniaxial yield stress `σ_y₀`.
    pub sigma_y0: f64,
    /// Isotropic hardening modulus `H` (slope dσ_y/dε̄ᵖ).
    pub hard: f64,
}

/// Per-integration-point plastic state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlasticState {
    /// Accumulated plastic strain (Voigt 6-vector, engineering shear).
    pub plastic_strain: [f64; 6],
    /// Accumulated equivalent (von Mises) plastic strain `ε̄ᵖ`.
    pub eq_plastic_strain: f64,
}

impl Default for PlasticState {
    fn default() -> Self {
        PlasticState {
            plastic_strain: [0.0; 6],
            eq_plastic_strain: 0.0,
        }
    }
}

impl PlasticMaterial {
    /// Construct a von Mises material with isotropic hardening.
    pub fn new(young: f64, poisson: f64, sigma_y0: f64, hard: f64) -> Self {
        PlasticMaterial {
            young,
            poisson,
            sigma_y0: sigma_y0,
            hard,
        }
    }

    /// Lamé parameters `(λ, μ)`.
    fn lame(&self) -> (f64, f64) {
        let mu = self.young / (2.0 * (1.0 + self.poisson));
        let lam = self.young * self.poisson / ((1.0 + self.poisson) * (1.0 - 2.0 * self.poisson));
        (lam, mu)
    }

    /// Shear modulus `μ`.
    fn shear(&self) -> f64 {
        self.young / (2.0 * (1.0 + self.poisson))
    }

    /// Isotropic elastic tangent `C` (6×6 Voigt, engineering shear).
    fn elastic_matrix(&self) -> [[f64; 6]; 6] {
        let (lam, mu) = self.lame();
        let mut c = [[0.0; 6]; 6];
        // Diagonal normal block.
        c[0][0] = lam + 2.0 * mu;
        c[1][1] = lam + 2.0 * mu;
        c[2][2] = lam + 2.0 * mu;
        // Off-diagonal normal coupling.
        c[0][1] = lam;
        c[0][2] = lam;
        c[1][0] = lam;
        c[1][2] = lam;
        c[2][0] = lam;
        c[2][1] = lam;
        // Shear block (engineering shear → `C₃₃ = μ`).
        c[3][3] = mu;
        c[4][4] = mu;
        c[5][5] = mu;
        c
    }

    /// Deviatoric part of a Voigt stress (shear components unchanged).
    fn deviatoric(sigma: &[f64; 6]) -> [f64; 6] {
        let p = (sigma[0] + sigma[1] + sigma[2]) / 3.0;
        [
            sigma[0] - p,
            sigma[1] - p,
            sigma[2] - p,
            sigma[3],
            sigma[4],
            sigma[5],
        ]
    }

    /// von Mises equivalent stress of a Voigt stress.
    pub fn von_mises(&self, sigma: &[f64; 6]) -> f64 {
        let s = Self::deviatoric(sigma);
        // s:s (tensor) = Σ s_ii² + 2 Σ s_ij², with engineering shear in s[3..].
        let ss = s[0] * s[0] + s[1] * s[1] + s[2] * s[2]
            + 2.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]);
        (1.5 * ss).sqrt()
    }

    /// Elastic-predictor / plastic-corrector return mapping for a trial stress.
    ///
    /// If the trial stress lies inside the (hardened) yield surface the state is
    /// unchanged and the stress is returned as-is. Otherwise the stress is
    /// projected back onto the yield surface along the von Mises flow
    /// direction and the plastic state is updated.
    pub fn return_map(&self, sigma_trial: &[f64; 6], prev: &PlasticState) -> ([f64; 6], PlasticState) {
        let q = self.von_mises(sigma_trial);
        let yield_stress = self.sigma_y0 + self.hard * prev.eq_plastic_strain;
        let phi = q - yield_stress;
        if phi <= 0.0 {
            return (*sigma_trial, *prev);
        }

        let g = self.shear();
        let s = Self::deviatoric(sigma_trial);
        let ss = s[0] * s[0] + s[1] * s[1] + s[2] * s[2]
            + 2.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]);
        let norm = ss.sqrt();
        let dg = phi / (3.0 * g + self.hard);
        let nscale = (1.5_f64).sqrt() / norm;

        let mut sigma = [0.0_f64; 6];
        let mut deps_p = [0.0_f64; 6];
        for i in 0..6 {
            let n_i = nscale * s[i];
            sigma[i] = sigma_trial[i] - 2.0 * g * dg * n_i;
            deps_p[i] = dg * n_i;
        }

        let mut state = *prev;
        for i in 0..6 {
            state.plastic_strain[i] += deps_p[i];
        }
        state.eq_plastic_strain += dg;
        (sigma, state)
    }

    /// Strain-driven update: elastic strain = `total − plastic`, trial stress
    /// via the elastic tangent, then [`PlasticMaterial::return_map`].
    pub fn update(&self, total_strain: &[f64; 6], prev: &PlasticState) -> ([f64; 6], PlasticState) {
        let c = self.elastic_matrix();
        let mut eps_e = [0.0_f64; 6];
        for i in 0..6 {
            eps_e[i] = total_strain[i] - prev.plastic_strain[i];
        }
        let sigma_trial = matvec(&c, &eps_e);
        self.return_map(&sigma_trial, prev)
    }
}

/// 6×6 · 6 matrix-vector product (row-major).
fn matvec(a: &[[f64; 6]; 6], x: &[f64; 6]) -> [f64; 6] {
    let mut y = [0.0_f64; 6];
    for i in 0..6 {
        let mut s = 0.0;
        for j in 0..6 {
            s += a[i][j] * x[j];
        }
        y[i] = s;
    }
    y
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn elastic_regime_is_linear_and_state_preserving() {
        let mat = PlasticMaterial::new(200e9, 0.3, 250e6, 1e9);
        // Trial stress well below yield (uniaxial 100 MPa < 250 MPa).
        let trial = [100e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let (sigma, state) = mat.return_map(&trial, &PlasticState::default());
        assert_eq!(sigma, trial);
        assert_eq!(state, PlasticState::default());
    }

    #[test]
    fn uniaxial_tension_saturates_and_accumulates_plastic_strain() {
        let mat = PlasticMaterial::new(200e9, 0.3, 250e6, 1e9);
        // Trial stress past first yield (400 MPa > 250 MPa).
        let trial = [400e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        let (sigma, state) = mat.return_map(&trial, &PlasticState::default());
        // Yield-surface consistency: von Mises == hardened yield stress.
        let q = mat.von_mises(&sigma);
        let yield_stress = mat.sigma_y0 + mat.hard * state.eq_plastic_strain;
        assert!((q - yield_stress).abs() < 1.0, "q={q}, yield={yield_stress}");
        // Stress is relieved below the trial.
        assert!(sigma[0] < trial[0]);
        // Plastic strain accumulated.
        assert!(state.eq_plastic_strain > 0.0);
        // Uniaxial direction: axial plastic strain positive, lateral negative.
        assert!(state.plastic_strain[0] > 0.0);
        assert!(state.plastic_strain[1] < 0.0);
        assert!(state.plastic_strain[2] < 0.0);
    }

    #[test]
    fn pure_shear_yields_at_sqrt3_times_yield() {
        let mat = PlasticMaterial::new(200e9, 0.3, 250e6, 1e9);
        // Pure shear τ; von Mises = √3 |τ|, so yield when τ = 250/√3 MPa.
        let tau = 200e6; // > 250/√3 ≈ 144 MPa
        let trial = [0.0, 0.0, 0.0, tau, 0.0, 0.0];
        let (sigma, state) = mat.return_map(&trial, &PlasticState::default());
        let q = mat.von_mises(&sigma);
        let yield_stress = mat.sigma_y0 + mat.hard * state.eq_plastic_strain;
        assert!((q - yield_stress).abs() < 1.0);
        assert!(state.eq_plastic_strain > 0.0);
        // Shear stress relieved.
        assert!((sigma[3]).abs() < tau.abs());
    }

    #[test]
    fn cyclic_hardening_surfaces_and_unloading_is_elastic() {
        let mat = PlasticMaterial::new(200e9, 0.3, 250e6, 1e9);
        // Load beyond yield, then unload within the hardened surface.
        let (s1, st1) = mat.return_map(&[400e6, 0.0, 0.0, 0.0, 0.0, 0.0], &PlasticState::default());
        let yield_after = mat.sigma_y0 + mat.hard * st1.eq_plastic_strain;
        assert!(mat.von_mises(&s1) <= yield_after + 1.0);

        // Unload to a smaller trial stress (well inside the hardened surface)
        // → no further plastic flow.
        let (s2, st2) = mat.return_map(&[100e6, 0.0, 0.0, 0.0, 0.0, 0.0], &st1);
        assert_eq!(s2, [100e6, 0.0, 0.0, 0.0, 0.0, 0.0]);
        assert_eq!(st2, st1);
    }

    #[test]
    fn strain_driven_matches_elastic_for_small_strain() {
        let mat = PlasticMaterial::new(200e9, 0.3, 250e6, 1e9);
        // Uniaxial *stress* of 100 MPa corresponds (elastically) to strain
        // ε = [σ/E, -νσ/E, -νσ/E, 0,0,0]; below yield → pure elastic.
        let s = 100e6;
        let e_strain = [s / mat.young, -mat.poisson * s / mat.young, -mat.poisson * s / mat.young, 0.0, 0.0, 0.0];
        let (sigma, state) = mat.update(&e_strain, &PlasticState::default());
        assert!((sigma[0] - s).abs() < 1.0);
        assert_eq!(state, PlasticState::default());
    }
}
