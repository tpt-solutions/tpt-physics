//! The D2Q9 lattice definition: velocity set, weights, and equilibrium.
//!
//! All constants use the standard lattice-unit convention `cs² = 1/3`, `Δx =
//! Δt = 1`.

/// The 9-velocity D2Q9 lattice.
pub struct D2Q9;

impl D2Q9 {
    /// Squared lattice sound speed `cs² = 1/3`.
    pub const CS2: f64 = 1.0 / 3.0;

    /// The nine discrete velocities `(ex, ey)`, ordered
    /// `0:rest, 1:E, 2:N, 3:W, 4:S, 5:NE, 6:NW, 7:SW, 8:SE`.
    pub const E: [(i32, i32); 9] = [
        (0, 0),
        (1, 0),
        (0, 1),
        (-1, 0),
        (0, -1),
        (1, 1),
        (-1, 1),
        (-1, -1),
        (1, -1),
    ];

    /// The nine lattice weights.
    pub const W: [f64; 9] = [
        4.0 / 9.0,
        1.0 / 9.0,
        1.0 / 9.0,
        1.0 / 9.0,
        1.0 / 9.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0,
    ];

    /// Index of the opposite velocity for each direction (used for bounce-back).
    pub const OPP: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];

    /// Second-order equilibrium distribution `f_i^eq(ρ, u)`.
    ///
    /// `f_i^eq = w_i ρ (1 + 3 e·u + 4.5 (e·u)² − 1.5 |u|²)`.
    #[inline]
    pub fn equilibrium(rho: f64, u: [f64; 2]) -> [f64; 9] {
        let u2 = u[0] * u[0] + u[1] * u[1];
        let mut feq = [0.0_f64; 9];
        for (i, slot) in feq.iter_mut().enumerate() {
            let (ex, ey) = D2Q9::E[i];
            let eu = ex as f64 * u[0] + ey as f64 * u[1];
            *slot = D2Q9::W[i] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * u2);
        }
        feq
    }
}

/// Velocity-set description shared by every lattice, so solvers and tests can
/// be written against the *model* rather than the concrete [`D2Q9`] type.
///
/// Velocities are padded to three components (`[ex, ey, ez]`, `ez = 0` in 2-D)
/// so a future D3Q19 / D3Q27 fits the same trait. [`crate::Lbm2D`] is
/// currently specialised to [`D2Q9`]; this trait is the seam for generalising
/// it, and it already lets the isotropy conditions be checked uniformly.
pub trait LatticeModel {
    /// Spatial dimension `D`.
    const D: usize;
    /// Number of discrete velocities `Q`.
    const Q: usize;
    /// Squared lattice sound speed `cs²`.
    const CS2: f64;

    /// The `Q` discrete velocities `[ex, ey, ez]`.
    fn velocities() -> &'static [[i32; 3]];
    /// The `Q` quadrature weights.
    fn weights() -> &'static [f64];
    /// For each velocity, the index of its opposite (`e_opp = −e`).
    fn opposites() -> &'static [usize];

    /// Second-order equilibrium `f_k^eq(ρ, u)` for velocity `k`; `u` holds
    /// `D` components.
    fn equilibrium_k(k: usize, rho: f64, u: &[f64]) -> f64 {
        let e = Self::velocities()[k];
        let mut eu = 0.0;
        let mut u2 = 0.0;
        for (d, &ud) in u.iter().enumerate().take(Self::D) {
            eu += e[d] as f64 * ud;
            u2 += ud * ud;
        }
        let c = 1.0 / Self::CS2;
        Self::weights()[k] * rho * (1.0 + c * eu + 0.5 * c * c * eu * eu - 0.5 * c * u2)
    }
}

const D2Q9_E3: [[i32; 3]; 9] = [
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [-1, 0, 0],
    [0, -1, 0],
    [1, 1, 0],
    [-1, 1, 0],
    [-1, -1, 0],
    [1, -1, 0],
];

impl LatticeModel for D2Q9 {
    const D: usize = 2;
    const Q: usize = 9;
    const CS2: f64 = D2Q9::CS2;

    fn velocities() -> &'static [[i32; 3]] {
        &D2Q9_E3
    }
    fn weights() -> &'static [f64] {
        &D2Q9::W
    }
    fn opposites() -> &'static [usize] {
        &D2Q9::OPP
    }
}

/// The 19-velocity three-dimensional lattice (rest, 6 face, 12 edge
/// neighbours). Provided as lattice data behind [`LatticeModel`]; there is no
/// 3-D solver yet.
pub struct D3Q19;

const D3Q19_E: [[i32; 3]; 19] = [
    [0, 0, 0],
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [0, 0, 1],
    [0, 0, -1],
    [1, 1, 0],
    [-1, -1, 0],
    [1, -1, 0],
    [-1, 1, 0],
    [1, 0, 1],
    [-1, 0, -1],
    [1, 0, -1],
    [-1, 0, 1],
    [0, 1, 1],
    [0, -1, -1],
    [0, 1, -1],
    [0, -1, 1],
];

const D3Q19_W: [f64; 19] = [
    1.0 / 3.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

/// Opposite-velocity indices: entries come in `(k, k+1)` pairs after rest.
const D3Q19_OPP: [usize; 19] = [
    0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17,
];

impl LatticeModel for D3Q19 {
    const D: usize = 3;
    const Q: usize = 19;
    const CS2: f64 = 1.0 / 3.0;

    fn velocities() -> &'static [[i32; 3]] {
        &D3Q19_E
    }
    fn weights() -> &'static [f64] {
        &D3Q19_W
    }
    fn opposites() -> &'static [usize] {
        &D3Q19_OPP
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The isotropy conditions every valid lattice must satisfy:
    /// `Σw = 1`, `Σw e = 0`, `Σw e_a e_b = cs² δ_ab`, and `OPP` negates `e`.
    fn check_isotropy<L: LatticeModel>() {
        let (e, w, opp) = (L::velocities(), L::weights(), L::opposites());
        assert_eq!((e.len(), w.len(), opp.len()), (L::Q, L::Q, L::Q));
        assert!((w.iter().sum::<f64>() - 1.0).abs() < 1e-14);
        for a in 0..L::D {
            let m1: f64 = (0..L::Q).map(|k| w[k] * e[k][a] as f64).sum();
            assert!(m1.abs() < 1e-14, "first moment {a}: {m1}");
            for b in 0..L::D {
                let m2: f64 = (0..L::Q).map(|k| w[k] * (e[k][a] * e[k][b]) as f64).sum();
                let want = if a == b { L::CS2 } else { 0.0 };
                assert!((m2 - want).abs() < 1e-14, "second moment {a}{b}: {m2}");
            }
        }
        for (k, ek) in e.iter().enumerate() {
            for (a, &x) in ek.iter().enumerate() {
                assert_eq!(e[opp[k]][a], -x);
            }
        }
    }

    #[test]
    fn d2q9_is_isotropic() {
        check_isotropy::<D2Q9>();
    }

    #[test]
    fn d3q19_is_isotropic() {
        check_isotropy::<D3Q19>();
    }

    #[test]
    fn d3q19_equilibrium_conserves_mass_and_momentum() {
        let (rho, u) = (1.1, [0.05, -0.02, 0.03]);
        let (mut m0, mut m1) = (0.0, [0.0; 3]);
        for k in 0..D3Q19::Q {
            let f = D3Q19::equilibrium_k(k, rho, &u);
            m0 += f;
            for (a, m) in m1.iter_mut().enumerate() {
                *m += f * D3Q19::velocities()[k][a] as f64;
            }
        }
        assert!((m0 - rho).abs() < 1e-12);
        for a in 0..3 {
            assert!((m1[a] - rho * u[a]).abs() < 1e-12, "axis {a}");
        }
    }

    #[test]
    fn trait_equilibrium_matches_inherent() {
        let (rho, u) = (1.2, [0.07, -0.04]);
        let feq = D2Q9::equilibrium(rho, u);
        for (k, f) in feq.iter().enumerate() {
            let g = <D2Q9 as LatticeModel>::equilibrium_k(k, rho, &u);
            assert!((f - g).abs() < 1e-14, "k={k}: {f} vs {g}");
        }
    }
}
