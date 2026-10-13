//! Lid-driven-cavity reference data (Ghia, Ghia & Shin, 1982) and helpers to
//! compare an [`Lbm2D`] cavity run against it.
//!
//! Ghia et al. tabulate the steady incompressible solution of the unit-square
//! cavity (lid moving in `+x` at `U = 1`) on a 129×129 grid: the `u`-velocity
//! along the vertical centreline `x = ½` and the `v`-velocity along the
//! horizontal centreline `y = ½`.
//!
//! # Geometry convention
//!
//! A cavity built by [`Lbm2D::cavity`] has solid wall *nodes* on the lattice
//! border, so with half-way bounce-back the physical walls sit half a cell
//! inside them. The fluid occupies `L = n − 2` cells and a node `i` (for
//! `1 ≤ i ≤ n − 2`) is at the physical coordinate `(i − ½) / L`. The helpers
//! here use that mapping, so profiles can be sampled at the exact Ghia
//! coordinates by linear interpolation (`wall` values are the no-slip / lid
//! speeds at coordinate 0 and 1).

use crate::Lbm2D;

/// One reference profile: `(coordinate, velocity / U_lid)` pairs, sorted by
/// descending coordinate as in the original tables.
pub type Profile = &'static [(f64, f64)];

/// `u / U` along the vertical centreline `x = ½` vs `y`, Re = 100.
pub const U_RE100: Profile = &[
    (1.0000, 1.00000),
    (0.9766, 0.84123),
    (0.9688, 0.78871),
    (0.9609, 0.73722),
    (0.9531, 0.68717),
    (0.8516, 0.23151),
    (0.7344, 0.00332),
    (0.6172, -0.13641),
    (0.5000, -0.20581),
    (0.4531, -0.21090),
    (0.2813, -0.15662),
    (0.1719, -0.10150),
    (0.1016, -0.06434),
    (0.0703, -0.04775),
    (0.0625, -0.04192),
    (0.0547, -0.03717),
    (0.0000, 0.00000),
];

/// `v / U` along the horizontal centreline `y = ½` vs `x`, Re = 100.
pub const V_RE100: Profile = &[
    (1.0000, 0.00000),
    (0.9688, -0.05906),
    (0.9609, -0.07391),
    (0.9531, -0.08864),
    (0.9453, -0.10313),
    (0.9063, -0.16914),
    (0.8594, -0.22445),
    (0.8047, -0.24533),
    (0.5000, 0.05454),
    (0.2344, 0.17527),
    (0.2266, 0.17507),
    (0.1563, 0.16077),
    (0.0938, 0.12317),
    (0.0781, 0.10890),
    (0.0703, 0.10091),
    (0.0625, 0.09233),
    (0.0000, 0.00000),
];

/// `u / U` along the vertical centreline `x = ½` vs `y`, Re = 400.
pub const U_RE400: Profile = &[
    (1.0000, 1.00000),
    (0.9766, 0.75837),
    (0.9688, 0.68439),
    (0.9609, 0.61756),
    (0.9531, 0.55892),
    (0.8516, 0.29093),
    (0.7344, 0.16256),
    (0.6172, 0.02135),
    (0.5000, -0.11477),
    (0.4531, -0.17119),
    (0.2813, -0.32726),
    (0.1719, -0.24299),
    (0.1016, -0.14612),
    (0.0703, -0.10338),
    (0.0625, -0.09266),
    (0.0547, -0.08186),
    (0.0000, 0.00000),
];

/// `v / U` along the horizontal centreline `y = ½` vs `x`, Re = 400.
pub const V_RE400: Profile = &[
    (1.0000, 0.00000),
    (0.9688, -0.12146),
    (0.9609, -0.15663),
    (0.9531, -0.19254),
    (0.9453, -0.22847),
    (0.9063, -0.33714),
    (0.8594, -0.44993),
    (0.8047, -0.38598),
    (0.5000, 0.05186),
    (0.2344, 0.30174),
    (0.2266, 0.30203),
    (0.1563, 0.28124),
    (0.0938, 0.22965),
    (0.0781, 0.20920),
    (0.0703, 0.19713),
    (0.0625, 0.18360),
    (0.0000, 0.00000),
];

/// Piecewise-linear value of a quantity along one lattice line at the
/// physical coordinate `p ∈ [0, 1]`. `get(i)` is the node value for the fluid
/// nodes `1 ..= n-2`; `w0` / `w1` are the wall values at `p = 0` / `p = 1`
/// (half a cell outside the first / last fluid node).
fn line_value(n: usize, p: f64, get: impl Fn(usize) -> f64, w0: f64, w1: f64) -> f64 {
    let l = (n - 2) as f64;
    // Node `i` sits at s = i; the walls sit at s = ½ and s = n − 3/2.
    let s = (p * l + 0.5).clamp(0.5, n as f64 - 1.5);
    if s <= 1.0 {
        let t = (s - 0.5) / 0.5;
        w0 + (get(1) - w0) * t
    } else if s >= (n - 2) as f64 {
        let t = (s - (n - 2) as f64) / 0.5;
        get(n - 2) + (w1 - get(n - 2)) * t
    } else {
        let i0 = s.floor() as usize;
        let t = s - i0 as f64;
        get(i0) * (1.0 - t) + get(i0 + 1) * t
    }
}

/// Bilinear sample of a cavity field at physical `(x, y) ∈ [0, 1]²`.
/// `field` selects the lattice array; `wall` is its value on the three
/// stationary walls and `lid` on the top wall.
fn sample(sim: &Lbm2D, field: &[f64], x: f64, y: f64, wall: f64, lid: f64) -> f64 {
    let (nx, ny) = (sim.nx, sim.ny);
    // Interpolate along x on each of the two bracketing fluid rows, then
    // along y (the wall values enter only through the y interpolation; the
    // side-wall values enter through the x interpolation).
    let row = |iy: usize| line_value(nx, x, |ix| field[iy * nx + ix], wall, wall);
    line_value(ny, y, row, wall, lid)
}

/// `u / U_lid` of a cavity run at physical `(x, y)`.
pub fn u_at(sim: &Lbm2D, u_lid: f64, x: f64, y: f64) -> f64 {
    sample(sim, &sim.ux, x, y, 0.0, u_lid) / u_lid
}

/// `v / U_lid` of a cavity run at physical `(x, y)`.
pub fn v_at(sim: &Lbm2D, u_lid: f64, x: f64, y: f64) -> f64 {
    sample(sim, &sim.uy, x, y, 0.0, 0.0) / u_lid
}

/// Error of a cavity run against a Ghia profile.
#[derive(Debug, Clone, Copy)]
pub struct ProfileError {
    /// Root-mean-square difference over the interior reference points.
    pub rms: f64,
    /// Largest absolute difference over the interior reference points.
    pub max: f64,
}

fn compare(reference: Profile, mut model: impl FnMut(f64) -> f64) -> ProfileError {
    let mut sum = 0.0;
    let mut max = 0.0_f64;
    let mut count = 0usize;
    for &(p, r) in reference {
        // The end points are the boundary conditions themselves.
        if p <= 0.0 || p >= 1.0 {
            continue;
        }
        let d = model(p) - r;
        sum += d * d;
        max = max.max(d.abs());
        count += 1;
    }
    ProfileError {
        rms: (sum / count as f64).sqrt(),
        max,
    }
}

/// Error of the vertical-centreline `u` profile against `reference` (use
/// [`U_RE100`] or [`U_RE400`]).
pub fn u_centerline_error(sim: &Lbm2D, u_lid: f64, reference: Profile) -> ProfileError {
    compare(reference, |y| u_at(sim, u_lid, 0.5, y))
}

/// Error of the horizontal-centreline `v` profile against `reference` (use
/// [`V_RE100`] or [`V_RE400`]).
pub fn v_centerline_error(sim: &Lbm2D, u_lid: f64, reference: Profile) -> ProfileError {
    compare(reference, |x| v_at(sim, u_lid, x, 0.5))
}
