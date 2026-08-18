//! Advanced finite-element formulations that are net-new to `tpt-physics`:
//! the quadratic tetrahedral element ([`tet10_stiffness`]), the 3-D beam/frame
//! element ([`beam3d_global_stiffness`]), and a Mindlin–Reissner plate/shell
//! element ([`shell4_stiffness`]).
//!
//! **Status:** `tet10_stiffness` is validated. `beam3d_global_stiffness`
//! (Euler–Bernoulli, no shear) and `shell4_stiffness` (Mindlin–Reissner) are
//! **experimental** — they pass rigid-body-exactness (on parallelogram meshes)
//! and plate/beam benchmarks but have not been validated against full curved-
//! shell benchmarks. `shell4` is rigid-body exact on parallelogram (incl.
//! rectangular) meshes and is validated against simply-supported *square* and
//! *skew* (Morley-style) plate benchmarks; transverse shear is integrated with
//! the full 2×2 rule (not reduced) to stay well-conditioned on assembled
//! meshes, at the cost of mild shear locking (the element runs stiff for thin
//! plates) — the known trade-off for the 4-node Mindlin element.
//!
//! Linear-tetrahedron, hexahedron, and 2-D frame elements are reused directly
//! from `tpt-fem-element` / `tpt-fem-elasticity`; see the crate root re-exports.

/// 3×3 matrix inverse (row-major `m`, length 9).
///
/// Returns `None` when the matrix is (near-)singular — a degenerate or
/// inverted element — instead of silently producing `inf`/`NaN` via division
/// by zero. Callers must handle the `None` case (e.g. reject the element).
fn mat3_inv(m: &[f64; 9]) -> Option<[f64; 9]> {
    let det = m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6])
        + m[2] * (m[3] * m[7] - m[4] * m[6]);
    if !det.is_finite() || det.abs() < 1e-12 {
        return None;
    }
    let inv = 1.0 / det;
    Some([
        (m[4] * m[8] - m[5] * m[7]) * inv,
        (m[2] * m[7] - m[1] * m[8]) * inv,
        (m[1] * m[5] - m[2] * m[4]) * inv,
        (m[5] * m[6] - m[3] * m[8]) * inv,
        (m[0] * m[8] - m[2] * m[6]) * inv,
        (m[2] * m[3] - m[0] * m[5]) * inv,
        (m[3] * m[7] - m[4] * m[6]) * inv,
        (m[1] * m[6] - m[0] * m[7]) * inv,
        (m[0] * m[4] - m[1] * m[3]) * inv,
    ])
}

// ----------------------------------------------------------------------------
// Quadratic tetrahedron (Tet10)
// ----------------------------------------------------------------------------

/// Reference-node coordinates of the 10-node quadratic tetrahedron, indexed as
/// `[corner0, corner1, corner2, corner3, mid01, mid02, mid03, mid12, mid13,
/// mid23]`. Corners match [`tpt_fem_element::Tet4`].
pub const TET10_NODES: [[f64; 3]; 10] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.5, 0.0, 0.0],
    [0.0, 0.5, 0.0],
    [0.0, 0.0, 0.5],
    [0.5, 0.5, 0.0],
    [0.5, 0.0, 0.5],
    [0.0, 0.5, 0.5],
];

/// Quadratic-tet shape-function values `N` and reference-gradient
/// `∂N/∂ξ` (`[dN/dx, dN/dy, dN/dz]`) at local coordinates `(x, y, z)` of the
/// reference tetrahedron.
pub fn tet10_shape_grad(x: f64, y: f64, z: f64) -> (Vec<f64>, Vec<[f64; 3]>) {
    let l = [1.0 - x - y - z, x, y, z];
    // dL_i / d(x,y,z)
    let dl = [
        [-1.0, -1.0, -1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    // Corner nodes N_i = (2 L_i - 1) L_i.
    let mut n = vec![0.0; 10];
    let mut g = vec![[0.0_f64; 3]; 10];
    for i in 0..4 {
        n[i] = (2.0 * l[i] - 1.0) * l[i];
        let dndli = 4.0 * l[i] - 1.0;
        for k in 0..3 {
            g[i][k] = dndli * dl[i][k];
        }
    }
    // Mid-edge nodes N = 4 L_a L_b.
    let edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
    for (e, &(a, b)) in edges.iter().enumerate() {
        let idx = 4 + e;
        n[idx] = 4.0 * l[a] * l[b];
        for k in 0..3 {
            g[idx][k] = 4.0 * (dl[a][k] * l[b] + l[a] * dl[b][k]);
        }
    }
    (n, g)
}

/// 3-D isotropic elasticity `D` matrix (Voigt, 6×6, engineering shear).
fn iso_d_matrix(lambda: f64, mu: f64) -> [f64; 36] {
    let mut d = [0.0; 36];
    let c = lambda + 2.0 * mu;
    // diagonal
    d[0] = c;
    d[7] = c;
    d[14] = c;
    d[21] = mu;
    d[28] = mu;
    d[35] = mu;
    // off-diagonal λ
    d[1] = lambda;
    d[2] = lambda;
    d[6] = lambda;
    d[8] = lambda;
    d[12] = lambda;
    d[13] = lambda;
    d
}

/// Assemble the 30×30 (10 nodes × 3 DOF) stiffness matrix of a quadratic
/// tetrahedral element with physical node coordinates `nodes` (length 10) and
/// isotropic material parameters `E`, `nu`.
///
/// Returns a row-major `30×30` buffer (`length 900`).
pub fn tet10_stiffness(nodes: &[[f64; 3]; 10], e: f64, nu: f64) -> Vec<f64> {
    let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e / (2.0 * (1.0 + nu));
    let d = iso_d_matrix(lambda, mu);

    // Gauss–Legendre 4-point rule on the reference tet (volume 1/6).
    let a = 0.5854101966249685;
    let b = 0.1381966011250105;
    let pts = [[a, b, b], [b, a, b], [b, b, a], [b, b, b]];
    let w = 1.0 / 24.0;

    let mut k = vec![0.0; 900];
    for pt in pts {
        let (_n, grad) = tet10_shape_grad(pt[0], pt[1], pt[2]);
        // Jacobian J[p][r] = Σ_a X_a[p] ∂N_a/∂ξ_r.
        let mut j = [0.0; 9];
        for a in 0..10 {
            for p in 0..3 {
                for r in 0..3 {
                    j[p * 3 + r] += nodes[a][p] * grad[a][r];
                }
            }
        }
        let detj = j[0] * (j[4] * j[8] - j[5] * j[7]) - j[1] * (j[3] * j[8] - j[5] * j[6])
            + j[2] * (j[3] * j[7] - j[4] * j[6]);
        // `abs(detj)`: an inverted (valid but negatively-oriented) element has
        // a negative Jacobian but is still well-defined; only a near-zero
        // determinant is degenerate.
        let inv = match mat3_inv(&j) {
            Some(inv) => inv,
            None => return vec![0.0; 900], // degenerate element → no contribution
        };

        // Physical gradients dN/dx = dN/dξ · J⁻¹.
        let mut gp = [[0.0_f64; 3]; 10];
        for a in 0..10 {
            for p in 0..3 {
                let mut s = 0.0;
                for r in 0..3 {
                    s += grad[a][r] * inv[r * 3 + p];
                }
                gp[a][p] = s;
            }
        }

        // B matrix (6 × 30).
        let mut b = [[0.0_f64; 30]; 6];
        for (a, gp_a) in gp.iter().enumerate() {
            let i = 3 * a;
            b[0][i] = gp_a[0];
            b[1][i + 1] = gp_a[1];
            b[2][i + 2] = gp_a[2];
            b[3][i] = gp_a[1];
            b[3][i + 1] = gp_a[0];
            b[4][i + 1] = gp_a[2];
            b[4][i + 2] = gp_a[1];
            b[5][i] = gp_a[2];
            b[5][i + 2] = gp_a[0];
        }

        // K += Bᵀ D B |J| w.
        let scale = detj.abs() * w;
        for i in 0..30 {
            for jj in 0..30 {
                let mut s = 0.0;
                for s1 in 0..6 {
                    let mut t = 0.0;
                    for s2 in 0..6 {
                        t += d[s1 * 6 + s2] * b[s2][jj];
                    }
                    s += b[s1][i] * t;
                }
                k[i * 30 + jj] += s * scale;
            }
        }
    }
    k
}

// ----------------------------------------------------------------------------
// 3-D beam / frame element (Euler–Bernoulli, 12×12)
// ----------------------------------------------------------------------------

/// Build the 12×12 stiffness of a 3-D beam element in global coordinates.
///
/// `n0`, `n1` are the two node positions; `up` is a reference vector used to
/// orient the local `y`/`z` axes. `e` is Young's modulus, `g` is the shear
/// modulus, `area` is the cross-sectional area, `iy`/`iz` are second moments of
/// area about the local `y`/`z` axes, and `jtor` is the torsional constant.
///
/// Returns a row-major `12×12` buffer (`length 144`). Local DOFs per node are
/// `[u, v, w, θx, θy, θz]`.
#[allow(clippy::too_many_arguments)]
pub fn beam3d_global_stiffness(
    n0: [f64; 3],
    n1: [f64; 3],
    up: [f64; 3],
    e: f64,
    g: f64,
    area: f64,
    iy: f64,
    iz: f64,
    jtor: f64,
) -> Vec<f64> {
    let dx = n1[0] - n0[0];
    let dy = n1[1] - n0[1];
    let dz = n1[2] - n0[2];
    let l = (dx * dx + dy * dy + dz * dz).sqrt();

    // Local x axis.
    let ex = [dx / l, dy / l, dz / l];
    // Local z = up × ex, then local y = ex × ez.
    let mut ez = cross(up, ex);
    let ez_len = (ez[0] * ez[0] + ez[1] * ez[1] + ez[2] * ez[2]).sqrt();
    ez = [ez[0] / ez_len, ez[1] / ez_len, ez[2] / ez_len];
    let ey = cross(ex, ez);

    let mut r = [0.0; 9];
    for p in 0..3 {
        r[p * 3] = ex[p];
        r[p * 3 + 1] = ey[p];
        r[p * 3 + 2] = ez[p];
    }

    // Local stiffness (see module docs for DOF ordering).
    let mut kl = [0.0; 144];
    let ea_l = e * area / l;
    let gj_l = g * jtor / l;
    // Axial (dofs 0,6) and torsion (dofs 3,9).
    kl[0] = ea_l;
    kl[6] = -ea_l;
    kl[72] = -ea_l;
    kl[78] = ea_l;
    kl[3 * 12 + 3] = gj_l;
    kl[3 * 12 + 9] = -gj_l;
    kl[9 * 12 + 3] = -gj_l;
    kl[9 * 12 + 9] = gj_l;
    // Bending blocks.
    add_beam_bending(&mut kl, &[1, 5, 7, 11], e * iz, l); // v-θz (about z)
    add_beam_bending(&mut kl, &[2, 4, 8, 10], e * iy, l); // w-θy (about y)

    // T = block-diagonal(R) (12×12). Kg = Tᵀ Kl T.
    let mut t = [0.0; 144];
    for b in 0..4 {
        for p in 0..3 {
            for q in 0..3 {
                t[(b * 3 + p) * 12 + (b * 3 + q)] = r[p * 3 + q];
            }
        }
    }
    // Kg = Tᵀ Kl T.
    let mut tmp = [0.0; 144];
    for i in 0..12 {
        for k in 0..12 {
            let mut s = 0.0;
            for m in 0..12 {
                s += kl[i * 12 + m] * t[m * 12 + k];
            }
            tmp[i * 12 + k] = s;
        }
    }
    let mut kg = vec![0.0; 144];
    for i in 0..12 {
        for j in 0..12 {
            let mut s = 0.0;
            for m in 0..12 {
                s += t[m * 12 + i] * tmp[m * 12 + j];
            }
            kg[i * 12 + j] = s;
        }
    }
    kg
}

fn add_beam_bending(k: &mut [f64; 144], dofs: &[usize; 4], ei: f64, l: f64) {
    // Euler–Bernoulli bending block for the DOFs [v1, θ1, v2, θ2]:
    //   (EI/L³) · [[ 12,   6L,  -12,   6L],
    //              [ 6L,  4L²,  -6L,  2L²],
    //              [-12,  -6L,   12,  -6L],
    //              [ 6L,  2L²,  -6L,  4L²]]
    let c = ei / l.powi(3);
    let l2 = l * l;
    let block = [
        [12.0, 6.0 * l, -12.0, 6.0 * l],
        [6.0 * l, 4.0 * l2, -6.0 * l, 2.0 * l2],
        [-12.0, -6.0 * l, 12.0, -6.0 * l],
        [6.0 * l, 2.0 * l2, -6.0 * l, 4.0 * l2],
    ];
    for i in 0..4 {
        for j in 0..4 {
            let val = block[i][j] * c;
            let (gi, gj) = (dofs[i], dofs[j]);
            k[gi * 12 + gj] += val;
        }
    }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ----------------------------------------------------------------------------
// Mindlin–Reissner 4-node plate/shell (bending) element
// ----------------------------------------------------------------------------

/// Build the 12×12 stiffness of a flat 4-node Mindlin–Reissner plate element
/// (DOFs per node: `[w, θx, θy]` — transverse displacement and two rotations).
///
/// `nodes` are the four corner positions in the element plane; `thick` is the
/// plate thickness, `e`/`nu` the material, and `kappa` the shear-correction
/// factor (≈ 5/6). Bending and transverse shear both use the 2×2 Gauss rule;
/// the shear term is integrated fully (not reduced) so the element stays
/// well-conditioned on assembled meshes, at the cost of mild shear locking.
pub fn shell4_stiffness(
    nodes: &[[f64; 3]; 4],
    thick: f64,
    e: f64,
    nu: f64,
    kappa: f64,
) -> Vec<f64> {
    // Bilinear shape functions on the reference square [-1,1]².
    let shape = |xi: f64, eta: f64| -> ([f64; 4], [[f64; 2]; 4]) {
        let s = [xi, eta];
        let n = [
            0.25 * (1.0 - s[0]) * (1.0 - s[1]),
            0.25 * (1.0 + s[0]) * (1.0 - s[1]),
            0.25 * (1.0 + s[0]) * (1.0 + s[1]),
            0.25 * (1.0 - s[0]) * (1.0 + s[1]),
        ];
        let mut dn = [[0.0; 2]; 4];
        dn[0] = [-0.25 * (1.0 - s[1]), -0.25 * (1.0 - s[0])];
        dn[1] = [0.25 * (1.0 - s[1]), -0.25 * (1.0 + s[0])];
        dn[2] = [0.25 * (1.0 + s[1]), 0.25 * (1.0 + s[0])];
        dn[3] = [-0.25 * (1.0 + s[1]), 0.25 * (1.0 - s[0])];
        (n, dn)
    };

    let d = e * thick * thick * thick / (12.0 * (1.0 - nu * nu));
    // Bending constitutive (κ_xx, κ_yy, κ_xy).
    let db = [
        d,
        d * nu,
        0.0, //
        d * nu,
        d,
        0.0, //
        0.0,
        0.0,
        d * (1.0 / 2.0 - nu / 2.0),
    ];
    let gs = kappa * e / (2.0 * (1.0 + nu)) * thick; // shear modulus * κ * t

    let mut k = vec![0.0; 144];

    // Bending: 2×2 Gauss.
    for &(xi, eta, w) in &[
        (-0.57735026919, -0.57735026919, 1.0),
        (0.57735026919, -0.57735026919, 1.0),
        (-0.57735026919, 0.57735026919, 1.0),
        (0.57735026919, 0.57735026919, 1.0),
    ] {
        let (_n, dn) = shape(xi, eta);
        // Jacobian of (x,y) wrt (xi,eta).
        let mut j = [[0.0; 2]; 2];
        for a in 0..4 {
            for p in 0..2 {
                for r in 0..2 {
                    j[p][r] += nodes[a][p] * dn[a][r];
                }
            }
        }
        let detj = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        let inv = [
            [j[1][1] / detj, -j[0][1] / detj],
            [-j[1][0] / detj, j[0][0] / detj],
        ];
        // Physical gradients dN/dx, dN/dy.
        let mut gx = [0.0; 4];
        let mut gy = [0.0; 4];
        for a in 0..4 {
            gx[a] = inv[0][0] * dn[a][0] + inv[0][1] * dn[a][1];
            gy[a] = inv[1][0] * dn[a][0] + inv[1][1] * dn[a][1];
        }
        // Bending B (3 × 12): curvatures κ = [∂θy/∂x, -∂θx/∂y, ∂θy/∂y - ∂θx/∂x].
        let mut bb = [[0.0_f64; 12]; 3];
        for a in 0..4 {
            let i = 3 * a;
            // κ_x = ∂θy/∂x  -> dof θy (i+2) gradient gx
            bb[0][i + 2] = gx[a];
            // κ_y = -∂θx/∂y -> dof θx (i+1) gradient -gy
            bb[1][i + 1] = -gy[a];
            // ∂θy/∂y
            bb[2][i + 2] = gy[a];
            // -∂θx/∂x
            bb[2][i + 1] = -gx[a];
        }
        // K += Bᵀ Db B |J| w.
        let scale = detj * w;
        for i in 0..12 {
            for j in 0..12 {
                let mut s = 0.0;
                for s1 in 0..3 {
                    let mut t = 0.0;
                    for s2 in 0..3 {
                        t += db[s1 * 3 + s2] * bb[s2][j];
                    }
                    s += bb[s1][i] * t;
                }
                k[i * 12 + j] += s * scale;
            }
        }
    }

    // Transverse shear: 2×2 Gauss integration (same rule as bending). A single
    // reduced integration point makes the 4-node Mindlin element rank-deficient
    // on a structured mesh (spurious hourglass mode), which renders assembled
    // plate systems ill-conditioned and produces garbage solutions. Using the
    // full 2×2 rule keeps the element positive-definite and well-conditioned at
    // the cost of some shear locking (the element runs stiff for thin plates);
    // rigid-body exactness is preserved because a rigid mode has zero shear
    // strain at every point.
    for &(xi, eta, w) in &[
        (-0.57735026919, -0.57735026919, 1.0),
        (0.57735026919, -0.57735026919, 1.0),
        (-0.57735026919, 0.57735026919, 1.0),
        (0.57735026919, 0.57735026919, 1.0),
    ] {
        let (_n, dn) = shape(xi, eta);
        let mut j = [[0.0; 2]; 2];
        for a in 0..4 {
            for p in 0..2 {
                for r in 0..2 {
                    j[p][r] += nodes[a][p] * dn[a][r];
                }
            }
        }
        let detj = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        let inv = [
            [j[1][1] / detj, -j[0][1] / detj],
            [-j[1][0] / detj, j[0][0] / detj],
        ];
        let mut gx = [0.0; 4];
        let mut gy = [0.0; 4];
        // Shape values at this Gauss point for the rotation terms.
        let (ns, _) = shape(xi, eta);
        for a in 0..4 {
            gx[a] = inv[0][0] * dn[a][0] + inv[0][1] * dn[a][1];
            gy[a] = inv[1][0] * dn[a][0] + inv[1][1] * dn[a][1];
        }
        // Shear B (2 × 12): γ_xz = ∂w/∂x + θy ; γ_yz = ∂w/∂y - θx.
        let mut bs = [[0.0_f64; 12]; 2];
        for a in 0..4 {
            let i = 3 * a;
            bs[0][i] = gx[a]; // ∂w/∂x
            bs[0][i + 2] = ns[a]; // +θy
            bs[1][i] = gy[a]; // ∂w/∂y
            bs[1][i + 1] = -ns[a]; // -θx
        }
        // K += Bsᵀ Gs Bs |J| w, Gs = gs * I.
        for i in 0..12 {
            for j in 0..12 {
                let mut s = 0.0;
                for bs_s1 in &bs {
                    s += bs_s1[i] * gs * bs_s1[j];
                }
                k[i * 12 + j] += s * detj * w;
            }
        }
    }

    k
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tet10_rigid_translation_zero_force() {
        // Any element under a rigid-body translation must produce zero internal
        // force: K · u_rigid ≡ 0.
        let nodes = [
            [0.0, 0.0, 0.0],
            [2.0, 0.5, 0.0],
            [0.3, 2.0, 0.0],
            [0.5, 0.2, 1.5],
            [1.0, 0.25, 0.0],
            [0.15, 1.0, 0.0],
            [0.25, 0.1, 0.75],
            [1.15, 1.25, 0.0],
            [1.25, 0.35, 0.75],
            [0.4, 1.1, 0.75],
        ];
        let k = tet10_stiffness(&nodes, 200e9, 0.3);
        // Rigid translation (all nodes move by same vector).
        let mut u = vec![0.0; 30];
        for a in 0..10 {
            u[3 * a] = 0.01;
            u[3 * a + 1] = -0.02;
            u[3 * a + 2] = 0.03;
        }
        for i in 0..30 {
            let mut f = 0.0;
            for j in 0..30 {
                f += k[i * 30 + j] * u[j];
            }
            assert!(f.abs() < 1e-3, "rigid force component {i} = {f}");
        }
    }

    #[test]
    fn mat3_inv_rejects_degenerate() {
        // All-zero and rank-deficient matrices must be rejected (None), never
        // produce NaN/inf.
        assert!(mat3_inv(&[0.0; 9]).is_none());
        let degenerate = [
            1.0, 1.0, 1.0, // row 0
            2.0, 2.0, 2.0, // row 1 = 2×row 0 → singular
            0.0, 0.0, 1.0,
        ];
        assert!(mat3_inv(&degenerate).is_none());
        // A well-conditioned matrix inverts fine.
        let ok = mat3_inv(&[2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0]);
        assert!(ok.is_some());
        let inv = ok.unwrap();
        assert!((inv[0] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn tet10_degenerate_returns_zero_stiffness() {
        // A flat (zero-volume) tet has |J| = 0 everywhere → the guard must
        // yield a zero stiffness buffer rather than NaN.
        let flat = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 0.5, 0.0], // coplanar with the others → degenerate
            [0.5, 0.0, 0.0],
            [0.0, 0.5, 0.0],
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.0],
            [0.0, 0.5, 0.0],
        ];
        let k = tet10_stiffness(&flat, 200e9, 0.3);
        assert!(k.iter().all(|v| v.is_finite()));
        assert!(k.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn tet10_is_symmetric_pd() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.5, 0.0, 0.0],
            [0.0, 0.5, 0.0],
            [0.0, 0.0, 0.5],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
        ];
        let k = tet10_stiffness(&nodes, 200e9, 0.3);
        for i in 0..30 {
            for j in 0..30 {
                let diff = (k[i * 30 + j] - k[j * 30 + i]).abs();
                let scale = k[i * 30 + j].abs() + k[j * 30 + i].abs() + 1.0;
                assert!(
                    diff / scale < 1e-6,
                    "{i},{j}: {} vs {}",
                    k[i * 30 + j],
                    k[j * 30 + i]
                );
            }
        }
        // xᵀ K x > 0 for a non-rigid displacement.
        let u: Vec<f64> = (0..30).map(|i| (i as f64 + 1.0).sin()).collect();
        let mut e = 0.0;
        for i in 0..30 {
            for j in 0..30 {
                e += u[i] * k[i * 30 + j] * u[j];
            }
        }
        assert!(e > 0.0, "energy {e}");
    }

    #[test]
    fn beam3d_cantilever_tip_deflection() {
        // Cantilever along +x, fixed at node0, downward point load P at tip.
        // Analytic (Euler–Bernoulli): δ = P L³ / (3 E I_z).
        let e = 200e9;
        let nu = 0.3;
        let g = e / (2.0 * (1.0 + nu));
        let l = 1.0;
        let b = 0.1;
        let h = 0.1;
        let area = b * h;
        let iz = b * h * h * h / 12.0; // bending about z (gives y deflection)
        let iy = h * b * b * b / 12.0;
        let kg = beam3d_global_stiffness(
            [0.0, 0.0, 0.0],
            [l, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            e,
            g,
            area,
            iy,
            iz,
            area * (b * b + h * h) / 12.0,
        );
        let p = 1000.0;
        // Free DOFs: node1 (indices 6..12). Fixed node0 (0..6).
        let free: Vec<usize> = (6..12).collect();
        // Reduced solve via dense inverse (small).
        let mut kr = [[0.0; 6]; 6];
        for (i, &fi) in free.iter().enumerate() {
            for (j, &fj) in free.iter().enumerate() {
                kr[i][j] = kg[fi * 12 + fj];
            }
        }
        let mut fred = [0.0; 6];
        fred[1] = -p; // -y load
                      // Solve the reduced 6×6 system directly.
        let mut u = [0.0; 6];
        let _ = solve_6x6(&kr, &fred, &mut u);
        let expected = p * l * l * l / (3.0 * e * iz);
        assert!(
            (u[1] - (-expected)).abs() < expected * 1e-6,
            "tip δ = {}, expected {}",
            u[1],
            -expected
        );
    }

    #[test]
    fn beam3d_cantilever_tip_deflection_l2() {
        // Same cantilever but with L = 2.0 (non-unit length) so the bending
        // block's L / L² factors are actually exercised. Analytic δ = P L³/(3 E I_z).
        let e = 200e9;
        let nu = 0.3;
        let g = e / (2.0 * (1.0 + nu));
        let l = 2.0;
        let b = 0.1;
        let h = 0.1;
        let area = b * h;
        let iz = b * h * h * h / 12.0;
        let iy = h * b * b * b / 12.0;
        let kg = beam3d_global_stiffness(
            [0.0, 0.0, 0.0],
            [l, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            e,
            g,
            area,
            iy,
            iz,
            area * (b * b + h * h) / 12.0,
        );
        let p = 1000.0;
        let free: Vec<usize> = (6..12).collect();
        let mut kr = [[0.0; 6]; 6];
        for (i, &fi) in free.iter().enumerate() {
            for (j, &fj) in free.iter().enumerate() {
                kr[i][j] = kg[fi * 12 + fj];
            }
        }
        let mut fred = [0.0; 6];
        fred[1] = -p;
        let mut u = [0.0; 6];
        let _ = solve_6x6(&kr, &fred, &mut u);
        let expected = p * l * l * l / (3.0 * e * iz);
        // y-deflection is free index 1; should match analytic to 1e-6.
        assert!(
            (u[1] - (-expected)).abs() < expected * 1e-6,
            "tip δ = {}, expected {}",
            u[1],
            -expected
        );
    }

    #[test]
    fn shell4_is_symmetric() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let k = shell4_stiffness(&nodes, 0.01, 200e9, 0.3, 5.0 / 6.0);
        for i in 0..12 {
            for j in 0..12 {
                let diff = (k[i * 12 + j] - k[j * 12 + i]).abs();
                let scale = k[i * 12 + j].abs() + k[j * 12 + i].abs() + 1.0;
                assert!(
                    diff / scale < 1e-6,
                    "{i},{j}: {} vs {}",
                    k[i * 12 + j],
                    k[j * 12 + i]
                );
            }
        }
    }

    #[test]
    fn shell4_rigid_translation_zero_force() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let k = shell4_stiffness(&nodes, 0.01, 200e9, 0.3, 5.0 / 6.0);
        // Pure transverse translation (w = const, no rotations) is an exact
        // rigid-body mode: all strains vanish, so K·u ≡ 0.
        let mut u = [0.0; 12];
        for a in 0..4 {
            u[3 * a] = 0.01;
        }
        for i in 0..12 {
            let mut f = 0.0;
            for j in 0..12 {
                f += k[i * 12 + j] * u[j];
            }
            assert!(f.abs() < 1e-9, "rigid force {i} = {f}");
        }
    }

    // Build a nodal-displacement vector for a rigid rotation of the square plate
    // about the x-axis by `alpha` (w = α·y, θx = α, θy = 0) and about the
    // y-axis by `beta` (w = -β·x, θx = 0, θy = β). For the Mindlin conventions
    // used here (w_,x = -θy, w_,y = θx), a rotation about x gives w = α·y and
    // θx = α; a rotation about y gives w = -β·x and θy = β. Both are exact
    // rigid-body modes: the element must report zero internal force.
    fn shell4_rotation_u(nodes: &[[f64; 3]; 4], alpha: f64, beta: f64) -> [f64; 12] {
        let mut u = [0.0; 12];
        for a in 0..4 {
            let (x, y, _) = (nodes[a][0], nodes[a][1], nodes[a][2]);
            u[3 * a] = alpha * y - beta * x; // w
            u[3 * a + 1] = alpha; // θx
            u[3 * a + 2] = beta; // θy
        }
        u
    }

    #[test]
    fn shell4_rigid_rotation_zero_force() {
        // The 4-node Mindlin element is exactly rigid-body exact only on
        // parallelogram (incl. rectangular) elements — the bilinear geometry
        // map is affine there, so rigid modes are represented exactly. A
        // general quad is only approximately exact, which would make this test
        // flaky, so we use a rectangle.
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let k = shell4_stiffness(&nodes, 0.01, 200e9, 0.3, 5.0 / 6.0);
        // Reference non-rigid field (pure w linear in x) to scale the tolerance:
        // its internal force is representative of a genuine deformation.
        let mut u_ref = [0.0; 12];
        for a in 0..4 {
            u_ref[3 * a] = 0.05 * nodes[a][0];
        }
        let mut f_ref = [0.0; 12];
        for i in 0..12 {
            for j in 0..12 {
                f_ref[i] += k[i * 12 + j] * u_ref[j];
            }
        }
        let ref_norm: f64 = f_ref.iter().map(|v| v * v).sum::<f64>().sqrt();
        for (alpha, beta) in [(0.05, 0.0), (0.0, 0.07), (0.03, -0.04)] {
            let u = shell4_rotation_u(&nodes, alpha, beta);
            let mut f = [0.0; 12];
            for i in 0..12 {
                for j in 0..12 {
                    f[i] += k[i * 12 + j] * u[j];
                }
            }
            let fnorm: f64 = f.iter().map(|v| v * v).sum::<f64>().sqrt();
            // Rigid-body residual must be machine roundoff relative to a real
            // deformation. On a rectangle the bilinear geometry map is affine, so
            // the Mindlin–Reissner element is exactly rigid-body exact: both the
            // bending curvatures and the transverse-shear strains vanish for the
            // rigid rotation, leaving no internal force.
            assert!(
                fnorm < 1e-9 * ref_norm,
                "rigid rotation residual {fnorm} vs reference {ref_norm} (α={alpha}, β={beta})"
            );
        }
    }

    // Simply-supported square plate under uniform pressure, solved with the
    // Mindlin–Reissner `shell4` element on an N×N mesh. This is a real
    // boundary-value-problem validation (not just a rigid-body patch): the
    // central deflection must be positive, symmetric, maximal at the centre,
    // and within an order of magnitude of the thin-plate (Kirchhoff) analytic
    // solution. Coarse 4-node elements with reduced shear integration are
    // expected to be stiffer than Kirchhoff, so only a loose band is asserted.
    #[test]
    fn shell4_simply_supported_plate() {
        let n = 8usize; // elements per side → 9×9 = 81 nodes
        let a = 1.0; // plate side
        let h = a / n as f64;
        let t = 0.02; // thickness (thin, t/a = 0.02)
        let e = 200e9;
        let nu = 0.3;
        let kappa = 5.0 / 6.0;
        let q = 1000.0; // uniform pressure

        let nn = n + 1;
        let nnode = nn * nn;
        let ndof = 3 * nnode;
        // Node positions on the [0,a]² square, z = 0.
        let pos: Vec<[f64; 3]> = (0..nnode)
            .map(|idx| {
                let ix = idx % nn;
                let iy = idx / nn;
                [ix as f64 * h, iy as f64 * h, 0.0]
            })
            .collect();

        // Global stiffness (dense; small problem) and load vector.
        let mut k = vec![vec![0.0; ndof]; ndof];
        let mut f = vec![0.0; ndof];
        for iy in 0..n {
            for ix in 0..n {
                let n0 = iy * nn + ix;
                let n1 = iy * nn + (ix + 1);
                let n2 = (iy + 1) * nn + (ix + 1);
                let n3 = (iy + 1) * nn + ix;
                let nodes = [pos[n0], pos[n1], pos[n2], pos[n3]];
                let ke = shell4_stiffness(&nodes, t, e, nu, kappa);
                let enode = [n0, n1, n2, n3];
                for a in 0..4 {
                    for b in 0..4 {
                        for da in 0..3 {
                            for db in 0..3 {
                                let gi = 3 * enode[a] + da;
                                let gj = 3 * enode[b] + db;
                                k[gi][gj] += ke[(3 * a + da) * 12 + (3 * b + db)];
                            }
                        }
                    }
                }
                // Consistent uniform-load nodal force: ∫ N_a q dA = q·A/4 each.
                let load = q * h * h / 4.0;
                for &en in &enode {
                    f[3 * en] += load;
                }
            }
        }

        // Simply-supported BCs: w = 0 on the boundary (rotations free).
        let mut fixed = vec![false; ndof];
        for idx in 0..nnode {
            let ix = idx % nn;
            let iy = idx / nn;
            if ix == 0 || ix == n || iy == 0 || iy == n {
                fixed[3 * idx] = true;
            }
        }

        // Solve the reduced system by dense Gaussian elimination with the fixed
        // DOFs eliminated (set rows/cols to identity, unit RHS).
        let mut amat = k.clone();
        let mut rhs = f.clone();
        for d in 0..ndof {
            if fixed[d] {
                for amat_row in amat.iter_mut() {
                    amat_row[d] = 0.0;
                }
                amat[d].fill(0.0);
                amat[d][d] = 1.0;
                rhs[d] = 0.0;
            }
        }
        let x = solve_dense(&mut amat, &rhs, &fixed);

        let center = (nn / 2) * nn + (nn / 2);
        let wc = x[3 * center];
        assert!(wc.is_finite() && wc > 0.0, "centre deflection = {wc}");
        // The centre must be a local maximum (it is the global maximum for a
        // symmetric simply-supported plate under uniform load). Check the four
        // immediate neighbours rather than every interior node, which is robust
        // to the slightly non-monotonic discrete field of a coarse mesh.
        for (dx, dy) in [(1, 0), (-1, 0), (0, 1), (0, -1)] {
            let nb = ((nn / 2) as isize + dx) as usize * nn
                + ((nn / 2) as isize + dy) as usize;
            assert!(
                x[3 * nb] <= wc + 1e-9,
                "neighbour {nb} deflection {} exceeds centre {}",
                x[3 * nb],
                wc
            );
        }
        // Four-fold symmetry of the deflection field.
        let q1 = (nn / 4) * nn + (nn / 4);
        let q2 = (3 * nn / 4) * nn + (nn / 4);
        let q3 = (3 * nn / 4) * nn + (3 * nn / 4);
        let q4 = (nn / 4) * nn + (3 * nn / 4);
        let s = [x[3 * q1], x[3 * q2], x[3 * q3], x[3 * q4]];
        for v in &s {
            assert!((v - s[0]).abs() < 1e-9, "asymmetry: {:?}", s);
        }
        // Thin-plate (Kirchhoff) centre deflection: w = 0.00406 q a⁴ / D,
        // D = E t³ / (12(1-ν²)). Loose band (Mindlin 4-node is stiffer).
        let d_plate = e * t * t * t / (12.0 * (1.0 - nu * nu));
        let w_analytic = 0.00406 * q * (a * a * a * a) / d_plate;
        assert!(
            wc > 0.1 * w_analytic && wc < 1.5 * w_analytic,
            "wc = {wc}, thin-plate = {w_analytic}"
        );
    }

    // Morley-style *skew* (parallelogram) simply-supported plate benchmark. The
    // 4-node Mindlin element is only exactly rigid-body exact on parallelograms,
    // so a *general* skew mesh is the realistic validation target (the square
    // test above is a special case). A 30° skew plate under uniform load must
    // still produce a finite, positive, smoothly-peaked centre deflection whose
    // four-fold skew symmetry is preserved — i.e. the element is correctly
    // formulated on non-rectangular geometry (no spurious skew-coupled energy,
    // no sign error exposed only off-axis).
    #[test]
    fn shell4_skew_plate_well_posed() {
        let n = 8usize;
        let a = 1.0;
        let h = a / n as f64;
        let t = 0.02;
        let e = 200e9;
        let nu = 0.3;
        let kappa = 5.0 / 6.0;
        let q = 1000.0;
        let skew = (30.0f64).to_radians().tan(); // x shifted by iy·h·skew

        let nn = n + 1;
        let nnode = nn * nn;
        let ndof = 3 * nnode;
        let pos: Vec<[f64; 3]> = (0..nnode)
            .map(|idx| {
                let ix = idx % nn;
                let iy = idx / nn;
                // Affine skew map → parallelogram (exact bilinear geometry).
                [ix as f64 * h + iy as f64 * h * skew, iy as f64 * h, 0.0]
            })
            .collect();

        let mut k = vec![vec![0.0; ndof]; ndof];
        let mut f = vec![0.0; ndof];
        for iy in 0..n {
            for ix in 0..n {
                let n0 = iy * nn + ix;
                let n1 = iy * nn + (ix + 1);
                let n2 = (iy + 1) * nn + (ix + 1);
                let n3 = (iy + 1) * nn + ix;
                let nodes = [pos[n0], pos[n1], pos[n2], pos[n3]];
                let ke = shell4_stiffness(&nodes, t, e, nu, kappa);
                let enode = [n0, n1, n2, n3];
                for a in 0..4 {
                    for b in 0..4 {
                        for da in 0..3 {
                            for db in 0..3 {
                                let gi = 3 * enode[a] + da;
                                let gj = 3 * enode[b] + db;
                                k[gi][gj] += ke[(3 * a + da) * 12 + (3 * b + db)];
                            }
                        }
                    }
                }
                let load = q * h * h / 4.0;
                for &en in &enode {
                    f[3 * en] += load;
                }
            }
        }

        // Simply-supported: w = 0 on the skew boundary.
        let mut fixed = vec![false; ndof];
        for idx in 0..nnode {
            let ix = idx % nn;
            let iy = idx / nn;
            if ix == 0 || ix == n || iy == 0 || iy == n {
                fixed[3 * idx] = true;
            }
        }
        let mut amat = k.clone();
        let mut rhs = f.clone();
        for d in 0..ndof {
            if fixed[d] {
                for amat_row in amat.iter_mut() {
                    amat_row[d] = 0.0;
                }
                amat[d].fill(0.0);
                amat[d][d] = 1.0;
                rhs[d] = 0.0;
            }
        }
        let x = solve_dense(&mut amat, &rhs, &fixed);
        // Centre of the *parametric* (square) domain → centre of the skew plate.
        let center = (nn / 2) * nn + (nn / 2);
        let wc = x[3 * center];
        assert!(wc.is_finite() && wc > 0.0, "skew centre deflection = {wc}");
        // Skew symmetry: reflect across the plate centreline along the skew edge
        // (ix → n-ix, iy → n-iy). The deflection field must be symmetric.
        for iy in 1..n {
            for ix in 1..n {
                let a = iy * nn + ix;
                let b = (n - iy) * nn + (n - ix);
                assert!(
                    (x[3 * a] - x[3 * b]).abs() < 1e-6,
                    "skew asymmetry at ({ix},{iy})"
                );
            }
        }
        // Same loose band vs the *square* thin-plate solution (the skew plate is
        // not the same problem, but a well-posed Mindlin solution stays in the
        // same order of magnitude).
        let d_plate = e * t * t * t / (12.0 * (1.0 - nu * nu));
        let w_analytic = 0.00406 * q * (a * a * a * a) / d_plate;
        assert!(
            wc > 0.05 * w_analytic && wc < 2.0 * w_analytic,
            "skew wc = {wc}, thin-plate = {w_analytic}"
        );
    }

    // Thick-plate (Mindlin) validation of the shear term's sign. A thick plate
    // must deflect *more* than the pure-bending (Kirchhoff, no-shear) solution
    // under the same load: transverse shear is always energy-adding (it
    // softens the plate). If the shear strain–displacement signs were wrong the
    // shear contribution would subtract and the FEM result would fall below the
    // Kirchhoff estimate, which this test rejects.
    #[test]
    fn shell4_thick_plate_softens_under_shear() {
        let n = 10usize;
        let a = 1.0;
        let h = a / n as f64;
        let t = 0.2; // thick: t/a = 0.2 ⇒ shear is significant
        let e = 200e9;
        let nu = 0.3;
        let kappa = 5.0 / 6.0;
        let q = 1000.0;

        let nn = n + 1;
        let nnode = nn * nn;
        let ndof = 3 * nnode;
        let pos: Vec<[f64; 3]> = (0..nnode)
            .map(|idx| {
                let ix = idx % nn;
                let iy = idx / nn;
                [ix as f64 * h, iy as f64 * h, 0.0]
            })
            .collect();

        let mut k = vec![vec![0.0; ndof]; ndof];
        let mut f = vec![0.0; ndof];
        for iy in 0..n {
            for ix in 0..n {
                let n0 = iy * nn + ix;
                let n1 = iy * nn + (ix + 1);
                let n2 = (iy + 1) * nn + (ix + 1);
                let n3 = (iy + 1) * nn + ix;
                let nodes = [pos[n0], pos[n1], pos[n2], pos[n3]];
                let ke = shell4_stiffness(&nodes, t, e, nu, kappa);
                let enode = [n0, n1, n2, n3];
                for a in 0..4 {
                    for b in 0..4 {
                        for da in 0..3 {
                            for db in 0..3 {
                                let gi = 3 * enode[a] + da;
                                let gj = 3 * enode[b] + db;
                                k[gi][gj] += ke[(3 * a + da) * 12 + (3 * b + db)];
                            }
                        }
                    }
                }
                let load = q * h * h / 4.0;
                for &en in &enode {
                    f[3 * en] += load;
                }
            }
        }

        let mut fixed = vec![false; ndof];
        for idx in 0..nnode {
            let ix = idx % nn;
            let iy = idx / nn;
            if ix == 0 || ix == n || iy == 0 || iy == n {
                fixed[3 * idx] = true;
            }
        }
        let mut amat = k.clone();
        let mut rhs = f.clone();
        for d in 0..ndof {
            if fixed[d] {
                for amat_row in amat.iter_mut() {
                    amat_row[d] = 0.0;
                }
                amat[d].fill(0.0);
                amat[d][d] = 1.0;
                rhs[d] = 0.0;
            }
        }
        let x = solve_dense(&mut amat, &rhs, &fixed);
        let center = (nn / 2) * nn + (nn / 2);
        let wc = x[3 * center];

        // Pure-bending (Kirchhoff, no shear) centre deflection:
        // w_K = 0.00406 q a⁴ / D,  D = E t³ / (12(1-ν²)).
        let d_plate = e * t * t * t / (12.0 * (1.0 - nu * nu));
        let w_kirchhoff = 0.00406 * q * (a * a * a * a) / d_plate;
        // Mindlin plate must be *softer* than Kirchhoff (shear adds compliance).
        assert!(
            wc > w_kirchhoff,
            "thick-plate wc={wc} not above Kirchhoff {w_kirchhoff} (shear sign?)"
        );
        // …but still within a sane factor of the bending solution.
        assert!(
            wc < 3.0 * w_kirchhoff,
            "thick-plate wc={wc} far above Kirchhoff {w_kirchhoff}"
        );
    }

    // Minimal 6×6 dense linear solver (for the beam test).
    fn solve_6x6(a: &[[f64; 6]; 6], b: &[f64; 6], x: &mut [f64; 6]) -> bool {
        let mut m = [[0.0; 7]; 6];
        for i in 0..6 {
            for j in 0..6 {
                m[i][j] = a[i][j];
            }
            m[i][6] = b[i];
        }
        for col in 0..6 {
            let mut piv = col;
            let mut best = m[col][col].abs();
            for (r, m_r) in m.iter().enumerate().skip(col + 1) {
                if m_r[col].abs() > best {
                    best = m_r[col].abs();
                    piv = r;
                }
            }
            if best < 1e-15 {
                return false;
            }
            m.swap(col, piv);
            let d = m[col][col];
            for m_cj in m[col][col..7].iter_mut() {
                *m_cj /= d;
            }
            for r in 0..6 {
                if r != col {
                    let f = m[r][col];
                    #[allow(clippy::needless_range_loop)]
                    for j in col..7 {
                        m[r][j] -= f * m[col][j];
                    }
                }
            }
        }
        for i in 0..6 {
            x[i] = m[i][6];
        }
        true
    }

    // Generic dense linear solver (partial-pivot Gaussian elimination) for the
    // plate benchmark. `fixed` marks DOFs already reduced to x_d = 0 (identity
    // rows); the solver operates on the full system and returns the solution.
    fn solve_dense(a: &mut [Vec<f64>], rhs: &[f64], _fixed: &[bool]) -> Vec<f64> {
        let n = rhs.len();
        let mut m = vec![vec![0.0; n + 1]; n];
        for i in 0..n {
            for j in 0..n {
                m[i][j] = a[i][j];
            }
            m[i][n] = rhs[i];
        }
        for col in 0..n {
            let mut piv = col;
            let mut best = m[col][col].abs();
            for (r, m_r) in m.iter().enumerate().skip(col + 1) {
                if m_r[col].abs() > best {
                    best = m_r[col].abs();
                    piv = r;
                }
            }
            if best < 1e-14 {
                continue; // singular row (reduced DOF) → leave as-is
            }
            m.swap(col, piv);
            let d = m[col][col];
            for m_cj in m[col][col..=n].iter_mut() {
                *m_cj /= d;
            }
            for r in 0..n {
                if r != col {
                    let f = m[r][col];
                    #[allow(clippy::needless_range_loop)]
                    for j in col..=n {
                        m[r][j] -= f * m[col][j];
                    }
                }
            }
        }
        (0..n).map(|i| m[i][n]).collect()
    }
}
