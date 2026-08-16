//! Advanced finite-element formulations that are net-new to `tpt-physics`:
//! the quadratic tetrahedral element ([`tet10_stiffness`]), the 3-D beam/frame
//! element ([`beam3d_global_stiffness`]), and a Mindlin–Reissner plate/shell
//! element ([`shell4_stiffness`]).
//!
//! Linear-tetrahedron, hexahedron, and 2-D frame elements are reused directly
//! from `tpt-fem-element` / `tpt-fem-elasticity`; see the crate root re-exports.

/// 3×3 matrix inverse (row-major `m`, length 9). Panics on singular input.
fn mat3_inv(m: &[f64; 9]) -> [f64; 9] {
    let det = m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6])
        + m[2] * (m[3] * m[7] - m[4] * m[6]);
    let inv = 1.0 / det;
    [
        (m[4] * m[8] - m[5] * m[7]) * inv,
        (m[2] * m[7] - m[1] * m[8]) * inv,
        (m[1] * m[5] - m[2] * m[4]) * inv,
        (m[5] * m[6] - m[3] * m[8]) * inv,
        (m[0] * m[8] - m[2] * m[6]) * inv,
        (m[2] * m[3] - m[0] * m[5]) * inv,
        (m[3] * m[7] - m[4] * m[6]) * inv,
        (m[1] * m[6] - m[0] * m[7]) * inv,
        (m[0] * m[4] - m[1] * m[3]) * inv,
    ]
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
        let inv = mat3_inv(&j);

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
        for a in 0..10 {
            let i = 3 * a;
            b[0][i] = gp[a][0];
            b[1][i + 1] = gp[a][1];
            b[2][i + 2] = gp[a][2];
            b[3][i] = gp[a][1];
            b[3][i + 1] = gp[a][0];
            b[4][i + 1] = gp[a][2];
            b[4][i + 2] = gp[a][1];
            b[5][i] = gp[a][2];
            b[5][i + 2] = gp[a][0];
        }

        // K += Bᵀ D B |J| w.
        let scale = detj * w;
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
        r[p * 3 + 0] = ex[p];
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
    let _ = add_beam_bending(&mut kl, &[1, 5, 7, 11], e * iz, l); // v-θz (about z)
    let _ = add_beam_bending(&mut kl, &[2, 4, 8, 10], e * iy, l); // w-θy (about y)

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
/// factor (≈ 5/6). Bending uses a 2×2 Gauss rule; transverse shear uses a
/// single reduced-integration point to suppress shear locking.
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
        dn[0] = [-0.25 * (1.0 - s[1]), -0.25 * (1.0 - s[1])];
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

    // Transverse shear: reduced 1-point integration (xi=eta=0).
    {
        let (n, dn) = shape(0.0, 0.0);
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
        for a in 0..4 {
            gx[a] = inv[0][0] * dn[a][0] + inv[0][1] * dn[a][1];
            gy[a] = inv[1][0] * dn[a][0] + inv[1][1] * dn[a][1];
        }
        // Shear B (2 × 12): γ_xz = ∂w/∂x - θy(x) ; γ_yz = ∂w/∂y + θx(x),
        // where θy(x) = Σ N_a θy_a (so the rotation coefficient is the shape N_a).
        let mut bs = [[0.0_f64; 12]; 2];
        for a in 0..4 {
            let i = 3 * a;
            bs[0][i] = gx[a]; // ∂w/∂x
            bs[0][i + 2] = -n[a]; // -θy
            bs[1][i] = gy[a]; // ∂w/∂y
            bs[1][i + 1] = n[a]; // +θx
        }
        // K += Bsᵀ Gs Bs |J|, Gs = gs * I.
        for i in 0..12 {
            for j in 0..12 {
                let mut s = 0.0;
                for s1 in 0..2 {
                    s += bs[s1][i] * gs * bs[s1][j];
                }
                k[i * 12 + j] += s * detj;
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
            for r in (col + 1)..6 {
                if m[r][col].abs() > best {
                    best = m[r][col].abs();
                    piv = r;
                }
            }
            if best < 1e-15 {
                return false;
            }
            m.swap(col, piv);
            let d = m[col][col];
            for j in col..7 {
                m[col][j] /= d;
            }
            for r in 0..6 {
                if r != col {
                    let f = m[r][col];
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
}
