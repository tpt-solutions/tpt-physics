//! Geometric-nonlinear continuum element (Total Lagrangian, St. Venant–
//! Kirchhoff). This is the net-new large-deformation framework: `tpt-fem-solve`
//! provides Newton–Raphson continuation, but only against a hand-written truss
//! residual — a general continuum tangent across real elements does not exist
//! there. Here we provide the consistent internal-force and tangent-stiffness
//! operators for the linear tetrahedral element (constant deformation gradient,
//! exact for Tet4). The same Total-Lagrangian pattern extends to Hex8/beam by
//! supplying their reference gradients and quadrature.

/// Reference-coordinate gradients `∂N_a/∂ξ` of the 4-node linear tetrahedron
/// (constant), indexed `[node][ξ-component]`.
const TET4_GRAD_XI: [[f64; 3]; 4] = [
    [-1.0, -1.0, -1.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];

/// 3×3 matrix multiply `a · b` (row-major `[[f64;3];3]`).
fn mm(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut s = 0.0;
            for k in 0..3 {
                s += a[i][k] * b[k][j];
            }
            c[i][j] = s;
        }
    }
    c
}

/// 3×3 inverse; panics on singular.
fn inv3(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    let inv = 1.0 / det;
    [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv,
        ],
    ]
}

/// Deformation gradient `F = ∂x/∂X` for a linear tetrahedral element.
pub fn deformation_gradient(ref_coords: &[[f64; 3]; 4], cur_coords: &[[f64; 3]; 4]) -> [[f64; 3]; 3] {
    let mut jref = [[0.0; 3]; 3];
    let mut jcur = [[0.0; 3]; 3];
    for a in 0..4 {
        for p in 0..3 {
            for r in 0..3 {
                jref[p][r] += ref_coords[a][p] * TET4_GRAD_XI[a][r];
                jcur[p][r] += cur_coords[a][p] * TET4_GRAD_XI[a][r];
            }
        }
    }
    mm(&jcur, &inv3(&jref))
}

/// Reference configuration: material gradients `∂N_a/∂X` (length-3 per node)
/// and the reference volume `V0`.
pub(crate) fn reference_basis(ref_coords: &[[f64; 3]; 4]) -> ([[f64; 3]; 4], f64) {
    let mut jref = [[0.0; 3]; 3];
    for a in 0..4 {
        for p in 0..3 {
            for r in 0..3 {
                jref[p][r] += ref_coords[a][p] * TET4_GRAD_XI[a][r];
            }
        }
    }
    let det = jref[0][0] * (jref[1][1] * jref[2][2] - jref[1][2] * jref[2][1])
        - jref[0][1] * (jref[1][0] * jref[2][2] - jref[1][2] * jref[2][0])
        + jref[0][2] * (jref[1][0] * jref[2][1] - jref[1][1] * jref[2][0]);
    let v0 = det.abs() / 6.0;
    let inv = inv3(&jref);
    let mut gx = [[0.0; 3]; 4];
    for a in 0..4 {
        for r in 0..3 {
            let mut s = 0.0;
            for s2 in 0..3 {
                s += TET4_GRAD_XI[a][s2] * inv[s2][r];
            }
            gx[a][r] = s;
        }
    }
    (gx, v0)
}

/// Green–Lagrange strain (6-vector, engineering) from `F`.
fn green_lagrange(f: &[[f64; 3]; 3]) -> [f64; 6] {
    // E = 0.5 (FᵀF - I).
    let ftf = mm(&transpose(f), f);
    [
        0.5 * (ftf[0][0] - 1.0),
        0.5 * (ftf[1][1] - 1.0),
        0.5 * (ftf[2][2] - 1.0),
        ftf[0][1], // 2E_xy
        ftf[1][2], // 2E_yz
        ftf[0][2], // 2E_xz
    ]
}

fn transpose(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut t = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = m[j][i];
        }
    }
    t
}

/// 2nd Piola–Kirchhoff stress (6-vector) for St. Venant–Kirchhoff material.
fn second_pk_stress(e: &[f64; 6], lambda: f64, mu: f64) -> [f64; 6] {
    let tr = e[0] + e[1] + e[2];
    [
        lambda * tr + 2.0 * mu * e[0],
        lambda * tr + 2.0 * mu * e[1],
        lambda * tr + 2.0 * mu * e[2],
        mu * e[3],
        mu * e[4],
        mu * e[5],
    ]
}

/// Iso elastic tangent `C` (6×6, Voigt) for St. Venant–Kirchhoff.
fn iso_tangent(lambda: f64, mu: f64) -> [f64; 36] {
    let mut c = [0.0; 36];
    let k = lambda + 2.0 * mu;
    c[0] = k; c[7] = k; c[14] = k;
    c[21] = mu; c[28] = mu; c[35] = mu;
    c[1] = lambda; c[2] = lambda;
    c[6] = lambda; c[8] = lambda;
    c[12] = lambda; c[13] = lambda;
    c
}

/// Internal nodal force vector (length 12) of a linear tetrahedral element
/// under the given current configuration, using a St. Venant–Kirchhoff
/// material (`lambda`, `mu`). Returns `f = ∫ Bᵀ S dV` (the elastic restoring
/// force; subtract from the external load in a residual).
pub fn tet4_internal_force(
    ref_coords: &[[f64; 3]; 4],
    cur_coords: &[[f64; 3]; 4],
    lambda: f64,
    mu: f64,
) -> [f64; 12] {
    let (gx, v0) = reference_basis(ref_coords);
    let f = deformation_gradient(ref_coords, cur_coords);
    let e = green_lagrange(&f);
    let s = second_pk_stress(&e, lambda, mu);
    // Stress tensor (symmetric) from Voigt.
    let sh = [
        [s[0], s[3], s[5]],
        [s[3], s[1], s[4]],
        [s[5], s[4], s[2]],
    ];
    let mut force = [0.0; 12];
    for a in 0..4 {
        for i in 0..3 {
            let mut ssum = 0.0;
            for r in 0..3 {
                ssum += gx[a][r] * sh[r][i];
            }
            force[3 * a + i] = v0 * ssum;
        }
    }
    force
}

/// Consistent tangent stiffness (length 144, row-major 12×12) of a linear
/// tetrahedral element: material part `Bᵀ C B` plus the geometric (initial
/// stress) part built from the 2nd P-K stress.
pub fn tet4_tangent(
    ref_coords: &[[f64; 3]; 4],
    cur_coords: &[[f64; 3]; 4],
    lambda: f64,
    mu: f64,
) -> Vec<f64> {
    let (gx, v0) = reference_basis(ref_coords);
    let f = deformation_gradient(ref_coords, cur_coords);
    let e = green_lagrange(&f);
    let s = second_pk_stress(&e, lambda, mu);
    let c = iso_tangent(lambda, mu);
    let sh = [
        [s[0], s[3], s[5]],
        [s[3], s[1], s[4]],
        [s[5], s[4], s[2]],
    ];

    let b_node = |a: usize, strain: usize, comp: usize| -> f64 {
        // B mapping node `a` displacement component `comp` → strain `strain`.
        // Ordering: εxx,εyy,εzz,γxy,γyz,γxz; comp: 0=x,1=y,2=z.
        match (strain, comp) {
            (0, 0) => gx[a][0],
            (1, 1) => gx[a][1],
            (2, 2) => gx[a][2],
            (3, 0) => gx[a][1],
            (3, 1) => gx[a][0],
            (4, 1) => gx[a][2],
            (4, 2) => gx[a][1],
            (5, 0) => gx[a][2],
            (5, 2) => gx[a][0],
            _ => 0.0,
        }
    };

    let mut k = vec![0.0; 144];
    for a in 0..4 {
        for b in 0..4 {
            for i in 0..3 {
                for j in 0..3 {
                    // Material part.
                    let mut mat = 0.0;
                    for s1 in 0..6 {
                        let mut t = 0.0;
                        for s2 in 0..6 {
                            t += c[s1 * 6 + s2] * b_node(b, s2, j);
                        }
                        mat += b_node(a, s1, i) * t;
                    }
                    // Geometric (initial-stress) part: g_aᵀ S g_b.
                    let mut geo = 0.0;
                    for r in 0..3 {
                        for l in 0..3 {
                            geo += gx[a][r] * sh[r][l] * gx[b][l];
                        }
                    }
                    k[(3 * a + i) * 12 + (3 * b + j)] += v0 * (mat + geo);
                }
            }
        }
    }
    k
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lame(e: f64, nu: f64) -> (f64, f64) {
        let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let mu = e / (2.0 * (1.0 + nu));
        (lambda, mu)
    }

    #[test]
    fn zero_deformation_zero_force() {
        let r = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let (la, mu) = lame(200e9, 0.3);
        let f = tet4_internal_force(&r, &r, la, mu);
        assert!(f.iter().all(|v| v.abs() < 1e-9));
    }

    #[test]
    fn internal_force_is_self_equilibrated() {
        // Any (non-rigid) deformation yields a nodal force sum of zero.
        let r = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let c = [[0.1, 0.2, 0.0], [1.1, -0.1, 0.3], [-0.2, 0.9, 0.1], [0.05, 0.05, 1.1]];
        let (la, mu) = lame(200e9, 0.3);
        let f = tet4_internal_force(&r, &c, la, mu);
        let mut s = [0.0; 3];
        for a in 0..4 {
            for i in 0..3 {
                s[i] += f[3 * a + i];
            }
        }
        assert!(s.iter().all(|v| v.abs() < 1e-3), "sum = {:?}", s);
    }

    #[test]
    fn rigorous_rotation_gives_zero_force() {
        // A pure rotation has F orthonormal → E = 0 → S = 0 → f = 0.
        let r = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let ang = 0.3_f64;
        let rot = [
            [1.0, 0.0, 0.0],
            [0.0, ang.cos(), -ang.sin()],
            [0.0, ang.sin(), ang.cos()],
        ];
        let c: Vec<[f64; 3]> = r.iter().map(|p| [rot[0][0]*p[0]+rot[0][1]*p[1]+rot[0][2]*p[2], rot[1][0]*p[0]+rot[1][1]*p[1]+rot[1][2]*p[2], rot[2][0]*p[0]+rot[2][1]*p[1]+rot[2][2]*p[2]]).collect();
        let (la, mu) = lame(200e9, 0.3);
        let f = tet4_internal_force(&r, &c.try_into().unwrap(), la, mu);
        assert!(f.iter().all(|v| v.abs() < 1e-2), "rotation force {:?}", f);
    }

    #[test]
    fn tangent_linearizes_internal_force() {
        // At the undeformed state the tangent equals the linear elastic
        // stiffness, so f(cur) ≈ K(cur-ref) for a small displacement.
        let r = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let (la, mu) = lame(200e9, 0.3);
        let k = tet4_tangent(&r, &r, la, mu);
        let mut u = [0.0; 12];
        for i in 0..12 {
            u[i] = (i as f64 * 0.013 - 0.05) * 1e-3; // ~1e-4 strain scale
        }
        let mut cur = r;
        for a in 0..4 {
            for i in 0..3 {
                cur[a][i] = r[a][i] + u[3 * a + i];
            }
        }
        let f = tet4_internal_force(&r, &cur, la, mu);
        for i in 0..12 {
            let mut ku = 0.0;
            for j in 0..12 {
                ku += k[i * 12 + j] * u[j];
            }
            let rel = (ku - f[i]).abs() / (f[i].abs().max(1.0));
            assert!(rel < 1e-2, "i={i} K·u={ku} f={} rel={rel}", f[i]);
        }
    }

    #[test]
    fn tangent_is_symmetric() {
        let r = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let c = [[0.1, 0.05, 0.0], [1.05, -0.02, 0.1], [-0.05, 0.95, 0.05], [0.02, 0.03, 1.05]];
        let (la, mu) = lame(200e9, 0.3);
        let k = tet4_tangent(&r, &c, la, mu);
        for i in 0..12 {
            for j in 0..12 {
                assert!((k[i * 12 + j] - k[j * 12 + i]).abs() < 1e-6);
            }
        }
    }
}
