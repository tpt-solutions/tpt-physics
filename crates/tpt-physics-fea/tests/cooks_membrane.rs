//! FEA validation: elastic–plastic Cook's membrane (J2 plasticity).
//!
//! A skewed cantilever ("Cook's membrane") is the standard benchmark for
//! nonlinear solid-mechanics codes. We solve it as a **plane-strain** problem
//! with the project's own [`PlasticMaterial`] (von Mises / J2 with linear
//! isotropic hardening) integrated by incremental return mapping at every
//! Gauss point of a bilinear (Q4) mesh. The validation checks the qualitative
//! behaviour any correct elastic–plastic solver must show, independent of the
//! 3D-printed spacer:
//!
//! 1. **Plasticity increases compliance** — under the same load the elastic
//!    solution (yield stress set artificially high) is *less* compliant than
//!    the elastic–plastic solution (finite yield stress).
//! 2. **Hardening stiffens** — a higher hardening modulus gives a *smaller*
//!    tip displacement than a lower one.
//!
//! This exercises [`tpt_physics_fea::plasticity`] through a real finite-element
//! assembly rather than the point-wise unit tests in that module.

use tpt_physics_fea::plasticity::{PlasticMaterial, PlasticState};

const N: usize = 6; // elements per side
const E: f64 = 1.0;
const NU: f64 = 0.3;
const P_TOTAL: f64 = 0.5; // total downward load on the right edge

/// Bilinear map of unit-square `(ξ, η)` to the skewed Cook's-membrane quad
/// with corners `(0,0)`, `(48,44)`, `(48,60)`, `(0,44)`.
fn node_pos(i: usize, j: usize) -> [f64; 2] {
    let xi = i as f64 / N as f64;
    let eta = j as f64 / N as f64;
    let x = 48.0 * xi;
    let y = 44.0 * xi + 44.0 * eta - 28.0 * xi * eta;
    [x, y]
}

/// Node id for grid indices `(i, j)` with `i,j ∈ [0, N]`.
fn nid(i: usize, j: usize) -> usize {
    i * (N + 1) + j
}

/// Q4 shape-function derivatives w.r.t. `(ξ, η)` at a Gauss point.
fn shape_derivs(xi: f64, eta: f64) -> [[f64; 2]; 4] {
    let xim2 = 1.0 - xi;
    let xip = 1.0 + xi;
    let etam = 1.0 - eta;
    let etap = 1.0 + eta;
    let dndxi = [-0.25 * etam, 0.25 * etam, 0.25 * etap, -0.25 * etap];
    let dndeta = [-0.25 * xim2, -0.25 * xip, 0.25 * xip, 0.25 * xim2];
    [
        [dndxi[0], dndeta[0]],
        [dndxi[1], dndeta[1]],
        [dndxi[2], dndeta[2]],
        [dndxi[3], dndeta[3]],
    ]
}

/// Assemble the 6×8 strain–displacement matrix `B` (plane strain: εzz = 0)
/// and the Jacobian determinant at a Gauss point of one element.
fn bmatrix(nodes: &[[f64; 2]; 4], xi: f64, eta: f64) -> ([[f64; 8]; 6], f64) {
    let d = shape_derivs(xi, eta);
    let mut j = [[0.0_f64; 2]; 2];
    for k in 0..4 {
        j[0][0] += d[k][0] * nodes[k][0];
        j[0][1] += d[k][0] * nodes[k][1];
        j[1][0] += d[k][1] * nodes[k][0];
        j[1][1] += d[k][1] * nodes[k][1];
    }
    let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
    let inv = [
        [j[1][1] / det, -j[0][1] / det],
        [-j[1][0] / det, j[0][0] / det],
    ];
    let mut b = [[0.0_f64; 8]; 6];
    for k in 0..4 {
        let dx = inv[0][0] * d[k][0] + inv[0][1] * d[k][1];
        let dy = inv[1][0] * d[k][0] + inv[1][1] * d[k][1];
        let c = k * 2;
        b[0][c] = dx; // εxx
        b[1][c + 1] = dy; // εyy
        b[3][c] = dy; // γxy
        b[3][c + 1] = dx;
    }
    (b, det)
}

/// Solve `A x = b` by Gaussian elimination with partial pivoting (`A` is
/// `n × n`, row-major, consumed).
fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for i in 0..n {
        let mut piv = i;
        let mut best = a[i][i].abs();
        for r in (i + 1)..n {
            if a[r][i].abs() > best {
                best = a[r][i].abs();
                piv = r;
            }
        }
        a.swap(i, piv);
        b.swap(i, piv);
        for r in (i + 1)..n {
            let f = a[r][i] / a[i][i];
            for c in i..n {
                a[r][c] -= f * a[i][c];
            }
            b[r] -= f * b[i];
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for c in (i + 1)..n {
            s -= a[i][c] * x[c];
        }
        x[i] = s / a[i][i];
    }
    x
}

/// Run the Cook's-membrane analysis with the given yield stress and hardening.
/// Returns the vertical tip displacement (node `(N, N)`).
fn cooks_membrane(sigma_y0: f64, hard: f64) -> f64 {
    let mat = PlasticMaterial::new(E, NU, sigma_y0, hard);
    let c = mat.elastic_matrix();
    let nnodes = (N + 1) * (N + 1);
    let ndof = 2 * nnodes;

    let mut pos = vec![[0.0_f64; 2]; nnodes];
    for i in 0..=N {
        for j in 0..=N {
            pos[nid(i, j)] = node_pos(i, j);
        }
    }
    let mut elems = Vec::new();
    for i in 0..N {
        for j in 0..N {
            elems.push([
                nid(i, j),
                nid(i + 1, j),
                nid(i + 1, j + 1),
                nid(i, j + 1),
            ]);
        }
    }

    // Constrained dofs: left edge (i == 0) fixed in u and v.
    let mut constrained = vec![false; ndof];
    for j in 0..=N {
        let id = nid(0, j);
        constrained[2 * id] = true;
        constrained[2 * id + 1] = true;
    }

    // Precompute Gauss-point B-matrices.
    let gp = [-1.0 / 3.0_f64.sqrt(), 1.0 / 3.0_f64.sqrt()];
    let mut gauss_elem = Vec::new();
    let mut gauss_b = Vec::new();
    let mut gauss_det = Vec::new();
    for (e, el) in elems.iter().enumerate() {
        let nodes = [pos[el[0]], pos[el[1]], pos[el[2]], pos[el[3]]];
        for &xi in &gp {
            for &eta in &gp {
                let (b, det) = bmatrix(&nodes, xi, eta);
                gauss_elem.push(e);
                gauss_b.push(b);
                gauss_det.push(det);
            }
        }
    }
    let ng = gauss_elem.len();

    // Elastic global tangent K (constant), reduced to free dofs.
    let mut k_full = vec![vec![0.0_f64; ndof]; ndof];
    for gi in 0..ng {
        let e = gauss_elem[gi];
        let bmat = &gauss_b[gi];
        let det = gauss_det[gi];
        let el = elems[e];
        let mut ldof = [0usize; 8];
        for m in 0..4 {
            ldof[2 * m] = 2 * el[m];
            ldof[2 * m + 1] = 2 * el[m] + 1;
        }
        for a in 0..6 {
            let mut cbq = [0.0_f64; 8];
            for q in 0..8 {
                let mut s = 0.0;
                for bb in 0..6 {
                    s += c[a][bb] * bmat[bb][q];
                }
                cbq[q] = s;
            }
            for p in 0..8 {
                let bp = bmat[a][p];
                if bp == 0.0 {
                    continue;
                }
                for q in 0..8 {
                    k_full[ldof[p]][ldof[q]] += bp * cbq[q] * det;
                }
            }
        }
    }

    // Reduced (free-dof) system.
    let mut red = vec![usize::MAX; ndof];
    let mut free = Vec::new();
    for d in 0..ndof {
        if !constrained[d] {
            red[d] = free.len();
            free.push(d);
        }
    }
    let nf = free.len();
    let mut k_red = vec![vec![0.0_f64; nf]; nf];
    for i in 0..ndof {
        if constrained[i] {
            continue;
        }
        for j in 0..ndof {
            if constrained[j] {
                continue;
            }
            k_red[red[i]][red[j]] = k_full[i][j];
        }
    }

    // External load: downward (negative y) on the right edge, equal share.
    let mut f_ext = vec![0.0_f64; nf];
    for j in 0..=N {
        let node = nid(N, j);
        let dof = 2 * node + 1;
        if !constrained[dof] {
            f_ext[red[dof]] += -P_TOTAL / (N + 1) as f64;
        }
    }

    // Incremental load stepping with Picard iteration. Plasticity is integrated
    // incrementally: each converged step applies return mapping to the *strain
    // increment* from the previous step (trial = σ_prev + C·Δε).
    let n_inc = 40;
    let max_iter = 30;
    let relax = 0.8;
    let mut u_full = vec![0.0_f64; ndof];
    let mut u_red = vec![0.0_f64; nf];
    let mut eps_prev = vec![[0.0_f64; 6]; ng];
    let mut sig_prev = vec![[0.0_f64; 6]; ng];
    let mut state_prev: Vec<PlasticState> = vec![PlasticState::default(); ng];

    for inc in 0..n_inc {
        let scale = (inc + 1) as f64 / n_inc as f64;
        let target: Vec<f64> = f_ext.iter().map(|v| v * scale).collect();

        let mut last_sig = sig_prev.clone();
        let mut last_state = state_prev.clone();
        let mut last_eps = eps_prev.clone();

        for _ in 0..max_iter {
            let mut f_int = vec![0.0_f64; nf];
            let mut new_sig = vec![[0.0_f64; 6]; ng];
            let mut new_state = vec![PlasticState::default(); ng];
            let mut new_eps = vec![[0.0_f64; 6]; ng];
            for gi in 0..ng {
                let e = gauss_elem[gi];
                let bmat = &gauss_b[gi];
                let det = gauss_det[gi];
                let el = elems[e];
                let mut ue = [0.0_f64; 8];
                for m in 0..4 {
                    ue[2 * m] = u_full[2 * el[m]];
                    ue[2 * m + 1] = u_full[2 * el[m] + 1];
                }
                let mut eps = [0.0_f64; 6];
                for a in 0..6 {
                    let mut s = 0.0;
                    for p in 0..8 {
                        s += bmat[a][p] * ue[p];
                    }
                    eps[a] = s;
                }
                // Trial stress for this increment.
                let mut trial = [0.0_f64; 6];
                for a in 0..6 {
                    let mut s = 0.0;
                    for bb in 0..6 {
                        s += c[a][bb] * (eps[bb] - eps_prev[gi][bb]);
                    }
                    trial[a] = sig_prev[gi][a] + s;
                }
                let (sig, st) = mat.return_map(&trial, &state_prev[gi]);
                // Internal force contribution B^T σ det.
                for a in 0..6 {
                    for p in 0..8 {
                        let val = bmat[a][p] * sig[a] * det;
                        if val == 0.0 {
                            continue;
                        }
                        let gdof = ldof_global(elems[e], p);
                        if !constrained[gdof] {
                            f_int[red[gdof]] += val;
                        }
                    }
                }
                new_sig[gi] = sig;
                new_state[gi] = st;
                new_eps[gi] = eps;
            }

            // Picard update.
            let mut rhs = vec![0.0_f64; nf];
            let mut norm = 0.0;
            for i in 0..nf {
                rhs[i] = target[i] - f_int[i];
                norm += rhs[i] * rhs[i];
            }
            norm = norm.sqrt();
            let du = solve(k_red.clone(), rhs);
            for i in 0..nf {
                u_red[i] += relax * du[i];
            }
            for d in 0..ndof {
                if !constrained[d] {
                    u_full[d] = u_red[red[d]];
                }
            }

            last_sig = new_sig;
            last_state = new_state;
            last_eps = new_eps;

            let tnorm: f64 = target.iter().map(|v| v * v).sum::<f64>().sqrt();
            if norm < 1e-6 * (1.0 + tnorm) {
                break;
            }
        }

        // Commit the converged incremental state.
        sig_prev = last_sig;
        state_prev = last_state;
        eps_prev = last_eps;
    }

    u_full[2 * nid(N, N) + 1]
}

#[inline]
fn ldof_global(el: [usize; 4], p: usize) -> usize {
    let m = p / 2;
    let comp = p % 2;
    2 * el[m] + comp
}

#[test]
fn elastic_plastic_cooks_membrane() {
    let u_elastic = cooks_membrane(1e6, 1e6); // effectively elastic (no yield)
    let u_plastic_soft = cooks_membrane(0.1, 0.2); // yields, low hardening
    let u_plastic_hard = cooks_membrane(0.1, 0.8); // yields, higher hardening

    assert!(
        u_elastic.is_finite() && u_plastic_soft.is_finite() && u_plastic_hard.is_finite(),
        "non-finite tip displacement"
    );
    // Plasticity makes the structure more compliant than the elastic solution.
    assert!(
        u_plastic_soft.abs() > u_elastic.abs(),
        "expected plastic compliance: |u_plastic| {} vs |u_elastic| {}",
        u_plastic_soft.abs(),
        u_elastic.abs()
    );
    // Hardening stiffens the response.
    assert!(
        u_plastic_hard.abs() < u_plastic_soft.abs(),
        "expected hardening to reduce displacement: {} vs {}",
        u_plastic_hard.abs(),
        u_plastic_soft.abs()
    );

    eprintln!(
        "Cook's membrane tip v: elastic={:.3}, plastic(H=0.5)={:.3}, plastic(H=2.0)={:.3}",
        u_elastic, u_plastic_soft, u_plastic_hard
    );
}
