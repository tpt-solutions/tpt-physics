//! Hello-world for `tpt-physics-fea`: a clamped-free (cantilever) 3-D beam
//! under a tip load, solved with the net-new Euler–Bernoulli frame element.
//!
//! Run with: `cargo run --example beam --release`

use tpt_physics_fea::elements::beam3d_global_stiffness;

fn main() {
    // A 1 m steel cantilever, 10 mm square cross-section, tip load P = 1 kN.
    let e = 200e9;
    let nu = 0.3;
    let g = e / (2.0 * (1.0 + nu));
    let l = 1.0;
    let b = 0.01;
    let h = 0.01;
    let area = b * h;
    let iz = b * h * h * h / 12.0; // bending about z (gives y deflection)
    let iy = h * b * b * b / 12.0;
    let jtor = area * (b * b + h * h) / 12.0;

    // Node 0 fixed at origin, node 1 free at x = l.
    let kg = beam3d_global_stiffness(
        [0.0, 0.0, 0.0],
        [l, 0.0, 0.0],
        [0.0, 0.0, 1.0], // up vector → orients local y/z
        e,
        g,
        area,
        iy,
        iz,
        jtor,
    );

    // Free DOFs are node 1 (indices 6..12); node 0 is fully clamped (0..6).
    let free: Vec<usize> = (6..12).collect();
    let p = 1000.0;
    let mut kf = [[0.0_f64; 6]; 6];
    let mut ff = [0.0_f64; 6];
    for (i, &fi) in free.iter().enumerate() {
        for (j, &fj) in free.iter().enumerate() {
            kf[i][j] = kg[fi * 12 + fj];
        }
    }
    ff[1] = -p; // downward (−y) tip load

    let mut uf = [0.0_f64; 6];
    solve_6x6(&kf, &ff, &mut uf);

    // Analytic tip deflection (Euler–Bernoulli): δ = P L³ / (3 E I_z).
    let expected = p * l * l * l / (3.0 * e * iz);
    let tip_dy = uf[1];

    println!("Cantilever beam (L={l} m, P={p} N, steel):");
    println!("  analytic tip deflection : {expected:.3e} m");
    println!("  FEA tip deflection     : {tip_dy:.3e} m");
    println!(
        "  relative error         : {:.2e}",
        (tip_dy - (-expected)).abs() / expected
    );
}

/// Minimal dense 6×6 linear solver (partial-pivot Gaussian elimination).
fn solve_6x6(a: &[[f64; 6]; 6], b: &[f64; 6], x: &mut [f64; 6]) -> bool {
    let mut m = [[0.0_f64; 7]; 6];
    for i in 0..6 {
        for j in 0..6 {
            m[i][j] = a[i][j];
        }
        m[i][6] = b[i];
    }
    for c in 0..6 {
        let mut piv = c;
        let mut best = m[c][c].abs();
        for r in (c + 1)..6 {
            if m[r][c].abs() > best {
                best = m[r][c].abs();
                piv = r;
            }
        }
        if best < 1e-15 {
            return false;
        }
        m.swap(c, piv);
        let d = m[c][c];
        for j in c..7 {
            m[c][j] /= d;
        }
        for r in 0..6 {
            if r != c {
                let f = m[r][c];
                for j in c..7 {
                    m[r][j] -= f * m[c][j];
                }
            }
        }
    }
    for i in 0..6 {
        x[i] = m[i][6];
    }
    true
}
