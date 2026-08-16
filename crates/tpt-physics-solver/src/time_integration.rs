//! Time-integration schemes.
//!
//! Both transient (dynamic) schemes are net-new: `tpt-fem-eigen` and
//! `tpt-fem-elasticity::solve_modal` only do frequency-domain modal analysis,
//! so no time-stepping solver exists in either sibling repo.

use crate::cg::cg;
use crate::error::SolverError;
use crate::linalg::LinearOperator;

/// Classic 4th-order Runge–Kutta step for a first-order ODE `y' = f(t, y)`.
///
/// Returns `(t + dt, y_new)`.
pub fn rk4<F>(f: F, t: f64, y: &[f64], dt: f64) -> (f64, Vec<f64>)
where
    F: Fn(f64, &[f64], &mut [f64]),
{
    let n = y.len();
    let mut k1 = vec![0.0; n];
    f(t, y, &mut k1);

    let mut tmp = vec![0.0; n];
    for i in 0..n {
        tmp[i] = y[i] + 0.5 * dt * k1[i];
    }
    let mut k2 = vec![0.0; n];
    f(t + 0.5 * dt, &tmp, &mut k2);

    for i in 0..n {
        tmp[i] = y[i] + 0.5 * dt * k2[i];
    }
    let mut k3 = vec![0.0; n];
    f(t + 0.5 * dt, &tmp, &mut k3);

    for i in 0..n {
        tmp[i] = y[i] + dt * k3[i];
    }
    let mut k4 = vec![0.0; n];
    f(t + dt, &tmp, &mut k4);

    let mut ynew = vec![0.0; n];
    for i in 0..n {
        ynew[i] = y[i] + (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    (t + dt, ynew)
}

/// Effective stiffness operator `K_eff = a0·M + a1·C + K` used by Newmark.
struct Effective<'a> {
    m: &'a dyn LinearOperator,
    c: &'a dyn LinearOperator,
    k: &'a dyn LinearOperator,
    a0: f64,
    a1: f64,
    n: usize,
}

impl LinearOperator for Effective<'_> {
    fn nrows(&self) -> usize {
        self.n
    }
    fn ncols(&self) -> usize {
        self.n
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        let mut tmp = vec![0.0; self.n];
        self.k.apply(x, y);
        self.m.apply(x, &mut tmp);
        for i in 0..self.n {
            y[i] += self.a0 * tmp[i];
        }
        self.c.apply(x, &mut tmp);
        for i in 0..self.n {
            y[i] += self.a1 * tmp[i];
        }
    }
}

/// Newmark-beta integrator for second-order structural dynamics
/// `M a + C v + K u = f(t)`.
///
/// With `β = 1/4, γ = 1/2` the scheme is unconditionally stable and
/// second-order accurate. The effective system at each step is solved with the
/// in-crate Conjugate Gradient solver (no direct-factorisation dependency).
pub struct NewmarkBeta {
    m: Box<dyn LinearOperator>,
    c: Box<dyn LinearOperator>,
    k: Box<dyn LinearOperator>,
    dt: f64,
    beta: f64,
    gamma: f64,
    n: usize,
    /// Current displacement `uₙ`.
    pub u: Vec<f64>,
    /// Current velocity `vₙ`.
    pub v: Vec<f64>,
    /// Current acceleration `aₙ`.
    pub a: Vec<f64>,
    /// Linear-solver tolerance.
    pub tol: f64,
    /// Linear-solver iteration cap.
    pub max_iter: usize,
}

impl NewmarkBeta {
    /// Build an integrator. `m`, `c`, `k` are the mass, damping, and stiffness
    /// operators (all `n×n`). `u0`, `v0` are the initial displacement and
    /// velocity; `a0` is computed from the initial equilibrium
    /// `a0 = M⁻¹ (f0 - C v0 - K u0)` if `None` is supplied.
    pub fn new(
        m: Box<dyn LinearOperator>,
        c: Box<dyn LinearOperator>,
        k: Box<dyn LinearOperator>,
        dt: f64,
        beta: f64,
        gamma: f64,
        u0: Vec<f64>,
        v0: Vec<f64>,
        a0: Option<Vec<f64>>,
    ) -> Result<Self, SolverError> {
        let n = m.nrows();
        if m.ncols() != n || c.nrows() != n || c.ncols() != n || k.nrows() != n || k.ncols() != n {
            return Err(SolverError::NotSquare {
                nrows: n,
                ncols: k.ncols(),
            });
        }
        let a0v = match a0 {
            Some(a) => a,
            None => {
                // a0 = M⁻¹ (f0 - C v0 - K u0) with f0 = 0  =>  a0 = -M⁻¹(C v0 + K u0).
                let mut kv = vec![0.0; n];
                k.apply(&u0, &mut kv);
                let mut cv = vec![0.0; n];
                c.apply(&v0, &mut cv);
                for i in 0..n {
                    kv[i] += cv[i];
                }
                // Solve M a0 = -kv via CG.
                let (a, _) = cg(
                    &*m,
                    &kv.iter().map(|x| -x).collect::<Vec<_>>(),
                    None,
                    1e-10,
                    200,
                )?;
                a
            }
        };
        Ok(NewmarkBeta {
            m,
            c,
            k,
            dt,
            beta,
            gamma,
            n,
            u: u0,
            v: v0,
            a: a0v,
            tol: 1e-10,
            max_iter: 500,
        })
    }

    /// Advance one step with the external force `f_next` at `t + dt`.
    pub fn step(&mut self, f_next: &[f64]) -> Result<(), SolverError> {
        let dt = self.dt;
        let beta = self.beta;
        let gamma = self.gamma;

        // Newmark constants (Bathe).
        let a0 = 1.0 / (beta * dt * dt);
        let a1 = gamma / (beta * dt);
        let a2 = 1.0 / (beta * dt);
        let a3 = 1.0 / (2.0 * beta) - 1.0;
        let a4 = gamma / beta - 1.0;
        let a5 = dt * (gamma / (2.0 * beta) - 1.0);

        // RHS: f_eff = f_next + M(a0 u + a2 v + a3 a) + C(a1 u + a4 v + a5 a).
        let mut mu = vec![0.0; self.n];
        let mut mc = vec![0.0; self.n];
        {
            let mut t1 = vec![0.0; self.n];
            for i in 0..self.n {
                t1[i] = a0 * self.u[i] + a2 * self.v[i] + a3 * self.a[i];
            }
            self.m.apply(&t1, &mut mu);
            let mut t2 = vec![0.0; self.n];
            for i in 0..self.n {
                t2[i] = a1 * self.u[i] + a4 * self.v[i] + a5 * self.a[i];
            }
            self.c.apply(&t2, &mut mc);
        }
        let mut feff = vec![0.0; self.n];
        for i in 0..self.n {
            feff[i] = f_next[i] + mu[i] + mc[i];
        }

        let eff = Effective {
            m: &*self.m,
            c: &*self.c,
            k: &*self.k,
            a0,
            a1,
            n: self.n,
        };
        let (u_next, _) = cg(&eff, &feff, None, self.tol, self.max_iter)?;

        let mut a_next = vec![0.0; self.n];
        for i in 0..self.n {
            a_next[i] = a0 * (u_next[i] - self.u[i]) - a2 * self.v[i] - a3 * self.a[i];
        }
        let mut v_next = vec![0.0; self.n];
        for i in 0..self.n {
            v_next[i] = self.v[i] + dt * ((1.0 - gamma) * self.a[i] + gamma * a_next[i]);
        }

        self.u = u_next;
        self.v = v_next;
        self.a = a_next;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::{csr_from_dense, DenseMat};

    #[test]
    fn rk4_oscillator() {
        // x'' = -x, state [x, v]; analytic period 2π.
        let f = |_t: f64, y: &[f64], out: &mut [f64]| {
            out[0] = y[1];
            out[1] = -y[0];
        };
        let mut t = 0.0;
        let mut y = vec![1.0, 0.0];
        let dt = 0.01;
        for _ in 0..628 {
            // ~ one period
            let (nt, ny) = rk4(&f, t, &y, dt);
            t = nt;
            y = ny;
        }
        assert!((y[0] - 1.0).abs() < 5e-3, "x at T = {}", y[0]);
        assert!((y[1] - 0.0).abs() < 5e-3, "v at T = {}", y[1]);
    }

    #[test]
    fn newmark_undamped_oscillator() {
        // m x'' + k x = 0, m = k = 1, x(0)=1, v(0)=0 => x(t)=cos(t).
        let m = Box::new(DenseMat::from_row_major(1, 1, vec![1.0]));
        let c = Box::new(DenseMat::from_row_major(1, 1, vec![0.0]));
        let k = Box::new(DenseMat::from_row_major(1, 1, vec![1.0]));
        let dt = 0.01;
        let mut integ = NewmarkBeta::new(m, c, k, dt, 0.25, 0.5, vec![1.0], vec![0.0], None)
            .expect("newmark build");
        let steps = 628; // ≈ 2π
        for _ in 0..steps {
            integ.step(&[0.0]).expect("step");
        }
        assert!((integ.u[0] - 1.0).abs() < 1e-2, "x(2π) = {}", integ.u[0]);
    }

    #[test]
    fn newmark_force_response() {
        // m x'' + k x = F (constant), undamped. With m = k = 1, x(0) = v(0) = 0,
        // the analytic response is x(t) = 1 - cos(t) (steady state 1 plus the
        // homogeneous oscillation). Verify the integrator tracks it.
        let m = Box::new(DenseMat::from_row_major(1, 1, vec![1.0]));
        let c = Box::new(DenseMat::from_row_major(1, 1, vec![0.0]));
        let k = Box::new(DenseMat::from_row_major(1, 1, vec![1.0]));
        let dt = 0.02;
        let steps = 2000;
        let mut integ =
            NewmarkBeta::new(m, c, k, dt, 0.25, 0.5, vec![0.0], vec![0.0], None).unwrap();
        for _ in 0..steps {
            integ.step(&[1.0]).unwrap();
        }
        let t = steps as f64 * dt;
        let analytic = 1.0 - t.cos();
        assert!(
            (integ.u[0] - analytic).abs() < 2e-2,
            "x({t}) = {}, analytic {}",
            integ.u[0],
            analytic
        );
    }

    #[test]
    fn newmark_rejects_mismatched_dims() {
        let m = Box::new(csr_from_dense(2, 2, &[1.0, 0.0, 0.0, 1.0]));
        let c = Box::new(csr_from_dense(2, 2, &[0.0; 4]));
        let k = Box::new(csr_from_dense(1, 1, &[1.0])); // wrong size
        let r = NewmarkBeta::new(m, c, k, 0.01, 0.25, 0.5, vec![0.0; 2], vec![0.0; 2], None);
        assert!(matches!(r, Err(SolverError::NotSquare { .. })));
    }
}
