//! Wall-load mapping in flux coordinates: closed-form per-cell accumulation
//! of a caller birth-rate density field onto a caller wall surface sharing
//! the same flux coordinates, weighted by the equilibrium Jacobian.
//!
//! Reads no files and solves nothing: the caller supplies the radial mesh,
//! the angular grid, the per-voxel birth-rate densities (from any source
//! model), and the per-voxel Jacobian values (evaluated with
//! [`crate::wout::fourier_jacobian`] or [`crate::wout::Wout::jacobian_grid`]);
//! this module is the accumulation math plus gates. There is no transport
//! here — no field-line tracing, no shadowing, no FEM/thermal step, no CAD:
//! each wall node sees the radial line integral at its own angles.
//!
//! Coordinates: `s` is the normalized toroidal-flux label (`s = 0` on axis,
//! `s = 1` at the edge), `θ ∈ [0, 2π)` poloidal, `ζ` toroidal over one
//! field period (`[0, 2π/nfp)`, the file's toroidal convention — the caller
//! scales totals for full-torus work exactly as for the J1–J3 helpers).
//! The wall surface sits at the outer radial edge and shares the `(θ, ζ)`
//! grid, so one wall node exists per angular voxel column.
//!
//! Equation labels:
//!
//! - (W1) accumulation: `q[j][k] = Σ_i S[i][j][k] · J[i][j][k] · Δs[i]`,
//!   with `S` the birth-rate density, `J = √g` the Jacobian, and `Δs[i]`
//!   the `i`-th radial cell width. Voxels are radial-major:
//!   `index = (i * ntheta + j) * nzeta + k`.
//! - (W2) one-field-period total: `Q = Σ_jk q[j][k] · Δθ · Δζ`, which
//!   equals the voxel sum `Σ_ijk S·J·Δs·Δθ·Δζ` by construction
//!   (discrete conservation to roundoff).
//! - (W3) axisymmetric limit: uniform `S`, uniform `J = J0` give the
//!   spatially constant analytic wall flux `q = S · J0 · (s_hi − s_lo)`
//!   over the radial span.

use crate::error::{Error, Result};

/// Wall loads over one field period: one radial line integral (W1) per
/// wall node plus the angle-weighted total (W2).
#[derive(Debug, Clone, PartialEq)]
pub struct WallLoad {
    /// Poloidal node count (outer over `θ`).
    pub ntheta: usize,
    /// Toroidal node count over one field period (inner over `ζ`).
    pub nzeta: usize,
    /// Per-node loads `loads[j * nzeta + k]` (same units as the caller
    /// birth-rate density times length: the `Δs` weight is dimensionless).
    pub loads: Vec<f64>,
    /// One-field-period total `Q` (W2): `Σ loads · Δθ · Δζ`.
    pub total: f64,
}

impl WallLoad {
    /// Borrow the load at wall node `(j, k)` (loud on out-of-range nodes).
    pub fn at(&self, j: usize, k: usize) -> Result<f64> {
        if j >= self.ntheta {
            return Err(Error::OutOfRange {
                what: "wall theta node",
                index: j,
                len: self.ntheta,
            });
        }
        if k >= self.nzeta {
            return Err(Error::OutOfRange {
                what: "wall zeta node",
                index: k,
                len: self.nzeta,
            });
        }
        Ok(self.loads[j * self.nzeta + k])
    }

    /// Split into nested rows (outer over `θ`, inner over `ζ`).
    pub fn to_nested(&self) -> Vec<Vec<f64>> {
        self.loads.chunks(self.nzeta).map(<[f64]>::to_vec).collect()
    }
}

/// Map a birth-rate density field onto the wall (W1–W2).
///
/// `s_edges` holds the `nr + 1` strictly increasing radial edges,
/// `ntheta`/`nzeta` the angular grid counts, `nfp` the field-period count
/// (positive; sets the one-field-period `Δζ = 2π/(nfp·nzeta)` behind the
/// total), and `birth`/`jacobian` the per-voxel birth-rate densities and
/// `√g` values in radial-major order (`(i * ntheta + j) * nzeta + k`).
/// Length mismatches are [`Error::BadShape`]; empty grids, a non-positive
/// `nfp`, non-increasing or non-finite edges, non-finite or negative
/// densities, and non-finite or negative Jacobians are [`Error::BadValue`].
/// A zero Jacobian is allowed (null-volume voxel); a negative one is not.
pub fn wall_load(
    s_edges: &[f64],
    ntheta: usize,
    nzeta: usize,
    nfp: u32,
    birth: &[f64],
    jacobian: &[f64],
) -> Result<WallLoad> {
    if ntheta == 0 {
        return Err(Error::BadValue("ntheta must be positive".into()));
    }
    if nzeta == 0 {
        return Err(Error::BadValue("nzeta must be positive".into()));
    }
    if nfp == 0 {
        return Err(Error::BadValue("nfp must be positive".into()));
    }
    if s_edges.len() < 2 {
        return Err(Error::BadShape {
            what: "s_edges",
            expected: 2,
            got: s_edges.len(),
        });
    }
    for pair in s_edges.windows(2) {
        if !pair[0].is_finite() || !pair[1].is_finite() {
            return Err(Error::BadValue("non-finite s edge".into()));
        }
        if pair[1] <= pair[0] {
            return Err(Error::BadValue(
                "s_edges must be strictly increasing".into(),
            ));
        }
    }
    let nr = s_edges.len() - 1;
    let n_vox = nr * ntheta * nzeta;
    if birth.len() != n_vox {
        return Err(Error::BadShape {
            what: "birth",
            expected: n_vox,
            got: birth.len(),
        });
    }
    if jacobian.len() != n_vox {
        return Err(Error::BadShape {
            what: "jacobian",
            expected: n_vox,
            got: jacobian.len(),
        });
    }
    for &s in birth {
        if !s.is_finite() {
            return Err(Error::BadValue("non-finite birth-rate density".into()));
        }
        if s < 0.0 {
            return Err(Error::BadValue("negative birth-rate density".into()));
        }
    }
    for &j in jacobian {
        if !j.is_finite() {
            return Err(Error::BadValue("non-finite Jacobian".into()));
        }
        if j < 0.0 {
            return Err(Error::BadValue("negative Jacobian".into()));
        }
    }

    let mut loads = vec![0.0; ntheta * nzeta];
    for i in 0..nr {
        let ds = s_edges[i + 1] - s_edges[i];
        for j in 0..ntheta {
            for k in 0..nzeta {
                let vox = (i * ntheta + j) * nzeta + k;
                loads[j * nzeta + k] += birth[vox] * jacobian[vox] * ds;
            }
        }
    }
    for &q in &loads {
        if !q.is_finite() {
            return Err(Error::BadValue("non-finite wall load".into()));
        }
    }
    let dtheta = 2.0 * std::f64::consts::PI / ntheta as f64;
    let dzeta = 2.0 * std::f64::consts::PI / (nzeta as f64 * f64::from(nfp));
    let total = loads.iter().sum::<f64>() * dtheta * dzeta;
    if !total.is_finite() {
        return Err(Error::BadValue("non-finite wall total".into()));
    }
    Ok(WallLoad {
        ntheta,
        nzeta,
        loads,
        total,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-12;

    fn close(got: f64, want: f64) -> bool {
        (got - want).abs() <= TOL * want.abs().max(1.0)
    }

    #[test]
    fn axisymmetric_uniform_matches_analytic_wall_flux() {
        // (W3): S = 2, J0 = 3 over s in [0, 1] give q = 6 everywhere;
        // Q = 6 * 12 nodes * (pi/2) * (pi/3) = 12 * pi^2 (nfp = 2).
        let birth = vec![2.0; 2 * 4 * 3];
        let jac = vec![3.0; 2 * 4 * 3];
        let wall = wall_load(&[0.0, 0.5, 1.0], 4, 3, 2, &birth, &jac).unwrap();
        assert_eq!((wall.ntheta, wall.nzeta), (4, 3));
        assert!(wall.loads.iter().all(|&q| close(q, 6.0)));
        let want_total = 12.0 * std::f64::consts::PI.powi(2);
        assert!(
            close(wall.total, want_total),
            "total {} vs hand value {want_total}",
            wall.total
        );
        // Nested view and node access agree.
        assert_eq!(wall.to_nested().len(), 4);
        assert_eq!(wall.at(3, 2).unwrap(), wall.loads[3 * 3 + 2]);
    }

    #[test]
    fn cosine_jacobian_hand_vector() {
        // One radial cell, S = 1, J = 1 + 0.5*cos(2θ) sampled at the four
        // cardinal poloidal nodes (nzeta = 1): hand loads [1.5, 0.5, 1.5, 0.5].
        let jac = vec![1.5, 0.5, 1.5, 0.5];
        let wall = wall_load(&[0.0, 1.0], 4, 1, 1, &[1.0; 4], &jac).unwrap();
        for (got, want) in wall.loads.iter().zip([1.5, 0.5, 1.5, 0.5]) {
            assert!(close(*got, want), "got {got} want {want}");
        }
        // Radial weighting: S = 4 below s = 0.5, S = 2 above, J = 1.
        let birth = vec![4.0, 4.0, 2.0, 2.0];
        let wall = wall_load(&[0.0, 0.5, 1.0], 2, 1, 1, &birth, &[1.0; 4]).unwrap();
        assert!(wall.loads.iter().all(|&q| close(q, 3.0)));
    }

    #[test]
    fn total_conserves_discrete_births() {
        // (W2): the angle-weighted total equals the direct voxel sum.
        let s_edges = vec![0.0, 0.25, 1.0];
        let birth = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let jac = vec![1.0, 0.5, 2.0, 1.5, 1.0, 1.0, 0.25, 2.0];
        let wall = wall_load(&s_edges, 2, 2, 3, &birth, &jac).unwrap();
        let dtheta = 2.0 * std::f64::consts::PI / 2.0;
        let dzeta = 2.0 * std::f64::consts::PI / (2.0 * 3.0);
        let mut want = 0.0;
        for i in 0..2 {
            let ds = s_edges[i + 1] - s_edges[i];
            for v in 0..4 {
                want += birth[i * 4 + v] * jac[i * 4 + v] * ds * dtheta * dzeta;
            }
        }
        assert!(close(wall.total, want));
        assert!(wall.at(2, 0).is_err());
        assert!(wall.at(0, 2).is_err());
    }

    #[test]
    fn malformed_inputs_fail_loudly() {
        let good_birth = vec![1.0; 4];
        let good_jac = vec![1.0; 4];
        // Empty / degenerate grids and counts.
        assert!(wall_load(&[], 2, 2, 1, &good_birth, &good_jac).is_err());
        assert!(wall_load(&[0.0], 2, 2, 1, &good_birth, &good_jac).is_err());
        assert!(wall_load(&[0.0, 1.0], 0, 2, 1, &good_birth, &good_jac).is_err());
        assert!(wall_load(&[0.0, 1.0], 2, 0, 1, &good_birth, &good_jac).is_err());
        assert!(wall_load(&[0.0, 1.0], 2, 2, 0, &good_birth, &good_jac).is_err());
        // Non-increasing, repeated, and non-finite edges.
        assert!(wall_load(&[0.0, 1.0, 0.5], 2, 1, 1, &good_birth, &good_jac).is_err());
        assert!(wall_load(&[0.0, 0.0, 1.0], 2, 1, 1, &good_birth, &good_jac).is_err());
        assert!(wall_load(&[0.0, f64::NAN], 2, 1, 1, &[1.0; 2], &[1.0; 2]).is_err());
        // Shape mismatches name the offending input.
        let err = wall_load(&[0.0, 1.0], 2, 2, 1, &[1.0; 3], &good_jac).unwrap_err();
        assert_eq!(
            err,
            Error::BadShape {
                what: "birth",
                expected: 4,
                got: 3
            }
        );
        assert!(wall_load(&[0.0, 1.0], 2, 2, 1, &good_birth, &[1.0; 5]).is_err());
        // Non-finite and negative densities / Jacobians.
        let mut bad = good_birth.clone();
        bad[1] = f64::NAN;
        assert!(wall_load(&[0.0, 1.0], 2, 2, 1, &bad, &good_jac).is_err());
        let mut bad = good_birth.clone();
        bad[1] = -1.0;
        assert!(wall_load(&[0.0, 1.0], 2, 2, 1, &bad, &good_jac).is_err());
        let mut bad = good_jac.clone();
        bad[0] = f64::INFINITY;
        assert!(wall_load(&[0.0, 1.0], 2, 2, 1, &good_birth, &bad).is_err());
        let mut bad = good_jac;
        bad[2] = -0.5;
        assert!(wall_load(&[0.0, 1.0], 2, 2, 1, &good_birth, &bad).is_err());
        // Zero Jacobians are allowed (null-volume voxel).
        let wall = wall_load(&[0.0, 1.0], 2, 2, 1, &good_birth, &[0.0; 4]).unwrap();
        assert!(wall.loads.iter().all(|&q| q == 0.0));
        assert_eq!(wall.total, 0.0);
    }
}
