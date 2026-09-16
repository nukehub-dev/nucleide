//! 1D cell-centred finite-volume + theta-stepping solver for (T1–T2).
//!
//! The slab `0 ≤ x ≤ L` carries `cells` finite volumes (`dx = L / cells`)
//! holding the mobile concentration `c_m` plus one trapped population per
//! trap species. Diffusion assembles the tridiagonal rate matrix `A` with
//! `dc/dt = A c + b` (`b` folds Dirichlet/Sieverts/Henry face values and
//! the volumetric source); each implicit theta step
//!
//! ```text
//! (I − dt·θ·A) cⁿ⁺¹ = (I + dt·(1−θ)·A) cⁿ + dt·b − Σⱼ(ctⱼ* − ctⱼⁿ)
//! ```
//!
//! solves through [`nucleide_linalg::tridiag`] (Thomas, O(N)) while the
//! per-cell trap populations update by the exact backward-Euler map of
//! (T2) at the latest mobile iterate:
//!
//! ```text
//! ctⱼⁿ⁺¹ = (ctⱼⁿ + dt·kⱼ·cⁿ⁺¹·Nⱼ) / (1 + dt·(kⱼ·cⁿ⁺¹ + pⱼ)),
//! ```
//!
//! which is positivity-preserving and bounded by `Nⱼ`. Mobile and traps
//! couple by Picard iteration to `rtol`/`atol` (the kinetics `solve.rs`
//! pattern: options struct with `validate()`, exact stepping onto output
//! times, a hard step budget); trap-free systems take a single solve per
//! step. Temperature is the caller-supplied steady profile from
//! [`TransportParams`](crate::params::TransportParams) — no heat solve in
//! v1. Recombination ends (`J = K_r c²`) close per implicit step (G6) by
//! the same affine face-response trick as the steady state: for frozen
//! faces the step is linear, so one base Thomas solve plus one sensitivity
//! column per recombination end pins the face-adjacent response exactly
//! and the faces close in closed form (one end) or Newton iteration with
//! the analytic Jacobian (two ends). The face iterate lives inside the
//! trap Picard loop (base solve → face close → trap map per pass) and
//! shares `rtol`/`atol`; a face Newton that exhausts its cap fails the
//! step with [`Error::NotConverged`] (hard fail, no dt retry).
//!
//! Closed-form helpers pin the algebraic gates without a solve:
//! [`time_lag`] (G2-lag), [`breakthrough_ratio`] (G2 series),
//! [`equilibrium_trapped`] (T2-eq, G3b), [`effective_diffusivity`] (G3a),
//! [`irreversible_fill`] (G3c), [`sieverts_concentration`] (G4).

use crate::bc::Boundary;
use crate::error::Error;
use crate::params::TransportParams;

/// Theta-method selector: `1/2` (Crank–Nicolson, default) or `1` (backward
/// Euler).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Theta {
    /// Second-order A-stable Crank–Nicolson (default).
    #[default]
    CrankNicolson,
    /// First-order L-stable backward Euler (extra damping of fast trap
    /// transients; larger phase error on the breakthrough tail).
    BackwardEuler,
}

impl Theta {
    /// Theta value (`1/2` or `1`).
    pub fn value(self) -> f64 {
        match self {
            Theta::CrankNicolson => 0.5,
            Theta::BackwardEuler => 1.0,
        }
    }
}

/// Solver tolerances and budgets.
#[derive(Debug, Clone, PartialEq)]
pub struct SolverOptions {
    /// Implicit integration method (theta value).
    pub theta: Theta,
    /// Relative tolerance for the per-step mobile/trap/face Picard coupling
    /// (recombination faces share these knobs; there are no dedicated face
    /// tolerances).
    pub rtol: f64,
    /// Absolute tolerance (in concentration units) for the Picard coupling.
    pub atol: f64,
    /// Smallest allowed internal step \[s\].
    pub dt_min: f64,
    /// Largest allowed internal step \[s\].
    pub dt_max: f64,
    /// Maximum accepted internal steps over the whole grid.
    pub max_steps: usize,
}

impl Default for SolverOptions {
    fn default() -> Self {
        Self {
            theta: Theta::CrankNicolson,
            rtol: 1e-9,
            atol: 1e-12,
            dt_min: 1e-14,
            dt_max: f64::INFINITY,
            max_steps: 1_000_000,
        }
    }
}

impl SolverOptions {
    /// Validate tolerances and budgets (finite, positive, ordered).
    pub fn validate(&self) -> Result<(), Error> {
        if !self.rtol.is_finite() || self.rtol <= 0.0 {
            return Err(Error::BadOption("rtol must be finite and > 0"));
        }
        if !self.atol.is_finite() || self.atol <= 0.0 {
            return Err(Error::BadOption("atol must be finite and > 0"));
        }
        if !self.dt_min.is_finite() || self.dt_min <= 0.0 {
            return Err(Error::BadOption("dt_min must be finite and > 0"));
        }
        if self.dt_max <= 0.0 || self.dt_max < self.dt_min {
            return Err(Error::BadOption("need 0 < dt_min <= dt_max"));
        }
        if self.max_steps == 0 {
            return Err(Error::BadOption("max_steps must be > 0"));
        }
        Ok(())
    }
}

/// Output time grid: strictly increasing times `> 0`.
#[derive(Debug, Clone, PartialEq)]
pub struct TimeGrid {
    /// Output times \[s\].
    pub times: Vec<f64>,
}

impl TimeGrid {
    /// Validate: at least one time, all finite and `> 0`, strictly
    /// increasing.
    pub fn new(times: Vec<f64>) -> Result<Self, Error> {
        if times.is_empty() {
            return Err(Error::BadGrid("time grid must not be empty"));
        }
        let mut prev = 0.0_f64;
        for t in &times {
            if !t.is_finite() || *t <= prev {
                return Err(Error::BadGrid(
                    "times must be finite, > 0, strictly increasing",
                ));
            }
            prev = *t;
        }
        Ok(Self { times })
    }
}

/// Initial state at `t = 0`: mobile profile plus per-cell trapped loads.
#[derive(Debug, Clone, PartialEq)]
pub struct InitialState {
    /// Mobile concentration per cell \[mol/m³\] (`>= 0`, finite).
    pub mobile: Vec<f64>,
    /// Trapped concentrations `[cell][trap]` \[mol/m³\] (`0 <= ct <= N_j`).
    pub trapped: Vec<Vec<f64>>,
}

impl InitialState {
    /// Zero initial state for the params geometry (mobile and traps zero).
    pub fn zeros(params: &TransportParams) -> Self {
        Self {
            mobile: vec![0.0; params.cells],
            trapped: vec![vec![0.0; params.traps.len()]; params.cells],
        }
    }

    /// Uniform mobile value with per-species uniform trapped loads.
    pub fn uniform(
        params: &TransportParams,
        mobile_value: f64,
        trapped_values: Vec<f64>,
    ) -> Result<Self, Error> {
        if trapped_values.len() != params.traps.len() {
            return Err(Error::BadState(
                "trapped values length must match trap count",
            ));
        }
        Self::new(
            params,
            vec![mobile_value; params.cells],
            vec![trapped_values; params.cells],
        )
    }

    /// Validate raw profiles against the params geometry: mobile length
    /// matches the cell count (finite, `>= 0`); trapped is
    /// `[cell][trap]`-shaped (finite, `0 <= ct <= N_j`).
    pub fn new(
        params: &TransportParams,
        mobile: Vec<f64>,
        trapped: Vec<Vec<f64>>,
    ) -> Result<Self, Error> {
        if mobile.len() != params.cells {
            return Err(Error::BadState("mobile length must match cell count"));
        }
        if mobile.iter().any(|c| !c.is_finite() || *c < 0.0) {
            return Err(Error::BadState("mobile must be finite and >= 0"));
        }
        if trapped.len() != params.cells {
            return Err(Error::BadState(
                "trapped outer length must match cell count",
            ));
        }
        for (cell, row) in trapped.iter().enumerate() {
            if row.len() != params.traps.len() {
                return Err(Error::BadState(
                    "trapped inner length must match trap count",
                ));
            }
            for (j, ct) in row.iter().enumerate() {
                let n = params.traps[j].site_density;
                if !ct.is_finite() || *ct < 0.0 || *ct > n {
                    return Err(Error::BadState("trapped loads must satisfy 0 <= ct <= N_j"));
                }
            }
            let _ = cell;
        }
        Ok(Self { mobile, trapped })
    }
}

/// Full transient: one mobile/trapped row plus outward surface fluxes per
/// output time (the `t = 0` state is the caller-supplied
/// [`InitialState`], echoed here as `initial`).
#[derive(Debug, Clone, PartialEq)]
pub struct Solution {
    /// Output times \[s\] (echo of the grid).
    pub times: Vec<f64>,
    /// Mobile concentration at each output time (`[time][cell]`).
    pub mobile: Vec<Vec<f64>>,
    /// Trapped concentrations (`[time][cell][trap]`).
    pub trapped: Vec<Vec<Vec<f64>>>,
    /// Outward flux at `x = 0` \[mol/m²/s\] (positive leaves the slab).
    pub flux_left: Vec<f64>,
    /// Outward flux at `x = L` \[mol/m²/s\] (positive leaves the slab).
    pub flux_right: Vec<f64>,
    /// Cell width `dx` \[m\] (for inventory integrals).
    pub dx: f64,
    /// The `t = 0` state.
    pub initial: InitialState,
}

impl Solution {
    /// Mobile inventory per unit area at each output time \[mol/m²\].
    pub fn inventory_mobile(&self) -> Vec<f64> {
        self.mobile
            .iter()
            .map(|row| row.iter().sum::<f64>() * self.dx)
            .collect()
    }

    /// Trapped inventory per unit area at each output time \[mol/m²\].
    pub fn inventory_trapped(&self) -> Vec<f64> {
        self.trapped
            .iter()
            .map(|rows| rows.iter().map(|row| row.iter().sum::<f64>()).sum::<f64>() * self.dx)
            .collect()
    }

    /// Total (mobile + trapped) inventory per unit area \[mol/m²\].
    pub fn inventory_total(&self) -> Vec<f64> {
        let m = self.inventory_mobile();
        let t = self.inventory_trapped();
        m.iter().zip(&t).map(|(a, b)| a + b).collect()
    }
}

/// Trap-free-style steady state: mobile profile, Langmuir trapped loads,
/// boundary fluxes, and inventories.
///
/// At steady state `dct/dt = 0`, so every trap sits on its Langmuir
/// isotherm (T2-eq) and the mobile profile solves the trap-free system
/// `A c + b = 0` (G1/G4); the trapped loads follow pointwise (G3b).
#[derive(Debug, Clone, PartialEq)]
pub struct SteadyState {
    /// Cell-centre positions \[m\].
    pub centres: Vec<f64>,
    /// Mobile concentration per cell \[mol/m³\].
    pub mobile: Vec<f64>,
    /// Trapped concentrations `[cell][trap]` \[mol/m³\].
    pub trapped: Vec<Vec<f64>>,
    /// Outward flux at `x = 0` \[mol/m²/s\].
    pub flux_left: f64,
    /// Outward flux at `x = L` \[mol/m²/s\].
    pub flux_right: f64,
    /// Mobile inventory per unit area \[mol/m²\].
    pub inventory_mobile: f64,
    /// Trapped inventory per unit area \[mol/m²\].
    pub inventory_trapped: f64,
}

// ---------------------------------------------------------------------------
// Closed-form helpers (algebraic gates)
// ---------------------------------------------------------------------------

/// Permeation time lag `t_lag = L²/6D` \[s\] (G2-lag).
///
/// Rejects non-finite or non-positive inputs.
pub fn time_lag(length: f64, diffusivity: f64) -> Result<f64, Error> {
    if !length.is_finite() || length <= 0.0 {
        return Err(Error::BadData("length must be finite and > 0"));
    }
    if !diffusivity.is_finite() || diffusivity <= 0.0 {
        return Err(Error::BadData("diffusivity must be finite and > 0"));
    }
    Ok(length * length / (6.0 * diffusivity))
}

/// Normalized outlet flux `J(L,t)/J_ss` of the G2 permeation transient
/// (trap-free slab, `c(0,t) = c0`, `c(L,t) = 0`, `c(x,0) = 0`):
///
/// ```text
/// J/Jss = 1 + 2 Σ_{n≥1} (−1)ⁿ exp(−D n² π² t / L²).
/// ```
///
/// The series truncates once `|term| < 1e-18` (the fixture convention).
/// Rejects non-finite/non-positive `D`/`L` and non-finite/non-positive `t`.
pub fn breakthrough_ratio(diffusivity: f64, length: f64, t: f64) -> Result<f64, Error> {
    if !diffusivity.is_finite() || diffusivity <= 0.0 {
        return Err(Error::BadData("diffusivity must be finite and > 0"));
    }
    if !length.is_finite() || length <= 0.0 {
        return Err(Error::BadData("length must be finite and > 0"));
    }
    if !t.is_finite() || t <= 0.0 {
        return Err(Error::BadData("time must be finite and > 0"));
    }
    let alpha = diffusivity * std::f64::consts::PI.powi(2) * t / length.powi(2);
    let mut sum = 0.0;
    let mut n = 1_u32;
    loop {
        let term = 2.0 * (-1.0_f64).powi(n as i32) * (-alpha * (n as f64).powi(2)).exp();
        sum += term;
        if term.abs() < 1e-18 {
            break;
        }
        n += 1;
        debug_assert!(n < 10_000_000, "G2 series failed to converge");
    }
    Ok(1.0 + sum)
}

/// Langmuir equilibrium load `c_t = N K c / (1 + K c)` \[mol/m³\] (T2-eq).
///
/// Rejects non-finite/negative inputs.
pub fn equilibrium_trapped(
    site_density: f64,
    equilibrium_constant: f64,
    c_mobile: f64,
) -> Result<f64, Error> {
    if !site_density.is_finite() || site_density < 0.0 {
        return Err(Error::BadData("site density must be finite and >= 0"));
    }
    if !equilibrium_constant.is_finite() || equilibrium_constant < 0.0 {
        return Err(Error::BadData(
            "equilibrium constant must be finite and >= 0",
        ));
    }
    if !c_mobile.is_finite() || c_mobile < 0.0 {
        return Err(Error::BadData(
            "mobile concentration must be finite and >= 0",
        ));
    }
    let kc = equilibrium_constant * c_mobile;
    Ok(site_density * kc / (1.0 + kc))
}

/// Oriani effective diffusivity `D_eff = D / (1 + K N)` \[m²/s\] (G3a):
/// the low-occupancy (`K c ≪ 1`) local-equilibrium limit of (T1–T2).
///
/// Rejects non-finite/negative inputs.
pub fn effective_diffusivity(
    diffusivity: f64,
    equilibrium_constant: f64,
    site_density: f64,
) -> Result<f64, Error> {
    if !diffusivity.is_finite() || diffusivity < 0.0 {
        return Err(Error::BadData("diffusivity must be finite and >= 0"));
    }
    if !equilibrium_constant.is_finite() || equilibrium_constant < 0.0 {
        return Err(Error::BadData(
            "equilibrium constant must be finite and >= 0",
        ));
    }
    if !site_density.is_finite() || site_density < 0.0 {
        return Err(Error::BadData("site density must be finite and >= 0"));
    }
    Ok(diffusivity / (1.0 + equilibrium_constant * site_density))
}

/// Irreversible-trap fill `c_t(t) = N (1 − e^{−k c t})` \[mol/m³\] (G3c):
/// (T2) at `p → 0` with mobile `c` held fixed and `c_t(0) = 0`.
///
/// Rejects non-finite/negative inputs.
pub fn irreversible_fill(
    rate_k: f64,
    c_mobile: f64,
    site_density: f64,
    t: f64,
) -> Result<f64, Error> {
    if !rate_k.is_finite() || rate_k < 0.0 {
        return Err(Error::BadData("trapping rate must be finite and >= 0"));
    }
    if !c_mobile.is_finite() || c_mobile < 0.0 {
        return Err(Error::BadData(
            "mobile concentration must be finite and >= 0",
        ));
    }
    if !site_density.is_finite() || site_density < 0.0 {
        return Err(Error::BadData("site density must be finite and >= 0"));
    }
    if !t.is_finite() || t < 0.0 {
        return Err(Error::BadData("time must be finite and >= 0"));
    }
    Ok(site_density * (1.0 - (-rate_k * c_mobile * t).exp()))
}

/// Sieverts surface concentration `c = K_S sqrt(p)` \[mol/m³\] (G4).
///
/// Rejects non-finite/negative inputs.
pub fn sieverts_concentration(solubility: f64, pressure: f64) -> Result<f64, Error> {
    if !solubility.is_finite() || solubility < 0.0 {
        return Err(Error::BadData("solubility must be finite and >= 0"));
    }
    if !pressure.is_finite() || pressure < 0.0 {
        return Err(Error::BadData("pressure must be finite and >= 0"));
    }
    Ok(solubility * pressure.sqrt())
}

/// Recombination rate `K_r = kr0 * exp(-e_r / R / temp)` \[m⁴/mol/s\] (G5
/// surface Arrhenius helper; same validation as
/// [`crate::params::arrhenius`]).
pub fn recombination_rate_arrhenius(kr0: f64, e_r: f64, temp: f64) -> Result<f64, Error> {
    crate::params::arrhenius(kr0, e_r, temp)
}

// ---------------------------------------------------------------------------
// Discrete operator
// ---------------------------------------------------------------------------

/// Tridiagonal rate system `dc/dt = A c + rhs` with `A` in
/// sub/diag/sup form and face diffusivities for flux evaluation.
pub(crate) struct DiffusionSystem {
    pub(crate) sub: Vec<f64>,
    pub(crate) diag: Vec<f64>,
    pub(crate) sup: Vec<f64>,
    pub(crate) rhs: Vec<f64>,
    pub(crate) face_d: Vec<f64>,
}

/// Face concentration imposed by an equilibrium end; `None` for zero flux.
pub(crate) fn face_concentration(bc: &Boundary) -> Option<f64> {
    bc.surface_concentration()
}

/// Assemble the cell-centred finite-volume diffusion operator.
///
/// Interior faces use the arithmetic-mean diffusivity; boundary faces the
/// adjacent cell value. A Dirichlet/Sieverts/Henry face of value `cf` with
/// conductance `g = 2 D/dx` contributes `−g/dx` to the diagonal and
/// `g cf/dx` to the right-hand side; a zero-flux face contributes nothing.
fn assemble(
    params: &TransportParams,
    left: &Boundary,
    right: &Boundary,
) -> Result<DiffusionSystem, Error> {
    let n = params.cells;
    let dx = params.dx();
    let d_cell = params.diffusivities()?;
    let mut face_d = vec![0.0; n + 1];
    face_d[0] = d_cell[0];
    face_d[n] = d_cell[n - 1];
    for i in 1..n {
        face_d[i] = 0.5 * (d_cell[i - 1] + d_cell[i]);
    }
    let mut sub = vec![0.0; n.saturating_sub(1)];
    let mut diag = vec![0.0; n];
    let mut sup = vec![0.0; n.saturating_sub(1)];
    let mut rhs = vec![0.0; n];
    for i in 0..n {
        // West coupling.
        let (gw, cw) = if i > 0 {
            (face_d[i] / dx, None)
        } else {
            match face_concentration(left) {
                Some(cf) => (2.0 * face_d[0] / dx, Some(cf)),
                None => (0.0, None),
            }
        };
        // East coupling.
        let (ge, ce) = if i + 1 < n {
            (face_d[i + 1] / dx, None)
        } else {
            match face_concentration(right) {
                Some(cf) => (2.0 * face_d[n] / dx, Some(cf)),
                None => (0.0, None),
            }
        };
        diag[i] = -(gw + ge) / dx;
        if i > 0 {
            sub[i - 1] = gw / dx;
        }
        if i + 1 < n {
            sup[i] = ge / dx;
        }
        let mut b = params.source_at(i);
        if let Some(cf) = cw {
            b += gw * cf / dx;
        }
        if let Some(cf) = ce {
            b += ge * cf / dx;
        }
        rhs[i] = b;
    }
    Ok(DiffusionSystem {
        sub,
        diag,
        sup,
        rhs,
        face_d,
    })
}

/// Outward surface flux \[mol/m²/s\] from the boundary-adjacent cell value
/// (`D (c_cell − c_face) / (dx/2)` for equilibrium ends, `0` for zero flux).
pub(crate) fn outward_flux(bc: &Boundary, c_cell: f64, d_face: f64, dx: f64) -> f64 {
    match face_concentration(bc) {
        Some(cf) => d_face * (c_cell - cf) / (dx * 0.5),
        None => 0.0,
    }
}

/// Exact backward-Euler map of (T2) at fixed mobile `c` over `dt`:
///
/// ```text
/// ctⱼⁿᵉʷ = (ctⱼᵒˡᵈ + dt·kⱼ·c·Nⱼ) / (1 + dt·(kⱼ·c + pⱼ)).
/// ```
pub(crate) fn trap_update(ct_old: f64, c: f64, k: f64, p: f64, site_density: f64, dt: f64) -> f64 {
    (ct_old + dt * k * c * site_density) / (1.0 + dt * (k * c + p))
}

pub(crate) fn tridiag_err(e: nucleide_linalg::TridiagError) -> Error {
    Error::Tridiag(e.to_string())
}

/// Cap on the per-step trap-coupling Picard iterations.
pub(crate) const MAX_PICARD: usize = 100;

// ---------------------------------------------------------------------------
// Steady state (G1/G3b/G4)
// ---------------------------------------------------------------------------

/// Cap on the two-face G5 Newton iterations.
pub(crate) const G5_NEWTON_MAX: usize = 50;
/// Relative tolerance on the G5 face residuals.
pub(crate) const G5_RTOL: f64 = 1e-9;
/// Absolute tolerance on the G5 face residuals.
pub(crate) const G5_ATOL: f64 = 1e-12;

/// Steady state with at least one recombination end (G5).
///
/// The mobile profile is affine in the recombination face values (the
/// steady system is linear for frozen faces, and traps evaluate pointwise
/// afterwards), so `1 + n_rec` Thomas solves pin the face-adjacent
/// response exactly: one basis profile with all recombination faces at
/// zero plus one unit-perturbation column per recombination end. A single
/// end then closes in closed form from
/// `K_r cf² + g (1 − b) cf − g a = 0`; two ends close by Newton iteration
/// on the face pair with the analytic Jacobian (quadratic convergence).
/// Fails with [`Error::NotConverged`] only if the face Newton exhausts
/// its cap.
fn steady_recombination(
    params: &TransportParams,
    left: &Boundary,
    right: &Boundary,
    kr_left: Option<f64>,
    kr_right: Option<f64>,
) -> Result<SteadyState, Error> {
    let n = params.cells;
    let dx = params.dx();
    let d_cell = params.diffusivities()?;
    // Recombination ends in (end, rate, conductance) form (`0` = left).
    let ends: Vec<(usize, f64, f64)> = [
        kr_left.map(|kr| (0, kr, 2.0 * d_cell[0] / dx)),
        kr_right.map(|kr| (1, kr, 2.0 * d_cell[n - 1] / dx)),
    ]
    .into_iter()
    .flatten()
    .collect();
    let solve_faces = |face: [f64; 2]| -> Result<Vec<f64>, Error> {
        let bl = kr_left.map_or_else(|| left.clone(), |_| Boundary::Dirichlet(face[0]));
        let br = kr_right.map_or_else(|| right.clone(), |_| Boundary::Dirichlet(face[1]));
        linear_steady(&assemble(params, &bl, &br)?)
    };
    // Basis: faces at zero, then one unit column per recombination end.
    // Response rows hold (a, b0, b1) with c_adj = a + b0*cfL + b1*cfR.
    let adj = |e: usize| if e == 0 { 0 } else { n - 1 };
    let base = solve_faces([0.0, 0.0])?;
    let mut resp = [[0.0_f64; 3]; 2];
    for e in 0..2 {
        resp[e][0] = base[adj(e)];
    }
    for (j, &(e, _, _)) in ends.iter().enumerate() {
        let mut face = [0.0, 0.0];
        face[e] = 1.0;
        let col = solve_faces(face)?;
        for r in 0..2 {
            resp[r][1 + j] = col[adj(r)] - base[adj(r)];
        }
    }
    // Face values from the recombination quadratics.
    let mut face = [0.0_f64; 2];
    match ends.as_slice() {
        [(e, kr, g)] => {
            let (a, b) = (resp[*e][0], resp[*e][1]);
            // K_r cf² + g(1−b) cf − g a = 0, positive root (a ≥ 0, b < 1).
            let lin = g * (1.0 - b);
            let disc = lin.mul_add(lin, 4.0 * kr * g * a.max(0.0));
            face[*e] = (disc.sqrt() - lin) / (2.0 * kr);
        }
        [(e0, kr0, g0), (e1, kr1, g1)] => {
            let (a0, b00, b01) = (resp[*e0][0], resp[*e0][1], resp[*e0][2]);
            let (a1, b10, b11) = (resp[*e1][0], resp[*e1][1], resp[*e1][2]);
            let (mut x0, mut x1) = (0.0_f64, 0.0_f64);
            let mut ok = false;
            for _ in 0..G5_NEWTON_MAX {
                let c0 = a0 + b00 * x0 + b01 * x1;
                let c1 = a1 + b10 * x0 + b11 * x1;
                let f0 = kr0 * x0 * x0 - g0 * (c0 - x0);
                let f1 = kr1 * x1 * x1 - g1 * (c1 - x1);
                if f0.abs() <= G5_ATOL + G5_RTOL * (kr0 * x0 * x0).abs()
                    && f1.abs() <= G5_ATOL + G5_RTOL * (kr1 * x1 * x1).abs()
                {
                    ok = true;
                    break;
                }
                let j00 = 2.0 * kr0 * x0 - g0 * (b00 - 1.0);
                let j01 = -g0 * b01;
                let j10 = -g1 * b10;
                let j11 = 2.0 * kr1 * x1 - g1 * (b11 - 1.0);
                let det = j00 * j11 - j01 * j10;
                if !det.is_finite() || det == 0.0 {
                    break;
                }
                x0 = (x0 - (j11 * f0 - j01 * f1) / det).max(0.0);
                x1 = (x1 - (j00 * f1 - j10 * f0) / det).max(0.0);
            }
            if !ok {
                return Err(Error::NotConverged);
            }
            face[*e0] = x0;
            face[*e1] = x1;
        }
        _ => {
            // Unreachable-in-practice: callers dispatch here only with at
            // least one recombination end. Loud error, never a panic.
            return Err(Error::BadBoundary(
                "steady_recombination needs a recombination end",
            ));
        }
    }
    let mobile = solve_faces(face)?;
    let bl = kr_left.map_or_else(|| left.clone(), |_| Boundary::Dirichlet(face[0]));
    let br = kr_right.map_or_else(|| right.clone(), |_| Boundary::Dirichlet(face[1]));
    let sys = assemble(params, &bl, &br)?;
    let flux_left = kr_left.map_or_else(
        || outward_flux(left, mobile[0], sys.face_d[0], dx),
        |kr| kr * face[0] * face[0],
    );
    let flux_right = kr_right.map_or_else(
        || outward_flux(right, mobile[n - 1], sys.face_d[n], dx),
        |kr| kr * face[1] * face[1],
    );
    finish_steady(params, mobile, flux_left, flux_right)
}

/// Langmuir trapped loads, inventories, and the [`SteadyState`] bundle for
/// a converged mobile profile with known outward fluxes.
fn finish_steady(
    params: &TransportParams,
    mobile: Vec<f64>,
    flux_left: f64,
    flux_right: f64,
) -> Result<SteadyState, Error> {
    let n = params.cells;
    let dx = params.dx();
    let mut trapped = vec![vec![0.0; params.traps.len()]; n];
    for i in 0..n {
        let rates = params.trap_rates_at(i)?;
        for (j, (k, p)) in rates.iter().enumerate() {
            let keq = k / p;
            trapped[i][j] = equilibrium_trapped(params.traps[j].site_density, keq, mobile[i])?;
        }
    }
    let inventory_mobile: f64 = mobile.iter().sum::<f64>() * dx;
    let inventory_trapped: f64 = trapped
        .iter()
        .map(|row| row.iter().sum::<f64>())
        .sum::<f64>()
        * dx;
    Ok(SteadyState {
        centres: params.cell_centres(),
        flux_left,
        flux_right,
        mobile,
        trapped,
        inventory_mobile,
        inventory_trapped,
    })
}

/// Solve one linear steady system `-A c = rhs` through the Thomas path.
pub(crate) fn linear_steady(sys: &DiffusionSystem) -> Result<Vec<f64>, Error> {
    let neg_sub: Vec<f64> = sys.sub.iter().map(|v| -v).collect();
    let neg_diag: Vec<f64> = sys.diag.iter().map(|v| -v).collect();
    let neg_sup: Vec<f64> = sys.sup.iter().map(|v| -v).collect();
    nucleide_linalg::tridiag::solve(&neg_sub, &neg_diag, &neg_sup, &sys.rhs).map_err(tridiag_err)
}

/// Trap-free-style steady state of (T1–T2).
///
/// Solves `A c + b = 0` through [`nucleide_linalg::tridiag`] (exact to
/// roundoff for the G1/G4 linear profiles) and evaluates the Langmuir
/// isotherm (T2-eq) pointwise for the trapped loads (G3b). A recombination
/// end (G5) closes through the exact face-response construction
/// ([`steady_recombination`]); the transient closes them per step (G6,
/// [`solve_recombination`]).
pub fn steady_state(
    params: &TransportParams,
    left: &Boundary,
    right: &Boundary,
) -> Result<SteadyState, Error> {
    let kr_left = left.recombination_rate();
    let kr_right = right.recombination_rate();
    if kr_left.is_some() || kr_right.is_some() {
        return steady_recombination(params, left, right, kr_left, kr_right);
    }
    let sys = assemble(params, left, right)?;
    let mobile = linear_steady(&sys)?;
    let dx = params.dx();
    let n = params.cells;
    let flux_left = outward_flux(left, mobile[0], sys.face_d[0], dx);
    let flux_right = outward_flux(right, mobile[n - 1], sys.face_d[n], dx);
    finish_steady(params, mobile, flux_left, flux_right)
}

// ---------------------------------------------------------------------------
// Transient (G2 + invariants)
// ---------------------------------------------------------------------------

/// Solve the (T1–T2) transient over `grid` from `initial`.
///
/// Between consecutive event times (`t = 0` plus every output time) the
/// stepper takes `ceil(span / dt_max)` theta steps of equal width (each
/// `>= dt_min`, else [`Error::BadOption`]); output rows land exactly on
/// the grid times. Mobile and traps couple by Picard iteration to
/// `rtol`/`atol` (trap-free steps take one solve). A recombination end
/// (G6) closes per step through the affine face-response construction
/// ([`solve_recombination`]), fused into the same Picard loop with the
/// same tolerances; outward fluxes on recombination ends report the
/// converged `K_r cf²`. Exceeding `max_steps` returns
/// [`Error::StepBudget`]; any Picard or face-Newton failure returns
/// [`Error::NotConverged`].
pub fn solve(
    params: &TransportParams,
    left: &Boundary,
    right: &Boundary,
    grid: &TimeGrid,
    initial: &InitialState,
    opts: &SolverOptions,
) -> Result<Solution, Error> {
    opts.validate()?;
    if initial.mobile.len() != params.cells
        || initial.trapped.len() != params.cells
        || initial
            .trapped
            .iter()
            .any(|row| row.len() != params.traps.len())
    {
        return Err(Error::BadState(
            "initial state shape must match params geometry",
        ));
    }
    let kr_left = left.recombination_rate();
    let kr_right = right.recombination_rate();
    if kr_left.is_some() || kr_right.is_some() {
        return solve_recombination(params, left, right, grid, initial, opts, kr_left, kr_right);
    }
    let sys = assemble(params, left, right)?;
    let n = params.cells;
    let ntraps = params.traps.len();
    let dx = params.dx();
    let theta = opts.theta.value();

    let mut c = initial.mobile.clone();
    let mut ct = initial.trapped.clone();

    let mut mobile_out = Vec::with_capacity(grid.times.len());
    let mut trapped_out = Vec::with_capacity(grid.times.len());
    let mut flux_left = Vec::with_capacity(grid.times.len());
    let mut flux_right = Vec::with_capacity(grid.times.len());

    let mut t_prev = 0.0_f64;
    let mut steps = 0_usize;
    // Per-cell trap rates are time-independent (steady T profile): cache.
    let mut rates: Vec<Vec<(f64, f64)>> = Vec::with_capacity(n);
    for i in 0..n {
        rates.push(params.trap_rates_at(i)?);
    }

    for &t_out in &grid.times {
        let span = t_out - t_prev;
        let n_sub = ((span / opts.dt_max).ceil() as usize).max(1);
        let dt = span / n_sub as f64;
        if dt < opts.dt_min {
            return Err(Error::BadOption("output spacing needs a step below dt_min"));
        }
        for _ in 0..n_sub {
            // Theta-step matrices: M = I − dt·θ·A (strictly dominant, so
            // the Thomas pivots never vanish for well-posed inputs).
            let m_sub: Vec<f64> = sys.sub.iter().map(|v| -dt * theta * v).collect();
            let m_diag: Vec<f64> = sys.diag.iter().map(|v| 1.0 - dt * theta * v).collect();
            let m_sup: Vec<f64> = sys.sup.iter().map(|v| -dt * theta * v).collect();
            // Explicit half: e = (I + dt·(1−θ)·A) c + dt·b.
            let mut e = vec![0.0; n];
            for i in 0..n {
                let mut a_c = sys.diag[i] * c[i];
                if i > 0 {
                    a_c += sys.sub[i - 1] * c[i - 1];
                }
                if i + 1 < n {
                    a_c += sys.sup[i] * c[i + 1];
                }
                e[i] = c[i] + dt * (1.0 - theta) * a_c + dt * sys.rhs[i];
            }
            if ntraps == 0 {
                c = nucleide_linalg::tridiag::solve(&m_sub, &m_diag, &m_sup, &e)
                    .map_err(tridiag_err)?;
            } else {
                // Picard: trap sink from the latest trapped iterate, then
                // the exact per-cell backward-Euler trap map at the new
                // mobile iterate, until consecutive iterates agree to
                // rtol/atol (the first pass always runs predict + correct).
                let mut c_iter = c.clone();
                let mut ct_iter = ct.clone();
                let mut converged = false;
                for _ in 0..MAX_PICARD {
                    let mut rhs = e.clone();
                    for i in 0..n {
                        let old_total: f64 = ct[i].iter().sum();
                        let star_total: f64 = ct_iter[i].iter().sum();
                        rhs[i] -= star_total - old_total;
                    }
                    let c_next = nucleide_linalg::tridiag::solve(&m_sub, &m_diag, &m_sup, &rhs)
                        .map_err(tridiag_err)?;
                    let mut ct_next = ct.clone();
                    for i in 0..n {
                        for j in 0..ntraps {
                            let (k, p) = rates[i][j];
                            ct_next[i][j] = trap_update(
                                ct[i][j],
                                c_next[i],
                                k,
                                p,
                                params.traps[j].site_density,
                                dt,
                            );
                        }
                    }
                    let mut err = 0.0_f64;
                    for i in 0..n {
                        let scale = opts.atol + opts.rtol * c_next[i].abs().max(c_iter[i].abs());
                        err = err.max((c_next[i] - c_iter[i]).abs() / scale);
                        for j in 0..ntraps {
                            let scale_t = opts.atol
                                + opts.rtol * ct_next[i][j].abs().max(ct_iter[i][j].abs());
                            err = err.max((ct_next[i][j] - ct_iter[i][j]).abs() / scale_t);
                        }
                    }
                    c_iter = c_next;
                    ct_iter = ct_next;
                    if err <= 1.0 {
                        converged = true;
                        break;
                    }
                }
                if !converged {
                    return Err(Error::NotConverged);
                }
                c = c_iter;
                ct = ct_iter;
            }
            steps += 1;
            if steps > opts.max_steps {
                return Err(Error::StepBudget(opts.max_steps));
            }
        }
        t_prev = t_out;
        mobile_out.push(c.clone());
        trapped_out.push(ct.clone());
        flux_left.push(outward_flux(left, c[0], sys.face_d[0], dx));
        flux_right.push(outward_flux(right, c[n - 1], sys.face_d[n], dx));
    }

    Ok(Solution {
        times: grid.times.clone(),
        mobile: mobile_out,
        trapped: trapped_out,
        flux_left,
        flux_right,
        dx,
        initial: initial.clone(),
    })
}

// ---------------------------------------------------------------------------
// Transient recombination (G6)
// ---------------------------------------------------------------------------

/// Positive root of the pointwise face quadratic `K_r cf² + g cf − g·a = 0`
/// for a frozen adjacent-cell value `a` (the `a` clamp plus the outer
/// clamp is the face clamp keeping `cf ≥ 0`). The rationalized quotient
/// `2·g·a / (√disc + g)` keeps the K_r → 0 limit accurate (the naive
/// `(√disc − g)/(2·K_r)` cancels catastrophically as K_r shrinks).
pub(crate) fn face_from_adjacent(kr: f64, g: f64, adj: f64) -> f64 {
    let a = adj.max(0.0);
    let disc = g.mul_add(g, 4.0 * kr * g * a);
    (2.0 * g * a / (disc.sqrt() + g)).max(0.0)
}

/// Closed form for one affine end `c_adj = a + b·cf`:
/// `K_r cf² + g(1−b) cf − g·a = 0`, positive root (mirrors the G5 single
/// end), rationalized for the same K_r → 0 accuracy as
/// [`face_from_adjacent`].
pub(crate) fn face_closed_single(kr: f64, g: f64, a: f64, b: f64) -> f64 {
    let a = a.max(0.0);
    let lin = g * (1.0 - b);
    let disc = lin.mul_add(lin, 4.0 * kr * g * a);
    (2.0 * g * a / (disc.sqrt() + lin)).max(0.0)
}

/// Newton iteration on the two-face system with the analytic Jacobian
/// (mirrors the G5 pair solve) from `start`, to the shared `opts`
/// tolerances on the flux residuals. Fails with [`Error::NotConverged`]
/// when the cap is exhausted or the Jacobian goes singular/non-finite
/// (hard fail: the step is rejected, no dt retry).
#[allow(clippy::too_many_arguments)]
pub(crate) fn face_newton_pair(
    kr0: f64,
    g0: f64,
    a0: f64,
    b00: f64,
    b01: f64,
    kr1: f64,
    g1: f64,
    a1: f64,
    b10: f64,
    b11: f64,
    start: [f64; 2],
    opts: &SolverOptions,
) -> Result<[f64; 2], Error> {
    let (mut x0, mut x1) = (start[0].max(0.0), start[1].max(0.0));
    for _ in 0..G5_NEWTON_MAX {
        let c0 = a0 + b00 * x0 + b01 * x1;
        let c1 = a1 + b10 * x0 + b11 * x1;
        let f0 = kr0 * x0 * x0 - g0 * (c0 - x0);
        let f1 = kr1 * x1 * x1 - g1 * (c1 - x1);
        if f0.abs() <= opts.atol + opts.rtol * (kr0 * x0 * x0).abs()
            && f1.abs() <= opts.atol + opts.rtol * (kr1 * x1 * x1).abs()
        {
            return Ok([x0, x1]);
        }
        let j00 = 2.0 * kr0 * x0 - g0 * (b00 - 1.0);
        let j01 = -g0 * b01;
        let j10 = -g1 * b10;
        let j11 = 2.0 * kr1 * x1 - g1 * (b11 - 1.0);
        let det = j00 * j11 - j01 * j10;
        if !det.is_finite() || det == 0.0 {
            break;
        }
        x0 = (x0 - (j11 * f0 - j01 * f1) / det).max(0.0);
        x1 = (x1 - (j00 * f1 - j10 * f0) / det).max(0.0);
    }
    Err(Error::NotConverged)
}

/// Close the recombination faces for a base profile `base` (the mobile
/// solve with zeroed faces) given the sensitivity columns; `start`
/// warm-starts the pair Newton. Returns the closed faces (`None` on a
/// non-recombination end).
#[allow(clippy::too_many_arguments)]
pub(crate) fn close_faces(
    base: &[f64],
    sens_left: Option<&Vec<f64>>,
    sens_right: Option<&Vec<f64>>,
    kr_left: Option<f64>,
    kr_right: Option<f64>,
    g_left: f64,
    g_right: f64,
    start: [f64; 2],
    opts: &SolverOptions,
) -> Result<[Option<f64>; 2], Error> {
    match (kr_left, kr_right) {
        (Some(kr), None) => {
            let p = sens_left.expect("left sensitivity built for a left end");
            Ok([Some(face_closed_single(kr, g_left, base[0], p[0])), None])
        }
        (None, Some(kr)) => {
            let p = sens_right.expect("right sensitivity built for a right end");
            let n = base.len();
            Ok([
                None,
                Some(face_closed_single(kr, g_right, base[n - 1], p[n - 1])),
            ])
        }
        (Some(kr0), Some(kr1)) => {
            let pl = sens_left.expect("left sensitivity built for a left end");
            let pr = sens_right.expect("right sensitivity built for a right end");
            let n = base.len();
            let pair = face_newton_pair(
                kr0,
                g_left,
                base[0],
                pl[0],
                pr[0],
                kr1,
                g_right,
                base[n - 1],
                pl[n - 1],
                pr[n - 1],
                start,
                opts,
            )?;
            Ok([Some(pair[0]), Some(pair[1])])
        }
        (None, None) => {
            // Unreachable-in-practice: callers dispatch here only with at
            // least one recombination end. Loud error, never a panic.
            Err(Error::BadBoundary("close_faces needs a recombination end"))
        }
    }
}

/// Mobile profile for closed faces: `base + P·f` (no extra Thomas solve).
pub(crate) fn apply_faces(
    mut base: Vec<f64>,
    sens_left: Option<&Vec<f64>>,
    sens_right: Option<&Vec<f64>>,
    faces: [Option<f64>; 2],
) -> Vec<f64> {
    if let (Some(p), Some(f)) = (sens_left, faces[0]) {
        for (row, col) in base.iter_mut().zip(p.iter()) {
            *row += col * f;
        }
    }
    if let (Some(p), Some(f)) = (sens_right, faces[1]) {
        for (row, col) in base.iter_mut().zip(p.iter()) {
            *row += col * f;
        }
    }
    base
}

/// Transient with at least one recombination end (G6).
///
/// There is no closed form for the recombination transient (the face
/// quadratic couples to the time derivative), so the per-step face system
/// closes numerically: for frozen traps the θ-step is linear in the
/// mobile profile with the faces as sources (`M c = r + Q f`), hence
/// affine (`c = base + P f`) with the sensitivity columns `P = M⁻¹ Q`
/// built once per `dt` span (`n_rec` Thomas solves; `M` only depends on
/// `dt`). Trap-free steps are then exact in `1 + n_rec` solves; steps
/// with traps fuse the face iterate into the trap Picard loop (base solve
/// → face close → trap map per pass, one base solve per pass) to the
/// shared `rtol`/`atol`. Outward fluxes on recombination ends report the
/// converged `K_r cf²`, against which the discrete mass balance closes to
/// roundoff. Any face-Newton or Picard failure is [`Error::NotConverged`].
#[allow(clippy::too_many_arguments)]
fn solve_recombination(
    params: &TransportParams,
    left: &Boundary,
    right: &Boundary,
    grid: &TimeGrid,
    initial: &InitialState,
    opts: &SolverOptions,
    kr_left: Option<f64>,
    kr_right: Option<f64>,
) -> Result<Solution, Error> {
    let n = params.cells;
    let ntraps = params.traps.len();
    let dx = params.dx();
    let theta = opts.theta.value();
    let d_cell = params.diffusivities()?;
    let g_left = 2.0 * d_cell[0] / dx;
    let g_right = 2.0 * d_cell[n - 1] / dx;
    // Base operator: recombination ends as Dirichlet(0), keeping the
    // conductance in the operator while the face source drops out (it
    // re-enters through the sensitivity columns and the explicit old-face
    // term below).
    let bl0 = kr_left.map_or_else(|| left.clone(), |_| Boundary::Dirichlet(0.0));
    let br0 = kr_right.map_or_else(|| right.clone(), |_| Boundary::Dirichlet(0.0));
    let sys = assemble(params, &bl0, &br0)?;

    // Per-cell trap rates are time-independent (steady T profile): cache.
    let mut rates: Vec<Vec<(f64, f64)>> = Vec::with_capacity(n);
    for i in 0..n {
        rates.push(params.trap_rates_at(i)?);
    }

    let mut c = initial.mobile.clone();
    let mut ct = initial.trapped.clone();
    let mut cf_left = kr_left.map(|kr| face_from_adjacent(kr, g_left, c[0]));
    let mut cf_right = kr_right.map(|kr| face_from_adjacent(kr, g_right, c[n - 1]));

    let mut mobile_out = Vec::with_capacity(grid.times.len());
    let mut trapped_out = Vec::with_capacity(grid.times.len());
    let mut flux_left = Vec::with_capacity(grid.times.len());
    let mut flux_right = Vec::with_capacity(grid.times.len());

    let mut t_prev = 0.0_f64;
    let mut steps = 0_usize;

    for &t_out in &grid.times {
        let span = t_out - t_prev;
        let n_sub = ((span / opts.dt_max).ceil() as usize).max(1);
        let dt = span / n_sub as f64;
        if dt < opts.dt_min {
            return Err(Error::BadOption("output spacing needs a step below dt_min"));
        }
        // Theta-step operator (fixed across the span) plus one sensitivity
        // column per recombination end.
        let m_sub: Vec<f64> = sys.sub.iter().map(|v| -dt * theta * v).collect();
        let m_diag: Vec<f64> = sys.diag.iter().map(|v| 1.0 - dt * theta * v).collect();
        let m_sup: Vec<f64> = sys.sup.iter().map(|v| -dt * theta * v).collect();
        let sens_left = kr_left
            .map(|_| {
                let mut v = vec![0.0; n];
                v[0] = dt * theta * g_left / dx;
                nucleide_linalg::tridiag::solve(&m_sub, &m_diag, &m_sup, &v).map_err(tridiag_err)
            })
            .transpose()?;
        let sens_right = kr_right
            .map(|_| {
                let mut v = vec![0.0; n];
                v[n - 1] = dt * theta * g_right / dx;
                nucleide_linalg::tridiag::solve(&m_sub, &m_diag, &m_sup, &v).map_err(tridiag_err)
            })
            .transpose()?;
        for _ in 0..n_sub {
            // Explicit half on the frozen base system plus the old-face
            // sources.
            let mut e = vec![0.0; n];
            for i in 0..n {
                let mut a_c = sys.diag[i] * c[i];
                if i > 0 {
                    a_c += sys.sub[i - 1] * c[i - 1];
                }
                if i + 1 < n {
                    a_c += sys.sup[i] * c[i + 1];
                }
                e[i] = c[i] + dt * (1.0 - theta) * a_c + dt * sys.rhs[i];
            }
            if let Some(cf) = cf_left {
                e[0] += dt * (1.0 - theta) * g_left / dx * cf;
            }
            if let Some(cf) = cf_right {
                e[n - 1] += dt * (1.0 - theta) * g_right / dx * cf;
            }
            if ntraps == 0 {
                let base = nucleide_linalg::tridiag::solve(&m_sub, &m_diag, &m_sup, &e)
                    .map_err(tridiag_err)?;
                let faces = close_faces(
                    &base,
                    sens_left.as_ref(),
                    sens_right.as_ref(),
                    kr_left,
                    kr_right,
                    g_left,
                    g_right,
                    [cf_left.unwrap_or(0.0), cf_right.unwrap_or(0.0)],
                    opts,
                )?;
                cf_left = faces[0];
                cf_right = faces[1];
                c = apply_faces(base, sens_left.as_ref(), sens_right.as_ref(), faces);
            } else {
                // Fused Picard: the face iterate is re-closed from every
                // base solve and the trap map runs at the new mobile
                // iterate; consecutive (mobile, trap, face) iterates agree
                // to rtol/atol (the first pass always runs).
                let mut c_iter = c.clone();
                let mut ct_iter = ct.clone();
                let mut fl_iter = cf_left;
                let mut fr_iter = cf_right;
                let mut converged = false;
                for _ in 0..MAX_PICARD {
                    let mut rhs = e.clone();
                    for i in 0..n {
                        let old_total: f64 = ct[i].iter().sum();
                        let star_total: f64 = ct_iter[i].iter().sum();
                        rhs[i] -= star_total - old_total;
                    }
                    let base = nucleide_linalg::tridiag::solve(&m_sub, &m_diag, &m_sup, &rhs)
                        .map_err(tridiag_err)?;
                    let faces = close_faces(
                        &base,
                        sens_left.as_ref(),
                        sens_right.as_ref(),
                        kr_left,
                        kr_right,
                        g_left,
                        g_right,
                        [fl_iter.unwrap_or(0.0), fr_iter.unwrap_or(0.0)],
                        opts,
                    )?;
                    let c_next = apply_faces(base, sens_left.as_ref(), sens_right.as_ref(), faces);
                    let mut ct_next = ct.clone();
                    for i in 0..n {
                        for j in 0..ntraps {
                            let (k, p) = rates[i][j];
                            ct_next[i][j] = trap_update(
                                ct[i][j],
                                c_next[i],
                                k,
                                p,
                                params.traps[j].site_density,
                                dt,
                            );
                        }
                    }
                    let mut err = 0.0_f64;
                    for i in 0..n {
                        let scale = opts.atol + opts.rtol * c_next[i].abs().max(c_iter[i].abs());
                        err = err.max((c_next[i] - c_iter[i]).abs() / scale);
                        for j in 0..ntraps {
                            let scale_t = opts.atol
                                + opts.rtol * ct_next[i][j].abs().max(ct_iter[i][j].abs());
                            err = err.max((ct_next[i][j] - ct_iter[i][j]).abs() / scale_t);
                        }
                    }
                    if let (Some(nf), Some(prev)) = (faces[0], fl_iter) {
                        let scale = opts.atol + opts.rtol * nf.abs().max(prev.abs());
                        err = err.max((nf - prev).abs() / scale);
                    }
                    if let (Some(nf), Some(prev)) = (faces[1], fr_iter) {
                        let scale = opts.atol + opts.rtol * nf.abs().max(prev.abs());
                        err = err.max((nf - prev).abs() / scale);
                    }
                    c_iter = c_next;
                    ct_iter = ct_next;
                    fl_iter = faces[0];
                    fr_iter = faces[1];
                    if err <= 1.0 {
                        converged = true;
                        break;
                    }
                }
                if !converged {
                    return Err(Error::NotConverged);
                }
                c = c_iter;
                ct = ct_iter;
                cf_left = fl_iter;
                cf_right = fr_iter;
            }
            steps += 1;
            if steps > opts.max_steps {
                return Err(Error::StepBudget(opts.max_steps));
            }
        }
        t_prev = t_out;
        mobile_out.push(c.clone());
        trapped_out.push(ct.clone());
        flux_left.push(cf_left.map_or_else(
            || outward_flux(left, c[0], sys.face_d[0], dx),
            |cf| kr_left.expect("rate present for a closed left face") * cf * cf,
        ));
        flux_right.push(cf_right.map_or_else(
            || outward_flux(right, c[n - 1], sys.face_d[n], dx),
            |cf| kr_right.expect("rate present for a closed right face") * cf * cf,
        ));
    }

    Ok(Solution {
        times: grid.times.clone(),
        mobile: mobile_out,
        trapped: trapped_out,
        flux_left,
        flux_right,
        dx,
        initial: initial.clone(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::tests::no_traps;

    fn tight(theta: Theta) -> SolverOptions {
        SolverOptions {
            theta,
            rtol: 1e-10,
            atol: 1e-14,
            dt_min: 1e-14,
            dt_max: 0.1,
            max_steps: 10_000_000,
        }
    }

    #[test]
    fn helpers_match_closed_forms() {
        // G2-lag: L^2/6D for the G1/G2 slab.
        assert!((time_lag(1e-3, 1e-9).unwrap() - 166.666_666_666_666_66).abs() < 1e-9);
        // G2 series endpoints: pre-breakthrough ~0, late ~1.
        assert!(breakthrough_ratio(1e-9, 1e-3, 1.0).unwrap().abs() < 1e-9);
        assert!(
            (breakthrough_ratio(1e-9, 1e-3, 1000.0).unwrap() - 0.999_896_553_627_592_4).abs()
                < 1e-12
        );
        // G3a: D_eff = D/4 at KN = 3.
        assert_eq!(effective_diffusivity(1e-9, 1.0, 3.0).unwrap(), 2.5e-10);
        // T2-eq: N K c/(1 + K c).
        assert!((equilibrium_trapped(2.0, 1.0, 0.5).unwrap() - 2.0 / 3.0).abs() < 1e-15);
        // G3c: N(1 − e^{−kct}) at t = 0 and late.
        assert_eq!(irreversible_fill(0.05, 1.0, 2.0, 0.0).unwrap(), 0.0);
        assert!((irreversible_fill(0.05, 1.0, 2.0, 1000.0).unwrap() - 2.0).abs() < 1e-12);
        // G4: K_S sqrt(p).
        assert_eq!(sieverts_concentration(2.0, 16.0).unwrap(), 8.0);
        assert!(time_lag(0.0, 1e-9).is_err());
        assert!(time_lag(1e-3, -1.0).is_err());
        assert!(breakthrough_ratio(0.0, 1e-3, 1.0).is_err());
        assert!(breakthrough_ratio(1e-9, 0.0, 1.0).is_err());
        assert!(breakthrough_ratio(1e-9, 1e-3, 0.0).is_err());
        assert!(breakthrough_ratio(1e-9, 1e-3, -1.0).is_err());
        assert!(effective_diffusivity(-1.0, 1.0, 3.0).is_err());
        assert!(effective_diffusivity(1e-9, -1.0, 3.0).is_err());
        assert!(effective_diffusivity(1e-9, 1.0, -3.0).is_err());
        assert!(equilibrium_trapped(-1.0, 1.0, 0.5).is_err());
        assert!(equilibrium_trapped(2.0, -1.0, 0.5).is_err());
        assert!(equilibrium_trapped(2.0, 1.0, -0.5).is_err());
        assert!(irreversible_fill(0.05, 1.0, 2.0, -1.0).is_err());
        assert!(irreversible_fill(-0.05, 1.0, 2.0, 1.0).is_err());
        assert!(irreversible_fill(0.05, -1.0, 2.0, 1.0).is_err());
        assert!(irreversible_fill(0.05, 1.0, -2.0, 1.0).is_err());
        assert!(sieverts_concentration(-1.0, 16.0).is_err());
        assert!(sieverts_concentration(2.0, -16.0).is_err());
    }

    #[test]
    fn steady_dirichlet_is_linear() {
        // G1 smoke: coarse grid recovers the linear profile to roundoff.
        let p = no_traps(1e-3, 8, 1e-9);
        let s = steady_state(
            &p,
            &Boundary::dirichlet(1.0).unwrap(),
            &Boundary::dirichlet(0.0).unwrap(),
        )
        .unwrap();
        for (i, c) in s.mobile.iter().enumerate() {
            let x = (i as f64 + 0.5) / 8.0;
            assert!((c - (1.0 - x)).abs() < 1e-12, "cell {i}: {c}");
        }
        assert!((s.flux_left + 1e-6).abs() < 1e-15);
        assert!((s.flux_right - 1e-6).abs() < 1e-15);
        assert!((s.inventory_mobile - 5e-4).abs() < 1e-15);
    }

    #[test]
    fn transient_holds_steady_state() {
        // Invariant: starting on the G1 profile with matching ends stays put.
        let p = no_traps(1e-3, 16, 1e-9);
        let s = steady_state(
            &p,
            &Boundary::dirichlet(1.0).unwrap(),
            &Boundary::dirichlet(0.0).unwrap(),
        )
        .unwrap();
        let init = InitialState::new(&p, s.mobile.clone(), vec![vec![]; 16]).unwrap();
        let sol = solve(
            &p,
            &Boundary::dirichlet(1.0).unwrap(),
            &Boundary::dirichlet(0.0).unwrap(),
            &TimeGrid::new(vec![10.0, 100.0]).unwrap(),
            &init,
            &tight(Theta::CrankNicolson),
        )
        .unwrap();
        for row in &sol.mobile {
            for (got, want) in row.iter().zip(&s.mobile) {
                assert!((got - want).abs() / want.max(1e-300) < 1e-9);
            }
        }
    }

    #[test]
    fn pure_neumann_steady_has_no_solution() {
        // Zero-flux ends with no source: −A is singular, so the Thomas
        // path surfaces its elimination error instead of a profile.
        let p = no_traps(1e-3, 8, 1e-9);
        assert!(steady_state(&p, &Boundary::ZeroFlux, &Boundary::ZeroFlux).is_err());
    }

    #[test]
    fn dirichlet_zeroflux_is_flat() {
        // ZeroFlux right end with uniform Dirichlet left: flat profile.
        let p = no_traps(1e-3, 8, 1e-9);
        let s = steady_state(&p, &Boundary::dirichlet(2.0).unwrap(), &Boundary::ZeroFlux).unwrap();
        for c in &s.mobile {
            assert!((c - 2.0).abs() < 1e-12);
        }
        assert!(s.flux_left.abs() < 1e-15);
        assert!(s.flux_right.abs() < 1e-15);
    }

    #[test]
    fn recombination_boundary_constructors() {
        // The Arrhenius constructor flows through the same validation.
        assert!(Boundary::recombination_arrhenius(1.0, 0.0, 500.0).is_ok());
        assert_eq!(
            Boundary::recombination_arrhenius(1.0, 0.0, 500.0).unwrap(),
            Boundary::recombination(1.0).unwrap()
        );
        assert!(Boundary::recombination_arrhenius(-1.0, 0.0, 500.0).is_err());
        assert!(Boundary::recombination_arrhenius(1.0, 0.0, 0.0).is_err());
    }

    #[test]
    fn g5a_dirichlet_recombination_closed_form() {
        // G5a: trap-free, S=0, c(0)=c0, J(L)=K_r c_s². D=1e-9, L=1e-3,
        // c0=1.0, K_r=1e-6 gives c_s=(√5−1)/2≈0.6180339887498949,
        // J=K_r c_s²≈3.819660112501052e-7, I=(c0+c_s)L/2≈8.090169943744474e-4.
        let p = no_traps(1e-3, 512, 1e-9);
        let s = steady_state(
            &p,
            &Boundary::dirichlet(1.0).unwrap(),
            &Boundary::recombination(1e-6).unwrap(),
        )
        .unwrap();
        let cs = 0.618_033_988_749_894_9;
        for (i, c) in s.mobile.iter().enumerate() {
            let x = (i as f64 + 0.5) / 512.0;
            assert!((c - (1.0 + (cs - 1.0) * x)).abs() < 1e-12, "cell {i}: {c}");
        }
        assert!((s.flux_right - 3.819_660_112_501_052e-7).abs() < 1e-12 * 3.82e-7 + 1e-18);
        assert!((s.flux_left + s.flux_right).abs() < 1e-12);
        assert!((s.inventory_mobile - 8.090_169_943_744_474e-4).abs() < 1e-12);
        assert!(s.mobile.iter().all(|&c| c >= 0.0));
    }

    #[test]
    fn g5b_large_rate_recovers_dirichlet() {
        // G5b: K_r→∞ recovers the G1 Dirichlet end (c_s→0, J→D c0/L).
        let p = no_traps(1e-3, 64, 1e-9);
        let s = steady_state(
            &p,
            &Boundary::dirichlet(1.0).unwrap(),
            &Boundary::recombination(1e12).unwrap(),
        )
        .unwrap();
        assert!((s.flux_right - 1e-6).abs() / 1e-6 < 1e-6);
        assert!((s.inventory_mobile - 5e-4).abs() / 5e-4 < 1e-6);
    }

    /// G5a reference numbers shared by the G6a gates (D=1e-9, L=1e-3,
    /// c0=1.0, K_r=1e-6): c_s=(√5−1)/2, J=K_r c_s².
    const G6_CS: f64 = 0.618_033_988_749_894_9;
    const G6_JSS: f64 = 3.819_660_112_501_052e-7;

    #[test]
    fn g6a_late_time_asymptote_to_g5a() {
        // G6a: the Dirichlet + recombination transient lands on the G5a
        // steady flux within 1e-6 relative at late time (same grid, so the
        // spatial error cancels and only the time asymptote is gated). The
        // approach rate is set by the slowest mode of the problem
        // linearized about the steady face value c_s: λ₁ = D y₁²/L² with
        // y₁ the first positive root of y·cos y + (2·K_r·c_s·L/D)·sin y = 0
        // (here 2·K_r·c_s·L/D = 1.236, y₁ ≈ 2.10, τ₁ = 1/λ₁ ≈ 1.4·t_lag),
        // so t = 60·t_lag ≈ 44·τ₁ leaves e⁻⁴⁴ ≈ 8e-20 of the slowest mode
        // and the 1e-6 gate is a pure asymptote check. The steps are kept
        // resolved (dt_max = 1) because the t = 0 corner kink excites
        // Crank–Nicolson fine-mode ringing that only damps as |g| → 1 on
        // modes with λ·dt ≫ 2 and would otherwise still sit in the two
        // face-adjacent cells (visible in flux_left) at this time.
        let lag = time_lag(1e-3, 1e-9).unwrap();
        let p = no_traps(1e-3, 128, 1e-9);
        let left = Boundary::dirichlet(1.0).unwrap();
        let right = Boundary::recombination(1e-6).unwrap();
        let steady = steady_state(&p, &left, &right).unwrap();
        assert!((steady.flux_right - G6_JSS).abs() / G6_JSS < 1e-9);
        let sol = solve(
            &p,
            &left,
            &right,
            &TimeGrid::new(vec![60.0 * lag]).unwrap(),
            &InitialState::zeros(&p),
            &SolverOptions {
                dt_max: 1.0,
                rtol: 1e-10,
                atol: 1e-14,
                ..Default::default()
            },
        )
        .unwrap();
        assert!((sol.flux_right[0] - G6_JSS).abs() / G6_JSS < 1e-6);
        assert!((sol.flux_left[0] + sol.flux_right[0]).abs() / G6_JSS < 1e-6);
        for (i, &c) in sol.mobile[0].iter().enumerate() {
            let x = (i as f64 + 0.5) / 128.0;
            assert!(
                (c - (1.0 + (G6_CS - 1.0) * x)).abs() < 1e-9,
                "cell {i}: {c}"
            );
        }
        assert!(sol.mobile[0].iter().all(|&c| c >= 0.0));
        assert!(sol.flux_right[0] >= 0.0);
    }

    #[test]
    fn g6b_large_rate_recovers_dirichlet_transient() {
        // G6b: K_r→∞ recovers the G2 Dirichlet transient (continuity in
        // 1/K_r, mirroring G5b; K_r→0 is the zero-flux limit instead).
        let p = no_traps(1e-3, 64, 1e-9);
        let left = Boundary::dirichlet(1.0).unwrap();
        let times = vec![100.0, 316.227_766_016_837_96, 1000.0];
        let opts = SolverOptions {
            dt_max: 1.0,
            rtol: 1e-10,
            atol: 1e-14,
            ..Default::default()
        };
        let dir = solve(
            &p,
            &left,
            &Boundary::dirichlet(0.0).unwrap(),
            &TimeGrid::new(times.clone()).unwrap(),
            &InitialState::zeros(&p),
            &opts,
        )
        .unwrap();
        let rec = solve(
            &p,
            &left,
            &Boundary::recombination(1e12).unwrap(),
            &TimeGrid::new(times).unwrap(),
            &InitialState::zeros(&p),
            &opts,
        )
        .unwrap();
        for (d, r) in dir.flux_right.iter().zip(&rec.flux_right) {
            assert!((r - d).abs() <= 1e-6 * d.abs().max(1e-12), "{r} vs {d}");
        }
        for row in &rec.mobile {
            assert!(row.iter().all(|&c| c >= 0.0));
        }
    }

    #[test]
    fn g6c_mass_balance_with_kr_flux() {
        // G6c: trap-free Dirichlet + recombination, one step per output
        // interval, so the discrete theta balance
        // Iᵏ−Iᵏ⁻¹ = dt·[θRᵏ + (1−θ)Rᵏ⁻¹] with R = −(F_L+F_R) closes against
        // the reported K_r cf² flux to roundoff. Initial faces are exact:
        // zero mobile gives cf = 0 on the recombination end, and the
        // t = 0 Dirichlet inflow is 2D(c_cell − cf)/dx by hand.
        let p = no_traps(1e-3, 16, 1e-9);
        let dx = p.dx();
        let left = Boundary::dirichlet(1.0).unwrap();
        let right = Boundary::recombination(1e-6).unwrap();
        let opts = SolverOptions {
            dt_max: 20.0,
            rtol: 1e-10,
            atol: 1e-14,
            ..Default::default()
        };
        let sol = solve(
            &p,
            &left,
            &right,
            &TimeGrid::new(vec![20.0, 40.0]).unwrap(),
            &InitialState::zeros(&p),
            &opts,
        )
        .unwrap();
        let inv: Vec<f64> = sol
            .mobile
            .iter()
            .map(|row| row.iter().sum::<f64>() * dx)
            .collect();
        // t = 0 rates: recombination face at zero, Dirichlet inflow exact.
        let (mut r_prev, mut i_prev) = (2.0 * 1e-9 * 1.0 / dx, 0.0_f64);
        for (k, dt) in [20.0, 20.0].iter().enumerate() {
            let r_now = -(sol.flux_left[k] + sol.flux_right[k]);
            let want = i_prev + dt * (0.5 * r_prev + 0.5 * r_now);
            assert!(
                (inv[k] - want).abs() < 1e-12 * want.abs().max(1e-300),
                "row {k}: {} vs {want}",
                inv[k]
            );
            r_prev = r_now;
            i_prev = inv[k];
        }
    }

    /// Mobile profile of the G6d case at `t_end` with uniform steps of
    /// `dt_max` (spans divide evenly, so `dt == dt_max`).
    fn g6d_profile(theta: Theta, dt_max: f64, t_end: f64) -> Vec<f64> {
        let p = no_traps(1e-3, 64, 1e-9);
        solve(
            &p,
            &Boundary::dirichlet(1.0).unwrap(),
            &Boundary::recombination(1e-6).unwrap(),
            &TimeGrid::new(vec![t_end]).unwrap(),
            &InitialState::zeros(&p),
            &SolverOptions {
                theta,
                dt_max,
                rtol: 1e-12,
                atol: 1e-15,
                ..Default::default()
            },
        )
        .unwrap()
        .mobile
        .pop()
        .expect("one output row")
    }

    fn max_diff(a: &[f64], b: &[f64]) -> f64 {
        a.iter()
            .zip(b.iter())
            .map(|(x, y)| (x - y).abs())
            .fold(0.0_f64, f64::max)
    }

    #[test]
    fn g6d_dt_halving_order() {
        // G6d: dt-halving on the fixed recombination case shows the
        // θ-method order (2 for Crank–Nicolson, 1 for backward Euler) on
        // two successive halvings within a pinned band. The t = 0 corner
        // kink (zero initial state against c(0) = 1) excites
        // Crank–Nicolson fine-mode ringing that decays as |g| → 1 on
        // modes with λ·dt ≫ 2, so for coarse dt the profile error at late
        // time is ringing, not truncation error; the order is therefore
        // pinned in the resolved band (t = 1200 s ≈ 7·t_lag, dt ≤ 3 s),
        // where both measured orders sit at ≈2.0 (CN) and ≈1.0 (BE).
        let t_end = 1200.0;
        let dt = [3.0_f64, 1.5, 0.75, 0.375];
        let orders = |theta: Theta| {
            let p0 = g6d_profile(theta, dt[0], t_end);
            let p1 = g6d_profile(theta, dt[1], t_end);
            let p2 = g6d_profile(theta, dt[2], t_end);
            let p3 = g6d_profile(theta, dt[3], t_end);
            [
                (max_diff(&p0, &p1) / max_diff(&p1, &p2)).log2(),
                (max_diff(&p1, &p2) / max_diff(&p2, &p3)).log2(),
            ]
        };
        for order in orders(Theta::CrankNicolson) {
            assert!((1.5..=2.5).contains(&order), "Crank–Nicolson order {order}");
        }
        for order in orders(Theta::BackwardEuler) {
            assert!((0.7..=1.3).contains(&order), "backward Euler order {order}");
        }
    }

    // ------------------------------------------------------------------
    // G6e: independent method-of-lines cross-check
    // ------------------------------------------------------------------

    /// Right-hand side of the independent MoL cross-check solver (G6e):
    /// node-centred second-order central differences on a packed state
    /// `y = [c(0..=n), ct blocks]` — mobile concentrations at nodes
    /// `x_i = i·h`, `i = 0..=n`, first (stride-1 Laplacian), then `nt`
    /// McNabb–Foster loads per node. Deliberately disjoint from the
    /// production path (cell-centred FV, implicit θ-stepping, Thomas
    /// solves): the recombination end closes through a ghost node — `-D
    /// ∂c/∂x = K_r c²` differenced centrally at `x = L` gives `c[n+1] =
    /// c[n−1] − 2h K_r c[n]²/D`, an O(h²) Robin closure consistent with the
    /// O(h²) central Laplacian (exact on the G5a linear steady profile) —
    /// traps evolve pointwise by `dct/dt = k c (N − ct) − p ct` with the
    /// mobile sink `−Σ dct/dt`, and the Dirichlet node `i = 0` is pinned
    /// (`dy[0] = 0`; its traps still equilibrate against `c0`).
    #[allow(clippy::too_many_arguments)]
    fn mol_rhs(
        y: &[f64],
        dy: &mut [f64],
        n: usize,
        nt: usize,
        d: f64,
        h: f64,
        kr: f64,
        c0: f64,
        traps: &[(f64, f64, f64)],
    ) {
        let r = 1.0 / (h * h);
        for i in 1..n {
            dy[i] = d * (y[i - 1] - 2.0 * y[i] + y[i + 1]) * r;
        }
        let ghost = y[n - 1] - 2.0 * h * kr * y[n] * y[n] / d;
        dy[n] = d * (y[n - 1] - 2.0 * y[n] + ghost) * r;
        dy[0] = 0.0;
        if !traps.is_empty() {
            let m = n + 1;
            for i in 0..=n {
                let c = if i == 0 { c0 } else { y[i] };
                let mut sink = 0.0;
                for (j, &(k, p, nj)) in traps.iter().enumerate() {
                    let idx = m + i * nt + j;
                    let dct = k * c * (nj - y[idx]) - p * y[idx];
                    dy[idx] = dct;
                    sink += dct;
                }
                if i > 0 {
                    dy[i] -= sink;
                }
            }
        }
    }

    /// One classical RK4 step of `dt` on the packed MoL state (`k`, `tmp`
    /// are scratch vectors of the state length).
    #[allow(clippy::too_many_arguments)]
    fn mol_rk4(
        y: &mut [f64],
        k: &mut [Vec<f64>; 4],
        tmp: &mut [f64],
        dt: f64,
        n: usize,
        nt: usize,
        d: f64,
        h: f64,
        kr: f64,
        c0: f64,
        traps: &[(f64, f64, f64)],
    ) {
        let m = y.len();
        mol_rhs(y, &mut k[0], n, nt, d, h, kr, c0, traps);
        for i in 0..m {
            tmp[i] = y[i] + 0.5 * dt * k[0][i];
        }
        mol_rhs(tmp, &mut k[1], n, nt, d, h, kr, c0, traps);
        for i in 0..m {
            tmp[i] = y[i] + 0.5 * dt * k[1][i];
        }
        mol_rhs(tmp, &mut k[2], n, nt, d, h, kr, c0, traps);
        for i in 0..m {
            tmp[i] = y[i] + dt * k[2][i];
        }
        mol_rhs(tmp, &mut k[3], n, nt, d, h, kr, c0, traps);
        for i in 0..m {
            y[i] += dt / 6.0 * (k[0][i] + 2.0 * (k[1][i] + k[2][i]) + k[3][i]);
        }
    }

    /// G6e MoL reference rows at `t_out` for the Dirichlet(1) +
    /// recombination(`kr`) slab from a clean start, optionally carrying
    /// trap species `(k, p, N)`, on `n` intervals at `dt = safety·h²/D`
    /// (each step lands exactly on the output time).
    #[allow(clippy::too_many_arguments)]
    fn mol_trajectory(
        length: f64,
        n: usize,
        d: f64,
        kr: f64,
        traps: &[(f64, f64, f64)],
        t_out: &[f64],
        safety: f64,
    ) -> Vec<Vec<f64>> {
        let h = length / n as f64;
        let dt = safety * h * h / d;
        let mut y = vec![0.0; (n + 1) * (1 + traps.len())];
        y[0] = 1.0;
        let mut k = [
            vec![0.0; y.len()],
            vec![0.0; y.len()],
            vec![0.0; y.len()],
            vec![0.0; y.len()],
        ];
        let mut tmp = vec![0.0; y.len()];
        let mut rows = Vec::with_capacity(t_out.len());
        let mut t = 0.0;
        for &target in t_out {
            while t < target {
                let step = dt.min(target - t);
                mol_rk4(
                    &mut y,
                    &mut k,
                    &mut tmp,
                    step,
                    n,
                    traps.len(),
                    d,
                    h,
                    kr,
                    1.0,
                    traps,
                );
                t += step;
            }
            rows.push(y.clone());
        }
        rows
    }

    /// Linear interpolation of the mobile profile (`trap = None`) or one
    /// trap block at `x` between nodes.
    fn mol_sample(row: &[f64], n: usize, nt: usize, trap: Option<usize>, h: f64, x: f64) -> f64 {
        let (base, stride) = match trap {
            None => (0, 1),
            Some(j) => (n + 1 + j, nt),
        };
        let u = x / h;
        let j = u.floor() as usize;
        let w = u - j as f64;
        row[base + j * stride] * (1.0 - w) + row[base + (j + 1) * stride] * w
    }

    /// Relative max-norm of `a − b` against `max |a|`.
    fn rel_max_diff(a: &[f64], b: &[f64]) -> f64 {
        let scale = a.iter().fold(0.0_f64, |m, &v| m.max(v.abs())).max(1e-300);
        max_diff(a, b) / scale
    }

    /// G6e production rows of the Dirichlet(1) + recombination(1e-6) case
    /// at `times` (mobile, trapped, outlet flux), chained at `times[0]`:
    /// the first segment takes resolved steps (dt_max = 0.1 s) so the
    /// Crank–Nicolson fine-mode ringing excited by the t = 0 kink (|g| → 1
    /// on stiff modes, decay time ≈ λ_max dt²/4 ≈ 2.6 s here) is dead well
    /// before the first comparison time, and the rest runs at `dt_b` from
    /// the smooth segment-end state (no kink, so nothing to ring). `dt_b`
    /// is per case: 1 s for the trap-free slab (CN phase error ≈ 2e-6
    /// there), 0.25 s for the trap-coupled slab (the mobile/trap splitting
    /// error of the fused Picard grows with dt, visible at 5e-4 on the
    /// trapped profile at dt = 1 s).
    #[allow(clippy::type_complexity)]
    fn g6e_production(
        params: &TransportParams,
        times: &[f64],
        dt_b: f64,
    ) -> (Vec<Vec<f64>>, Vec<Vec<Vec<f64>>>, Vec<f64>) {
        let left = Boundary::dirichlet(1.0).unwrap();
        let right = Boundary::recombination(1e-6).unwrap();
        let seg_a = solve(
            params,
            &left,
            &right,
            &TimeGrid::new(vec![times[0]]).unwrap(),
            &InitialState::zeros(params),
            &SolverOptions {
                dt_max: 0.1,
                rtol: 1e-9,
                atol: 1e-13,
                ..Default::default()
            },
        )
        .unwrap();
        let last = seg_a.times.len() - 1;
        let init = InitialState::new(
            params,
            seg_a.mobile[last].clone(),
            seg_a.trapped[last].clone(),
        )
        .unwrap();
        let rest: Vec<f64> = times[1..].iter().map(|&t| t - times[0]).collect();
        let seg_b = solve(
            params,
            &left,
            &right,
            &TimeGrid::new(rest).unwrap(),
            &init,
            &SolverOptions {
                dt_max: dt_b,
                rtol: 1e-9,
                atol: 1e-13,
                ..Default::default()
            },
        )
        .unwrap();
        let mut rows = vec![seg_a.mobile[last].clone()];
        rows.extend(seg_b.mobile.iter().cloned());
        let mut trapped = vec![seg_a.trapped[last].clone()];
        trapped.extend(seg_b.trapped.iter().cloned());
        let mut fluxes = vec![seg_a.flux_right[last]];
        fluxes.extend_from_slice(&seg_b.flux_right);
        (rows, trapped, fluxes)
    }

    #[test]
    fn g6e_independent_mol_crosscheck() {
        // G6e: G6a–G6d are self-referential (the same code on refined steps
        // and grids), so a consistent spatial-discretization or face-formula
        // bug could pass them all. This gate cross-checks the mid-transient
        // trajectory against an independent solver — method of lines with
        // second-order central differences on a node-centred grid and a
        // hand-rolled explicit RK4 — sharing no discretization, face
        // closure, or time integrator with the production kernel. Cases at
        // t = 0.5/1/2 t_lag on the G6a slab (D = 1e-9, L = 1e-3,
        // Dirichlet(1) upstream, recombination(1e-6) downstream, clean
        // start): trap-free, plus a trap-coupled variant (k, p, N) =
        // (0.05, 0.01, 2) at 500 K — the McNabb–Foster map is a pointwise
        // ODE, so the MoL solver carries it inside the same RK4 stages at
        // negligible extra stiffness (k·c, p ≪ D/h²), which is why the
        // trap-coupled case is in this gate rather than skipped. Production
        // uses
        // 512 cells; the MoL uses 256 intervals at dt = 0.45 h²/D (inside
        // the h²/2D parabolic stability limit). Both methods and both face
        // closures are O(h²), so the cross-method difference is expected at
        // the 1e-5 level relative on the profile max-norm; measured worst
        // values are ≈2.4e-5 (mobile profile), ≈1.2e-4 (trapped profile),
        // and ≈2e-4 relative on the outlet flux where it is significant,
        // pinned below at ~5x headroom. The early-time outlet fluxes are
        // tiny (≤ 1e-8 mol/m²/s), so the flux pin is relative with a floor
        // at 1e-3 of the G5a steady flux — a face-formula bug shows up at
        // O(1) relative error late in the transient, where the flux is
        // significant.
        let lag = time_lag(1e-3, 1e-9).unwrap();
        let times = vec![0.5 * lag, lag, 2.0 * lag];
        let mol_n = 256;
        let h_mol = 1e-3 / mol_n as f64;
        let mol_free = mol_trajectory(1e-3, mol_n, 1e-9, 1e-6, &[], &times, 0.45);
        let mol_trap = mol_trajectory(1e-3, mol_n, 1e-9, 1e-6, &[(0.05, 0.01, 2.0)], &times, 0.45);

        let p = no_traps(1e-3, 512, 1e-9);
        let pt = TransportParams::new(
            1e-3,
            512,
            1e-9,
            0.0,
            vec![crate::TrapSpec::new(0.05, 0.0, 0.01, 0.0, 2.0).unwrap()],
            vec![500.0],
            vec![],
        )
        .unwrap();
        let (rows_free, _, flux_free) = g6e_production(&p, &times, 1.0);
        let (rows_trap, trapped_trap, flux_trap) = g6e_production(&pt, &times, 0.25);

        let mut worst_profile = 0.0_f64;
        let mut worst_trapped = 0.0_f64;
        for k in 0..times.len() {
            let x = |i: usize| (i as f64 + 0.5) * 1e-3 / 512.0;
            // Trap-free mobile profile at the 512 cell centres.
            let sampled: Vec<f64> = (0..512)
                .map(|i| mol_sample(&mol_free[k], mol_n, 0, None, h_mol, x(i)))
                .collect();
            worst_profile = worst_profile.max(rel_max_diff(&rows_free[k], &sampled));
            // Trap-coupled mobile and trapped profiles.
            let sampled_t: Vec<f64> = (0..512)
                .map(|i| mol_sample(&mol_trap[k], mol_n, 1, None, h_mol, x(i)))
                .collect();
            worst_profile = worst_profile.max(rel_max_diff(&rows_trap[k], &sampled_t));
            let sampled_ct: Vec<f64> = (0..512)
                .map(|i| mol_sample(&mol_trap[k], mol_n, 1, Some(0), h_mol, x(i)))
                .collect();
            let prod_ct: Vec<f64> = trapped_trap[k].iter().map(|r| r[0]).collect();
            worst_trapped = worst_trapped.max(rel_max_diff(&prod_ct, &sampled_ct));
            // Outlet fluxes (K_r cf² on both sides; production reports the
            // converged face value, the MoL the nodal mobile value at the
            // boundary node x = L).
            for (prod, mol_row) in [(&flux_free, &mol_free[k]), (&flux_trap, &mol_trap[k])] {
                let mol_flux = 1e-6 * mol_row[mol_n] * mol_row[mol_n];
                assert!(
                    (prod[k] - mol_flux).abs() <= 5e-3 * prod[k].abs().max(1e-3 * G6_JSS),
                    "t = {}: flux {} vs {mol_flux}",
                    times[k],
                    prod[k]
                );
            }
        }
        assert!(
            worst_profile <= 1.5e-4,
            "cross-method mobile profile error {worst_profile}"
        );
        assert!(
            worst_trapped <= 6e-4,
            "cross-method trapped profile error {worst_trapped}"
        );
    }

    #[test]
    fn g6_positivity_with_face_clamp() {
        // Positivity under the face clamp (cf ≥ 0): one-end and two-end
        // cases stay nonneg on a coarse grid.
        let p = no_traps(1e-3, 8, 1e-9);
        let opts = SolverOptions {
            dt_max: 1.0,
            ..Default::default()
        };
        let one = solve(
            &p,
            &Boundary::dirichlet(1.0).unwrap(),
            &Boundary::recombination(1e-6).unwrap(),
            &TimeGrid::new(vec![1.0, 10.0, 100.0]).unwrap(),
            &InitialState::zeros(&p),
            &opts,
        )
        .unwrap();
        for row in &one.mobile {
            assert!(row.iter().all(|&c| c >= 0.0));
        }
        assert!(one.flux_right.iter().all(|&f| f >= 0.0));
        // Two recombination ends draining a uniform load (exercises the
        // pair Newton): symmetric drain, nonneg throughout.
        let two = solve(
            &p,
            &Boundary::recombination(1e-6).unwrap(),
            &Boundary::recombination(2e-6).unwrap(),
            &TimeGrid::new(vec![1.0, 10.0, 100.0]).unwrap(),
            &InitialState::uniform(&p, 1.0, vec![]).unwrap(),
            &opts,
        )
        .unwrap();
        for row in &two.mobile {
            assert!(row.iter().all(|&c| c >= 0.0));
        }
        assert!(two.flux_left.iter().all(|&f| f >= 0.0));
        assert!(two.flux_right.iter().all(|&f| f >= 0.0));
    }

    #[test]
    fn recombination_face_failure_is_not_converged() {
        // Failure surface: stiff traps (the MAX_PICARD 2-cycle) behind a
        // stiff recombination end with a large dt_max fails the fused
        // face/trap Picard loop as the existing NotConverged (no new
        // variant, no dt retry).
        let stiff = TransportParams::new(
            1e-3,
            8,
            1e-9,
            0.0,
            vec![crate::TrapSpec::new(0.5, 0.0, 1e-12, 0.0, 3.0).unwrap()],
            vec![500.0],
            vec![],
        )
        .unwrap();
        let init_s = InitialState::uniform(&stiff, 1.0, vec![0.0]).unwrap();
        assert_eq!(
            solve(
                &stiff,
                &Boundary::ZeroFlux,
                &Boundary::recombination(1e-3).unwrap(),
                &TimeGrid::new(vec![1.0]).unwrap(),
                &init_s,
                &SolverOptions {
                    dt_max: 1.0,
                    rtol: 1e-12,
                    atol: 1e-15,
                    ..Default::default()
                },
            ),
            Err(Error::NotConverged)
        );
    }

    #[test]
    fn rejects_bad_options_and_states() {
        let p = no_traps(1e-3, 8, 1e-9);
        assert!(TimeGrid::new(vec![]).is_err());
        assert!(TimeGrid::new(vec![1.0, 1.0]).is_err());
        assert!(TimeGrid::new(vec![-1.0]).is_err());
        assert!(InitialState::new(&p, vec![0.0; 7], vec![vec![]; 8]).is_err());
        assert!(InitialState::new(&p, vec![-1.0; 8], vec![vec![]; 8]).is_err());
        assert!(InitialState::new(&p, vec![f64::NAN; 8], vec![vec![]; 8]).is_err());
        // Trapped-shape branches need a trap-bearing geometry.
        let tp = TransportParams::new(
            1e-3,
            8,
            1e-9,
            0.0,
            vec![crate::TrapSpec::new(0.05, 0.0, 0.01, 0.0, 2.0).unwrap()],
            vec![500.0],
            vec![],
        )
        .unwrap();
        assert!(InitialState::new(&tp, vec![1.0; 8], vec![vec![0.0]; 7]).is_err());
        assert!(InitialState::new(&tp, vec![1.0; 8], vec![vec![]; 8]).is_err());
        assert!(InitialState::new(&tp, vec![1.0; 8], vec![vec![3.0]; 8]).is_err());
        assert!(InitialState::new(&tp, vec![1.0; 8], vec![vec![-0.5]; 8]).is_err());
        assert!(InitialState::uniform(&tp, f64::NAN, vec![0.0]).is_err());
        assert!(SolverOptions {
            rtol: -1.0,
            ..Default::default()
        }
        .validate()
        .is_err());
        assert!(SolverOptions {
            max_steps: 0,
            ..Default::default()
        }
        .validate()
        .is_err());
        assert!(SolverOptions {
            dt_min: 1.0,
            dt_max: 0.5,
            ..Default::default()
        }
        .validate()
        .is_err());
        // A dt_min above the output spacing is rejected, not stalled on.
        let init = InitialState::zeros(&p);
        let grid = TimeGrid::new(vec![0.5]).unwrap();
        let opts = SolverOptions {
            dt_min: 1.0,
            ..Default::default()
        };
        assert!(solve(
            &p,
            &Boundary::dirichlet(1.0).unwrap(),
            &Boundary::dirichlet(0.0).unwrap(),
            &grid,
            &init,
            &opts,
        )
        .is_err());
        // A hand-built state whose trapped rows miss the trap count is
        // rejected at solve time even though it type-checks.
        let bad_trapped = InitialState {
            mobile: vec![0.0; 8],
            trapped: vec![vec![]; 8],
        };
        let one_trap = TransportParams::new(
            1e-3,
            8,
            1e-9,
            0.0,
            vec![crate::TrapSpec::new(0.05, 0.0, 0.01, 0.0, 2.0).unwrap()],
            vec![500.0],
            vec![],
        )
        .unwrap();
        assert!(InitialState::uniform(&one_trap, 1.0, vec![]).is_err());
        assert!(solve(
            &one_trap,
            &Boundary::ZeroFlux,
            &Boundary::ZeroFlux,
            &TimeGrid::new(vec![1.0]).unwrap(),
            &bad_trapped,
            &SolverOptions::default(),
        )
        .is_err());
    }

    #[test]
    fn step_budget_and_picard_failure_surface() {
        let p = no_traps(1e-3, 8, 1e-9);
        let init = InitialState::zeros(&p);
        // Four internal steps needed, one allowed.
        let budgeted = SolverOptions {
            dt_max: 0.5,
            max_steps: 1,
            ..Default::default()
        };
        assert_eq!(
            solve(
                &p,
                &Boundary::dirichlet(1.0).unwrap(),
                &Boundary::dirichlet(0.0).unwrap(),
                &TimeGrid::new(vec![1.0, 2.0]).unwrap(),
                &init,
                &budgeted,
            ),
            Err(Error::StepBudget(1))
        );
        // A bounded Picard 2-cycle defeats the iteration cap: with
        // dt·k·N = 1.5 and dt·k = 0.5 the trap map is f(c) = 1.5c/(1+0.5c)
        // and the sealed uniform slab iterates h(c) = 1 − f(c), whose exact
        // {0, 1} 2-cycle never meets rtol/atol.
        let stiff = TransportParams::new(
            1e-3,
            8,
            1e-9,
            0.0,
            vec![crate::TrapSpec::new(0.5, 0.0, 1e-12, 0.0, 3.0).unwrap()],
            vec![500.0],
            vec![],
        )
        .unwrap();
        let init_s = InitialState::uniform(&stiff, 1.0, vec![0.0]).unwrap();
        assert_eq!(
            solve(
                &stiff,
                &Boundary::ZeroFlux,
                &Boundary::ZeroFlux,
                &TimeGrid::new(vec![1.0]).unwrap(),
                &init_s,
                &SolverOptions {
                    dt_max: 1.0,
                    rtol: 1e-12,
                    atol: 1e-15,
                    ..Default::default()
                },
            ),
            Err(Error::NotConverged)
        );
    }
}
