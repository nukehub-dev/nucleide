//! Multi-layer series stacks with linear interface conditions (G7–G9).
//!
//! A [`LayerStack`] chains `N ≥ 1` slabs — e.g. a W/Cu/CuCrZr first-wall
//! stack — each with its own thickness, diffusivity, solubility, trap
//! species, temperature, and source. Inside every layer the mobile balance
//! (T1) and the McNabb–Foster kinetics (T2) hold with the layer's own
//! coefficients; temperature stays caller-supplied per layer (no heat solve).
//!
//! Internal interfaces carry a linear local-equilibrium law: the Sieverts
//! potential `u = c_m / K_S` \[Pa¹ᐟ²\] or the Henry potential `u = c_m / K_H`
//! \[Pa\] (the layer's [`LayerSpec::solubility`] plays whichever constant the
//! adjacent interface's [`Interface`] law declares) and the flux
//! `J = −D ∂c_m/∂x` are continuous. Recombination internal interface laws
//! remain loud [`Error::UnsupportedInterface`] errors (recorded limitation).
//! Because both supported conditions are linear in the mobile concentration,
//! every interior face flux
//!
//! ```text
//! J_f = (c_i / K_i − c_{i+1} / K_{i+1}) / R_f,
//! R_f = dx_i / (2 Φ_i) + dx_{i+1} / (2 Φ_{i+1}),   Φ = D · K,
//! ```
//!
//! folds directly into the tridiagonal θ-step matrix — the interface needs
//! no Newton machinery of its own (strictly simpler than the recombination
//! faces of G5/G6, which are nonlinear in the face value). The outer ends
//! reuse the landed [`Boundary`] taxonomy unchanged, including recombination
//! ends through the same affine face-response construction (G5/G6), here
//! built on the layered step matrix. A one-layer stack dispatches to the
//! landed single-slab kernel, so `N = 1` reproduces it exactly (G7c).

use crate::bc::Boundary;
use crate::error::Error;
use crate::params::{arrhenius, TransportParams, TrapSpec};
use crate::solve::{
    self, apply_faces, close_faces, face_closed_single, face_concentration, face_from_adjacent,
    face_newton_pair, linear_steady, outward_flux, trap_update, DiffusionSystem, InitialState,
    SolverOptions, TimeGrid, G5_ATOL, G5_RTOL, MAX_PICARD,
};

/// Internal interface condition between adjacent layers.
///
/// The taxonomy mirrors [`crate::bc::Boundary`]. The linear laws
/// [`Sieverts`] and [`Henry`] are supported: the local-equilibrium condition
/// `u = c_m / K` continuous (Sieverts potential with `K = K_S`, Henry
/// potential with `K = K_H`; the layer [`LayerSpec::solubility`] carries the
/// matching constant) with continuous flux. [`Recombination`] internal
/// interfaces are rejected loudly by [`LayerStack::new`] with
/// [`Error::UnsupportedInterface`] (recorded limitation, not a silent pass).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Interface {
    /// Sieverts local equilibrium: `c_m / K_S` continuous, flux continuous.
    Sieverts,
    /// Linear (Henry) local equilibrium: `c_m / K_H` continuous, flux
    /// continuous. Shares the Sieverts resistance form (the layer
    /// [`LayerSpec::solubility`] plays `K_H`), so it folds into the
    /// θ-step matrix exactly like [`Interface::Sieverts`].
    Henry,
    /// Recombination interface law `J = K_r c²` — not supported (loud
    /// error): the nonlinear face would need an unproven per-interface
    /// Newton construction (recorded for a later cycle).
    Recombination,
}

impl Interface {
    /// Support check: the linear laws pass; recombination fails loudly.
    fn validate(&self) -> Result<(), Error> {
        match self {
            Interface::Sieverts | Interface::Henry => Ok(()),
            Interface::Recombination => Err(Error::UnsupportedInterface(
                "internal recombination interfaces are not supported (Sieverts and Henry only)",
            )),
        }
    }
}

/// One layer of a [`LayerStack`]: thickness, cell count, transport data.
///
/// All coefficients are caller data with the landed Arrhenius stance
/// (T-Arr): `D(T) = d0·exp(−e_d/RT)`; the solubility `K_S`
/// \[mol/m³/Pa¹ᐟ²\] is a constant defining the layer's Sieverts potential
/// `u = c_m / K_S` at internal interfaces (outer ends carry their own
/// surface laws). Trap species, temperature, and source follow
/// [`TransportParams`]: one value (uniform) or one per cell of this layer.
#[derive(Debug, Clone, PartialEq)]
pub struct LayerSpec {
    /// Layer thickness \[m\] (`> 0`).
    pub thickness: f64,
    /// Finite-volume cells in this layer (`>= 1`).
    pub cells: usize,
    /// Diffusivity pre-factor `D_0` \[m²/s\] (`> 0`).
    pub d0: f64,
    /// Diffusivity activation energy `E_D` \[J/mol\] (`>= 0`).
    pub e_d: f64,
    /// Sieverts solubility `K_S` \[mol/m³/Pa¹ᐟ²\] or Henry constant `K_H`
    /// \[mol/m³/Pa\] (`> 0`), per the law declared at the adjacent internal
    /// interfaces: the layer's potential `u = c_m / K` (Sieverts potential
    /// under [`Interface::Sieverts`], Henry potential under
    /// [`Interface::Henry`]) at internal interfaces.
    pub solubility: f64,
    /// Trap species of this layer (possibly empty: pure Fickian diffusion).
    pub traps: Vec<TrapSpec>,
    /// Temperature \[K\]: one value (uniform) or one per cell of this layer.
    pub temperature: Vec<f64>,
    /// Volumetric source `S` \[mol/m³/s\]: empty (zero), one value, or one
    /// per cell of this layer.
    pub source: Vec<f64>,
}

impl LayerSpec {
    /// Validate one layer: positive finite thickness, at least one cell,
    /// finite positive diffusivity/solubility data, finite non-negative
    /// activation energy, and uniform-or-per-cell temperature/source
    /// profiles matching this layer's cell count.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        thickness: f64,
        cells: usize,
        d0: f64,
        e_d: f64,
        solubility: f64,
        traps: Vec<TrapSpec>,
        temperature: Vec<f64>,
        source: Vec<f64>,
    ) -> Result<Self, Error> {
        if !thickness.is_finite() || thickness <= 0.0 {
            return Err(Error::BadData("layer thickness must be finite and > 0"));
        }
        if cells < 1 {
            return Err(Error::BadGrid("need at least 1 cell per layer"));
        }
        if !d0.is_finite() || d0 <= 0.0 {
            return Err(Error::BadData("layer D0 must be finite and > 0"));
        }
        if !e_d.is_finite() || e_d < 0.0 {
            return Err(Error::BadData("layer E_D must be finite and >= 0"));
        }
        if !solubility.is_finite() || solubility <= 0.0 {
            return Err(Error::BadData("layer solubility must be finite and > 0"));
        }
        if temperature.len() != 1 && temperature.len() != cells {
            return Err(Error::BadData(
                "layer temperature must be uniform (1 value) or per-cell",
            ));
        }
        if temperature.iter().any(|t| !t.is_finite() || *t <= 0.0) {
            return Err(Error::BadData("layer temperatures must be finite and > 0"));
        }
        if !source.is_empty() && source.len() != 1 && source.len() != cells {
            return Err(Error::BadData(
                "layer source must be empty (zero), uniform (1 value), or per-cell",
            ));
        }
        if source.iter().any(|s| !s.is_finite()) {
            return Err(Error::BadData("layer source values must be finite"));
        }
        Ok(Self {
            thickness,
            cells,
            d0,
            e_d,
            solubility,
            traps,
            temperature,
            source,
        })
    }
}

/// A series stack of [`LayerSpec`] layers with one [`Interface`] per gap.
///
/// Cell indices run across the whole stack (layer 0 first); per-layer
/// coefficients expand to per-cell arrays on demand. Construction validates
/// every layer and rejects recombination interfaces loudly (see
/// [`Interface`]).
#[derive(Debug, Clone, PartialEq)]
pub struct LayerStack {
    /// Layers from the left (`x = 0`) face to the right (`x = L`) face.
    pub layers: Vec<LayerSpec>,
    /// Interface condition between `layers[k]` and `layers[k + 1]`
    /// (length `layers.len() − 1`; Sieverts or Henry — see [`Interface`]).
    pub interfaces: Vec<Interface>,
    /// Cached first-cell index of each layer (length `layers.len()`).
    /// Derived from `layers` at construction so per-cell queries stay
    /// O(log layers) instead of re-summing offsets on every call.
    starts: Vec<usize>,
}

impl LayerStack {
    /// Validate a stack: at least one layer, one interface per gap, every
    /// interface supported, and a total cell count that fits in `usize`.
    pub fn new(layers: Vec<LayerSpec>, interfaces: Vec<Interface>) -> Result<Self, Error> {
        if layers.is_empty() {
            return Err(Error::BadGrid("a layer stack needs at least one layer"));
        }
        if interfaces.len() + 1 != layers.len() {
            return Err(Error::BadGrid(
                "need exactly one interface between adjacent layers",
            ));
        }
        for interface in &interfaces {
            interface.validate()?;
        }
        // Validate the total cell count fits in usize, caching each
        // layer's first-cell offset as we go (see `starts`).
        let mut starts = Vec::with_capacity(layers.len());
        layers.iter().try_fold(0_usize, |total, layer| {
            starts.push(total);
            total
                .checked_add(layer.cells)
                .ok_or(Error::BadGrid("total cell count overflows usize"))
        })?;
        Ok(Self {
            layers,
            interfaces,
            starts,
        })
    }

    /// Total cell count across all layers.
    pub fn total_cells(&self) -> usize {
        self.layers.iter().map(|l| l.cells).sum()
    }

    /// Total stack length `Σ thickness` \[m\].
    pub fn total_length(&self) -> f64 {
        self.layers.iter().map(|l| l.thickness).sum()
    }

    /// Index of the layer owning cell `i`, or `None` when `i` is outside the
    /// stack's total cell count. Binary search over the cached `starts`
    /// offsets: O(log layers), safe to call per cell inside stepper loops.
    fn layer_of(&self, i: usize) -> Option<usize> {
        let k = self.starts.partition_point(|&s| s <= i).checked_sub(1)?;
        (i < self.starts[k] + self.layers[k].cells).then_some(k)
    }

    /// Loud [`Error::BadCellIndex`] for a public per-cell query.
    fn bad_index(&self, i: usize) -> Error {
        Error::BadCellIndex {
            index: i,
            total: self.total_cells(),
        }
    }

    /// Per-cell widths `dx` \[m\].
    pub fn cell_widths(&self) -> Vec<f64> {
        let mut dx = Vec::with_capacity(self.total_cells());
        for layer in &self.layers {
            dx.extend(std::iter::repeat_n(
                layer.thickness / layer.cells as f64,
                layer.cells,
            ));
        }
        dx
    }

    /// Cell-centre positions across the stack \[m\].
    pub fn cell_centres(&self) -> Vec<f64> {
        let mut centres = Vec::with_capacity(self.total_cells());
        let mut x0 = 0.0_f64;
        for layer in &self.layers {
            let dx = layer.thickness / layer.cells as f64;
            centres.extend((0..layer.cells).map(|i| x0 + (i as f64 + 0.5) * dx));
            x0 += layer.thickness;
        }
        centres
    }

    /// Temperature at cell `i` \[K\] (uniform broadcasts). A cell index
    /// outside the stack is a loud [`Error::BadCellIndex`].
    pub fn temp_at(&self, i: usize) -> Result<f64, Error> {
        let k = self.layer_of(i).ok_or_else(|| self.bad_index(i))?;
        let layer = &self.layers[k];
        if layer.temperature.len() == 1 {
            Ok(layer.temperature[0])
        } else {
            Ok(layer.temperature[i - self.layer_start(k)])
        }
    }

    /// First cell index of layer `k` (cached at construction).
    fn layer_start(&self, k: usize) -> usize {
        self.starts[k]
    }

    /// Source at cell `i` \[mol/m³/s\] (empty means zero, uniform broadcasts).
    /// A cell index outside the stack is a loud [`Error::BadCellIndex`].
    pub fn source_at(&self, i: usize) -> Result<f64, Error> {
        let k = self.layer_of(i).ok_or_else(|| self.bad_index(i))?;
        let layer = &self.layers[k];
        if layer.source.is_empty() {
            Ok(0.0)
        } else if layer.source.len() == 1 {
            Ok(layer.source[0])
        } else {
            Ok(layer.source[i - self.layer_start(k)])
        }
    }

    /// Diffusivity at cell `i` \[m²/s\] (Arrhenius in the cell temperature).
    pub fn diffusivity_at(&self, i: usize) -> Result<f64, Error> {
        let layer = &self.layers[self.layer_of(i).ok_or_else(|| self.bad_index(i))?];
        arrhenius(layer.d0, layer.e_d, self.temp_at(i)?)
    }

    /// Diffusivity at every cell \[m²/s\].
    pub fn diffusivities(&self) -> Result<Vec<f64>, Error> {
        (0..self.total_cells())
            .map(|i| self.diffusivity_at(i))
            .collect()
    }

    /// The layer's interface constant `K` at cell `i` — `K_S`
    /// \[mol/m³/Pa¹ᐟ²\] under a Sieverts law, `K_H` \[mol/m³/Pa\] under a
    /// Henry law (see [`Interface`]). A cell index outside the stack is a
    /// loud [`Error::BadCellIndex`].
    pub fn solubility_at(&self, i: usize) -> Result<f64, Error> {
        Ok(self.layers[self.layer_of(i).ok_or_else(|| self.bad_index(i))?].solubility)
    }

    /// Trap `(k, p)` rates at cell `i` for every species of the owning layer.
    pub fn trap_rates_at(&self, i: usize) -> Result<Vec<(f64, f64)>, Error> {
        let t = self.temp_at(i)?;
        self.layers[self.layer_of(i).ok_or_else(|| self.bad_index(i))?]
            .traps
            .iter()
            .map(|trap| trap.rates(t))
            .collect()
    }

    /// Zero initial state for the stack geometry (mobile and traps zero).
    pub fn zero_state(&self) -> InitialState {
        let mut trapped = Vec::with_capacity(self.total_cells());
        for layer in &self.layers {
            trapped.extend(std::iter::repeat_n(
                vec![0.0_f64; layer.traps.len()],
                layer.cells,
            ));
        }
        InitialState {
            mobile: vec![0.0; self.total_cells()],
            trapped,
        }
    }

    /// Validate an initial state against the stack geometry: mobile length,
    /// per-cell trapped rows matching the owning layer's trap count, and the
    /// landed sign/finiteness/site bounds.
    fn check_state(&self, initial: &InitialState) -> Result<(), Error> {
        let n = self.total_cells();
        if initial.mobile.len() != n {
            return Err(Error::BadState(
                "initial mobile length must match total cell count",
            ));
        }
        if initial.mobile.iter().any(|c| !c.is_finite() || *c < 0.0) {
            return Err(Error::BadState("initial mobile must be finite and >= 0"));
        }
        if initial.trapped.len() != n {
            return Err(Error::BadState(
                "initial trapped outer length must match total cell count",
            ));
        }
        for i in 0..n {
            let traps =
                &self.layers[self.layer_of(i).expect("i < total_cells by construction")].traps;
            let row = &initial.trapped[i];
            if row.len() != traps.len() {
                return Err(Error::BadState(
                    "initial trapped inner length must match the owning layer's trap count",
                ));
            }
            for (j, ct) in row.iter().enumerate() {
                let max = traps[j].site_density;
                if !ct.is_finite() || *ct < 0.0 || *ct > max {
                    return Err(Error::BadState(
                        "initial trapped loads must satisfy 0 <= ct <= N_j",
                    ));
                }
            }
        }
        Ok(())
    }
}

/// Steady state of a [`LayerStack`]: per-cell mobile profile, Langmuir
/// trapped loads on the owning layer's isotherm, boundary fluxes, and
/// per-cell geometry for inventory integrals. Mirrors
/// [`crate::solve::SteadyState`] with per-cell widths.
#[derive(Debug, Clone, PartialEq)]
pub struct LayeredSteadyState {
    /// Cell-centre positions across the stack \[m\].
    pub centres: Vec<f64>,
    /// Per-cell widths \[m\].
    pub dx: Vec<f64>,
    /// Mobile concentration per cell \[mol/m³\].
    pub mobile: Vec<f64>,
    /// Trapped concentrations `[cell][trap]` \[mol/m³\].
    pub trapped: Vec<Vec<f64>>,
    /// Outward flux at `x = 0` \[mol/m²/s\] (positive leaves the stack).
    pub flux_left: f64,
    /// Outward flux at `x = L` \[mol/m²/s\] (positive leaves the stack).
    pub flux_right: f64,
    /// Mobile inventory per unit area \[mol/m²\].
    pub inventory_mobile: f64,
    /// Trapped inventory per unit area \[mol/m²\].
    pub inventory_trapped: f64,
}

/// Full transient of a [`LayerStack`]: one mobile/trapped row plus outward
/// surface fluxes per output time. Mirrors [`crate::solve::Solution`] with
/// per-cell widths.
#[derive(Debug, Clone, PartialEq)]
pub struct LayeredSolution {
    /// Output times \[s\] (echo of the grid).
    pub times: Vec<f64>,
    /// Mobile concentration at each output time (`[time][cell]`).
    pub mobile: Vec<Vec<f64>>,
    /// Trapped concentrations (`[time][cell][trap]`).
    pub trapped: Vec<Vec<Vec<f64>>>,
    /// Outward flux at `x = 0` \[mol/m²/s\] (positive leaves the stack).
    pub flux_left: Vec<f64>,
    /// Outward flux at `x = L` \[mol/m²/s\] (positive leaves the stack).
    pub flux_right: Vec<f64>,
    /// Per-cell widths \[m\] (for inventory integrals).
    pub dx: Vec<f64>,
    /// The `t = 0` state.
    pub initial: InitialState,
}

impl LayeredSolution {
    /// Mobile inventory per unit area at each output time \[mol/m²\].
    pub fn inventory_mobile(&self) -> Vec<f64> {
        self.mobile
            .iter()
            .map(|row| row.iter().zip(&self.dx).map(|(c, w)| c * w).sum())
            .collect()
    }

    /// Trapped inventory per unit area at each output time \[mol/m²\].
    pub fn inventory_trapped(&self) -> Vec<f64> {
        self.trapped
            .iter()
            .map(|rows| {
                rows.iter()
                    .zip(&self.dx)
                    .map(|(row, w)| row.iter().sum::<f64>() * w)
                    .sum()
            })
            .collect()
    }

    /// Total (mobile + trapped) inventory per unit area \[mol/m²\].
    pub fn inventory_total(&self) -> Vec<f64> {
        let m = self.inventory_mobile();
        let t = self.inventory_trapped();
        m.iter().zip(&t).map(|(a, b)| a + b).collect()
    }
}

// ---------------------------------------------------------------------------
// Layered discrete operator
// ---------------------------------------------------------------------------

/// Layered rate system `dc/dt = A c + rhs` plus the per-cell geometry and
/// face resistances needed for flux evaluation. `A` is stored in
/// sub/diag/sup form exactly like the landed [`DiffusionSystem`]; the extra
/// fields carry the interface data.
struct LayeredSystem {
    sub: Vec<f64>,
    diag: Vec<f64>,
    sup: Vec<f64>,
    rhs: Vec<f64>,
    /// Per-cell diffusivities (boundary flux evaluation).
    d_cell: Vec<f64>,
    /// Per-cell widths.
    dx: Vec<f64>,
}

impl LayeredSystem {
    /// Borrow the tridiagonal rate system alone (for [`linear_steady`]).
    fn rate_system(&self) -> DiffusionSystem {
        DiffusionSystem {
            sub: self.sub.clone(),
            diag: self.diag.clone(),
            sup: self.sup.clone(),
            rhs: self.rhs.clone(),
            face_d: self.d_cell.clone(),
        }
    }
}

/// Assemble the layered cell-centred finite-volume diffusion operator.
///
/// Interior faces carry the linear interface flux — Sieverts or Henry, the
/// same resistance form with the layer [`LayerSpec::solubility`] playing
/// `K_S` or `K_H` per the adjacent [`Interface`] (also *within* one layer,
/// where it reduces to the landed arithmetic-mean face for uniform `D`):
/// `J_f = (c_i/K_i − c_{i+1}/K_{i+1})/R_f`. Boundary faces reuse the landed
/// half-cell conductance `2 D_cell/dx` against the equilibrium face
/// concentration. The mobile concentration stays the unknown, so the face
/// fluxes are linear and the matrix tridiagonal.
fn assemble_layers(
    stack: &LayerStack,
    left: &Boundary,
    right: &Boundary,
) -> Result<LayeredSystem, Error> {
    let n = stack.total_cells();
    let dx = stack.cell_widths();
    let d_cell = stack.diffusivities()?;
    let sol: Vec<f64> = (0..n)
        .map(|i| stack.solubility_at(i))
        .collect::<std::result::Result<Vec<_>, Error>>()?;
    let phi: Vec<f64> = d_cell.iter().zip(&sol).map(|(d, s)| d * s).collect();
    let mut face_r = vec![0.0; n + 1];
    face_r[0] = dx[0] / (2.0 * phi[0]);
    face_r[n] = dx[n - 1] / (2.0 * phi[n - 1]);
    for i in 1..n {
        face_r[i] = dx[i - 1] / (2.0 * phi[i - 1]) + dx[i] / (2.0 * phi[i]);
    }
    let mut sub = vec![0.0; n.saturating_sub(1)];
    let mut diag = vec![0.0; n];
    let mut sup = vec![0.0; n.saturating_sub(1)];
    let mut rhs = vec![0.0; n];
    for i in 0..n {
        // West coupling: interior face (Sieverts form) or landed boundary.
        if i > 0 {
            let rw = face_r[i];
            diag[i] -= 1.0 / (rw * dx[i] * sol[i]);
            sub[i - 1] = 1.0 / (rw * dx[i] * sol[i - 1]);
        } else if let Some(cf) = face_concentration(left) {
            diag[i] -= 2.0 * d_cell[i] / (dx[i] * dx[i]);
            rhs[i] += 2.0 * d_cell[i] * cf / (dx[i] * dx[i]);
        }
        // East coupling.
        if i + 1 < n {
            let re = face_r[i + 1];
            diag[i] -= 1.0 / (re * dx[i] * sol[i]);
            sup[i] = 1.0 / (re * dx[i] * sol[i + 1]);
        } else if let Some(cf) = face_concentration(right) {
            diag[i] -= 2.0 * d_cell[i] / (dx[i] * dx[i]);
            rhs[i] += 2.0 * d_cell[i] * cf / (dx[i] * dx[i]);
        }
        rhs[i] += stack.source_at(i)?;
    }
    Ok(LayeredSystem {
        sub,
        diag,
        sup,
        rhs,
        d_cell,
        dx,
    })
}

/// Layered steady bundle: Langmuir trapped loads on each layer's isotherm
/// plus inventories for a converged mobile profile.
fn finish_layered_steady(
    stack: &LayerStack,
    sys: &LayeredSystem,
    mobile: Vec<f64>,
    flux_left: f64,
    flux_right: f64,
) -> Result<LayeredSteadyState, Error> {
    let n = stack.total_cells();
    let mut trapped = Vec::with_capacity(n);
    for (i, &c) in mobile.iter().enumerate() {
        let rates = stack.trap_rates_at(i)?;
        let mut row = Vec::with_capacity(rates.len());
        for (j, (k, p)) in rates.iter().enumerate() {
            let traps =
                &stack.layers[stack.layer_of(i).expect("i < total_cells by construction")].traps;
            row.push(solve_langmuir(traps[j].site_density, *k, *p, c)?);
        }
        trapped.push(row);
    }
    let inventory_mobile: f64 = mobile.iter().zip(&sys.dx).map(|(c, w)| c * w).sum();
    let inventory_trapped: f64 = trapped
        .iter()
        .zip(&sys.dx)
        .map(|(row, w)| row.iter().sum::<f64>() * w)
        .sum();
    Ok(LayeredSteadyState {
        centres: stack.cell_centres(),
        dx: sys.dx.clone(),
        mobile,
        trapped,
        flux_left,
        flux_right,
        inventory_mobile,
        inventory_trapped,
    })
}

/// Langmuir load `N K c/(1 + K c)` with `K = k/p` (T2-eq).
fn solve_langmuir(site_density: f64, k: f64, p: f64, c: f64) -> Result<f64, Error> {
    let keq = k / p;
    solve::equilibrium_trapped(site_density, keq, c)
}

/// Map a single-layer stack onto the landed [`TransportParams`] (the G7c
/// dispatch keeps `N = 1` on the landed kernel exactly).
fn landed_params(layer: &LayerSpec) -> Result<TransportParams, Error> {
    TransportParams::new(
        layer.thickness,
        layer.cells,
        layer.d0,
        layer.e_d,
        layer.traps.clone(),
        layer.temperature.clone(),
        layer.source.clone(),
    )
}

// ---------------------------------------------------------------------------
// Steady state (G7)
// ---------------------------------------------------------------------------

/// Steady state of a [`LayerStack`] with at least one recombination end.
///
/// Mirrors [`solve::steady_state`]'s G5 face-response construction on the
/// layered operator: for frozen faces the steady system is linear, so one
/// basis profile plus one unit-perturbation column per recombination end
/// pins the face-adjacent response exactly; one end closes in closed form,
/// two by the analytic-Jacobian Newton (shared helpers, G5 tolerances).
fn steady_layers_recombination(
    stack: &LayerStack,
    left: &Boundary,
    right: &Boundary,
    kr_left: Option<f64>,
    kr_right: Option<f64>,
) -> Result<LayeredSteadyState, Error> {
    let n = stack.total_cells();
    let dx = stack.cell_widths();
    let d_cell = stack.diffusivities()?;
    let ends: Vec<(usize, f64, f64)> = [
        kr_left.map(|kr| (0, kr, 2.0 * d_cell[0] / dx[0])),
        kr_right.map(|kr| (1, kr, 2.0 * d_cell[n - 1] / dx[n - 1])),
    ]
    .into_iter()
    .flatten()
    .collect();
    let solve_faces = |face: [f64; 2]| -> Result<Vec<f64>, Error> {
        let bl = kr_left.map_or_else(|| left.clone(), |_| Boundary::Dirichlet(face[0]));
        let br = kr_right.map_or_else(|| right.clone(), |_| Boundary::Dirichlet(face[1]));
        linear_steady(&assemble_layers(stack, &bl, &br)?.rate_system())
    };
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
    let mut face = [0.0_f64; 2];
    match ends.as_slice() {
        [(e, kr, g)] => {
            face[*e] = face_closed_single(*kr, *g, resp[*e][0], resp[*e][1]);
        }
        [(e0, kr0, g0), (e1, kr1, g1)] => {
            let (a0, b00, b01) = (resp[*e0][0], resp[*e0][1], resp[*e0][2]);
            let (a1, b10, b11) = (resp[*e1][0], resp[*e1][1], resp[*e1][2]);
            let opts = SolverOptions {
                rtol: G5_RTOL,
                atol: G5_ATOL,
                ..Default::default()
            };
            let pair = face_newton_pair(
                *kr0,
                *g0,
                a0,
                b00,
                b01,
                *kr1,
                *g1,
                a1,
                b10,
                b11,
                [0.0, 0.0],
                &opts,
            )?;
            face[*e0] = pair[0];
            face[*e1] = pair[1];
        }
        _ => {
            // Unreachable-in-practice: callers dispatch here only with at
            // least one recombination end. Loud error, never a panic.
            return Err(Error::BadBoundary(
                "steady_layers_recombination needs a recombination end",
            ));
        }
    }
    let mobile = solve_faces(face)?;
    let bl = kr_left.map_or_else(|| left.clone(), |_| Boundary::Dirichlet(face[0]));
    let br = kr_right.map_or_else(|| right.clone(), |_| Boundary::Dirichlet(face[1]));
    let sys = assemble_layers(stack, &bl, &br)?;
    let flux_left = kr_left.map_or_else(
        || outward_flux(left, mobile[0], sys.d_cell[0], sys.dx[0]),
        |kr| kr * face[0] * face[0],
    );
    let flux_right = kr_right.map_or_else(
        || outward_flux(right, mobile[n - 1], sys.d_cell[n - 1], sys.dx[n - 1]),
        |kr| kr * face[1] * face[1],
    );
    finish_layered_steady(stack, &sys, mobile, flux_left, flux_right)
}

/// Trap-free-style steady state of a [`LayerStack`] (G7).
///
/// Solves the layered linear system through [`nucleide_linalg::tridiag`] and
/// evaluates each layer's Langmuir isotherm pointwise. A one-layer stack
/// dispatches to [`solve::steady_state`] and reproduces it exactly (G7c);
/// recombination outer ends close through the layered face-response
/// construction ([`steady_layers_recombination`], the landed G5 machinery on
/// the extended interface system).
pub fn steady_layers(
    stack: &LayerStack,
    left: &Boundary,
    right: &Boundary,
) -> Result<LayeredSteadyState, Error> {
    if stack.layers.len() == 1 {
        let params = landed_params(&stack.layers[0])?;
        let s = solve::steady_state(&params, left, right)?;
        let dx = vec![params.dx(); params.cells];
        return Ok(LayeredSteadyState {
            centres: s.centres,
            dx,
            mobile: s.mobile,
            trapped: s.trapped,
            flux_left: s.flux_left,
            flux_right: s.flux_right,
            inventory_mobile: s.inventory_mobile,
            inventory_trapped: s.inventory_trapped,
        });
    }
    let kr_left = left.recombination_rate();
    let kr_right = right.recombination_rate();
    if kr_left.is_some() || kr_right.is_some() {
        return steady_layers_recombination(stack, left, right, kr_left, kr_right);
    }
    let sys = assemble_layers(stack, left, right)?;
    let mobile = linear_steady(&sys.rate_system())?;
    let flux_left = outward_flux(left, mobile[0], sys.d_cell[0], sys.dx[0]);
    let n = stack.total_cells();
    let flux_right = outward_flux(right, mobile[n - 1], sys.d_cell[n - 1], sys.dx[n - 1]);
    finish_layered_steady(stack, &sys, mobile, flux_left, flux_right)
}

// ---------------------------------------------------------------------------
// Transient (G8)
// ---------------------------------------------------------------------------

/// Solve the layered (T1–T2) transient over `grid` from `initial`.
///
/// The θ-stepper mirrors [`solve::solve`]: diffusion implicit through the
/// shared `linalg::tridiag` Thomas solve, traps by the exact per-cell
/// backward-Euler map with Picard coupling to `rtol`/`atol`. The linear
/// interface fluxes (Sieverts/Henry) sit inside the step matrix and no
/// interface iteration exists; recombination outer ends (G6 machinery) close
/// per step through the affine face-response construction on the layered
/// matrix, fused into the same Picard loop. A one-layer stack dispatches to
/// [`solve::solve`] exactly (G7c).
pub fn solve_layers(
    stack: &LayerStack,
    left: &Boundary,
    right: &Boundary,
    grid: &TimeGrid,
    initial: &InitialState,
    opts: &SolverOptions,
) -> Result<LayeredSolution, Error> {
    opts.validate()?;
    stack.check_state(initial)?;
    if stack.layers.len() == 1 {
        let params = landed_params(&stack.layers[0])?;
        let sol = solve::solve(&params, left, right, grid, initial, opts)?;
        let n = stack.total_cells();
        return Ok(LayeredSolution {
            times: sol.times,
            mobile: sol.mobile,
            trapped: sol.trapped,
            flux_left: sol.flux_left,
            flux_right: sol.flux_right,
            dx: vec![sol.dx; n],
            initial: sol.initial,
        });
    }
    let kr_left = left.recombination_rate();
    let kr_right = right.recombination_rate();
    if kr_left.is_some() || kr_right.is_some() {
        return solve_layers_recombination(
            stack, left, right, grid, initial, opts, kr_left, kr_right,
        );
    }
    let sys = assemble_layers(stack, left, right)?;
    let n = stack.total_cells();
    let ntraps_per_cell: Vec<usize> = (0..n)
        .map(|i| {
            stack.layers[stack.layer_of(i).expect("i < total_cells by construction")]
                .traps
                .len()
        })
        .collect();
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
    let rates: Vec<Vec<(f64, f64)>> = (0..n)
        .map(|i| stack.trap_rates_at(i))
        .collect::<std::result::Result<Vec<_>, Error>>()?;

    for &t_out in &grid.times {
        let span = t_out - t_prev;
        let n_sub = ((span / opts.dt_max).ceil() as usize).max(1);
        let dt = span / n_sub as f64;
        if dt < opts.dt_min {
            return Err(Error::BadOption("output spacing needs a step below dt_min"));
        }
        // dt is constant across the substeps of this span, so build the
        // theta-step matrix once per span (matches the recombination paths).
        let m_sub: Vec<f64> = sys.sub.iter().map(|v| -dt * theta * v).collect();
        let m_diag: Vec<f64> = sys.diag.iter().map(|v| 1.0 - dt * theta * v).collect();
        let m_sup: Vec<f64> = sys.sup.iter().map(|v| -dt * theta * v).collect();
        for _ in 0..n_sub {
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
            if ntraps_per_cell.iter().all(|&m| m == 0) {
                c = nucleide_linalg::tridiag::solve(&m_sub, &m_diag, &m_sup, &e)
                    .map_err(crate::solve::tridiag_err)?;
            } else {
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
                        .map_err(crate::solve::tridiag_err)?;
                    let mut ct_next = ct.clone();
                    for i in 0..n {
                        for j in 0..ct_next[i].len() {
                            let (k, p) = rates[i][j];
                            ct_next[i][j] = trap_update(
                                ct[i][j],
                                c_next[i],
                                k,
                                p,
                                stack.layers
                                    [stack.layer_of(i).expect("i < total_cells by construction")]
                                .traps[j]
                                    .site_density,
                                dt,
                            );
                        }
                    }
                    let mut err = 0.0_f64;
                    for i in 0..n {
                        let scale = opts.atol + opts.rtol * c_next[i].abs().max(c_iter[i].abs());
                        err = err.max((c_next[i] - c_iter[i]).abs() / scale);
                        for j in 0..ct_next[i].len() {
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
        flux_left.push(outward_flux(left, c[0], sys.d_cell[0], sys.dx[0]));
        flux_right.push(outward_flux(
            right,
            c[n - 1],
            sys.d_cell[n - 1],
            sys.dx[n - 1],
        ));
    }

    Ok(LayeredSolution {
        times: grid.times.clone(),
        mobile: mobile_out,
        trapped: trapped_out,
        flux_left,
        flux_right,
        dx: sys.dx.clone(),
        initial: initial.clone(),
    })
}

/// Layered transient with at least one recombination end.
///
/// Mirrors [`solve::solve_recombination`]: the recombination ends become
/// Dirichlet(0) faces in the base layered operator (the conductance stays in
/// the matrix; the face source re-enters through the sensitivity columns and
/// the explicit old-face term), one sensitivity column `P = M⁻¹ Q` per
/// recombination end is built per `dt` span, and the faces close per step in
/// closed form (one end) or by the analytic-Jacobian Newton (two ends),
/// fused into the trap Picard loop with the shared `rtol`/`atol`.
#[allow(clippy::too_many_arguments)]
fn solve_layers_recombination(
    stack: &LayerStack,
    left: &Boundary,
    right: &Boundary,
    grid: &TimeGrid,
    initial: &InitialState,
    opts: &SolverOptions,
    kr_left: Option<f64>,
    kr_right: Option<f64>,
) -> Result<LayeredSolution, Error> {
    let n = stack.total_cells();
    let dx = stack.cell_widths();
    let d_cell = stack.diffusivities()?;
    let theta = opts.theta.value();
    let g_left = 2.0 * d_cell[0] / dx[0];
    let g_right = 2.0 * d_cell[n - 1] / dx[n - 1];
    let bl0 = kr_left.map_or_else(|| left.clone(), |_| Boundary::Dirichlet(0.0));
    let br0 = kr_right.map_or_else(|| right.clone(), |_| Boundary::Dirichlet(0.0));
    let sys = assemble_layers(stack, &bl0, &br0)?;

    let rates: Vec<Vec<(f64, f64)>> = (0..n)
        .map(|i| stack.trap_rates_at(i))
        .collect::<std::result::Result<Vec<_>, Error>>()?;
    let any_traps = rates.iter().any(|r| !r.is_empty());

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
        let m_sub: Vec<f64> = sys.sub.iter().map(|v| -dt * theta * v).collect();
        let m_diag: Vec<f64> = sys.diag.iter().map(|v| 1.0 - dt * theta * v).collect();
        let m_sup: Vec<f64> = sys.sup.iter().map(|v| -dt * theta * v).collect();
        let sens_left = kr_left
            .map(|_| {
                let mut v = vec![0.0; n];
                v[0] = dt * theta * g_left / dx[0];
                nucleide_linalg::tridiag::solve(&m_sub, &m_diag, &m_sup, &v)
                    .map_err(crate::solve::tridiag_err)
            })
            .transpose()?;
        let sens_right = kr_right
            .map(|_| {
                let mut v = vec![0.0; n];
                v[n - 1] = dt * theta * g_right / dx[n - 1];
                nucleide_linalg::tridiag::solve(&m_sub, &m_diag, &m_sup, &v)
                    .map_err(crate::solve::tridiag_err)
            })
            .transpose()?;
        for _ in 0..n_sub {
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
                e[0] += dt * (1.0 - theta) * g_left / dx[0] * cf;
            }
            if let Some(cf) = cf_right {
                e[n - 1] += dt * (1.0 - theta) * g_right / dx[n - 1] * cf;
            }
            if !any_traps {
                let base = nucleide_linalg::tridiag::solve(&m_sub, &m_diag, &m_sup, &e)
                    .map_err(crate::solve::tridiag_err)?;
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
                        .map_err(crate::solve::tridiag_err)?;
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
                        for j in 0..ct_next[i].len() {
                            let (k, p) = rates[i][j];
                            ct_next[i][j] = trap_update(
                                ct[i][j],
                                c_next[i],
                                k,
                                p,
                                stack.layers
                                    [stack.layer_of(i).expect("i < total_cells by construction")]
                                .traps[j]
                                    .site_density,
                                dt,
                            );
                        }
                    }
                    let mut err = 0.0_f64;
                    for i in 0..n {
                        let scale = opts.atol + opts.rtol * c_next[i].abs().max(c_iter[i].abs());
                        err = err.max((c_next[i] - c_iter[i]).abs() / scale);
                        for j in 0..ct_next[i].len() {
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
            || outward_flux(left, c[0], sys.d_cell[0], sys.dx[0]),
            |cf| kr_left.expect("rate present for a closed left face") * cf * cf,
        ));
        flux_right.push(cf_right.map_or_else(
            || outward_flux(right, c[n - 1], sys.d_cell[n - 1], sys.dx[n - 1]),
            |cf| kr_right.expect("rate present for a closed right face") * cf * cf,
        ));
    }

    Ok(LayeredSolution {
        times: grid.times.clone(),
        mobile: mobile_out,
        trapped: trapped_out,
        flux_left,
        flux_right,
        dx: sys.dx.clone(),
        initial: initial.clone(),
    })
}

// ---------------------------------------------------------------------------
// Gates: multi-layer series stacks (G7/G8)
// ---------------------------------------------------------------------------
//
// Provenance: every number below is synthetic — hand-built layer stacks with
// hand-derived closed forms (series-resistance network over per-layer
// permeabilities Φ = D·K_S; recombination outlet closes the face quadratic
// against the same resistance). No evaluated-library data, consistent with
// the landed-kernel analytic-gate stance.

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::tests::no_traps;

    /// Gate stack A (2-layer): L₁ = L₂ = 5e-4 m, D₁ = 1e-9, D₂ = 5e-10 m²/s,
    /// K_S₁ = 2.0, K_S₂ = 0.5 mol/m³/Pa¹ᐟ², 128 + 128 cells, 500 K.
    /// Φ₁ = 2e-9, Φ₂ = 2.5e-10 mol/m/s/Pa¹ᐟ².
    fn gate_stack_a() -> LayerStack {
        LayerStack::new(
            vec![
                LayerSpec::new(5e-4, 128, 1e-9, 0.0, 2.0, vec![], vec![500.0], vec![]).unwrap(),
                LayerSpec::new(5e-4, 128, 5e-10, 0.0, 0.5, vec![], vec![500.0], vec![]).unwrap(),
            ],
            vec![Interface::Sieverts],
        )
        .unwrap()
    }

    /// Series resistance of stack A: L₁/Φ₁ + L₂/Φ₂ = 2.5e5 + 2.0e6.
    const STACK_A_R: f64 = 2.25e6;

    #[test]
    fn layer_offsets_resolve_every_cell() {
        let stack = gate_stack_a();
        assert_eq!(stack.layer_start(0), 0);
        assert_eq!(stack.layer_start(1), 128);
        for i in 0..128 {
            assert_eq!(stack.layer_of(i), Some(0));
        }
        for i in 128..256 {
            assert_eq!(stack.layer_of(i), Some(1));
        }
        assert_eq!(stack.layer_of(256), None);
        assert!(matches!(
            stack.temp_at(256),
            Err(Error::BadCellIndex {
                index: 256,
                total: 256
            })
        ));
    }

    /// Closed-form steady flux for Dirichlet(c0)|Dirichlet(cN): J = Δu/R.
    fn series_flux(stack: &LayerStack, c0: f64, cn: f64, left: &Boundary, right: &Boundary) -> f64 {
        let _ = (left, right);
        let u0 = c0 / stack.layers[0].solubility;
        let un = cn / stack.layers[stack.layers.len() - 1].solubility;
        let r: f64 = stack
            .layers
            .iter()
            .map(|l| l.thickness / (l.d0 * l.solubility))
            .sum();
        (u0 - un) / r
    }

    /// Interior face fluxes of a converged layered steady profile: face `i`
    /// (between cells `i−1`, `i`) carries `J = (u_{i−1} − u_i)/R_f`.
    fn interior_face_fluxes(stack: &LayerStack, s: &LayeredSteadyState) -> Vec<f64> {
        let n = stack.total_cells();
        let d = stack.diffusivities().unwrap();
        let sol: Vec<f64> = (0..n).map(|i| stack.solubility_at(i).unwrap()).collect();
        let dx = stack.cell_widths();
        let phi: Vec<f64> = d.iter().zip(&sol).map(|(a, b)| a * b).collect();
        (1..n)
            .map(|i| {
                let r = dx[i - 1] / (2.0 * phi[i - 1]) + dx[i] / (2.0 * phi[i]);
                (s.mobile[i - 1] / sol[i - 1] - s.mobile[i] / sol[i]) / r
            })
            .collect()
    }

    /// Sieverts potential at interface `iface` from each side of the face,
    /// given the converged face flux `j`: the adjacent cell centre sits half
    /// a cell-width away, so `u_face = u_adj ∓ j·R_half`.
    fn interface_u_sides(
        stack: &LayerStack,
        s: &LayeredSteadyState,
        iface: usize,
        j: f64,
    ) -> (f64, f64) {
        let k = iface;
        let il = stack.layer_start(k) + stack.layers[k].cells - 1;
        let ir = il + 1;
        let d = stack.diffusivities().unwrap();
        let s_l = stack.layers[k].solubility;
        let s_r = stack.layers[k + 1].solubility;
        let dx = stack.cell_widths();
        let r_half_l = dx[il] / (2.0 * d[il] * s_l);
        let r_half_r = dx[ir] / (2.0 * d[ir] * s_r);
        let u_left = s.mobile[il] / s_l - j * r_half_l;
        let u_right = s.mobile[ir] / s_r + j * r_half_r;
        (u_left, u_right)
    }

    #[test]
    fn g7a_two_layer_series_resistance_steady() {
        // G7a: 2-layer stack, Dirichlet(1.0)|Dirichlet(0.0). Closed form
        // u₀ = 0.5, R = 2.25e6 → J = 2.2222...e-7; interface potential
        // u₁ = u₀ − J·L₁/Φ₁ = 0.4444...; the discrete profile is the exact
        // piecewise-linear u-profile at the nodes, so 1e-12 pins are pure
        // roundoff (same tolerance class as G5a).
        let stack = gate_stack_a();
        let s = steady_layers(
            &stack,
            &Boundary::dirichlet(1.0).unwrap(),
            &Boundary::dirichlet(0.0).unwrap(),
        )
        .unwrap();
        let j = series_flux(
            &stack,
            1.0,
            0.0,
            &Boundary::dirichlet(1.0).unwrap(),
            &Boundary::dirichlet(0.0).unwrap(),
        );
        assert!(
            (s.flux_right - j).abs() <= 1e-12 * j + 1e-18,
            "{} vs {j}",
            s.flux_right
        );
        assert!((s.flux_left + j).abs() <= 1e-12 * j + 1e-18);
        // Interface flux continuity: every interior face carries the same J
        // to roundoff (the recombination-free steady rows enforce this
        // exactly; the pin guards the assembly path).
        for (k, jf) in interior_face_fluxes(&stack, &s).iter().enumerate() {
            assert!(
                (jf - j).abs() <= 1e-12 * j + 1e-18,
                "face {}: {jf} vs {j}",
                k + 1
            );
        }
        // Profile: piecewise-linear u(x) sampled at the cell centres.
        let u0 = 1.0 / 2.0;
        let u1 = u0 - j * 5e-4 / (1e-9 * 2.0);
        for (i, &c) in s.mobile.iter().enumerate() {
            let x = s.centres[i];
            let u_exact = if x <= 5e-4 {
                u0 - j * x / (1e-9 * 2.0)
            } else {
                u1 - j * (x - 5e-4) / (5e-10 * 0.5)
            };
            let c_exact = if x <= 5e-4 {
                2.0 * u_exact
            } else {
                0.5 * u_exact
            };
            assert!((c - c_exact).abs() <= 1e-12, "cell {i}: {c} vs {c_exact}");
        }
        // Sieverts interface condition: u at the interface from each side
        // (half-cell-corrected) is continuous across the jump in c and hits
        // the closed-form interface potential.
        let (u_left, u_right) = interface_u_sides(&stack, &s, 0, j);
        assert!((u_left - u_right).abs() <= 1e-12, "{u_left} vs {u_right}");
        assert!((u_left - u1).abs() <= 1e-12);
        assert!(
            (s.inventory_mobile - s.mobile.iter().zip(&s.dx).map(|(c, w)| c * w).sum::<f64>())
                .abs()
                < 1e-18
        );
    }

    #[test]
    fn g7b_three_layer_interface_flux_continuity() {
        // G7b: 3-layer stack (N > 2 exercises the per-gap bookkeeping):
        // L = 3e-4 m each, D = 1e-9/4e-10/2.5e-10 m²/s, K_S = 1.0/0.8/0.6,
        // Dirichlet(1.2)|Dirichlet(0.1). Closed form over Φ = D·K_S.
        let stack = LayerStack::new(
            vec![
                LayerSpec::new(3e-4, 64, 1e-9, 0.0, 1.0, vec![], vec![500.0], vec![]).unwrap(),
                LayerSpec::new(3e-4, 64, 4e-10, 0.0, 0.8, vec![], vec![500.0], vec![]).unwrap(),
                LayerSpec::new(3e-4, 64, 2.5e-10, 0.0, 0.6, vec![], vec![500.0], vec![]).unwrap(),
            ],
            vec![Interface::Sieverts, Interface::Sieverts],
        )
        .unwrap();
        let s = steady_layers(
            &stack,
            &Boundary::dirichlet(1.2).unwrap(),
            &Boundary::dirichlet(0.1).unwrap(),
        )
        .unwrap();
        let j = series_flux(
            &stack,
            1.2,
            0.1,
            &Boundary::dirichlet(1.2).unwrap(),
            &Boundary::dirichlet(0.1).unwrap(),
        );
        assert!(
            (s.flux_right - j).abs() <= 1e-12 * j + 1e-18,
            "{} vs {j}",
            s.flux_right
        );
        assert!((s.flux_left + j).abs() <= 1e-12 * j + 1e-18);
        for (k, jf) in interior_face_fluxes(&stack, &s).iter().enumerate() {
            assert!(
                (jf - j).abs() <= 1e-12 * j + 1e-18,
                "face {}: {jf} vs {j}",
                k + 1
            );
        }
        // u continuity at both interfaces, from each side (half-cell
        // corrected), vs the closed form u(xₖ) = u₀ − J·Σ_{j≤k} L_j/Φ_j.
        let u0 = 1.2;
        let u1 = u0 - j * 3e-4 / (1e-9 * 1.0);
        let u2 = u1 - j * 3e-4 / (4e-10 * 0.8);
        for (iface, want) in [(0_usize, u1), (1, u2)] {
            let (ul, ur) = interface_u_sides(&stack, &s, iface, j);
            assert!((ul - ur).abs() <= 1e-12, "interface {iface}: {ul} vs {ur}");
            assert!(
                (ul - want).abs() <= 1e-12,
                "interface {iface}: {ul} vs {want}"
            );
        }
    }

    #[test]
    fn g7c_single_layer_recovers_landed_kernel_exactly() {
        // G7c (regression anchor): a 1-layer stack dispatches to the landed
        // single-slab kernel, so the results are identical by construction —
        // asserted field-by-field, steady and transient.
        let layer = LayerSpec::new(1e-3, 64, 1e-9, 0.0, 3.0, vec![], vec![500.0], vec![]).unwrap();
        let stack = LayerStack::new(vec![layer.clone()], vec![]).unwrap();
        let params = landed_params(&layer).unwrap();
        let left = Boundary::dirichlet(1.0).unwrap();
        let right = Boundary::dirichlet(0.0).unwrap();
        let landed = solve::steady_state(&params, &left, &right).unwrap();
        let layered = steady_layers(&stack, &left, &right).unwrap();
        assert_eq!(landed.mobile, layered.mobile);
        assert_eq!(landed.trapped, layered.trapped);
        assert_eq!(landed.flux_left, layered.flux_left);
        assert_eq!(landed.flux_right, layered.flux_right);
        assert_eq!(landed.inventory_mobile, layered.inventory_mobile);
        assert_eq!(landed.inventory_trapped, layered.inventory_trapped);
        assert_eq!(landed.centres, layered.centres);
        // Transient, including a trap species so the per-layer trap
        // bookkeeping is on the dispatch path too.
        let layer_t = LayerSpec::new(
            1e-3,
            64,
            1e-9,
            0.0,
            3.0,
            vec![TrapSpec::new(0.05, 0.0, 0.01, 0.0, 2.0).unwrap()],
            vec![500.0],
            vec![],
        )
        .unwrap();
        let stack_t = LayerStack::new(vec![layer_t.clone()], vec![]).unwrap();
        let params_t = landed_params(&layer_t).unwrap();
        let grid = TimeGrid::new(vec![10.0, 100.0]).unwrap();
        let init = stack_t.zero_state();
        let opts = SolverOptions {
            rtol: 1e-10,
            atol: 1e-14,
            dt_max: 1.0,
            ..Default::default()
        };
        let landed_t = solve::solve(&params_t, &left, &right, &grid, &init, &opts).unwrap();
        let layered_t = solve_layers(&stack_t, &left, &right, &grid, &init, &opts).unwrap();
        assert_eq!(landed_t.mobile, layered_t.mobile);
        assert_eq!(landed_t.trapped, layered_t.trapped);
        assert_eq!(landed_t.flux_left, layered_t.flux_left);
        assert_eq!(landed_t.flux_right, layered_t.flux_right);
        assert_eq!(landed_t.times, layered_t.times);
        assert_eq!(vec![landed_t.dx; 64], layered_t.dx);
    }

    #[test]
    fn g7d_property_continuous_stack_matches_landed() {
        // G7d: a 2-layer stack with identical D and K_S in both layers is a
        // plain slab split at a transparent interface (c/K_S continuous with
        // equal K_S is c continuous; the face resistance halves to the
        // landed arithmetic mean). The layered path itself (not the 1-layer
        // dispatch) must then match the landed kernel to roundoff.
        let stack = LayerStack::new(
            vec![
                LayerSpec::new(5e-4, 32, 1e-9, 0.0, 1.0, vec![], vec![500.0], vec![]).unwrap(),
                LayerSpec::new(5e-4, 32, 1e-9, 0.0, 1.0, vec![], vec![500.0], vec![]).unwrap(),
            ],
            vec![Interface::Sieverts],
        )
        .unwrap();
        let slab = no_traps(1e-3, 64, 1e-9);
        let left = Boundary::dirichlet(1.0).unwrap();
        let right = Boundary::dirichlet(0.0).unwrap();
        let landed = solve::steady_state(&slab, &left, &right).unwrap();
        let layered = steady_layers(&stack, &left, &right).unwrap();
        assert!((layered.flux_right - landed.flux_right).abs() <= 1e-12 * landed.flux_right);
        for (a, b) in layered.mobile.iter().zip(&landed.mobile) {
            assert!((a - b).abs() <= 1e-12, "{a} vs {b}");
        }
        let grid = TimeGrid::new(vec![100.0, 1000.0]).unwrap();
        let init = stack.zero_state();
        let opts = SolverOptions {
            rtol: 1e-10,
            atol: 1e-14,
            dt_max: 1.0,
            ..Default::default()
        };
        let landed_t = solve::solve(&slab, &left, &right, &grid, &init, &opts).unwrap();
        let layered_t = solve_layers(&stack, &left, &right, &grid, &init, &opts).unwrap();
        for (a, b) in layered_t.flux_right.iter().zip(&landed_t.flux_right) {
            assert!((a - b).abs() <= 1e-9 * b.abs().max(1e-12), "{a} vs {b}");
        }
    }

    #[test]
    fn g7e_layered_steady_recombination_outer_end() {
        // G7e: recombination outer ends are reused on the layered system.
        // Stack A with a recombination outlet closes in closed form: the
        // stack is the series resistance R between u₀ and the outlet face
        // potential u_L, and J = K_r·(K_S₂·u_L)². Solve
        // K_r·K_S₂²·R·u_L² + u_L − u₀ = 0 for u_L (positive root), then
        // J = K_r·(K_S₂·u_L)² — the same quadratic structure as G5a with
        // the slab resistance replaced by the stack resistance.
        let stack = gate_stack_a();
        let kr = 1e-6;
        let r = STACK_A_R;
        let u0 = 0.5;
        let a = kr * 0.5_f64.powi(2) * r;
        let u_l = (-1.0 + (1.0 + 4.0 * a * u0).sqrt()) / (2.0 * a);
        let j = kr * (0.5 * u_l).powi(2);
        let s = steady_layers(
            &stack,
            &Boundary::dirichlet(1.0).unwrap(),
            &Boundary::recombination(kr).unwrap(),
        )
        .unwrap();
        assert!(
            (s.flux_right - j).abs() <= 1e-12 * j + 1e-18,
            "{} vs {j}",
            s.flux_right
        );
        assert!((s.flux_left + j).abs() <= 1e-12 * j + 1e-18);
        // Outlet face concentration: c_f = K_S₂·u_L = 0.5·u_L.
        let cf = 0.5 * u_l;
        // The last cell value sits half a cell uphill: c_adj = c_f + K_S₂·J·R_half
        // (the u-drop across the half cell scales by the solubility in c).
        let r_half = (5e-4 / 128.0) / (2.0 * 5e-10 * 0.5);
        let c_adj = cf + 0.5 * j * r_half;
        assert!(
            (s.mobile[255] - c_adj).abs() <= 1e-9 * c_adj,
            "{}",
            s.mobile[255]
        );
        for &c in &s.mobile {
            assert!(c >= 0.0);
        }
    }

    /// G8a gate stack: same transport data as stack A on a 32 + 32 grid, so
    /// the finest-mode scale matches the landed G6d band (t = 1200 s,
    /// dt ≤ 3 s): Crank–Nicolson ringing from the t = 0 corner kink is dead
    /// well before the pin time at these dt.
    fn gate_stack_a_coarse() -> LayerStack {
        LayerStack::new(
            vec![
                LayerSpec::new(5e-4, 32, 1e-9, 0.0, 2.0, vec![], vec![500.0], vec![]).unwrap(),
                LayerSpec::new(5e-4, 32, 5e-10, 0.0, 0.5, vec![], vec![500.0], vec![]).unwrap(),
            ],
            vec![Interface::Sieverts],
        )
        .unwrap()
    }

    /// Mobile profile of a gate stack at `t_end` with uniform steps of `dt_max`.
    fn stack_a_profile(theta: crate::solve::Theta, dt_max: f64, t_end: f64) -> Vec<f64> {
        let stack = gate_stack_a_coarse();
        solve_layers(
            &stack,
            &Boundary::dirichlet(1.0).unwrap(),
            &Boundary::dirichlet(0.0).unwrap(),
            &TimeGrid::new(vec![t_end]).unwrap(),
            &stack.zero_state(),
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
    fn g8a_dt_halving_order_preserved() {
        // G8a: dt-halving on the fixed 2-layer stack shows the θ-method
        // order (≈2 Crank–Nicolson, ≈1 backward Euler) on two successive
        // halvings, pinned in the resolved band exactly like G6d (the t = 0
        // corner kink excites Crank–Nicolson fine-mode ringing, so the band
        // is measured where truncation error dominates). The interface
        // fluxes sit inside the step matrix, so the orders must match the
        // landed single-slab bands.
        let t_end = 1200.0;
        let dt = [3.0_f64, 1.5, 0.75, 0.375];
        let orders = |theta: crate::solve::Theta| {
            let p0 = stack_a_profile(theta, dt[0], t_end);
            let p1 = stack_a_profile(theta, dt[1], t_end);
            let p2 = stack_a_profile(theta, dt[2], t_end);
            let p3 = stack_a_profile(theta, dt[3], t_end);
            [
                (max_diff(&p0, &p1) / max_diff(&p1, &p2)).log2(),
                (max_diff(&p1, &p2) / max_diff(&p2, &p3)).log2(),
            ]
        };
        for order in orders(crate::solve::Theta::CrankNicolson) {
            assert!((1.5..=2.5).contains(&order), "Crank–Nicolson order {order}");
        }
        for order in orders(crate::solve::Theta::BackwardEuler) {
            assert!((0.7..=1.3).contains(&order), "backward Euler order {order}");
        }
    }

    /// Chained layered transient rows at `t_end`: resolved steps (dt = 0.5)
    /// to `t_a` kill the t = 0 corner-kink ringing, then coarse steps carry
    /// the smooth state to `t_end` (the discrete steady state is a fixed
    /// point of the θ-step, so only the decayed transient modes matter).
    fn stack_a_chained(
        right: &Boundary,
        t_a: f64,
        t_end: f64,
        dt_late: f64,
    ) -> (Vec<f64>, f64, f64) {
        let stack = gate_stack_a();
        let left = Boundary::dirichlet(1.0).unwrap();
        let seg_a = solve_layers(
            &stack,
            &left,
            right,
            &TimeGrid::new(vec![t_a]).unwrap(),
            &stack.zero_state(),
            &SolverOptions {
                dt_max: 0.5,
                rtol: 1e-10,
                atol: 1e-14,
                ..Default::default()
            },
        )
        .unwrap();
        let init = InitialState {
            mobile: seg_a.mobile[0].clone(),
            trapped: seg_a.trapped[0].clone(),
        };
        let seg_b = solve_layers(
            &stack,
            &left,
            right,
            &TimeGrid::new(vec![t_end - t_a]).unwrap(),
            &init,
            &SolverOptions {
                dt_max: dt_late,
                rtol: 1e-10,
                atol: 1e-14,
                ..Default::default()
            },
        )
        .unwrap();
        (
            seg_b.mobile[0].clone(),
            seg_b.flux_left[0],
            seg_b.flux_right[0],
        )
    }

    #[test]
    fn g8b_late_time_asymptote_to_layered_steady() {
        // G8b: the 2-layer transient from a clean slab lands on the layered
        // steady state at t = 60·maxₖ(Lₖ²/6Dₖ) = 5000 s (the stack lag
        // scale; the coupled slowest mode decays faster than the slowest
        // single layer, so e⁻⁶⁰ leaves only the asymptote). Both a
        // Dirichlet and a recombination outlet are exercised (the latter
        // closes the layered face construction per step).
        let t_end = 5000.0;
        // Dirichlet outlet: flux against the series closed form.
        let stack = gate_stack_a();
        let j_dir = 0.5 / STACK_A_R;
        let (prof_d, fl_d, fr_d) =
            stack_a_chained(&Boundary::dirichlet(0.0).unwrap(), 300.0, t_end, 5.0);
        assert!((fr_d - j_dir).abs() <= 1e-6 * j_dir, "{fr_d} vs {j_dir}");
        assert!((fl_d + fr_d).abs() <= 1e-6 * j_dir);
        let steady = steady_layers(
            &stack,
            &Boundary::dirichlet(1.0).unwrap(),
            &Boundary::dirichlet(0.0).unwrap(),
        )
        .unwrap();
        assert!(max_diff(&prof_d, &steady.mobile) <= 1e-9);
        // Recombination outlet: flux against the G7e closed form and the
        // profile against the layered steady recombination state.
        let kr = 1e-6;
        let a = kr * 0.5_f64.powi(2) * STACK_A_R;
        let u_l = (-1.0 + (1.0 + 4.0 * a * 0.5).sqrt()) / (2.0 * a);
        let j_rec = kr * (0.5 * u_l).powi(2);
        let right = Boundary::recombination(kr).unwrap();
        let (prof_r, fl_r, fr_r) = stack_a_chained(&right, 300.0, t_end, 5.0);
        assert!((fr_r - j_rec).abs() <= 1e-6 * j_rec, "{fr_r} vs {j_rec}");
        assert!((fl_r + fr_r).abs() <= 1e-6 * j_rec);
        let steady_r = steady_layers(&stack, &Boundary::dirichlet(1.0).unwrap(), &right).unwrap();
        assert!(max_diff(&prof_r, &steady_r.mobile) <= 1e-9);
        assert!(prof_d.iter().all(|&c| c >= 0.0));
        assert!(prof_r.iter().all(|&c| c >= 0.0));
    }

    #[test]
    fn g8c_theta_mass_balance_layered() {
        // G8c: trap-free 2-layer stack, one θ-step per output interval, so
        // the discrete balance Iᵏ − Iᵏ⁻¹ = dt·[θRᵏ + (1−θ)Rᵏ⁻¹] with
        // R = −(F_L + F_R) closes against the reported boundary fluxes to
        // roundoff. Internal interface fluxes cancel in the sum (they are
        // internal), so the landed G6c identity carries over unchanged.
        let stack = gate_stack_a();
        let dx = stack.cell_widths();
        let left = Boundary::dirichlet(1.0).unwrap();
        let right = Boundary::dirichlet(0.0).unwrap();
        let opts = SolverOptions {
            dt_max: 20.0,
            rtol: 1e-10,
            atol: 1e-14,
            ..Default::default()
        };
        let sol = solve_layers(
            &stack,
            &left,
            &right,
            &TimeGrid::new(vec![20.0, 40.0]).unwrap(),
            &stack.zero_state(),
            &opts,
        )
        .unwrap();
        let inv: Vec<f64> = sol
            .mobile
            .iter()
            .map(|row| row.iter().zip(&dx).map(|(c, w)| c * w).sum::<f64>())
            .collect();
        let r_now = |k: usize| -(sol.flux_left[k] + sol.flux_right[k]);
        let (mut r_prev, mut i_prev) = (2.0 * 1e-9 * 1.0 / dx[0], 0.0_f64);
        for (k, dt) in [20.0, 20.0].iter().enumerate() {
            let want = i_prev + dt * (0.5 * r_prev + 0.5 * r_now(k));
            assert!(
                (inv[k] - want).abs() < 1e-12 * want.abs().max(1e-300),
                "row {k}: {} vs {want}",
                inv[k]
            );
            r_prev = r_now(k);
            i_prev = inv[k];
        }
    }

    #[test]
    fn g8_positivity_and_per_layer_traps() {
        // Positivity on the layered transient with *different* trap species
        // per layer (layer 1 trap-free, layer 2 one species): the per-cell
        // trap bookkeeping and the interface coupling stay nonneg.
        let stack = LayerStack::new(
            vec![
                LayerSpec::new(5e-4, 32, 1e-9, 0.0, 2.0, vec![], vec![500.0], vec![]).unwrap(),
                LayerSpec::new(
                    5e-4,
                    32,
                    5e-10,
                    0.0,
                    0.5,
                    vec![TrapSpec::new(0.05, 0.0, 0.01, 0.0, 2.0).unwrap()],
                    vec![500.0],
                    vec![],
                )
                .unwrap(),
            ],
            vec![Interface::Sieverts],
        )
        .unwrap();
        let init = InitialState {
            mobile: vec![0.0; 64],
            trapped: {
                let mut v: Vec<Vec<f64>> = Vec::new();
                v.extend(std::iter::repeat_n(vec![], 32));
                v.extend(std::iter::repeat_n(vec![0.0], 32));
                v
            },
        };
        let sol = solve_layers(
            &stack,
            &Boundary::dirichlet(1.0).unwrap(),
            &Boundary::dirichlet(0.0).unwrap(),
            &TimeGrid::new(vec![10.0, 100.0, 1000.0]).unwrap(),
            &init,
            &SolverOptions {
                dt_max: 1.0,
                ..Default::default()
            },
        )
        .unwrap();
        for row in &sol.mobile {
            assert!(row.iter().all(|&c| c >= 0.0));
        }
        for rows in &sol.trapped {
            for row in rows {
                assert!(row.iter().all(|&ct| (0.0..=2.0).contains(&ct)));
            }
        }
        assert!(sol.flux_right.iter().all(|&f| f >= 0.0));
        // Inventories are monotone non-decreasing while filling from the
        // left end (inflow ≥ 0, outflow ≥ 0, trapping is a sink).
        let inv = sol.inventory_total();
        for w in inv.windows(2) {
            assert!(w[1] >= w[0] - 1e-12);
        }
    }

    // -----------------------------------------------------------------------
    // Gates: Henry internal interfaces (G9)
    // -----------------------------------------------------------------------
    //
    // Provenance: synthetic stacks with hand-derived series-resistance
    // oracles — the same resistance form as G7/G8, with the layer
    // `solubility` playing `K_H` (Henry gap) or `K_S` (Sieverts gap) per the
    // adjacent interface's law. No evaluated-library data, consistent with
    // the analytic-gate stance.

    /// Henry gate stack (2-layer): L₁ = L₂ = 4e-4 m, D₁ = 2e-9,
    /// D₂ = 5e-10 m²/s, K_H₁ = 1.5, K_H₂ = 0.75 mol/m³/Pa, 96 + 96 cells,
    /// 500 K. Φ₁ = 3e-9, Φ₂ = 3.75e-10 mol/m/s/Pa.
    fn gate_stack_henry() -> LayerStack {
        LayerStack::new(
            vec![
                LayerSpec::new(4e-4, 96, 2e-9, 0.0, 1.5, vec![], vec![500.0], vec![]).unwrap(),
                LayerSpec::new(4e-4, 96, 5e-10, 0.0, 0.75, vec![], vec![500.0], vec![]).unwrap(),
            ],
            vec![Interface::Henry],
        )
        .unwrap()
    }

    #[test]
    fn g9a_henry_series_resistance_steady() {
        // G9a: 2-layer stack joined by a Henry interface, Dirichlet(1.2)|
        // Dirichlet(0.2). Closed form u₀ = 1.2/1.5 = 0.8, u_L = 0.2/0.75,
        // R = 1.2e6 → J = 4.4444...e-7; interface potential
        // u₁ = u₀ − J·L₁/Φ₁ = 0.7407...; the same (L2) series-resistance
        // form as G7a with K_H in place of K_S, pinned at G5-class
        // tolerances.
        let stack = gate_stack_henry();
        let left = Boundary::dirichlet(1.2).unwrap();
        let right = Boundary::dirichlet(0.2).unwrap();
        let s = steady_layers(&stack, &left, &right).unwrap();
        let j = series_flux(&stack, 1.2, 0.2, &left, &right);
        assert!(
            (s.flux_right - j).abs() <= 1e-12 * j + 1e-18,
            "{} vs {j}",
            s.flux_right
        );
        assert!((s.flux_left + j).abs() <= 1e-12 * j + 1e-18);
        // Interface flux continuity: every interior face — including the
        // Henry gap — carries the same J to roundoff (the resistance form
        // is law-agnostic, so the G7a pin carries over unchanged).
        for (k, jf) in interior_face_fluxes(&stack, &s).iter().enumerate() {
            assert!(
                (jf - j).abs() <= 1e-12 * j + 1e-18,
                "face {}: {jf} vs {j}",
                k + 1
            );
        }
        // Profile: piecewise-linear Henry-potential u(x) at the cell centres.
        let u0 = 1.2 / 1.5;
        let u1 = u0 - j * 4e-4 / (2e-9 * 1.5);
        for (i, &c) in s.mobile.iter().enumerate() {
            let x = s.centres[i];
            let u_exact = if x <= 4e-4 {
                u0 - j * x / (2e-9 * 1.5)
            } else {
                u1 - j * (x - 4e-4) / (5e-10 * 0.75)
            };
            let c_exact = if x <= 4e-4 {
                1.5 * u_exact
            } else {
                0.75 * u_exact
            };
            assert!((c - c_exact).abs() <= 1e-12, "cell {i}: {c} vs {c_exact}");
        }
        // Henry interface condition: the potential u = c/K_H at the
        // interface from each side (half-cell-corrected) is continuous
        // across the jump in c and hits the closed-form interface potential.
        let (u_left, u_right) = interface_u_sides(&stack, &s, 0, j);
        assert!((u_left - u_right).abs() <= 1e-12, "{u_left} vs {u_right}");
        assert!((u_left - u1).abs() <= 1e-12);
        assert!(s.mobile.iter().all(|&c| c >= 0.0));
    }

    /// Mixed-law gate stack (3-layer): L = 3e-4 m each, D = 1e-9 / 4e-10 /
    /// 2.5e-10 m²/s, K = 1.0 / 0.8 / 0.6 (K_S at the Sieverts gap, K_H at
    /// the Henry gap), 64 cells per layer. Φ = 1e-9, 3.2e-10, 1.5e-10
    /// mol/m/s/Pa.
    fn gate_stack_mixed() -> LayerStack {
        LayerStack::new(
            vec![
                LayerSpec::new(3e-4, 64, 1e-9, 0.0, 1.0, vec![], vec![500.0], vec![]).unwrap(),
                LayerSpec::new(3e-4, 64, 4e-10, 0.0, 0.8, vec![], vec![500.0], vec![]).unwrap(),
                LayerSpec::new(3e-4, 64, 2.5e-10, 0.0, 0.6, vec![], vec![500.0], vec![]).unwrap(),
            ],
            vec![Interface::Sieverts, Interface::Henry],
        )
        .unwrap()
    }

    #[test]
    fn g9b_mixed_interface_laws_flux_continuity() {
        // G9b: one interface of each law on the same stack. Both are linear
        // local-equilibrium conditions with the same resistance form, so the
        // (L2) closed form holds over the mixed stack and every interior
        // face — across the Sieverts gap *and* the Henry gap — carries the
        // same flux to roundoff.
        let stack = gate_stack_mixed();
        let left = Boundary::dirichlet(1.0).unwrap();
        let right = Boundary::dirichlet(0.0).unwrap();
        let s = steady_layers(&stack, &left, &right).unwrap();
        let j = series_flux(&stack, 1.0, 0.0, &left, &right);
        assert!(
            (s.flux_right - j).abs() <= 1e-12 * j + 1e-18,
            "{} vs {j}",
            s.flux_right
        );
        assert!((s.flux_left + j).abs() <= 1e-12 * j + 1e-18);
        for (k, jf) in interior_face_fluxes(&stack, &s).iter().enumerate() {
            assert!(
                (jf - j).abs() <= 1e-12 * j + 1e-18,
                "face {}: {jf} vs {j}",
                k + 1
            );
        }
        // u continuity at both interfaces, from each side (half-cell
        // corrected), vs the closed form u(xₖ) = u₀ − J·Σ_{j≤k} L_j/Φ_j —
        // the potential law switches with the interface, the resistance
        // form does not.
        let u0 = 1.0_f64;
        let u1 = u0 - j * 3e-4 / (1e-9 * 1.0);
        let u2 = u1 - j * 3e-4 / (4e-10 * 0.8);
        for (iface, want) in [(0_usize, u1), (1, u2)] {
            let (ul, ur) = interface_u_sides(&stack, &s, iface, j);
            assert!((ul - ur).abs() <= 1e-12, "interface {iface}: {ul} vs {ur}");
            assert!(
                (ul - want).abs() <= 1e-12,
                "interface {iface}: {ul} vs {want}"
            );
        }
        assert!(s.mobile.iter().all(|&c| c >= 0.0));
    }

    #[test]
    fn loud_errors_and_validation() {
        // Recombination internal interfaces stay loud named errors; the
        // linear Henry law constructs fine (G9).
        let layer = |d: f64, s: f64| {
            LayerSpec::new(5e-4, 8, d, 0.0, s, vec![], vec![500.0], vec![]).unwrap()
        };
        assert!(LayerStack::new(
            vec![layer(1e-9, 1.0), layer(1e-9, 1.0)],
            vec![Interface::Henry]
        )
        .is_ok());
        let err = LayerStack::new(
            vec![layer(1e-9, 1.0), layer(1e-9, 1.0)],
            vec![Interface::Recombination],
        )
        .unwrap_err();
        assert!(matches!(err, Error::UnsupportedInterface(_)), "{err}");
        // Shape errors.
        assert!(LayerStack::new(vec![], vec![]).is_err());
        assert!(LayerStack::new(vec![layer(1e-9, 1.0)], vec![Interface::Sieverts]).is_err());
        assert!(LayerStack::new(
            vec![layer(1e-9, 1.0), layer(1e-9, 1.0)],
            vec![Interface::Sieverts, Interface::Sieverts],
        )
        .is_err());
        // Layer-spec validation mirrors the landed params stance.
        assert!(LayerSpec::new(0.0, 8, 1e-9, 0.0, 1.0, vec![], vec![500.0], vec![]).is_err());
        assert!(LayerSpec::new(5e-4, 0, 1e-9, 0.0, 1.0, vec![], vec![500.0], vec![]).is_err());
        assert!(LayerSpec::new(5e-4, 8, 0.0, 0.0, 1.0, vec![], vec![500.0], vec![]).is_err());
        assert!(LayerSpec::new(5e-4, 8, 1e-9, 0.0, 0.0, vec![], vec![500.0], vec![]).is_err());
        assert!(
            LayerSpec::new(5e-4, 8, 1e-9, 0.0, 1.0, vec![], vec![500.0, 600.0], vec![]).is_err()
        );
        // Initial-state validation against the owning layer's trap count.
        let stack = LayerStack::new(
            vec![
                LayerSpec::new(5e-4, 8, 1e-9, 0.0, 1.0, vec![], vec![500.0], vec![]).unwrap(),
                LayerSpec::new(
                    5e-4,
                    8,
                    1e-9,
                    0.0,
                    1.0,
                    vec![TrapSpec::new(0.05, 0.0, 0.01, 0.0, 2.0).unwrap()],
                    vec![500.0],
                    vec![],
                )
                .unwrap(),
            ],
            vec![Interface::Sieverts],
        )
        .unwrap();
        let left = Boundary::dirichlet(1.0).unwrap();
        let right = Boundary::dirichlet(0.0).unwrap();
        let good = InitialState {
            mobile: vec![0.0; 16],
            trapped: {
                let mut v: Vec<Vec<f64>> = Vec::new();
                v.extend(std::iter::repeat_n(vec![], 8));
                v.extend(std::iter::repeat_n(vec![0.0], 8));
                v
            },
        };
        let grid = TimeGrid::new(vec![1.0]).unwrap();
        let opts = SolverOptions::default();
        assert!(solve_layers(&stack, &left, &right, &grid, &good, &opts).is_ok());
        let bad_mobile = InitialState {
            mobile: vec![0.0; 15],
            ..good.clone()
        };
        assert!(solve_layers(&stack, &left, &right, &grid, &bad_mobile, &opts).is_err());
        let mut bad_trapped = good.clone();
        bad_trapped.trapped[15] = vec![];
        assert!(solve_layers(&stack, &left, &right, &grid, &bad_trapped, &opts).is_err());
        let mut bad_load = good;
        bad_load.trapped[15] = vec![3.0];
        assert!(solve_layers(&stack, &left, &right, &grid, &bad_load, &opts).is_err());
        // dt_min above the output spacing is rejected, not stalled on.
        let opts = SolverOptions {
            dt_min: 2.0,
            ..Default::default()
        };
        assert!(solve_layers(&stack, &left, &right, &grid, &stack.zero_state(), &opts).is_err());
    }

    #[test]
    fn per_cell_getters_reject_out_of_range_index() {
        // The audit-confirmed panic path: `stack.temp_at(99)` on a 16-cell
        // stack (via the public getters) must be a named error, never a
        // panic.
        let stack = LayerStack::new(
            vec![
                LayerSpec::new(5e-4, 8, 1e-9, 0.0, 1.0, vec![], vec![500.0], vec![1.0]).unwrap(),
                LayerSpec::new(5e-4, 8, 1e-9, 0.0, 2.0, vec![], vec![600.0], vec![]).unwrap(),
            ],
            vec![Interface::Sieverts],
        )
        .unwrap();
        let want = Error::BadCellIndex {
            index: 99,
            total: 16,
        };
        assert_eq!(stack.temp_at(99).unwrap_err(), want);
        assert_eq!(stack.source_at(99).unwrap_err(), want);
        assert_eq!(stack.diffusivity_at(99).unwrap_err(), want);
        assert_eq!(stack.solubility_at(99).unwrap_err(), want);
        assert_eq!(stack.trap_rates_at(99).unwrap_err(), want);
        // In-range cells keep working, including across the layer boundary.
        assert_eq!(stack.temp_at(7).unwrap(), 500.0);
        assert_eq!(stack.temp_at(8).unwrap(), 600.0);
        assert_eq!(stack.source_at(0).unwrap(), 1.0);
        assert_eq!(stack.source_at(8).unwrap(), 0.0);
        assert_eq!(stack.solubility_at(15).unwrap(), 2.0);
    }
}
