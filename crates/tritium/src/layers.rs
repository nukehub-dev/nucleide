//! Multi-layer series stacks with linear interface conditions (G7–G9).
//!
//! A [`LayerStack`] chains `N ≥ 1` slabs — e.g. a W/Cu/CuCrZr first-wall
//! stack — each with its own thickness, diffusivity, solubility, trap
//! species, temperature, and source. Inside every layer the mobile balance
//! (T1) and the McNabb–Foster kinetics (T2) hold with the layer's own
//! coefficients; temperature stays caller-supplied per layer (no heat solve).
//!
//! Internal interfaces carry either a linear local-equilibrium law — the
//! Sieverts potential `u = c_m / K_S` \[Pa¹ᐟ²\] or the Henry potential
//! `u = c_m / K_H` \[Pa\] (the layer's [`LayerSpec::solubility`] plays
//! whichever constant the adjacent interface's [`Interface`] law declares)
//! — or the vented-sink recombination law (R-S, G10/G11): a single face
//! concentration `x ≥ 0` with half-cell fluxes `J_L = g_L (c_L − x)` and
//! `J_R = g_R (x − c_R)` (`g = 2D/dx` per side, no solubility involved) and
//! the desorption sink `K_r x²` venting from the face, so
//! `J_L − J_R = K_r x²`. The flux-continuous product law (`J = K_r x₀x₁`)
//! stays permanently rejected (spike outcome: spurious insulated root,
//! symmetry breaking, blocked transient) — there is no spelling for it, so
//! it cannot be constructed.
//! Because the supported linear conditions are linear in the mobile
//! concentration, every linear interior face flux
//!
//! ```text
//! J_f = (c_i / K_i − c_{i+1} / K_{i+1}) / R_f,
//! R_f = dx_i / (2 Φ_i) + dx_{i+1} / (2 Φ_{i+1}),   Φ = D · K,
//! ```
//!
//! folds directly into the tridiagonal θ-step matrix — the linear interfaces
//! need no Newton machinery of their own (strictly simpler than the
//! recombination faces of G5/G6, which are nonlinear in the face value). A
//! vented-sink gap instead CUTs the matrix (block-diagonal; each block sees
//! a Dirichlet cut face): the adjacent-cell values stay affine in the face
//! value, so the scalar residual `R(x) = P − Q·x − K_r·x² = 0` closes in
//! closed form (`x = 2P/(√(Q²+4K_rP)+Q)`; `P > 0` required, a loud error
//! otherwise). The outer ends reuse the landed [`Boundary`] taxonomy
//! unchanged, including recombination ends through the same affine
//! face-response construction (G5/G6), here built on the layered step matrix. A one-layer stack dispatches to the
//! landed single-slab kernel, so `N = 1` reproduces it exactly (G7c).

use crate::bc::Boundary;
use crate::error::Error;
use crate::params::{arrhenius, TransportParams, TrapSpec};
use crate::solve::{
    self, apply_faces, close_faces, face_closed_single, face_concentration, face_from_adjacent,
    face_newton_pair, linear_steady, outward_flux, trap_update, DiffusionSystem, InitialState,
    SolverOptions, TimeGrid, G5_ATOL, G5_NEWTON_MAX, G5_RTOL, MAX_PICARD,
};

/// Internal interface condition between adjacent layers.
///
/// The taxonomy mirrors [`crate::bc::Boundary`]. The linear laws
/// [`Sieverts`] and [`Henry`] are supported: the local-equilibrium condition
/// `u = c_m / K` continuous (Sieverts potential with `K = K_S`, Henry
/// potential with `K = K_H`; the layer [`LayerSpec::solubility`] carries the
/// matching constant) with continuous flux. [`Recombination`] is the
/// vented-sink law (R-S, G10/G11): a single face concentration `x ≥ 0` with
/// half-cell fluxes `J_L = g_L·(c_L − x)`, `J_R = g_R·(x − c_R)` (rightward
/// positive, `g = 2D/dx` per side) and the desorption sink `K_r·x²` venting
/// from the face, closing `J_L − J_R = K_r·x²` in closed form on the CUT
/// (block-diagonal) matrix. The `K_r → 0` limit is a continuous-concentration
/// joint, not a Sieverts law (no `K` ratio enters); a non-positive permeation
/// drive `P ≤ 0` has no nonnegative root and fails loudly at solve time.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Interface {
    /// Sieverts local equilibrium: `c_m / K_S` continuous, flux continuous.
    Sieverts,
    /// Linear (Henry) local equilibrium: `c_m / K_H` continuous, flux
    /// continuous. Shares the Sieverts resistance form (the layer
    /// [`LayerSpec::solubility`] plays `K_H`), so it folds into the
    /// θ-step matrix exactly like [`Interface::Sieverts`].
    Henry,
    /// Vented-sink recombination law (R-S): the gap desorbs `K_r·x²`
    /// \[mol/m²/s\] from the single face concentration `x` \[mol/m³\].
    Recombination {
        /// Desorption rate `K_r` \[m⁴/mol/s\] (`> 0`, finite).
        rate: f64,
    },
}

impl Interface {
    /// Build a vented-sink recombination interface (finite, positive rate).
    pub fn recombination(rate: f64) -> Result<Self, Error> {
        if !rate.is_finite() || rate <= 0.0 {
            return Err(Error::BadData(
                "recombination interface rate must be finite and > 0",
            ));
        }
        Ok(Interface::Recombination { rate })
    }

    /// Desorption rate `K_r` for a recombination interface; `None` for the
    /// linear laws.
    pub fn recombination_rate(&self) -> Option<f64> {
        match self {
            Interface::Recombination { rate } => Some(*rate),
            _ => None,
        }
    }

    /// Support check: every law passes once its data validates (the rate
    /// check above for recombination).
    fn validate(&self) -> Result<(), Error> {
        match self {
            Interface::Sieverts | Interface::Henry => Ok(()),
            Interface::Recombination { rate } => Self::recombination(*rate).map(|_| ()),
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
/// every layer and every interface rate (see [`Interface`]).
#[derive(Debug, Clone, PartialEq)]
pub struct LayerStack {
    /// Layers from the left (`x = 0`) face to the right (`x = L`) face.
    pub layers: Vec<LayerSpec>,
    /// Interface condition between `layers[k]` and `layers[k + 1]`
    /// (length `layers.len() − 1`; Sieverts, Henry, or vented-sink
    /// recombination — see [`Interface`]).
    pub interfaces: Vec<Interface>,
    /// Cached first-cell index of each layer (length `layers.len()`).
    /// Derived from `layers` at construction so per-cell queries stay
    /// O(log layers) instead of re-summing offsets on every call.
    starts: Vec<usize>,
}

impl LayerStack {
    /// Validate a stack: at least one layer, one interface per gap, every
    /// interface rate valid, and a total cell count that fits in `usize`.
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
    /// Vented-sink face concentration per stack gap \[mol/m³\]: `Some(x)`
    /// at recombination gaps, `None` at linear (Sieverts/Henry) gaps.
    pub interface_faces: Vec<Option<f64>>,
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
    /// Vented-sink face concentration per output time per stack gap
    /// (`[time][gap]`, `Some(x)` at recombination gaps, `None` at linear
    /// gaps) for closing the discrete balance with the desorption term.
    pub interface_faces: Vec<Vec<Option<f64>>>,
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
        self.rate_system_with(self.rhs.clone())
    }

    /// Borrow the tridiagonal rate system with a replacement right-hand
    /// side (for the vented-sink sensitivity columns, which share the CUT
    /// matrix with a unit face source).
    fn rate_system_with(&self, rhs: Vec<f64>) -> DiffusionSystem {
        DiffusionSystem {
            sub: self.sub.clone(),
            diag: self.diag.clone(),
            sup: self.sup.clone(),
            rhs,
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
///
/// Faces flagged in `cut` (vented-sink recombination gaps) are CUT instead:
/// the cross-coupling drops out (block-diagonal) and each side keeps only
/// its half-cell conductance `g = 2 D_cell/dx` on the diagonal — the face
/// source re-enters through the sensitivity columns, never through `rhs`.
fn assemble_layers_inner(
    stack: &LayerStack,
    left: &Boundary,
    right: &Boundary,
    cut: &[bool],
) -> Result<LayeredSystem, Error> {
    let n = stack.total_cells();
    debug_assert_eq!(cut.len(), n + 1);
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
        // West coupling: CUT face (vented-sink Dirichlet cut), interior
        // face (Sieverts form), or landed boundary.
        if i > 0 && cut[i] {
            diag[i] -= 2.0 * d_cell[i] / (dx[i] * dx[i]);
        } else if i > 0 {
            let rw = face_r[i];
            diag[i] -= 1.0 / (rw * dx[i] * sol[i]);
            sub[i - 1] = 1.0 / (rw * dx[i] * sol[i - 1]);
        } else if let Some(cf) = face_concentration(left) {
            diag[i] -= 2.0 * d_cell[i] / (dx[i] * dx[i]);
            rhs[i] += 2.0 * d_cell[i] * cf / (dx[i] * dx[i]);
        }
        // East coupling.
        if i + 1 < n && cut[i + 1] {
            diag[i] -= 2.0 * d_cell[i] / (dx[i] * dx[i]);
        } else if i + 1 < n {
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

/// Assemble the layered cell-centred finite-volume diffusion operator (no
/// CUT faces; see [`assemble_layers_inner`]).
fn assemble_layers(
    stack: &LayerStack,
    left: &Boundary,
    right: &Boundary,
) -> Result<LayeredSystem, Error> {
    assemble_layers_inner(stack, left, right, &vec![false; stack.total_cells() + 1])
}

/// Layered steady bundle: Langmuir trapped loads on each layer's isotherm
/// plus inventories for a converged mobile profile.
fn finish_layered_steady(
    stack: &LayerStack,
    sys: &LayeredSystem,
    mobile: Vec<f64>,
    flux_left: f64,
    flux_right: f64,
    interface_faces: Vec<Option<f64>>,
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
        interface_faces,
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
    let gaps = stack.interfaces.len();
    finish_layered_steady(stack, &sys, mobile, flux_left, flux_right, vec![None; gaps])
}

/// Trap-free-style steady state of a [`LayerStack`] (G7/G10).
///
/// Solves the layered linear system through [`nucleide_linalg::tridiag`] and
/// evaluates each layer's Langmuir isotherm pointwise. A one-layer stack
/// dispatches to [`solve::steady_state`] and reproduces it exactly (G7c);
/// recombination outer ends close through the layered face-response
/// construction ([`steady_layers_recombination`], the landed G5 machinery on
/// the extended interface system); stacks with vented-sink recombination
/// internal interfaces close through the CUT-matrix construction
/// ([`steady_layers_vented`], G10).
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
            interface_faces: vec![],
        });
    }
    if has_vent_gaps(stack) {
        return steady_layers_vented(stack, left, right);
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
    let gaps = stack.interfaces.len();
    finish_layered_steady(stack, &sys, mobile, flux_left, flux_right, vec![None; gaps])
}

// ---------------------------------------------------------------------------
// Vented-sink recombination internal interfaces (G10/G11)
// ---------------------------------------------------------------------------
//
// Law (R-S): one interface unknown `x ≥ 0` (single face concentration). The
// bulk matrix is CUT at the gap (block-diagonal; each block sees a Dirichlet
// cut face, so a single gap's cross-sensitivities are exactly zero).
// Adjacent-cell values are affine (`cL = aL + bL·x`, `cR = aR + bR·x` from
// one base plus one column solve per block); with half-cell conductances
// `g = 2D/dx` per side (rightward positive) the residual
// `R(x) = J_L − J_R − K_r·x² = P − Q·x − K_r·x² = 0`
// (`P = g_L·a_L + g_R·a_R`, `Q = g_L·(1−b_L) + g_R·(1−b_R) > 0`) closes in
// closed form `x = 2P/(√(Q²+4K_rP)+Q)`. Stacks with several vent faces (more
// gaps, or recombination outer ends alongside) couple through shared blocks
// and close by Newton iteration with the analytic Jacobian instead — the
// same affine construction, the G5/G6 pair-Newton shape. The θ-stepper
// coupling mirrors the G6 outer-end machinery with the faces moved interior:
// per `dt` span the cut-matrix sensitivity columns, per step the explicit
// vector carrying the old-face source plus the FULL time-independent
// outer-face source, the close fused into the trap Picard loop in the same
// slot as G6. `P ≤ 0` (reverse or zero drive) has no nonnegative root and
// fails loudly; `K_r → 0` is a continuous-concentration joint (not
// Sieverts — no solubility ratio enters the gap).

/// Whether any stack gap carries the vented-sink recombination law.
fn has_vent_gaps(stack: &LayerStack) -> bool {
    stack
        .interfaces
        .iter()
        .any(|i| i.recombination_rate().is_some())
}

/// One vented-sink internal gap: the two adjacent cells, their half-cell
/// conductances, and the desorption rate.
struct VentGap {
    gap: usize,
    left_cell: usize,
    right_cell: usize,
    g_left: f64,
    g_right: f64,
    kr: f64,
}

/// Collect the vented-sink gaps in gap order (non-empty on the vented path).
fn collect_vent_gaps(stack: &LayerStack, dx: &[f64], d: &[f64]) -> Vec<VentGap> {
    let mut gaps = Vec::new();
    for (k, interface) in stack.interfaces.iter().enumerate() {
        if let Interface::Recombination { rate } = interface {
            let right_cell = stack.layer_start(k + 1);
            let left_cell = right_cell - 1;
            gaps.push(VentGap {
                gap: k,
                left_cell,
                right_cell,
                g_left: 2.0 * d[left_cell] / dx[left_cell],
                g_right: 2.0 * d[right_cell] / dx[right_cell],
                kr: *rate,
            });
        }
    }
    gaps
}

/// One unknown of the vented face system: a recombination outer end or an
/// internal vented-sink gap. Order in `faces`: outer-left (if any), then the
/// gaps in gap order, then outer-right (if any).
enum VentFace {
    OuterLeft {
        kr: f64,
        g: f64,
    },
    Gap {
        gap: usize,
        kr: f64,
        g_left: f64,
        g_right: f64,
        left_cell: usize,
        right_cell: usize,
    },
    OuterRight {
        kr: f64,
        g: f64,
    },
}

impl VentFace {
    fn kr(&self) -> f64 {
        match self {
            VentFace::OuterLeft { kr, .. }
            | VentFace::Gap { kr, .. }
            | VentFace::OuterRight { kr, .. } => *kr,
        }
    }

    /// Operator source entries `(cell, g/dx)` coupling this face value into
    /// the CUT system (the face source `Q·y`, scaled by `dt·θ` per step).
    fn sources(&self, dx: &[f64], n: usize) -> Vec<(usize, f64)> {
        match self {
            VentFace::OuterLeft { g, .. } => vec![(0, g / dx[0])],
            VentFace::Gap {
                g_left,
                g_right,
                left_cell,
                right_cell,
                ..
            } => vec![
                (*left_cell, g_left / dx[*left_cell]),
                (*right_cell, g_right / dx[*right_cell]),
            ],
            VentFace::OuterRight { g, .. } => vec![(n - 1, g / dx[n - 1])],
        }
    }
}

/// Adjacent-cell response of one vent face: each coupled side answers
/// `c(y) = a + Σ_l b[l]·y_l` (outer faces couple one side, gaps two).
struct VentResponse {
    a: Vec<f64>,
    b: Vec<Vec<f64>>,
}

/// CUT assembly shared by the vented steady state and transient: the faces,
/// the CUT operator (recombination outer ends as Dirichlet(0) faces,
/// vented gaps CUT), and the face positions of the outer ends.
struct VentAssembly {
    faces: Vec<VentFace>,
    sys: LayeredSystem,
    outer_left: Option<usize>,
    outer_right: Option<usize>,
}

fn vent_assemble(
    stack: &LayerStack,
    left: &Boundary,
    right: &Boundary,
) -> Result<VentAssembly, Error> {
    let n = stack.total_cells();
    let dx = stack.cell_widths();
    let d = stack.diffusivities()?;
    let gaps = collect_vent_gaps(stack, &dx, &d);
    debug_assert!(!gaps.is_empty());
    let mut faces = Vec::new();
    let mut outer_left = None;
    let mut outer_right = None;
    if let Some(kr) = left.recombination_rate() {
        outer_left = Some(faces.len());
        faces.push(VentFace::OuterLeft {
            kr,
            g: 2.0 * d[0] / dx[0],
        });
    }
    for gap in &gaps {
        faces.push(VentFace::Gap {
            gap: gap.gap,
            kr: gap.kr,
            g_left: gap.g_left,
            g_right: gap.g_right,
            left_cell: gap.left_cell,
            right_cell: gap.right_cell,
        });
    }
    if let Some(kr) = right.recombination_rate() {
        outer_right = Some(faces.len());
        faces.push(VentFace::OuterRight {
            kr,
            g: 2.0 * d[n - 1] / dx[n - 1],
        });
    }
    let bl0 = left
        .recombination_rate()
        .map_or_else(|| left.clone(), |_| Boundary::Dirichlet(0.0));
    let br0 = right
        .recombination_rate()
        .map_or_else(|| right.clone(), |_| Boundary::Dirichlet(0.0));
    let mut cut = vec![false; n + 1];
    for gap in &gaps {
        cut[stack.layer_start(gap.gap + 1)] = true;
    }
    let sys = assemble_layers_inner(stack, &bl0, &br0, &cut)?;
    Ok(VentAssembly {
        faces,
        sys,
        outer_left,
        outer_right,
    })
}

/// Adjacent-cell response rows for every vent face from a base profile and
/// the sensitivity columns (`c = base + Σ col·y`).
fn vent_responses(base: &[f64], cols: &[Vec<f64>], faces: &[VentFace]) -> Vec<VentResponse> {
    faces
        .iter()
        .map(|face| match face {
            VentFace::OuterLeft { .. } => VentResponse {
                a: vec![base[0]],
                b: vec![cols.iter().map(|col| col[0]).collect()],
            },
            VentFace::Gap {
                left_cell,
                right_cell,
                ..
            } => VentResponse {
                a: vec![base[*left_cell], base[*right_cell]],
                b: vec![
                    cols.iter().map(|col| col[*left_cell]).collect(),
                    cols.iter().map(|col| col[*right_cell]).collect(),
                ],
            },
            VentFace::OuterRight { .. } => {
                let last = base.len() - 1;
                VentResponse {
                    a: vec![base[last]],
                    b: vec![cols.iter().map(|col| col[last]).collect()],
                }
            }
        })
        .collect()
}

/// Residual of vent face `j`: outer ends `K_r·y² − g·(c − y)`, gaps
/// `g_L·(c_L − y) − g_R·(y − c_R) − K_r·y²`.
fn vent_residual(face: &VentFace, resp: &VentResponse, y: &[f64], j: usize) -> f64 {
    let c =
        |side: usize| resp.a[side] + resp.b[side].iter().zip(y).map(|(b, v)| b * v).sum::<f64>();
    match face {
        VentFace::OuterLeft { kr, g } | VentFace::OuterRight { kr, g } => {
            kr * y[j] * y[j] - g * (c(0) - y[j])
        }
        VentFace::Gap {
            kr,
            g_left,
            g_right,
            ..
        } => g_left * (c(0) - y[j]) - g_right * (y[j] - c(1)) - kr * y[j] * y[j],
    }
}

/// Analytic Jacobian entry `∂f_j/∂y_l` of the vent face system.
fn vent_jacobian(face: &VentFace, resp: &VentResponse, y: &[f64], j: usize, l: usize) -> f64 {
    let kron = f64::from(j == l);
    match face {
        VentFace::OuterLeft { kr, g } | VentFace::OuterRight { kr, g } => {
            -g * (resp.b[0][l] - kron) + 2.0 * kr * y[j] * kron
        }
        VentFace::Gap {
            kr,
            g_left,
            g_right,
            ..
        } => {
            g_left * resp.b[0][l] + g_right * resp.b[1][l]
                - (g_left + g_right) * kron
                - 2.0 * kr * y[j] * kron
        }
    }
}

/// Closed form for a lone vent face (`R(y) = P − Q·y − K_r·y² = 0`,
/// rationalized positive root). `P ≤ 0` has no nonnegative root and fails
/// loudly — reverse or zero drive is not supported on vented gaps (steady
/// construction; per-step transient closes use the clamped form below).
fn vent_closed_single(face: &VentFace, resp: &VentResponse) -> Result<f64, Error> {
    let (p, q) = vent_pq(face, resp, false);
    if p <= 0.0 {
        return Err(Error::BadData(
            "vented-sink interface has non-positive permeation drive (no nonnegative face root without forward drive)",
        ));
    }
    Ok(vent_closed_pq(face.kr(), p, q))
}

/// Closed form for a lone vent face on a transient step: the G6 clamping
/// stance — negative adjacent responses (Crank–Nicolson ringing or
/// intermediate Picard iterates, never a converged physical state) clamp
/// to zero instead of failing the step, and zero drive closes at `x = 0`.
fn vent_closed_single_step(face: &VentFace, resp: &VentResponse) -> f64 {
    let (p, q) = vent_pq(face, resp, true);
    vent_closed_pq(face.kr(), p, q)
}

/// Drive and slope of the lone-face residual `R(y) = P − Q·y − K_r·y²`,
/// optionally clamping negative adjacents to zero first.
fn vent_pq(face: &VentFace, resp: &VentResponse, clamp: bool) -> (f64, f64) {
    let adj = |a: f64| if clamp { a.max(0.0) } else { a };
    match face {
        VentFace::OuterLeft { g, .. } | VentFace::OuterRight { g, .. } => {
            (g * adj(resp.a[0]), g * (1.0 - resp.b[0][0]))
        }
        VentFace::Gap {
            g_left, g_right, ..
        } => (
            g_left * adj(resp.a[0]) + g_right * adj(resp.a[1]),
            g_left * (1.0 - resp.b[0][0]) + g_right * (1.0 - resp.b[1][0]),
        ),
    }
}

/// Rationalized positive root `2P/(√(Q²+4K_rP)+Q)` (`P = 0` closes at 0).
fn vent_closed_pq(kr: f64, p: f64, q: f64) -> f64 {
    let p = p.max(0.0);
    let disc = q.mul_add(q, 4.0 * kr * p);
    (2.0 * p / (disc.sqrt() + q)).max(0.0)
}

/// Newton iteration on the coupled vent face system with the analytic
/// Jacobian from `start` (the G5 pair-Newton shape generalized), to `rtol` /
/// `atol` on the flux residuals scaled by `|K_r·y²|`. Fails with
/// [`Error::NotConverged`] when the cap is exhausted or the Jacobian goes
/// singular/non-finite (hard fail, never a partial face).
fn vent_newton(
    faces: &[VentFace],
    resps: &[VentResponse],
    start: Vec<f64>,
    rtol: f64,
    atol: f64,
) -> Result<Vec<f64>, Error> {
    let m = faces.len();
    let mut y: Vec<f64> = start.into_iter().map(|v| v.max(0.0)).collect();
    for _ in 0..G5_NEWTON_MAX {
        let f: Vec<f64> = faces
            .iter()
            .zip(resps)
            .enumerate()
            .map(|(j, (face, resp))| vent_residual(face, resp, &y, j))
            .collect();
        if faces
            .iter()
            .zip(&y)
            .zip(&f)
            .all(|((face, &yj), &fj)| fj.abs() <= atol + rtol * (face.kr() * yj * yj).abs())
        {
            return Ok(y);
        }
        let jac: Vec<Vec<f64>> = faces
            .iter()
            .zip(resps)
            .enumerate()
            .map(|(j, (face, resp))| {
                (0..m)
                    .map(|l| vent_jacobian(face, resp, &y, j, l))
                    .collect()
            })
            .collect();
        let Some(step) = solve_dense(&jac, &f) else {
            break;
        };
        for (yj, &s) in y.iter_mut().zip(&step) {
            *yj = (*yj - s).max(0.0);
        }
    }
    Err(Error::NotConverged)
}

/// Small dense solve by Gaussian elimination with partial pivoting (`None`
/// on a singular or non-finite pivot — the vent Newton treats it as a hard
/// face failure).
fn solve_dense(mat: &[Vec<f64>], rhs: &[f64]) -> Option<Vec<f64>> {
    let m = rhs.len();
    let mut a: Vec<Vec<f64>> = mat.to_vec();
    let mut b: Vec<f64> = rhs.to_vec();
    for col in 0..m {
        let mut piv = col;
        for row in col + 1..m {
            if a[row][col].abs() > a[piv][col].abs() {
                piv = row;
            }
        }
        if !a[piv][col].is_finite() || a[piv][col] == 0.0 {
            return None;
        }
        a.swap(col, piv);
        b.swap(col, piv);
        for row in col + 1..m {
            let f = a[row][col] / a[col][col];
            if !f.is_finite() {
                return None;
            }
            // Clone the short pivot segment once (m is tiny) so the update
            // reads through an iterator instead of a range loop.
            let pivot: Vec<f64> = a[col][col..m].to_vec();
            for (ark, ack) in a[row][col..m].iter_mut().zip(&pivot) {
                *ark -= f * ack;
            }
            b[row] -= f * b[col];
        }
    }
    let mut x = vec![0.0; m];
    for row in (0..m).rev() {
        let mut s = b[row];
        for k in row + 1..m {
            s -= a[row][k] * x[k];
        }
        if !a[row][row].is_finite() || a[row][row] == 0.0 {
            return None;
        }
        x[row] = s / a[row][row];
        if !x[row].is_finite() {
            return None;
        }
    }
    Some(x)
}

/// Close the vent faces for the steady state (G5 tolerances from zero; a
/// lone gap closes in closed form with the loud non-positive-drive error).
fn close_vent_faces_steady(
    base: &[f64],
    cols: &[Vec<f64>],
    faces: &[VentFace],
) -> Result<Vec<f64>, Error> {
    let resps = vent_responses(base, cols, faces);
    if faces.len() == 1 {
        return Ok(vec![vent_closed_single(&faces[0], &resps[0])?]);
    }
    vent_newton(faces, &resps, vec![0.0; faces.len()], G5_RTOL, G5_ATOL)
}

/// Close the vent faces for one transient step (caller tolerances,
/// warm-started; a lone gap closes in the clamped form, coupled faces by
/// Newton — a non-converging face system is [`Error::NotConverged`]).
fn close_vent_faces_step(
    base: &[f64],
    cols: &[Vec<f64>],
    faces: &[VentFace],
    start: Vec<f64>,
    rtol: f64,
    atol: f64,
) -> Result<Vec<f64>, Error> {
    let resps = vent_responses(base, cols, faces);
    if faces.len() == 1 {
        return Ok(vec![vent_closed_single_step(&faces[0], &resps[0])]);
    }
    vent_newton(faces, &resps, start, rtol, atol)
}

/// Vent-face profile update: `base + Σ col·y` (no extra Thomas solve).
fn apply_vent_faces(mut base: Vec<f64>, cols: &[Vec<f64>], y: &[f64]) -> Vec<f64> {
    for (col, &yj) in cols.iter().zip(y) {
        for (row, &s) in base.iter_mut().zip(col.iter()) {
            *row += s * yj;
        }
    }
    base
}

/// Steady state of a [`LayerStack`] with vented-sink recombination gaps
/// (G10).
///
/// The CUT operator is linear for frozen faces, so one basis profile plus
/// one unit-source column per vent face pins the adjacent response exactly;
/// a lone gap closes in closed form, coupled faces by the analytic-Jacobian
/// Newton (G5 tolerances). Traps evaluate pointwise on the Langmuir
/// isotherm afterwards, exactly like the linear steady path.
fn steady_layers_vented(
    stack: &LayerStack,
    left: &Boundary,
    right: &Boundary,
) -> Result<LayeredSteadyState, Error> {
    let n = stack.total_cells();
    let asm = vent_assemble(stack, left, right)?;
    let m = asm.faces.len();
    let dx = asm.sys.dx.clone();
    let base = linear_steady(&asm.sys.rate_system())?;
    let mut cols = Vec::with_capacity(m);
    for face in &asm.faces {
        let mut rhs = vec![0.0; n];
        for (cell, s) in face.sources(&dx, n) {
            rhs[cell] += s;
        }
        cols.push(linear_steady(&asm.sys.rate_system_with(rhs))?);
    }
    let y = close_vent_faces_steady(&base, &cols, &asm.faces)?;
    let mobile = apply_vent_faces(base, &cols, &y);
    let kr_left = left.recombination_rate();
    let kr_right = right.recombination_rate();
    let flux_left = kr_left.map_or_else(
        || outward_flux(left, mobile[0], asm.sys.d_cell[0], asm.sys.dx[0]),
        |kr| kr * y[asm.outer_left.expect("rate present for a closed left face")].powi(2),
    );
    let flux_right = kr_right.map_or_else(
        || {
            outward_flux(
                right,
                mobile[n - 1],
                asm.sys.d_cell[n - 1],
                asm.sys.dx[n - 1],
            )
        },
        |kr| {
            kr * y[asm
                .outer_right
                .expect("rate present for a closed right face")]
            .powi(2)
        },
    );
    let mut interface_faces = vec![None; stack.interfaces.len()];
    for (j, face) in asm.faces.iter().enumerate() {
        if let VentFace::Gap { gap, .. } = face {
            interface_faces[*gap] = Some(y[j]);
        }
    }
    finish_layered_steady(
        stack,
        &asm.sys,
        mobile,
        flux_left,
        flux_right,
        interface_faces,
    )
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
/// matrix, fused into the same Picard loop; stacks with vented-sink
/// recombination gaps close per step through the CUT-matrix construction
/// ([`solve_layers_vented`], G11). A one-layer stack dispatches to
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
        let rows = sol.times.len();
        return Ok(LayeredSolution {
            times: sol.times,
            mobile: sol.mobile,
            trapped: sol.trapped,
            flux_left: sol.flux_left,
            flux_right: sol.flux_right,
            dx: vec![sol.dx; n],
            initial: sol.initial,
            interface_faces: vec![vec![]; rows],
        });
    }
    if has_vent_gaps(stack) {
        return solve_layers_vented(stack, left, right, grid, initial, opts);
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

    let gaps = stack.interfaces.len();
    Ok(LayeredSolution {
        times: grid.times.clone(),
        mobile: mobile_out,
        trapped: trapped_out,
        flux_left,
        flux_right,
        dx: sys.dx.clone(),
        initial: initial.clone(),
        interface_faces: vec![vec![None; gaps]; grid.times.len()],
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
        interface_faces: vec![vec![None; stack.interfaces.len()]; grid.times.len()],
    })
}

/// Layered transient with vented-sink recombination gaps (G11).
///
/// Mirrors [`solve_layers_recombination`] with the faces moved interior: the
/// CUT operator carries the half-cell conductance of every vent face (the
/// face source re-enters through the sensitivity columns and the explicit
/// old-face term), one sensitivity column `S = M⁻¹·dt·θ·q` per vent face is
/// built per `dt` span, and the faces close per step in closed form (a lone
/// gap) or by the analytic-Jacobian Newton (coupled faces, warm-started from
/// the previous step), fused into the trap Picard loop with the shared
/// `rtol`/`atol`. Outward fluxes on recombination ends report the converged
/// `K_r·y²`; the discrete balance closes against the reported boundary
/// fluxes *plus* the desorption sinks `Σ K_r·x²` from the recorded gap
/// faces. Any face or Picard failure is [`Error::NotConverged`]; a
/// non-positive step drive is the loud vent-drive error.
#[allow(clippy::too_many_lines)]
fn solve_layers_vented(
    stack: &LayerStack,
    left: &Boundary,
    right: &Boundary,
    grid: &TimeGrid,
    initial: &InitialState,
    opts: &SolverOptions,
) -> Result<LayeredSolution, Error> {
    let n = stack.total_cells();
    let asm = vent_assemble(stack, left, right)?;
    let m = asm.faces.len();
    let theta = opts.theta.value();
    let dx = asm.sys.dx.clone();
    let sources: Vec<Vec<(usize, f64)>> =
        asm.faces.iter().map(|face| face.sources(&dx, n)).collect();

    let rates: Vec<Vec<(f64, f64)>> = (0..n)
        .map(|i| stack.trap_rates_at(i))
        .collect::<std::result::Result<Vec<_>, Error>>()?;
    let any_traps = rates.iter().any(|r| !r.is_empty());

    let mut c = initial.mobile.clone();
    let mut ct = initial.trapped.clone();
    let mut y = vec![0.0; m];

    let mut mobile_out = Vec::with_capacity(grid.times.len());
    let mut trapped_out = Vec::with_capacity(grid.times.len());
    let mut flux_left = Vec::with_capacity(grid.times.len());
    let mut flux_right = Vec::with_capacity(grid.times.len());
    let mut faces_out: Vec<Vec<Option<f64>>> = Vec::with_capacity(grid.times.len());

    let mut t_prev = 0.0_f64;
    let mut steps = 0_usize;

    for &t_out in &grid.times {
        let span = t_out - t_prev;
        let n_sub = ((span / opts.dt_max).ceil() as usize).max(1);
        let dt = span / n_sub as f64;
        if dt < opts.dt_min {
            return Err(Error::BadOption("output spacing needs a step below dt_min"));
        }
        let m_sub: Vec<f64> = asm.sys.sub.iter().map(|v| -dt * theta * v).collect();
        let m_diag: Vec<f64> = asm.sys.diag.iter().map(|v| 1.0 - dt * theta * v).collect();
        let m_sup: Vec<f64> = asm.sys.sup.iter().map(|v| -dt * theta * v).collect();
        let mut sens = Vec::with_capacity(m);
        for src in &sources {
            let mut v = vec![0.0; n];
            for (cell, s) in src {
                v[*cell] += dt * theta * s;
            }
            sens.push(
                nucleide_linalg::tridiag::solve(&m_sub, &m_diag, &m_sup, &v)
                    .map_err(crate::solve::tridiag_err)?,
            );
        }
        for _ in 0..n_sub {
            let mut e = vec![0.0; n];
            for i in 0..n {
                let mut a_c = asm.sys.diag[i] * c[i];
                if i > 0 {
                    a_c += asm.sys.sub[i - 1] * c[i - 1];
                }
                if i + 1 < n {
                    a_c += asm.sys.sup[i] * c[i + 1];
                }
                e[i] = c[i] + dt * (1.0 - theta) * a_c + dt * asm.sys.rhs[i];
            }
            for (src, &yj) in sources.iter().zip(&y) {
                for (cell, s) in src {
                    e[*cell] += dt * (1.0 - theta) * s * yj;
                }
            }
            if !any_traps {
                let base = nucleide_linalg::tridiag::solve(&m_sub, &m_diag, &m_sup, &e)
                    .map_err(crate::solve::tridiag_err)?;
                y = close_vent_faces_step(&base, &sens, &asm.faces, y, opts.rtol, opts.atol)?;
                c = apply_vent_faces(base, &sens, &y);
            } else {
                let mut c_iter = c.clone();
                let mut ct_iter = ct.clone();
                let mut y_iter = y.clone();
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
                    let y_next = close_vent_faces_step(
                        &base,
                        &sens,
                        &asm.faces,
                        y_iter.clone(),
                        opts.rtol,
                        opts.atol,
                    )?;
                    let c_next = apply_vent_faces(base, &sens, &y_next);
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
                    for jj in 0..asm.faces.len() {
                        let (nf, prev) = (y_next[jj], y_iter[jj]);
                        let scale = opts.atol + opts.rtol * nf.abs().max(prev.abs());
                        err = err.max((nf - prev).abs() / scale);
                    }
                    c_iter = c_next;
                    ct_iter = ct_next;
                    y_iter = y_next;
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
                y = y_iter;
            }
            steps += 1;
            if steps > opts.max_steps {
                return Err(Error::StepBudget(opts.max_steps));
            }
        }
        t_prev = t_out;
        mobile_out.push(c.clone());
        trapped_out.push(ct.clone());
        flux_left.push(asm.outer_left.map_or_else(
            || outward_flux(left, c[0], asm.sys.d_cell[0], asm.sys.dx[0]),
            |idx| {
                left.recombination_rate()
                    .expect("rate present for a closed left face")
                    * y[idx].powi(2)
            },
        ));
        flux_right.push(asm.outer_right.map_or_else(
            || outward_flux(right, c[n - 1], asm.sys.d_cell[n - 1], asm.sys.dx[n - 1]),
            |idx| {
                right
                    .recombination_rate()
                    .expect("rate present for a closed right face")
                    * y[idx].powi(2)
            },
        ));
        let mut row = vec![None; stack.interfaces.len()];
        for (j, face) in asm.faces.iter().enumerate() {
            if let VentFace::Gap { gap, .. } = face {
                row[*gap] = Some(y[j]);
            }
        }
        faces_out.push(row);
    }

    Ok(LayeredSolution {
        times: grid.times.clone(),
        mobile: mobile_out,
        trapped: trapped_out,
        flux_left,
        flux_right,
        dx: asm.sys.dx.clone(),
        initial: initial.clone(),
        interface_faces: faces_out,
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
        // The vented-sink recombination law constructs fine (G10); its rate
        // validates like a recombination end (finite, positive).
        let layer = |d: f64, s: f64| {
            LayerSpec::new(5e-4, 8, d, 0.0, s, vec![], vec![500.0], vec![]).unwrap()
        };
        assert!(LayerStack::new(
            vec![layer(1e-9, 1.0), layer(1e-9, 1.0)],
            vec![Interface::Henry]
        )
        .is_ok());
        assert!(LayerStack::new(
            vec![layer(1e-9, 1.0), layer(1e-9, 1.0)],
            vec![Interface::recombination(1e-7).unwrap()]
        )
        .is_ok());
        assert!(Interface::recombination(0.0).is_err());
        assert!(Interface::recombination(f64::INFINITY).is_err());
        assert!(Interface::recombination(f64::NAN).is_err());
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

    // -----------------------------------------------------------------------
    // Gates: vented-sink recombination internal interfaces (G10/G11)
    // -----------------------------------------------------------------------
    //
    // Provenance: synthetic stacks with hand-derived closed forms. A
    // trap-free Dirichlet-bounded block is nodally exact for linear profiles
    // (the FV face fluxes reproduce `D·Δc/L` exactly), so the discrete vent
    // residual is the continuum quadratic `R(x) = P − Q·x − K_r·x²` to
    // roundoff and the face value below is a genuine hand oracle. No
    // evaluated-library data, consistent with the analytic-gate stance.

    /// Vented gate stack R (2-layer): L₁ = L₂ = 5e-4 m, D₁ = 1e-9,
    /// D₂ = 5e-10 m²/s, K₁ = 2.0, K₂ = 0.5 (the solubilities play no role
    /// at the vented gap — the face is a concentration, not a potential),
    /// 64 + 64 cells, 500 K, one vented-sink gap.
    fn gate_stack_rec(kr: f64) -> LayerStack {
        LayerStack::new(
            vec![
                LayerSpec::new(5e-4, 64, 1e-9, 0.0, 2.0, vec![], vec![500.0], vec![]).unwrap(),
                LayerSpec::new(5e-4, 64, 5e-10, 0.0, 0.5, vec![], vec![500.0], vec![]).unwrap(),
            ],
            vec![Interface::recombination(kr).unwrap()],
        )
        .unwrap()
    }

    /// Continuum vented-sink root for a 2-layer Dirichlet(c0)|Dirichlet(cN)
    /// stack: `R(x) = D₁(c₀−x)/L₁ − D₂(x−cN)/L₂ − K_r·x² = 0`.
    fn vent_closed_2layer(d1: f64, l1: f64, d2: f64, l2: f64, c0: f64, cn: f64, kr: f64) -> f64 {
        let p = d1 * c0 / l1 + d2 * cn / l2;
        let q = d1 / l1 + d2 / l2;
        2.0 * p / ((q * q + 4.0 * kr * p).sqrt() + q)
    }

    /// Gap fluxes `(J_L, J_R)` (rightward positive) of a converged vented
    /// steady profile from the reported face value.
    fn vent_gap_fluxes(
        stack: &LayerStack,
        s: &LayeredSteadyState,
        gap: usize,
        x: f64,
    ) -> (f64, f64) {
        let ir = stack.layer_start(gap + 1);
        let il = ir - 1;
        let d = stack.diffusivities().unwrap();
        let dx = stack.cell_widths();
        let g_l = 2.0 * d[il] / dx[il];
        let g_r = 2.0 * d[ir] / dx[ir];
        (g_l * (s.mobile[il] - x), g_r * (x - s.mobile[ir]))
    }

    #[test]
    fn g10a_two_layer_vented_sink_closed_form() {
        // G10a: 2-layer stack with one vented gap, Dirichlet(1.0)|
        // Dirichlet(0.0), K_r = 3e-7. The discrete profile is the exact
        // block-linear continuum profile at the nodes, so the face value,
        // fluxes, and profile pin at G5-class tolerances (1e-12 relative
        // with a 1e-18 absolute floor).
        let kr = 3e-7;
        let stack = gate_stack_rec(kr);
        let s = steady_layers(
            &stack,
            &Boundary::dirichlet(1.0).unwrap(),
            &Boundary::dirichlet(0.0).unwrap(),
        )
        .unwrap();
        let x_star = vent_closed_2layer(1e-9, 5e-4, 5e-10, 5e-4, 1.0, 0.0, kr);
        let x = s.interface_faces[0].expect("vented gap reports its face");
        assert!((x - x_star).abs() <= 1e-12 * x_star, "{x} vs {x_star}");
        let j = 5e-10 * x_star / 5e-4;
        assert!(
            (s.flux_right - j).abs() <= 1e-12 * j + 1e-18,
            "{} vs {j}",
            s.flux_right
        );
        // Outer balance carries the desorption sink: F_L + F_R + K_r·x² = 0.
        assert!((s.flux_left + s.flux_right + kr * x * x).abs() <= 1e-12 * j + 1e-18);
        // Gap residual to the pinned tolerance.
        let (jl, jr) = vent_gap_fluxes(&stack, &s, 0, x);
        assert!(
            (jl - jr - kr * x * x).abs() <= 1e-12 * j + 1e-18,
            "{jl} vs {jr}"
        );
        // Block-linear profile at the cell centres.
        for (i, &c) in s.mobile.iter().enumerate() {
            let xc = s.centres[i];
            let want = if xc <= 5e-4 {
                1.0 + (x_star - 1.0) * xc / 5e-4
            } else {
                x_star * (1.0 - (xc - 5e-4) / 5e-4)
            };
            assert!((c - want).abs() <= 1e-12, "cell {i}: {c} vs {want}");
        }
        assert!(s.mobile.iter().all(|&c| c >= 0.0));
    }

    #[test]
    fn g10b_vent_rate_limits() {
        // K_r → ∞ pins the left block to its Dirichlet(0)-terminated
        // profile (flux D₁·c₀/L₁ = 2e-6, pinned to 4 digits).
        let s = steady_layers(
            &gate_stack_rec(100.0),
            &Boundary::dirichlet(1.0).unwrap(),
            &Boundary::dirichlet(0.0).unwrap(),
        )
        .unwrap();
        assert!((s.flux_left + 2e-6).abs() <= 1e-3 * 2e-6);
        // K_r → 0 is the continuous-concentration joint: identical layers
        // recover the uncut landed slab (not a Sieverts jump — there is no
        // K ratio at a vented gap).
        let tiny = LayerStack::new(
            vec![
                LayerSpec::new(5e-4, 32, 1e-9, 0.0, 1.0, vec![], vec![500.0], vec![]).unwrap(),
                LayerSpec::new(5e-4, 32, 1e-9, 0.0, 1.0, vec![], vec![500.0], vec![]).unwrap(),
            ],
            vec![Interface::recombination(1e-24).unwrap()],
        )
        .unwrap();
        let s = steady_layers(
            &tiny,
            &Boundary::dirichlet(1.0).unwrap(),
            &Boundary::dirichlet(0.0).unwrap(),
        )
        .unwrap();
        let slab = no_traps(1e-3, 64, 1e-9);
        let landed = solve::steady_state(
            &slab,
            &Boundary::dirichlet(1.0).unwrap(),
            &Boundary::dirichlet(0.0).unwrap(),
        )
        .unwrap();
        assert!(max_diff(&s.mobile, &landed.mobile) <= 1e-12);
        // Same limit on contrasting solubilities: c stays joint-continuous
        // across the gap (the Sieverts stack jumps by ~0.67 there).
        let s = steady_layers(
            &gate_stack_rec(1e-24),
            &Boundary::dirichlet(1.0).unwrap(),
            &Boundary::dirichlet(0.0).unwrap(),
        )
        .unwrap();
        assert!((s.mobile[63] - s.mobile[64]).abs() <= 1e-2);
    }

    /// Mixed-law gate stack (4-layer): L = 3e-4 m each,
    /// D = 1e-9 / 4e-10 / 2.5e-10 / 1e-9 m²/s,
    /// K = 1.0 / 0.8 / 0.6 / 1.2, 32 cells per layer, gaps
    /// [Sieverts, Henry, vented-sink]. The CUT partitions the stack into a
    /// linear 3-layer super-block (resistance R_L) and a 1-layer block
    /// (R_R); with the face in concentration units the residual is the
    /// scalar quadratic `R(x) = u₀/R_L − x/(K₂R_L) − x/(K₃R_R) − K_r·x²`.
    fn gate_stack_mixed_vent(kr: f64) -> LayerStack {
        LayerStack::new(
            vec![
                LayerSpec::new(3e-4, 32, 1e-9, 0.0, 1.0, vec![], vec![500.0], vec![]).unwrap(),
                LayerSpec::new(3e-4, 32, 4e-10, 0.0, 0.8, vec![], vec![500.0], vec![]).unwrap(),
                LayerSpec::new(3e-4, 32, 2.5e-10, 0.0, 0.6, vec![], vec![500.0], vec![]).unwrap(),
                LayerSpec::new(3e-4, 32, 1e-9, 0.0, 1.2, vec![], vec![500.0], vec![]).unwrap(),
            ],
            vec![
                Interface::Sieverts,
                Interface::Henry,
                Interface::recombination(kr).unwrap(),
            ],
        )
        .unwrap()
    }

    #[test]
    fn g10c_mixed_law_stack_with_vented_gap() {
        // G10c: one interface of each law on the same stack. The linear gaps
        // hold flux continuity to roundoff; the vented gap closes the scalar
        // quadratic against the super-block resistances at 1e-12.
        let kr = 5e-7;
        let stack = gate_stack_mixed_vent(kr);
        let s = steady_layers(
            &stack,
            &Boundary::dirichlet(1.0).unwrap(),
            &Boundary::dirichlet(0.0).unwrap(),
        )
        .unwrap();
        let r_l = 3e-4 / (1e-9 * 1.0) + 3e-4 / (4e-10 * 0.8) + 3e-4 / (2.5e-10 * 0.6);
        let r_r = 3e-4 / (1e-9 * 1.2);
        let p = 1.0 / r_l;
        let q = 1.0 / (0.6 * r_l) + 1.0 / (1.2 * r_r);
        let x_star = 2.0 * p / ((q * q + 4.0 * kr * p).sqrt() + q);
        let x = s.interface_faces[2].expect("vented gap reports its face");
        assert!(s.interface_faces[0].is_none() && s.interface_faces[1].is_none());
        assert!((x - x_star).abs() <= 1e-12 * x_star, "{x} vs {x_star}");
        let jl = (1.0 - x_star / 0.6) / r_l;
        let jr = (x_star / 1.2) / r_r;
        assert!(
            (s.flux_right - jr).abs() <= 1e-12 * jr + 1e-18,
            "{} vs {jr}",
            s.flux_right
        );
        assert!((s.flux_left + jl).abs() <= 1e-12 * jl + 1e-18);
        assert!((s.flux_left + s.flux_right + kr * x * x).abs() <= 1e-12 * jl + 1e-18);
        let (gl, gr) = vent_gap_fluxes(&stack, &s, 2, x);
        assert!(
            (gl - gr - kr * x * x).abs() <= 1e-12 * jl + 1e-18,
            "{gl} vs {gr}"
        );
        // Flux continuity across the Sieverts gap *and* the Henry gap: every
        // linear interior face carries its super-block flux to roundoff (the
        // vented-gap face itself carries the jump, so it is skipped).
        let rec_face = stack.layer_start(3);
        for (k, jf) in interior_face_fluxes(&stack, &s).iter().enumerate() {
            if k + 1 == rec_face {
                continue;
            }
            let want = if k + 1 < rec_face { jl } else { jr };
            assert!(
                (jf - want).abs() <= 1e-12 * want.abs() + 1e-18,
                "face {}: {jf} vs {want}",
                k + 1
            );
        }
        // Potential continuity at both linear gaps from each side.
        let (u0l, u0r) = interface_u_sides(&stack, &s, 0, jl);
        assert!((u0l - u0r).abs() <= 1e-12);
        let (u1l, u1r) = interface_u_sides(&stack, &s, 1, jl);
        assert!((u1l - u1r).abs() <= 1e-12);
        assert!(s.mobile.iter().all(|&c| c >= 0.0));
    }

    #[test]
    fn g10d_two_vented_gaps_coupled_newton() {
        // G10d: two vented gaps share the middle block, so the faces couple
        // and close by Newton (no scalar oracle — residual, outer balance
        // with both desorption sinks, and positivity pin the solve).
        let stack = LayerStack::new(
            vec![
                LayerSpec::new(5e-4, 32, 1e-9, 0.0, 2.0, vec![], vec![500.0], vec![]).unwrap(),
                LayerSpec::new(5e-4, 32, 4e-10, 0.0, 0.8, vec![], vec![500.0], vec![]).unwrap(),
                LayerSpec::new(5e-4, 32, 2.5e-10, 0.0, 0.6, vec![], vec![500.0], vec![]).unwrap(),
            ],
            vec![
                Interface::recombination(3e-7).unwrap(),
                Interface::recombination(5e-7).unwrap(),
            ],
        )
        .unwrap();
        let s = steady_layers(
            &stack,
            &Boundary::dirichlet(1.0).unwrap(),
            &Boundary::dirichlet(0.0).unwrap(),
        )
        .unwrap();
        let x0 = s.interface_faces[0].expect("gap 0 reports its face");
        let x1 = s.interface_faces[1].expect("gap 1 reports its face");
        assert!(x0 > 0.0 && x1 > 0.0);
        // Faces order downhill with the drive.
        assert!(x0 > x1);
        // Residuals pin the Newton contract (G5_ATOL + G5_RTOL·K_r·x²,
        // plus reassociation noise from this independent recompute); the
        // outer balance then closes modulo the two residuals.
        let scale = s.flux_right.abs().max(1e-300);
        let mut res = 0.0_f64;
        for (gap, &x, kr) in [(0_usize, &x0, 3e-7), (1_usize, &x1, 5e-7)] {
            let (jl, jr) = vent_gap_fluxes(&stack, &s, gap, x);
            let r = jl - jr - kr * x * x;
            assert!(
                r.abs() <= G5_ATOL + G5_RTOL * kr * x * x + 1e-16,
                "gap {gap}: {jl} vs {jr}"
            );
            res += r.abs();
        }
        let desorb = 3e-7 * x0 * x0 + 5e-7 * x1 * x1;
        assert!((s.flux_left + s.flux_right + desorb).abs() <= res + 1e-12 * scale + 1e-18);
        assert!(s.mobile.iter().all(|&c| c >= 0.0));
    }

    #[test]
    fn g10e_vented_gap_with_recombination_outer_end() {
        // G10e: a vented gap plus a recombination outer end couple through
        // the end block and close by Newton — the internal residual, the
        // outer balance with the internal sink, and positivity pin it.
        let stack = gate_stack_rec(3e-7);
        let kr_out = 1e-6;
        let s = steady_layers(
            &stack,
            &Boundary::dirichlet(1.0).unwrap(),
            &Boundary::recombination(kr_out).unwrap(),
        )
        .unwrap();
        let x = s.interface_faces[0].expect("vented gap reports its face");
        assert!(x > 0.0);
        let scale = s.flux_right.abs().max(1e-300);
        let (jl, jr) = vent_gap_fluxes(&stack, &s, 0, x);
        let r = jl - jr - 3e-7 * x * x;
        assert!(
            r.abs() <= G5_ATOL + G5_RTOL * 3e-7 * x * x + 1e-16,
            "{jl} vs {jr}"
        );
        // The outer recombination face carries its own Newton contract
        // (G5_ATOL + G5_RTOL·K_r·cf²), which joins the balance bound.
        assert!(
            (s.flux_left + s.flux_right + 3e-7 * x * x).abs()
                <= r.abs() + G5_ATOL + G5_RTOL * s.flux_right.abs() + 1e-12 * scale + 1e-18
        );
        assert!(s.mobile.iter().all(|&c| c >= 0.0));
    }

    #[test]
    fn g10f_vented_steady_with_traps_on_isotherm() {
        // Traps evaluate on the owning layer's Langmuir isotherm behind the
        // vented close, exactly like the linear steady path.
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
            vec![Interface::recombination(3e-7).unwrap()],
        )
        .unwrap();
        let s = steady_layers(
            &stack,
            &Boundary::dirichlet(1.0).unwrap(),
            &Boundary::dirichlet(0.0).unwrap(),
        )
        .unwrap();
        for (i, (&c, row)) in s.mobile.iter().zip(&s.trapped).enumerate() {
            let k = stack.layer_of(i).expect("cell owned");
            for (j, &ct) in row.iter().enumerate() {
                let trap = &stack.layers[k].traps[j];
                let (kk, pp) = trap.rates(500.0).unwrap();
                let want = solve::equilibrium_trapped(trap.site_density, kk / pp, c).unwrap();
                assert!((ct - want).abs() <= 1e-12 * want.max(1e-300));
            }
        }
        assert_eq!(s.trapped[0].len(), 0);
        assert_eq!(s.trapped[32].len(), 1);
        let x = s.interface_faces[0].expect("vented gap reports its face");
        let (jl, jr) = vent_gap_fluxes(&stack, &s, 0, x);
        let scale = s.flux_right.abs().max(1e-300);
        assert!((jl - jr - 3e-7 * x * x).abs() <= 1e-12 * scale + 1e-18);
    }

    #[test]
    fn g10g_nonpositive_drive_is_a_loud_error() {
        // Zero drive (both ends at zero, no source) gives P = 0: no
        // nonnegative vent root exists, so the steady state fails loudly
        // instead of returning a silent zero. The transient from the zero
        // state instead follows the G6 clamping stance per step and returns
        // the zero trajectory.
        let stack = gate_stack_rec(3e-7);
        let zero = Boundary::dirichlet(0.0).unwrap();
        let err = steady_layers(&stack, &zero, &zero).unwrap_err();
        assert!(matches!(err, Error::BadData(_)), "{err}");
        let grid = TimeGrid::new(vec![10.0]).unwrap();
        let opts = SolverOptions::default();
        let sol = solve_layers(&stack, &zero, &zero, &grid, &stack.zero_state(), &opts).unwrap();
        assert!(sol.mobile[0].iter().all(|&c| c == 0.0));
        assert_eq!(sol.flux_left, vec![0.0]);
        assert_eq!(sol.flux_right, vec![0.0]);
        assert_eq!(sol.interface_faces[0], vec![Some(0.0)]);
    }

    /// Vented gate stack, coarse (32 + 32 cells) for the transient gates.
    fn gate_stack_rec_coarse(kr: f64) -> LayerStack {
        LayerStack::new(
            vec![
                LayerSpec::new(5e-4, 32, 1e-9, 0.0, 2.0, vec![], vec![500.0], vec![]).unwrap(),
                LayerSpec::new(5e-4, 32, 5e-10, 0.0, 0.5, vec![], vec![500.0], vec![]).unwrap(),
            ],
            vec![Interface::recombination(kr).unwrap()],
        )
        .unwrap()
    }

    #[test]
    fn g11a_theta_balance_with_desorption() {
        // G11a: trap-free vented stack, one θ-step per output interval — the
        // discrete balance Iᵏ − Iᵏ⁻¹ = dt·[θRᵏ + (1−θ)Rᵏ⁻¹] with
        // R = −(F_L + F_R) − K_r·x² closes against the reported boundary
        // fluxes plus the desorption sink to roundoff.
        let kr = 3e-7;
        let stack = gate_stack_rec_coarse(kr);
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
        let x: Vec<f64> = sol
            .interface_faces
            .iter()
            .map(|row| row[0].unwrap())
            .collect();
        assert!(x.iter().all(|&v| v >= 0.0));
        let r_now = |k: usize| -(sol.flux_left[k] + sol.flux_right[k]) - kr * x[k] * x[k];
        // t = 0 from the zero state: F_L = −2D₁c₀/dx₀, F_R = 0, x = 0.
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
        for row in &sol.mobile {
            assert!(row.iter().all(|&c| c >= 0.0));
        }
    }

    /// Vented mobile profile at `t_end` with uniform steps of `dt_max`.
    fn vent_profile(theta: crate::solve::Theta, dt_max: f64, t_end: f64) -> Vec<f64> {
        let stack = gate_stack_rec(3e-7);
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

    #[test]
    fn g11b_dt_halving_order_preserved() {
        // G11b: dt-halving on the vented stack shows the θ-method order in
        // the resolved band at t = 200 s. Crank–Nicolson pins the spike's
        // (8, 4, 2) pair (measured 1.89): finer pairs compare a
        // ringing-dominated solution against a ringing-free one (the G6d
        // phenomenon — coarse-dt late-time errors are ringing, not
        // truncation), so only the spike pair is gated. Backward Euler is
        // L-stable and pins both pairs.
        let t_end = 200.0;
        let cn: Vec<Vec<f64>> = [8.0_f64, 4.0, 2.0]
            .iter()
            .map(|&dt| vent_profile(crate::solve::Theta::CrankNicolson, dt, t_end))
            .collect();
        let order = (max_diff(&cn[0], &cn[1]) / max_diff(&cn[1], &cn[2])).log2();
        assert!((1.5..=2.5).contains(&order), "Crank–Nicolson order {order}");
        let be: Vec<Vec<f64>> = [8.0_f64, 4.0, 2.0, 1.0]
            .iter()
            .map(|&dt| vent_profile(crate::solve::Theta::BackwardEuler, dt, t_end))
            .collect();
        for w in be.windows(3) {
            let order = (max_diff(&w[0], &w[1]) / max_diff(&w[1], &w[2])).log2();
            assert!((0.7..=1.3).contains(&order), "backward Euler order {order}");
        }
    }

    #[test]
    fn g11c_late_time_asymptote_to_vented_steady() {
        // G11c: the vented transient from a clean slab lands on the vented
        // steady state at t = 20000 s (the interface mode decays slower than
        // the linear stack lag, so the run is ~240 linear lags).
        let kr = 3e-7;
        let stack = gate_stack_rec(kr);
        let left = Boundary::dirichlet(1.0).unwrap();
        let right = Boundary::dirichlet(0.0).unwrap();
        let seg_a = solve_layers(
            &stack,
            &left,
            &right,
            &TimeGrid::new(vec![300.0]).unwrap(),
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
            &right,
            &TimeGrid::new(vec![19700.0]).unwrap(),
            &init,
            &SolverOptions {
                dt_max: 5.0,
                rtol: 1e-10,
                atol: 1e-14,
                ..Default::default()
            },
        )
        .unwrap();
        let steady = steady_layers(&stack, &left, &right).unwrap();
        assert!((seg_b.flux_right[0] - steady.flux_right).abs() <= 1e-6 * steady.flux_right.abs());
        let rate = seg_b.interface_faces[0][0].unwrap().powi(2) * kr;
        assert!(
            (seg_b.flux_left[0] + seg_b.flux_right[0] + rate).abs()
                <= 1e-6 * steady.flux_right.abs()
        );
        assert!(max_diff(&seg_b.mobile[0], &steady.mobile) <= 1e-9);
        assert!(seg_b.mobile[0].iter().all(|&c| c >= 0.0));
    }

    #[test]
    fn g11d_trap_coupled_mass_balance_with_desorption() {
        // G11d (build-cycle gate — the unprototyped Picard fusion): the
        // trap-coupled vented transient closes the discrete balance over
        // mobile + trapped inventory against the boundary fluxes plus the
        // desorption sink. Traps sit in the driven layer so they engage on
        // the gate interval, and the run is backward Euler: its explicit
        // vector stays nonneg structurally, keeping the early-front regime
        // out of the face clamp (Crank–Nicolson explicit wiggles near a
        // steep front would trip it at O(1e-13)).
        let kr = 3e-7;
        let stack = LayerStack::new(
            vec![
                LayerSpec::new(
                    5e-4,
                    16,
                    1e-9,
                    0.0,
                    2.0,
                    vec![TrapSpec::new(0.05, 0.0, 0.01, 0.0, 2.0).unwrap()],
                    vec![500.0],
                    vec![],
                )
                .unwrap(),
                LayerSpec::new(5e-4, 16, 5e-10, 0.0, 0.5, vec![], vec![500.0], vec![]).unwrap(),
            ],
            vec![Interface::recombination(kr).unwrap()],
        )
        .unwrap();
        let dx = stack.cell_widths();
        let left = Boundary::dirichlet(1.0).unwrap();
        let right = Boundary::dirichlet(0.0).unwrap();
        let opts = SolverOptions {
            theta: crate::solve::Theta::BackwardEuler,
            dt_max: 2.0,
            rtol: 1e-10,
            atol: 1e-14,
            ..Default::default()
        };
        let sol = solve_layers(
            &stack,
            &left,
            &right,
            &TimeGrid::new(vec![2.0, 4.0]).unwrap(),
            &stack.zero_state(),
            &opts,
        )
        .unwrap();
        let total: Vec<f64> = sol
            .mobile
            .iter()
            .zip(&sol.trapped)
            .map(|(mrow, trows)| {
                mrow.iter().zip(&dx).map(|(c, w)| c * w).sum::<f64>()
                    + trows
                        .iter()
                        .zip(&dx)
                        .map(|(row, w)| row.iter().sum::<f64>() * w)
                        .sum::<f64>()
            })
            .collect();
        let x: Vec<f64> = sol
            .interface_faces
            .iter()
            .map(|row| row[0].unwrap())
            .collect();
        let r_now = |k: usize| -(sol.flux_left[k] + sol.flux_right[k]) - kr * x[k] * x[k];
        // Backward Euler: Iᵏ − Iᵏ⁻¹ = dt·Rᵏ.
        let mut i_prev = 0.0_f64;
        for (k, dt) in [2.0, 2.0].iter().enumerate() {
            let want = i_prev + dt * r_now(k);
            assert!(
                (total[k] - want).abs() < 1e-9 * want.abs().max(1e-300),
                "row {k}: {} vs {want}",
                total[k]
            );
            i_prev = total[k];
        }
        for row in &sol.mobile {
            assert!(row.iter().all(|&c| c >= 0.0));
        }
        // The driven-layer traps genuinely engaged (not a vacuous pass).
        let mut ct_max = 0.0_f64;
        for rows in &sol.trapped {
            for row in rows.iter().take(16) {
                assert!(row.iter().all(|&ct| (0.0..=2.0).contains(&ct)));
                ct_max = ct_max.max(row[0]);
            }
        }
        assert!(ct_max > 1e-6);
    }
}
