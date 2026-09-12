//! The [`Material`] composition model: construction, conversions, mixing,
//! arithmetic, and (de)serialization.

use std::collections::BTreeMap;
use std::ops::{Add, Div, Mul, Sub};

use nucleide_nuclei::NuclideId;
use serde::de::{Deserialize, Deserializer};
use serde::ser::{Serialize, SerializeStruct, Serializer};

use crate::Error;

/// True for values that cannot serve as a positive total mass.
fn not_positive(v: f64) -> bool {
    v.is_nan() || v <= 0.0
}

/// True for values that cannot serve as a mixing fraction.
fn is_negative(v: f64) -> bool {
    v.is_nan() || v < 0.0
}

/// Source of per-nuclide atomic masses in g/mol.
///
/// Atomic-mass-dependent operations ([`Material::from_atom_frac`] and
/// [`Material::atom_fractions`]) are generic over this trait so the material
/// crate never depends on the nuclear-data tables directly. Integrating the
/// real tables later is a single `impl MassProvider for nucleide_nuclei::data::...`.
pub trait MassProvider {
    /// Atomic mass of the nuclide identified by raw `nucid`
    /// (`(Z*1000 + A)*10_000 + state`), or `None` if unknown.
    fn mass(&self, nucid: u32) -> Option<f64>;
}

/// A [`MassProvider`] that knows no masses.
///
/// Useful as an explicit placeholder; every lookup returns `None`, so
/// mass-dependent conversions fail with [`crate::Error::MissingMass`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct NoMasses;

impl MassProvider for NoMasses {
    fn mass(&self, _nucid: u32) -> Option<f64> {
        None
    }
}

/// [`MassProvider`] backed by the AME2020 tables in `nucleide_nuclei::data`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Ame2020;

impl MassProvider for Ame2020 {
    fn mass(&self, nucid: u32) -> Option<f64> {
        nucleide_nuclei::data::atomic_mass(nucid)
    }
}

/// A nuclear material: nuclide masses plus optional density and metadata.
///
/// The composition stores absolute masses per nuclide, in grams by
/// convention; only relative amounts matter for fraction-based consumers,
/// which normalize on demand. Density is deliberately separate from the
/// composition: it is a property of the physical stream and is
/// not scaled or combined by the arithmetic operators except where noted.
///
/// Combining two materials clears density and metadata (a mixture has no
/// single density); scalar scaling preserves them.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Material {
    /// Stored masses (grams) keyed by nuclide.
    pub comp: BTreeMap<NuclideId, f64>,
    density: Option<f64>,
    metadata: Option<serde_json::Value>,
}

impl Material {
    /// An empty material.
    pub fn new() -> Self {
        Self::default()
    }

    /// Build a material from atom counts/fractions, converting to masses via
    /// `m_i = n_i * M_i` with atomic masses from `masses`.
    ///
    /// Entries with zero atom count are skipped. Fails with
    /// [`crate::Error::MissingMass`] if any nonzero entry lacks a known
    /// atomic mass.
    pub fn from_atom_frac(
        atoms: &[(NuclideId, f64)],
        masses: &impl MassProvider,
        density: Option<f64>,
    ) -> crate::Result<Self> {
        let mut mat = Self {
            density,
            ..Self::default()
        };
        for &(id, atoms) in atoms {
            if atoms == 0.0 {
                continue;
            }
            let am = masses
                .mass(id.nucid())
                .ok_or(crate::Error::MissingMass(id))?;
            mat.comp.insert(id, am * atoms);
        }
        Ok(mat)
    }

    /// Add `mass` grams of `id`, accumulating when already present.
    pub fn add_nuclide(&mut self, id: NuclideId, mass: f64) {
        *self.comp.entry(id).or_insert(0.0) += mass;
    }

    /// Remove a nuclide, returning its stored mass if present.
    pub fn remove_nuclide(&mut self, id: NuclideId) -> Option<f64> {
        self.comp.remove(&id)
    }

    /// Drop the entire composition (density and metadata are kept).
    pub fn clear(&mut self) {
        self.comp.clear();
    }

    /// Total stored mass in grams.
    pub fn mass(&self) -> f64 {
        self.comp.values().sum()
    }

    /// Mass density previously set on this material, if any.
    pub fn density(&self) -> Option<f64> {
        self.density
    }

    /// Set (or unset) the mass density.
    pub fn set_density(&mut self, density: Option<f64>) {
        self.density = density;
    }

    /// Free-form metadata attached to this material.
    pub fn metadata(&self) -> Option<&serde_json::Value> {
        self.metadata.as_ref()
    }

    /// Replace the free-form metadata.
    pub fn set_metadata(&mut self, metadata: Option<serde_json::Value>) {
        self.metadata = metadata;
    }

    /// Normalized weight fractions; they sum to one.
    pub fn weight_fractions(&self) -> crate::Result<BTreeMap<NuclideId, f64>> {
        let total = self.mass();
        if not_positive(total) {
            return Err(crate::Error::Degenerate);
        }
        Ok(self.comp.iter().map(|(&id, &m)| (id, m / total)).collect())
    }

    /// Normalized atom fractions; they sum to one.
    ///
    /// Each nuclide contributes moles proportional to `mass / M`; atomic
    /// masses come from `masses`.
    pub fn atom_fractions(
        &self,
        masses: &impl MassProvider,
    ) -> crate::Result<BTreeMap<NuclideId, f64>> {
        let mut moles = BTreeMap::new();
        let mut total = 0.0;
        for (&id, &m) in &self.comp {
            let am = masses
                .mass(id.nucid())
                .ok_or(crate::Error::MissingMass(id))?;
            let n = m / am;
            moles.insert(id, n);
            total += n;
        }
        if not_positive(total) {
            return Err(crate::Error::Degenerate);
        }
        Ok(moles.into_iter().map(|(id, n)| (id, n / total)).collect())
    }

    /// Mix streams weighted by relative mass amounts.
    ///
    /// Fractions need not sum to one; they are relative weights of each
    /// stream's full mass.
    pub fn mix_by_mass(parts: &[(&Material, f64)]) -> crate::Result<Self> {
        let mut out = Self::new();
        for &(mat, frac) in parts {
            if is_negative(frac) {
                return Err(crate::Error::NegativeFraction(frac));
            }
            for (&id, &m) in &mat.comp {
                out.add_nuclide(id, frac * m);
            }
        }
        if not_positive(out.mass()) {
            return Err(crate::Error::Degenerate);
        }
        Ok(out)
    }

    /// Mix streams weighted by relative volumes, converting each stream's
    /// contribution through its own density (`m = v * rho`). Every input
    /// must have a positive density set.
    pub fn mix_by_volume(parts: &[(&Material, f64)]) -> crate::Result<Self> {
        let mut out = Self::new();
        for &(mat, vol) in parts {
            if is_negative(vol) {
                return Err(crate::Error::NegativeFraction(vol));
            }
            match mat.density() {
                Some(rho) if rho > 0.0 => {
                    for (&id, &m) in &mat.comp {
                        out.add_nuclide(id, vol * rho * m / mat.mass());
                    }
                }
                _ => return Err(crate::Error::MissingDensity),
            }
        }
        if not_positive(out.mass()) {
            return Err(crate::Error::Degenerate);
        }
        Ok(out)
    }

    /// Split this material into product and tails streams by per-nuclide
    /// separation efficiency.
    ///
    /// Each listed nuclide sends the fraction `eff` of its stored mass to
    /// the product stream and `1 - eff` to the tails stream; nuclides absent
    /// from `effs` send nothing to product (`eff = 0`). Mass is conserved
    /// per nuclide: `product + tails == self` up to floating-point rounding.
    /// Efficiencies must be finite values in `[0, 1]` (else
    /// [`crate::Error::InvalidEfficiency`]); a repeated nuclide keeps its
    /// last-listed efficiency. Both outputs clear density and metadata (a
    /// split stream has no single density), and nuclides with exactly zero
    /// mass on a side are dropped from that side.
    pub fn separate(&self, effs: &[(NuclideId, f64)]) -> crate::Result<(Self, Self)> {
        let mut table = BTreeMap::new();
        for &(id, eff) in effs {
            if !eff.is_finite() || eff < 0.0 || eff > 1.0 {
                return Err(crate::Error::InvalidEfficiency(eff));
            }
            table.insert(id, eff);
        }
        let mut product = Self::new();
        let mut tails = Self::new();
        for (&id, &m) in &self.comp {
            let eff = table.get(&id).copied().unwrap_or(0.0);
            let p = m * eff;
            let t = m - p;
            if p != 0.0 {
                product.comp.insert(id, p);
            }
            if t != 0.0 {
                tails.comp.insert(id, t);
            }
        }
        Ok((product, tails))
    }

    /// Blend streams at fixed ratios with explicit normalization.
    ///
    /// Ratios are relative target proportions: they are normalized by their
    /// sum (`w_i = r_i / Σr`) and the output is the weighted average
    /// `Σ w_i · mat_i` (density and metadata cleared, as for the arithmetic
    /// operators). Unlike the cycamore mixer this never falls back to a
    /// silent uniform split: an empty slice or an all-zero (or non-finite)
    /// ratio sum fails with [`crate::Error::Degenerate`], and any negative
    /// or non-finite ratio fails with [`crate::Error::NegativeFraction`].
    pub fn blend(parts: &[(&Material, f64)]) -> crate::Result<Self> {
        if parts.is_empty() {
            return Err(crate::Error::Degenerate);
        }
        let mut sum = 0.0;
        for &(_, ratio) in parts {
            if !ratio.is_finite() || ratio < 0.0 {
                return Err(crate::Error::NegativeFraction(ratio));
            }
            sum += ratio;
        }
        if !(sum > 0.0 && sum.is_finite()) {
            return Err(crate::Error::Degenerate);
        }
        let mut out = Self::new();
        for &(mat, ratio) in parts {
            let w = ratio / sum;
            for (&id, &m) in &mat.comp {
                out.add_nuclide(id, w * m);
            }
        }
        if not_positive(out.mass()) {
            return Err(crate::Error::Degenerate);
        }
        Ok(out)
    }

    /// Scale all stored masses by `factor`, keeping density and metadata.
    fn scaled(&self, factor: f64) -> Self {
        Self {
            comp: self.comp.iter().map(|(&id, &m)| (id, m * factor)).collect(),
            density: self.density,
            metadata: self.metadata.clone(),
        }
    }
}

impl Add for Material {
    type Output = Material;

    /// Combine two mass streams: per-nuclide masses add. Density and
    /// metadata are cleared on the mixture. Nuclides whose combined mass is
    /// exactly zero are dropped.
    fn add(self, rhs: Material) -> Material {
        let mut comp = self.comp;
        for (id, m) in rhs.comp {
            *comp.entry(id).or_insert(0.0) += m;
        }
        comp.retain(|_, m| *m != 0.0);
        Material {
            comp,
            density: None,
            metadata: None,
        }
    }
}

impl Sub for Material {
    type Output = Material;

    /// Remove a mass stream: per-nuclide masses subtract. Density and
    /// metadata are cleared. Nuclides whose combined mass is exactly zero
    /// are dropped.
    fn sub(self, rhs: Material) -> Material {
        let mut comp = self.comp;
        for (id, m) in rhs.comp {
            *comp.entry(id).or_insert(0.0) -= m;
        }
        comp.retain(|_, m| *m != 0.0);
        Material {
            comp,
            density: None,
            metadata: None,
        }
    }
}

impl Mul<f64> for Material {
    type Output = Material;

    /// Scale every stored mass by `rhs` (density unchanged).
    fn mul(self, rhs: f64) -> Material {
        self.scaled(rhs)
    }
}

impl Div<f64> for Material {
    type Output = Material;

    /// Divide every stored mass by `rhs` (density unchanged).
    ///
    /// # Panics
    /// If `rhs` is zero.
    fn div(self, rhs: f64) -> Material {
        assert!(rhs != 0.0, "cannot divide a material mass by zero");
        self.scaled(1.0 / rhs)
    }
}

impl Add<f64> for Material {
    type Output = Material;

    /// Raise the total mass by `rhs` grams while preserving relative
    /// composition. Density is unchanged.
    ///
    /// # Panics
    /// If the current mass is non-positive or the new total would be.
    fn add(self, rhs: f64) -> Material {
        let total = self.mass();
        let new_total = total + rhs;
        assert!(
            total > 0.0 && new_total > 0.0,
            "cannot add {rhs} g to a material of {total} g"
        );
        self.scaled(new_total / total)
    }
}

impl Sub<f64> for Material {
    type Output = Material;

    /// Lower the total mass by `rhs` grams while preserving relative
    /// composition.
    ///
    /// # Panics
    /// If the current mass is non-positive or the remainder would be.
    fn sub(self, rhs: f64) -> Material {
        let total = self.mass();
        let new_total = total - rhs;
        assert!(
            total > 0.0 && new_total > 0.0,
            "cannot subtract {rhs} g from a material of {total} g"
        );
        self.scaled(new_total / total)
    }
}

impl Serialize for Material {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let comp: BTreeMap<String, f64> =
            self.comp.iter().map(|(id, m)| (id.to_name(), *m)).collect();
        let mut state = serializer.serialize_struct("Material", 3)?;
        state.serialize_field("comp", &comp)?;
        state.serialize_field("density", &self.density)?;
        state.serialize_field("metadata", &self.metadata)?;
        state.end()
    }
}

impl<'de> Deserialize<'de> for Material {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(serde::Deserialize)]
        struct RawMaterial {
            comp: BTreeMap<String, f64>,
            density: Option<f64>,
            metadata: Option<serde_json::Value>,
        }

        let raw = RawMaterial::deserialize(deserializer)?;
        let mut comp = BTreeMap::new();
        for (name, mass) in raw.comp {
            let id = NuclideId::from_name(&name).map_err(|source| {
                serde::de::Error::custom(Error::BadNuclide {
                    name: name.clone(),
                    source,
                })
            })?;
            comp.insert(id, mass);
        }
        Ok(Material {
            comp,
            density: raw.density,
            metadata: raw.metadata,
        })
    }
}

// ---------------------------------------------------------------------------
// Radioanalytics
//
// Activity (and derived specific activity) computed from stored masses via a
// [`DecayProvider`] plus the existing [`MassProvider`]. Decay heat
// ([`Material::decay_heat`]) additionally needs mean recoverable decay
// energies ([`DecayEnergyProvider`], backed by the `decay_energy.tsv` table
// in `nucleide_nuclei::data`).
//
// Dose per gram ([`Material::dose_per_g`]) follows PyNE's
// `Material::dose_per_g` using the EPA/DOE/GENII ingestion/inhalation/air-soil
// factors in `dose_factors.tsv` ([`DoseProvider`], backed by
// [`DoseFactors`]/`nucleide_nuclei::data::DoseData`). Screening-level only —
// not for safety decisions (upstream HNF-5636/PyNE disclaimer).
// ---------------------------------------------------------------------------

/// Avogadro constant, atoms per mole (exact, 2019 SI).
pub const AVOGADRO: f64 = 6.022_140_76e23;

/// One unified atomic mass unit in grams (2022 CODATA).
pub const GRAMS_PER_U: f64 = 1.660_539_068_92e-24;

/// One MeV in joules (exact, 2019 SI: 1 eV = 1.602176634e-19 J).
pub const MEV_TO_JOULES: f64 = 1.602_176_634e-13;

/// Curies per becquerel (PyNE `Ci_per_Bq`, `src/data.cpp`).
pub const CI_PER_BQ: f64 = 2.702_702_7e-11;

/// Picocuries per becquerel (PyNE `pCi_per_Bq`, `Material::dose_per_g`).
pub const PCI_PER_BQ: f64 = 27.027_027;

/// Dose pathway and source types (re-exported from `nucleide-nuclei` so
/// analytics call sites need one import).
pub use nucleide_nuclei::data::{DosePathway, DoseSource};

/// Source of per-nuclide decay constants λ in inverse seconds.
///
/// Like [`MassProvider`], injected as a trait so analytics never hard-depend
/// on decay data availability.
pub trait DecayProvider {
    /// Decay constant of the nuclide identified by raw `nucid`, or `None`
    /// if unknown (stable nuclides included).
    fn decay_constant(&self, nucid: u32) -> Option<f64>;
}

/// A [`DecayProvider`] that knows no decays.
///
/// Every lookup returns `None`, so activities fail explicitly with
/// [`AnalyticsError::MissingDecay`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct NoDecay;

impl DecayProvider for NoDecay {
    fn decay_constant(&self, _nucid: u32) -> Option<f64> {
        None
    }
}

/// [`DecayProvider`] backed by the ENDF/B-VIII.0 half-life table in
/// `nucleide_nuclei::data`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ChainDecays;

impl DecayProvider for ChainDecays {
    fn decay_constant(&self, nucid: u32) -> Option<f64> {
        nucleide_nuclei::data::decay_constant(nucid)
    }
}

/// Provider bundle for radioanalytic quantities.
///
/// Both providers are needed at once — atom numbers come from masses,
/// activities from decay constants — so they travel together:
///
/// ```
/// use nucleide_material::{Analytics, Ame2020, ChainDecays};
/// # let mut mat = nucleide_material::Material::new();
/// # let co = nucleide_nuclei::NuclideId::from_name("Co60").unwrap();
/// # mat.add_nuclide(co, 1e-6);
/// let an = Analytics { masses: &Ame2020, decays: &ChainDecays };
/// let a = mat.activity(&an).unwrap();
/// ```
pub struct Analytics<'a> {
    /// Atomic masses (u) for gram → atom conversion.
    pub masses: &'a dyn MassProvider,
    /// Decay constants λ (1/s).
    pub decays: &'a dyn DecayProvider,
}

impl std::fmt::Debug for Analytics<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Analytics").finish_non_exhaustive()
    }
}

/// Errors from radioanalytics beyond the shared [`enum@crate::Error`] set.
#[derive(Debug, Error)]
pub enum AnalyticsError {
    /// No decay data was available for a requested nuclide.
    ///
    /// Reserved for provider-level absence. The built-in analytics paths
    /// below never construct this anymore: a known atomic mass with no
    /// decay constant is a stable nuclide (λ = 0) and contributes exactly
    /// 0.0 instead of erroring, so only a genuinely unknown nuclide (no
    /// mass data, [`crate::Error::MissingMass`]) can still fail. The
    /// variant is retained for API compatibility.
    #[error("no decay data available for nuclide `{0}`")]
    MissingDecay(NuclideId),
    /// No mean decay energy was available for a requested nuclide.
    #[error("no decay energy available for nuclide `{0}`")]
    MissingEnergy(NuclideId),
    /// No dose factor was available for a requested nuclide/pathway/source.
    #[error("no dose factor available for nuclide `{0}`")]
    MissingDose(NuclideId),
    /// An underlying composition failure (missing mass, degenerate total).
    #[error(transparent)]
    Core(#[from] crate::Error),
}

/// Source of per-nuclide mean recoverable decay energies in MeV per decay.
///
/// Like [`MassProvider`], injected as a trait so analytics never hard-depend
/// on decay-energy data availability. Kept separate from [`DecayProvider`]
/// (which supplies decay constants) so depletion callers can mix sources;
/// `nucleide_nuclei::data::DecayData` implements this trait, giving Stream A
/// a single provider for both without any material↔depletion coupling
/// (material never depends on depletion; nuclei never depends on material).
pub trait DecayEnergyProvider {
    /// Mean recoverable energy per decay of the nuclide identified by raw
    /// `nucid`, in MeV, or `None` if unknown (stable nuclides included).
    fn decay_energy_mev(&self, nucid: u32) -> Option<f64>;
}

/// A [`DecayEnergyProvider`] that knows no decay energies.
///
/// Every lookup returns `None`, so heat calculations fail explicitly with
/// [`AnalyticsError::MissingEnergy`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct NoDecayEnergies;

impl DecayEnergyProvider for NoDecayEnergies {
    fn decay_energy_mev(&self, _nucid: u32) -> Option<f64> {
        None
    }
}

/// [`DecayEnergyProvider`] backed by the ENDF/B-VII.1 prompt-decay-energy
/// table in `nucleide_nuclei::data` (mean-field evaluation values, generated
/// by `scripts/gen-nuclear-data.py` — see the table docs before quoting
/// heat numbers).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct DecayEnergies;

impl DecayEnergyProvider for DecayEnergies {
    fn decay_energy_mev(&self, nucid: u32) -> Option<f64> {
        nucleide_nuclei::data::decay_energy_mev(nucid)
    }
}

impl DecayEnergyProvider for nucleide_nuclei::data::DecayData {
    fn decay_energy_mev(&self, nucid: u32) -> Option<f64> {
        nucleide_nuclei::data::decay_energy_mev(nucid)
    }
}

/// Source of per-nuclide dose factors (raw table values).
///
/// Like [`MassProvider`], injected as a trait so analytics never hard-depend
/// on dose-table availability. Kept separate from [`DecayProvider`]/
/// [`DecayEnergyProvider`] so callers can mix sources; the blanket impl for
/// `nucleide_nuclei::data::DoseData` gives a single provider with no
/// material↔nuclei circularity (material depends on nuclei, never the reverse).
pub trait DoseProvider {
    /// Raw dose factor for the nuclide identified by raw `nucid`, or `None`
    /// if the nuclide has no row for this `pathway`/`source`.
    ///
    /// Implementations backed by [`crate::DoseFactors`] return the stored
    /// `-1` sentinel for GENII/DOE air (PyNE missing-air convention);
    /// [`Material::dose_per_g`] treats negative factors as missing and fails
    /// with [`AnalyticsError::MissingDose`].
    fn dose_factor(&self, nucid: u32, pathway: DosePathway, source: DoseSource) -> Option<f64>;
}

/// A [`DoseProvider`] that knows no dose factors.
///
/// Every lookup returns `None`, so dose calculations fail explicitly with
/// [`AnalyticsError::MissingDose`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct NoDoses;

impl DoseProvider for NoDoses {
    fn dose_factor(&self, _nucid: u32, _pathway: DosePathway, _source: DoseSource) -> Option<f64> {
        None
    }
}

/// [`DoseProvider`] backed by the HNF-5636/PyNE dose-factor table in
/// `nucleide_nuclei::data` (generated by `scripts/gen-nuclear-data.py` —
/// see the table docs; screening-level only, not for safety decisions).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct DoseFactors;

impl DoseProvider for DoseFactors {
    fn dose_factor(&self, nucid: u32, pathway: DosePathway, source: DoseSource) -> Option<f64> {
        nucleide_nuclei::data::dose_factor(nucid, pathway, source)
    }
}

impl DoseProvider for nucleide_nuclei::data::DoseData {
    fn dose_factor(&self, nucid: u32, pathway: DosePathway, source: DoseSource) -> Option<f64> {
        nucleide_nuclei::data::dose_factor(nucid, pathway, source)
    }
}

impl Material {
    /// Activity `A = λ·N` per nuclide, in becquerels.
    ///
    /// Atom counts follow from stored masses through `masses`
    /// (`N = m / (M · u)` with `u = 1.66053906892e-24 g`) and decay
    /// constants through `decays`. Stable-as-zero: a known atomic mass
    /// with no decay constant is a stable nuclide (λ = 0, mirroring the
    /// chain rule where `None` decay → 0.0) and contributes exactly 0.0.
    /// Fails with [`AnalyticsError::Core`] wrapping
    /// [`crate::Error::MissingMass`] when an atomic mass is unknown, so
    /// genuinely unknown nuclides never collapse to silent zeros.
    pub fn activity(
        &self,
        analytics: &Analytics<'_>,
    ) -> std::result::Result<BTreeMap<NuclideId, f64>, AnalyticsError> {
        let mut out = BTreeMap::new();
        for (&id, &grams) in &self.comp {
            let mass_u = analytics
                .masses
                .mass(id.nucid())
                .ok_or(crate::Error::MissingMass(id))?;
            let lambda = analytics.decays.decay_constant(id.nucid()).unwrap_or(0.0);
            let atoms = grams / (mass_u * GRAMS_PER_U);
            out.insert(id, lambda * atoms);
        }
        Ok(out)
    }

    /// Specific activity of the whole material, in Bq/g: total activity
    /// divided by total stored mass. Fails with
    /// [`AnalyticsError::Core`](`crate::Error::Degenerate`) for empty or
    /// non-positive materials; otherwise identical error behavior to
    /// [`Material::activity`].
    pub fn specific_activity(&self, analytics: &Analytics<'_>) -> Result<f64, AnalyticsError> {
        let total_mass = self.mass();
        if not_positive(total_mass) {
            return Err(crate::Error::Degenerate.into());
        }
        let mut total_activity = 0.0;
        for value in self.activity(analytics)?.values() {
            total_activity += value;
        }
        Ok(total_activity / total_mass)
    }

    /// Decay heat per nuclide, in watts: `P_i = A_i · E_i`.
    ///
    /// Activities come from [`Material::activity`] (masses via `analytics`,
    /// decay constants via `analytics.decays`); mean recoverable energies
    /// per decay come from `energies` in MeV, converted with
    /// [`MEV_TO_JOULES`]. Energies are screening-level placeholders (see
    /// [`DecayEnergies`]), so heat numbers are order-of-magnitude checks,
    /// not calorimetry.
    ///
    /// Fails with [`AnalyticsError::MissingEnergy`] for radioactive
    /// nuclides without a decay-energy row; stable nuclides (zero
    /// activity, hence `P = A·E = 0` regardless of `E`) skip the energy
    /// lookup and contribute exactly 0.0. Otherwise identical error
    /// behavior to [`Material::activity`].
    pub fn decay_heat(
        &self,
        analytics: &Analytics<'_>,
        energies: &impl DecayEnergyProvider,
    ) -> Result<BTreeMap<NuclideId, f64>, AnalyticsError> {
        let activities = self.activity(analytics)?;
        let mut out = BTreeMap::new();
        for (&id, &activity_bq) in &activities {
            // Exact zero by construction (λ = 0 or zero stored mass); no
            // energy row exists for stable nuclides, and none is needed.
            if activity_bq == 0.0 {
                out.insert(id, 0.0);
                continue;
            }
            let mev = energies
                .decay_energy_mev(id.nucid())
                .ok_or(AnalyticsError::MissingEnergy(id))?;
            out.insert(id, activity_bq * mev * MEV_TO_JOULES);
        }
        Ok(out)
    }

    /// Total decay heat of the whole material, in watts: the sum of
    /// [`Material::decay_heat`]. Same error behavior.
    pub fn total_decay_heat(
        &self,
        analytics: &Analytics<'_>,
        energies: &impl DecayEnergyProvider,
    ) -> Result<f64, AnalyticsError> {
        let mut total = 0.0;
        for value in self.decay_heat(analytics, energies)?.values() {
            total += value;
        }
        Ok(total)
    }

    /// Dose per gram per nuclide, mirroring PyNE `Material::dose_per_g`.
    ///
    /// For weight fraction `w_i = m_i / m_tot`:
    ///
    /// ```text
    /// dose_i = Ci_per_Bq · N_A · w_i · λ_i · DF_i / M_i   (air/soil)
    /// dose_i = pCi_per_Bq · N_A · w_i · λ_i · DF_i / M_i  (ingest/inhale)
    /// ```
    ///
    /// with [`CI_PER_BQ`] = 2.7027027e-11, [`PCI_PER_BQ`] = 27.027027,
    /// `N_A` = [`AVOGADRO`] (PyNE uses 6.0221415e23; the difference is <0.1 ppm),
    /// `λ` from `analytics.decays`, `M` (g/mol) from `analytics.masses`, and
    /// `DF` from `doses`. Units follow the table: air `mrem/h per g per m^3`,
    /// soil `mrem/h per g per m^2`, ingest/inhale `mrem per g`. The returned
    /// map holds each nuclide's per-gram contribution; sum for the total.
    ///
    /// Screening-level only — not for safety decisions. Stable nuclides
    /// (known mass, λ = 0 or absent) contribute exactly 0.0 and skip the
    /// `doses` lookup entirely, so no [`AnalyticsError::MissingDose`] is
    /// raised for them. Fails with [`AnalyticsError::MissingDose`] when a
    /// radioactive nuclide's factor is absent or negative (`-1` GENII/DOE
    /// air sentinel); otherwise identical error behavior to
    /// [`Material::activity`] (plus `Degenerate` for empty materials).
    pub fn dose_per_g(
        &self,
        analytics: &Analytics<'_>,
        doses: &impl DoseProvider,
        pathway: DosePathway,
        source: DoseSource,
    ) -> Result<BTreeMap<NuclideId, f64>, AnalyticsError> {
        let total_mass = self.mass();
        if not_positive(total_mass) {
            return Err(crate::Error::Degenerate.into());
        }
        let per_bq = match pathway {
            DosePathway::Air | DosePathway::Soil => CI_PER_BQ,
            DosePathway::Ingest | DosePathway::Inhale => PCI_PER_BQ,
        };
        let mut out = BTreeMap::new();
        for (&id, &grams) in &self.comp {
            let mass_u = analytics
                .masses
                .mass(id.nucid())
                .ok_or(crate::Error::MissingMass(id))?;
            let lambda = analytics.decays.decay_constant(id.nucid()).unwrap_or(0.0);
            if lambda <= 0.0 {
                out.insert(id, 0.0);
                continue;
            }
            let df = doses
                .dose_factor(id.nucid(), pathway, source)
                .filter(|v| v.is_finite() && *v >= 0.0)
                .ok_or(AnalyticsError::MissingDose(id))?;
            let w = grams / total_mass;
            out.insert(id, per_bq * AVOGADRO * w * lambda * df / mass_u);
        }
        Ok(out)
    }

    /// Total dose per gram of the whole material: the sum of
    /// [`Material::dose_per_g`]. Same units and error behavior.
    pub fn total_dose_per_g(
        &self,
        analytics: &Analytics<'_>,
        doses: &impl DoseProvider,
        pathway: DosePathway,
        source: DoseSource,
    ) -> Result<f64, AnalyticsError> {
        let mut total = 0.0;
        for value in self.dose_per_g(analytics, doses, pathway, source)?.values() {
            total += value;
        }
        Ok(total)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn id(name: &str) -> NuclideId {
        NuclideId::from_name(name).unwrap()
    }

    fn close(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-12, "{a} != {b}");
    }

    /// Round-number atomic masses so mixture math stays hand-checkable.
    struct Table(BTreeMap<u32, f64>);

    impl Table {
        fn new(pairs: &[(&str, f64)]) -> Self {
            Self(
                pairs
                    .iter()
                    .map(|&(name, m)| (id(name).nucid(), m))
                    .collect(),
            )
        }
    }

    impl MassProvider for Table {
        fn mass(&self, nucid: u32) -> Option<f64> {
            self.0.get(&nucid).copied()
        }
    }

    fn water_table() -> Table {
        Table::new(&[("H1", 1.0), ("O16", 16.0)])
    }

    #[test]
    fn empty_material_has_zero_mass() {
        let mat = Material::new();
        close(mat.mass(), 0.0);
        assert!(mat.comp.is_empty());
        assert_eq!(mat.density(), None);
    }

    #[test]
    fn add_nuclide_accumulates_and_remove_returns_mass() {
        let mut mat = Material::new();
        let u5 = id("U235");
        mat.add_nuclide(u5, 10.0);
        mat.add_nuclide(u5, 5.0);
        close(mat.mass(), 15.0);
        close(mat.remove_nuclide(u5).unwrap(), 15.0);
        assert_eq!(mat.remove_nuclide(u5), None);
    }

    #[test]
    fn clear_drops_composition_only() {
        let mut mat = Material::new();
        mat.add_nuclide(id("U235"), 3.0);
        mat.add_nuclide(id("U238"), 1.0);
        mat.set_density(Some(19.1));
        mat.clear();
        assert!(mat.comp.is_empty());
        assert_eq!(mat.density(), Some(19.1));
    }

    #[test]
    fn from_atom_frac_water_hand_computed() {
        let mat = Material::from_atom_frac(
            &[(id("H1"), 2.0), (id("O16"), 1.0)],
            &water_table(),
            Some(1.0),
        )
        .unwrap();

        close(mat.comp[&id("H1")], 2.0);
        close(mat.comp[&id("O16")], 16.0);
        close(mat.mass(), 18.0);

        let wf = mat.weight_fractions().unwrap();
        close(wf[&id("H1")], 1.0 / 9.0);
        close(wf[&id("O16")], 8.0 / 9.0);

        let af = mat.atom_fractions(&water_table()).unwrap();
        close(af[&id("H1")], 2.0 / 3.0);
        close(af[&id("O16")], 1.0 / 3.0);
    }

    #[test]
    fn from_atom_frac_skips_zero_counts_and_sets_density() {
        let mat =
            Material::from_atom_frac(&[(id("H1"), 0.0), (id("O16"), 1.0)], &water_table(), None)
                .unwrap();
        assert!(!mat.comp.contains_key(&id("H1")));
        assert!(mat.comp.contains_key(&id("O16")));
        assert_eq!(mat.density(), None);
    }

    #[test]
    fn from_atom_frac_without_masses_errors() {
        let err = Material::from_atom_frac(&[(id("U235"), 1.0)], &NoMasses, None).unwrap_err();
        assert!(matches!(err, Error::MissingMass(_)));
    }

    #[test]
    fn weight_fractions_normalize_to_one() {
        let mut mat = Material::new();
        mat.add_nuclide(id("U235"), 19.0);
        mat.add_nuclide(id("U238"), 1.0);
        let wf = mat.weight_fractions().unwrap();
        close(wf[&id("U235")], 0.95);
        close(wf[&id("U238")], 0.05);
        close(wf.values().sum(), 1.0);
    }

    #[test]
    fn weight_fractions_of_empty_material_error() {
        assert!(matches!(
            Material::new().weight_fractions(),
            Err(Error::Degenerate)
        ));
    }

    #[test]
    fn atom_fractions_missing_mass_errors() {
        let mut mat = Material::new();
        mat.add_nuclide(id("U235"), 1.0);
        assert!(matches!(
            mat.atom_fractions(&NoMasses),
            Err(Error::MissingMass(_))
        ));
    }

    #[test]
    fn adding_materials_mixes_by_mass() {
        let mut fuel = Material::new();
        fuel.add_nuclide(id("U235"), 3.0);
        fuel.set_density(Some(19.0));

        let mut matrix = Material::new();
        matrix.add_nuclide(id("U238"), 1.0);
        matrix.set_density(Some(10.0));

        let mixed = fuel + matrix;
        close(mixed.mass(), 4.0);
        let wf = mixed.weight_fractions().unwrap();
        close(wf[&id("U235")], 0.75);
        close(wf[&id("U238")], 0.25);
        assert_eq!(mixed.density(), None, "mixtures have no single density");
    }

    #[test]
    fn subtracting_materials_removes_stream() {
        let mut a = Material::new();
        a.add_nuclide(id("U235"), 3.0);
        a.add_nuclide(id("U238"), 1.0);
        let mut b = Material::new();
        b.add_nuclide(id("U238"), 1.0);

        let rest = a - b;
        assert_eq!(rest.comp.len(), 1);
        close(rest.comp[&id("U235")], 3.0);
    }

    #[test]
    fn scalar_mul_div_scale_masses_and_keep_density() {
        let mut mat = Material::new();
        mat.add_nuclide(id("U235"), 3.0);
        mat.add_nuclide(id("U238"), 1.0);
        mat.set_density(Some(19.1));

        let doubled = mat.clone() * 2.0;
        close(doubled.mass(), 8.0);
        close(doubled.comp[&id("U235")], 6.0);
        assert_eq!(doubled.density(), Some(19.1));

        let quartered = doubled / 4.0;
        close(quartered.mass(), 2.0);
        close(quartered.comp[&id("U238")], 0.5);
    }

    #[test]
    fn scalar_add_sub_shift_total_mass_proportionally() {
        let mut mat = Material::new();
        mat.add_nuclide(id("U235"), 2.0);
        mat.set_density(Some(19.1));

        let grown = mat.clone() + 1.0;
        close(grown.mass(), 3.0);
        close(grown.comp[&id("U235")], 3.0);

        let shrunk = grown - 1.0;
        close(shrunk.mass(), 2.0);
        close(shrunk.comp[&id("U235")], 2.0);
        assert_eq!(shrunk.density(), Some(19.1));
    }

    #[test]
    #[should_panic(expected = "divide")]
    fn divide_by_zero_panics() {
        let _ = Material::new() / 0.0;
    }

    #[test]
    #[should_panic(expected = "cannot add")]
    fn scalar_add_to_zero_mass_panics() {
        let _ = Material::new() + 5.0;
    }

    #[test]
    #[should_panic(expected = "cannot subtract")]
    fn scalar_sub_below_zero_panics() {
        let mut mat = Material::new();
        mat.add_nuclide(id("U235"), 1.0);
        let _ = mat - 2.0;
    }

    #[test]
    fn mix_by_mass_weights_full_streams() {
        let mut a = Material::new();
        a.add_nuclide(id("U235"), 1.0);
        a.add_nuclide(id("Pu239"), 1.0);
        let mut b = Material::new();
        b.add_nuclide(id("U238"), 1.0);

        let mixed = Material::mix_by_mass(&[(&a, 1.0), (&b, 2.0)]).unwrap();
        // Stream a carries 2 g (U235 + Pu239) at weight 1; stream b carries
        // 1 g of U238 at weight 2.
        close(mixed.mass(), 4.0);
        let wf = mixed.weight_fractions().unwrap();
        close(wf[&id("U235")], 0.25);
        close(wf[&id("Pu239")], 0.25);
        close(wf[&id("U238")], 0.5);
    }

    #[test]
    fn mix_by_volume_converts_through_densities() {
        let mut heavy = Material::new();
        heavy.add_nuclide(id("U238"), 1.0);
        heavy.set_density(Some(10.0));
        let mut light = Material::new();
        light.add_nuclide(id("H1"), 1.0);
        light.set_density(Some(2.0));

        // 1 volume unit at rho=10 plus 1.5 units at rho=2:
        let mixed = Material::mix_by_volume(&[(&heavy, 1.0), (&light, 1.5)]).unwrap();
        close(mixed.comp[&id("U238")], 10.0);
        close(mixed.comp[&id("H1")], 3.0);
    }

    #[test]
    fn mix_by_volume_requires_density() {
        let mut mat = Material::new();
        mat.add_nuclide(id("U235"), 1.0);
        assert!(matches!(
            Material::mix_by_volume(&[(&mat, 1.0)]),
            Err(Error::MissingDensity)
        ));
    }

    #[test]
    fn negative_mix_fraction_rejected() {
        let mut mat = Material::new();
        mat.add_nuclide(id("U235"), 1.0);
        assert!(matches!(
            Material::mix_by_mass(&[(&mat, -1.0)]),
            Err(Error::NegativeFraction(_))
        ));
    }

    /// Feed for the separation tests: U235 10 g, U238 90 g, Pu239 1 g,
    /// Pu240 2 g, Am241 3 g, Am242 2.8 g (108.8 g total).
    fn sep_feed() -> Material {
        let mut mat = Material::new();
        mat.add_nuclide(id("U235"), 10.0);
        mat.add_nuclide(id("U238"), 90.0);
        mat.add_nuclide(id("Pu239"), 1.0);
        mat.add_nuclide(id("Pu240"), 2.0);
        mat.add_nuclide(id("Am241"), 3.0);
        mat.add_nuclide(id("Am242"), 2.8);
        mat
    }

    #[test]
    fn separate_splits_by_efficiency_and_conserves_mass() {
        // Element shorthands expanded per nuclide: U at 0.7, Pu at 0.4,
        // Am241 at 0.4; unlisted Am242 goes entirely to tails.
        let feed = sep_feed();
        let effs = [
            (id("U235"), 0.7),
            (id("U238"), 0.7),
            (id("Pu239"), 0.4),
            (id("Pu240"), 0.4),
            (id("Am241"), 0.4),
        ];
        let (product, tails) = feed.separate(&effs).unwrap();

        // Hand-computed product masses (g).
        close(product.comp[&id("U235")], 7.0);
        close(product.comp[&id("U238")], 63.0);
        close(product.comp[&id("Pu239")], 0.4);
        close(product.comp[&id("Pu240")], 0.8);
        close(product.comp[&id("Am241")], 1.2);
        assert!(!product.comp.contains_key(&id("Am242")));
        close(product.mass(), 72.4);

        // Hand-computed tails masses (g).
        close(tails.comp[&id("U235")], 3.0);
        close(tails.comp[&id("U238")], 27.0);
        close(tails.comp[&id("Pu239")], 0.6);
        close(tails.comp[&id("Pu240")], 1.2);
        close(tails.comp[&id("Am241")], 1.8);
        close(tails.comp[&id("Am242")], 2.8);
        close(tails.mass(), 36.4);

        // Per-nuclide conservation: product + tails == feed.
        for (&nuc, &m) in &feed.comp {
            let p = product.comp.get(&nuc).copied().unwrap_or(0.0);
            let t = tails.comp.get(&nuc).copied().unwrap_or(0.0);
            close(p + t, m);
        }
        close(product.mass() + tails.mass(), feed.mass());
        assert_eq!(product.density(), None);
        assert_eq!(tails.density(), None);
    }

    #[test]
    fn separate_edge_efficiencies_route_wholly() {
        let feed = sep_feed();
        // eff 1 sends everything to product; eff 0 sends all to tails.
        let (all_product, no_tails) = feed
            .separate(&[
                (id("U235"), 1.0),
                (id("U238"), 1.0),
                (id("Pu239"), 1.0),
                (id("Pu240"), 1.0),
                (id("Am241"), 1.0),
                (id("Am242"), 1.0),
            ])
            .unwrap();
        close(all_product.mass(), feed.mass());
        assert!(no_tails.comp.is_empty());

        let (no_product, all_tails) = feed.separate(&[]).unwrap();
        assert!(no_product.comp.is_empty());
        close(all_tails.mass(), feed.mass());
    }

    #[test]
    fn separate_rejects_out_of_range_efficiencies() {
        let feed = sep_feed();
        for bad in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
            assert!(
                matches!(
                    feed.separate(&[(id("U235"), bad)]),
                    Err(Error::InvalidEfficiency(_))
                ),
                "efficiency {bad} must be rejected"
            );
        }
    }

    #[test]
    fn blend_normalizes_fixed_ratios() {
        let mut a = Material::new();
        a.add_nuclide(id("U235"), 1.0);
        a.add_nuclide(id("Pu239"), 1.0);
        let mut b = Material::new();
        b.add_nuclide(id("U238"), 1.0);

        // Ratios [1, 2] normalize to [1/3, 2/3]: the output is the
        // weighted average (1/3)*a + (2/3)*b, total 4/3 g.
        let out = Material::blend(&[(&a, 1.0), (&b, 2.0)]).unwrap();
        close(out.comp[&id("U235")], 1.0 / 3.0);
        close(out.comp[&id("Pu239")], 1.0 / 3.0);
        close(out.comp[&id("U238")], 2.0 / 3.0);
        close(out.mass(), 4.0 / 3.0);
        let wf = out.weight_fractions().unwrap();
        close(wf[&id("U235")], 0.25);
        close(wf[&id("Pu239")], 0.25);
        close(wf[&id("U238")], 0.5);

        // Equal ratios [2, 2] give the plain mean: total 1.5 g.
        let half = Material::blend(&[(&a, 2.0), (&b, 2.0)]).unwrap();
        close(half.comp[&id("U235")], 0.5);
        close(half.comp[&id("Pu239")], 0.5);
        close(half.comp[&id("U238")], 0.5);
        close(half.mass(), 1.5);
        assert_eq!(half.density(), None);
    }

    #[test]
    fn blend_rejects_degenerate_and_negative_recipes() {
        let mut a = Material::new();
        a.add_nuclide(id("U235"), 1.0);
        // Empty and all-zero recipes are degenerate (no silent 1/N split).
        assert!(matches!(Material::blend(&[]), Err(Error::Degenerate)));
        assert!(matches!(
            Material::blend(&[(&a, 0.0)]),
            Err(Error::Degenerate)
        ));
        // Negative, NaN, and infinite ratios are rejected outright.
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(
                matches!(
                    Material::blend(&[(&a, bad)]),
                    Err(Error::NegativeFraction(_))
                ),
                "ratio {bad} must be rejected"
            );
        }
    }

    #[test]
    fn json_round_trip_preserves_everything() {
        let mut mat = Material::new();
        mat.add_nuclide(id("U235"), 19.0);
        mat.add_nuclide(id("Am242_m1"), 1.0);
        mat.set_density(Some(19.1));
        mat.set_metadata(Some(serde_json::json!({"enrichment": 0.03})));

        let text = serde_json::to_string(&mat).unwrap();
        let parsed: Material = serde_json::from_str(&text).unwrap();
        assert_eq!(parsed, mat);
    }

    #[test]
    fn json_uses_gnds_names_as_keys() {
        let mut mat = Material::new();
        mat.add_nuclide(id("U235"), 1.0);
        let text = serde_json::to_string(&mat).unwrap();
        assert!(
            text.contains("\"comp\":{\"U235\":1.0}"),
            "unexpected serialization: {text}"
        );
    }

    #[test]
    fn json_rejects_unknown_nuclide_names() {
        let err = serde_json::from_str::<Material>(
            r#"{"comp":{"Notanuclide":1.0},"density":null,"metadata":null}"#,
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("invalid nuclide name `Notanuclide`"), "{err}");
    }
}

#[cfg(test)]
mod radio_tests {
    use super::*;
    use std::f64::consts::LN_2;

    fn nid(name: &str) -> NuclideId {
        NuclideId::from_name(name).unwrap()
    }

    #[test]
    fn activity_of_one_gram_co60_matches_hand_calculation() {
        let mut mat = Material::new();
        mat.add_nuclide(nid("Co60"), 1.0);

        let analytics = Analytics {
            masses: &Ame2020,
            decays: &ChainDecays,
        };
        let activity = mat.activity(&analytics).unwrap();
        let co60 = nid("Co60");

        // λ from the half-life table, independently of ChainDecays.
        let lambda = LN_2 / nucleide_nuclei::data::half_life(co60.nucid()).unwrap();
        assert_eq!(activity.keys().next().copied(), Some(co60));
        assert_eq!(
            ChainDecays.decay_constant(co60.nucid()),
            Some(lambda),
            "ChainDecays must be ln(2)/t_half of the tabulated half-life"
        );
        // N = m / (M · u) atoms for 1 g.
        let mass_u = nucleide_nuclei::data::atomic_mass(co60.nucid()).unwrap();
        let expected = lambda * (1.0 / (mass_u * GRAMS_PER_U));
        assert!((activity[&co60] - expected).abs() / expected < 1e-12);
    }

    #[test]
    fn specific_activity_is_activity_per_gram_in_becquerels() {
        // 5 g of Cs137: specific activity must equal total activity / 5.
        let mut mat = Material::new();
        mat.add_nuclide(nid("Cs137"), 5.0);

        let analytics = Analytics {
            masses: &Ame2020,
            decays: &ChainDecays,
        };
        let total: f64 = mat.activity(&analytics).unwrap().values().sum();
        let spec = mat.specific_activity(&analytics).unwrap();
        assert!((spec - total / 5.0).abs() < 1e-6 * spec.abs());
        // Cs137 specific activity is ~3.2 TBq/g; sanity-band the units.
        assert!(spec > 1e12 && spec < 1e14, "{spec} Bq/g");
    }

    #[test]
    fn chain_decays_lambda_is_ln2_over_tabulated_half_life() {
        let nucid = nid("Co60").nucid();
        let lambda = ChainDecays.decay_constant(nucid).unwrap();
        let t_half = nucleide_nuclei::data::half_life(nucid).unwrap();
        assert!((lambda - LN_2 / t_half).abs() < 1e-18);
        // Stable Fe56 has no tabulated decay data.
        assert_eq!(ChainDecays.decay_constant(nid("Fe56").nucid()), None);
    }

    #[test]
    fn no_decay_provider_treats_known_masses_as_stable_zero() {
        // Provider-level absence resolves through the same stable-as-zero
        // rule: a known mass with no decay constant contributes 0.0 rather
        // than raising MissingDecay.
        let mut mat = Material::new();
        mat.add_nuclide(nid("Co60"), 1.0);

        let analytics = Analytics {
            masses: &Ame2020,
            decays: &NoDecay,
        };
        let activity = mat.activity(&analytics).unwrap();
        assert_eq!(activity[&nid("Co60")], 0.0);
        assert_eq!(mat.specific_activity(&analytics).unwrap(), 0.0);
    }

    #[test]
    fn activity_needs_masses_and_nonempty_materials() {
        let mut mat = Material::new();
        mat.add_nuclide(nid("Co60"), 1.0);
        let no_masses = Analytics {
            masses: &NoMasses,
            decays: &ChainDecays,
        };
        assert!(matches!(
            mat.activity(&no_masses),
            Err(AnalyticsError::Core(crate::Error::MissingMass(_)))
        ));

        let empty = Analytics {
            masses: &Ame2020,
            decays: &ChainDecays,
        };
        assert!(matches!(
            Material::new().specific_activity(&empty),
            Err(AnalyticsError::Core(crate::Error::Degenerate))
        ));
    }

    #[test]
    fn decay_heat_of_one_gram_co60_matches_hand_calculation() {
        let mut mat = Material::new();
        mat.add_nuclide(nid("Co60"), 1.0);

        let analytics = Analytics {
            masses: &Ame2020,
            decays: &ChainDecays,
        };
        let heat = mat.decay_heat(&analytics, &DecayEnergies).unwrap();
        let co60 = nid("Co60");

        // P = A * E with E from the ENDF/B-VII.1 decay-energy table.
        let activity: f64 = mat.activity(&analytics).unwrap()[&co60];
        let e_mev = DecayEnergies.decay_energy_mev(co60.nucid()).unwrap();
        let expected = activity * e_mev * MEV_TO_JOULES;
        assert!((heat[&co60] - expected).abs() / expected < 1e-12);
        // 1 g Co60 is ~40 TBq * ~4.2e-13 J ≈ ~17 W; sanity-band the units.
        assert!(heat[&co60] > 5.0 && heat[&co60] < 50.0, "{}", heat[&co60]);

        let total = mat.total_decay_heat(&analytics, &DecayEnergies).unwrap();
        assert!((total - expected).abs() / expected < 1e-12);
    }

    #[test]
    fn nuclei_decay_data_serves_as_energy_provider() {
        // Stream A wiring: a single DecayData covers λ and MeV with no
        // material<->depletion coupling.
        let provider = nucleide_nuclei::data::DecayData;
        assert_eq!(
            DecayEnergies.decay_energy_mev(nid("Cs137").nucid()),
            provider.decay_energy_mev(nid("Cs137").nucid())
        );

        let mut mat = Material::new();
        mat.add_nuclide(nid("Cs137"), 2.0);
        let analytics = Analytics {
            masses: &Ame2020,
            decays: &ChainDecays,
        };
        let via_facade = mat.total_decay_heat(&analytics, &provider).unwrap();
        let via_struct = mat.total_decay_heat(&analytics, &DecayEnergies).unwrap();
        assert!((via_facade - via_struct).abs() < 1e-18);
    }

    #[test]
    fn decay_heat_missing_energy_errors() {
        // Cf237 has a half-life (activity computes) but no ENDF decay tape,
        // hence no decay-energy row.
        let mut mat = Material::new();
        mat.add_nuclide(nid("Cf237"), 1.0);
        let analytics = Analytics {
            masses: &Ame2020,
            decays: &ChainDecays,
        };
        match mat.decay_heat(&analytics, &DecayEnergies).unwrap_err() {
            AnalyticsError::MissingEnergy(id) => assert_eq!(id, nid("Cf237")),
            other => panic!("{other:?}"),
        }
        // The explicit no-data provider fails on the first nuclide too.
        let mut co = Material::new();
        co.add_nuclide(nid("Co60"), 1.0);
        match co.decay_heat(&analytics, &NoDecayEnergies).unwrap_err() {
            AnalyticsError::MissingEnergy(id) => assert_eq!(id, nid("Co60")),
            other => panic!("{other:?}"),
        }
    }

    /// Synthetic mass provider with small inline maps (never vendored data).
    struct MassTable(BTreeMap<u32, f64>);

    impl MassProvider for MassTable {
        fn mass(&self, nucid: u32) -> Option<f64> {
            self.0.get(&nucid).copied()
        }
    }

    /// Synthetic dose provider with small inline maps (never vendored data).
    struct DoseTable {
        factors: BTreeMap<(u32, DosePathway, DoseSource), f64>,
    }

    impl DoseTable {
        fn new(pairs: &[(&str, DosePathway, DoseSource, f64)]) -> Self {
            Self {
                factors: pairs
                    .iter()
                    .map(|&(name, p, s, v)| (nid(name).nucid(), p, s, v))
                    .map(|(n, p, s, v)| ((n, p, s), v))
                    .collect(),
            }
        }
    }

    impl DoseProvider for DoseTable {
        fn dose_factor(&self, nucid: u32, pathway: DosePathway, source: DoseSource) -> Option<f64> {
            self.factors.get(&(nucid, pathway, source)).copied()
        }
    }

    struct ConstDecays(BTreeMap<u32, f64>);

    impl DecayProvider for ConstDecays {
        fn decay_constant(&self, nucid: u32) -> Option<f64> {
            self.0.get(&nucid).copied()
        }
    }

    #[test]
    fn dose_per_g_matches_pyne_equation_with_synthetic_data() {
        use DosePathway as P;
        // Hand-checkable masses/lambdas; DFs are synthetic (not PyNE values).
        let masses = MassTable(
            [(nid("H1").nucid(), 1.0), (nid("Co60").nucid(), 60.0)]
                .into_iter()
                .collect(),
        );
        let decays = ConstDecays(
            [(nid("H1").nucid(), 0.1), (nid("Co60").nucid(), 0.2)]
                .into_iter()
                .collect(),
        );
        let analytics = Analytics {
            masses: &masses,
            decays: &decays,
        };
        let doses = DoseTable::new(&[
            ("H1", P::Ingest, DoseSource::Epa, 2.0),
            ("Co60", P::Ingest, DoseSource::Epa, 3.0),
            ("H1", P::Air, DoseSource::Epa, 4.0),
            ("Co60", P::Air, DoseSource::Epa, 5.0),
        ]);
        let mut mat = Material::new();
        mat.add_nuclide(nid("H1"), 1.0);
        mat.add_nuclide(nid("Co60"), 3.0);
        // w_H=0.25, w_Co=0.75.
        let ingest = mat
            .dose_per_g(&analytics, &doses, P::Ingest, DoseSource::Epa)
            .unwrap();
        let e_h = PCI_PER_BQ * AVOGADRO * 0.25 * 0.1 * 2.0 / 1.0;
        let e_co = PCI_PER_BQ * AVOGADRO * 0.75 * 0.2 * 3.0 / 60.0;
        assert!((ingest[&nid("H1")] - e_h).abs() / e_h < 1e-12);
        assert!((ingest[&nid("Co60")] - e_co).abs() / e_co < 1e-12);
        let total = mat
            .total_dose_per_g(&analytics, &doses, P::Ingest, DoseSource::Epa)
            .unwrap();
        assert!((total - (e_h + e_co)).abs() / total < 1e-12);
        // Air uses Ci_per_Bq instead of pCi_per_Bq.
        let air = mat
            .dose_per_g(&analytics, &doses, P::Air, DoseSource::Epa)
            .unwrap();
        let e_air = CI_PER_BQ * AVOGADRO * 0.25 * 0.1 * 4.0 / 1.0;
        assert!((air[&nid("H1")] - e_air).abs() / e_air < 1e-12);
    }

    #[test]
    fn dose_per_g_of_one_gram_co60_matches_hand_calculation() {
        use DosePathway as P;
        use DoseSource as S;
        let mut mat = Material::new();
        mat.add_nuclide(nid("Co60"), 1.0);
        let analytics = Analytics {
            masses: &Ame2020,
            decays: &ChainDecays,
        };
        let per_nuc = mat
            .dose_per_g(&analytics, &DoseFactors, P::Ingest, S::Epa)
            .unwrap();
        let co60 = nid("Co60");
        let mass_u = nucleide_nuclei::data::atomic_mass(co60.nucid()).unwrap();
        let lambda = LN_2 / nucleide_nuclei::data::half_life(co60.nucid()).unwrap();
        let df = DoseFactors
            .dose_factor(co60.nucid(), P::Ingest, S::Epa)
            .unwrap();
        // Raw table factor is the Co-60 ingest EPA spot (2.69e-05 mrem/pCi).
        assert!((df - 2.69e-05).abs() / 2.69e-05 < 1e-9);
        let expected = PCI_PER_BQ * AVOGADRO * 1.0 * lambda * df / mass_u;
        assert!((per_nuc[&co60] - expected).abs() / expected < 1e-12);
        let total = mat
            .total_dose_per_g(&analytics, &DoseFactors, P::Ingest, S::Epa)
            .unwrap();
        assert!((total - expected).abs() / expected < 1e-12);
    }

    #[test]
    fn dose_per_g_missing_dose_errors() {
        use DosePathway as P;
        use DoseSource as S;
        let analytics = Analytics {
            masses: &Ame2020,
            decays: &ChainDecays,
        };
        // Fe56 is stable with no dose row: stable-as-zero skips the dose
        // provider, so it contributes exactly 0.0 instead of erroring.
        let mut fe = Material::new();
        fe.add_nuclide(nid("Fe56"), 1.0);
        let fe_dose = fe
            .dose_per_g(&analytics, &DoseFactors, P::Ingest, S::Epa)
            .unwrap();
        assert_eq!(fe_dose[&nid("Fe56")], 0.0);
        assert_eq!(
            fe.total_dose_per_g(&analytics, &DoseFactors, P::Ingest, S::Epa)
                .unwrap(),
            0.0
        );
        // GENII air is a -1 sentinel, treated as missing.
        let mut h3 = Material::new();
        h3.add_nuclide(nid("H3"), 1.0);
        match h3
            .dose_per_g(&analytics, &DoseFactors, P::Air, S::Genii)
            .unwrap_err()
        {
            AnalyticsError::MissingDose(id) => assert_eq!(id, nid("H3")),
            other => panic!("{other:?}"),
        }
        // Explicit no-data provider fails too.
        let mut co = Material::new();
        co.add_nuclide(nid("Co60"), 1.0);
        match co
            .dose_per_g(&analytics, &NoDoses, P::Ingest, S::Epa)
            .unwrap_err()
        {
            AnalyticsError::MissingDose(id) => assert_eq!(id, nid("Co60")),
            other => panic!("{other:?}"),
        }
        // Empty materials are degenerate.
        match Material::new()
            .total_dose_per_g(&analytics, &DoseFactors, P::Ingest, S::Epa)
            .unwrap_err()
        {
            AnalyticsError::Core(crate::Error::Degenerate) => {}
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn stable_water_contributes_exact_zeros() {
        use DosePathway as P;
        use DoseSource as S;
        // H2O: every member has an AME2020 mass but no decay data, so all
        // three observables resolve to exactly 0.0 on every pathway/source.
        let mut mat = Material::new();
        mat.add_nuclide(nid("H1"), 2.0);
        mat.add_nuclide(nid("O16"), 16.0);
        let analytics = Analytics {
            masses: &Ame2020,
            decays: &ChainDecays,
        };

        let activity = mat.activity(&analytics).unwrap();
        assert_eq!(activity[&nid("H1")], 0.0);
        assert_eq!(activity[&nid("O16")], 0.0);
        assert_eq!(mat.specific_activity(&analytics).unwrap(), 0.0);

        let heat = mat.decay_heat(&analytics, &DecayEnergies).unwrap();
        assert_eq!(heat[&nid("H1")], 0.0);
        assert_eq!(heat[&nid("O16")], 0.0);
        assert_eq!(
            mat.total_decay_heat(&analytics, &DecayEnergies).unwrap(),
            0.0
        );

        for pathway in [P::Air, P::Soil, P::Ingest, P::Inhale] {
            for source in [S::Epa, S::Doe, S::Genii] {
                let dose = mat
                    .dose_per_g(&analytics, &DoseFactors, pathway, source)
                    .unwrap();
                assert_eq!(dose[&nid("H1")], 0.0, "{pathway:?}/{source:?}");
                assert_eq!(dose[&nid("O16")], 0.0, "{pathway:?}/{source:?}");
                assert_eq!(
                    mat.total_dose_per_g(&analytics, &DoseFactors, pathway, source)
                        .unwrap(),
                    0.0,
                    "{pathway:?}/{source:?}"
                );
            }
        }
    }

    #[test]
    fn mixed_stable_plus_radioactive_matches_radioactive_only() {
        use DosePathway as P;
        use DoseSource as S;
        // Activity and heat are absolute per nuclide: the U235 entries are
        // identical with or without the stable diluent, which adds 0.0.
        let mut mixed = Material::new();
        mixed.add_nuclide(nid("U235"), 1.0);
        mixed.add_nuclide(nid("H1"), 1.0);
        let mut pure = Material::new();
        pure.add_nuclide(nid("U235"), 1.0);
        let analytics = Analytics {
            masses: &Ame2020,
            decays: &ChainDecays,
        };

        let mixed_act = mixed.activity(&analytics).unwrap();
        let pure_act = pure.activity(&analytics).unwrap();
        assert_eq!(mixed_act[&nid("U235")], pure_act[&nid("U235")]);
        assert!(pure_act[&nid("U235")] > 0.0);
        assert_eq!(mixed_act[&nid("H1")], 0.0);

        let mixed_heat = mixed.decay_heat(&analytics, &DecayEnergies).unwrap();
        let pure_heat = pure.decay_heat(&analytics, &DecayEnergies).unwrap();
        assert_eq!(mixed_heat[&nid("U235")], pure_heat[&nid("U235")]);
        assert!(pure_heat[&nid("U235")] > 0.0);
        assert_eq!(mixed_heat[&nid("H1")], 0.0);

        // Dose is per gram, so the U235 entry scales by its weight fraction
        // (1/2 here) while H1 contributes exactly 0.0.
        let mixed_dose = mixed
            .dose_per_g(&analytics, &DoseFactors, P::Ingest, S::Epa)
            .unwrap();
        let pure_dose = pure
            .dose_per_g(&analytics, &DoseFactors, P::Ingest, S::Epa)
            .unwrap();
        assert_eq!(mixed_dose[&nid("H1")], 0.0);
        assert_eq!(mixed_dose[&nid("U235")], pure_dose[&nid("U235")] * 0.5);
    }

    #[test]
    fn unknown_nuclide_without_mass_still_errors() {
        use DosePathway as P;
        use DoseSource as S;
        // Og296 parses (Z = 118) but has no AME2020 row anywhere: the
        // guard rail against silent zeros for genuinely unknown nuclides.
        let og = nid("Og296");
        assert_eq!(
            nucleide_nuclei::data::atomic_mass(og.nucid()),
            None,
            "Og296 must stay absent from the mass table"
        );
        let mut mat = Material::new();
        mat.add_nuclide(og, 1.0);
        let analytics = Analytics {
            masses: &Ame2020,
            decays: &ChainDecays,
        };

        match mat.activity(&analytics).unwrap_err() {
            AnalyticsError::Core(crate::Error::MissingMass(id)) => assert_eq!(id, og),
            other => panic!("{other:?}"),
        }
        match mat.decay_heat(&analytics, &DecayEnergies).unwrap_err() {
            AnalyticsError::Core(crate::Error::MissingMass(id)) => assert_eq!(id, og),
            other => panic!("{other:?}"),
        }
        match mat
            .dose_per_g(&analytics, &DoseFactors, P::Ingest, S::Epa)
            .unwrap_err()
        {
            AnalyticsError::Core(crate::Error::MissingMass(id)) => assert_eq!(id, og),
            other => panic!("{other:?}"),
        }
        let msg = crate::Error::MissingMass(og).to_string();
        assert!(msg.contains("Og296"), "{msg}");
    }

    #[test]
    fn nuclei_dose_data_serves_as_dose_provider() {
        use DosePathway as P;
        use DoseSource as S;
        let provider = nucleide_nuclei::data::DoseData;
        assert_eq!(
            DoseFactors.dose_factor(nid("Cs137").nucid(), P::Inhale, S::Epa),
            provider.dose_factor(nid("Cs137").nucid(), P::Inhale, S::Epa)
        );
        let mut mat = Material::new();
        mat.add_nuclide(nid("Cs137"), 2.0);
        let analytics = Analytics {
            masses: &Ame2020,
            decays: &ChainDecays,
        };
        let via_facade = mat
            .total_dose_per_g(&analytics, &provider, P::Inhale, S::Epa)
            .unwrap();
        let via_struct = mat
            .total_dose_per_g(&analytics, &DoseFactors, P::Inhale, S::Epa)
            .unwrap();
        assert!((via_facade - via_struct).abs() < 1e-18);
    }
}

#[cfg(test)]
mod ame_tests {
    use super::*;

    #[test]
    fn ame2020_provider_resolves_water() {
        // H2O from atom fractions with real masses
        let m = Material::from_atom_frac(
            &[
                (nucleide_nuclei::NuclideId::from_name("H1").unwrap(), 2.0),
                (nucleide_nuclei::NuclideId::from_name("O16").unwrap(), 1.0),
            ],
            &Ame2020,
            Some(1.0),
        )
        .unwrap();
        let af = m.atom_fractions(&Ame2020).unwrap();
        assert!(
            (af[&nucleide_nuclei::NuclideId::from_name("H1").unwrap()] - 2.0 / 3.0).abs() < 1e-12
        );
    }
}
