//! The [`Material`] composition model: construction, conversions, mixing,
//! arithmetic, and (de)serialization.

use std::collections::BTreeMap;
use std::ops::{Add, Div, Mul, Sub};

use nuclei::NuclideId;
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
/// real tables later is a single `impl MassProvider for nuclei::data::...`.
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

/// [`MassProvider`] backed by the AME2020 tables in `nuclei::data`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Ame2020;

impl MassProvider for Ame2020 {
    fn mass(&self, nucid: u32) -> Option<f64> {
        nuclei::data::atomic_mass(nucid)
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
// [`DecayProvider`] plus the existing [`MassProvider`]. Proper decay heat
// needs per-branch decay energies (mean beta/gamma/alpha energy release),
// which no table currently provides; add `decay_heat` once a decay-energy
// table lands in `nuclei::data` — the plumbing here is exactly this module's.
// ---------------------------------------------------------------------------

/// Avogadro constant, atoms per mole (exact, 2019 SI).
pub const AVOGADRO: f64 = 6.022_140_76e23;

/// One unified atomic mass unit in grams (2022 CODATA).
pub const GRAMS_PER_U: f64 = 1.660_539_068_92e-24;

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
/// `nuclei::data`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ChainDecays;

impl DecayProvider for ChainDecays {
    fn decay_constant(&self, nucid: u32) -> Option<f64> {
        nuclei::data::decay_constant(nucid)
    }
}

/// Provider bundle for radioanalytic quantities.
///
/// Both providers are needed at once — atom numbers come from masses,
/// activities from decay constants — so they travel together:
///
/// ```
/// use material::{Analytics, Ame2020, ChainDecays};
/// # let mut mat = material::Material::new();
/// # let co = nuclei::NuclideId::from_name("Co60").unwrap();
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

/// Errors from radioanalytics beyond the shared [`crate::Error`] set.
#[derive(Debug, Error)]
pub enum AnalyticsError {
    /// No decay data was available for a requested nuclide.
    #[error("no decay data available for nuclide `{0}`")]
    MissingDecay(NuclideId),
    /// An underlying composition failure (missing mass, degenerate total).
    #[error(transparent)]
    Core(#[from] crate::Error),
}

impl Material {
    /// Activity `A = λ·N` per nuclide, in becquerels.
    ///
    /// Atom counts follow from stored masses through `masses`
    /// (`N = m / (M · u)` with `u = 1.66053906892e-24 g`) and decay
    /// constants through `decays`. Fails with
    /// [`AnalyticsError::MissingDecay`] for nuclides without decay data and
    /// [`AnalyticsError::Core`] wrapping [`crate::Error::MissingMass`] when
    /// an atomic mass is unknown.
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
            let lambda = analytics
                .decays
                .decay_constant(id.nucid())
                .ok_or(AnalyticsError::MissingDecay(id))?;
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
        let lambda = LN_2 / nuclei::data::half_life(co60.nucid()).unwrap();
        assert_eq!(activity.keys().next().copied(), Some(co60));
        assert_eq!(
            ChainDecays.decay_constant(co60.nucid()),
            Some(lambda),
            "ChainDecays must be ln(2)/t_half of the tabulated half-life"
        );
        // N = m / (M · u) atoms for 1 g.
        let mass_u = nuclei::data::atomic_mass(co60.nucid()).unwrap();
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
        let t_half = nuclei::data::half_life(nucid).unwrap();
        assert!((lambda - LN_2 / t_half).abs() < 1e-18);
        // Stable Fe56 has no tabulated decay data.
        assert_eq!(ChainDecays.decay_constant(nid("Fe56").nucid()), None);
    }

    #[test]
    fn no_decay_provider_yields_missing_decay_error() {
        let mut mat = Material::new();
        mat.add_nuclide(nid("Co60"), 1.0);

        let analytics = Analytics {
            masses: &Ame2020,
            decays: &NoDecay,
        };
        match mat.activity(&analytics).unwrap_err() {
            AnalyticsError::MissingDecay(id) => assert_eq!(id, nid("Co60")),
            other => panic!("{other:?}"),
        }
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
}

#[cfg(test)]
mod ame_tests {
    use super::*;

    #[test]
    fn ame2020_provider_resolves_water() {
        // H2O from atom fractions with real masses
        let m = Material::from_atom_frac(
            &[
                (nuclei::NuclideId::from_name("H1").unwrap(), 2.0),
                (nuclei::NuclideId::from_name("O16").unwrap(), 1.0),
            ],
            &Ame2020,
            Some(1.0),
        )
        .unwrap();
        let af = m.atom_fractions(&Ame2020).unwrap();
        assert!((af[&nuclei::NuclideId::from_name("H1").unwrap()] - 2.0 / 3.0).abs() < 1e-12);
    }
}
