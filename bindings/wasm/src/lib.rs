//! Browser-facing WASM bindings for Nucleide.
//!
//! This crate is a thin facade over the workspace crates. It exposes the same
//! core capabilities as `nucleide._internal` but through a `wasm-bindgen` JS
//! API so tutorials can run live in the browser without a Python runtime.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use nucleide_nuclei::NuclideId;

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

fn js_err(e: impl std::fmt::Display) -> JsValue {
    JsValue::from_str(&e.to_string())
}

fn to_js<T: Serialize>(value: &T) -> Result<JsValue, JsValue> {
    // json_compatible serializes Rust maps as plain JS objects; the default
    // serializer emits JS `Map`s, which `Object.entries` in the UI cannot read.
    value
        .serialize(&serde_wasm_bindgen::Serializer::json_compatible())
        .map_err(js_err)
}

// ---------------------------------------------------------------------------
// Nuclide
// ---------------------------------------------------------------------------

/// A nuclide identifier with cross-code naming conversions.
#[wasm_bindgen]
pub struct WasmNuclide {
    inner: NuclideId,
}

#[derive(Serialize)]
struct NuclideInfo {
    name: String,
    nucid: u32,
    z: u32,
    a: u32,
    state: u32,
    zzaaam: u32,
    zaid: u32,
    zzllaaam: String,
    serpent: String,
    nist: String,
    cinder: u32,
    alara: String,
    sza: u32,
    mass: Option<f64>,
    abundance: Option<f64>,
}

#[wasm_bindgen]
impl WasmNuclide {
    /// Parse a nuclide name in any accepted dialect (e.g. `"Pu-241"`,
    /// `"241Pu"`, `"Ba137m"`).
    #[wasm_bindgen(constructor)]
    pub fn new(name: &str) -> Result<WasmNuclide, JsValue> {
        NuclideId::from_name(name)
            .map(|inner| WasmNuclide { inner })
            .map_err(js_err)
    }

    /// Build from a ZZAAAM integer (`942410` → Pu-241 ground state).
    #[wasm_bindgen(js_name = fromZzaaam)]
    pub fn from_zzaaam(v: u32) -> Result<WasmNuclide, JsValue> {
        NuclideId::from_zzaaam(v)
            .map(|inner| WasmNuclide { inner })
            .map_err(js_err)
    }

    /// Build from an HDF5/NNDC nucid string (`"942410"` → Pu-241 ground
    /// state).
    #[wasm_bindgen(js_name = fromNucid)]
    pub fn from_nucid(s: &str) -> Result<WasmNuclide, JsValue> {
        s.parse::<NuclideId>()
            .map(|inner| WasmNuclide { inner })
            .map_err(js_err)
    }

    /// Return all common identifiers and nuclear data as one JS object.
    #[wasm_bindgen(js_name = toObject)]
    pub fn to_object(&self) -> Result<JsValue, JsValue> {
        to_js(&NuclideInfo {
            name: self.inner.to_name(),
            nucid: self.inner.nucid(),
            z: self.inner.z(),
            a: self.inner.a(),
            state: self.inner.state(),
            zzaaam: self.inner.zzaaam(),
            zaid: nucleide_nuclei::dialects::to_zaid(self.inner),
            zzllaaam: nucleide_nuclei::dialects::zzllaaam(self.inner),
            serpent: nucleide_nuclei::dialects::serpent(self.inner),
            nist: nucleide_nuclei::dialects::nist(self.inner),
            cinder: nucleide_nuclei::dialects::to_cinder(self.inner),
            alara: nucleide_nuclei::dialects::alara(self.inner),
            sza: nucleide_nuclei::dialects::to_sza(self.inner),
            mass: nucleide_nuclei::data::atomic_mass(self.inner.nucid()),
            abundance: nucleide_nuclei::data::natural_abundance(self.inner.nucid()),
        })
    }

    /// Canonical GNDS name (e.g. `"Pu241"`).
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.inner.to_name()
    }

    /// Nucid integer (`Z*10_000 + A*10 + state`).
    #[wasm_bindgen(getter)]
    pub fn nucid(&self) -> u32 {
        self.inner.nucid()
    }

    /// Proton number Z.
    #[wasm_bindgen(getter)]
    pub fn z(&self) -> u32 {
        self.inner.z()
    }

    /// Mass number A.
    #[wasm_bindgen(getter)]
    pub fn a(&self) -> u32 {
        self.inner.a()
    }

    /// Metastable state index (0 = ground state).
    #[wasm_bindgen(getter)]
    pub fn state(&self) -> u32 {
        self.inner.state()
    }

    /// ZZAAAM integer (`Z*10_000 + A*10 + state`).
    #[wasm_bindgen(getter)]
    pub fn zzaaam(&self) -> u32 {
        self.inner.zzaaam()
    }

    /// MCNP ZAID integer (`922350` → U-235; metastables add `300 + 100*S`).
    #[wasm_bindgen(getter)]
    pub fn zaid(&self) -> u32 {
        nucleide_nuclei::dialects::to_zaid(self.inner)
    }

    /// ZZLLAAAM string `"ZZ-LL-AAAM"` with a lowercase isomer letter
    /// (e.g. `"95-Am-242m"`).
    #[wasm_bindgen(getter)]
    pub fn zzllaaam(&self) -> String {
        nucleide_nuclei::dialects::zzllaaam(self.inner)
    }

    /// Serpent-style name (`"Ll-AAA"`, e.g. `"Am-242m"`).
    #[wasm_bindgen(getter)]
    pub fn serpent(&self) -> String {
        nucleide_nuclei::dialects::serpent(self.inner)
    }

    /// NIST-style name: mass number then symbol, state dropped
    /// (e.g. `"242Am"`).
    #[wasm_bindgen(getter)]
    pub fn nist(&self) -> String {
        nucleide_nuclei::dialects::nist(self.inner)
    }

    /// Cinder `AAAZZZM` integer (`A*10_000 + Z*10 + state`).
    #[wasm_bindgen(getter)]
    pub fn cinder(&self) -> u32 {
        nucleide_nuclei::dialects::to_cinder(self.inner)
    }

    /// ALARA name: lowercase `"ll:AAA"` (e.g. `"pu:239"`).
    #[wasm_bindgen(getter)]
    pub fn alara(&self) -> String {
        nucleide_nuclei::dialects::alara(self.inner)
    }

    /// SZA integer (`state*10^6 + Z*10^3 + A`).
    #[wasm_bindgen(getter)]
    pub fn sza(&self) -> u32 {
        nucleide_nuclei::dialects::to_sza(self.inner)
    }

    /// Atomic mass [amu] from AME2020, or null when not tabulated.
    #[wasm_bindgen(getter)]
    pub fn mass(&self) -> Option<f64> {
        nucleide_nuclei::data::atomic_mass(self.inner.nucid())
    }

    /// Natural abundance [atom fraction], or null when not tabulated.
    #[wasm_bindgen(getter)]
    pub fn abundance(&self) -> Option<f64> {
        nucleide_nuclei::data::natural_abundance(self.inner.nucid())
    }

    /// FLUKA isotope name (e.g. `"235-U"`), erroring when the nuclide has no
    /// entry in the vendored FLUKA table.
    #[wasm_bindgen]
    pub fn fluka(&self) -> Result<String, JsValue> {
        nucleide_nuclei::dialects::id_to_fluka(self.inner)
            .map(String::from)
            .map_err(js_err)
    }
}

// ---------------------------------------------------------------------------
// Nuclear data helpers
// ---------------------------------------------------------------------------

fn resolve_nucid(key: &str) -> Result<NuclideId, JsValue> {
    key.parse::<NuclideId>().map_err(js_err)
}

/// Atomic mass [amu] for a nuclide name in any accepted dialect, or null
/// when not tabulated.
#[wasm_bindgen]
pub fn atomic_mass(key: &str) -> Result<Option<f64>, JsValue> {
    let id = resolve_nucid(key)?;
    Ok(nucleide_nuclei::data::atomic_mass(id.nucid()))
}

/// Natural abundance [atom fraction] for a nuclide name in any accepted
/// dialect, or null when not tabulated.
#[wasm_bindgen]
pub fn natural_abundance(key: &str) -> Result<Option<f64>, JsValue> {
    let id = resolve_nucid(key)?;
    Ok(nucleide_nuclei::data::natural_abundance(id.nucid()))
}

/// Half-life [s] for a nuclide name in any accepted dialect, or null for
/// stable/un-tabulated nuclides.
#[wasm_bindgen]
pub fn half_life(key: &str) -> Result<Option<f64>, JsValue> {
    let id = resolve_nucid(key)?;
    Ok(nucleide_nuclei::data::half_life(id.nucid()))
}

/// Decay constant [1/s] (`ln 2 / half-life`) for a nuclide name in any
/// accepted dialect, or null for stable/un-tabulated nuclides.
#[wasm_bindgen]
pub fn decay_constant(key: &str) -> Result<Option<f64>, JsValue> {
    let id = resolve_nucid(key)?;
    Ok(nucleide_nuclei::data::half_life(id.nucid()).map(|t| std::f64::consts::LN_2 / t))
}

/// One evaluated decay branch, JSON-shaped for the table-level mirror.
#[derive(serde::Serialize)]
struct JsDecayBranch {
    progeny: String,
    branching_fraction: f64,
    mode: String,
}

/// Evaluated decay branches for a nuclide (empty when the nuclide has no
/// branch rows, e.g. stable nuclides).
#[wasm_bindgen]
pub fn decay_branches(key: &str) -> Result<JsValue, JsValue> {
    let id = resolve_nucid(key)?;
    let branches: Vec<JsDecayBranch> = nucleide_nuclei::data::decay_branches(id.nucid())
        .unwrap_or_default()
        .into_iter()
        .map(|b| JsDecayBranch {
            progeny: nucleide_nuclei::NuclideId::from_nucid(b.progeny).to_name(),
            branching_fraction: b.branching_fraction,
            mode: b.mode.as_str().to_string(),
        })
        .collect();
    to_js(&branches)
}

/// Branching fraction from parent to progeny (GNDS names), if tabulated.
///
/// Named apart from the chain-scoped `branching_fraction` below.
#[wasm_bindgen]
pub fn decay_branch_fraction(parent: &str, progeny: &str) -> Result<Option<f64>, JsValue> {
    let p = resolve_nucid(parent)?;
    let d = resolve_nucid(progeny)?;
    Ok(nucleide_nuclei::data::branching_fraction(
        p.nucid(),
        d.nucid(),
    ))
}

/// Normalize a nuclide name in any accepted dialect to canonical GNDS form.
///
/// Accepts symbol-first (`Pu241`, `Pu-241`, `Ba137m`), mass-first (`241Pu`,
/// `40K`), and isomer suffix letters (`Ir-192n` → second isomer).
#[wasm_bindgen]
pub fn normalize_nuclide(name: &str) -> Result<String, JsValue> {
    nucleide_nuclei::dialects::normalize_nuclide_name(name)
        .map(|id| id.to_name())
        .map_err(js_err)
}

/// Q-value [MeV] of neutron capture for a nuclide name, or null when not
/// tabulated.
#[wasm_bindgen]
pub fn q_value_capture(key: &str) -> Result<Option<f64>, JsValue> {
    let id = resolve_nucid(key)?;
    Ok(nucleide_nuclei::data::q_value_neutron_capture(id.nucid()))
}

/// Q-value [MeV] of alpha decay for a nuclide name, or null when not
/// tabulated.
#[wasm_bindgen]
pub fn q_value_alpha(key: &str) -> Result<Option<f64>, JsValue> {
    let id = resolve_nucid(key)?;
    Ok(nucleide_nuclei::data::q_value_alpha(id.nucid()))
}

// ---------------------------------------------------------------------------
// Material
// ---------------------------------------------------------------------------

/// A nuclear material built from a chemical formula.
#[wasm_bindgen]
pub struct WasmMaterial {
    inner: nucleide_material::Material,
}

#[derive(Deserialize)]
struct MaterialPart {
    formula: String,
    fraction: f64,
}

#[wasm_bindgen]
impl WasmMaterial {
    /// Build a material from a formula such as "H2O" or "UO2" using natural
    /// abundances and AME2020 atomic masses.
    #[wasm_bindgen(constructor)]
    pub fn from_formula(formula: &str) -> Result<WasmMaterial, JsValue> {
        let mat = nucleide_material::Material::from_formula(
            formula,
            &nucleide_material::Ame2020,
            &nucleide_material::NaturalAbundances,
            None,
        )
        .map_err(js_err)?;
        Ok(WasmMaterial { inner: mat })
    }

    /// Build a material from atom fractions supplied as `{ "U235": 0.0072, ... }`.
    #[wasm_bindgen(js_name = fromAtomFrac)]
    pub fn from_atom_frac(atoms: JsValue) -> Result<WasmMaterial, JsValue> {
        let map: BTreeMap<String, f64> = serde_wasm_bindgen::from_value(atoms).map_err(js_err)?;
        let atoms: Vec<(NuclideId, f64)> = map
            .into_iter()
            .map(|(k, v)| Ok((k.parse::<NuclideId>().map_err(js_err)?, v)))
            .collect::<Result<Vec<_>, JsValue>>()?;
        let mat =
            nucleide_material::Material::from_atom_frac(&atoms, &nucleide_material::Ame2020, None)
                .map_err(js_err)?;
        Ok(WasmMaterial { inner: mat })
    }

    /// Mix several formula materials by relative mass fractions.
    /// `parts` is `[{ formula: "UO2", fraction: 0.9 }, { formula: "H2O", fraction: 0.1 }]`.
    #[wasm_bindgen(js_name = mixByMass)]
    pub fn mix_by_mass(parts: JsValue) -> Result<WasmMaterial, JsValue> {
        let parts: Vec<MaterialPart> = serde_wasm_bindgen::from_value(parts).map_err(js_err)?;
        let materials: Vec<WasmMaterial> = parts
            .iter()
            .map(|p| WasmMaterial::from_formula(&p.formula))
            .collect::<Result<Vec<_>, _>>()?;
        let refs: Vec<(&nucleide_material::Material, f64)> = materials
            .iter()
            .zip(parts.iter().map(|p| p.fraction))
            .map(|(m, f)| (&m.inner, f))
            .collect();
        let mixed = nucleide_material::Material::mix_by_mass(&refs).map_err(js_err)?;
        Ok(WasmMaterial { inner: mixed })
    }

    /// Total stored mass of the composition [g].
    #[wasm_bindgen(getter)]
    pub fn mass(&self) -> f64 {
        self.inner.mass()
    }

    /// Mass density [g/cm3] set on this material, if any.
    #[wasm_bindgen(getter)]
    pub fn density(&self) -> Option<f64> {
        self.inner.density()
    }

    /// `{ "U235": 0.0072, "U238": 0.9928, ... }` weight fractions.
    #[wasm_bindgen(js_name = weightFractions)]
    pub fn weight_fractions(&self) -> Result<JsValue, JsValue> {
        let frac = self.inner.weight_fractions().map_err(js_err)?;
        let map: BTreeMap<String, f64> =
            frac.into_iter().map(|(id, v)| (id.to_name(), v)).collect();
        to_js(&map)
    }

    /// `{ "U235": 0.0072, "U238": 0.9928, ... }` atom fractions.
    #[wasm_bindgen(js_name = atomFractions)]
    pub fn atom_fractions(&self) -> Result<JsValue, JsValue> {
        let frac = self
            .inner
            .atom_fractions(&nucleide_material::Ame2020)
            .map_err(js_err)?;
        let map: BTreeMap<String, f64> =
            frac.into_iter().map(|(id, v)| (id.to_name(), v)).collect();
        to_js(&map)
    }

    /// Serialize this material to OpenMC-style XML.
    #[wasm_bindgen(js_name = toXml)]
    pub fn to_xml(&self, name: &str, density: f64) -> Result<String, JsValue> {
        self.inner.to_xml(name, density, "g/cm3").map_err(js_err)
    }
}

// ---------------------------------------------------------------------------
// Materials Compendium
// ---------------------------------------------------------------------------

/// The DOE/PNNL Materials Compendium parsed from its JSON distribution.
#[wasm_bindgen]
pub struct WasmMaterialsCompendium {
    inner: nucleide_material::MaterialsLibrary,
}

#[derive(Serialize)]
struct CompendiumEntryInfo {
    name: String,
    acronym: Vec<String>,
    mat_num: u32,
    density: f64,
    atom_density: f64,
    source: String,
    comment: Vec<String>,
    weight_fractions: BTreeMap<String, f64>,
}

#[wasm_bindgen]
impl WasmMaterialsCompendium {
    /// Parse the compendium from its JSON text.
    #[wasm_bindgen(js_name = fromJson)]
    pub fn from_json(text: &str) -> Result<WasmMaterialsCompendium, JsValue> {
        let inner = nucleide_material::MaterialsLibrary::from_json(text).map_err(js_err)?;
        Ok(WasmMaterialsCompendium { inner })
    }

    /// Number of materials in the compendium.
    #[wasm_bindgen(getter)]
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// Whether the compendium holds no materials.
    #[wasm_bindgen(getter)]
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// All material names.
    pub fn names(&self) -> Result<JsValue, JsValue> {
        to_js(&self.inner.names())
    }

    /// One entry with its isotope-level weight fractions keyed by nuclide name.
    pub fn get(&self, name: &str) -> Result<JsValue, JsValue> {
        let entry = self
            .inner
            .get(name)
            .ok_or_else(|| js_err(format!("no compendium material named '{name}'")))?;
        let weight_fractions = entry
            .weight_fractions()
            .into_iter()
            .map(|(zaid, frac)| {
                let key = nucleide_nuclei::dialects::from_zaid(zaid)
                    .map(|id| id.to_name())
                    .unwrap_or_else(|_| zaid.to_string());
                (key, frac)
            })
            .collect();
        let info = CompendiumEntryInfo {
            name: entry.name.clone(),
            acronym: entry.acronym.clone(),
            mat_num: entry.mat_num,
            density: entry.density,
            atom_density: entry.atom_density,
            source: entry.source.clone(),
            comment: entry.comment.clone(),
            weight_fractions,
        };
        to_js(&info)
    }
}

// ---------------------------------------------------------------------------
// Enrichment cascade
// ---------------------------------------------------------------------------

/// An enrichment cascade in the M* parameterization, holding the feed,
/// product, and tails streams alongside the solved stage counts.
#[wasm_bindgen]
pub struct WasmCascade {
    inner: nucleide_enrichment::Cascade,
}

#[derive(Deserialize)]
struct CascadeConfig {
    alpha: f64,
    #[serde(rename = "Mstar")]
    mstar: f64,
    #[serde(rename = "enrichingKey")]
    enriching_key: String,
    #[serde(rename = "strippingKey")]
    stripping_key: String,
    #[serde(rename = "N")]
    n: f64,
    #[serde(rename = "M")]
    m: f64,
    #[serde(rename = "feedAssay")]
    feed_assay: f64,
    #[serde(rename = "productAssay")]
    product_assay: f64,
    #[serde(rename = "tailsAssay")]
    tails_assay: f64,
    feed: Option<BTreeMap<String, f64>>,
}

#[derive(Serialize)]
struct CascadeResult {
    alpha: f64,
    #[serde(rename = "Mstar")]
    mstar: f64,
    #[serde(rename = "feedAssay")]
    feed_assay: f64,
    #[serde(rename = "productAssay")]
    product_assay: f64,
    #[serde(rename = "tailsAssay")]
    tails_assay: f64,
    #[serde(rename = "stagesEnriching")]
    stages_enriching: f64,
    #[serde(rename = "stagesStripping")]
    stages_stripping: f64,
    #[serde(rename = "swuPerFeed")]
    swu_per_feed: f64,
    #[serde(rename = "swuPerProduct")]
    swu_per_product: f64,
    #[serde(rename = "productPerFeed")]
    product_per_feed: f64,
    #[serde(rename = "tailsPerFeed")]
    tails_per_feed: f64,
    feed: BTreeMap<String, f64>,
    product: BTreeMap<String, f64>,
    tails: BTreeMap<String, f64>,
}

#[derive(Serialize)]
struct StagePointJson {
    stage: u32,
    #[serde(rename = "assayJ")]
    assay_j: f64,
}

impl WasmCascade {
    fn to_result(&self) -> Result<JsValue, JsValue> {
        to_js(&CascadeResult {
            alpha: self.inner.alpha,
            mstar: self.inner.Mstar,
            feed_assay: self.inner.x_feed_j,
            product_assay: self.inner.x_prod_j,
            tails_assay: self.inner.x_tail_j,
            stages_enriching: self.inner.N,
            stages_stripping: self.inner.M,
            swu_per_feed: self.inner.swu_per_feed,
            swu_per_product: self.inner.swu_per_prod,
            product_per_feed: nucleide_enrichment::prod_per_feed(
                self.inner.x_feed_j,
                self.inner.x_prod_j,
                self.inner.x_tail_j,
            ),
            tails_per_feed: nucleide_enrichment::tail_per_feed(
                self.inner.x_feed_j,
                self.inner.x_prod_j,
                self.inner.x_tail_j,
            ),
            feed: stream_to_map(&self.inner.mat_feed),
            product: stream_to_map(&self.inner.mat_prod),
            tails: stream_to_map(&self.inner.mat_tail),
        })
    }
}

fn stream_to_map(stream: &nucleide_enrichment::Stream) -> BTreeMap<String, f64> {
    stream
        .comp
        .iter()
        .map(|(id, v)| (id.to_name(), *v))
        .collect()
}

#[wasm_bindgen]
impl WasmCascade {
    /// Default uranium cascade (0.72 % feed, 5 % product, 0.25 % tails).
    #[wasm_bindgen(js_name = defaultUranium)]
    pub fn default_uranium() -> WasmCascade {
        WasmCascade {
            inner: nucleide_enrichment::default_uranium_cascade(),
        }
    }

    /// Build a cascade from a configuration object.
    #[wasm_bindgen(constructor)]
    pub fn new(config: JsValue) -> Result<WasmCascade, JsValue> {
        let cfg: CascadeConfig = serde_wasm_bindgen::from_value(config).map_err(js_err)?;
        let j = cfg.enriching_key.parse::<NuclideId>().map_err(js_err)?;
        let k = cfg.stripping_key.parse::<NuclideId>().map_err(js_err)?;

        let mat_feed = if let Some(feed) = cfg.feed {
            nucleide_enrichment::Stream::from_comp(
                feed.into_iter()
                    .map(|(k, v)| Ok((k.parse::<NuclideId>().map_err(js_err)?, v)))
                    .collect::<Result<BTreeMap<_, _>, JsValue>>()?,
            )
        } else {
            nucleide_enrichment::default_uranium_cascade().mat_feed
        };

        let mut cascade = nucleide_enrichment::Cascade {
            alpha: cfg.alpha,
            Mstar: cfg.mstar,
            j,
            k,
            N: cfg.n,
            M: cfg.m,
            x_feed_j: cfg.feed_assay,
            x_prod_j: cfg.product_assay,
            x_tail_j: cfg.tails_assay,
            mat_feed,
            ..Default::default()
        };
        cascade.reset_xjs();
        Ok(WasmCascade { inner: cascade })
    }

    /// Solve the cascade in place.
    pub fn solve(&mut self) -> Result<(), JsValue> {
        self.inner = nucleide_enrichment::solve_numeric(
            &self.inner,
            nucleide_enrichment::DEFAULT_TOLERANCE,
            nucleide_enrichment::DEFAULT_MAX_ITER,
        )
        .map_err(js_err)?;
        Ok(())
    }

    /// Optimize Mstar to minimize total flow, then solve.
    #[wasm_bindgen(js_name = solveMulticomponent)]
    pub fn solve_multicomponent(&mut self) -> Result<(), JsValue> {
        self.inner = nucleide_enrichment::multicomponent(
            &self.inner,
            nucleide_enrichment::DEFAULT_TOLERANCE,
            nucleide_enrichment::DEFAULT_MAX_ITER,
        )
        .map_err(js_err)?;
        Ok(())
    }

    /// Return the full cascade state as a JS object.
    #[wasm_bindgen(js_name = toObject)]
    pub fn to_object(&self) -> Result<JsValue, JsValue> {
        self.to_result()
    }

    /// Stage separation factor α.
    #[wasm_bindgen(getter)]
    pub fn alpha(&self) -> f64 {
        self.inner.alpha
    }

    /// Feed assay of the enriching key (atom fraction).
    #[wasm_bindgen(getter, js_name = feedAssay)]
    pub fn feed_assay(&self) -> f64 {
        self.inner.x_feed_j
    }

    /// Product assay of the enriching key (atom fraction).
    #[wasm_bindgen(getter, js_name = productAssay)]
    pub fn product_assay(&self) -> f64 {
        self.inner.x_prod_j
    }

    /// Tails assay of the enriching key (atom fraction).
    #[wasm_bindgen(getter, js_name = tailsAssay)]
    pub fn tails_assay(&self) -> f64 {
        self.inner.x_tail_j
    }

    /// Number of enriching stages N.
    #[wasm_bindgen(getter, js_name = stagesEnriching)]
    pub fn stages_enriching(&self) -> f64 {
        self.inner.N
    }

    /// Number of stripping stages M.
    #[wasm_bindgen(getter, js_name = stagesStripping)]
    pub fn stages_stripping(&self) -> f64 {
        self.inner.M
    }

    /// Separative work per unit feed.
    #[wasm_bindgen(getter, js_name = swuPerFeed)]
    pub fn swu_per_feed(&self) -> f64 {
        self.inner.swu_per_feed
    }

    /// Separative work per unit product.
    #[wasm_bindgen(getter, js_name = swuPerProduct)]
    pub fn swu_per_product(&self) -> f64 {
        self.inner.swu_per_prod
    }

    /// Per-stage assay of the enriching key through the ideal cascade.
    #[wasm_bindgen(js_name = stageProfile)]
    pub fn stage_profile(&self) -> Result<JsValue, JsValue> {
        let profile = self.inner.stage_profile().map_err(js_err)?;
        let out: Vec<StagePointJson> = profile
            .into_iter()
            .map(|p| StagePointJson {
                stage: p.stage,
                assay_j: p.assay_j,
            })
            .collect();
        to_js(&out)
    }
}

// ---------------------------------------------------------------------------
// Depletion
// ---------------------------------------------------------------------------

/// A depletion chain: nuclides with their decay branches and neutron-induced
/// transmutation reactions.
#[wasm_bindgen]
pub struct WasmChain {
    inner: std::sync::Arc<nucleide_depletion::Chain>,
}

#[wasm_bindgen]
impl WasmChain {
    /// Parse a depletion-chain XML document from a string.
    #[wasm_bindgen(js_name = fromXml)]
    pub fn from_xml(xml: &str) -> Result<WasmChain, JsValue> {
        nucleide_depletion::Chain::from_xml(xml)
            .map(|inner| WasmChain {
                inner: std::sync::Arc::new(inner),
            })
            .map_err(js_err)
    }

    /// Nuclide names in chain (matrix) order.
    #[wasm_bindgen(getter)]
    pub fn nuclides(&self) -> Result<JsValue, JsValue> {
        let names: Vec<String> = self.inner.nuclides.iter().map(|n| n.name.clone()).collect();
        to_js(&names)
    }
}

/// Run one depletion step (CRAM or the analytic Bateman fast path).
///
/// `n0` is a JS object mapping nuclide names to atom counts. `rates` is a JS
/// object mapping `"Name:reaction"` strings to one-group rates [1/s]. `order`
/// is 16 or 48. `method` is an optional solver spelling (`"cram16"`,
/// `"cram48"`, `"bateman"`, `"bateman_hp"`, default `"cram48"`); an
/// explicitly non-default `method` overrides `order`. Bateman arms fall back
/// to CRAM-48 on non-decay systems.
#[wasm_bindgen]
pub fn deplete(
    chain: &WasmChain,
    n0: JsValue,
    dt: f64,
    rates: JsValue,
    order: u8,
    method: Option<String>,
) -> Result<JsValue, JsValue> {
    let n0: BTreeMap<String, f64> = serde_wasm_bindgen::from_value(n0).map_err(js_err)?;
    let rates: BTreeMap<String, f64> = serde_wasm_bindgen::from_value(rates).map_err(js_err)?;

    let method = resolve_method(order, method.as_deref())?;
    let reaction_rates = parse_reaction_rates(&chain.inner, &rates)?;

    let sys = nucleide_depletion::DepletionSystem::build((*chain.inner).clone(), &reaction_rates)
        .map_err(js_err)?;
    let result = nucleide_depletion::deplete_with_method(&sys, method, &n0, dt).map_err(js_err)?;
    to_js(&result.atoms)
}

fn parse_cram_order(order: u8) -> Result<nucleide_depletion::Order, JsValue> {
    match order {
        16 => Ok(nucleide_depletion::Order::Order16),
        48 => Ok(nucleide_depletion::Order::Order48),
        other => Err(js_err(format!(
            "unsupported CRAM order {other} (supported: 16, 48)"
        ))),
    }
}

/// Parse a solver `method` spelling; mirrors the Python `parse_method`.
fn parse_method(name: &str) -> Result<nucleide_depletion::Method, JsValue> {
    name.parse().map_err(js_err)
}

/// Resolve legacy `order` plus optional `method` into a core `Method`
/// (non-default `method` wins; absent/default defers to `order`).
fn resolve_method(order: u8, method: Option<&str>) -> Result<nucleide_depletion::Method, JsValue> {
    match method {
        None => parse_cram_order(order).map(nucleide_depletion::Method::Cram),
        Some(name) => {
            let parsed = parse_method(name)?;
            if parsed == nucleide_depletion::Method::default_cram() {
                parse_cram_order(order).map(nucleide_depletion::Method::Cram)
            } else {
                Ok(parsed)
            }
        }
    }
}

fn parse_integrator(name: &str) -> Result<nucleide_depletion::Integrator, JsValue> {
    if name.eq_ignore_ascii_case("predictor") {
        Ok(nucleide_depletion::Integrator::Predictor)
    } else if name.eq_ignore_ascii_case("cecm") {
        Ok(nucleide_depletion::Integrator::Cecm)
    } else if name.eq_ignore_ascii_case("cf4") {
        Ok(nucleide_depletion::Integrator::Cf4)
    } else {
        Err(js_err(format!(
            "unsupported integrator `{name}` (supported: predictor, cecm, cf4)"
        )))
    }
}

fn parse_reaction_rates(
    chain: &nucleide_depletion::Chain,
    rates: &BTreeMap<String, f64>,
) -> Result<nucleide_depletion::ReactionRates, JsValue> {
    let mut reaction_rates = nucleide_depletion::ReactionRates::new();
    for (key, v) in rates {
        let (nuc, rx) = key
            .split_once(':')
            .ok_or_else(|| js_err(format!("rate key `{key}` must be `Name:reaction`")))?;
        let idx = chain
            .index_of(nuc)
            .ok_or_else(|| js_err(format!("rate for unknown nuclide `{nuc}`")))?;
        reaction_rates
            .entry(idx)
            .or_default()
            .insert(rx.to_string(), *v);
    }
    Ok(reaction_rates)
}

#[derive(Serialize)]
struct DepleteSeriesResult {
    times: Vec<f64>,
    atoms: Vec<BTreeMap<String, f64>>,
    activity: Vec<BTreeMap<String, f64>>,
    decay_heat: Vec<BTreeMap<String, f64>>,
}

/// Run a multi-step depletion series with activity and decay heat.
///
/// `n0` maps nuclide names to atom counts; `dts` is the list of step lengths
/// [s]; `rates` maps `"Name:reaction"` strings to one-group rates [1/s] and
/// applies unchanged to every step (one `Step` per `dt` with the same rates).
/// `integrator` is `"predictor"`, `"cecm"`, or `"cf4"`; `order` is 16 or 48;
/// `method` is an optional solver spelling (`"cram16"`, `"cram48"`,
/// `"bateman"`, `"bateman_hp"`) that overrides `order` when non-default.
/// Bateman steps with live rates fall back to CRAM-48.
///
/// Returns `{ times, atoms, activity, decay_heat }` where `times` holds the
/// cumulative nodes starting at `t = 0` (so every list has `dts.len() + 1`
/// entries and row 0 echoes the initial state) and the other three are
/// per-node maps keyed by nuclide name (`activity` in Bq, `decay_heat` in W).
#[wasm_bindgen(js_name = depleteSeries)]
pub fn deplete_series(
    chain: &WasmChain,
    n0: JsValue,
    dts: JsValue,
    rates: JsValue,
    integrator: &str,
    order: u8,
    method: Option<String>,
) -> Result<JsValue, JsValue> {
    let n0: BTreeMap<String, f64> = serde_wasm_bindgen::from_value(n0).map_err(js_err)?;
    let dts: Vec<f64> = serde_wasm_bindgen::from_value(dts).map_err(js_err)?;
    let rates: BTreeMap<String, f64> = serde_wasm_bindgen::from_value(rates).map_err(js_err)?;

    let method = resolve_method(order, method.as_deref())?;
    let integrator = parse_integrator(integrator)?;
    let reaction_rates = parse_reaction_rates(&chain.inner, &rates)?;

    let mut n0_vec = vec![0.0; chain.inner.nuclides.len()];
    for (name, value) in &n0 {
        let idx = chain
            .inner
            .index_of(name)
            .ok_or_else(|| js_err(format!("n0 has unknown nuclide `{name}`")))?;
        n0_vec[idx] = *value;
    }

    let steps: Vec<nucleide_depletion::Step> = dts
        .iter()
        .map(|dt| nucleide_depletion::Step::new(*dt, reaction_rates.clone()))
        .collect();
    let sys = nucleide_depletion::DepletionSystem::build((*chain.inner).clone(), &reaction_rates)
        .map_err(js_err)?;
    let series =
        nucleide_depletion::integrate_with_method(&sys, &n0_vec, &steps, integrator, method)
            .map_err(js_err)?;

    let names: Vec<String> = chain
        .inner
        .nuclides
        .iter()
        .map(|nuc| nuc.name.clone())
        .collect();
    let key_rows = |rows: &[Vec<f64>]| -> Vec<BTreeMap<String, f64>> {
        rows.iter()
            .map(|row| {
                names
                    .iter()
                    .cloned()
                    .zip(row.iter().copied())
                    .collect::<BTreeMap<String, f64>>()
            })
            .collect()
    };
    to_js(&DepleteSeriesResult {
        times: series.times,
        atoms: key_rows(&series.atoms),
        activity: key_rows(&series.activity),
        decay_heat: key_rows(&series.decay_heat),
    })
}

// ---------------------------------------------------------------------------
// MCNP I/O
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct McnpMaterialJson {
    number: u32,
    fractions: BTreeMap<String, f64>,
    #[serde(rename = "fractionType")]
    fraction_type: String,
    density: Option<f64>,
    comments: Vec<String>,
}

/// Parse MCNP material cards (`mX`) from deck text into a JSON summary.
#[wasm_bindgen(js_name = parseMcnpMaterials)]
pub fn parse_mcnp_materials(text: &str) -> Result<JsValue, JsValue> {
    let mats = nucleide_mcnp_io::inp::materials_from_str(text).map_err(js_err)?;
    let out: Vec<McnpMaterialJson> = mats
        .into_iter()
        .map(|m| McnpMaterialJson {
            number: m.number,
            fractions: m
                .fractions
                .into_iter()
                .map(|(id, v)| (id.to_name(), v))
                .collect(),
            fraction_type: match m.fraction_type {
                nucleide_mcnp_io::inp::FracKind::Atom => "atom".to_string(),
                nucleide_mcnp_io::inp::FracKind::Mass => "mass".to_string(),
            },
            density: m.density,
            comments: m.comments,
        })
        .collect();
    to_js(&out)
}

#[derive(Serialize)]
struct XsdirTableJson {
    name: String,
    zaid: String,
    #[serde(rename = "serpentType")]
    serpent_type: Option<String>,
    awr: f64,
    filename: String,
    #[serde(rename = "fileType")]
    file_type: i64,
    temperature: Option<f64>,
    metastable: Option<bool>,
}

#[derive(Serialize)]
struct XsdirSummary {
    datapath: Option<String>,
    #[serde(rename = "awrCount")]
    awr_count: usize,
    #[serde(rename = "tableCount")]
    table_count: usize,
    tables: Vec<XsdirTableJson>,
}

/// Parse an MCNP `xsdir` cross-section directory into a JSON summary.
#[wasm_bindgen(js_name = parseXsdir)]
pub fn parse_xsdir(text: &str) -> Result<JsValue, JsValue> {
    let xsdir = nucleide_mcnp_io::xsdir::Xsdir::parse(text).map_err(js_err)?;
    let tables: Vec<XsdirTableJson> = xsdir
        .tables
        .into_iter()
        .map(|t| XsdirTableJson {
            name: t.name.clone(),
            zaid: t.zaid().to_string(),
            serpent_type: t.serpent_type().map(|v| v.to_string()),
            awr: t.awr,
            filename: t.filename.clone(),
            file_type: t.filetype,
            temperature: t.temperature,
            metastable: t.metastable(),
        })
        .collect();
    to_js(&XsdirSummary {
        datapath: xsdir.datapath,
        awr_count: xsdir.awr.len(),
        table_count: tables.len(),
        tables,
    })
}

#[derive(Serialize)]
struct MeshTallySummary {
    #[serde(rename = "tallyNumber")]
    tally_number: u32,
    particle: char,
    #[serde(rename = "doseResponse")]
    dose_response: bool,
    dims: Vec<usize>,
    #[serde(rename = "numVes")]
    num_ves: usize,
    #[serde(rename = "numEGroups")]
    num_e_groups: usize,
    #[serde(rename = "xBounds")]
    x_bounds: Vec<f64>,
    #[serde(rename = "yBounds")]
    y_bounds: Vec<f64>,
    #[serde(rename = "zBounds")]
    z_bounds: Vec<f64>,
    #[serde(rename = "eBounds")]
    e_bounds: Vec<f64>,
    result: Vec<Vec<f64>>,
    #[serde(rename = "relError")]
    rel_error: Vec<Vec<f64>>,
    #[serde(rename = "totalResult")]
    total_result: Vec<f64>,
    #[serde(rename = "totalRelError")]
    total_rel_error: Vec<f64>,
}

#[derive(Serialize)]
struct MeshtalSummary {
    version: String,
    title: String,
    histories: u64,
    #[serde(rename = "tallyCount")]
    tally_count: usize,
    tallies: BTreeMap<String, MeshTallySummary>,
}

/// Parse an MCNP `meshtal` file into a JSON summary of its mesh tallies.
#[wasm_bindgen(js_name = parseMeshtal)]
pub fn parse_meshtal(text: &str) -> Result<JsValue, JsValue> {
    let meshtal = nucleide_mcnp_io::meshtal::Meshtal::parse(text).map_err(js_err)?;
    let tallies: BTreeMap<String, MeshTallySummary> = meshtal
        .tallies
        .into_iter()
        .map(|(num, t)| {
            let summary = MeshTallySummary {
                tally_number: t.tally_number,
                particle: t.particle.letter(),
                dose_response: t.dose_response,
                dims: t.dims().to_vec(),
                num_ves: t.num_ves(),
                num_e_groups: t.num_e_groups(),
                x_bounds: t.x_bounds,
                y_bounds: t.y_bounds,
                z_bounds: t.z_bounds,
                e_bounds: t.e_bounds,
                result: t.result,
                rel_error: t.rel_error,
                total_result: t.total_result,
                total_rel_error: t.total_rel_error,
            };
            (num.to_string(), summary)
        })
        .collect();
    to_js(&MeshtalSummary {
        version: meshtal.version,
        title: meshtal.title,
        histories: meshtal.histories,
        tally_count: tallies.len(),
        tallies,
    })
}

#[derive(Serialize)]
struct WwinpSummary {
    ni: u32,
    nr: u32,
    ne: Vec<u32>,
    nf: [u32; 3],
    origin: [f64; 3],
    nc: [u32; 3],
    bounds: Vec<Vec<f64>>,
    e: Vec<Vec<f64>>,
    ww: Vec<Vec<Vec<f64>>>,
}

/// Parse an MCNP `WWINP` weight-window mesh into a JSON summary.
#[wasm_bindgen(js_name = parseWwinp)]
pub fn parse_wwinp(text: &str) -> Result<JsValue, JsValue> {
    let wwinp = nucleide_mcnp_io::wwinp::Wwinp::parse(text).map_err(js_err)?;
    to_js(&WwinpSummary {
        ni: wwinp.ni,
        nr: wwinp.nr,
        ne: wwinp.ne,
        nf: wwinp.nf,
        origin: wwinp.origin,
        nc: wwinp.nc,
        bounds: wwinp.bounds,
        e: wwinp.e,
        ww: wwinp.ww,
    })
}

// ---------------------------------------------------------------------------
// Variance reduction
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct MagicSummary {
    #[serde(rename = "lowerBoundsWw")]
    lower_bounds_ww: Vec<f64>,
    #[serde(rename = "groupsPerVe")]
    groups_per_ve: usize,
    #[serde(rename = "scaleFactors")]
    scale_factors: Vec<f64>,
    #[serde(rename = "eUpperBounds")]
    e_upper_bounds: Vec<f64>,
    #[serde(rename = "wwTagName")]
    ww_tag_name: String,
    #[serde(rename = "eUpperBoundsTagName")]
    e_upper_bounds_tag_name: String,
}

/// Compute MAGIC-method weight-window lower bounds from one mesh tally of a
/// meshtal file.
///
/// `tally_number` selects the tally within `meshtal_text`; `selection` is
/// `"total"` or `"perGroup"`.
#[wasm_bindgen(js_name = magicBounds)]
pub fn magic_bounds(
    meshtal_text: &str,
    tally_number: u32,
    selection: &str,
    tolerance: f64,
    null_value: f64,
) -> Result<JsValue, JsValue> {
    let meshtal = nucleide_mcnp_io::meshtal::Meshtal::parse(meshtal_text).map_err(js_err)?;
    let tally = meshtal
        .tallies
        .get(&tally_number)
        .ok_or_else(|| js_err(format!("tally {tally_number} not found")))?;
    let sel = match selection {
        "total" => nucleide_vr_tools::magic::MagicSelection::Total,
        "perGroup" => nucleide_vr_tools::magic::MagicSelection::PerGroup,
        _ => return Err(js_err("selection must be 'total' or 'perGroup'")),
    };
    let out = nucleide_vr_tools::magic::magic_with(
        tally,
        sel,
        nucleide_vr_tools::magic::MagicParams {
            tolerance,
            null_value,
        },
    )
    .map_err(js_err)?;
    to_js(&MagicSummary {
        lower_bounds_ww: out.lower_bounds_ww,
        groups_per_ve: out.groups_per_ve,
        scale_factors: out.scale_factors,
        e_upper_bounds: out.e_upper_bounds,
        ww_tag_name: out.ww_tag_name,
        e_upper_bounds_tag_name: out.e_upper_bounds_tag_name,
    })
}

/// Sample an index from a discrete probability mass function using the alias
/// method; `r1` and `r2` are independent uniform randoms in `[0, 1)`.
#[wasm_bindgen(js_name = aliasTableSample)]
pub fn alias_table_sample(pdf: Vec<f64>, r1: f64, r2: f64) -> Result<usize, JsValue> {
    let table = nucleide_vr_tools::sampling::AliasTable::new(&pdf).map_err(js_err)?;
    Ok(table.sample(r1, r2))
}

#[derive(Serialize)]
struct SampledVoxelSummary {
    index: usize,
    i: usize,
    j: usize,
    k: usize,
    weight: f64,
}

/// Sample a source voxel from one mesh tally of a meshtal file.
///
/// `tally_number` selects the tally within `meshtal_text`; `mode` is
/// `"analog"` (sample tally values) or `"uniform"` (sample voxel volumes);
/// `r1` and `r2` are independent uniform randoms in `[0, 1)`.
#[wasm_bindgen(js_name = meshSourceSample)]
pub fn mesh_source_sample(
    meshtal_text: &str,
    tally_number: u32,
    mode: &str,
    r1: f64,
    r2: f64,
) -> Result<JsValue, JsValue> {
    let meshtal = nucleide_mcnp_io::meshtal::Meshtal::parse(meshtal_text).map_err(js_err)?;
    let tally = meshtal
        .tallies
        .get(&tally_number)
        .ok_or_else(|| js_err(format!("tally {tally_number} not found")))?;
    let mode = match mode {
        "analog" => nucleide_vr_tools::sampling::Mode::Analog,
        "uniform" => nucleide_vr_tools::sampling::Mode::Uniform,
        _ => return Err(js_err("mode must be 'analog' or 'uniform'")),
    };
    let sampler =
        nucleide_vr_tools::sampling::MeshSourceSampler::new(tally, mode, None).map_err(js_err)?;
    let sample = sampler.sample(r1, r2);
    to_js(&SampledVoxelSummary {
        index: sample.index,
        i: sample.i,
        j: sample.j,
        k: sample.k,
        weight: sample.weight,
    })
}

// ---------------------------------------------------------------------------
// Activation (ALARA deck/output, FISPACT-II output, R2S workflow)
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct AlaraMixtureJson {
    name: String,
    entries: Vec<String>,
}

#[derive(Serialize)]
struct AlaraFluxJson {
    name: String,
    file: String,
    scale: f64,
    skip: usize,
    format: String,
}

#[derive(Serialize)]
struct AlaraScheduleJson {
    name: String,
    items: Vec<Vec<String>>,
}

#[derive(Serialize)]
struct AlaraDeckSummary {
    block_kinds: Vec<String>,
    mixtures: Vec<AlaraMixtureJson>,
    fluxes: Vec<AlaraFluxJson>,
    cooling_times_s: Vec<f64>,
    schedules: Vec<AlaraScheduleJson>,
}

fn mixture_entry_text(entry: &nucleide_alara_io::deck::MixtureEntry) -> String {
    use nucleide_alara_io::deck::MixtureEntry as E;
    match entry {
        E::Material {
            name,
            rel_density,
            vol_fraction,
        } => format!("material {name} {rel_density} {vol_fraction}"),
        E::Element {
            symbol,
            rel_density,
            vol_fraction,
        } => format!("element {symbol} {rel_density} {vol_fraction}"),
        E::Like {
            mixture,
            rel_density,
        } => format!("like {mixture} {rel_density}"),
        E::Target { target_kind, name } => format!("target {target_kind} {name}"),
    }
}

/// Parse an ALARA input deck into a JSON summary.
#[wasm_bindgen(js_name = parseAlaraDeck)]
pub fn parse_alara_deck(text: &str) -> Result<JsValue, JsValue> {
    let deck = nucleide_alara_io::AlaraDeck::parse(text).map_err(js_err)?;
    to_js(&AlaraDeckSummary {
        block_kinds: deck
            .block_kinds()
            .iter()
            .map(|kind| kind.to_string())
            .collect(),
        mixtures: deck
            .mixtures
            .iter()
            .map(|mix| AlaraMixtureJson {
                name: mix.name.clone(),
                entries: mix.entries.iter().map(mixture_entry_text).collect(),
            })
            .collect(),
        fluxes: deck
            .fluxes
            .iter()
            .map(|flux| AlaraFluxJson {
                name: flux.name.clone(),
                file: flux.file.clone(),
                scale: flux.scale,
                skip: flux.skip,
                format: flux.format.clone(),
            })
            .collect(),
        cooling_times_s: deck
            .cooling
            .as_ref()
            .map(|cooling| cooling.times_s.clone())
            .unwrap_or_default(),
        schedules: deck
            .schedules
            .iter()
            .map(|schedule| AlaraScheduleJson {
                name: schedule.name.clone(),
                items: schedule
                    .items
                    .iter()
                    .map(|item| item.tokens.clone())
                    .collect(),
            })
            .collect(),
    })
}

#[derive(Serialize)]
struct ResponseRowJson {
    time_s: f64,
    time_label: String,
    nuclide: String,
    half_life_s: f64,
    run_lbl: String,
    block: String,
    block_name: String,
    block_num: i64,
    variable: String,
    var_unit: String,
    value: f64,
}

fn response_row_json(row: &nucleide_alara_io::output::ResponseRow) -> ResponseRowJson {
    ResponseRowJson {
        time_s: row.time_s,
        time_label: row.time_label.clone(),
        nuclide: row.nuclide.clone(),
        half_life_s: row.half_life_s,
        run_lbl: row.run_lbl.clone(),
        block: row.block.as_str().to_string(),
        block_name: row.block_name.clone(),
        block_num: row.block_num,
        variable: row.variable.as_str().to_string(),
        var_unit: row.var_unit.clone(),
        value: row.value,
    }
}

fn distinct_sorted(values: impl Iterator<Item = String>) -> Vec<String> {
    let set: std::collections::BTreeSet<String> = values.collect();
    set.into_iter().collect()
}

#[derive(Serialize)]
struct AlaraOutputSummary {
    rows: Vec<ResponseRowJson>,
    variables: Vec<String>,
    blocks: Vec<String>,
}

/// Parse an ALARA activation-output listing into a JSON summary.
#[wasm_bindgen(js_name = parseAlaraOutput)]
pub fn parse_alara_output(text: &str, run_lbl: &str) -> Result<JsValue, JsValue> {
    let frame = nucleide_alara_io::output::ResponseFrame::parse(text, run_lbl).map_err(js_err)?;
    to_js(&AlaraOutputSummary {
        rows: frame.rows.iter().map(response_row_json).collect(),
        variables: distinct_sorted(
            frame
                .rows
                .iter()
                .map(|row| row.variable.as_str().to_string()),
        ),
        blocks: distinct_sorted(frame.rows.iter().map(|row| row.block.as_str().to_string())),
    })
}

#[derive(Serialize)]
struct FispactOutputSummary {
    rows: Vec<ResponseRowJson>,
    variables: Vec<String>,
}

/// Parse a FISPACT-II inventory listing into a JSON summary.
#[wasm_bindgen(js_name = parseFispactOutput)]
pub fn parse_fispact_output(text: &str, run_lbl: &str) -> Result<JsValue, JsValue> {
    let frame = nucleide_fispact_io::parse_to_frame(text, run_lbl).map_err(js_err)?;
    to_js(&FispactOutputSummary {
        rows: frame.rows.iter().map(response_row_json).collect(),
        variables: distinct_sorted(
            frame
                .rows
                .iter()
                .map(|row| row.variable.as_str().to_string()),
        ),
    })
}

#[derive(Serialize)]
struct R2sStepJson {
    zone: String,
    flux: String,
}

#[derive(Serialize)]
struct R2sSummary {
    steps: Vec<R2sStepJson>,
    cooling_s: Vec<f64>,
    top_schedule: String,
    total_s: f64,
}

/// Pick the flux block for `zone`: broadcast a lone flux, else match by name.
///
/// Mirrors `nucleide_r2s::R2sWorkflow::from_deck`, which is not depended on here
/// because the `r2s` crate enables the `depletion`/`rayon` feature that this
/// `wasm32-unknown-unknown` build keeps disabled.
fn r2s_resolve_flux(
    fluxes: &[nucleide_alara_io::deck::FluxDef],
    zone: &str,
) -> Result<String, String> {
    if fluxes.len() == 1 {
        return Ok(fluxes[0].name.clone());
    }
    fluxes
        .iter()
        .find(|flux| flux.name == zone)
        .map(|flux| flux.name.clone())
        .ok_or_else(|| {
            format!(
                "zone `{zone}` matches no flux block ({} flux blocks, no name match)",
                fluxes.len()
            )
        })
}

/// Top-schedule discovery mirroring `nucleide_r2s::R2sWorkflow::from_deck`.
fn r2s_top_schedule(deck: &nucleide_alara_io::deck::AlaraDeck) -> Result<String, String> {
    if deck.schedules.is_empty() {
        return Err("deck defines no schedules".to_string());
    }
    if deck.schedules.len() == 1 {
        return Ok(deck.schedules[0].name.clone());
    }
    let referenced: std::collections::HashSet<&str> = deck
        .schedules
        .iter()
        .flat_map(|schedule| schedule.items.iter())
        .filter(|item| item.tokens.len() == 4)
        .map(|item| item.tokens[0].as_str())
        .collect();
    let tops: Vec<&str> = deck
        .schedules
        .iter()
        .map(|schedule| schedule.name.as_str())
        .filter(|name| !referenced.contains(name))
        .collect();
    match tops.as_slice() {
        [top] => Ok((*top).to_string()),
        [] => Err(
            "no top-level schedule: every schedule is referenced (possible recursion)".to_string(),
        ),
        _ => Err(format!(
            "ambiguous top-level schedules: {}",
            tops.join(", ")
        )),
    }
}

/// Convert one raw deck schedule item (4- or 6-token form) for expansion.
fn r2s_sched_item(tokens: &[String], line: usize) -> Result<nucleide_alara_io::SchedItem, String> {
    match tokens {
        [op_text, op_unit, flux, history, delay_text, delay_unit] => {
            let op: f64 = op_text
                .parse()
                .map_err(|_| format!("line {line}: expected operating time, found `{op_text}`"))?;
            let delay: f64 = delay_text
                .parse()
                .map_err(|_| format!("line {line}: expected delay, found `{delay_text}`"))?;
            Ok(nucleide_alara_io::SchedItem::Pulse {
                op_time_s: nucleide_alara_io::parse_time_to_seconds(op, op_unit)
                    .map_err(|e| e.to_string())?,
                flux: flux.clone(),
                history: history.clone(),
                delay_s: nucleide_alara_io::parse_time_to_seconds(delay, delay_unit)
                    .map_err(|e| e.to_string())?,
            })
        }
        [name, history, delay_text, delay_unit] => {
            let delay: f64 = delay_text
                .parse()
                .map_err(|_| format!("line {line}: expected delay, found `{delay_text}`"))?;
            Ok(nucleide_alara_io::SchedItem::SubSchedule {
                name: name.clone(),
                history: history.clone(),
                delay_s: nucleide_alara_io::parse_time_to_seconds(delay, delay_unit)
                    .map_err(|e| e.to_string())?,
            })
        }
        _ => Err(format!(
            "line {line}: expected 4- or 6-token schedule item, found {}",
            tokens.join(" ")
        )),
    }
}

/// Derive an R2S workflow summary from an ALARA deck.
#[wasm_bindgen(js_name = r2sFromDeck)]
pub fn r2s_from_deck(text: &str) -> Result<JsValue, JsValue> {
    let deck = nucleide_alara_io::AlaraDeck::parse(text).map_err(js_err)?;
    let loading = deck
        .mat_loading
        .as_ref()
        .ok_or_else(|| js_err("deck defines no `mat_loading` block"))?;
    let zones: Vec<&str> = loading
        .entries
        .iter()
        .filter(|entry| !entry.mixture.eq_ignore_ascii_case("void"))
        .map(|entry| entry.zone.as_str())
        .collect();
    if zones.is_empty() {
        return Err(js_err("deck defines no non-`void` zones in `mat_loading`"));
    }
    if deck.fluxes.is_empty() {
        return Err(js_err("deck defines no `flux` blocks"));
    }
    let steps: Vec<R2sStepJson> = zones
        .iter()
        .map(|zone| {
            Ok(R2sStepJson {
                zone: (*zone).to_string(),
                flux: r2s_resolve_flux(&deck.fluxes, zone).map_err(js_err)?,
            })
        })
        .collect::<Result<Vec<_>, JsValue>>()?;
    let cooling_s = deck
        .cooling
        .as_ref()
        .map(|cooling| cooling.times_s.clone())
        .unwrap_or_default();
    let top_schedule = r2s_top_schedule(&deck).map_err(js_err)?;

    let mut schedules = Vec::with_capacity(deck.schedules.len());
    for raw in &deck.schedules {
        let mut items = Vec::with_capacity(raw.items.len());
        for entry in &raw.items {
            items.push(r2s_sched_item(&entry.tokens, entry.line).map_err(js_err)?);
        }
        schedules.push(nucleide_alara_io::schedule::ScheduleDef {
            name: raw.name.clone(),
            items,
        });
    }
    let histories: Vec<nucleide_alara_io::schedule::PulseHistory> = deck
        .pulse_histories
        .iter()
        .map(|history| nucleide_alara_io::schedule::PulseHistory {
            name: history.name.clone(),
            levels: history
                .levels
                .iter()
                .map(|level| nucleide_alara_io::schedule::PulseLevel {
                    count: level.pulses,
                    delay_s: level.delay_s,
                })
                .collect(),
        })
        .collect();
    let flat = nucleide_alara_io::schedule::expand_from(&top_schedule, &schedules, &histories)
        .map_err(js_err)?;
    to_js(&R2sSummary {
        steps,
        cooling_s,
        top_schedule,
        total_s: nucleide_alara_io::total_time(&flat),
    })
}

// ---------------------------------------------------------------------------
// Serpent output (_res / _dep / _det)
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct SerpentVariableJson {
    name: String,
    kind: &'static str,
    shape: String,
    value: Option<serde_json::Value>,
}

fn serpent_scalar_json(value: &nucleide_serpent_io::Value) -> serde_json::Value {
    match value {
        nucleide_serpent_io::Value::Num(x) => serde_json::json!(x),
        nucleide_serpent_io::Value::Str(s) => serde_json::json!(s),
    }
}

fn serpent_variable_json(name: &str, entry: &nucleide_serpent_io::Entry) -> SerpentVariableJson {
    match entry {
        nucleide_serpent_io::Entry::Scalar(v) => SerpentVariableJson {
            name: name.to_string(),
            kind: "scalar",
            shape: String::new(),
            value: Some(serpent_scalar_json(v)),
        },
        nucleide_serpent_io::Entry::Vector(vals) => SerpentVariableJson {
            name: name.to_string(),
            kind: "vector",
            shape: format!("[{}]", vals.len()),
            value: None,
        },
        nucleide_serpent_io::Entry::Matrix(m) => SerpentVariableJson {
            name: name.to_string(),
            kind: "matrix",
            shape: format!("[{}×{}]", m.rows(), m.cols()),
            value: None,
        },
    }
}

fn serpent_variables(table: &nucleide_serpent_io::Table) -> Vec<SerpentVariableJson> {
    table
        .iter()
        .map(|(name, entry)| serpent_variable_json(name, entry))
        .collect()
}

#[derive(Serialize)]
struct SerpentResSummary {
    variable_count: usize,
    version: Option<String>,
    title: Option<String>,
    keff: Option<Vec<f64>>,
    variables: Vec<SerpentVariableJson>,
}

/// Parse a Serpent `_res.m` results file into a JSON summary.
#[wasm_bindgen(js_name = parseSerpentRes)]
pub fn parse_serpent_res(text: &str) -> Result<JsValue, JsValue> {
    let table = nucleide_serpent_io::parse_res(text).map_err(js_err)?;
    let version = table
        .get_vec_str("VERSION")
        .ok()
        .and_then(|v| v.first().cloned());
    let title = table
        .get_vec_str("TITLE")
        .ok()
        .and_then(|v| v.first().map(|s| s.trim().to_string()));
    let keff = table
        .get_matrix("IMP_KEFF")
        .ok()
        .and_then(|m| m.row_f64(0).ok())
        .map(|row| row.iter().take(2).copied().collect());
    to_js(&SerpentResSummary {
        variable_count: table.len(),
        version,
        title,
        keff,
        variables: serpent_variables(&table),
    })
}

#[derive(Serialize)]
struct SerpentDepSummary {
    variable_count: usize,
    nuclides: Vec<String>,
    zai: Vec<f64>,
    variables: Vec<SerpentVariableJson>,
}

/// Parse a Serpent `_dep.m` depletion file into a JSON summary.
#[wasm_bindgen(js_name = parseSerpentDep)]
pub fn parse_serpent_dep(text: &str) -> Result<JsValue, JsValue> {
    let table = nucleide_serpent_io::parse_dep(text).map_err(js_err)?;
    let nuclides = table
        .get_vec_str("NAMES")
        .map(|names| names.iter().map(|s| s.trim().to_string()).collect())
        .unwrap_or_default();
    let zai = table.get_vec_f64("ZAI").unwrap_or_default();
    to_js(&SerpentDepSummary {
        variable_count: table.len(),
        nuclides,
        zai,
        variables: serpent_variables(&table),
    })
}

#[derive(Serialize)]
struct SerpentDetSummary {
    variable_count: usize,
    detectors: Vec<String>,
    variables: Vec<SerpentVariableJson>,
}

/// Parse a Serpent `_det.m` detector file into a JSON summary.
#[wasm_bindgen(js_name = parseSerpentDet)]
pub fn parse_serpent_det(text: &str) -> Result<JsValue, JsValue> {
    let table = nucleide_serpent_io::parse_det(text).map_err(js_err)?;
    let detectors = table
        .iter()
        .filter(|(name, entry)| {
            name.starts_with("DET") && matches!(entry, nucleide_serpent_io::Entry::Matrix(_))
        })
        .map(|(name, _)| name.clone())
        .collect();
    to_js(&SerpentDetSummary {
        variable_count: table.len(),
        detectors,
        variables: serpent_variables(&table),
    })
}

// ---------------------------------------------------------------------------
// FLUKA USRBIN
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct UsrbinTallyJson {
    name: String,
    particle: String,
    coord_sys: String,
    dims: [usize; 3],
    x_bounds: Vec<f64>,
    y_bounds: Vec<f64>,
    z_bounds: Vec<f64>,
    part_data: Vec<f64>,
    error_data: Vec<f64>,
}

#[derive(Serialize)]
struct UsrbinSummary {
    tally_count: usize,
    tallies: Vec<UsrbinTallyJson>,
}

/// Parse FLUKA `.lis` USRBIN tallies into a JSON summary.
#[wasm_bindgen(js_name = parseUsrbin)]
pub fn parse_usrbin(text: &str) -> Result<JsValue, JsValue> {
    let tallies = nucleide_fluka_io::usrbin::parse_usrbin(text).map_err(js_err)?;
    to_js(&UsrbinSummary {
        tally_count: tallies.len(),
        tallies: tallies
            .iter()
            .map(|t| UsrbinTallyJson {
                name: t.name.clone(),
                particle: t.particle.clone(),
                coord_sys: t.coord_sys.to_string(),
                dims: t.dims(),
                x_bounds: t.x_bounds.clone(),
                y_bounds: t.y_bounds.clone(),
                z_bounds: t.z_bounds.clone(),
                part_data: t.part_data.clone(),
                error_data: t.error_data.clone(),
            })
            .collect(),
    })
}

// ---------------------------------------------------------------------------
// ORIGEN tapes (TAPE5 input echo / TAPE6 inventory / TAPE9 decay constants)
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct OrigenTape5StepJson {
    flux: f64,
    days: f64,
}

#[derive(Serialize)]
struct OrigenGramsJson {
    nuclide: String,
    grams: f64,
}

#[derive(Serialize)]
struct OrigenTape5MaterialJson {
    name: String,
    entries: Vec<OrigenGramsJson>,
}

#[derive(Serialize)]
struct OrigenTape5Summary {
    titles: Vec<String>,
    steps: Vec<OrigenTape5StepJson>,
    materials: Vec<OrigenTape5MaterialJson>,
}

/// Parse an ORIGEN TAPE5 input echo into a JSON summary.
#[wasm_bindgen(js_name = parseOrigenTape5)]
pub fn parse_origen_tape5(text: &str) -> Result<JsValue, JsValue> {
    let tape = nucleide_origen_io::Tape5::parse(text).map_err(js_err)?;
    to_js(&OrigenTape5Summary {
        titles: tape.titles,
        steps: tape
            .irradiation_steps
            .iter()
            .map(|s| OrigenTape5StepJson {
                flux: s.flux,
                days: s.days,
            })
            .collect(),
        materials: tape
            .materials
            .iter()
            .map(|m| OrigenTape5MaterialJson {
                name: m.name.clone(),
                entries: m
                    .grams
                    .iter()
                    .map(|(nuclide, grams)| OrigenGramsJson {
                        nuclide: nuclide.clone(),
                        grams: *grams,
                    })
                    .collect(),
            })
            .collect(),
    })
}

#[derive(Serialize)]
struct OrigenTape6RecordJson {
    nuclide: String,
    grams: f64,
    activity_bq: f64,
}

#[derive(Serialize)]
struct OrigenTape6Summary {
    total_activity_bq: f64,
    records: Vec<OrigenTape6RecordJson>,
}

/// Parse an ORIGEN TAPE6 output inventory into a JSON summary.
#[wasm_bindgen(js_name = parseOrigenTape6)]
pub fn parse_origen_tape6(text: &str) -> Result<JsValue, JsValue> {
    let tape = nucleide_origen_io::Tape6::parse(text).map_err(js_err)?;
    to_js(&OrigenTape6Summary {
        total_activity_bq: tape.total_activity(),
        records: tape
            .records
            .iter()
            .map(|r| OrigenTape6RecordJson {
                nuclide: r.nuclide.clone(),
                grams: r.grams,
                activity_bq: r.activity_bq,
            })
            .collect(),
    })
}

#[derive(Serialize)]
struct OrigenTape9EntryJson {
    nuclide: String,
    decay_const: f64,
}

#[derive(Serialize)]
struct OrigenTape9Summary {
    entries: Vec<OrigenTape9EntryJson>,
}

/// Parse an ORIGEN TAPE9 decay-constant table into a JSON summary.
#[wasm_bindgen(js_name = parseOrigenTape9)]
pub fn parse_origen_tape9(text: &str) -> Result<JsValue, JsValue> {
    let entries = nucleide_origen_io::Tape9Entry::parse(text).map_err(js_err)?;
    to_js(&OrigenTape9Summary {
        entries: entries
            .iter()
            .map(|e| OrigenTape9EntryJson {
                nuclide: e.nuclide.clone(),
                decay_const: e.decay_const,
            })
            .collect(),
    })
}

// ---------------------------------------------------------------------------
// Deterministic transport (CCCC ISOTXS)
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct IsotxsNuclideJson {
    label: String,
    zaid: String,
    groups: usize,
    total_xs: Vec<f64>,
}

#[derive(Serialize)]
struct IsotxsSummary {
    nuclides: Vec<IsotxsNuclideJson>,
    groups: usize,
}

/// Parse an ISOTXS multigroup library into a JSON summary.
#[wasm_bindgen(js_name = parseIsotxs)]
pub fn parse_isotxs(text: &str) -> Result<JsValue, JsValue> {
    let lib = nucleide_cccc_io::IsotxsLib::parse(text).map_err(js_err)?;
    to_js(&IsotxsSummary {
        groups: lib.nuclides.first().map(|n| n.groups).unwrap_or(0),
        nuclides: lib
            .nuclides
            .iter()
            .map(|n| IsotxsNuclideJson {
                label: n.label.clone(),
                zaid: n.zaid.clone(),
                groups: n.groups,
                total_xs: n.total_xs.clone(),
            })
            .collect(),
    })
}

// ---------------------------------------------------------------------------
// MCNP full-deck problem
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct DeckCellJson {
    num: u32,
    mat: u32,
    #[serde(skip_serializing_if = "Option::is_none")]
    dens: Option<f64>,
    geom: String,
    params: Vec<String>,
}

#[derive(Serialize)]
struct DeckSurfJson {
    num: u32,
    reflecting: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    transform: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    periodic: Option<u32>,
    kind: String,
    coeffs: Vec<f64>,
}

#[derive(Serialize)]
struct ModeJson {
    particles: Vec<String>,
}

#[derive(Serialize)]
struct TransformJson {
    number: u32,
    displacement: [f64; 3],
    rotation: Vec<f64>,
    #[serde(rename = "inDegrees")]
    in_degrees: bool,
    #[serde(rename = "mainToAux")]
    main_to_aux: bool,
    hidden: bool,
}

#[derive(Serialize)]
struct UniverseJson {
    number: u32,
    cells: Vec<u32>,
    #[serde(rename = "notTruncated")]
    not_truncated: Vec<u32>,
}

#[derive(Serialize)]
struct LatticeJson {
    cell: u32,
    lattice: u8,
}

#[derive(Serialize)]
struct FillJson {
    cell: u32,
    kind: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    universe: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(rename = "minIndex")]
    min_index: Option<[i32; 3]>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(rename = "maxIndex")]
    max_index: Option<[i32; 3]>,
    #[serde(skip_serializing_if = "Option::is_none")]
    universes: Option<Vec<Option<u32>>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    transform: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(rename = "hiddenTransform")]
    hidden_transform: Option<Vec<f64>>,
    #[serde(rename = "inDegrees")]
    in_degrees: bool,
}

#[derive(Serialize)]
struct ImportanceJson {
    cell: u32,
    particle: String,
    value: f64,
}

#[derive(Serialize)]
struct VolumeJson {
    cell: u32,
    volume: f64,
}

#[derive(Serialize)]
struct TallyJson {
    number: u32,
    #[serde(rename = "type")]
    tally_type: u8,
    particles: Vec<String>,
    entries: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    fm: Option<Vec<String>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(rename = "eBins")]
    e_bins: Option<Vec<String>>,
}

/// A parsed MCNP input deck with format-preserving write-back.
///
/// Thin facade over [`nucleide_mcnp_io::problem::DeckProblem`]: text in via
/// [`parse_deck`](nucleide_mcnp_io::problem::parse_deck), text out via
/// [`dumps`](nucleide_mcnp_io::problem::DeckProblem::dumps). There is no
/// file-reading entry point (no filesystem in the browser).
#[wasm_bindgen]
pub struct WasmDeckProblem {
    inner: nucleide_mcnp_io::problem::DeckProblem,
}

#[wasm_bindgen]
impl WasmDeckProblem {
    /// Parse a deck from its text.
    #[wasm_bindgen(js_name = fromText)]
    pub fn from_text(text: &str) -> Result<WasmDeckProblem, JsValue> {
        nucleide_mcnp_io::problem::parse_deck(text)
            .map(|inner| WasmDeckProblem { inner })
            .map_err(js_err)
    }

    /// Serialize back to MCNP input text (byte-identical when unedited).
    pub fn dumps(&self) -> String {
        self.inner.dumps()
    }

    /// Message (first) line.
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.message.clone()
    }

    /// Title card (second line).
    #[wasm_bindgen(getter)]
    pub fn title(&self) -> String {
        self.inner.title.clone()
    }

    /// Cell cards (`dens` is absent for void cells).
    pub fn cells(&self) -> Result<JsValue, JsValue> {
        let out: Vec<DeckCellJson> = self
            .inner
            .cells
            .iter()
            .map(|c| DeckCellJson {
                num: c.num,
                mat: c.mat,
                dens: c.dens,
                geom: c.geom.render(),
                params: c.params.clone(),
            })
            .collect();
        to_js(&out)
    }

    /// Surface cards (`transform`/`periodic` are absent when not present).
    pub fn surfs(&self) -> Result<JsValue, JsValue> {
        let out: Vec<DeckSurfJson> = self
            .inner
            .surfs
            .iter()
            .map(|s| DeckSurfJson {
                num: s.num,
                reflecting: s.reflecting,
                transform: s.transform,
                periodic: s.periodic,
                kind: s.kind.keyword().to_string(),
                coeffs: s.coeffs.clone(),
            })
            .collect();
        to_js(&out)
    }

    /// Material numbers in file order.
    #[wasm_bindgen(js_name = materialNumbers)]
    pub fn material_numbers(&self) -> Vec<u32> {
        self.inner.materials.iter().map(|m| m.number).collect()
    }

    /// Data-card names in file order (`MODE`, `M1`, `KCODE`, ...).
    #[wasm_bindgen(js_name = dataNames)]
    pub fn data_names(&self) -> Vec<String> {
        self.inner.data.iter().map(|d| d.name.clone()).collect()
    }

    /// Typed `MODE` card.
    pub fn mode(&self) -> Result<JsValue, JsValue> {
        let mode = self.inner.mode().map_err(js_err)?;
        to_js(&ModeJson {
            particles: mode.particles,
        })
    }

    /// Typed `TRn` cards in file order.
    pub fn transforms(&self) -> Result<JsValue, JsValue> {
        let out: Vec<TransformJson> = self
            .inner
            .transforms()
            .map_err(js_err)?
            .into_iter()
            .map(|t| TransformJson {
                number: t.number,
                displacement: t.displacement,
                rotation: t.rotation,
                in_degrees: t.is_in_degrees,
                main_to_aux: t.is_main_to_aux,
                hidden: t.hidden,
            })
            .collect();
        to_js(&out)
    }

    /// Auto-created universes (from cell/data `U` usage, including 0).
    pub fn universes(&self) -> Result<JsValue, JsValue> {
        let out: Vec<UniverseJson> = self
            .inner
            .universes()
            .map_err(js_err)?
            .into_iter()
            .map(|u| UniverseJson {
                number: u.number,
                cells: u.cells,
                not_truncated: u.not_truncated,
            })
            .collect();
        to_js(&out)
    }

    /// Cell `LAT` assignments in file order.
    pub fn lattices(&self) -> Result<JsValue, JsValue> {
        let out: Vec<LatticeJson> = self
            .inner
            .lattices()
            .map_err(js_err)?
            .into_iter()
            .map(|l| LatticeJson {
                cell: l.cell,
                lattice: l.lattice,
            })
            .collect();
        to_js(&out)
    }

    /// Cell `FILL` assignments in file order (`kind` is `single` or
    /// `matrix`; matrix empties are `null`).
    pub fn fills(&self) -> Result<JsValue, JsValue> {
        use nucleide_mcnp_io::semantic::{FillTarget, FillTransform};
        let out: Vec<FillJson> = self
            .inner
            .fills()
            .map_err(js_err)?
            .into_iter()
            .map(|f| {
                let mut json = FillJson {
                    cell: f.cell,
                    kind: String::new(),
                    universe: None,
                    min_index: None,
                    max_index: None,
                    universes: None,
                    transform: None,
                    hidden_transform: None,
                    in_degrees: f.in_degrees,
                };
                match &f.target {
                    FillTarget::Single(u) => {
                        json.kind = "single".to_string();
                        json.universe = Some(*u);
                    }
                    FillTarget::Matrix {
                        min_index,
                        max_index,
                        universes,
                    } => {
                        json.kind = "matrix".to_string();
                        json.min_index = Some(*min_index);
                        json.max_index = Some(*max_index);
                        json.universes = Some(universes.clone());
                    }
                }
                match &f.transform {
                    None => {}
                    Some(FillTransform::Reference(n)) => {
                        json.transform = Some(*n);
                    }
                    Some(FillTransform::Hidden(t)) => {
                        let mut coords: Vec<f64> = t.displacement.to_vec();
                        coords.extend(t.rotation.iter().copied());
                        json.hidden_transform = Some(coords);
                    }
                }
                json
            })
            .collect();
        to_js(&out)
    }

    /// Cell importance entries in file order.
    pub fn importances(&self) -> Result<JsValue, JsValue> {
        let out: Vec<ImportanceJson> = self
            .inner
            .importances()
            .map_err(js_err)?
            .into_iter()
            .map(|v| ImportanceJson {
                cell: v.cell,
                particle: v.particle,
                value: v.value,
            })
            .collect();
        to_js(&out)
    }

    /// Manual cell volumes in file order.
    pub fn volumes(&self) -> Result<JsValue, JsValue> {
        let out: Vec<VolumeJson> = self
            .inner
            .volumes()
            .map_err(js_err)?
            .into_iter()
            .map(|v| VolumeJson {
                cell: v.cell,
                volume: v.volume,
            })
            .collect();
        to_js(&out)
    }

    /// Typed tallies (`Fn` with grouped `FMn`/`En`) in number order.
    pub fn tallies(&self) -> Result<JsValue, JsValue> {
        let out: Vec<TallyJson> = self
            .inner
            .tallies()
            .map_err(js_err)?
            .into_iter()
            .map(|t| TallyJson {
                number: t.number,
                tally_type: t.tally_type,
                particles: t.particles,
                entries: t.entries,
                fm: t.fm,
                e_bins: t.e_bins,
            })
            .collect();
        to_js(&out)
    }

    /// Cell inventory: `cell -> material` map in file order.
    #[wasm_bindgen(js_name = cellInventory)]
    pub fn cell_inventory(&self) -> Result<JsValue, JsValue> {
        // Keys as strings: serde-wasm-bindgen maps need string keys to
        // become plain JS objects (JS object keys are strings anyway).
        let inv: BTreeMap<String, u32> = nucleide_mcnp_io::problem::cell_inventory(&self.inner)
            .into_iter()
            .map(|(cell, mat)| (cell.to_string(), mat))
            .collect();
        to_js(&inv)
    }

    /// Validate every L3 semantic rule (throws on error).
    pub fn validate(&self) -> Result<(), JsValue> {
        self.inner.validate().map_err(js_err)
    }

    /// Non-fatal validation notes (particle/mode mismatches, not errors).
    #[wasm_bindgen(js_name = validationNotes)]
    pub fn validation_notes(&self) -> Vec<String> {
        self.inner.validation_notes()
    }

    /// Set a cell's density, re-rendering that card canonically.
    #[wasm_bindgen(js_name = setCellDensity)]
    pub fn set_cell_density(&mut self, cell: u32, dens: f64) -> Result<(), JsValue> {
        self.inner.set_cell_density(cell, dens).map_err(js_err)
    }

    /// Set a cell's material number, re-rendering that card canonically.
    #[wasm_bindgen(js_name = setCellMaterial)]
    pub fn set_cell_material(&mut self, cell: u32, mat: u32) -> Result<(), JsValue> {
        self.inner.set_cell_material(cell, mat).map_err(js_err)
    }

    /// Set the `MODE` card particles, re-rendering that card canonically
    /// (appending one when absent).
    #[wasm_bindgen(js_name = setMode)]
    pub fn set_mode(&mut self, particles: Vec<String>) -> Result<(), JsValue> {
        self.inner.set_mode(particles).map_err(js_err)
    }

    /// Set a cell's universe (`notTruncated` writes `U=-n`).
    #[wasm_bindgen(js_name = setCellUniverse)]
    pub fn set_cell_universe(
        &mut self,
        cell: u32,
        universe: u32,
        not_truncated: Option<bool>,
    ) -> Result<(), JsValue> {
        self.inner
            .set_cell_universe(cell, universe, not_truncated.unwrap_or(false))
            .map_err(js_err)
    }

    /// Set (`1`/`2`) or clear (`undefined`) a cell's lattice.
    #[wasm_bindgen(js_name = setCellLattice)]
    pub fn set_cell_lattice(&mut self, cell: u32, lattice: Option<u8>) -> Result<(), JsValue> {
        self.inner.set_cell_lattice(cell, lattice).map_err(js_err)
    }

    /// Set a cell's fill to a single universe.
    #[wasm_bindgen(js_name = setCellFill)]
    pub fn set_cell_fill(&mut self, cell: u32, universe: u32) -> Result<(), JsValue> {
        self.inner.set_cell_fill(cell, universe).map_err(js_err)
    }
}

// ---------------------------------------------------------------------------
// Decay inventory
// ---------------------------------------------------------------------------

fn inventory_system(
    chain: &nucleide_depletion::Chain,
    rates: &nucleide_depletion::ReactionRates,
) -> Result<nucleide_depletion::DepletionSystem, JsValue> {
    nucleide_depletion::DepletionSystem::build(chain.clone(), rates).map_err(js_err)
}

fn inventory_atoms(
    chain: &nucleide_depletion::Chain,
    atoms: &BTreeMap<String, f64>,
) -> Result<Vec<f64>, JsValue> {
    let mut vec = vec![0.0; chain.len()];
    for (name, value) in atoms {
        let idx = chain
            .index_of(name)
            .ok_or_else(|| js_err(format!("unknown nuclide `{name}` for this chain")))?;
        vec[idx] = *value;
    }
    Ok(vec)
}

/// A unit-aware decay inventory over a depletion chain.
///
/// Thin facade over [`nucleide_depletion::DecayInventory`]: quantities convert
/// through [`QuantityUnit`](nucleide_depletion::QuantityUnit) on entry and on
/// exit, and [`decay`](WasmInventory::decay) routes through the core series
/// driver (predictor over one step) so reaction rates and the solver
/// `method` stay honored exactly like the Python `Inventory.decay`.
#[wasm_bindgen]
pub struct WasmInventory {
    chain: std::sync::Arc<nucleide_depletion::Chain>,
    inner: nucleide_depletion::DecayInventory,
}

#[wasm_bindgen]
impl WasmInventory {
    /// Build from quantities in `units` (default `"atoms"`; `Bq`/`Ci`
    /// activity, `g`/`kg` mass, `mol`, ... — see `QuantityUnit`).
    #[wasm_bindgen(constructor)]
    pub fn new(
        chain: &WasmChain,
        comp: JsValue,
        units: Option<String>,
    ) -> Result<WasmInventory, JsValue> {
        let comp: BTreeMap<String, f64> = serde_wasm_bindgen::from_value(comp).map_err(js_err)?;
        let unit = units
            .as_deref()
            .unwrap_or("atoms")
            .parse::<nucleide_depletion::QuantityUnit>()
            .map_err(js_err)?;
        let sys = inventory_system(&chain.inner, &nucleide_depletion::ReactionRates::new())?;
        let inner =
            nucleide_depletion::DecayInventory::from_units(&comp, unit, &sys).map_err(js_err)?;
        Ok(WasmInventory {
            chain: chain.inner.clone(),
            inner,
        })
    }

    /// Parse [`to_csv`](WasmInventory::to_csv) output back into an inventory
    /// over `chain` (names outside the chain are an error).
    #[wasm_bindgen(js_name = fromCsv)]
    pub fn from_csv(chain: &WasmChain, text: &str) -> Result<WasmInventory, JsValue> {
        let inner = nucleide_depletion::DecayInventory::from_csv(text).map_err(js_err)?;
        for name in inner.atoms.keys() {
            if chain.inner.index_of(name).is_none() {
                return Err(js_err(format!("unknown nuclide `{name}` for this chain")));
            }
        }
        Ok(WasmInventory {
            chain: chain.inner.clone(),
            inner,
        })
    }

    /// Atom counts by nuclide name.
    pub fn numbers(&self) -> Result<JsValue, JsValue> {
        to_js(&self.inner.numbers())
    }

    /// Decay over `dt` in `time_unit` (default `"s"`); optional one-group
    /// `rates` (`"Name:reaction"` keys), CRAM `order` (default 48), and
    /// solver `method` (`"cram16"`, `"cram48"`, `"bateman"`, `"bateman_hp"`,
    /// default `"cram48"` — an explicitly non-default `method` overrides
    /// `order`). Unlike the decay-only core, this honors `rates`; a Bateman
    /// `method` with live rates falls back to CRAM-48.
    pub fn decay(
        &self,
        dt: f64,
        time_unit: Option<String>,
        rates: JsValue,
        order: Option<u8>,
        method: Option<String>,
    ) -> Result<WasmInventory, JsValue> {
        let method = resolve_method(order.unwrap_or(48), method.as_deref())?;
        let unit = nucleide_depletion::time_unit_from_str(time_unit.as_deref().unwrap_or("s"))
            .map_err(js_err)?;
        let rate_map: BTreeMap<String, f64> = if rates.is_undefined() || rates.is_null() {
            BTreeMap::new()
        } else {
            serde_wasm_bindgen::from_value(rates).map_err(js_err)?
        };
        let reaction_rates = parse_reaction_rates(&self.chain, &rate_map)?;
        let template = inventory_system(&self.chain, &reaction_rates)?;
        let seconds = dt * unit.as_seconds();
        let steps = vec![nucleide_depletion::Step::new(
            seconds,
            reaction_rates.clone(),
        )];
        let series = nucleide_depletion::integrate_with_method(
            &template,
            &inventory_atoms(&self.chain, &self.inner.atoms)?,
            &steps,
            nucleide_depletion::Integrator::Predictor,
            method,
        )
        .map_err(js_err)?;
        let names: Vec<String> = self.chain.nuclides.iter().map(|n| n.name.clone()).collect();
        let atoms = names
            .iter()
            .zip(series.atoms.last().cloned().unwrap_or_default())
            .map(|(n, v)| (n.clone(), v))
            .collect();
        Ok(WasmInventory {
            chain: self.chain.clone(),
            inner: nucleide_depletion::DecayInventory { atoms },
        })
    }

    /// Activity per nuclide in `units`.
    pub fn activities(&self, units: &str) -> Result<JsValue, JsValue> {
        let unit = units
            .parse::<nucleide_depletion::QuantityUnit>()
            .map_err(js_err)?;
        let sys = inventory_system(&self.chain, &nucleide_depletion::ReactionRates::new())?;
        to_js(&self.inner.activities(&sys, unit).map_err(js_err)?)
    }

    /// Mass per nuclide in `units`.
    pub fn masses(&self, units: &str) -> Result<JsValue, JsValue> {
        let unit = units
            .parse::<nucleide_depletion::QuantityUnit>()
            .map_err(js_err)?;
        to_js(&self.inner.masses(unit).map_err(js_err)?)
    }

    /// Moles per nuclide in `units`.
    pub fn moles(&self, units: &str) -> Result<JsValue, JsValue> {
        let unit = units
            .parse::<nucleide_depletion::QuantityUnit>()
            .map_err(js_err)?;
        to_js(&self.inner.moles(unit).map_err(js_err)?)
    }

    /// Fraction of total activity per nuclide.
    #[wasm_bindgen(js_name = activityFractions)]
    pub fn activity_fractions(&self) -> Result<JsValue, JsValue> {
        let sys = inventory_system(&self.chain, &nucleide_depletion::ReactionRates::new())?;
        to_js(&self.inner.activity_fractions(&sys).map_err(js_err)?)
    }

    /// Fraction of total mass per nuclide.
    #[wasm_bindgen(js_name = massFractions)]
    pub fn mass_fractions(&self) -> Result<JsValue, JsValue> {
        to_js(&self.inner.mass_fractions().map_err(js_err)?)
    }

    /// Fraction of total atoms per nuclide (mole fractions).
    #[wasm_bindgen(js_name = moleFractions)]
    pub fn mole_fractions(&self) -> Result<JsValue, JsValue> {
        to_js(&self.inner.mole_fractions())
    }

    /// Human-readable half-lives (`"3.2 d"`, `"stable"`, `"unknown"`).
    #[wasm_bindgen(js_name = halfLivesReadable)]
    pub fn half_lives_readable(&self) -> Result<JsValue, JsValue> {
        to_js(&self.inner.half_lives_readable())
    }

    /// Add two inventories (atom counts sum).
    pub fn add(&self, other: &WasmInventory) -> WasmInventory {
        WasmInventory {
            chain: self.chain.clone(),
            inner: self.inner.add(&other.inner),
        }
    }

    /// Subtract (clamped at zero).
    pub fn sub(&self, other: &WasmInventory) -> WasmInventory {
        WasmInventory {
            chain: self.chain.clone(),
            inner: self.inner.sub(&other.inner),
        }
    }

    /// Scale all counts by `s`.
    pub fn mul(&self, s: f64) -> WasmInventory {
        WasmInventory {
            chain: self.chain.clone(),
            inner: self.inner.mul(s),
        }
    }

    /// Divide all counts by `s`.
    pub fn div(&self, s: f64) -> WasmInventory {
        WasmInventory {
            chain: self.chain.clone(),
            inner: self.inner.div(s),
        }
    }

    /// Serialize as `nuclide,quantity,unit` CSV rows (atom counts).
    #[wasm_bindgen(js_name = toCsv)]
    pub fn to_csv(&self) -> String {
        self.inner.to_csv()
    }
}

/// Time-integrated decays per nuclide over one step, keyed by name.
///
/// Diagonal Bateman integral over `[0, dt]` (stable nuclides report `0.0`);
/// optional one-group `rates` select the depletion system.
#[wasm_bindgen(js_name = cumulativeDecays)]
pub fn cumulative_decays(
    chain: &WasmChain,
    n0: JsValue,
    dt: f64,
    rates: JsValue,
) -> Result<JsValue, JsValue> {
    let n0: BTreeMap<String, f64> = serde_wasm_bindgen::from_value(n0).map_err(js_err)?;
    let rate_map: BTreeMap<String, f64> = if rates.is_undefined() || rates.is_null() {
        BTreeMap::new()
    } else {
        serde_wasm_bindgen::from_value(rates).map_err(js_err)?
    };
    let sys = inventory_system(
        &chain.inner,
        &parse_reaction_rates(&chain.inner, &rate_map)?,
    )?;
    let vec = inventory_atoms(&chain.inner, &n0)?;
    let out = nucleide_depletion::cumulative_decays(&sys, &vec, dt).map_err(js_err)?;
    let keyed: BTreeMap<String, f64> = chain
        .inner
        .nuclides
        .iter()
        .zip(out)
        .map(|(nuc, v)| (nuc.name.clone(), v))
        .collect();
    to_js(&keyed)
}

/// `(child, branching_ratio, decay_mode)` triples for a chain nuclide
/// (empty when the name is unknown).
pub fn progeny(chain: &WasmChain, name: &str) -> Result<JsValue, JsValue> {
    to_js(&nucleide_depletion::progeny(&chain.inner, name))
}

/// Branching fraction from parent to child (`undefined` when absent).
#[wasm_bindgen(js_name = branchingFraction)]
pub fn branching_fraction(
    chain: &WasmChain,
    parent: &str,
    child: &str,
) -> Result<Option<f64>, JsValue> {
    Ok(nucleide_depletion::branching_fraction(
        &chain.inner,
        parent,
        child,
    ))
}

/// Decay-mode label from parent to child (`undefined` when absent).
#[wasm_bindgen(js_name = decayMode)]
pub fn decay_mode(chain: &WasmChain, parent: &str, child: &str) -> Result<Option<String>, JsValue> {
    Ok(nucleide_depletion::decay_mode(&chain.inner, parent, child))
}

/// `(parent, child, branching_ratio, decay_mode)` edges of a chain.
#[wasm_bindgen(js_name = chainEdges)]
pub fn chain_edges(chain: &WasmChain) -> Result<JsValue, JsValue> {
    to_js(&nucleide_depletion::chain_edges(&chain.inner))
}

// ---------------------------------------------------------------------------
// Emission (five-dialect cards + mass drift)
// ---------------------------------------------------------------------------

/// Optional per-dialect overrides for the `emit*` functions (`density` stays
/// a top-level argument, like the Python `emit_cards` signature).
#[derive(Deserialize, Default)]
struct EmitOptsJson {
    #[serde(default, rename = "mcnpNumber", alias = "mcnp_number")]
    mcnp_number: Option<u32>,
    #[serde(default, rename = "xsSuffix", alias = "xs_suffix")]
    xs_suffix: Option<String>,
    #[serde(default, rename = "serpentLib", alias = "serpent_lib")]
    serpent_lib: Option<String>,
    #[serde(default, rename = "flukaFid", alias = "fluka_fid")]
    fluka_fid: Option<u32>,
    #[serde(default, rename = "partisnZone", alias = "partisn_zone")]
    partisn_zone: Option<u32>,
}

#[derive(Serialize)]
struct DroppedJson {
    nuclide: String,
    mass: f64,
    reason: String,
}

#[derive(Serialize)]
struct DriftRowJson {
    code: String,
    #[serde(rename = "massIn")]
    mass_in: f64,
    #[serde(rename = "massOut")]
    mass_out: f64,
    #[serde(rename = "relDrift")]
    rel_drift: f64,
    dropped: Vec<DroppedJson>,
    reparsed: bool,
}

fn comp_to_emit_material(comp: JsValue) -> Result<nucleide_material::Material, JsValue> {
    let map: BTreeMap<String, f64> = serde_wasm_bindgen::from_value(comp).map_err(js_err)?;
    let mut mat = nucleide_material::Material::new();
    for (name, grams) in &map {
        let id = name
            .parse::<NuclideId>()
            .map_err(|e| js_err(format!("`{name}`: {e}")))?;
        mat.add_nuclide(id, *grams);
    }
    Ok(mat)
}

fn emit_options(name: &str, opts: JsValue) -> Result<nucleide_emit::EmitOptions, JsValue> {
    let parsed: EmitOptsJson = if opts.is_undefined() || opts.is_null() {
        EmitOptsJson::default()
    } else {
        serde_wasm_bindgen::from_value(opts).map_err(js_err)?
    };
    let mut out = nucleide_emit::EmitOptions::new(name);
    if let Some(v) = parsed.mcnp_number {
        out.mcnp_number = v;
    }
    if let Some(v) = parsed.xs_suffix {
        out.xs_suffix = v;
    }
    if let Some(v) = parsed.serpent_lib {
        out.serpent_lib = v;
    }
    if let Some(v) = parsed.fluka_fid {
        out.fluka_fid = v;
    }
    if let Some(v) = parsed.partisn_zone {
        out.partisn_zone = v;
    }
    Ok(out)
}

fn emit_cards_json(emitted: &[nucleide_emit::Emitted]) -> BTreeMap<String, String> {
    emitted
        .iter()
        .map(|e| (e.code.to_string(), e.text.clone()))
        .collect()
}

fn drift_rows_json(table: &nucleide_emit::DriftTable) -> Vec<DriftRowJson> {
    table
        .rows
        .iter()
        .map(|r| DriftRowJson {
            code: r.code.to_string(),
            mass_in: r.mass_in,
            mass_out: r.mass_out,
            rel_drift: r.rel_drift,
            dropped: r
                .dropped
                .iter()
                .map(|d| DroppedJson {
                    nuclide: d.id.to_name(),
                    mass: d.mass,
                    reason: d.reason.clone(),
                })
                .collect(),
            reparsed: r.reparsed,
        })
        .collect()
}

/// Emit one composition through all five code dialects (MCNP, Serpent, FLUKA,
/// ALARA, PARTISN). Returns `{code: card_text}`.
///
/// `comp` maps GNDS nuclide names to grams; `density` is the mass density
/// [g/cm³] for dialects that need one (Serpent/FLUKA/PARTISN throw
/// `MissingDensity` without it); `opts` optionally overrides `mcnpNumber`
/// (default 1), `xsSuffix` (default `"80c"`), `serpentLib` (default `"03c"`),
/// `flukaFid` (default 1), and `partisnZone` (default 1).
#[wasm_bindgen(js_name = emitCards)]
pub fn emit_cards(
    comp: JsValue,
    name: &str,
    density: Option<f64>,
    opts: JsValue,
) -> Result<JsValue, JsValue> {
    let mut mat = comp_to_emit_material(comp)?;
    mat.set_density(density);
    let options = emit_options(name, opts)?;
    let (emitted, _) = nucleide_emit::emit_drift(&mat, &options).map_err(js_err)?;
    to_js(&emit_cards_json(&emitted))
}

/// Mass-drift report for one composition across all five code dialects.
///
/// Same inputs as [`emit_cards`]; returns
/// `[{code, massIn, massOut, relDrift, dropped: [{nuclide, mass, reason}], reparsed}]`.
#[wasm_bindgen(js_name = emitDriftTable)]
pub fn emit_drift_table(
    comp: JsValue,
    name: &str,
    density: Option<f64>,
    opts: JsValue,
) -> Result<JsValue, JsValue> {
    let mut mat = comp_to_emit_material(comp)?;
    mat.set_density(density);
    let options = emit_options(name, opts)?;
    let (_, table) = nucleide_emit::emit_drift(&mat, &options).map_err(js_err)?;
    to_js(&drift_rows_json(&table))
}

/// Emit one ARMI-keyed composition through all five code dialects (MCNP,
/// Serpent, FLUKA, ALARA, PARTISN). Returns `{code: card_text}`.
///
/// `comp` maps ARMI nuclide keys (`nU235`, `92235`, `U-2355`, ...) to grams;
/// keys resolve via `nucleide_emit::armi::from_armi_mass_fracs`, so elemental
/// keys, bare `AM242`, and negative/non-finite masses throw. `density` is the
/// hot mass density [g/cm³] for dialects that need one.
#[wasm_bindgen(js_name = emitArmiCards)]
pub fn emit_armi_cards(
    comp: JsValue,
    name: &str,
    density: Option<f64>,
    opts: JsValue,
) -> Result<JsValue, JsValue> {
    let map: BTreeMap<String, f64> = serde_wasm_bindgen::from_value(comp).map_err(js_err)?;
    let mat = nucleide_emit::armi::from_armi_mass_fracs(map, density).map_err(js_err)?;
    let options = emit_options(name, opts)?;
    let (emitted, _) = nucleide_emit::emit_drift(&mat, &options).map_err(js_err)?;
    to_js(&emit_cards_json(&emitted))
}

/// Mass-drift report for one ARMI-keyed composition across all five code
/// dialects.
///
/// Same inputs as [`emit_armi_cards`]; returns
/// `[{code, massIn, massOut, relDrift, dropped: [{nuclide, mass, reason}], reparsed}]`.
#[wasm_bindgen(js_name = emitArmiDriftTable)]
pub fn emit_armi_drift_table(
    comp: JsValue,
    name: &str,
    density: Option<f64>,
    opts: JsValue,
) -> Result<JsValue, JsValue> {
    let map: BTreeMap<String, f64> = serde_wasm_bindgen::from_value(comp).map_err(js_err)?;
    let mat = nucleide_emit::armi::from_armi_mass_fracs(map, density).map_err(js_err)?;
    let options = emit_options(name, opts)?;
    let (_, table) = nucleide_emit::emit_drift(&mat, &options).map_err(js_err)?;
    to_js(&drift_rows_json(&table))
}

// ---------------------------------------------------------------------------
// Dose factors (screening-level only)
// ---------------------------------------------------------------------------

fn parse_dose_pathway(s: &str) -> Result<nucleide_nuclei::data::DosePathway, JsValue> {
    nucleide_nuclei::data::DosePathway::parse(s).ok_or_else(|| {
        js_err(format!(
            "unknown dose pathway `{s}` (supported: air, soil, ingest, inhale)"
        ))
    })
}

fn parse_dose_source(s: &str) -> Result<nucleide_nuclei::data::DoseSource, JsValue> {
    nucleide_nuclei::data::DoseSource::parse(s).ok_or_else(|| {
        js_err(format!(
            "unknown dose source `{s}` (supported: EPA, DOE, GENII)"
        ))
    })
}

/// Raw dose factor for a nuclide name, pathway, and source.
///
/// `pathway` is one of `air`/`soil`/`ingest`/`inhale` (`ext_air`/`ext_soil`
/// aliases accepted); `source` is one of `EPA`/`DOE`/`GENII` (default `EPA`).
/// Returns `undefined` when the nuclide has no row; GENII/DOE air resolve to
/// `-1.0` (PyNE missing-air sentinel). Screening-level only — not for safety
/// decisions.
#[wasm_bindgen(js_name = doseFactor)]
pub fn dose_factor(
    name: &str,
    pathway: &str,
    source: Option<String>,
) -> Result<Option<f64>, JsValue> {
    let id = name
        .parse::<NuclideId>()
        .map_err(|e| js_err(format!("`{name}`: {e}")))?;
    let p = parse_dose_pathway(pathway)?;
    let s = parse_dose_source(source.as_deref().unwrap_or("EPA"))?;
    Ok(nucleide_nuclei::data::dose_factor(id.nucid(), p, s))
}

/// Total dose per gram of a composition (`{nuclide name: grams}`).
///
/// Thin wrapper over `Material::total_dose_per_g` (AME2020 masses, ENDF/B-VIII.0
/// decay constants, HNF-5636/PyNE dose factors). `pathway` is one of
/// `air`/`soil`/`ingest`/`inhale`; `source` is `EPA`/`DOE`/`GENII` (default
/// `EPA`). Units follow the table: air `mrem/h per g per m^3`, soil
/// `mrem/h per g per m^2`, ingest/inhale `mrem per g`. Screening-level only —
/// not for safety decisions. Stable nuclides (known mass, no decay constant)
/// contribute 0 without a dose-factor lookup. Throws when a nuclide lacks
/// mass data, or a radioactive nuclide lacks dose data (including `-1`
/// GENII/DOE air sentinels).
#[wasm_bindgen(js_name = dosePerGram)]
pub fn dose_per_gram(comp: JsValue, pathway: &str, source: Option<String>) -> Result<f64, JsValue> {
    let mat = comp_to_emit_material(comp)?;
    let analytics = nucleide_material::Analytics {
        masses: &nucleide_material::Ame2020,
        decays: &nucleide_material::ChainDecays,
    };
    let p = parse_dose_pathway(pathway)?;
    let s = parse_dose_source(source.as_deref().unwrap_or("EPA"))?;
    mat.total_dose_per_g(&analytics, &nucleide_material::DoseFactors, p, s)
        .map_err(js_err)
}

// ---------------------------------------------------------------------------
// R2S snapshot (local port — no `nucleide-r2s` dependency)
// ---------------------------------------------------------------------------
//
// `nucleide-r2s` depends on `nucleide-depletion` with default features, which
// would re-enable Rayon through Cargo feature unification and break this
// `wasm32-unknown-unknown` build (the crate offers no `default-features = false`
// switch). The thin snapshot→deck adapter from `nucleide_r2s::snapshot` is
// therefore re-implemented here over `nucleide-alara-io` deck types plus the
// nuclei ARMI bridge — the same no-`r2s`-dep precedent as `r2s_from_deck`
// above. Key rules mirror the `nucleide-emit` ARMI-input rule (elemental keys
// rejected with "expand first", bare `AM242` rejected, explicit `AM242G`
// ground accepted).

/// Geometry string emitted for snapshot decks (volumes method still needs a
/// `geometry` block for parser consumers; the value carries no mesh axes).
const SNAPSHOT_GEOMETRY: &str = "rectangular";
/// Default schedule name synthesized when no `scheduleText` is given.
const SNAPSHOT_SCHEDULE: &str = "snap_schedule";
/// Default pulsing-history name synthesized when no `scheduleText` is given.
const SNAPSHOT_HISTORY: &str = "snap_once";

#[derive(Deserialize)]
#[allow(dead_code)] // `material`/`xs_type` are informational round-trip metadata, never emitted.
struct SnapshotZoneJson {
    id: String,
    #[serde(rename = "volumeCm3", alias = "volume_cm3")]
    volume_cm3: f64,
    #[serde(default, rename = "zbottomCm", alias = "zbottom_cm")]
    zbottom_cm: Option<f64>,
    #[serde(default, rename = "ztopCm", alias = "ztop_cm")]
    ztop_cm: Option<f64>,
    #[serde(default)]
    material: Option<String>,
    #[serde(default, rename = "xsType", alias = "xs_type")]
    xs_type: Option<String>,
    #[serde(default, rename = "temperatureC", alias = "temperature_C")]
    temperature_c: Option<f64>,
    composition: BTreeMap<String, f64>,
    #[serde(default)]
    flux: Option<String>,
}

#[derive(Deserialize)]
struct SnapshotFluxJson {
    name: String,
    file: String,
    scale: f64,
}

#[derive(Deserialize)]
struct SnapshotInputJson {
    zones: Vec<SnapshotZoneJson>,
    #[serde(rename = "fluxDefs", alias = "flux_defs")]
    flux_defs: Vec<SnapshotFluxJson>,
    #[serde(rename = "coolingS", alias = "cooling_s")]
    cooling_s: Vec<f64>,
    #[serde(default, rename = "scheduleText", alias = "schedule_text")]
    schedule_text: Option<String>,
    #[serde(default)]
    output: Option<String>,
}

#[derive(Serialize)]
struct SnapshotBundleJson {
    workflow: R2sSummary,
    deck: String,
    decks: Vec<String>,
}

/// Mixture name for a zone (`mix_<zone>`; the prefix keeps zones named
/// `void` distinct from the `void` mixture).
fn snap_mixture_name(zone: &str) -> String {
    format!("mix_{zone}")
}

/// Record a raw block so the canonical writer replays the typed views in deck
/// order (synthesized blocks use line 0, like `emit_decks` does).
fn snap_push_block(deck: &mut nucleide_alara_io::deck::AlaraDeck, kind: &str, body: Vec<String>) {
    deck.blocks.push(nucleide_alara_io::deck::RawBlock {
        kind: kind.to_string(),
        line: 0,
        body,
    });
}

/// Zone ids are opaque deck tokens: non-empty with no whitespace.
fn snap_check_zone_id(zone: &str) -> Result<(), String> {
    if zone.trim().is_empty() {
        return Err("snapshot zone id is empty".to_string());
    }
    if zone.chars().any(char::is_whitespace) {
        return Err(format!(
            "snapshot zone `{zone}` contains whitespace (deck tokens must not)"
        ));
    }
    Ok(())
}

/// Volumes must be finite and positive (ALARA interval volumes in cm³).
fn snap_check_volume(zone: &str, volume: f64) -> Result<(), String> {
    if !volume.is_finite() || volume <= 0.0 {
        return Err(format!(
            "snapshot zone `{zone}` has non-positive non-finite volume {volume} (cm³)"
        ));
    }
    Ok(())
}

/// Axial extents are informational but must be sane when both are dumped.
fn snap_check_extent(zone: &SnapshotZoneJson) -> Result<(), String> {
    for (label, value) in [
        ("zbottom_cm", zone.zbottom_cm),
        ("ztop_cm", zone.ztop_cm),
        ("temperature_C", zone.temperature_c),
    ] {
        if let Some(v) = value {
            if !v.is_finite() {
                return Err(format!(
                    "snapshot zone `{}` has non-finite {label} {v}",
                    zone.id
                ));
            }
        }
    }
    if let (Some(bottom), Some(top)) = (zone.zbottom_cm, zone.ztop_cm) {
        if top < bottom {
            return Err(format!(
                "snapshot zone `{}` has ztop_cm {top} below zbottom_cm {bottom}",
                zone.id
            ));
        }
    }
    Ok(())
}

/// Flux names are deck tokens referencing schedule items: non-empty, no
/// whitespace.
fn snap_check_flux_name(name: &str) -> Result<(), String> {
    if name.trim().is_empty() {
        return Err("snapshot flux name is empty".to_string());
    }
    if name.chars().any(char::is_whitespace) {
        return Err(format!(
            "snapshot flux `{name}` contains whitespace (deck tokens must not)"
        ));
    }
    Ok(())
}

/// Canonical element symbol when `t` is a bare elemental key (`ZR` → `Zr`),
/// else `None`. Only pure-letter input qualifies, so no nuclide key with
/// mass digits, tags, or separators can collide.
fn snap_elemental_symbol(t: &str) -> Option<String> {
    if t.is_empty() || !t.bytes().all(|b| b.is_ascii_alphabetic()) {
        return None;
    }
    let mut chars = t.chars();
    let mut canon = String::with_capacity(t.len());
    if let Some(first) = chars.next() {
        canon.extend(first.to_uppercase());
    }
    canon.extend(chars.flat_map(|c| c.to_lowercase()));
    if nucleide_nuclei::element_z(&canon).is_some() {
        Some(canon)
    } else {
        None
    }
}

/// Strip one leading ARMI database `n`/`N` when followed by a letter
/// (`nZr` → `Zr`), mirroring the bridge so prefixed elementals are caught.
fn snap_strip_db_prefix(t: &str) -> Option<&str> {
    let mut chars = t.chars();
    let first = chars.next()?;
    let second = chars.next()?;
    if (first == 'n' || first == 'N') && second.is_ascii_alphabetic() {
        Some(&t[1..])
    } else {
        None
    }
}

/// Uppercased key with the database prefix stripped (`nAm242` → `AM242`),
/// for the bare-`AM242` disambiguation check.
fn snap_bare_core(t: &str) -> Option<String> {
    if t.is_empty() {
        return None;
    }
    let u = t.to_ascii_uppercase();
    let core = match u.strip_prefix('N') {
        Some(rest) if rest.starts_with(|c: char| c.is_ascii_alphabetic()) => rest,
        _ => u.as_str(),
    };
    Some(core.to_string())
}

/// Resolve one ARMI-side key: v1 rejections first, then the nuclei bridge.
/// Mirrors the `nucleide-emit` ARMI-input rule.
fn snap_key_to_nucid(key: &str) -> Result<NuclideId, String> {
    let t = key.trim();
    if let Some(sym) =
        snap_elemental_symbol(t).or_else(|| snap_strip_db_prefix(t).and_then(snap_elemental_symbol))
    {
        return Err(format!(
            "snapshot key `{t}` is elemental (`{sym}`): pass post-expansion nuclide \
             number densities; the adapter never reimplements \
             expandElementalMassFracsToNuclides — expand first"
        ));
    }
    match snap_bare_core(t).as_deref() {
        Some("AM242") => {
            return Err(format!(
                "snapshot key `{t}` is ambiguous bare `AM242` (ARMI means the m-state): \
                 pass `AM242M` for Am-242m or `AM242G` for ground explicitly"
            ));
        }
        Some("AM242G") => {
            return NuclideId::new(95, 242, 0).map_err(|e| format!("snapshot key `{t}`: {e}"));
        }
        _ => {}
    }
    nucleide_nuclei::armi::armi_name_to_nucid(t).map_err(|e| format!("snapshot key `{t}`: {e}"))
}

/// Resolve one zone's composition to canonical ids, summing duplicate keys
/// (e.g. `U235` plus `nU235`) the way a volume homogenization would.
/// Densities must be finite and non-negative; empty input means void.
fn snap_resolve_composition(
    zone: &str,
    pairs: &BTreeMap<String, f64>,
) -> Result<BTreeMap<NuclideId, f64>, String> {
    let mut ids = BTreeMap::new();
    for (key, ndens) in pairs {
        let id = snap_key_to_nucid(key)?;
        if !ndens.is_finite() || *ndens < 0.0 {
            return Err(format!(
                "snapshot zone `{zone}` key `{key}` has number density {ndens} \
                 (expected finite atoms/barn-cm >= 0)"
            ));
        }
        *ids.entry(id).or_insert(0.0) += *ndens;
    }
    Ok(ids)
}

/// Merge caller-supplied irradiation-history blocks into the deck.
///
/// Only `schedule` and `pulsehistory` blocks may cross this boundary; any
/// other block kind (geometry, fluxes, cooling, …) is a caller error naming
/// the kind, keeping the dict-in schema versionless.
fn snap_merge_schedule_text(
    deck: &mut nucleide_alara_io::deck::AlaraDeck,
    text: &str,
) -> Result<(), String> {
    use nucleide_alara_io::deck::AlaraDeck;
    let fragment = AlaraDeck::parse(&format!("geometry {SNAPSHOT_GEOMETRY}\n{text}"))
        .map_err(|e| e.to_string())?;
    let mut bad: std::collections::BTreeSet<&str> = std::collections::BTreeSet::new();
    let mut geometries = 0;
    for block in &fragment.blocks {
        match block.kind.as_str() {
            "geometry" => geometries += 1,
            "schedule" | "pulsehistory" => {}
            other => {
                bad.insert(other);
            }
        }
    }
    if geometries > 1 {
        bad.insert("geometry");
    }
    if !bad.is_empty() {
        let kinds: Vec<&str> = bad.into_iter().collect();
        return Err(format!(
            "snapshot schedule_text must hold schedule/pulsehistory blocks only, \
             found {}",
            kinds.join(", ")
        ));
    }
    for block in fragment.blocks {
        if block.kind == "schedule" || block.kind == "pulsehistory" {
            deck.blocks
                .push(nucleide_alara_io::deck::RawBlock { line: 0, ..block });
        }
    }
    deck.schedules = fragment.schedules;
    deck.pulse_histories = fragment.pulse_histories;
    Ok(())
}

/// Synthesize one 1-day-per-flux schedule plus a single-pulse history, so a
/// snapshot without `scheduleText` still derives a workflow and expands.
fn snap_synthesize_schedule(
    deck: &mut nucleide_alara_io::deck::AlaraDeck,
    flux_defs: &[SnapshotFluxJson],
) -> Result<(), String> {
    use nucleide_alara_io::deck::{PulseHistory, PulseLevel, ScheduleDef, ScheduleItemRef};
    if flux_defs.is_empty() {
        return Err("snapshot defines no flux definitions".to_string());
    }
    let items = flux_defs
        .iter()
        .map(|flux| ScheduleItemRef {
            tokens: vec![
                "1".to_string(),
                "d".to_string(),
                flux.name.clone(),
                SNAPSHOT_HISTORY.to_string(),
                "0".to_string(),
                "s".to_string(),
            ],
            line: 0,
        })
        .collect::<Vec<_>>();
    snap_push_block(
        deck,
        "schedule",
        items.iter().map(|item| item.tokens.join(" ")).collect(),
    );
    deck.schedules.push(ScheduleDef {
        name: SNAPSHOT_SCHEDULE.to_string(),
        items,
        line: 0,
    });
    snap_push_block(deck, "pulsehistory", vec!["1 0 s".to_string()]);
    deck.pulse_histories.push(PulseHistory {
        name: SNAPSHOT_HISTORY.to_string(),
        levels: vec![PulseLevel {
            pulses: 1,
            delay_s: 0.0,
        }],
        line: 0,
    });
    Ok(())
}

/// Build a validated [`AlaraDeck`](nucleide_alara_io::deck::AlaraDeck)
/// template from a versionless snapshot.
///
/// Runs deck validation before returning, so dangling mixture, volume, or
/// schedule references surface here. A missing `cooling` block (empty
/// `coolingS`) is *not* a deck error; it fails later in workflow validation,
/// mirroring parsed decks.
fn snap_deck_from_snapshot(
    input: &SnapshotInputJson,
) -> Result<nucleide_alara_io::deck::AlaraDeck, String> {
    use nucleide_alara_io::deck::{
        AlaraDeck, Cooling, FluxDef, Geometry, MatLoading, MatLoadingEntry, Mixture, MixtureEntry,
        OutputDef, VolumeEntry, Volumes,
    };
    if input.zones.is_empty() {
        return Err("snapshot defines no zones".to_string());
    }
    if input.flux_defs.is_empty() {
        return Err("snapshot defines no flux definitions".to_string());
    }
    let mut seen_zones = std::collections::BTreeSet::new();
    for zone in &input.zones {
        snap_check_zone_id(&zone.id)?;
        if !seen_zones.insert(zone.id.as_str()) {
            return Err(format!(
                "snapshot zone `{}` appears more than once",
                zone.id
            ));
        }
        snap_check_volume(&zone.id, zone.volume_cm3)?;
        snap_check_extent(zone)?;
    }
    let mut seen_fluxes = std::collections::BTreeSet::new();
    for flux in &input.flux_defs {
        snap_check_flux_name(&flux.name)?;
        if !seen_fluxes.insert(flux.name.as_str()) {
            return Err(format!(
                "snapshot flux `{}` appears more than once",
                flux.name
            ));
        }
        if flux.file.trim().is_empty() {
            return Err(format!("snapshot flux `{}` names an empty file", flux.name));
        }
        if !flux.scale.is_finite() {
            return Err(format!(
                "snapshot flux `{}` has non-finite scale {}",
                flux.name, flux.scale
            ));
        }
    }
    for time in &input.cooling_s {
        if !time.is_finite() || *time < 0.0 {
            return Err(format!(
                "snapshot cooling time `{time}` is not a non-negative finite value"
            ));
        }
    }

    let mut deck = AlaraDeck::default();
    snap_push_block(&mut deck, "geometry", vec![SNAPSHOT_GEOMETRY.to_string()]);
    deck.geometry = Some(Geometry {
        kind: SNAPSHOT_GEOMETRY.to_string(),
        line: 0,
    });

    let mut volume_entries = Vec::with_capacity(input.zones.len());
    let mut loading_entries = Vec::with_capacity(input.zones.len());
    for zone in &input.zones {
        volume_entries.push(VolumeEntry {
            volume: zone.volume_cm3,
            zone: zone.id.clone(),
        });
        let resolved = snap_resolve_composition(&zone.id, &zone.composition)?;
        if resolved.is_empty() {
            loading_entries.push(MatLoadingEntry {
                zone: zone.id.clone(),
                mixture: "void".to_string(),
            });
        } else {
            let mixture = snap_mixture_name(&zone.id);
            let entries = resolved
                .into_iter()
                .map(|(id, ndens)| MixtureEntry::Element {
                    symbol: id.to_name(),
                    rel_density: 1.0,
                    vol_fraction: ndens,
                })
                .collect::<Vec<_>>();
            snap_push_block(
                &mut deck,
                "mixture",
                entries
                    .iter()
                    .map(|entry| match entry {
                        MixtureEntry::Element {
                            symbol,
                            rel_density,
                            vol_fraction,
                        } => Ok(format!("element {symbol} {rel_density} {vol_fraction}")),
                        other => Err(format!(
                            "internal error: snapshot mixture entry is not an element ({other:?})"
                        )),
                    })
                    .collect::<Result<Vec<_>, String>>()?,
            );
            deck.mixtures.push(Mixture {
                name: mixture.clone(),
                entries,
                line: 0,
            });
            loading_entries.push(MatLoadingEntry {
                zone: zone.id.clone(),
                mixture,
            });
        }
    }
    snap_push_block(
        &mut deck,
        "volumes",
        volume_entries
            .iter()
            .map(|entry| format!("{} {}", entry.volume, entry.zone))
            .collect(),
    );
    deck.volumes = Some(Volumes {
        entries: volume_entries,
        line: 0,
    });
    snap_push_block(
        &mut deck,
        "mat_loading",
        loading_entries
            .iter()
            .map(|entry| format!("{} {}", entry.zone, entry.mixture))
            .collect(),
    );
    deck.mat_loading = Some(MatLoading {
        entries: loading_entries,
        line: 0,
    });

    for flux in &input.flux_defs {
        snap_push_block(
            &mut deck,
            "flux",
            vec![format!(
                "{} {} {} 0 default",
                flux.name, flux.file, flux.scale
            )],
        );
        deck.fluxes.push(FluxDef {
            name: flux.name.clone(),
            file: flux.file.clone(),
            scale: flux.scale,
            skip: 0,
            format: "default".to_string(),
            line: 0,
        });
    }

    match input.schedule_text.as_deref() {
        Some(text) if !text.trim().is_empty() => snap_merge_schedule_text(&mut deck, text)?,
        _ => snap_synthesize_schedule(&mut deck, &input.flux_defs)?,
    }

    if !input.cooling_s.is_empty() {
        snap_push_block(
            &mut deck,
            "cooling",
            input
                .cooling_s
                .iter()
                .map(|time| format!("{time} s"))
                .collect(),
        );
        deck.cooling = Some(Cooling {
            times_s: input.cooling_s.clone(),
            line: 0,
        });
    }

    if let Some(resolution) = input.output.as_deref() {
        let kind = resolution.trim().to_ascii_lowercase();
        if !["interval", "zone", "mixture"].contains(&kind.as_str()) {
            return Err(format!(
                "snapshot output `{resolution}` is not an ALARA resolution \
                 (interval, zone, or mixture)"
            ));
        }
        snap_push_block(&mut deck, "output", vec!["number_density".to_string()]);
        deck.outputs.push(OutputDef {
            resolution: kind,
            entries: vec!["number_density".to_string()],
            line: 0,
        });
    }

    deck.validate().map_err(|e| e.to_string())?;
    Ok(deck)
}

/// Cross-check derived steps against the template deck: every step zone must
/// appear in `mat_loading`, every step flux must be a defined `flux` block,
/// and both the workflow (`cooling_s`) and the deck (`cooling` block) must
/// carry a non-empty cooling history.
fn snap_validate_against(
    deck: &nucleide_alara_io::deck::AlaraDeck,
    steps: &[R2sStepJson],
    cooling_s: &[f64],
) -> Result<(), String> {
    if steps.is_empty() {
        return Err("empty R2S workflow".to_string());
    }
    for step in steps {
        if step.zone.trim().is_empty() {
            return Err("R2S step names an empty zone".to_string());
        }
        if step.flux.trim().is_empty() {
            return Err(format!(
                "R2S step for zone `{}` names an empty flux",
                step.zone
            ));
        }
    }
    let loading = deck
        .mat_loading
        .as_ref()
        .ok_or_else(|| "deck defines no `mat_loading` block".to_string())?;
    let zones: std::collections::BTreeSet<&str> = loading
        .entries
        .iter()
        .map(|entry| entry.zone.as_str())
        .collect();
    for step in steps {
        if !zones.contains(step.zone.as_str()) {
            return Err(format!(
                "workflow zone `{}` is not in deck `mat_loading`",
                step.zone
            ));
        }
    }
    let fluxes: std::collections::BTreeSet<&str> =
        deck.fluxes.iter().map(|flux| flux.name.as_str()).collect();
    for step in steps {
        if !fluxes.contains(step.flux.as_str()) {
            return Err(format!(
                "workflow flux `{}` for zone `{}` is not a defined deck flux",
                step.flux, step.zone
            ));
        }
    }
    if cooling_s.is_empty() {
        return Err("workflow defines no cooling times".to_string());
    }
    match &deck.cooling {
        Some(cooling) if !cooling.times_s.is_empty() => Ok(()),
        _ => Err("deck defines no cooling times".to_string()),
    }
}

/// Emit one deck clone per workflow step with `solve_zones` narrowed to that
/// step's zone (typed field plus the raw block the canonical writer replays).
fn snap_emit_decks(
    template: &nucleide_alara_io::deck::AlaraDeck,
    zones: &[String],
) -> Vec<nucleide_alara_io::deck::AlaraDeck> {
    use nucleide_alara_io::deck::{RawBlock, SolveZones};
    zones
        .iter()
        .map(|zone| {
            let mut deck = template.clone();
            let line = deck.solve_zones.as_ref().map(|z| z.line).unwrap_or(0);
            deck.solve_zones = Some(SolveZones {
                zones: vec![zone.clone()],
                line,
            });
            let mut found = false;
            for raw in deck
                .blocks
                .iter_mut()
                .filter(|block| block.kind == "solve_zones")
            {
                raw.body = vec![zone.clone()];
                found = true;
            }
            if !found {
                deck.blocks.push(RawBlock {
                    kind: "solve_zones".to_string(),
                    line,
                    body: vec![zone.clone()],
                });
            }
            deck
        })
        .collect()
}

/// Build an R2S workflow bundle from a versionless snapshot object.
///
/// `snapshot` mirrors the Python `r2s_from_snapshot` dict: `zones` (list of
/// `{id, volumeCm3, composition: {ARMI-name: ndens}}` with optional
/// `zbottomCm`/`ztopCm`/`material`/`xsType`/`temperatureC`/`flux`),
/// `fluxDefs` (list of `{name, file, scale}`), `coolingS` (list of seconds),
/// plus optional `scheduleText` and `output`. Snake-case aliases
/// (`volume_cm3`, `flux_defs`, `cooling_s`, `schedule_text`, ...) are accepted
/// for every camelCase key.
///
/// Returns `{workflow, deck, decks}`: the workflow summary (same shape as
/// [`r2s_from_deck`]), the canonical template deck text, and one canonical
/// deck text per step. Composition keys follow the emit ARMI-input rule
/// (post-expansion nuclide keys; elemental keys, bare `AM242`, and unknown
/// names throw); densities are atoms/barn-cm. Empty `coolingS` throws via
/// workflow validation.
#[wasm_bindgen(js_name = r2sFromSnapshot)]
pub fn r2s_from_snapshot(snapshot: JsValue) -> Result<JsValue, JsValue> {
    let input: SnapshotInputJson = serde_wasm_bindgen::from_value(snapshot).map_err(js_err)?;
    let deck = snap_deck_from_snapshot(&input).map_err(js_err)?;

    let loading = deck
        .mat_loading
        .as_ref()
        .ok_or_else(|| js_err("deck defines no `mat_loading` block"))?;
    let zones: Vec<&str> = loading
        .entries
        .iter()
        .filter(|entry| !entry.mixture.eq_ignore_ascii_case("void"))
        .map(|entry| entry.zone.as_str())
        .collect();
    if zones.is_empty() {
        return Err(js_err("deck defines no non-`void` zones in `mat_loading`"));
    }
    if deck.fluxes.is_empty() {
        return Err(js_err("deck defines no `flux` blocks"));
    }
    let mut steps: Vec<R2sStepJson> = zones
        .iter()
        .map(|zone| {
            Ok(R2sStepJson {
                zone: (*zone).to_string(),
                flux: r2s_resolve_flux(&deck.fluxes, zone).map_err(js_err)?,
            })
        })
        .collect::<Result<Vec<_>, JsValue>>()?;

    let pins: BTreeMap<&str, &str> = input
        .zones
        .iter()
        .filter_map(|zone| zone.flux.as_deref().map(|flux| (zone.id.as_str(), flux)))
        .collect();
    if !pins.is_empty() {
        let defined: std::collections::BTreeSet<&str> =
            deck.fluxes.iter().map(|flux| flux.name.as_str()).collect();
        for step in &mut steps {
            if let Some(pinned) = pins.get(step.zone.as_str()) {
                if !defined.contains(*pinned) {
                    return Err(js_err(format!(
                        "snapshot zone `{}` pins undefined flux `{pinned}`",
                        step.zone
                    )));
                }
                step.flux = (*pinned).to_string();
            }
        }
    }

    let cooling_s = deck
        .cooling
        .as_ref()
        .map(|cooling| cooling.times_s.clone())
        .unwrap_or_default();
    let top_schedule = r2s_top_schedule(&deck).map_err(js_err)?;

    let mut schedules = Vec::with_capacity(deck.schedules.len());
    for raw in &deck.schedules {
        let mut items = Vec::with_capacity(raw.items.len());
        for entry in &raw.items {
            items.push(r2s_sched_item(&entry.tokens, entry.line).map_err(js_err)?);
        }
        schedules.push(nucleide_alara_io::schedule::ScheduleDef {
            name: raw.name.clone(),
            items,
        });
    }
    let histories: Vec<nucleide_alara_io::schedule::PulseHistory> = deck
        .pulse_histories
        .iter()
        .map(|history| nucleide_alara_io::schedule::PulseHistory {
            name: history.name.clone(),
            levels: history
                .levels
                .iter()
                .map(|level| nucleide_alara_io::schedule::PulseLevel {
                    count: level.pulses,
                    delay_s: level.delay_s,
                })
                .collect(),
        })
        .collect();
    let flat = nucleide_alara_io::schedule::expand_from(&top_schedule, &schedules, &histories)
        .map_err(js_err)?;

    snap_validate_against(&deck, &steps, &cooling_s).map_err(js_err)?;
    let step_zones: Vec<String> = steps.iter().map(|s| s.zone.clone()).collect();
    let decks = snap_emit_decks(&deck, &step_zones);

    to_js(&SnapshotBundleJson {
        workflow: R2sSummary {
            steps,
            cooling_s,
            top_schedule,
            total_s: nucleide_alara_io::total_time(&flat),
        },
        deck: deck.to_string(),
        decks: decks.iter().map(ToString::to_string).collect(),
    })
}

// ---------------------------------------------------------------------------
// Point kinetics (step-reactivity transient + prompt jump)
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct KineticsTransientResult {
    times: Vec<f64>,
    n: Vec<f64>,
    #[serde(rename = "promptJump")]
    prompt_jump: Option<f64>,
    #[serde(rename = "betaTotal")]
    beta_total: f64,
}

fn check_finite_vec(values: &[f64], label: &str) -> Result<(), JsValue> {
    if values.iter().all(|v| v.is_finite()) {
        Ok(())
    } else {
        Err(js_err(format!("{label} must all be finite")))
    }
}

/// Solve a step-reactivity point-kinetics transient and report the
/// prompt-jump estimate alongside the `n(t)` series.
///
/// Thin facade over `nucleide-kinetics`: `betas`/`lambdas`/`lambda_gen`
/// build [`KineticParams`](nucleide_kinetics::KineticParams), the step levels
/// build a `Step` [`Reactivity`](nucleide_kinetics::Reactivity), `times`
/// build a [`TimeGrid`](nucleide_kinetics::TimeGrid), precursors default to
/// the equilibrium (E2) populations for `n0`, and the solve runs with default
/// [`SolverOptions`](nucleide_kinetics::SolverOptions). Returns
/// `{ times, n, promptJump, betaTotal }`; `promptJump` is `null` when the
/// post-step level is at/past prompt critical (`rho_final >= beta`), where
/// the (E4) formula has no solution but the solver still runs.
#[wasm_bindgen(js_name = kineticsTransient)]
#[allow(clippy::too_many_arguments)] // thin JS facade: one scalar per solver input, by design
pub fn kinetics_transient(
    betas: Vec<f64>,
    lambdas: Vec<f64>,
    lambda_gen: f64,
    t_step: f64,
    rho_init: f64,
    rho_final: f64,
    times: Vec<f64>,
    n0: f64,
) -> Result<JsValue, JsValue> {
    check_finite_vec(&betas, "betas")?;
    check_finite_vec(&lambdas, "lambdas")?;
    check_finite_vec(&times, "times")?;
    for (value, label) in [
        (lambda_gen, "lambdaGen"),
        (t_step, "tStep"),
        (rho_init, "rhoInit"),
        (rho_final, "rhoFinal"),
        (n0, "n0"),
    ] {
        if !value.is_finite() {
            return Err(js_err(format!("{label} must be finite")));
        }
    }
    let params =
        nucleide_kinetics::KineticParams::new(betas, lambdas, lambda_gen).map_err(js_err)?;
    let rho = nucleide_kinetics::Reactivity::Step {
        t_step,
        rho_init,
        rho_final,
    };
    let grid = nucleide_kinetics::TimeGrid::new(times).map_err(js_err)?;
    let state = nucleide_kinetics::State::new(&params, n0, None).map_err(js_err)?;
    let sol = nucleide_kinetics::solve(
        &params,
        &rho,
        &grid,
        &state,
        &nucleide_kinetics::SolverOptions::default(),
    )
    .map_err(js_err)?;
    let prompt_jump =
        nucleide_kinetics::prompt_jump(n0, rho_init, rho_final, params.beta_total()).ok();
    to_js(&KineticsTransientResult {
        times: sol.times,
        n: sol.n,
        prompt_jump,
        beta_total: params.beta_total(),
    })
}

// ---------------------------------------------------------------------------
// Spectroscopy (smoothing + peak counting)
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct SpectroscopySmoothResult {
    smoothed: Vec<f64>,
    gross: f64,
    background: f64,
    net: f64,
}

/// Smooth a counts vector and count a peak window on it.
///
/// Thin facade over `nucleide-spectroscopy`: `method` selects the smoothing
/// pass (`"rect3"`, `"rect5"`, `"rect7"` → E1 rectangular with that width,
/// `"five-point"` → E2); `c1`/`c2` delimit the peak window positionally
/// (channels labelled `0..N-1`, E4 half-open gross, E3 `m == 1` background,
/// E5 net). Returns `{ smoothed, gross, background, net }`.
#[wasm_bindgen(js_name = spectroscopySmooth)]
pub fn spectroscopy_smooth(
    counts: Vec<f64>,
    method: &str,
    c1: i32,
    c2: i32,
) -> Result<JsValue, JsValue> {
    check_finite_vec(&counts, "counts")?;
    let smoothed = match method {
        "rect3" => nucleide_spectroscopy::rect_smooth(&counts, 3).map_err(js_err)?,
        "rect5" => nucleide_spectroscopy::rect_smooth(&counts, 5).map_err(js_err)?,
        "rect7" => nucleide_spectroscopy::rect_smooth(&counts, 7).map_err(js_err)?,
        "five-point" => nucleide_spectroscopy::five_point_smooth(&counts).map_err(js_err)?,
        _ => {
            return Err(js_err(format!(
                "unknown smoothing method `{method}` (supported: rect3, rect5, rect7, five-point)"
            )));
        }
    };
    let channels: Vec<f64> = (0..counts.len()).map(|c| c as f64).collect();
    let (c1, c2) = (i64::from(c1), i64::from(c2));
    let gross = nucleide_spectroscopy::gross_count(&counts, &channels, c1, c2).map_err(js_err)?;
    let background =
        nucleide_spectroscopy::calc_bg(&counts, &channels, c1, c2, 1).map_err(js_err)?;
    to_js(&SpectroscopySmoothResult {
        smoothed,
        gross,
        background,
        net: gross - background,
    })
}
