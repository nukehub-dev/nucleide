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

#[derive(Serialize)]
struct MctalTallySummary {
    number: u32,
    #[serde(rename = "particleType")]
    particle_type: i32,
    #[serde(rename = "detectorType")]
    detector_type: Option<i32>,
    /// Effective bin counts per card (`f`, `d`, `u`, `s`, `m`, `c`, `e`, `t`).
    bins: Vec<usize>,
    /// Parsed `(value, rel_error)` pair count.
    pairs: usize,
    /// Sum of tally values (errors excluded).
    total: f64,
}

#[derive(Serialize)]
struct MctalSummary {
    #[serde(rename = "codeName")]
    code_name: String,
    #[serde(rename = "codeVersion")]
    code_version: String,
    #[serde(rename = "nHistories")]
    n_histories: u64,
    #[serde(rename = "tallyNums")]
    tally_nums: Vec<u32>,
    npert: Option<String>,
    tallies: Vec<MctalTallySummary>,
    #[serde(rename = "nCycles")]
    n_cycles: usize,
}

/// Parse an MCNP `MCTAL` file into a JSON summary of its header, standard
/// tally bodies, and kcode cycle count.
#[wasm_bindgen(js_name = parseMctal)]
pub fn parse_mctal(text: &str) -> Result<JsValue, JsValue> {
    let mctal = nucleide_mcnp_io::mctal::Mctal::parse(text).map_err(js_err)?;
    let tallies = mctal
        .tallies
        .iter()
        .map(|t| MctalTallySummary {
            number: t.number,
            particle_type: t.particle_type,
            detector_type: t.detector_type,
            bins: vec![
                t.f.bins(),
                t.d.bins(),
                t.u.bins(),
                t.s.bins(),
                t.m.bins(),
                t.c.bins(),
                t.e.bins(),
                t.t.bins(),
            ],
            pairs: t.vals.len(),
            total: t.total_val(),
        })
        .collect();
    to_js(&MctalSummary {
        code_name: mctal.code_name,
        code_version: mctal.code_version,
        n_histories: mctal.n_histories,
        tally_nums: mctal.tally_nums,
        npert: mctal.npert,
        tallies,
        n_cycles: mctal.n_cycles,
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
    keff_history: Option<Vec<[f64; 2]>>,
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
    let keff_matrix = table.get_matrix("IMP_KEFF").ok();
    let keff = keff_matrix
        .and_then(|m| m.row_f64(0).ok())
        .map(|row| row.iter().take(2).copied().collect());
    // One `[mean, err]` row per burnup block: the full IMP_KEFF matrix, not
    // just the first row. Rows with a single column carry no uncertainty.
    let keff_history = keff_matrix.map(|m| {
        (0..m.rows())
            .filter_map(|r| m.row_f64(r).ok())
            .filter_map(|row| {
                row.first()
                    .copied()
                    .map(|mean| [mean, row.get(1).copied().unwrap_or(0.0)])
            })
            .collect::<Vec<[f64; 2]>>()
    });
    to_js(&SerpentResSummary {
        variable_count: table.len(),
        version,
        title,
        keff,
        keff_history: keff_history.filter(|h| !h.is_empty()),
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
struct SerpentDetSpectrum {
    name: String,
    energy_mid: Vec<f64>,
    values: Vec<f64>,
    errors: Vec<f64>,
}

/// One spectrum per detector value matrix. Serpent 1 rows carry ten bin
/// indices then the tally value, its relative error, and the history count;
/// Serpent 2 rows carry eleven bin indices (or ten without time bins) then
/// the tally value and its relative error. Serpent 1 detectors are recognized
/// by their `<name>_VALS`/`<name>_EBINS` bin-count scalars, mirroring the
/// crate's reshape rules. The energy midpoint column of the matching
/// `DET<name>E` grid pairs with the rows 1:1; grids that do not align
/// (multi-axis detectors) leave `energy_mid` empty.
fn serpent_detector_spectra(table: &nucleide_serpent_io::Table) -> Vec<SerpentDetSpectrum> {
    let mut spectra = Vec::new();
    for (name, entry) in table.iter() {
        let nucleide_serpent_io::Entry::Matrix(m) = entry else {
            continue;
        };
        if !name.starts_with("DET") || m.cols() < 12 {
            continue;
        }
        let serpent1 = table.contains_key(format!("{name}_VALS").as_str())
            || table.contains_key(format!("{name}_EBINS").as_str());
        let (value_col, error_col) = if serpent1 {
            (m.cols() - 3, m.cols() - 2)
        } else {
            (m.cols() - 2, m.cols() - 1)
        };
        let mut values = Vec::with_capacity(m.rows());
        let mut errors = Vec::with_capacity(m.rows());
        let mut numeric = true;
        for r in 0..m.rows() {
            match m.row_f64(r) {
                Ok(row) => {
                    values.push(row[value_col]);
                    errors.push(row[error_col]);
                }
                Err(_) => {
                    numeric = false;
                    break;
                }
            }
        }
        if !numeric {
            continue;
        }
        let energy_mid = table
            .get_matrix(format!("{name}E"))
            .ok()
            .filter(|g| g.cols() == 3 && g.rows() == m.rows())
            .map(|g| {
                (0..g.rows())
                    .filter_map(|r| g.row_f64(r).ok())
                    .filter_map(|row| row.get(2).copied())
                    .collect()
            })
            .unwrap_or_default();
        spectra.push(SerpentDetSpectrum {
            name: name.clone(),
            energy_mid,
            values,
            errors,
        });
    }
    spectra
}

#[derive(Serialize)]
struct SerpentDetSummary {
    variable_count: usize,
    detectors: Vec<String>,
    spectra: Vec<SerpentDetSpectrum>,
    variables: Vec<SerpentVariableJson>,
}

/// Parse a Serpent `_det.m` detector file into a JSON summary.
///
/// `detectors` lists detector value matrices only (DET-prefixed, ≥12
/// columns): bin-grid matrices such as `DET<name>E` are inputs to the
/// detectors, not detectors themselves.
#[wasm_bindgen(js_name = parseSerpentDet)]
pub fn parse_serpent_det(text: &str) -> Result<JsValue, JsValue> {
    let table = nucleide_serpent_io::parse_det(text).map_err(js_err)?;
    let detectors = table
        .iter()
        .filter(|(name, entry)| {
            name.starts_with("DET")
                && matches!(entry, nucleide_serpent_io::Entry::Matrix(m) if m.cols() >= 12)
        })
        .map(|(name, _)| name.clone())
        .collect();
    to_js(&SerpentDetSummary {
        variable_count: table.len(),
        detectors,
        spectra: serpent_detector_spectra(&table),
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
// PARTISN writer (structured deck dicts, exact keys)
// ---------------------------------------------------------------------------

/// One PARTISN zone of the structured deck dict.
///
/// Exact keys (no aliases): `id`, `material`, `isotxs_labels`, `density`.
/// Mirrors the Python `partisn_render`/`partisn_validate` zone shape.
#[derive(Deserialize)]
struct PartisnZoneJson {
    id: u32,
    material: String,
    isotxs_labels: Vec<String>,
    density: f64,
}

/// Minimal PARTISN deck of the structured deck dict.
///
/// Exact keys (no aliases): `title`, `dim`, `zones`, `source` (optional,
/// `null`/absent omits the SOURCE card). Mirrors the Python
/// `partisn_render`/`partisn_validate` deck shape.
#[derive(Deserialize)]
struct PartisnDeckJson {
    title: String,
    dim: u8,
    zones: Vec<PartisnZoneJson>,
    #[serde(default)]
    source: Option<String>,
}

fn partisn_deck_from_json(deck: JsValue) -> Result<nucleide_cccc_io::PartisnDeck, JsValue> {
    let parsed: PartisnDeckJson = serde_wasm_bindgen::from_value(deck).map_err(js_err)?;
    Ok(nucleide_cccc_io::PartisnDeck {
        title: parsed.title,
        dim: parsed.dim,
        zones: parsed
            .zones
            .into_iter()
            .map(|z| nucleide_cccc_io::partisn::PartisnZone {
                id: z.id,
                material: z.material,
                isotxs_labels: z.isotxs_labels,
                density: z.density,
            })
            .collect(),
        source: parsed.source,
    })
}

/// Render a structured PARTISN deck dict to PARTISN input text.
///
/// `deck` uses the exact keys `title`, `dim` (1|2|3), `zones` (each with
/// exact keys `id`, `material`, `isotxs_labels`, `density`), and optional
/// `source` (`null`/absent omits the SOURCE card).
///
/// An empty `zones` list renders by design (TITLE/DIM/END with no ZONE
/// cards); callers must validate separately via [`partisn_validate`].
#[wasm_bindgen(js_name = partisnRender)]
pub fn partisn_render(deck: JsValue) -> Result<String, JsValue> {
    Ok(partisn_deck_from_json(deck)?.render())
}

/// Validate a structured PARTISN deck dict against ISOTXS library text.
///
/// Same exact deck keys as [`partisn_render`]. Throws when `zones` is empty,
/// when `dim` is not 1/2/3, or when a zone names an ISOTXS label absent from
/// the library.
#[wasm_bindgen(js_name = partisnValidate)]
pub fn partisn_validate(deck: JsValue, isotxs_text: &str) -> Result<(), JsValue> {
    let rust_deck = partisn_deck_from_json(deck)?;
    if rust_deck.zones.is_empty() {
        return Err(js_err(
            "PARTISN deck has no zones (empty zones list is vacuous)",
        ));
    }
    let lib = nucleide_cccc_io::IsotxsLib::parse(isotxs_text).map_err(js_err)?;
    rust_deck.validate(&lib).map_err(js_err)
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
// R2S voxel tags (local port — no `nucleide-r2s` dependency)
// ---------------------------------------------------------------------------
//
// `nucleide-r2s` depends on `nucleide-depletion` with default features, which
// would re-enable Rayon through Cargo feature unification and break this
// `wasm32-unknown-unknown` build (the same reason the snapshot adapter above
// is re-implemented here). The per-voxel tag math from `nucleide_r2s::tags`
// (`VoxelTags` + `tag_zone_totals` + `split_zone_totals` +
// `photon_groups_at` + `sum_group_strengths`,
// `crates/r2s/src/tags.rs:25-197`) is therefore copy-ported here over
// `nucleide-alara-io` photon types only, drifting with the owner crate by
// design. Only the ALARA photon types cross this boundary
// (`PhotonSource`/`PhotonGroup` in `photon_groups_at`/`sum_group_strengths`);
// the zone-source rows are a local stand-in for
// `nucleide_r2s::photon::ZonePhotonSource` (only `total()` is needed).
// Errors surface as strings via `js_err`; dict inputs arrive via
// `serde_wasm_bindgen` like the snapshot adapter above.

/// Browser-demo cap on voxel counts (mirrors the RTFLUX values cap).
const MAX_VOXEL_TAGS: usize = 200;

/// Per-voxel photon-source tags over `n_voxels` voxels.
///
/// Local port of `nucleide_r2s::tags::VoxelTags`; see the section header.
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
struct VoxelTags {
    /// Zone count the indices in `zone_of_voxel` refer to.
    n_zones: usize,
    /// Zone index per voxel (length `n_voxels`, each `< n_zones`).
    zone_of_voxel: Vec<usize>,
    /// Total source strength per voxel (length `n_voxels`).
    source_strength: Vec<f64>,
    /// Decay (cooling) time per voxel in seconds (length `n_voxels`).
    decay_time_s: Vec<f64>,
}

impl VoxelTags {
    /// Voxel count.
    fn n_voxels(&self) -> usize {
        self.zone_of_voxel.len()
    }

    /// Build tags from parallel arrays, checking lengths, zone-index
    /// bounds, and finiteness (non-finite strengths or times are rejected;
    /// sign conventions stay caller-side).
    fn new(
        n_zones: usize,
        zone_of_voxel: Vec<usize>,
        source_strength: Vec<f64>,
        decay_time_s: Vec<f64>,
    ) -> Result<Self, String> {
        let n = zone_of_voxel.len();
        if source_strength.len() != n || decay_time_s.len() != n {
            return Err(format!(
                "voxel tag arrays must share one length, found zone_of_voxel={n}, \
                 source_strength={}, decay_time_s={}",
                source_strength.len(),
                decay_time_s.len()
            ));
        }
        if let Some(bad) = zone_of_voxel.iter().position(|z| *z >= n_zones) {
            return Err(format!(
                "voxel {bad} names zone {} of only {n_zones}",
                zone_of_voxel[bad]
            ));
        }
        if source_strength.iter().any(|v| !v.is_finite()) {
            return Err("voxel source strengths must be finite".to_string());
        }
        if decay_time_s.iter().any(|v| !v.is_finite()) {
            return Err("voxel decay times must be finite".to_string());
        }
        Ok(VoxelTags {
            n_zones,
            zone_of_voxel,
            source_strength,
            decay_time_s,
        })
    }

    /// Sum of all voxel strengths.
    fn total_strength(&self) -> f64 {
        self.source_strength.iter().sum()
    }

    /// Sum of strengths over voxels naming `zone`.
    fn zone_total(&self, zone: usize) -> f64 {
        self.zone_of_voxel
            .iter()
            .zip(&self.source_strength)
            .filter(|(z, _)| **z == zone)
            .map(|(_, v)| *v)
            .sum()
    }
}

/// Zone photon-source row (local stand-in for
/// `nucleide_r2s::photon::ZonePhotonSource`; only `total()` feeds the tag
/// math below).
#[derive(Debug, Clone, Deserialize)]
struct VoxelZoneSource {
    #[allow(dead_code)] // positional tag math never names zones; kept for the ported shape.
    zone: String,
    groups: Vec<f64>,
}

impl VoxelZoneSource {
    /// Total photon strength.
    fn total(&self) -> f64 {
        self.groups.iter().sum()
    }
}

/// Copy each zone total onto every voxel of that zone (tag-as-attribute:
/// voxel strengths in a zone sum to `count * total`, not `total`).
/// Decay times are `0.0` (shutdown); voxel count is `zone_of_voxel.len()`.
fn voxel_tag_zone_totals(
    zones: &[VoxelZoneSource],
    zone_of_voxel: &[usize],
) -> Result<VoxelTags, String> {
    let totals: Vec<f64> = zones.iter().map(VoxelZoneSource::total).collect();
    let strengths: Vec<f64> = zone_of_voxel
        .iter()
        .map(|z| {
            totals
                .get(*z)
                .copied()
                .ok_or_else(|| format!("voxel names zone {z} of only {} zones", totals.len()))
        })
        .collect::<Result<_, _>>()?;
    VoxelTags::new(
        zones.len(),
        zone_of_voxel.to_vec(),
        strengths,
        vec![0.0; zone_of_voxel.len()],
    )
}

/// Distribute each zone total conservatively over its voxels (voxel =
/// `total / voxel count in zone`; zones with no voxels contribute nothing).
/// The voxel strengths in a zone sum back to the zone total exactly when
/// the division is exact, and [`VoxelTags::total_strength`] equals the sum
/// of zone totals over zones owning at least one voxel.
fn voxel_split_zone_totals(
    zones: &[VoxelZoneSource],
    zone_of_voxel: &[usize],
) -> Result<VoxelTags, String> {
    let totals: Vec<f64> = zones.iter().map(VoxelZoneSource::total).collect();
    let mut counts = vec![0usize; zones.len()];
    for z in zone_of_voxel {
        counts
            .get_mut(*z)
            .ok_or_else(|| format!("voxel names zone {z} of only {} zones", zones.len()))?;
        counts[*z] += 1;
    }
    let strengths: Vec<f64> = zone_of_voxel
        .iter()
        .map(|z| totals[*z] / counts[*z] as f64)
        .collect();
    VoxelTags::new(
        zones.len(),
        zone_of_voxel.to_vec(),
        strengths,
        vec![0.0; zone_of_voxel.len()],
    )
}

/// Select `.photonSrc` group spectra for caller-named `nuclides` at exactly
/// `time_s` seconds (exact match, the same convention the R2S shutdown
/// assembly uses for `0.0`). Unknown nuclides select nothing; `TOTAL`
/// aggregates are ordinary rows — pass `"TOTAL"` to select them. No
/// rescaling is applied: strengths keep the file's own normalization and
/// stay caller-side to interpret.
fn voxel_photon_groups_at<'a>(
    photon: &'a nucleide_alara_io::photon::PhotonSource,
    nuclides: &[&str],
    time_s: f64,
) -> Vec<&'a nucleide_alara_io::photon::PhotonGroup> {
    photon
        .groups
        .iter()
        .filter(|g| g.time_s == time_s && nuclides.contains(&g.nuclide.as_str()))
        .collect()
}

/// Add selected group spectra element-wise (ALARA group order preserved).
/// Empty selection yields an empty spectrum; ragged group counts are an
/// error. The sums conserve the input total exactly up to float rounding:
/// `sums.iter().sum()` equals the sum of the inputs' totals.
fn voxel_sum_group_strengths(
    groups: &[&nucleide_alara_io::photon::PhotonGroup],
) -> Result<Vec<f64>, String> {
    let mut sums: Vec<f64> = Vec::new();
    for g in groups {
        if sums.is_empty() {
            sums.clone_from(&g.strengths);
        } else {
            if sums.len() != g.strengths.len() {
                return Err(format!(
                    "photon group spectra have ragged group counts ({} vs {})",
                    sums.len(),
                    g.strengths.len()
                ));
            }
            for (s, v) in sums.iter_mut().zip(&g.strengths) {
                *s += *v;
            }
        }
    }
    Ok(sums)
}

#[derive(Deserialize)]
struct VoxelTagsInputJson {
    totals: Vec<f64>,
    #[serde(rename = "zoneOfVoxel", alias = "zone_of_voxel")]
    zone_of_voxel: Vec<usize>,
    #[serde(default)]
    split: bool,
}

#[derive(Serialize)]
struct VoxelTagsResult {
    n_zones: usize,
    n_voxels: usize,
    zone_of_voxel: Vec<usize>,
    source_strength: Vec<f64>,
    decay_time_s: Vec<f64>,
    zone_totals: Vec<f64>,
    total: f64,
}

/// Tag per-voxel source strengths from per-zone totals.
///
/// Dict-in `{totals, zoneOfVoxel, split?}` (snake_case aliases accepted):
/// `totals` carries one total source strength per zone; with `split=false`
/// every voxel copies its zone total (tag-as-attribute), with `split=true`
/// each zone total is divided conservatively over its voxels (mirrors the
/// Python `r2s_tag_zone_strength` facade). Voxel counts above
/// [`MAX_VOXEL_TAGS`] are rejected for the browser demo. Returns
/// `{n_zones, n_voxels, zone_of_voxel, source_strength, decay_time_s,
/// zone_totals, total}` with shutdown (`0.0`) decay times.
#[wasm_bindgen(js_name = voxelTagsFromTotals)]
pub fn voxel_tags_from_totals(input: JsValue) -> Result<JsValue, JsValue> {
    let parsed: VoxelTagsInputJson = serde_wasm_bindgen::from_value(input).map_err(js_err)?;
    if parsed.zone_of_voxel.is_empty() {
        return Err(js_err(
            "zone_of_voxel is empty (tagging zero voxels is vacuous)",
        ));
    }
    if parsed.zone_of_voxel.len() > MAX_VOXEL_TAGS {
        return Err(js_err(format!(
            "voxel count {} exceeds the demo cap of {MAX_VOXEL_TAGS}",
            parsed.zone_of_voxel.len()
        )));
    }
    let zones: Vec<VoxelZoneSource> = parsed
        .totals
        .into_iter()
        .enumerate()
        .map(|(i, total)| VoxelZoneSource {
            zone: format!("zone{i}"),
            groups: if total == 0.0 {
                Vec::new()
            } else {
                vec![total]
            },
        })
        .collect();
    let tags = if parsed.split {
        voxel_split_zone_totals(&zones, &parsed.zone_of_voxel)
    } else {
        voxel_tag_zone_totals(&zones, &parsed.zone_of_voxel)
    }
    .map_err(js_err)?;
    let total = tags.total_strength();
    let n_voxels = tags.n_voxels();
    let zone_totals: Vec<f64> = (0..tags.n_zones).map(|z| tags.zone_total(z)).collect();
    to_js(&VoxelTagsResult {
        n_zones: tags.n_zones,
        n_voxels,
        zone_of_voxel: tags.zone_of_voxel,
        source_strength: tags.source_strength,
        decay_time_s: tags.decay_time_s,
        zone_totals,
        total,
    })
}

#[derive(Deserialize)]
struct VoxelPhotonInputJson {
    #[serde(rename = "photonText", alias = "photon_text")]
    photon_text: String,
    #[serde(default)]
    nuclides: Vec<String>,
    #[serde(rename = "timeS", alias = "time_s", default)]
    time_s: f64,
}

#[derive(Serialize)]
struct VoxelPhotonGroupJson {
    nuclide: String,
    time_s: f64,
    strengths: Vec<f64>,
}

#[derive(Serialize)]
struct VoxelPhotonResult {
    groups: Vec<VoxelPhotonGroupJson>,
    sums: Vec<f64>,
    total: f64,
}

/// Select and sum `.photonSrc` group spectra for `nuclides` at `time_s`.
///
/// Dict-in `{photonText, nuclides, timeS}` (snake_case aliases accepted):
/// parses ALARA photon-source text, keeps rows matching the named nuclides
/// at exactly `time_s` seconds (shutdown `0.0`), and adds them element-wise
/// in ALARA group order (mirrors the Python `r2s_photon_group_sums`
/// facade). Returns `{groups, sums, total}`; no rescaling is applied.
#[wasm_bindgen(js_name = voxelPhotonSums)]
pub fn voxel_photon_sums(input: JsValue) -> Result<JsValue, JsValue> {
    let parsed: VoxelPhotonInputJson = serde_wasm_bindgen::from_value(input).map_err(js_err)?;
    if !parsed.time_s.is_finite() {
        return Err(js_err(format!(
            "timeS must be finite (got {})",
            parsed.time_s
        )));
    }
    if parsed.nuclides.is_empty() {
        return Err(js_err("nuclides is empty (no rows can match)"));
    }
    let source =
        nucleide_alara_io::photon::PhotonSource::from_str(&parsed.photon_text).map_err(js_err)?;
    let names: Vec<&str> = parsed.nuclides.iter().map(String::as_str).collect();
    let at = voxel_photon_groups_at(&source, &names, parsed.time_s);
    let sums = voxel_sum_group_strengths(&at).map_err(js_err)?;
    to_js(&VoxelPhotonResult {
        groups: at
            .iter()
            .map(|g| VoxelPhotonGroupJson {
                nuclide: g.nuclide.clone(),
                time_s: g.time_s,
                strengths: g.strengths.clone(),
            })
            .collect(),
        total: sums.iter().sum(),
        sums,
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
// Tritium transport (permeation breakthrough curve)
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct TritiumBreakthroughResult {
    times: Vec<f64>,
    #[serde(rename = "fluxOverJss")]
    flux_over_jss: Vec<f64>,
    #[serde(rename = "tLag")]
    t_lag: f64,
    #[serde(rename = "jss")]
    jss: f64,
}

/// Solve a trap-free permeation transient and report the normalized outlet
/// flux alongside the time lag.
///
/// Thin facade over `nucleide-tritium`: the slab `0 ≤ x ≤ L` starts empty
/// with Dirichlet ends `c(0) = c0`, `c(L) = 0` (the G2 setup) at uniform
/// diffusivity `diffusivity`, and the solve runs on a fixed 200-cell grid
/// with default [`SolverOptions`](nucleide_tritium::SolverOptions) and an
/// internal step capped at `(L²/D)/2000`. Returns `{ times, fluxOverJss,
/// tLag, jss }` with `tLag = L²/6D` (G2-lag) and `jss = D·c0/L` (G1).
#[wasm_bindgen(js_name = tritiumBreakthrough)]
pub fn tritium_breakthrough(
    length: f64,
    diffusivity: f64,
    c0: f64,
    times: Vec<f64>,
) -> Result<JsValue, JsValue> {
    check_finite_vec(&times, "times")?;
    for (value, label) in [(length, "length"), (diffusivity, "diffusivity"), (c0, "c0")] {
        if !value.is_finite() {
            return Err(js_err(format!("{label} must be finite")));
        }
    }
    let params = nucleide_tritium::TransportParams::new(
        length,
        200,
        diffusivity,
        0.0,
        vec![],
        vec![500.0],
        vec![],
    )
    .map_err(js_err)?;
    let left = nucleide_tritium::Boundary::dirichlet(c0).map_err(js_err)?;
    let right = nucleide_tritium::Boundary::dirichlet(0.0).map_err(js_err)?;
    let grid = nucleide_tritium::TimeGrid::new(times).map_err(js_err)?;
    let initial = nucleide_tritium::InitialState::zeros(&params);
    let diffusive = length * length / diffusivity;
    let sol = nucleide_tritium::solve(
        &params,
        &left,
        &right,
        &grid,
        &initial,
        &nucleide_tritium::SolverOptions {
            dt_max: diffusive / 2000.0,
            ..Default::default()
        },
    )
    .map_err(js_err)?;
    let jss = diffusivity * c0 / length;
    let flux_over_jss: Vec<f64> = sol.flux_right.iter().map(|f| f / jss).collect();
    to_js(&TritiumBreakthroughResult {
        times: sol.times,
        flux_over_jss,
        t_lag: nucleide_tritium::time_lag(length, diffusivity).map_err(js_err)?,
        jss,
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

// ---------------------------------------------------------------------------
// UQ-lite sampling (seeded MVN over caller-supplied blocks + moments)
// ---------------------------------------------------------------------------

/// Demo cap on draw count: the browser slice stays cheap (small-n only).
const MAX_UQ_SAMPLES: usize = 5000;

#[derive(Serialize)]
struct UqSampleResult {
    samples: Vec<Vec<f64>>,
    #[serde(rename = "sampleMean")]
    sample_mean: Vec<f64>,
    #[serde(rename = "sampleCov")]
    sample_cov: Vec<Vec<f64>>,
    method: String,
    #[serde(rename = "minEigen")]
    min_eigen: Option<f64>,
    #[serde(rename = "maxEigen")]
    max_eigen: Option<f64>,
}

/// Draw seeded multivariate-normal samples over a caller-supplied covariance
/// block and report the sample moments alongside the inputs.
///
/// Thin facade over `nucleide-linalg` `sample` (`sample_mvn` +
/// `sample_mean` + `sample_cov`; the decay perturbers stay Python-side and
/// fission yields stay named-open): `mean` is a plain array, `cov` a nested
/// array, `n` the draw count (capped at [`MAX_UQ_SAMPLES`] for the browser),
/// `seed` a plain number (integer-valued, passed through as `u64`).
/// Finiteness and shape validation come from the crate's errors, mapped to
/// `JsValue`. Returns `{ samples, sampleMean, sampleCov, method, minEigen,
/// maxEigen }`; `method` is `"cholesky"` or `"eigen_clip"` (with the
/// unclipped extremes; `null` on the Cholesky path).
#[wasm_bindgen(js_name = uqSample)]
pub fn uq_sample(mean: Vec<f64>, cov: JsValue, n: usize, seed: f64) -> Result<JsValue, JsValue> {
    let cov: Vec<Vec<f64>> = serde_wasm_bindgen::from_value(cov).map_err(js_err)?;
    if !seed.is_finite() || seed < 0.0 || seed.fract() != 0.0 {
        return Err(js_err(format!(
            "seed must be a finite non-negative integer (got {seed})"
        )));
    }
    if seed >= 9_007_199_254_740_992.0 {
        return Err(js_err(format!(
            "seed must be below 2^53 for exact f64 integer precision (got {seed})"
        )));
    }
    if n < 2 {
        return Err(js_err(format!(
            "n = {n} needs n >= 2 (sample moments need at least 2 draws; \
             the kernel supports n = 1 via Rust/Python)"
        )));
    }
    if n > MAX_UQ_SAMPLES {
        return Err(js_err(format!(
            "n = {n} exceeds the demo cap of {MAX_UQ_SAMPLES} draws"
        )));
    }
    let set = nucleide_linalg::sample::sample_mvn(&mean, &cov, n, seed as u64).map_err(js_err)?;
    let sample_mean = nucleide_linalg::sample::sample_mean(&set.samples).map_err(js_err)?;
    let sample_cov = nucleide_linalg::sample::sample_cov(&set.samples).map_err(js_err)?;
    let (method, min_eigen, max_eigen) = match &set.method {
        nucleide_linalg::sample::FactorMethod::Cholesky => ("cholesky".to_string(), None, None),
        nucleide_linalg::sample::FactorMethod::EigenClip {
            min_eigen,
            max_eigen,
        } => ("eigen_clip".to_string(), Some(*min_eigen), Some(*max_eigen)),
    };
    to_js(&UqSampleResult {
        samples: set.samples,
        sample_mean,
        sample_cov,
        method,
        min_eigen,
        max_eigen,
    })
}

/// Draw seeded Latin-hypercube samples over a caller-supplied covariance
/// block and report the sample moments alongside the inputs.
///
/// Thin facade over `nucleide-linalg` `sample` (`sample_lhs` +
/// `sample_mean` + `sample_cov`): stratified uniforms (one jittered draw
/// per stratum per dimension) through the hand-rolled inverse-normal CDF,
/// then the shared factor path. `mean` is a plain array, `cov` a nested
/// array, `n` the draw count (capped at [`MAX_UQ_SAMPLES`] for the browser),
/// `seed` a plain number (integer-valued, passed through as `u64` — no
/// BigInt types cross this boundary). Returns `{ samples, sampleMean,
/// sampleCov, method, minEigen, maxEigen }` like [`uq_sample`].
#[wasm_bindgen(js_name = sampleLhs)]
pub fn sample_lhs(mean: Vec<f64>, cov: JsValue, n: usize, seed: f64) -> Result<JsValue, JsValue> {
    let cov: Vec<Vec<f64>> = serde_wasm_bindgen::from_value(cov).map_err(js_err)?;
    if !seed.is_finite() || seed < 0.0 || seed.fract() != 0.0 {
        return Err(js_err(format!(
            "seed must be a finite non-negative integer (got {seed})"
        )));
    }
    if seed >= 9_007_199_254_740_992.0 {
        return Err(js_err(format!(
            "seed must be below 2^53 for exact f64 integer precision (got {seed})"
        )));
    }
    if n < 2 {
        return Err(js_err(format!(
            "n = {n} needs n >= 2 (sample moments need at least 2 draws; \
             the kernel supports n = 1 via Rust/Python)"
        )));
    }
    if n > MAX_UQ_SAMPLES {
        return Err(js_err(format!(
            "n = {n} exceeds the demo cap of {MAX_UQ_SAMPLES} draws"
        )));
    }
    let set = nucleide_linalg::sample::sample_lhs(&mean, &cov, n, seed as u64).map_err(js_err)?;
    let sample_mean = nucleide_linalg::sample::sample_mean(&set.samples).map_err(js_err)?;
    let sample_cov = nucleide_linalg::sample::sample_cov(&set.samples).map_err(js_err)?;
    let (method, min_eigen, max_eigen) = match &set.method {
        nucleide_linalg::sample::FactorMethod::Cholesky => ("cholesky".to_string(), None, None),
        nucleide_linalg::sample::FactorMethod::EigenClip {
            min_eigen,
            max_eigen,
        } => ("eigen_clip".to_string(), Some(*min_eigen), Some(*max_eigen)),
    };
    to_js(&UqSampleResult {
        samples: set.samples,
        sample_mean,
        sample_cov,
        method,
        min_eigen,
        max_eigen,
    })
}

// ---------------------------------------------------------------------------
// MCPL particle lists (bytes-based; no filesystem in the browser)
// ---------------------------------------------------------------------------
//
// Spike note (A1): `nucleide-mcpl-io` pulls `thiserror` + `flate2` (pure-Rust
// miniz_oxide/crc32fast, no `std::fs`) + `nucleide-mcnp-io` (already a WASM
// dep). `cargo check -p nucleide-wasm --target wasm32-unknown-unknown`
// passes. Gzip is detected by magic bytes (`1f 8b`), never by suffix, and
// the `open`/`write_to_path` path-based APIs are reference only (Python
// shape at `bindings/python/src/lib.rs`): WASM uses `from_bytes` /
// `encode_file` / `ssw2mcpl_bytes` plus `SurfSrc::from_bytes` and the SSW
// `write_to(&mut Vec<u8>, …)` writer form.

/// Demo cap on surfaced particles: the browser slice stays cheap.
const MAX_MCPL_PARTICLES: usize = 200;

fn maybe_gunzip_bytes(data: &[u8]) -> Result<Vec<u8>, JsValue> {
    if data.len() >= 2 && data[0] == 0x1f && data[1] == 0x8b {
        use std::io::Read as _;
        let mut dec = flate2::read::GzDecoder::new(data);
        let mut out = Vec::new();
        dec.read_to_end(&mut out)
            .map_err(|e| js_err(format!("gzip decode: {e}")))?;
        Ok(out)
    } else {
        Ok(data.to_vec())
    }
}

fn gzip_bytes(data: &[u8]) -> Result<Vec<u8>, JsValue> {
    use std::io::Write as _;
    let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    enc.write_all(data)
        .map_err(|e| js_err(format!("gzip encode: {e}")))?;
    enc.finish()
        .map_err(|e| js_err(format!("gzip finish: {e}")))
}

#[derive(Serialize)]
struct McplParticleJson {
    ekin: f64,
    position: [f64; 3],
    direction: [f64; 3],
    time: f64,
    weight: f64,
    pdgcode: i32,
    userflags: u32,
}

#[derive(Serialize)]
struct McplBlobJson {
    key: String,
    len: usize,
}

#[derive(Serialize)]
struct McplSummary {
    version: u16,
    #[serde(rename = "nparticles")]
    nparticles: u64,
    srcname: String,
    comments: Vec<String>,
    #[serde(rename = "hasUserflags")]
    has_userflags: bool,
    #[serde(rename = "hasPolarisation")]
    has_polarisation: bool,
    #[serde(rename = "doublePrec")]
    double_prec: bool,
    #[serde(rename = "universalPdgcode")]
    universal_pdgcode: Option<i32>,
    #[serde(rename = "universalWeight")]
    universal_weight: Option<f64>,
    blobs: Vec<McplBlobJson>,
    particles: Vec<McplParticleJson>,
    truncated: bool,
}

fn mcpl_summary_of(file: &nucleide_mcpl_io::McplFile) -> Result<McplSummary, JsValue> {
    let ps = file.particles().map_err(js_err)?;
    let truncated = ps.len() > MAX_MCPL_PARTICLES;
    let particles: Vec<McplParticleJson> = ps
        .iter()
        .take(MAX_MCPL_PARTICLES)
        .map(|p| McplParticleJson {
            ekin: p.ekin,
            position: p.position,
            direction: p.direction,
            time: p.time,
            weight: p.weight,
            pdgcode: p.pdgcode,
            userflags: p.userflags,
        })
        .collect();
    Ok(McplSummary {
        version: file.header.version,
        nparticles: file.header.nparticles,
        srcname: file.header.srcname.clone(),
        comments: file.header.comments.clone(),
        has_userflags: file.header.has_userflags,
        has_polarisation: file.header.has_polarisation,
        double_prec: file.header.double_prec,
        universal_pdgcode: file.header.universal_pdgcode,
        universal_weight: file.header.universal_weight,
        blobs: file
            .header
            .blobs
            .iter()
            .map(|b| McplBlobJson {
                key: b.key.clone(),
                len: b.data.len(),
            })
            .collect(),
        particles,
        truncated,
    })
}

/// Parse MCPL bytes (`Uint8Array`) into a capped JSON summary.
///
/// Gzip is sniffed by magic bytes, never by suffix. Particle lists are
/// capped at [`MAX_MCPL_PARTICLES`] (`truncated: true` when capped).
#[wasm_bindgen(js_name = readMcpl)]
pub fn read_mcpl(bytes: js_sys::Uint8Array) -> Result<JsValue, JsValue> {
    let raw = bytes.to_vec();
    let data = maybe_gunzip_bytes(&raw)?;
    let file = nucleide_mcpl_io::McplFile::from_bytes(data).map_err(js_err)?;
    to_js(&mcpl_summary_of(&file)?)
}

#[derive(Deserialize)]
struct McplHeaderJson {
    #[serde(default = "default_srcname")]
    srcname: String,
    #[serde(default)]
    comments: Vec<String>,
    #[serde(default, rename = "hasUserflags", alias = "has_userflags")]
    has_userflags: bool,
    #[serde(default, rename = "hasPolarisation", alias = "has_polarisation")]
    has_polarisation: bool,
    #[serde(default, rename = "doublePrec", alias = "double_prec")]
    double_prec: bool,
    #[serde(default, rename = "universalPdgcode", alias = "universal_pdgcode")]
    universal_pdgcode: Option<i32>,
    #[serde(default, rename = "universalWeight", alias = "universal_weight")]
    universal_weight: Option<f64>,
    #[serde(default)]
    gzip: bool,
}

fn default_srcname() -> String {
    "nucleide-wasm".to_string()
}

#[derive(Deserialize)]
struct McplParticleIn {
    #[serde(default)]
    ekin: f64,
    #[serde(default)]
    position: Option<[f64; 3]>,
    #[serde(default)]
    direction: Option<[f64; 3]>,
    #[serde(default)]
    time: f64,
    #[serde(default = "default_weight")]
    weight: f64,
    #[serde(default = "default_pdg")]
    pdgcode: i32,
    #[serde(default)]
    userflags: u32,
}

fn default_weight() -> f64 {
    1.0
}
fn default_pdg() -> i32 {
    2112
}

/// Encode MCPL bytes from a header object plus particle rows.
///
/// Returns a `Uint8Array` (gzip-compressed when `header.gzip` is set).
#[wasm_bindgen(js_name = writeMcpl)]
pub fn write_mcpl(header: JsValue, particles: JsValue) -> Result<js_sys::Uint8Array, JsValue> {
    let h: McplHeaderJson = serde_wasm_bindgen::from_value(header).map_err(js_err)?;
    let rows: Vec<McplParticleIn> = serde_wasm_bindgen::from_value(particles).map_err(js_err)?;
    if rows.len() > MAX_MCPL_PARTICLES {
        return Err(js_err(format!(
            "particle count {} exceeds the demo cap of {MAX_MCPL_PARTICLES}",
            rows.len()
        )));
    }
    let ps: Vec<nucleide_mcpl_io::Particle> = rows
        .into_iter()
        .map(|r| nucleide_mcpl_io::Particle {
            ekin: r.ekin,
            polarisation: [0.0; 3],
            position: r.position.unwrap_or([0.0; 3]),
            direction: r.direction.unwrap_or([0.0, 0.0, 1.0]),
            time: r.time,
            weight: r.weight,
            pdgcode: r.pdgcode,
            userflags: r.userflags,
        })
        .collect();
    let header = nucleide_mcpl_io::Header {
        has_userflags: h.has_userflags,
        has_polarisation: h.has_polarisation,
        double_prec: h.double_prec,
        universal_pdgcode: h.universal_pdgcode,
        universal_weight: h.universal_weight,
        srcname: h.srcname,
        comments: h.comments,
        blobs: Vec::new(),
        nparticles: ps.len() as u64,
        ..nucleide_mcpl_io::Header::default()
    };
    let bytes = nucleide_mcpl_io::encode_file(&header, &ps).map_err(js_err)?;
    let out = if h.gzip { gzip_bytes(&bytes)? } else { bytes };
    Ok(js_sys::Uint8Array::from(out.as_slice()))
}

#[derive(Deserialize, Default)]
struct Ssw2McplOptsJson {
    #[serde(default, rename = "doublePrec", alias = "double_prec")]
    double_prec: bool,
    #[serde(
        default = "default_true",
        rename = "surfToUserflags",
        alias = "surf_to_userflags"
    )]
    surf_to_userflags: bool,
    #[serde(default)]
    gzip: bool,
    #[serde(default)]
    srcname: Option<String>,
    #[serde(default)]
    comments: Vec<String>,
}

fn default_true() -> bool {
    true
}

/// Convert SSW bytes to MCPL bytes (neutron/gamma-only v1).
///
/// `sswBytes` is a `Uint8Array` of the SSW file; `surfs`/`kinds` pair each
/// track (`kinds` holds `"neutron"`/`"gamma"` only — v1 rejects other PDG
/// codes loudly). Returns MCPL file bytes as a `Uint8Array`.
#[wasm_bindgen(js_name = ssw2mcpl)]
pub fn ssw2mcpl(
    ssw_bytes: js_sys::Uint8Array,
    surfs: Vec<u32>,
    kinds: Vec<String>,
    options: JsValue,
) -> Result<js_sys::Uint8Array, JsValue> {
    let opts: Ssw2McplOptsJson = if options.is_undefined() || options.is_null() {
        Ssw2McplOptsJson::default()
    } else {
        serde_wasm_bindgen::from_value(options).map_err(js_err)?
    };
    let ssw = nucleide_mcnp_io::surfsrc::SurfSrc::from_bytes(ssw_bytes.to_vec()).map_err(js_err)?;
    let raw = ssw.read_tracklist().map_err(js_err)?;
    if raw.len() != surfs.len() || raw.len() != kinds.len() {
        return Err(js_err(format!(
            "ssw2mcpl: SSW holds {} tracks but got {} surfs and {} kinds (one surf+kind per track required)",
            raw.len(),
            surfs.len(),
            kinds.len()
        )));
    }
    let mut tracks = Vec::with_capacity(raw.len());
    for (i, ((t, surf), kind)) in raw.iter().zip(surfs).zip(kinds.iter()).enumerate() {
        let kind = nucleide_mcpl_io::ssw::SswParticleKind::parse(kind).ok_or_else(|| {
            js_err(format!(
                "track {i} kind `{kind}` unknown (expected \"neutron\" or \"gamma\")"
            ))
        })?;
        // Demo scope is neutron/gamma-only v1; other table kinds are
        // named-open here (the crate converts the full 5-kind table).
        if !matches!(
            kind,
            nucleide_mcpl_io::ssw::SswParticleKind::Neutron
                | nucleide_mcpl_io::ssw::SswParticleKind::Gamma
        ) {
            return Err(js_err(format!(
                "track {i} kind `{}` is named-open in this demo (neutron/gamma only)",
                kind.as_str()
            )));
        }
        tracks.push(nucleide_mcpl_io::ssw::SswTrack {
            ekin: t.erg,
            time_shakes: t.tme,
            position: [t.x, t.y, t.z],
            direction: [t.u, t.v, t.cs],
            weight: t.wgt,
            surf,
            kind,
        });
    }
    let options = nucleide_mcpl_io::ssw::Ssw2McplOptions {
        double_prec: opts.double_prec,
        surf_to_userflags: opts.surf_to_userflags,
        gzip: opts.gzip,
        deck_blob: None,
        srcname: opts.srcname.unwrap_or_else(|| "ssw2mcpl".to_string()),
        comments: opts.comments,
        polarisation: None,
        universal_pdg: false,
        universal_weight: false,
    };
    let bytes = nucleide_mcpl_io::ssw::ssw2mcpl_bytes(&tracks, &options).map_err(js_err)?;
    Ok(js_sys::Uint8Array::from(bytes.as_slice()))
}

/// Convert MCPL bytes back to SSW bytes against a reference SSW header.
///
/// `surface` overrides every track's surface id (`[1, 999999]`); without it
/// each particle's `userflags` supplies the id. Returns SSW file bytes.
#[wasm_bindgen(js_name = mcpl2ssw)]
pub fn mcpl2ssw(
    mcpl_bytes: js_sys::Uint8Array,
    reference_ssw_bytes: js_sys::Uint8Array,
    surface: Option<u32>,
) -> Result<js_sys::Uint8Array, JsValue> {
    let mcpl_data = maybe_gunzip_bytes(&mcpl_bytes.to_vec())?;
    let mcpl = nucleide_mcpl_io::McplFile::from_bytes(mcpl_data).map_err(js_err)?;
    let particles = mcpl.particles().map_err(js_err)?;
    if particles.len() > MAX_MCPL_PARTICLES {
        return Err(js_err(format!(
            "particle count {} exceeds the demo cap of {MAX_MCPL_PARTICLES}",
            particles.len()
        )));
    }
    let reference = nucleide_mcnp_io::surfsrc::SurfSrc::from_bytes(reference_ssw_bytes.to_vec())
        .map_err(js_err)?;
    let options = nucleide_mcpl_io::ssw::Mcpl2SswOptions {
        surface,
        ..Default::default()
    };
    let (header, tracks) =
        nucleide_mcpl_io::ssw::mcpl2ssw(&particles, &reference.header, &options).map_err(js_err)?;
    let mut out = Vec::new();
    nucleide_mcnp_io::surfsrc::write_to(&mut out, &header, &tracks).map_err(js_err)?;
    Ok(js_sys::Uint8Array::from(out.as_slice()))
}

// ---------------------------------------------------------------------------
// Material scalar bundle (separator / blender / CUSUM)
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct SeparateResult {
    product: BTreeMap<String, f64>,
    tails: BTreeMap<String, f64>,
}

fn comp_to_material(comp: JsValue) -> Result<nucleide_material::Material, JsValue> {
    comp_to_emit_material(comp)
}

fn material_to_comp(mat: &nucleide_material::Material) -> Result<BTreeMap<String, f64>, JsValue> {
    let mut out = BTreeMap::new();
    for (id, grams) in &mat.comp {
        out.insert(id.to_name(), *grams);
    }
    Ok(out)
}

/// Split one composition by per-nuclide product efficiencies.
///
/// `comp` maps GNDS names to grams; `effs` maps GNDS names to `∈ [0, 1]`.
/// Returns `{ product, tails }` (grams by nuclide name).
#[wasm_bindgen(js_name = materialSeparate)]
pub fn material_separate(comp: JsValue, effs: JsValue) -> Result<JsValue, JsValue> {
    let mat = comp_to_material(comp)?;
    let eff_map: BTreeMap<String, f64> = serde_wasm_bindgen::from_value(effs).map_err(js_err)?;
    let effs: Vec<(NuclideId, f64)> = eff_map
        .into_iter()
        .map(|(k, v)| Ok((k.parse::<NuclideId>().map_err(js_err)?, v)))
        .collect::<Result<Vec<_>, JsValue>>()?;
    let (product, tails) = mat.separate(&effs).map_err(js_err)?;
    to_js(&SeparateResult {
        product: material_to_comp(&product)?,
        tails: material_to_comp(&tails)?,
    })
}

/// Blend several compositions by fixed ratios.
///
/// `parts` is `[{ comp: {nuclide: grams}, ratio: number }]`; ratios are
/// normalized. Returns the blended composition (grams by nuclide name).
#[wasm_bindgen(js_name = materialBlend)]
pub fn material_blend(parts: JsValue) -> Result<JsValue, JsValue> {
    #[derive(Deserialize)]
    struct BlendPart {
        comp: BTreeMap<String, f64>,
        ratio: f64,
    }
    let parts: Vec<BlendPart> = serde_wasm_bindgen::from_value(parts).map_err(js_err)?;
    let mats: Vec<nucleide_material::Material> = parts
        .iter()
        .map(|p| {
            let mut mat = nucleide_material::Material::new();
            for (name, grams) in &p.comp {
                let id = name
                    .parse::<NuclideId>()
                    .map_err(|e| js_err(format!("`{name}`: {e}")))?;
                mat.add_nuclide(id, *grams);
            }
            Ok(mat)
        })
        .collect::<Result<Vec<_>, JsValue>>()?;
    let refs: Vec<(&nucleide_material::Material, f64)> =
        mats.iter().zip(parts.iter().map(|p| p.ratio)).collect();
    let out = nucleide_material::Material::blend(&refs).map_err(js_err)?;
    to_js(&material_to_comp(&out)?)
}

#[derive(Serialize)]
struct CusumResult {
    alarmed: bool,
    statistic: f64,
    mean: f64,
    std: f64,
    count: usize,
}

/// Run a one-sided upper Page CUSUM over a scalar series.
///
/// Thin facade over `nucleide-material` `Cusum` (Welford statistics):
/// feeds `series` in order with the crate defaults (`k = 0.5`, `h = 4.0`,
/// `startup = 10`) unless overridden. Returns
/// `{ alarmed, statistic, mean, std, count }`.
#[wasm_bindgen(js_name = cusumDetect)]
pub fn cusum_detect(
    series: Vec<f64>,
    k: Option<f64>,
    h: Option<f64>,
    startup: Option<usize>,
) -> Result<JsValue, JsValue> {
    let mut cusum =
        nucleide_material::Cusum::new(k.unwrap_or(0.5), h.unwrap_or(4.0), startup.unwrap_or(10))
            .map_err(js_err)?;
    let mut alarmed = false;
    for x in series {
        alarmed = cusum.update(x);
    }
    to_js(&CusumResult {
        alarmed,
        statistic: cusum.statistic(),
        mean: cusum.mean(),
        std: cusum.std(),
        count: cusum.count(),
    })
}

// ---------------------------------------------------------------------------
// Spectroscopy scalar bundle (TSV / calib / SPE readers)
// ---------------------------------------------------------------------------

/// Parse interchange-TSV decay lines into `[[energyMeV, intensity]]`.
///
/// Thin facade over `parse_lines_tsv` (comments/`#`, blanks skipped).
#[wasm_bindgen(js_name = parseLinesTsv)]
pub fn parse_lines_tsv(text: &str) -> Result<JsValue, JsValue> {
    let rows = nucleide_spectroscopy::parse_lines_tsv(text).map_err(js_err)?;
    to_js(&rows)
}

/// Quadratic energy bins (E6): `ebin[ch] = a0 + a1*ch + a2*ch^2`.
#[wasm_bindgen(js_name = energyBins)]
pub fn energy_bins(channels: Vec<f64>, fit: Vec<f64>) -> Result<JsValue, JsValue> {
    to_js(&nucleide_spectroscopy::energy_bins(&channels, &fit).map_err(js_err)?)
}

/// Detector efficiency at `energyMev` (E7, `effFit` 1 or 2).
#[wasm_bindgen(js_name = detectorEfficiency)]
pub fn detector_efficiency(energy_mev: f64, coeff: Vec<f64>, eff_fit: i32) -> Result<f64, JsValue> {
    nucleide_spectroscopy::detector_efficiency(energy_mev, &coeff, i64::from(eff_fit))
        .map_err(js_err)
}

#[derive(Serialize)]
struct SpeSummary {
    spec_name: String,
    channels: usize,
    #[serde(rename = "startChan")]
    start_chan: i64,
    #[serde(rename = "liveTime")]
    live_time: f64,
    #[serde(rename = "realTime")]
    real_time: f64,
    #[serde(rename = "detId")]
    det_id: String,
    #[serde(rename = "energyFit")]
    energy_fit: Vec<f64>,
    counts: Vec<f64>,
    ebins: Vec<f64>,
    truncated: bool,
}

fn spe_summary_of(spec: &nucleide_spectroscopy::GammaSpectrum) -> SpeSummary {
    let n = spec.spectrum.counts.len();
    let truncated = n > MAX_MCPL_PARTICLES;
    SpeSummary {
        spec_name: spec.spectrum.spec_name.clone(),
        channels: n,
        start_chan: spec.spectrum.start_chan_num,
        live_time: spec.live_time,
        real_time: spec.real_time,
        det_id: spec.det_id.clone(),
        energy_fit: spec.calib_e_fit.clone(),
        counts: spec
            .spectrum
            .counts
            .iter()
            .take(MAX_MCPL_PARTICLES)
            .copied()
            .collect(),
        ebins: spec
            .spectrum
            .ebin
            .iter()
            .take(MAX_MCPL_PARTICLES)
            .copied()
            .collect(),
        truncated,
    }
}

/// Parse a dollar-format `.spe` file into a capped JSON summary.
///
/// `GammaSpectrum` has no `Serialize`; this follows the `IsotxsSummary`
/// precedent (summary struct, capped lists).
#[wasm_bindgen(js_name = parseDollarSpe)]
pub fn parse_dollar_spe(text: &str) -> Result<JsValue, JsValue> {
    let spec = nucleide_spectroscopy::parse_dollar_spe(text, "").map_err(js_err)?;
    to_js(&spe_summary_of(&spec))
}

/// Parse a plain-format `.spe` file into a capped JSON summary.
#[wasm_bindgen(js_name = parsePlainSpe)]
pub fn parse_plain_spe(text: &str) -> Result<JsValue, JsValue> {
    let spec = nucleide_spectroscopy::parse_plain_spe(text, "").map_err(js_err)?;
    to_js(&spe_summary_of(&spec))
}

// ---------------------------------------------------------------------------
// Kinetics scalar leftovers (inhour / stable period / prompt jump)
// ---------------------------------------------------------------------------

fn kinetic_params(
    betas: Vec<f64>,
    lambdas: Vec<f64>,
    lambda_gen: f64,
) -> Result<nucleide_kinetics::KineticParams, JsValue> {
    nucleide_kinetics::KineticParams::new(betas, lambdas, lambda_gen).map_err(js_err)
}

/// Inhour right-hand side `rho(omega)` (E3) [Δk].
#[wasm_bindgen(js_name = inhourRho)]
pub fn inhour_rho(
    betas: Vec<f64>,
    lambdas: Vec<f64>,
    lambda_gen: f64,
    omega: f64,
) -> Result<f64, JsValue> {
    let params = kinetic_params(betas, lambdas, lambda_gen)?;
    nucleide_kinetics::rho_of_omega(&params, omega).map_err(js_err)
}

/// Asymptotic stable period `T = 1/omega` [s] for `0 < rho < beta`.
#[wasm_bindgen(js_name = stablePeriod)]
pub fn stable_period(
    betas: Vec<f64>,
    lambdas: Vec<f64>,
    lambda_gen: f64,
    rho: f64,
) -> Result<f64, JsValue> {
    let params = kinetic_params(betas, lambdas, lambda_gen)?;
    nucleide_kinetics::stable_period(&params, rho).map_err(js_err)
}

/// Prompt-jump factor `n_after = n_before * beta / (beta - rho_after)`.
///
/// Thin facade over the E4 formula; errors at/past prompt critical.
#[wasm_bindgen(js_name = promptJump)]
pub fn prompt_jump(
    n_before: f64,
    rho_before: f64,
    rho_after: f64,
    beta_total: f64,
) -> Result<f64, JsValue> {
    nucleide_kinetics::prompt_jump(n_before, rho_before, rho_after, beta_total).map_err(js_err)
}

// ---------------------------------------------------------------------------
// Deterministic scalar (RTFLUX; PARTISN stays RECORD)
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct RtfluxSummary {
    kind: String,
    groups: usize,
    npoints: usize,
    values: Vec<f64>,
    truncated: bool,
}

/// Parse an RTFLUX/ATFLUX/RZFLUX flux file into a capped JSON summary.
///
/// Mirrors the Python `kind`-switch (`rtflux`|`atflux`|`rzflux`); PARTISN
/// deck writing stays RECORD (Python-only).
#[wasm_bindgen(js_name = parseRtflux)]
pub fn parse_rtflux(text: &str, kind: &str) -> Result<JsValue, JsValue> {
    let want = match kind.to_ascii_lowercase().as_str() {
        "rtflux" => nucleide_cccc_io::rtflux::FluxKind::Rtflux,
        "atflux" => nucleide_cccc_io::rtflux::FluxKind::Atflux,
        "rzflux" => nucleide_cccc_io::rtflux::FluxKind::Rzflux,
        other => {
            return Err(js_err(format!(
                "kind must be rtflux|atflux|rzflux, got `{other}`"
            )));
        }
    };
    let flux = nucleide_cccc_io::FluxFile::parse(want, text).map_err(js_err)?;
    let npoints = flux.values.len() / flux.groups.max(1);
    let truncated = flux.values.len() > MAX_MCPL_PARTICLES;
    to_js(&RtfluxSummary {
        kind: flux.kind.keyword().to_string(),
        groups: flux.groups,
        npoints,
        values: flux.values.into_iter().take(MAX_MCPL_PARTICLES).collect(),
        truncated,
    })
}

// ---------------------------------------------------------------------------
// Damage (NRT/arc-dpa folds + closed-form point evaluations)
// ---------------------------------------------------------------------------

/// NRT-dpa from one spectral fold: `seconds · 1e-24 · Σ_g flux[g]·response[g]`
/// with `flux` the per-group integrated flux [n/cm²/s] and `response` the
/// group dpa cross sections [barns] over `bounds` (`G + 1` MeV boundaries).
#[wasm_bindgen(js_name = damageNrtDpa)]
pub fn damage_nrt_dpa(
    flux: Vec<f64>,
    response: Vec<f64>,
    bounds: Vec<f64>,
    seconds: f64,
) -> Result<f64, JsValue> {
    nucleide_damage::nrt_dpa(&flux, &response, &bounds, seconds).map_err(js_err)
}

/// arc-dpa: same fold with arc-corrected dpa cross sections.
#[wasm_bindgen(js_name = damageArcDpa)]
pub fn damage_arc_dpa(
    flux: Vec<f64>,
    response: Vec<f64>,
    bounds: Vec<f64>,
    seconds: f64,
) -> Result<f64, JsValue> {
    nucleide_damage::arc_dpa(&flux, &response, &bounds, seconds).map_err(js_err)
}

/// Gas production [appm]: `seconds · 1e-18 · Σ_g flux[g]·response[g]` with the
/// caller's gas-production cross sections [barns].
#[wasm_bindgen(js_name = damageGasAppm)]
pub fn damage_gas_appm(
    flux: Vec<f64>,
    response: Vec<f64>,
    bounds: Vec<f64>,
    seconds: f64,
) -> Result<f64, JsValue> {
    nucleide_damage::gas_appm(&flux, &response, &bounds, seconds).map_err(js_err)
}

/// He/dpa ratio [appm per dpa] from one flux fold of the He production and
/// damage cross sections; a zero damage fold is a loud error, never `inf`.
#[wasm_bindgen(js_name = damageHeDpaRatio)]
pub fn damage_he_dpa_ratio(
    flux: Vec<f64>,
    he_response: Vec<f64>,
    damage_response: Vec<f64>,
    bounds: Vec<f64>,
    seconds: f64,
) -> Result<f64, JsValue> {
    nucleide_damage::he_dpa_ratio(&flux, &he_response, &damage_response, &bounds, seconds)
        .map_err(js_err)
}

/// Lindhard damage energy `T_dam = T·P(ε)` [eV] for a self-recoil `target`
/// (Robinson fit; energies in eV throughout).
#[wasm_bindgen(js_name = damageEnergy)]
pub fn damage_energy(t_ev: f64, target: &str) -> Result<f64, JsValue> {
    let id = target.parse::<NuclideId>().map_err(js_err)?;
    nucleide_damage::damage_energy(t_ev, &id, &id).map_err(js_err)
}

/// NRT displacement count `N_d(T)` for a self-recoil `target` with threshold
/// displacement energy `ed_ev` [eV].
#[wasm_bindgen(js_name = nrtDisplacements)]
pub fn nrt_displacements(t_ev: f64, ed_ev: f64, target: &str) -> Result<f64, JsValue> {
    let id = target.parse::<NuclideId>().map_err(js_err)?;
    nucleide_damage::nrt_displacements(t_ev, ed_ev, &id).map_err(js_err)
}

/// arc-dpa efficiency `ξ(T_dam)` (Nordlund et al. 2018 Eq. (7)) at damage
/// energy `t_dam_ev` [eV] with material constants `b_arc` (`< 0`) and `c_arc`
/// (`∈ (0, 1)`).
#[wasm_bindgen(js_name = arcEfficiency)]
pub fn arc_efficiency(t_dam_ev: f64, ed_ev: f64, b_arc: f64, c_arc: f64) -> Result<f64, JsValue> {
    let params = nucleide_damage::ArcParams::new(b_arc, c_arc).map_err(js_err)?;
    nucleide_damage::arc_efficiency(t_dam_ev, ed_ev, &params).map_err(js_err)
}

// ---------------------------------------------------------------------------
// Fusion neutron sources (ring/point/parametric: moments, sampling, cards)
// ---------------------------------------------------------------------------

/// Demo cap on sampled source particles: the browser slice stays cheap.
const MAX_FUSION_SAMPLES: usize = 5000;

fn parse_fusion_reaction(s: &str) -> Result<nucleide_plasma_source::FusionReaction, JsValue> {
    match s.to_ascii_lowercase().replace([' ', '-', '_'], "").as_str() {
        "dt" | "td" => Ok(nucleide_plasma_source::FusionReaction::Dt),
        "dd" => Ok(nucleide_plasma_source::FusionReaction::Dd),
        _ => Err(js_err(format!(
            "unknown fusion reaction `{s}` (supported: dt, dd)"
        ))),
    }
}

fn check_seed(seed: f64) -> Result<u64, JsValue> {
    if !seed.is_finite() || seed < 0.0 || seed.fract() != 0.0 {
        return Err(js_err(format!(
            "seed must be a finite non-negative integer (got {seed})"
        )));
    }
    if seed >= 9_007_199_254_740_992.0 {
        return Err(js_err(format!(
            "seed must be below 2^53 for exact f64 integer precision (got {seed})"
        )));
    }
    Ok(seed as u64)
}

fn check_sample_count(n: usize) -> Result<(), JsValue> {
    if n == 0 {
        return Err(js_err("n must be >= 1"));
    }
    if n > MAX_FUSION_SAMPLES {
        return Err(js_err(format!(
            "n = {n} exceeds the demo cap of {MAX_FUSION_SAMPLES} particles"
        )));
    }
    Ok(())
}

#[derive(Serialize)]
struct FusionSpectrumMoments {
    reaction: String,
    #[serde(rename = "nominalMeV")]
    nominal_mev: f64,
    #[serde(rename = "meanMeV")]
    mean_mev: f64,
    #[serde(rename = "sigmaMeV")]
    sigma_mev: f64,
    mono: bool,
}

/// Ballabio spectrum moments of `reaction` (`"dt"` / `"dd"`) at ion
/// temperature `ti_kev` [keV]: the nominal line, the shifted Gaussian mean,
/// and the Brysk width (`sigma` zero at `T_i = 0`, the monoenergetic line).
#[wasm_bindgen(js_name = fusionSpectrumMoments)]
pub fn fusion_spectrum_moments(reaction: &str, ti_kev: f64) -> Result<JsValue, JsValue> {
    let reaction = parse_fusion_reaction(reaction)?;
    let (mean_mev, sigma_mev) = reaction.moments_mev(ti_kev).map_err(js_err)?;
    to_js(&FusionSpectrumMoments {
        reaction: reaction.label().to_string(),
        nominal_mev: reaction.nominal_energy_mev(),
        mean_mev,
        sigma_mev,
        mono: sigma_mev == 0.0,
    })
}

/// Thermonuclear reactivity ⟨σv⟩ [m³/s] of `reaction` (`"dt"` / `"dd"`) at ion
/// temperature `ti_kev` [keV] (Bosch & Hale 1992 fit; zero at `T_i = 0`).
#[wasm_bindgen(js_name = fusionReactivity)]
pub fn fusion_reactivity(reaction: &str, ti_kev: f64) -> Result<f64, JsValue> {
    parse_fusion_reaction(reaction)?
        .reactivity_m3_per_s(ti_kev)
        .map_err(js_err)
}

#[derive(Deserialize)]
struct FusionRingSpecJson {
    #[serde(rename = "radiusCm", alias = "radius_cm")]
    radius_cm: f64,
    #[serde(rename = "heightCm", alias = "height_cm")]
    height_cm: f64,
    reaction: String,
    #[serde(rename = "tiKev", alias = "ti_kev")]
    ti_kev: f64,
}

#[derive(Deserialize)]
struct FusionPointSpecJson {
    #[serde(rename = "xCm", alias = "x_cm")]
    x_cm: f64,
    #[serde(rename = "yCm", alias = "y_cm")]
    y_cm: f64,
    #[serde(rename = "zCm", alias = "z_cm")]
    z_cm: f64,
    reaction: String,
    #[serde(rename = "tiKev", alias = "ti_kev")]
    ti_kev: f64,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct FusionParametricSpecJson {
    // Unknown keys are rejected loudly: silently dropping a `fuel_mixture`,
    // `species_temperatures`, or `sector` key would return full-torus
    // single-fuel results for a caller that asked for something else.
    #[serde(rename = "majorRadiusCm", alias = "major_radius_cm")]
    major_radius_cm: f64,
    #[serde(rename = "minorRadiusCm", alias = "minor_radius_cm")]
    minor_radius_cm: f64,
    elongation: f64,
    triangularity: f64,
    #[serde(rename = "shafranovFactorCm", alias = "shafranov_factor_cm")]
    shafranov_factor_cm: f64,
    mode: String,
    fuel: String,
    #[serde(rename = "centreDensityM3", alias = "centre_density_m3")]
    centre_density_m3: f64,
    #[serde(rename = "densityPeaking", alias = "density_peaking")]
    density_peaking: f64,
    #[serde(rename = "pedestalDensityM3", alias = "pedestal_density_m3")]
    pedestal_density_m3: f64,
    #[serde(rename = "separatrixDensityM3", alias = "separatrix_density_m3")]
    separatrix_density_m3: f64,
    #[serde(rename = "centreTempKev", alias = "centre_temp_kev")]
    centre_temp_kev: f64,
    #[serde(rename = "tempPeaking", alias = "temp_peaking")]
    temp_peaking: f64,
    #[serde(rename = "tempBeta", alias = "temp_beta")]
    temp_beta: f64,
    #[serde(rename = "pedestalTempKev", alias = "pedestal_temp_kev")]
    pedestal_temp_kev: f64,
    #[serde(rename = "separatrixTempKev", alias = "separatrix_temp_kev")]
    separatrix_temp_kev: f64,
    #[serde(rename = "pedestalRadiusCm", alias = "pedestal_radius_cm")]
    pedestal_radius_cm: f64,
    #[serde(rename = "fuelDeuterium", alias = "fuel_deuterium")]
    fuel_deuterium: Option<f64>,
    #[serde(rename = "fuelTritium", alias = "fuel_tritium")]
    fuel_tritium: Option<f64>,
    #[serde(rename = "tailFraction", alias = "tail_fraction")]
    tail_fraction: Option<f64>,
    #[serde(rename = "tailTempKev", alias = "tail_temp_kev")]
    tail_temp_kev: Option<f64>,
    #[serde(rename = "speciesDeuteriumKev", alias = "species_deuterium_kev")]
    species_deuterium_kev: Option<f64>,
    #[serde(rename = "speciesTritiumKev", alias = "species_tritium_kev")]
    species_tritium_kev: Option<f64>,
}

/// Read one string field off a raw JS object.
fn js_get_string(obj: &JsValue, key: &str) -> Result<String, JsValue> {
    js_sys::Reflect::get(obj, &JsValue::from_str(key))
        .ok()
        .and_then(|v| v.as_string())
        .ok_or_else(|| js_err(format!("spec is missing string field `{key}`")))
}

/// Read one finite-number field off a raw JS object.
fn js_get_f64(obj: &JsValue, key: &str) -> Result<f64, JsValue> {
    let v = js_sys::Reflect::get(obj, &JsValue::from_str(key))
        .ok()
        .and_then(|v| v.as_f64())
        .ok_or_else(|| js_err(format!("spec is missing numeric field `{key}`")))?;
    if !v.is_finite() {
        return Err(js_err(format!(
            "spec field `{key}` must be finite (got {v})"
        )));
    }
    Ok(v)
}

fn ring_config(spec: &JsValue) -> Result<nucleide_plasma_source::PlasmaSourceConfig, JsValue> {
    let parsed: FusionRingSpecJson =
        serde_wasm_bindgen::from_value(spec.clone()).map_err(js_err)?;
    Ok(nucleide_plasma_source::PlasmaSourceConfig::ring(
        parsed.radius_cm,
        parsed.height_cm,
        parse_fusion_reaction(&parsed.reaction)?,
        parsed.ti_kev,
    ))
}

fn point_config(spec: &JsValue) -> Result<nucleide_plasma_source::PlasmaSourceConfig, JsValue> {
    let parsed: FusionPointSpecJson =
        serde_wasm_bindgen::from_value(spec.clone()).map_err(js_err)?;
    Ok(nucleide_plasma_source::PlasmaSourceConfig::point(
        parsed.x_cm,
        parsed.y_cm,
        parsed.z_cm,
        parse_fusion_reaction(&parsed.reaction)?,
        parsed.ti_kev,
    ))
}

fn parametric_config(
    spec: &JsValue,
) -> Result<nucleide_plasma_source::ParametricPlasmaConfig, JsValue> {
    // The sampling envelope (`kind` dispatch, `n`/`seed` counts) rides on the
    // same object `sample_fusion_source` matches on, but the spec struct
    // denies unknown fields (so a stray `fuel_mixture`/`species_temperatures`/
    // `sector` key fails loudly instead of silently sampling the wrong
    // plasma). Strip just the envelope keys before struct parsing.
    let body: js_sys::Object = spec.clone().into();
    for key in ["kind", "n", "seed"] {
        let _ = js_sys::Reflect::delete_property(&body, &JsValue::from_str(key));
    }
    let parsed: FusionParametricSpecJson =
        serde_wasm_bindgen::from_value(JsValue::from(body)).map_err(js_err)?;
    use nucleide_plasma_source::{
        DensityProfile, MillerGeometry, ParametricPlasmaConfig, ProfileMode, TemperatureProfile,
    };
    let mode = ProfileMode::parse(&parsed.mode).map_err(js_err)?;
    let fuel_mixture = match (parsed.fuel_deuterium, parsed.fuel_tritium) {
        (Some(d), Some(t)) => Some(nucleide_plasma_source::FuelMixture::new(d, t).map_err(js_err)?),
        (None, None) => None,
        _ => {
            return Err(js_err(
                "fuelDeuterium and fuelTritium must be given together",
            ));
        }
    };
    let tail = match (parsed.tail_fraction, parsed.tail_temp_kev) {
        (Some(f), Some(t)) => {
            if fuel_mixture.is_none() {
                return Err(js_err(
                    "tailFraction/tailTempKev need a D/T fuel mixture (fuelDeuterium + fuelTritium)",
                ));
            }
            Some(nucleide_plasma_source::DeuteriumTail::new(f, t).map_err(js_err)?)
        }
        (None, None) => None,
        _ => {
            return Err(js_err(
                "tailFraction and tailTempKev must be given together",
            ));
        }
    };
    // Per-species ion temperatures compose with mixtures and with
    // single-fuel configs (single-fuel D-T reacts at the mass-weighted
    // T_DT, single-fuel D-D at T_D); the pair entries validate loudly in
    // the core.
    let species_temperatures = match (parsed.species_deuterium_kev, parsed.species_tritium_kev) {
        (Some(d), Some(t)) => {
            Some(nucleide_plasma_source::SpeciesIonTemperatures::new(d, t).map_err(js_err)?)
        }
        (None, None) => None,
        _ => {
            return Err(js_err(
                "speciesDeuteriumKev and speciesTritiumKev must be given together",
            ));
        }
    };
    Ok(ParametricPlasmaConfig {
        geometry: MillerGeometry {
            major_radius_cm: parsed.major_radius_cm,
            minor_radius_cm: parsed.minor_radius_cm,
            elongation: parsed.elongation,
            triangularity: parsed.triangularity,
            shafranov_factor_cm: parsed.shafranov_factor_cm,
        },
        mode,
        ion_density: DensityProfile {
            centre_m3: parsed.centre_density_m3,
            peaking_factor: parsed.density_peaking,
            pedestal_m3: parsed.pedestal_density_m3,
            separatrix_m3: parsed.separatrix_density_m3,
        },
        ion_temperature: TemperatureProfile {
            centre_kev: parsed.centre_temp_kev,
            peaking_factor: parsed.temp_peaking,
            beta: parsed.temp_beta,
            pedestal_kev: parsed.pedestal_temp_kev,
            separatrix_kev: parsed.separatrix_temp_kev,
        },
        pedestal_radius_cm: parsed.pedestal_radius_cm,
        fuel: parse_fusion_reaction(&parsed.fuel)?,
        fuel_mixture,
        species_temperatures,
        tail,
        sector: None,
        weight: 1.0,
    })
}

#[derive(Serialize)]
struct FusionParticleJson {
    #[serde(rename = "positionCm")]
    position_cm: [f64; 3],
    direction: [f64; 3],
    #[serde(rename = "energyMeV")]
    energy_mev: f64,
    weight: f64,
}

#[derive(Serialize)]
struct FusionSampleResult {
    kind: String,
    count: usize,
    particles: Vec<FusionParticleJson>,
}

fn particle_json(particles: &[nucleide_plasma_source::Particle]) -> Vec<FusionParticleJson> {
    particles
        .iter()
        .map(|p| FusionParticleJson {
            position_cm: p.position_cm,
            direction: p.direction,
            energy_mev: p.energy_mev,
            weight: p.weight,
        })
        .collect()
}

/// Sample `n` seeded fusion-source particles.
///
/// `spec` is a discriminated union on `kind`: `"ring"` takes `radiusCm`,
/// `heightCm`, `reaction` (`"dt"`/`"dd"`), `tiKev`; `"point"` takes `xCm`,
/// `yCm`, `zCm` instead of the ring geometry; `"parametric"` takes the Miller
/// geometry (`majorRadiusCm`, `minorRadiusCm`, `elongation`, `triangularity`,
/// `shafranovFactorCm`), the confinement `mode` (`"L"`/`"H"`/`"A"`), the
/// `fuel`, and the Fausser profile parameters (`centreDensityM3`,
/// `densityPeaking`, `pedestalDensityM3`, `separatrixDensityM3`,
/// `centreTempKev`, `tempPeaking`, `tempBeta`, `pedestalTempKev`,
/// `separatrixTempKev`, `pedestalRadiusCm`). All kinds take `n` (capped at
/// [`MAX_FUSION_SAMPLES`]) and an integer `seed`. Returns `{kind, count,
/// particles}` with one `{positionCm, direction, energyMeV, weight}` row per
/// particle (lengths in cm, energies in MeV).
#[wasm_bindgen(js_name = sampleFusionSource)]
pub fn sample_fusion_source(spec: JsValue) -> Result<JsValue, JsValue> {
    let kind = js_get_string(&spec, "kind")?.to_ascii_lowercase();
    let (n, seed) = sample_args(&spec)?;
    let result = match kind.as_str() {
        "ring" => {
            let mut sampler = nucleide_plasma_source::SourceSampler::new(ring_config(&spec)?, seed)
                .map_err(js_err)?;
            let particles = sampler.sample_n(n);
            FusionSampleResult {
                kind: "ring".to_string(),
                count: particles.len(),
                particles: particle_json(&particles),
            }
        }
        "point" => {
            let mut sampler =
                nucleide_plasma_source::SourceSampler::new(point_config(&spec)?, seed)
                    .map_err(js_err)?;
            let particles = sampler.sample_n(n);
            FusionSampleResult {
                kind: "point".to_string(),
                count: particles.len(),
                particles: particle_json(&particles),
            }
        }
        "parametric" => {
            let mut sampler =
                nucleide_plasma_source::ParametricSampler::new(parametric_config(&spec)?, seed)
                    .map_err(js_err)?;
            let particles = sampler.sample_n(n);
            FusionSampleResult {
                kind: "parametric".to_string(),
                count: particles.len(),
                particles: particle_json(&particles),
            }
        }
        other => {
            return Err(js_err(format!(
                "unknown fusion source kind `{other}` (supported: ring, point, parametric)"
            )));
        }
    };
    to_js(&result)
}

/// Read the shared `n`/`seed` scalars off the raw spec object.
fn sample_args(spec: &JsValue) -> Result<(usize, u64), JsValue> {
    let n = js_get_f64(spec, "n")?;
    check_sample_count(n as usize)?;
    let seed = check_seed(js_get_f64(spec, "seed")?)?;
    Ok((n as usize, seed))
}

#[derive(Serialize)]
struct FusionCardsResult {
    mcnp: String,
    serpent: String,
}

/// Emit MCNP `SDEF` and Serpent `src` cards for one fusion source.
///
/// `spec` matches [`sample_fusion_source`] (the `n`/`seed` fields are
/// ignored); optional `nBins` (default 21) sets the Gaussian tabulation bin
/// count and `mcnpVersion` (default 6) selects the `PAR=` designator dialect
/// for the `SDEF` card. Returns `{mcnp, serpent}` card texts.
#[wasm_bindgen(js_name = emitFusionSourceCards)]
pub fn emit_fusion_source_cards(spec: JsValue) -> Result<JsValue, JsValue> {
    let get_opt = |key: &str| -> Option<f64> {
        js_sys::Reflect::get(&spec, &JsValue::from_str(key))
            .ok()
            .and_then(|v| v.as_f64())
    };
    let n_bins = get_opt("nBins").unwrap_or(21.0) as usize;
    let version = get_opt("mcnpVersion").unwrap_or(6.0) as u32;
    let kind = js_get_string(&spec, "kind")?.to_ascii_lowercase();
    let (mcnp, serpent) = match kind.as_str() {
        "ring" | "point" => {
            let config = if kind == "ring" {
                ring_config(&spec)?
            } else {
                point_config(&spec)?
            };
            (
                nucleide_plasma_source::emit_sdef(&config, version, n_bins).map_err(js_err)?,
                nucleide_plasma_source::emit_serpent(&config, n_bins).map_err(js_err)?,
            )
        }
        "parametric" => {
            let config = parametric_config(&spec)?;
            (
                nucleide_plasma_source::emit_sdef_parametric(&config, version, n_bins)
                    .map_err(js_err)?,
                nucleide_plasma_source::emit_serpent_parametric(&config, n_bins).map_err(js_err)?,
            )
        }
        other => {
            return Err(js_err(format!(
                "unknown fusion source kind `{other}` (supported: ring, point, parametric)"
            )));
        }
    };
    to_js(&FusionCardsResult {
        mcnp: mcnp.text,
        serpent: serpent.text,
    })
}

// ---------------------------------------------------------------------------
// Spectrum unfolding (forward fold + SAND-II solve)
// ---------------------------------------------------------------------------

/// Forward operator: fold `spectrum` (one value per energy group) through the
/// `response` matrix (one row per detector) into calculated rates per
/// detector.
#[wasm_bindgen(js_name = unfoldForwardFold)]
pub fn unfold_forward_fold(response: JsValue, spectrum: Vec<f64>) -> Result<JsValue, JsValue> {
    let response: Vec<Vec<f64>> = serde_wasm_bindgen::from_value(response).map_err(js_err)?;
    to_js(&nucleide_unfold::forward_fold(&response, &spectrum).map_err(js_err)?)
}

#[derive(Serialize)]
struct SandiiSolutionJson {
    spectrum: Vec<f64>,
    rates: Vec<f64>,
    #[serde(rename = "rateFactors")]
    rate_factors: Vec<f64>,
    iterations: usize,
    tolerance: f64,
    #[serde(rename = "maxRelChange")]
    max_rel_change: f64,
}

/// Run the SAND-II iterative adjustment to convergence.
///
/// `response` is a nested array (one row per detector, one value per energy
/// group), `rates` the measured rate per detector, and `guess` one strictly
/// positive starting value per group. `tolerance` (default 1e-3) is the
/// per-group relative-change convergence bound and `maxIterations` (default
/// 200) the adjustment cap — exhausting it is a loud error, never a partial
/// spectrum. Returns `{spectrum, rates, rateFactors, iterations, tolerance,
/// maxRelChange}`.
#[wasm_bindgen(js_name = sandiiSolve)]
pub fn sandii_solve(
    response: JsValue,
    rates: Vec<f64>,
    guess: Vec<f64>,
    tolerance: Option<f64>,
    max_iterations: Option<f64>,
) -> Result<JsValue, JsValue> {
    let response: Vec<Vec<f64>> = serde_wasm_bindgen::from_value(response).map_err(js_err)?;
    let tolerance = tolerance.unwrap_or(nucleide_unfold::DEFAULT_TOLERANCE);
    let max_iterations = max_iterations
        .unwrap_or(nucleide_unfold::DEFAULT_MAX_ITERATIONS as f64)
        .round() as usize;
    let solution =
        nucleide_unfold::sandii::unfold(&response, &rates, &guess, tolerance, max_iterations)
            .map_err(js_err)?;
    to_js(&SandiiSolutionJson {
        spectrum: solution.spectrum,
        rates: solution.rates,
        rate_factors: solution.rate_factors,
        iterations: solution.iterations,
        tolerance: solution.tolerance,
        max_rel_change: solution.max_rel_change,
    })
}

// ---------------------------------------------------------------------------
// Clearance screening (EU 2013/59/Euratom Annex VII Table A default)
// ---------------------------------------------------------------------------

/// The vendored EU 2013/59/Euratom Annex VII Table A default: activity
/// concentrations for solid material in Bq/g. Returns `{count, entries}` with
/// one `{nuclide, limitBqG}` row per table entry in canonical nucid order.
#[wasm_bindgen(js_name = euClearanceTable)]
pub fn eu_clearance_table() -> Result<JsValue, JsValue> {
    #[derive(Serialize)]
    struct ClearanceEntryJson {
        nuclide: String,
        #[serde(rename = "limitBqG")]
        limit_bq_g: f64,
    }
    let table = nucleide_alara_io::ClearanceTable::eu_annex_vii();
    let entries: Vec<ClearanceEntryJson> = table
        .iter()
        .map(|(id, limit)| ClearanceEntryJson {
            nuclide: id.to_name(),
            limit_bq_g: limit,
        })
        .collect();
    #[derive(Serialize)]
    struct ClearanceTableJson {
        count: usize,
        entries: Vec<ClearanceEntryJson>,
    }
    to_js(&ClearanceTableJson {
        count: entries.len(),
        entries,
    })
}

fn clearance_inventory(comp: JsValue) -> Result<Vec<(NuclideId, f64)>, JsValue> {
    let map: BTreeMap<String, f64> = serde_wasm_bindgen::from_value(comp).map_err(js_err)?;
    map.into_iter()
        .map(|(name, activity)| {
            let id = name
                .parse::<NuclideId>()
                .map_err(|e| js_err(format!("`{name}`: {e}")))?;
            Ok((id, activity))
        })
        .collect()
}

/// Clearance index `CI = Σ_i A_i / CL_i` of a `{nuclide: Bq/g}` inventory
/// against the default EU Annex VII Table A (every inventory nuclide needs a
/// table entry; activities must be finite and non-negative).
#[wasm_bindgen(js_name = clearanceIndex)]
pub fn clearance_index(inventory: JsValue) -> Result<f64, JsValue> {
    let inventory = clearance_inventory(inventory)?;
    let table = nucleide_alara_io::ClearanceTable::eu_annex_vii();
    nucleide_alara_io::clearance::clearance_index(&inventory, &table).map_err(js_err)
}

/// Sum-of-fractions screening of a `{nuclide: Bq/g}` inventory against the
/// default EU Annex VII Table A: `<= 1` satisfies the screening criterion
/// (boundary included). Returns `{sum, class: "satisfied" | "exceeded",
/// maxFraction, maxNuclide}` with the dominant contributor (`maxNuclide` is
/// `null` for an empty inventory).
#[wasm_bindgen(js_name = clearanceSumOfFractions)]
pub fn clearance_sum_of_fractions(inventory: JsValue) -> Result<JsValue, JsValue> {
    #[derive(Serialize)]
    struct SumOfFractionsJson {
        sum: f64,
        class: String,
        #[serde(rename = "maxFraction")]
        max_fraction: f64,
        #[serde(rename = "maxNuclide")]
        max_nuclide: Option<String>,
    }
    let inventory = clearance_inventory(inventory)?;
    let table = nucleide_alara_io::ClearanceTable::eu_annex_vii();
    let out = nucleide_alara_io::clearance::sum_of_fractions(&inventory, &table).map_err(js_err)?;
    to_js(&SumOfFractionsJson {
        sum: out.sum,
        class: out.class.to_string(),
        max_fraction: out.max_fraction,
        max_nuclide: out.max_nuclide.map(|id| id.to_name()),
    })
}

// ---------------------------------------------------------------------------
// Heating, lattice, blanket, coil, sublet, and equilibrium facades
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct EcrhScalarsJson {
    #[serde(rename = "gyrofrequencyGhz")]
    gyrofrequency_ghz: f64,
    #[serde(rename = "resonantFieldT")]
    resonant_field_t: f64,
    #[serde(rename = "relativisticFieldT")]
    relativistic_field_t: Option<f64>,
    #[serde(rename = "o1CutoffDensityM3")]
    o1_cutoff_density_m3: f64,
    #[serde(rename = "x1CutoffDensityM3")]
    x1_cutoff_density_m3: Option<f64>,
}

fn check_harmonic(harmonic: f64) -> Result<u32, JsValue> {
    if !harmonic.is_finite() || harmonic.fract() != 0.0 || harmonic < 1.0 || harmonic > 10.0 {
        return Err(js_err(format!(
            "harmonic must be an integer in [1, 10] (got {harmonic})"
        )));
    }
    Ok(harmonic as u32)
}

/// Closed-form ECRH scalar quantities (ECRH-1–ECRH-5).
///
/// `frequencyGhz` is the wave frequency in GHz, `harmonic` the cyclotron
/// harmonic (integer 1–10), `bT` the local field in tesla, and `teKev` the
/// optional electron temperature in keV (relativistic shift). Returns
/// `{gyrofrequencyGhz, resonantFieldT, relativisticFieldT,
/// o1CutoffDensityM3, x1CutoffDensityM3}` (densities in m⁻³; the X1 cut-off
/// is `null` where the wave is evanescent, `f ≤ f_ce`).
#[wasm_bindgen(js_name = ecrhScalars)]
pub fn ecrh_scalars(
    frequency_ghz: f64,
    harmonic: f64,
    b_t: f64,
    te_kev: Option<f64>,
) -> Result<JsValue, JsValue> {
    use nucleide_plasma_source::ecrh;
    let n = check_harmonic(harmonic)?;
    let cold = ecrh::cold_resonant_field_t(frequency_ghz, n).map_err(js_err)?;
    let rel = te_kev
        .map(ecrh::relativistic_gamma)
        .transpose()
        .map_err(js_err)?
        .map(|g| cold * g);
    to_js(&EcrhScalarsJson {
        gyrofrequency_ghz: ecrh::gyrofrequency_ghz(b_t).map_err(js_err)?,
        resonant_field_t: cold,
        relativistic_field_t: rel,
        o1_cutoff_density_m3: ecrh::o1_cutoff_density_m3(frequency_ghz).map_err(js_err)?,
        x1_cutoff_density_m3: ecrh::x1_cutoff_density_m3(frequency_ghz, b_t).map_err(js_err)?,
    })
}

/// Whole-beamline ECRH accessibility report.
///
/// `sM`/`bT`/`neM3` are the caller beamline (positions in metres, fields in
/// tesla, densities in m⁻³); `frequencyGhz`/`harmonic` select the wave and
/// `teKev` the optional electron temperature. Returns
/// `{resonantFieldT, relativisticFieldT, o1CutoffDensityM3, resonanceM,
/// o1CutoffM, x1CutoffM}` with located positions in metres.
#[wasm_bindgen(js_name = ecrhAccessibility)]
pub fn ecrh_accessibility(
    s_m: Vec<f64>,
    b_t: Vec<f64>,
    ne_m3: Vec<f64>,
    frequency_ghz: f64,
    harmonic: f64,
    te_kev: Option<f64>,
) -> Result<JsValue, JsValue> {
    use nucleide_plasma_source::ecrh;
    let n = check_harmonic(harmonic)?;
    let out = ecrh::accessibility(&s_m, &b_t, &ne_m3, frequency_ghz, n, te_kev).map_err(js_err)?;
    #[derive(Serialize)]
    struct AccessibilityJson {
        #[serde(rename = "resonantFieldT")]
        resonant_field_t: f64,
        #[serde(rename = "relativisticFieldT")]
        relativistic_field_t: Option<f64>,
        #[serde(rename = "o1CutoffDensityM3")]
        o1_cutoff_density_m3: f64,
        #[serde(rename = "resonanceM")]
        resonance_m: Vec<f64>,
        #[serde(rename = "o1CutoffM")]
        o1_cutoff_m: Vec<f64>,
        #[serde(rename = "x1CutoffM")]
        x1_cutoff_m: Vec<f64>,
    }
    to_js(&AccessibilityJson {
        resonant_field_t: out.resonant_field_t,
        relativistic_field_t: out.relativistic_field_t,
        o1_cutoff_density_m3: out.o1_cutoff_density_m3,
        resonance_m: out.resonance_m,
        o1_cutoff_m: out.o1_cutoff_m,
        x1_cutoff_m: out.x1_cutoff_m,
    })
}

#[derive(Deserialize)]
struct LatticePointJson {
    #[serde(rename = "positionCm")]
    position_cm: [f64; 3],
    rate: f64,
    #[serde(rename = "tiKev")]
    ti_kev: f64,
}

/// Sample `n` seeded particles from an arbitrary-3D birth-rate lattice.
///
/// `spec` takes `points` (one `{positionCm: [x, y, z], rate, tiKev}` row per
/// node, lengths in cm, temperatures in keV), `reaction` (`"dt"`/`"dd"`),
/// optional `fieldPeriods`/`baseAngle` symmetry, `n` (capped at
/// [`MAX_FUSION_SAMPLES`]) and an integer `seed`. Returns `{kind:
/// "lattice", totalStrength, count, particles}` with one `{positionCm,
/// direction, energyMeV, weight}` row per particle.
#[wasm_bindgen(js_name = sampleLatticeSource)]
pub fn sample_lattice_source(spec: JsValue) -> Result<JsValue, JsValue> {
    use nucleide_plasma_source::lattice::{LatticePoint, LatticeSourceConfig, LatticeSymmetry};
    #[derive(Deserialize)]
    struct LatticeSpecJson {
        points: Vec<LatticePointJson>,
        reaction: String,
        #[serde(rename = "fieldPeriods")]
        field_periods: Option<f64>,
        #[serde(rename = "baseAngle")]
        base_angle: Option<f64>,
    }
    let parsed: LatticeSpecJson = serde_wasm_bindgen::from_value(spec.clone()).map_err(js_err)?;
    if parsed.points.is_empty() {
        return Err(js_err("lattice needs at least one point"));
    }
    let points: Vec<LatticePoint> = parsed
        .points
        .iter()
        .map(|p| LatticePoint {
            position_cm: p.position_cm,
            rate: p.rate,
            ion_temperature_kev: p.ti_kev,
        })
        .collect();
    let symmetry = match parsed.field_periods {
        Some(n) => {
            if !n.is_finite() || n.fract() != 0.0 || n < 1.0 {
                return Err(js_err(format!(
                    "fieldPeriods must be an integer >= 1 (got {n})"
                )));
            }
            Some(LatticeSymmetry::new(n as u32, parsed.base_angle.unwrap_or(0.0)).map_err(js_err)?)
        }
        None => None,
    };
    let config = LatticeSourceConfig {
        points,
        reaction: parse_fusion_reaction(&parsed.reaction)?,
        fuel_mixture: None,
        species_temperatures: None,
        symmetry,
        weight: 1.0,
    };
    let total_strength = config.total_strength().map_err(js_err)?;
    let (n, seed) = sample_args(&spec)?;
    let mut sampler =
        nucleide_plasma_source::lattice::LatticeSampler::new(config, seed).map_err(js_err)?;
    let particles = sampler.sample_n(n);
    #[derive(Serialize)]
    struct LatticeSampleResult {
        kind: String,
        #[serde(rename = "totalStrength")]
        total_strength: f64,
        count: usize,
        particles: Vec<FusionParticleJson>,
    }
    to_js(&LatticeSampleResult {
        kind: "lattice".to_string(),
        total_strength,
        count: particles.len(),
        particles: particle_json(&particles),
    })
}

/// TBR / blanket bookkeeping scalars (B1–B4 plus burn accounting).
///
/// `tritonsBred`/`sourceNeutrons` are the caller tallies, `portFractions`
/// the per-port coverage losses in `[0, 1)`, `fusionPowerMw` the fusion
/// power in MW, and `requiredTbr` the optional requirement threshold.
/// Returns `{rawTbr, effectiveTbr, margin, meets, burnGPerDay,
/// surplusGPerDay}` (`meets` is `null` without a threshold).
#[wasm_bindgen(js_name = tbrScalars)]
pub fn tbr_scalars(
    tritons_bred: f64,
    source_neutrons: f64,
    port_fractions: Vec<f64>,
    fusion_power_mw: f64,
    required_tbr: Option<f64>,
) -> Result<JsValue, JsValue> {
    use nucleide_blanket as blanket;
    let raw = blanket::tbr_from_tallies(tritons_bred, source_neutrons).map_err(js_err)?;
    let effective = blanket::apply_port_penalty(raw, &port_fractions).map_err(js_err)?;
    let meets = required_tbr
        .map(|r| blanket::meets_requirement(effective, r))
        .transpose()
        .map_err(js_err)?;
    #[derive(Serialize)]
    struct TbrScalarsJson {
        #[serde(rename = "rawTbr")]
        raw_tbr: f64,
        #[serde(rename = "effectiveTbr")]
        effective_tbr: f64,
        margin: f64,
        meets: Option<bool>,
        #[serde(rename = "burnGPerDay")]
        burn_g_per_day: f64,
        #[serde(rename = "surplusGPerDay")]
        surplus_g_per_day: f64,
    }
    to_js(&TbrScalarsJson {
        raw_tbr: raw,
        effective_tbr: effective,
        margin: blanket::breeding_margin(effective).map_err(js_err)?,
        meets,
        burn_g_per_day: blanket::tritium_burn_g_per_day(fusion_power_mw).map_err(js_err)?,
        surplus_g_per_day: blanket::net_surplus_g_per_day(effective, fusion_power_mw)
            .map_err(js_err)?,
    })
}

/// Coil fast-flux / fluence over caller groups.
///
/// `flux` holds the per-group fluxes, `bounds` the `G+1` MeV group edges,
/// `thresholdMev` the fast threshold, and `seconds` the hold time. Returns
/// `{fastFlux, fastFluence}`.
#[wasm_bindgen(js_name = coilFastFlux)]
pub fn coil_fast_flux(
    flux: Vec<f64>,
    bounds: Vec<f64>,
    threshold_mev: f64,
    seconds: f64,
) -> Result<JsValue, JsValue> {
    use nucleide_damage::coil;
    #[derive(Serialize)]
    struct CoilFluxJson {
        #[serde(rename = "fastFlux")]
        fast_flux: f64,
        #[serde(rename = "fastFluence")]
        fast_fluence: f64,
    }
    to_js(&CoilFluxJson {
        fast_flux: coil::fast_flux(&flux, &bounds, threshold_mev).map_err(js_err)?,
        fast_fluence: coil::fast_fluence(&flux, &bounds, threshold_mev, seconds).map_err(js_err)?,
    })
}

/// Weakest-link coil lifetime from caller limit/rate tables.
///
/// `limits`/`rates` share units per channel; zero-rate channels never age.
/// Returns `{seconds, limiting}` (`limiting` is `null` when nothing ages).
#[wasm_bindgen(js_name = coilLifetime)]
pub fn coil_lifetime(limits: Vec<f64>, rates: Vec<f64>) -> Result<JsValue, JsValue> {
    use nucleide_damage::coil;
    let out = coil::coil_lifetime(&limits, &rates).map_err(js_err)?;
    #[derive(Serialize)]
    struct CoilLifetimeJson {
        seconds: f64,
        limiting: Option<usize>,
    }
    to_js(&CoilLifetimeJson {
        seconds: out.seconds,
        limiting: out.limiting,
    })
}

#[derive(Deserialize)]
struct SubletActivityEntryJson {
    nuclide: String,
    #[serde(rename = "activityBq")]
    activity_bq: f64,
    irt: u8,
    #[serde(rename = "alphaFrac")]
    alpha_frac: Option<f64>,
}

/// Sublet S1 total activity with the α/β/γ split (Bq).
///
/// `entries` holds one `{nuclide, activityBq, irt, alphaFrac}` row per
/// inventory nuclide (`alphaFrac` only for IRT 12, 13, 15). Returns
/// `{totalBq, alphaBq, betaBq, gammaBq, exTritiumBq}`.
#[wasm_bindgen(js_name = subletActivity)]
pub fn sublet_activity(entries: JsValue) -> Result<JsValue, JsValue> {
    use nucleide_alara_io::sublet::{total_activity, ActivityEntry};
    let rows: Vec<SubletActivityEntryJson> =
        serde_wasm_bindgen::from_value(entries).map_err(js_err)?;
    let parsed: Vec<ActivityEntry> = rows
        .iter()
        .map(|r| {
            Ok::<_, JsValue>(ActivityEntry {
                nuclide: r
                    .nuclide
                    .parse::<NuclideId>()
                    .map_err(|e| js_err(format!("`{}`: {e}", r.nuclide)))?,
                activity_bq: r.activity_bq,
                irt: r.irt,
                alpha_frac: r.alpha_frac,
            })
        })
        .collect::<Result<_, _>>()?;
    let out = total_activity(&parsed).map_err(js_err)?;
    #[derive(Serialize)]
    struct SubletActivityJson {
        #[serde(rename = "totalBq")]
        total_bq: f64,
        #[serde(rename = "alphaBq")]
        alpha_bq: f64,
        #[serde(rename = "betaBq")]
        beta_bq: f64,
        #[serde(rename = "gammaBq")]
        gamma_bq: f64,
        #[serde(rename = "exTritiumBq")]
        ex_tritium_bq: f64,
    }
    to_js(&SubletActivityJson {
        total_bq: out.total_bq,
        alpha_bq: out.alpha_bq,
        beta_bq: out.beta_bq,
        gamma_bq: out.gamma_bq,
        ex_tritium_bq: out.ex_tritium_bq,
    })
}

#[derive(Deserialize)]
struct SubletHeatEntryJson {
    nuclide: String,
    #[serde(rename = "activityBq")]
    activity_bq: f64,
    #[serde(rename = "eAlphaEv")]
    e_alpha_ev: f64,
    #[serde(rename = "eBetaEv")]
    e_beta_ev: f64,
    #[serde(rename = "eGammaEv")]
    e_gamma_ev: f64,
}

/// Sublet S2 decay heat with the α/β/γ split (kW).
///
/// `entries` holds one `{nuclide, activityBq, eAlphaEv, eBetaEv, eGammaEv}`
/// row per inventory nuclide (caller decay energies in eV). Returns
/// `{alphaKw, betaKw, gammaKw, totalKw, exTritiumKw}`.
#[wasm_bindgen(js_name = subletDecayHeat)]
pub fn sublet_decay_heat(entries: JsValue) -> Result<JsValue, JsValue> {
    use nucleide_alara_io::sublet::{decay_heat, DecayHeatEntry};
    let rows: Vec<SubletHeatEntryJson> = serde_wasm_bindgen::from_value(entries).map_err(js_err)?;
    let parsed: Vec<DecayHeatEntry> = rows
        .iter()
        .map(|r| {
            Ok::<_, JsValue>(DecayHeatEntry {
                nuclide: r
                    .nuclide
                    .parse::<NuclideId>()
                    .map_err(|e| js_err(format!("`{}`: {e}", r.nuclide)))?,
                activity_bq: r.activity_bq,
                e_alpha_ev: r.e_alpha_ev,
                e_beta_ev: r.e_beta_ev,
                e_gamma_ev: r.e_gamma_ev,
            })
        })
        .collect::<Result<_, _>>()?;
    let out = decay_heat(&parsed).map_err(js_err)?;
    #[derive(Serialize)]
    struct SubletHeatJson {
        #[serde(rename = "alphaKw")]
        alpha_kw: f64,
        #[serde(rename = "betaKw")]
        beta_kw: f64,
        #[serde(rename = "gammaKw")]
        gamma_kw: f64,
        #[serde(rename = "totalKw")]
        total_kw: f64,
        #[serde(rename = "exTritiumKw")]
        ex_tritium_kw: f64,
    }
    to_js(&SubletHeatJson {
        alpha_kw: out.alpha_kw,
        beta_kw: out.beta_kw,
        gamma_kw: out.gamma_kw,
        total_kw: out.total_kw,
        ex_tritium_kw: out.ex_tritium_kw,
    })
}

#[derive(Deserialize)]
struct SubletHazardEntryJson {
    nuclide: String,
    #[serde(rename = "activityBq")]
    activity_bq: f64,
    #[serde(rename = "coeffSvPerBq")]
    coeff_sv_per_bq: f64,
}

fn sublet_hazard_entries(entries: JsValue) -> Result<Vec<nucleide_alara_io::HazardEntry>, JsValue> {
    let rows: Vec<SubletHazardEntryJson> =
        serde_wasm_bindgen::from_value(entries).map_err(js_err)?;
    rows.iter()
        .map(|r| {
            Ok::<_, JsValue>(nucleide_alara_io::HazardEntry {
                nuclide: r
                    .nuclide
                    .parse::<NuclideId>()
                    .map_err(|e| js_err(format!("`{}`: {e}", r.nuclide)))?,
                activity_bq: r.activity_bq,
                coeff_sv_per_bq: r.coeff_sv_per_bq,
            })
        })
        .collect::<Result<_, _>>()
}

/// Sublet S4 ingestion hazard `ΣAi·e^ing_i` (Sv, 50-year committed).
///
/// `entries` holds one `{nuclide, activityBq, coeffSvPerBq}` row per
/// inventory nuclide (caller ingestion coefficients in Sv/Bq, never
/// vendored). Returns `{totalSv, exTritiumSv}`.
#[wasm_bindgen(js_name = subletIngestionHazard)]
pub fn sublet_ingestion_hazard(entries: JsValue) -> Result<JsValue, JsValue> {
    let parsed = sublet_hazard_entries(entries)?;
    let out = nucleide_alara_io::sublet::ingestion_hazard(&parsed).map_err(js_err)?;
    #[derive(Serialize)]
    struct SubletHazardJson {
        #[serde(rename = "totalSv")]
        total_sv: f64,
        #[serde(rename = "exTritiumSv")]
        ex_tritium_sv: f64,
    }
    to_js(&SubletHazardJson {
        total_sv: out.total_sv,
        ex_tritium_sv: out.ex_tritium_sv,
    })
}

/// Sublet S5 inhalation hazard `ΣAi·e^inh_i` (Sv, 50-year committed).
///
/// Same contract as [`sublet_ingestion_hazard`] with the caller inhalation
/// coefficients. Returns `{totalSv, exTritiumSv}`.
#[wasm_bindgen(js_name = subletInhalationHazard)]
pub fn sublet_inhalation_hazard(entries: JsValue) -> Result<JsValue, JsValue> {
    let parsed = sublet_hazard_entries(entries)?;
    let out = nucleide_alara_io::sublet::inhalation_hazard(&parsed).map_err(js_err)?;
    #[derive(Serialize)]
    struct SubletHazardJson {
        #[serde(rename = "totalSv")]
        total_sv: f64,
        #[serde(rename = "exTritiumSv")]
        ex_tritium_sv: f64,
    }
    to_js(&SubletHazardJson {
        total_sv: out.total_sv,
        ex_tritium_sv: out.ex_tritium_sv,
    })
}

#[derive(Deserialize)]
struct SubletTransportEntryJson {
    nuclide: String,
    #[serde(rename = "activityBq")]
    activity_bq: f64,
    #[serde(rename = "a2Tbq")]
    a2_tbq: f64,
}

/// Sublet S6 transport ratio `ΣAi/(A2,i·C2)` with the effective A2 (TBq).
///
/// `entries` holds one `{nuclide, activityBq, a2Tbq}` row per inventory
/// nuclide (caller transport limits in TBq, never vendored). Returns
/// `{ratio, totalBq, effectiveA2Tbq}` (the effective A2 is defined by the
/// ratio, exactly).
#[wasm_bindgen(js_name = subletTransportRatio)]
pub fn sublet_transport_ratio(entries: JsValue) -> Result<JsValue, JsValue> {
    let rows: Vec<SubletTransportEntryJson> =
        serde_wasm_bindgen::from_value(entries).map_err(js_err)?;
    let parsed: Vec<nucleide_alara_io::TransportEntry> = rows
        .iter()
        .map(|r| {
            Ok::<_, JsValue>(nucleide_alara_io::TransportEntry {
                nuclide: r
                    .nuclide
                    .parse::<NuclideId>()
                    .map_err(|e| js_err(format!("`{}`: {e}", r.nuclide)))?,
                activity_bq: r.activity_bq,
                a2_tbq: r.a2_tbq,
            })
        })
        .collect::<Result<_, _>>()?;
    let out = nucleide_alara_io::sublet::transport_ratio(&parsed).map_err(js_err)?;
    #[derive(Serialize)]
    struct SubletTransportJson {
        ratio: f64,
        #[serde(rename = "totalBq")]
        total_bq: f64,
        #[serde(rename = "effectiveA2Tbq")]
        effective_a2_tbq: f64,
    }
    to_js(&SubletTransportJson {
        ratio: out.ratio,
        total_bq: out.total_bq,
        effective_a2_tbq: out.effective_a2_tbq,
    })
}

#[derive(Deserialize)]
struct SubletIaeaEntryJson {
    nuclide: String,
    #[serde(rename = "activityBq")]
    activity_bq: f64,
    #[serde(rename = "limitBqPerKg")]
    limit_bq_per_kg: f64,
}

/// Sublet S7 IAEA clearance index `ΣAi/(Mtot·Li)` (dimensionless).
///
/// `totalMassKg` is the total mass `Mtot` (finite, `> 0`); `entries` holds
/// one `{nuclide, activityBq, limitBqPerKg}` row per inventory nuclide
/// (caller IAEA levels in Bq/kg, never vendored). Returns `{index,
/// clearanceClass ("satisfied" at `index <= 1`, boundary included, else
/// "exceeded"), maxFraction, maxNuclide}` (`maxNuclide` is `null` for an
/// empty inventory). Screening arithmetic only, never a compliance decision.
#[wasm_bindgen(js_name = subletIaeaClearance)]
pub fn sublet_iaea_clearance(total_mass_kg: f64, entries: JsValue) -> Result<JsValue, JsValue> {
    let rows: Vec<SubletIaeaEntryJson> = serde_wasm_bindgen::from_value(entries).map_err(js_err)?;
    let parsed: Vec<nucleide_alara_io::IaeaEntry> = rows
        .iter()
        .map(|r| {
            Ok::<_, JsValue>(nucleide_alara_io::IaeaEntry {
                nuclide: r
                    .nuclide
                    .parse::<NuclideId>()
                    .map_err(|e| js_err(format!("`{}`: {e}", r.nuclide)))?,
                activity_bq: r.activity_bq,
                limit_bq_per_kg: r.limit_bq_per_kg,
            })
        })
        .collect::<Result<_, _>>()?;
    let out =
        nucleide_alara_io::sublet::iaea_clearance_index(total_mass_kg, &parsed).map_err(js_err)?;
    #[derive(Serialize)]
    struct SubletIaeaJson {
        index: f64,
        #[serde(rename = "clearanceClass")]
        clearance_class: String,
        #[serde(rename = "maxFraction")]
        max_fraction: f64,
        #[serde(rename = "maxNuclide")]
        max_nuclide: Option<String>,
    }
    to_js(&SubletIaeaJson {
        index: out.index,
        clearance_class: out.class.to_string(),
        max_fraction: out.max_fraction,
        max_nuclide: out.max_nuclide.map(|id| id.to_name()),
    })
}

#[derive(Deserialize)]
struct SubletDoseGroupJson {
    intensity: f64,
    #[serde(rename = "muAir")]
    mu_air: f64,
    mu: f64,
}

fn sublet_dose_groups(groups: JsValue) -> Result<Vec<nucleide_alara_io::DoseGroup>, JsValue> {
    let rows: Vec<SubletDoseGroupJson> = serde_wasm_bindgen::from_value(groups).map_err(js_err)?;
    Ok(rows
        .iter()
        .map(|r| nucleide_alara_io::DoseGroup {
            intensity: r.intensity,
            mu_air: r.mu_air,
            mu: r.mu,
        })
        .collect())
}

/// Sublet S3 slab dose `C·B/2·Σ μa/μm·Sγ` (Sv/h; `B = 2`, `C = 3.6e9·|e|`).
///
/// `activityBqPerKg` is the caller specific activity; `groups` holds one
/// `{intensity, muAir, mu}` row per gamma group (caller attenuation,
/// never vendored; the slab ratio divides by `mu`). Returns `{doseSvPerH}`.
#[wasm_bindgen(js_name = subletDoseSlab)]
pub fn sublet_dose_slab(activity_bq_per_kg: f64, groups: JsValue) -> Result<JsValue, JsValue> {
    let parsed = sublet_dose_groups(groups)?;
    let out = nucleide_alara_io::dose::dose_slab(activity_bq_per_kg, &parsed).map_err(js_err)?;
    #[derive(Serialize)]
    struct SubletDoseJson {
        #[serde(rename = "doseSvPerH")]
        dose_sv_per_h: f64,
    }
    to_js(&SubletDoseJson {
        dose_sv_per_h: out.dose_sv_per_h,
    })
}

/// Sublet S3 point dose `C·Σ μa/(4πr²)·e^(−μr)·mₛ·Sγ` (Sv/h).
///
/// Same group contract as [`sublet_dose_slab`], plus `sourceMassKg` (`mₛ`)
/// and `distanceM` (`r`). A finite `r` below 0.3 m clamps to 0.3 m and the
/// returned dict reports it loudly via `clamped: true` with the
/// `distanceUsedM` actually used (never silent). Returns `{doseSvPerH,
/// distanceUsedM, clamped}`.
#[wasm_bindgen(js_name = subletDosePoint)]
pub fn sublet_dose_point(
    activity_bq_per_kg: f64,
    source_mass_kg: f64,
    distance_m: f64,
    groups: JsValue,
) -> Result<JsValue, JsValue> {
    let parsed = sublet_dose_groups(groups)?;
    let out = nucleide_alara_io::dose::dose_point(
        activity_bq_per_kg,
        source_mass_kg,
        distance_m,
        &parsed,
    )
    .map_err(js_err)?;
    #[derive(Serialize)]
    struct SubletPointDoseJson {
        #[serde(rename = "doseSvPerH")]
        dose_sv_per_h: f64,
        #[serde(rename = "distanceUsedM")]
        distance_used_m: f64,
        clamped: bool,
    }
    to_js(&SubletPointDoseJson {
        dose_sv_per_h: out.dose_sv_per_h,
        distance_used_m: out.distance_used_m,
        clamped: out.clamped,
    })
}

/// Sublet S3 mixture fold `μm(Eᵢ) = Σⱼ fⱼ·μmⱼ(Eᵢ)`.
///
/// `fractions` weights every group of `elementMus` (fractions finite,
/// `>= 0`, summing to 1). Returns one mixture coefficient per group.
#[wasm_bindgen(js_name = subletDoseMixtureMu)]
pub fn sublet_dose_mixture_mu(
    fractions: Vec<f64>,
    element_mus: JsValue,
) -> Result<JsValue, JsValue> {
    let rows: Vec<Vec<f64>> = serde_wasm_bindgen::from_value(element_mus).map_err(js_err)?;
    let out = nucleide_alara_io::dose::mixture_mu(&fractions, &rows).map_err(js_err)?;
    to_js(&out)
}

/// He/dpa ratio UQ: per-draw ratios over the seeded joint-MVN block.
///
/// `mean`/`cov` describe relative perturbations of the stacked `[flux,
/// heResponse, damageResponse]` vector (dimension `3G`). Each draw refolds
/// He and dpa and forms the ratio per draw; returns `{metric, nominal,
/// mean, std, expected, analyticStd, k, n, seed, passed}` with the
/// second-order bias-corrected expectation and the first-order
/// delta-propagated std. A draw at non-positive dpa fails loudly.
#[wasm_bindgen(js_name = damageHeDpaRatioUq)]
#[allow(clippy::too_many_arguments)]
pub fn damage_he_dpa_ratio_uq(
    flux: Vec<f64>,
    he_response: Vec<f64>,
    damage_response: Vec<f64>,
    bounds: Vec<f64>,
    seconds: f64,
    mean: Vec<f64>,
    cov: JsValue,
    n: usize,
    seed: f64,
    k: f64,
) -> Result<JsValue, JsValue> {
    let cov: Vec<Vec<f64>> = serde_wasm_bindgen::from_value(cov).map_err(js_err)?;
    let out = nucleide_damage::he_dpa_ratio_uq(
        &flux,
        &he_response,
        &damage_response,
        &bounds,
        seconds,
        &mean,
        &cov,
        n,
        check_seed(seed)?,
        k,
    )
    .map_err(js_err)?;
    #[derive(Serialize)]
    struct RatioUqJson {
        metric: String,
        nominal: f64,
        mean: f64,
        std: f64,
        expected: f64,
        #[serde(rename = "analyticStd")]
        analytic_std: f64,
        q16: f64,
        q50: f64,
        q84: f64,
        #[serde(rename = "expectedQ16")]
        expected_q16: f64,
        #[serde(rename = "expectedQ50")]
        expected_q50: f64,
        #[serde(rename = "expectedQ84")]
        expected_q84: f64,
        #[serde(rename = "quantilesPassed")]
        quantiles_passed: bool,
        k: f64,
        n: usize,
        seed: u64,
        passed: bool,
    }
    to_js(&RatioUqJson {
        metric: out.metric.name().to_string(),
        nominal: out.nominal,
        mean: out.mean,
        std: out.std,
        expected: out.expected,
        analytic_std: out.analytic_std,
        q16: out.q16,
        q50: out.q50,
        q84: out.q84,
        expected_q16: out.expected_q16,
        expected_q50: out.expected_q50,
        expected_q84: out.expected_q84,
        quantiles_passed: out.quantiles_passed,
        k: out.k,
        n: out.n,
        seed: out.seed,
        passed: out.passed,
    })
}

fn indata_value_json(value: &nucleide_equilib_io::IndataValue) -> serde_json::Value {
    match value {
        nucleide_equilib_io::IndataValue::Int(v) => serde_json::json!({"Int": v}),
        nucleide_equilib_io::IndataValue::Float(v) => serde_json::json!({"Float": v}),
        nucleide_equilib_io::IndataValue::Bool(v) => serde_json::json!({"Bool": v}),
        nucleide_equilib_io::IndataValue::Str(v) => serde_json::json!({"Str": v}),
    }
}

/// Parse a STELLOPT-style `&INDATA ... /` block.
///
/// Returns `{scalars, indexed}` where `scalars` maps upper-cased names to
/// `{Int|Float|Bool|Str: value}` and `indexed` maps names to
/// `[{index: [...], value}]` rows.
#[wasm_bindgen(js_name = parseIndata)]
pub fn parse_indata(text: &str) -> Result<JsValue, JsValue> {
    let parsed = nucleide_equilib_io::parse_indata(text).map_err(js_err)?;
    #[derive(Serialize)]
    struct IndexedRowJson {
        index: Vec<i64>,
        value: serde_json::Value,
    }
    #[derive(Serialize)]
    struct IndataJson {
        scalars: std::collections::BTreeMap<String, serde_json::Value>,
        indexed: std::collections::BTreeMap<String, Vec<IndexedRowJson>>,
    }
    to_js(&IndataJson {
        scalars: parsed
            .scalars
            .iter()
            .map(|(k, v)| (k.clone(), indata_value_json(v)))
            .collect(),
        indexed: parsed
            .indexed
            .iter()
            .map(|(k, rows)| {
                (
                    k.clone(),
                    rows.iter()
                        .map(|(idx, v)| IndexedRowJson {
                            index: (*idx).clone(),
                            value: indata_value_json(v),
                        })
                        .collect(),
                )
            })
            .collect(),
    })
}

/// Map a birth-rate density field onto the wall (W1–W2).
///
/// `sEdges` holds the `nr + 1` radial edges, `ntheta`/`nzeta` the angular
/// grid counts, `nfp` the field-period count, and `birth`/`jacobian` the
/// per-voxel densities and `√g` values in radial-major order. Returns
/// `{ntheta, nzeta, loads, total}` with nested per-node loads and the
/// one-field-period total.
#[wasm_bindgen(js_name = wallLoad)]
pub fn wall_load(
    s_edges: Vec<f64>,
    ntheta: f64,
    nzeta: f64,
    nfp: f64,
    birth: Vec<f64>,
    jacobian: Vec<f64>,
) -> Result<JsValue, JsValue> {
    for (label, v) in [("ntheta", ntheta), ("nzeta", nzeta), ("nfp", nfp)] {
        if !v.is_finite() || v.fract() != 0.0 || v < 1.0 {
            return Err(js_err(format!("{label} must be an integer >= 1 (got {v})")));
        }
    }
    let out = nucleide_equilib_io::wall_load(
        &s_edges,
        ntheta as usize,
        nzeta as usize,
        nfp as u32,
        &birth,
        &jacobian,
    )
    .map_err(js_err)?;
    #[derive(Serialize)]
    struct WallLoadJson {
        ntheta: usize,
        nzeta: usize,
        loads: Vec<Vec<f64>>,
        total: f64,
    }
    to_js(&WallLoadJson {
        ntheta: out.ntheta,
        nzeta: out.nzeta,
        loads: out.to_nested(),
        total: out.total,
    })
}
