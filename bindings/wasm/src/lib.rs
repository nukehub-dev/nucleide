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
    #[wasm_bindgen(constructor)]
    pub fn new(name: &str) -> Result<WasmNuclide, JsValue> {
        NuclideId::from_name(name)
            .map(|inner| WasmNuclide { inner })
            .map_err(js_err)
    }

    #[wasm_bindgen(js_name = fromZzaaam)]
    pub fn from_zzaaam(v: u32) -> Result<WasmNuclide, JsValue> {
        NuclideId::from_zzaaam(v)
            .map(|inner| WasmNuclide { inner })
            .map_err(js_err)
    }

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

    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.inner.to_name()
    }

    #[wasm_bindgen(getter)]
    pub fn nucid(&self) -> u32 {
        self.inner.nucid()
    }

    #[wasm_bindgen(getter)]
    pub fn z(&self) -> u32 {
        self.inner.z()
    }

    #[wasm_bindgen(getter)]
    pub fn a(&self) -> u32 {
        self.inner.a()
    }

    #[wasm_bindgen(getter)]
    pub fn state(&self) -> u32 {
        self.inner.state()
    }

    #[wasm_bindgen(getter)]
    pub fn zzaaam(&self) -> u32 {
        self.inner.zzaaam()
    }

    #[wasm_bindgen(getter)]
    pub fn zaid(&self) -> u32 {
        nucleide_nuclei::dialects::to_zaid(self.inner)
    }

    #[wasm_bindgen(getter)]
    pub fn zzllaaam(&self) -> String {
        nucleide_nuclei::dialects::zzllaaam(self.inner)
    }

    #[wasm_bindgen(getter)]
    pub fn serpent(&self) -> String {
        nucleide_nuclei::dialects::serpent(self.inner)
    }

    #[wasm_bindgen(getter)]
    pub fn nist(&self) -> String {
        nucleide_nuclei::dialects::nist(self.inner)
    }

    #[wasm_bindgen(getter)]
    pub fn cinder(&self) -> u32 {
        nucleide_nuclei::dialects::to_cinder(self.inner)
    }

    #[wasm_bindgen(getter)]
    pub fn alara(&self) -> String {
        nucleide_nuclei::dialects::alara(self.inner)
    }

    #[wasm_bindgen(getter)]
    pub fn sza(&self) -> u32 {
        nucleide_nuclei::dialects::to_sza(self.inner)
    }

    #[wasm_bindgen(getter)]
    pub fn mass(&self) -> Option<f64> {
        nucleide_nuclei::data::atomic_mass(self.inner.nucid())
    }

    #[wasm_bindgen(getter)]
    pub fn abundance(&self) -> Option<f64> {
        nucleide_nuclei::data::natural_abundance(self.inner.nucid())
    }

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

#[wasm_bindgen]
pub fn atomic_mass(key: &str) -> Result<Option<f64>, JsValue> {
    let id = resolve_nucid(key)?;
    Ok(nucleide_nuclei::data::atomic_mass(id.nucid()))
}

#[wasm_bindgen]
pub fn natural_abundance(key: &str) -> Result<Option<f64>, JsValue> {
    let id = resolve_nucid(key)?;
    Ok(nucleide_nuclei::data::natural_abundance(id.nucid()))
}

#[wasm_bindgen]
pub fn half_life(key: &str) -> Result<Option<f64>, JsValue> {
    let id = resolve_nucid(key)?;
    Ok(nucleide_nuclei::data::half_life(id.nucid()))
}

#[wasm_bindgen]
pub fn decay_constant(key: &str) -> Result<Option<f64>, JsValue> {
    let id = resolve_nucid(key)?;
    Ok(nucleide_nuclei::data::half_life(id.nucid()).map(|t| std::f64::consts::LN_2 / t))
}

#[wasm_bindgen]
pub fn q_value_capture(key: &str) -> Result<Option<f64>, JsValue> {
    let id = resolve_nucid(key)?;
    Ok(nucleide_nuclei::data::q_value_neutron_capture(id.nucid()))
}

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

    #[wasm_bindgen(getter)]
    pub fn mass(&self) -> f64 {
        self.inner.mass()
    }

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

    #[wasm_bindgen(getter)]
    pub fn len(&self) -> usize {
        self.inner.len()
    }

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

    #[wasm_bindgen(getter)]
    pub fn alpha(&self) -> f64 {
        self.inner.alpha
    }

    #[wasm_bindgen(getter, js_name = feedAssay)]
    pub fn feed_assay(&self) -> f64 {
        self.inner.x_feed_j
    }

    #[wasm_bindgen(getter, js_name = productAssay)]
    pub fn product_assay(&self) -> f64 {
        self.inner.x_prod_j
    }

    #[wasm_bindgen(getter, js_name = tailsAssay)]
    pub fn tails_assay(&self) -> f64 {
        self.inner.x_tail_j
    }

    #[wasm_bindgen(getter, js_name = stagesEnriching)]
    pub fn stages_enriching(&self) -> f64 {
        self.inner.N
    }

    #[wasm_bindgen(getter, js_name = stagesStripping)]
    pub fn stages_stripping(&self) -> f64 {
        self.inner.M
    }

    #[wasm_bindgen(getter, js_name = swuPerFeed)]
    pub fn swu_per_feed(&self) -> f64 {
        self.inner.swu_per_feed
    }

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

    #[wasm_bindgen(getter)]
    pub fn nuclides(&self) -> Result<JsValue, JsValue> {
        let names: Vec<String> = self.inner.nuclides.iter().map(|n| n.name.clone()).collect();
        to_js(&names)
    }
}

/// Run one CRAM depletion step.
///
/// `n0` is a JS object mapping nuclide names to atom counts. `rates` is a JS
/// object mapping `"Name:reaction"` strings to one-group rates [1/s]. `order`
/// is 16 or 48.
#[wasm_bindgen]
pub fn deplete(
    chain: &WasmChain,
    n0: JsValue,
    dt: f64,
    rates: JsValue,
    order: u8,
) -> Result<JsValue, JsValue> {
    let n0: BTreeMap<String, f64> = serde_wasm_bindgen::from_value(n0).map_err(js_err)?;
    let rates: BTreeMap<String, f64> = serde_wasm_bindgen::from_value(rates).map_err(js_err)?;

    let order = match order {
        16 => nucleide_depletion::Order::Order16,
        48 => nucleide_depletion::Order::Order48,
        other => return Err(js_err(format!("unsupported CRAM order {other}"))),
    };

    let mut reaction_rates = nucleide_depletion::ReactionRates::new();
    for (key, v) in rates {
        let (nuc, rx) = key
            .split_once(':')
            .ok_or_else(|| js_err(format!("rate key `{key}` must be `Name:reaction`")))?;
        let idx = chain
            .inner
            .index_of(nuc)
            .ok_or_else(|| js_err(format!("rate for unknown nuclide `{nuc}`")))?;
        reaction_rates
            .entry(idx)
            .or_default()
            .insert(rx.to_string(), v);
    }

    let sys = nucleide_depletion::DepletionSystem::build((*chain.inner).clone(), &reaction_rates)
        .map_err(js_err)?;
    let result = nucleide_depletion::deplete(&sys, order, &n0, dt).map_err(js_err)?;
    to_js(&result.atoms)
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
