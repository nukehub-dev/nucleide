//! DOE/PNNL Materials Compendium ingestion (Revision 2, 411 materials).
//!
//! Parses the official `MaterialsCompendium.json` of the DOE/PNNL
//! Materials Compendium, Revision 2 (companion license vendored under
//! `fixtures/data/`) into typed entries convertible to [`Material`]s with
//! full isotopic weight fractions.
//!
//! Schema notes (verified against the dataset):
//! - top level `{siteVersion: String, data: [entry; 411]}`;
//! - every entry carries `Elements[].Isotopes[]` — no degenerate cases;
//! - `Isotopes[].WeightFraction` is the material-level mass fraction
//!   (element-level fractions live in `IsotopicWeightFraction`);
//! - `ZAID` is numeric (`1001`, `95242`-style metastables absent here);
//! - names and MatNum values are unique across all 411 entries.

use std::collections::BTreeMap;
use std::path::Path;

use serde::Deserialize;

use crate::material::Material;
use nuclei::NuclideId;

/// Errors from compendium loading.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    Io(String),
    Json(String),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::Io(m) => write!(f, "io error: {m}"),
            Error::Json(m) => write!(f, "compendium JSON error: {m}"),
        }
    }
}
impl std::error::Error for Error {}

/// Accept either a JSON string or an array of strings (the upstream dataset
/// is inconsistent: most entries use arrays, some a single bare string).
fn string_or_vec<'de, D>(de: D) -> Result<Vec<String>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    struct V;
    impl<'de2> serde::de::Visitor<'de2> for V {
        type Value = Vec<String>;
        fn expecting(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            f.write_str("string or list of strings")
        }
        fn visit_str<E: serde::de::Error>(self, v: &str) -> Result<Self::Value, E> {
            Ok(vec![v.to_string()])
        }
        fn visit_seq<S: serde::de::SeqAccess<'de2>>(
            self,
            mut seq: S,
        ) -> Result<Self::Value, S::Error> {
            let mut out = Vec::new();
            while let Some(s) = seq.next_element::<String>()? {
                out.push(s);
            }
            Ok(out)
        }
    }
    de.deserialize_any(V)
}

/// One isotope row of an element inside a compendium material.
#[derive(Debug, Clone, Deserialize)]
pub struct CompendiumIsotope {
    #[serde(rename = "Isotope")]
    pub isotope: String,
    /// MCNP-style ZAID; serialized as a JSON string upstream ("1001").
    #[serde(rename = "ZAID")]
    pub zaid: String,
    #[serde(rename = "Abundance")]
    pub abundance: f64,
    /// Mass fraction within the parent element.
    #[serde(rename = "IsotopicWeightFraction")]
    pub isotopic_weight_fraction: f64,
    /// Mass fraction of the whole material.
    #[serde(rename = "WeightFraction")]
    pub weight_fraction: f64,
}

/// One elemental constituent.
#[derive(Debug, Clone, Deserialize)]
pub struct CompendiumElement {
    #[serde(rename = "Element")]
    pub element: String,
    #[serde(rename = "AtomFraction")]
    pub atom_fraction: f64,
    #[serde(rename = "Isotopes")]
    pub isotopes: Vec<CompendiumIsotope>,
}

/// One compendium material entry.
#[derive(Debug, Clone, Deserialize)]
pub struct CompendiumEntry {
    #[serde(rename = "Name")]
    pub name: String,
    #[serde(rename = "Acronym", default, deserialize_with = "string_or_vec")]
    pub acronym: Vec<String>,
    #[serde(rename = "MatNum")]
    pub mat_num: u32,
    /// Nominal density [g/cm³].
    #[serde(rename = "Density")]
    pub density: f64,
    #[serde(rename = "MaterialAtomDensity", default)]
    pub atom_density: f64,
    #[serde(rename = "Source", default)]
    pub source: String,
    #[serde(rename = "Comment", default, deserialize_with = "string_or_vec")]
    pub comment: Vec<String>,
    #[serde(rename = "Elements")]
    pub elements: Vec<CompendiumElement>,
}

impl CompendiumEntry {
    /// All (isotope ZAID → material-level weight fraction) pairs.
    pub fn weight_fractions(&self) -> BTreeMap<u32, f64> {
        let mut out = BTreeMap::new();
        for el in &self.elements {
            for iso in &el.isotopes {
                let zaid: u32 = match iso.zaid.parse() {
                    Ok(v) => v,
                    Err(_) => continue,
                };
                *out.entry(zaid).or_insert(0.0) += iso.weight_fraction;
            }
        }
        out
    }

    /// Convert to a [`Material`] whose composition holds relative masses
    /// equal to the isotopic weight fractions (normalized to 1 g total).
    ///
    /// Density and provenance are attached as metadata; set the real mass or
    /// density separately when building transport inputs.
    pub fn to_material(&self) -> Result<Material, Error> {
        use nuclei::dialects;
        let mut mat = Material::new();
        for (zaid, wf) in self.weight_fractions() {
            if wf <= 0.0 {
                continue;
            }
            // Natural-element zaids (AAA == 0) cannot map to a single
            // nuclide; keep them as placeholder ground-state ids like the
            // mcnp-io inp convention.
            let id = dialects::from_zaid(zaid)
                .unwrap_or_else(|_| NuclideId::from_nucid((zaid / 1000) * 10_000_000));
            mat.add_nuclide(id, wf);
        }
        mat.set_metadata(Some(serde_json::json!({
            "source": "DOE/PNNL Materials Compendium Rev.2",
            "acronym": self.acronym,
            "mat_num": self.mat_num,
            "density_g_cm3": self.density,
            "atom_density": self.atom_density,
            "reference": self.source,
        })));
        Ok(mat)
    }
}

/// The full parsed compendium with fast lookups.
#[derive(Debug, Clone)]
pub struct MaterialsLibrary {
    pub site_version: String,
    pub entries: Vec<CompendiumEntry>,
    by_name: BTreeMap<String, usize>,
    by_matnum: BTreeMap<u32, usize>,
}

impl MaterialsLibrary {
    /// Parse compendium JSON text.
    pub fn from_json(text: &str) -> Result<Self, Error> {
        #[derive(Deserialize)]
        struct Top {
            #[serde(rename = "siteVersion")]
            site_version: String,
            data: Vec<CompendiumEntry>,
        }
        let top: Top = serde_json::from_str(text).map_err(|e| Error::Json(e.to_string()))?;
        let mut lib = Self {
            site_version: top.site_version,
            entries: top.data,
            by_name: BTreeMap::new(),
            by_matnum: BTreeMap::new(),
        };
        for (i, e) in lib.entries.iter().enumerate() {
            lib.by_name.insert(e.name.to_ascii_lowercase(), i);
            lib.by_matnum.insert(e.mat_num, i);
        }
        Ok(lib)
    }

    /// Read and parse a compendium JSON file.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self, Error> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path)
            .map_err(|e| Error::Io(format!("{}: {}", path.display(), e)))?;
        Self::from_json(&text)
    }

    /// Number of materials.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Case-insensitive lookup by display name ("Air (dry, near sea level)").
    pub fn get(&self, name: &str) -> Option<&CompendiumEntry> {
        self.by_name
            .get(&name.to_ascii_lowercase())
            .map(|&i| &lib_entries(self)[i])
    }

    /// Lookup by compendium material number (1..=411).
    pub fn get_by_matnum(&self, num: u32) -> Option<&CompendiumEntry> {
        self.by_matnum.get(&num).map(|&i| &lib_entries(self)[i])
    }

    /// All display names in file order.
    pub fn names(&self) -> Vec<&str> {
        self.entries.iter().map(|e| e.name.as_str()).collect()
    }
}

// Borrow helper keeping the index maps and storage in one struct without
// split borrows surfacing in the public API.
fn lib_entries(lib: &MaterialsLibrary) -> &Vec<CompendiumEntry> {
    &lib.entries
}

#[cfg(test)]
mod tests {
    use super::*;

    fn library() -> MaterialsLibrary {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../fixtures/data/MaterialsCompendium.json"
        );
        MaterialsLibrary::from_file(path).unwrap()
    }

    #[test]
    fn loads_all_411_materials() {
        let lib = library();
        assert_eq!(lib.len(), 411);
        assert!(!lib.is_empty());
        assert!(!lib.site_version.is_empty());
    }

    #[test]
    fn names_unique_and_ordered() {
        let lib = library();
        let names = lib.names();
        let unique: std::collections::BTreeSet<&str> = names.iter().copied().collect();
        assert_eq!(unique.len(), 411);
        // Dataset order: B-110 Bone Equivalent Plastic is entry zero.
        assert_eq!(names[0], "Bone Equivalent Plastic, B-110");
    }

    #[test]
    fn lookup_by_name_case_insensitive() {
        let lib = library();
        let air = lib.get("air (DRY, near sea level)").expect("air present");
        assert_eq!(air.mat_num, 4); // early compendium entry
        assert!(air.density > 0.001);
    }

    #[test]
    fn lookup_by_matnum() {
        let lib = library();
        assert!(lib.get_by_matnum(1).is_some());
        assert!(lib.get_by_matnum(412).is_none());
    }

    #[test]
    fn weight_fractions_sum_near_one() {
        let lib = library();
        for entry in &lib.entries {
            let total: f64 = entry.weight_fractions().values().sum();
            assert!(
                (total - 1.0).abs() < 1e-2,
                "material `{}` fractions sum to {total}",
                entry.name
            );
        }
    }

    #[test]
    fn bone_plastic_h1_spot_value() {
        // From the dataset: B-110 Bone Equivalent Plastic, H1 wf = 0.035491
        let lib = library();
        let bone = lib.get("Bone Equivalent Plastic, B-110").unwrap();
        let wf = bone.weight_fractions();
        assert!((wf[&1_001] - 0.035491).abs() < 1e-6);
    }

    #[test]
    fn to_material_builds_named_composition_with_metadata() {
        let lib = library();
        let entry = lib.get("Acetone").unwrap();
        let mat = entry.to_material().unwrap();
        assert!(!mat.comp.is_empty());
        let meta = mat.metadata().expect("metadata attached");
        assert_eq!(meta["mat_num"], entry.mat_num);
        // Fractions are relative masses; atom_fractions must resolve via AME2020.
        let af = mat.atom_fractions(&crate::Ame2020).unwrap();
        assert_eq!(af.len(), mat.comp.len());
    }

    #[test]
    fn missing_file_errors() {
        let err = MaterialsLibrary::from_file("/nonexistent/compendium.json");
        assert!(matches!(err, Err(Error::Io(_))));
    }
}
