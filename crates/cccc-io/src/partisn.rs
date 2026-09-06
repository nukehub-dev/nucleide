//! PARTISN discrete-ordinates input writer.
//!
//! Emits a minimal PARTISN input in this documented subset:
//!
//! ```text
//! TITLE <title>
//! DIM <1|2|3>
//! ZONE <id> MATERIAL <material> DENSITY <density> ISOTXS <label>...
//! ... one ZONE card per zone ...
//! [SOURCE <source>]
//! END
//! ```
//!
//! Validation rules ([`PartisnDeck::validate`]):
//!
//! - `dim` must be 1, 2, or 3; anything else fails.
//! - Every label in every zone's `isotxs_labels` must exist in the referenced
//!   [`IsotxsLib`](crate::isotxs::IsotxsLib); a dangling label fails with a
//!   [`crate::error::Error::CrossRef`] error naming the zone and the label.

use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};
use crate::isotxs::IsotxsLib;

/// One PARTISN material zone with its ISOTXS nuclide mapping.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PartisnZone {
    /// Zone identifier.
    pub id: u32,
    /// Material name assigned to the zone.
    pub material: String,
    /// ISOTXS nuclide labels mapped into this zone.
    pub isotxs_labels: Vec<String>,
    /// Zone material density.
    pub density: f64,
}

/// Minimal PARTISN discrete-ordinates input deck.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct PartisnDeck {
    /// Problem title (emitted on the title card).
    pub title: String,
    /// Problem dimensionality: 1, 2, or 3.
    pub dim: u8,
    /// Material zones with ISOTXS nuclide mappings.
    pub zones: Vec<PartisnZone>,
    /// Optional isotropic volumetric source description.
    pub source: Option<String>,
}

impl PartisnDeck {
    /// Render the deck to PARTISN input text.
    pub fn render(&self) -> String {
        let mut out = String::new();
        out.push_str(&format!("TITLE {}\n", self.title));
        out.push_str(&format!("DIM {}\n", self.dim));
        for zone in &self.zones {
            out.push_str(&format!(
                "ZONE {} MATERIAL {} DENSITY {} ISOTXS {}\n",
                zone.id,
                zone.material,
                zone.density,
                zone.isotxs_labels.join(" ")
            ));
        }
        if let Some(source) = &self.source {
            out.push_str(&format!("SOURCE {source}\n"));
        }
        out.push_str("END\n");
        out
    }

    /// Check the deck against an ISOTXS library.
    ///
    /// Fails when `dim` is not 1, 2, or 3, or when a zone names an ISOTXS
    /// label absent from `isotxs` (reported as [`Error::CrossRef`]).
    pub fn validate(&self, isotxs: &IsotxsLib) -> Result<()> {
        if !(1..=3).contains(&self.dim) {
            return Err(Error::Parse {
                line: 1,
                msg: format!(
                    "PARTISN DIM {} out of range (expected 1, 2, or 3)",
                    self.dim
                ),
            });
        }
        for zone in &self.zones {
            for label in &zone.isotxs_labels {
                if isotxs.find(label).is_none() {
                    return Err(Error::CrossRef(format!(
                        "zone {} references unknown ISOTXS nuclide `{label}`",
                        zone.id
                    )));
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture_lib() -> IsotxsLib {
        IsotxsLib::from_file(format!(
            "{}/../../fixtures/cccc/isotxs_sample",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap()
    }

    fn sample_deck() -> PartisnDeck {
        PartisnDeck {
            title: "synthetic slab".to_string(),
            dim: 1,
            zones: vec![
                PartisnZone {
                    id: 1,
                    material: "fuel".to_string(),
                    isotxs_labels: vec!["U235".to_string()],
                    density: 10.0,
                },
                PartisnZone {
                    id: 2,
                    material: "blanket".to_string(),
                    isotxs_labels: vec!["PU239".to_string()],
                    density: 5.0,
                },
            ],
            source: Some("isotropic".to_string()),
        }
    }

    #[test]
    fn renders_title_dim_zones_source_and_end() {
        let text = sample_deck().render();
        assert!(text.contains("synthetic slab"));
        assert!(text.contains("DIM 1"));
        assert!(text.contains("ZONE 1"));
        assert!(text.contains("MATERIAL fuel"));
        assert!(text.contains("U235"));
        assert!(text.contains("ZONE 2"));
        assert!(text.contains("PU239"));
        assert!(text.contains("ISOTXS"));
        assert!(text.contains("SOURCE isotropic"));
        assert!(text.ends_with("END\n"));
    }

    #[test]
    fn renders_without_source_omits_source_card() {
        let mut deck = sample_deck();
        deck.source = None;
        let text = deck.render();
        assert!(!text.contains("SOURCE"));
        assert!(text.ends_with("END\n"));
    }

    #[test]
    fn validates_against_isotxs_fixture() {
        sample_deck().validate(&fixture_lib()).unwrap();
    }

    #[test]
    fn dangling_label_is_crossref() {
        let mut deck = sample_deck();
        deck.zones[0].isotxs_labels.push("U238".to_string());
        let err = deck.validate(&fixture_lib()).unwrap_err();
        match err {
            Error::CrossRef(msg) => {
                assert!(msg.contains("U238"));
                assert!(msg.contains('1'));
            }
            err => panic!("expected CrossRef error, got {err:?}"),
        }
    }

    #[test]
    fn bad_dim_is_rejected() {
        let lib = fixture_lib();
        for dim in [0, 4, 255] {
            let mut deck = sample_deck();
            deck.dim = dim;
            assert!(
                deck.validate(&lib).is_err(),
                "dim {dim} should not validate"
            );
        }
    }
}
