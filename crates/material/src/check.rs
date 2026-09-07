//! Composition label-collision and conservation checks.
//!
//! Pure functions over [`Material`]; no I/O. [`check_labels`] truncates
//! every member's cross-code labels to DIF3D/MC2-style field widths and
//! reports the truncated forms claimed by more than one nuclide.
//! [`audit`] validates mass-fraction normalization, sign, key uniqueness,
//! and atomic-mass availability.

use std::collections::{BTreeMap, BTreeSet};

use nucleide_nuclei::{armi, dialects, NuclideId};

use crate::{MassProvider, Material};

/// Default truncation widths probed by [`check_labels`]: 6 (DIF3D label
/// limit) and 8 (MC2-3 label limit).
pub const DEFAULT_WIDTHS: [usize; 2] = [6, 8];

/// Tolerance for the mass-fraction-sum check in [`audit`].
const FRACTION_TOL: f64 = 1e-9;

/// One truncated label form claimed by more than one nuclide.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Collision {
    /// The shared truncated form.
    pub truncated: String,
    /// The field width it was truncated to.
    pub width: usize,
    /// The colliding nuclides, sorted.
    pub members: Vec<NuclideId>,
}

/// The class of a conservation problem found by [`audit`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AuditKind {
    /// Normalized mass fractions do not sum to 1.0 (or no total exists).
    FractionsDontSum,
    /// A stored mass is negative or NaN.
    NegativeMass,
    /// The same nuclide id appears more than once.
    DuplicateNuclide,
    /// No atomic mass is available for a nuclide.
    UnknownMass,
}

/// One conservation problem found by [`audit`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AuditIssue {
    /// The class of problem.
    pub kind: AuditKind,
    /// Human-readable detail (nuclide, offending value).
    pub detail: String,
}

/// Cross-code labels emitted for one nuclide: GNDS name, MCNP ZAID,
/// Serpent, ALARA, and the ARMI database label.
fn cross_code_labels(id: NuclideId) -> [String; 5] {
    [
        id.to_name(),
        dialects::to_zaid(id).to_string(),
        dialects::serpent(id),
        dialects::alara(id),
        armi::nucid_to_armi_label(id),
    ]
}

/// Truncate every member's cross-code labels to each of `widths` and report
/// every truncated form claimed by more than one nuclide.
///
/// Results are ordered by `(width, truncated)`. Width `0` is skipped: it
/// maps every label to the empty string.
pub fn check_labels(mat: &Material, widths: &[usize]) -> Vec<Collision> {
    let mut claims: BTreeMap<(usize, String), BTreeSet<NuclideId>> = BTreeMap::new();
    for &id in mat.comp.keys() {
        for label in cross_code_labels(id) {
            for &width in widths {
                if width == 0 {
                    continue;
                }
                let truncated: String = label.chars().take(width).collect();
                claims.entry((width, truncated)).or_default().insert(id);
            }
        }
    }
    claims
        .into_iter()
        .filter(|(_, members)| members.len() > 1)
        .map(|((width, truncated), members)| Collision {
            truncated,
            width,
            members: members.into_iter().collect(),
        })
        .collect()
}

/// Validate mass-fraction normalization, mass signs, key uniqueness, and
/// atomic-mass availability.
///
/// - [`AuditKind::FractionsDontSum`]: the total mass is not positive-finite
///   (empty, zero, NaN, or infinite total, so no fractions exist) or the
///   normalized fractions sum outside `1.0 ± 1e-9`.
/// - [`AuditKind::NegativeMass`]: one issue per nuclide with a negative or
///   NaN stored mass.
/// - [`AuditKind::DuplicateNuclide`]: defensive only — `comp` is keyed by
///   [`NuclideId`], so duplicate keys cannot occur in memory; the variant
///   exists for file/FFI consumers sharing [`AuditKind`].
/// - [`AuditKind::UnknownMass`]: one issue per nuclide missing from
///   `masses`.
pub fn audit(mat: &Material, masses: &impl MassProvider) -> Vec<AuditIssue> {
    let mut issues = Vec::new();
    for (&id, &mass) in &mat.comp {
        if mass.is_nan() || mass < 0.0 {
            issues.push(AuditIssue {
                kind: AuditKind::NegativeMass,
                detail: format!("non-negative mass expected for {id}: {mass}"),
            });
        }
    }
    {
        let mut seen = BTreeSet::new();
        for &id in mat.comp.keys() {
            if !seen.insert(id.nucid()) {
                issues.push(AuditIssue {
                    kind: AuditKind::DuplicateNuclide,
                    detail: format!("duplicate nuclide {id}"),
                });
            }
        }
    }
    let total: f64 = mat.comp.values().sum();
    if !total.is_finite() || total <= 0.0 {
        issues.push(AuditIssue {
            kind: AuditKind::FractionsDontSum,
            detail: format!("no mass fractions: total mass is {total}"),
        });
    } else {
        let sum: f64 = mat.comp.values().map(|mass| mass / total).sum();
        if (sum - 1.0).abs() > FRACTION_TOL {
            issues.push(AuditIssue {
                kind: AuditKind::FractionsDontSum,
                detail: format!("mass fractions sum to {sum}, expected 1.0"),
            });
        }
    }
    for &id in mat.comp.keys() {
        if masses.mass(id.nucid()).is_none() {
            issues.push(AuditIssue {
                kind: AuditKind::UnknownMass,
                detail: format!("no atomic mass available for {id}"),
            });
        }
    }
    issues
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    fn id(name: &str) -> NuclideId {
        NuclideId::from_name(name).unwrap()
    }

    struct Table(HashMap<u32, f64>);

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

    fn full_table() -> Table {
        Table::new(&[
            ("U235", 235.0),
            ("U238", 238.0),
            ("Pu239", 239.0),
            ("Am242_m1", 242.0),
            ("Am242_m2", 242.0),
        ])
    }

    #[test]
    fn default_widths_cover_dif3d_and_mc2_limits() {
        assert_eq!(DEFAULT_WIDTHS, [6, 8]);
    }

    #[test]
    fn isomer_pair_collides_at_width_six() {
        let mut mat = Material::new();
        mat.add_nuclide(id("Am242_m1"), 1.0);
        mat.add_nuclide(id("Am242_m2"), 1.0);

        let hits = check_labels(&mat, &[6]);
        let gnds = hits
            .iter()
            .find(|c| c.truncated == "Am242_" && c.width == 6)
            .expect("Am242_m1/Am242_m2 share 6-char GNDS prefix Am242_");
        assert_eq!(gnds.members, vec![id("Am242_m1"), id("Am242_m2")]);
    }

    #[test]
    fn collisions_cover_default_widths_and_stay_ordered() {
        let mut mat = Material::new();
        mat.add_nuclide(id("Am242_m1"), 1.0);
        mat.add_nuclide(id("Am242_m2"), 1.0);

        let hits = check_labels(&mat, &DEFAULT_WIDTHS);
        assert!(!hits.is_empty());
        assert!(hits.iter().all(|c| c.members.len() > 1));
        let keys: Vec<(usize, &str)> = hits
            .iter()
            .map(|c| (c.width, c.truncated.as_str()))
            .collect();
        let mut sorted = keys.clone();
        sorted.sort();
        assert_eq!(keys, sorted, "collisions ordered by (width, truncated)");
        // ALARA drops state, so both isomers are literally "am:242".
        assert!(hits.iter().any(|c| c.truncated == "am:242"));
    }

    #[test]
    fn distinct_nuclides_do_not_collide() {
        let mut mat = Material::new();
        mat.add_nuclide(id("U235"), 19.0);
        mat.add_nuclide(id("U238"), 1.0);
        assert_eq!(check_labels(&mat, &DEFAULT_WIDTHS), vec![]);
        assert_eq!(check_labels(&Material::new(), &DEFAULT_WIDTHS), vec![]);

        let mut single = Material::new();
        single.add_nuclide(id("U235"), 1.0);
        assert_eq!(check_labels(&single, &[6, 8]), vec![]);
        assert_eq!(check_labels(&single, &[0]), vec![]);
    }

    #[test]
    fn clean_material_audits_empty() {
        let mut mat = Material::new();
        mat.add_nuclide(id("U235"), 19.0);
        mat.add_nuclide(id("U238"), 1.0);
        assert_eq!(audit(&mat, &full_table()), vec![]);
    }

    #[test]
    fn negative_and_nan_masses_are_flagged() {
        let mut mat = Material::new();
        mat.add_nuclide(id("U235"), -1.0);
        mat.add_nuclide(id("U238"), 2.0);
        let issues = audit(&mat, &full_table());
        assert_eq!(issues.len(), 1);
        assert_eq!(issues[0].kind, AuditKind::NegativeMass);
        assert!(issues[0].detail.contains("U235"));

        let mut nan = Material::new();
        nan.add_nuclide(id("Pu239"), f64::NAN);
        nan.add_nuclide(id("U235"), 1.0);
        let kinds: Vec<AuditKind> = audit(&nan, &full_table()).iter().map(|i| i.kind).collect();
        assert!(kinds.contains(&AuditKind::NegativeMass));
        assert!(kinds.contains(&AuditKind::FractionsDontSum));
    }

    #[test]
    fn empty_material_fractions_do_not_sum() {
        let issues = audit(&Material::new(), &full_table());
        assert_eq!(issues.len(), 1);
        assert_eq!(issues[0].kind, AuditKind::FractionsDontSum);
    }

    #[test]
    fn missing_masses_are_flagged() {
        let mut mat = Material::new();
        mat.add_nuclide(id("U235"), 1.0);
        mat.add_nuclide(id("Pu239"), 1.0);
        let partial = Table::new(&[("U235", 235.0)]);
        let issues = audit(&mat, &partial);
        assert_eq!(issues.len(), 1);
        assert_eq!(issues[0].kind, AuditKind::UnknownMass);
        assert!(issues[0].detail.contains("Pu239"));
    }

    #[test]
    fn duplicate_variant_unreachable_for_map_backed_material() {
        // `comp` is a BTreeMap: re-adding accumulates instead of duplicating.
        let mut mat = Material::new();
        mat.add_nuclide(id("U235"), 1.0);
        mat.add_nuclide(id("U235"), 2.0);
        assert_eq!(mat.comp.len(), 1);
        assert!(audit(&mat, &full_table())
            .iter()
            .all(|i| i.kind != AuditKind::DuplicateNuclide));
    }
}
