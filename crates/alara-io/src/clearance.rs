//! Clearance / waste-classification analytics over parsed inventories.
//!
//! Pure arithmetic on activation-code inventories: the **clearance index**
//!
//! ```text
//! CI = sum_i A_i / CL_i
//! ```
//!
//! over per-nuclide activities `A_i` and clearance levels `CL_i`, plus the
//! **sum-of-fractions rule** (IAEA Safety Standards Series No. RS-G-1.7 §5,
//! referenced by designation only): a material satisfies the clearance
//! screening criterion when `sum_i A_i / CL_i <= 1` (boundary included). The
//! methodology follows Sublet et al., *Nuclear Data Sheets* **139** (2017) 77
//! (FISPACT-II radiological indices); further indices from that paper stay
//! recorded, not implemented.
//!
//! This module is **screening arithmetic, not a compliance decision**: real
//! clearance requires the governing regulatory table, material bookkeeping,
//! and the national transposition of the underlying directive — see
//! [`ClearanceTable::eu_annex_vii`] for the provenance and limits of the
//! vendored default.
//!
//! Two further vendored sets accompany the EU default: the Spanish CSN
//! conditional-clearance tables for NORM landfill disposal (Tables 1–3 for
//! inert / non-hazardous / hazardous landfills, per NORM material nature,
//! with the Table 4 chain keys expanded to per-member entries — see
//! [`ClearanceTable::es_conditional_inert`]). The caller selects the table
//! explicitly; there is no cross-table logic and no "most permissive wins".
//!
//! Unit discipline is the caller's: `A_i` and `CL_i` must carry the same
//! basis. The vendored [`ClearanceTable::eu_annex_vii`] table is an
//! *activity-concentration* table (Bq/g), so inventories compared against it
//! must be massic activities (Bq/g), not total activities (Bq). Nuclide keys
//! are canonical [`NuclideId`]s resolved through the shared dialect
//! machinery — no second naming convention is introduced here.
//!
//! Inventories normally come from parsed [`crate::output::ResponseFrame`]
//! rows (see [`inventory_from_frame`]); any other source can supply the
//! same `&[(NuclideId, f64)]` shape directly.

use std::collections::BTreeMap;

use nucleide_nuclei::NuclideId;

use crate::error::{Error, Result};
use crate::output::{ResponseFrame, ResponseVar};

/// Vendored EU 2013/59/Euratom Annex VII Table A transcription (Bq/g).
///
/// Transcribed from the official legal text (EUR-Lex CELEX:32013L0059,
/// OJ L 13, 17.1.2014, Annex VII Table A, "Activity concentrations and
/// activities per unit mass of radionuclides ... for solid materials"),
/// accessed 2026-09-15. EU legal text is reusable with attribution under
/// Commission Implementing Decision 2011/833/EU. Columns: `GNDS name`,
/// `limit Bq/g`. Plain published facts (nuclide, number); never vendor IAEA
/// tables (RS-G-1.7, GSG-17) — reference them by designation only.
const EU_ANNEX_VII_A_TSV: &str = include_str!("data/eu_annex_vii_a.tsv");

/// Vendored Spanish CSN natural decay-chain definitions (Tabla 4).
///
/// Transcribed from the Consejo de Seguridad Nuclear draft technical opinion
/// `CSN/PDT/AICD/TGE/2503/02` (`TGE/VAR/2025/1`), "Propuesta de dictamen
/// técnico para el informe favorable de niveles de desclasificación para la
/// gestión en vertedero de material radiactivo de origen natural (NORM)",
/// hosted on csn.es and accessed 2026-09-16; official regulatory text under
/// RD 1029/2022 and RD 1217/2024 (RINR), both transposing Directive
/// 2013/59/Euratom. Columns: `chain key`, `member1,member2,...` in source
/// order. The source's branching annotations (e.g. `Pa-234 (0.3%)`) are
/// secular-equilibrium composition notes, not part of the names; every
/// listed isotope is a member.
const ES_NORM_CHAINS_TSV: &str = include_str!("data/es_norm_chains.tsv");

/// Vendored Spanish CSN conditional NORM clearance levels for a landfill of
/// inert waste (Tabla 1, "Vertedero de residuos inertes"), Bq/g.
///
/// Same source as [`ES_NORM_CHAINS_TSV`] (its Tabla 1), transcribed with
/// attribution. Columns: `chain key`, then the source columns in source
/// order `rocas` (rocks), `cenizas` (ashes), `arenas` (sands), `escorias`
/// (slags), `gas/petroleo` (oil & gas NORM waste). Chain keys expand per
/// [`ES_NORM_CHAINS_TSV`]: the key's level applies to each member in secular
/// equilibrium, the EU Part-2 precedent.
const ES_CSN_INERT_TSV: &str = include_str!("data/es_csn_inert.tsv");

/// Vendored Spanish CSN conditional NORM clearance levels for a landfill of
/// non-hazardous waste (Tabla 2, "Vertedero de residuos no peligrosos"),
/// Bq/g. Same shape and source as [`ES_CSN_INERT_TSV`] (its Tabla 2).
const ES_CSN_NON_HAZARDOUS_TSV: &str = include_str!("data/es_csn_non_hazardous.tsv");

/// Vendored Spanish CSN conditional NORM clearance levels for a landfill of
/// hazardous waste (Tabla 3, "Vertedero de residuos peligrosos"), Bq/g.
/// Same shape and source as [`ES_CSN_INERT_TSV`] (its Tabla 3).
const ES_CSN_HAZARDOUS_TSV: &str = include_str!("data/es_csn_hazardous.tsv");

/// NORM waste material natures of the Spanish CSN landfill clearance columns.
///
/// The CSN Tables 1–3 give one clearance level per chain key for each
/// material nature a NORM waste can take; the caller picks the column that
/// matches the waste under screening. Variant order matches the source
/// column order (ROCAS, CENIZAS, ARENAS, ESCORIAS, GAS/PETROLEO).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EsNormMaterial {
    /// `ROCAS` — rocks.
    Rocks,
    /// `CENIZAS` — ashes.
    Ashes,
    /// `ARENAS` — sands.
    Sands,
    /// `ESCORIAS` — slags.
    Slags,
    /// `GAS/PETROLEO` — oil & gas industry NORM waste.
    OilGas,
}

impl EsNormMaterial {
    /// All material natures in source column order (index == TSV column).
    const ALL: [Self; 5] = [
        Self::Rocks,
        Self::Ashes,
        Self::Sands,
        Self::Slags,
        Self::OilGas,
    ];

    /// Column index of this material nature in the embedded CSN table TSVs.
    fn column(self) -> usize {
        match self {
            Self::Rocks => 0,
            Self::Ashes => 1,
            Self::Sands => 2,
            Self::Slags => 3,
            Self::OilGas => 4,
        }
    }
}

/// Per-nuclide clearance-level table used by [`clearance_index`] and
/// [`sum_of_fractions`].
///
/// Keys are canonical [`NuclideId`]s; values are the clearance levels in the
/// unit basis the caller is working in (Bq total, Bq/g, ...). Build one with
/// [`ClearanceTable::insert`] or take the vendored
/// [`ClearanceTable::eu_annex_vii`] default.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ClearanceTable {
    limits: BTreeMap<NuclideId, f64>,
}

impl ClearanceTable {
    /// Empty caller-supplied table.
    pub fn new() -> Self {
        Self::default()
    }

    /// Insert or replace one nuclide's clearance level.
    ///
    /// `limit` must be finite and strictly positive. Errors:
    /// [`Error::BadClearanceValue`] on `limit <= 0.0` or non-finite input.
    pub fn insert(&mut self, nuclide: NuclideId, limit: f64) -> Result<()> {
        if !limit.is_finite() || limit <= 0.0 {
            return Err(Error::BadClearanceValue {
                nuclide: nuclide.to_name(),
                msg: format!("clearance limit must be finite and > 0, got {limit}"),
            });
        }
        self.limits.insert(nuclide, limit);
        Ok(())
    }

    /// Clearance level for `nuclide`, or `None` when the table carries none.
    pub fn get(&self, nuclide: NuclideId) -> Option<f64> {
        self.limits.get(&nuclide).copied()
    }

    /// Number of table entries.
    pub fn len(&self) -> usize {
        self.limits.len()
    }

    /// True when the table has no entries.
    pub fn is_empty(&self) -> bool {
        self.limits.is_empty()
    }

    /// Iterate entries in canonical nucid order.
    pub fn iter(&self) -> impl Iterator<Item = (NuclideId, f64)> + '_ {
        self.limits.iter().map(|(nuc, limit)| (*nuc, *limit))
    }

    /// Vendored default: EU 2013/59/Euratom Annex VII Table A, activity
    /// concentrations for solid material in Bq/g.
    ///
    /// Provenance: official legal text transcribed 2026-09-15 from EUR-Lex
    /// CELEX:32013L0059 (Council Directive 2013/59/Euratom of 5 December
    /// 2013, OJ L 13, 17.1.2014, p. 1), Annex VII Table A, consolidated
    /// table version; reusable with attribution per Commission Implementing
    /// Decision (EU) 2011/833/EU. The transcription is a plain-facts table
    /// (`name`, `Bq/g`); see the module docs for the unit-basis contract.
    ///
    /// The table is embedded at build time and parsed lazily on first call.
    /// The transcription is committed data (like the `nuclei` static tables)
    /// pinned by the `eu_table_row_count` test, so this constructor treats a
    /// corrupt transcription as unreachable-in-practice and panics with the
    /// offending line; fallible callers use [`Self::try_eu_annex_vii`].
    pub fn eu_annex_vii() -> Self {
        static TABLE: std::sync::OnceLock<ClearanceTable> = std::sync::OnceLock::new();
        TABLE
            .get_or_init(|| {
                Self::try_eu_annex_vii().unwrap_or_else(|e| panic!("eu_annex_vii_a.tsv: {e}"))
            })
            .clone()
    }

    /// Fallible parse of the embedded EU Annex VII Table A transcription.
    ///
    /// Same data as [`Self::eu_annex_vii`] without the lazy cache: every
    /// malformed line (missing tab, unknown nuclide, bad limit) is a loud
    /// [`Error::Parse`] carrying the 1-based TSV line number.
    pub fn try_eu_annex_vii() -> Result<Self> {
        let mut table = ClearanceTable::new();
        for (index, line) in EU_ANNEX_VII_A_TSV.lines().enumerate() {
            let line_no = index + 1;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let (name, value_text) = trimmed.split_once('\t').ok_or_else(|| Error::Parse {
                line: line_no,
                msg: format!("expected `name<TAB>Bq/g`, found `{trimmed}`"),
            })?;
            let nuclide = NuclideId::from_name(name.trim()).map_err(|_| Error::Parse {
                line: line_no,
                msg: format!("unknown nuclide `{name}`"),
            })?;
            let limit: f64 = value_text.trim().parse().map_err(|_| Error::Parse {
                line: line_no,
                msg: format!("bad limit `{value_text}`"),
            })?;
            table.insert(nuclide, limit).map_err(|e| Error::Parse {
                line: line_no,
                msg: format!("{e}"),
            })?;
        }
        Ok(table)
    }

    /// Vendored Spanish CSN table for a landfill of **inert** waste
    /// (Tabla 1), activity concentrations in Bq/g.
    ///
    /// Selects the per-chain-key levels for one NORM material nature
    /// ([`EsNormMaterial`]) and expands the Table 4 chain keys to per-member
    /// entries at the parent value. Provenance and flattening contract: see
    /// [`Self::try_es_conditional_inert`]; same `try_` + cache split as
    /// [`Self::eu_annex_vii`].
    pub fn es_conditional_inert(material: EsNormMaterial) -> Self {
        static TABLES: std::sync::OnceLock<[ClearanceTable; 5]> = std::sync::OnceLock::new();
        TABLES.get_or_init(|| Self::parse_es_family(ES_CSN_INERT_TSV, "es_csn_inert.tsv"))
            [material.column()]
        .clone()
    }

    /// Fallible parse of the embedded inert-landfill table (Tabla 1).
    ///
    /// Same data as [`Self::es_conditional_inert`] without the lazy cache:
    /// every malformed line (wrong field count, unknown chain key, bad
    /// limit, unresolvable chain member) is a loud [`Error::Parse`] carrying
    /// the 1-based TSV line number.
    pub fn try_es_conditional_inert(material: EsNormMaterial) -> Result<Self> {
        Self::try_es_table(ES_CSN_INERT_TSV, "es_csn_inert.tsv", material)
    }

    /// Vendored Spanish CSN table for a landfill of **non-hazardous** waste
    /// (Tabla 2), activity concentrations in Bq/g.
    ///
    /// Same contract as [`Self::es_conditional_inert`].
    pub fn es_conditional_non_hazardous(material: EsNormMaterial) -> Self {
        static TABLES: std::sync::OnceLock<[ClearanceTable; 5]> = std::sync::OnceLock::new();
        TABLES.get_or_init(|| {
            Self::parse_es_family(ES_CSN_NON_HAZARDOUS_TSV, "es_csn_non_hazardous.tsv")
        })[material.column()]
        .clone()
    }

    /// Fallible parse of the embedded non-hazardous-landfill table (Tabla 2).
    ///
    /// Same contract as [`Self::try_es_conditional_inert`].
    pub fn try_es_conditional_non_hazardous(material: EsNormMaterial) -> Result<Self> {
        Self::try_es_table(
            ES_CSN_NON_HAZARDOUS_TSV,
            "es_csn_non_hazardous.tsv",
            material,
        )
    }

    /// Vendored Spanish CSN table for a landfill of **hazardous** waste
    /// (Tabla 3), activity concentrations in Bq/g.
    ///
    /// Same contract as [`Self::es_conditional_inert`].
    pub fn es_conditional_hazardous(material: EsNormMaterial) -> Self {
        static TABLES: std::sync::OnceLock<[ClearanceTable; 5]> = std::sync::OnceLock::new();
        TABLES.get_or_init(|| Self::parse_es_family(ES_CSN_HAZARDOUS_TSV, "es_csn_hazardous.tsv"))
            [material.column()]
        .clone()
    }

    /// Fallible parse of the embedded hazardous-landfill table (Tabla 3).
    ///
    /// Same contract as [`Self::try_es_conditional_inert`].
    pub fn try_es_conditional_hazardous(material: EsNormMaterial) -> Result<Self> {
        Self::try_es_table(ES_CSN_HAZARDOUS_TSV, "es_csn_hazardous.tsv", material)
    }

    /// Parse one embedded CSN landfill table TSV for one material nature.
    ///
    /// Rows are chain keys (validated against the Tabla 4 transcription) with
    /// five per-material columns in source order; the selected column's level
    /// expands to every chain member at the parent value (the EU Part-2
    /// precedent). Rows apply in source order: a nuclide claimed by several
    /// keys keeps the **last** claiming key's value, i.e. the most specific
    /// subchain listed for it (the source lists each secular-equilibrium
    /// chain before its subchains). Screening against an
    /// equilibrium-chain characterization (CSN unit rule over chain
    /// precursors) is a caller-side table build, not this flattening.
    fn try_es_table(tsv: &str, file: &str, material: EsNormMaterial) -> Result<Self> {
        const COLUMNS: usize = 5;
        let chains = Self::parse_es_norm_chains()?;
        let column = material.column();
        let mut table = ClearanceTable::new();
        let mut seen_keys = std::collections::BTreeSet::new();
        for (index, line) in tsv.lines().enumerate() {
            let line_no = index + 1;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = trimmed.split('\t').collect();
            if fields.len() != COLUMNS + 1 {
                return Err(Error::Parse {
                    line: line_no,
                    msg: format!(
                        "{file}: expected `chain key` plus {COLUMNS} material columns, found {} fields",
                        fields.len()
                    ),
                });
            }
            let key = fields[0].trim();
            if !seen_keys.insert(key.to_string()) {
                return Err(Error::Parse {
                    line: line_no,
                    msg: format!("{file}: duplicate chain key `{key}`"),
                });
            }
            let members = chains.get(key).ok_or_else(|| Error::Parse {
                line: line_no,
                msg: format!("{file}: unknown chain key `{key}`"),
            })?;
            let limit: f64 = fields[1 + column]
                .trim()
                .parse()
                .map_err(|_| Error::Parse {
                    line: line_no,
                    msg: format!("{}: bad limit `{}`", file, fields[1 + column]),
                })?;
            for member in members {
                table.insert(*member, limit).map_err(|e| Error::Parse {
                    line: line_no,
                    msg: format!("{file}: {e}"),
                })?;
            }
        }
        Ok(table)
    }

    /// Parse the embedded Tabla 4 transcription into chain key -> members.
    ///
    /// Members resolve through the shared [`nucleide_nuclei`] dialect
    /// machinery (dashed spellings such as `Pa-234m` parse directly); a
    /// member that does not resolve is a loud [`Error::Parse`] carrying the
    /// 1-based TSV line number.
    fn parse_es_norm_chains() -> Result<BTreeMap<String, Vec<NuclideId>>> {
        let mut chains = BTreeMap::new();
        for (index, line) in ES_NORM_CHAINS_TSV.lines().enumerate() {
            let line_no = index + 1;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let (key, members_text) = trimmed.split_once('\t').ok_or_else(|| Error::Parse {
                line: line_no,
                msg: format!(
                    "es_norm_chains.tsv: expected `chain key<TAB>members`, found `{trimmed}`"
                ),
            })?;
            let mut members = Vec::new();
            for name in members_text.split(',') {
                let name = name.trim();
                let nuclide = NuclideId::from_name(name).map_err(|_| Error::Parse {
                    line: line_no,
                    msg: format!("es_norm_chains.tsv: unknown chain member `{name}`"),
                })?;
                members.push(nuclide);
            }
            if members.is_empty() {
                return Err(Error::Parse {
                    line: line_no,
                    msg: format!("es_norm_chains.tsv: chain key `{key}` has no members"),
                });
            }
            if chains.insert(key.to_string(), members).is_some() {
                return Err(Error::Parse {
                    line: line_no,
                    msg: format!("es_norm_chains.tsv: duplicate chain key `{key}`"),
                });
            }
        }
        Ok(chains)
    }

    /// Parse all five material columns of one CSN landfill table (cached
    /// constructor backend; panics with the offending line on a corrupt
    /// transcription, unreachable-in-practice per the row-count tests).
    fn parse_es_family(tsv: &'static str, file: &'static str) -> [ClearanceTable; 5] {
        std::array::from_fn(|column| {
            let material = EsNormMaterial::ALL[column];
            Self::try_es_table(tsv, file, material).unwrap_or_else(|e| panic!("{file}: {e}"))
        })
    }
}

/// Sum-of-fractions screening outcome (RS-G-1.7 §5, referenced by designation).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClearanceClass {
    /// `sum_i A_i / CL_i <= 1`: the inventory satisfies the sum-of-fractions
    /// screening criterion. The boundary `== 1` lands here ("does not
    /// exceed"), matching the exemption/clearance wording.
    Satisfied,
    /// `sum_i A_i / CL_i > 1`: the inventory does not satisfy the criterion.
    Exceeded,
}

impl std::fmt::Display for ClearanceClass {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Satisfied => "satisfied",
            Self::Exceeded => "exceeded",
        })
    }
}

/// Result of [`sum_of_fractions`]: the fraction sum, its screening class,
/// and the single dominant contributor (largest per-nuclide fraction).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SumOfFractions {
    /// `sum_i A_i / CL_i` over the inventory.
    pub sum: f64,
    /// Screening classification of `sum`.
    pub class: ClearanceClass,
    /// Largest per-nuclide fraction `A_i / CL_i`.
    pub max_fraction: f64,
    /// Nuclide carrying [`Self::max_fraction`]; `None` for an empty inventory
    /// (where [`Self::max_fraction`] is `0.0`).
    pub max_nuclide: Option<NuclideId>,
}

/// Clearance index `CI = sum_i A_i / CL_i` over an inventory.
///
/// `inventory` holds `(nuclide, activity)` pairs; `table` carries the
/// clearance levels in the same unit basis (see the module contract). Every
/// inventory nuclide must have a table entry and every activity must be
/// finite and `>= 0` — violations are loud [`Error`]s, never silent skips.
/// Duplicate nuclides are the caller's problem to avoid; they sum like any
/// other pair.
pub fn clearance_index(inventory: &[(NuclideId, f64)], table: &ClearanceTable) -> Result<f64> {
    Ok(fraction_sum(inventory, table)?.sum)
}

/// Sum-of-fractions screening: the fraction sum, its class (`<= 1` passes,
/// boundary included), and the dominant nuclide.
///
/// Same input contract as [`clearance_index`].
pub fn sum_of_fractions(
    inventory: &[(NuclideId, f64)],
    table: &ClearanceTable,
) -> Result<SumOfFractions> {
    fraction_sum(inventory, table)
}

/// Shared accumulation for [`clearance_index`] / [`sum_of_fractions`].
fn fraction_sum(inventory: &[(NuclideId, f64)], table: &ClearanceTable) -> Result<SumOfFractions> {
    let mut sum = 0.0f64;
    let mut max_fraction = 0.0f64;
    let mut max_nuclide = None;
    for (nuclide, activity) in inventory {
        if !activity.is_finite() || *activity < 0.0 {
            return Err(Error::BadClearanceValue {
                nuclide: nuclide.to_name(),
                msg: format!("activity must be finite and >= 0, got {activity}"),
            });
        }
        let limit = table
            .get(*nuclide)
            .ok_or_else(|| Error::MissingClearanceLimit {
                nuclide: nuclide.to_name(),
                table_len: table.len(),
            })?;
        let fraction = activity / limit;
        sum += fraction;
        if fraction > max_fraction {
            max_fraction = fraction;
            max_nuclide = Some(*nuclide);
        }
    }
    // Finite inputs can still overflow the sum (huge activity over a tiny
    // limit): an infinite total is a loud error, never `Ok(inf)`.
    if !sum.is_finite() {
        return Err(Error::BadClearanceValue {
            nuclide: "total".to_string(),
            msg: format!(
                "fraction sum overflowed to {sum} over {} entries",
                inventory.len()
            ),
        });
    }
    // `inventory` empty: no dominant nuclide (`max_fraction` is 0.0).
    let class = if sum <= 1.0 {
        ClearanceClass::Satisfied
    } else {
        ClearanceClass::Exceeded
    };
    Ok(SumOfFractions {
        sum,
        class,
        max_fraction,
        max_nuclide,
    })
}

/// Extract `(nuclide, activity)` inventory pairs from one cooling time of a
/// parsed [`ResponseFrame`].
///
/// Keeps [`ResponseVar::SpecificActivity`] rows at `time_s` and drops `total`
/// aggregate rows. Rows whose nuclide name does not resolve through the
/// shared dialect machinery are a loud [`Error::Parse`] (parser-validated
/// frames cannot produce them; hand-built frames can).
///
/// Unit note: the returned activities are in the frame's own activity unit
/// (ALARA `Bq/cm3` vs FISPACT-II `Bq`, see `var_unit`); the caller must pick
/// a [`ClearanceTable`] with the matching basis.
pub fn inventory_from_frame(frame: &ResponseFrame, time_s: f64) -> Result<Vec<(NuclideId, f64)>> {
    let mut pairs = Vec::new();
    for row in &frame.rows {
        if row.variable != ResponseVar::SpecificActivity || row.time_s != time_s {
            continue;
        }
        if row.is_total() {
            continue;
        }
        let nuclide = row.nuclide_id().map_err(|e| Error::Parse {
            line: 0,
            msg: format!("frame row nuclide `{}` does not resolve: {e}", row.nuclide),
        })?;
        pairs.push((nuclide, row.value));
    }
    Ok(pairs)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Hand-built three-nuclide table (synthetic, round numbers).
    fn toy_table() -> ClearanceTable {
        let mut table = ClearanceTable::new();
        table
            .insert(NuclideId::from_name("Co60").unwrap(), 10.0)
            .unwrap();
        table
            .insert(NuclideId::from_name("H3").unwrap(), 5.0)
            .unwrap();
        table
            .insert(NuclideId::from_name("Fe55").unwrap(), 2.0)
            .unwrap();
        table
    }

    #[test]
    fn hand_computed_ci_vectors_at_exact_equality() {
        let table = toy_table();
        // Single nuclide at its own limit: CI == 1 exactly.
        let inv = [(NuclideId::from_name("Co60").unwrap(), 10.0)];
        assert_eq!(clearance_index(&inv, &table).unwrap(), 1.0);
        // Two nuclides at half limit each: 0.5 + 0.5 == 1 exactly.
        let inv = [
            (NuclideId::from_name("Co60").unwrap(), 5.0),
            (NuclideId::from_name("H3").unwrap(), 2.5),
        ];
        assert_eq!(clearance_index(&inv, &table).unwrap(), 1.0);
        // Mixed vector: 10/10 + 5/5 + 4/2 == 4 exactly.
        let inv = [
            (NuclideId::from_name("Co60").unwrap(), 10.0),
            (NuclideId::from_name("H3").unwrap(), 5.0),
            (NuclideId::from_name("Fe55").unwrap(), 4.0),
        ];
        assert_eq!(clearance_index(&inv, &table).unwrap(), 4.0);
        // Zero activity contributes nothing.
        let inv = [(NuclideId::from_name("Fe55").unwrap(), 0.0)];
        assert_eq!(clearance_index(&inv, &table).unwrap(), 0.0);
        // Empty inventory: CI == 0.
        assert_eq!(clearance_index(&[], &table).unwrap(), 0.0);
    }

    #[test]
    fn sum_of_fractions_boundary_probes_both_sides() {
        let table = toy_table();
        // Exactly == 1 on the boundary: satisfied (boundary included).
        let inv = [(NuclideId::from_name("Co60").unwrap(), 10.0)];
        let out = sum_of_fractions(&inv, &table).unwrap();
        assert_eq!(out.sum, 1.0);
        assert_eq!(out.class, ClearanceClass::Satisfied);
        assert_eq!(out.max_fraction, 1.0);
        assert_eq!(out.max_nuclide, Some(NuclideId::from_name("Co60").unwrap()));

        // One ulp above the boundary: exceeded.
        let inv = [(NuclideId::from_name("Co60").unwrap(), 10.000000000000002)];
        let out = sum_of_fractions(&inv, &table).unwrap();
        assert!(out.sum > 1.0);
        assert_eq!(out.class, ClearanceClass::Exceeded);

        // One ulp below the boundary: satisfied.
        let inv = [(NuclideId::from_name("Co60").unwrap(), 9.999999999999998)];
        let out = sum_of_fractions(&inv, &table).unwrap();
        assert!(out.sum < 1.0);
        assert_eq!(out.class, ClearanceClass::Satisfied);

        // Fraction accumulation both sides: 0.4 + 0.5 = 0.9 vs 0.6 + 0.5 = 1.1.
        let below = [
            (NuclideId::from_name("Co60").unwrap(), 4.0),
            (NuclideId::from_name("H3").unwrap(), 2.5),
        ];
        assert_eq!(
            sum_of_fractions(&below, &table).unwrap().class,
            ClearanceClass::Satisfied
        );
        let above = [
            (NuclideId::from_name("Co60").unwrap(), 6.0),
            (NuclideId::from_name("H3").unwrap(), 2.5),
        ];
        let out = sum_of_fractions(&above, &table).unwrap();
        assert_eq!(out.class, ClearanceClass::Exceeded);
        assert_eq!(out.max_nuclide, Some(NuclideId::from_name("Co60").unwrap()));
        assert!(out.max_fraction > out.sum - out.max_fraction);
    }

    #[test]
    fn eu_table_loads_with_pinned_entry_count() {
        // Pins the committed transcription: any edit to the TSV must update
        // this count deliberately, which keeps `eu_annex_vii()` honest about
        // its unreachable-in-practice panic.
        let table = ClearanceTable::try_eu_annex_vii().unwrap();
        assert_eq!(table.len(), 260);
        assert_eq!(ClearanceTable::eu_annex_vii().len(), 260);
        // Spot check: Co-60 Table A value is 0.1 Bq/g.
        let co60 = table.get(NuclideId::from_name("Co60").unwrap()).unwrap();
        assert_eq!(co60, 0.1);
    }

    #[test]
    fn overflowing_fraction_sum_is_a_loud_error() {
        // Finite inputs whose sum overflows: loud, never Ok(inf).
        // f64::MAX over a 1e-308 limit overflows the single fraction to inf.
        let mut table = ClearanceTable::new();
        table
            .insert(NuclideId::from_name("Co60").unwrap(), 1e-308)
            .unwrap();
        let inv = [(NuclideId::from_name("Co60").unwrap(), f64::MAX)];
        match clearance_index(&inv, &table) {
            Err(Error::BadClearanceValue { nuclide, msg }) => {
                assert_eq!(nuclide, "total");
                assert!(msg.contains("overflowed"), "msg was `{msg}`");
            }
            other => panic!("expected overflow error, got {other:?}"),
        }
    }

    #[test]
    fn missing_limit_and_bad_values_are_loud_errors() {
        let table = toy_table();
        let unknown = [(NuclideId::from_name("Mn54").unwrap(), 1.0)];
        match clearance_index(&unknown, &table) {
            Err(Error::MissingClearanceLimit { nuclide, table_len }) => {
                assert_eq!(nuclide, "Mn54");
                assert_eq!(table_len, 3);
            }
            other => panic!("expected MissingClearanceLimit, got {other:?}"),
        }
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            let inv = [(NuclideId::from_name("Co60").unwrap(), bad)];
            assert!(
                matches!(
                    clearance_index(&inv, &table),
                    Err(Error::BadClearanceValue { .. })
                ),
                "activity {bad} should fail"
            );
        }
        let mut bad_table = ClearanceTable::new();
        for bad in [0.0, -2.0, f64::NAN] {
            assert!(
                matches!(
                    bad_table.insert(NuclideId::from_name("Co60").unwrap(), bad),
                    Err(Error::BadClearanceValue { .. })
                ),
                "limit {bad} should fail"
            );
        }
        assert!(bad_table.is_empty());
        bad_table
            .insert(NuclideId::from_name("Co60").unwrap(), 1.0)
            .unwrap();
        assert_eq!(bad_table.len(), 1);
    }

    #[test]
    fn frame_inventory_extracts_specific_activity_at_one_time() {
        let text = "*** Specific Activity [Bq/cm3] ***\n\
                     Interval #1 (Zone: inner) :\n\
                     isotope  t_1/2(s)   shutdown      1 d\n\
                     =====\n\
                     co-60 \t1.6636e+08  4.0000e+01  2.0000e+01\n\
                     h-3 \t3.8881e+08   1.0000e+01  5.0000e+00\n\
                     =====\n                     total   0           5.0000e+01  2.5000e+01\n";
        let frame = ResponseFrame::parse(text, "r").unwrap();
        let pairs = inventory_from_frame(&frame, 86_400.0).unwrap();
        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0], (NuclideId::from_name("Co60").unwrap(), 20.0));
        assert_eq!(pairs[1], (NuclideId::from_name("H3").unwrap(), 5.0));
        // Shutdown time: other values.
        let pairs = inventory_from_frame(&frame, 0.0).unwrap();
        assert_eq!(pairs[0].1, 40.0);
        // A time with no rows yields an empty (valid) inventory.
        assert!(inventory_from_frame(&frame, -1.0).unwrap().is_empty());
    }

    #[test]
    fn eu_annex_vii_default_table_loads_and_covers_classic_clearance_nuclides() {
        let table = ClearanceTable::eu_annex_vii();
        assert!(!table.is_empty());
        // Spot values transcribed from the directive (Bq/g, solid materials):
        // H-3 -> 100, C-14 -> 1, Co-60 -> 0.1, Cs-137 -> 0.1, Sr-90 -> 1,
        // Pu-239 -> 0.1; Part 2: U-238 -> 1, Th-232 -> 1, K-40 -> 10.
        for (name, expected) in [
            ("H3", 100.0),
            ("C14", 1.0),
            ("Co60", 0.1),
            ("Cs137", 0.1),
            ("Sr90", 1.0),
            ("Pu239", 0.1),
            ("U238", 1.0),
            ("Th232", 1.0),
            ("K40", 10.0),
        ] {
            let nuc = NuclideId::from_name(name).unwrap();
            assert_eq!(
                table.get(nuc),
                Some(expected),
                "{name} limit should be {expected} Bq/g"
            );
        }
        // All entries positive; canonical-keyed.
        assert!(table.iter().all(|(_, limit)| limit > 0.0));
        assert_eq!(table, ClearanceTable::eu_annex_vii());
    }

    /// Data rows (non-comment, non-blank) of an embedded TSV.
    fn tsv_data_rows(tsv: &str) -> usize {
        tsv.lines()
            .filter(|line| {
                let trimmed = line.trim();
                !trimmed.is_empty() && !trimmed.starts_with('#')
            })
            .count()
    }

    #[test]
    fn es_chains_transcribed_with_pinned_row_count_and_members() {
        // Pins the committed Tabla 4 transcription (15 chain keys).
        assert_eq!(tsv_data_rows(ES_NORM_CHAINS_TSV), 15);
        let chains = ClearanceTable::parse_es_norm_chains().unwrap();
        assert_eq!(chains.len(), 15);
        // Member counts straight from the source table.
        assert_eq!(chains["U-238sec"].len(), 15);
        assert_eq!(chains["U-nat"].len(), 7);
        assert_eq!(chains["Ra-226+"].len(), 6);
        assert_eq!(chains["U-235sec"].len(), 13);
        assert_eq!(chains["Ac-227+"].len(), 10);
        assert_eq!(chains["Th-232sec"].len(), 11);
        assert_eq!(chains["Th-228+"].len(), 8);
        assert_eq!(chains["K-40"].len(), 1);
        // Metastable and short-lived members resolve through the shared
        // dialect machinery (dashed source spellings).
        for name in ["Pa-234m", "Rn-222", "Po-218", "Fr-223", "Tl-207", "Po-211"] {
            let nuc = NuclideId::from_name(name).unwrap();
            assert!(
                chains.values().flatten().any(|member| *member == nuc),
                "{name} should appear as a chain member"
            );
        }
    }

    #[test]
    fn es_tables_load_with_pinned_entry_counts() {
        // Each table TSV pins 15 chain-key rows; the flattened union of the
        // Tabla 4 members is 40 nuclides for every landfill type/material.
        assert_eq!(tsv_data_rows(ES_CSN_INERT_TSV), 15);
        assert_eq!(tsv_data_rows(ES_CSN_NON_HAZARDOUS_TSV), 15);
        assert_eq!(tsv_data_rows(ES_CSN_HAZARDOUS_TSV), 15);
        for table in [
            ClearanceTable::try_es_conditional_inert(EsNormMaterial::Rocks).unwrap(),
            ClearanceTable::try_es_conditional_non_hazardous(EsNormMaterial::Ashes).unwrap(),
            ClearanceTable::try_es_conditional_hazardous(EsNormMaterial::OilGas).unwrap(),
        ] {
            assert_eq!(table.len(), 40);
            assert!(table.iter().all(|(_, limit)| limit > 0.0));
        }
        // Cached constructors agree with the fallible ones.
        for material in EsNormMaterial::ALL {
            assert_eq!(
                ClearanceTable::es_conditional_inert(material),
                ClearanceTable::try_es_conditional_inert(material).unwrap()
            );
            assert_eq!(
                ClearanceTable::es_conditional_non_hazardous(material),
                ClearanceTable::try_es_conditional_non_hazardous(material).unwrap()
            );
            assert_eq!(
                ClearanceTable::es_conditional_hazardous(material),
                ClearanceTable::try_es_conditional_hazardous(material).unwrap()
            );
        }
    }

    #[test]
    fn es_table_rejects_non_positive_limits_and_duplicate_keys() {
        // A negative or non-finite limit must fail loudly at parse (never a
        // silent screening depression); a duplicated chain-key row is
        // transcription corruption, not the documented last-wins overlap
        // across distinct keys.
        let base = ES_CSN_INERT_TSV
            .lines()
            .find(|l| !l.trim().is_empty() && !l.trim().starts_with('#'))
            .unwrap();
        let neg = ES_CSN_INERT_TSV.replacen(base, "K-40\t-1\t-1\t-1\t-1\t-1", 1);
        assert!(ClearanceTable::try_es_table(&neg, "es_neg.tsv", EsNormMaterial::Rocks).is_err());
        let dup = format!("{ES_CSN_INERT_TSV}\n{base}\n");
        assert!(ClearanceTable::try_es_table(&dup, "es_dup.tsv", EsNormMaterial::Rocks).is_err());
    }

    #[test]
    fn es_transcription_spots_per_table_and_material() {
        // Hand-read from the CSN PDF (Tablas 1-3), flattened with the
        // documented last-claiming-key rule; every value below is a cell of
        // the source matrix reached through chain expansion.
        let inert_rocas = ClearanceTable::es_conditional_inert(EsNormMaterial::Rocks);
        for (name, expected) in [
            ("U-238", 10.0),   // U-nat
            ("Pa-234m", 10.0), // U-nat
            ("U-234", 10.0),   // U-nat
            ("Th-230", 10.0),  // Th-230
            ("Ra-226", 10.0),  // Ra-226+
            ("Bi-214", 10.0),  // Ra-226+
            ("Pb-210", 10.0),  // Pb-210+
            ("Bi-210", 10.0),  // Pb-210+
            ("Po-210", 5.0),   // Po-210
            ("U-235", 10.0),   // U-235+
            ("Th-231", 10.0),  // U-235+
            ("Pa-231", 10.0),  // Pa-231
            ("Ac-227", 5.0),   // Ac-227+
            ("Tl-207", 5.0),   // Ac-227+
            ("Th-232", 5.0),   // Th-232
            ("Ra-228", 10.0),  // Ra-228+
            ("Ac-228", 10.0),  // Ra-228+
            ("Th-228", 5.0),   // Th-228+
            ("Tl-208", 5.0),   // Th-228+
            ("K-40", 10.0),    // K-40
        ] {
            let nuc = NuclideId::from_name(name).unwrap();
            assert_eq!(inert_rocas.get(nuc), Some(expected), "inert/rocas {name}");
        }

        // Column differentiation: the GAS/PETROLEO column of the same table.
        let inert_gas = ClearanceTable::es_conditional_inert(EsNormMaterial::OilGas);
        for (name, expected) in [
            ("U-238", 500.0),
            ("U-234", 500.0),
            ("Th-230", 500.0),
            ("Ra-226", 50.0),
            ("Pb-210", 100.0),
            ("Po-210", 50.0),
            ("U-235", 100.0),
            ("Pa-231", 500.0),
            ("Ac-227", 50.0),
            ("Th-232", 100.0),
            ("Ra-228", 50.0),
            ("Th-228", 10.0),
            ("K-40", 100.0),
        ] {
            let nuc = NuclideId::from_name(name).unwrap();
            assert_eq!(inert_gas.get(nuc), Some(expected), "inert/gas {name}");
        }

        // Landfill differentiation on shared material columns.
        let non_haz_rocas = ClearanceTable::es_conditional_non_hazardous(EsNormMaterial::Rocks);
        assert_eq!(
            non_haz_rocas.get(NuclideId::from_name("Po-210").unwrap()),
            Some(10.0)
        );
        assert_eq!(
            non_haz_rocas.get(NuclideId::from_name("Pb-210").unwrap()),
            Some(10.0)
        );
        let non_haz_gas = ClearanceTable::es_conditional_non_hazardous(EsNormMaterial::OilGas);
        assert_eq!(
            non_haz_gas.get(NuclideId::from_name("Pb-210").unwrap()),
            Some(500.0)
        );
        assert_eq!(
            non_haz_gas.get(NuclideId::from_name("Po-210").unwrap()),
            Some(100.0)
        );
        let haz_rocas = ClearanceTable::es_conditional_hazardous(EsNormMaterial::Rocks);
        assert_eq!(
            haz_rocas.get(NuclideId::from_name("Po-210").unwrap()),
            Some(100.0)
        );
        assert_eq!(
            haz_rocas.get(NuclideId::from_name("K-40").unwrap()),
            Some(50.0)
        );
        let haz_cenizas = ClearanceTable::es_conditional_hazardous(EsNormMaterial::Ashes);
        assert_eq!(
            haz_cenizas.get(NuclideId::from_name("Po-210").unwrap()),
            Some(50.0)
        );
        let haz_gas = ClearanceTable::es_conditional_hazardous(EsNormMaterial::OilGas);
        for (name, expected) in [("U-238", 500.0), ("Po-210", 500.0), ("K-40", 500.0)] {
            let nuc = NuclideId::from_name(name).unwrap();
            assert_eq!(haz_gas.get(nuc), Some(expected), "hazardous/gas {name}");
        }
    }

    #[test]
    fn es_chain_expansion_places_parent_value_on_every_member() {
        // Replays the committed TSVs row by row and asserts the flattened
        // table equals the documented expansion: each chain key's level on
        // each of its members, last claiming key winning.
        let chains = ClearanceTable::parse_es_norm_chains().unwrap();
        for (tsv, file) in [
            (ES_CSN_INERT_TSV, "es_csn_inert.tsv"),
            (ES_CSN_NON_HAZARDOUS_TSV, "es_csn_non_hazardous.tsv"),
            (ES_CSN_HAZARDOUS_TSV, "es_csn_hazardous.tsv"),
        ] {
            for material in EsNormMaterial::ALL {
                let table = ClearanceTable::try_es_table(tsv, file, material).unwrap();
                let mut expected: BTreeMap<NuclideId, f64> = BTreeMap::new();
                for line in tsv.lines().filter(|l| !l.trim().starts_with('#')) {
                    let trimmed = line.trim();
                    if trimmed.is_empty() {
                        continue;
                    }
                    let fields: Vec<&str> = trimmed.split('\t').collect();
                    let limit: f64 = fields[1 + material.column()].parse().unwrap();
                    for member in &chains[fields[0].trim()] {
                        expected.insert(*member, limit);
                    }
                }
                assert_eq!(
                    table.iter().collect::<Vec<_>>(),
                    expected.iter().map(|(n, v)| (*n, *v)).collect::<Vec<_>>(),
                    "{file} column {} mismatch",
                    material.column()
                );
                // Every member of every chain key carries its key's level
                // unless a later-listed key claims it (checked above).
                assert_eq!(expected.len(), 40, "{file}: union of Tabla 4 members");
            }
        }
    }

    #[test]
    fn es_hand_computed_screening_vectors_at_exact_equality() {
        // Inert landfill, rocks: 5/10 + 5/10 + 2.5/5 + 5/10 == 2 exactly.
        let table = ClearanceTable::es_conditional_inert(EsNormMaterial::Rocks);
        let inv = [
            (NuclideId::from_name("U-238").unwrap(), 5.0),
            (NuclideId::from_name("Th-230").unwrap(), 5.0),
            (NuclideId::from_name("Po-210").unwrap(), 2.5),
            (NuclideId::from_name("K-40").unwrap(), 5.0),
        ];
        assert_eq!(clearance_index(&inv, &table).unwrap(), 2.0);
        // Non-hazardous landfill, sands: 5/10 + 5/10 + 5/10 + 10/10 == 2.5.
        let table = ClearanceTable::es_conditional_non_hazardous(EsNormMaterial::Sands);
        let inv = [
            (NuclideId::from_name("Pb-210").unwrap(), 5.0),
            (NuclideId::from_name("Bi-210").unwrap(), 5.0),
            (NuclideId::from_name("K-40").unwrap(), 5.0),
            (NuclideId::from_name("Ra-228").unwrap(), 10.0),
        ];
        assert_eq!(clearance_index(&inv, &table).unwrap(), 2.5);
        // Hazardous landfill, oil & gas: 250/500 + 25/50 + 250/500 + 250/500
        // == 2 exactly.
        let table = ClearanceTable::es_conditional_hazardous(EsNormMaterial::OilGas);
        let inv = [
            (NuclideId::from_name("U-238").unwrap(), 250.0),
            (NuclideId::from_name("Ra-226").unwrap(), 25.0),
            (NuclideId::from_name("Po-210").unwrap(), 250.0),
            (NuclideId::from_name("K-40").unwrap(), 250.0),
        ];
        assert_eq!(clearance_index(&inv, &table).unwrap(), 2.0);
        // Sum-of-fractions shares the arithmetic: half the inert vector.
        let table = ClearanceTable::es_conditional_inert(EsNormMaterial::Rocks);
        let inv = [(NuclideId::from_name("Po-210").unwrap(), 2.5)];
        let out = sum_of_fractions(&inv, &table).unwrap();
        assert_eq!(out.sum, 0.5);
        assert_eq!(out.class, ClearanceClass::Satisfied);
        assert_eq!(
            out.max_nuclide,
            Some(NuclideId::from_name("Po-210").unwrap())
        );
    }
}
