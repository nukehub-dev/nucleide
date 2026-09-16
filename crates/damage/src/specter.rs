//! Opt-in vendored SPECTER fallback table: Table VII displacement-damage
//! cross sections plus the Table II `E_d` column.
//!
//! The fold kernels ([`crate::fold`]) take caller-supplied response slices
//! and never reach for this table implicitly. A caller that has no
//! NJOY/SPECTER-class pipeline of its own can instead build a one-group
//! response slice from the vendored values explicitly:
//!
//! ```text
//! response = [SpecterTable::for_spectrum(SpecterSpectrum::Hfir)
//!                 .dpa_xs_barns("Fe").unwrap()]
//! dpa = nrt_dpa(&[fluence], &response, &[0.0, 20.0], 1.0)
//! ```
//!
//! ## What is vendored (and what is not)
//!
//! The embedded `data/specter_table_vii.tsv` transcribes Table VII of
//! Greenwood & Smither, ANL/FPP/TM-197 (January 1985 — a US-government
//! work, public domain): spectrum-averaged *damage-energy* cross sections
//! in keV-barns for 24 elements (Be through Pb) across 7 spectra, plus the
//! Table II Lindhard-cutoff threshold `E_d` (eV) per element. Displacement
//! cross sections in barns follow the report's own rule,
//! `barns = keV-b × 0.8 / (2 × E_d[keV])`, with `E_d` from the vendored
//! column or a value of the caller's choice ([`SpecterTable::dpa_xs_barns`]
//! vs [`SpecterTable::dpa_xs_barns_with_ed`]).
//!
//! Scope decisions, recorded: **displacement XS only** (no gas-production
//! columns — the smaller honest step); **Table VII, not Appendix A** (the
//! group-wise Appendix A printouts are image-only pages in the report scan
//! and cannot meet the transcription bar, so they stay out); **24 elements,
//! not the 41 Table I entries** (Table VII is the only cleanly transcribable
//! summary table — the same narrowing precedent as documented project-wide).
//!
//! ## Unit and spectrum basis
//!
//! Stored values are damage-energy cross sections (keV-b), *not*
//! displacement cross sections: the `0.8/2E_d` conversion is applied on
//! access. All seven columns are whole-spectrum averages (the report's
//! footnote: "averages over the entire spectrum"), so a vendored value is a
//! one-group response for its spectrum — bounds document the fold but never
//! enter the sum, per the [`crate::fold`] convention. Column semantics:
//!
//! - `thermal`: thermal (n,gamma) damage only (`T_GAM × σ_0` at 2200 m/s),
//!   not a full-spectrum average — do not fold it with a thermal flux as if
//!   it were one.
//! - `fission`: 235U fission spectrum average.
//! - `14mev`: average of the SPECTER energy group from 14–15 MeV.
//! - `hfir` / `ebr2` / `fftf` / `fusion`: HFIR (ORNL, PTP position), EBR-II
//!   (ANL-W, Row 2), FFTF (HEDL, MOTA), and the UWMAK first-wall fusion
//!   spectrum.
//!
//! The HFIR column reproduces the validation oracle's Table VI spots
//! (Fe/Ti/Cu spectrum-averaged dpa cross sections) at Table VII print
//! precision (3 significant figures): the table and the oracle read the
//! same report pages. Element keys are the report's symbols (`"Fe"`,
//! case-sensitive); `"Ag"` is natural silver and `"W"` natural tungsten
//! (their Table II isotopes share one `E_d` each — see the TSV header).

use std::collections::BTreeMap;
use std::sync::OnceLock;

use crate::error::{Error, Result};

/// Embedded Table VII transcription (see the module docs for provenance).
const SPECTER_TABLE_VII_TSV: &str = include_str!("data/specter_table_vii.tsv");

/// Number of spectrum columns in the TSV (thermal .. fusion).
const N_SPECTRA: usize = 7;

/// Vendored SPECTER Table VII spectra, in TSV column order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SpecterSpectrum {
    /// Thermal (n,gamma) damage only (`T_GAM × σ_0`); not a full-spectrum
    /// average.
    Thermal,
    /// 235U fission spectrum average.
    Fission,
    /// Average of the SPECTER energy group from 14–15 MeV.
    FourteenMeV,
    /// HFIR (ORNL, PTP position) whole-spectrum average.
    Hfir,
    /// EBR-II (ANL-W, Row 2) whole-spectrum average.
    EbrII,
    /// FFTF (HEDL, MOTA) whole-spectrum average.
    Fftf,
    /// UWMAK first-wall fusion spectrum whole-spectrum average.
    Fusion,
}

impl SpecterSpectrum {
    /// All seven spectra in TSV column order.
    pub const ALL: [Self; N_SPECTRA] = [
        Self::Thermal,
        Self::Fission,
        Self::FourteenMeV,
        Self::Hfir,
        Self::EbrII,
        Self::Fftf,
        Self::Fusion,
    ];

    /// Column index of this spectrum in the embedded TSV.
    fn column(self) -> usize {
        match self {
            Self::Thermal => 0,
            Self::Fission => 1,
            Self::FourteenMeV => 2,
            Self::Hfir => 3,
            Self::EbrII => 4,
            Self::Fftf => 5,
            Self::Fusion => 6,
        }
    }

    /// Canonical lowercase name of this spectrum (the TSV column key).
    pub fn name(self) -> &'static str {
        match self {
            Self::Thermal => "thermal",
            Self::Fission => "fission",
            Self::FourteenMeV => "14mev",
            Self::Hfir => "hfir",
            Self::EbrII => "ebr2",
            Self::Fftf => "fftf",
            Self::Fusion => "fusion",
        }
    }

    /// Parse a spectrum name (case-insensitive; separators ignored).
    ///
    /// Accepts the canonical [`Self::name`] spellings plus `ebr-ii` and
    /// `14_mev`/`14 mev` for the digit-led columns. Anything else is a loud
    /// [`Error::UnknownSpecterSpectrum`].
    pub fn parse(name: &str) -> Result<Self> {
        let key: String = name
            .chars()
            .filter(|c| !c.is_whitespace() && *c != '-' && *c != '_')
            .collect::<String>()
            .to_ascii_lowercase();
        for spectrum in Self::ALL {
            if key == spectrum.name() {
                return Ok(spectrum);
            }
        }
        // `ebr2` canonical; the report spells the reactor EBR-II.
        if key == "ebrii" {
            return Ok(Self::EbrII);
        }
        Err(Error::UnknownSpecterSpectrum {
            spectrum: name.to_string(),
        })
    }
}

impl std::fmt::Display for SpecterSpectrum {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.name())
    }
}

/// One element's vendored row for a single spectrum.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpecterEntry {
    /// Damage-energy cross section in keV-b, verbatim from Table VII.
    pub damage_energy_kev_b: f64,
    /// Lindhard-cutoff threshold displacement energy in eV, from Table II.
    pub ed_ev: f64,
}

impl SpecterEntry {
    /// Displacement cross section in barns via the report's own rule with
    /// the vendored `E_d`: `keV-b × 0.8 / (2 × E_d[keV])`.
    pub fn dpa_xs_barns(&self) -> f64 {
        self.damage_energy_kev_b * 0.8 / (2.0 * self.ed_ev / 1000.0)
    }

    /// Displacement cross section in barns with a caller-chosen `E_d`
    /// (eV) — the report's "or a value of your choice".
    ///
    /// `ed_ev` must be finite and strictly positive (loud
    /// [`Error::NonFinite`]/[`Error::NonPositive`], reusing the
    /// crate-wide quantity checks).
    pub fn dpa_xs_barns_with_ed(&self, ed_ev: f64) -> Result<f64> {
        if !ed_ev.is_finite() {
            return Err(Error::NonFinite("ed_ev"));
        }
        if ed_ev <= 0.0 {
            return Err(Error::NonPositive("ed_ev"));
        }
        Ok(self.damage_energy_kev_b * 0.8 / (2.0 * ed_ev / 1000.0))
    }
}

/// One spectrum's vendored SPECTER fallback table: element symbol
/// (`"Fe"`, case-sensitive) to [`SpecterEntry`].
///
/// Build one with [`SpecterTable::for_spectrum`] (cached) or
/// [`SpecterTable::try_for_spectrum`] (fallible) — the same `try_` + cache
/// split as the vendored clearance tables. The table is an opt-in caller
/// input to the [`crate::fold`] kernels, never consulted implicitly.
#[derive(Debug, Clone, PartialEq)]
pub struct SpecterTable {
    spectrum: SpecterSpectrum,
    entries: BTreeMap<String, SpecterEntry>,
}

impl SpecterTable {
    /// Vendored fallback table for `spectrum`, parsed lazily on first call.
    ///
    /// The transcription is committed data pinned by row-count tests, so a
    /// corrupt transcription is treated as unreachable-in-practice and
    /// panics with the offending line; fallible callers use
    /// [`Self::try_for_spectrum`].
    pub fn for_spectrum(spectrum: SpecterSpectrum) -> Self {
        static TABLES: OnceLock<[SpecterTable; N_SPECTRA]> = OnceLock::new();
        TABLES.get_or_init(|| {
            std::array::from_fn(|column| {
                let spectrum = SpecterSpectrum::ALL[column];
                Self::try_for_spectrum(spectrum)
                    .unwrap_or_else(|e| panic!("specter_table_vii.tsv: {e}"))
            })
        })[spectrum.column()]
        .clone()
    }

    /// Fallible parse of the embedded Table VII transcription for
    /// `spectrum`.
    ///
    /// Same data as [`Self::for_spectrum`] without the lazy cache: every
    /// malformed line (wrong field count, bad value, bad `E_d`, duplicate
    /// element) is a loud [`Error::Parse`] carrying the 1-based TSV line
    /// number.
    pub fn try_for_spectrum(spectrum: SpecterSpectrum) -> Result<Self> {
        let library = parse_library_tsv(SPECTER_TABLE_VII_TSV)?;
        let column = spectrum.column();
        let entries = library
            .into_iter()
            .map(|(element, row)| {
                (
                    element,
                    SpecterEntry {
                        damage_energy_kev_b: row.values[column],
                        ed_ev: row.ed_ev,
                    },
                )
            })
            .collect();
        Ok(Self { spectrum, entries })
    }

    /// The spectrum this table projects.
    pub fn spectrum(&self) -> SpecterSpectrum {
        self.spectrum
    }

    /// Number of element entries (24 for the committed transcription).
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// True when the table has no entries.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// True when the table carries `element`.
    pub fn contains(&self, element: &str) -> bool {
        self.entries.contains_key(element)
    }

    /// The vendored entry for `element`, or `None` when outside the
    /// 24-element Table VII set.
    pub fn get(&self, element: &str) -> Option<SpecterEntry> {
        self.entries.get(element).copied()
    }

    /// Verbatim Table VII damage-energy cross section (keV-b) for
    /// `element` at this table's spectrum, or `None` for unknown elements.
    pub fn damage_energy_kev_b(&self, element: &str) -> Option<f64> {
        self.get(element).map(|entry| entry.damage_energy_kev_b)
    }

    /// Vendored Table II `E_d` (eV) for `element`, or `None` for unknown
    /// elements.
    pub fn ed_ev(&self, element: &str) -> Option<f64> {
        self.get(element).map(|entry| entry.ed_ev)
    }

    /// Displacement cross section in barns for `element` at this table's
    /// spectrum via the vendored `E_d`, or `None` for unknown elements.
    pub fn dpa_xs_barns(&self, element: &str) -> Option<f64> {
        self.get(element).map(|entry| entry.dpa_xs_barns())
    }

    /// Displacement cross section in barns for `element` with a
    /// caller-chosen `E_d` (eV).
    ///
    /// Unknown elements are a loud [`Error::UnknownSpecterElement`]; a bad
    /// `E_d` is [`Error::NonFinite`]/[`Error::NonPositive`] (see
    /// [`SpecterEntry::dpa_xs_barns_with_ed`]).
    pub fn dpa_xs_barns_with_ed(&self, element: &str, ed_ev: f64) -> Result<f64> {
        let entry = self
            .get(element)
            .ok_or_else(|| Error::UnknownSpecterElement {
                element: element.to_string(),
            })?;
        entry.dpa_xs_barns_with_ed(ed_ev)
    }

    /// Element symbols in lexicographic order (BTreeMap canonical order,
    /// like the vendored clearance tables — not report order).
    pub fn elements(&self) -> Vec<String> {
        self.entries.keys().cloned().collect()
    }

    /// Iterate entries in lexicographic element-symbol order.
    pub fn iter(&self) -> impl Iterator<Item = (&String, &SpecterEntry)> {
        self.entries.iter()
    }
}

/// Vendored Table II `E_d` (eV) for `element`.
///
/// Convenience over [`SpecterTable::ed_ev`] with a loud
/// [`Error::UnknownSpecterElement`] instead of `None`; the value is
/// identical in every spectrum view.
pub fn specter_ed_ev(element: &str) -> Result<f64> {
    SpecterTable::for_spectrum(SpecterSpectrum::Hfir)
        .ed_ev(element)
        .ok_or_else(|| Error::UnknownSpecterElement {
            element: element.to_string(),
        })
}

/// One parsed TSV row: all seven spectra plus `E_d`.
#[derive(Debug, Clone, Copy, PartialEq)]
struct FullRow {
    /// Damage-energy cross sections in keV-b, TSV column order.
    values: [f64; N_SPECTRA],
    /// Table II `E_d` in eV.
    ed_ev: f64,
}

/// Parse Table VII TSV text into element rows (1-based line numbers in
/// errors). The embedded transcription goes through here via
/// [`SpecterTable::try_for_spectrum`]; tests drive it with synthetic text.
fn parse_library_tsv(text: &str) -> Result<BTreeMap<String, FullRow>> {
    const FIELDS: usize = N_SPECTRA + 2;
    let mut library = BTreeMap::new();
    for (index, line) in text.lines().enumerate() {
        let line_no = index + 1;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() != FIELDS {
            return Err(Error::Parse {
                line: line_no,
                msg: format!(
                    "specter_table_vii.tsv: expected `element` plus {N_SPECTRA} \
                     spectrum columns plus `ed_ev` ({FIELDS} fields), found {}",
                    fields.len()
                ),
            });
        }
        let element = fields[0].trim().to_string();
        if element.is_empty() {
            return Err(Error::Parse {
                line: line_no,
                msg: "specter_table_vii.tsv: empty element symbol".to_string(),
            });
        }
        let mut values = [0.0f64; N_SPECTRA];
        for (column, text) in fields[1..=N_SPECTRA].iter().enumerate() {
            let value: f64 = text.trim().parse().map_err(|_| Error::Parse {
                line: line_no,
                msg: format!("specter_table_vii.tsv: bad keV-b value `{text}`"),
            })?;
            if !value.is_finite() || value <= 0.0 {
                return Err(Error::Parse {
                    line: line_no,
                    msg: format!(
                        "specter_table_vii.tsv: keV-b value must be finite and > 0, got `{text}`"
                    ),
                });
            }
            values[column] = value;
        }
        let ed_ev: f64 = fields[N_SPECTRA + 1]
            .trim()
            .parse()
            .map_err(|_| Error::Parse {
                line: line_no,
                msg: format!(
                    "specter_table_vii.tsv: bad ed_ev `{}`",
                    fields[N_SPECTRA + 1]
                ),
            })?;
        if !ed_ev.is_finite() || ed_ev <= 0.0 {
            return Err(Error::Parse {
                line: line_no,
                msg: format!(
                    "specter_table_vii.tsv: ed_ev must be finite and > 0, got `{}`",
                    fields[N_SPECTRA + 1]
                ),
            });
        }
        if library
            .insert(element.clone(), FullRow { values, ed_ev })
            .is_some()
        {
            return Err(Error::Parse {
                line: line_no,
                msg: format!("specter_table_vii.tsv: duplicate element `{element}`"),
            });
        }
    }
    Ok(library)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fold::nrt_dpa;

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
    fn table_vii_loads_with_pinned_row_count() {
        // Pins the committed transcription: 24 element rows (Be through Pb,
        // report order). Any edit to the TSV must update this count
        // deliberately, which keeps `for_spectrum` honest about its
        // unreachable-in-practice panic.
        assert_eq!(tsv_data_rows(SPECTER_TABLE_VII_TSV), 24);
        for spectrum in SpecterSpectrum::ALL {
            let table = SpecterTable::try_for_spectrum(spectrum).unwrap();
            assert_eq!(table.len(), 24);
            assert_eq!(SpecterTable::for_spectrum(spectrum), table);
            assert_eq!(table.spectrum(), spectrum);
            assert!(!table.is_empty());
        }
        // Lexicographic (BTreeMap canonical) order, like the clearance tables.
        let table = SpecterTable::for_spectrum(SpecterSpectrum::Hfir);
        let elements = table.elements();
        assert_eq!(elements.first().map(String::as_str), Some("Ag"));
        assert_eq!(elements.last().map(String::as_str), Some("Zr"));
        assert_eq!(elements.len(), 24);
        assert!(table.contains("Fe"));
        assert!(!table.contains("U"));
        assert!(!table.contains("U"));
    }

    #[test]
    fn transcription_spots_across_all_columns() {
        // Hand-read from the report scan (Table VII, report p. -31-):
        // (element, spectrum, keV-b). Covers every column and both ends of
        // the element range, including the trailing-point integers.
        let table = |s| SpecterTable::for_spectrum(s);
        for (element, spectrum, expected) in [
            ("Be", SpecterSpectrum::Thermal, 0.010),
            ("C", SpecterSpectrum::Thermal, 0.002),
            ("Co", SpecterSpectrum::Thermal, 13.38),
            ("Ag", SpecterSpectrum::Thermal, 10.06),
            ("Pb", SpecterSpectrum::Thermal, 0.021),
            ("Nb", SpecterSpectrum::Thermal, 0.128),
            ("Au", SpecterSpectrum::Fission, 50.2),
            ("V", SpecterSpectrum::Fission, 101.0),
            ("Ni", SpecterSpectrum::FourteenMeV, 300.0),
            ("Cu", SpecterSpectrum::FourteenMeV, 296.0),
            ("Na", SpecterSpectrum::FourteenMeV, 140.0),
            ("Fe", SpecterSpectrum::Hfir, 19.1),
            ("Ti", SpecterSpectrum::Hfir, 21.9),
            ("Cu", SpecterSpectrum::Hfir, 18.6),
            ("Ag", SpecterSpectrum::Hfir, 23.0),
            ("Mo", SpecterSpectrum::EbrII, 53.6),
            ("Zr", SpecterSpectrum::EbrII, 54.3),
            ("Ta", SpecterSpectrum::Fftf, 16.9),
            ("W", SpecterSpectrum::Fftf, 16.5),
            ("Cr", SpecterSpectrum::Fusion, 104.0),
            ("Ni", SpecterSpectrum::Fusion, 109.4),
            ("W", SpecterSpectrum::Fusion, 68.0),
            ("Pb", SpecterSpectrum::Fusion, 71.4),
        ] {
            assert_eq!(
                table(spectrum).damage_energy_kev_b(element),
                Some(expected),
                "{element} {spectrum}"
            );
        }
        // Table II E_d spots (report p. -14-): Be 31, Al 27, Mo 60,
        // Ag 60 (both isotopes), Ta 53, W 90 (all four isotopes),
        // Au 30, Pb 25, Fe 40.
        let hfir = SpecterTable::for_spectrum(SpecterSpectrum::Hfir);
        for (element, expected) in [
            ("Be", 31.0),
            ("C", 31.0),
            ("Na", 25.0),
            ("Al", 27.0),
            ("Fe", 40.0),
            ("Mo", 60.0),
            ("Ag", 60.0),
            ("Ta", 53.0),
            ("W", 90.0),
            ("Au", 30.0),
            ("Pb", 25.0),
        ] {
            assert_eq!(hfir.ed_ev(element), Some(expected), "{element} E_d");
        }
        // Unknown elements read None, never a default.
        assert_eq!(hfir.damage_energy_kev_b("U"), None);
        assert_eq!(hfir.ed_ev("H"), None);
        assert_eq!(hfir.dpa_xs_barns("Fe235"), None);
    }

    #[test]
    fn hfir_column_reproduces_oracle_spots_at_table_print_precision() {
        // The table must reproduce its own oracle: Table VII HFIR keV-b
        // through the report's 0.8/2E_d rule against the validation
        // harness's Table VI spectrum-averaged dpa cross sections
        // (Fe 191.18, Ti 218.73, Cu 186.49 barns). Table VII prints 3
        // significant figures (half-ulp up to ~2.7e-3 on these cells), so
        // the gate sits at 5e-3 relative — the finest agreement the
        // transcribed inputs allow.
        let table = SpecterTable::for_spectrum(SpecterSpectrum::Hfir);
        for (element, printed_xs) in [("Fe", 191.18), ("Ti", 218.73), ("Cu", 186.49)] {
            let got = table.dpa_xs_barns(element).unwrap();
            let err = (got - printed_xs).abs() / printed_xs;
            assert!(err < 5e-3, "{element}: got {got}, want {printed_xs}");
        }
        // Same check through the one-group fluence fold the oracle runs:
        // HFIR-CTR32 fluence 4.78373E+22 n/cm2, seconds = 1.
        let fluence = vec![4.78373e22];
        let bounds = vec![0.0, 20.0];
        for (element, printed_dpa) in [("Fe", 9.1455), ("Ti", 10.464), ("Cu", 8.9212)] {
            let xs = table.dpa_xs_barns(element).unwrap();
            let got = nrt_dpa(&fluence, &[xs], &bounds, 1.0).unwrap();
            let err = (got - printed_dpa).abs() / printed_dpa;
            assert!(err < 5e-3, "{element}: got {got}, want {printed_dpa}");
        }
    }

    #[test]
    fn vendored_fold_equals_caller_slice_fold_exactly() {
        // The table is just another caller input: folding a vendored value
        // and folding the same floats supplied by the caller agree bit for
        // bit, and both equal the hand product.
        let table = SpecterTable::for_spectrum(SpecterSpectrum::Fusion);
        let xs = table.dpa_xs_barns("Fe").unwrap();
        assert_eq!(xs, 101.8 * 0.8 / (2.0 * 40.0 / 1000.0));
        let bounds = vec![0.0, 20.0];
        // One-group fusion fold of the vendored Fe value.
        let vendored = nrt_dpa(&[3.0e14], &[xs], &bounds, 1.0).unwrap();
        let caller = nrt_dpa(&[3.0e14], &[xs], &bounds, 1.0).unwrap();
        assert_eq!(vendored, caller);
        // Hand product with the fold's own association
        // (scale * seconds * acc): exact to the ulp.
        assert_eq!(vendored, 1.0e-24 * 1.0 * (3.0e14 * xs));
        // Caller-chosen E_d ("a value of your choice"): E_d = 30 eV scales
        // the 40 eV default by exactly 4/3.
        let custom = table.dpa_xs_barns_with_ed("Fe", 30.0).unwrap();
        assert_eq!(custom, xs * 40.0 / 30.0);
        // Multi-element response vector folds to the hand sum.
        let elements = ["Ti", "Fe", "Cu"];
        let response: Vec<f64> = elements
            .iter()
            .map(|el| table.dpa_xs_barns(el).unwrap())
            .collect();
        let flux3 = vec![1.0e14, 2.0e14, 4.0e14];
        let hand: f64 = flux3.iter().zip(response.iter()).map(|(f, r)| f * r).sum();
        assert_eq!(
            nrt_dpa(&flux3, &response, &[0.0, 1.0, 10.0, 20.0], 2.0).unwrap(),
            1.0e-24 * 2.0 * hand
        );
    }

    #[test]
    fn spectrum_names_parse_case_insensitively() {
        for spectrum in SpecterSpectrum::ALL {
            assert_eq!(SpecterSpectrum::parse(spectrum.name()), Ok(spectrum));
            assert_eq!(
                SpecterSpectrum::parse(&spectrum.name().to_ascii_uppercase()),
                Ok(spectrum)
            );
        }
        assert_eq!(SpecterSpectrum::parse("ebr-ii"), Ok(SpecterSpectrum::EbrII));
        assert_eq!(
            SpecterSpectrum::parse("14_mev"),
            Ok(SpecterSpectrum::FourteenMeV)
        );
        assert_eq!(SpecterSpectrum::parse("HFIR"), Ok(SpecterSpectrum::Hfir));
        match SpecterSpectrum::parse("pwr") {
            Err(Error::UnknownSpecterSpectrum { spectrum }) => assert_eq!(spectrum, "pwr"),
            other => panic!("expected UnknownSpecterSpectrum, got {other:?}"),
        }
        assert_eq!(SpecterSpectrum::Hfir.to_string(), "hfir");
    }

    #[test]
    fn malformed_rows_and_bad_keys_are_loud_errors() {
        // Wrong field count.
        match parse_library_tsv("Fe\t1.0\t2.0\n") {
            Err(Error::Parse { line, msg }) => {
                assert_eq!(line, 1);
                assert!(msg.contains("9 fields"), "msg was `{msg}`");
            }
            other => panic!("expected Parse, got {other:?}"),
        }
        // Bad keV-b value.
        let bad_value = "Fe\t1.01\t84.4\t290.0\t19.1\t46.3\t27.3\tnan\t40\n";
        assert!(matches!(
            parse_library_tsv(bad_value),
            Err(Error::Parse { .. })
        ));
        // Non-positive keV-b value.
        let zero_value = "Fe\t1.01\t84.4\t290.0\t19.1\t46.3\t27.3\t0.0\t40\n";
        assert!(matches!(
            parse_library_tsv(zero_value),
            Err(Error::Parse { .. })
        ));
        // Bad E_d.
        let bad_ed = "Fe\t1.01\t84.4\t290.0\t19.1\t46.3\t27.3\t101.8\t0\n";
        assert!(matches!(
            parse_library_tsv(bad_ed),
            Err(Error::Parse { .. })
        ));
        // Duplicate element.
        let row = "Fe\t1.01\t84.4\t290.0\t19.1\t46.3\t27.3\t101.8\t40\n";
        let doubled = format!("{row}{row}");
        match parse_library_tsv(&doubled) {
            Err(Error::Parse { line, msg }) => {
                assert_eq!(line, 2);
                assert!(msg.contains("duplicate"), "msg was `{msg}`");
            }
            other => panic!("expected Parse, got {other:?}"),
        }
        // Unknown element with a caller E_d.
        let table = SpecterTable::for_spectrum(SpecterSpectrum::Hfir);
        match table.dpa_xs_barns_with_ed("U", 40.0) {
            Err(Error::UnknownSpecterElement { element }) => assert_eq!(element, "U"),
            other => panic!("expected UnknownSpecterElement, got {other:?}"),
        }
        // Bad caller E_d values reuse the crate-wide quantity checks.
        assert!(matches!(
            table.dpa_xs_barns_with_ed("Fe", 0.0),
            Err(Error::NonPositive("ed_ev"))
        ));
        assert!(matches!(
            table.dpa_xs_barns_with_ed("Fe", f64::NAN),
            Err(Error::NonFinite("ed_ev"))
        ));
    }
}
