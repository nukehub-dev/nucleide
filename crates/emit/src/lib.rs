//! Single-material emission to legacy transport-code cards.
//!
//! [`emit_all`] renders one [`Material`] through five code dialects — MCNP,
//! Serpent, FLUKA, ALARA, PARTISN — and [`drift_table`] reports how much mass
//! survives each translation. Every emitter is pure glue over the workspace
//! `*-io` crates plus [`nucleide_nuclei::dialects`]; no new physics lives here.
//!
//! # Drift semantics
//!
//! Each emitter returns the mass it could represent per nuclide
//! ([`Emitted::accounted`]) plus what it had to drop ([`Emitted::dropped`],
//! e.g. a nuclide with no FLUKA name). [`DriftRow::rel_drift`] is
//! `(mass_in - mass_out) / mass_in`. [`Emitted::reparsed`] tells whether the
//! emitted text was machine-verified by feeding it back through that code's
//! reader (MCNP and ALARA only — Serpent, FLUKA, and PARTISN have no material
//! readers in this workspace, so their drift is analytic).
//!
//! # Example
//!
//! ```rust
//! use nucleide_emit::{EmitOptions, emit_drift};
//! use nucleide_material::Material;
//! use nucleide_nuclei::NuclideId;
//!
//! let mut mat = Material::new();
//! mat.add_nuclide(NuclideId::from_name("U235").unwrap(), 5.0);
//! mat.add_nuclide(NuclideId::from_name("U238").unwrap(), 95.0);
//! let opts = EmitOptions::new("leu").with_density(10.0);
//! let (emitted, table) = emit_drift(&mat, &opts).unwrap();
//! assert_eq!(emitted.len(), 5);
//! assert!(table.worst_rel_drift() < 1e-9);
//! ```

pub mod alara;
pub mod fluka;
pub mod mcnp;
pub mod partisn;
pub mod serpent;

use nucleide_material::Material;
use nucleide_nuclei::NuclideId;

/// One of the five supported transport-code dialects.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Code {
    Mcnp,
    Serpent,
    Fluka,
    Alara,
    Partisn,
}

impl Code {
    /// All codes in emission order.
    pub fn all() -> [Code; 5] {
        [
            Code::Mcnp,
            Code::Serpent,
            Code::Fluka,
            Code::Alara,
            Code::Partisn,
        ]
    }
}

impl std::fmt::Display for Code {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Code::Mcnp => write!(f, "MCNP"),
            Code::Serpent => write!(f, "Serpent"),
            Code::Fluka => write!(f, "FLUKA"),
            Code::Alara => write!(f, "ALARA"),
            Code::Partisn => write!(f, "PARTISN"),
        }
    }
}

/// Shared emission settings.
#[derive(Debug, Clone)]
pub struct EmitOptions {
    /// Material/card name used by every dialect.
    pub name: String,
    /// MCNP `m<number>` card number.
    pub mcnp_number: u32,
    /// MCNP cross-section library suffix (`80c` renders `92235.80c`).
    pub xs_suffix: String,
    /// Mass density [g/cm³] for dialects that need one (Serpent, FLUKA,
    /// PARTISN). Falls back to [`Material::density`]; errors when absent.
    pub density: Option<f64>,
    /// FLUKA material index number.
    pub fluka_fid: u32,
    /// PARTISN zone id for the single emitted zone.
    pub partisn_zone: u32,
}

impl EmitOptions {
    /// Options with conventional defaults (`m1`, `80c`, FLUKA fid 1, zone 1).
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            mcnp_number: 1,
            xs_suffix: "80c".to_string(),
            density: None,
            fluka_fid: 1,
            partisn_zone: 1,
        }
    }

    /// Override the mass density [g/cm³].
    pub fn with_density(mut self, density: f64) -> Self {
        self.density = Some(density);
        self
    }

    fn density_for(&self, mat: &Material) -> Result<f64> {
        self.density
            .or_else(|| mat.density())
            .ok_or(Error::MissingDensity)
    }
}

/// One nuclide the emitter could not represent.
#[derive(Debug, Clone, PartialEq)]
pub struct Dropped {
    /// Nuclide that was skipped.
    pub id: NuclideId,
    /// Mass [g] left out of the emitted cards.
    pub mass: f64,
    /// Machine-readable reason (e.g. `"no-fluka-name"`, `"no-atomic-mass"`).
    pub reason: String,
}

/// Cards emitted for one dialect.
#[derive(Debug, Clone)]
pub struct Emitted {
    /// Dialect these cards belong to.
    pub code: Code,
    /// Card text, ready to paste into a deck.
    pub text: String,
    /// Mass [g] accounted per emitted nuclide, in [`Material`] order.
    pub accounted: Vec<(NuclideId, f64)>,
    /// Nuclides skipped with reasons.
    pub dropped: Vec<Dropped>,
    /// Whether `text` was verified by re-parsing with the code's own reader.
    /// Only MCNP and ALARA have material readers in this workspace.
    pub reparsed: bool,
}

impl Emitted {
    /// Mass [g] represented by these cards.
    pub fn mass_out(&self) -> f64 {
        self.accounted.iter().map(|(_, m)| m).sum()
    }
}

/// One row of the mass-drift report.
#[derive(Debug, Clone)]
pub struct DriftRow {
    /// Dialect this row covers.
    pub code: Code,
    /// Input mass [g].
    pub mass_in: f64,
    /// Mass represented in the emitted cards [g].
    pub mass_out: f64,
    /// `(mass_in - mass_out) / mass_in`; zero for lossless emission.
    pub rel_drift: f64,
    /// Nuclides skipped with reasons.
    pub dropped: Vec<Dropped>,
    /// Whether the emitted text was re-parse verified.
    pub reparsed: bool,
}

/// Mass conservation across all five dialects.
#[derive(Debug, Clone)]
pub struct DriftTable {
    /// Material name from [`EmitOptions::name`].
    pub name: String,
    /// One row per dialect, in [`Code::all`] order.
    pub rows: Vec<DriftRow>,
}

impl DriftTable {
    /// Largest relative drift over all dialects.
    pub fn worst_rel_drift(&self) -> f64 {
        self.rows
            .iter()
            .map(|r| r.rel_drift.abs())
            .fold(0.0, f64::max)
    }
}

/// Errors raised while emitting cards.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Material is empty or its masses sum to a non-positive value.
    #[error("material is empty or its masses sum to a non-positive value")]
    Degenerate,
    /// A dialect needs a mass density but neither options nor material has one.
    #[error("emission requires a mass density but none was set")]
    MissingDensity,
    /// The emitted text failed to re-parse with the code's own reader.
    #[error("emitted {code} text failed to re-parse: {detail}")]
    Reparse {
        /// Dialect whose text did not survive its reader.
        code: Code,
        /// Reader error.
        detail: String,
    },
    /// Material-layer failure (mass tables, fractions).
    #[error(transparent)]
    Material(#[from] nucleide_material::Error),
    /// Nuclide-dialect failure.
    #[error(transparent)]
    Nuclei(#[from] nucleide_nuclei::Error),
    /// MCNP reader failure during re-parse verification.
    #[error(transparent)]
    Mcnp(#[from] nucleide_mcnp_io::inp::Error),
    /// ALARA reader failure during re-parse verification.
    #[error(transparent)]
    Alara(#[from] nucleide_alara_io::Error),
    /// FLUKA naming failure.
    #[error(transparent)]
    Fluka(#[from] nucleide_fluka_io::material::Error),
}

/// Crate-local result alias.
pub type Result<T> = std::result::Result<T, Error>;

/// Emit one [`Material`] through all five dialects, in [`Code::all`] order.
pub fn emit_all(mat: &Material, opts: &EmitOptions) -> Result<Vec<Emitted>> {
    if mat.mass() <= 0.0 || !mat.mass().is_finite() {
        return Err(Error::Degenerate);
    }
    Ok(vec![
        mcnp::emit_mcnp(mat, opts)?,
        serpent::emit_serpent(mat, opts)?,
        fluka::emit_fluka(mat, opts)?,
        alara::emit_alara(mat, opts)?,
        partisn::emit_partisn(mat, opts)?,
    ])
}

/// Build the mass-drift report for already-emitted cards.
pub fn drift_table(mat: &Material, emitted: &[Emitted]) -> DriftTable {
    let mass_in = mat.mass();
    let rows = emitted
        .iter()
        .map(|e| {
            let mass_out = e.mass_out();
            DriftRow {
                code: e.code,
                mass_in,
                mass_out,
                rel_drift: if mass_in == 0.0 {
                    0.0
                } else {
                    (mass_in - mass_out) / mass_in
                },
                dropped: e.dropped.clone(),
                reparsed: e.reparsed,
            }
        })
        .collect();
    DriftTable {
        name: String::new(),
        rows,
    }
}

/// Emit through all dialects and report mass drift in one call.
pub fn emit_drift(mat: &Material, opts: &EmitOptions) -> Result<(Vec<Emitted>, DriftTable)> {
    let emitted = emit_all(mat, opts)?;
    let mut table = drift_table(mat, &emitted);
    table.name = opts.name.clone();
    Ok((emitted, table))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn metal() -> Material {
        // Uranium metal: every dialect represents both nuclides, so the
        // full table is lossless. (Light isotopes such as H1/O16 have no
        // FLUKA isotope-table entry and exercise the dropped path instead.)
        let mut mat = Material::new();
        mat.add_nuclide(NuclideId::from_name("U235").unwrap(), 5.0);
        mat.add_nuclide(NuclideId::from_name("U238").unwrap(), 95.0);
        mat
    }

    #[test]
    fn emit_all_covers_five_codes_lossless() {
        let opts = EmitOptions::new("umetal").with_density(19.1);
        let (emitted, table) = emit_drift(&metal(), &opts).unwrap();
        assert_eq!(emitted.len(), 5);
        assert_eq!(
            emitted.iter().map(|e| e.code).collect::<Vec<_>>(),
            Code::all()
        );
        assert_eq!(table.name, "umetal");
        assert_eq!(table.rows.len(), 5);
        for row in &table.rows {
            assert!((row.mass_in - 100.0).abs() < 1e-12, "{:?}", row.code);
            assert!((row.mass_out - 100.0).abs() < 1e-9, "{:?}", row.code);
            assert!(row.rel_drift.abs() < 1e-9, "{:?}", row.code);
            assert!(row.dropped.is_empty(), "{:?}", row.code);
        }
        assert!(table.worst_rel_drift() < 1e-9);
        let reparsed: Vec<Code> = emitted
            .iter()
            .filter(|e| e.reparsed)
            .map(|e| e.code)
            .collect();
        assert_eq!(reparsed, vec![Code::Mcnp, Code::Alara]);
    }

    #[test]
    fn drift_reports_fluka_loss() {
        // O16 has no isotope entry in the vendored FLUKA table (H1 maps to
        // HYDROG-1), so 80 of 100 g drift away on that row while the other
        // four stay lossless.
        let mut mat = Material::new();
        mat.add_nuclide(NuclideId::from_name("H1").unwrap(), 20.0);
        mat.add_nuclide(NuclideId::from_name("O16").unwrap(), 80.0);
        let opts = EmitOptions::new("water").with_density(1.0);
        let (_, table) = emit_drift(&mat, &opts).unwrap();
        let fluka = table.rows.iter().find(|r| r.code == Code::Fluka).unwrap();
        assert!((fluka.rel_drift - 0.8).abs() < 1e-12);
        assert_eq!(fluka.dropped.len(), 1);
        assert_eq!(fluka.dropped[0].id, NuclideId::from_name("O16").unwrap());
        assert!((table.worst_rel_drift() - 0.8).abs() < 1e-12);
        for row in table.rows.iter().filter(|r| r.code != Code::Fluka) {
            assert!(row.rel_drift.abs() < 1e-9, "{:?}", row.code);
        }
    }

    #[test]
    fn code_display_names() {
        assert_eq!(
            Code::all()
                .iter()
                .map(|c| c.to_string())
                .collect::<Vec<_>>(),
            vec!["MCNP", "Serpent", "FLUKA", "ALARA", "PARTISN"]
        );
    }

    #[test]
    fn empty_material_is_degenerate() {
        let opts = EmitOptions::new("void").with_density(1.0);
        assert!(matches!(
            emit_all(&Material::new(), &opts),
            Err(Error::Degenerate)
        ));
    }
}
