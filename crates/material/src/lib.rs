//! Materials: compositions, mixing, and serialization for nuclear engineering.
//!
//! A [`Material`] is a map from [`NuclideId`] to a stored mass in grams,
//! plus an optional mass density (g/cm3 by convention) and free-form JSON
//! metadata. Density lives outside the composition: it is a property of
//! the stream rather than part of the composition itself.
//!
//! Atomic-mass-dependent conversions ([`Material::from_atom_frac`] and
//! [`Material::atom_fractions`]) take a [`MassProvider`] so the material
//! crate stays independent of the nuclear-data tables; wire up
//! `nuclei::data` once it lands.
//!
//! ```
//! use material::{MassProvider, Material, NoMasses};
//! use nuclei::NuclideId;
//!
//! let mut mat = Material::new();
//! mat.add_nuclide(NuclideId::from_name("U235").unwrap(), 19.0);
//! mat.add_nuclide(NuclideId::from_name("U238").unwrap(), 1.0);
//! assert_eq!(mat.mass(), 20.0);
//!
//! // Atom conversions need atomic masses:
//! let atoms = mat.atom_fractions(&NoMasses);
//! assert!(atoms.is_err());
//! ```

mod compendium;
mod expansion;
mod material;
mod xml;

pub use compendium::{
    CompendiumElement, CompendiumEntry, CompendiumIsotope, Error as CompendiumError,
    MaterialsLibrary,
};
pub use expansion::{
    parse_formula, AbundanceProvider, FormulaError, FormulaResult, NaturalAbundances, NoAbundances,
};
pub use material::{
    Ame2020, Analytics, AnalyticsError, ChainDecays, DecayProvider, MassProvider, Material,
    NoDecay, NoMasses, AVOGADRO, GRAMS_PER_U,
};
pub use xml::MaterialsDoc;

use nuclei::NuclideId;
use thiserror::Error;

/// Result alias for the material crate.
pub type Result<T> = std::result::Result<T, Error>;

/// Errors produced by material construction, conversion, and export.
#[derive(Debug, Error)]
pub enum Error {
    /// A nuclide name could not be parsed.
    #[error("invalid nuclide name `{name}`")]
    BadNuclide {
        /// The rejected name.
        name: String,
        /// Underlying parsing error from the nuclei crate.
        #[source]
        source: nuclei::Error,
    },
    /// An atomic mass was required but not supplied.
    #[error("no atomic mass available for nuclide `{0}`")]
    MissingMass(NuclideId),
    /// The composition is empty or its masses sum to a non-positive value.
    #[error("material is empty or its masses sum to a non-positive value")]
    Degenerate,
    /// A volume-based operation hit a material without a density.
    #[error("operation requires a mass density but none was set")]
    MissingDensity,
    /// A mixing fraction was negative.
    #[error("negative mixing fraction `{0}`")]
    NegativeFraction(f64),
    /// Writing the XML document failed.
    #[error(transparent)]
    Write(#[from] std::io::Error),
}
