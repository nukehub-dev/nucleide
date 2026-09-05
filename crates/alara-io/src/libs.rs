//! ALARA material, element, and waste-disposal-rating support libraries.
//!
//! The `material_lib` / `element_lib` deck blocks point at the text files
//! parsed here; the older [`LibSpec`] nickname list is kept for
//! compatibility. Element identifiers keep their file spelling verbatim
//! (including enriched `li:90`-style suffixes); [`split_enriched`] and
//! [`canonical_element_z`] help callers interpret them via `nuclei`
//! without ever panicking.

use serde::{Deserialize, Serialize};
use std::path::Path;

use crate::error::{Error, Result};

/// One ALARA data-library reference.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct LibEntry {
    /// Library nickname, as referenced by `material_lib` blocks.
    pub name: String,
    /// Filesystem path or identifier of the library data.
    pub path: String,
}

/// The set of libraries a deck references.
///
/// Binary `.lib` contents stay out of scope; this type only tracks
/// nicknames and paths.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct LibSpec {
    /// Known libraries in deck order.
    pub libraries: Vec<LibEntry>,
}

impl LibSpec {
    /// Empty library set.
    pub fn new() -> Self {
        Self::default()
    }

    /// Register one library reference.
    pub fn push(&mut self, name: &str, path: &str) {
        self.libraries.push(LibEntry {
            name: name.to_string(),
            path: path.to_string(),
        });
    }

    /// Parse `NAME PATH` lines from in-memory text.
    ///
    /// Blank lines and `#` comments are skipped. Empty input is an error.
    pub fn parse(text: &str) -> Result<Self> {
        if text.trim().is_empty() {
            return Err(Error::Parse {
                line: 1,
                msg: "empty library list".to_string(),
            });
        }
        let mut spec = Self::new();
        for (index, raw) in text.lines().enumerate() {
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut words = line.split_whitespace();
            match (words.next(), words.next(), words.next()) {
                (Some(name), Some(path), None) => spec.push(name, path),
                _ => {
                    return Err(Error::Parse {
                        line: index + 1,
                        msg: format!("expected `NAME PATH`, found `{line}`"),
                    });
                }
            }
        }
        Ok(spec)
    }

    /// Look up a library by nickname.
    pub fn find(&self, name: &str) -> Option<&LibEntry> {
        self.libraries.iter().find(|entry| entry.name == name)
    }

    /// Resolve a nickname, reporting dangling references as [`Error::CrossRef`].
    pub fn resolve(&self, name: &str) -> Result<&LibEntry> {
        self.find(name)
            .ok_or_else(|| Error::CrossRef(format!("undefined library `{name}`")))
    }

    /// Number of registered libraries.
    pub fn len(&self) -> usize {
        self.libraries.len()
    }

    /// True when no libraries are registered.
    pub fn is_empty(&self) -> bool {
        self.libraries.is_empty()
    }
}

/// Split an ALARA element identifier into its base chemical symbol and
/// optional enrichment modifier: `li:90` yields `("li", Some("90"))` while
/// plain `li` yields `("li", None)`.
pub fn split_enriched(symbol: &str) -> (&str, Option<&str>) {
    match symbol.split_once(':') {
        Some((base, tag)) => (base, Some(tag)),
        None => (symbol, None),
    }
}

/// Atomic number for a free-form element symbol via `nuclei`.
///
/// ALARA libraries spell symbols lowercase (`li`, `mn:56`); the base symbol
/// (before any `:` modifier) is canonicalized before lookup. Returns `None`
/// for unrecognized symbols instead of panicking.
pub fn canonical_element_z(symbol: &str) -> Option<u32> {
    let (base, _) = split_enriched(symbol);
    if base.is_empty() {
        return None;
    }
    let mut chars = base.chars();
    let mut canonical = String::with_capacity(base.len());
    if let Some(first) = chars.next() {
        canonical.extend(first.to_uppercase());
    }
    canonical.push_str(&chars.as_str().to_lowercase());
    nuclei::element_z(&canonical).filter(|z| *z > 0)
}

/// Collect non-blank, non-`#`-comment lines with their 1-based line numbers.
fn logical_lines(text: &str) -> Vec<(usize, &str)> {
    text.lines()
        .enumerate()
        .filter_map(|(index, raw)| {
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') {
                None
            } else {
                Some((index + 1, line))
            }
        })
        .collect()
}

/// One constituent of a [`MaterialEntry`]: `<element> <weight-fraction-%> <Z>`.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct MatEntry {
    /// Element identifier spelled as in the library file (e.g. `li`).
    pub element: String,
    /// Weight fraction in percent.
    pub weight_frac: f64,
    /// Atomic number as listed in the file.
    pub z: u32,
}

impl MatEntry {
    /// Base chemical symbol (before any `:` enrichment modifier).
    pub fn base_symbol(&self) -> &str {
        split_enriched(&self.element).0
    }

    /// Atomic number for the base symbol via `nuclei`, if recognized.
    pub fn validated_z(&self) -> Option<u32> {
        canonical_element_z(&self.element)
    }
}

/// One material definition: `<name> <density> <n>` plus `n` [`MatEntry`] lines.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct MaterialEntry {
    /// Material identifier (no whitespace), e.g. `FLiBe`.
    pub name: String,
    /// Theoretical density in g/cm3.
    pub density: f64,
    /// Constituent elements in file order.
    pub entries: Vec<MatEntry>,
}

impl MaterialEntry {
    /// Look up a constituent by its file-spelled identifier.
    pub fn find(&self, element: &str) -> Option<&MatEntry> {
        self.entries.iter().find(|entry| entry.element == element)
    }

    /// Number of constituent elements.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// True when no constituents are stored.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

/// Parsed ALARA material library: a sequence of [`MaterialEntry`] blocks.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct MatLib {
    /// Material definitions in file order.
    pub materials: Vec<MaterialEntry>,
}

impl MatLib {
    /// Parse a material library from in-memory text in a single pass.
    ///
    /// Blank lines and `#` comment lines are skipped; every other line must
    /// belong to a `<name> <density> <n>` header or one of its `n`
    /// `<element> <weight-%> <Z>` entries. Errors carry 1-based lines.
    pub fn parse(text: &str) -> Result<Self> {
        if text.trim().is_empty() {
            return Err(Error::Parse {
                line: 1,
                msg: "empty material library".to_string(),
            });
        }
        let lines = logical_lines(text);
        if lines.is_empty() {
            return Err(Error::Parse {
                line: 1,
                msg: "empty material library".to_string(),
            });
        }
        let mut materials: Vec<MaterialEntry> = Vec::new();
        let mut cursor = 0;
        while cursor < lines.len() {
            let (line_no, header) = lines[cursor];
            let words: Vec<&str> = header.split_whitespace().collect();
            if words.len() != 3 {
                return Err(Error::Parse {
                    line: line_no,
                    msg: format!("expected `<name> <density> <n>`, found `{header}`"),
                });
            }
            let density: f64 = words[1].parse().map_err(|_| Error::Parse {
                line: line_no,
                msg: format!("expected material density, found `{}`", words[1]),
            })?;
            let count: usize = words[2].parse().map_err(|_| Error::Parse {
                line: line_no,
                msg: format!("expected element count, found `{}`", words[2]),
            })?;
            cursor += 1;
            if lines.len() - cursor < count {
                let at = lines.last().map_or(line_no, |(n, _)| *n);
                return Err(Error::Parse {
                    line: at,
                    msg: format!(
                        "material `{}` declares {count} elements, found {}",
                        words[0],
                        lines.len() - cursor
                    ),
                });
            }
            let mut entries = Vec::with_capacity(count.min(1024));
            for _ in 0..count {
                let (entry_line, entry_text) = lines[cursor];
                cursor += 1;
                let parts: Vec<&str> = entry_text.split_whitespace().collect();
                if parts.len() != 3 {
                    return Err(Error::Parse {
                        line: entry_line,
                        msg: format!("expected `<element> <weight-%> <Z>`, found `{entry_text}`"),
                    });
                }
                let weight_frac: f64 = parts[1].parse().map_err(|_| Error::Parse {
                    line: entry_line,
                    msg: format!("expected weight fraction, found `{}`", parts[1]),
                })?;
                let z: u32 = parts[2].parse().map_err(|_| Error::Parse {
                    line: entry_line,
                    msg: format!("expected atomic number, found `{}`", parts[2]),
                })?;
                entries.push(MatEntry {
                    element: parts[0].to_string(),
                    weight_frac,
                    z,
                });
            }
            materials.push(MaterialEntry {
                name: words[0].to_string(),
                density,
                entries,
            });
        }
        Ok(Self { materials })
    }

    /// Read and parse a material library from disk.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        let text = std::fs::read_to_string(path.as_ref())?;
        Self::parse(&text)
    }

    /// Look up a material by identifier (exact match).
    pub fn find(&self, name: &str) -> Option<&MaterialEntry> {
        self.materials.iter().find(|entry| entry.name == name)
    }

    /// Number of material definitions.
    pub fn len(&self) -> usize {
        self.materials.len()
    }

    /// True when no materials are stored.
    pub fn is_empty(&self) -> bool {
        self.materials.is_empty()
    }
}

impl std::str::FromStr for MatLib {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        Self::parse(s)
    }
}

/// One isotope of an [`ElementEntry`]: `<mass-number> <abundance-%>`.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct IsotopeAbund {
    /// Mass number.
    pub mass: u32,
    /// Atomic abundance in percent.
    pub abund_pct: f64,
}

/// One element definition: `<symbol> <atomic-weight> <Z> <density> <n>`
/// plus `n` [`IsotopeAbund`] lines.
///
/// `symbol` keeps its file spelling verbatim, including enriched
/// `li:90`-style suffixes; use [`ElementEntry::base_symbol`] and
/// [`ElementEntry::enrichment_tag`] to interpret it.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct ElementEntry {
    /// Element identifier as spelled in the file (e.g. `li`, `li:90`).
    pub symbol: String,
    /// Molar mass.
    pub atomic_weight: f64,
    /// Atomic number as listed in the file.
    pub z: u32,
    /// Theoretical density in g/cm3.
    pub density: f64,
    /// Isotopic abundances in file order.
    pub isotopes: Vec<IsotopeAbund>,
}

impl ElementEntry {
    /// Base chemical symbol (before any `:` enrichment modifier).
    pub fn base_symbol(&self) -> &str {
        split_enriched(&self.symbol).0
    }

    /// Enrichment modifier after `:`, if any (e.g. `Some("90")`).
    pub fn enrichment_tag(&self) -> Option<&str> {
        split_enriched(&self.symbol).1
    }

    /// Atomic number for the base symbol via `nuclei`, if recognized.
    pub fn validated_z(&self) -> Option<u32> {
        canonical_element_z(&self.symbol)
    }

    /// Look up an isotope by mass number.
    pub fn find(&self, mass: u32) -> Option<&IsotopeAbund> {
        self.isotopes.iter().find(|iso| iso.mass == mass)
    }

    /// Number of isotopes stored.
    pub fn len(&self) -> usize {
        self.isotopes.len()
    }

    /// True when no isotopes are stored.
    pub fn is_empty(&self) -> bool {
        self.isotopes.is_empty()
    }
}

/// Parsed ALARA element library: a sequence of [`ElementEntry`] blocks.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct EleLib {
    /// Element definitions in file order.
    pub elements: Vec<ElementEntry>,
}

impl EleLib {
    /// Parse an element library from in-memory text in a single pass.
    ///
    /// Blank lines and `#` comment lines are skipped; every other line must
    /// belong to a `<symbol> <atomic-weight> <Z> <density> <n>` header or one
    /// of its `n` `<mass> <abundance-%>` entries. Errors carry 1-based lines.
    pub fn parse(text: &str) -> Result<Self> {
        if text.trim().is_empty() {
            return Err(Error::Parse {
                line: 1,
                msg: "empty element library".to_string(),
            });
        }
        let lines = logical_lines(text);
        if lines.is_empty() {
            return Err(Error::Parse {
                line: 1,
                msg: "empty element library".to_string(),
            });
        }
        let mut elements: Vec<ElementEntry> = Vec::new();
        let mut cursor = 0;
        while cursor < lines.len() {
            let (line_no, header) = lines[cursor];
            let words: Vec<&str> = header.split_whitespace().collect();
            if words.len() != 5 {
                return Err(Error::Parse {
                    line: line_no,
                    msg: format!(
                        "expected `<symbol> <atomic-weight> <Z> <density> <n>`, found `{header}`"
                    ),
                });
            }
            let atomic_weight: f64 = words[1].parse().map_err(|_| Error::Parse {
                line: line_no,
                msg: format!("expected atomic weight, found `{}`", words[1]),
            })?;
            let z: u32 = words[2].parse().map_err(|_| Error::Parse {
                line: line_no,
                msg: format!("expected atomic number, found `{}`", words[2]),
            })?;
            let density: f64 = words[3].parse().map_err(|_| Error::Parse {
                line: line_no,
                msg: format!("expected density, found `{}`", words[3]),
            })?;
            let count: usize = words[4].parse().map_err(|_| Error::Parse {
                line: line_no,
                msg: format!("expected isotope count, found `{}`", words[4]),
            })?;
            cursor += 1;
            if lines.len() - cursor < count {
                let at = lines.last().map_or(line_no, |(n, _)| *n);
                return Err(Error::Parse {
                    line: at,
                    msg: format!(
                        "element `{}` declares {count} isotopes, found {}",
                        words[0],
                        lines.len() - cursor
                    ),
                });
            }
            let mut isotopes = Vec::with_capacity(count.min(1024));
            for _ in 0..count {
                let (iso_line, iso_text) = lines[cursor];
                cursor += 1;
                let parts: Vec<&str> = iso_text.split_whitespace().collect();
                if parts.len() != 2 {
                    return Err(Error::Parse {
                        line: iso_line,
                        msg: format!("expected `<mass> <abundance-%>`, found `{iso_text}`"),
                    });
                }
                let mass: u32 = parts[0].parse().map_err(|_| Error::Parse {
                    line: iso_line,
                    msg: format!("expected mass number, found `{}`", parts[0]),
                })?;
                let abund_pct: f64 = parts[1].parse().map_err(|_| Error::Parse {
                    line: iso_line,
                    msg: format!("expected abundance, found `{}`", parts[1]),
                })?;
                isotopes.push(IsotopeAbund { mass, abund_pct });
            }
            elements.push(ElementEntry {
                symbol: words[0].to_string(),
                atomic_weight,
                z,
                density,
                isotopes,
            });
        }
        Ok(Self { elements })
    }

    /// Read and parse an element library from disk.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        let text = std::fs::read_to_string(path.as_ref())?;
        Self::parse(&text)
    }

    /// Look up an element by its file-spelled symbol (exact match).
    pub fn find(&self, symbol: &str) -> Option<&ElementEntry> {
        self.elements.iter().find(|entry| entry.symbol == symbol)
    }

    /// Number of element definitions.
    pub fn len(&self) -> usize {
        self.elements.len()
    }

    /// True when no elements are stored.
    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }
}

impl std::str::FromStr for EleLib {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        Self::parse(s)
    }
}

/// One waste-disposal-rating / clearance-index limit: `<zaid-key> <limit>`.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct WdrLimit {
    /// Isotope key spelled as in the file (kza `ZZAAAM`, e.g. `60140`).
    pub zaid_key: String,
    /// Disposal limit (specific activity in the file's units).
    pub limit: f64,
}

/// Parsed ALARA waste-disposal-rating / clearance-index file (NRCA/NRCC).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct WdrLib {
    /// Library name (from [`WdrLib::parse`] or the source file stem).
    pub name: String,
    /// Limits in file order.
    pub limits: Vec<WdrLimit>,
}

impl WdrLib {
    /// Parse `<zaid-key> <limit>` lines from in-memory text in a single pass.
    ///
    /// Blank lines and `#` comment lines are skipped. Errors carry 1-based
    /// line numbers.
    pub fn parse(name: &str, text: &str) -> Result<Self> {
        if text.trim().is_empty() {
            return Err(Error::Parse {
                line: 1,
                msg: format!("empty waste-disposal library `{name}`"),
            });
        }
        let mut limits: Vec<WdrLimit> = Vec::new();
        for (index, raw) in text.lines().enumerate() {
            let line_no = index + 1;
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() != 2 {
                return Err(Error::Parse {
                    line: line_no,
                    msg: format!("expected `<zaid-key> <limit>`, found `{line}`"),
                });
            }
            let limit: f64 = parts[1].parse().map_err(|_| Error::Parse {
                line: line_no,
                msg: format!("expected disposal limit, found `{}`", parts[1]),
            })?;
            limits.push(WdrLimit {
                zaid_key: parts[0].to_string(),
                limit,
            });
        }
        if limits.is_empty() {
            return Err(Error::Parse {
                line: 1,
                msg: format!("empty waste-disposal library `{name}`"),
            });
        }
        Ok(Self {
            name: name.to_string(),
            limits,
        })
    }

    /// Read and parse a waste-disposal file from disk, naming it after the
    /// file stem (e.g. `NRCA`).
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        let text = std::fs::read_to_string(path.as_ref())?;
        let name = path
            .as_ref()
            .file_stem()
            .and_then(|stem| stem.to_str())
            .unwrap_or_default();
        Self::parse(name, &text)
    }

    /// Look up a limit by its file-spelled isotope key (exact match).
    pub fn find(&self, zaid_key: &str) -> Option<&WdrLimit> {
        self.limits.iter().find(|entry| entry.zaid_key == zaid_key)
    }

    /// Number of limits stored.
    pub fn len(&self) -> usize {
        self.limits.len()
    }

    /// True when no limits are stored.
    pub fn is_empty(&self) -> bool {
        self.limits.is_empty()
    }
}

impl std::str::FromStr for WdrLib {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        Self::parse("", s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    fn fixture(name: &str) -> String {
        let path = format!(
            "{}/../../fixtures/alara/libs/{name}",
            env!("CARGO_MANIFEST_DIR")
        );
        std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("read fixture {path}: {error}"))
    }

    #[test]
    fn parse_registers_libraries() {
        let spec = LibSpec::parse("liba /data/a\n# comment\nlibb /data/b\n").unwrap();
        assert_eq!(spec.len(), 2);
        assert!(!spec.is_empty());
        assert_eq!(spec.find("liba").unwrap().path, "/data/a");
        assert_eq!(spec.resolve("libb").unwrap().path, "/data/b");
    }

    #[test]
    fn parse_rejects_malformed_lines() {
        assert!(matches!(
            LibSpec::parse("onlyone\n"),
            Err(Error::Parse { line: 1, .. })
        ));
        assert!(matches!(
            LibSpec::parse("a b c\n"),
            Err(Error::Parse { .. })
        ));
    }

    #[test]
    fn resolve_reports_dangling_references() {
        let spec = LibSpec::new();
        assert!(spec.is_empty());
        assert!(matches!(spec.resolve("missing"), Err(Error::CrossRef(_))));
    }

    #[test]
    fn split_enriched_separates_base_and_tag() {
        assert_eq!(split_enriched("li"), ("li", None));
        assert_eq!(split_enriched("li:90"), ("li", Some("90")));
        assert_eq!(split_enriched("mn:56"), ("mn", Some("56")));
    }

    #[test]
    fn canonical_element_z_validates_without_panicking() {
        assert_eq!(canonical_element_z("li"), Some(3));
        assert_eq!(canonical_element_z("Li"), Some(3));
        assert_eq!(canonical_element_z("li:90"), Some(3));
        assert_eq!(canonical_element_z("mn:56"), Some(25));
        assert_eq!(canonical_element_z("xx"), None);
        assert_eq!(canonical_element_z(""), None);
    }

    #[test]
    fn matlib_parses_blocks_with_counts() {
        let lib = MatLib::parse("WATER 1.0 2\nh 11.111 1\no 88.889 8\n").unwrap();
        assert_eq!(lib.len(), 1);
        let water = lib.find("WATER").unwrap();
        assert_eq!(water.density, 1.0);
        assert_eq!(water.len(), 2);
        assert_eq!(
            water.find("h").unwrap(),
            &MatEntry {
                element: "h".to_string(),
                weight_frac: 11.111,
                z: 1,
            }
        );
        assert_eq!(water.find("o").unwrap().validated_z(), Some(8));
        assert!(lib.find("MISSING").is_none());
    }

    #[test]
    fn matlib_reports_bad_floats_with_line_numbers() {
        let error = MatLib::parse("WATER 1.0 1\nh nope 1\n").unwrap_err();
        assert!(matches!(error, Error::Parse { line: 2, .. }));
        let error = MatLib::parse("WATER nope 1\nh 1.0 1\n").unwrap_err();
        assert!(matches!(error, Error::Parse { line: 1, .. }));
    }

    #[test]
    fn matlib_reports_truncated_blocks() {
        let error = MatLib::parse("WATER 1.0 2\nh 11.111 1\n").unwrap_err();
        assert!(matches!(error, Error::Parse { .. }));
    }

    #[test]
    fn matlib_rejects_empty_and_malformed_headers() {
        assert!(matches!(MatLib::parse("  \n"), Err(Error::Parse { .. })));
        assert!(matches!(
            MatLib::parse("# only a comment\n"),
            Err(Error::Parse { .. })
        ));
        assert!(matches!(
            MatLib::parse("WATER 1.0\n"),
            Err(Error::Parse { line: 1, .. })
        ));
        assert!(matches!(
            MatLib::from_str("WATER 1.0 2\nh 11.111 1\no 88.889 8\n")
                .unwrap()
                .len(),
            1
        ));
    }

    #[test]
    fn elelib_parses_enriched_symbols_verbatim() {
        let lib = EleLib::parse("li:90 6.11521 3 0.53 2\n6 90.0\n7 10.0\n").unwrap();
        assert_eq!(lib.len(), 1);
        let li = lib.find("li:90").unwrap();
        assert_eq!(li.symbol, "li:90");
        assert_eq!(li.base_symbol(), "li");
        assert_eq!(li.enrichment_tag(), Some("90"));
        assert_eq!(li.z, 3);
        assert_eq!(li.validated_z(), Some(3));
        assert_eq!(li.len(), 2);
        assert_eq!(li.find(6).unwrap().abund_pct, 90.0);
        let natural = EleLib::parse("li 6.941 3 0.53 1\n7 100.0\n").unwrap();
        assert_eq!(natural.find("li").unwrap().enrichment_tag(), None);
    }

    #[test]
    fn elelib_reports_bad_entries_with_line_numbers() {
        let error = EleLib::parse("li 6.941 three 0.53 1\n7 100.0\n").unwrap_err();
        assert!(matches!(error, Error::Parse { line: 1, .. }));
        let error = EleLib::parse("li 6.941 3 0.53 1\n7 nope\n").unwrap_err();
        assert!(matches!(error, Error::Parse { line: 2, .. }));
        assert!(matches!(EleLib::parse(""), Err(Error::Parse { .. })));
    }

    #[test]
    fn wdrlib_parses_key_limit_pairs() {
        let lib = WdrLib::parse("test", "10030 4.0e+02\n60140 8.0e+00\n").unwrap();
        assert_eq!(lib.name, "test");
        assert_eq!(lib.len(), 2);
        assert!(!lib.is_empty());
        assert_eq!(lib.find("10030").unwrap().limit, 400.0);
        assert!(lib.find("99999").is_none());
        assert!(matches!(
            WdrLib::from_str("10030 4.0e+02\n").unwrap().len(),
            1
        ));
    }

    #[test]
    fn wdrlib_rejects_malformed_lines() {
        assert!(matches!(
            WdrLib::parse("t", "10030\n"),
            Err(Error::Parse { line: 1, .. })
        ));
        assert!(matches!(
            WdrLib::parse("t", "10030 nope\n"),
            Err(Error::Parse { line: 1, .. })
        ));
        assert!(matches!(
            WdrLib::parse("t", "  \n"),
            Err(Error::Parse { .. })
        ));
    }

    #[test]
    fn fixture_matlib_holds_expected_materials() {
        let lib = MatLib::from_file(format!(
            "{}/../../fixtures/alara/libs/sampleMatlib",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap();
        assert_eq!(lib.len(), 8);
        for name in ["Li", "FLiBe", "WATER", "C1020", "CONC", "b4c"] {
            assert!(lib.find(name).is_some(), "missing material {name}");
        }
        let li = lib.find("Li").unwrap();
        assert_eq!(li.density, 1.0);
        assert_eq!(li.len(), 11);
        assert_eq!(li.find("li").unwrap().weight_frac, 99.901);
        let water = lib.find("WATER").unwrap();
        assert_eq!(water.len(), 2);
        let conc = lib.find("CONC").unwrap();
        assert_eq!(conc.len(), 13);
        let steel = lib.find("SS-304").unwrap();
        assert_eq!(steel.len(), 11);
        assert_eq!(steel.find("fe").unwrap().weight_frac, 70.578);
    }

    #[test]
    fn fixture_elelib_holds_expected_elements() {
        let lib = EleLib::from_file(format!(
            "{}/../../fixtures/alara/libs/myElelib",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap();
        assert_eq!(lib.len(), 86);
        let li90 = lib.find("li:90").unwrap();
        assert_eq!(li90.symbol, "li:90");
        assert_eq!(li90.base_symbol(), "li");
        assert_eq!(li90.enrichment_tag(), Some("90"));
        assert_eq!(li90.z, 3);
        assert!((li90.atomic_weight - 6.11521).abs() < 1e-5);
        assert_eq!(li90.len(), 2);
        assert_eq!(li90.find(6).unwrap().abund_pct, 90.0);
        assert_eq!(li90.find(7).unwrap().abund_pct, 10.0);
        assert!(lib.find("mn:56").is_some());
        assert!(lib.find("fe:56").is_some());
        let h = lib.find("h").unwrap();
        assert_eq!(h.enrichment_tag(), None);
        assert_eq!(h.validated_z(), Some(1));
    }

    #[test]
    fn fixture_wdr_files_hold_expected_limits() {
        let nrca = WdrLib::from_file(format!(
            "{}/../../fixtures/alara/libs/NRCA",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap();
        assert_eq!(nrca.name, "NRCA");
        assert_eq!(nrca.len(), 10);
        assert_eq!(nrca.find("10030").unwrap().limit, 400.0);
        assert_eq!(nrca.find("551370").unwrap().limit, 10.0);

        let nrcc = WdrLib::from_file(format!(
            "{}/../../fixtures/alara/libs/NRCC",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap();
        assert_eq!(nrcc.name, "NRCC");
        assert_eq!(nrcc.len(), 9);
        assert_eq!(nrcc.find("60140").unwrap().limit, 80.0);
    }

    #[test]
    fn fixture_text_matches_parse_from_string() {
        assert_eq!(MatLib::parse(&fixture("sampleMatlib")).unwrap().len(), 8);
        assert_eq!(EleLib::parse(&fixture("myElelib")).unwrap().len(), 86);
        assert_eq!(WdrLib::parse("NRCA", &fixture("NRCA")).unwrap().len(), 10);
        assert_eq!(WdrLib::parse("NRCC", &fixture("NRCC")).unwrap().len(), 9);
    }
}
