//! Depletion chain handling: depletion-chain XML format.

use std::collections::BTreeMap;
use std::fmt;
use std::path::Path;

/// One decay mode (e.g. beta → Xe135, branching ratio 1.0).
#[derive(Debug, Clone, PartialEq)]
pub struct DecayMode {
    /// ENDF decay type ("beta", "alpha", "sf", ...).
    pub kind: String,
    /// Daughter nuclide name.
    pub target: String,
    pub branching_ratio: f64,
}

/// One transmutation reaction channel.
#[derive(Debug, Clone, PartialEq)]
pub struct Reaction {
    /// Reaction name as written in chain files ("(n,gamma)", "fission", "(n,2n)", ...).
    pub kind: String,
    /// Daughter nuclide; `None` means pure loss ("Nothing").
    pub target: Option<String>,
    /// Q value [eV].
    pub q: f64,
}

/// Fission product yields at one incident energy.
#[derive(Debug, Clone, PartialEq)]
pub struct FissionYields {
    pub energy: f64,
    /// Product nuclide → independent yield fraction.
    pub products: BTreeMap<String, f64>,
}

/// A single chain nuclide with all its transmutation channels.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ChainNuclide {
    pub name: String,
    /// Half-life [s]; `None` for stable nuclides.
    pub half_life: Option<f64>,
    /// Decay energy [eV].
    pub decay_energy: f64,
    pub decay_modes: Vec<DecayMode>,
    pub reactions: Vec<Reaction>,
    pub neutron_fission_yields: Vec<FissionYields>,
}

impl ChainNuclide {
    /// Decay constant λ = ln(2)/T½ [1/s]; stable nuclides give 0.
    pub fn decay_constant(&self) -> f64 {
        self.half_life
            .map(|t| std::f64::consts::LN_2 / t)
            .unwrap_or(0.0)
    }
}

/// Errors from chain parsing and matrix construction.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    Io(String),
    Xml(String),
    /// Reference to a nuclide absent from the chain.
    UnknownNuclide {
        name: String,
        context: &'static str,
    },
    BadStructure(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(m) => write!(f, "io error: {m}"),
            Error::Xml(m) => write!(f, "XML error: {m}"),
            Error::UnknownNuclide { name, context } => {
                write!(f, "unknown nuclide `{name}` ({context})")
            }
            Error::BadStructure(m) => write!(f, "malformed chain: {m}"),
        }
    }
}
impl std::error::Error for Error {}

/// Parsed depletion chain in topological file order.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Chain {
    pub nuclides: Vec<ChainNuclide>,
    index: BTreeMap<String, usize>,
}

impl Chain {
    /// Read a chain XML file.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self, Error> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path)
            .map_err(|e| Error::Io(format!("{}: {}", path.display(), e)))?;
        Chain::from_xml(&text)
    }

    /// Parse chain XML text (depletion-chain format).
    pub fn from_xml(text: &str) -> Result<Self, Error> {
        let doc = roxmltree::Document::parse(text).map_err(|e| Error::Xml(e.to_string()))?;
        let root = doc.root_element();
        if !root.has_tag_name("depletion_chain") {
            return Err(Error::BadStructure(format!(
                "expected <depletion_chain> root, found <{}>",
                root.tag_name().name()
            )));
        }

        let mut nuclides = Vec::new();
        for n_el in root.children().filter(|c| c.is_element()) {
            if n_el.tag_name().name() != "nuclide" {
                continue;
            }
            let mut nuc = ChainNuclide {
                name: attr(n_el, "name")?.to_string(),
                ..Default::default()
            };
            if let Some(t) = opt_attr(n_el, "half_life") {
                nuc.half_life = Some(parse_f64("half_life", t)?);
            }
            if let Some(t) = opt_attr(n_el, "decay_energy") {
                nuc.decay_energy = parse_f64("decay_energy", t)?;
            }
            for child in n_el.children().filter(|c| c.is_element()) {
                match child.tag_name().name() {
                    "decay" => {
                        nuc.decay_modes.push(DecayMode {
                            kind: attr(child, "type")?.to_string(),
                            target: attr(child, "target")?.to_string(),
                            branching_ratio: parse_f64(
                                "branching_ratio",
                                opt_attr(child, "branching_ratio").unwrap_or("1.0"),
                            )?,
                        });
                    }
                    "reaction" => {
                        let target = opt_attr(child, "target").and_then(|t| {
                            if t.eq_ignore_ascii_case("nothing") {
                                None
                            } else {
                                Some(t.to_string())
                            }
                        });
                        nuc.reactions.push(Reaction {
                            kind: attr(child, "type")?.to_string(),
                            target,
                            q: parse_f64("Q", opt_attr(child, "Q").unwrap_or("0.0"))?,
                        });
                    }
                    "neutron_fission_yields" => {
                        nuc.neutron_fission_yields = parse_fission_yields(child)?;
                    }
                    _ => {} // source / other metadata ignored
                }
            }
            nuclides.push(nuc);
        }

        if nuclides.is_empty() {
            return Err(Error::BadStructure("chain has no nuclides".into()));
        }

        Self::from_nuclides(nuclides)
    }

    /// Build a chain directly from nuclides (index computed here).
    pub fn from_nuclides(nuclides: Vec<ChainNuclide>) -> Result<Self, Error> {
        let mut chain = Chain {
            nuclides,
            index: BTreeMap::new(),
        };
        for (i, n) in chain.nuclides.iter().enumerate() {
            if chain.index.insert(n.name.clone(), i).is_some() {
                return Err(Error::BadStructure(format!(
                    "duplicate nuclide `{}`",
                    n.name
                )));
            }
        }
        Ok(chain)
    }

    /// Index of a nuclide by name.
    pub fn index_of(&self, name: &str) -> Option<usize> {
        self.index.get(name).copied()
    }

    /// Number of nuclides (matrix dimension).
    pub fn len(&self) -> usize {
        self.nuclides.len()
    }

    pub fn is_empty(&self) -> bool {
        self.nuclides.is_empty()
    }
}

fn parse_fission_yields(el: roxmltree::Node) -> Result<Vec<FissionYields>, Error> {
    // energies list + one <fission_yields energy="..."> block per energy
    let energies_text = el
        .children()
        .find(|c| c.is_element() && c.tag_name().name() == "energies")
        .and_then(|c| c.text())
        .ok_or_else(|| Error::BadStructure("fission yields missing <energies>".into()))?;
    let _energies: Vec<f64> = energies_text
        .split_whitespace()
        .map(|t| parse_f64("fy energy", t))
        .collect::<Result<_, _>>()?;

    let mut out = Vec::new();
    for fy in el
        .children()
        .filter(|c| c.is_element() && c.tag_name().name() == "fission_yields")
    {
        let energy = parse_f64("fy energy", fy.attribute("energy").unwrap_or_default())?;
        let products = fy
            .children()
            .find(|c| c.is_element() && c.tag_name().name() == "products")
            .and_then(|c| c.text())
            .unwrap_or_default()
            .split_whitespace()
            .map(String::from)
            .collect::<Vec<_>>();
        let data = fy
            .children()
            .find(|c| c.is_element() && c.tag_name().name() == "data")
            .and_then(|c| c.text())
            .unwrap_or_default()
            .split_whitespace()
            .map(|t| parse_f64("fy data", t))
            .collect::<Result<Vec<_>, _>>()?;
        if products.len() != data.len() {
            return Err(Error::BadStructure(format!(
                "fission yields products/data length mismatch: {} vs {}",
                products.len(),
                data.len()
            )));
        }
        out.push(FissionYields {
            energy,
            products: products.into_iter().zip(data).collect(),
        });
    }
    Ok(out)
}

fn attr<'a>(el: roxmltree::Node<'a, 'a>, name: &str) -> Result<&'a str, Error> {
    el.attribute(name).ok_or_else(|| {
        Error::Xml(format!(
            "<{}> missing attribute `{}`",
            el.tag_name().name(),
            name
        ))
    })
}

fn opt_attr<'a>(el: roxmltree::Node<'a, 'a>, name: &str) -> Option<&'a str> {
    el.attribute(name)
}

fn parse_f64(context: &'static str, text: &str) -> Result<f64, Error> {
    text.trim()
        .parse::<f64>()
        .map_err(|_| Error::Xml(format!("bad float for {context}: `{text}`")))
}
