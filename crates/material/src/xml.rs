//! Materials XML export.
//!
//! Emits `<material>` elements carrying a `name` attribute, a `<density>`
//! child with `value`/`units` attributes (units default to `g/cm3`), and one
//! `<nuclide name="U235" wo="..."/>` child per component. Weight fractions
//! are written through the `wo` attribute; atom fractions would use `ao`.
//! Nuclides are emitted in ascending nucid order and subnormal values are
//! clamped to zero before writing.

use std::io::Write;

use quick_xml::events::{BytesDecl, BytesEnd, BytesStart, Event};
use quick_xml::Writer;

use crate::{Error, Material};

/// Smallest positive normal f64; smaller magnitudes are clamped to zero so
/// downstream readers never see a subnormal.
const SMALLEST_NORMAL: f64 = 2.225_073_858_507_201_4e-308;

/// Default density units.
pub const DEFAULT_DENSITY_UNITS: &str = "g/cm3";

/// Format a float the way Python's `str(float)` would, with subnormals
/// clamped to zero.
fn fmt_num(v: f64) -> String {
    let v = if v != 0.0 && v.abs() < SMALLEST_NORMAL {
        0.0
    } else {
        v
    };
    let s = format!("{v}");
    if v.is_finite() && !s.contains(['.', 'e', 'E']) {
        format!("{s}.0")
    } else {
        s
    }
}

impl Material {
    /// Serialize this material as a `<material>` XML fragment.
    ///
    /// Components are written as weight fractions (`wo` attributes). The
    /// density and its units are taken from the arguments rather than from
    /// [`Material::density`], matching the free-standing export style.
    pub fn to_xml(&self, name: &str, density: f64, units: &str) -> crate::Result<String> {
        let mut writer = Writer::new_with_indent(Vec::<u8>::new(), b' ', 2);
        write_material(&mut writer, self, name, density, units)?;
        Ok(String::from_utf8_lossy(&writer.into_inner()).into_owned())
    }
}

fn write_material<W: Write>(
    writer: &mut Writer<W>,
    mat: &Material,
    name: &str,
    density: f64,
    units: &str,
) -> crate::Result<()> {
    let fractions = mat.weight_fractions()?;

    let mut root = BytesStart::new("material");
    root.push_attribute(("name", name));
    writer.write_event(Event::Start(root))?;

    let mut den = BytesStart::new("density");
    let value = fmt_num(density);
    den.push_attribute(("value", value.as_str()));
    den.push_attribute(("units", units));
    writer.write_event(Event::Empty(den))?;

    for (id, wf) in fractions {
        let mut nuc = BytesStart::new("nuclide");
        let nuc_name = id.to_name();
        let wo = fmt_num(wf);
        nuc.push_attribute(("name", nuc_name.as_str()));
        nuc.push_attribute(("wo", wo.as_str()));
        writer.write_event(Event::Empty(nuc))?;
    }

    writer.write_event(Event::End(BytesEnd::new("material")))?;
    Ok(())
}

/// A `<materials>` document bundling named [`Material`]s, optionally pointing
/// at a cross-sections file via the root `cross_sections` attribute — the
/// container shape of a materials collection document.
#[derive(Debug, Clone, Default)]
pub struct MaterialsDoc {
    cross_sections: Option<String>,
    entries: Vec<(String, Material)>,
}

impl MaterialsDoc {
    /// An empty document.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the `cross_sections` attribute on the root element.
    pub fn cross_sections(mut self, path: impl Into<String>) -> Self {
        self.cross_sections = Some(path.into());
        self
    }

    /// Append a named material.
    pub fn push(mut self, name: impl Into<String>, material: Material) -> Self {
        self.entries.push((name.into(), material));
        self
    }

    /// Number of materials in the document.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// True when the document holds no materials.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Serialize the full `<materials>` document, including the XML
    /// declaration. Each material uses its stored density (in
    /// [`DEFAULT_DENSITY_UNITS`]); missing densities are an error.
    pub fn to_xml(&self) -> crate::Result<String> {
        let mut writer = Writer::new_with_indent(Vec::<u8>::new(), b' ', 2);
        writer.write_event(Event::Decl(BytesDecl::new("1.0", Some("UTF-8"), None)))?;

        let mut root = BytesStart::new("materials");
        if let Some(path) = &self.cross_sections {
            root.push_attribute(("cross_sections", path.as_str()));
        }
        writer.write_event(Event::Start(root))?;

        for (name, mat) in &self.entries {
            let rho = mat.density().ok_or(Error::MissingDensity)?;
            write_material(&mut writer, mat, name, rho, DEFAULT_DENSITY_UNITS)?;
        }

        writer.write_event(Event::End(BytesEnd::new("materials")))?;
        Ok(String::from_utf8_lossy(&writer.into_inner()).into_owned())
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use quick_xml::events::Event;
    use quick_xml::Reader;

    use super::*;

    fn id(name: &str) -> nuclei::NuclideId {
        nuclei::NuclideId::from_name(name).unwrap()
    }

    fn uranium() -> Material {
        let mut mat = Material::new();
        mat.add_nuclide(id("U235"), 19.0);
        mat.add_nuclide(id("U238"), 1.0);
        mat
    }

    /// Collect (tag, attributes) pairs from an XML document.
    fn elements(xml: &str) -> Vec<(String, BTreeMap<String, String>)> {
        let mut reader = Reader::from_str(xml);
        reader.config_mut().trim_text(true);
        let mut out = Vec::new();
        let mut event = reader.read_event().unwrap();
        while event != Event::Eof {
            if let Event::Empty(bs) | Event::Start(bs) = &event {
                let attrs: BTreeMap<String, String> = bs
                    .attributes()
                    .map(|a| {
                        let a = a.unwrap();
                        (
                            String::from_utf8_lossy(a.key.as_ref()).into_owned(),
                            a.decode_and_unescape_value(reader.decoder())
                                .unwrap()
                                .into_owned(),
                        )
                    })
                    .collect();
                out.push((
                    String::from_utf8_lossy(bs.name().as_ref()).into_owned(),
                    attrs,
                ));
            }
            event = reader.read_event().unwrap();
        }
        out
    }

    fn map(pairs: &[(&str, &str)]) -> BTreeMap<String, String> {
        pairs
            .iter()
            .map(|&(k, v)| (k.to_string(), v.to_string()))
            .collect()
    }

    #[test]
    fn fragment_shape() {
        let xml = uranium().to_xml("uo2", 19.1, "g/cm3").unwrap();
        assert_eq!(
            xml,
            "<material name=\"uo2\">\n  \
             <density value=\"19.1\" units=\"g/cm3\"/>\n  \
             <nuclide name=\"U235\" wo=\"0.95\"/>\n  \
             <nuclide name=\"U238\" wo=\"0.05\"/>\n\
             </material>"
        );
    }

    #[test]
    fn fragment_attributes_are_well_formed() {
        let xml = uranium().to_xml("u", 10.0, "g/cm3").unwrap();
        let elems = elements(&xml);

        assert_eq!(elems[0], ("material".to_string(), map(&[("name", "u")])));
        assert_eq!(
            elems[1],
            (
                "density".to_string(),
                map(&[("value", "10.0"), ("units", "g/cm3")])
            )
        );
        for (tag, attrs) in &elems[2..] {
            assert_eq!(tag, "nuclide");
            assert_eq!(attrs.len(), 2);
            assert!(attrs.contains_key("name"));
            assert!(attrs.contains_key("wo"), "weight fraction attr missing");
        }
    }

    #[test]
    fn fragment_of_empty_material_errors() {
        let err = Material::new().to_xml("void", 1.0, "g/cm3").unwrap_err();
        assert!(matches!(err, Error::Degenerate));
    }

    #[test]
    fn subnormal_values_are_clamped_to_zero() {
        assert_eq!(fmt_num(1e-320), "0.0");
        assert_eq!(fmt_num(-0.0), "-0.0");
        assert_eq!(fmt_num(1.0), "1.0");
        assert_eq!(fmt_num(0.95), "0.95");
        assert_eq!(fmt_num(f64::INFINITY), "inf");
    }

    #[test]
    fn materials_doc_emits_cross_sections_root() {
        let mut water = Material::new();
        water.add_nuclide(id("H1"), 2.0);
        water.set_density(Some(0.998));

        let doc = MaterialsDoc::new()
            .cross_sections("/data/cross_sections.xml")
            .push("water", water);

        let xml = doc.to_xml().unwrap();
        assert_eq!(
            xml,
            "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
             <materials cross_sections=\"/data/cross_sections.xml\">\n  \
             <material name=\"water\">\n    \
             <density value=\"0.998\" units=\"g/cm3\"/>\n    \
             <nuclide name=\"H1\" wo=\"1.0\"/>\n  \
             </material>\n\
             </materials>"
        );
    }

    #[test]
    fn materials_doc_requires_density_per_material() {
        let bare = uranium(); // no density set
        let err = MaterialsDoc::new().push("u", bare).to_xml().unwrap_err();
        assert!(matches!(err, Error::MissingDensity));
    }
}
