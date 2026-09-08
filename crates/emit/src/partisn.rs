//! PARTISN single-zone deck emission.
//!
//! The material becomes one [`PartisnDeck`] zone carrying an ISOTXS label
//! per nuclide (uppercase canonical names, matching the fixture style).
//! PARTISN has no deck reader in this workspace — and [`validate`] needs an
//! [`IsotxsLib`] the emitter does not have — so drift is analytic
//! ([`Emitted::reparsed`] is always false).
//!
//! [`validate`]: nucleide_cccc_io::PartisnDeck::validate
//! [`IsotxsLib`]: nucleide_cccc_io::IsotxsLib

use nucleide_cccc_io::{partisn::PartisnZone, PartisnDeck};
use nucleide_material::Material;

use crate::{Code, EmitOptions, Emitted, Error, Result};

/// Emit `mat` as a single-zone PARTISN deck.
pub fn emit_partisn(mat: &Material, opts: &EmitOptions) -> Result<Emitted> {
    let density = opts.density_for(mat)?;
    let fracs = mat.weight_fractions().map_err(|_| Error::Degenerate)?;
    let mut accounted = Vec::with_capacity(fracs.len());
    let mut labels = Vec::with_capacity(fracs.len());
    for id in fracs.keys() {
        labels.push(id.to_name().to_ascii_uppercase());
        accounted.push((*id, mat.comp.get(id).copied().unwrap_or(0.0)));
    }
    let deck = PartisnDeck {
        title: opts.name.clone(),
        dim: 1,
        zones: vec![PartisnZone {
            id: opts.partisn_zone,
            material: opts.name.clone(),
            isotxs_labels: labels,
            density,
        }],
        source: None,
    };
    Ok(Emitted {
        code: Code::Partisn,
        text: deck.render(),
        accounted,
        dropped: Vec::new(),
        reparsed: false,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use nucleide_nuclei::NuclideId;

    #[test]
    fn partisn_deck_shape() {
        let mut mat = Material::new();
        mat.add_nuclide(NuclideId::from_name("U235").unwrap(), 1.0);
        let opts = EmitOptions::new("fuel").with_density(10.0);
        let out = emit_partisn(&mat, &opts).unwrap();
        assert!(out.text.contains("TITLE fuel"));
        assert!(out
            .text
            .contains("ZONE 1 MATERIAL fuel DENSITY 10 ISOTXS U235"));
        assert!(out.text.ends_with("END\n"));
        assert!(!out.reparsed);
        assert!((out.mass_out() - 1.0).abs() < 1e-12);
    }
}
