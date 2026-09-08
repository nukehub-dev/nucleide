//! FLUKA `COMPOUND` card emission.
//!
//! The whole material becomes one compound via [`compound_str`] with mass
//! fractions. Nuclides without a FLUKA name or atomic mass are reported in
//! [`Emitted::dropped`] instead of failing the emission. FLUKA has no
//! material reader in this workspace, so drift is analytic
//! ([`Emitted::reparsed`] is always false).

use nucleide_fluka_io::material::{compound_str, Component, FlukaNuc, FracType};

use crate::{Code, Dropped, EmitOptions, Emitted, Error, Result};
use nucleide_material::Material;

/// Emit `mat` as FLUKA `MATERIAL` + `COMPOUND` cards.
pub fn emit_fluka(mat: &Material, opts: &EmitOptions) -> Result<Emitted> {
    let density = opts.density_for(mat)?;
    let fracs = mat.weight_fractions().map_err(|_| Error::Degenerate)?;
    let mut accounted = Vec::with_capacity(fracs.len());
    let mut dropped = Vec::new();
    let mut components = Vec::with_capacity(fracs.len());
    for (id, w) in &fracs {
        let mass = mat.comp.get(id).copied().unwrap_or(0.0);
        let nuc = FlukaNuc::from(*id);
        if let Err(e) = nuc.fluka_name() {
            dropped.push(Dropped {
                id: *id,
                mass,
                reason: format!("no-fluka-name: {e}"),
            });
            continue;
        }
        if nuc.atomic_mass().is_none() {
            dropped.push(Dropped {
                id: *id,
                mass,
                reason: "no-atomic-mass".to_string(),
            });
            continue;
        }
        components.push(Component::new(*id, *w));
        accounted.push((*id, mass));
    }
    let text = compound_str(
        opts.fluka_fid,
        &opts.name,
        density,
        FracType::Mass,
        &components,
    )
    .map_err(Error::Fluka)?;
    Ok(Emitted {
        code: Code::Fluka,
        text,
        accounted,
        dropped,
        reparsed: false,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use nucleide_nuclei::NuclideId;

    #[test]
    fn fluka_compound_shape() {
        let mut mat = Material::new();
        mat.add_nuclide(NuclideId::from_name("U235").unwrap(), 4.0);
        mat.add_nuclide(NuclideId::from_name("U238").unwrap(), 96.0);
        let opts = EmitOptions::new("URANIUM").with_density(19.1);
        let out = emit_fluka(&mat, &opts).unwrap();
        assert!(out.text.contains("COMPOUND"));
        assert!(out.text.contains("235-U"));
        assert!(out.text.contains("238-U"));
        assert!(out.dropped.is_empty());
        assert!((out.mass_out() - 100.0).abs() < 1e-12);
    }

    #[test]
    fn fluka_drops_nameless_nuclide() {
        // O16 has no isotope entry in the vendored FLUKA table (only
        // natural OXYGEN, which never matches a NuclideId).
        let mut mat = Material::new();
        mat.add_nuclide(NuclideId::from_name("U235").unwrap(), 1.0);
        mat.add_nuclide(NuclideId::from_name("O16").unwrap(), 2.0);
        let opts = EmitOptions::new("UOX").with_density(10.0);
        let out = emit_fluka(&mat, &opts).unwrap();
        assert_eq!(out.dropped.len(), 1);
        assert_eq!(out.dropped[0].id, NuclideId::from_name("O16").unwrap());
        assert!(out.dropped[0].reason.starts_with("no-fluka-name"));
        assert!((out.mass_out() - 1.0).abs() < 1e-12);
    }
}
