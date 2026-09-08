//! Serpent `mat` card emission.
//!
//! Mass density is emitted negative and mass fractions are emitted negative,
//! per the Serpent convention (negative `dens` selects g/cm³, negative
//! fractions select mass fractions; Serpent normalises automatically).
//! Nuclides use `{zaid}.{lib}` identifiers (e.g. `92235.03c`), matching the
//! official examples — bare names do not resolve against an `acelib` table.
//! One nuclide per line; Serpent cards need no continuation markers. Serpent
//! has no material reader in this workspace, so drift is analytic
//! ([`Emitted::reparsed`] is always false).

use nucleide_material::Material;
use nucleide_nuclei::dialects::to_zaid;

use crate::{Code, EmitOptions, Emitted, Error, Result};

/// Emit `mat` as a Serpent `mat` card.
pub fn emit_serpent(mat: &Material, opts: &EmitOptions) -> Result<Emitted> {
    let density = opts.density_for(mat)?;
    let fracs = mat.weight_fractions().map_err(|_| Error::Degenerate)?;
    let mut accounted = Vec::with_capacity(fracs.len());
    let mut text = format!("mat {} {}\n", opts.name, -density);
    for (id, w) in &fracs {
        text.push_str(&format!("  {}.{} {}\n", to_zaid(*id), opts.serpent_lib, -w));
        accounted.push((*id, mat.comp.get(id).copied().unwrap_or(0.0)));
    }
    Ok(Emitted {
        code: Code::Serpent,
        text,
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
    fn serpent_card_shape() {
        let mut mat = Material::new();
        mat.add_nuclide(NuclideId::from_name("U235").unwrap(), 5.0);
        mat.add_nuclide(NuclideId::from_name("U238").unwrap(), 95.0);
        let opts = EmitOptions::new("leu").with_density(10.0);
        let out = emit_serpent(&mat, &opts).unwrap();
        assert!(out.text.starts_with("mat leu -10\n"));
        assert!(out.text.contains("92235.03c -0.05"));
        assert!(out.text.contains("92238.03c -0.95"));
        assert!(!out.reparsed);
        assert!((out.mass_out() - 100.0).abs() < 1e-12);
    }

    #[test]
    fn serpent_needs_density() {
        let mut mat = Material::new();
        mat.add_nuclide(NuclideId::from_name("U235").unwrap(), 1.0);
        let opts = EmitOptions::new("leu");
        assert!(matches!(
            emit_serpent(&mat, &opts),
            Err(Error::MissingDensity)
        ));
    }
}
