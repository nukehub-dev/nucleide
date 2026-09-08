//! ALARA `mixture` block emission.
//!
//! Each nuclide becomes an `element <name> 1.0 <massfrac>` entry using the
//! canonical [`NuclideId`] name (which the deck reader parses back). The
//! block text matches [`AlaraDeck`] display formatting and is verified by
//! re-parsing the fragment with [`AlaraDeck::parse`].

use nucleide_alara_io::AlaraDeck;
use nucleide_material::Material;

use crate::{Code, EmitOptions, Emitted, Error, Result};

/// Emit `mat` as an ALARA `mixture` block.
pub fn emit_alara(mat: &Material, opts: &EmitOptions) -> Result<Emitted> {
    let fracs = mat.weight_fractions().map_err(|_| Error::Degenerate)?;
    let mut accounted = Vec::with_capacity(fracs.len());
    let mut text = format!("mixture {}\n", opts.name);
    for (id, w) in &fracs {
        text.push_str(&format!("element {} 1 {w}\n", id.to_name()));
        accounted.push((*id, mat.comp.get(id).copied().unwrap_or(0.0)));
    }
    text.push_str("end\n");

    let reparsed = verify_round_trip(&text, opts, &fracs)?;

    Ok(Emitted {
        code: Code::Alara,
        text,
        accounted,
        dropped: Vec::new(),
        reparsed,
    })
}

fn verify_round_trip(
    text: &str,
    opts: &EmitOptions,
    fracs: &std::collections::BTreeMap<nucleide_nuclei::NuclideId, f64>,
) -> Result<bool> {
    let deck = AlaraDeck::parse(text).map_err(|e| Error::Reparse {
        code: Code::Alara,
        detail: e.to_string(),
    })?;
    let mix = deck
        .find_mixture(&opts.name)
        .ok_or_else(|| Error::Reparse {
            code: Code::Alara,
            detail: format!("mixture {} missing after re-parse", opts.name),
        })?;
    if mix.entries.len() != fracs.len() {
        return Ok(false);
    }
    for entry in &mix.entries {
        let (symbol, vol) = match entry {
            nucleide_alara_io::deck::MixtureEntry::Element {
                symbol,
                vol_fraction,
                ..
            } => (symbol, *vol_fraction),
            _ => return Ok(false),
        };
        let id = nucleide_nuclei::NuclideId::from_name(symbol).map_err(|e| Error::Reparse {
            code: Code::Alara,
            detail: e.to_string(),
        })?;
        let want = fracs.get(&id).copied().unwrap_or(0.0);
        if (vol - want).abs() > 1e-9 * want.abs().max(1e-300) {
            return Ok(false);
        }
    }
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nucleide_nuclei::NuclideId;

    #[test]
    fn alara_mixture_shape_and_round_trip() {
        let mut mat = Material::new();
        mat.add_nuclide(NuclideId::from_name("Fe56").unwrap(), 3.0);
        mat.add_nuclide(NuclideId::from_name("C12").unwrap(), 1.0);
        let opts = EmitOptions::new("steelish");
        let out = emit_alara(&mat, &opts).unwrap();
        assert!(out.text.starts_with("mixture steelish\n"));
        assert!(out.text.contains("element Fe56 1 0.75"));
        assert!(out.text.contains("element C12 1 0.25"));
        assert!(out.text.ends_with("end\n"));
        assert!(out.reparsed);
        assert!((out.mass_out() - 4.0).abs() < 1e-12);
    }
}
