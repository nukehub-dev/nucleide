//! MCNP `m<number>` card emission.
//!
//! Mass fractions are emitted negative ([`FracKind::Mass`] convention) with
//! `{zaid}.{xs_suffix}` identifiers from [`to_zaid`]. Four pairs per line;
//! continuation lines start with five blanks. The text is verified by
//! re-parsing with [`materials_from_str`] and comparing fractions.

use nucleide_material::Material;
use nucleide_mcnp_io::inp::{materials_from_str, FracKind};
use nucleide_nuclei::{dialects::to_zaid, NuclideId};

use crate::{Code, EmitOptions, Emitted, Error, Result};

/// Pairs of `(zaid, fraction)` per card line.
const PAIRS_PER_LINE: usize = 4;

/// Emit `mat` as an MCNP material card.
pub fn emit_mcnp(mat: &Material, opts: &EmitOptions) -> Result<Emitted> {
    let fracs = mat.weight_fractions().map_err(|_| Error::Degenerate)?;
    let mut accounted = Vec::with_capacity(fracs.len());
    let mut pairs = Vec::with_capacity(fracs.len());
    for (id, w) in &fracs {
        // Natural-element placeholders (a == 0) emit as Z000, which MCNP
        // accepts; the reader keeps them as placeholders.
        pairs.push((to_zaid(*id), -*w));
        accounted.push((*id, mat.comp.get(id).copied().unwrap_or(0.0)));
    }

    let mut text = format!("m{}", opts.mcnp_number);
    for (i, chunk) in pairs.chunks(PAIRS_PER_LINE).enumerate() {
        if i > 0 {
            text.push_str("\n     ");
        }
        for (zaid, frac) in chunk {
            text.push_str(&format!(" {zaid}.{} {frac}", opts.xs_suffix));
        }
    }
    text.push('\n');

    // Re-parse verification: every emitted nuclide must come back with its
    // fraction within tolerance.
    let reparsed = verify_round_trip(&text, opts, &fracs)?;

    Ok(Emitted {
        code: Code::Mcnp,
        text,
        accounted,
        dropped: Vec::new(),
        reparsed,
    })
}

fn verify_round_trip(
    text: &str,
    opts: &EmitOptions,
    fracs: &std::collections::BTreeMap<NuclideId, f64>,
) -> Result<bool> {
    let mats = materials_from_str(text).map_err(|e| Error::Reparse {
        code: Code::Mcnp,
        detail: e.to_string(),
    })?;
    let back = mats
        .iter()
        .find(|m| m.number == opts.mcnp_number)
        .ok_or_else(|| Error::Reparse {
            code: Code::Mcnp,
            detail: format!("m{} missing after re-parse", opts.mcnp_number),
        })?;
    if back.fraction_type != FracKind::Mass {
        return Ok(false);
    }
    // Recovered fractions are positive magnitudes; compare against |emitted|.
    let mut recovered: std::collections::BTreeMap<NuclideId, f64> =
        std::collections::BTreeMap::new();
    for (id, f) in &back.fractions {
        *recovered.entry(*id).or_insert(0.0) += f.abs();
    }
    for (id, w) in fracs {
        // Natural-element placeholders round-trip by Z only.
        let got = if id.a() == 0 {
            recovered
                .iter()
                .filter(|(k, _)| k.z() == id.z())
                .map(|(_, v)| v)
                .sum()
        } else {
            recovered.get(id).copied().unwrap_or(0.0)
        };
        if (got - w).abs() > 1e-9 * w.abs().max(1e-300) {
            return Ok(false);
        }
    }
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;

    pub(crate) fn uo2() -> Material {
        let mut mat = Material::new();
        mat.add_nuclide(NuclideId::from_name("U235").unwrap(), 5.0);
        mat.add_nuclide(NuclideId::from_name("U238").unwrap(), 95.0);
        mat.add_nuclide(NuclideId::from_name("O16").unwrap(), 100.0);
        mat
    }

    #[test]
    fn mcnp_card_shape_and_round_trip() {
        let opts = EmitOptions::new("uo2");
        let out = emit_mcnp(&uo2(), &opts).unwrap();
        assert!(out.text.starts_with("m1 "));
        assert!(out.text.contains("92235.80c -0.025"));
        assert!(out.text.contains("8016.80c -0.5"));
        assert!(out.reparsed);
        assert!(out.dropped.is_empty());
        assert!((out.mass_out() - 200.0).abs() < 1e-12);
    }

    #[test]
    fn mcnp_empty_is_degenerate() {
        let opts = EmitOptions::new("void");
        assert!(matches!(
            emit_mcnp(&Material::new(), &opts),
            Err(Error::Degenerate)
        ));
    }
}
