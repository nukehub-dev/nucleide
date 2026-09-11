//! ARMI blueprint mass fractions → [`Material`].
//!
//! One-way bridge from ARMI-side `material.massFrac` dicts (post-expansion
//! nuclide keys) into the emission pipeline. Keys resolve through the
//! already-landed [`armi_name_to_nucid`]
//! bridge, so every documented precedence there applies here unchanged:
//! `nU235` database names, bare `U235`/`AM242M`, MCNP ZAIDs-as-strings
//! (`92235`), AAAZZZS ids (`2350920`), and unambiguous MC2-3 labels
//! (`U-2355`, `U235_7`, `PU2397`). Truncated isomers, letter-tagged V.2
//! labels, naturals, and lumped/dummy nuclides (`DUMP1`, `DUMP2`, `LFP*`)
//! fail in the bridge and surface here as [`Error::Nuclei`] naming the key.
//!
//! # v1 caller-side rules
//!
//! The emitter never reimplements ARMI-side resolution; the caller does:
//!
//! - **Elemental keys rejected.** `ZR`, `FE`, `NA`, `O`, `C`, `W`, … (any
//!   case, with or without the `n` database prefix) fail with an "expand
//!   first" error. The caller passes ARMI's post-expansion
//!   `material.massFrac`, never raw elemental fractions.
//! - **Number fractions / number densities out of scope.** Values are masses
//!   in grams; the caller converts first. Any positive masses are accepted
//!   (emission normalizes via weight fractions); negative or non-finite
//!   values fail naming the key.
//! - **Enrichment shorthands, `balance`, temperatures out of scope.**
//!   `U235_wt_frac`/`TD_frac`, the `balance` keyword, and `Tinput`/`Thot`
//!   are resolved by the caller, which passes the hot density in g/cm³.
//! - **Bare `AM242` rejected.** ARMI uses the bare alias for the m-state
//!   while GNDS parsing would silently yield ground, so guessing either way
//!   risks a silent isomer flip. Pass `AM242M` for Am-242m explicitly; the
//!   ground state is `AM242G` (accepted here as an explicit alias for
//!   Am-242g, alongside the GNDS `Am242_m0` and ZAID `95642` forms).
//!
//! Density is set exactly like the Python `emit_drift_inner` helper
//! (`set_density(density)`, unconditionally): `None` leaves the material
//! without a density, so Serpent/FLUKA/PARTISN emission still fails with
//! [`Error::MissingDensity`] while MCNP/ALARA are unaffected.

use nucleide_material::Material;
use nucleide_nuclei::{armi::armi_name_to_nucid, NuclideId};

use crate::{Error, Result};

/// Build a [`Material`] from ARMI-side mass fractions.
///
/// `pairs` maps ARMI nuclide keys to masses in grams (any positive values;
/// emission normalizes); `density` is the hot mass density in g/cm³, set
/// unconditionally like the Python `emit_drift_inner` helper. See the module
/// docs for the v1 caller-side rules (elemental keys, `AM242`, units).
pub fn from_armi_mass_fracs<I, S>(pairs: I, density: Option<f64>) -> Result<Material>
where
    I: IntoIterator<Item = (S, f64)>,
    S: AsRef<str>,
{
    let mut mat = Material::new();
    for (key, mass) in pairs {
        let id = armi_key_to_nucid(key.as_ref())?;
        if !mass.is_finite() || mass < 0.0 {
            return Err(Error::ArmiKey {
                key: key.as_ref().trim().to_string(),
                reason: format!("mass {mass} is not a non-negative finite value"),
            });
        }
        mat.add_nuclide(id, mass);
    }
    mat.set_density(density);
    Ok(mat)
}

/// Resolve one ARMI-side key: v1 rejections first, then the nuclei bridge.
fn armi_key_to_nucid(key: &str) -> Result<NuclideId> {
    let t = key.trim();
    if let Some(sym) = elemental_symbol(t).or_else(|| strip_db_prefix(t).and_then(elemental_symbol))
    {
        return Err(Error::ArmiKey {
            key: t.to_string(),
            reason: format!(
                "elemental key `{sym}`: pass ARMI post-expansion nuclide mass fractions \
                 (material.massFrac); the emitter never reimplements \
                 expandElementalMassFracsToNuclides — expand first"
            ),
        });
    }
    match bare_core(t).as_deref() {
        Some("AM242") => {
            return Err(Error::ArmiKey {
                key: t.to_string(),
                reason: "ambiguous bare `AM242` (ARMI means the m-state): pass `AM242M` for \
                         Am-242m or `AM242G` for ground explicitly"
                    .to_string(),
            });
        }
        // Explicit ground alias (the bridge has no `G` suffix form).
        Some("AM242G") => return Ok(NuclideId::new(95, 242, 0)?),
        _ => {}
    }
    // Name the key on bridge failures (transparent `Error::Nuclei`).
    armi_name_to_nucid(t)
        .map_err(|e| Error::Nuclei(nucleide_nuclei::Error::BadNumber(format!("{t}: {e}"))))
}

/// Canonical element symbol when `t` is a bare elemental key (`ZR` → `Zr`),
/// else `None`. Only pure-letter input qualifies, so no nuclide key with
/// mass digits, tags, or separators can collide.
fn elemental_symbol(t: &str) -> Option<String> {
    if t.is_empty() || !t.bytes().all(|b| b.is_ascii_alphabetic()) {
        return None;
    }
    let mut chars = t.chars();
    let mut canon = String::with_capacity(t.len());
    if let Some(first) = chars.next() {
        canon.extend(first.to_uppercase());
    }
    canon.extend(chars.flat_map(|c| c.to_lowercase()));
    if nucleide_nuclei::element_z(&canon).is_some() {
        Some(canon)
    } else {
        None
    }
}

/// Strip one leading ARMI database `n`/`N` when followed by a letter
/// (`nZr` → `Zr`), mirroring the bridge so prefixed elementals are caught.
fn strip_db_prefix(t: &str) -> Option<&str> {
    let mut chars = t.chars();
    let first = chars.next()?;
    let second = chars.next()?;
    if (first == 'n' || first == 'N') && second.is_ascii_alphabetic() {
        Some(&t[1..])
    } else {
        None
    }
}

/// Uppercased key with the database prefix stripped (`nAm242` → `AM242`),
/// for the bare-`AM242` disambiguation check.
fn bare_core(t: &str) -> Option<String> {
    if t.is_empty() {
        return None;
    }
    let u = t.to_ascii_uppercase();
    let core = match u.strip_prefix('N') {
        Some(rest) if rest.starts_with(|c: char| c.is_ascii_alphabetic()) => rest,
        _ => u.as_str(),
    };
    Some(core.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{emit_all, Code, EmitOptions};

    fn gnds_metal() -> Material {
        let mut mat = Material::new();
        mat.add_nuclide(NuclideId::from_name("U235").unwrap(), 5.0);
        mat.add_nuclide(NuclideId::from_name("U238").unwrap(), 95.0);
        mat.set_density(Some(19.1));
        mat
    }

    #[test]
    fn armi_matches_gnds_cards_with_goldens() {
        // ARMI database names emit byte-identical cards to the GNDS baseline.
        let armi = from_armi_mass_fracs([("nU235", 5.0), ("nU238", 95.0)], Some(19.1)).unwrap();
        let gnds = gnds_metal();
        let opts = EmitOptions::new("umetal");
        let armi_out = emit_all(&armi, &opts).unwrap();
        let gnds_out = emit_all(&gnds, &opts).unwrap();
        assert_eq!(armi_out.len(), 5);
        for (a, g) in armi_out.iter().zip(gnds_out.iter()) {
            assert_eq!(a.code, g.code);
            assert_eq!(a.text, g.text, "{:?}", a.code);
        }
        let mcnp = &armi_out[Code::Mcnp as usize].text;
        assert_eq!(mcnp, "m1 92235.80c -0.05 92238.80c -0.95\n");
        let serpent = &armi_out[Code::Serpent as usize].text;
        assert_eq!(
            serpent,
            "mat umetal -19.1\n  92235.03c -0.05\n  92238.03c -0.95\n"
        );
    }

    #[test]
    fn isomer_mcc3_and_zaid_keys_resolve() {
        let mat = from_armi_mass_fracs(
            [
                ("nAm242m", 2.0),
                ("U-2355", 4.0),
                ("U235_7", 1.0),
                ("PU2397", 3.0),
                ("92238", 5.0),
                ("2350920", 1.0),
                ("1001", 2.0),
            ],
            Some(10.0),
        )
        .unwrap();
        let get = |name: &str| {
            mat.comp
                .get(&NuclideId::from_name(name).unwrap())
                .copied()
                .unwrap_or(0.0)
        };
        assert!((get("Am242_m1") - 2.0).abs() < 1e-12);
        assert!((get("U235") - 6.0).abs() < 1e-12);
        assert!((get("Pu239") - 3.0).abs() < 1e-12);
        assert!((get("U238") - 5.0).abs() < 1e-12);
        assert!((get("H1") - 2.0).abs() < 1e-12);
        let opts = EmitOptions::new("mix");
        let out = emit_all(&mat, &opts).unwrap();
        let mcnp = &out[Code::Mcnp as usize].text;
        for zaid in [
            "95242.80c",
            "92235.80c",
            "94239.80c",
            "92238.80c",
            "1001.80c",
        ] {
            assert!(mcnp.contains(zaid), "{mcnp}");
        }
    }

    #[test]
    fn garbage_keys_fail_naming_the_key() {
        for bad in ["nXx999", "DUMP1", "DUMP2", "LFP38", "NOTANUCLIDE", ""] {
            let err = from_armi_mass_fracs([(bad, 1.0)], Some(1.0)).unwrap_err();
            assert!(matches!(err, Error::Nuclei(_)), "{bad}: {err}");
            assert!(err.to_string().contains(bad), "{bad}: {err}");
        }
    }

    #[test]
    fn bare_am242_rejected_but_explicit_states_accepted() {
        for bare in ["AM242", "am242", "nAm242", "NAm242"] {
            let err = from_armi_mass_fracs([(bare, 1.0)], Some(1.0)).unwrap_err();
            assert!(matches!(err, Error::ArmiKey { .. }), "{bare}: {err}");
            assert!(err.to_string().contains("AM242M"), "{bare}: {err}");
        }
        let m = from_armi_mass_fracs([("AM242M", 1.0)], None).unwrap();
        assert!(m
            .comp
            .contains_key(&NuclideId::from_name("Am242_m1").unwrap()));
        for g in ["AM242G", "nAm242G"] {
            let m = from_armi_mass_fracs([(g, 1.0)], None).unwrap();
            let id = NuclideId::new(95, 242, 0).unwrap();
            assert!(m.comp.contains_key(&id), "{g}");
            assert_eq!(m.comp.len(), 1);
        }
    }

    #[test]
    fn elemental_keys_rejected_with_expand_first() {
        for el in [
            "ZR", "FE", "NA", "O", "C", "W", "Zr", "zr", "nZr", "U", "XE",
        ] {
            let err = from_armi_mass_fracs([(el, 1.0)], Some(1.0)).unwrap_err();
            assert!(matches!(err, Error::ArmiKey { .. }), "{el}: {err}");
            assert!(err.to_string().contains("expand first"), "{el}: {err}");
            assert!(err.to_string().contains(el.trim()), "{el}: {err}");
        }
    }

    #[test]
    fn negative_and_non_finite_masses_rejected() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1.0] {
            let err = from_armi_mass_fracs([("nU235", bad)], Some(1.0)).unwrap_err();
            assert!(matches!(err, Error::ArmiKey { .. }), "{bad}: {err}");
            assert!(err.to_string().contains("nU235"), "{bad}: {err}");
        }
    }

    #[test]
    fn density_flows_to_material_like_emit_drift_inner() {
        let with = from_armi_mass_fracs([("nU235", 1.0)], Some(19.1)).unwrap();
        assert_eq!(with.density(), Some(19.1));
        let without = from_armi_mass_fracs([("nU235", 1.0)], None).unwrap();
        assert_eq!(without.density(), None);
        // None leaves the existing MissingDensity behavior for Serpent while
        // MCNP (density-free) still emits.
        let opts = EmitOptions::new("umetal");
        assert!(matches!(
            crate::serpent::emit_serpent(&without, &opts),
            Err(Error::MissingDensity)
        ));
        assert!(crate::mcnp::emit_mcnp(&without, &opts).unwrap().reparsed);
    }
}
