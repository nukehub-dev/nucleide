//! ARMI nuclide-dialect bridge: database names, MCNP ZAIDs-as-strings,
//! MC2-3-style labels, and AAAZZZS ids.
//!
//! This module translates the label conventions used by
//! `armi.nucDirectory.nuclideBases` into [`NuclideId`] without vendoring
//! ARMI's tables:
//!
//! - ARMI database names: `nU235`, `nPu239`, `nAm242m` (a leading `n` plus
//!   the capitalized ARMI name). Bare ARMI/Capitalized names (`U235`,
//!   `AM242M`) are accepted as well.
//! - MCNP ZAIDs-as-strings: `92235`, `95242` (Am-242m), `95642` (Am-242
//!   ground) via [`dialects::from_zaid`].
//! - AAAZZZS ids: `2350920` (U-235), numerically identical to the Cinder
//!   encoding, via [`dialects::from_cinder`].
//! - MC2-3-style labels, only where unambiguous:
//!   - ENDF/B-V.2 dash form: `U-2355`, `B-10 5` (symbol, dash, mass digits,
//!     single `5`/`7` library tag).
//!   - ENDF/B-VII.x underscore form: `U235_7`, `FE56_7`, `H1___7`,
//!     `B10__7`.
//!   - Compact 6-character form: `PU2397`, `AM2427`, `CM2445`, `U2355`
//!     (letters plus mass digits plus a single `5`/`7` tag).
//!   - Full-mass trailing-`M` V.2 ids (`AM242M`, `TE129M`); these coincide
//!     with the ARMI isomer names and resolve identically.
//!
//! Out of scope (rejected as ambiguous): truncated isomer labels whose mass
//! digits were shortened to fit the 6-character ISOTXS field (`AM42M7`,
//! `TE29M7`); V.2 labels with letter version tags (`LI-7 V`, `PU239V`);
//! natural-element rows (`URANIUM`, `C____7`); lumped/dummy nuclides.
//!
//! Precedence inside [`armi_name_to_nucid`]: all-digit strings read as MCNP
//! ZAID first, then AAAZZZS (so `10010` is Ne-10-as-ZAID, not H-1-as-AAAZZZS).
//! Otherwise an MC2-3 interpretation wins over [`NuclideId::from_name`], so a
//! trailing `_7` is always the ENDF/B-VII library tag, never metastable
//! state 7 (write state 7 GNDS-style as `U235_m7`). A leading `n` is only
//! stripped when the direct parse fails, so nitrogen (`N15`) never loses
//! its symbol.

use crate::dialects;
use crate::{element_symbol, Error, NuclideId};

/// Library tags accepted on MC2-3-style labels (`5` = ENDF/B-V.2,
/// `7` = ENDF/B-VII.x).
const LIB_TAGS: [char; 2] = ['5', '7'];

/// Parse any ARMI-side label into a [`NuclideId`].
///
/// Accepts ARMI database names (`nU235`), bare ARMI/GNDS names (`U235`,
/// `AM242M`), MCNP ZAIDs-as-strings (`92235`), MC2-3-style labels
/// (`U-2355`, `U235_7`), and AAAZZZS ids (`2350920`). See the module docs
/// for precedence and the unambiguous MC2-3 subset.
pub fn armi_name_to_nucid(s: &str) -> Result<NuclideId, Error> {
    let t = s.trim();
    if t.is_empty() {
        return Err(Error::MissingMassNumber(s.to_string()));
    }
    if t.bytes().all(|b| b.is_ascii_digit()) {
        return digit_forms(t);
    }
    if let Ok(id) = mcc3_to_nucid(t) {
        return Ok(id);
    }
    if let Ok(id) = NuclideId::from_name(t) {
        return Ok(id);
    }
    if let Some(stripped) = strip_db_prefix(t) {
        if let Ok(id) = mcc3_to_nucid(stripped) {
            return Ok(id);
        }
        return NuclideId::from_name(stripped);
    }
    NuclideId::from_name(t)
}

/// Emit the ARMI database name for `id` (`U235` → `nU235`,
/// Am-242m → `nAm242m`), matching ARMI's `f"n{name.capitalize()}"`.
///
/// Total over every raw id (one invalid-id policy shared with
/// [`NuclideId::to_name`](crate::NuclideId::to_name)): validated ids keep the
/// exact historical spelling; raw ids outside the validated `(Z, A, state)`
/// domain render the `to_name` diagnostic fallback with the same `n` prefix
/// and capitalization, so label pipelines (`emit`, `r2s` snapshot,
/// `material` checks) can never panic on an unchecked integer.
pub fn nucid_to_armi_label(id: NuclideId) -> String {
    let name = match (id.is_valid(), element_symbol(id.z())) {
        // `is_valid` pins `1 <= Z <= 118`, so the symbol lookup succeeds;
        // the pair match keeps this total even if the table ever shrinks.
        (true, Some(sym)) => match id.state() {
            0 => format!("{}{}", sym.to_ascii_uppercase(), id.a()),
            1 => format!("{}{}M", sym.to_ascii_uppercase(), id.a()),
            st => format!("{}{}M{st}", sym.to_ascii_uppercase(), id.a()),
        },
        _ => id.to_name(),
    };
    let mut chars = name.chars();
    let mut out = String::with_capacity(name.len() + 1);
    out.push('n');
    if let Some(first) = chars.next() {
        out.extend(first.to_uppercase());
    }
    out.extend(chars.flat_map(|c| c.to_lowercase()));
    out
}

/// Parse an MC2-3-style label into a [`NuclideId`].
///
/// Accepts the dash (`U-2355`), underscore (`U235_7`), and compact
/// (`PU2397`) tagged forms plus full-mass trailing-`M` V.2 ids (`AM242M`);
/// every result is a ground state except the trailing-`M` forms (state 1).
/// Truncated isomer labels (`AM42M7`), letter-tagged V.2 labels (`LI-7 V`),
/// and anything else outside the unambiguous subset fail.
pub fn mcc3_to_nucid(s: &str) -> Result<NuclideId, Error> {
    let t = s.trim();
    if t.is_empty() {
        return Err(Error::MissingMassNumber(s.to_string()));
    }
    let u = t.to_ascii_uppercase();
    if let Some((sym, tail)) = u.split_once('-') {
        return mcc3_dash(s, sym.trim(), tail);
    }
    if u.contains('_') {
        return mcc3_underscore(s, &u);
    }
    if let Some(result) = mcc3_compact(&u) {
        return result;
    }
    if u.ends_with('M') {
        return NuclideId::from_name(&u);
    }
    Err(Error::BadNumber(s.to_string()))
}

/// All-digit input: MCNP ZAID first, then AAAZZZS (Cinder encoding).
fn digit_forms(t: &str) -> Result<NuclideId, Error> {
    let v: u32 = t.parse().map_err(|_| Error::BadNumber(t.to_string()))?;
    if let Ok(id) = dialects::from_zaid(v) {
        return Ok(id);
    }
    if let Ok(id) = dialects::from_cinder(v) {
        return Ok(id);
    }
    Err(Error::BadNumber(t.to_string()))
}

/// Strip one leading ARMI database `n`/`N` when followed by a letter
/// (`nU235` → `U235`). The caller tries the direct parse first so real
/// `N`-element names (`N15`) are never mangled.
fn strip_db_prefix(t: &str) -> Option<&str> {
    let mut chars = t.chars();
    let first = chars.next()?;
    let second = chars.next()?;
    if (first == 'n' || first == 'N') && second.is_ascii_alphabetic() {
        // `first` is ASCII, so byte index 1 is a char boundary.
        Some(&t[1..])
    } else {
        None
    }
}

/// Atomic number for a unilateral (uppercased) element symbol.
fn z_of_upper(sym: &str) -> Option<u32> {
    if sym.is_empty() {
        return None;
    }
    let mut chars = sym.chars();
    let mut canonical = String::with_capacity(sym.len());
    if let Some(first) = chars.next() {
        canonical.extend(first.to_uppercase());
    }
    canonical.extend(chars.flat_map(|c| c.to_lowercase()));
    crate::element_z(&canonical)
}

/// V.2 dash form: `U-2355`, `B-10 5`, `H-2  5`.
fn mcc3_dash(orig: &str, sym: &str, tail: &str) -> Result<NuclideId, Error> {
    if tail.chars().any(|c| c.is_ascii_alphabetic()) {
        return Err(Error::BadNumber(orig.to_string()));
    }
    let digits: String = tail.chars().filter(|c| c.is_ascii_digit()).collect();
    if digits.is_empty() {
        return Err(Error::MissingMassNumber(orig.to_string()));
    }
    if digits.len() < 2 || digits.len() > 4 {
        return Err(Error::BadNumber(orig.to_string()));
    }
    let (mass_str, tag) = digits.split_at(digits.len() - 1);
    let Some(tag_char) = tag.chars().next() else {
        return Err(Error::BadNumber(orig.to_string()));
    };
    if !LIB_TAGS.contains(&tag_char) {
        return Err(Error::BadNumber(orig.to_string()));
    }
    let z = z_of_upper(sym).ok_or_else(|| Error::UnknownElement(sym.to_string()))?;
    let a: u32 = mass_str
        .parse()
        .map_err(|_| Error::BadNumber(mass_str.to_string()))?;
    NuclideId::new(z, a, 0)
}

/// VII.x underscore form: `U235_7`, `FE56_7`, `H1___7`, `B10__7`.
fn mcc3_underscore(orig: &str, u: &str) -> Result<NuclideId, Error> {
    if u.chars().any(|c| c == 'M') {
        // Truncated isomers (`U23M_7`) and GNDS `m` forms (`AM242_M1`)
        // are not underscore MC2-3 labels.
        return Err(Error::BadNumber(orig.to_string()));
    }
    let bytes = u.as_bytes();
    let mut i = 0;
    while i < bytes.len() && bytes[i].is_ascii_alphabetic() {
        i += 1;
    }
    let sym = &u[..i];
    let mut j = i;
    while j < bytes.len() && bytes[j].is_ascii_digit() {
        j += 1;
    }
    let mass_str = &u[i..j];
    let rest = &u[j..];
    let mut rest_chars = rest.chars();
    let tag = rest_chars.next_back();
    let ok = (1..=2).contains(&sym.len())
        && !mass_str.is_empty()
        && mass_str.len() <= 3
        && rest.len() >= 2
        && rest_chars.all(|c| c == '_')
        && tag.is_some_and(|c| LIB_TAGS.contains(&c));
    if !ok {
        return Err(Error::BadNumber(orig.to_string()));
    }
    let z = z_of_upper(sym).ok_or_else(|| Error::UnknownElement(sym.to_string()))?;
    let a: u32 = mass_str
        .parse()
        .map_err(|_| Error::BadNumber(mass_str.to_string()))?;
    NuclideId::new(z, a, 0)
}

/// Compact tagged form: `PU2397`, `AM2427`, `CM2445`, `U2355`.
///
/// Returns `None` when the shape does not match so the caller can fall
/// through to the trailing-`M` rule and the final error.
fn mcc3_compact(u: &str) -> Option<Result<NuclideId, Error>> {
    let split = u.find(|c: char| !c.is_ascii_alphabetic())?;
    let (sym, digits) = u.split_at(split);
    if !(1..=2).contains(&sym.len()) || digits.len() != 4 {
        return None;
    }
    if !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let (mass_str, tag) = digits.split_at(3);
    let tag_char = tag.chars().next()?;
    if !LIB_TAGS.contains(&tag_char) {
        return None;
    }
    Some((|| {
        let z = z_of_upper(sym).ok_or_else(|| Error::UnknownElement(sym.to_string()))?;
        let a: u32 = mass_str
            .parse()
            .map_err(|_| Error::BadNumber(mass_str.to_string()))?;
        NuclideId::new(z, a, 0)
    })())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn nid(z: u32, a: u32, s: u32) -> NuclideId {
        NuclideId::new(z, a, s).unwrap()
    }

    #[test]
    fn parses_db_and_bare_armi_names() {
        assert_eq!(armi_name_to_nucid("nU235").unwrap(), nid(92, 235, 0));
        assert_eq!(armi_name_to_nucid("U235").unwrap(), nid(92, 235, 0));
        assert_eq!(armi_name_to_nucid("NU235").unwrap(), nid(92, 235, 0));
        assert_eq!(armi_name_to_nucid("nU238").unwrap(), nid(92, 238, 0));
        assert_eq!(armi_name_to_nucid("nPu239").unwrap(), nid(94, 239, 0));
        assert_eq!(armi_name_to_nucid("PU239").unwrap(), nid(94, 239, 0));
    }

    #[test]
    fn parses_armi_isomer_names() {
        assert_eq!(armi_name_to_nucid("AM242M").unwrap(), nid(95, 242, 1));
        assert_eq!(armi_name_to_nucid("nAm242m").unwrap(), nid(95, 242, 1));
        assert_eq!(armi_name_to_nucid("AM242M2").unwrap(), nid(95, 242, 2));
    }

    #[test]
    fn db_prefix_never_eats_nitrogen() {
        assert_eq!(armi_name_to_nucid("N15").unwrap(), nid(7, 15, 0));
        assert_eq!(armi_name_to_nucid("nN15").unwrap(), nid(7, 15, 0));
    }

    #[test]
    fn parses_zaid_strings() {
        assert_eq!(armi_name_to_nucid("92235").unwrap(), nid(92, 235, 0));
        assert_eq!(armi_name_to_nucid("94239").unwrap(), nid(94, 239, 0));
        assert_eq!(armi_name_to_nucid("95242").unwrap(), nid(95, 242, 1));
        assert_eq!(armi_name_to_nucid("95642").unwrap(), nid(95, 242, 0));
        assert_eq!(armi_name_to_nucid("1001").unwrap(), nid(1, 1, 0));
    }

    #[test]
    fn parses_aaazzzs_ids() {
        assert_eq!(armi_name_to_nucid("2350920").unwrap(), nid(92, 235, 0));
        assert_eq!(armi_name_to_nucid("2390940").unwrap(), nid(94, 239, 0));
        assert_eq!(armi_name_to_nucid("2420951").unwrap(), nid(95, 242, 1));
    }

    #[test]
    fn parses_mcc3_forms() {
        assert_eq!(mcc3_to_nucid("U-2355").unwrap(), nid(92, 235, 0));
        assert_eq!(mcc3_to_nucid("U-2365").unwrap(), nid(92, 236, 0));
        assert_eq!(mcc3_to_nucid("B-10 5").unwrap(), nid(5, 10, 0));
        assert_eq!(mcc3_to_nucid("H-2  5").unwrap(), nid(1, 2, 0));
        assert_eq!(mcc3_to_nucid("U235_7").unwrap(), nid(92, 235, 0));
        assert_eq!(mcc3_to_nucid("FE56_7").unwrap(), nid(26, 56, 0));
        assert_eq!(mcc3_to_nucid("H1___7").unwrap(), nid(1, 1, 0));
        assert_eq!(mcc3_to_nucid("B10__7").unwrap(), nid(5, 10, 0));
        assert_eq!(mcc3_to_nucid("LI6__7").unwrap(), nid(3, 6, 0));
        assert_eq!(mcc3_to_nucid("PU2397").unwrap(), nid(94, 239, 0));
        assert_eq!(mcc3_to_nucid("AM2427").unwrap(), nid(95, 242, 0));
        assert_eq!(mcc3_to_nucid("CM2445").unwrap(), nid(96, 244, 0));
        assert_eq!(mcc3_to_nucid("U2355").unwrap(), nid(92, 235, 0));
        assert_eq!(mcc3_to_nucid("AM242M").unwrap(), nid(95, 242, 1));
        assert_eq!(mcc3_to_nucid("TE129M").unwrap(), nid(52, 129, 1));
        // Lowercase MCC3 input is accepted.
        assert_eq!(mcc3_to_nucid("u235_7").unwrap(), nid(92, 235, 0));
        assert_eq!(mcc3_to_nucid("u-2355").unwrap(), nid(92, 235, 0));
        // And through the full bridge.
        assert_eq!(armi_name_to_nucid("U-2355").unwrap(), nid(92, 235, 0));
        assert_eq!(armi_name_to_nucid("U235_7").unwrap(), nid(92, 235, 0));
        assert_eq!(armi_name_to_nucid("nU235_7").unwrap(), nid(92, 235, 0));
    }

    #[test]
    fn emits_armi_db_names() {
        assert_eq!(nucid_to_armi_label(nid(92, 235, 0)), "nU235");
        assert_eq!(nucid_to_armi_label(nid(92, 238, 0)), "nU238");
        assert_eq!(nucid_to_armi_label(nid(94, 239, 0)), "nPu239");
        assert_eq!(nucid_to_armi_label(nid(95, 242, 1)), "nAm242m");
        assert_eq!(nucid_to_armi_label(nid(95, 242, 2)), "nAm242m2");
        assert_eq!(nucid_to_armi_label(nid(1, 1, 0)), "nH1");
    }

    #[test]
    fn invalid_raw_ids_get_fallback_labels() {
        // Formerly `.expect()` panics on the element-symbol lookup; now the
        // `to_name` diagnostic with the ARMI `n`-prefix treatment.
        assert_eq!(nucid_to_armi_label(NuclideId::from_nucid(0)), "nZ0a0[m0]");
        assert_eq!(
            nucid_to_armi_label(NuclideId::from_nucid(920_050_000)),
            "nZ92a5[m0]"
        );
        assert_eq!(
            nucid_to_armi_label(NuclideId::from_nucid(u32::MAX)),
            "nZ429a496[m5]"
        );
    }

    #[test]
    fn db_labels_round_trip() {
        for id in [
            nid(92, 235, 0),
            nid(92, 238, 0),
            nid(94, 239, 0),
            nid(95, 242, 0),
            nid(95, 242, 1),
            nid(1, 1, 0),
            nid(26, 56, 0),
        ] {
            assert_eq!(armi_name_to_nucid(&nucid_to_armi_label(id)).unwrap(), id);
        }
    }

    #[test]
    fn rejects_garbage() {
        for bad in [
            "",
            "   ",
            "n",
            "NOTANUCLIDE",
            "???",
            "AM42M7",
            "TE29M7",
            "LI-7 V",
            "PU239V",
            "U-235X5",
            "nXx999",
            "99999999",
        ] {
            assert!(armi_name_to_nucid(bad).is_err(), "{bad}");
        }
        // Truncated isomers and letter-tagged V.2 labels fail directly too.
        assert!(mcc3_to_nucid("AM42M7").is_err());
        assert!(mcc3_to_nucid("TE29M7").is_err());
        assert!(mcc3_to_nucid("LI-7 V").is_err());
        assert!(mcc3_to_nucid("U235_8").is_err());
        assert!(mcc3_to_nucid("U235_").is_err());
        assert!(mcc3_to_nucid("U235").is_err());
        assert!(mcc3_to_nucid("Am242_m1").is_err());
    }

    #[test]
    fn serpent_forms_still_resolve_through_bridge() {
        // `U-235` is Serpent, not MCC3: the MCC3 dash read (U-23, invalid)
        // fails and the GNDS fallback recovers U-235.
        assert_eq!(armi_name_to_nucid("U-235").unwrap(), nid(92, 235, 0));
        assert!(!mcc3_to_nucid("U-235").unwrap_err().to_string().is_empty());
    }
}
