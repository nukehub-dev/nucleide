//! SDEF decay-source cards (E9) from caller-supplied decay lines.
//!
//! The upstream point-source emitter is monoenergetic: it renders
//! `SDEF POS=...`, an optional `VEC=... DIR=1` line, `ERG=<E>`, `WGT=<w>`,
//! and `PAR=<designator>`. This module generalizes it to a caller-supplied
//! line list: every energy/intensity pair is an input (the E8 caller-constants
//! posture applied to N lines), intensities are normalized to probabilities
//! (E9), and the card text follows the upstream field order with the
//! monoenergetic `ERG=<E>` replaced by a discrete-energy distribution
//! `ERG=D1` plus paired `SI1`/`SP1` cards when more than one line survives.
//!
//! # Verified surface syntax, unverified sampling semantics
//!
//! The emitted shapes are verified surface syntax: `SDEF key=value` fields,
//! 5-space continuations, `ERG=Dn` distribution references, and the
//! `SIn L` / `SPn D` card pair (one entry per line) all round-trip through a
//! third-party MCNP input parser (MontePy's parser test suite). What is
//! **not** verified from local files is the MCNP-side sampling semantics:
//! what the `L`/`D` flags select beyond "discrete list + per-entry values",
//! the defaults when `SPn` is omitted, zero/negative probability handling,
//! whether MCNP renormalizes, and any bins-vs-lines distinction for energy.
//! Those semantics are the caller's responsibility against the MCNP manual;
//! this renderer only guarantees the card shapes above and that the emitted
//! probabilities are normalized to sum to 1.0 (E9) before 6-significant-digit
//! card formatting. Re-checked against that parser suite: it exercises only
//! the surface shapes above, so the semantics stay unverified.

use nucleide_nuclei::particles::ParticleId;

use crate::Error;

/// MCNP fixed-format card width; longer entry lists wrap onto continuation
/// lines.
const CARD_WIDTH: usize = 80;

/// Caller-supplied point-source fields for an SDEF card, mirroring the
/// upstream monoenergetic `PointSource` constructor (position, direction,
/// weight, particle). The energy comes from the decay line list instead.
#[derive(Debug, Clone, PartialEq)]
pub struct PointSource {
    /// Source position x \[cm\].
    pub x: f64,
    /// Source position y \[cm\].
    pub y: f64,
    /// Source position z \[cm\].
    pub z: f64,
    /// Direction vector u component (emit `VEC`/`DIR=1` only when nonzero).
    pub u: f64,
    /// Direction vector v component.
    pub v: f64,
    /// Direction vector w component.
    pub w: f64,
    /// Particle weight for `WGT=`; must be finite and positive.
    pub weight: f64,
    /// Particle species for `PAR=` (designator from `nucleide-nuclei`).
    pub particle: ParticleId,
}

impl Default for PointSource {
    /// Upstream constructor defaults: origin, no direction, unit weight,
    /// neutron.
    fn default() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            u: 0.0,
            v: 0.0,
            w: 0.0,
            weight: 1.0,
            particle: ParticleId::Neutron,
        }
    }
}

/// Normalize caller-supplied decay lines into `(energy, probability)` bins
/// (E9): merge duplicate energies by summing intensities, sort ascending,
/// and divide by the total intensity.
///
/// Specified behavior (no upstream to inherit — upstream stops at the
/// monoenergetic card):
///
/// - empty input is an error ([`Error::EmptyLines`]);
/// - non-finite energies or intensities are an error
///   ([`Error::NonFiniteLine`]);
/// - negative energies ([`Error::NegativeEnergy`]) and negative intensities
///   ([`Error::NegativeIntensity`]) are errors;
/// - zero-intensity lines are dropped (an `SP1` entry with zero probability
///   has unverified MCNP semantics); if none survive, [`Error::NoEmission`];
/// - a line with positive intensity can still normalize to a `0` probability
///   by `f64` underflow when intensities span extreme ranges (e.g. `1e-308`
///   against `1e308`); such bins are kept and emitted as literal `0` entries
///   on the `SP1 D` card — only pre-normalization zero-intensity lines are
///   dropped;
/// - intensities that overflow to a non-finite total are
///   [`Error::NonFiniteLine`]; a total that underflows to zero is
///   [`Error::NoEmission`];
/// - duplicate energies merge by summing, then bins are sorted ascending, so
///   the result is independent of input order;
/// - the returned probabilities sum to 1.0 in exact arithmetic.
pub fn normalize_decay_lines(lines: &[(f64, f64)]) -> Result<Vec<(f64, f64)>, Error> {
    if lines.is_empty() {
        return Err(Error::EmptyLines);
    }
    for &(e, i) in lines {
        if !e.is_finite() || !i.is_finite() {
            return Err(Error::NonFiniteLine(if e.is_finite() { i } else { e }));
        }
        if e < 0.0 {
            return Err(Error::NegativeEnergy(e));
        }
        if i < 0.0 {
            return Err(Error::NegativeIntensity(i));
        }
    }
    let mut merged: Vec<(f64, f64)> = lines.iter().copied().filter(|&(_, i)| i > 0.0).collect();
    if merged.is_empty() {
        return Err(Error::NoEmission);
    }
    merged.sort_by(|a, b| a.0.partial_cmp(&b.0).expect("energies are finite"));
    let mut deduped: Vec<(f64, f64)> = Vec::with_capacity(merged.len());
    for (e, i) in merged {
        match deduped.last_mut() {
            Some(last) if last.0 == e => last.1 += i,
            _ => deduped.push((e, i)),
        }
    }
    let total: f64 = deduped.iter().map(|l| l.1).sum();
    if !total.is_finite() {
        return Err(Error::NonFiniteLine(total));
    }
    if total <= 0.0 {
        return Err(Error::NoEmission);
    }
    Ok(deduped.iter().map(|&(e, i)| (e, i / total)).collect())
}

/// Render a decay-source SDEF card plus, for multi-line sources, the paired
/// `SI1`/`SP1` distribution cards (E9).
///
/// Returns the normalized bins (as in [`normalize_decay_lines`]) and the card
/// text. `version` selects the `PAR=` designator dialect: 5 via
/// [`ParticleId::mcnp`], 6 via [`ParticleId::mcnp6`]; any other value is an
/// error, and a particle with no designator for that version is an error.
///
/// Field order follows the upstream monoenergetic emitter: `SDEF POS=` is
/// always first; `VEC=<u> <v> <w> DIR=1` joins on a 5-space continuation
/// line only when some direction component is nonzero (comparison-based, so
/// `-0.0` counts as zero); then `ERG=`, `WGT=`, `PAR=`. A single surviving
/// line keeps the upstream inline `ERG=<E>`; multiple lines emit `ERG=D1`
/// with the energies on `SI1 L` and the probabilities on `SP1 D` (surface
/// syntax verified against an MCNP input parser; sampling semantics are the
/// caller's responsibility — see the module docs).
pub fn sdef_card(
    lines: &[(f64, f64)],
    source: &PointSource,
    version: u32,
) -> Result<(Vec<(f64, f64)>, String), Error> {
    let bins = normalize_decay_lines(lines)?;
    let text = render_card(&bins, source, version)?;
    Ok((bins, text))
}

fn render_card(bins: &[(f64, f64)], source: &PointSource, version: u32) -> Result<String, Error> {
    for (field, value) in [
        ("x", source.x),
        ("y", source.y),
        ("z", source.z),
        ("u", source.u),
        ("v", source.v),
        ("w", source.w),
    ] {
        if !value.is_finite() {
            return Err(Error::NonFiniteSourceField(field));
        }
    }
    if !source.weight.is_finite() {
        return Err(Error::NonFiniteSourceField("weight"));
    }
    if source.weight <= 0.0 {
        return Err(Error::NonPositiveWeight(source.weight));
    }
    let par = match version {
        5 => source.particle.mcnp(),
        6 => source.particle.mcnp6(),
        _ => return Err(Error::UnsupportedMcnpVersion(version)),
    }
    .ok_or(Error::ParticleNotScorable(source.particle.name(), version))?;

    let mut out = format!(
        "SDEF POS={} {} {}",
        fmt_g6(source.x),
        fmt_g6(source.y),
        fmt_g6(source.z)
    );
    if source.u != 0.0 || source.v != 0.0 || source.w != 0.0 {
        out.push_str(&format!(
            "\n     VEC={} {} {} DIR=1",
            fmt_g6(source.u),
            fmt_g6(source.v),
            fmt_g6(source.w)
        ));
    }
    if bins.len() == 1 {
        out.push_str(&format!("\n     ERG={}", fmt_g6(bins[0].0)));
    } else {
        // Multi-line: discrete-distribution form (surface syntax verified;
        // sampling semantics are the caller's responsibility).
        out.push_str("\n     ERG=D1");
    }
    out.push_str(&format!("\n     WGT={}", fmt_g6(source.weight)));
    out.push_str(&format!("\n     PAR={par}"));
    if bins.len() > 1 {
        let energies: Vec<String> = bins.iter().map(|b| fmt_g6(b.0)).collect();
        let probs: Vec<String> = bins.iter().map(|b| fmt_g6(b.1)).collect();
        push_list(&mut out, "SI1 L", &energies);
        push_list(&mut out, "SP1 D", &probs);
    }
    Ok(out)
}

/// Append a data card (`head` plus single-space-separated entries), wrapping
/// onto 5-space continuation lines so no line exceeds [`CARD_WIDTH`]. The
/// wrap indent matches the SDEF continuation shape; wrapping itself is a
/// renderer convention (entry order is unchanged).
fn push_list(out: &mut String, head: &str, entries: &[String]) {
    out.push('\n');
    out.push_str(head);
    let mut col = head.len();
    for e in entries {
        let need = 1 + e.len();
        if col + need > CARD_WIDTH {
            out.push_str("\n     ");
            col = 5;
        }
        out.push(' ');
        out.push_str(e);
        col += need;
    }
}

/// Format like C++ `std::ostream <<` default formatting (`%g` with precision
/// 6, trailing zeros stripped): what the upstream C++ emitter produces for
/// the plain `double` fields it streams into the card text.
fn fmt_g6(x: f64) -> String {
    debug_assert!(x.is_finite());
    if x == 0.0 {
        return if x.is_sign_negative() {
            "-0".into()
        } else {
            "0".into()
        };
    }
    const P: i32 = 6;
    let sci = format!("{:.*e}", (P - 1) as usize, x);
    let (mant, exp) = sci.split_once('e').expect("scientific format has e");
    let exp: i32 = exp.parse().expect("scientific exponent is an integer");
    let mant = trim_zeros(mant);
    if !(-4..P).contains(&exp) {
        format!("{mant}e{exp:+03}")
    } else {
        let decimals = (P - 1 - exp).max(0) as usize;
        trim_zeros(&format!("{x:.decimals$}"))
    }
}

/// Strip trailing fractional zeros (`"1.40000"` -> `"1.4"`, `"14.000"` ->
/// `"14"`); integer strings pass through.
fn trim_zeros(s: &str) -> String {
    if s.contains('.') {
        let t = s.trim_end_matches('0');
        t.strip_suffix('.').unwrap_or(t).to_string()
    } else {
        s.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const CONT: &str = "\n     ";

    fn beam() -> PointSource {
        PointSource {
            x: 1.0,
            y: 2.0,
            z: 3.0,
            u: 0.0,
            v: 0.0,
            w: 1.0,
            weight: 0.5,
            particle: ParticleId::Photon,
        }
    }

    fn lines() -> Vec<(f64, f64)> {
        vec![(1.17, 1.0), (0.662, 2.0), (1.33, 1.0)]
    }

    // --- C++ defaultfloat (precision 6) number formatting ----------------

    #[test]
    fn fmt_g6_matches_cpp_defaultfloat() {
        assert_eq!(fmt_g6(14.0), "14");
        assert_eq!(fmt_g6(0.662), "0.662");
        assert_eq!(fmt_g6(0.85), "0.85");
        assert_eq!(fmt_g6(1.0 / 3.0), "0.333333");
        assert_eq!(fmt_g6(100.0), "100");
        assert_eq!(fmt_g6(0.0001), "0.0001");
        assert_eq!(fmt_g6(-1.5), "-1.5");
        assert_eq!(fmt_g6(0.0), "0");
        assert_eq!(fmt_g6(-0.0), "-0");
        // Scientific branch, C++ exponent padding included.
        assert_eq!(fmt_g6(1e-5), "1e-05");
        assert_eq!(fmt_g6(2.5e-10), "2.5e-10");
        assert_eq!(fmt_g6(1234567.0), "1.23457e+06");
        assert_eq!(fmt_g6(999999.9), "1e+06");
        assert_eq!(fmt_g6(123456.0), "123456");
    }

    // --- byte-exact card goldens ------------------------------------------

    #[test]
    fn single_line_isotropic_golden() {
        let (bins, card) = sdef_card(&[(0.662, 2.0)], &PointSource::default(), 5).unwrap();
        assert_eq!(bins, vec![(0.662, 1.0)]);
        assert_eq!(
            card,
            "SDEF POS=0 0 0\
             \n     ERG=0.662\
             \n     WGT=1\
             \n     PAR=n"
        );
    }

    #[test]
    fn single_line_beam_golden() {
        let (bins, card) = sdef_card(&[(0.662, 2.0)], &beam(), 5).unwrap();
        assert_eq!(bins, vec![(0.662, 1.0)]);
        assert_eq!(
            card,
            "SDEF POS=1 2 3\
             \n     VEC=0 0 1 DIR=1\
             \n     ERG=0.662\
             \n     WGT=0.5\
             \n     PAR=p"
        );
    }

    #[test]
    fn multi_line_distribution_golden() {
        let (bins, card) = sdef_card(&lines(), &PointSource::default(), 5).unwrap();
        assert_eq!(bins, vec![(0.662, 0.5), (1.17, 0.25), (1.33, 0.25)]);
        assert_eq!(
            card,
            "SDEF POS=0 0 0\
             \n     ERG=D1\
             \n     WGT=1\
             \n     PAR=n\
             \nSI1 L 0.662 1.17 1.33\
             \nSP1 D 0.5 0.25 0.25"
        );
    }

    #[test]
    fn zero_intensity_line_drops_to_single_line_form() {
        let (bins, card) =
            sdef_card(&[(0.662, 2.0), (1.17, 0.0)], &PointSource::default(), 5).unwrap();
        assert_eq!(bins, vec![(0.662, 1.0)]);
        assert!(card.contains("\n     ERG=0.662\n"));
        assert!(!card.contains("SI1"));
    }

    // --- normalization arithmetic (E9) ------------------------------------

    #[test]
    fn merges_duplicates_sums_and_sorts() {
        // Unsorted, with a duplicate that must merge before normalizing.
        let bins = normalize_decay_lines(&[(1.0, 2.0), (0.5, 1.0), (1.0, 1.0)]).unwrap();
        assert_eq!(bins, vec![(0.5, 0.25), (1.0, 0.75)]);
        let sum: f64 = bins.iter().map(|b| b.1).sum();
        assert!((sum - 1.0).abs() < 1e-15);
    }

    #[test]
    fn normalization_is_order_independent() {
        let a = normalize_decay_lines(&lines()).unwrap();
        let mut shuffled = lines();
        shuffled.reverse();
        let b = normalize_decay_lines(&shuffled).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn three_way_normalization_arithmetic() {
        let bins = normalize_decay_lines(&[(0.662, 2.0), (1.17, 1.0), (1.33, 1.0)]).unwrap();
        assert_eq!(bins[0], (0.662, 0.5));
        assert!((bins[1].1 - 0.25).abs() < 1e-15);
        assert!((bins[2].1 - 0.25).abs() < 1e-15);
        assert!(bins.windows(2).all(|w| w[0].0 < w[1].0));
    }

    // --- edge cases: specified, not inherited ------------------------------

    #[test]
    fn empty_lines_rejected() {
        assert_eq!(normalize_decay_lines(&[]), Err(Error::EmptyLines));
    }

    #[test]
    fn all_zero_intensities_rejected() {
        assert_eq!(
            normalize_decay_lines(&[(0.662, 0.0)]),
            Err(Error::NoEmission)
        );
    }

    #[test]
    fn negative_intensity_rejected() {
        assert_eq!(
            normalize_decay_lines(&[(0.662, -1.0)]),
            Err(Error::NegativeIntensity(-1.0))
        );
    }

    #[test]
    fn negative_energy_rejected() {
        assert_eq!(
            normalize_decay_lines(&[(-0.5, 1.0)]),
            Err(Error::NegativeEnergy(-0.5))
        );
    }

    #[test]
    fn non_finite_values_rejected() {
        assert!(
            matches!(normalize_decay_lines(&[(f64::NAN, 1.0)]), Err(Error::NonFiniteLine(v)) if v.is_nan())
        );
        assert_eq!(
            normalize_decay_lines(&[(0.662, f64::INFINITY)]),
            Err(Error::NonFiniteLine(f64::INFINITY))
        );
    }

    #[test]
    fn overflowing_total_rejected() {
        assert_eq!(
            normalize_decay_lines(&[(0.662, 1e308), (1.17, 1e308)]),
            Err(Error::NonFiniteLine(f64::INFINITY))
        );
    }

    // --- point-source validation + particle/version mapping ----------------

    #[test]
    fn bad_version_rejected() {
        assert_eq!(
            sdef_card(&[(0.662, 1.0)], &PointSource::default(), 4),
            Err(Error::UnsupportedMcnpVersion(4))
        );
    }

    #[test]
    fn unscorable_particle_rejected() {
        // Proton gains an MCNP designator only in version 6.
        let proton = PointSource {
            particle: ParticleId::Proton,
            ..PointSource::default()
        };
        assert_eq!(
            sdef_card(&[(0.662, 1.0)], &proton, 5),
            Err(Error::ParticleNotScorable("Proton", 5))
        );
        let (_, card) = sdef_card(&[(0.662, 1.0)], &proton, 6).unwrap();
        assert!(card.ends_with("\n     PAR=h"));
    }

    #[test]
    fn bad_weight_rejected() {
        for weight in [0.0, -1.0, f64::NAN] {
            let src = PointSource {
                weight,
                ..PointSource::default()
            };
            assert!(
                sdef_card(&[(0.662, 1.0)], &src, 5).is_err(),
                "weight {weight}"
            );
        }
    }

    #[test]
    fn non_finite_source_field_rejected() {
        let src = PointSource {
            z: f64::NAN,
            ..PointSource::default()
        };
        assert_eq!(
            sdef_card(&[(0.662, 1.0)], &src, 5),
            Err(Error::NonFiniteSourceField("z"))
        );
    }

    #[test]
    fn negative_zero_direction_counts_as_isotropic() {
        let src = PointSource {
            u: -0.0,
            ..PointSource::default()
        };
        let (_, card) = sdef_card(&[(0.662, 1.0)], &src, 5).unwrap();
        assert!(!card.contains("VEC="));
    }

    // --- card wrapping ------------------------------------------------------

    #[test]
    fn long_lists_wrap_under_card_width() {
        let lines: Vec<(f64, f64)> = (1..=12).map(|i| (0.1 * i as f64, i as f64)).collect();
        let (_, card) = sdef_card(&lines, &PointSource::default(), 5).unwrap();
        for line in card.lines() {
            assert!(line.len() <= CARD_WIDTH, "line over {CARD_WIDTH}: {line:?}");
        }
        let cont_count = card.matches(CONT.trim_start_matches('\n')).count();
        assert!(cont_count >= 4, "SDEF fields + wrapped SI/SP continuations");
        // Continuation lines carry the 5-space indent.
        let wrapped = card
            .lines()
            .filter(|l| l.starts_with("      0.102564"))
            .count();
        assert_eq!(
            wrapped, 1,
            "wrapped SP1 continuation with 5-space indent + separator"
        );
        // Entry order survives wrapping: probabilities stay ascending.
        let sp: Vec<f64> = card
            .split("SP1 D ")
            .nth(1)
            .unwrap()
            .split_whitespace()
            .filter_map(|t| t.parse().ok())
            .collect();
        assert_eq!(sp.len(), 12);
        assert!(sp.windows(2).all(|w| w[0] < w[1]));
    }
}
