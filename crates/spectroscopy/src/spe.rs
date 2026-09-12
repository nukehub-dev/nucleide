//! Text `.spe` readers: dollar format and plain format.
//!
//! The dollar reader takes files whose first line is `$SPEC_ID:` (anything
//! else is an [`Error::UnsupportedFormat`]); the plain reader rejects that
//! magic. Pinned quirks, all reproduced from the upstream readers:
//!
//! * Duplicate tags: first occurrence wins (dollar); later keys overwrite
//!   (plain).
//! * Missing dollar tags raise [`Error::MissingTag`]; unknown keys are
//!   ignored by both readers (`$ROI:`, `$PRESETS:`, `$ENER_FIT:` are never
//!   read).
//! * `$MEAS_TIM:` stores live time first, real time second.
//! * `$DATA:` gives `start last`; the channel count is `last + 1`, followed
//!   by that many one-float-per-line data lines.
//! * `$MCA_CAL:`/`$SHAPE_CAL:` triplets sit two lines below the tag (the
//!   line between holds the count); only the first three tokens parse, so a
//!   trailing `keV` unit is ignored.
//! * Plain `Energy Fit:`/`FWHM Fit:` values sit at single-space-split
//!   indices 0/2/4 (double-space separated values); the acquisition time
//!   re-joins three colon parts.
//! * Dollar channel labels are positional `0..len` even when
//!   `start_chan_num` is nonzero; E3–E5 index counts positionally too.

use crate::{spectrum::GammaSpectrum, Error};

/// Find the first line exactly equal to `tag`.
fn tag_index(lines: &[&str], tag: &'static str) -> Result<usize, Error> {
    lines
        .iter()
        .position(|l| *l == tag)
        .ok_or(Error::MissingTag(tag))
}

/// Parse one whitespace-separated float, trimming first.
fn parse_f64(text: &str, what: &'static str) -> Result<f64, Error> {
    text.trim()
        .parse::<f64>()
        .map_err(|_| Error::MalformedValue(what))
}

/// Parse one trimmed integer.
fn parse_int(text: &str, what: &'static str) -> Result<i64, Error> {
    text.trim()
        .parse::<i64>()
        .map_err(|_| Error::MalformedValue(what))
}

/// Read a dollar-format `.spe` file already loaded as text.
///
/// `file_name` is recorded verbatim on the spectrum (empty for pure-text
/// parsing).
pub fn parse_dollar_spe(text: &str, file_name: &str) -> Result<GammaSpectrum, Error> {
    let lines: Vec<&str> = text.lines().collect();
    if lines.first() != Some(&"$SPEC_ID:") {
        return Err(Error::UnsupportedFormat("expected a $SPEC_ID: first line"));
    }
    let after = |tag: &'static str, offset: usize| -> Result<&str, Error> {
        let i = tag_index(&lines, tag)?;
        lines.get(i + offset).copied().ok_or(Error::TruncatedData)
    };

    let mut spec = GammaSpectrum::new();
    spec.file_name = file_name.to_string();
    spec.spectrum.spec_name = after("$SPEC_ID:", 1)?.trim().to_string();

    let rem1: Vec<&str> = after("$SPEC_REM:", 1)?.split(' ').collect();
    spec.det_id = rem1.get(1).copied().unwrap_or("").trim().to_string();
    if spec.det_id.is_empty() {
        return Err(Error::MalformedValue("DET# id"));
    }
    let rem2: Vec<&str> = after("$SPEC_REM:", 2)?.split(' ').collect();
    spec.det_descp = rem2.get(1).copied().unwrap_or("").trim().to_string();
    if spec.det_descp.is_empty() {
        return Err(Error::MalformedValue("DETDESC# description"));
    }

    let date: Vec<&str> = after("$DATE_MEA:", 1)?.split(' ').collect();
    spec.start_date = date.first().copied().unwrap_or("").to_string();
    spec.start_time = date.get(1).copied().unwrap_or("").to_string();
    if spec.start_date.is_empty() || spec.start_time.is_empty() {
        return Err(Error::MalformedValue("DATE_MEA date/time"));
    }

    // Live time first, real time second.
    let times: Vec<&str> = after("$MEAS_TIM:", 1)?.split(' ').collect();
    if times.len() < 2 {
        return Err(Error::MalformedValue("MEAS_TIM live/real"));
    }
    spec.live_time = parse_f64(times[0], "MEAS_TIM live")?;
    spec.real_time = parse_f64(times[1], "MEAS_TIM real")?;

    let data: Vec<&str> = after("$DATA:", 1)?.split(' ').collect();
    if data.len() < 2 {
        return Err(Error::MalformedValue("DATA start/last"));
    }
    spec.spectrum.start_chan_num = parse_int(data[0], "DATA start")?;
    let last = parse_int(data[1], "DATA last")?;
    if last < spec.spectrum.start_chan_num {
        return Err(Error::MalformedValue("DATA last below start"));
    }
    // Channel count is `last + 1`, read positionally (labels stay 0-based).
    let n = (last + 1) as usize;
    spec.spectrum.num_channels = n;
    let data_at = tag_index(&lines, "$DATA:")?;
    if lines.len() < data_at + 2 + n {
        return Err(Error::TruncatedData);
    }
    for line in lines.iter().skip(data_at + 2).take(n) {
        spec.spectrum.counts.push(parse_f64(line, "DATA counts")?);
    }

    // Triplet two lines below the tag; extra tokens (e.g. `keV`) ignored.
    for (tag, slot) in [
        ("$MCA_CAL:", &mut spec.calib_e_fit),
        ("$SHAPE_CAL:", &mut spec.calib_fwhm_fit),
    ] {
        let coeffs: Vec<&str> = after(tag, 2)?.split(' ').collect();
        if coeffs.len() < 3 {
            return Err(Error::MalformedValue("calibration triplet"));
        }
        slot.push(parse_f64(coeffs[0], "calibration a0")?);
        slot.push(parse_f64(coeffs[1], "calibration a1")?);
        slot.push(parse_f64(coeffs[2], "calibration a2")?);
    }

    spec.spectrum.channels = (0..n).map(|c| c as f64).collect();
    spec.calc_ebins()?;
    Ok(spec)
}

/// Read a plain-format `.spe` file already loaded as text.
///
/// `key: value` header lines plus `channel: counts` data lines after the
/// `SPECTRUM` marker. Unknown keys are ignored.
pub fn parse_plain_spe(text: &str, file_name: &str) -> Result<GammaSpectrum, Error> {
    let lines: Vec<&str> = text.lines().collect();
    if lines.first() == Some(&"$SPEC_ID:") {
        return Err(Error::UnsupportedFormat(
            "dollar format is not supported by this function",
        ));
    }
    let mut spec = GammaSpectrum::new();
    spec.file_name = file_name.to_string();
    let mut inspec = false;
    for item in &lines {
        let parts: Vec<&str> = item.split(':').collect();
        if inspec && parts.len() > 1 {
            spec.spectrum.channels.push(parse_f64(parts[0], "channel")?);
            spec.spectrum.counts.push(parse_f64(parts[1], "counts")?);
        }
        match parts[0] {
            "Spectrum name" if parts.len() > 1 => {
                spec.spectrum.spec_name = parts[1].trim().to_string();
            }
            "Detector ID" if parts.len() > 1 => {
                spec.det_id = parts[1].trim().to_string();
            }
            "Detector description" if parts.len() > 1 => {
                spec.det_descp = parts[1].trim().to_string();
            }
            "Real Time" if parts.len() > 1 => {
                spec.real_time = parse_f64(parts[1], "Real Time")?;
            }
            "Live Time" if parts.len() > 1 => {
                spec.live_time = parse_f64(parts[1], "Live Time")?;
            }
            "Acquisition start date" if parts.len() > 1 => {
                spec.start_date = parts[1].trim().to_string();
            }
            "Acquisition start time" if parts.len() > 3 => {
                spec.start_time = format!("{}:{}:{}", parts[1].trim(), parts[2], parts[3].trim());
            }
            "Starting channel number" if parts.len() > 1 => {
                spec.spectrum.start_chan_num = parse_int(parts[1], "Starting channel number")?;
            }
            "Number of channels" if parts.len() > 1 => {
                let n = parse_int(parts[1], "Number of channels")?;
                if n < 0 {
                    return Err(Error::MalformedValue("Number of channels"));
                }
                spec.spectrum.num_channels = n as usize;
            }
            "Energy Fit" if parts.len() > 1 => {
                // Double-space separated values: single-space split keeps
                // the gaps, so the three coefficients sit at 0/2/4.
                let vals: Vec<&str> = parts[1].trim().split(' ').collect();
                if vals.len() < 5 {
                    return Err(Error::MalformedValue("Energy Fit coefficients"));
                }
                spec.calib_e_fit = vec![
                    parse_f64(vals[0], "Energy Fit a0")?,
                    parse_f64(vals[2], "Energy Fit a1")?,
                    parse_f64(vals[4], "Energy Fit a2")?,
                ];
            }
            "FWHM Fit" if parts.len() > 1 => {
                let vals: Vec<&str> = parts[1].trim().split(' ').collect();
                if vals.len() < 5 {
                    return Err(Error::MalformedValue("FWHM Fit coefficients"));
                }
                spec.calib_fwhm_fit = vec![
                    parse_f64(vals[0], "FWHM Fit a0")?,
                    parse_f64(vals[2], "FWHM Fit a1")?,
                    parse_f64(vals[4], "FWHM Fit a2")?,
                ];
            }
            "SPECTRUM" => inspec = true,
            _ => {}
        }
    }
    spec.calc_ebins()?;
    Ok(spec)
}

#[cfg(test)]
mod tests {
    use super::*;

    const DOLLAR: &str = "\
$SPEC_ID:
SYNTH
$SPEC_REM:
DET# 7
DETDESC# SYNDET
$DATE_MEA:
02/03/2026 09:15:00
$MEAS_TIM:
90 100
$DATA:
0 3
1
2
4
8
$ROI:
0
$PRESETS:
None
0
0
$ENER_FIT:
1.5 2.0
$MCA_CAL:
3
1.500000E+000 2.000000E+000 5.000000E-001 keV
$SHAPE_CAL:
3
7.000000E-001 5.000000E-004 2.000000E-007
";

    #[test]
    fn dollar_pins_header_quirks() {
        let s = parse_dollar_spe(DOLLAR, "syn.spe").unwrap();
        assert_eq!(s.spectrum.spec_name, "SYNTH");
        assert_eq!(s.det_id, "7");
        assert_eq!(s.det_descp, "SYNDET");
        assert_eq!(s.start_date, "02/03/2026");
        assert_eq!(s.start_time, "09:15:00");
        // MEAS_TIM is live-first: real = 100, live = 90.
        assert_eq!(
            (s.real_time, s.live_time, s.dead_time()),
            (100.0, 90.0, 10.0)
        );
        assert_eq!(s.spectrum.start_chan_num, 0);
        assert_eq!(s.spectrum.num_channels, 4);
        assert_eq!(s.spectrum.counts, vec![1.0, 2.0, 4.0, 8.0]);
        assert_eq!(s.spectrum.channels, vec![0.0, 1.0, 2.0, 3.0]);
        assert_eq!(s.calib_e_fit, vec![1.5, 2.0, 0.5]);
        assert_eq!(s.calib_fwhm_fit, vec![0.7, 0.0005, 0.0000002]);
        // E6 on positional channels: ebin[2] = 1.5 + 4 + 2 = 7.5.
        assert_eq!(s.spectrum.ebin, vec![1.5, 4.0, 7.5, 12.0]);
    }

    #[test]
    fn fwhm_triplet_is_carried_never_evaluated() {
        // Gap record: the upstream analysis subset parses the
        // `$SHAPE_CAL:`/`FWHM Fit` triplet and stores it, but defines no
        // evaluation routine over it — so neither does this crate. Pin the
        // carried-without-evaluation posture: mutating the triplet leaves
        // the E6 energy bins untouched.
        let mut s = parse_dollar_spe(DOLLAR, "").unwrap();
        let before = s.spectrum.ebin.clone();
        s.calib_fwhm_fit = vec![0.0, 0.0, 0.0];
        s.calc_ebins().unwrap();
        assert_eq!(s.spectrum.ebin, before);
    }

    #[test]
    fn dollar_rejects_plain_magic_and_missing_tags() {
        assert_eq!(
            parse_dollar_spe("Spectrum name: x\n", ""),
            Err(Error::UnsupportedFormat("expected a $SPEC_ID: first line"))
        );
        let missing_cal = DOLLAR.replace(
            "$MCA_CAL:\n3\n1.500000E+000 2.000000E+000 5.000000E-001 keV\n",
            "",
        );
        assert_eq!(
            parse_dollar_spe(&missing_cal, ""),
            Err(Error::MissingTag("$MCA_CAL:"))
        );
    }

    #[test]
    fn dollar_first_tag_wins_and_start_stays_positional() {
        // Duplicate $DATE_MEA:: the first occurrence wins.
        let dup = DOLLAR.replacen(
            "$DATE_MEA:\n02/03/2026 09:15:00",
            "$DATE_MEA:\n01/01/2000 00:00:00\n$DATE_MEA:\n02/03/2026 09:15:00",
            1,
        );
        // Note: inserting a second tag line shifts nothing else; the first
        // $DATE_MEA: still resolves to the 2000 date.
        let s = parse_dollar_spe(&dup, "").unwrap();
        assert_eq!(s.start_date, "01/01/2000");

        // Nonzero start: labels stay 0-based, counts stay positional.
        let shifted = DOLLAR.replace("$DATA:\n0 3", "$DATA:\n2 2");
        let s = parse_dollar_spe(&shifted, "").unwrap();
        assert_eq!(s.spectrum.start_chan_num, 2);
        assert_eq!(s.spectrum.num_channels, 3);
        assert_eq!(s.spectrum.channels, vec![0.0, 1.0, 2.0]);
        assert_eq!(s.spectrum.counts, vec![1.0, 2.0, 4.0]);
    }

    const PLAIN: &str = "\
Spectrum name:  SYNTH-PLAIN
Detector ID:  7
Detector description: SYNDET
Mystery Key:  ignored entirely
Real Time:  100.5
Live Time:  90.25
Acquisition start date:  03-Feb-2026
Acquisition start time:  09:15:00
Starting channel number:  0
Number of channels:  4
Energy Fit:  1.500000E+000  2.000000E+000  5.000000E-001
FWHM Fit:  7.000000E-001  5.000000E-004  2.000000E-007
SPECTRUM

     0:    1.00000E+000
     1:    2.00000E+000
     2:    4.00000E+000
     3:    8.00000E+000
";

    #[test]
    fn plain_pins_keys_and_ignores_unknown() {
        let s = parse_plain_spe(PLAIN, "syn.spe").unwrap();
        assert_eq!(s.spectrum.spec_name, "SYNTH-PLAIN");
        assert_eq!(s.det_id, "7");
        assert_eq!(s.det_descp, "SYNDET");
        assert_eq!((s.real_time, s.live_time), (100.5, 90.25));
        assert_eq!(s.start_date, "03-Feb-2026");
        assert_eq!(s.start_time, "09:15:00");
        assert_eq!(s.spectrum.channels, vec![0.0, 1.0, 2.0, 3.0]);
        assert_eq!(s.spectrum.counts, vec![1.0, 2.0, 4.0, 8.0]);
        assert_eq!(s.calib_e_fit, vec![1.5, 2.0, 0.5]);
        assert_eq!(s.spectrum.ebin, vec![1.5, 4.0, 7.5, 12.0]);
    }

    #[test]
    fn dollar_rejects_short_headers_and_truncated_data() {
        for (mutant, fragment) in [
            (DOLLAR.replace("DET# 7", "DET#"), "DET# id"),
            (
                DOLLAR.replace("DETDESC# SYNDET", "DETDESC#"),
                "DETDESC# description",
            ),
            (
                DOLLAR.replace("02/03/2026 09:15:00", "02/03/2026"),
                "DATE_MEA date/time",
            ),
            (DOLLAR.replace("90 100", "90"), "MEAS_TIM live/real"),
            (DOLLAR.replace("0 3", "0"), "DATA start/last"),
            (DOLLAR.replace("0 3", "3 2"), "DATA last below start"),
            (
                DOLLAR.replace("1.500000E+000 2.000000E+000 5.000000E-001 keV", "1.5 2.0"),
                "calibration triplet",
            ),
        ] {
            let err = parse_dollar_spe(&mutant, "").unwrap_err();
            assert!(err.to_string().contains(fragment), "{fragment}: {err}");
        }
        // Fewer count lines than the declared channel count.
        let truncated = DOLLAR.split("4\n8\n").next().unwrap().to_string() + "4\n";
        assert_eq!(parse_dollar_spe(&truncated, ""), Err(Error::TruncatedData));
    }

    #[test]
    fn non_numeric_values_name_their_field() {
        for (mutant, fragment) in [
            (DOLLAR.replace("90 100", "XX 100"), "MEAS_TIM live"),
            (DOLLAR.replace("0 3", "X 3"), "DATA start"),
            (DOLLAR.replacen("\n4\n", "\nXX\n", 1), "DATA counts"),
            (
                DOLLAR.replace("1.500000E+000 2.000000E+000 5.000000E-001 keV", "X Y Z"),
                "calibration a0",
            ),
        ] {
            let err = parse_dollar_spe(&mutant, "").unwrap_err();
            assert!(err.to_string().contains(fragment), "{fragment}: {err}");
        }
        for (mutant, fragment) in [
            (
                PLAIN.replace("     0:    1.00000E+000", "     X:    1.00000E+000"),
                "channel",
            ),
            (
                PLAIN.replace("     0:    1.00000E+000", "     0:    X"),
                "counts",
            ),
            (
                PLAIN.replace("Real Time:  100.5", "Real Time:  X"),
                "Real Time",
            ),
            (
                PLAIN.replace("Starting channel number:  0", "Starting channel number:  X"),
                "Starting channel number",
            ),
            (
                PLAIN.replace("Number of channels:  4", "Number of channels:  X"),
                "Number of channels",
            ),
            (
                PLAIN.replace(
                    "Energy Fit:  1.500000E+000  2.000000E+000  5.000000E-001",
                    "Energy Fit:  X  2.0  0.5",
                ),
                "Energy Fit a0",
            ),
        ] {
            let err = parse_plain_spe(&mutant, "").unwrap_err();
            assert!(err.to_string().contains(fragment), "{fragment}: {err}");
        }
    }

    #[test]
    fn plain_rejects_bad_counts_and_short_fits() {
        assert_eq!(
            parse_plain_spe(
                &PLAIN.replace("Number of channels:  4", "Number of channels:  -1"),
                ""
            ),
            Err(Error::MalformedValue("Number of channels"))
        );
        assert_eq!(
            parse_plain_spe(
                &PLAIN.replace(
                    "Energy Fit:  1.500000E+000  2.000000E+000  5.000000E-001",
                    "Energy Fit:  1.0"
                ),
                ""
            ),
            Err(Error::MalformedValue("Energy Fit coefficients"))
        );
        assert_eq!(
            parse_plain_spe(
                &PLAIN.replace(
                    "FWHM Fit:  7.000000E-001  5.000000E-004  2.000000E-007",
                    "FWHM Fit:  1.0"
                ),
                ""
            ),
            Err(Error::MalformedValue("FWHM Fit coefficients"))
        );
    }

    #[test]
    fn plain_rejects_dollar_magic() {
        assert_eq!(
            parse_plain_spe("$SPEC_ID:\nx\n", ""),
            Err(Error::UnsupportedFormat(
                "dollar format is not supported by this function"
            ))
        );
    }
}
