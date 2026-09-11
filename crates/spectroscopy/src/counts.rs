//! Peak counting: background (E3), gross counts (E4), net counts (E5).
//!
//! Channel bounds `c1`/`c2` index the counts vector **positionally** (the
//! upstream quirk: a nonzero `start_chan_num` never shifts the slice).
//! Indexing follows Python-slice semantics — out-of-range stops clamp to
//! the vector ends — while the `c2 <= max(channels)` range check runs first.

use crate::Error;

/// Resolve a possibly-negative index the way a Python slice does.
fn resolve(idx: i64, len: i64) -> usize {
    let v = if idx < 0 { len + idx } else { idx };
    v.clamp(0, len) as usize
}

/// Sum `counts` over the half-open positional range `[start, stop)`.
///
/// Like a Python slice, an inverted range (possible when a negative start
/// resolves past the stop, e.g. `counts[-2:0]`) sums to zero instead of
/// panicking.
fn slice_sum(counts: &[f64], start: i64, stop: i64) -> f64 {
    let n = counts.len() as i64;
    let (s, e) = (resolve(start, n), resolve(stop, n));
    if s >= e {
        return 0.0;
    }
    counts[s..e].iter().sum()
}

/// Shared E3/E4 range validation, in upstream check order.
fn check_range(channels: &[f64], c1: i64, c2: i64) -> Result<(), Error> {
    if c1 > c2 {
        return Err(Error::BadChannelRange { c1, c2 });
    }
    if c1 < 0 {
        return Err(Error::NegativeChannel(c1));
    }
    let max = channels
        .iter()
        .copied()
        .reduce(|a, b| if a >= b { a } else { b })
        .ok_or(Error::EmptyChannels)?;
    if c2 as f64 > max {
        return Err(Error::ChannelOutOfRange { c2, max });
    }
    Ok(())
}

/// Background under a peak (E3, `m == 1` only).
///
/// `low = sum(counts[c1-2..c1])`, `high = sum(counts[c2..c2+2])`, and
/// `bg = (low + high) * (c2 - c1 + 1) / 6`. Any other `m` is an error:
/// no Compton/step model exists upstream.
pub fn calc_bg(counts: &[f64], channels: &[f64], c1: i64, c2: i64, m: i64) -> Result<f64, Error> {
    check_range(channels, c1, c2)?;
    if m != 1 {
        return Err(Error::BadBackgroundMethod(m));
    }
    let low = slice_sum(counts, c1 - 2, c1);
    let high = slice_sum(counts, c2, c2 + 2);
    Ok((low + high) * (c2 - c1 + 1) as f64 / 6.0)
}

/// Total counts between two channels (E4, half-open: excludes `c2`).
pub fn gross_count(counts: &[f64], channels: &[f64], c1: i64, c2: i64) -> Result<f64, Error> {
    check_range(channels, c1, c2)?;
    Ok(slice_sum(counts, c1, c2))
}

/// Net counts: gross minus background (E5).
pub fn net_counts(
    counts: &[f64],
    channels: &[f64],
    c1: i64,
    c2: i64,
    m: i64,
) -> Result<f64, Error> {
    Ok(gross_count(counts, channels, c1, c2)? - calc_bg(counts, channels, c1, c2, m)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn channels(n: usize) -> Vec<f64> {
        (0..n).map(|c| c as f64).collect()
    }

    // Hand-computed on counts [10,10,1,2,3,4,5,6,20,20,20], c1=2, c2=6:
    // low = 10+10 = 20, high = 5+6 = 11, bg = 31*5/6 = 155/6;
    // gross = 1+2+3+4 = 10; net = 10 - 155/6 = -95/6.
    const COUNTS: [f64; 11] = [10.0, 10.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 20.0, 20.0, 20.0];

    #[test]
    fn bg_gross_net_hand_values() {
        let ch = channels(11);
        assert_eq!(calc_bg(&COUNTS, &ch, 2, 6, 1).unwrap(), 155.0 / 6.0);
        assert_eq!(gross_count(&COUNTS, &ch, 2, 6).unwrap(), 10.0);
        // net = 10 - 155/6: one ulp away from -95/6, so compare loosely.
        let net = net_counts(&COUNTS, &ch, 2, 6, 1).unwrap();
        assert!((net - (-95.0 / 6.0)).abs() < 1e-12, "net = {net}");
    }

    #[test]
    fn range_errors_match_upstream_order() {
        let ch = channels(11);
        assert_eq!(
            calc_bg(&COUNTS, &ch, 6, 2, 1),
            Err(Error::BadChannelRange { c1: 6, c2: 2 })
        );
        assert_eq!(
            calc_bg(&COUNTS, &ch, -1, 2, 1),
            Err(Error::NegativeChannel(-1))
        );
        assert_eq!(
            calc_bg(&COUNTS, &ch, 2, 11, 1),
            Err(Error::ChannelOutOfRange { c2: 11, max: 10.0 })
        );
        assert_eq!(
            calc_bg(&COUNTS, &ch, 2, 6, 2),
            Err(Error::BadBackgroundMethod(2))
        );
        assert_eq!(
            gross_count(&COUNTS, &ch, 6, 2),
            Err(Error::BadChannelRange { c1: 6, c2: 2 })
        );
        assert_eq!(gross_count(&COUNTS, &[], 0, 0), Err(Error::EmptyChannels));
    }

    #[test]
    fn positional_indexing_ignores_label_offset() {
        // Channels labelled 5.. but indexed positionally (upstream quirk).
        let labels = vec![5.0, 6.0, 7.0, 8.0];
        let counts = vec![1.0, 2.0, 4.0, 8.0];
        assert_eq!(gross_count(&counts, &labels, 0, 2).unwrap(), 3.0);
    }

    #[test]
    fn small_c1_matches_python_slice_edges() {
        // counts[-2:0] and counts[-1:1] are empty upstream: low = 0.
        // c1=0: bg = 11*7/6; c1=1: bg = 11*6/6.
        let ch = channels(11);
        assert_eq!(calc_bg(&COUNTS, &ch, 0, 6, 1).unwrap(), 77.0 / 6.0);
        assert_eq!(calc_bg(&COUNTS, &ch, 1, 6, 1).unwrap(), 11.0);
    }
}
