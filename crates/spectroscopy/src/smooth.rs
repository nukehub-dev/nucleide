//! Spectrum smoothing: rectangular (E1) and five-point (E2) passes.
//!
//! Both passes copy the edge channels that have no full neighbourhood and
//! smooth the interior. Inputs and outputs are plain count vectors; the
//! spectrum name bookkeeping (`+ " smoothed"`) is caller-side.

use crate::Error;

/// Rectangular smoothing (E1).
///
/// With `ext = (m - 1) / 2`, `smooth[i] = sum(counts[i-ext..=i+ext]) / m`
/// for `ext <= i < N - ext`; the first and last `ext` channels are copied.
/// `m` must be odd and at least 3.
pub fn rect_smooth(counts: &[f64], m: usize) -> Result<Vec<f64>, Error> {
    if m < 3 {
        return Err(Error::SmoothWidthTooSmall(m));
    }
    if m % 2 == 0 {
        return Err(Error::SmoothWidthEven(m));
    }
    if counts.is_empty() {
        return Err(Error::EmptyCounts);
    }
    let ext = (m - 1) / 2;
    let n = counts.len();
    let mut smooth = Vec::with_capacity(n);
    // Edges: copy the first `ext` channels (fewer when the spectrum is
    // shorter than the window; the interior loop below then stays empty).
    for v in counts.iter().take(ext.min(n)) {
        smooth.push(*v);
    }
    // Interior: full-window mean while the window fits.
    let mut i = ext;
    while i + ext < n {
        let sum: f64 = counts[i - ext..=i + ext].iter().sum();
        smooth.push(sum / m as f64);
        i += 1;
    }
    // Trailing edge: copy the rest.
    for v in counts.iter().skip(i) {
        smooth.push(*v);
    }
    debug_assert_eq!(smooth.len(), n);
    Ok(smooth)
}

/// Five-point smoothing (E2).
///
/// `smooth[i] = (c[i-2] + c[i+2] + 2*c[i+1] + 2*c[i-1] + 3*c[i]) / 9` for
/// `2 <= i < N - 2`; the first and last two channels are copied. The
/// weighting is the low-statistics recommendation cited by the upstream
/// module (Phillips 1978); the citation is carried as text only.
pub fn five_point_smooth(counts: &[f64]) -> Result<Vec<f64>, Error> {
    let n = counts.len();
    if n < 4 {
        return Err(Error::TooFewChannels(n));
    }
    let mut smooth = Vec::with_capacity(n);
    smooth.push(counts[0]);
    smooth.push(counts[1]);
    let mut i = 2;
    while i + 2 < n {
        smooth.push(
            (counts[i - 2]
                + counts[i + 2]
                + 2.0 * counts[i + 1]
                + 2.0 * counts[i - 1]
                + 3.0 * counts[i])
                / 9.0,
        );
        i += 1;
    }
    smooth.push(counts[i]);
    smooth.push(counts[i + 1]);
    debug_assert_eq!(smooth.len(), n);
    Ok(smooth)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Hand-computed: m=3 on [1,3,2,6,4] gives
    // [1, (1+3+2)/3, (3+2+6)/3, (2+6+4)/3, 4].
    #[test]
    fn rect_m3_nonlinear() {
        let out = rect_smooth(&[1.0, 3.0, 2.0, 6.0, 4.0], 3).unwrap();
        assert_eq!(out, vec![1.0, 2.0, 11.0 / 3.0, 4.0, 4.0]);
    }

    // Hand-computed: m=5, ext=2 on [1,2,4,8,16,32,64] gives
    // [1,2,31/5,62/5,124/5,32,64].
    #[test]
    fn rect_m5_powers_of_two() {
        let out = rect_smooth(&[1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0], 5).unwrap();
        assert_eq!(out, vec![1.0, 2.0, 6.2, 12.4, 24.8, 32.0, 64.0]);
    }

    #[test]
    fn rect_rejects_bad_widths() {
        assert_eq!(
            rect_smooth(&[1.0, 2.0, 3.0], 2),
            Err(Error::SmoothWidthTooSmall(2))
        );
        assert_eq!(
            rect_smooth(&[1.0, 2.0, 3.0], 4),
            Err(Error::SmoothWidthEven(4))
        );
        assert_eq!(rect_smooth(&[], 3), Err(Error::EmptyCounts));
    }

    // Hand-computed on [2,5,1,6,3,8,4]:
    // i=2: (2+3+12+10+3)/9 = 30/9; i=3: (5+8+6+2+18)/9 = 39/9;
    // i=4: (1+4+16+12+9)/9 = 42/9.
    #[test]
    fn five_point_nonlinear() {
        let out = five_point_smooth(&[2.0, 5.0, 1.0, 6.0, 3.0, 8.0, 4.0]).unwrap();
        assert_eq!(
            out,
            vec![2.0, 5.0, 30.0 / 9.0, 39.0 / 9.0, 42.0 / 9.0, 8.0, 4.0]
        );
    }

    #[test]
    fn five_point_copies_edges_and_len4() {
        let out = five_point_smooth(&[7.0, 1.0, 5.0, 3.0]).unwrap();
        assert_eq!(out, vec![7.0, 1.0, 5.0, 3.0]);
        assert_eq!(
            five_point_smooth(&[1.0, 2.0, 3.0]),
            Err(Error::TooFewChannels(3))
        );
    }
}
