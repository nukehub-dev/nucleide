//! X-ray line algebra (E8) over caller-supplied atomic constants.
//!
//! The upstream routine reads fluorescence yields, shell ratios, and line
//! energies from its HDF5 atomic table; that table is not vendored here
//! (license unverified), so every constant arrives via [`AtomicData`].
//! Conversion inputs use [`Option`]: `None` plays the upstream NaN-sentinel
//! role (conversion absent). There is no combined upstream Python binding
//! for this routine (only a material method), so this explicit function is
//! the documented Nucleide surface.

/// Caller-supplied atomic constants for one element.
///
/// Yields and ratios are fractions; energies are in keV.
#[derive(Debug, Clone, PartialEq)]
pub struct AtomicData {
    /// K-shell fluorescence yield (fraction).
    pub k_shell_fluor: f64,
    /// L-shell fluorescence yield (fraction).
    pub l_shell_fluor: f64,
    /// Probability a K-shell hole is filled from the L shell (fraction).
    pub prob: f64,
    /// K-beta to K-alpha fluorescence ratio.
    pub kb_to_ka: f64,
    /// K-alpha-2 to K-alpha-1 fluorescence ratio.
    pub ka2_to_ka1: f64,
    /// K-alpha-1 line energy \[keV\].
    pub ka1_en_kev: f64,
    /// K-alpha-2 line energy \[keV\].
    pub ka2_en_kev: f64,
    /// K-beta line energy \[keV\].
    pub kb_en_kev: f64,
    /// L line energy \[keV\].
    pub l_en_kev: f64,
}

/// One X-ray line: energy \[keV\] plus intensity (photons per conversion).
#[derive(Debug, Clone, PartialEq)]
pub struct XrayLine {
    /// Line energy \[keV\].
    pub energy_kev: f64,
    /// Line intensity.
    pub intensity: f64,
}

/// X-ray lines for the given K/L conversion coefficients (E8).
///
/// With `k_conv`: `xk = k_fluor*k_conv`, `xka = xk/(1+kb_to_ka)`,
/// `xka1 = xka/(1+ka2_to_ka1)`, `xka2 = xka-xka1`, `xkb = xk-xka`, and —
/// when `l_conv` is also present — `xl = (l_conv+k_conv*prob)*l_fluor`.
/// With only `l_conv`: `xl = l_conv*l_fluor`. With neither, all intensities
/// are zero. Output order is always Ka1, Ka2, Kb, L.
pub fn xray_lines(data: &AtomicData, k_conv: Option<f64>, l_conv: Option<f64>) -> [XrayLine; 4] {
    let (mut xka1, mut xka2, mut xkb, mut xl) = (0.0, 0.0, 0.0, 0.0);
    if let Some(k) = k_conv {
        let xk = data.k_shell_fluor * k;
        let xka = xk / (1.0 + data.kb_to_ka);
        xka1 = xka / (1.0 + data.ka2_to_ka1);
        xka2 = xka - xka1;
        xkb = xk - xka;
        if let Some(l) = l_conv {
            xl = (l + k * data.prob) * data.l_shell_fluor;
        }
    } else if let Some(l) = l_conv {
        xl = l * data.l_shell_fluor;
    }
    [
        XrayLine {
            energy_kev: data.ka1_en_kev,
            intensity: xka1,
        },
        XrayLine {
            energy_kev: data.ka2_en_kev,
            intensity: xka2,
        },
        XrayLine {
            energy_kev: data.kb_en_kev,
            intensity: xkb,
        },
        XrayLine {
            energy_kev: data.l_en_kev,
            intensity: xl,
        },
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn data() -> AtomicData {
        AtomicData {
            k_shell_fluor: 0.9,
            l_shell_fluor: 0.4,
            prob: 0.8,
            kb_to_ka: 0.2,
            ka2_to_ka1: 0.5,
            ka1_en_kev: 10.0,
            ka2_en_kev: 20.0,
            kb_en_kev: 30.0,
            l_en_kev: 40.0,
        }
    }

    // Hand-computed with k=2, l=3: xk=1.8, xka=1.5, xka1=1.0, xka2=0.5,
    // xkb=0.3, xl=(3+1.6)*0.4=1.84 (float rounding: compare loosely).
    #[test]
    fn both_conversions() {
        let out = xray_lines(&data(), Some(2.0), Some(3.0));
        let got: Vec<(f64, f64)> = out.iter().map(|l| (l.energy_kev, l.intensity)).collect();
        let want = [(10.0, 1.0), (20.0, 0.5), (30.0, 0.3), (40.0, 1.84)];
        for ((ge, gi), (we, wi)) in got.iter().zip(want.iter()) {
            assert_eq!(ge, we);
            assert!((gi - wi).abs() < 1e-12, "{gi} vs {wi}");
        }
    }

    #[test]
    fn l_only_and_neither() {
        let l_only = xray_lines(&data(), None, Some(3.0));
        assert!((l_only[3].intensity - 1.2).abs() < 1e-12);
        assert!(l_only[..3].iter().all(|l| l.intensity == 0.0));

        let neither = xray_lines(&data(), None, None);
        assert!(neither.iter().all(|l| l.intensity == 0.0));
        // Energies still echo the caller constants.
        assert_eq!(neither[0].energy_kev, 10.0);
    }

    #[test]
    fn k_only_leaves_l_dark() {
        let out = xray_lines(&data(), Some(2.0), None);
        assert_eq!(out[0].intensity, 1.0);
        assert_eq!(out[3].intensity, 0.0);
    }
}
