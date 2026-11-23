//! Single-line activity estimates from caller-supplied peak measurements (E10).

use crate::{Error, Result};

/// Activity and counting standard uncertainty for one gamma line.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GammaLineActivity {
    /// Signed live-time-averaged activity estimate in Bq.
    pub activity_bq: f64,
    /// Standard uncertainty from the supplied peak-area uncertainty, in Bq.
    pub counting_uncertainty_bq: f64,
}

/// Estimate activity from a background-subtracted full-energy peak area.
///
/// Supply a finite signed area and a finite nonnegative standard uncertainty
/// in counts, positive live time in seconds, and the absolute peak efficiency
/// and emission probability as fractions in `(0, 1]`. Corrections for the
/// measurement and uncertainty in the exposure factors remain caller inputs.
pub fn gamma_line_activity(
    net_counts: f64,
    net_counts_uncertainty: f64,
    live_time_s: f64,
    peak_efficiency: f64,
    emission_probability: f64,
) -> Result<GammaLineActivity> {
    if !net_counts.is_finite() {
        return Err(Error::BadActivityInput("net_counts must be finite"));
    }
    if !net_counts_uncertainty.is_finite() || net_counts_uncertainty < 0.0 {
        return Err(Error::BadActivityInput(
            "net_counts_uncertainty must be finite and nonnegative",
        ));
    }
    if !live_time_s.is_finite() || live_time_s <= 0.0 {
        return Err(Error::BadActivityInput(
            "live_time_s must be finite and positive",
        ));
    }
    for (value, message) in [
        (peak_efficiency, "peak_efficiency must be in (0, 1]"),
        (
            emission_probability,
            "emission_probability must be in (0, 1]",
        ),
    ] {
        if !(value > 0.0 && value <= 1.0) {
            return Err(Error::BadActivityInput(message));
        }
    }

    let time_efficiency = live_time_s * peak_efficiency;
    let exposure = time_efficiency * emission_probability;
    // Recover the low parts of the two products before dividing.
    let error = live_time_s.mul_add(peak_efficiency, -time_efficiency);
    let exposure = exposure
        + error.mul_add(
            emission_probability,
            time_efficiency.mul_add(emission_probability, -exposure),
        );
    let to_activity = |counts: f64| {
        if counts == 0.0 {
            return Ok(counts);
        }
        let activity = if exposure.is_normal() {
            counts / exposure
        } else {
            // Preserve representable ratios when the exposure itself is subnormal.
            let (time, time_power) = binary_parts(live_time_s);
            let (efficiency, efficiency_power) = binary_parts(peak_efficiency);
            let (probability, probability_power) = binary_parts(emission_probability);
            let (area, area_power) = binary_parts(counts.abs());
            let (value, power) = binary_parts(area / (time * efficiency * probability));
            let power = power + area_power - time_power - efficiency_power - probability_power;
            if !(-1075..=1023).contains(&power) {
                return Err(Error::ActivityOutOfRange);
            }
            let normal_power = power.max(-1022);
            (value * 2.0_f64.powi(normal_power) * 2.0_f64.powi(power - normal_power))
                .copysign(counts)
        };
        if activity.is_finite() && activity != 0.0 {
            Ok(activity)
        } else {
            Err(Error::ActivityOutOfRange)
        }
    };
    Ok(GammaLineActivity {
        activity_bq: to_activity(net_counts)?,
        counting_uncertainty_bq: to_activity(net_counts_uncertainty)?,
    })
}

// Express a positive finite value as a significand in [1, 2) times 2^power.
fn binary_parts(value: f64) -> (f64, i32) {
    let (value, adjustment) = if value.is_subnormal() {
        (value * 2.0_f64.powi(54), -54)
    } else {
        (value, 0)
    };
    let bits = value.to_bits();
    let significand = f64::from_bits((bits & ((1_u64 << 52) - 1)) | 1.0_f64.to_bits());
    (significand, (bits >> 52) as i32 - 1023 + adjustment)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn computes_signed_activity_and_counting_uncertainty() {
        for (counts, expected) in [(120.0, 40.0), (0.0, 0.0), (-12.0, -4.0)] {
            let result = gamma_line_activity(counts, 12.0, 60.0, 0.1, 0.5).unwrap();
            assert_eq!(result.activity_bq, expected);
            assert_eq!(result.counting_uncertainty_bq, 4.0);
        }
    }

    #[test]
    fn preserves_finite_ratios_even_when_exposure_underflows() {
        let least_positive = f64::from_bits(1);
        let result =
            gamma_line_activity(-least_positive, least_positive, least_positive, 0.5, 0.5).unwrap();
        assert_eq!(result.activity_bq, -4.0);
        assert_eq!(result.counting_uncertainty_bq, 4.0);
        assert_eq!(
            gamma_line_activity(1.0, 0.0, least_positive, 1.0, 1.0),
            Err(Error::ActivityOutOfRange)
        );
        assert_eq!(
            gamma_line_activity(least_positive, 0.0, f64::MAX, 1.0, 1.0),
            Err(Error::ActivityOutOfRange)
        );
    }

    #[test]
    fn exposure_rounding_preserves_activity_near_the_upper_limit() {
        let result =
            gamma_line_activity(f64::MAX, 0.0, 7.0, 0.65, 0.219_780_219_780_219_78).unwrap();
        assert_eq!(result.activity_bq, f64::MAX);
    }
}
