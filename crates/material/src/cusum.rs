//! One-sided upper Page CUSUM change detector over a scalar stream.
//!
//! Dependency-free: the in-control mean and variance are tracked online
//! with Welford's recurrence, and the CUSUM statistic follows Page's
//! cumulative-sum rule with reference shift and alarm threshold scaled by
//! the running standard deviation. No simulator concepts (time, agents)
//! are involved: this is a pure function of the observed sequence.
//!
//! Update order per observation `x` (after skipping non-finite inputs,
//! which leave the state untouched):
//!
//! ```text
//! n += 1
//! mean += (x - mean) / n            (Welford)
//! M2   += (x - mean_old) * (x - mean)
//! std   = sqrt(M2 / (n - 1))        (sample std; 0 for n < 2)
//! S     = max(0, S + (x - mean) - k * std)
//! alarm = n > startup && S > h * std
//! ```

/// One-sided upper Page CUSUM detector with Welford running statistics.
///
/// `ref_shift_k` is the reference shift in units of the running standard
/// deviation (the smallest sustained mean shift worth flagging),
/// `alarm_h` is the alarm threshold in the same units, and `startup` is
/// the number of initial observations during which alarms are suppressed
/// while the running statistics settle.
#[derive(Debug, Clone, PartialEq)]
pub struct Cusum {
    ref_shift_k: f64,
    alarm_h: f64,
    startup: usize,
    count: usize,
    mean: f64,
    m2: f64,
    statistic: f64,
    alarmed: bool,
}

impl Cusum {
    /// Build a detector; fails when `ref_shift_k` is negative or
    /// non-finite, `alarm_h` is non-positive or non-finite.
    pub fn new(ref_shift_k: f64, alarm_h: f64, startup: usize) -> crate::Result<Self> {
        if !ref_shift_k.is_finite() || ref_shift_k < 0.0 {
            return Err(crate::Error::InvalidCusum(format!(
                "reference shift k must be finite and >= 0, got {ref_shift_k}"
            )));
        }
        if !alarm_h.is_finite() || alarm_h <= 0.0 {
            return Err(crate::Error::InvalidCusum(format!(
                "alarm threshold h must be finite and > 0, got {alarm_h}"
            )));
        }
        Ok(Self {
            ref_shift_k,
            alarm_h,
            startup,
            count: 0,
            mean: 0.0,
            m2: 0.0,
            statistic: 0.0,
            alarmed: false,
        })
    }

    /// Number of observations consumed (non-finite inputs do not count).
    pub fn count(&self) -> usize {
        self.count
    }

    /// Running mean of the observations seen so far (0 with no data).
    pub fn mean(&self) -> f64 {
        self.mean
    }

    /// Running sample variance (`M2 / (n - 1)`; 0 with fewer than 2 points).
    pub fn variance(&self) -> f64 {
        if self.count >= 2 {
            self.m2 / (self.count as f64 - 1.0)
        } else {
            0.0
        }
    }

    /// Running sample standard deviation.
    pub fn std(&self) -> f64 {
        self.variance().sqrt()
    }

    /// Current CUSUM statistic `S >= 0`.
    pub fn statistic(&self) -> f64 {
        self.statistic
    }

    /// Whether the detector is currently alarmed.
    ///
    /// Live (not latched): re-evaluated on every [`Cusum::update`] as
    /// `count > startup && statistic > h * std`.
    pub fn status(&self) -> bool {
        self.alarmed
    }

    /// Feed one observation; returns the resulting alarm [`Cusum::status`].
    ///
    /// Non-finite inputs are ignored (no state change) so a corrupt sample
    /// cannot poison the running statistics.
    pub fn update(&mut self, x: f64) -> bool {
        if !x.is_finite() {
            return self.alarmed;
        }
        self.count += 1;
        let n = self.count as f64;
        let old_mean = self.mean;
        self.mean += (x - old_mean) / n;
        self.m2 += (x - old_mean) * (x - self.mean);
        let std = self.std();
        let shifted = x - self.mean - self.ref_shift_k * std;
        self.statistic = (self.statistic + shifted).max(0.0);
        self.alarmed = self.count > self.startup && self.statistic > self.alarm_h * std;
        self.alarmed
    }

    /// Drop all observations; tuning parameters are kept.
    pub fn reset(&mut self) {
        self.count = 0;
        self.mean = 0.0;
        self.m2 = 0.0;
        self.statistic = 0.0;
        self.alarmed = false;
    }
}

impl Default for Cusum {
    /// Canonical safeguards tuning: `k = 0.5`, `h = 4.0`, `startup = 10`.
    fn default() -> Self {
        Self {
            ref_shift_k: 0.5,
            alarm_h: 4.0,
            startup: 10,
            count: 0,
            mean: 0.0,
            m2: 0.0,
            statistic: 0.0,
            alarmed: false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) {
        assert!(
            (a - b).abs() <= tol,
            "{a} != {b} within absolute tolerance {tol}"
        );
    }

    #[test]
    fn default_tuning_matches_canonical_values() {
        let cusum = Cusum::default();
        assert_eq!(cusum.ref_shift_k, 0.5);
        assert_eq!(cusum.alarm_h, 4.0);
        assert_eq!(cusum.startup, 10);
        assert_eq!(cusum.count(), 0);
        assert_eq!(cusum.statistic(), 0.0);
        assert!(!cusum.status());
        // `new` with the same tuning agrees with `default`.
        assert_eq!(Cusum::new(0.5, 4.0, 10).unwrap(), cusum);
    }

    #[test]
    fn invalid_parameters_rejected() {
        assert!(matches!(
            Cusum::new(f64::NAN, 4.0, 10),
            Err(crate::Error::InvalidCusum(_))
        ));
        assert!(matches!(
            Cusum::new(-0.1, 4.0, 10),
            Err(crate::Error::InvalidCusum(_))
        ));
        assert!(matches!(
            Cusum::new(0.5, 0.0, 10),
            Err(crate::Error::InvalidCusum(_))
        ));
        assert!(matches!(
            Cusum::new(0.5, f64::INFINITY, 10),
            Err(crate::Error::InvalidCusum(_))
        ));
    }

    #[test]
    fn welford_tracks_naive_mean_and_sample_variance() {
        // Hand-picked sequence; oracle moments recomputed naively below.
        let xs = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let mut cusum = Cusum::new(0.5, 4.0, 100).unwrap();
        let mut seen: Vec<f64> = Vec::new();
        for &x in &xs {
            seen.push(x);
            cusum.update(x);
            let n = seen.len() as f64;
            let mean = seen.iter().sum::<f64>() / n;
            let var = if seen.len() >= 2 {
                seen.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0)
            } else {
                0.0
            };
            close(cusum.mean(), mean, 1e-12);
            close(cusum.variance(), var, 1e-12);
            close(cusum.std(), var.sqrt(), 1e-12);
            assert_eq!(cusum.count(), seen.len());
        }
        // Naive oracle totals for the full sequence: sum 40, mean 5.0,
        // M2 = 9 + 1 + 1 + 1 + 0 + 0 + 4 + 16 = 32.0.
        close(cusum.mean(), 5.0, 1e-12);
        close(cusum.variance(), 32.0 / 7.0, 1e-12);
    }

    #[test]
    fn constant_stream_never_alarms() {
        let mut cusum = Cusum::default();
        for _ in 0..50 {
            assert!(!cusum.update(1.0));
        }
        assert_eq!(cusum.statistic(), 0.0);
        assert_eq!(cusum.mean(), 1.0);
        assert_eq!(cusum.variance(), 0.0);
        assert!(!cusum.status());
    }

    #[test]
    fn step_change_alarms_after_changepoint() {
        // Ten in-control points at 1.0, then a sustained step to 2.0.
        let mut cusum = Cusum::default();
        for _ in 0..10 {
            assert!(!cusum.update(1.0), "baseline must stay quiet");
        }
        assert_eq!(cusum.statistic(), 0.0);
        // First two post-change points build the statistic without firing.
        // At n = 11: mean = 12/11, var = 1/11,
        //   S = 2 - 12/11 - 0.5*sqrt(1/11) = 0.7583352368020273.
        // At n = 12: mean = 14/12, var = 2/11,
        //   S = S11 + (2 - 14/12) - 0.5*sqrt(2/11) = 1.39704383409498.
        assert!(!cusum.update(2.0));
        close(cusum.statistic(), 0.7583352368020273, 1e-12);
        assert!(!cusum.update(2.0));
        close(cusum.statistic(), 1.39704383409498, 1e-12);
        // The shift persists, so the detector must fire on the very next
        // point (n = 13, S = 1.947010098498992 > 4*std) and stay alarmed.
        assert!(cusum.update(2.0));
        close(cusum.statistic(), 1.947010098498992, 1e-12);
        assert!(cusum.status());
        assert_eq!(cusum.count(), 13);
    }

    #[test]
    fn reset_clears_state_but_keeps_tuning() {
        let mut cusum = Cusum::default();
        for _ in 0..10 {
            cusum.update(1.0);
        }
        for _ in 0..10 {
            cusum.update(2.0);
        }
        assert!(cusum.status());
        cusum.reset();
        assert_eq!(cusum.count(), 0);
        assert_eq!(cusum.mean(), 0.0);
        assert_eq!(cusum.variance(), 0.0);
        assert_eq!(cusum.statistic(), 0.0);
        assert!(!cusum.status());
        assert_eq!(cusum, Cusum::default());
    }

    #[test]
    fn non_finite_inputs_leave_state_untouched() {
        let mut cusum = Cusum::default();
        cusum.update(1.0);
        let snapshot = cusum.clone();
        assert_eq!(cusum.update(f64::NAN), snapshot.status());
        assert_eq!(cusum.update(f64::INFINITY), snapshot.status());
        assert_eq!(cusum, snapshot);
    }
}
