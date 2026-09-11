//! Prescribed-reactivity insertions in Δk.
//!
//! Taxonomy (take from the upstream PyRK `reactivity_insertion` module):
//! constant (the default no-insertion object), step (Heaviside),
//! impulse/pulse (finite-width box), and ramp (linear rise between two
//! times, then a final level). `Polyline` (piecewise-linear knots) is a
//! Nucleide extension for table-driven schedules such as the
//! `fixtures/kinetics/ramp_table.csv` input; it reduces to the ramp form on
//! two knots. Left behind: everything reactivity-feedback, decay-heat, and
//! input-plumbing related.
//!
//! Conventions mirror the upstream classes: step/impulse/ramp comparisons
//! are right-continuous at knots (a step at `t_step` reads `rho_final` for
//! `t >= t_step`), and the ramp interpolates linearly between its start and
//! end with an independent post-ramp level `rho_final`.

use crate::error::Error;

/// Prescribed reactivity schedule `rho(t)` \[Δk\].
#[derive(Debug, Clone, PartialEq)]
pub enum Reactivity {
    /// `rho(t) = rho` (the no-insertion default).
    Constant {
        /// Reactivity \[Δk\].
        rho: f64,
    },
    /// Heaviside step: `rho_init` for `t < t_step`, else `rho_final`.
    Step {
        /// Step time \[s\].
        t_step: f64,
        /// Reactivity before the step \[Δk\].
        rho_init: f64,
        /// Reactivity at/after the step \[Δk\].
        rho_final: f64,
    },
    /// Finite-width pulse: `rho_max` on `[t_start, t_end]`, else `rho_init`.
    Impulse {
        /// Pulse start \[s\].
        t_start: f64,
        /// Pulse end \[s\].
        t_end: f64,
        /// Baseline reactivity \[Δk\].
        rho_init: f64,
        /// Pulse reactivity \[Δk\].
        rho_max: f64,
    },
    /// Linear rise from `rho_init` at `t_start` toward `rho_rise` at `t_end`,
    /// then `rho_final` afterwards (the post-ramp level is independent, as
    /// in the upstream ramp class).
    Ramp {
        /// Ramp start \[s\].
        t_start: f64,
        /// Ramp end \[s\].
        t_end: f64,
        /// Reactivity before the ramp \[Δk\].
        rho_init: f64,
        /// Reactivity reached at `t_end` along the ramp \[Δk\].
        rho_rise: f64,
        /// Reactivity after the ramp \[Δk\].
        rho_final: f64,
    },
    /// Piecewise-linear interpolation through `(times, values)` knots;
    /// constant extension outside the knot span. Nucleide extension (the
    /// upstream taxonomy has no table-driven insertion).
    Polyline {
        /// Knot times \[s\], strictly increasing, length `>= 2`.
        times: Vec<f64>,
        /// Knot reactivities \[Δk\], same length as `times`.
        values: Vec<f64>,
    },
}

fn finite(x: f64) -> bool {
    x.is_finite()
}

impl Reactivity {
    /// Zero-reactivity (no-insertion) schedule.
    pub fn zero() -> Self {
        Self::Constant { rho: 0.0 }
    }

    /// Validate schedule parameters.
    ///
    /// Rejects non-finite reactivities/times, negative times, degenerate
    /// intervals (`t_end <= t_start`, `t_step < 0`), and polylines with
    /// fewer than two knots, length mismatch, or non-increasing times.
    /// Reactivity magnitudes are *not* bounded (super-prompt-critical steps
    /// are representable; [`crate::prompt_jump`] reports them as errors
    /// only where its formula has no solution).
    pub fn validate(&self) -> Result<(), Error> {
        match self {
            Self::Constant { rho } => {
                if !finite(*rho) {
                    return Err(Error::BadReactivity("rho must be finite"));
                }
            }
            Self::Step {
                t_step,
                rho_init,
                rho_final,
            } => {
                if !finite(*t_step) || *t_step < 0.0 {
                    return Err(Error::BadReactivity("t_step must be finite and >= 0"));
                }
                if !finite(*rho_init) || !finite(*rho_final) {
                    return Err(Error::BadReactivity("step levels must be finite"));
                }
            }
            Self::Impulse {
                t_start,
                t_end,
                rho_init,
                rho_max,
            } => {
                if !finite(*t_start) || !finite(*t_end) || *t_start < 0.0 || *t_end <= *t_start {
                    return Err(Error::BadReactivity("impulse needs 0 <= t_start < t_end"));
                }
                if !finite(*rho_init) || !finite(*rho_max) {
                    return Err(Error::BadReactivity("impulse levels must be finite"));
                }
            }
            Self::Ramp {
                t_start,
                t_end,
                rho_init,
                rho_rise,
                rho_final,
            } => {
                if !finite(*t_start) || !finite(*t_end) || *t_start < 0.0 || *t_end <= *t_start {
                    return Err(Error::BadReactivity("ramp needs 0 <= t_start < t_end"));
                }
                if !finite(*rho_init) || !finite(*rho_rise) || !finite(*rho_final) {
                    return Err(Error::BadReactivity("ramp levels must be finite"));
                }
            }
            Self::Polyline { times, values } => {
                if times.len() < 2 || times.len() != values.len() {
                    return Err(Error::BadReactivity(
                        "polyline needs >= 2 knots with matching values",
                    ));
                }
                let mut prev = f64::NEG_INFINITY;
                for (t, v) in times.iter().zip(values) {
                    if !finite(*t) || !finite(*v) {
                        return Err(Error::BadReactivity("polyline knots must be finite"));
                    }
                    if *t <= prev {
                        return Err(Error::BadReactivity(
                            "polyline times must be strictly increasing",
                        ));
                    }
                    prev = *t;
                }
                if times[0] < 0.0 {
                    return Err(Error::BadReactivity("polyline times must be >= 0"));
                }
            }
        }
        Ok(())
    }

    /// Evaluate `rho(t)` \[Δk\] with right-continuous knot semantics.
    ///
    /// # Panics
    ///
    /// Panics when the schedule holds no evaluable knots (e.g. an
    /// unvalidated empty `Polyline` indexes an empty vector). Validate once
    /// up front with [`Self::validate`] ([`crate::solve()`] already does this).
    pub fn eval(&self, t: f64) -> f64 {
        match self {
            Self::Constant { rho } => *rho,
            Self::Step {
                t_step,
                rho_init,
                rho_final,
            } => {
                if t < *t_step {
                    *rho_init
                } else {
                    *rho_final
                }
            }
            Self::Impulse {
                t_start,
                t_end,
                rho_init,
                rho_max,
            } => {
                if t < *t_start || t > *t_end {
                    *rho_init
                } else {
                    *rho_max
                }
            }
            Self::Ramp {
                t_start,
                t_end,
                rho_init,
                rho_rise,
                rho_final,
            } => {
                if t < *t_start {
                    *rho_init
                } else if t <= *t_end {
                    rho_init + (rho_rise - rho_init) * (t - t_start) / (t_end - t_start)
                } else {
                    *rho_final
                }
            }
            Self::Polyline { times, values } => {
                if t <= times[0] {
                    return values[0];
                }
                if t >= times[times.len() - 1] {
                    return values[values.len() - 1];
                }
                let k = times.partition_point(|tk| *tk <= t);
                let t0 = times[k - 1];
                let t1 = times[k];
                let v0 = values[k - 1];
                let v1 = values[k];
                v0 + (v1 - v0) * (t - t0) / (t1 - t0)
            }
        }
    }

    /// Knot times where `rho(t)` is discontinuous or kinked.
    ///
    /// The solver steps exactly onto these so Heaviside edges and
    /// ramp/polyline kinks never fall inside an integration step.
    pub fn knots(&self) -> Vec<f64> {
        match self {
            Self::Constant { .. } => vec![],
            Self::Step { t_step, .. } => vec![*t_step],
            Self::Impulse { t_start, t_end, .. } => vec![*t_start, *t_end],
            Self::Ramp { t_start, t_end, .. } => vec![*t_start, *t_end],
            Self::Polyline { times, .. } => times.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constant_zero() {
        let r = Reactivity::zero();
        r.validate().unwrap();
        assert_eq!(r.eval(0.0), 0.0);
        assert_eq!(r.eval(1e3), 0.0);
        assert!(r.knots().is_empty());
    }

    #[test]
    fn step_is_right_continuous() {
        // Upstream StepReactivityInsertion: t < t_step -> init, else final.
        let r = Reactivity::Step {
            t_step: 1.0,
            rho_init: 0.0,
            rho_final: 0.002,
        };
        r.validate().unwrap();
        assert_eq!(r.eval(0.999), 0.0);
        assert_eq!(r.eval(1.0), 0.002);
        assert_eq!(r.eval(1.001), 0.002);
        assert_eq!(r.knots(), vec![1.0]);
    }

    #[test]
    fn impulse_box() {
        // Upstream ImpulseReactivityInsertion: max on [t_start, t_end].
        let r = Reactivity::Impulse {
            t_start: 1.0,
            t_end: 2.0,
            rho_init: 0.0,
            rho_max: 0.003,
        };
        r.validate().unwrap();
        assert_eq!(r.eval(0.5), 0.0);
        assert_eq!(r.eval(1.0), 0.003);
        assert_eq!(r.eval(2.0), 0.003);
        assert_eq!(r.eval(2.5), 0.0);
    }

    #[test]
    fn ramp_slope_and_final() {
        // Upstream RampReactivityInsertion: linear between, rho_final after.
        let r = Reactivity::Ramp {
            t_start: 1.0,
            t_end: 3.0,
            rho_init: 0.0,
            rho_rise: 0.004,
            rho_final: 0.003,
        };
        r.validate().unwrap();
        assert_eq!(r.eval(0.5), 0.0);
        assert!((r.eval(2.0) - 0.002).abs() < 1e-15);
        assert_eq!(r.eval(3.0), 0.004);
        assert_eq!(r.eval(3.5), 0.003);
    }

    #[test]
    fn polyline_interp_and_clamp() {
        let r = Reactivity::Polyline {
            times: vec![0.0, 1.0, 3.0],
            values: vec![0.0, 0.002, 0.002],
        };
        r.validate().unwrap();
        assert_eq!(r.eval(-1.0), 0.0);
        assert!((r.eval(0.5) - 0.001).abs() < 1e-15);
        assert_eq!(r.eval(2.0), 0.002);
        assert_eq!(r.eval(99.0), 0.002);
    }

    #[test]
    fn rejects_bad_schedules() {
        assert!(Reactivity::Constant { rho: f64::NAN }.validate().is_err());
        assert!(Reactivity::Step {
            t_step: -1.0,
            rho_init: 0.0,
            rho_final: 0.1
        }
        .validate()
        .is_err());
        assert!(Reactivity::Impulse {
            t_start: 2.0,
            t_end: 2.0,
            rho_init: 0.0,
            rho_max: 0.1
        }
        .validate()
        .is_err());
        assert!(Reactivity::Ramp {
            t_start: 1.0,
            t_end: 1.0,
            rho_init: 0.0,
            rho_rise: 0.1,
            rho_final: 0.1
        }
        .validate()
        .is_err());
        assert!(Reactivity::Polyline {
            times: vec![0.0],
            values: vec![0.0]
        }
        .validate()
        .is_err());
        assert!(Reactivity::Polyline {
            times: vec![0.0, 1.0],
            values: vec![0.0]
        }
        .validate()
        .is_err());
        assert!(Reactivity::Polyline {
            times: vec![1.0, 0.0],
            values: vec![0.0, 0.1]
        }
        .validate()
        .is_err());
    }

    #[test]
    fn rejects_nonfinite_levels_and_knots() {
        assert!(Reactivity::Step {
            t_step: 1.0,
            rho_init: f64::NAN,
            rho_final: 0.1
        }
        .validate()
        .is_err());
        assert!(Reactivity::Impulse {
            t_start: 1.0,
            t_end: 2.0,
            rho_init: 0.0,
            rho_max: f64::INFINITY
        }
        .validate()
        .is_err());
        assert!(Reactivity::Ramp {
            t_start: 1.0,
            t_end: 2.0,
            rho_init: 0.0,
            rho_rise: f64::NAN,
            rho_final: 0.1
        }
        .validate()
        .is_err());
        assert!(Reactivity::Polyline {
            times: vec![0.0, f64::NAN],
            values: vec![0.0, 0.1]
        }
        .validate()
        .is_err());
        assert!(Reactivity::Polyline {
            times: vec![-1.0, 1.0],
            values: vec![0.0, 0.1]
        }
        .validate()
        .is_err());
    }

    #[test]
    fn knots_cover_every_variant() {
        assert!(
            Reactivity::Impulse {
                t_start: 1.0,
                t_end: 2.0,
                rho_init: 0.0,
                rho_max: 0.1
            }
            .knots()
                == vec![1.0, 2.0]
        );
        assert!(
            Reactivity::Ramp {
                t_start: 1.0,
                t_end: 3.0,
                rho_init: 0.0,
                rho_rise: 0.1,
                rho_final: 0.1
            }
            .knots()
                == vec![1.0, 3.0]
        );
        assert!(
            Reactivity::Polyline {
                times: vec![0.0, 1.0],
                values: vec![0.0, 0.1]
            }
            .knots()
                == vec![0.0, 1.0]
        );
    }
}
