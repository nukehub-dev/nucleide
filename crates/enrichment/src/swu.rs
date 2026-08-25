//! Separative work unit (SWU) analytics: the value (separation potential)
//! function and per-stream SWU ratios.
//!
//! Closed-form separative-work helpers.
//! (`value_func`, `swu_per_feed`, `swu_per_prod`, `swu_per_tail`). The value
//! function is the standard Dirac separation potential
//!
//! ```text
//! V(x) = (2x - 1) ln(x / (1 - x))
//! ```
//!
//! (upstream's docstring writes `log(x/(x-1))` but its code — and its test
//! suite — use `log(x/(1-x))`; the code convention is ported here).

use crate::cascade::{
    feed_per_prod, feed_per_tail, prod_per_feed, prod_per_tail, tail_per_feed, tail_per_prod,
};

/// Value or separation potential of an assay `x`
/// (separative value function): `V(x) = (2x - 1) ln(x/(1 - x))`.
///
/// Zero at `x = 0.5`, positive elsewhere on `(0, 1)`, and symmetric under
/// `x -> 1 - x`: `V(1-x) = (1-2x) ln((1-x)/x) = V(x)`.
pub fn value_func(x: f64) -> f64 {
    (2.0 * x - 1.0) * (x / (1.0 - x)).ln()
}

/// SWU per unit mass of feed material for assays `x_feed`, `x_prod`,
/// `x_tail`:
///
/// ```text
/// P/F * V(xP) + T/F * V(xT) - V(xF)
/// ```
pub fn swu_per_feed(x_feed: f64, x_prod: f64, x_tail: f64) -> f64 {
    prod_per_feed(x_feed, x_prod, x_tail) * value_func(x_prod)
        + tail_per_feed(x_feed, x_prod, x_tail) * value_func(x_tail)
        - value_func(x_feed)
}

/// SWU per unit mass of product material:
///
/// ```text
/// V(xP) + T/P * V(xT) - F/P * V(xF)
/// ```
pub fn swu_per_prod(x_feed: f64, x_prod: f64, x_tail: f64) -> f64 {
    value_func(x_prod) + tail_per_prod(x_feed, x_prod, x_tail) * value_func(x_tail)
        - feed_per_prod(x_feed, x_prod, x_tail) * value_func(x_feed)
}

/// SWU per unit mass of tails material:
///
/// ```text
/// P/T * V(xP) + V(xT) - F/T * V(xF)
/// ```
pub fn swu_per_tail(x_feed: f64, x_prod: f64, x_tail: f64) -> f64 {
    prod_per_tail(x_feed, x_prod, x_tail) * value_func(x_prod) + value_func(x_tail)
        - feed_per_tail(x_feed, x_prod, x_tail) * value_func(x_feed)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// pytest.approx defaults: absolute 1e-12 OR relative 1e-6.
    fn approx(obs: f64, exp: f64) -> bool {
        (obs - exp).abs() <= 1e-12 + 1e-6 * exp.abs()
    }

    fn assert_rel_close(obs: f64, exp: f64, rel: f64, what: &str) {
        assert!(
            (obs / exp - 1.0).abs() < rel,
            "{what}: observed {obs}, expected {exp} (rel {rel})"
        );
    }

    #[test]
    fn value_func_definition() {
        // Port of test_value.
        let x = 0.0072_f64;
        let exp = (2.0 * x - 1.0) * (x / (1.0 - x)).ln();
        let obs = value_func(x);
        assert!(approx(obs, exp), "value_func({x}): {obs} vs {exp}");
    }

    #[test]
    fn value_func_symmetry_and_zero() {
        assert_eq!(value_func(0.5), 0.0);
        let x = 0.11_f64;
        assert!((value_func(x) - value_func(1.0 - x)).abs() < 1e-12);
        assert!(value_func(0.75) > 0.0 && value_func(0.25) > 0.0);
    }

    #[test]
    fn swu_from_all_streams_agree() {
        // Port of test_swu: 15.1596 kg feed <-> 1.5 kg product <->
        // 13.6596 kg tails all require ~11.765 kg SWU.
        let (xf, xp, xt) = (0.0072_f64, 0.05_f64, 0.0025_f64);
        let (feed, prod, tails) = (15.1596_f64, 1.5_f64, 13.6596_f64);
        let exp = 11765.0 / 1e3;

        assert_rel_close(feed * swu_per_feed(xf, xp, xt), exp, 1e-4, "SWU from feed");
        assert_rel_close(
            prod * swu_per_prod(xf, xp, xt),
            exp,
            1e-4,
            "SWU from product",
        );
        assert_rel_close(
            tails * swu_per_tail(xf, xp, xt),
            exp,
            1e-4,
            "SWU from tails",
        );
    }

    #[test]
    fn swu_per_stream_identities() {
        // The three ratios are consistent views of one separative work:
        // SWU/F = (P/F) * SWU/P and SWU/F = (T/F) * SWU/T.
        let (xf, xp, xt) = (0.0072_f64, 0.05_f64, 0.0025_f64);
        let ppf = crate::prod_per_feed(xf, xp, xt);
        let tpf = crate::tail_per_feed(xf, xp, xt);
        let swuf = swu_per_feed(xf, xp, xt);
        assert!((swuf - ppf * swu_per_prod(xf, xp, xt)).abs() < 1e-12);
        assert!((swuf - tpf * swu_per_tail(xf, xp, xt)).abs() < 1e-12);
        // Enrichment costs work: positive for a realistic assay ladder.
        assert!(swuf > 0.0);
    }
}
