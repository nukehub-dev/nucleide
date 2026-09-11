//! Fixture-pinned oracle gates (`fixtures/spectroscopy/`).
//!
//! These tests replay the hand-built synthetic fixtures through the public
//! API; they fail if the implementation drifts from the recorded
//! hand-computed values.

use nucleide_spectroscopy::{calc_bg, gross_count, net_counts};

fn fixture(name: &str) -> serde_json::Value {
    let path = format!(
        "{}/../../fixtures/spectroscopy/{name}",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(&path).expect("fixture readable");
    serde_json::from_str(&text).expect("fixture parses")
}

fn vec_of(v: &serde_json::Value) -> Vec<f64> {
    serde_json::from_value(v.clone()).unwrap()
}

#[test]
fn smooth_oracle_replay() {
    let fix = fixture("smooth_oracle.json");
    let counts = vec_of(&fix["counts"]);

    let rect = &fix["rect_smooth_m5"];
    let m = rect["m"].as_u64().unwrap() as usize;
    let want = vec_of(&rect["expected"]);
    let got = nucleide_spectroscopy::rect_smooth(&counts, m).unwrap();
    assert_eq!(got.len(), want.len());
    for (g, w) in got.iter().zip(want.iter()) {
        assert!((g - w).abs() < 1e-12, "rect {g} vs {w}");
    }

    let want = vec_of(&fix["five_point"]["expected"]);
    let got = nucleide_spectroscopy::five_point_smooth(&counts).unwrap();
    assert_eq!(got.len(), want.len());
    for (g, w) in got.iter().zip(want.iter()) {
        assert!((g - w).abs() < 1e-12, "five-point {g} vs {w}");
    }
}

#[test]
fn counts_oracle_replay() {
    let fix = fixture("smooth_oracle.json");
    let counts = vec_of(&fix["counts"]);
    let channels: Vec<f64> = (0..counts.len()).map(|c| c as f64).collect();
    let c = &fix["counts_oracle"];
    let (c1, c2, m) = (
        c["c1"].as_i64().unwrap(),
        c["c2"].as_i64().unwrap(),
        c["m"].as_i64().unwrap(),
    );
    let rel = |g: f64, w: f64| (g - w).abs() / w.abs().max(1e-30);
    assert!(
        rel(
            calc_bg(&counts, &channels, c1, c2, m).unwrap(),
            c["bg"].as_f64().unwrap()
        ) < 1e-12
    );
    assert!(
        rel(
            gross_count(&counts, &channels, c1, c2).unwrap(),
            c["gross"].as_f64().unwrap()
        ) < 1e-12
    );
    assert!(
        rel(
            net_counts(&counts, &channels, c1, c2, m).unwrap(),
            c["net"].as_f64().unwrap()
        ) < 1e-12
    );
}

#[test]
fn spe_fixtures_cross_format_counts_match() {
    let dir = format!("{}/../../fixtures/spectroscopy", env!("CARGO_MANIFEST_DIR"));
    let dollar = std::fs::read_to_string(format!("{dir}/dollar_min.spe")).unwrap();
    let plain = std::fs::read_to_string(format!("{dir}/plain_min.spe")).unwrap();
    let d = nucleide_spectroscopy::parse_dollar_spe(&dollar, "dollar_min.spe").unwrap();
    let p = nucleide_spectroscopy::parse_plain_spe(&plain, "plain_min.spe").unwrap();
    assert_eq!(d.spectrum.counts, p.spectrum.counts);
    assert_eq!(d.spectrum.counts.len(), d.spectrum.num_channels);
    assert_eq!(p.spectrum.counts.len(), p.spectrum.num_channels);
}
