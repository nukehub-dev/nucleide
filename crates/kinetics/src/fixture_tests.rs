//! Fixture-pinned oracle gates (`fixtures/kinetics/`).
//!
//! These tests replay the hand-built synthetic fixtures through the public
//! API; they fail if the implementation drifts from the recorded
//! closed-form values.

use crate::{KineticParams, Reactivity, Solution, SolverOptions, State, TimeGrid};

fn fixture(name: &str) -> serde_json::Value {
    let path = format!(
        "{}/../../fixtures/kinetics/{name}",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(&path).expect("fixture readable");
    serde_json::from_str(&text).expect("fixture parses")
}

fn params_of(v: &serde_json::Value) -> KineticParams {
    let p = &v["params"];
    let betas: Vec<f64> = serde_json::from_value(p["betas"].clone()).unwrap();
    let lambdas: Vec<f64> = serde_json::from_value(p["lambdas"].clone()).unwrap();
    let gen: f64 = serde_json::from_value(p["Lambda"].clone()).unwrap();
    KineticParams::new(betas, lambdas, gen).unwrap()
}

fn tight() -> SolverOptions {
    SolverOptions {
        rtol: 1e-10,
        atol: 1e-14,
        ..Default::default()
    }
}

#[test]
fn step_oracle_initial_rate_and_jump_factor() {
    let v = fixture("step_oracle.json");
    let p = params_of(&v);
    // Algebraic gates @1e-12 against the recorded closed-form values.
    let rc = Reactivity::Constant {
        rho: v["step"]["rho_final"].as_f64().unwrap(),
    };
    let s = State::new(&p, v["n0"].as_f64().unwrap(), None).unwrap();
    let rate = crate::solve::initial_rate(&p, &rc, &s);
    let want_rate = v["initial_rate"].as_f64().unwrap();
    assert!((rate - want_rate).abs() / want_rate < 1e-12);
    let pj = crate::prompt_jump(
        1.0,
        v["step"]["rho_init"].as_f64().unwrap(),
        v["step"]["rho_final"].as_f64().unwrap(),
        p.beta_total(),
    )
    .unwrap();
    let want_pj = v["prompt_jump_factor"].as_f64().unwrap();
    assert!((pj - want_pj).abs() / want_pj < 1e-12);
}

#[test]
fn step_oracle_transient_matches_closed_form() {
    let v = fixture("step_oracle.json");
    let p = params_of(&v);
    let rho = Reactivity::Step {
        t_step: v["step"]["t_step"].as_f64().unwrap(),
        rho_init: v["step"]["rho_init"].as_f64().unwrap(),
        rho_final: v["step"]["rho_final"].as_f64().unwrap(),
    };
    let n0 = v["n0"].as_f64().unwrap();
    let s = State::new(&p, n0, None).unwrap();
    let ts: Vec<f64> = serde_json::from_value(v["analytic_since_step"]["t"].clone()).unwrap();
    let want: Vec<f64> = serde_json::from_value(v["analytic_since_step"]["n"].clone()).unwrap();
    let t_step = v["step"]["t_step"].as_f64().unwrap();
    let grid_times: Vec<f64> = ts.iter().map(|t| t_step + t).collect();
    let sol: Solution =
        crate::solve(&p, &rho, &TimeGrid::new(grid_times).unwrap(), &s, &tight()).unwrap();
    for (n, w) in sol.n.iter().zip(&want) {
        assert!((n - w).abs() / w < 1e-6, "{n} vs {w}");
    }
}

#[test]
fn inhour_fixture_period_and_residual() {
    let v = fixture("inhour_check.json");
    let p = params_of(&v);
    let rho = v["rho"].as_f64().unwrap();
    let t = crate::stable_period(&p, rho).unwrap();
    let want = v["stable_period_s"].as_f64().unwrap();
    assert!((t - want).abs() / want < 1e-9, "{t} vs {want}");
    let res = crate::inhour_residual(&p, rho, 1.0 / t).unwrap().abs();
    assert!(res < 1e-10, "residual {res}");
}

#[test]
fn ramp_table_polyline_stays_physical() {
    // Input-only fixture: identical history feeds the Rust solver here and
    // the PyRK cross-check in validation/.
    let path = format!(
        "{}/../../fixtures/kinetics/ramp_table.csv",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(&path).expect("ramp table readable");
    let mut times = vec![];
    let mut values = vec![];
    for (k, line) in text.lines().enumerate() {
        if k == 0 {
            assert_eq!(line.trim(), "t_s,rho_dk");
            continue;
        }
        let (a, b) = line.split_once(',').expect("two columns");
        times.push(a.trim().parse::<f64>().unwrap());
        values.push(b.trim().parse::<f64>().unwrap());
    }
    assert!(times.len() >= 10);
    let p = KineticParams::new(
        vec![0.00021, 0.00141, 0.00127, 0.00255, 0.00074, 0.00032],
        vec![0.01, 0.03, 0.1, 0.3, 1.0, 3.0],
        1e-5,
    )
    .unwrap();
    let rho = Reactivity::Polyline { times, values };
    rho.validate().unwrap();
    let s = State::new(&p, 1.0, None).unwrap();
    let grid_times: Vec<f64> = (1..=40).map(|k| 0.25 * k as f64).collect();
    let sol = crate::solve(&p, &rho, &TimeGrid::new(grid_times).unwrap(), &s, &tight()).unwrap();
    for (n, c) in sol.n.iter().zip(&sol.c) {
        assert!(*n > 0.0);
        assert!(c.iter().all(|x| *x > 0.0));
    }
    for w in sol.n.windows(2) {
        // Tolerance-aware monotonicity: stationary segments (rho = 0) hold
        // only to the stepper tolerance, so admit sub-tolerance dips.
        assert!(w[1] >= w[0] - 1e-9, "ramp-up-and-hold must not decrease");
    }
    assert!(*sol.n.last().unwrap() > 1.0);
}
