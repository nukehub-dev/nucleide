"""Python-side tests for the point-kinetics core (O1-O4 gates + input errors)."""

import json
import math
from pathlib import Path

import pytest

import nucleide.kinetics as kin

FIX = Path(__file__).parent.parent / "fixtures" / "kinetics"

BETAS = [0.00021, 0.00141, 0.00127, 0.00255, 0.00074, 0.00032]
LAMBDAS = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0]
LAMBDA = 1e-5
BETA = sum(BETAS)


def test_o1_initial_rate() -> None:
    oracle = json.loads((FIX / "step_oracle.json").read_text())
    rate = kin.initial_rate([0.0065], [0.08], 1e-4, {"kind": "constant", "rho": 0.002}, 1.0)
    assert rate == pytest.approx(oracle["initial_rate"], rel=1e-12)


def test_o2_prompt_jump_plateau() -> None:
    oracle = json.loads((FIX / "step_oracle.json").read_text())
    pj = kin.prompt_jump(1.0, 0.0, 0.002, 0.0065)
    assert pj == pytest.approx(oracle["prompt_jump_factor"], rel=1e-12)
    # 6-group plateau: prompt transient (Lambda/(beta-rho) ~ 3 ms) is over
    # at 0.05 s while precursors are still frozen.
    out = kin.solve(
        BETAS,
        LAMBDAS,
        LAMBDA,
        {"kind": "step", "t_step": 1.0, "rho_init": 0.0, "rho_final": 0.003},
        [1.05],
        1.0,
        rtol=1e-10,
        atol=1e-14,
    )
    want = kin.prompt_jump(1.0, 0.0, 0.003, BETA)
    assert out["n"][0] == pytest.approx(want, rel=0.02)


def test_o3_one_group_closed_form() -> None:
    oracle = json.loads((FIX / "step_oracle.json").read_text())
    step = oracle["step"]
    ts = [step["t_step"] + t for t in oracle["analytic_since_step"]["t"]]
    out = kin.solve(
        oracle["params"]["betas"],
        oracle["params"]["lambdas"],
        oracle["params"]["Lambda"],
        {"kind": "step", **step},
        ts,
        oracle["n0"],
        rtol=1e-10,
        atol=1e-14,
    )
    for n, w in zip(out["n"], oracle["analytic_since_step"]["n"], strict=True):
        assert n == pytest.approx(w, rel=1e-6)


def test_o4_tail_matches_stable_period() -> None:
    check = json.loads((FIX / "inhour_check.json").read_text())
    period = kin.stable_period(
        check["params"]["betas"],
        check["params"]["lambdas"],
        check["params"]["Lambda"],
        check["rho"],
    )
    assert period == pytest.approx(check["stable_period_s"], rel=1e-9)
    out = kin.solve(
        check["params"]["betas"],
        check["params"]["lambdas"],
        check["params"]["Lambda"],
        {"kind": "constant", "rho": check["rho"]},
        [50.0, 100.0],
        1.0,
        rtol=1e-10,
        atol=1e-14,
    )
    slope = math.log(out["n"][1] / out["n"][0]) / 50.0
    assert slope == pytest.approx(1.0 / period, rel=0.02)


def test_equilibrium_and_inhour_helpers() -> None:
    c0 = kin.equilibrium(BETAS, LAMBDAS, LAMBDA, 2.0)
    assert len(c0) == 6
    assert c0[0] == pytest.approx(0.00021 / (0.01 * LAMBDA) * 2.0, rel=1e-12)
    assert kin.inhour_rho(BETAS, LAMBDAS, LAMBDA, 0.0) == pytest.approx(0.0, abs=0.0)
    rho = kin.inhour_rho(BETAS, LAMBDAS, LAMBDA, 0.05605439092656976)
    assert rho == pytest.approx(0.002, rel=1e-9)


def test_insertion_kinds_and_errors() -> None:
    n0 = kin.solve(BETAS, LAMBDAS, LAMBDA, {"kind": "constant", "rho": 0.0}, [1.0], 1.0)
    assert n0["n"][0] == pytest.approx(1.0, rel=1e-9)
    imp = kin.solve(
        BETAS,
        LAMBDAS,
        LAMBDA,
        {"kind": "impulse", "t_start": 1.0, "t_end": 2.0, "rho_init": 0.0, "rho_max": 0.002},
        [0.5, 1.5, 3.0],
        1.0,
    )
    assert imp["n"][0] == pytest.approx(1.0, rel=1e-6)
    assert imp["n"][1] > 1.0
    ramp = kin.solve(
        BETAS,
        LAMBDAS,
        LAMBDA,
        {
            "kind": "ramp",
            "t_start": 1.0,
            "t_end": 3.0,
            "rho_init": 0.0,
            "rho_rise": 0.002,
            "rho_final": 0.002,
        },
        [0.5, 2.0, 5.0],
        1.0,
    )
    assert ramp["n"][0] == pytest.approx(1.0, rel=1e-6)
    assert ramp["n"][2] > ramp["n"][1] > 1.0
    with pytest.raises(ValueError):
        kin.solve(BETAS, LAMBDAS, LAMBDA, {"kind": "nope"}, [1.0], 1.0)
    with pytest.raises(ValueError):
        kin.solve(BETAS, LAMBDAS, LAMBDA, {"kind": "constant", "rho": 0.0}, [1.0, 1.0], 1.0)
    with pytest.raises(ValueError):
        kin.solve([0.0065], [0.08], 0.0, {"kind": "constant", "rho": 0.0}, [1.0], 1.0)
    with pytest.raises(ValueError):
        kin.prompt_jump(1.0, 0.0, BETA, BETA)
    with pytest.raises(ValueError):
        kin.stable_period(BETAS, LAMBDAS, LAMBDA, 0.0)
    with pytest.raises(ValueError):
        kin.solve(
            BETAS, LAMBDAS, LAMBDA, {"kind": "constant", "rho": 0.0}, [1.0], 1.0, method="rk4"
        )
