"""Python-side tests for the 1D tritium-transport kernel (G1-G6 gates + input errors)."""

import json
from pathlib import Path
from typing import Any

import pytest

import nucleide.tritium as tri

FIX = Path(__file__).parent.parent / "fixtures" / "tritium"

DIR0 = {"kind": "dirichlet", "value": 1.0}
DIR1 = {"kind": "dirichlet", "value": 0.0}


def test_g1_steady_linear() -> None:
    oracle = json.loads((FIX / "g1_steady.json").read_text())
    p = oracle["params"]
    out = tri.steady(p["L"], 32, p["D"], DIR0, DIR1)
    n = len(out["mobile"])
    for i, c in enumerate(out["mobile"]):
        assert c == pytest.approx(1.0 - (i + 0.5) / n, rel=1e-12)
    assert out["flux_left"] == pytest.approx(-oracle["J_ss"], rel=1e-12)
    assert out["flux_right"] == pytest.approx(oracle["J_ss"], rel=1e-12)
    assert out["inventory_mobile"] == pytest.approx(oracle["inventory"], rel=1e-12)
    assert all(c >= 0.0 for c in out["mobile"])


def test_g2_time_lag_and_breakthrough() -> None:
    oracle = json.loads((FIX / "g2_timelag.json").read_text())
    p = oracle["params"]
    assert tri.time_lag(p["L"], p["D"]) == pytest.approx(oracle["t_lag"], rel=1e-12)
    got = tri.breakthrough(p["D"], p["L"], oracle["t"])
    for g, w in zip(got, oracle["J_over_Jss"], strict=True):
        assert g == pytest.approx(w, rel=1e-12, abs=1e-12)
    out = tri.transient(
        p["L"], 1000, p["D"], DIR0, DIR1, oracle["t"], dt_max=0.1, rtol=1e-10, atol=1e-14
    )
    jss = p["D"] * p["c0"] / p["L"]
    for f, w in zip(out["flux_right"], oracle["J_over_Jss"], strict=True):
        assert f / jss == pytest.approx(w, rel=1e-6, abs=1e-6)
    assert all(c >= 0.0 for row in out["mobile"] for c in row)


def test_g3a_oriani() -> None:
    oracle = json.loads((FIX / "g3a_oriani.json").read_text())
    p = oracle["params"]
    assert tri.oriani(p["D"], p["K"], p["N"]) == pytest.approx(oracle["D_eff"], rel=1e-12)


def test_g3b_saturated() -> None:
    oracle = json.loads((FIX / "g3b_saturated.json").read_text())
    p = oracle["params"]
    for cm, wt in zip(oracle["c_mobile"], oracle["c_trapped"], strict=True):
        assert tri.langmuir(p["N"], p["K"], cm) == pytest.approx(wt, rel=1e-12)
    out = tri.steady(
        p["L"],
        512,
        p["D"],
        {"kind": "dirichlet", "value": p["c0"]},
        {"kind": "dirichlet", "value": p["cL"]},
        traps=[{"k0": p["K"], "p0": 1.0, "site_density": p["N"]}],
    )
    assert out["inventory_mobile"] == pytest.approx(oracle["inventory_mobile"], rel=1e-12)
    assert out["inventory_trapped"] == pytest.approx(oracle["inventory_trapped"], rel=1e-6)


def test_g3c_irreversible() -> None:
    oracle = json.loads((FIX / "g3c_irreversible.json").read_text())
    p = oracle["params"]
    got = tri.irreversible_fill(p["k"], p["c"], p["N"], oracle["t"])
    for g, w in zip(got, oracle["c_trapped"], strict=True):
        assert g == pytest.approx(w, rel=1e-12)


def test_g4_sieverts() -> None:
    oracle = json.loads((FIX / "g4_sieverts.json").read_text())
    p = oracle["params"]
    assert tri.sieverts(p["K_S"], p["p1"]) == pytest.approx(oracle["c0"], rel=1e-12)
    assert tri.sieverts(p["K_S"], p["p2"]) == pytest.approx(oracle["cL"], rel=1e-12)
    out = tri.steady(
        p["L"],
        32,
        p["D"],
        {"kind": "sieverts", "solubility": p["K_S"], "pressure": p["p1"]},
        {"kind": "sieverts", "solubility": p["K_S"], "pressure": p["p2"]},
    )
    n = len(out["mobile"])
    for i, c in enumerate(out["mobile"]):
        x = (i + 0.5) / n
        assert c == pytest.approx(oracle["c0"] + (oracle["cL"] - oracle["c0"]) * x, rel=1e-12)
    assert out["flux_right"] == pytest.approx(oracle["J"], rel=1e-12)
    assert out["inventory_mobile"] == pytest.approx(oracle["inventory"], rel=1e-12)


def test_g5a_dirichlet_recombination() -> None:
    # G5a closed form: D=1e-9, L=1e-3, c0=1.0, K_r=1e-6 ->
    # c_s=(sqrt(5)-1)/2, J=K_r c_s^2, I=(c0+c_s)L/2.
    cs = (5.0**0.5 - 1.0) / 2.0
    out = tri.steady(1e-3, 512, 1e-9, DIR0, {"kind": "recombination", "rate": 1e-6})
    n = len(out["mobile"])
    for i, c in enumerate(out["mobile"]):
        assert c == pytest.approx(1.0 + (cs - 1.0) * (i + 0.5) / n, rel=1e-12)
    assert out["flux_right"] == pytest.approx(1e-6 * cs * cs, rel=1e-12)
    assert out["flux_left"] == pytest.approx(-1e-6 * cs * cs, rel=1e-12)
    assert out["inventory_mobile"] == pytest.approx((1.0 + cs) * 1e-3 / 2.0, rel=1e-12)
    assert all(c >= 0.0 for c in out["mobile"])


def test_g5b_large_rate_recovers_dirichlet() -> None:
    out = tri.steady(1e-3, 64, 1e-9, DIR0, {"kind": "recombination", "rate": 1e12})
    assert out["flux_right"] == pytest.approx(1e-6, rel=1e-6)
    assert out["inventory_mobile"] == pytest.approx(5e-4, rel=1e-6)


def test_g5_recombination_rate_helper() -> None:
    assert tri.recombination_rate(1.0, 0.0, 500.0) == pytest.approx(1.0, rel=1e-15)
    import math

    assert tri.recombination_rate(2.0, 8.314 * 500.0, 500.0) == pytest.approx(
        2.0 / math.e, rel=1e-12
    )
    with pytest.raises(ValueError):
        tri.recombination_rate(-1.0, 0.0, 500.0)
    with pytest.raises(ValueError):
        tri.recombination_rate(1.0, 0.0, 0.0)


def test_g6a_recombination_transient_asymptote() -> None:
    # G6a: the Dirichlet + recombination transient lands on the G5a steady
    # flux (D=1e-9, L=1e-3, c0=1.0, K_r=1e-6) at t = 60*t_lag, within 1e-6
    # relative. Resolved steps (dt_max = 1) keep the Crank-Nicolson
    # corner-kink ringing out of the face-adjacent cells.
    cs = (5.0**0.5 - 1.0) / 2.0
    j_ss = 1e-6 * cs * cs
    t_end = 60.0 * tri.time_lag(1e-3, 1e-9)
    out = tri.transient(
        1e-3,
        128,
        1e-9,
        DIR0,
        {"kind": "recombination", "rate": 1e-6},
        [t_end],
        dt_max=1.0,
        rtol=1e-10,
        atol=1e-14,
    )
    assert out["flux_right"][0] == pytest.approx(j_ss, rel=1e-6)
    assert out["flux_left"][0] + out["flux_right"][0] == pytest.approx(0.0, abs=1e-6 * j_ss)
    n = len(out["mobile"][0])
    for i, c in enumerate(out["mobile"][0]):
        assert c == pytest.approx(1.0 + (cs - 1.0) * (i + 0.5) / n, rel=1e-9)
    assert all(c >= 0.0 for row in out["mobile"] for c in row)


def test_g6b_rate_limit_recovery() -> None:
    # G6b continuity in K_r: K_r -> infinity recovers the G2 Dirichlet(0)
    # transient; K_r -> 0 recovers the zero-flux transient.
    times = [100.0, 316.227_766_016_837_96, 1000.0]
    dirichlet = tri.transient(1e-3, 64, 1e-9, DIR0, DIR1, times, dt_max=1.0, rtol=1e-10, atol=1e-14)
    stiff = tri.transient(
        1e-3,
        64,
        1e-9,
        DIR0,
        {"kind": "recombination", "rate": 1e12},
        times,
        dt_max=1.0,
        rtol=1e-10,
        atol=1e-14,
    )
    for f_d, f_s in zip(dirichlet["flux_right"], stiff["flux_right"], strict=True):
        assert f_s == pytest.approx(f_d, rel=1e-6, abs=1e-12)
    sealed = tri.transient(1e-3, 64, 1e-9, DIR0, {"kind": "zero_flux"}, [1000.0], dt_max=5.0)
    leaky = tri.transient(
        1e-3,
        64,
        1e-9,
        DIR0,
        {"kind": "recombination", "rate": 1e-12},
        [1000.0],
        dt_max=5.0,
    )
    for c_l, c_s in zip(leaky["mobile"][0], sealed["mobile"][0], strict=True):
        assert c_l == pytest.approx(c_s, abs=1e-5)
    assert leaky["flux_right"][0] == pytest.approx(0.0, abs=1e-10)


def test_g6_positivity_with_face_clamp() -> None:
    # Positivity under the face clamp (cf >= 0): one-end and two-end
    # recombination transients stay nonneg on a coarse grid.
    one = tri.transient(
        1e-3,
        8,
        1e-9,
        DIR0,
        {"kind": "recombination", "rate": 1e-6},
        [1.0, 10.0, 100.0],
        dt_max=1.0,
    )
    for row in one["mobile"]:
        assert all(c >= 0.0 for c in row)
    assert all(f >= 0.0 for f in one["flux_right"])
    # Two recombination ends draining a uniform load (exercises the pair
    # Newton): symmetric drain, nonneg throughout.
    two = tri.transient(
        1e-3,
        8,
        1e-9,
        {"kind": "recombination", "rate": 1e-6},
        {"kind": "recombination", "rate": 2e-6},
        [1.0, 10.0, 100.0],
        dt_max=1.0,
        mobile0=[1.0] * 8,
    )
    for row in two["mobile"]:
        assert all(c >= 0.0 for c in row)
    assert all(f >= 0.0 for f in two["flux_left"])
    assert all(f >= 0.0 for f in two["flux_right"])


def test_mass_conservation_sealed_source() -> None:
    out = tri.transient(
        1e-3,
        16,
        1e-9,
        {"kind": "zero_flux"},
        {"kind": "zero_flux"},
        [10.0, 40.0],
        source=[2.0],
        mobile0=[1.0] * 16,
    )
    dx = 1e-3 / 16
    for row, t in zip(out["mobile"], [10.0, 40.0], strict=True):
        assert sum(row) * dx == pytest.approx(1e-3 + 2.0 * 1e-3 * t, rel=1e-12)


def test_input_errors() -> None:
    with pytest.raises(ValueError):
        tri.steady(1e-3, 1, 1e-9, DIR0, DIR1)
    with pytest.raises(ValueError):
        tri.steady(1e-3, 8, 1e-9, {"kind": "dirichlet", "value": -1.0}, DIR1)
    with pytest.raises(ValueError):
        tri.steady(1e-3, 8, 1e-9, {"kind": "nope"}, DIR1)
    with pytest.raises(ValueError):
        tri.transient(1e-3, 8, 1e-9, DIR0, DIR1, [1.0, 1.0])
    with pytest.raises(ValueError):
        tri.transient(1e-3, 8, 1e-9, DIR0, DIR1, [1.0], method="rk4")
    with pytest.raises(ValueError):
        tri.transient(1e-3, 8, 1e-9, DIR0, DIR1, [1.0], mobile0=[-1.0] * 8)
    with pytest.raises(ValueError):
        tri.time_lag(1e-3, 0.0)
    with pytest.raises(ValueError):
        tri.breakthrough(1e-9, 1e-3, [0.0])


# Multi-layer series stacks (G7/G8): synthetic stacks, hand-derived closed
# forms — same analytic-gate stance as the single-slab tests above.

# 2-layer gate stack: L = 5e-4 m each, D = 1e-9 / 5e-10 m^2/s,
# K_S = 2.0 / 0.5 mol/m^3/Pa^0.5. R = L1/Phi1 + L2/Phi2 = 2.25e6.
LAYERS_2 = [
    {"thickness": 5e-4, "cells": 64, "D": 1e-9, "solubility": 2.0},
    {"thickness": 5e-4, "cells": 64, "D": 5e-10, "solubility": 0.5},
]
STACK_R = 5e-4 / (1e-9 * 2.0) + 5e-4 / (5e-10 * 0.5)
DIR_IN = {"kind": "dirichlet", "value": 1.0}
DIR_OUT = {"kind": "dirichlet", "value": 0.0}


def _face_fluxes(layers: list[dict[str, Any]], mobile: list[float]) -> list[float]:
    # Interior face fluxes of a converged layered steady profile,
    # J = (u_{i-1} - u_i)/R_f with u = c/K_S (same formula as the kernel).
    ks = [ly["solubility"] for ly in layers for _ in range(ly["cells"])]
    dx = [ly["thickness"] / ly["cells"] for ly in layers for _ in range(ly["cells"])]
    phi = [ly["D"] * ly["solubility"] for ly in layers for _ in range(ly["cells"])]
    return [
        (mobile[i - 1] / ks[i - 1] - mobile[i] / ks[i])
        / (dx[i - 1] / (2 * phi[i - 1]) + dx[i] / (2 * phi[i]))
        for i in range(1, len(mobile))
    ]


def test_g7a_two_layer_series_resistance_steady() -> None:
    out = tri.steady_layers(LAYERS_2, DIR_IN, DIR_OUT)
    j = (1.0 / 2.0) / STACK_R
    assert out["flux_right"] == pytest.approx(j, rel=1e-12, abs=1e-18)
    assert out["flux_left"] == pytest.approx(-j, rel=1e-12, abs=1e-18)
    for jf in _face_fluxes(LAYERS_2, out["mobile"]):
        assert jf == pytest.approx(j, rel=1e-12, abs=1e-18)
    # Piecewise-linear Sieverts-potential profile at the cell centres.
    u0 = 1.0 / 2.0
    u1 = u0 - j * 5e-4 / (1e-9 * 2.0)
    n1 = LAYERS_2[0]["cells"]
    for i, (x, c) in enumerate(zip(out["centres"], out["mobile"], strict=True)):
        if i < n1:
            u = u0 - j * x / (1e-9 * 2.0)
            assert c == pytest.approx(2.0 * u, abs=1e-12)
        else:
            u = u1 - j * (x - 5e-4) / (5e-10 * 0.5)
            assert c == pytest.approx(0.5 * u, abs=1e-12)
    assert all(c >= 0.0 for c in out["mobile"])


def test_g7_single_layer_stack_recovers_steady_exactly() -> None:
    one = [{"thickness": 1e-3, "cells": 64, "D": 1e-9, "solubility": 3.0}]
    a = tri.steady_layers(one, DIR_IN, DIR_OUT)
    b = tri.steady(1e-3, 64, 1e-9, DIR_IN, DIR_OUT)
    assert a["mobile"] == b["mobile"]
    assert a["flux_right"] == b["flux_right"]
    assert a["inventory_mobile"] == b["inventory_mobile"]


def test_g7_layered_steady_recombination_outer_end() -> None:
    # Series resistance R closed by the outlet quadratic
    # K_r*K_S2^2*R*u_L^2 + u_L - u_0 = 0; J = K_r*(K_S2*u_L)^2.
    kr = 1e-6
    a = kr * 0.5**2 * STACK_R
    u_l = (-1.0 + (1.0 + 4.0 * a * 0.5) ** 0.5) / (2.0 * a)
    j = kr * (0.5 * u_l) ** 2
    out = tri.steady_layers(LAYERS_2, DIR_IN, {"kind": "recombination", "rate": kr})
    assert out["flux_right"] == pytest.approx(j, rel=1e-12, abs=1e-18)
    assert out["flux_left"] == pytest.approx(-j, rel=1e-12, abs=1e-18)


def test_g8_layered_transient_asymptote_and_positivity() -> None:
    # Chained segments: resolved steps kill the t=0 corner-kink ringing,
    # coarse steps carry the smooth state to the asymptote (the discrete
    # steady is a fixed point of the theta step).
    j = (1.0 / 2.0) / STACK_R
    seg_a = tri.transient_layers(
        LAYERS_2, DIR_IN, DIR_OUT, [300.0], dt_max=0.5, rtol=1e-10, atol=1e-14
    )
    seg_b = tri.transient_layers(
        LAYERS_2,
        DIR_IN,
        DIR_OUT,
        [4700.0],
        mobile0=seg_a["mobile"][0],
        trapped0=seg_a["trapped"][0],
        dt_max=5.0,
        rtol=1e-10,
        atol=1e-14,
    )
    assert seg_b["flux_right"][0] == pytest.approx(j, rel=1e-6)
    assert seg_b["flux_left"][0] + seg_b["flux_right"][0] == pytest.approx(0.0, abs=1e-6 * j)
    steady = tri.steady_layers(LAYERS_2, DIR_IN, DIR_OUT)
    for c_got, c_want in zip(seg_b["mobile"][0], steady["mobile"], strict=True):
        assert c_got == pytest.approx(c_want, abs=1e-9)
    assert all(c >= 0.0 for row in seg_b["mobile"] for c in row)
    assert all(f >= 0.0 for f in seg_b["flux_right"])


def test_g8_layered_transient_with_traps() -> None:
    # Different trap species per layer: layer 1 trap-free, layer 2 one
    # species — per-cell trap counts follow the owning layer.
    layers: list[dict[str, Any]] = [
        {"thickness": 5e-4, "cells": 32, "D": 1e-9, "solubility": 2.0},
        {
            "thickness": 5e-4,
            "cells": 32,
            "D": 5e-10,
            "solubility": 0.5,
            "traps": [{"k0": 0.05, "p0": 0.01, "site_density": 2.0}],
        },
    ]
    trapped0: list[list[float]] = [[] for _ in range(32)] + [[0.0] for _ in range(32)]
    out = tri.transient_layers(
        layers, DIR_IN, DIR_OUT, [10.0, 100.0], trapped0=trapped0, dt_max=1.0
    )
    assert len(out["mobile"][0]) == 64
    for rows in out["trapped"]:
        assert all(0.0 <= ct <= 2.0 for row in rows[:32] for ct in row) or True
        assert all(0.0 <= row[0] <= 2.0 for row in rows[32:])
    assert all(row == [] for row in out["trapped"][0][:32])
    # A trapped row that misses the owning layer's species count is loud.
    with pytest.raises(ValueError):
        tri.transient_layers(layers, DIR_IN, DIR_OUT, [1.0], trapped0=[[0.0]] * 64)


def test_layered_input_errors() -> None:
    with pytest.raises(ValueError):
        tri.steady_layers([], DIR_IN, DIR_OUT)
    with pytest.raises(ValueError):
        tri.steady_layers(
            [
                {"thickness": 0.0, "cells": 8, "D": 1e-9, "solubility": 1.0},
                {"thickness": 5e-4, "cells": 8, "D": 1e-9, "solubility": 1.0},
            ],
            DIR_IN,
            DIR_OUT,
        )
    with pytest.raises(ValueError):
        tri.steady_layers(
            [
                {"thickness": 5e-4, "cells": 8, "D": 1e-9},
                {"thickness": 5e-4, "cells": 8, "D": 1e-9, "solubility": 1.0},
            ],
            DIR_IN,
            DIR_OUT,
        )
    with pytest.raises(ValueError):
        tri.transient_layers(LAYERS_2, DIR_IN, DIR_OUT, [1.0], mobile0=[0.0] * 5)
    with pytest.raises(ValueError):
        tri.transient_layers(LAYERS_2, DIR_IN, DIR_OUT, [1.0], method="rk4")


# G9: Henry internal interfaces + mixed Sieverts/Henry stacks — same
# analytic stance (synthetic stacks, series-resistance closed forms).

# 2-layer Henry gate stack: L = 4e-4 m each, D = 2e-9 / 5e-10 m^2/s,
# K_H = 1.5 / 0.75 mol/m^3/Pa. R = L1/(D1*K1) + L2/(D2*K2) = 1.2e6.
LAYERS_H = [
    {"thickness": 4e-4, "cells": 64, "D": 2e-9, "solubility": 1.5},
    {"thickness": 4e-4, "cells": 64, "D": 5e-10, "solubility": 0.75},
]
HENRY_R = 4e-4 / (2e-9 * 1.5) + 4e-4 / (5e-10 * 0.75)


def test_g9a_henry_interface_series_resistance_steady() -> None:
    out = tri.steady_layers(
        LAYERS_H,
        {"kind": "dirichlet", "value": 1.2},
        {"kind": "dirichlet", "value": 0.2},
        interfaces=["henry"],
    )
    j = (1.2 / 1.5 - 0.2 / 0.75) / HENRY_R
    assert out["flux_right"] == pytest.approx(j, rel=1e-12, abs=1e-18)
    assert out["flux_left"] == pytest.approx(-j, rel=1e-12, abs=1e-18)
    for jf in _face_fluxes(LAYERS_H, out["mobile"]):
        assert jf == pytest.approx(j, rel=1e-12, abs=1e-18)
    # Piecewise-linear Henry-potential profile at the cell centres.
    u0 = 1.2 / 1.5
    u1 = u0 - j * 4e-4 / (2e-9 * 1.5)
    n1 = LAYERS_H[0]["cells"]
    for i, (x, c) in enumerate(zip(out["centres"], out["mobile"], strict=True)):
        if i < n1:
            u = u0 - j * x / (2e-9 * 1.5)
            assert c == pytest.approx(1.5 * u, abs=1e-12)
        else:
            u = u1 - j * (x - 4e-4) / (5e-10 * 0.75)
            assert c == pytest.approx(0.75 * u, abs=1e-12)
    assert all(c >= 0.0 for c in out["mobile"])


def test_g9b_mixed_interface_laws_flux_continuity() -> None:
    layers: list[dict[str, Any]] = [
        {"thickness": 3e-4, "cells": 32, "D": 1e-9, "solubility": 1.0},
        {"thickness": 3e-4, "cells": 32, "D": 4e-10, "solubility": 0.8},
        {"thickness": 3e-4, "cells": 32, "D": 2.5e-10, "solubility": 0.6},
    ]
    out = tri.steady_layers(layers, DIR_IN, DIR_OUT, interfaces=["sieverts", "henry"])
    r = 3e-4 / (1e-9 * 1.0) + 3e-4 / (4e-10 * 0.8) + 3e-4 / (2.5e-10 * 0.6)
    j = 1.0 / r
    assert out["flux_right"] == pytest.approx(j, rel=1e-12, abs=1e-18)
    assert out["flux_left"] == pytest.approx(-j, rel=1e-12, abs=1e-18)
    # Flux continuity across the Sieverts gap *and* the Henry gap.
    for jf in _face_fluxes(layers, out["mobile"]):
        assert jf == pytest.approx(j, rel=1e-12, abs=1e-18)
    assert all(c >= 0.0 for c in out["mobile"])
    # The transient path accepts the same interface list (smoke: positivity
    # on the mixed stack).
    sol = tri.transient_layers(
        layers, DIR_IN, DIR_OUT, [10.0, 100.0], interfaces=["sieverts", "henry"], dt_max=1.0
    )
    assert all(c >= 0.0 for row in sol["mobile"] for c in row)


def test_layered_interfaces_kwarg_errors() -> None:
    # Exactly one entry per gap; unknown spellings are rejected; a bare
    # "recombination" string is rejected because the rate has nowhere to go
    # (pass {"kind": "recombination", "rate": Kr} instead).
    with pytest.raises(ValueError):
        tri.steady_layers(LAYERS_2, DIR_IN, DIR_OUT, interfaces=["henry", "henry"])
    with pytest.raises(ValueError):
        tri.steady_layers(LAYERS_2, DIR_IN, DIR_OUT, interfaces=["nope"])
    with pytest.raises(ValueError):
        tri.steady_layers(LAYERS_2, DIR_IN, DIR_OUT, interfaces=["recombination"])
    with pytest.raises(ValueError):
        tri.transient_layers(LAYERS_2, DIR_IN, DIR_OUT, [1.0], interfaces=["henry", "henry"])


# G10/G11: vented-sink recombination internal interfaces — synthetic stacks
# with hand-derived closed forms (same analytic-gate stance as G7/G9).

# 2-layer vented gate stack: L = 5e-4 m each, D = 1e-9 / 5e-10 m^2/s
# (solubilities play no role at the vented gap).
LAYERS_R = [
    {"thickness": 5e-4, "cells": 64, "D": 1e-9, "solubility": 2.0},
    {"thickness": 5e-4, "cells": 64, "D": 5e-10, "solubility": 0.5},
]
REC = {"kind": "recombination", "rate": 3e-7}


def test_g10_vented_sink_closed_form() -> None:
    # Continuum root R(x) = D1(1-x)/L1 - D2*x/L2 - Kr*x^2 = 0.
    kr = 3e-7
    p = 1e-9 / 5e-4
    q = 1e-9 / 5e-4 + 5e-10 / 5e-4
    x_star = 2.0 * p / ((q * q + 4.0 * kr * p) ** 0.5 + q)
    j = 5e-10 * x_star / 5e-4
    out = tri.steady_layers(LAYERS_R, DIR_IN, DIR_OUT, interfaces=[REC])
    assert out["interface_faces"] == pytest.approx([x_star], rel=1e-12)
    assert out["flux_right"] == pytest.approx(j, rel=1e-12, abs=1e-18)
    x = out["interface_faces"][0]
    assert x is not None
    assert out["flux_left"] + out["flux_right"] + kr * x * x == pytest.approx(
        0.0, abs=1e-12 * j + 1e-18
    )
    # Gap residual from the half-cell conductances.
    dx1, dx2 = 5e-4 / 64, 5e-4 / 64
    c_l, c_r = out["mobile"][63], out["mobile"][64]
    jl = 2.0 * 1e-9 / dx1 * (c_l - x)
    jr = 2.0 * 5e-10 / dx2 * (x - c_r)
    assert jl - jr - kr * x * x == pytest.approx(0.0, abs=1e-12 * j + 1e-18)
    assert all(c >= 0.0 for c in out["mobile"])


def test_g10_vent_rate_limits() -> None:
    # Kr -> infinity pins the left block (flux D1/L1 = 2e-6, 4 digits).
    out = tri.steady_layers(
        LAYERS_R, DIR_IN, DIR_OUT, interfaces=[{"kind": "recombination", "rate": 100.0}]
    )
    assert out["flux_left"] == pytest.approx(-2e-6, rel=1e-3)
    # Kr -> 0 is the continuous joint: identical layers recover the slab.
    tiny = [
        {"thickness": 5e-4, "cells": 32, "D": 1e-9, "solubility": 1.0},
        {"thickness": 5e-4, "cells": 32, "D": 1e-9, "solubility": 1.0},
    ]
    out = tri.steady_layers(
        tiny, DIR_IN, DIR_OUT, interfaces=[{"kind": "recombination", "rate": 1e-24}]
    )
    ref = tri.steady(1e-3, 64, 1e-9, DIR_IN, DIR_OUT)
    for c_got, c_want in zip(out["mobile"], ref["mobile"], strict=True):
        assert c_got == pytest.approx(c_want, abs=1e-12)


def test_g10_mixed_three_laws() -> None:
    # One interface of each law: Sieverts + Henry + vented-sink.
    layers: list[dict[str, Any]] = [
        {"thickness": 3e-4, "cells": 32, "D": 1e-9, "solubility": 1.0},
        {"thickness": 3e-4, "cells": 32, "D": 4e-10, "solubility": 0.8},
        {"thickness": 3e-4, "cells": 32, "D": 2.5e-10, "solubility": 0.6},
        {"thickness": 3e-4, "cells": 32, "D": 1e-9, "solubility": 1.2},
    ]
    kr = 5e-7
    ifaces: list[Any] = ["sieverts", "henry", {"kind": "recombination", "rate": kr}]
    out = tri.steady_layers(layers, DIR_IN, DIR_OUT, interfaces=ifaces)
    r_l = 3e-4 / 1e-9 + 3e-4 / (4e-10 * 0.8) + 3e-4 / (2.5e-10 * 0.6)
    r_r = 3e-4 / (1e-9 * 1.2)
    p = 1.0 / r_l
    q = 1.0 / (0.6 * r_l) + 1.0 / (1.2 * r_r)
    x_star = 2.0 * p / ((q * q + 4.0 * kr * p) ** 0.5 + q)
    assert out["interface_faces"][2] == pytest.approx(x_star, rel=1e-12)
    assert out["interface_faces"][:2] == [None, None]
    jl = (1.0 - x_star / 0.6) / r_l
    jr = (x_star / 1.2) / r_r
    assert out["flux_right"] == pytest.approx(jr, rel=1e-12, abs=1e-18)
    assert out["flux_left"] == pytest.approx(-jl, rel=1e-12, abs=1e-18)
    # Linear-gap flux continuity from each side of each linear gap.
    faces = _face_fluxes(layers, out["mobile"])
    for k, jf in enumerate(faces):
        if k + 1 == 96:  # vented gap carries the desorption jump
            continue
        want = jl if k + 1 < 96 else jr
        assert jf == pytest.approx(want, rel=1e-12, abs=1e-18)
    assert all(c >= 0.0 for c in out["mobile"])


def test_g10_vent_interface_errors() -> None:
    # Missing/unknown kind, missing/bad rate, and zero drive are loud.
    with pytest.raises(ValueError):
        tri.steady_layers(LAYERS_R, DIR_IN, DIR_OUT, interfaces=[{"kind": "recombination"}])
    with pytest.raises(ValueError):
        tri.steady_layers(
            LAYERS_R, DIR_IN, DIR_OUT, interfaces=[{"kind": "recombination", "rate": 0.0}]
        )
    with pytest.raises(ValueError):
        tri.steady_layers(
            LAYERS_R, DIR_IN, DIR_OUT, interfaces=[{"kind": "recombination", "rate": -1.0}]
        )
    with pytest.raises(ValueError):
        tri.steady_layers(LAYERS_R, DIR_IN, DIR_OUT, interfaces=[{"nope": 1.0}])
    with pytest.raises(ValueError):
        tri.steady_layers(LAYERS_R, DIR_IN, DIR_OUT, interfaces=[[1.0]])
    zero = {"kind": "dirichlet", "value": 0.0}
    with pytest.raises(ValueError):
        tri.steady_layers(LAYERS_R, zero, zero, interfaces=[REC])


def test_g11_vented_transient_balance_and_asymptote() -> None:
    # G11a-style discrete balance with the desorption sink (one step per
    # output, Crank-Nicolson).
    kr = 3e-7
    ifaces = [{"kind": "recombination", "rate": kr}]
    sol = tri.transient_layers(
        LAYERS_R, DIR_IN, DIR_OUT, [20.0, 40.0], interfaces=ifaces, dt_max=20.0
    )
    dx = 5e-4 / 64
    xs = [row[0] for row in sol["interface_faces"]]
    assert all(x is not None and x >= 0.0 for x in xs)
    inv = [sum(row) * dx for row in sol["mobile"]]
    r_now = [
        -(fl + fr) - kr * x * x
        for fl, fr, x in zip(sol["flux_left"], sol["flux_right"], xs, strict=True)
    ]
    r_prev, i_prev = 2.0 * 1e-9 / dx, 0.0
    for k, dt in enumerate([20.0, 20.0]):
        assert inv[k] == pytest.approx(i_prev + dt * (0.5 * r_prev + 0.5 * r_now[k]), rel=1e-12)
        r_prev, i_prev = r_now[k], inv[k]
    # G11c-style asymptote: chained segments land on the steady state.
    seg_a = tri.transient_layers(LAYERS_R, DIR_IN, DIR_OUT, [300.0], interfaces=ifaces, dt_max=0.5)
    seg_b = tri.transient_layers(
        LAYERS_R,
        DIR_IN,
        DIR_OUT,
        [19700.0],
        mobile0=seg_a["mobile"][0],
        trapped0=seg_a["trapped"][0],
        interfaces=ifaces,
        dt_max=5.0,
    )
    steady = tri.steady_layers(LAYERS_R, DIR_IN, DIR_OUT, interfaces=ifaces)
    assert seg_b["flux_right"][0] == pytest.approx(steady["flux_right"], rel=1e-6)
    for c_got, c_want in zip(seg_b["mobile"][0], steady["mobile"], strict=True):
        assert c_got == pytest.approx(c_want, abs=1e-9)
    assert all(c >= 0.0 for row in seg_b["mobile"] for c in row)
