"""Neutron spectrum unfolding cross-check (`nucleide.unfold` vs analytic gates).

Analytic checks always run. The correlated-fit numerical comparison uses
SciPy when installed:

1. Analytic gates U1-U4 (SAND-II): forward-fold against a hand matvec,
   exact-guess fixed point, determined-system recovery from a biased guess at
   pinned tolerances, and underdetermined rate reproduction + error reduction.
2. IRDFF-II analytical benchmark-shape probes U5: the analytical fields
   named in Trkov et al., Nucl. Data Sheets 163 (2020) 1 (open access,
   arXiv:1909.03336) — a 293.6 K thermal Maxwellian, a pure 1/E field, and a
   25 keV Maxwellian — reused as published facts with citation. The IRDFF-II
   tabulated group spectra are IAEA-copyright data files and are deliberately
   NOT used; provenance stays caller-supplied.
3. STAYSL-class least-squares gates U9-U13 (Perey, ORNL/TM-6062, 1977):
   determined recovery with caller sigmas as weights, underdetermined rate
   reproduction, IRDFF-II shape recovery, the cross-method fixed-point check
   (a SAND-II fixed point is also a least-squares fixed point — both methods
   consume the same fold), and the convergence/validation contract.
4. GRAVEL chi-square-weighted gates U14-U18 (Matzke, PTB-N-19, 1994):
   determined recovery with `N^2/sigma^2` measurement weights,
   underdetermined rate reproduction under counting-statistics sigmas,
   IRDFF-II shape recovery, the cross-method fixed-point checks (at a fixed
   point every rate ratio is 1, so the weights are irrelevant — SAND-II,
   STAYSL-class, and GRAVEL fixed points coincide where the formulations
   agree), and the convergence/validation contract.
5. MAXED maximum-entropy gates U19-U23 (Reginatto & Goldhagen, Health Phys.
   77 (1999) 579): determined recovery with `1/sigma^2` chi-square weights,
   underdetermined rate reproduction under counting-statistics sigmas,
   IRDFF-II shape recovery with group bounds bracketing each shape's
   support (the solve-based house rule — trailing groups at the 1e-30
   positivity floor would otherwise read a spurious relative error of
   exactly 1.0), the cross-method fixed-point checks (at a fixed point the
   dual residual is zero, so MAXED holds SAND-II fixed points and vice
   versa), and the convergence/validation contract including the
   chi-square-target hard fail.
6. Correlated STAYSL-class gates U24-U27: diagonal equivalence, analytic
   interior and bounded optima, detector permutation and unit scaling,
   covariance validation, and comparison with a separate bounded solver.

All response matrices are hand-built synthetics (Gaussian bumps over
log-spaced groups), never evaluated data.
"""

from __future__ import annotations

import math
import sys
from collections.abc import Callable

from common import Report, fmt

import nucleide.unfold as unfold

FAILURES = 0


def _check(ok: bool, label: str) -> str:
    global FAILURES
    if not ok:
        FAILURES += 1
        print(f"FAIL: {label}", file=sys.stderr)
    return "PASS" if ok else "FAIL"


def _synthetic_response(n_det: int, n_groups: int, width: float) -> list[list[float]]:
    """Hand-built synthetic response: Gaussian bumps + a small uniform tail."""
    rows = []
    for i in range(n_det):
        center = i * (n_groups - 1) / max(n_det - 1, 1)
        rows.append([math.exp(-(((j - center) / width) ** 2)) + 1e-3 for j in range(n_groups)])
    return rows


def _log_midpoints(n_groups: int, lo: float, hi: float) -> list[float]:
    step = math.log(hi / lo) / n_groups
    return [math.exp(math.log(lo) + (j + 0.5) * step) for j in range(n_groups)]


def _maxwellian_plus_inv_e(midpoints: list[float], kt: float) -> list[float]:
    return [e * math.exp(-e / kt) + 1e-6 / e for e in midpoints]


def _max_rel_err(a: list[float], b: list[float]) -> float:
    return max(abs(x - y) / y for x, y in zip(a, b, strict=True))


def analytic_gates() -> tuple[list[list[str]], list[str]]:
    """U1-U4 on synthetic forward folds; returns (gate rows, prose notes)."""
    rows: list[list[str]] = []
    notes: list[str] = []

    response = _synthetic_response(3, 4, 1.5)
    spectrum = [1.0, 2.0, 3.0, 4.0]
    rates = unfold.forward_fold(response, spectrum)
    worst = max(
        abs(r - sum(x * p for x, p in zip(row, spectrum, strict=True)))
        for r, row in zip(rates, response, strict=True)
    )
    rows.append(["U1 forward fold", fmt(worst), "< 1e-12", _check(worst < 1e-12, "U1 fold")])

    response = _synthetic_response(4, 12, 2.0)
    midpoints = _log_midpoints(12, 1e-6, 10.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    out = unfold.sandii(response, rates, truth)
    ok = (
        out["iterations"] == 1
        and out["max_rel_change"] == 0.0
        and out["spectrum"] == truth
        and max(abs(f - 1.0) for f in out["rate_factors"]) < 1e-12
    )
    rows.append(["U2 exact fixed point", "-", "1 no-op adjust", _check(ok, "U2 fixed point")])
    notes.append("U2: an exact guess is a fixed point (one confirming no-op adjustment).")

    response = _synthetic_response(6, 6, 0.9)
    midpoints = _log_midpoints(6, 1e-6, 10.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    guess = [3.0 * p for p in truth]
    out = unfold.sandii(response, rates, guess, tolerance=1e-12, max_iterations=10_000)
    worst = _max_rel_err(out["spectrum"], truth)
    rate_err = max(abs(f - 1.0) for f in out["rate_factors"])
    ok = worst < 1e-6 and rate_err < 1e-9
    rows.append(
        ["U3 determined recovery", fmt(worst), "< 1e-6", _check(worst < 1e-6, "U3 recovery")]
    )
    rows.append(
        [
            "U3 rate factors",
            fmt(rate_err),
            "< 1e-9",
            _check(rate_err < 1e-9, "U3 rate factors"),
        ]
    )
    notes.append(
        f"U3 determined 6x6 recovery from a 3x-biased guess: {out['iterations']} adjustments."
    )

    n_groups = 24
    response = _synthetic_response(6, n_groups, 3.0)
    midpoints = _log_midpoints(n_groups, 1e-6, 12.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    # Non-uniform bias: part of it lies in the response null space and must
    # keep the guess; the recoverable part must still converge.
    guess = [p * (1.3 + 0.7 * (j % 5) / 4.0) for j, p in enumerate(truth)]
    out = unfold.sandii(response, rates, guess, tolerance=1e-9, max_iterations=50_000)
    rate_err = max(abs(f - 1.0) for f in out["rate_factors"])
    recovered_err = _max_rel_err(out["spectrum"], truth)
    guess_err = _max_rel_err(guess, truth)
    rows.append(
        [
            "U4 underdetermined rate repro",
            fmt(rate_err),
            "< 1e-6",
            _check(rate_err < 1e-6, "U4 rate repro"),
        ]
    )
    rows.append(
        [
            "U4 spectrum error reduction",
            fmt(recovered_err),
            "< guess err / 2",
            _check(recovered_err < guess_err / 2.0, "U4 error reduction"),
        ]
    )
    notes.append(
        f"U4 underdetermined case (6 detectors, {n_groups} groups): rates reproduced; "
        f"spectrum error {recovered_err:.3e} vs guess error {guess_err:.3e} — the "
        "recoverable bias component converges, the null-space part keeps the guess."
    )
    return rows, notes


def irdff_ii_probes() -> tuple[list[list[str]], list[str]]:
    """U5 IRDFF-II analytical benchmark-shape recovery (published facts)."""
    rows: list[list[str]] = []
    notes: list[str] = []
    n_groups = 6
    midpoints = _log_midpoints(n_groups, 1e-9, 20.0)
    thermal_kt = 8.617333262e-5 * 293.6  # 293.6 K in MeV
    shapes = {
        "thermal Maxwellian 293.6 K": [e * math.exp(-e / thermal_kt) + 1e-30 for e in midpoints],
        "pure 1/E": [1.0 / e + 1e-30 for e in midpoints],
        "Maxwellian 25 keV": [math.sqrt(e) * math.exp(-e / 2.5e-2) + 1e-30 for e in midpoints],
    }
    for name, shape in shapes.items():
        response = _synthetic_response(n_groups, n_groups, 0.8)
        rates = unfold.forward_fold(response, shape)
        guess = [5.0 * p for p in shape]
        out = unfold.sandii(response, rates, guess, tolerance=1e-11, max_iterations=100_000)
        worst = _max_rel_err(out["spectrum"], shape)
        rows.append([f"U5 {name}", fmt(worst), "< 1e-5", _check(worst < 1e-5, f"U5 {name}")])
    notes.append(
        "U5 shapes are the IRDFF-II analytical benchmark fields named in "
        "Trkov et al., Nucl. Data Sheets 163 (2020) 1 (arXiv:1909.03336), "
        "re-derived from their closed forms; the IAEA-copyright tabulated "
        "spectra are never used."
    )
    return rows, notes


def contract_gates() -> tuple[list[list[str]], list[str]]:
    """U6-U8 convergence and validation contract probes."""
    rows: list[list[str]] = []
    notes: list[str] = []

    response = _synthetic_response(4, 8, 1.5)
    midpoints = _log_midpoints(8, 1e-6, 10.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    guess = [2.0 * p for p in truth]
    try:
        unfold.sandii(response, rates, guess, tolerance=1e-12, max_iterations=1)
        ok = False
    except ValueError:
        ok = True
    rows.append(["U6 cap is a hard fail", "-", "ValueError", _check(ok, "U6 hard fail")])

    out = unfold.sandii([[1.0, 0.0, 0.0]], [4.0], [2.0, 7.0, 3.0], tolerance=1e-12)
    ok = out["spectrum"][0] == 4.0 and out["spectrum"][1:] == [7.0, 3.0]
    rows.append(["U7 unconstrained keeps guess", "-", "bit-exact", _check(ok, "U7 unconstrained")])

    failures: list[str] = []

    def _expect(label: str, fn: Callable[[], object], needle: str) -> None:
        try:
            fn()
            failures.append(label)
        except ValueError as exc:
            if needle not in str(exc):
                failures.append(f"{label} (wrong message: {exc})")

    _expect("rates length", lambda: unfold.sandii(response, [1.0], guess), "shape mismatch")
    _expect(
        "negative response",
        lambda: unfold.sandii([[1.0, -1.0], [1.0, 1.0]], [1.0, 1.0], [1.0, 1.0]),
        "negative entry",
    )
    _expect(
        "zero guess",
        lambda: unfold.sandii([[1.0, 1.0]], [1.0], [1.0, 0.0]),
        "strictly positive",
    )
    _expect(
        "unreachable rate",
        lambda: unfold.sandii([[0.0, 0.0]], [1.0], [1.0, 1.0]),
        "unreachable",
    )
    _expect(
        "fold shape",
        lambda: unfold.forward_fold([[1.0, 1.0]], [1.0]),
        "shape mismatch",
    )
    rows.append(["U8 named input errors", "-", "5 cases", _check(not failures, "U8 errors")])
    if failures:
        notes.append(f"U8 missing/weak rejections: {', '.join(failures)}.")
    notes.append(
        "U6-U8: non-convergence raises (never a partial spectrum), unconstrained "
        "groups keep the guess bit-exactly, and every malformed input names its cause."
    )
    return rows, notes


def staysl_gates() -> tuple[list[list[str]], list[str]]:
    """U9-U13 STAYSL-class least-squares gates (Perey, ORNL/TM-6062, 1977)."""
    rows: list[list[str]] = []
    notes: list[str] = []

    response = _synthetic_response(6, 6, 0.9)
    midpoints = _log_midpoints(6, 1e-6, 10.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    sigmas = [1.0] * len(rates)
    guess = [3.0 * p for p in truth]
    # Caller sigmas pinned as weights: a redundant sloppy reading must lose.
    sloppy = [2.0, 2.2]
    out_w = unfold.staysl(
        [[1.0, 0.0], [1.0, 0.0]],
        sloppy,
        [0.01, 1.0],
        [1.5, 7.0],
        tolerance=1e-12,
        damping=1e-8,
    )
    ok_w = abs(out_w["spectrum"][0] - 2.0) < 0.02 and out_w["spectrum"][1] == 7.0
    rows.append(["U9 caller sigmas weight", "-", "precise wins", _check(ok_w, "U9 weights")])
    out = unfold.staysl(response, rates, sigmas, guess, tolerance=1e-9, damping=1e-12)
    worst = _max_rel_err(out["spectrum"], truth)
    rate_err = max(abs(f - 1.0) for f in out["rate_factors"])
    ok = worst < 1e-6 and rate_err < 1e-6
    rows.append(
        ["U9 determined recovery", fmt(worst), "< 1e-6", _check(worst < 1e-6, "U9 recovery")]
    )
    rows.append(
        [
            "U9 rate factors",
            fmt(rate_err),
            "< 1e-6",
            _check(rate_err < 1e-6, "U9 rate factors"),
        ]
    )
    notes.append(
        f"U9 determined 6x6 recovery from a 3x-biased guess: {out['iterations']} cycles; "
        "caller-supplied per-detector sigmas are the least-squares weights "
        "(a redundant reading with 100x the sigma loses to the precise one)."
    )

    n_groups = 24
    response = _synthetic_response(6, n_groups, 3.0)
    midpoints = _log_midpoints(n_groups, 1e-6, 12.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    sigmas = [1.0] * len(rates)
    guess = [p * (1.3 + 0.7 * (j % 5) / 4.0) for j, p in enumerate(truth)]
    out = unfold.staysl(response, rates, sigmas, guess, tolerance=1e-9, damping=1e-6)
    rate_err = max(abs(f - 1.0) for f in out["rate_factors"])
    rows.append(
        [
            "U10 underdetermined rate repro",
            fmt(rate_err),
            "< 1e-4",
            _check(rate_err < 1e-4, "U10 rate repro"),
        ]
    )
    notes.append(
        f"U10 underdetermined case (6 detectors, {n_groups} groups): the least-squares "
        f"fixed point interpolates the measured rates ({out['iterations']} cycles). Unlike "
        "SAND-II's multiplicative spreading it does not promise a spectrum closer "
        "than the guess — weakly constrained group directions take their values from "
        "the fit, not the prior — so the gate pins rate reproduction only."
    )

    thermal_kt = 8.617333262e-5 * 293.6  # 293.6 K in MeV
    shapes = {
        "thermal Maxwellian 293.6 K": [
            e * math.exp(-e / thermal_kt) + 1e-30 for e in _log_midpoints(18, 1e-9, 1e-3)
        ],
        "pure 1/E": [1.0 / e + 1e-30 for e in _log_midpoints(18, 1e-6, 10.0)],
        "Maxwellian 25 keV": [
            math.sqrt(e) * math.exp(-e / 2.5e-2) + 1e-30 for e in _log_midpoints(18, 1e-3, 0.3)
        ],
    }
    for name, shape in shapes.items():
        response = _synthetic_response(18, 18, 0.8)
        rates = unfold.forward_fold(response, shape)
        sigmas = [1.0] * len(rates)
        guess = [5.0 * p for p in shape]
        out = unfold.staysl(response, rates, sigmas, guess, tolerance=1e-9, damping=1e-10)
        worst = _max_rel_err(out["spectrum"], shape)
        rows.append([f"U11 {name}", fmt(worst), "< 1e-5", _check(worst < 1e-5, f"U11 {name}")])
    notes.append(
        "U11 shapes are the IRDFF-II analytical benchmark fields named in "
        "Trkov et al., Nucl. Data Sheets 163 (2020) 1 (arXiv:1909.03336), "
        "re-derived from their closed forms with group bounds bracketing each "
        "shape's support; the IAEA-copyright tabulated spectra are never used."
    )

    response = _synthetic_response(6, 6, 1.0)
    midpoints = _log_midpoints(6, 1e-6, 10.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    guess = [2.5 * p for p in truth]
    sandii_out = unfold.sandii(response, rates, guess, tolerance=1e-12, max_iterations=100_000)
    sigmas = [1.0] * len(rates)
    out = unfold.staysl(
        response,
        rates,
        sigmas,
        sandii_out["spectrum"],
        tolerance=1e-4,
        damping=1e-6,
    )
    worst = _max_rel_err(out["spectrum"], sandii_out["spectrum"])
    ok = out["iterations"] == 1 and out["max_rel_change"] < 1e-4 and worst < 1e-5
    rows.append(
        [
            "U12 SAND-II fp is a LS fp",
            fmt(worst),
            "< 1e-5",
            _check(ok, "U12 cross fixed point"),
        ]
    )
    notes.append(
        "U12: both methods consume the same fold — once SAND-II has converged, "
        "the anchored least-squares cycle about its spectrum is a no-op, so a "
        "SAND-II fixed point is a least-squares fixed point too."
    )

    failures: list[str] = []

    def _expect(label: str, fn: Callable[[], object], needle: str) -> None:
        try:
            fn()
            failures.append(label)
        except ValueError as exc:
            if needle not in str(exc):
                failures.append(f"{label} (wrong message: {exc})")

    response = _synthetic_response(4, 8, 1.5)
    midpoints = _log_midpoints(8, 1e-6, 10.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    sigmas = [1.0] * len(rates)
    guess = [2.0 * p for p in truth]
    _expect(
        "cap exhaustion",
        lambda: unfold.staysl(response, rates, sigmas, guess, tolerance=1e-12, max_iterations=1),
        "did not converge",
    )
    _expect("sigmas length", lambda: unfold.staysl(response, rates, [1.0], guess), "shape mismatch")
    _expect(
        "zero sigma",
        lambda: unfold.staysl(response, rates, [0.0] * len(rates), guess),
        "sigma must be > 0",
    )
    _expect(
        "unreachable rate",
        lambda: unfold.staysl([[0.0, 0.0, 0.0]], [1.0], [1.0], guess),
        "unreachable",
    )
    rows.append(["U13 named contract errors", "-", "4 cases", _check(not failures, "U13 errors")])
    if failures:
        notes.append(f"U13 missing/weak rejections: {', '.join(failures)}.")
    notes.append(
        "U13: non-convergence raises (never a partial spectrum), malformed sigmas "
        "name their cause, and unreachable rates are loud."
    )
    return rows, notes


def gravel_gates() -> tuple[list[list[str]], list[str]]:
    """U14-U18 GRAVEL chi-square-weighted gates (Matzke, PTB-N-19, 1994)."""
    rows: list[list[str]] = []
    notes: list[str] = []

    response = _synthetic_response(6, 6, 0.9)
    midpoints = _log_midpoints(6, 1e-6, 10.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    sigmas = [1.0] * len(rates)
    guess = [3.0 * p for p in truth]
    # Caller sigmas pinned as weights: a redundant sloppy reading must lose.
    out_w = unfold.gravel(
        [[1.0, 0.0], [1.0, 0.0]],
        [2.0, 2.2],
        [0.01, 1.0],
        [1.5, 7.0],
        tolerance=1e-12,
    )
    ok_w = abs(out_w["spectrum"][0] - 2.0) < 0.02 and out_w["spectrum"][1] == 7.0
    rows.append(["U14 caller sigmas weight", "-", "precise wins", _check(ok_w, "U14 weights")])
    out = unfold.gravel(response, rates, sigmas, guess, tolerance=1e-12, max_iterations=10_000)
    worst = _max_rel_err(out["spectrum"], truth)
    rate_err = max(abs(f - 1.0) for f in out["rate_factors"])
    rows.append(
        ["U14 determined recovery", fmt(worst), "< 1e-6", _check(worst < 1e-6, "U14 recovery")]
    )
    rows.append(
        [
            "U14 rate factors",
            fmt(rate_err),
            "< 1e-9",
            _check(rate_err < 1e-9, "U14 rate factors"),
        ]
    )
    notes.append(
        f"U14 determined 6x6 recovery from a 3x-biased guess: {out['iterations']} adjustments; "
        "caller-supplied per-detector sigmas are the chi-square weights "
        "(a redundant reading with 100x the sigma loses to the precise one)."
    )

    n_groups = 24
    response = _synthetic_response(6, n_groups, 3.0)
    midpoints = _log_midpoints(n_groups, 1e-6, 12.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    sigmas = [math.sqrt(r) for r in rates]
    guess = [p * (1.3 + 0.7 * (j % 5) / 4.0) for j, p in enumerate(truth)]
    out = unfold.gravel(response, rates, sigmas, guess, tolerance=1e-9, max_iterations=50_000)
    rate_err = max(abs(f - 1.0) for f in out["rate_factors"])
    recovered_err = _max_rel_err(out["spectrum"], truth)
    guess_err = _max_rel_err(guess, truth)
    rows.append(
        [
            "U15 underdetermined rate repro",
            fmt(rate_err),
            "< 1e-6",
            _check(rate_err < 1e-6, "U15 rate repro"),
        ]
    )
    rows.append(
        [
            "U15 spectrum error reduction",
            fmt(recovered_err),
            "< guess err",
            _check(recovered_err < guess_err, "U15 error reduction"),
        ]
    )
    notes.append(
        f"U15 underdetermined case (6 detectors, {n_groups} groups) under "
        f"counting-statistics sigmas: rates reproduced ({out['iterations']} "
        f"adjustments), spectrum error {recovered_err:.3e} vs guess error "
        f"{guess_err:.3e}."
    )

    thermal_kt = 8.617333262e-5 * 293.6  # 293.6 K in MeV
    shapes = {
        "thermal Maxwellian 293.6 K": [
            e * math.exp(-e / thermal_kt) + 1e-30 for e in _log_midpoints(18, 1e-9, 1e-3)
        ],
        "pure 1/E": [1.0 / e + 1e-30 for e in _log_midpoints(18, 1e-6, 10.0)],
        "Maxwellian 25 keV": [
            math.sqrt(e) * math.exp(-e / 2.5e-2) + 1e-30 for e in _log_midpoints(18, 1e-3, 0.3)
        ],
    }
    for name, shape in shapes.items():
        response = _synthetic_response(18, 18, 0.8)
        rates = unfold.forward_fold(response, shape)
        sigmas = [1.0] * len(rates)
        guess = [5.0 * p for p in shape]
        out = unfold.gravel(response, rates, sigmas, guess, tolerance=1e-11, max_iterations=100_000)
        worst = _max_rel_err(out["spectrum"], shape)
        rows.append([f"U16 {name}", fmt(worst), "< 1e-5", _check(worst < 1e-5, f"U16 {name}")])
    notes.append(
        "U16 shapes are the IRDFF-II analytical benchmark fields named in "
        "Trkov et al., Nucl. Data Sheets 163 (2020) 1 (arXiv:1909.03336), "
        "re-derived from their closed forms with group bounds bracketing each "
        "shape's support; the IAEA-copyright tabulated spectra are never used."
    )

    response = _synthetic_response(6, 6, 1.0)
    midpoints = _log_midpoints(6, 1e-6, 10.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    guess = [2.5 * p for p in truth]
    sigmas = [1.0] * len(rates)
    sandii_out = unfold.sandii(response, rates, guess, tolerance=1e-12, max_iterations=100_000)
    out = unfold.gravel(response, rates, sigmas, sandii_out["spectrum"], tolerance=1e-4)
    worst = _max_rel_err(out["spectrum"], sandii_out["spectrum"])
    ok = out["iterations"] == 1 and out["max_rel_change"] < 1e-4 and worst < 1e-5
    rows.append(
        [
            "U17 SAND-II fp is a GRAVEL fp",
            fmt(worst),
            "< 1e-5",
            _check(ok, "U17 SAND-II cross fixed point"),
        ]
    )
    gravel_out = unfold.gravel(
        response, rates, sigmas, guess, tolerance=1e-12, max_iterations=100_000
    )
    back = unfold.sandii(response, rates, gravel_out["spectrum"], tolerance=1e-4)
    worst = _max_rel_err(back["spectrum"], gravel_out["spectrum"])
    ok = back["iterations"] == 1 and back["max_rel_change"] < 1e-4 and worst < 1e-5
    rows.append(
        [
            "U17 GRAVEL fp is a SAND-II fp",
            fmt(worst),
            "< 1e-5",
            _check(ok, "U17 GRAVEL cross fixed point"),
        ]
    )
    ls_out = unfold.staysl(
        response, rates, sigmas, gravel_out["spectrum"], tolerance=1e-4, damping=1e-6
    )
    worst = _max_rel_err(ls_out["spectrum"], gravel_out["spectrum"])
    ok = ls_out["iterations"] == 1 and ls_out["max_rel_change"] < 1e-4 and worst < 1e-5
    rows.append(
        [
            "U17 GRAVEL fp is a LS fp",
            fmt(worst),
            "< 1e-5",
            _check(ok, "U17 least-squares cross fixed point"),
        ]
    )
    notes.append(
        "U17: all three methods consume the same fold — at a fixed point every "
        "rate ratio is 1, so the GRAVEL weights are irrelevant and each "
        "method holds the others' fixed points after one confirming cycle."
    )

    failures: list[str] = []

    def _expect(label: str, fn: Callable[[], object], needle: str) -> None:
        try:
            fn()
            failures.append(label)
        except ValueError as exc:
            if needle not in str(exc):
                failures.append(f"{label} (wrong message: {exc})")

    response = _synthetic_response(4, 8, 1.5)
    midpoints = _log_midpoints(8, 1e-6, 10.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    sigmas = [1.0] * len(rates)
    guess = [2.0 * p for p in truth]
    _expect(
        "cap exhaustion",
        lambda: unfold.gravel(response, rates, sigmas, guess, tolerance=1e-12, max_iterations=1),
        "did not converge",
    )
    _expect("sigmas length", lambda: unfold.gravel(response, rates, [1.0], guess), "shape mismatch")
    _expect(
        "zero sigma",
        lambda: unfold.gravel(response, rates, [0.0] * len(rates), guess),
        "sigma must be > 0",
    )
    _expect(
        "non-finite sigma",
        lambda: unfold.gravel(response, rates, [1.0] * (len(rates) - 1) + [float("nan")], guess),
        "non-finite sigma",
    )
    _expect(
        "unreachable rate",
        lambda: unfold.gravel([[0.0, 0.0, 0.0]], [1.0], [1.0], guess),
        "unreachable",
    )
    rows.append(["U18 named contract errors", "-", "5 cases", _check(not failures, "U18 errors")])
    if failures:
        notes.append(f"U18 missing/weak rejections: {', '.join(failures)}.")
    notes.append(
        "U18: non-convergence raises (never a partial spectrum), malformed sigmas "
        "name their cause, and unreachable rates are loud."
    )
    return rows, notes


def maxed_gates() -> tuple[list[list[str]], list[str]]:
    """U19-U23 MAXED maximum-entropy gates (Reginatto & Goldhagen 1999)."""
    rows: list[list[str]] = []
    notes: list[str] = []

    response = _synthetic_response(6, 6, 0.9)
    midpoints = _log_midpoints(6, 1e-6, 10.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    sigmas = [1.0] * len(rates)
    guess = [3.0 * p for p in truth]
    # Caller sigmas pinned as weights: a redundant sloppy reading must lose.
    out_w = unfold.maxed(
        [[1.0, 0.0], [1.0, 0.0]],
        [2.0, 2.2],
        [0.01, 1.0],
        [1.5, 7.0],
        tolerance=1e-12,
    )
    ok_w = abs(out_w["spectrum"][0] - 2.0) < 0.02 and out_w["spectrum"][1] == 7.0
    rows.append(["U19 caller sigmas weight", "-", "precise wins", _check(ok_w, "U19 weights")])
    out = unfold.maxed(response, rates, sigmas, guess, tolerance=1e-12, max_iterations=10_000)
    worst = _max_rel_err(out["spectrum"], truth)
    rate_err = max(abs(f - 1.0) for f in out["rate_factors"])
    rows.append(
        ["U19 determined recovery", fmt(worst), "< 1e-6", _check(worst < 1e-6, "U19 recovery")]
    )
    rows.append(
        [
            "U19 rate factors",
            fmt(rate_err),
            "< 1e-9",
            _check(rate_err < 1e-9, "U19 rate factors"),
        ]
    )
    notes.append(
        f"U19 determined 6x6 recovery from a 3x-biased guess: {out['iterations']} adjustments; "
        "caller-supplied per-detector sigmas are the chi-square weights "
        "(a redundant reading with 100x the sigma loses to the precise one)."
    )

    n_groups = 24
    response = _synthetic_response(6, n_groups, 3.0)
    midpoints = _log_midpoints(n_groups, 1e-6, 12.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    sigmas = [math.sqrt(r) for r in rates]
    guess = [p * (1.3 + 0.7 * (j % 5) / 4.0) for j, p in enumerate(truth)]
    out = unfold.maxed(response, rates, sigmas, guess, tolerance=1e-9, max_iterations=50_000)
    rate_err = max(abs(f - 1.0) for f in out["rate_factors"])
    recovered_err = _max_rel_err(out["spectrum"], truth)
    guess_err = _max_rel_err(guess, truth)
    rows.append(
        [
            "U20 underdetermined rate repro",
            fmt(rate_err),
            "< 1e-6",
            _check(rate_err < 1e-6, "U20 rate repro"),
        ]
    )
    rows.append(
        [
            "U20 spectrum error reduction",
            fmt(recovered_err),
            "< guess err",
            _check(recovered_err < guess_err, "U20 error reduction"),
        ]
    )
    notes.append(
        f"U20 underdetermined case (6 detectors, {n_groups} groups) under "
        f"counting-statistics sigmas: rates reproduced ({out['iterations']} "
        f"adjustments), spectrum error {recovered_err:.3e} vs guess error "
        f"{guess_err:.3e}."
    )

    thermal_kt = 8.617333262e-5 * 293.6  # 293.6 K in MeV
    shapes = {
        "thermal Maxwellian 293.6 K": [
            e * math.exp(-e / thermal_kt) + 1e-30 for e in _log_midpoints(18, 1e-9, 1e-3)
        ],
        "pure 1/E": [1.0 / e + 1e-30 for e in _log_midpoints(18, 1e-6, 10.0)],
        "Maxwellian 25 keV": [
            math.sqrt(e) * math.exp(-e / 2.5e-2) + 1e-30 for e in _log_midpoints(18, 1e-3, 0.3)
        ],
    }
    for name, shape in shapes.items():
        response = _synthetic_response(18, 18, 0.8)
        rates = unfold.forward_fold(response, shape)
        sigmas = [1.0] * len(rates)
        guess = [5.0 * p for p in shape]
        out = unfold.maxed(response, rates, sigmas, guess, tolerance=1e-11, max_iterations=100_000)
        worst = _max_rel_err(out["spectrum"], shape)
        rows.append([f"U21 {name}", fmt(worst), "< 1e-5", _check(worst < 1e-5, f"U21 {name}")])
    notes.append(
        "U21 shapes are the IRDFF-II analytical benchmark fields named in "
        "Trkov et al., Nucl. Data Sheets 163 (2020) 1 (arXiv:1909.03336), "
        "re-derived from their closed forms with group bounds bracketing each "
        "shape's support; the IAEA-copyright tabulated spectra are never used."
    )

    response = _synthetic_response(6, 6, 1.0)
    midpoints = _log_midpoints(6, 1e-6, 10.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    guess = [2.5 * p for p in truth]
    sigmas = [1.0] * len(rates)
    sandii_out = unfold.sandii(response, rates, guess, tolerance=1e-12, max_iterations=100_000)
    out = unfold.maxed(response, rates, sigmas, sandii_out["spectrum"], tolerance=1e-4)
    worst = _max_rel_err(out["spectrum"], sandii_out["spectrum"])
    ok = out["iterations"] == 1 and out["max_rel_change"] < 1e-4 and worst < 1e-5
    rows.append(
        [
            "U22 SAND-II fp is a MAXED fp",
            fmt(worst),
            "< 1e-5",
            _check(ok, "U22 SAND-II cross fixed point"),
        ]
    )
    maxed_out = unfold.maxed(
        response, rates, sigmas, guess, tolerance=1e-12, max_iterations=100_000
    )
    back = unfold.sandii(response, rates, maxed_out["spectrum"], tolerance=1e-4)
    worst = _max_rel_err(back["spectrum"], maxed_out["spectrum"])
    ok = back["iterations"] == 1 and back["max_rel_change"] < 1e-4 and worst < 1e-5
    rows.append(
        [
            "U22 MAXED fp is a SAND-II fp",
            fmt(worst),
            "< 1e-5",
            _check(ok, "U22 MAXED cross fixed point"),
        ]
    )
    notes.append(
        "U22: both methods consume the same fold — at a fixed point the dual "
        "residual is zero, so MAXED holds the SAND-II point (and vice versa) "
        "after one confirming adjustment."
    )

    failures: list[str] = []

    def _expect(label: str, fn: Callable[[], object], needle: str) -> None:
        try:
            fn()
            failures.append(label)
        except ValueError as exc:
            if needle not in str(exc):
                failures.append(f"{label} (wrong message: {exc})")

    response = _synthetic_response(4, 8, 1.5)
    midpoints = _log_midpoints(8, 1e-6, 10.0)
    truth = _maxwellian_plus_inv_e(midpoints, 2.53e-5)
    rates = unfold.forward_fold(response, truth)
    sigmas = [1.0] * len(rates)
    guess = [2.0 * p for p in truth]
    _expect(
        "cap exhaustion",
        lambda: unfold.maxed(response, rates, sigmas, guess, tolerance=1e-12, max_iterations=1),
        "did not converge",
    )
    _expect(
        "chi-square target violation",
        lambda: unfold.maxed(
            [[1.0, 0.0], [1.0, 0.0]],
            [2.0, 2.1],
            [0.01, 0.01],
            [1.5, 7.0],
            target_chi2=1.0,
        ),
        "did not converge",
    )
    _expect("sigmas length", lambda: unfold.maxed(response, rates, [1.0], guess), "shape mismatch")
    _expect(
        "zero sigma",
        lambda: unfold.maxed(response, rates, [0.0] * len(rates), guess),
        "sigma must be > 0",
    )
    _expect(
        "non-finite sigma",
        lambda: unfold.maxed(response, rates, [1.0] * (len(rates) - 1) + [float("nan")], guess),
        "non-finite sigma",
    )
    _expect(
        "bad target",
        lambda: unfold.maxed(response, rates, sigmas, guess, target_chi2=-1.0),
        "target_chi2",
    )
    _expect(
        "unreachable rate",
        lambda: unfold.maxed([[0.0, 0.0, 0.0]], [1.0], [1.0], guess),
        "unreachable",
    )
    rows.append(["U23 named contract errors", "-", "7 cases", _check(not failures, "U23 errors")])
    if failures:
        notes.append(f"U23 missing/weak rejections: {', '.join(failures)}.")
    notes.append(
        "U23: non-convergence raises (never a partial spectrum) — both on "
        "cap exhaustion and when the chi-square bottoms out above an "
        "unreachable target on disagreeing redundant readings; malformed "
        "sigmas, targets, and unreachable rates name their cause."
    )
    return rows, notes


def correlated_staysl_gates() -> tuple[list[list[str]], list[str]]:
    """U24-U26 analytic and input-contract checks."""
    identity = [[1.0, 0.0], [0.0, 1.0]]
    rates = [3.0, 2.0]
    guess = [1.0, 1.0]
    covariance = [[1.0, 0.5], [0.5, 1.0]]
    fit = unfold.staysl_correlated(identity, rates, covariance, guess, damping=1.0)
    diagonal = unfold.staysl_correlated(identity, rates, identity, guess, damping=1.0)
    sigma = unfold.staysl(identity, rates, [1.0, 1.0], guess, damping=1.0)
    errors = [
        abs(actual - expected)
        for actual, expected in zip(fit["spectrum"], [29 / 15, 19 / 15], strict=True)
    ]
    errors.extend(
        abs(actual - expected)
        for actual, expected in zip(diagonal["spectrum"], sigma["spectrum"], strict=True)
    )
    errors.extend(
        [abs(fit["measurement_chi2"] - 268 / 225), abs(fit["damping_penalty"] - 212 / 225)]
    )
    error = max(errors)
    rows = [
        [
            "U24 interior fit and diagonal equivalence",
            fmt(error),
            "< 1e-12",
            _check(error < 1e-12, "U24 interior and diagonal fits"),
        ]
    ]

    response = [[0.0, 3.0, 3.0], [0.0, 1.0, 0.0], [3.0, 2.0, 0.0]]
    rates = [6.0, 4.0, 1.0]
    guess = [1.0, 1.0, 1.0]
    covariance = [[1.0, 0.75, 0.0], [0.75, 1.0, 0.0], [0.0, 0.0, 1.0]]
    fit = unfold.staysl_correlated(response, rates, covariance, guess, damping=1.0)
    permutation = [2, 0, 1]
    scale = [1e-60, 1e60, 3.0]
    permuted = unfold.staysl_correlated(
        [[scale[i] * value for value in response[i]] for i in permutation],
        [scale[i] * rates[i] for i in permutation],
        [[scale[i] * covariance[i][j] * scale[j] for j in permutation] for i in permutation],
        guess,
        damping=1.0,
    )
    error = max(
        abs(result - expected)
        for spectrum in (fit["spectrum"], permuted["spectrum"])
        for result, expected in zip(spectrum, [0.0, 151 / 141, 11 / 47], strict=True)
    )
    ok = error < 1e-12 and all(
        out["spectrum"][0] == 0.0 and out["optimality"] <= out["tolerance"]
        for out in (fit, permuted)
    )
    rows.append(
        ["U25 bounded optimum and rescaling", fmt(error), "< 1e-12", _check(ok, "U25 bound")]
    )
    failed = []
    for invalid, match in [
        ([[1.0, 0.0], [0.0]], "square matrix"),
        ([[1.0, 0.1], [0.9, 1.0]], "symmetric"),
        ([[1.0, 1.0], [1.0, 1.0]], "positive definite"),
        ([[1.0, 1 - 1e-13], [1 - 1e-13, 1.0]], "ill-conditioned"),
    ]:
        try:
            unfold.staysl_correlated(identity, [3.0, 2.0], invalid, [1.0, 1.0])
        except ValueError as exc:
            if match in str(exc):
                continue
        failed.append(match)
    try:
        unfold.staysl_correlated(response, rates, covariance, guess, damping=1.0, max_iterations=1)
    except ValueError as exc:
        if "did not converge" not in str(exc):
            failed.append("iteration cap")
    else:
        failed.append("iteration cap")
    rows.append(
        [
            "U26 covariance and iteration errors",
            ", ".join(failed) or "-",
            "5 cases",
            _check(not failed, "U26 inputs and iteration cap"),
        ]
    )
    return rows, [
        "U24 compares the interior solution and both objective terms to closed forms. "
        "U25 verifies a bounded optimum under detector permutation and changes of rate units."
    ]


def scipy_bounded_gate() -> tuple[list[list[str]], list[str]]:
    """U27 independent bounded least-squares reference on synthetic data."""
    try:
        import numpy as np
        import scipy
        from scipy.linalg import cholesky, solve_triangular
        from scipy.optimize import lsq_linear
    except ModuleNotFoundError as exc:
        if exc.name not in ("numpy", "scipy"):
            raise
        message = f"U27 skipped because {exc.name} is unavailable"
        print(f"SKIP: {message}", file=sys.stderr)
        return [["U27 SciPy bounded comparison", "-", "-", "SKIP"]], [f"{message}."]

    rng = np.random.default_rng(702)
    max_spectrum_error = max_objective_error = 0.0
    interior = bound = 0
    references_converged = True
    for _ in range(40):
        detectors, groups = rng.integers(2, 8, size=2)
        correlations = rng.normal(size=(detectors, detectors))
        covariance = correlations @ correlations.T + 0.15 * np.eye(detectors)
        response = rng.uniform(0.0, 4.0, size=(detectors, groups))
        rates = rng.uniform(0.0, 6.0, size=detectors)
        guess = rng.uniform(0.1, 3.0, size=groups)
        damping = 10.0 ** rng.uniform(-3.0, 1.0)
        factor = cholesky(covariance, lower=True)
        design = np.vstack(
            (solve_triangular(factor, response, lower=True), math.sqrt(damping) * np.eye(groups))
        )
        target = np.r_[solve_triangular(factor, rates, lower=True), math.sqrt(damping) * guess]
        reference = lsq_linear(design, target, bounds=(0.0, np.inf), method="bvls", tol=1e-12)
        references_converged &= reference.success
        fit = unfold.staysl_correlated(
            response.tolist(), rates.tolist(), covariance.tolist(), guess.tolist(), damping=damping
        )
        if any(value == 0.0 for value in fit["spectrum"]):
            bound += 1
        else:
            interior += 1
        spectrum_error = max(abs(np.asarray(fit["spectrum"]) - reference.x))
        reference_objective = 2.0 * reference.cost
        objective_error = abs(
            fit["measurement_chi2"] + fit["damping_penalty"] - reference_objective
        )
        max_spectrum_error = max(max_spectrum_error, spectrum_error)
        max_objective_error = max(
            max_objective_error, objective_error / max(1.0, reference_objective)
        )
    rows = [
        [
            "U27 SciPy spectrum max abs error",
            fmt(max_spectrum_error),
            "< 1e-8",
            _check(max_spectrum_error < 1e-8, "U27 spectrum"),
        ],
        [
            "U27 SciPy objective scaled abs error",
            fmt(max_objective_error),
            "< 1e-10",
            _check(max_objective_error < 1e-10, "U27 objective"),
        ],
        [
            "U27 converged references and active-set coverage",
            f"{interior} interior, {bound} on a bound",
            "both kinds",
            _check(references_converged and interior > 0 and bound > 0, "U27 case coverage"),
        ],
    ]
    return rows, [
        f"U27 uses SciPy {scipy.__version__} BVLS with its own Cholesky factorization "
        "on 40 seeded synthetic cases."
    ]


def main() -> int:
    report = Report("unfold", "Neutron spectrum unfolding (`unfold_vs_analytic.py`)")
    report.prose(
        "Validation for `nucleide.unfold`: analytic "
        "gates U1-U4 on synthetic forward folds (fold identity, exact fixed "
        "point, determined and underdetermined round-trips at pinned "
        "tolerances), IRDFF-II analytical benchmark-shape probes U5 "
        "(Trkov et al. 2020, published facts with citation), STAYSL-class "
        "least-squares gates U9-U13, GRAVEL chi-square-weighted gates "
        "U14-U18, MAXED maximum-entropy gates U19-U23, and correlated "
        "STAYSL-class gates U24-U26. The U27 numerical comparison runs when "
        "SciPy is installed."
    )
    rows1, notes1 = analytic_gates()
    for note in notes1:
        report.prose(note)
    report.table(["Gate", "Value", "Tol", "Status"], rows1)
    rows2, notes2 = irdff_ii_probes()
    for note in notes2:
        report.prose(note)
    report.table(["Gate", "Rel err", "Tol", "Status"], rows2)
    rows3, notes3 = contract_gates()
    for note in notes3:
        report.prose(note)
    report.table(["Gate", "Value", "Tol", "Status"], rows3)
    rows4, notes4 = staysl_gates()
    for note in notes4:
        report.prose(note)
    report.table(["Gate", "Value", "Tol", "Status"], rows4)
    rows5, notes5 = gravel_gates()
    for note in notes5:
        report.prose(note)
    report.table(["Gate", "Value", "Tol", "Status"], rows5)
    rows6, notes6 = maxed_gates()
    for note in notes6:
        report.prose(note)
    report.table(["Gate", "Value", "Tol", "Status"], rows6)
    for gates in (correlated_staysl_gates, scipy_bounded_gate):
        rows, notes = gates()
        for note in notes:
            report.prose(note)
        report.table(["Gate", "Value", "Tol", "Status"], rows)
    report.emit()

    if FAILURES:
        print(f"FAIL: {FAILURES} unfold check(s) failed", file=sys.stderr)
        return 1
    print("All unfold checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
