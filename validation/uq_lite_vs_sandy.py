"""UQ-lite sampling cross-check (`nucleide.uq` vs synthetic gates + SANDY).

Two tiers:

1. Synthetic gates (always run): U1/U2 covariance recovery on the
   hand-built `fixtures/uq/cov_2x2.json` / `cov_3x3.json` blocks (sample
   mean/covariance within k standard errors at the pinned seeds — the same
   statistical derivation as the Rust tests, never round-number guesses),
   U3 decay-branch deficit preservation + seed reproducibility, U4
   eigen-clip fallback on a rank-deficient block. No evaluated data: every
   input is a synthetic round value.
2. SANDY cross-check: the same nucleide draws wrapped in
   ``sandy.samples.Samples`` (rows = variables, columns = realizations) with
   ``get_mean``/``get_cov`` compared against ``nucleide.uq`` ``sample_mean``/
   ``sample_cov`` at 1e-9. SANDY is an optional oracle dependency (PyPI
   ``sandy`` 1.1.0, pure wheel, no NJOY needed for this path): if it — or
   pandas — cannot be imported, tier 2 is reported as SKIP with its reason
   (never silently). Tape-driven comparisons (``sandy.sampling`` CLI over
   ENDF files, ERRORR covariances) stay NJOY-gated skips: no tapes are
   vendored here by design.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

from common import Report, fmt, rel_diff

import nucleide.uq as uq

FAILURES = 0
REPO_ROOT = Path(__file__).resolve().parent.parent
FIX = REPO_ROOT / "fixtures" / "uq"


def _check(ok: bool, label: str) -> str:
    global FAILURES
    if not ok:
        FAILURES += 1
        print(f"FAIL: {label}", file=sys.stderr)
    return "PASS" if ok else "FAIL"


def _load(name: str) -> dict:
    return json.loads((FIX / name).read_text())


def _recovery_gate(fx: dict, label: str) -> tuple[list[str], str]:
    """One covariance-recovery gate; returns (table row, prose note)."""
    mean, cov, n, k, seed = fx["mean"], fx["cov"], fx["n"], fx["k"], fx["seed"]
    out = uq.sample_mvn(mean, cov, n, seed)
    dim = len(mean)
    sm = uq.sample_mean(out["samples"])
    sc = uq.sample_cov(out["samples"])
    worst = 0.0
    for i in range(dim):
        worst = max(worst, abs(sm[i] - mean[i]) / (math.sqrt(cov[i][i] / n) or 1.0))
        for j in range(dim):
            se = math.sqrt((cov[i][i] * cov[j][j] + cov[i][j] ** 2) / (n - 1))
            worst = max(worst, abs(sc[i][j] - cov[i][j]) / (se or 1.0))
    ok = out["method"] == "cholesky" and worst <= k
    note = f"{label}: method {out['method']}, worst moment error {worst:.3f} SE (k = {k})."
    return [label, f"<= {k} SE", f"{worst:.3f} SE", _check(ok, label)], note


def tier1() -> tuple[list[list[str]], list[str]]:
    """Synthetic gates U1-U4 (always run)."""
    rows: list[list[str]] = []
    notes: list[str] = []

    row, note = _recovery_gate(_load("cov_2x2.json"), "U1 2x2 recovery")
    rows.append(row)
    notes.append(note)
    row, note = _recovery_gate(_load("cov_3x3.json"), "U2 3x3 recovery")
    rows.append(row)
    notes.append(note)

    fx = _load("decay_perturb.json")
    branches = uq.perturb_branches(fx["base_branches"], fx["rel"])
    deficit_ok = abs(sum(branches) - sum(fx["base_branches"])) < 1e-12
    a = uq.sample_mvn([0.0], [[1.0]], 8, 42)
    b = uq.sample_mvn([0.0], [[1.0]], 8, 42)
    repro_ok = a["samples"] == b["samples"]
    ok = deficit_ok and repro_ok and all(v >= 0.0 for v in branches)
    notes.append(f"U3 deficit kept {sum(branches):.6f}; seed reproducibility {repro_ok}.")
    rows.append(
        [
            "U3 deficit+repro",
            "deficit kept, same seed identical",
            f"{sum(branches):.6f}",
            _check(ok, "U3 deficit+repro"),
        ]
    )

    out = uq.sample_mvn([0.0, 0.0], [[1.0, 1.0], [1.0, 1.0]], 20000, 20260913)
    ok = out["method"] == "eigen_clip"
    notes.append(f"U4 rank-1 block factorisation path: {out['method']}.")
    rows.append(["U4 eigen-clip path", "eigen_clip", out["method"], _check(ok, "U4 path")])
    return rows, notes


def tier2() -> tuple[list[list[str]], list[str], bool]:
    """Live SANDY moment cross-check; SKIP loudly when absent.

    Returns (rows, notes, skipped). When SANDY + pandas import, the rows are
    4-column gate rows; otherwise a 2-column SKIP table.
    """
    try:
        import pandas as pd
        import sandy
        from sandy.samples import Samples
    except Exception as exc:
        note = f"Tier 2 (SANDY cross-check) SKIPPED: {exc}"
        return [["SANDY Samples moments", "SKIP (see prose)"]], [note], True

    try:
        import importlib.metadata as md

        sandy_version = md.version("sandy")
    except Exception:
        sandy_version = "unknown"

    fx = _load("cov_2x2.json")
    out = uq.sample_mvn(fx["mean"], fx["cov"], fx["n"], fx["seed"])
    frame = pd.DataFrame(list(map(list, zip(*out["samples"], strict=True))))
    wrapped = Samples(frame)
    sandy_mean = list(wrapped.get_mean().values)
    sandy_cov = wrapped.get_cov().values.tolist()
    sm = uq.sample_mean(out["samples"])
    sc = uq.sample_cov(out["samples"])
    worst_mean = max(rel_diff(a, b) for a, b in zip(sandy_mean, sm, strict=True))
    worst_cov = max(
        rel_diff(a, b)
        for ra, rb in zip(sandy_cov, sc, strict=True)
        for a, b in zip(ra, rb, strict=True)
    )
    ok_mean = worst_mean < 1e-9
    ok_cov = worst_cov < 1e-9
    notes = [
        f"Tier 2 runs SANDY {sandy_version} Samples.get_mean/get_cov over the "
        "nucleide U1 draws (tape-free: no ENDF input, no NJOY).",
        "Tape-driven sandy.sampling/ERRORR comparisons stay NJOY-gated skips "
        "(no tapes vendored, by design).",
    ]
    rows = [
        ["SANDY mean", "< 1e-9", fmt(worst_mean), _check(ok_mean, "SANDY mean")],
        ["SANDY cov", "< 1e-9", fmt(worst_cov), _check(ok_cov, "SANDY cov")],
    ]
    _ = sandy
    return rows, notes, False


def main() -> int:
    report = Report("uq_lite", "UQ-lite sampling vs SANDY")
    rows1, notes1 = tier1()
    report.table(["Gate", "Expected", "Got", "Status"], rows1)
    for note in notes1:
        report.prose(note)
    rows2, notes2, skipped = tier2()
    if skipped:
        report.table(["Check", "Status"], rows2)
    else:
        report.table(["Check", "Expected", "Got", "Status"], rows2)
    for note in notes2:
        report.prose(note)
    report.emit()
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
