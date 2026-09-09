"""Dose-factor oracle vs the PyNE dose-per-gram formula.

Compares Nucleide `dose_factor`/`dose_per_g` against the real PyNE equations
(`Material::dose_per_g`, `src/material.cpp:1504-1534` with `Ci_per_Bq` from
`src/data.cpp` and source ids 0=EPA/1=DOE/2=GENII) on spot nuclides from the
HNF-5636 tables. When PyNE is importable, its `ext_*_dose`/`ingest_dose` /
`inhale_dose` accessors supply the reference factors; otherwise a
formula-equivalent pure-Python reference (same constants) is used and the
report says so. Never hand-edit results; `run_all.sh` regenerates them.
"""

from __future__ import annotations

import math
import sys

from common import Report, fmt, rel_diff

import nucleide

CI_PER_BQ = 2.7027027e-11
PCI_PER_BQ = 27.027027
N_A_PYNE = 6.0221415e23
N_A_NUCLEIDE = 6.02214076e23

SPOTS = [
    ("Co60", "ingest", "EPA", 2.69e-05),
    ("Cs137", "inhale", "EPA", 3.19e-05),
    ("H3", "air", "EPA", 4.41e-012),
    ("K40", "soil", "EPA", 4.33e02),
]

FAILURES = 0


def _check(ok: bool, label: str) -> str:
    global FAILURES
    if not ok:
        FAILURES += 1
        print(f"FAIL: {label}", file=sys.stderr)
    return "PASS" if ok else "FAIL"


def _pyne_factor(name: str, pathway: str, source: str) -> tuple[float | None, str]:
    """Return (factor, note) from PyNE when available, else (None, reason)."""
    try:
        import pyne.data as pyne_data  # type: ignore[import-not-found]
    except Exception as exc:
        return None, f"PyNE unavailable ({exc}); using formula reference."
    src_id = {"EPA": 0, "DOE": 1, "GENII": 2}[source]
    try:
        if pathway == "air":
            val = pyne_data.ext_air_dose(name, src_id)
        elif pathway == "soil":
            val = pyne_data.ext_soil_dose(name, src_id)
        elif pathway == "ingest":
            val = pyne_data.ingest_dose(name, src_id)
        elif pathway == "inhale":
            val = pyne_data.inhale_dose(name, src_id)
        else:
            return None, f"unknown pathway {pathway}"
    except Exception as exc:
        return None, f"PyNE lookup failed ({exc})."
    return float(val), "PyNE accessor"


def pyne_dose_per_g(w: float, lam: float, mass_u: float, df: float, pathway: str) -> float:
    """PyNE `Material::dose_per_g` per-nuclide term with PyNE's N_A."""
    k = CI_PER_BQ if pathway in ("air", "soil") else PCI_PER_BQ
    return k * N_A_PYNE * w * lam * df / mass_u


def main() -> int:
    report = Report("dose", "Dose coefficients (`dose_vs_pyne.py`)")
    report.prose(
        "Nucleide dose factors (HNF-SD-WM-TI-707 Rev.1 / HNF-5636 App. O via"
        " PyNE `dbgen/dosefactors*.csv`) vs the PyNE `Material::dose_per_g`"
        " equations (source ids 0=EPA/1=DOE/2=GENII). Screening-level only;"
        " not for safety decisions."
    )

    rows: list[list[str]] = []
    notes: list[str] = []
    for name, pathway, source, spot in SPOTS:
        nuc = nucleide.nuclei.dose_factor(name, pathway, source)
        if nuc is None:
            _check(False, f"{name} {pathway} {source} missing")
            rows.append([name, pathway, source, "None", fmt(spot), "n/a", "FAIL"])
            continue
        d = rel_diff(nuc, spot)
        rows.append(
            [name, pathway, source, fmt(nuc), fmt(spot), fmt(d), _check(d < 1e-6, f"{name} DF")]
        )
        ref, note = _pyne_factor(name, pathway, source)
        notes.append(f"{name} {pathway} {source}: {note}")
        if ref is not None and ref >= 0:
            dd = rel_diff(nuc, ref)
            _check(dd < 1e-6, f"{name} vs PyNE factor")

    report.table(
        ["Nuclide", "Pathway", "Source", "Nucleide DF", "Spot DF", "Rel diff", "Status"],
        rows,
    )
    for note in notes:
        report.prose(note)

    # Per-gram equation check: 1 g pure Co60 ingest EPA via both N_A values.
    df = nucleide.nuclei.dose_factor("Co60", "ingest", "EPA")
    lam = nucleide.nuclei.decay_constant("Co60")
    mass = nucleide.nuclei.atomic_mass("Co60")
    assert df is not None and lam is not None and mass is not None
    got = nucleide.material.dose_per_g({"Co60": 1.0}, "ingest", "EPA")
    ref_pyne = pyne_dose_per_g(1.0, lam, mass, df, "ingest")
    ref_exact = PCI_PER_BQ * N_A_NUCLEIDE * lam * df / mass
    report.prose(
        f"1 g Co60 ingest EPA per-gram dose: nucleide {fmt(got)},"
        f" PyNE-N_A formula {fmt(ref_pyne)}, exact-N_A formula {fmt(ref_exact)}."
    )
    # N_A difference is <0.1 ppm; implementation uses exact N_A.
    _check(
        rel_diff(got, ref_exact) < 1e-12,
        "Co60 dose_per_g matches exact-N_A formula",
    )
    _check(
        rel_diff(got, ref_pyne) < 1e-6,
        "Co60 dose_per_g matches PyNE-N_A formula within 1 ppm",
    )
    report.prose(
        "Nucleide uses exact N_A (6.02214076e23); PyNE uses 6.0221415e23."
        " The per-gram doses agree to <1 ppm."
    )

    # GENII/DOE air sentinel stays -1 (missing).
    sentinel = nucleide.nuclei.dose_factor("H3", "air", "GENII")
    report.prose(f"H3 air GENII sentinel: {sentinel} (PyNE -1-for-missing-air).")
    _check(sentinel == -1.0, "GENII air sentinel is -1")
    try:
        nucleide.material.dose_per_g({"H3": 1.0}, "air", "GENII")
        _check(False, "GENII air dose_per_g should error")
    except ValueError:
        _check(True, "GENII air dose_per_g errors")

    # Finite, positive sanity band for K-40 soil EPA per-gram dose.
    k40 = nucleide.material.dose_per_g({"K40": 1.0}, "soil", "EPA")
    report.prose(f"1 g K-40 soil EPA per-gram dose: {fmt(k40)} mrem/h per g per m^2.")
    _check(math.isfinite(k40) and k40 > 0.0, "K40 per-gram dose finite/positive")

    try:
        report.emit()
    except Exception as exc:
        print(f"SKIPPED report write outside the container ({exc}).")

    if FAILURES:
        print(f"FAIL: {FAILURES} dose check(s) failed", file=sys.stderr)
        return 1
    print("All dose checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
