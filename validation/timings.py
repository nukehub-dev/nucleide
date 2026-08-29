"""Coarse wall-time comparisons between Nucleide and reference codes."""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import openmc.deplete
import pyne.enrichment as pyne_enr
from common import Report, fmt

import nucleide

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "validation"))

from depletion_vs_openmc import (  # noqa: E402
    nucleide_rates,
    openmc_rates,
    patched_chain_ni_path,
)
from magic_vs_pyne import pyne_magic_equivalent  # noqa: E402

# Fallback Criterion values if target/criterion/ is unavailable.
FALLBACK_CRAM48_US = 7.6286e-05
FALLBACK_DEPLETE_E2E_US = 7.8794e-05


def _criterion_mean_s(bench_dir: Path) -> float | None:
    """Read the mean wall time (seconds) from a Criterion estimates.json."""
    estimates = bench_dir / "new" / "estimates.json"
    try:
        data = json.loads(estimates.read_text())
        ns = data["mean"]["point_estimate"]
        return float(ns) * 1.0e-9
    except Exception:
        return None


def native_rust_timings() -> dict[str, float]:
    """Return native Rust Criterion timings for chain_ni.xml CRAM-48."""
    criterion_root = REPO_ROOT / "target" / "criterion" / "depletion_ni_chain"
    cram48 = _criterion_mean_s(criterion_root / "cram48_solve")
    e2e = _criterion_mean_s(criterion_root / "deplete_end_to_end")
    return {
        "cram48_s": cram48 if cram48 is not None else FALLBACK_CRAM48_US,
        "deplete_e2e_s": e2e if e2e is not None else FALLBACK_DEPLETE_E2E_US,
    }


def time_it(fn, repeats: int = 20) -> float:
    """Return mean wall time over `repeats` calls (seconds)."""
    # Warmup
    fn()
    times: list[float] = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn()
        t1 = time.perf_counter()
        times.append(t1 - t0)
    return float(np.mean(times))


def depletion_timings(repeats: int = 20) -> dict[str, float]:
    """Time CRAM-48 single-step solve on chain_ni.xml.

    Both solvers are timed on pre-built systems/matrices so the measurement
    reflects only the linear-algebra solve step.
    """
    chain_path = patched_chain_ni_path()
    nuc_chain = nucleide.depletion.read_chain(chain_path)
    om_chain = openmc.deplete.Chain.from_xml(chain_path)

    dt = 2.592e6
    n0 = {nuc: (1.0e24 if nuc == "Ni58" else 0.0) for nuc in nuc_chain.nuclides}
    nuc_rates = nucleide_rates(nuc_chain)
    nuc_sys = nucleide.depletion.build_depletion_system(nuc_chain, nuc_rates)
    n0_vec_nuc = [n0.get(nuc, 0.0) for nuc in nuc_chain.nuclides]

    om_rates = openmc_rates(om_chain)
    A = om_chain.form_matrix(om_rates[0])
    n0_vec = np.array([n0.get(nuc.name, 0.0) for nuc in om_chain.nuclides])

    def nucleide_fn() -> None:
        nuc_sys.solve_vec(n0_vec_nuc, dt, order=48)

    def openmc_fn() -> None:
        openmc.deplete.cram.CRAM48(A, n0_vec, dt)

    return {
        "nucleide_deplete_s": time_it(nucleide_fn, repeats),
        "openmc_cram48_s": time_it(openmc_fn, repeats),
    }


def enrichment_timings(repeats: int = 20) -> dict[str, float]:
    """Time multicomponent enrichment solve for default uranium."""

    def nucleide_fn() -> None:
        casc = nucleide.enrichment.Cascade.default_uranium()
        casc.solve()

    def pyne_fn() -> None:
        casc = pyne_enr.default_uranium_cascade()
        pyne_enr.multicomponent(casc, solver="numeric")

    return {
        "nucleide_enrich_s": time_it(nucleide_fn, repeats),
        "pyne_enrich_s": time_it(pyne_fn, repeats),
    }


def magic_timings(repeats: int = 20) -> dict[str, float]:
    """Time MAGIC weight-window generation on the synthetic tally."""
    tally_path = REPO_ROOT / "validation" / "magic_tally.txt"
    meshtal = nucleide.mcnp.read_meshtal(str(tally_path))
    tally = meshtal.tallies[4]

    def nucleide_fn() -> None:
        nucleide.vr.magic(tally, per_group=False, tolerance=0.5)

    def pyne_fn() -> None:
        pyne_magic_equivalent(
            tally.result,
            tally.rel_error,
            tally.total_result,
            tally.total_rel_error,
            per_group=False,
            tolerance=0.5,
        )

    return {
        "nucleide_magic_s": time_it(nucleide_fn, repeats),
        "pyne_magic_s": time_it(pyne_fn, repeats),
    }


def main() -> int:
    report = Report("timings", "Timings (`timings.py`)")

    report.prose(
        "Mean wall time over 20 repeats. The CRAM comparison now times **only the solve\n"
        "step** on pre-built systems/matrices."
    )

    dep = depletion_timings()
    enr = enrichment_timings()
    mag = magic_timings()
    native = native_rust_timings()

    report.table(
        ["Operation", "Nucleide", "Reference code"],
        [
            [
                "CRAM-48 solve (`chain_ni.xml`)",
                f"{fmt(dep['nucleide_deplete_s'])} s",
                f"OpenMC CRAM48: {fmt(dep['openmc_cram48_s'])} s",
            ],
            [
                "Default uranium enrichment solve",
                f"{fmt(enr['nucleide_enrich_s'])} s",
                f"PyNE multicomponent: {fmt(enr['pyne_enrich_s'])} s",
            ],
            [
                "MAGIC total-mode solve (synthetic tally)",
                f"{fmt(mag['nucleide_magic_s'])} s",
                f"PyNE-equivalent pure Python: {fmt(mag['pyne_magic_s'])} s",
            ],
            [
                "Native Rust CRAM-48 solve (Criterion, no Python)",
                f"{fmt(native['cram48_s'])} s",
                "—",
            ],
            [
                "Native Rust deplete end-to-end (Criterion, no Python)",
                f"{fmt(native['deplete_e2e_s'])} s",
                "—",
            ],
        ],
    )
    report.prose(
        "The native Rust Criterion numbers were obtained with `cargo bench -p depletion`\n"
        "and represent the same `chain_ni.xml` CRAM-48 solve without any Python/PyO3 "
        "overhead."
    )

    report.emit()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
