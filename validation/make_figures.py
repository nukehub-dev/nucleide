"""Generate paper figures from the validation results.

- ``figures/timings.png``: bar chart (log scale) of Nucleide vs reference-code
  wall times, read from ``results/timings.json``.
- ``figures/depletion_agreement.png``: log-log scatter of Nucleide vs OpenMC
  per-nuclide final densities with identity line plus a relative-difference
  residual strip per panel (``chain_ni.xml`` and the CASL/VERA chain).
- ``figures/kinetics_transient.png``: semilog step response of the 1-group
  point-kinetics transient (Nucleide vs closed-form analytic) with the
  prompt-jump level marked, plus a per-node relative-error residual strip,
  read from ``results/kinetics.json``.
- ``figures/spectroscopy_overlay.png``: raw vs smoothed counts with the ROI
  (channels 2..5) gross/background/net totals annotated, read from
  ``results/spectroscopy.json``.

Deterministic: fixed style, no timestamps or randomness. Figures are
generated artifacts — regenerate with ``run_all.sh``, never hand-edit.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR / "results"
FIGURES_DIR = SCRIPT_DIR / "figures"

_FLOAT_RE = re.compile(r"[0-9]+\.[0-9]+e[+-][0-9]+")

plt.rcParams.update(
    {
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.axisbelow": True,  # grid below bars/points, never striped through them
    }
)

# PNG metadata is kept constant so reruns on the same environment produce
# byte-stable content.
PNG_METADATA = {"Software": "nucleide validation/make_figures.py"}

COLOR_NUCLEIDE = "#1f77b4"
COLOR_NATIVE = "#2ca02c"
COLOR_REFERENCE = "#d62728"


def _seconds(cell: str) -> float:
    """Extract the ``%e`` wall time from a timings table cell."""
    m = _FLOAT_RE.search(cell)
    if not m:
        raise ValueError(f"no wall time found in {cell!r}")
    return float(m.group(0))


def _human_time(t: float) -> str:
    """Format seconds as a compact µs/ms/s label for bar annotations."""
    if t >= 1.0:
        return f"{t:.3g} s"
    if t >= 1e-3:
        return f"{t * 1e3:.3g} ms"
    return f"{t * 1e6:.3g} µs"


def timings_figure() -> Path:
    """Bar chart of wall times from results/timings.json."""
    data = json.loads((RESULTS_DIR / "timings.json").read_text())
    table = next(s for s in data["sections"] if s["kind"] == "table")
    rows = {row[0]: row for row in table["rows"]}

    # Each group: (group label, [(bar kind, bar annotation, seconds), ...]).
    # Bar kinds: "nucleide", "native", "reference". The reference annotation
    # names the concrete code, since it differs per group.
    groups = [
        (
            "CRAM-48 solve\n(chain_ni.xml)",
            [
                ("nucleide", "", _seconds(rows["CRAM-48 solve (`chain_ni.xml`)"][1])),
                ("reference", "OpenMC", _seconds(rows["CRAM-48 solve (`chain_ni.xml`)"][2])),
                (
                    "native",
                    "",
                    _seconds(rows["Native Rust CRAM-48 solve (Criterion, no Python)"][1]),
                ),
            ],
        ),
        (
            "U enrichment solve\n(multicomponent)",
            [
                ("nucleide", "", _seconds(rows["Default uranium enrichment solve"][1])),
                ("reference", "PyNE", _seconds(rows["Default uranium enrichment solve"][2])),
            ],
        ),
        (
            "MAGIC solve\n(synthetic tally)",
            [
                ("nucleide", "", _seconds(rows["MAGIC total-mode solve (synthetic tally)"][1])),
                (
                    "reference",
                    "pure Python",
                    _seconds(rows["MAGIC total-mode solve (synthetic tally)"][2]),
                ),
            ],
        ),
    ]

    colors = {
        "nucleide": COLOR_NUCLEIDE,
        "native": COLOR_NATIVE,
        "reference": COLOR_REFERENCE,
    }
    legend_labels = {
        "nucleide": "Nucleide (Python)",
        "native": "Nucleide (native Rust)",
        "reference": "reference code",
    }

    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    width = 0.8 / 3
    seen: set[str] = set()
    ymin = min(value for _, bars in groups for _, _, value in bars) / 3
    for g, (_, bars) in enumerate(groups):
        for k, (kind, ref_name, value) in enumerate(bars):
            x = g + (k - (len(bars) - 1) / 2) * width * 1.1
            ax.bar(
                x,
                value,
                width=width,
                color=colors[kind],
                label=legend_labels[kind] if kind not in seen else None,
            )
            seen.add(kind)
            # Wall time above the bar; reference-code name inside the bar, or
            # above it when the bar is too short to hold vertical text.
            ax.annotate(
                _human_time(value),
                (x, value),
                xytext=(0, 2),
                textcoords="offset points",
                ha="center",
                fontsize=7.5,
            )
            if ref_name:
                short = np.log10(value / ymin) < 1.6
                ax.annotate(
                    ref_name,
                    (x, value),
                    xytext=(0, 14 if short else -2),
                    textcoords="offset points",
                    ha="center",
                    va="top" if not short else "bottom",
                    fontsize=7.5,
                    color="#888888" if short else "white",
                    rotation=0 if short else 90,
                )
    ax.set_yscale("log")
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([g[0] for g in groups])
    ax.set_ylabel("mean wall time [s]")
    ax.set_title("Wall time per solve (validation container, 20 repeats)")
    # Headroom so the legend and bar labels clear the tallest bars.
    top = max(value for _, bars in groups for _, _, value in bars) * 120
    ax.set_ylim(top=top)
    ax.legend(loc="upper left", framealpha=0.9)
    ax.grid(axis="y", which="major", alpha=0.4)
    ax.grid(axis="y", which="minor", alpha=0.15)

    out = FIGURES_DIR / "timings.png"
    fig.savefig(out, metadata=PNG_METADATA)
    plt.close(fig)
    return out


def depletion_agreement_figure() -> Path:
    """Identity scatter + relative-difference residuals for both chains."""
    sys.path.insert(0, str(SCRIPT_DIR))
    from depletion_casl_vs_openmc import run_casl_vectors
    from depletion_vs_openmc import run_chain_ni_vectors

    panels = [
        ("chain_ni.xml (Ni-58 activation)", run_chain_ni_vectors(), "depletion.json"),
        ("CASL/VERA chain (fresh UO2)", run_casl_vectors(), "depletion_casl.json"),
    ]

    fig = plt.figure(figsize=(8.6, 5.2))
    grid = fig.add_gridspec(2, 2, height_ratios=[3, 1], hspace=0.08, wspace=0.25)
    for col, (title, (nuc_out, om_out), json_name) in enumerate(panels):
        names = sorted(nuc_out)
        x = np.array([om_out.get(n, 0.0) for n in names])
        y = np.array([nuc_out[n] for n in names])
        mask = (x > 0.0) & (y > 0.0)
        x, y = x[mask], y[mask]
        # Same density threshold as the comparison scripts: exclude nuclides
        # whose density is negligible relative to the chain maximum.
        max_dens = float(max(x.max(), y.max()))
        keep = np.maximum(x, y) >= 1.0e-12 * max_dens
        x, y = x[keep], y[keep]
        rel = np.abs(y - x) / np.maximum(x, y)
        max_rel = float(rel.max())

        # Sanity: the figure's number must match the reported JSON value.
        report = json.loads((RESULTS_DIR / json_name).read_text())
        report_table = next(s for s in report["sections"] if s["kind"] == "table")
        reported = next(
            float(_FLOAT_RE.search(row[1]).group(0))
            for row in report_table["rows"]
            if "Max relative difference" in row[0]
        )
        if abs(max_rel - reported) > 0.01 * reported:
            raise ValueError(
                f"{json_name}: figure max rel {max_rel:.3e} != reported {reported:.3e}"
            )

        ax = fig.add_subplot(grid[0, col])
        lo = min(x.min(), y.min())
        hi = max(x.max(), y.max())
        ax.loglog([lo, hi], [lo, hi], "k-", lw=0.8, label="identity")
        ax.loglog(x, y, ".", ms=4, alpha=0.6, color=COLOR_NUCLEIDE, label=f"{mask.sum()} nuclides")
        ax.annotate(
            f"max rel diff = {max_rel:.1e}",
            (0.03, 0.97),
            xycoords="axes fraction",
            va="top",
            fontsize=8,
        )
        ax.set_title(title)
        ax.legend(loc="lower right")
        ax.grid(which="both", alpha=0.3)
        if col == 0:
            ax.set_ylabel("Nucleide final density [atoms]")
        ax.tick_params(labelbottom=False)

        axr = fig.add_subplot(grid[1, col], sharex=ax)
        axr.loglog(x, np.maximum(rel, 1e-18), ".", ms=3, alpha=0.6, color=COLOR_REFERENCE)
        axr.set_ylim(1e-17, 1e-13)
        axr.grid(which="both", alpha=0.3)
        axr.set_xlabel("OpenMC final density [atoms]")
        if col == 0:
            axr.set_ylabel("rel. diff.", fontsize=8)

    fig.suptitle("Per-nuclide final densities after one 30-day CRAM-48 step, Nucleide vs OpenMC")

    out = FIGURES_DIR / "depletion_agreement.png"
    fig.savefig(out, metadata=PNG_METADATA)
    plt.close(fig)
    return out


def _table_by_headers(data: dict, headers: list[str]) -> list[list[str]]:
    """Return the rows of the table section with exactly these headers."""
    for section in data["sections"]:
        if section["kind"] == "table" and section["headers"] == headers:
            return section["rows"]
    raise ValueError(f"no table with headers {headers!r}")


def _quantity(data: dict, label: str) -> float:
    """Return the float value of a ``["Quantity", "Value"]`` row by label."""
    for row in _table_by_headers(data, ["Quantity", "Value"]):
        if row[0] == label:
            return float(row[1])
    raise ValueError(f"no quantity row {label!r}")


def _gate_rel_err(data: dict, gate: str) -> float:
    """Return the ``Rel err`` float of a ``Gate`` table row by label."""
    for row in _table_by_headers(data, ["Gate", "Rel err", "Tol", "Status"]):
        if row[0] == gate:
            return float(row[1])
    raise ValueError(f"no gate row {gate!r}")


def kinetics_transient_figure() -> Path:
    """Step response of the 1-group transient (Nucleide vs analytic + residuals)."""
    data = json.loads((RESULTS_DIR / "kinetics.json").read_text())
    rows = _table_by_headers(data, ["t since step (s)", "Nucleide n", "Analytic n"])
    t = np.array([float(r[0]) for r in rows])
    nuc = np.array([float(r[1]) for r in rows])
    ana = np.array([float(r[2]) for r in rows])
    pj_level = _quantity(data, "Prompt-jump level (n0 = 1)")
    # Residuals are arithmetic on the two plotted series (same convention as
    # the depletion_agreement.py lower strips).
    rel = np.abs(nuc - ana) / np.abs(ana)

    # Sanity: the figure's number must match the reported JSON gate value.
    worst = float(np.max(rel))
    reported = _gate_rel_err(data, "O3 1-group transient")
    if abs(worst - reported) > 0.01 * max(reported, 1e-30):
        raise ValueError(f"kinetics.json: figure worst rel {worst:.3e} != {reported:.3e}")
    exp = int(np.floor(np.log10(worst)))
    mant = worst / 10.0**exp

    # Title parameters (rho step 0 -> 0.002, n0 = 1) are the recorded values
    # in fixtures/kinetics/step_oracle.json.
    fig = plt.figure(figsize=(6.4, 4.2))
    grid = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.08)
    ax = fig.add_subplot(grid[0])
    ax.semilogx(t, nuc, "o-", ms=4, color=COLOR_NUCLEIDE, label="Nucleide n")
    ax.semilogx(t, ana, "k--", lw=1.0, label="Analytic n")
    ax.axhline(pj_level, color=COLOR_REFERENCE, ls="--", lw=1.0, label="Prompt-jump level")
    ax.annotate(
        f"worst rel. err. vs analytic: ${mant:.1f}\\times10^{{{exp}}}$",
        (0.03, 0.95),
        xycoords="axes fraction",
        va="top",
        fontsize=8,
    )
    ax.set_ylabel("n(t)")
    ax.set_title("1-group step transient (ρ step 0 → 0.002, n₀ = 1): Nucleide vs analytic")
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), framealpha=0.9)
    ax.grid(which="both", alpha=0.3)
    ax.tick_params(labelbottom=False)

    axr = fig.add_subplot(grid[1], sharex=ax)
    axr.loglog(t, np.maximum(rel, 1e-18), ".", ms=5, color=COLOR_REFERENCE)
    axr.set_ylim(float(rel[rel > 0].min()) / 2, worst * 5)
    axr.set_xlabel("t since step [s]")
    axr.set_ylabel("rel. err.", fontsize=8)
    axr.grid(which="both", alpha=0.3)

    out = FIGURES_DIR / "kinetics_transient.png"
    fig.savefig(out, metadata=PNG_METADATA)
    plt.close(fig)
    return out


def spectroscopy_overlay_figure() -> Path:
    """Raw vs smoothed counts with the ROI-totals annotation."""
    data = json.loads((RESULTS_DIR / "spectroscopy.json").read_text())
    rows = _table_by_headers(
        data, ["Channel", "Raw counts", "Rect-smoothed (m=5)", "Five-point smoothed"]
    )
    channels = np.array([float(r[0]) for r in rows])
    raw = np.array([float(r[1]) for r in rows])
    rect = np.array([float(r[2]) for r in rows])
    five = np.array([float(r[3]) for r in rows])
    bg = _quantity(data, "Background level (E3, channels 2..5)")

    # Sanity: overlay must cover exactly the 7-channel oracle vector.
    if len(rows) != 7 or not np.all(channels == np.arange(7.0)):
        raise ValueError("spectroscopy.json: overlay is not the 7-channel oracle vector")
    for gate in ("E1 rect smooth", "E2 five-point smooth"):
        if _gate_rel_err(data, gate) != 0.0:
            raise ValueError(f"spectroscopy.json: {gate} is not an exact match")

    # ROI totals over channels 2..5 (half-open [2, 5), as in the E4 gate):
    # gross sums the plotted raw series, net subtracts the JSON background
    # total. These equal the recorded hand values checked by
    # spectroscopy_vs_pyne.py (gross 10.0, background 76/6, net -16/6).
    gross = float(raw[2] + raw[3] + raw[4])
    net = gross - bg

    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    ax.plot(channels, raw, "ks", ms=5, label="Raw counts")
    ax.plot(channels, rect, "o-", ms=4, color=COLOR_NUCLEIDE, label="Rect-smoothed (m=5)")
    ax.plot(channels, five, "^-", ms=4, color=COLOR_NATIVE, label="Five-point smoothed")
    ax.text(
        0.02,
        0.97,
        f"ROI totals (channels 2..5):\ngross = {gross:.1f}, "
        f"background = {bg:.4f},\nnet = {net:.4f}\n"
        "smoothers match hand values exactly",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.9},
    )
    ax.set_xlabel("channel")
    ax.set_ylabel("counts")
    ax.set_title("Smoothing overlay on the 7-channel oracle vector")
    ax.set_xticks(channels)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), framealpha=0.9)
    ax.grid(which="both", alpha=0.3)

    out = FIGURES_DIR / "spectroscopy_overlay.png"
    fig.savefig(out, metadata=PNG_METADATA)
    plt.close(fig)
    return out


def main() -> int:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    t = timings_figure()
    print(f"Wrote {t}")
    d = depletion_agreement_figure()
    print(f"Wrote {d}")
    k = kinetics_transient_figure()
    print(f"Wrote {k}")
    s = spectroscopy_overlay_figure()
    print(f"Wrote {s}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
