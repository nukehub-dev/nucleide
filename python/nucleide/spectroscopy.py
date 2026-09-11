"""Gamma-ray spectroscopy and measurement (backed by the `nucleide-spectroscopy` crate)."""

from typing import Any

from nucleide._internal import (
    spectroscopy_calc_bg,
    spectroscopy_detector_efficiency,
    spectroscopy_energy_bins,
    spectroscopy_five_point_smooth,
    spectroscopy_gross_count,
    spectroscopy_net_counts,
    spectroscopy_parse_dollar_spe,
    spectroscopy_parse_spe,
    spectroscopy_read_dollar_spe,
    spectroscopy_read_spe,
    spectroscopy_rect_smooth,
    spectroscopy_xray_lines,
)

__all__ = [
    "rect_smooth",
    "five_point_smooth",
    "calc_bg",
    "gross_count",
    "net_counts",
    "energy_bins",
    "detector_efficiency",
    "xray_lines",
    "parse_dollar_spe",
    "parse_spe",
    "read_dollar_spe",
    "read_spe",
]


def rect_smooth(counts: list[float], m: int) -> list[float]:
    """Rectangular smoothing: window mean over ``m`` (odd, >= 3), edges copied."""
    return spectroscopy_rect_smooth(counts, m)


def five_point_smooth(counts: list[float]) -> list[float]:
    """Five-point smoothing (low-statistics weights); first/last two copied."""
    return spectroscopy_five_point_smooth(counts)


def calc_bg(counts: list[float], channels: list[float], c1: int, c2: int, m: int) -> float:
    """Background under a peak (``m == 1`` only); bounds index ``counts`` positionally."""
    return spectroscopy_calc_bg(counts, channels, c1, c2, m)


def gross_count(counts: list[float], channels: list[float], c1: int, c2: int) -> float:
    """Gross counts between two channels, half-open (excludes ``c2``)."""
    return spectroscopy_gross_count(counts, channels, c1, c2)


def net_counts(counts: list[float], channels: list[float], c1: int, c2: int, m: int) -> float:
    """Net counts: ``gross_count`` minus ``calc_bg``."""
    return spectroscopy_net_counts(counts, channels, c1, c2, m)


def energy_bins(channels: list[float], calib_e_fit: list[float]) -> list[float]:
    """Energy per channel from the ``[a0, a1, a2]`` quadratic fit."""
    return spectroscopy_energy_bins(channels, calib_e_fit)


def detector_efficiency(energy_mev: float, eff_coeff: list[float], eff_fit: int = 1) -> float:
    """Detector efficiency at ``energy_mev`` [MeV] (``eff_fit`` 1 or 2)."""
    return spectroscopy_detector_efficiency(energy_mev, eff_coeff, eff_fit)


def xray_lines(
    atomic: dict[str, float],
    k_conv: float | None = None,
    l_conv: float | None = None,
) -> list[tuple[float, float]]:
    """X-ray lines as ``[(energy_kev, intensity)]`` (Ka1, Ka2, Kb, L).

    ``atomic`` carries the nine caller-supplied constants
    (``k_shell_fluor``, ``l_shell_fluor``, ``prob``, ``kb_to_ka``,
    ``ka2_to_ka1``, ``ka1_en_kev``, ``ka2_en_kev``, ``kb_en_kev``,
    ``l_en_kev``); no atomic table is vendored. ``None`` (or NaN, the
    upstream sentinel) marks a conversion absent. Upstream exposes no
    combined function for this routine — only a material method — so this
    explicit entry point is the documented Nucleide surface.
    """
    return spectroscopy_xray_lines(atomic, k_conv, l_conv)


def parse_dollar_spe(text: str) -> dict[str, Any]:
    """Parse dollar-format ``.spe`` text (first line must be ``$SPEC_ID:``)."""
    return spectroscopy_parse_dollar_spe(text)


def parse_spe(text: str) -> dict[str, Any]:
    """Parse plain-format ``.spe`` text (rejects the ``$SPEC_ID:`` magic)."""
    return spectroscopy_parse_spe(text)


def read_dollar_spe(path: str) -> dict[str, Any]:
    """Read a dollar-format ``.spe`` file."""
    return spectroscopy_read_dollar_spe(path)


def read_spe(path: str) -> dict[str, Any]:
    """Read a plain-format ``.spe`` file."""
    return spectroscopy_read_spe(path)
