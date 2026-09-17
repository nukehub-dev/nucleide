"""IAEA IRDFF-II dosimetry response pack (synthetic data, fully offline).

Nothing IAEA-copyright is touched here except through the live-fetch gate:
unit tests build synthetic zips/sections in the exact ``IRDFF-II.g725``
``MF=3`` layout and monkeypatch ``nucleide.data.IRDFF_SHA256`` to the
synthetic digest — ``fetch_irdff`` enforces the pinned hash on every call,
so the override is how a synthetic "release" is pinned.

``TestLiveFetch`` downloads the official IAEA zip (hash-verified by
``fetch_irdff`` itself) and parses the full v1+v2 registry through the
landed ``unfold`` iterators. It needs network access to the IAEA site, so
it raises a loud ``SKIP`` when the download fails; a hash mismatch still
raises loudly and never skips.
"""

import hashlib
import math
import zipfile
from pathlib import Path
from typing import Any

import pytest

import nucleide.data
import nucleide.unfold as unfold
from nucleide._internal import parse_irdff_g725
from nucleide.data import IRDFF_MEMBER

#: The v1 foil-activation subset, in registry order (pins the facade default).
V1_NAMES = [
    "au197_ng",
    "in115_ng",
    "u235_nf",
    "u238_nf",
    "fe56_np",
    "ni58_np",
    "al27_na",
    "na23_n2n",
]

#: The v2 extension set, in registry order (pins the expansion list).
V2_NAMES = [
    "f19_n2n",
    "b10_na",
    "mg24_np",
    "al27_np",
    "si28_np",
    "p31_np",
    "s32_np",
    "sc45_ng",
    "ti46_np",
    "ti47_np",
    "ti48_np",
    "mn55_n2n",
    "fe54_np",
    "co59_ng",
    "co59_np",
    "ni58_n2n",
    "cu63_na",
    "zn64_np",
    "in113_ng",
    "ta181_ng",
    "w186_ng",
    "th232_nf",
    "np237_nf",
    "pu239_nf",
    "bi209_n2n",
    "bi209_n3n",
]

#: Synthetic (MAT, MT, ZA) keys mirroring the real registry shape only.
SYNTH_KEYS = {
    "au197_ng": (7925, 102, 79197.0),
    "in115_ng": (4931, 102, 49115.0),
    "u235_nf": (9228, 18, 92235.0),
    "u238_nf": (9237, 18, 92238.0),
    "fe56_np": (2631, 103, 26056.0),
    "ni58_np": (2825, 103, 28058.0),
    "al27_na": (1325, 107, 13027.0),
    "na23_n2n": (1125, 16, 11023.0),
    "sc45_ng": (2125, 102, 21045.0),
    "co59_ng": (2725, 102, 27059.0),
    "ti46_np": (2225, 103, 22046.0),
    "mn55_n2n": (2525, 16, 25055.0),
    "pu239_nf": (9437, 18, 94239.0),
}

N_GROUPS = 725


def synth_bounds() -> list[float]:
    """Synthetic 726-boundary structure with the pinned shape (never IAEA data)."""
    lo, hi = 1e-5, 6e7
    step = math.log(hi / lo) / N_GROUPS
    return [lo * math.exp(step * g) for g in range(N_GROUPS + 1)]


def fmt11(x: float) -> str:
    """One eleven-column ENDF-style field (embedded exponents below 1.0)."""
    if x == 0.0:
        s = "0.0"
    elif abs(x) < 1.0:
        mantissa, exp = f"{x:.6e}".split("e")
        s = f"{mantissa}{'-' if exp.startswith('-') else '+'}{exp[1:].lstrip('0') or '0'}"
    elif abs(x) < 1e5:
        s = f"{x:.5f}"
    else:
        s = f"{x:.1f}"
    return f"{s:>11}"


def section_text(name: str, energies: list[float], values: list[float], seq0: int = 1) -> str:
    """One synthetic ``MF=3`` section in the exact distributed layout."""
    assert len(energies) == len(values)
    mat, mt, za = SYNTH_KEYS[name]
    zero = f"{0:>11}"
    lines = [
        f"{fmt11(za)}{fmt11(100.0)}{zero}{zero}{zero}{zero}{mat:>4} 3{mt:>3}{seq0:>5}",
        f"{fmt11(0.0)}{fmt11(0.0)}{zero}{zero}{1:>11}{len(energies):>11}"
        f"{mat:>4} 3{mt:>3}{seq0 + 1:>5}",
        f"{len(energies):>11}{1:>11}{mat:>4} 3{mt:>3}{seq0 + 2:>5}",
    ]
    for chunk in range(0, len(energies), 3):
        fields = ""
        for k in range(chunk, min(chunk + 3, len(energies))):
            fields += fmt11(energies[k]) + fmt11(values[k])
        fields += " " * (22 * (3 - (min(chunk + 3, len(energies)) - chunk)))
        lines.append(f"{fields}{mat:>4} 3{mt:>3}{seq0 + 3 + chunk // 3:>5}")
    return "\n".join(lines) + "\n"


def full_values(sigma: float) -> list[float]:
    """Constant full-range row values with the structural terminator slot."""
    return [sigma] * N_GROUPS + [0.0]


def synth_pack_text() -> str:
    """Synthetic pack: Au/Sc/Co full-range captures + Ti/Mn threshold rows."""
    bounds = synth_bounds()
    text = section_text("au197_ng", bounds, full_values(2.0))
    text += section_text("sc45_ng", bounds, full_values(3.0))
    text += section_text("co59_ng", bounds, full_values(1.5))
    ti_energies = bounds[500:]
    ti_values = [5.0] * (len(ti_energies) - 1) + [0.0]
    text += section_text("ti46_np", ti_energies, ti_values)
    mn_energies = bounds[589:]
    mn_values = [4.0] * (len(mn_energies) - 1) + [0.0]
    text += section_text("mn55_n2n", mn_energies, mn_values)
    return text


def make_pinned_zip(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Build a synthetic zip and pin it as the expected IRDFF-II download."""
    zip_path = tmp_path / "irdff_synthetic.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr(IRDFF_MEMBER, "synthetic")
    digest = hashlib.sha256(zip_path.read_bytes()).hexdigest()
    monkeypatch.setattr(nucleide.data, "IRDFF_SHA256", digest)
    return zip_path


def fetch_kwargs(tmp_path: Path, zip_path: Path) -> dict[str, Any]:
    return {"dest": tmp_path / "cache", "url": zip_path.as_uri()}


class TestFetch:
    def test_download_then_cache_hit(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        zip_path = make_pinned_zip(tmp_path, monkeypatch)
        got = nucleide.data.fetch_irdff(**fetch_kwargs(tmp_path, zip_path))
        assert got == str(tmp_path / "cache" / "irdff_synthetic.zip")
        # A verified cache is reused without touching the network: the second
        # call's URL is unreachable and must never be opened.
        again = nucleide.data.fetch_irdff(
            dest=tmp_path / "cache", url="http://127.0.0.1:59999/irdff_synthetic.zip"
        )
        assert again == got

    def test_hash_mismatch_names_hashes_and_cleans_up(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        zip_path = make_pinned_zip(tmp_path, monkeypatch)
        monkeypatch.setattr(nucleide.data, "IRDFF_SHA256", "0" * 64)
        with pytest.raises(RuntimeError, match="sha256 mismatch.*expected"):
            nucleide.data.fetch_irdff(**fetch_kwargs(tmp_path, zip_path))
        assert not (tmp_path / "cache" / "irdff_synthetic.zip").exists()

    def test_corrupt_cache_is_a_loud_mismatch(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        zip_path = make_pinned_zip(tmp_path, monkeypatch)
        dest = tmp_path / "cache"
        dest.mkdir()
        (dest / "irdff_synthetic.zip").write_bytes(b"corrupt")
        with pytest.raises(RuntimeError, match="sha256 mismatch"):
            nucleide.data.fetch_irdff(dest=dest, url=zip_path.as_uri())

    def test_network_failure_without_cache(self, tmp_path: Path) -> None:
        with pytest.raises(RuntimeError, match="failed to download.*No usable cached copy"):
            nucleide.data.fetch_irdff(
                dest=tmp_path / "cache",
                url="http://127.0.0.1:59999/irdff_synthetic.zip",
            )


class TestParseBinding:
    def test_explicit_v1_and_v2_names_in_request_order(self) -> None:
        pack = parse_irdff_g725(synth_pack_text(), ["mn55_n2n", "au197_ng", "sc45_ng", "ti46_np"])
        assert pack["reactions"] == ["mn55_n2n", "au197_ng", "sc45_ng", "ti46_np"]
        assert len(pack["groups"]) == N_GROUPS + 1
        assert all(len(row) == N_GROUPS for row in pack["response"])
        # Constant rows round-trip exactly; threshold rows are zero below group 500.
        assert pack["response"][1] == [2.0] * N_GROUPS
        assert pack["response"][3][:500] == [0.0] * 500
        assert pack["response"][3][500:] == [5.0] * (N_GROUPS - 500)

    def test_default_none_parses_the_v1_pack(self) -> None:
        bounds = synth_bounds()
        text = "".join(section_text(name, bounds, full_values(1.0)) for name in V1_NAMES)
        pack = parse_irdff_g725(text)
        assert pack["reactions"] == V1_NAMES
        assert all(row == [1.0] * N_GROUPS for row in pack["response"])

    def test_unknown_name_is_loud(self) -> None:
        with pytest.raises(ValueError, match="unknown IRDFF-II reaction `in115_nm`"):
            parse_irdff_g725(synth_pack_text(), ["in115_nm"])

    def test_missing_section_names_mat_mt(self) -> None:
        with pytest.raises(ValueError, match="no MF=3 section for pu239_nf"):
            parse_irdff_g725(synth_pack_text(), ["pu239_nf"])


class TestFolds:
    def test_synthetic_rows_fold_by_hand(self) -> None:
        pack = parse_irdff_g725(synth_pack_text(), ["au197_ng", "ti46_np", "mn55_n2n"])
        spectrum = [1.0] * N_GROUPS
        rates = unfold.forward_fold(pack["response"], spectrum)
        # Hand-checked: constant rows sum their support exactly.
        assert rates == [2.0 * N_GROUPS, 5.0 * (N_GROUPS - 500), 4.0 * (N_GROUPS - 589)]

    def test_exact_guess_is_a_fixed_point(self) -> None:
        pack = parse_irdff_g725(synth_pack_text(), ["au197_ng", "sc45_ng", "co59_ng"])
        truth = [1.0 + (g % 7) for g in range(N_GROUPS)]
        rates = unfold.forward_fold(pack["response"], truth)
        out = unfold.sandii(pack["response"], rates, list(truth))
        assert out["iterations"] == 1
        assert out["rate_factors"] == pytest.approx([1.0] * 3, rel=1e-12)


def _live_zip(tmp_path: Path) -> str:
    """Fetch the official IAEA zip, or loud-SKIP when offline (never on mismatch)."""
    try:
        return nucleide.data.fetch_irdff(dest=tmp_path / "cache")
    except RuntimeError as exc:
        if "failed to download" in str(exc):
            pytest.skip(f"IRDFF-II live-fetch gate needs network access: {exc}")
        raise


def _first_nonzero(row: list[float]) -> int:
    return next(i for i, v in enumerate(row) if v > 0.0)


class TestLiveFetch:
    def test_full_registry_structure_and_hand_checked_orderings(self, tmp_path: Path) -> None:
        zip_path = _live_zip(tmp_path)
        with zipfile.ZipFile(zip_path) as archive:
            text = archive.read(IRDFF_MEMBER).decode("utf-8", errors="replace")
        wanted = V1_NAMES + V2_NAMES
        pack = parse_irdff_g725(text, wanted)
        assert pack["reactions"] == wanted
        groups = pack["groups"]
        assert len(groups) == N_GROUPS + 1
        assert all(b < c for b, c in zip(groups, groups[1:], strict=True))
        assert groups[0] == pytest.approx(1e-5, rel=1e-9)
        assert groups[-1] == pytest.approx(6e7, rel=1e-9)
        rows = dict(zip(pack["reactions"], pack["response"], strict=True))
        for name, row in rows.items():
            assert len(row) == N_GROUPS, name
            assert all(math.isfinite(v) and v >= 0.0 for v in row), name
        # Captures and fission chambers respond at thermal energies;
        # threshold rows are exactly zero there.
        for name in ("au197_ng", "sc45_ng", "co59_ng", "u235_nf", "pu239_nf"):
            assert rows[name][0] > 0.0, name
        for name in ("fe56_np", "ni58_np", "al27_na", "na23_n2n", "ti46_np"):
            assert rows[name][0] == 0.0, name
        # Hand-checked first-response groups on the pinned file: higher
        # thresholds light up later, in this exact order.
        first = {name: _first_nonzero(rows[name]) for name in V2_NAMES}
        assert first["bi209_n2n"] == 559
        assert first["mn55_n2n"] == 589
        assert first["f19_n2n"] == 594
        assert first["ni58_n2n"] == 609
        assert first["bi209_n3n"] == 629
        assert first["bi209_n2n"] < first["mn55_n2n"] < first["ni58_n2n"]
        assert first["bi209_n2n"] < first["bi209_n3n"]

    def test_fetched_rows_drive_the_unfold_iterators(self, tmp_path: Path) -> None:
        zip_path = _live_zip(tmp_path)
        with zipfile.ZipFile(zip_path) as archive:
            text = archive.read(IRDFF_MEMBER).decode("utf-8", errors="replace")
        pack = parse_irdff_g725(text, V1_NAMES + V2_NAMES)
        lo, hi = 1e-5, 6e7
        step = math.log(hi / lo) / N_GROUPS
        midpoints = [lo * math.exp(step * (g + 0.5)) for g in range(N_GROUPS)]
        truth = [math.exp(-(((math.log(e) + 5.0) / 3.0) ** 2)) + 1e-9 for e in midpoints]
        rates = unfold.forward_fold(pack["response"], truth)
        assert all(r > 0.0 for r in rates)
        # Exact guess is a fixed point: fetched rows reproduce the fold.
        out = unfold.sandii(pack["response"], rates, list(truth))
        assert out["iterations"] == 1
        assert out["rate_factors"] == pytest.approx([1.0] * len(rates), rel=1e-12)
        # Biased guess still recovers the rates at the pinned tolerance.
        guess = [3.0 * p for p in truth]
        out = unfold.sandii(pack["response"], rates, guess, tolerance=1e-9, max_iterations=50_000)
        assert out["rate_factors"] == pytest.approx([1.0] * len(rates), rel=1e-4)
