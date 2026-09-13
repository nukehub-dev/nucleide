"""Python-side tests for the spectroscopy core (E1-E9 gates + format quirks)."""

import json
from pathlib import Path

import pytest

import nucleide.spectroscopy as sp

FIX = Path(__file__).parent.parent / "fixtures" / "spectroscopy"

COUNTS = [2.0, 5.0, 1.0, 6.0, 3.0, 8.0, 4.0]
CHANNELS = [float(c) for c in range(len(COUNTS))]

EFF_COEFF = [
    -2.81861504261204,
    -0.727352820018942,
    -0.0395798886481904,
    -0.0592305254664096,
    0.023772637347443,
    0.0325306475072671,
]

ATOMIC = {
    "k_shell_fluor": 0.9,
    "l_shell_fluor": 0.4,
    "prob": 0.8,
    "kb_to_ka": 0.2,
    "ka2_to_ka1": 0.5,
    "ka1_en_kev": 10.0,
    "ka2_en_kev": 20.0,
    "kb_en_kev": 30.0,
    "l_en_kev": 40.0,
}


def test_e1_rect_smooth() -> None:
    # Hand-computed m=5, ext=2: [2, 5, 17/5, 23/5, 22/5, 8, 4].
    assert sp.rect_smooth(COUNTS, 5) == pytest.approx([2.0, 5.0, 3.4, 4.6, 4.4, 8.0, 4.0])
    # Hand-computed m=3 on [1,3,2,6,4]: [1, 2, 11/3, 4, 4].
    assert sp.rect_smooth([1.0, 3.0, 2.0, 6.0, 4.0], 3) == pytest.approx(
        [1.0, 2.0, 11.0 / 3.0, 4.0, 4.0]
    )
    with pytest.raises(ValueError):
        sp.rect_smooth(COUNTS, 2)
    with pytest.raises(ValueError):
        sp.rect_smooth(COUNTS, 4)


def test_e2_five_point_smooth() -> None:
    # Hand-computed: [2, 5, 30/9, 39/9, 42/9, 8, 4].
    assert sp.five_point_smooth(COUNTS) == pytest.approx(
        [2.0, 5.0, 30.0 / 9.0, 39.0 / 9.0, 42.0 / 9.0, 8.0, 4.0]
    )


def test_e3_e4_e5_counts() -> None:
    # Hand-computed at c1=2, c2=5, m=1: low=7, high=12, bg=76/6.
    assert sp.calc_bg(COUNTS, CHANNELS, 2, 5, 1) == pytest.approx(76.0 / 6.0)
    # Half-open: counts[2:5] = 1+6+3.
    assert sp.gross_count(COUNTS, CHANNELS, 2, 5) == pytest.approx(10.0)
    assert sp.net_counts(COUNTS, CHANNELS, 2, 5, 1) == pytest.approx(-16.0 / 6.0)
    with pytest.raises(ValueError):
        sp.calc_bg(COUNTS, CHANNELS, 5, 2, 1)
    with pytest.raises(ValueError):
        sp.calc_bg(COUNTS, CHANNELS, -1, 2, 1)
    with pytest.raises(ValueError):
        sp.calc_bg(COUNTS, CHANNELS, 2, 7, 1)
    with pytest.raises(ValueError):
        sp.calc_bg(COUNTS, CHANNELS, 2, 5, 2)
    with pytest.raises(ValueError):
        sp.gross_count(COUNTS, CHANNELS, 5, 2)


def test_bg_small_c1_matches_python_slices() -> None:
    counts = [10.0, 10.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 20.0, 20.0, 20.0]
    channels = [float(c) for c in range(len(counts))]
    # counts[-2:0] and counts[-1:1] are empty upstream, so low = 0.
    assert sp.calc_bg(counts, channels, 0, 6, 1) == pytest.approx(77.0 / 6.0)
    assert sp.calc_bg(counts, channels, 1, 6, 1) == pytest.approx(11.0)


def test_e6_energy_bins() -> None:
    # Hand-computed fit [1.5, 2.0, 0.5] on [0, 1, 2].
    assert sp.energy_bins([0.0, 1.0, 2.0], [1.5, 2.0, 0.5]) == pytest.approx([1.5, 4.0, 7.5])


def test_e7_efficiency_golden() -> None:
    # Oracle golden: at E=1 MeV, ln E=0, eff=exp(a0)=0.059688551591347033.
    assert sp.detector_efficiency(1.0, EFF_COEFF, 1) == pytest.approx(
        0.059688551591347033, rel=1e-12
    )
    # Hand-computed fit2: coeff [0.5, -1.0] at E=2 gives exp(0)=1.
    assert sp.detector_efficiency(2.0, [0.5, -1.0], 2) == pytest.approx(1.0)
    with pytest.raises(ValueError):
        sp.detector_efficiency(1.0, EFF_COEFF, 10)


def test_e8_xray_lines() -> None:
    # Hand-computed with k=2, l=3: [(10,1), (20,0.5), (30,0.3), (40,1.84)].
    lines = sp.xray_lines(ATOMIC, k_conv=2.0, l_conv=3.0)
    assert [e for e, _ in lines] == pytest.approx([10.0, 20.0, 30.0, 40.0])
    assert [i for _, i in lines] == pytest.approx([1.0, 0.5, 0.3, 1.84])
    l_only = sp.xray_lines(ATOMIC, l_conv=3.0)
    assert [i for _, i in l_only] == pytest.approx([0.0, 0.0, 0.0, 1.2])
    # NaN plays the upstream absent-conversion sentinel.
    assert sp.xray_lines(ATOMIC, k_conv=float("nan"), l_conv=3.0) == l_only


def test_spe_readers_cross_format_counts_match() -> None:
    dollar = sp.read_dollar_spe(str(FIX / "dollar_min.spe"))
    plain = sp.read_spe(str(FIX / "plain_min.spe"))
    assert dollar["counts"] == plain["counts"]
    assert dollar["spec_name"] == "SYNTHETIC-DOLLAR"
    assert dollar["det_id"] == "7"
    assert dollar["det_descp"] == "SYNDET"
    assert dollar["real_time"] == pytest.approx(100.0)
    assert dollar["live_time"] == pytest.approx(90.0)
    assert dollar["start_date"] == "02/03/2026"
    assert dollar["start_time"] == "09:15:00"
    assert dollar["num_channels"] == 8
    assert dollar["calib_e_fit"] == pytest.approx([1.5, 2.0, 0.5])
    assert plain["spec_name"] == "SYNTHETIC-PLAIN"
    assert plain["real_time"] == pytest.approx(100.5)
    assert plain["live_time"] == pytest.approx(90.25)


def test_spe_magic_rejected_by_other_reader() -> None:
    dollar_text = (FIX / "dollar_min.spe").read_text()
    plain_text = (FIX / "plain_min.spe").read_text()
    with pytest.raises(ValueError):
        sp.parse_spe(dollar_text)
    with pytest.raises(ValueError):
        sp.parse_dollar_spe(plain_text)


def test_e9_sdef_fixture_goldens() -> None:
    # Byte-exact card text for the hand-built synthetic line lists in
    # sdef_oracle.json (provenance recorded in the fixture).
    fix = json.loads((FIX / "sdef_oracle.json").read_text())
    for case in ("single_isotropic", "single_beam", "multi_distribution", "wrapped"):
        c = fix[case]
        kwargs = dict(c["source"].items())
        bins, card = sp.sdef_decay_source(
            [(float(e), float(i)) for e, i in c["lines"]],
            version=c["version"],
            **kwargs,
        )
        assert [e for e, _ in bins] == pytest.approx([e for e, _ in c["expected_bins"]])
        assert [p for _, p in bins] == pytest.approx(
            [p for _, p in c["expected_bins"]], rel=1e-12, abs=1e-15
        )
        assert card == c["expected_card"]
        # Every emitted card line stays inside MCNP's 80-column fixed width.
        assert all(len(line) <= 80 for line in card.splitlines())


def test_e9_normalization_merges_duplicates_and_sorts() -> None:
    bins, card = sp.sdef_decay_source([(1.0, 2.0), (0.5, 1.0), (1.0, 1.0)])
    assert bins == pytest.approx([(0.5, 0.25), (1.0, 0.75)])
    assert "\nSI1 L 0.5 1\nSP1 D 0.25 0.75" in card


def test_e9_single_line_isotropic_and_beam() -> None:
    _, iso = sp.sdef_decay_source([(0.662, 2.0)])
    assert iso == "SDEF POS=0 0 0\n     ERG=0.662\n     WGT=1\n     PAR=n"
    _, beam = sp.sdef_decay_source([(0.662, 2.0)], z=3.0, w=1.0, weight=0.5, particle="Photon")
    assert beam == (
        "SDEF POS=0 0 3\n     VEC=0 0 1 DIR=1\n     ERG=0.662\n     WGT=0.5\n     PAR=p"
    )


def test_e9_particle_version_dialect() -> None:
    # Proton gains an MCNP designator only in version 6.
    with pytest.raises(ValueError):
        sp.sdef_decay_source([(0.662, 1.0)], particle="Proton", version=5)
    _, card = sp.sdef_decay_source([(0.662, 1.0)], particle="Proton", version=6)
    assert card.endswith("\n     PAR=h")
    with pytest.raises(ValueError):
        sp.sdef_decay_source([(0.662, 1.0)], particle="Waka waka")
    with pytest.raises(ValueError):
        sp.sdef_decay_source([(0.662, 1.0)], version=4)


def test_e9_rejected_inputs() -> None:
    for lines in [
        [],
        [(0.662, 0.0)],
        [(0.662, -1.0)],
        [(-0.5, 1.0)],
        [(float("nan"), 1.0)],
        [(0.662, float("inf"))],
    ]:
        with pytest.raises(ValueError):
            sp.sdef_decay_source(lines)
    with pytest.raises(ValueError):
        sp.sdef_decay_source([(0.662, 1.0)], weight=0.0)
    with pytest.raises(ValueError):
        sp.sdef_decay_source([(0.662, 1.0)], z=float("nan"))


def test_e9_zero_intensity_line_drops_to_single_line_form() -> None:
    bins, card = sp.sdef_decay_source([(0.662, 2.0), (1.17, 0.0)])
    assert bins == pytest.approx([(0.662, 1.0)])
    assert card == "SDEF POS=0 0 0\n     ERG=0.662\n     WGT=1\n     PAR=n"


def test_e9_lines_tsv_fixture_feeds_normalization() -> None:
    # decay_lines_sample.tsv rows (0.662x2 + 1.0, 1.170x1.0) merge to
    # (0.662 → 3/4, 1.17 → 1/4) through the same E9 path as list input.
    text = (FIX / "decay_lines_sample.tsv").read_text()
    assert sp.parse_lines_tsv(text) == pytest.approx([(0.662, 2.0), (1.17, 1.0), (0.662, 1.0)])
    assert sp.read_decay_lines(str(FIX / "decay_lines_sample.tsv")) == pytest.approx(
        [(0.662, 2.0), (1.17, 1.0), (0.662, 1.0)]
    )
    bins, _ = sp.sdef_decay_source(sp.parse_lines_tsv(text))
    assert bins == pytest.approx([(0.662, 0.75), (1.17, 0.25)])


def test_e9_lines_tsv_rejects_malformed_rows() -> None:
    with pytest.raises(ValueError):
        sp.parse_lines_tsv("")
    with pytest.raises(ValueError):
        sp.parse_lines_tsv("# only a comment\n")
    with pytest.raises(ValueError):
        sp.parse_lines_tsv("0.662\n")
    with pytest.raises(ValueError):
        sp.parse_lines_tsv("0.662 1.0 3.0\n")
    with pytest.raises(ValueError):
        sp.parse_lines_tsv("0.662 lots\n")


def test_e7_fit_fixture_closed_form() -> None:
    # Synthetic closed-form points (fixtures/spectroscopy/efficiency_fit.json):
    # every efficiency sits exactly on its E7 curve, so the fit recovers the
    # recorded coefficients within 1e-9 and re-evaluates within 1e-9.
    fix = json.loads((FIX / "efficiency_fit.json").read_text())
    for case in ("fit1_degree2", "fit2_degree1"):
        c = fix[case]
        got = sp.fit_efficiency(c["energies"], c["effs"], c["weights"], c["order"], c["fit"])
        assert got == pytest.approx(c["expected_coeff"], rel=1e-9, abs=1e-12)
        for e, v in zip(c["energies"], c["effs"], strict=True):
            assert sp.detector_efficiency(e, got, c["fit"]) == pytest.approx(v, rel=1e-9)


def test_e7_fit_zero_weight_skips_row() -> None:
    # eff = E pins the fit-1 line a0 + a1 ln E at [0, 1]; (3, 300) is skipped.
    got = sp.fit_efficiency([1.0, 2.0, 3.0], [1.0, 2.0, 300.0], [1.0, 1.0, 0.0], 1, 1)
    assert got == pytest.approx([0.0, 1.0], abs=1e-9)


def test_e7_fit_rejects() -> None:
    with pytest.raises(ValueError):
        sp.fit_efficiency([], [], [], 1, 1)
    with pytest.raises(ValueError):
        sp.fit_efficiency([0.5, 1.0], [0.1], [1.0, 1.0], 1, 1)
    with pytest.raises(ValueError):
        sp.fit_efficiency([0.5, 1.0], [0.1, 0.2], [1.0, 1.0], 1, 7)
    with pytest.raises(ValueError):
        sp.fit_efficiency([0.5, 1.0], [0.1, 0.2], [1.0, 1.0], 5, 1)
    with pytest.raises(ValueError):
        sp.fit_efficiency([0.0, 1.0], [0.1, 0.2], [1.0, 1.0], 1, 1)
    with pytest.raises(ValueError):
        sp.fit_efficiency([0.5, 1.0], [0.0, 0.2], [1.0, 1.0], 1, 1)
    with pytest.raises(ValueError):
        sp.fit_efficiency([0.5, 1.0], [0.1, 0.2], [1.0, -1.0], 1, 1)
    with pytest.raises(ValueError):
        sp.fit_efficiency([0.5, 1.0], [0.1, 0.2], [0.0, 0.0], 0, 1)
