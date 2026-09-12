"""Python-side tests for MCPL interchange + MCTAL tally bodies (synthetic only).

All particle records are hand-built axis vectors with closed-form packing
math; no upstream files are read. Gzip transport is covered semantically
(compressed bytes are encoder-dependent and never asserted).
"""

import gzip
import os
from pathlib import Path

import pytest

import nucleide.mcpl as mcpl
from nucleide.mcnp import read_mctal, read_ssw

FIXTURE = os.path.join(
    os.path.dirname(__file__),
    "..",
    "fixtures",
    "mcnp",
    "mctal",
    "synthetic_tally_bodies.mctal",
)

SSW_REF = os.path.join(
    os.path.dirname(__file__),
    "..",
    "fixtures",
    "mcpl",
    "ssw_conversion",
    "reference.w",
)

SSW_SURFS = [100, 200]
SSW_KINDS = ["neutron", "gamma"]

HEADER = {
    "srcname": "nucleide-test",
    "comments": ["synthetic probe"],
    "has_userflags": True,
    "has_polarisation": False,
    "double_prec": False,
    "universal_pdgcode": None,
    "universal_weight": None,
    "blobs": [],
}

PARTICLES = [
    {
        "ekin": 2.5,
        "polarisation": [0.0, 0.0, 0.0],
        "position": [1.0, -2.0, 0.5],
        "direction": [0.0, 0.0, 1.0],
        "time": 3.0,
        "weight": 1.0,
        "pdgcode": 2112,
        "userflags": 0,
    },
    {
        "ekin": 0.662,
        "polarisation": [0.0, 0.0, 0.0],
        "position": [0.0, 0.0, 0.0],
        "direction": [1.0, 0.0, 0.0],
        "time": 0.0,
        "weight": 0.5,
        "pdgcode": 22,
        "userflags": 7,
    },
]


def test_mcpl_round_trip_single(tmp_path: Path) -> None:
    path = str(tmp_path / "probe.mcpl")
    mcpl.write_mcpl(path, HEADER, PARTICLES)
    back = mcpl.read_mcpl(path)
    assert back.version == 3
    assert back.nparticles == 2
    assert back.srcname == "nucleide-test"
    assert back.comments == ["synthetic probe"]
    assert back.has_userflags is True
    got = back.particles()
    assert len(got) == 2
    assert got[0]["ekin"] == pytest.approx(2.5, abs=1e-6)
    assert got[0]["direction"] == pytest.approx([0.0, 0.0, 1.0], abs=1e-6)
    assert got[0]["pdgcode"] == 2112
    assert got[1]["pdgcode"] == 22
    assert got[1]["userflags"] == 7
    assert got[1]["weight"] == pytest.approx(0.5, abs=1e-6)


def test_mcpl_round_trip_double_universal(tmp_path: Path) -> None:
    header = {
        **HEADER,
        "double_prec": True,
        "has_polarisation": True,
        "universal_pdgcode": 2112,
        "universal_weight": 1.5,
    }
    path = str(tmp_path / "probe.mcpl")
    mcpl.write_mcpl(path, header, PARTICLES)
    back = mcpl.read_mcpl(path)
    assert back.double_prec is True
    assert back.universal_pdgcode == 2112
    assert back.universal_weight == pytest.approx(1.5)
    for p in back.particles():
        assert p["pdgcode"] == 2112
        assert p["weight"] == pytest.approx(1.5)


def test_mcpl_gzip_transparent(tmp_path: Path) -> None:
    path = str(tmp_path / "probe.mcpl.gz")
    mcpl.write_mcpl(path, HEADER, PARTICLES)
    # Raw bytes carry the gzip magic (encoder settings differ across
    # writers, so only the magic is asserted, never the full byte stream).
    with open(path, "rb") as fh:
        assert fh.read(2) == b"\x1f\x8b"
    # gzip.open reads through the layer transparently (decompressed magic).
    with gzip.open(path, "rb") as fh:
        assert fh.read(4) == b"MCPL"
    back = mcpl.read_mcpl(path)
    assert back.nparticles == 2
    assert back.particles()[0]["ekin"] == pytest.approx(2.5, abs=1e-6)


def test_mcpl_rejects_bad_inputs(tmp_path: Path) -> None:
    bad = [dict(PARTICLES[0], direction=[1.0, 1.0, 1.0])]
    with pytest.raises(ValueError):
        mcpl.write_mcpl(str(tmp_path / "x.mcpl"), HEADER, bad)
    bad_e = [dict(PARTICLES[0], ekin=-1.0)]
    with pytest.raises(ValueError):
        mcpl.write_mcpl(str(tmp_path / "y.mcpl"), HEADER, bad_e)


def test_mctal_bodies_fixture_closed_form() -> None:
    m = read_mctal(FIXTURE)
    assert m.tally_nums == [4, 14]
    bodies = {t["number"]: t for t in m.tallies}
    assert set(bodies) == {4, 14}
    t4 = bodies[4]
    assert t4["f"] == {"count": 2, "values": pytest.approx([1.0, 2.0])}
    assert t4["e"] == {"count": 2, "values": pytest.approx([0.5, 2.0])}
    assert [v for v, _ in t4["vals"]] == pytest.approx([11.0, 12.0, 21.0, 22.0])
    assert t4["total"] == pytest.approx(66.0)
    assert m.tally_vals_array(4).shape == (4, 2)
    t14 = bodies[14]
    assert [v for v, _ in t14["vals"]] == pytest.approx([7.0])
    with pytest.raises(ValueError):
        m.tally_vals_array(999)


def test_ssw2mcpl_fixture_closed_form(tmp_path: Path) -> None:
    out = str(tmp_path / "conv.mcpl")
    assert mcpl.ssw2mcpl(SSW_REF, out, SSW_SURFS, SSW_KINDS) == 2
    back = mcpl.read_mcpl(out)
    assert back.nparticles == 2
    assert back.srcname == "ssw2mcpl"
    assert back.has_userflags is True
    got = back.particles()
    # Track 1 (neutron): energy verbatim, 3.0e5 shakes -> 3.0 ms, surf 100.
    assert got[0]["ekin"] == pytest.approx(2.5)
    assert got[0]["time"] == pytest.approx(3.0)
    assert got[0]["pdgcode"] == 2112
    assert got[0]["userflags"] == 100
    assert got[0]["position"] == pytest.approx([1.0, -2.0, 0.5])
    assert got[0]["direction"] == pytest.approx([0.0, 0.0, 1.0])
    assert got[0]["weight"] == pytest.approx(1.0)
    # Track 2 (gamma): 0.662 MeV through single precision, surf 200.
    assert got[1]["ekin"] == pytest.approx(0.662, abs=1e-6)
    assert got[1]["time"] == pytest.approx(0.0)
    assert got[1]["pdgcode"] == 22
    assert got[1]["userflags"] == 200
    assert got[1]["direction"] == pytest.approx([1.0, 0.0, 0.0], abs=1e-6)


def test_ssw2mcpl_options(tmp_path: Path) -> None:
    out = str(tmp_path / "conv.mcpl")
    mcpl.ssw2mcpl(SSW_REF, out, SSW_SURFS, SSW_KINDS, {"surf_to_userflags": False})
    back = mcpl.read_mcpl(out)
    assert back.has_userflags is False
    assert [p["userflags"] for p in back.particles()] == [0, 0]

    mcpl.ssw2mcpl(
        SSW_REF,
        out,
        SSW_SURFS,
        SSW_KINDS,
        {
            "double_prec": True,
            "srcname": "probe",
            "comments": ["synthetic"],
            "deck_blob": ("ssw_deck", b"c synthetic deck"),
        },
    )
    back = mcpl.read_mcpl(out)
    assert back.double_prec is True
    assert back.srcname == "probe"
    assert back.comments == ["synthetic"]
    assert back.blobs == [("ssw_deck", b"c synthetic deck")]

    gz = str(tmp_path / "conv.mcpl.gz")
    mcpl.ssw2mcpl(SSW_REF, gz, SSW_SURFS, SSW_KINDS, {"gzip": True})
    with open(gz, "rb") as fh:
        assert fh.read(2) == b"\x1f\x8b"
    assert mcpl.read_mcpl(gz).nparticles == 2


def test_ssw2mcpl_rejects_bad_inputs(tmp_path: Path) -> None:
    out = str(tmp_path / "conv.mcpl")
    with pytest.raises(ValueError):
        mcpl.ssw2mcpl(SSW_REF, out, [100], SSW_KINDS)
    with pytest.raises(ValueError):
        mcpl.ssw2mcpl(SSW_REF, out, SSW_SURFS, ["neutron", "proton"])
    with pytest.raises(ValueError):
        mcpl.ssw2mcpl(SSW_REF, out, SSW_SURFS, SSW_KINDS, {"double_prec": "yes"})
    with pytest.raises(ValueError):
        mcpl.ssw2mcpl(SSW_REF, out, SSW_SURFS, SSW_KINDS, ["not-a-dict"])  # type: ignore[arg-type]


def test_mcpl2ssw_round_trip_closed_form(tmp_path: Path) -> None:
    probe = str(tmp_path / "probe.mcpl")
    assert mcpl.ssw2mcpl(SSW_REF, probe, SSW_SURFS, SSW_KINDS) == 2
    out = str(tmp_path / "back.w")
    assert mcpl.mcpl2ssw(probe, SSW_REF, out) == 2
    tracks = read_ssw(out).tracks()
    assert len(tracks) == 2
    # Neutron: 3.0 ms -> 3.0e5 shakes, geometry verbatim, no cosine forced.
    assert tracks[0]["erg"] == pytest.approx(2.5)
    assert tracks[0]["tme"] == pytest.approx(3.0e5, rel=1e-6)
    assert tracks[0]["wgt"] == pytest.approx(1.0)
    assert (tracks[0]["x"], tracks[0]["y"], tracks[0]["z"]) == pytest.approx((1.0, -2.0, 0.5))
    assert (tracks[0]["u"], tracks[0]["v"], tracks[0]["cs"]) == pytest.approx((0.0, 0.0, 1.0))
    # Gamma.
    assert tracks[1]["erg"] == pytest.approx(0.662, abs=1e-6)
    assert (tracks[1]["u"], tracks[1]["v"], tracks[1]["cs"]) == pytest.approx(
        (1.0, 0.0, 0.0), abs=1e-6
    )
    # Reference header passes through (counts patched to the converted
    # tally); re-parse is stable.
    ref = read_ssw(SSW_REF)
    again = read_ssw(out)
    assert (again.nrss, again.np1, again.orignp1) == (2, 2, -2)
    assert again.kod == ref.kod
    assert again.ver == ref.ver


def test_mcpl2ssw_surface_override(tmp_path: Path) -> None:
    probe = str(tmp_path / "probe.mcpl")
    mcpl.ssw2mcpl(SSW_REF, probe, SSW_SURFS, SSW_KINDS)
    out = str(tmp_path / "back.w")
    # Override wins over userflags; the write still succeeds (surface ids
    # live outside the track records in this file layout).
    assert mcpl.mcpl2ssw(probe, SSW_REF, out, surface=7) == 2
    assert len(read_ssw(out).tracks()) == 2
    with pytest.raises(ValueError):
        mcpl.mcpl2ssw(probe, SSW_REF, out, surface=1_000_000)


def test_mcpl2ssw_rejects_unsupported_pdg(tmp_path: Path) -> None:
    probe = str(tmp_path / "proton.mcpl")
    mcpl.write_mcpl(probe, HEADER, [dict(PARTICLES[0], pdgcode=2212)])
    with pytest.raises(ValueError):
        mcpl.mcpl2ssw(probe, SSW_REF, str(tmp_path / "back.w"))
