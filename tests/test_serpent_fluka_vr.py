"""Golden tests for Serpent/FLUKA readers, MAGIC, sampling, writers."""

from pathlib import Path

import pytest

import nucleide

FIX = Path(__file__).parent.parent / "fixtures"


class TestSerpent:
    def test_res(self) -> None:
        r = nucleide.serpent.read_serpent(str(FIX / "serpent" / "serp2_res.m"), "res")
        assert "ABS_KEFF" in r and "CONVERSION_RATIO" in r
        # Matrix entry [cycle][mean, stdev]
        keff = r["ABS_KEFF"]
        assert keff[0][0][0] == pytest.approx(1.01503, abs=1e-4)
        assert keff[0][0][1] == pytest.approx(0.00324, abs=1e-5)

    def test_dep(self) -> None:
        d = nucleide.serpent.read_serpent(str(FIX / "serpent" / "sample1_dep.m"), "dep")
        assert "ZAI" in d and "DAYS" in d
        assert len(d["ZAI"]) > 0

    def test_det(self) -> None:
        det = nucleide.serpent.read_serpent(str(FIX / "serpent" / "serp2_det.m"), "det")
        assert any("DET" in k for k in det)


class TestFluka:
    def test_usrbin_single(self) -> None:
        tallies = nucleide.fluka.read_usrbin(str(FIX / "fluka" / "fluka_usrbin_single.lis"))
        assert len(tallies) == 1
        t = tallies[0]
        nx, ny, nz = t.dims()
        assert nx * ny * nz == len(t.data)
        assert len(t.error) == len(t.data)

    def test_usrbin_multiple(self) -> None:
        tallies = nucleide.fluka.read_usrbin(str(FIX / "fluka" / "fluka_usrbin_multiple.lis"))
        assert len(tallies) > 1

    def test_usrbin_degenerate(self) -> None:
        tallies = nucleide.fluka.read_usrbin(str(FIX / "fluka" / "fluka_usrbin_degenerate.lis"))
        assert len(tallies) >= 1


class TestMagicAndSampling:
    def setup_method(self) -> None:
        self.meshtal = nucleide.mcnp.read_meshtal(
            str(FIX / "mcnp" / "meshtal" / "mcnp_meshtal_single_meshtal.txt")
        )
        self.tally = self.meshtal.tallies[4]

    def test_magic_total(self) -> None:
        out = nucleide.vr.magic(self.tally)
        nve = self.tally.num_ves()
        assert len(out.lower_bounds_ww) == nve
        assert out.groups_per_ve == 1
        # Null entries where relative error exceeded tolerance are 0.0
        vals = out.lower_bounds_ww
        assert any(v > 0 for v in vals)

    def test_magic_per_group(self) -> None:
        out = nucleide.vr.magic(self.tally, per_group=True)
        assert len(out.lower_bounds_ww) == self.tally.num_ves() * 3

    def test_alias_table_round_trip(self) -> None:
        pdf = [0.9, 0.05, 0.03, 0.02]
        table = nucleide.vr.AliasTable(pdf)
        assert len(table) == 4
        # Deterministic draws land in range
        for i in range(100):
            idx = table.sample(i / 128.0, (i * 7 % 97) / 97.0)
            assert 0 <= idx < 4

    def test_sampler_modes(self) -> None:
        analog = nucleide.vr.MeshSourceSampler(self.tally, "analog")
        uniform = nucleide.vr.MeshSourceSampler(self.tally, "uniform")
        s = analog.sample(0.3, 0.7)
        assert set(s) == {"index", "i", "j", "k", "weight"}
        assert s["weight"] == 1.0  # analog birth weight is unity
        u = uniform.sample(0.3, 0.7)
        assert u["weight"] > 0

    def test_user_mode(self) -> None:
        nve = self.tally.num_ves()
        sampler = nucleide.vr.MeshSourceSampler(self.tally, "user", user_pdf=[1.0] * nve)
        s = sampler.sample(0.5, 0.5)
        assert 0 <= s["index"] < nve


class TestWriters:
    def test_ssw_round_trip_bytes(self, tmp_path: Path) -> None:
        src = FIX / "mcnp" / "ssw" / "mcnp_surfsrc_onetrack.w"
        ssw = nucleide.mcnp.read_ssw(str(src))
        tracks = ssw.tracks()
        out = tmp_path / "roundtrip.w"
        nucleide.mcnp.write_ssw(ssw, str(out), tracks)
        assert out.read_bytes() == src.read_bytes()

    def test_mesh_to_geom_structure(self) -> None:
        xb = [-200.0, -66.67, 66.67, 200.0]
        yb = [-200.0, 200.0]
        zb = [-200.0, 200.0]
        nve = 3 * 1 * 1
        mats = [("water", 1.0)] * nve
        deck = nucleide.mcnp.mesh_to_geom(xb, yb, zb, mats, "test deck")
        assert deck.startswith("test deck")
        # surfaces as px planes + graveyard shell, per the mesh_to_geom oracle
        assert "px -200.0" in deck
        assert "0 -1:4:-5:6:-7:8" in deck
