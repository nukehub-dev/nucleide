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
        assert keff[0][0] == pytest.approx(1.01503, abs=1e-4)
        assert keff[0][1] == pytest.approx(0.00324, abs=1e-5)

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


class TestWeightWindowEmission:
    """MAGIC -> OpenMC settings.xml / Serpent WWINP emission round trips."""

    def setup_method(self) -> None:
        self.meshtal = nucleide.mcnp.read_meshtal(
            str(FIX / "mcnp" / "meshtal" / "mcnp_meshtal_single_meshtal.txt")
        )
        self.tally = self.meshtal.tallies[4]

    @staticmethod
    def xfastest(lower: list[float], dims: tuple[int, int, int], groups: int) -> list[float]:
        """Reorder ve-major z-fastest windows to the x-fastest group-outermost
        flat layout both emitters write."""
        nx, ny, nz = dims
        flat: list[float] = []
        for g in range(groups):
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        ve = (i * ny + j) * nz + k
                        flat.append(lower[ve * groups + g])
        return flat

    def test_openmc_fragment_recovers_magic_windows(self) -> None:
        import xml.etree.ElementTree as ET

        out = nucleide.vr.magic(self.tally, per_group=True)
        result = nucleide.vr.emit_openmc_weight_windows(self.tally, out, mesh_id=7, window_id=9)
        root = ET.fromstring(f"<settings>{result['xml']}</settings>")
        mesh = root.find("mesh")
        assert mesh is not None and mesh.get("id") == "7"
        dims = self.tally.dims()
        assert [int(v) for v in mesh.findtext("dimension", "").split()] == list(dims)
        ww = root.find("weight_windows")
        assert ww is not None and ww.get("id") == "9"
        assert ww.findtext("mesh") == "7"
        assert ww.findtext("particle_type") == "neutron"
        n_groups = out.groups_per_ve
        e_bounds = [float(v) for v in ww.findtext("energy_bounds", "").split()]
        assert e_bounds == [0.0] + [e * 1.0e6 for e in out.e_upper_bounds]
        nft = dims[0] * dims[1] * dims[2]
        lower = [float(v) for v in ww.findtext("lower_ww_bounds", "").split()]
        upper = [float(v) for v in ww.findtext("upper_ww_bounds", "").split()]
        assert len(lower) == len(upper) == nft * n_groups
        expected = self.xfastest(list(out.lower_bounds_ww), dims, n_groups)
        assert lower == pytest.approx(expected, abs=1e-12)
        assert upper == pytest.approx([5.0 * v for v in expected], abs=1e-12)
        assert float(ww.findtext("survival_ratio", "0")) == 3.0
        assert int(ww.findtext("max_split", "0")) == 10
        assert float(ww.findtext("weight_cutoff", "1")) == 1e-38
        assert any("upper-bounds-synthesized" in n for n in result["notes"])

    def test_serpent_wwin_recovers_magic_windows(self) -> None:
        out = nucleide.vr.magic(self.tally, per_group=True)
        result = nucleide.vr.emit_serpent_wwin(self.tally, out, name="ww1", file="m.wwd")
        assert result["card"] == 'wwin ww1 wf "m.wwd" 2'
        assert any("FMT=1" in n for n in result["notes"])

        # Minimal WWINP structural parse (header + block-3 window rows).
        toks = iter(result["text"].split())
        next(toks)  # if = 1
        next(toks)  # iv = 1
        ni = int(next(toks))
        nr = int(next(toks))
        assert (ni, nr) == (1, 10)
        ne = [int(next(toks)) for _ in range(1)]
        groups = ne[0]
        nf = [int(float(next(toks))) for _ in range(3)]
        origin = [float(next(toks)) for _ in range(3)]
        nc = [int(float(next(toks))) for _ in range(3)]
        next(toks)  # nwg
        dims = self.tally.dims()
        assert nf == list(dims)
        assert origin == [
            self.tally.x_bounds[0],
            self.tally.y_bounds[0],
            self.tally.z_bounds[0],
        ]
        for axis in range(3):
            for _ in range(3 * nc[axis] + 1):
                next(toks)
        energies = [float(next(toks)) for _ in range(groups)]
        assert energies == pytest.approx(list(out.e_upper_bounds), abs=1e-9)
        nft = nf[0] * nf[1] * nf[2]
        rows = [[float(next(toks)) for _ in range(nft)] for _ in range(groups)]
        expected = self.xfastest(list(out.lower_bounds_ww), dims, groups)
        for g in range(groups):
            assert rows[g] == pytest.approx(expected[g * nft : (g + 1) * nft], abs=1e-6)

    def test_serpent_total_mode_and_photon_header(self) -> None:
        out = nucleide.vr.magic(self.tally)
        result = nucleide.vr.emit_serpent_wwin(self.tally, out)
        toks = result["text"].split()
        assert int(toks[2]) == 1  # ni
        assert [int(v) for v in toks[4:5]] == [1]  # one energy group
        assert result["card"].endswith('wf "wwindows.wwd" 2')

    def test_emission_loud_errors(self) -> None:
        out = nucleide.vr.magic(self.tally)
        # A negative null_value nulls every cell, exercising the loud
        # negative-window rejection in both emitters.
        out_neg = nucleide.vr.magic_with(
            self.tally, selection="per_group", tolerance=-1.0, null_value=-1.0
        )
        assert any(v == -1.0 for v in out_neg.lower_bounds_ww)
        with pytest.raises(ValueError, match="negative"):
            nucleide.vr.emit_serpent_wwin(self.tally, out_neg)
        with pytest.raises(ValueError, match="negative"):
            nucleide.vr.emit_openmc_weight_windows(self.tally, out_neg)
        with pytest.raises(ValueError, match="survival_ratio"):
            nucleide.vr.emit_openmc_weight_windows(self.tally, out, survival_ratio=1.0)


class TestKdeSampler:
    def test_fit_draw_pdf(self) -> None:
        kde = nucleide.vr.KdeSampler([[1.0, 2.0], [3.0, 4.0]], [0.5, 2.0])
        assert kde.n_samples() == 2
        assert kde.bandwidths() == [0.5, 2.0]
        assert kde.draw(0.75, [1.0, -0.5]) == [3.5, 3.0]
        assert kde.pdf([3.0, 4.0]) > kde.pdf([5.0, 6.0]) > 0.0

    def test_adaptive_hand_vector(self) -> None:
        kde = nucleide.vr.KdeSampler([[0.0], [0.2], [0.4], [10.0]], "adaptive", k=1, scale=1.0)
        assert kde.n_samples() == 4
        sw = kde.sample_bandwidths()
        assert sw[0][0] == pytest.approx(0.2, abs=1e-12)
        assert sw[1][0] == pytest.approx(0.2, abs=1e-12)
        assert sw[2][0] == pytest.approx(0.2, abs=1e-12)
        assert sw[3][0] == pytest.approx(9.6, abs=1e-9)
        assert sw[0][0] != sw[3][0]
        base = kde.bandwidths()
        assert base[0] == pytest.approx(0.04079367910119179, abs=1e-12)
        assert base[3] == pytest.approx(1.9580965968572057, abs=1e-12)
        assert kde.pdf([0.0]) == pytest.approx(0.8746688805236177, abs=1e-9)
        got = kde.draw(0.5, [1.0])
        assert got[0] == 0.4 + sw[2][0] * 1.0
        assert got[0] == pytest.approx(0.6, abs=1e-12)

    def test_adaptive_rank_scale_and_defaults(self) -> None:
        rows = [[0.0], [0.2], [0.4], [10.0]]
        kde = nucleide.vr.KdeSampler(rows, "adaptive", k=2)
        sw = kde.sample_bandwidths()
        assert sw[0][0] == pytest.approx(0.4, abs=1e-12)
        assert sw[3][0] == pytest.approx(9.8, abs=1e-9)
        scaled = nucleide.vr.KdeSampler(rows, "adaptive", scale=2.0)
        assert scaled.sample_bandwidths()[0][0] == pytest.approx(0.4, abs=1e-12)
        # Defaults pin k = 1, scale = 1.0.
        defaulted = nucleide.vr.KdeSampler(rows, "adaptive")
        assert defaulted.sample_bandwidths()[0][0] == pytest.approx(0.2, abs=1e-12)

    def test_adaptive_2d_correlated_pairs(self) -> None:
        kde = nucleide.vr.KdeSampler(
            [[0.0, 0.0], [1.0, 1.0], [4.0, 4.0], [5.0, 5.0]], "adaptive", k=1, scale=1.0
        )
        sw = kde.sample_bandwidths()
        assert all(
            row[0] == pytest.approx(1.4142135623730951, abs=1e-12)
            and row[1] == pytest.approx(1.4142135623730951, abs=1e-12)
            for row in sw
        )
        assert kde.pdf([2.5, 2.5]) > kde.pdf([2.5, 5.0])
        got = kde.draw(0.0, [1.0, 1.0])
        assert got[0] == got[1]
        assert got[0] == pytest.approx(1.4142135623730951, abs=1e-12)

    def test_adaptive_loud_errors(self) -> None:
        rows = [[0.0], [0.2], [0.4], [10.0]]
        with pytest.raises(ValueError, match="k = 0"):
            nucleide.vr.KdeSampler(rows, "adaptive", k=0)
        with pytest.raises(ValueError):
            nucleide.vr.KdeSampler(rows, "adaptive", k=4)
        with pytest.raises(ValueError):
            nucleide.vr.KdeSampler([[1.0]], "adaptive", k=1)
        with pytest.raises(ValueError):
            nucleide.vr.KdeSampler(rows, "adaptive", scale=0.0)
        with pytest.raises(ValueError):
            nucleide.vr.KdeSampler(rows, "adaptive", scale=-1.0)
        with pytest.raises(ValueError):
            nucleide.vr.KdeSampler(rows, "adaptive", scale=float("nan"))
        with pytest.raises(ValueError):
            nucleide.vr.KdeSampler(rows, "adaptive", scale=float("inf"))
        with pytest.raises(ValueError, match="duplicate"):
            nucleide.vr.KdeSampler([[0.0], [0.0], [1.0]], "adaptive", k=1)
        with pytest.raises(ValueError, match="zero variance"):
            nucleide.vr.KdeSampler([[1.0, 0.0], [1.0, 1.0], [1.0, 2.0]], "adaptive", k=1)
        with pytest.raises(ValueError, match="k/scale only apply"):
            nucleide.vr.KdeSampler(rows, "silverman", k=1)

    def test_silverman_and_errors(self) -> None:
        kde = nucleide.vr.KdeSampler([[0.0], [1.0], [2.0], [3.0]])
        assert kde.bandwidths()[0] > 0.0
        with pytest.raises(ValueError, match="zero variance"):
            nucleide.vr.KdeSampler([[3.0]] * 8)
        with pytest.raises(ValueError):
            nucleide.vr.KdeSampler([[1.0, 2.0], [3.0]])
        with pytest.raises(ValueError):
            nucleide.vr.KdeSampler([[1.0]], "scott")
        kde = nucleide.vr.KdeSampler([[0.0], [2.0]])
        with pytest.raises(ValueError):
            kde.draw(1.0, [0.0])
        with pytest.raises(ValueError):
            kde.pdf([0.0, 0.0])


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
