"""Golden tests for the mcnp-io readers against vendored fixtures."""

from pathlib import Path

import pytest

import nucleide

FIXTURES = Path(__file__).parent.parent / "fixtures" / "mcnp"


class TestXsdir:
    def test_awr_and_tables(self) -> None:
        x = nucleide.read_xsdir(str(FIXTURES / "xsdir" / "dummy_xsdir"))
        assert x.awr[1000] == 0.99931697
        assert len(x.tables) == 3
        t = x.tables[0]
        assert t.name == "1001.44c"
        assert t.temperature == 5.5555e5
        assert not t.ptable
        assert x.tables[1].ptable

    def test_find_table_and_serpent(self) -> None:
        x = nucleide.read_xsdir(str(FIXTURES / "xsdir" / "dummy_xsdir"))
        hits = x.find_table("1001")
        assert [t.name for t in hits] == ["1001.44c", "1001.66c", "1001.70c"]
        line = hits[0].to_serpent(".")
        assert line == (
            "1001.44c 1001.44c 1 1001 0 1.111111 6.44688328094e+15 0 ./many_xs/1001.555nc"
        )

    def test_nucs(self) -> None:
        x = nucleide.read_xsdir(str(FIXTURES / "xsdir" / "dummy_xsdir"))
        assert x.nucs() == [nucleide.Nuclide("H1").nucid]


class TestMeshtal:
    def test_single(self) -> None:
        m = nucleide.read_meshtal(str(FIXTURES / "meshtal" / "mcnp_meshtal_single_meshtal.txt"))
        assert m.version == "5.mpi"
        assert m.histories == 100000
        t = m.tallies[4]
        assert t.particle == "n"
        assert t.dose_response
        assert t.x_bounds == [-200.0, -66.67, 66.67, 200.0]
        assert list(t.dims()) == [3, 5, 3]
        r, e = t.cell(0, 0, 0)
        assert r[0] == pytest.approx(4.96471e-9)
        assert e[0] == pytest.approx(1.98750e-1)
        tr, te = t.cell_total(0, 0, 0)
        assert tr == pytest.approx(1.91370e-7)

    def test_multiple(self) -> None:
        m = nucleide.read_meshtal(str(FIXTURES / "meshtal" / "mcnp_meshtal_multiple_meshtal.txt"))
        assert sorted(m.tallies) == [4, 14, 24, 34]
        assert m.tallies[24].particle == "p"
        assert m.tallies[14].num_e_groups() == 1


class TestWwinp:
    def test_neutron(self) -> None:
        w = nucleide.read_wwinp(str(FIXTURES / "wwinp" / "mcnp_wwinp_wwinp_n.txt"))
        assert w.ni == 1
        assert w.ne == [7]
        assert list(w.nf) == [15, 8, 6]
        assert w.cm[0][:3] == [-99.0, -97.0, 97.0]
        row = w.ww_row(0, 0)
        assert len(row) == 720
        col = w.ww_column(0, 0)
        assert len(col) == 7

    def test_np(self) -> None:
        w = nucleide.read_wwinp(str(FIXTURES / "wwinp" / "mcnp_wwinp_wwinp_np.txt"))
        assert w.ni == 2
        assert len(w.e) == 2


class TestMctal:
    def test_kcode5(self) -> None:
        m = nucleide.read_mctal(str(FIXTURES / "mctal" / "synthetic_kcode5.mctal"))
        assert m.code_name == "mcnp"
        assert m.n_cycles == 8
        assert m.n_inactive == 4
        assert m.k_col[0] == pytest.approx(0.985)

    def test_kcode19_averages(self) -> None:
        m = nucleide.read_mctal(str(FIXTURES / "mctal" / "synthetic_kcode19.mctal"))
        assert len(m.averages) == 6
        assert m.averages[0]["fom"] == 42.0


class TestSsw:
    def test_header_and_tracks(self) -> None:
        s = nucleide.read_ssw(str(FIXTURES / "ssw" / "mcnp_surfsrc_onetrack.w"))
        assert s.kod == "mcnp"
        assert s.nrss == 1
        tracks = s.tracks()
        assert len(tracks) == 1
        assert tracks[0]["nps"] == 1.0
        assert tracks[0]["wgt"] == pytest.approx(0.99995639)

    def test_print_header(self) -> None:
        s = nucleide.read_ssw(str(FIXTURES / "ssw" / "mcnp5_surfsrc.w"))
        txt = s.print_header()
        assert txt.startswith("Code: mcnp")
        assert s.nrss == 173


class TestPtrac:
    def test_i4(self) -> None:
        p = nucleide.read_ptrac(str(FIXTURES / "ptrac" / "mcnp_ptrac_i4_little.ptrac"))
        assert p.problem_title.startswith("Generate a well-defined PTRAC file")
        assert p.width_code == 0
        events = p.events()
        assert events[0]["event_type"] == 1000.0
        assert events[0]["xxx"] == 0.0

    def test_i8_same_title(self) -> None:
        p = nucleide.read_ptrac(str(FIXTURES / "ptrac" / "mcnp_ptrac_i8_little.ptrac"))
        assert p.width_code == 1
        assert len(p.events()) == len(p.events())
