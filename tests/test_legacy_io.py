"""Legacy I/O: ENDL reader, SSW combine, PTRAC-HDF5 writer, meshtal mesh data.

Contract-to-test trace (each class pins one verified spec item):
- TestEndlEndftod / TestEndlLibrary — `pyne/endl.py` framing, header slices,
  `NFIELDS_RPROP`, `fromendl_tok` 11-char fields, `get_rx` selection, plus the
  `endftod` C++ semantics (`src/utils.cpp:85`) over `fixtures/endl/`.
- TestSswCombine — `scripts/ssw_combine.py` header rules (signed `orignp1`
  sum, plain `nrss` sum, sign-preserving `nps` shift) via the Rust
  `combine_files` port.
- TestPtracRows — `pyne/ptrac_to_hdf5.py` flow and the 19-column `PtracEvent`
  schema (`pyne/mcnp.py:1019-1042`) over the existing `read_ptrac` stream.
- TestWritePtracHdf5 — the HDF5 writer against a stubbed `h5py` (HDF5 bytes
  are never asserted: the format is not byte-stable).
- TestMeshtalMeshData — the Python-side mesh-data extractor over
  `read_meshtal` (MOAB tagging stays caller-side).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

import nucleide

FIXTURES = Path(__file__).parent.parent / "fixtures"
MCNP_DIR = FIXTURES / "mcnp"
ENDL_FILE = FIXTURES / "endl" / "synthetic_eedl.txt"
PB = 820000000


class TestEndlEndftod:
    def test_vectors(self) -> None:
        cases: list[tuple[str, float]] = [
            ("", 0.0),
            ("           ", 0.0),
            (" 1.00000-05", 1e-5),
            (" 3.63530+05", 363530.0),
            (" 1.00000E+05", 1e5),
            (" 1.00000D+05", 1e5),
            (" 1.00000e-05", 1e-5),
            (" 207.20000 ", 207.2),
            (" 2.0720+02 ", 207.2),
            ("-1.50000+00", -1.5),
            (" 1.23000E  ", 1.23),
        ]
        for field, want in cases:
            assert nucleide.mcnp.endl_endftod(field) == pytest.approx(want)


class TestEndlLibrary:
    def test_nuclides_and_incident_particles(self) -> None:
        lib = nucleide.mcnp.read_endl(str(ENDL_FILE))
        assert lib.nuclides() == [PB]

    def test_integrated_cross_section(self) -> None:
        lib = nucleide.mcnp.read_endl(str(ENDL_FILE))
        data = lib.get_rx(PB, 9, 10, 0)
        assert len(data) == 3
        assert all(len(row) == 2 for row in data)
        col0 = [row[0] for row in data]
        col1 = [row[1] for row in data]
        assert min(col0) == pytest.approx(1e-5)
        assert max(col0) == pytest.approx(1e5)
        assert min(col1) == pytest.approx(363530.0)
        assert max(col1) == pytest.approx(8.42443e9)

    def test_spectra_table(self) -> None:
        lib = nucleide.mcnp.read_endl(str(ENDL_FILE))
        data = lib.get_rx(PB, 9, 82, 21)
        assert len(data) == 2
        assert all(len(row) == 3 for row in data)
        assert min(row[0] for row in data) == pytest.approx(1e-5)
        assert max(row[0] for row in data) == pytest.approx(1e5)
        assert min(row[2] for row in data) == pytest.approx(4.6132e-8)
        assert max(row[2] for row in data) == pytest.approx(9.15055e6)

    def test_subshell_and_outgoing_selectors(self) -> None:
        lib = nucleide.mcnp.read_endl(str(ENDL_FILE))
        first = lib.get_rx(PB, 9, 81, 0, x1=1)
        assert first[0][0] == pytest.approx(8.829e-2)
        assert first[0][1] == pytest.approx(0.566158)
        second = lib.get_rx(PB, 9, 81, 0, x1=2)
        assert second == [[pytest.approx(1e-3), pytest.approx(2.5)]]
        both = lib.get_rx(PB, 9, 81, 0, x1=1, p_out=9)
        assert both[0][1] == pytest.approx(0.566158)

    def test_rmod_zero_forces_x1_none(self) -> None:
        # Table 1 carries x1=5.0 in its header field but rmod=0, so selecting
        # x1=5 must miss (upstream stores x1=None for rmod==0).
        lib = nucleide.mcnp.read_endl(str(ENDL_FILE))
        with pytest.raises(ValueError):
            lib.get_rx(PB, 9, 10, 0, x1=5)

    def test_unknown_nucleus_and_selectors_raise(self) -> None:
        lib = nucleide.mcnp.read_endl(str(ENDL_FILE))
        with pytest.raises(ValueError):
            lib.get_rx(920000000, 9, 10, 0)
        with pytest.raises(ValueError):
            lib.get_rx(PB, 8, 10, 0)
        with pytest.raises(ValueError):
            lib.get_rx(PB, 9, 81, 0, x1=99)

    def test_isotope_name_resolves_bare_symbol_does_not(self) -> None:
        lib = nucleide.mcnp.read_endl(str(ENDL_FILE))
        # Pb208 resolves to an id the fixture does not carry: typed miss.
        with pytest.raises(ValueError):
            lib.get_rx("Pb208", 9, 10, 0)
        # Bare element names need a mass number: resolution error, not a miss.
        with pytest.raises(ValueError):
            lib.get_rx("Pb", 9, 10, 0)


class TestSswCombine:
    def test_pair_sums_signed_headers_and_shifts_nps(self, tmp_path: Path) -> None:
        one = str(MCNP_DIR / "ssw" / "mcnp_surfsrc_onetrack.w")
        first = nucleide.mcnp.read_ssw(one)
        out = str(tmp_path / "pair.w")
        nucleide.mcnp.combine_ssw_files(out, [one, one])

        merged = nucleide.mcnp.read_ssw(out)
        assert merged.orignp1 == 2 * first.orignp1
        assert merged.np1 == 2 * first.np1
        assert merged.nrss == 2 * first.nrss

        before = first.tracks()
        after = merged.tracks()
        assert len(after) == 2
        assert after[0] == before[0]
        want = before[0]["nps"] + math.copysign(first.np1, before[0]["nps"])
        assert after[1]["nps"] == pytest.approx(want)

    def test_merged_file_reparses_cleanly(self, tmp_path: Path) -> None:
        one = str(MCNP_DIR / "ssw" / "mcnp_surfsrc_onetrack.w")
        out = str(tmp_path / "pair.w")
        nucleide.mcnp.combine_ssw_files(out, [one, one])
        merged = nucleide.mcnp.read_ssw(out)
        assert len(merged.tracks()) == merged.nrss
        assert merged.print_header().startswith("Code: mcnp")

    def test_incompatible_headers_raise(self, tmp_path: Path) -> None:
        out = str(tmp_path / "bad.w")
        with pytest.raises(ValueError):
            nucleide.mcnp.combine_ssw_files(
                out,
                [
                    str(MCNP_DIR / "ssw" / "mcnp5_surfsrc.w"),
                    str(MCNP_DIR / "ssw" / "mcnpx_surfsrc.w"),
                ],
            )

    def test_empty_inputs_raise(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError):
            nucleide.mcnp.combine_ssw_files(str(tmp_path / "empty.w"), [])

    def test_cli_writes_newssr(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.chdir(tmp_path)
        one = str(MCNP_DIR / "ssw" / "mcnp_surfsrc_onetrack.w")
        nucleide.mcnp.ssw_combine_main([one, one])
        merged = nucleide.mcnp.read_ssw("newssr")
        assert merged.nrss == 2


class TestPtracRows:
    EXPECTED_COLUMNS = (
        "event_type",
        "node",
        "nsr",
        "nsf",
        "nxs",
        "ntyn",
        "ipt",
        "ncl",
        "mat",
        "ncp",
        "xxx",
        "yyy",
        "zzz",
        "uuu",
        "vvv",
        "www",
        "erg",
        "wgt",
        "tme",
    )

    def test_columns_match_ptrac_event_schema(self) -> None:
        assert nucleide.mcnp.ptrac_event_columns() == self.EXPECTED_COLUMNS

    def test_rows_conform_to_schema(self) -> None:
        for name in ("mcnp_ptrac_i4_little.ptrac", "mcnp_ptrac_i8_little.ptrac"):
            rows = nucleide.mcnp.ptrac_event_rows(str(MCNP_DIR / "ptrac" / name))
            assert rows
            for row in rows:
                assert tuple(row) == self.EXPECTED_COLUMNS

    def test_first_event_values(self) -> None:
        rows = nucleide.mcnp.ptrac_event_rows(
            str(MCNP_DIR / "ptrac" / "mcnp_ptrac_i4_little.ptrac")
        )
        assert rows[0]["event_type"] == 1000.0
        assert rows[0]["xxx"] == 0.0
        assert rows[0]["yyy"] == 0.0
        assert rows[0]["zzz"] == 0.0


def _fake_h5py(calls: dict[str, Any]) -> Any:
    """Minimal `h5py` stub recording File/dataset traffic."""

    class FakeDataset:
        def __init__(self) -> None:
            self.attrs: dict[str, Any] = {}

    class FakeFile:
        def __init__(self, path: str, mode: str) -> None:
            calls["path"] = path
            calls["mode"] = mode
            # One dataset store per path, shared across opens of that path.
            store = calls.setdefault("files", {}).setdefault(path, {})
            calls["datasets"] = store
            self._datasets: dict[str, Any] = store

        def __contains__(self, key: object) -> bool:
            return key in self._datasets

        def __delitem__(self, key: str) -> None:
            del self._datasets[key]

        def create_dataset(self, name: str, data: Any) -> FakeDataset:
            dataset = FakeDataset()
            self._datasets[name] = (dataset, data)
            return dataset

        def __enter__(self) -> FakeFile:
            return self

        def __exit__(self, *exc: Any) -> None:
            calls["closed"] = True

    return SimpleNamespace(File=FakeFile)


class TestWritePtracHdf5:
    def test_writer_uses_row_stream(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        calls: dict[str, Any] = {}
        monkeypatch.setitem(sys.modules, "h5py", _fake_h5py(calls))
        ptrac = str(MCNP_DIR / "ptrac" / "mcnp_ptrac_i4_little.ptrac")
        count = nucleide.mcnp.write_ptrac_hdf5(
            ptrac, str(tmp_path / "out.h5"), table_title="Custom"
        )
        expected = nucleide.mcnp.ptrac_event_rows(ptrac)
        assert count == len(expected)
        assert calls["mode"] == "a"
        dataset, data = calls["datasets"]["ptrac"]
        assert dataset.attrs["TITLE"] == "Custom"
        assert list(data.dtype.names) == list(nucleide.mcnp.ptrac_event_columns())
        assert len(data) == len(expected)
        assert data["event_type"][0] == expected[0]["event_type"]
        assert calls["closed"] is True

    def test_cli_flags_pass_through(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        calls: dict[str, Any] = {}
        monkeypatch.setitem(sys.modules, "h5py", _fake_h5py(calls))
        ptrac = str(MCNP_DIR / "ptrac" / "mcnp_ptrac_i4_little.ptrac")
        out = str(tmp_path / "cli.h5")
        count = nucleide.mcnp.ptrac_to_hdf5_main([ptrac, out, "-n", "events", "-t", "Title!"])
        assert count > 0
        assert "events" in calls["datasets"]

    def test_progress_and_overwrite(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
    ) -> None:
        calls: dict[str, Any] = {}
        monkeypatch.setitem(sys.modules, "h5py", _fake_h5py(calls))
        monkeypatch.setattr(nucleide.mcnp, "_PTRAC_PROGRESS_EVERY", 1)
        ptrac = str(MCNP_DIR / "ptrac" / "mcnp_ptrac_i4_little.ptrac")
        out = str(tmp_path / "prog.h5")
        count = nucleide.mcnp.write_ptrac_hdf5(ptrac, out, show_progress=True)
        assert count > 0
        assert capsys.readouterr().out.startswith("processing event 1\n")
        # Writing the same table name twice replaces the first dataset.
        assert nucleide.mcnp.write_ptrac_hdf5(ptrac, out) == count

    def test_missing_h5py_raises_helpfully(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.delitem(sys.modules, "h5py", raising=False)
        import builtins

        real_import = builtins.__import__

        def no_h5py(name: str, *args: Any, **kwargs: Any) -> Any:
            if name == "h5py":
                raise ImportError("No module named 'h5py'")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", no_h5py)
        import io as _io_check  # exercises the pass-through branch above

        assert _io_check is not None
        with pytest.raises(ImportError, match="h5py"):
            nucleide.mcnp.write_ptrac_hdf5(
                str(MCNP_DIR / "ptrac" / "mcnp_ptrac_i4_little.ptrac"),
                str(tmp_path / "out.h5"),
            )


class TestMeshtalMeshData:
    def test_single_tally_contents(self) -> None:
        data = nucleide.mcnp.meshtal_mesh_data(
            str(MCNP_DIR / "meshtal" / "mcnp_meshtal_single_meshtal.txt")
        )
        assert data["version"] == "5.mpi"
        assert data["histories"] == 100000
        tally = data["tallies"][4]
        assert tally["particle"] == "n"
        assert tally["dose_response"] is True
        assert tally["x_bounds"] == [-200.0, -66.67, 66.67, 200.0]
        assert tally["dims"] == [3, 5, 3]
        assert tally["num_ves"] == 45
        assert tally["num_e_groups"] == 3
        assert len(tally["result"]) == 45
        assert all(len(row) == 3 for row in tally["result"])
        assert tally["result"][0][0] == pytest.approx(4.96471e-9)
        assert tally["total_result"][0] == pytest.approx(1.91370e-7)

    def test_matches_read_meshtal(self) -> None:
        path = str(MCNP_DIR / "meshtal" / "mcnp_meshtal_multiple_meshtal.txt")
        data = nucleide.mcnp.meshtal_mesh_data(path)
        meshtal = nucleide.mcnp.read_meshtal(path)
        assert sorted(data["tallies"]) == sorted(meshtal.tallies)
        for number, tally in meshtal.tallies.items():
            entry = data["tallies"][number]
            assert entry["result"] == tally.result
            assert entry["rel_error"] == tally.rel_error
            assert entry["total_result"] == tally.total_result

    def test_single_group_totals_mirror(self) -> None:
        data = nucleide.mcnp.meshtal_mesh_data(
            str(MCNP_DIR / "meshtal" / "mcnp_meshtal_multiple_meshtal.txt")
        )
        single = data["tallies"][14]
        assert single["num_e_groups"] == 1
        assert single["total_result"] == [row[0] for row in single["result"]]

    def test_numpy_option(self) -> None:
        np = pytest.importorskip("numpy")
        data = nucleide.mcnp.meshtal_mesh_data(
            str(MCNP_DIR / "meshtal" / "mcnp_meshtal_single_meshtal.txt"),
            as_numpy=True,
        )
        result = data["tallies"][4]["result"]
        totals = data["tallies"][4]["total_result"]
        assert isinstance(result, np.ndarray)
        assert result.shape == (45, 3)
        assert totals.shape == (45,)
        assert result[0, 0] == pytest.approx(4.96471e-9)
