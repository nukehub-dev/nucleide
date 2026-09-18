"""Python-side tests for the equilibrium data readers (gates + input errors)."""

import math
import struct
from pathlib import Path

import pytest

import nucleide.equilib as eq

HDF5_SIG = bytes([0x89, 0x48, 0x44, 0x46, 0x0D, 0x0A, 0x1A, 0x0A])


def _pad(buf: bytearray, name: str) -> None:
    buf += struct.pack(">I", len(name)) + name.encode("ascii")
    while len(buf) % 4:
        buf.append(0)


def synthetic_wout(*, cdf2: bool = False, with_nyq: bool = True, drop: str = "") -> bytes:
    """Hand-packed synthetic `wout` (published var lists only, never solver output)."""
    dims = [("radius", 3), ("mn_mode", 2), ("mn_mode_nyq", 2)]

    def ints(vals: list[int]) -> bytes:
        return struct.pack(f">{len(vals)}i", *vals)

    def doubles(vals: list[float]) -> bytes:
        return struct.pack(f">{len(vals)}d", *vals)

    varlist = [
        ("nfp", [], 4, ints([3])),
        ("ns", [], 4, ints([3])),
        ("mpol", [], 4, ints([5])),
        ("ntor", [], 4, ints([2])),
        ("phiedge", [], 6, doubles([1.5])),
        ("volume_p", [], 6, doubles([30.0])),
        ("xm", [1], 6, doubles([0.0, 2.0])),
        ("xn", [1], 6, doubles([0.0, -1.0])),
        ("rmnc", [0, 1], 6, doubles([1.0, 0.1, 1.1, 0.1, 1.2, 0.1])),
        ("zmns", [0, 1], 6, doubles([0.0, 0.2, 0.0, 0.25, 0.0, 0.3])),
        ("lmns", [0, 1], 6, doubles([0.0, 0.01, 0.0, 0.02, 0.0, 0.03])),
        ("gmnc", [0, 2], 6, doubles([2.0, 0.0, 1.0, 0.5, 1.0, -1.0])),
        ("bmnc", [0, 1], 6, doubles([5.0, 0.0, 5.1, 0.0, 5.2, 0.0])),
    ]
    if with_nyq:
        varlist += [
            ("xm_nyq", [2], 6, doubles([0.0, 2.0])),
            ("xn_nyq", [2], 6, doubles([0.0, -1.0])),
        ]
    varlist = [v for v in varlist if v[0] != drop]
    head = bytearray(b"CDF\x02" if cdf2 else b"CDF\x01")
    head += struct.pack(">Q" if cdf2 else ">I", 0)
    head += struct.pack(">II", 10, len(dims))
    for name, length in dims:
        _pad(head, name)
        head += struct.pack(">Q" if cdf2 else ">I", length)
    head += struct.pack(">II", 0, 0)
    head += struct.pack(">II", 11, len(varlist))
    header_len = len(head)
    for name, dids, _code, _payload in varlist:
        header_len += 4 + len(name) + (-len(name) % 4)
        header_len += 4 + 4 * len(dids) + 8 + 4 + (16 if cdf2 else 8)
    begin = header_len
    offsets = []
    for _name, _dids, _code, payload in varlist:
        offsets.append(begin)
        begin += len(payload) + (-len(payload) % 4)
    for (name, dids, code, payload), off in zip(varlist, offsets, strict=True):
        _pad(head, name)
        head += struct.pack(">I", len(dids))
        for d in dids:
            head += struct.pack(">I", d)
        head += struct.pack(">II", 0, 0)
        if cdf2:
            head += struct.pack(">IQQ", code, len(payload), off)
        else:
            head += struct.pack(">III", code, len(payload), off)
    assert len(head) == header_len
    out = bytes(head)
    for (_name, _dids, _code, payload), _off in zip(varlist, offsets, strict=True):
        out += payload + b"\x00" * (-len(payload) % 4)
    return out


def _write(tmp_path: Path, name: str, blob: bytes) -> str:
    path = tmp_path / name
    path.write_bytes(blob)
    return str(path)


def test_probe_accepts_cdf1_and_cdf2(tmp_path: Path) -> None:
    assert eq.probe_variant(_write(tmp_path, "a.nc", synthetic_wout())) == "classic"
    assert eq.probe_variant(_write(tmp_path, "b.nc", synthetic_wout(cdf2=True))) == "64-bit offset"


def test_probe_rejects_wrong_variants_loudly(tmp_path: Path) -> None:
    bad_hdf5 = _write(tmp_path, "h.nc", HDF5_SIG)
    with pytest.raises(ValueError, match="facade-side"):
        eq.probe_variant(bad_hdf5)
    bad_cdf5 = _write(tmp_path, "c.nc", b"CDF\x05\x00\x00\x00\x00")
    with pytest.raises(ValueError, match="facade-side"):
        eq.probe_variant(bad_cdf5)
    for blob in (b"NOPE", b"CD", b"CDF", b"CDF\x03"):
        bad = _write(tmp_path, "g.nc", blob)
        with pytest.raises(ValueError):
            eq.probe_variant(bad)


def test_read_wout_round_trip(tmp_path: Path) -> None:
    w = eq.read_wout(_write(tmp_path, "w.nc", synthetic_wout()))
    assert w["variant"] == "classic"
    assert (w["nfp"], w["ns"], w["n_modes"], w["n_modes_nyq"]) == (3, 3, 2, 2)
    assert (w["mpol"], w["ntor"], w["phiedge"], w["volume_p"]) == (5, 2, 1.5, 30.0)
    assert w["xm"] == [0.0, 2.0]
    assert w["xn"] == [0.0, -1.0]
    assert w["xm_nyq"] == [0.0, 2.0]
    assert w["rmnc"] == {"rows": 3, "cols": 2, "values": [[1.0, 0.1], [1.1, 0.1], [1.2, 0.1]]}
    assert w["gmnc"]["values"] == [[2.0, 0.0], [1.0, 0.5], [1.0, -1.0]]
    assert w["fields"]["bmnc"]["values"][0] == [5.0, 0.0]
    # CDF-2 reads identically.
    w2 = eq.read_wout(_write(tmp_path, "w2.nc", synthetic_wout(cdf2=True)))
    assert w2["variant"] == "64-bit offset"
    assert w2["gmnc"] == w["gmnc"]


def test_read_wout_missing_variable_fails(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="gmnc"):
        eq.read_wout(_write(tmp_path, "m.nc", synthetic_wout(drop="gmnc")))


def test_jacobian_hand_vectors() -> None:
    # Pure (0,0) row: identically 2.0.
    assert eq.jacobian([0.0], [0.0], [2.0], 0.0, 0.0) == 2.0
    assert eq.jacobian([0.0], [0.0], [2.0], 1.3, 2.1) == 2.0
    # 1.0 + 0.5*cos(2θ + ζ).
    assert eq.jacobian([0.0, 2.0], [0.0, -1.0], [1.0, 0.5], 0.0, 0.0) == pytest.approx(1.5)
    assert eq.jacobian([0.0, 2.0], [0.0, -1.0], [1.0, 0.5], math.pi / 2, 0.0) == pytest.approx(0.5)
    assert eq.jacobian([0.0, 2.0], [0.0, -1.0], [1.0, 0.5], 0.0, math.pi) == pytest.approx(0.5)
    # Shape and non-finite guards.
    with pytest.raises(ValueError):
        eq.jacobian([0.0], [0.0, 1.0], [1.0], 0.0, 0.0)
    with pytest.raises(ValueError):
        eq.jacobian([0.0], [0.0], [float("nan")], 0.0, 0.0)


def test_jacobian_grid_contract() -> None:
    grid = eq.jacobian_grid([0.0], [0.0], [1.0], 3, 4, 5)
    assert len(grid) == 4
    assert all(len(row) == 5 for row in grid)
    assert all(v == 1.0 for row in grid for v in row)
    two = eq.jacobian_grid([0.0, 2.0], [0.0, -1.0], [1.0, 0.5], 3, 4, 1)
    assert two[0][0] == pytest.approx(1.5)
    assert two[1][0] == pytest.approx(0.5)
    with pytest.raises(ValueError):
        eq.jacobian_grid([0.0], [0.0], [1.0], 0, 4, 5)


INDATA = """
! Comment line.
&INDATA
  LFREEB = F
  MGRID_FILE = 'none'
  DELT = 0.9
  NFP = 3
  NCURR = 1
  MPOL = 5 NTOR = 4
  PMASS_TYPE = 'power_series'
  AM = 0.0 1.0D0 3*0.0
  AI(2) = 1.5E-3
  RBC(0,0) = 1.0
  RBC(1,0) = -2.5D-1
  ZBS(1,0) = 2.5e-1
/
"""


def test_indata_grammar_vectors() -> None:
    parsed = eq.parse_indata(INDATA)
    assert eq.indata_int(parsed, "nfp") == 3
    assert eq.indata_int(parsed, "NFP") == 3
    assert eq.indata_int(parsed, "mpol") == 5
    assert eq.indata_float(parsed, "delt") == pytest.approx(0.9)
    assert eq.indata_bool(parsed, "lfreeb") is False
    assert eq.indata_str(parsed, "mgrid_file") == "none"
    assert eq.indata_str(parsed, "pmass_type") == "power_series"
    assert eq.indata_series_1d(parsed, "am") == [(0, 0.0), (1, 1.0), (2, 0.0), (3, 0.0), (4, 0.0)]
    assert eq.indata_series_1d(parsed, "ai") == [(2, 1.5e-3)]
    assert eq.indata_table_2d(parsed, "rbc") == [((0, 0), 1.0), ((1, 0), -0.25)]
    assert eq.indata_table_2d(parsed, "zbs") == [((1, 0), 0.25)]
    assert eq.is_free_boundary(parsed) is False
    assert eq.ncurr_is_iprime(parsed) is True


def test_indata_rejections() -> None:
    with pytest.raises(ValueError, match="sliced"):
        eq.parse_indata("&INDATA\nRBC(0:4,2) = 1.0\n/\n")
    with pytest.raises(ValueError, match="INDATA"):
        eq.parse_indata("NFP = 3\n")
    with pytest.raises(ValueError, match="terminated"):
        eq.parse_indata("&INDATA\nNFP = 3\n")
    parsed = eq.parse_indata("&INDATA\nA = T\nB = 1.5\nC = 3.0\nD = 3.5\n/\n")
    with pytest.raises(ValueError, match="integer"):
        eq.indata_int(parsed, "a")
    with pytest.raises(ValueError, match="missing"):
        eq.indata_int(parsed, "zzz")
    assert eq.indata_int(parsed, "c") == 3
    with pytest.raises(ValueError, match="integer"):
        eq.indata_int(parsed, "b")
    with pytest.raises(ValueError, match="integer"):
        eq.indata_int(parsed, "d")
    with pytest.raises(ValueError, match="integer"):
        eq.indata_int(eq.parse_indata("&INDATA\nS = 'x'\n/\n"), "s")
    with pytest.raises(ValueError, match="logical"):
        eq.indata_bool(parsed, "b")
    with pytest.raises(ValueError, match="missing"):
        eq.indata_bool(parsed, "zzz")
    with pytest.raises(ValueError, match="float"):
        eq.indata_float(parsed, "a")
    with pytest.raises(ValueError, match="missing"):
        eq.indata_float(parsed, "zzz")
    assert eq.indata_float(parsed, "b") == pytest.approx(1.5)
    with pytest.raises(ValueError, match="string"):
        eq.indata_str(parsed, "a")
    with pytest.raises(ValueError, match="missing"):
        eq.indata_str(parsed, "zzz")
    with pytest.raises(ValueError, match="missing"):
        eq.indata_series_1d(parsed, "am")
    with pytest.raises(ValueError, match="missing"):
        eq.indata_table_2d(parsed, "rbc")
    series = eq.parse_indata("&INDATA\nAM(0) = 1.0\nBB(1) = T\nCC(1,2) = 1.0\n/\n")
    assert eq.indata_series_1d(series, "am") == [(0, 1.0)]
    with pytest.raises(ValueError, match="1-D float"):
        eq.indata_series_1d(series, "bb")
    with pytest.raises(ValueError, match="1-D float"):
        eq.indata_series_1d(series, "cc")
    with pytest.raises(ValueError, match="2-D float"):
        eq.indata_table_2d(series, "bb")
    with pytest.raises(ValueError, match="2-D float"):
        eq.indata_table_2d(series, "am")
    assert eq.indata_table_2d(series, "cc") == [((1, 2), 1.0)]
    with pytest.raises(ValueError, match="logical"):
        eq.is_free_boundary(eq.parse_indata("&INDATA\nLFREEB = 1\n/\n"))
    assert eq.is_free_boundary(eq.parse_indata("&INDATA\n/\n")) is False
    assert eq.is_free_boundary(eq.parse_indata("&INDATA\nLFREEB = .TRUE.\n/\n")) is True
    assert eq.ncurr_is_iprime(eq.parse_indata("&INDATA\nNCURR = 0\n/\n")) is False
    assert eq.ncurr_is_iprime(eq.parse_indata("&INDATA\n/\n")) is None


def test_wall_load_axisymmetric_matches_analytic() -> None:
    # (W3): S = 2, J0 = 3 over s in [0, 1] give q = 6 everywhere;
    # Q = 12 * pi^2 over one field period with nfp = 2.
    out = eq.wall_load([0.0, 0.5, 1.0], 4, 3, 2, [2.0] * 24, [3.0] * 24)
    assert (out["ntheta"], out["nzeta"]) == (4, 3)
    assert len(out["loads"]) == 4
    assert all(len(row) == 3 for row in out["loads"])
    assert all(v == pytest.approx(6.0) for row in out["loads"] for v in row)
    assert out["total"] == pytest.approx(12.0 * math.pi**2)


def test_wall_load_cosine_and_conservation() -> None:
    # J = 1 + 0.5*cos(2θ) at the cardinal nodes, S = 1, one radial cell.
    out = eq.wall_load([0.0, 1.0], 4, 1, 1, [1.0] * 4, [1.5, 0.5, 1.5, 0.5])
    assert [row[0] for row in out["loads"]] == pytest.approx([1.5, 0.5, 1.5, 0.5])
    # Conservation: total equals the direct voxel sum.
    s_edges = [0.0, 0.25, 1.0]
    birth = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    jac = [1.0, 0.5, 2.0, 1.5, 1.0, 1.0, 0.25, 2.0]
    out = eq.wall_load(s_edges, 2, 2, 3, birth, jac)
    dtheta = 2.0 * math.pi / 2.0
    dzeta = 2.0 * math.pi / (2.0 * 3.0)
    want = sum(
        birth[i * 4 + v] * jac[i * 4 + v] * (s_edges[i + 1] - s_edges[i]) * dtheta * dzeta
        for i in range(2)
        for v in range(4)
    )
    assert out["total"] == pytest.approx(want)


def test_wall_load_rejects_malformed() -> None:
    with pytest.raises(ValueError, match="strictly increasing"):
        eq.wall_load([0.0, 1.0, 0.5], 2, 1, 1, [1.0] * 4, [1.0] * 4)
    with pytest.raises(ValueError, match="positive"):
        eq.wall_load([0.0, 1.0], 0, 1, 1, [1.0] * 2, [1.0] * 2)
    with pytest.raises(ValueError, match="positive"):
        eq.wall_load([0.0, 1.0], 2, 1, 0, [1.0] * 2, [1.0] * 2)
    with pytest.raises(ValueError, match="birth"):
        eq.wall_load([0.0, 1.0], 2, 1, 1, [1.0], [1.0] * 2)
    with pytest.raises(ValueError, match="negative"):
        eq.wall_load([0.0, 1.0], 2, 1, 1, [-1.0, 1.0], [1.0] * 2)
    with pytest.raises(ValueError, match="[Jj]acobian"):
        eq.wall_load([0.0, 1.0], 2, 1, 1, [1.0] * 2, [1.0, -0.5])
