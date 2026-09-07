"""Python-side tests for the multi-step `deplete_series` driver (Stream C).

Uses only a synthetic 3-nuclide chain (A -> B -> C, pure decay) written to
``tmp_path`` — no licensed fixtures.
"""

import math
from pathlib import Path

import pytest

import nucleide

# A: lambda = 1e-6 /s, B: lambda = 1e-5 /s, C: stable.
T_HALF_A = math.log(2) / 1e-6
T_HALF_B = math.log(2) / 1e-5

CHAIN_XML = f"""<depletion_chain>
  <nuclide name="A" half_life="{T_HALF_A!r}">
    <decay type="beta" target="B" branching_ratio="1.0"/>
  </nuclide>
  <nuclide name="B" half_life="{T_HALF_B!r}">
    <decay type="beta" target="C" branching_ratio="1.0"/>
  </nuclide>
  <nuclide name="C"/>
</depletion_chain>
"""


def _write_chain(tmp_path: Path, text: str = CHAIN_XML) -> str:
    path = tmp_path / "chain_synth.xml"
    path.write_text(text)
    return str(path)


def test_series_length_and_times(tmp_path: Path) -> None:
    chain = nucleide.depletion.read_chain(_write_chain(tmp_path))
    dts = [1e5, 2e5, 5e4]
    out = nucleide.depletion.deplete_series(chain, {"A": 1e15}, dts)
    assert out["times"] == pytest.approx([1e5, 3e5, 3.5e5])
    assert len(out["atoms"]) == 3
    assert len(out["activity"]) == 3
    assert len(out["decay_heat"]) == 3
    assert set(out.keys()) == {"times", "atoms", "activity", "decay_heat"}


def test_series_conserves_atoms(tmp_path: Path) -> None:
    # Pure decay with no loss channels: total atoms conserved each step.
    chain = nucleide.depletion.read_chain(_write_chain(tmp_path))
    out = nucleide.depletion.deplete_series(chain, {"A": 1e15}, [1e5, 1e5, 1e5])
    for step in out["atoms"]:
        total = sum(step.values())
        assert total == pytest.approx(1e15, rel=1e-8)


def test_series_matches_single_step(tmp_path: Path) -> None:
    # The standalone driver threads single-step `deplete` calls forward.
    chain = nucleide.depletion.read_chain(_write_chain(tmp_path))
    n0 = {"A": 1e15, "B": 2e14}
    series = nucleide.depletion.deplete_series(chain, n0, [1e5])
    single = nucleide.depletion.deplete(chain, n0, 1e5, order=48)
    for name, value in single.items():
        assert series["atoms"][0][name] == pytest.approx(value, rel=1e-12)


def test_series_activity_spot(tmp_path: Path) -> None:
    chain = nucleide.depletion.read_chain(_write_chain(tmp_path))
    out = nucleide.depletion.deplete_series(chain, {"A": 1e15}, [1e5])
    lam_a = math.log(2) / T_HALF_A
    expected_a = 1e15 * math.exp(-lam_a * 1e5)
    assert out["atoms"][0]["A"] == pytest.approx(expected_a, rel=1e-8)
    assert out["activity"][0]["A"] == pytest.approx(lam_a * expected_a, rel=1e-8)
    # C is stable: contributes atoms but no activity.
    assert out["atoms"][0]["C"] > 0
    assert out["activity"][0]["C"] == 0.0
    # Synthetic A/B/C carry no decay energies: heat maps stay all-zero.
    assert set(out["decay_heat"][0]) == {"A", "B", "C"}
    assert all(v == 0.0 for v in out["decay_heat"][0].values())


def test_series_heat_uses_endf_table(tmp_path: Path) -> None:
    # A single Co60 nuclide: heat = A * E * J/MeV with E from MF8/MT457.
    chain_xml = (
        '<depletion_chain><nuclide name="Co60" '
        f'half_life="{60 * 60 * 24 * 365.25 * 5.2714!r}"/></depletion_chain>'
    )
    chain = nucleide.depletion.read_chain(_write_chain(tmp_path, chain_xml))
    out = nucleide.depletion.deplete_series(chain, {"Co60": 1e20}, [1e6])
    heat = out["decay_heat"][0]["Co60"]
    act = out["activity"][0]["Co60"]
    assert heat == pytest.approx(act * 2.60061 * 1.602e-13, rel=1e-6)
    assert heat > 0


def test_series_invalid_integrator_raises(tmp_path: Path) -> None:
    chain = nucleide.depletion.read_chain(_write_chain(tmp_path))
    with pytest.raises(ValueError, match="unsupported integrator"):
        nucleide.depletion.deplete_series(chain, {"A": 1.0}, [1.0], integrator="corrector")


def test_series_rates_list_length_mismatch_raises(tmp_path: Path) -> None:
    chain = nucleide.depletion.read_chain(_write_chain(tmp_path))
    with pytest.raises(ValueError, match="rates_list"):
        nucleide.depletion.deplete_series(chain, {"A": 1.0}, [1.0, 2.0], rates_list=[None])


def test_data_accessors_live(tmp_path: Path) -> None:
    xs = nucleide.nuclei.simple_xs("H1")
    assert xs is not None
    thermal, fast = xs
    assert thermal == pytest.approx(20.84, rel=0.05)
    assert fast == pytest.approx(0.687, rel=0.05)
    assert nucleide.nuclei.scattering_length("H1") == pytest.approx(-3.7406, abs=0.01)
    # Cs137 prompt mean excludes the Ba137_m1 daughter gamma.
    assert nucleide.nuclei.decay_energy("Cs137") == pytest.approx(0.1794, rel=0.05)
    assert nucleide.nuclei.decay_energy("Ba137_m1") == pytest.approx(0.6614, rel=0.05)
    # Documented coverage gaps stay graceful.
    assert nucleide.nuclei.simple_xs("Cs137") is None
    assert nucleide.nuclei.simple_xs("Og294") is None
    assert nucleide.nuclei.decay_energy("Fe56") is None
    # 1 g Co60 ≈ 17 W (evaluation-backed band).
    assert nucleide.material.decay_heat({"Co60": 1.0}) == pytest.approx(17.0, rel=0.3)
    with pytest.raises(ValueError):
        nucleide.nuclei.simple_xs("NotANuclide")
