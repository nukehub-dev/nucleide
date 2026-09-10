"""ENDF/B-VIII.0 decay branches + isomer masses (license-free decay-data pack)."""

import math

import pytest

import nucleide


def test_k40_two_branches_sum_to_one() -> None:
    branches = nucleide.nuclei.decay_branches("K40")
    assert len(branches) == 2
    by_progeny = {prog: (bf, mode) for prog, bf, mode in branches}
    bf_ca, mode_ca = by_progeny["Ca40"]
    assert bf_ca == pytest.approx(0.8914, rel=1e-9)
    assert mode_ca == "beta-"
    bf_ar, mode_ar = by_progeny["Ar40"]
    assert bf_ar == pytest.approx(0.1086, rel=1e-9)
    assert mode_ar == "ec/beta+"
    assert math.fsum(bf for _, bf, _ in branches) == pytest.approx(1.0, abs=1e-9)
    # The same VIII.0 checkout pins the half-life table entry.
    assert nucleide.nuclei.half_life("K40") == pytest.approx(3.93839e16, rel=1e-6)


def test_branch_spots() -> None:
    assert nucleide.nuclei.decay_branches("Es254") == [("Bk250", 1.0, "alpha")]
    assert nucleide.nuclei.decay_branches("Ba137_m1") == [("Ba137", 1.0, "IT")]
    he8 = {prog: (bf, mode) for prog, bf, mode in nucleide.nuclei.decay_branches("He8")}
    bf_li8, mode_li8 = he8["Li8"]
    assert bf_li8 == pytest.approx(0.84, rel=1e-9)
    assert mode_li8 == "beta-"
    bf_li7, mode_li7 = he8["Li7"]
    assert bf_li7 == pytest.approx(0.16, rel=1e-9)
    assert mode_li7 == "beta-"
    # Es254m mixes alpha/beta-/EC/IT; the SF branch is dropped.
    es_m1 = nucleide.nuclei.decay_branches("Es254_m1")
    assert len(es_m1) == 4
    assert all(mode != "sf" for _, _, mode in es_m1)
    assert nucleide.nuclei.decay_branch_fraction("Es254_m1", "Fm254") == pytest.approx(
        0.98, rel=1e-9
    )


def test_branch_absent_stays_graceful() -> None:
    assert nucleide.nuclei.decay_branches("Fe56") == []
    assert nucleide.nuclei.decay_branches("Te123") == []
    assert nucleide.nuclei.decay_branches("Ca46") == []
    assert nucleide.nuclei.decay_branch_fraction("K40", "K40") is None
    assert nucleide.nuclei.decay_branch_fraction("Fe56", "Fe56") is None
    with pytest.raises(ValueError):
        nucleide.nuclei.decay_branches("NotANuclide")


def test_isomer_masses() -> None:
    ground = nucleide.nuclei.atomic_mass("Ba137")
    assert ground is not None
    excited = nucleide.nuclei.atomic_mass("Ba137_m1")
    assert excited is not None
    # m = m_ground + 661659 eV / 931.49410242 MeV/u.
    assert excited == pytest.approx(ground + 0.661659 / 931.49410242, abs=1e-9)
    assert excited > ground
    # Te123 keeps its stable ground mass; its isomer carries ELIS.
    te = nucleide.nuclei.atomic_mass("Te123")
    assert te is not None
    te_m1 = nucleide.nuclei.atomic_mass("Te123_m1")
    assert te_m1 is not None
    assert te_m1 > te


def test_normalize_nuclide_name() -> None:
    assert nucleide.nuclei.normalize_nuclide("241Pu") == "Pu241"
    assert nucleide.nuclei.normalize_nuclide("Pu-241") == "Pu241"
    assert nucleide.nuclei.normalize_nuclide("Ba137m") == "Ba137_m1"
    assert nucleide.nuclei.normalize_nuclide("Ba-137m") == "Ba137_m1"
    assert nucleide.nuclei.normalize_nuclide("Ir-192n") == "Ir192_m2"
    assert nucleide.nuclei.normalize_nuclide("40K") == "K40"
    with pytest.raises(ValueError):
        nucleide.nuclei.normalize_nuclide("NotANuclide")
