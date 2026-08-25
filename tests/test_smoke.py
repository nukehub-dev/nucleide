"""Python-side tests for the naming/data core (run after `maturin develop`)."""

import nucleide


def test_version() -> None:
    assert nucleide.__version__ == "0.1.0"


def test_nuclide_uranium235() -> None:
    u5 = nucleide.Nuclide("U235")
    assert u5.name == "U235"
    assert u5.nucid == 922350000
    assert u5.zzaaam == 922350
    assert (u5.z, u5.a, u5.state) == (92, 235, 0)


def test_nuclide_metastable() -> None:
    am = nucleide.Nuclide("Am242_m1")
    assert am.name == "Am242_m1"
    assert am.nucid == 952420001
    assert am.state == 1


def test_nuclide_dialects() -> None:
    u5 = nucleide.Nuclide("U235")
    assert u5.zaid == 92235
    assert u5.zzllaaam == "92-U-235"
    assert u5.serpent == "U-235"
    assert u5.alara == "u:235"
    assert u5.cinder == 2350920
    assert nucleide.from_zaid(95242).name == "Am242_m1"


def test_fluka_round_trip() -> None:
    assert nucleide.Nuclide("U235").fluka() == "235-U"
    assert nucleide.Nuclide("H1").fluka() == "HYDROG-1"
    assert nucleide.Nuclide("Li7").fluka() == "LITHIU-7"


def test_atomic_masses() -> None:
    h1 = nucleide.atomic_mass("H1")
    u5 = nucleide.atomic_mass("U235")
    assert h1 is not None
    assert u5 is not None
    assert abs(h1 - 1.007825031898) < 1e-11
    assert abs(u5 - 235.043928117) < 1e-8
    assert nucleide.atomic_mass(999999999) is None
    assert nucleide.Nuclide("Am242_m1").mass is None  # metastables: ground-state table
    by_id = nucleide.atomic_mass(10010000)
    assert by_id is not None
    assert abs(by_id - 1.007825031898) < 1e-11


def test_natural_abundance() -> None:
    u5 = nucleide.natural_abundance("U235")
    o16 = nucleide.natural_abundance("O16")
    fe56 = nucleide.Nuclide("Fe56").abundance
    assert u5 is not None
    assert o16 is not None
    assert fe56 is not None
    assert abs(u5 - 0.007204) < 1e-9
    assert abs(o16 - 0.9976206) < 1e-9
    assert nucleide.natural_abundance("Cf255") is None
    assert abs(fe56 - 0.91754) < 1e-9


def test_particle() -> None:
    n = nucleide.Particle("neutron")
    assert n.name == "Neutron"
    assert nucleide.Particle("g").name == "Photon"
    assert nucleide.Particle("n").mcnp() == "n"
    assert nucleide.Particle("proton").geant4() == "proton"
    assert nucleide.Particle(2112).name == "Neutron"


def test_rxname() -> None:
    fid = nucleide.rxname_id("fission")
    aid = nucleide.rxname_id("absorption")
    # Hashed reaction ids; numeric strings resolve via MT lookup
    assert fid == nucleide.rxname_id("18")
    assert aid == nucleide.rxname_id("27")
    assert nucleide.rxname_mt(fid) == 18
    assert nucleide.rxname_name(fid) == "fission"
