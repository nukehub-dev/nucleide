"""Python-side tests for the naming/data core (run after `maturin develop`)."""

import nucleide


def test_version() -> None:
    assert nucleide.__version__ == "0.1.0"


def test_nuclide_uranium235() -> None:
    u5 = nucleide.nuclei.Nuclide("U235")
    assert u5.name == "U235"
    assert u5.nucid == 922350000
    assert u5.zzaaam == 922350
    assert (u5.z, u5.a, u5.state) == (92, 235, 0)


def test_nuclide_metastable() -> None:
    am = nucleide.nuclei.Nuclide("Am242_m1")
    assert am.name == "Am242_m1"
    assert am.nucid == 952420001
    assert am.state == 1


def test_nuclide_dialects() -> None:
    u5 = nucleide.nuclei.Nuclide("U235")
    assert u5.zaid == 92235
    assert u5.zzllaaam == "92-U-235"
    assert u5.serpent == "U-235"
    assert u5.alara == "u:235"
    assert u5.cinder == 2350920
    assert nucleide.nuclei.from_zaid(95242).name == "Am242_m1"


def test_fluka_round_trip() -> None:
    assert nucleide.nuclei.Nuclide("U235").fluka() == "235-U"
    assert nucleide.nuclei.Nuclide("H1").fluka() == "HYDROG-1"
    assert nucleide.nuclei.Nuclide("Li7").fluka() == "LITHIU-7"


def test_atomic_masses() -> None:
    h1 = nucleide.nuclei.atomic_mass("H1")
    u5 = nucleide.nuclei.atomic_mass("U235")
    assert h1 is not None
    assert u5 is not None
    assert abs(h1 - 1.007825031898) < 1e-11
    assert abs(u5 - 235.043928117) < 1e-8
    assert nucleide.nuclei.atomic_mass(999999999) is None
    assert nucleide.nuclei.Nuclide("Am242_m1").mass is None  # metastables: ground-state table
    by_id = nucleide.nuclei.atomic_mass(10010000)
    assert by_id is not None
    assert abs(by_id - 1.007825031898) < 1e-11


def test_natural_abundance() -> None:
    u5 = nucleide.nuclei.natural_abundance("U235")
    o16 = nucleide.nuclei.natural_abundance("O16")
    fe56 = nucleide.nuclei.Nuclide("Fe56").abundance
    assert u5 is not None
    assert o16 is not None
    assert fe56 is not None
    assert abs(u5 - 0.007204) < 1e-9
    assert abs(o16 - 0.9976206) < 1e-9
    assert nucleide.nuclei.natural_abundance("Cf255") is None
    assert abs(fe56 - 0.91754) < 1e-9


def test_particle() -> None:
    n = nucleide.nuclei.Particle("neutron")
    assert n.name == "Neutron"
    assert nucleide.nuclei.Particle("g").name == "Photon"
    assert nucleide.nuclei.Particle("n").mcnp() == "n"
    assert nucleide.nuclei.Particle("proton").geant4() == "proton"
    assert nucleide.nuclei.Particle(2112).name == "Neutron"


def test_rxname() -> None:
    fid = nucleide.nuclei.rxname_id("fission")
    aid = nucleide.nuclei.rxname_id("absorption")
    # Hashed reaction ids; numeric strings resolve via MT lookup
    assert fid == nucleide.nuclei.rxname_id("18")
    assert aid == nucleide.nuclei.rxname_id("27")
    assert nucleide.nuclei.rxname_mt(fid) == 18
    assert nucleide.nuclei.rxname_name(fid) == "fission"


def test_to_xml() -> None:
    comp = nucleide.material.from_formula("H2O")
    xml = nucleide.material.to_xml(comp, "water", 1.0, "g/cm3")
    assert '<material name="water"' in xml
    assert 'density value="1.0" units="g/cm3"' in xml
    assert "H1" in xml
    assert "O16" in xml
