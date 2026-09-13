"""Thin facade tests for the reader facade bundle + missing thin-facade paths.

Every wrapper here is a thin pyfn over existing Rust (no new math/data);
assertions use synthetic inline texts or committed fixtures only.
"""

from pathlib import Path

import pytest

import nucleide

FIX = Path(__file__).parent.parent / "fixtures"


class TestRxnameGraph:
    def test_id_name_mt_round_trip(self) -> None:
        fid = nucleide.nuclei.rxname_id("fission")
        assert nucleide.nuclei.rxname_name(fid) == "fission"
        assert nucleide.nuclei.rxname_mt(fid) == 18
        # Numeric MT strings resolve through the same registry.
        assert nucleide.nuclei.rxname_id("18") == fid
        assert nucleide.nuclei.rxname_id("absorption") == nucleide.nuclei.rxname_id("27")

    def test_label_doc_reaction(self) -> None:
        fid = nucleide.nuclei.rxname_id("fission")
        assert nucleide.nuclei.rxname_label(fid) == "(z,fission)"
        assert "fiss" in nucleide.nuclei.rxname_doc(fid).lower()
        row = nucleide.nuclei.rxname_reaction(fid)
        assert row is not None
        assert row["name"] == "fission" and row["mt"] == 18
        assert row["label"] == nucleide.nuclei.rxname_label(fid)
        assert nucleide.nuclei.rxname_reaction(0xFFFFFFFF) is None
        assert nucleide.nuclei.rxname_label(0xFFFFFFFF) == ""
        assert nucleide.nuclei.rxname_doc(0xFFFFFFFF) == ""

    def test_nucdelta_child_parent(self) -> None:
        u235 = nucleide.nuclei.Nuclide("U235").nucid
        u236 = nucleide.nuclei.Nuclide("U236").nucid
        absorp = nucleide.nuclei.rxname_id("absorption")
        assert nucleide.nuclei.rxname_id_from_nucdelta(u235, u236, "n") == absorp
        assert nucleide.nuclei.rxname_child("U235", "absorption") == "U236"
        assert nucleide.nuclei.rxname_child("U235", absorp, "n") == "U236"
        assert nucleide.nuclei.rxname_parent("U236", "absorption") == "U235"
        with pytest.raises(ValueError):
            nucleide.nuclei.rxname_child("U235", "absorption", projectile="bogus")


class TestParticleHelpers:
    def test_validity_helpers(self) -> None:
        assert nucleide.nuclei.particle_is_valid("n")
        assert nucleide.nuclei.particle_is_valid("neutron")
        assert nucleide.nuclei.particle_is_valid("U235")  # heavy ion counts
        assert not nucleide.nuclei.particle_is_valid("not-a-particle")
        assert nucleide.nuclei.particle_is_valid_pdc(2112)
        assert not nucleide.nuclei.particle_is_valid_pdc(999999)
        assert nucleide.nuclei.particle_is_hydrogen("Proton")
        assert nucleide.nuclei.particle_is_hydrogen("H1")
        assert not nucleide.nuclei.particle_is_hydrogen("U235")
        assert nucleide.nuclei.particle_is_heavy_ion("U235")
        assert not nucleide.nuclei.particle_is_heavy_ion("Proton")

    def test_particle_class_paths(self) -> None:
        assert nucleide.nuclei.Particle("n").name == "Neutron"
        assert nucleide.nuclei.Particle(2112).name == "Neutron"
        assert nucleide.nuclei.Particle("g").mcnp() == "p"


class TestNuclearDataPaths:
    def test_simple_xs_scattering_decay_energy(self) -> None:
        xs = nucleide.nuclei.simple_xs("H1")
        assert xs is not None
        thermal, fast = xs
        assert thermal == pytest.approx(20.8401, rel=1e-6)
        assert fast == pytest.approx(0.687144, rel=1e-6)
        assert nucleide.nuclei.scattering_length("H1") == pytest.approx(-3.7406, rel=1e-6)
        assert nucleide.nuclei.decay_energy("Co60") == pytest.approx(2.60061, rel=1e-6)
        # Prompt-only convention: Cs137 reports the parent share.
        assert nucleide.nuclei.decay_energy("Cs137") == pytest.approx(0.179448, rel=1e-4)


class TestDoseExtras:
    def test_f1_and_lung_model(self) -> None:
        assert nucleide.nuclei.dose_f1("H3") == pytest.approx(1.0)
        assert nucleide.nuclei.dose_f1("Co60") == pytest.approx(0.3)
        assert nucleide.nuclei.dose_lung_model("Co60") == "Y"
        assert nucleide.nuclei.dose_lung_model("H3") == "V"
        assert nucleide.nuclei.dose_f1("U235", source="DOE") is None or isinstance(
            nucleide.nuclei.dose_f1("U235", source="DOE"), float
        )


class TestMaterialMix:
    def test_mix_by_mass(self) -> None:
        out = nucleide.material.mix_by_mass([({"U235": 19.0}, 1.0), ({"U238": 1.0}, 1.0)])
        assert out == {"U235": 19.0, "U238": 1.0}
        with pytest.raises(ValueError):
            nucleide.material.mix_by_mass([({"U235": 1.0}, -1.0)])

    def test_mix_by_volume(self) -> None:
        out = nucleide.material.mix_by_volume(
            [({"U235": 19.0}, 1.0, 19.1), ({"U238": 1.0}, 1.0, 19.1)]
        )
        assert set(out) == {"U235", "U238"}
        with pytest.raises(ValueError):
            nucleide.material.mix_by_volume([({"U235": 1.0}, 1.0, 0.0)])

    def test_specific_activity(self) -> None:
        assert nucleide.material.specific_activity({"U235": 1.0}) == pytest.approx(
            79960.38, rel=1e-4
        )

    def test_materials_doc(self) -> None:
        xml = nucleide.material.materials_doc_to_xml(
            [("fuel", {"U235": 19.0}, 10.0)], cross_sections="xs.xml"
        )
        assert "<materials" in xml and "fuel" in xml and "xs.xml" in xml

    def test_expand_collapse(self) -> None:
        expanded = nucleide.material.expand_elements({"U": 20.0})
        assert set(expanded) == {"U234", "U235", "U238"}
        assert sum(expanded.values()) == pytest.approx(20.0, rel=1e-9)
        collapsed = nucleide.material.collapse_elements({"U235": 19.0, "U238": 1.0})
        assert collapsed == {"U": 20.0}


class TestFlukaStrings:
    def test_material_str(self) -> None:
        card = nucleide.fluka.fluka_material_str(1, "U235", 19.1)
        assert card.startswith("MATERIAL") and "235-U" in card
        # Builtin elements need no card (empty string, matching upstream).
        assert nucleide.fluka.fluka_material_str(1, "N", 1.0) == ""
        with pytest.raises(ValueError):
            nucleide.fluka.fluka_material_str(1, "not-a-nuclide", 1.0)

    def test_compound_str(self) -> None:
        card = nucleide.fluka.fluka_compound_str(
            1, "UOX", 10.0, "mass", [("U235", 0.5), ("U238", 0.5)]
        )
        assert "MATERIAL" in card and "COMPOUND" in card and "UOX" in card
        with pytest.raises(ValueError):
            nucleide.fluka.fluka_compound_str(1, "EMPTY", 1.0, "mass", [])

    def test_builtin_set(self) -> None:
        builtins = nucleide.fluka.fluka_builtin_set()
        assert len(builtins) == 37
        assert "WATER" in builtins and builtins == sorted(builtins)


class TestAlaraTotals:
    def test_validate_and_check_block(self) -> None:
        deck = (FIX / "alara" / "decks" / "sample2").read_text()
        nucleide.alara.alara_validate_deck(deck)
        nucleide.alara.alara_check_block("geometry", 1)
        with pytest.raises(ValueError):
            nucleide.alara.alara_check_block("bogus", 3)

    def test_flux_totals(self) -> None:
        assert nucleide.alara.alara_flux_total("f", "1.0 2.0 3.0") == pytest.approx(6.0)
        assert nucleide.alara.alara_flux_len("f", "1.0 2.0 3.0") == 3

    def test_output_totals(self) -> None:
        text = (FIX / "alara" / "output" / "sample2.out").read_text()
        total = nucleide.alara.alara_output_total_activity(text, "sample")
        assert total > 0
        totals = nucleide.alara.alara_output_totals(text, "sample")
        assert len(totals) > 0
        assert all(r["nuclide"] == "total" for r in totals)

    def test_photon_and_schedule_totals(self) -> None:
        assert nucleide.alara.alara_photon_total_strength("U235 1.0 d 1.0 2.0") == pytest.approx(
            3.0
        )
        deck = (FIX / "alara" / "decks" / "sample2").read_text()
        assert nucleide.alara.alara_schedule_total_time(deck) == pytest.approx(31557600.0)


class TestOrigenFind:
    TAPE6 = "U235 1.0 100.0\nPu239 2.0 200.0\n"
    TAPE9 = "U235 1e-9\nPu239 2e-9\n"

    def test_tape6_find_and_total(self) -> None:
        row = nucleide.origen.origen_tape6_find(self.TAPE6, "u235")
        assert row is not None and row["grams"] == pytest.approx(1.0)
        assert nucleide.origen.origen_tape6_total_activity(self.TAPE6) == pytest.approx(300.0)
        assert nucleide.origen.origen_tape6_find(self.TAPE6, "Cm244") is None

    def test_tape9_find(self) -> None:
        row = nucleide.origen.origen_tape9_find(self.TAPE9, "PU239")
        assert row is not None and row["decay_const"] == pytest.approx(2e-9)
        assert nucleide.origen.origen_tape9_find(self.TAPE9, "Cm244") is None


class TestCcccAccessors:
    def test_rtflux_accessors(self) -> None:
        text = (FIX / "cccc" / "rtflux_sample").read_text()
        assert nucleide.cccc.cccc_rtflux_npoints(text) == 2
        assert nucleide.cccc.cccc_rtflux_total(text) == pytest.approx(21.0)
        assert nucleide.cccc.cccc_rtflux_point(text, "rtflux", 0) == [1.0, 2.0, 3.0]
        assert nucleide.cccc.cccc_rtflux_point(text, "rtflux", 99) is None

    def test_isotxs_accessors(self) -> None:
        text = (FIX / "cccc" / "isotxs_sample").read_text()
        assert nucleide.cccc.cccc_isotxs_len(text) == 2
        row = nucleide.cccc.cccc_isotxs_find(text, "U235")
        assert row is not None and row["groups"] == 3
        assert nucleide.cccc.cccc_isotxs_find(text, "missing") is None


class TestFispactSniffer:
    def test_suffix_convention(self) -> None:
        assert nucleide.fispact.fispact_is_output("inventory.fis") is True
        assert nucleide.fispact.fispact_is_output("inventory.out") is False


class TestEnrichmentRatios:
    def test_mass_ratios_close(self) -> None:
        f, p, t = 0.0072, 0.05, 0.002
        prod_feed = nucleide.enrichment.prod_per_feed(f, p, t)
        assert prod_feed == pytest.approx((f - t) / (p - t))
        assert prod_feed + nucleide.enrichment.tail_per_feed(f, p, t) == pytest.approx(1.0)
        assert nucleide.enrichment.feed_per_prod(f, p, t) == pytest.approx(1.0 / prod_feed)
        assert nucleide.enrichment.tail_per_prod(f, p, t) == pytest.approx(
            nucleide.enrichment.tail_per_feed(f, p, t) / prod_feed
        )
        assert nucleide.enrichment.prod_per_tail(f, p, t) == pytest.approx(
            prod_feed / nucleide.enrichment.tail_per_feed(f, p, t)
        )
        assert nucleide.enrichment.feed_per_tail(f, p, t) == pytest.approx(
            1.0 / nucleide.enrichment.tail_per_feed(f, p, t)
        )

    def test_alphastar(self) -> None:
        assert nucleide.enrichment.alphastar_i(1.05, 236.0, 235.0) == pytest.approx(1.05)


class TestKineticsFromIfp:
    def test_valid_params(self) -> None:
        params = nucleide.kinetics.from_ifp([0.002, 0.004], 1e-4, [0.01, 0.1])
        assert params["groups"] == 2
        assert params["beta_total"] == pytest.approx(0.006)
        assert params["lambda_gen"] == pytest.approx(1e-4)

    def test_invalid_raises(self) -> None:
        with pytest.raises(ValueError):
            nucleide.kinetics.from_ifp([0.002], 1e-4, [0.01, 0.1])


class TestVrPaths:
    def setup_method(self) -> None:
        meshtal = nucleide.mcnp.read_meshtal(
            str(FIX / "mcnp" / "meshtal" / "mcnp_meshtal_single_meshtal.txt")
        )
        self.tally = meshtal.tallies[4]

    def test_magic_with_selections(self) -> None:
        total = nucleide.vr.magic_with(self.tally, "total", 0.5, 0.0)
        assert total.groups_per_ve == 1
        assert len(total.lower_bounds_ww) == self.tally.num_ves()
        grouped = nucleide.vr.magic_with(self.tally, "per_group")
        assert grouped.groups_per_ve == 3
        assert len(grouped.lower_bounds_ww) == self.tally.num_ves() * 3
        with pytest.raises(ValueError):
            nucleide.vr.magic_with(self.tally, "bogus")

    def test_sampler_accessors(self) -> None:
        sampler = nucleide.vr.MeshSourceSampler(self.tally, "analog")
        assert sampler.mode() == "analog"
        assert sampler.num_voxels() == self.tally.num_ves()
        assert sampler.table_len() == self.tally.num_ves()
        uniform = nucleide.vr.MeshSourceSampler(self.tally, "uniform")
        assert uniform.mode() == "uniform"
        user = nucleide.vr.MeshSourceSampler(
            self.tally, "user", user_pdf=[1.0] * self.tally.num_ves()
        )
        assert user.mode() == "user"


class TestMcplStatsum:
    def test_round_trip(self) -> None:
        comment = nucleide.mcpl.mcpl_statsum_comment("nps", 1.0)
        assert nucleide.mcpl.mcpl_statsum_validate(comment) == "nps"

    def test_bad_comment_raises(self) -> None:
        with pytest.raises(ValueError):
            nucleide.mcpl.mcpl_statsum_validate("not-a-stat-comment")


class TestMctalMeshTfc:
    def test_mesh_fixture(self) -> None:
        m = nucleide.mcnp.read_mctal(str(FIX / "mcnp" / "mctal" / "synthetic_mesh.mctal"))
        assert m.tally_nums == [4]
        assert m.tallies == []
        assert len(m.mesh_tallies) == 1
        body = m.mesh_tallies[0]
        assert body["detector_type"] == -1
        assert body["dims"] == [2, 1, 1]
        assert body["cora"] == [0.0, 5.0, 10.0]
        assert body["corb"] == [0.0, 10.0]
        assert body["corc"] == [0.0, 10.0]
        assert body["vals"] == [(11.0, 1.375), (12.0, 1.5)]
        assert body["total"] == pytest.approx(23.0)

    def test_tfc_variant_fixture(self) -> None:
        m = nucleide.mcnp.read_mctal(str(FIX / "mcnp" / "mctal" / "synthetic_tfc_variants.mctal"))
        assert len(m.tallies) == 1
        assert m.mesh_tallies == []
        tally = m.tallies[0]
        assert tally["e"]["variant"] == "t"
        assert tally["e"]["values"] == [0.5, 2.0]
        assert tally["t"]["flag"] == 0
        assert tally["vals"] == [(11.0, 1.375), (12.0, 1.5)]
        tfc = tally["tfc"]
        assert tfc is not None
        assert tfc["jtf"] == [1, 2, 3, 4, 5, 6, 7, 8, 9]
        assert len(tfc["rows"]) == 2
        assert tfc["rows"][0]["fom"] is None
        assert tfc["rows"][1]["fom"] == pytest.approx(3.5)
