"""Dose factors and dose-per-gram (HNF-5636/PyNE tables, real equations)."""

import math

import pytest

import nucleide


class TestDoseFactors:
    def test_spot_factors(self) -> None:
        assert nucleide.nuclei.dose_factor("Co60", "ingest", "EPA") == pytest.approx(
            2.69e-05, rel=1e-9
        )
        assert nucleide.nuclei.dose_factor("Cs137", "inhale", "EPA") == pytest.approx(
            3.19e-05, rel=1e-9
        )
        assert nucleide.nuclei.dose_factor("H3", "air", "EPA") == pytest.approx(4.41e-012, rel=1e-9)
        assert nucleide.nuclei.dose_factor("K40", "soil", "EPA") == pytest.approx(4.33e02, rel=1e-9)

    def test_sources_and_plus_d_fold(self) -> None:
        # Cs-137+D folds into Cs137; GENII/DOE give parallel evaluations.
        epa = nucleide.nuclei.dose_factor("Cs137", "inhale", "EPA")
        genii = nucleide.nuclei.dose_factor("Cs137", "inhale", "GENII")
        doe = nucleide.nuclei.dose_factor("Cs137", "inhale", "DOE")
        assert epa == pytest.approx(3.19e-05, rel=1e-9)
        assert genii == pytest.approx(2.98e-05, rel=1e-9)
        assert doe == pytest.approx(3.20e-05, rel=1e-9)

    def test_missing_air_is_minus_one(self) -> None:
        assert nucleide.nuclei.dose_factor("H3", "air", "GENII") == -1.0
        assert nucleide.nuclei.dose_factor("H3", "air", "DOE") == -1.0
        assert nucleide.nuclei.dose_factor("Fe56", "ingest", "EPA") is None

    def test_bad_inputs(self) -> None:
        with pytest.raises(ValueError):
            nucleide.nuclei.dose_factor("Co60", "nope", "EPA")
        with pytest.raises(ValueError):
            nucleide.nuclei.dose_factor("Co60", "ingest", "nope")
        with pytest.raises(ValueError):
            nucleide.nuclei.dose_factor("Notanuclide", "ingest", "EPA")


class TestDosePerG:
    def test_co60_ingest_matches_pyne_equation(self) -> None:
        # dose = pCi_per_Bq * N_A * w * lambda * DF / M for 1 g pure Co60.
        df = nucleide.nuclei.dose_factor("Co60", "ingest", "EPA")
        assert df is not None
        lam = nucleide.nuclei.decay_constant("Co60")
        mass = nucleide.nuclei.atomic_mass("Co60")
        assert lam is not None and mass is not None
        expected = 27.027027 * 6.02214076e23 * 1.0 * lam * df / mass
        got = nucleide.material.dose_per_g({"Co60": 1.0}, "ingest", "EPA")
        assert got == pytest.approx(expected, rel=1e-12)

    def test_h3_air_matches_pyne_equation(self) -> None:
        df = nucleide.nuclei.dose_factor("H3", "air", "EPA")
        assert df is not None
        lam = nucleide.nuclei.decay_constant("H3")
        mass = nucleide.nuclei.atomic_mass("H3")
        assert lam is not None and mass is not None
        expected = 2.7027027e-11 * 6.02214076e23 * 1.0 * lam * df / mass
        got = nucleide.material.dose_per_g({"H3": 2.0}, "air", "EPA")
        # Per-gram: 2 g pure H3 gives the same per-gram dose.
        assert got == pytest.approx(expected, rel=1e-12)

    def test_mixture_weights_by_mass_fraction(self) -> None:
        # 1 g Co60 + 3 g Cs137: total per-gram dose is the mass-weighted mean.
        co = nucleide.material.dose_per_g({"Co60": 1.0}, "ingest", "EPA")
        cs = nucleide.material.dose_per_g({"Cs137": 1.0}, "ingest", "EPA")
        mix = nucleide.material.dose_per_g({"Co60": 1.0, "Cs137": 3.0}, "ingest", "EPA")
        assert mix == pytest.approx((co + 3.0 * cs) / 4.0, rel=1e-12)

    def test_missing_dose_errors(self) -> None:
        # I135 decays but has no dose row (dose covers 93 folded nuclides).
        with pytest.raises(ValueError, match="no dose factor"):
            nucleide.material.dose_per_g({"I135": 1.0}, "ingest", "EPA")
        with pytest.raises(ValueError, match="no dose factor"):
            nucleide.material.dose_per_g({"H3": 1.0}, "air", "GENII")
        with pytest.raises(ValueError):
            nucleide.material.dose_per_g({"Co60": 1.0}, "nope", "EPA")

    def test_k40_soil_sanity_band(self) -> None:
        # 1 g K-40 soil EPA per-gram dose is large but finite; check magnitude.
        got = nucleide.material.dose_per_g({"K40": 1.0}, "soil", "EPA")
        assert math.isfinite(got) and got > 0.0
