#!/usr/bin/env python3
"""OpenMC 0.16 real-data normalization probe (no dose or transport claim).

Runs a small fixed-source neutron problem with an actual NNDC/OpenMC HDF5
cross-section library and checks the activation-to-transport normalization
contract against the transport result:

- A centered 14 MeV point source in a 10 cm radius water sphere.
- Two runs with the same seed and particle count, source strengths 1 and 10.
- The fixed-source tally must scale exactly by the source strength, so the
  per-source conversion ``raw * rate / volume`` for the unit-strength run
  must agree with the rate-scaled conversion ``raw / volume`` for the
  strength-10 run.
- A second material (iron) must produce a different tally, confirming that
  the HDF5 cross sections were actually loaded and used.
- When ``nucleide.openmc.normalize_flux`` is importable, the same conversion
  is checked through the package API.

This is a normalization/units and data-loading check. It is not a transport
benchmark, a shutdown-dose calculation, or an independent physics oracle.
The real nuclear data are supplied by the caller and never vendored or
downloaded by this script.

Usage:
    python validation/openmc_real_data_probe.py --data-dir /path/to/hdf5

The data directory must contain ``H1.h5``, ``O16.h5`` and ``Fe56.h5``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import tempfile
import traceback
from pathlib import Path
from typing import Any

EXPECTED_OPENMC = "0.16.0"
RADIUS_CM = 10.0
SOURCE_ENERGY_EV = 14.0e6
PARTICLES = 20_000
BATCHES = 5
SEED = 20261116
STRENGTHS = (1.0, 10.0)
REQUIRED_DATA = ("H1", "O16", "Fe56")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_library(openmc: Any, data_dir: Path, work: Path) -> Path:
    """Register the required HDF5 files and write a temporary cross_sections.xml."""
    missing = [name for name in REQUIRED_DATA if not (data_dir / f"{name}.h5").is_file()]
    if missing:
        raise FileNotFoundError(
            f"data directory {data_dir} is missing required files: "
            + ", ".join(f"{name}.h5" for name in missing)
        )
    library = openmc.data.DataLibrary()
    for name in REQUIRED_DATA:
        library.register_file(str(data_dir / f"{name}.h5"))
    path = work / "cross_sections.xml"
    library.export_to_xml(str(path))
    return path


def run_case(openmc: Any, work: Path, material_name: str, strength: float) -> dict[str, Any]:
    """Run one fixed-source sphere case and return the cell flux tally."""
    case_dir = work / f"{material_name}-s{strength:g}"
    case_dir.mkdir()
    sphere = openmc.Sphere(r=RADIUS_CM, boundary_type="vacuum")
    cell = openmc.Cell(region=-sphere)
    if material_name == "water":
        material = openmc.Material(name="water")
        material.add_nuclide("H1", 2.0)
        material.add_nuclide("O16", 1.0)
        material.set_density("g/cm3", 1.0)
    elif material_name == "iron":
        material = openmc.Material(name="iron")
        material.add_nuclide("Fe56", 1.0)
        material.set_density("g/cm3", 7.874)
    else:
        raise ValueError(f"unknown material {material_name}")
    cell.fill = material
    geometry = openmc.Geometry([cell])
    materials = openmc.Materials([material])
    source = openmc.IndependentSource(
        space=openmc.stats.Point((0.0, 0.0, 0.0)),
        energy=openmc.stats.Discrete([SOURCE_ENERGY_EV], [1.0]),
        strength=strength,
        particle="neutron",
    )
    settings = openmc.Settings(
        particles=PARTICLES,
        batches=BATCHES,
        run_mode="fixed source",
        seed=SEED,
        output={"summary": False, "tallies": False},
    )
    settings.source = [source]
    tally = openmc.Tally(name=f"{material_name} flux")
    tally.filters = [openmc.CellFilter(cell)]
    tally.scores = ["flux"]
    model = openmc.Model(
        geometry=geometry,
        materials=materials,
        settings=settings,
        tallies=openmc.Tallies([tally]),
    )
    model.export_to_xml(case_dir)
    openmc.run(output=False, cwd=str(case_dir))
    statepoint = openmc.StatePoint(str(case_dir / f"statepoint.{BATCHES}.h5"))
    try:
        data = statepoint.tallies[tally.id]
        mean = float(data.mean.ravel()[0])
        std_dev = float(data.std_dev.ravel()[0])
    finally:
        close = getattr(statepoint, "close", None)
        if callable(close):
            close()
    if not math.isfinite(mean) or mean <= 0.0:
        raise RuntimeError(
            f"{material_name} strength {strength}: nonpositive/nonfinite tally {mean}"
        )
    if not math.isfinite(std_dev):
        std_dev = None
    return {
        "material": material_name,
        "source_strength": strength,
        "mean_volume_integrated_track_length": mean,
        "std_dev": std_dev,
        "relative_error": (std_dev / mean if std_dev is not None else None),
    }


def package_check(
    water_unit: dict[str, Any], water_ten: dict[str, Any], volume: float
) -> dict[str, Any]:
    """Check the same conversion through the optional Nucleide package."""
    try:
        from nucleide import openmc as omc
    except Exception as error:  # pragma: no cover - optional environment
        return {"status": "SKIP", "reason": f"nucleide import failed: {error}"}
    if not hasattr(omc, "normalize_flux"):
        return {"status": "SKIP", "reason": "nucleide.openmc.normalize_flux is absent"}
    arrays_unit = {
        "flux": [[water_unit["mean_volume_integrated_track_length"]]],
        "rel_err": [[0.0]],
        "normalization": "unknown",
        "units": "unknown",
    }
    arrays_ten = {
        "flux": [[water_ten["mean_volume_integrated_track_length"]]],
        "rel_err": [[0.0]],
        "normalization": "unknown",
        "units": "unknown",
    }
    unit = omc.normalize_flux(
        arrays_unit,
        10.0,
        [volume],
        input_normalization="per-source-volume-integrated",
    )
    ten = omc.normalize_flux(
        arrays_ten,
        10.0,
        [volume],
        input_normalization="rate-scaled-volume-integrated",
    )
    unit_flux = float(unit["flux"][0][0])
    ten_flux = float(ten["flux"][0][0])
    return {
        "status": "PASS"
        if math.isclose(unit_flux, ten_flux, rel_tol=1e-12, abs_tol=0.0)
        else "FAIL",
        "per_source_flux": unit_flux,
        "rate_scaled_flux": ten_flux,
        "relative_difference": abs(unit_flux - ten_flux) / ten_flux,
        "units": unit["units"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        required=True,
        help="directory containing H1.h5, O16.h5 and Fe56.h5",
    )
    args = parser.parse_args()

    try:
        import openmc
    except Exception as error:
        print(json.dumps({"status": "SKIP", "reason": f"OpenMC import failed: {error}"}))
        return 2
    if openmc.__version__ != EXPECTED_OPENMC:
        print(
            json.dumps(
                {
                    "status": "SKIP",
                    "reason": f"OpenMC {openmc.__version__} != {EXPECTED_OPENMC}",
                }
            )
        )
        return 2
    data_dir = args.data_dir.expanduser().resolve()
    if not data_dir.is_dir():
        print(json.dumps({"status": "SKIP", "reason": f"data directory not found: {data_dir}"}))
        return 2

    volume = 4.0 / 3.0 * math.pi * RADIUS_CM**3
    with tempfile.TemporaryDirectory(prefix="openmc-real-data-") as tmp:
        work = Path(tmp)
        cross_sections = write_library(openmc, data_dir, work)
        previous = os.environ.get("OPENMC_CROSS_SECTIONS")
        os.environ["OPENMC_CROSS_SECTIONS"] = str(cross_sections)
        try:
            water_unit = run_case(openmc, work, "water", STRENGTHS[0])
            water_ten = run_case(openmc, work, "water", STRENGTHS[1])
            iron_ten = run_case(openmc, work, "iron", STRENGTHS[1])
        except Exception as error:
            print(
                json.dumps(
                    {
                        "status": "FAIL",
                        "error": f"{type(error).__name__}: {error}",
                        "traceback": traceback.format_exc(),
                    },
                    indent=2,
                )
            )
            return 1
        finally:
            if previous is None:
                os.environ.pop("OPENMC_CROSS_SECTIONS", None)
            else:
                os.environ["OPENMC_CROSS_SECTIONS"] = previous

    raw_unit = water_unit["mean_volume_integrated_track_length"]
    raw_ten = water_ten["mean_volume_integrated_track_length"]
    ratio = raw_ten / raw_unit
    per_source_flux = raw_unit * 10.0 / volume
    rate_scaled_flux = raw_ten / volume
    iron_flux = iron_ten["mean_volume_integrated_track_length"] / volume
    water_flux = rate_scaled_flux
    ratio_ok = math.isclose(ratio, 10.0, rel_tol=1e-12, abs_tol=0.0)
    conversion_ok = math.isclose(per_source_flux, rate_scaled_flux, rel_tol=1e-12, abs_tol=0.0)
    data_applied = not math.isclose(water_flux, iron_flux, rel_tol=1e-3, abs_tol=0.0)
    api = package_check(water_unit, water_ten, volume)
    api_ok = api.get("status") in ("PASS", "SKIP")
    checks = [
        {"name": "source_strength_scales_tally", "pass": ratio_ok, "ratio": ratio},
        {
            "name": "per_source_and_rate_scaled_flux_agree",
            "pass": conversion_ok,
            "per_source_flux": per_source_flux,
            "rate_scaled_flux": rate_scaled_flux,
            "relative_difference": abs(per_source_flux - rate_scaled_flux) / rate_scaled_flux,
        },
        {
            "name": "real_cross_sections_change_transport",
            "pass": data_applied,
            "water_flux": water_flux,
            "iron_flux": iron_flux,
        },
        {
            "name": "nucleide_package_normalize_flux",
            "pass": api_ok,
            "status": api.get("status"),
            "detail": api,
        },
    ]
    passed = all(check["pass"] for check in checks)
    summary = {
        "probe": "openmc_real_data_probe",
        "status": "PASS" if passed else "FAIL",
        "openmc_version": openmc.__version__,
        "data_dir": str(data_dir),
        "required_data": {
            name: {
                "path": str(data_dir / f"{name}.h5"),
                "sha256": sha256_file(data_dir / f"{name}.h5"),
                "size_bytes": (data_dir / f"{name}.h5").stat().st_size,
            }
            for name in REQUIRED_DATA
        },
        "case": {
            "geometry": f"vacuum sphere r={RADIUS_CM} cm",
            "materials": ["water (H1/O16, 1.0 g/cm3)", "iron (Fe56, 7.874 g/cm3)"],
            "source": (
                f"centered 14 MeV point source, {PARTICLES} particles x "
                f"{BATCHES} batches, seed {SEED}"
            ),
            "volume_cm3": volume,
        },
        "water_unit": water_unit,
        "water_ten": water_ten,
        "iron_ten": iron_ten,
        "checks": checks,
        "limitations": [
            "This is a normalization/units and data-loading check, not an independent "
            "transport or shutdown-dose validation.",
            "The sphere and source are synthetic; no evaluated activation or photon "
            "source data are used.",
            "The real HDF5 files are caller-supplied and never vendored or downloaded.",
        ],
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
