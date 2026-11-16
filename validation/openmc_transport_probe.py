#!/usr/bin/env python3
"""OpenMC 0.16 caller-data transport normalization probe (no dose claim).

Runs small fixed-source neutron and photon problems with an actual
NNDC/OpenMC HDF5 library and checks the activation-to-transport
normalization contract against transport:

- Water and iron spheres (10 cm radius) with a centered point source.
- Two runs with the same seed and particle count, source strengths 1 and 10.
- The fixed-source tally must scale by source strength, so the per-source
  conversion ``raw * rate / volume`` for the unit-strength run must agree with
  the rate-scaled conversion ``raw / volume`` for the strength-10 run.
- The iron sphere must produce a different flux than water, confirming that
  the loaded cross sections affect transport.
- When ``nucleide.openmc.normalize_flux`` is importable, the same conversion
  is checked through the package API.

Both particle types are exercised because the handoff source is photons and
the normalization convention is shared. This is a normalization/units and
data-loading check. It is not a transport benchmark, a shutdown-dose
calculation, or an independent physics oracle. The nuclear data are
supplied by the caller and never vendored or downloaded.

The data directory must contain ``H1.h5``, ``O16.h5`` and ``Fe56.h5`` plus
``photon/H.h5``, ``photon/O.h5`` and ``photon/Fe.h5``.

Usage:
    python validation/openmc_transport_probe.py --data-dir /path/to/hdf5
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
PARTICLES = 20_000
BATCHES = 5
SEED = 20261116
STRENGTHS = (1.0, 10.0)
NEUTRON_DATA = ("H1", "O16", "Fe56")
PHOTON_DATA = ("H", "O", "Fe")
NEUTRON_ENERGY_EV = 14.0e6
PHOTON_ENERGY_EV = 2.0e6
PARTICLES_TO_RUN = ("photon", "neutron")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_library(data_dir: Path, work: Path) -> Path:
    """Write a temporary cross_sections.xml with neutron and photon libraries."""
    entries: list[str] = []
    for name in NEUTRON_DATA:
        path = data_dir / f"{name}.h5"
        if not path.is_file():
            raise FileNotFoundError(f"missing neutron data file: {name}.h5")
        entries.append(f'  <library materials="{name}" path="{path}" type="neutron"/>')
    for name in PHOTON_DATA:
        path = data_dir / "photon" / f"{name}.h5"
        if not path.is_file():
            raise FileNotFoundError(f"missing photon data file: photon/{name}.h5")
        entries.append(f'  <library materials="{name}" path="{path}" type="photon"/>')
    xml = '<?xml version="1.0"?>\n<cross_sections>\n' + "\n".join(entries) + "\n</cross_sections>\n"
    cross_sections = work / "cross_sections.xml"
    cross_sections.write_text(xml, encoding="utf-8")
    return cross_sections


def run_case(
    openmc: Any,
    work: Path,
    material_name: str,
    particle: str,
    strength: float,
) -> dict[str, Any]:
    """Run one fixed-source sphere case and return the cell flux tally."""
    if particle not in PARTICLES_TO_RUN:
        raise ValueError(f"unknown particle {particle}")
    energy = PHOTON_ENERGY_EV if particle == "photon" else NEUTRON_ENERGY_EV
    case_dir = work / f"{particle}-{material_name}-s{strength:g}"
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
        energy=openmc.stats.Discrete([energy], [1.0]),
        strength=strength,
        particle=particle,
    )
    settings = openmc.Settings(
        particles=PARTICLES,
        batches=BATCHES,
        run_mode="fixed source",
        seed=SEED,
        output={"summary": False, "tallies": False},
    )
    settings.source = [source]
    if particle == "photon":
        settings.photon_transport = True
    tally = openmc.Tally(name=f"{particle}-{material_name} flux")
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
            f"{particle} {material_name} strength {strength}: nonpositive/nonfinite tally {mean}"
        )
    if not math.isfinite(std_dev):
        std_dev = None
    return {
        "particle": particle,
        "material": material_name,
        "source_strength": strength,
        "source_energy_ev": energy,
        "mean_volume_integrated_track_length": mean,
        "std_dev": std_dev,
        "relative_error": (std_dev / mean if std_dev is not None else None),
    }


def package_check(unit: dict[str, Any], ten: dict[str, Any], volume: float) -> dict[str, Any]:
    """Check the same conversion through the optional Nucleide package."""
    try:
        from nucleide import openmc as omc
    except Exception as error:  # pragma: no cover - optional environment
        return {"status": "SKIP", "reason": f"nucleide import failed: {error}"}
    if not hasattr(omc, "normalize_flux"):
        return {"status": "SKIP", "reason": "nucleide.openmc.normalize_flux is absent"}
    arrays_unit = {
        "flux": [[unit["mean_volume_integrated_track_length"]]],
        "rel_err": [[0.0]],
        "normalization": "unknown",
        "units": "unknown",
    }
    arrays_ten = {
        "flux": [[ten["mean_volume_integrated_track_length"]]],
        "rel_err": [[0.0]],
        "normalization": "unknown",
        "units": "unknown",
    }
    per_source = omc.normalize_flux(
        arrays_unit,
        10.0,
        [volume],
        input_normalization="per-source-volume-integrated",
    )
    rate_scaled = omc.normalize_flux(
        arrays_ten,
        10.0,
        [volume],
        input_normalization="rate-scaled-volume-integrated",
    )
    unit_flux = float(per_source["flux"][0][0])
    ten_flux = float(rate_scaled["flux"][0][0])
    return {
        "status": "PASS"
        if math.isclose(unit_flux, ten_flux, rel_tol=1e-12, abs_tol=0.0)
        else "FAIL",
        "per_source_flux": unit_flux,
        "rate_scaled_flux": ten_flux,
        "relative_difference": abs(unit_flux - ten_flux) / ten_flux,
        "units": per_source["units"],
    }


def conversion_checks(
    particle: str,
    unit: dict[str, Any],
    ten: dict[str, Any],
    iron: dict[str, Any],
    volume: float,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Build the scaling, conversion and data-loading checks for one particle."""
    raw_unit = unit["mean_volume_integrated_track_length"]
    raw_ten = ten["mean_volume_integrated_track_length"]
    ratio = raw_ten / raw_unit
    per_source_flux = raw_unit * 10.0 / volume
    rate_scaled_flux = raw_ten / volume
    iron_flux = iron["mean_volume_integrated_track_length"] / volume
    checks = [
        {
            "name": f"{particle}_source_strength_scales_tally",
            "pass": math.isclose(ratio, 10.0, rel_tol=1e-12, abs_tol=0.0),
            "ratio": ratio,
        },
        {
            "name": f"{particle}_per_source_and_rate_scaled_flux_agree",
            "pass": math.isclose(per_source_flux, rate_scaled_flux, rel_tol=1e-12, abs_tol=0.0),
            "per_source_flux": per_source_flux,
            "rate_scaled_flux": rate_scaled_flux,
            "relative_difference": abs(per_source_flux - rate_scaled_flux) / rate_scaled_flux,
        },
        {
            "name": f"{particle}_cross_sections_change_transport",
            "pass": not math.isclose(rate_scaled_flux, iron_flux, rel_tol=1e-3, abs_tol=0.0),
            "water_flux": rate_scaled_flux,
            "iron_flux": iron_flux,
        },
    ]
    return checks, {
        "water_unit": unit,
        "water_ten": ten,
        "iron_ten": iron,
        "per_source_flux": per_source_flux,
        "rate_scaled_flux": rate_scaled_flux,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        required=True,
        help="directory containing H1.h5, O16.h5, Fe56.h5 and photon/{H,O,Fe}.h5",
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
        print(f"data directory not found: {data_dir}", file=sys.stderr)
        print(json.dumps({"status": "SKIP", "reason": "data directory not found"}))
        return 2

    volume = 4.0 / 3.0 * math.pi * RADIUS_CM**3
    with tempfile.TemporaryDirectory(prefix="openmc-transport-") as tmp:
        work = Path(tmp)
        try:
            cross_sections = write_library(data_dir, work)
        except Exception as error:
            print(
                json.dumps(
                    {
                        "status": "SKIP",
                        "reason": f"{type(error).__name__}: {error}",
                    }
                )
            )
            return 2
        previous = os.environ.get("OPENMC_CROSS_SECTIONS")
        os.environ["OPENMC_CROSS_SECTIONS"] = str(cross_sections)
        try:
            results: dict[str, dict[str, Any]] = {}
            for particle in PARTICLES_TO_RUN:
                unit = run_case(openmc, work, "water", particle, STRENGTHS[0])
                ten = run_case(openmc, work, "water", particle, STRENGTHS[1])
                iron = run_case(openmc, work, "iron", particle, STRENGTHS[1])
                results[particle] = {"unit": unit, "ten": ten, "iron": iron}
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

    checks: list[dict[str, Any]] = []
    particle_blocks: dict[str, Any] = {}
    for particle in PARTICLES_TO_RUN:
        block = results[particle]
        particle_checks, summary = conversion_checks(
            particle, block["unit"], block["ten"], block["iron"], volume
        )
        checks.extend(particle_checks)
        particle_blocks[particle] = summary
    api = package_check(results["photon"]["unit"], results["photon"]["ten"], volume)
    api_ok = api.get("status") in ("PASS", "SKIP")
    checks.append(
        {
            "name": "nucleide_package_normalize_flux",
            "pass": api_ok,
            "status": api.get("status"),
            "detail": api,
        }
    )
    passed = all(check["pass"] for check in checks)
    summary = {
        "probe": "openmc_transport_probe",
        "status": "PASS" if passed else "FAIL",
        "openmc_version": openmc.__version__,
        "data_dir": "<caller-supplied>",
        "required_data": {
            "neutron": {
                name: {
                    "relative_path": f"{name}.h5",
                    "sha256": sha256_file(data_dir / f"{name}.h5"),
                    "size_bytes": (data_dir / f"{name}.h5").stat().st_size,
                }
                for name in NEUTRON_DATA
            },
            "photon": {
                name: {
                    "relative_path": f"photon/{name}.h5",
                    "sha256": sha256_file(data_dir / "photon" / f"{name}.h5"),
                    "size_bytes": (data_dir / "photon" / f"{name}.h5").stat().st_size,
                }
                for name in PHOTON_DATA
            },
        },
        "case": {
            "geometry": f"vacuum sphere r={RADIUS_CM} cm",
            "materials": [
                "water (H1/O16, 1.0 g/cm3)",
                "iron (Fe56, 7.874 g/cm3)",
            ],
            "sources": [
                f"photon {PHOTON_ENERGY_EV / 1e6:g} MeV point, centered",
                f"neutron {NEUTRON_ENERGY_EV / 1e6:g} MeV point, centered",
            ],
            "histories": f"{PARTICLES} particles x {BATCHES} batches, seed {SEED}",
            "volume_cm3": volume,
        },
        "particles": particle_blocks,
        "checks": checks,
        "limitations": [
            "This is a normalization/units and data-loading check, not an "
            "independent transport or shutdown-dose validation.",
            "The sphere and source are synthetic; no evaluated activation data are used.",
            "The HDF5 files are caller-supplied and never vendored or downloaded.",
        ],
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
