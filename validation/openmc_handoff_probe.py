#!/usr/bin/env python3
"""OpenMC 0.16.0 activation-to-transport handoff probe (no transport run).

This probe answers seven questions with the OpenMC package installed in the
validation image:

1. Is the exact expected OpenMC version installed?
2. How does OpenMC 0.16.0 normalize regular tally results, and what role does
   ``SourceBase.strength`` / ``MeshSource.strength`` play? A vacuum-sphere
   fixed-source run at strengths 1 and 10 checks the observed scale directly
   without interacting material or evaluated data.
3. Does a C-order ``[nx][ny][nz]`` source array reach the expected
   ``(i, j, k)`` mesh element when passed to ``openmc.MeshSource``?
4. Does an OpenMC histogram distribution interpret ``p`` as a differential
   density (per unit energy), not as an integrated bin strength?
5. Do mesh id, source count, strengths, energy edges and source ordering
   survive ``Settings.export_to_xml`` / ``Settings.from_xml``?
6. Does real OpenMC source sampling reproduce the exact joint voxel/group
   probabilities, using a fixed seed and a predeclared tolerance?

The script runs no interacting transport calculation. Source sampling is
isolated in a child process and uses ``openmc.Model.sample_external_source``,
which starts OpenMC in volume-calculation mode. A second child runs a vacuum
sphere with a centered point source and a cell flux tally at two source
strengths; this checks the volume-integrated, source-rate-scaled normalization
without cross sections or material interactions. OpenMC still requires a
cross-sections XML file to initialize; the children write a placeholder
containing one nonexistent library entry. No evaluated data are read or
downloaded, and only the vacuum track-length tally is produced. The
placeholder is an initialization shim, not physics input.

Output is a single JSON document on stdout.  Human-readable failures and
OpenMC warnings go to stderr.  Exit status is 0 only when every check passes.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import inspect
import json
import math
import os
import re
import subprocess
import sys
import tempfile
import traceback
import warnings
from pathlib import Path
from typing import Any

EXPECTED_OPENMC = "0.16.0"
PROBE_NAME = "openmc_handoff_probe"

# Fixed, predeclared sampling contract.
PRN_SEED = 20261116
JOINT_SAMPLES = 200_000
JOINT_TOLERANCE = 0.01
ORDER_SAMPLES = 120_000
ZERO_SAMPLES = 20_000

# Void fixed-source normalization check: an isotropic point source at the
# center of a vacuum sphere has exactly one track of length R per source
# particle, so the volume-integrated flux tally mean is R * total strength.
VOID_RADIUS = 10.0
VOID_STRENGTHS = (1.0, 10.0)
VOID_PARTICLES = 2_000
VOID_SEED = 20261116

# Two-voxel, two-group scenario from the transport handoff plan.
JOINT_EDGES = (1.0, 2.0, 10.0)
JOINT_VOXEL_STRENGTHS = (100.0, 300.0)
JOINT_INTEGRATED = ((1.0, 0.0), (0.0, 1.0))
JOINT_EXACT = ((0.25, 0.0), (0.0, 0.75))

# 2 x 3 x 2 ordering mesh.  Every cell gets a unique strength and a unique
# non-overlapping energy band so sampled energies reveal the selected cell.
ORDER_DIMENSION = (2, 3, 2)
ORDER_LOWER = (0.0, 0.0, 0.0)
ORDER_UPPER = (2.0, 3.0, 2.0)
ORDER_BAND_WIDTH = 10.0
ORDER_BAND_LOW = 1.0
ORDER_BAND_SPAN = 4.0

# Irregular three-bin histogram with a zero-strength middle bin.
HIST_EDGES = (1.0, 2.0, 5.0, 10.0)
HIST_INTEGRATED = (0.25, 0.0, 0.75)


def histogram_density(edges: tuple[float, ...], integrated: tuple[float, ...]) -> list[float]:
    """Convert integrated bin strengths to differential densities per unit x."""
    if len(integrated) != len(edges) - 1:
        raise ValueError("integrated strengths must match the number of bins")
    return [
        strength / (hi - lo)
        for strength, lo, hi in zip(integrated, edges[:-1], edges[1:], strict=True)
    ]


HIST_DENSITIES = tuple(histogram_density(HIST_EDGES, HIST_INTEGRATED))


def sha256_file(path: Path) -> str:
    """Hash a file without loading the whole file at once."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def find_openmc_source_dir() -> tuple[Path | None, str]:
    """Find the OpenMC source tree that the installed package was built from."""
    candidates: list[tuple[Path, str]] = []
    env_dir = os.environ.get("OPENMC_SOURCE_DIR")
    if env_dir:
        candidates.append((Path(env_dir), "OPENMC_SOURCE_DIR"))
    try:
        distribution = importlib.metadata.distribution("openmc")
        direct_url = json.loads(distribution.read_text("direct_url.json"))
        url = direct_url.get("url", "")
        if url.startswith("file://"):
            candidates.append((Path(url.removeprefix("file://")), "direct_url.json"))
    except Exception:
        pass
    candidates.extend(
        [
            (Path("/opt/openmc"), "default /opt/openmc"),
            (Path("/usr/local/src/openmc"), "default /usr/local/src/openmc"),
        ]
    )
    for candidate, origin in candidates:
        if (candidate / "src/tallies/tally.cpp").is_file():
            return candidate.resolve(), origin
    return None, "not found"


def find_pattern(text: str, pattern: str, context: int = 2) -> dict[str, Any]:
    """Locate a regex match and return its line number and a small context window."""
    match = re.search(pattern, text, flags=re.MULTILINE)
    if match is None:
        return {"found": False, "pattern": pattern}
    start_line = text.count("\n", 0, match.start()) + 1
    lines = text.splitlines()
    lo = max(1, start_line - context)
    hi = min(len(lines), start_line + context)
    return {
        "found": True,
        "pattern": pattern,
        "line": start_line,
        "match": match.group(0).strip(),
        "context": "\n".join(lines[lo - 1 : hi]),
    }


def collect_source_provenance(source_dir: Path | None) -> dict[str, Any]:
    """Inspect the installed OpenMC C++/docs source for normalization and sampling rules."""
    if source_dir is None:
        return {"available": False, "reason": "OpenMC source tree not found"}

    files = {
        "tally_cpp": source_dir / "src/tallies/tally.cpp",
        "settings_cpp": source_dir / "src/settings.cpp",
        "distribution_cpp": source_dir / "src/distribution.cpp",
        "source_cpp": source_dir / "src/source.cpp",
        "simulation_cpp": source_dir / "src/simulation.cpp",
        "tallies_rst": source_dir / "docs/source/usersguide/tallies.rst",
        "settings_rst": source_dir / "docs/source/usersguide/settings.rst",
        "io_settings_rst": source_dir / "docs/source/io_formats/settings.rst",
    }
    missing = [str(path) for path in files.values() if not path.is_file()]
    if missing:
        return {"available": False, "reason": "missing source files", "missing": missing}

    patterns: dict[str, list[tuple[str, str]]] = {
        "tally_cpp": [
            ("total_source", r"double total_source = 0\.0;"),
            (
                "contributing_particles",
                r"double contributing_particles = settings::reduce_tallies",
            ),
            (
                "regular_norm",
                r"double norm =\s*total_source / "
                r"\(contributing_particles \* settings::gen_per_batch\);",
            ),
            ("regular_accumulate", r"double val = results_\(i, j, TallyResult::VALUE\) \* norm;"),
            (
                "global_total_weight",
                r"double val = gt\(i, TallyResult::VALUE\) / simulation::total_weight;",
            ),
        ],
        "settings_cpp": [
            ("source_strengths_push", r"source_strengths\.push_back\(s->strength\(\)\);"),
            (
                "external_sources_probability",
                r"model::external_sources_probability\.assign\(source_strengths\);",
            ),
        ],
        "distribution_cpp": [
            (
                "discrete_integral",
                r"integral_ = std::accumulate\(prob_\.begin\(\), prob_\.end\(\), 0\.0\);",
            ),
            ("discrete_normalize", r"p_i /= integral_;"),
        ],
        "source_cpp": [
            (
                "mesh_strengths",
                r"strengths\.push_back\(sources_\.back\(\)->strength\(\)\);",
            ),
            (
                "mesh_spatial",
                r"space_ = std::make_unique<MeshSpatial>\(mesh_idx, strengths\);",
            ),
            (
                "mesh_sample_element",
                r"int32_t element = space_->sample_element_index\(seed\);",
            ),
        ],
        "simulation_cpp": [
            ("total_weight_reset", r"simulation::total_weight = 0\.0;"),
            ("total_weight_accum", r"simulation::total_weight \+= p\.wgt\(\);"),
        ],
        "tallies_rst": [
            (
                "per_source_particle",
                r"all tally scores are normalized per\s+source particle simulated",
            ),
            (
                "strength_as_unit_conversion",
                r"the tally results just need to be multipled by 10\\?\s*"
                r":sup:`4`\. This can either be done manually or using the\s+"
                r":attr:`openmc\.SourceBase\.strength` attribute\.",
            ),
        ],
        "settings_rst": [
            (
                "relative_strength",
                r"attribute that indicates the relative strength of a source distribution",
            ),
            (
                "uniform_sampling_weight",
                r"applying the source strength as a weight on the sampled source particles",
            ),
        ],
        "io_settings_rst": [
            (
                "mesh_element_relative_strength",
                r"relative source strength of each mesh element",
            ),
        ],
    }

    text_by_name: dict[str, str] = {}
    provenance: dict[str, Any] = {"available": True, "source_dir": str(source_dir)}
    for name, path in files.items():
        text = path.read_text(encoding="utf-8", errors="replace")
        text_by_name[name] = text
        provenance[name] = {
            "path": str(path),
            "sha256": sha256_file(path),
            "patterns": {label: find_pattern(text, pattern) for label, pattern in patterns[name]},
        }

    # Compute the exact regular-tally normalization context without relying on
    # the literal string "total_weight" appearing anywhere in the file.
    tally_text = text_by_name["tally_cpp"]
    regular_norm_match = re.search(
        r"double norm =\s*total_source / "
        r"\(contributing_particles \* settings::gen_per_batch\);",
        tally_text,
    )
    regular_norm_window = ""
    if regular_norm_match is not None:
        start = max(0, regular_norm_match.start() - 800)
        end = min(len(tally_text), regular_norm_match.end() + 800)
        regular_norm_window = tally_text[start:end]
    provenance["derived"] = {
        "regular_norm_uses_simulation_total_weight": "simulation::total_weight"
        in regular_norm_window,
        "regular_norm_uses_n_particles": "contributing_particles" in regular_norm_window
        and "settings::gen_per_batch" in regular_norm_window
        and "settings::n_particles"
        in provenance["tally_cpp"]["patterns"]["contributing_particles"].get("context", ""),
        "regular_norm_scales_by_total_source": "total_source" in regular_norm_window
        and "external_sources_probability" in text_by_name["tally_cpp"],
        "mesh_source_strength_is_sum": "sum(s.strength for s in self.sources)"
        in inspect.getsource(__import__("openmc").MeshSource.strength.fget),
    }
    return provenance


def collect_python_provenance(openmc: Any) -> dict[str, Any]:
    """Record the exact installed Python API lines that define the bridge contract."""
    targets = {
        "mesh_source_sources_setter": openmc.MeshSource.sources.fset,
        "mesh_source_strength_getter": openmc.MeshSource.strength.fget,
        "mesh_source_set_total_strength": openmc.MeshSource.set_total_strength,
        "source_base_strength_setter": openmc.SourceBase.strength.fset,
        "regular_mesh_indices": openmc.RegularMesh.indices.fget,
        "regular_mesh_get_indices_at_coords": openmc.RegularMesh.get_indices_at_coords,
        "tabular_cdf": openmc.stats.Tabular.cdf,
        "tabular_integral": openmc.stats.Tabular.integral,
    }
    result: dict[str, Any] = {}
    for name, target in targets.items():
        try:
            lines, start = inspect.getsourcelines(target)
            result[name] = {
                "file": inspect.getsourcefile(target),
                "line": start,
                "snippet": "".join(lines),
            }
        except (OSError, TypeError) as error:
            result[name] = {"error": str(error)}
    result["histogram_class_available"] = hasattr(openmc.stats, "Histogram")
    result["mesh_source_doc"] = inspect.getdoc(openmc.MeshSource)
    return result


def build_ordering_mesh(openmc: Any, mesh_id: int) -> tuple[Any, list[list[list[Any]]], Any]:
    """Build a 2 x 3 x 2 mesh with a unique strength and energy band per cell."""
    nx, ny, nz = ORDER_DIMENSION
    mesh = openmc.RegularMesh()
    mesh.id = mesh_id
    mesh.dimension = ORDER_DIMENSION
    mesh.lower_left = ORDER_LOWER
    mesh.upper_right = ORDER_UPPER

    sources_c: list[list[list[Any]]] = []
    for i in range(nx):
        plane: list[list[Any]] = []
        for j in range(ny):
            column: list[Any] = []
            for k in range(nz):
                c_index = (i * ny + j) * nz + k
                low = ORDER_BAND_LOW + ORDER_BAND_WIDTH * c_index
                energy = openmc.stats.Uniform(low, low + ORDER_BAND_SPAN)
                source = openmc.IndependentSource(
                    energy=energy,
                    strength=float(c_index + 1),
                    particle="photon",
                )
                column.append(source)
            plane.append(column)
        sources_c.append(plane)

    mesh_source = openmc.MeshSource(mesh, sources_c)
    return mesh, sources_c, mesh_source


def build_joint_mesh(openmc: Any, mesh_id: int) -> tuple[Any, Any]:
    """Build the two-voxel, two-group joint scenario with zero cross bins."""
    mesh = openmc.RegularMesh()
    mesh.id = mesh_id
    mesh.dimension = (2, 1, 1)
    mesh.lower_left = (0.0, 0.0, 0.0)
    mesh.upper_right = (2.0, 1.0, 1.0)

    sources = []
    for voxel, strength in enumerate(JOINT_VOXEL_STRENGTHS):
        densities = histogram_density(JOINT_EDGES, JOINT_INTEGRATED[voxel])
        energy = openmc.stats.Tabular(
            list(JOINT_EDGES),
            densities,
            interpolation="histogram",
        )
        sources.append(
            openmc.IndependentSource(energy=energy, strength=strength, particle="photon")
        )
    return mesh, openmc.MeshSource(mesh, sources)


def build_histogram_case(openmc: Any) -> dict[str, Any]:
    """Verify histogram density semantics, a zero bin, and a zero-strength voxel."""
    distribution = openmc.stats.Tabular(
        list(HIST_EDGES),
        list(HIST_DENSITIES),
        interpolation="histogram",
    )
    wrong_control = openmc.stats.Tabular(
        list(HIST_EDGES),
        list(HIST_INTEGRATED),
        interpolation="histogram",
    )
    all_zero = openmc.stats.Tabular(
        list(HIST_EDGES),
        [0.0, 0.0, 0.0],
        interpolation="histogram",
    )

    midpoints = [0.5 * (lo + hi) for lo, hi in zip(HIST_EDGES[:-1], HIST_EDGES[1:], strict=True)]
    widths = [hi - lo for lo, hi in zip(HIST_EDGES[:-1], HIST_EDGES[1:], strict=True)]
    cdf = distribution.cdf().tolist()
    integral = float(distribution.integral())
    zero_density = float(distribution.p[1])

    zero_strength_source = openmc.IndependentSource(
        energy=all_zero, strength=0.0, particle="photon"
    )
    mesh = openmc.RegularMesh()
    mesh.id = 19003
    mesh.dimension = (2, 1, 1)
    mesh.lower_left = (0.0, 0.0, 0.0)
    mesh.upper_right = (2.0, 1.0, 1.0)
    nonzero_source = openmc.IndependentSource(energy=distribution, strength=1.0, particle="photon")
    zero_voxel_mesh_source = openmc.MeshSource(mesh, [nonzero_source, zero_strength_source])

    return {
        "histogram_class_available": hasattr(openmc.stats, "Histogram"),
        "supported_equivalent": "openmc.stats.Tabular(..., interpolation='histogram')",
        "edges_eV": list(HIST_EDGES),
        "bin_widths_eV": widths,
        "integrated_strengths": list(HIST_INTEGRATED),
        "expected_density_pdf_per_eV": list(HIST_DENSITIES),
        "readback_x": distribution.x.tolist(),
        "readback_p": distribution.p.tolist(),
        "interpolation": distribution.interpolation,
        "cdf": cdf,
        "integral": integral,
        "evaluate_at_bin_midpoints": [float(distribution.evaluate(x)) for x in midpoints],
        "zero_bin": {
            "index": 1,
            "width_eV": widths[1],
            "integrated_strength": HIST_INTEGRATED[1],
            "density": zero_density,
            "evaluate_midpoint": float(distribution.evaluate(midpoints[1])),
            "cdf_flat": cdf[1] == cdf[2],
        },
        "wrong_control_integral": float(wrong_control.integral()),
        "wrong_control_note": (
            "Passing integrated strengths as histogram densities gives integral "
            f"{float(wrong_control.integral()):.6g}, not 1"
        ),
        "all_zero_distribution": {
            "integral": float(all_zero.integral()),
            "cdf": all_zero.cdf().tolist(),
            "pdf_valid": float(all_zero.integral()) > 0.0,
            "note": "A zero distribution is data, not a samplable PDF",
        },
        "zero_strength_voxel": {
            "source_strength": float(zero_strength_source.strength),
            "mesh_source_total_strength": float(zero_voxel_mesh_source.strength),
            "element_strengths": [
                float(source.strength) for source in zero_voxel_mesh_source.sources
            ],
            "relative_probabilities": [
                float(source.strength) / float(zero_voxel_mesh_source.strength)
                for source in zero_voxel_mesh_source.sources
            ],
        },
    }


def check_ordering_objects(
    openmc: Any, mesh: Any, sources_c: list[list[list[Any]]], mesh_source: Any
) -> dict[str, Any]:
    """Verify object-level Fortran-order flattening and 0-based coordinate lookup."""
    nx, ny, nz = ORDER_DIMENSION
    flat_sources = mesh_source.sources
    mapping_ok = True
    strengths_f_order: list[float] = [0.0] * int(mesh_source.sources.size)
    bands_f_order: list[float] = [0.0] * int(mesh_source.sources.size)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                f_index = i + nx * j + nx * ny * k
                source = sources_c[i][j][k]
                if flat_sources[f_index] is not source:
                    mapping_ok = False
                strengths_f_order[f_index] = float(flat_sources[f_index].strength)
                bands_f_order[f_index] = float(flat_sources[f_index].energy.a)

    expected_indices = [
        (i + 1, j + 1, k + 1) for k in range(nz) for j in range(ny) for i in range(nx)
    ]
    indices_match = list(mesh.indices) == expected_indices

    spacing = [
        (hi - lo) / n for lo, hi, n in zip(ORDER_LOWER, ORDER_UPPER, ORDER_DIMENSION, strict=True)
    ]
    coordinate_lookup_ok = True
    coordinate_lookup_examples: list[dict[str, Any]] = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                center = (
                    ORDER_LOWER[0] + (i + 0.5) * spacing[0],
                    ORDER_LOWER[1] + (j + 0.5) * spacing[1],
                    ORDER_LOWER[2] + (k + 0.5) * spacing[2],
                )
                actual = mesh.get_indices_at_coords(center)
                expected = (i, j, k)
                if actual != expected:
                    coordinate_lookup_ok = False
                if len(coordinate_lookup_examples) < 4:
                    coordinate_lookup_examples.append(
                        {"center": center, "get_indices_at_coords": actual}
                    )

    return {
        "mesh_id": int(mesh.id),
        "dimension": list(mesh.dimension),
        "n_elements": int(mesh.n_elements),
        "source_count": int(flat_sources.size),
        "source_strengths_f_order": strengths_f_order,
        "source_energy_lower_edges_f_order": bands_f_order,
        "mesh_source_total_strength": float(mesh_source.strength),
        "expected_total_strength": float(sum(range(1, nx * ny * nz + 1))),
        "object_mapping_fortran_order_ok": mapping_ok,
        "regular_mesh_indices_are_1_based_fortran": indices_match,
        "get_indices_at_coords_is_0_based": coordinate_lookup_ok,
        "coordinate_lookup_examples": coordinate_lookup_examples,
        "mapping_rule": {
            "openmc_internal_flat_index": "i + nx*j + nx*ny*k (Fortran order, 0-based)",
            "nucleide_c_order_flat_index": "(i*ny + j)*nz + k (C order, 0-based)",
            "pass_c_order_array": (
                "reshape the C-ordered values to (nx, ny, nz) and pass that array; "
                "OpenMC ravels it with order='F'"
            ),
            "pass_flat_list": ("openmc_sources[i + nx*j + nx*ny*k] = c_sources[(i*ny + j)*nz + k]"),
        },
    }


def check_mesh_source_strength_normalization(
    openmc: Any, mesh: Any, mesh_source: Any
) -> dict[str, Any]:
    """Show that MeshSource strength is the sum and that normalizing scales children."""
    clone_sources = []
    for source in mesh_source.sources:
        clone_sources.append(
            openmc.IndependentSource(
                energy=source.energy,
                strength=float(source.strength),
                particle="photon",
            )
        )
    clone = openmc.MeshSource(mesh, clone_sources)
    before = [float(source.strength) for source in clone.sources]
    clone.normalize_source_strengths()
    after = [float(source.strength) for source in clone.sources]
    ratios_before = [value / sum(before) for value in before]
    ratios_after = [value / sum(after) for value in after]
    return {
        "mesh_source_total_strength": float(mesh_source.strength),
        "sum_of_element_strengths": float(sum(source.strength for source in mesh_source.sources)),
        "before_normalization": before,
        "after_normalize_source_strengths": after,
        "after_total_strength": float(clone.strength),
        "relative_selection_preserved": all(
            math.isclose(a, b, rel_tol=0.0, abs_tol=1e-12)
            for a, b in zip(ratios_before, ratios_after, strict=True)
        ),
        "note": (
            "normalize_source_strengths() makes the total 1.0 (per-source-particle "
            "tally convention); element ratios are unchanged"
        ),
    }


def read_settings_silencing_mesh_id_warning(openmc: Any, path: Path) -> Any:
    """Read settings while ignoring the expected duplicate mesh-id warning."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", openmc.IDWarning)
        return openmc.Settings.from_xml(path)


def check_xml_roundtrip(
    openmc: Any, ordering_mesh: Any, ordering_ms: Any, joint_mesh: Any, joint_ms: Any
) -> dict[str, Any]:
    """Export both MeshSources through settings.xml, read them back, and compare."""
    with tempfile.TemporaryDirectory(prefix="openmc-handoff-xml-") as tmp:
        directory = Path(tmp)
        settings_path = directory / "settings.xml"
        settings = openmc.Settings(
            particles=1,
            batches=1,
            run_mode="fixed source",
            output={"summary": False},
        )
        settings.source = [ordering_ms, joint_ms]
        settings.export_to_xml(settings_path)
        xml_text = settings_path.read_text(encoding="utf-8")

        restored = read_settings_silencing_mesh_id_warning(openmc, settings_path)
        by_mesh_id = {
            int(source.mesh.id): source
            for source in restored.source
            if isinstance(source, openmc.MeshSource)
        }
        ordering_restored = by_mesh_id.get(int(ordering_mesh.id))
        joint_restored = by_mesh_id.get(int(joint_mesh.id))

        ordering_sources_ok = ordering_restored is not None
        ordering_mesh_ok = ordering_restored is not None
        ordering_strengths_ok = ordering_restored is not None
        ordering_energies_ok = ordering_restored is not None
        if ordering_restored is not None:
            ordering_mesh_ok = (
                tuple(ordering_restored.mesh.dimension) == ORDER_DIMENSION
                and tuple(ordering_restored.mesh.lower_left) == ORDER_LOWER
                and tuple(ordering_restored.mesh.upper_right) == ORDER_UPPER
                and int(ordering_restored.mesh.id) == int(ordering_mesh.id)
            )
            ordering_sources_ok = int(ordering_restored.sources.size) == 12
            strengths = [float(source.strength) for source in ordering_restored.sources]
            original_strengths = [float(source.strength) for source in ordering_ms.sources]
            ordering_strengths_ok = strengths == original_strengths
            energies = []
            for source in ordering_restored.sources:
                if not isinstance(source.energy, openmc.stats.Uniform):
                    ordering_energies_ok = False
                    break
                energies.append((float(source.energy.a), float(source.energy.b)))
            original_energies = [
                (float(source.energy.a), float(source.energy.b)) for source in ordering_ms.sources
            ]
            ordering_energies_ok = ordering_energies_ok and energies == original_energies

        joint_mesh_ok = joint_restored is not None
        joint_strengths_ok = joint_restored is not None
        joint_edges_ok = joint_restored is not None
        if joint_restored is not None:
            joint_mesh_ok = (
                tuple(joint_restored.mesh.dimension) == (2, 1, 1)
                and tuple(joint_restored.mesh.lower_left) == (0.0, 0.0, 0.0)
                and tuple(joint_restored.mesh.upper_right) == (2.0, 1.0, 1.0)
                and int(joint_restored.mesh.id) == int(joint_mesh.id)
            )
            joint_strengths_ok = [float(source.strength) for source in joint_restored.sources] == [
                float(value) for value in JOINT_VOXEL_STRENGTHS
            ]
            joint_edges_ok = True
            for voxel, source in enumerate(joint_restored.sources):
                if not isinstance(source.energy, openmc.stats.Tabular):
                    joint_edges_ok = False
                    break
                if source.energy.interpolation != "histogram":
                    joint_edges_ok = False
                    break
                expected = histogram_density(JOINT_EDGES, JOINT_INTEGRATED[voxel])
                if source.energy.x.tolist() != list(JOINT_EDGES):
                    joint_edges_ok = False
                    break
                if source.energy.p.tolist() != expected:
                    joint_edges_ok = False
                    break

        # The <source type="mesh"> strength attribute is derived from its child
        # sources.  Change it on a temporary copy and confirm it is ignored.
        mesh_strength_attribute_ignored = None
        mesh_strength_pattern = (
            rf'(<source type="mesh" strength=")[0-9.]+(" mesh="{int(ordering_mesh.id)}">)'
        )
        modified_text, replacements = re.subn(
            mesh_strength_pattern,
            r"\g<1>999999.0\g<2>",
            xml_text,
            count=1,
        )
        if replacements == 1:
            modified_path = directory / "modified-settings.xml"
            modified_path.write_text(modified_text, encoding="utf-8")
            modified = read_settings_silencing_mesh_id_warning(openmc, modified_path)
            modified_by_mesh_id = {
                int(source.mesh.id): source
                for source in modified.source
                if isinstance(source, openmc.MeshSource)
            }
            modified_ordering = modified_by_mesh_id.get(int(ordering_mesh.id))
            if modified_ordering is not None:
                mesh_strength_attribute_ignored = math.isclose(
                    float(modified_ordering.strength),
                    float(ordering_ms.strength),
                    rel_tol=0.0,
                    abs_tol=1e-12,
                )

        return {
            "route": "openmc.Settings.export_to_xml() -> openmc.Settings.from_xml()",
            "settings_xml_sha256": sha256_file(settings_path),
            "mesh_elements_in_xml": len(re.findall(r"<mesh id=", xml_text)),
            "mesh_sources_in_xml": len(re.findall(r'<source type="mesh"', xml_text)),
            "restored_source_count": len(restored.source),
            "source_order_in_xml": (
                "OpenMC internal Fortran order (MeshSource.sources after its "
                "order='F' ravel), not the caller's C order"
            ),
            "source_order_preserved": bool(ordering_strengths_ok and ordering_energies_ok),
            "ordering": {
                "mesh_ok": ordering_mesh_ok,
                "source_count_ok": ordering_sources_ok,
                "strengths_ok": ordering_strengths_ok,
                "energy_edges_ok": ordering_energies_ok,
            },
            "joint": {
                "mesh_ok": joint_mesh_ok,
                "source_count_ok": joint_restored is not None
                and int(joint_restored.sources.size) == 2,
                "strengths_ok": joint_strengths_ok,
                "energy_edges_ok": joint_edges_ok,
            },
            "mesh_source_strength_attribute_ignored": mesh_strength_attribute_ignored,
            "roundtrip_ok": all(
                [
                    ordering_mesh_ok,
                    ordering_sources_ok,
                    ordering_strengths_ok,
                    ordering_energies_ok,
                    joint_mesh_ok,
                    joint_restored is not None,
                    joint_strengths_ok,
                    joint_edges_ok,
                ]
            ),
        }


def write_cross_sections_placeholder(directory: Path) -> Path:
    """Write the smallest cross_sections.xml that lets OpenMC initialize."""
    path = directory / "cross_sections_placeholder.xml"
    path.write_text(
        '<?xml version="1.0"?>\n'
        "<cross_sections>\n"
        '  <library materials="H1" path="/tmp/nonexistent-openmc-handoff.h5" '
        'type="neutron"/>\n'
        "</cross_sections>\n",
        encoding="utf-8",
    )
    return path


def redirect_c_stdout_to_stderr() -> int:
    """Duplicate fd 1 onto stderr so OpenMC C++ output cannot corrupt JSON."""
    sys.stdout.flush()
    saved = os.dup(1)
    os.dup2(2, 1)
    return saved


def restore_c_stdout(saved: int) -> None:
    """Restore the original stdout file descriptor."""
    sys.stdout.flush()
    os.dup2(saved, 1)
    os.close(saved)


def sample_with_openmc(
    openmc: Any, mesh: Any, mesh_source: Any, n_samples: int, seed: int, sphere_radius: float
) -> Any:
    """Sample source sites through the installed OpenMC library, without transport."""
    settings = openmc.Settings(
        particles=1,
        batches=1,
        run_mode="fixed source",
        output={"summary": False},
    )
    settings.source = [mesh_source]
    sphere = openmc.Sphere(r=sphere_radius, boundary_type="vacuum")
    cell = openmc.Cell(region=-sphere)
    geometry = openmc.Geometry([cell])
    model = openmc.Model(settings=settings, geometry=geometry)
    saved = redirect_c_stdout_to_stderr()
    try:
        particles = model.sample_external_source(
            n_samples=n_samples,
            prn_seed=seed,
            as_array=True,
        )
    finally:
        restore_c_stdout(saved)
    return particles


def sample_joint(openmc: Any) -> dict[str, Any]:
    """Sample the two-voxel/two-group scenario and compare joint frequencies."""
    import numpy as np

    mesh, mesh_source = build_joint_mesh(openmc, mesh_id=19101)
    particles = sample_with_openmc(
        openmc,
        mesh,
        mesh_source,
        n_samples=JOINT_SAMPLES,
        seed=PRN_SEED,
        sphere_radius=100.0,
    )
    positions = np.asarray(particles["r"])
    indices = np.asarray([mesh.get_indices_at_coords(position) for position in positions])
    voxel = indices[:, 0]
    energies = np.asarray(particles["E"])
    group = np.searchsorted(np.asarray(JOINT_EDGES), energies, side="right") - 1

    joint = np.zeros((2, 2), dtype=np.int64)
    np.add.at(joint, (voxel, group), 1)
    empirical = joint / float(JOINT_SAMPLES)
    exact = np.asarray(JOINT_EXACT)
    max_abs_error = float(np.max(np.abs(empirical - exact)))
    zero_cross_count = int(joint[0, 1] + joint[1, 0])
    weights = np.asarray(particles["wgt"])
    inside = bool(
        np.all(voxel >= 0)
        and np.all(voxel < 2)
        and np.all(group >= 0)
        and np.all(group < 2)
        and np.all(np.isfinite(energies))
    )
    weights_unity = bool(np.all(weights == 1.0))
    passed = max_abs_error <= JOINT_TOLERANCE and zero_cross_count == 0 and weights_unity and inside
    return {
        "n_samples": JOINT_SAMPLES,
        "seed": PRN_SEED,
        "tolerance_abs": JOINT_TOLERANCE,
        "exact_joint": exact.tolist(),
        "counts": joint.tolist(),
        "empirical_joint": empirical.tolist(),
        "max_abs_error": max_abs_error,
        "zero_cross_count": zero_cross_count,
        "inside_mesh": inside,
        "weights": {
            "min": float(np.min(weights)),
            "max": float(np.max(weights)),
            "unique": np.unique(weights).tolist(),
            "all_unity": weights_unity,
        },
        "particle_codes": np.unique(np.asarray(particles["particle"])).tolist(),
        "mesh_total_strength": float(mesh_source.strength),
        "relative_selection_probabilities": [
            float(source.strength) / float(mesh_source.strength) for source in mesh_source.sources
        ],
        "pass": bool(passed),
    }


def sample_ordering(openmc: Any) -> dict[str, Any]:
    """Sample the 2 x 3 x 2 ordering mesh and check the C-to-F mapping in practice."""
    import numpy as np

    mesh, _, mesh_source = build_ordering_mesh(openmc, mesh_id=19102)
    particles = sample_with_openmc(
        openmc,
        mesh,
        mesh_source,
        n_samples=ORDER_SAMPLES,
        seed=PRN_SEED,
        sphere_radius=100.0,
    )
    positions = np.asarray(particles["r"])
    indices = np.asarray([mesh.get_indices_at_coords(position) for position in positions])
    i = indices[:, 0]
    j = indices[:, 1]
    k = indices[:, 2]
    c_index = (i * ORDER_DIMENSION[1] + j) * ORDER_DIMENSION[2] + k
    f_index = i + ORDER_DIMENSION[0] * j + ORDER_DIMENSION[0] * ORDER_DIMENSION[1] * k
    energies = np.asarray(particles["E"])
    band = np.floor((energies - ORDER_BAND_LOW) / ORDER_BAND_WIDTH).astype(np.int64)
    band_in_range = bool(np.all(band >= 0) and np.all(band < 12))
    mismatch_count = int(np.count_nonzero(band != c_index))
    counts_c = np.bincount(c_index, minlength=12)
    counts_f = np.bincount(f_index, minlength=12)
    weights = np.asarray(particles["wgt"])
    all_cells_seen = bool(np.all(counts_c > 0))
    weights_unity = bool(np.all(weights == 1.0))
    passed = (
        band_in_range
        and mismatch_count == 0
        and all_cells_seen
        and weights_unity
        and bool(np.all(i >= 0))
        and bool(np.all(i < 2))
        and bool(np.all(j >= 0))
        and bool(np.all(j < 3))
        and bool(np.all(k >= 0))
        and bool(np.all(k < 2))
    )
    return {
        "n_samples": ORDER_SAMPLES,
        "seed": PRN_SEED,
        "mismatch_count": mismatch_count,
        "band_in_range": band_in_range,
        "all_12_cells_seen": all_cells_seen,
        "counts_by_c_order_index": counts_c.tolist(),
        "counts_by_openmc_f_order_index": counts_f.tolist(),
        "weights_unity": weights_unity,
        "pass": bool(passed),
    }


def sample_zero_voxel(openmc: Any) -> dict[str, Any]:
    """Sample a mesh with one nonzero and one zero-strength element."""
    import numpy as np

    mesh = openmc.RegularMesh()
    mesh.id = 19103
    mesh.dimension = (2, 1, 1)
    mesh.lower_left = (0.0, 0.0, 0.0)
    mesh.upper_right = (2.0, 1.0, 1.0)
    sources = [
        openmc.IndependentSource(
            energy=openmc.stats.Uniform(1.0, 2.0), strength=1.0, particle="photon"
        ),
        openmc.IndependentSource(
            energy=openmc.stats.Uniform(3.0, 4.0), strength=0.0, particle="photon"
        ),
    ]
    mesh_source = openmc.MeshSource(mesh, sources)
    particles = sample_with_openmc(
        openmc,
        mesh,
        mesh_source,
        n_samples=ZERO_SAMPLES,
        seed=PRN_SEED,
        sphere_radius=100.0,
    )
    positions = np.asarray(particles["r"])
    indices = np.asarray([mesh.get_indices_at_coords(position) for position in positions])
    voxel = indices[:, 0]
    counts = np.bincount(voxel, minlength=2)
    weights = np.asarray(particles["wgt"])
    passed = (
        int(counts[0]) == ZERO_SAMPLES
        and int(counts[1]) == 0
        and bool(np.all(weights == 1.0))
        and float(mesh_source.strength) == 1.0
    )
    return {
        "n_samples": ZERO_SAMPLES,
        "seed": PRN_SEED,
        "counts": counts.tolist(),
        "element_strengths": [float(source.strength) for source in mesh_source.sources],
        "mesh_total_strength": float(mesh_source.strength),
        "weights": {
            "min": float(np.min(weights)),
            "max": float(np.max(weights)),
            "all_unity": bool(np.all(weights == 1.0)),
        },
        "pass": bool(passed),
    }


def run_void_tally_worker() -> int:
    """Child-process entry point: run a void fixed-source tally normalization check.

    The model is a vacuum sphere with an isotropic point source at its center.
    Every particle travels exactly one radius before leaving, so the
    volume-integrated track-length flux tally mean is ``radius * total source
    strength``. Running at strengths 1 and 10 distinguishes a per-source-
    particle tally (unchanged except by statistics) from a rate-scaled tally
    (ten times larger), and distinguishes an integrated flux score from a
    volume-averaged one. No evaluated data or interacting material is used.
    """
    try:
        import openmc
    except Exception as error:  # pragma: no cover - only for broken environments
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

    result: dict[str, Any] = {
        "radius_cm": VOID_RADIUS,
        "n_particles": VOID_PARTICLES,
        "seed": VOID_SEED,
        "strengths": list(VOID_STRENGTHS),
        "runs": [],
        "transport_run": "void-track-length-only",
    }
    with tempfile.TemporaryDirectory(prefix="openmc-void-tally-") as tmp:
        directory = Path(tmp)
        cross_sections = write_cross_sections_placeholder(directory)
        previous = os.environ.get("OPENMC_CROSS_SECTIONS")
        os.environ["OPENMC_CROSS_SECTIONS"] = str(cross_sections)
        try:
            for strength in VOID_STRENGTHS:
                model_dir = directory / f"strength-{strength:g}"
                model_dir.mkdir()
                sphere = openmc.Sphere(r=VOID_RADIUS, boundary_type="vacuum")
                cell = openmc.Cell(region=-sphere)
                geometry = openmc.Geometry([cell])
                source = openmc.IndependentSource(
                    space=openmc.stats.Point((0.0, 0.0, 0.0)),
                    energy=openmc.stats.Discrete([14.0e6], [1.0]),
                    strength=strength,
                    particle="neutron",
                )
                settings = openmc.Settings(
                    particles=VOID_PARTICLES,
                    batches=1,
                    run_mode="fixed source",
                    output={"summary": False, "tallies": False},
                )
                settings.seed = VOID_SEED
                settings.source = [source]
                tally = openmc.Tally(name="void flux")
                tally.filters = [openmc.CellFilter(cell)]
                tally.scores = ["flux"]
                model = openmc.Model(
                    geometry=geometry,
                    materials=openmc.Materials([]),
                    settings=settings,
                    tallies=openmc.Tallies([tally]),
                )
                model.export_to_xml(model_dir)
                saved = redirect_c_stdout_to_stderr()
                try:
                    openmc.run(output=False, cwd=str(model_dir))
                finally:
                    restore_c_stdout(saved)
                statepoint = openmc.StatePoint(str(model_dir / "statepoint.1.h5"))
                try:
                    data = statepoint.tallies[tally.id]
                    mean = float(data.mean.ravel()[0])
                    std_dev = float(data.std_dev.ravel()[0])
                    if not math.isfinite(std_dev):
                        # A zero-variance void tally reports NaN in 0.16.0;
                        # keep the JSON document strict.
                        std_dev = None
                finally:
                    close = getattr(statepoint, "close", None)
                    if callable(close):
                        close()
                result["runs"].append(
                    {
                        "source_strength": strength,
                        "mean": mean,
                        "std_dev": std_dev,
                        "expected_integrated_track_length": VOID_RADIUS * strength,
                        "expected_if_volume_averaged": (
                            VOID_RADIUS * strength / (4.0 / 3.0 * math.pi * VOID_RADIUS**3)
                        ),
                    }
                )
        finally:
            if previous is None:
                os.environ.pop("OPENMC_CROSS_SECTIONS", None)
            else:
                os.environ["OPENMC_CROSS_SECTIONS"] = previous
    print(json.dumps(result, sort_keys=True))
    return 0


def run_sample_worker() -> int:
    """Child-process entry point: sample sources, print one JSON document."""
    try:
        import openmc
    except Exception as error:  # pragma: no cover - only for broken environments
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

    with tempfile.TemporaryDirectory(prefix="openmc-handoff-sample-") as tmp:
        directory = Path(tmp)
        cross_sections = write_cross_sections_placeholder(directory)
        previous = os.environ.get("OPENMC_CROSS_SECTIONS")
        os.environ["OPENMC_CROSS_SECTIONS"] = str(cross_sections)
        try:
            result: dict[str, Any] = {
                "cross_sections_placeholder": {
                    "path": str(cross_sections),
                    "sha256": sha256_file(cross_sections),
                    "library_path": "/tmp/nonexistent-openmc-handoff.h5",
                    "library_path_exists": Path("/tmp/nonexistent-openmc-handoff.h5").exists(),
                    "note": (
                        "Initialization shim only. OpenMC runs in volume-calculation "
                        "mode; no cross-section data are read and no transport run occurs."
                    ),
                },
                "transport_run": False,
            }
            result["joint"] = sample_joint(openmc)
            result["ordering"] = sample_ordering(openmc)
            result["zero_voxel"] = sample_zero_voxel(openmc)
        finally:
            if previous is None:
                os.environ.pop("OPENMC_CROSS_SECTIONS", None)
            else:
                os.environ["OPENMC_CROSS_SECTIONS"] = previous
    print(json.dumps(result, sort_keys=True))
    return 0


def record_check(checks: dict[str, Any], name: str, passed: bool, **details: Any) -> None:
    """Record a named boolean check with supporting detail."""
    checks[name] = {"pass": bool(passed), **details}


def build_summary(openmc: Any) -> dict[str, Any]:
    """Run all Python-level checks and the isolated sampling worker."""
    checks: dict[str, Any] = {}
    summary: dict[str, Any] = {
        "probe": PROBE_NAME,
        "status": "FAIL",
        "openmc_version": openmc.__version__,
        "python_version": sys.version.split()[0],
        "openmc_package_file": str(Path(openmc.__file__).resolve()),
        "transport_run": False,
        "void_track_length_check": True,
        "limitations": [
            "No evaluated cross-section or photon data are installed. The only "
            "transport run is a vacuum track-length tally normalization check; no "
            "interacting transport calculation or physical tally result is produced "
            "or claimed.",
            "Source sampling verifies the OpenMC source API and joint probabilities only.",
            "The placeholder cross_sections.xml is an initialization shim, not physics input.",
            "The 2 x 3 x 2 and two-voxel scenarios are synthetic; they are not an ACTINV case.",
            "An all-zero MeshSource was not sampled: OpenMC 0.16.0 DiscreteIndex::normalize "
            "divides by integral_ with no zero guard. A single zero-strength voxel inside a "
            "nonzero source is verified; the adapter must reject an all-zero total strength "
            "itself.",
        ],
    }

    try:
        metadata_version = importlib.metadata.version("openmc")
    except importlib.metadata.PackageNotFoundError:
        metadata_version = "not installed as a distribution"
    version_ok = openmc.__version__ == EXPECTED_OPENMC and metadata_version == EXPECTED_OPENMC
    record_check(
        checks,
        "openmc_version_exactly_0_16_0",
        version_ok,
        runtime_version=openmc.__version__,
        distribution_version=metadata_version,
        expected=EXPECTED_OPENMC,
    )
    summary["openmc_distribution_version"] = metadata_version

    source_dir, source_origin = find_openmc_source_dir()
    summary["openmc_source"] = {
        "path": str(source_dir) if source_dir is not None else None,
        "origin": source_origin,
    }
    provenance = collect_source_provenance(source_dir)
    summary["normalization"] = {
        "inspected_source": provenance,
        "inspected_python": collect_python_provenance(openmc),
    }
    record_check(
        checks,
        "openmc_source_tree_available",
        bool(provenance.get("available")),
        source_dir=str(source_dir) if source_dir is not None else None,
        reason=provenance.get("reason"),
    )

    derived = provenance.get("derived", {}) if provenance.get("available") else {}
    per_source_particle_docs = (
        provenance.get("available", False)
        and provenance["tallies_rst"]["patterns"]["per_source_particle"]["found"]
    )
    strength_unit_conversion_docs = (
        provenance.get("available", False)
        and provenance["tallies_rst"]["patterns"]["strength_as_unit_conversion"]["found"]
    )
    relative_strength_docs = (
        provenance.get("available", False)
        and provenance["settings_rst"]["patterns"]["relative_strength"]["found"]
    )
    regular_norm_found = (
        provenance.get("available", False)
        and provenance["tally_cpp"]["patterns"]["regular_norm"]["found"]
    )
    global_total_weight_found = (
        provenance.get("available", False)
        and provenance["tally_cpp"]["patterns"]["global_total_weight"]["found"]
    )
    regular_uses_total_weight = bool(
        derived.get("regular_norm_uses_simulation_total_weight", False)
    )
    regular_uses_n_particles = bool(derived.get("regular_norm_uses_n_particles", False))
    regular_scales_by_total_source = bool(derived.get("regular_norm_scales_by_total_source", False))
    mesh_strength_sum = bool(derived.get("mesh_source_strength_is_sum", False))
    discrete_selection = (
        provenance.get("available", False)
        and provenance["distribution_cpp"]["patterns"]["discrete_integral"]["found"]
        and provenance["distribution_cpp"]["patterns"]["discrete_normalize"]["found"]
    )
    mesh_element_selection = (
        provenance.get("available", False)
        and provenance["source_cpp"]["patterns"]["mesh_sample_element"]["found"]
        and provenance["source_cpp"]["patterns"]["mesh_strengths"]["found"]
    )

    record_check(
        checks,
        "tally_results_use_per_source_particle_denominator",
        bool(per_source_particle_docs and regular_norm_found and regular_uses_n_particles),
        expected=True,
        observed=bool(per_source_particle_docs and regular_uses_n_particles),
        normalization_expression=(
            "norm = total_source / (contributing_particles * settings::gen_per_batch); "
            "contributing_particles = settings::n_particles"
        ),
        detail=(
            "OpenMC 0.16 docs and Tally::accumulate confirm a per-source-particle "
            "denominator. The numerator is the total source strength, so the result is "
            "per-source-particle x total source strength."
        ),
    )
    record_check(
        checks,
        "regular_tallies_do_not_use_simulation_total_weight_denominator",
        not regular_uses_total_weight,
        expected=True,
        observed=not regular_uses_total_weight,
        detail=(
            "Regular tallies in 0.16.0 use settings::n_particles, not "
            "simulation::total_weight. simulation::total_weight is used for global "
            "tallies (k/leakage) and MPI reduction. In default unweighted sampling "
            "total_weight equals n_particles numerically, but it is not the regular-"
            "tally denominator."
        ),
    )
    record_check(
        checks,
        "source_strength_serves_as_relative_selection_weight",
        bool(relative_strength_docs and discrete_selection and mesh_element_selection),
        expected=True,
        observed=bool(relative_strength_docs and discrete_selection and mesh_element_selection),
        detail=(
            "Per-element strength is used as a relative selection weight: DiscreteIndex "
            "normalizes it, and MeshSource samples an element by relative strength."
        ),
    )
    record_check(
        checks,
        "source_strength_is_relative_weight_and_tally_scale",
        bool(
            relative_strength_docs
            and regular_scales_by_total_source
            and strength_unit_conversion_docs
        ),
        expected=True,
        observed=bool(
            relative_strength_docs
            and regular_scales_by_total_source
            and strength_unit_conversion_docs
        ),
        detail=(
            "Source strength is both a relative selection weight and the tally scale: "
            "the sum of source strengths is total_source in Tally::accumulate, and the "
            "docs describe using SourceBase.strength to convert per-source-particle "
            "tallies to absolute units."
        ),
    )
    record_check(
        checks,
        "source_strength_sum_scales_tally_results",
        bool(regular_scales_by_total_source and strength_unit_conversion_docs),
        expected=True,
        observed=bool(regular_scales_by_total_source and strength_unit_conversion_docs),
        detail=(
            "Observed behavior: total_source = external_sources_probability.integral() "
            "= sum of source strengths; MeshSource.strength = sum of element strengths."
        ),
    )
    record_check(
        checks,
        "global_tallies_use_simulation_total_weight",
        bool(global_total_weight_found),
        expected=True,
        observed=bool(global_total_weight_found),
        detail="Global tally accumulation divides by simulation::total_weight.",
    )
    record_check(
        checks,
        "mesh_source_strength_is_sum_of_element_strengths",
        mesh_strength_sum,
        expected=True,
        observed=mesh_strength_sum,
        detail="openmc.MeshSource.strength returns sum(s.strength for s in self.sources).",
    )

    ordering_mesh, sources_c, ordering_ms = build_ordering_mesh(openmc, mesh_id=19001)
    joint_mesh, joint_ms = build_joint_mesh(openmc, mesh_id=19002)
    ordering = check_ordering_objects(openmc, ordering_mesh, sources_c, ordering_ms)
    summary["mapping"] = ordering
    record_check(
        checks,
        "c_order_array_maps_to_expected_openmc_mesh_elements",
        bool(ordering["object_mapping_fortran_order_ok"]),
        expected=True,
        observed=ordering["object_mapping_fortran_order_ok"],
        mapping_rule=ordering["mapping_rule"],
    )
    record_check(
        checks,
        "regular_mesh_indices_fortran_order_1_based",
        bool(ordering["regular_mesh_indices_are_1_based_fortran"]),
        expected=True,
        observed=ordering["regular_mesh_indices_are_1_based_fortran"],
    )
    record_check(
        checks,
        "get_indices_at_coords_zero_based",
        bool(ordering["get_indices_at_coords_is_0_based"]),
        expected=True,
        observed=ordering["get_indices_at_coords_is_0_based"],
    )
    summary["normalization"]["mesh_source_strength"] = check_mesh_source_strength_normalization(
        openmc, joint_mesh, joint_ms
    )
    summary["normalization"]["observed"] = {
        "regular_tally_denominator": (
            "settings::n_particles * settings::gen_per_batch "
            "(settings::work_per_rank when reduce_tallies is false)"
        ),
        "regular_tally_scale": (
            "external_sources_probability.integral(), the sum of source strengths"
        ),
        "simulation_total_weight_role": (
            "global tallies and MPI reduction; accumulated from p.wgt() for each "
            "source particle and numerically equal to n_particles in default "
            "unweighted sampling"
        ),
        "mesh_source_total_strength": float(joint_ms.strength),
        "mesh_source_element_strengths": [float(source.strength) for source in joint_ms.sources],
        "mesh_source_element_selection_probabilities": [
            float(source.strength) / float(joint_ms.strength) for source in joint_ms.sources
        ],
        "per_source_particle_conversion": (
            "If element strengths sum to 1.0, raw tallies are per source particle and "
            "physical_flux = raw * source_rate_s / volume_cm3. If element strengths are "
            "absolute photons/s, their sum is already applied by OpenMC and multiplying "
            "by source_rate_s again would double-count it."
        ),
    }

    histogram = build_histogram_case(openmc)
    summary["histogram"] = histogram
    density_ok = (
        math.isclose(histogram["integral"], 1.0, rel_tol=0.0, abs_tol=1e-12)
        and all(
            math.isclose(actual, expected, rel_tol=0.0, abs_tol=1e-12)
            for actual, expected in zip(
                histogram["readback_p"],
                histogram["expected_density_pdf_per_eV"],
                strict=True,
            )
        )
        and histogram["interpolation"] == "histogram"
        and math.isclose(histogram["wrong_control_integral"], 4.0, rel_tol=0.0, abs_tol=1e-12)
    )
    record_check(
        checks,
        "histogram_p_is_differential_density_per_eV",
        bool(density_ok),
        expected=True,
        observed=bool(density_ok),
        integral=histogram["integral"],
        wrong_control_integral=histogram["wrong_control_integral"],
    )
    zero_bin_ok = (
        histogram["zero_bin"]["integrated_strength"] == 0.0
        and histogram["zero_bin"]["density"] == 0.0
        and histogram["zero_bin"]["evaluate_midpoint"] == 0.0
        and bool(histogram["zero_bin"]["cdf_flat"])
    )
    record_check(
        checks,
        "zero_bin_has_zero_density_and_zero_integral_contribution",
        bool(zero_bin_ok),
        expected=True,
        observed=bool(zero_bin_ok),
    )
    zero_voxel_ok = (
        histogram["zero_strength_voxel"]["source_strength"] == 0.0
        and histogram["zero_strength_voxel"]["mesh_source_total_strength"] == 1.0
        and histogram["zero_strength_voxel"]["relative_probabilities"] == [1.0, 0.0]
        and not histogram["all_zero_distribution"]["pdf_valid"]
    )
    record_check(
        checks,
        "zero_strength_voxel_is_data_but_not_samplable_pdf",
        bool(zero_voxel_ok),
        expected=True,
        observed=bool(zero_voxel_ok),
    )
    record_check(
        checks,
        "openmc_stats_Histogram_class_absent_use_Tabular",
        not histogram["histogram_class_available"],
        expected=False,
        observed=histogram["histogram_class_available"],
        supported_equivalent=histogram["supported_equivalent"],
        note="The brief's openmc.stats.Histogram name does not exist in 0.16.0.",
    )

    roundtrip = check_xml_roundtrip(openmc, ordering_mesh, ordering_ms, joint_mesh, joint_ms)
    summary["xml_roundtrip"] = roundtrip
    record_check(
        checks,
        "settings_xml_roundtrip_preserves_mesh_sources",
        bool(roundtrip["roundtrip_ok"]),
        expected=True,
        observed=bool(roundtrip["roundtrip_ok"]),
        route=roundtrip["route"],
    )
    record_check(
        checks,
        "mesh_source_strength_attribute_is_derived_not_authoritative",
        roundtrip["mesh_source_strength_attribute_ignored"] is True,
        expected=True,
        observed=roundtrip["mesh_source_strength_attribute_ignored"],
        note=(
            'Changing the <source type="mesh" strength=...> attribute on a temporary '
            "XML copy does not change the read-back MeshSource.strength."
        ),
    )

    worker_command = [sys.executable, str(Path(__file__).resolve()), "--sample-worker"]
    worker: dict[str, Any] = {}
    worker_error = None
    try:
        completed = subprocess.run(
            worker_command,
            capture_output=True,
            text=True,
            timeout=180,
            check=False,
        )
        summary["sampling_worker"] = {
            "command": worker_command,
            "returncode": completed.returncode,
            "stdout_tail": completed.stdout[-4000:],
            "stderr_tail": completed.stderr[-4000:],
        }
        if completed.returncode == 0:
            worker = json.loads(completed.stdout)
        else:
            worker_error = f"sampling worker exited {completed.returncode}"
    except (OSError, subprocess.SubprocessError, json.JSONDecodeError) as error:
        worker_error = f"sampling worker failed: {type(error).__name__}: {error}"
        summary["sampling_worker"] = {"command": worker_command, "error": worker_error}

    summary["sampling"] = worker
    worker_ok = bool(worker) and worker_error is None
    record_check(
        checks,
        "source_sampling_worker_succeeded",
        worker_ok,
        expected=True,
        observed=worker_ok,
        error=worker_error,
    )
    joint_sampling = worker.get("joint", {}) if worker_ok else {}
    joint_sampling_ok = bool(joint_sampling.get("pass", False)) if worker_ok else False
    record_check(
        checks,
        "sampled_joint_voxel_group_probabilities_match",
        joint_sampling_ok,
        expected=True,
        observed=joint_sampling_ok,
        exact_joint=joint_sampling.get("exact_joint"),
        empirical_joint=joint_sampling.get("empirical_joint"),
        counts=joint_sampling.get("counts"),
        tolerance_abs=JOINT_TOLERANCE,
        n_samples=JOINT_SAMPLES,
        seed=PRN_SEED,
    )
    zero_cross_ok = worker_ok and int(joint_sampling.get("zero_cross_count", -1)) == 0
    record_check(
        checks,
        "sampled_zero_cross_bins_are_never_selected",
        zero_cross_ok,
        expected=True,
        observed=zero_cross_ok,
        zero_cross_count=joint_sampling.get("zero_cross_count"),
    )
    weights_ok = worker_ok and bool(joint_sampling.get("weights", {}).get("all_unity", False))
    record_check(
        checks,
        "default_sampling_particle_weights_are_unity",
        weights_ok,
        expected=True,
        observed=weights_ok,
        weights=joint_sampling.get("weights"),
        detail="Source strength selects a voxel; it is not applied as a particle weight.",
    )
    ordering_sampling = worker.get("ordering", {}) if worker_ok else {}
    ordering_sampling_ok = bool(ordering_sampling.get("pass", False)) if worker_ok else False
    record_check(
        checks,
        "sampled_2x3x2_cells_match_c_order_band_mapping",
        ordering_sampling_ok,
        expected=True,
        observed=ordering_sampling_ok,
        mismatch_count=ordering_sampling.get("mismatch_count"),
        counts_by_c_order_index=ordering_sampling.get("counts_by_c_order_index"),
        n_samples=ORDER_SAMPLES,
        seed=PRN_SEED,
    )
    zero_voxel_sampling = worker.get("zero_voxel", {}) if worker_ok else {}
    zero_voxel_sampling_ok = bool(zero_voxel_sampling.get("pass", False)) if worker_ok else False
    record_check(
        checks,
        "sampled_zero_strength_voxel_is_never_selected",
        zero_voxel_sampling_ok,
        expected=True,
        observed=zero_voxel_sampling_ok,
        counts=zero_voxel_sampling.get("counts"),
        n_samples=ZERO_SAMPLES,
        seed=PRN_SEED,
    )

    void_command = [sys.executable, str(Path(__file__).resolve()), "--void-tally-worker"]
    void_worker: dict[str, Any] = {}
    void_error = None
    try:
        completed = subprocess.run(
            void_command,
            capture_output=True,
            text=True,
            timeout=180,
            check=False,
        )
        summary["void_tally_worker"] = {
            "command": void_command,
            "returncode": completed.returncode,
            "stdout_tail": completed.stdout[-4000:],
            "stderr_tail": completed.stderr[-4000:],
        }
        if completed.returncode == 0:
            void_worker = json.loads(completed.stdout)
        else:
            void_error = f"void tally worker exited {completed.returncode}"
    except (OSError, subprocess.SubprocessError, json.JSONDecodeError) as error:
        void_error = f"void tally worker failed: {type(error).__name__}: {error}"
        summary["void_tally_worker"] = {"command": void_command, "error": void_error}

    summary["void_tally"] = void_worker
    void_runs = void_worker.get("runs", []) if void_worker else []
    void_exact = len(void_runs) == len(VOID_STRENGTHS) and all(
        math.isclose(
            float(run["mean"]),
            float(run["expected_integrated_track_length"]),
            rel_tol=1e-9,
            abs_tol=1e-12,
        )
        and not math.isclose(
            float(run["mean"]),
            float(run["expected_if_volume_averaged"]),
            rel_tol=1e-3,
            abs_tol=0.0,
        )
        for run in void_runs
    )
    ratio_ok = len(void_runs) == 2 and float(void_runs[1]["mean"]) == 10.0 * float(
        void_runs[0]["mean"]
    )
    record_check(
        checks,
        "void_fixed_source_flux_tally_is_integrated_and_rate_scaled",
        bool(void_error is None and void_exact and ratio_ok),
        expected=True,
        observed=bool(void_error is None and void_exact and ratio_ok),
        runs=void_runs,
        expected_integrated_track_length="radius * total source strength",
        expected_if_volume_averaged="radius * strength / sphere volume",
        error=void_error,
        detail=(
            "A vacuum sphere with a centered isotropic point source has exactly one "
            "track of length R per particle. The tally mean is R * strength, not "
            "R * strength / volume, and the 10x strength run is 10x the unit run."
        ),
    )

    summary["api_mismatches"] = [
        {
            "api": "openmc.stats.Histogram",
            "status": "absent in 0.16.0",
            "use_instead": "openmc.stats.Tabular(x, p, interpolation='histogram')",
            "detail": "Tabular p values are densities per unit x, not integrated bin strengths.",
        },
        {
            "api": "regular tally normalization",
            "status": "brief premise mismatch",
            "detail": (
                "Regular tallies divide by settings::n_particles and multiply by the "
                "sum of source strengths; they do not divide by simulation::total_weight."
            ),
        },
        {
            "api": "SourceBase.strength / MeshSource.strength",
            "status": "brief premise mismatch",
            "detail": (
                "Element strengths are relative selection weights, but their sum is the "
                "tally normalization scale. A mesh source with absolute photons/s "
                "strengths therefore produces tallies already multiplied by the total rate."
            ),
        },
        {
            "api": "openmc.MeshSource.sample",
            "status": "no such method",
            "detail": (
                "Use openmc.Model.sample_external_source or openmc.lib.sample_external_source. "
                "OpenMC must be initialized; source sites must lie inside the geometry."
            ),
        },
        {
            "api": "RegularMesh index bases",
            "status": "mixed",
            "detail": (
                "RegularMesh.indices yields 1-based (i, j, k) tuples in Fortran order; "
                "RegularMesh.get_indices_at_coords returns 0-based indices."
            ),
        },
        {
            "api": "MeshSource.sources ravel order",
            "status": "Fortran order",
            "detail": (
                "A multidimensional source array is flattened with order='F'. The "
                "bridge's C-order voxel index is (i*ny + j)*nz + k."
            ),
        },
        {
            "api": "MeshSource XML strength attribute",
            "status": "derived",
            "detail": (
                'The <source type="mesh" strength=...> attribute is the sum of child '
                "source strengths; from_xml_element rebuilds the sum from children."
            ),
        },
        {
            "api": "all-zero MeshSource",
            "status": "no explicit rejection in the inspected source path",
            "detail": (
                "DiscreteIndex::normalize divides by integral_ with no zero-integral guard. "
                "A zero element inside a nonzero source is valid, but all-zero total "
                "strength must be rejected by the adapter before export or sampling. It "
                "was not sampled here."
            ),
        },
    ]

    failed = [name for name, value in checks.items() if not value["pass"]]
    summary["checks"] = checks
    summary["failed_checks"] = failed
    summary["status"] = "PASS" if not failed else "FAIL"
    summary["exit_code"] = 0 if not failed else 1
    return summary


def main() -> int:
    """Validate the environment, run the probe, and print one JSON summary."""
    try:
        import openmc
    except Exception as error:
        reason = f"OpenMC import failed: {type(error).__name__}: {error}"
        print(f"SKIP: {reason}", file=sys.stderr)
        print(json.dumps({"probe": PROBE_NAME, "status": "SKIP", "reason": reason}))
        return 2

    if openmc.__version__ != EXPECTED_OPENMC:
        reason = f"OpenMC {openmc.__version__} is installed; expected exactly {EXPECTED_OPENMC}"
        print(f"SKIP: {reason}", file=sys.stderr)
        print(json.dumps({"probe": PROBE_NAME, "status": "SKIP", "reason": reason}))
        return 2

    try:
        summary = build_summary(openmc)
    except Exception as error:
        summary = {
            "probe": PROBE_NAME,
            "status": "FAIL",
            "openmc_version": getattr(openmc, "__version__", None),
            "error": f"{type(error).__name__}: {error}",
            "traceback": traceback.format_exc(),
            "exit_code": 1,
        }
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 1

    print(json.dumps(summary, indent=2, sort_keys=True))
    return int(summary["exit_code"])


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--sample-worker":
        sys.exit(run_sample_worker())
    if len(sys.argv) == 2 and sys.argv[1] == "--void-tally-worker":
        sys.exit(run_void_tally_worker())
    sys.exit(main())
