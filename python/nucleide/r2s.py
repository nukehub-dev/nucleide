"""Rigorous two-step (R2S) shutdown-dose-rate orchestration (backed by the `nucleide-r2s` crate).

Combines neutron flux meshes, ALARA activation decks, decay photon sources,
and mesh source sampling into one reproducible workflow description.
Transport and activation solving stay inside their respective codes; this
module only builds decks, links zones to spectra, and assembles photon
sources.

``tag_zone_strength(..., split=True, voxel_volumes_cm3=[...])`` splits
integrated zone strengths (photons/s) at uniform within-zone density. Positive
finite voxel volumes are caller-supplied; every positive zone must be mapped.
Outputs remain integrated strengths and shutdown time metadata stays 0.0.
Omitting volumes keeps the existing equal-count split; volumes require split.

Photon-source approximation: :func:`r2s_assemble` sums shutdown
``SpecificActivity`` over a zone and splits the total uniformly over the
requested energy groups. The uniform split preserves only the total shutdown
strength; real decay photons follow the nuclide- and energy-dependent lines
in ALARA ``.photonSrc`` spectra.

Transport handoff: :func:`photon_source` builds and validates an explicit
regular-mesh source at one cooling time (ascending eV edges, integrated
photons/s per voxel/group, C-order voxels, positive volumes). The returned
dict carries checked totals and a bounded structured report with stable issue
codes. :func:`mixed_voxel_photon_source` combines spectra from multiple
activation zones using caller-supplied intersection volumes.
:func:`photon_source_exported` restores sparse group exports with
explicit dense indices; :func:`source_probabilities` returns the joint
voxel/group distribution (never factored marginals),
:func:`source_energy_probabilities` returns per-voxel spectra and
:func:`source_density` applies each volume exactly once. :func:`select_cooling`
matches offsets exactly. OpenMC export and manifests live in
:mod:`nucleide.openmc`.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, TypedDict

from nucleide._internal import (
    r2s_assemble,
    r2s_expand,
    r2s_expand_sweep,
    r2s_from_deck,
    r2s_from_snapshot,
    r2s_mixed_voxel_photon_source,
    r2s_photon_group_sums,
    r2s_photon_source,
    r2s_photon_source_exported,
    r2s_select_cooling,
    r2s_snapshot_inventory,
    r2s_source_density,
    r2s_source_energy_probabilities,
    r2s_source_probabilities,
    r2s_source_report,
    r2s_tag_zone_strength,
    r2s_validate,
)

# Short aliases for the workflow-builder summaries.
from_deck = r2s_from_deck
from_snapshot = r2s_from_snapshot
validate = r2s_validate
expand = r2s_expand
assemble = r2s_assemble
tag_zone_strength = r2s_tag_zone_strength
photon_group_sums = r2s_photon_group_sums
snapshot_inventory = r2s_snapshot_inventory
expand_sweep = r2s_expand_sweep


def photon_source(
    mesh_shape: Sequence[int],
    lower_left_cm: Sequence[float],
    upper_right_cm: Sequence[float],
    voxel_volumes_cm3: Sequence[float],
    edges_ev: Sequence[float],
    cooling_s: float,
    strengths_photons_per_s: Sequence[Sequence[float]],
) -> dict[str, Any]:
    """Build and validate an integrated photon source at one cooling time.

    The regular mesh is ``(nx, ny, nz)`` with C-order voxel index
    ``(i*ny + j)*nz + k``; ``edges_ev`` are ascending group edges in eV and
    ``strengths_photons_per_s[v][g]`` are integrated photons/s per voxel and
    group. ``voxel_volumes_cm3`` is parallel to the voxel order and is used
    only when converting to a density (never multiplied into an integrated
    strength twice). The returned dict carries the validated fields, checked
    totals and the structured report, including the explicit modelling
    assumptions (uniform position within a voxel, uniform density per eV in a
    group). An all-zero source stays valid data but cannot be sampled.
    """
    return r2s_photon_source(
        list(mesh_shape),
        list(lower_left_cm),
        list(upper_right_cm),
        list(voxel_volumes_cm3),
        list(edges_ev),
        float(cooling_s),
        [list(row) for row in strengths_photons_per_s],
    )


class MixedVoxelZoneSpectrum(TypedDict):
    """One zone's declared volume, grid, cooling offset and integrated spectrum."""

    volume_cm3: float
    edges_ev: Sequence[float]
    cooling_s: float
    strengths_photons_per_s: Sequence[float]


def mixed_voxel_photon_source(
    mesh_shape: Sequence[int],
    lower_left_cm: Sequence[float],
    upper_right_cm: Sequence[float],
    voxel_volumes_cm3: Sequence[float],
    edges_ev: Sequence[float],
    cooling_s: float,
    zones: Sequence[MixedVoxelZoneSpectrum],
    overlaps: Sequence[tuple[int, int, float]],
) -> dict[str, Any]:
    """Allocate a photon source from the intersection volumes of disjoint zones.

    Each ``(voxel, zone, overlap_cm3)`` row has a unique pair of zero-based
    indices. Voxels are in flat C-order. Each zone's eV edges and cooling
    offset must match the selected source exactly. Emitting zones must have
    mapped volumes within a relative tolerance of 1e-12 of their declared
    volumes. The full voxel volumes must match the regular mesh geometry
    within the same tolerance. Reject overfill beyond that tolerance.

    Return the regular-mesh photon-source dict with an additional
    ``allocation`` section for zone coverage, residual voxel volumes and
    input and allocated group strengths. Unassigned volume has no source
    before homogenization. The exported distribution is uniform in position
    within each whole voxel and uniform per eV within each group.
    """
    return r2s_mixed_voxel_photon_source(
        list(mesh_shape),
        list(lower_left_cm),
        list(upper_right_cm),
        list(voxel_volumes_cm3),
        list(edges_ev),
        float(cooling_s),
        [
            (
                float(zone["volume_cm3"]),
                list(zone["edges_ev"]),
                float(zone["cooling_s"]),
                list(zone["strengths_photons_per_s"]),
            )
            for zone in zones
        ],
        list(overlaps),
    )


def photon_source_exported(
    mesh_shape: Sequence[int],
    lower_left_cm: Sequence[float],
    upper_right_cm: Sequence[float],
    voxel_volumes_cm3: Sequence[float],
    full_edges_ev: Sequence[float],
    cooling_s: float,
    exported_group_indices: Sequence[int],
    exported_strengths_photons_per_s: Sequence[Sequence[float]],
) -> dict[str, Any]:
    """Rebuild a dense source from exported groups plus explicit indices.

    ``exported_group_indices`` names the dense group positions of the
    exported columns and must be strictly increasing, unique and in range.
    Omitted groups are filled with exactly zero; centroids are never accepted
    as bin edges and missing/duplicate indices are loud errors.
    """
    return r2s_photon_source_exported(
        list(mesh_shape),
        list(lower_left_cm),
        list(upper_right_cm),
        list(voxel_volumes_cm3),
        list(full_edges_ev),
        float(cooling_s),
        list(exported_group_indices),
        [list(row) for row in exported_strengths_photons_per_s],
    )


def source_report(source: dict[str, Any]) -> dict[str, Any]:
    """Return the structured, bounded validation report for a source dict.

    Issue codes are stable, deterministic, and classified by severity
    (``error`` / ``warning`` / ``assumption``). The report is an inspection
    surface: constructors and exporters still raise on invalid input rather
    than silently repairing it.
    """
    return r2s_source_report(*_source_args(source))


def source_probabilities(source: dict[str, Any]) -> list[list[float]]:
    """Joint voxel/group sampling probabilities (never factor the marginals)."""
    return r2s_source_probabilities(*_source_args(source))


def source_energy_probabilities(source: dict[str, Any]) -> list[list[float]]:
    """Per-voxel energy probabilities, with all-zero rows for zero voxels."""
    return r2s_source_energy_probabilities(*_source_args(source))


def source_density(source: dict[str, Any]) -> list[list[float]]:
    """Integrated photons/(s cm3), applying each voxel volume exactly once."""
    return r2s_source_density(*_source_args(source))


def select_cooling(available_s: Sequence[float], requested_s: Sequence[float]) -> list[int]:
    """Exact cooling-offset selection: no interpolation or nearest substitution."""
    return r2s_select_cooling([float(t) for t in available_s], [float(t) for t in requested_s])


def _source_args(
    source: dict[str, Any],
) -> tuple[
    list[int],
    list[float],
    list[float],
    list[float],
    list[float],
    float,
    list[list[float]],
]:
    """Unpack a source dict into the shared binding's positional fields."""
    if not isinstance(source, dict):
        raise TypeError("photon source must be the dict returned by nucleide.r2s.photon_source")
    return (
        list(source["mesh_shape"]),
        list(source["lower_left_cm"]),
        list(source["upper_right_cm"]),
        list(source["voxel_volumes_cm3"]),
        list(source["edges_ev"]),
        float(source["cooling_s"]),
        [list(row) for row in source["strengths_photons_per_s"]],
    )


__all__ = [
    "r2s_from_deck",
    "r2s_from_snapshot",
    "r2s_validate",
    "r2s_expand",
    "r2s_assemble",
    "r2s_tag_zone_strength",
    "r2s_photon_group_sums",
    "r2s_snapshot_inventory",
    "r2s_expand_sweep",
    "from_deck",
    "from_snapshot",
    "validate",
    "expand",
    "assemble",
    "tag_zone_strength",
    "photon_group_sums",
    "snapshot_inventory",
    "expand_sweep",
    "photon_source",
    "MixedVoxelZoneSpectrum",
    "mixed_voxel_photon_source",
    "photon_source_exported",
    "source_report",
    "source_probabilities",
    "source_energy_probabilities",
    "source_density",
    "select_cooling",
]
