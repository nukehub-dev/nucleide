//! Activation-to-transport photon-source handoff for R2S.
//!
//! This module carries the activation-side results into a transport-ready
//! source without silently reinterpreting them:
//!
//! - [`normalize_flux`] converts a declared flux tally into a physical flux
//!   using an explicit [`TallyNormalization`] contract; per-source and
//!   source-rate-scaled variants never multiply the rate twice.
//! - [`PhotonSource`] is the dense structured-mesh photon source: voxel
//!   volumes, energy-bin edges, one cooling time, and integrated
//!   photons/s per voxel and group. [`PhotonSource::from_exported_groups`]
//!   expands a sparse group export back onto the full energy grid with
//!   zeros in omitted bins.
//! - [`select_cooling_indices`] selects cooling times by exact equality
//!   only; no interpolation or nearest-neighbour matching happens here.
//! - [`SourceReport`] validates a [`PhotonSourceDraft`] and explains every
//!   problem with a stable code, severity, location, message, and remedy.
//!   [`PhotonSource::from_draft`] rejects a draft exactly when the report
//!   contains an [`IssueSeverity::Error`].
//!
//! # Flat voxel order
//!
//! A voxel is identified by `(i, j, k)` with flat index
//! `(i * ny + j) * nz + k` (C order), where `mesh_shape = [nx, ny, nz]`.
//!
//! # Energy bins
//!
//! `edges_ev` holds bin edges, not centroids: `G + 1` strictly ascending
//! finite values whose first entry is `>= 0`. Centroids are never accepted
//! or converted into edges.
//!
//! # Issue codes
//!
//! Codes are stable, kebab-case, and reported in deterministic check order
//! (mesh, bounds, volumes, edges, cooling, strengths, totals, warnings,
//! assumptions). [`SOURCE_ISSUE_LIMIT`] caps the stored examples;
//! [`SourceReport::issue_count`] keeps the true count.
//!
//! Errors:
//!
//! - `source-mesh-degenerate`: a mesh dimension is zero.
//! - `source-mesh-overflow`: the mesh voxel product overflows `usize`.
//! - `source-bounds-non-finite`: a mesh-bound component is not finite.
//! - `source-bounds-not-ascending`: `upper_right_cm <= lower_left_cm` on an
//!   axis.
//! - `source-volume-count-mismatch`: volume count differs from the mesh
//!   voxel product.
//! - `source-volume-not-finite`: a voxel volume is not finite.
//! - `source-volume-not-positive`: a voxel volume is zero or negative.
//! - `source-edges-too-short`: fewer than two energy edges.
//! - `source-edge-not-finite`: an energy edge is not finite.
//! - `source-edges-not-ascending`: an energy edge does not exceed the
//!   previous edge.
//! - `source-edge-negative`: the first energy edge is negative.
//! - `source-cooling-not-finite`: the cooling time is not finite.
//! - `source-cooling-negative`: the cooling time is negative.
//! - `source-strength-count-mismatch`: strength-row count differs from
//!   the mesh voxel product.
//! - `source-group-count-mismatch`: a voxel row width differs from the
//!   energy-group count.
//! - `source-strength-not-finite`: a group strength is not finite.
//! - `source-strength-negative`: a group strength is negative.
//! - `source-total-not-finite`: compensated totals overflow.
//!
//! Warnings:
//!
//! - `source-all-zero`: every strength is zero, so sampling has no
//!   probability mass.
//!
//! Assumptions:
//!
//! - `source-assumption-uniform-within-bin`: emission has uniform density
//!   per eV inside each group (not uniform per group).
//! - `source-assumption-uniform-within-voxel`: emission is uniform in
//!   position inside each voxel.

use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};

/// Maximum number of [`SourceIssue`] examples stored in a [`SourceReport`].
pub const SOURCE_ISSUE_LIMIT: usize = 100;

/// How a flux tally is normalized before conversion to physical flux.
///
/// OpenMC fixed-source tallies are normalized by the number of source
/// particles **and** scaled by the sum of the configured source strengths
/// (`settings::run_mode == FIXED_SOURCE` in OpenMC 0.16's
/// `src/tallies/tally.cpp`). A source with strength `1.0` therefore yields a
/// per-source-particle tally, while a source configured in absolute
/// particles/s or photons/s yields a source-rate-scaled tally. The state
/// must be declared, because multiplying a rate-scaled tally by the rate a
/// second time is wrong.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TallyNormalization {
    /// Track-length / volume-integrated tally per source particle (OpenMC's
    /// default flux tally when the configured source strengths sum to 1).
    VolumeIntegratedPerSource,
    /// Volume-integrated tally already scaled by the configured total source
    /// strength (OpenMC source strengths sum to the physical rate).
    VolumeIntegratedRateScaled,
    /// Volume-averaged flux per source particle (e.g. random-ray
    /// volume-normalized tallies with normalized source weights).
    VolumeAveragedPerSource,
    /// Volume-averaged tally already scaled by the configured total source
    /// strength.
    VolumeAveragedRateScaled,
    /// Already physical flux; converting again is an error.
    Physical,
}

/// Convert a flux tally into physical flux using a declared normalization.
///
/// `tally[v][g]` is the tally in voxel `v` and group `g`; `volumes_cm3[v]`
/// is that voxel's volume in cm³. The tally must be nonempty and
/// rectangular, one volume is required per tally row, source rates must be
/// finite and `>= 0`, and every tally value must be finite and `>= 0`.
///
/// The per-source variants multiply by `source_rate_s`:
/// [`TallyNormalization::VolumeIntegratedPerSource`] computes
/// `value * source_rate_s / volume`; the volume-averaged variant computes
/// `value * source_rate_s`. The rate-scaled variants skip that multiplication
/// because OpenMC already applied the configured total source strength:
/// volume-integrated divides by volume only, volume-averaged returns the
/// values unchanged. [`TallyNormalization::Physical`] is rejected because
/// normalizing twice would silently corrupt the units. Shape, group order,
/// zeros, and caller inputs are preserved; non-finite intermediates or
/// results are errors.
pub fn normalize_flux(
    tally: &[Vec<f64>],
    source_rate_s: f64,
    volumes_cm3: &[f64],
    input: TallyNormalization,
) -> Result<Vec<Vec<f64>>> {
    if matches!(input, TallyNormalization::Physical) {
        return Err(Error::Invalid(
            "flux tally is already normalized; do not normalize twice".to_string(),
        ));
    }
    if !source_rate_s.is_finite() || source_rate_s < 0.0 {
        return Err(Error::Invalid(
            "source rate must be finite and nonnegative".to_string(),
        ));
    }
    let Some(first) = tally.first() else {
        return Err(Error::Invalid("flux tally is empty".to_string()));
    };
    if first.is_empty() {
        return Err(Error::Invalid(
            "flux tally rows must be nonempty".to_string(),
        ));
    }
    for (voxel, row) in tally.iter().enumerate() {
        if row.len() != first.len() {
            return Err(Error::Invalid(format!(
                "flux tally is ragged: voxel {voxel} has {} groups but voxel 0 has {}",
                row.len(),
                first.len()
            )));
        }
    }
    if volumes_cm3.len() != tally.len() {
        return Err(Error::Invalid(format!(
            "flux tally has {} voxel rows but {} volumes were supplied",
            tally.len(),
            volumes_cm3.len()
        )));
    }
    for (voxel, &volume) in volumes_cm3.iter().enumerate() {
        if !volume.is_finite() || volume <= 0.0 {
            return Err(Error::Invalid(format!(
                "voxel {voxel} volume must be finite and positive"
            )));
        }
    }
    for (voxel, row) in tally.iter().enumerate() {
        for (group, &value) in row.iter().enumerate() {
            if !value.is_finite() || value < 0.0 {
                return Err(Error::Invalid(format!(
                    "voxel {voxel} group {group} tally value must be finite and nonnegative"
                )));
            }
        }
    }

    let rate_scaled = matches!(
        input,
        TallyNormalization::VolumeIntegratedRateScaled
            | TallyNormalization::VolumeAveragedRateScaled
    );
    let mut flux = Vec::with_capacity(tally.len());
    for (voxel, row) in tally.iter().enumerate() {
        let mut converted = Vec::with_capacity(row.len());
        for (group, &value) in row.iter().enumerate() {
            let scaled = if rate_scaled {
                value
            } else {
                let scaled = value * source_rate_s;
                if !scaled.is_finite() {
                    return Err(Error::Invalid(format!(
                        "flux normalization overflowed at voxel {voxel} group {group}"
                    )));
                }
                scaled
            };
            let normalized = match input {
                TallyNormalization::VolumeIntegratedPerSource
                | TallyNormalization::VolumeIntegratedRateScaled => scaled / volumes_cm3[voxel],
                TallyNormalization::VolumeAveragedPerSource
                | TallyNormalization::VolumeAveragedRateScaled => scaled,
                TallyNormalization::Physical => unreachable!("physical input rejected above"),
            };
            if !normalized.is_finite() {
                return Err(Error::Invalid(format!(
                    "flux normalization overflowed at voxel {voxel} group {group}"
                )));
            }
            converted.push(normalized);
        }
        flux.push(converted);
    }
    Ok(flux)
}

/// Dense structured-mesh photon source handed from activation to transport.
///
/// Fields are public so callers can assemble and inspect the contract, but
/// every accessor revalidates the full layout through [`SourceReport`], so a
/// mutated source cannot silently produce totals or sampling probabilities
/// from inconsistent data.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PhotonSource {
    /// Mesh dimensions `[nx, ny, nz]`; every dimension is at least 1.
    /// The flat voxel index is `(i * ny + j) * nz + k` (C order).
    pub mesh_shape: [usize; 3],
    /// Lower-left mesh corner in cm (finite).
    pub lower_left_cm: [f64; 3],
    /// Upper-right mesh corner in cm (finite, strictly above the lower
    /// corner on every axis).
    pub upper_right_cm: [f64; 3],
    /// Voxel volumes in cm³, one per voxel in flat order (finite and
    /// positive).
    pub voxel_volumes_cm3: Vec<f64>,
    /// Energy-bin edges in eV: `G + 1` strictly ascending, finite values
    /// whose first entry is `>= 0`. These are edges, never centroids.
    pub edges_ev: Vec<f64>,
    /// Cooling time in seconds after shutdown (finite and `>= 0`).
    pub cooling_s: f64,
    /// Integrated group strengths in photons/s, indexed `[voxel][group]`
    /// (finite and `>= 0`, one row per voxel, `G` entries per row).
    pub strengths_photons_per_s: Vec<Vec<f64>>,
}

impl PhotonSource {
    /// Build and validate a dense photon source.
    ///
    /// Equivalent to [`Self::from_draft`] on a freshly built
    /// [`PhotonSourceDraft`]: the constructor and the report share one set
    /// of checks.
    pub fn new(
        mesh_shape: [usize; 3],
        lower_left_cm: [f64; 3],
        upper_right_cm: [f64; 3],
        voxel_volumes_cm3: Vec<f64>,
        edges_ev: Vec<f64>,
        cooling_s: f64,
        strengths_photons_per_s: Vec<Vec<f64>>,
    ) -> Result<Self> {
        let draft = PhotonSourceDraft::new(
            mesh_shape,
            lower_left_cm,
            upper_right_cm,
            voxel_volumes_cm3,
            edges_ev,
            cooling_s,
            strengths_photons_per_s,
        );
        Self::from_draft(&draft)
    }

    /// Build a dense photon source from an unvalidated draft.
    ///
    /// Returns [`Error::Invalid`] when the shared validation report contains
    /// any [`IssueSeverity::Error`]; warnings such as
    /// `source-all-zero` do not block construction.
    pub fn from_draft(draft: &PhotonSourceDraft) -> Result<Self> {
        let report = SourceReport::from_draft(draft);
        if !report.is_valid() {
            return Err(validation_error(&report));
        }
        if report.totals.is_none() {
            return Err(Error::Invalid(
                "photon source totals are not computable; check the mesh, edges, volumes, and \
                 strengths"
                    .to_string(),
            ));
        }
        Ok(Self {
            mesh_shape: draft.mesh_shape,
            lower_left_cm: draft.lower_left_cm,
            upper_right_cm: draft.upper_right_cm,
            voxel_volumes_cm3: draft.voxel_volumes_cm3.clone(),
            edges_ev: draft.edges_ev.clone(),
            cooling_s: draft.cooling_s,
            strengths_photons_per_s: draft.strengths_photons_per_s.clone(),
        })
    }

    /// Number of voxel rows in the strength layout.
    ///
    /// Validation requires this to equal `nx * ny * nz`.
    pub fn n_voxels(&self) -> usize {
        self.strengths_photons_per_s.len()
    }

    /// Number of energy groups (`edges_ev.len() - 1`, or zero for an
    /// invalid edge list).
    pub fn n_groups(&self) -> usize {
        self.edges_ev.len().saturating_sub(1)
    }

    /// Compensated total photon strength in photons/s.
    ///
    /// Uses Neumaier compensation and rejects non-finite cells or
    /// aggregates. An all-zero source returns `0.0` even though it cannot be
    /// sampled.
    pub fn checked_total_photons_per_s(&self) -> Result<f64> {
        Ok(self.checked_totals()?.total_photons_per_s)
    }

    /// Compensated per-voxel totals in photons/s in flat voxel order.
    pub fn voxel_totals(&self) -> Result<Vec<f64>> {
        Ok(self.checked_totals()?.voxel_totals_photons_per_s)
    }

    /// Compensated per-group totals in photons/s in energy-group order.
    pub fn group_totals(&self) -> Result<Vec<f64>> {
        Ok(self.checked_totals()?.group_totals_photons_per_s)
    }

    /// Photon density per voxel and group in photons/s/cm³.
    ///
    /// Each integrated strength is divided by its voxel volume exactly once.
    /// Zero-strength voxels stay zero; non-finite densities are errors.
    pub fn density_photons_per_s_cm3(&self) -> Result<Vec<Vec<f64>>> {
        self.checked_totals()?;
        let mut densities = Vec::with_capacity(self.strengths_photons_per_s.len());
        for (voxel, (row, &volume)) in self
            .strengths_photons_per_s
            .iter()
            .zip(&self.voxel_volumes_cm3)
            .enumerate()
        {
            let mut voxel_densities = Vec::with_capacity(row.len());
            for (group, &value) in row.iter().enumerate() {
                let density = value / volume;
                if !density.is_finite() {
                    return Err(Error::Invalid(format!(
                        "photon density at voxel {voxel} group {group} is not finite"
                    )));
                }
                voxel_densities.push(density);
            }
            densities.push(voxel_densities);
        }
        Ok(densities)
    }

    /// Joint voxel/group sampling probabilities (`strength / total`).
    ///
    /// Returns shape `[voxel][group]`; this is the joint distribution over
    /// voxels and groups. Do not factor it into independent spatial and
    /// energy marginals. A non-finite total or a total `<= 0` (including an
    /// all-zero source) is an error.
    pub fn joint_probabilities(&self) -> Result<Vec<Vec<f64>>> {
        let totals = self.checked_totals()?;
        let total = totals.total_photons_per_s;
        if !total.is_finite() || total <= 0.0 {
            return Err(Error::Invalid(format!(
                "photon source total {total} photons/s cannot be sampled; joint probabilities \
                 require a finite positive total"
            )));
        }
        let mut probabilities = Vec::with_capacity(self.strengths_photons_per_s.len());
        for (voxel, row) in self.strengths_photons_per_s.iter().enumerate() {
            let mut voxel_probabilities = Vec::with_capacity(row.len());
            for (group, &value) in row.iter().enumerate() {
                let probability = value / total;
                if !probability.is_finite() {
                    return Err(Error::Invalid(format!(
                        "joint probability at voxel {voxel} group {group} is not finite"
                    )));
                }
                voxel_probabilities.push(probability);
            }
            probabilities.push(voxel_probabilities);
        }
        Ok(probabilities)
    }

    /// Per-voxel energy probabilities (`strength / voxel total`).
    ///
    /// Returns shape `[voxel][group]`. Each positive-total row sums to one
    /// and describes only the energy distribution inside that voxel; a
    /// zero-source voxel yields an all-zero row (it carries no sampling
    /// distribution). Use [`Self::joint_probabilities`] for the joint
    /// voxel/group distribution. An all-zero source is not an error here
    /// because zero rows are meaningful data; callers that need a sampling
    /// distribution must reject a zero total themselves.
    pub fn energy_probabilities(&self) -> Result<Vec<Vec<f64>>> {
        let mut probabilities = Vec::with_capacity(self.strengths_photons_per_s.len());
        for (voxel, row) in self.strengths_photons_per_s.iter().enumerate() {
            let total: f64 = row.iter().sum();
            if !total.is_finite() {
                return Err(Error::Invalid(format!(
                    "voxel {voxel} total strength is not finite"
                )));
            }
            let mut voxel_probabilities = Vec::with_capacity(row.len());
            for (group, &value) in row.iter().enumerate() {
                let probability = if total > 0.0 { value / total } else { 0.0 };
                if !probability.is_finite() {
                    return Err(Error::Invalid(format!(
                        "energy probability at voxel {voxel} group {group} is not finite"
                    )));
                }
                voxel_probabilities.push(probability);
            }
            probabilities.push(voxel_probabilities);
        }
        Ok(probabilities)
    }

    /// Expand a sparse group export back onto the full energy grid.
    ///
    /// `full_edges_ev` describes the complete `G`-group grid;
    /// `exported_group_indices` lists the exported bins in strictly
    /// increasing, unique order; `exported_strengths[v][k]` is voxel `v`'s
    /// strength in exported bin `k`. Omitted groups are filled with zeros.
    ///
    /// Empty, duplicate, out-of-range, or unsorted indices are rejected,
    /// and every exported row must have exactly one value per index.
    /// `full_edges_ev` must be bin edges; centroids are never accepted or
    /// converted.
    #[allow(clippy::too_many_arguments)]
    pub fn from_exported_groups(
        full_edges_ev: &[f64],
        mesh_shape: [usize; 3],
        lower_left_cm: [f64; 3],
        upper_right_cm: [f64; 3],
        voxel_volumes_cm3: Vec<f64>,
        cooling_s: f64,
        exported_group_indices: &[usize],
        exported_strengths: &[Vec<f64>],
    ) -> Result<Self> {
        let full_groups = full_edges_ev
            .len()
            .checked_sub(1)
            .filter(|groups| *groups > 0)
            .ok_or_else(|| {
                Error::Invalid("full energy grid must contain at least two edges".to_string())
            })?;
        if exported_group_indices.is_empty() {
            return Err(Error::Invalid(
                "exported group indices must be nonempty".to_string(),
            ));
        }
        let mut previous = None;
        for &index in exported_group_indices {
            if index >= full_groups {
                return Err(Error::Invalid(format!(
                    "exported group index {index} is outside the {full_groups}-group full grid"
                )));
            }
            if let Some(previous) = previous {
                if index == previous {
                    return Err(Error::Invalid(format!(
                        "exported group indices contain duplicate {index}"
                    )));
                }
                if index < previous {
                    return Err(Error::Invalid(
                        "exported group indices must be strictly increasing".to_string(),
                    ));
                }
            }
            previous = Some(index);
        }
        let width = exported_group_indices.len();
        if let Some((voxel, row)) = exported_strengths
            .iter()
            .enumerate()
            .find(|(_, row)| row.len() != width)
        {
            return Err(Error::Invalid(format!(
                "exported strength row {voxel} has {} groups but {width} indices were supplied",
                row.len()
            )));
        }

        let mut dense = vec![vec![0.0; full_groups]; exported_strengths.len()];
        for (voxel, row) in exported_strengths.iter().enumerate() {
            for (column, &index) in exported_group_indices.iter().enumerate() {
                dense[voxel][index] = row[column];
            }
        }
        Self::new(
            mesh_shape,
            lower_left_cm,
            upper_right_cm,
            voxel_volumes_cm3,
            full_edges_ev.to_vec(),
            cooling_s,
            dense,
        )
    }

    fn checked_totals(&self) -> Result<SourceTotals> {
        let report = inspect_source(
            self.mesh_shape,
            self.lower_left_cm,
            self.upper_right_cm,
            &self.voxel_volumes_cm3,
            &self.edges_ev,
            self.cooling_s,
            &self.strengths_photons_per_s,
        );
        if !report.is_valid() {
            return Err(validation_error(&report));
        }
        report
            .totals
            .ok_or_else(|| Error::Invalid("photon source totals are not computable".to_string()))
    }
}

/// Unvalidated photon-source draft.
///
/// Carries the same fields as [`PhotonSource`] but imposes no invariants.
/// Run [`SourceReport::from_draft`] to obtain structured diagnostics, or
/// [`PhotonSource::from_draft`] to validate and construct.
#[derive(Debug, Clone, Default)]
pub struct PhotonSourceDraft {
    /// Mesh dimensions `[nx, ny, nz]`; may be degenerate.
    pub mesh_shape: [usize; 3],
    /// Lower-left mesh corner in cm; may be non-finite.
    pub lower_left_cm: [f64; 3],
    /// Upper-right mesh corner in cm; may be non-finite or inverted.
    pub upper_right_cm: [f64; 3],
    /// Candidate voxel volumes in cm³.
    pub voxel_volumes_cm3: Vec<f64>,
    /// Candidate energy-bin edges in eV.
    pub edges_ev: Vec<f64>,
    /// Candidate cooling time in seconds.
    pub cooling_s: f64,
    /// Candidate integrated group strengths in photons/s.
    pub strengths_photons_per_s: Vec<Vec<f64>>,
}

impl PhotonSourceDraft {
    /// Build an unvalidated draft from the photon-source fields.
    pub fn new(
        mesh_shape: [usize; 3],
        lower_left_cm: [f64; 3],
        upper_right_cm: [f64; 3],
        voxel_volumes_cm3: Vec<f64>,
        edges_ev: Vec<f64>,
        cooling_s: f64,
        strengths_photons_per_s: Vec<Vec<f64>>,
    ) -> Self {
        Self {
            mesh_shape,
            lower_left_cm,
            upper_right_cm,
            voxel_volumes_cm3,
            edges_ev,
            cooling_s,
            strengths_photons_per_s,
        }
    }
}

/// Severity of one [`SourceIssue`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IssueSeverity {
    /// The source cannot be used or constructed.
    Error,
    /// The source is constructible but a sampling limitation is worth
    /// surfacing (for example an all-zero total).
    Warning,
    /// A modelling assumption implied by the source representation.
    Assumption,
}

/// One structured diagnostic from [`SourceReport`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourceIssue {
    /// Stable kebab-case issue code (see the module-level list).
    pub code: &'static str,
    /// How the issue affects construction or sampling.
    pub severity: IssueSeverity,
    /// Deterministic location string such as `"mesh"`, `"voxel 3"`, or
    /// `"group 4"`.
    pub location: String,
    /// Human-readable description of the problem.
    pub message: String,
    /// Suggested corrective action; assumptions carry a "none" remedy.
    pub remedy: String,
}

/// Compensated totals computed for a source draft.
#[derive(Debug, Clone, PartialEq)]
pub struct SourceTotals {
    /// Compensated total over every voxel and group in photons/s.
    pub total_photons_per_s: f64,
    /// Compensated per-voxel totals in flat voxel order.
    pub voxel_totals_photons_per_s: Vec<f64>,
    /// Compensated per-group totals in energy-group order.
    pub group_totals_photons_per_s: Vec<f64>,
}

/// Structured validation report for a [`PhotonSourceDraft`].
#[derive(Debug, Clone, PartialEq)]
pub struct SourceReport {
    /// Stored diagnostics in deterministic check order, capped at
    /// [`SOURCE_ISSUE_LIMIT`].
    pub issues: Vec<SourceIssue>,
    /// Total diagnostics found before capping.
    pub issue_count: usize,
    /// Whether [`Self::issues`] omits diagnostics beyond the limit.
    pub truncated: bool,
    /// Compensated totals, or `None` when they are not computable.
    pub totals: Option<SourceTotals>,
}

impl SourceReport {
    /// Build a validation report for `draft`.
    ///
    /// The report never repairs the draft. Checks run in a fixed order:
    /// mesh, bounds, volumes, edges, cooling, strengths, totals, warnings,
    /// then assumptions.
    pub fn from_draft(draft: &PhotonSourceDraft) -> SourceReport {
        inspect_source(
            draft.mesh_shape,
            draft.lower_left_cm,
            draft.upper_right_cm,
            &draft.voxel_volumes_cm3,
            &draft.edges_ev,
            draft.cooling_s,
            &draft.strengths_photons_per_s,
        )
    }

    /// Whether the report contains no [`IssueSeverity::Error`].
    ///
    /// Warnings and assumptions do not make a source invalid.
    pub fn is_valid(&self) -> bool {
        !self
            .issues
            .iter()
            .any(|issue| issue.severity == IssueSeverity::Error)
    }
}

/// Select the index of each requested cooling time in `available_s`.
///
/// Matching is exact `==` only: no interpolation and no nearest-neighbour
/// snapping. Requested order is preserved, so repeated available times
/// resolve to the first occurrence. Available and requested times must be
/// finite and `>= 0`; duplicate requested offsets and missing offsets are
/// errors, and a missing-offset error names both the requested value and
/// the available times.
pub fn select_cooling_indices(available_s: &[f64], requested_s: &[f64]) -> Result<Vec<usize>> {
    for (index, &value) in available_s.iter().enumerate() {
        if !value.is_finite() || value < 0.0 {
            return Err(Error::Invalid(format!(
                "available cooling time {index} must be finite and nonnegative (found {value})"
            )));
        }
    }
    for (index, &value) in requested_s.iter().enumerate() {
        if !value.is_finite() || value < 0.0 {
            return Err(Error::Invalid(format!(
                "requested cooling time {index} must be finite and nonnegative (found {value})"
            )));
        }
    }
    for (index, &value) in requested_s.iter().enumerate() {
        if requested_s[..index].contains(&value) {
            return Err(Error::Invalid(format!(
                "requested cooling time {value} s appears more than once"
            )));
        }
    }

    let mut selected = Vec::with_capacity(requested_s.len());
    for &value in requested_s {
        let index = available_s
            .iter()
            .position(|&other| other == value)
            .ok_or_else(|| {
                let available = available_s
                    .iter()
                    .map(|time| time.to_string())
                    .collect::<Vec<_>>()
                    .join(", ");
                Error::Invalid(format!(
                    "requested cooling time {value} s is not among available times [{available}]"
                ))
            })?;
        selected.push(index);
    }
    Ok(selected)
}

/// Validate one draft and build its report.
///
/// This is the single interpretation shared by [`SourceReport::from_draft`]
/// and every [`PhotonSource`] accessor.
fn inspect_source(
    mesh_shape: [usize; 3],
    lower_left_cm: [f64; 3],
    upper_right_cm: [f64; 3],
    voxel_volumes_cm3: &[f64],
    edges_ev: &[f64],
    cooling_s: f64,
    strengths_photons_per_s: &[Vec<f64>],
) -> SourceReport {
    let mut issues = IssueBuilder::default();

    let mut mesh_ok = true;
    for (axis, &dimension) in mesh_shape.iter().enumerate() {
        if dimension == 0 {
            mesh_ok = false;
            issues.error(
                "source-mesh-degenerate",
                "mesh",
                format!("mesh dimension {axis} must be at least 1 (found {dimension})"),
                "give every mesh dimension at least one voxel",
            );
        }
    }
    let product = mesh_shape
        .iter()
        .try_fold(1usize, |product, &dimension| product.checked_mul(dimension));
    let product = match product {
        Some(product) => Some(product),
        None => {
            mesh_ok = false;
            issues.error(
                "source-mesh-overflow",
                "mesh",
                format!(
                    "mesh shape {}x{}x{} overflows the voxel count",
                    mesh_shape[0], mesh_shape[1], mesh_shape[2]
                ),
                "use a smaller mesh; the voxel count must fit in usize",
            );
            None
        }
    };

    let mut bounds_ok = true;
    for (axis, (&lower, &upper)) in lower_left_cm.iter().zip(upper_right_cm.iter()).enumerate() {
        if !lower.is_finite() {
            bounds_ok = false;
            issues.error(
                "source-bounds-non-finite",
                format!("lower_left_cm[{axis}]"),
                format!("lower mesh bound on axis {axis} is not finite"),
                "supply finite mesh bounds in cm",
            );
        }
        if !upper.is_finite() {
            bounds_ok = false;
            issues.error(
                "source-bounds-non-finite",
                format!("upper_right_cm[{axis}]"),
                format!("upper mesh bound on axis {axis} is not finite"),
                "supply finite mesh bounds in cm",
            );
        }
    }
    for (axis, (&lower, &upper)) in lower_left_cm.iter().zip(upper_right_cm.iter()).enumerate() {
        if lower.is_finite() && upper.is_finite() && upper <= lower {
            bounds_ok = false;
            issues.error(
                "source-bounds-not-ascending",
                format!("bounds[{axis}]"),
                format!("upper bound {upper} cm must exceed lower bound {lower} cm on axis {axis}"),
                "set upper_right_cm above lower_left_cm on every axis",
            );
        }
    }

    let mut volume_count_ok = false;
    if let Some(product) = product {
        if voxel_volumes_cm3.len() == product {
            volume_count_ok = true;
        } else {
            issues.error(
                "source-volume-count-mismatch",
                "voxel_volumes_cm3",
                format!(
                    "mesh has {product} voxels but {} volumes were supplied",
                    voxel_volumes_cm3.len()
                ),
                "supply one volume per mesh voxel",
            );
        }
    }
    let mut volumes_ok = volume_count_ok;
    for (voxel, &volume) in voxel_volumes_cm3.iter().enumerate() {
        if !volume.is_finite() {
            volumes_ok = false;
            issues.error(
                "source-volume-not-finite",
                format!("voxel {voxel}"),
                "voxel volume is not finite",
                "supply finite voxel volumes in cm3",
            );
        } else if volume <= 0.0 {
            volumes_ok = false;
            issues.error(
                "source-volume-not-positive",
                format!("voxel {voxel}"),
                format!("voxel volume {volume} cm3 must be positive"),
                "use positive voxel volumes; zero-volume voxels cannot carry density",
            );
        }
    }

    let edges_len_ok = edges_ev.len() >= 2;
    let mut edges_ok = edges_len_ok;
    if !edges_len_ok {
        issues.error(
            "source-edges-too-short",
            "edges_ev",
            format!(
                "energy grid has {} edges; at least 2 are required",
                edges_ev.len()
            ),
            "supply at least one energy bin (two edges) in eV",
        );
    }
    for (edge_index, &edge) in edges_ev.iter().enumerate() {
        if !edge.is_finite() {
            edges_ok = false;
            issues.error(
                "source-edge-not-finite",
                format!("edge {edge_index}"),
                format!("energy edge {edge_index} is not finite"),
                "supply finite energy edges in eV",
            );
        }
    }
    if edges_len_ok {
        for (edge_index, window) in edges_ev.windows(2).enumerate() {
            let (previous, edge) = (window[0], window[1]);
            let next = edge_index + 1;
            if previous.is_finite() && edge.is_finite() && edge <= previous {
                edges_ok = false;
                issues.error(
                    "source-edges-not-ascending",
                    format!("edge {next}"),
                    format!(
                        "energy edge {next} ({edge} eV) must exceed edge {edge_index} \
                         ({previous} eV)"
                    ),
                    "sort energy edges strictly ascending",
                );
            }
        }
        if edges_ev[0].is_finite() && edges_ev[0] < 0.0 {
            edges_ok = false;
            issues.error(
                "source-edge-negative",
                "edge 0",
                format!("first energy edge {} eV must be >= 0", edges_ev[0]),
                "start the energy grid at 0 eV or above",
            );
        }
    }

    if !cooling_s.is_finite() {
        issues.error(
            "source-cooling-not-finite",
            "cooling_s",
            "cooling time is not finite",
            "supply a finite cooling time in seconds",
        );
    } else if cooling_s < 0.0 {
        issues.error(
            "source-cooling-negative",
            "cooling_s",
            format!("cooling time {cooling_s} s must be >= 0"),
            "supply cooling_s >= 0 (seconds after shutdown)",
        );
    }

    let mut rows_ok = false;
    if let Some(product) = product {
        if strengths_photons_per_s.len() == product {
            rows_ok = true;
        } else {
            issues.error(
                "source-strength-count-mismatch",
                "strengths_photons_per_s",
                format!(
                    "mesh has {product} voxels but {} strength rows were supplied",
                    strengths_photons_per_s.len()
                ),
                "supply one strength row per mesh voxel",
            );
        }
    }

    let group_count = if edges_len_ok {
        Some(edges_ev.len() - 1)
    } else {
        None
    };
    let mut widths_ok = group_count.is_some();
    if let Some(group_count) = group_count {
        for (voxel, row) in strengths_photons_per_s.iter().enumerate() {
            if row.len() != group_count {
                widths_ok = false;
                issues.error(
                    "source-group-count-mismatch",
                    format!("voxel {voxel}"),
                    format!(
                        "voxel {voxel} has {} group strengths; expected {group_count}",
                        row.len()
                    ),
                    "pad or trim each voxel row to the energy-group count",
                );
            }
        }
    }

    let mut values_finite = true;
    let mut any_negative = false;
    for (voxel, row) in strengths_photons_per_s.iter().enumerate() {
        for (group, &value) in row.iter().enumerate() {
            if !value.is_finite() {
                values_finite = false;
                issues.error(
                    "source-strength-not-finite",
                    format!("voxel {voxel} group {group}"),
                    "strength is not finite",
                    "replace non-finite strengths with finite photons/s",
                );
            } else if value < 0.0 {
                any_negative = true;
                issues.error(
                    "source-strength-negative",
                    format!("voxel {voxel} group {group}"),
                    format!("strength {value} photons/s must be >= 0"),
                    "replace negative strengths with photons/s >= 0",
                );
            }
        }
    }

    let mut totals = None;
    if mesh_ok && rows_ok && widths_ok && values_finite {
        if let Some(group_count) = group_count {
            let mut voxel_totals = Vec::with_capacity(strengths_photons_per_s.len());
            let mut group_sums = vec![CompensatedSum::default(); group_count];
            let mut total_sum = CompensatedSum::default();
            let mut aggregates_finite = true;
            for row in strengths_photons_per_s {
                let mut voxel_sum = CompensatedSum::default();
                for (group, &value) in row.iter().enumerate() {
                    voxel_sum.add(value);
                    group_sums[group].add(value);
                    total_sum.add(value);
                }
                let voxel_total = voxel_sum.total();
                aggregates_finite &= voxel_total.is_finite();
                voxel_totals.push(voxel_total);
            }
            let group_totals: Vec<f64> = group_sums.iter().map(CompensatedSum::total).collect();
            let total = total_sum.total();
            aggregates_finite &=
                total.is_finite() && group_totals.iter().all(|value| value.is_finite());
            if aggregates_finite {
                if !any_negative && total == 0.0 {
                    issues.warning(
                        "source-all-zero",
                        "totals",
                        "every group strength is zero; sampling has no probability mass",
                        "provide a nonzero source or skip sampling",
                    );
                }
                totals = Some(SourceTotals {
                    total_photons_per_s: total,
                    voxel_totals_photons_per_s: voxel_totals,
                    group_totals_photons_per_s: group_totals,
                });
            } else {
                issues.error(
                    "source-total-not-finite",
                    "totals",
                    "compensated photon total is not finite",
                    "reduce the source magnitude; the aggregate must be finite",
                );
            }
        }
    }

    let geometry_ok = mesh_ok && bounds_ok && volume_count_ok && volumes_ok;
    if geometry_ok && edges_ok {
        issues.assumption(
            "source-assumption-uniform-within-bin",
            "edges_ev",
            "photon emission is uniform in energy within each group",
            "none; sampling treats each energy bin as uniform",
        );
        issues.assumption(
            "source-assumption-uniform-within-voxel",
            "mesh",
            "photon emission is uniform in position within each voxel",
            "none; sampling treats each voxel as uniform",
        );
    }

    issues.finish(totals)
}

/// Neumaier-compensated accumulator.
#[derive(Clone, Copy, Default)]
struct CompensatedSum {
    sum: f64,
    correction: f64,
}

impl CompensatedSum {
    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
    }

    fn total(&self) -> f64 {
        self.sum + self.correction
    }
}

/// Issue collector that caps stored examples but counts every issue.
#[derive(Default)]
struct IssueBuilder {
    issues: Vec<SourceIssue>,
    issue_count: usize,
}

impl IssueBuilder {
    fn push(&mut self, issue: SourceIssue) {
        self.issue_count += 1;
        if self.issues.len() < SOURCE_ISSUE_LIMIT {
            self.issues.push(issue);
        }
    }

    fn error(
        &mut self,
        code: &'static str,
        location: impl Into<String>,
        message: impl Into<String>,
        remedy: impl Into<String>,
    ) {
        self.push(SourceIssue {
            code,
            severity: IssueSeverity::Error,
            location: location.into(),
            message: message.into(),
            remedy: remedy.into(),
        });
    }

    fn warning(
        &mut self,
        code: &'static str,
        location: impl Into<String>,
        message: impl Into<String>,
        remedy: impl Into<String>,
    ) {
        self.push(SourceIssue {
            code,
            severity: IssueSeverity::Warning,
            location: location.into(),
            message: message.into(),
            remedy: remedy.into(),
        });
    }

    fn assumption(
        &mut self,
        code: &'static str,
        location: impl Into<String>,
        message: impl Into<String>,
        remedy: impl Into<String>,
    ) {
        self.push(SourceIssue {
            code,
            severity: IssueSeverity::Assumption,
            location: location.into(),
            message: message.into(),
            remedy: remedy.into(),
        });
    }

    fn finish(self, totals: Option<SourceTotals>) -> SourceReport {
        let truncated = self.issue_count > self.issues.len();
        SourceReport {
            issues: self.issues,
            issue_count: self.issue_count,
            truncated,
            totals,
        }
    }
}

/// First-error summary used by construction and accessor failures.
fn validation_error(report: &SourceReport) -> Error {
    match report
        .issues
        .iter()
        .find(|issue| issue.severity == IssueSeverity::Error)
    {
        Some(issue) => {
            let more = report.issue_count.saturating_sub(1);
            if more == 0 {
                Error::Invalid(format!(
                    "invalid photon source: {} at {}: {}",
                    issue.code, issue.location, issue.message
                ))
            } else {
                Error::Invalid(format!(
                    "invalid photon source: {} at {}: {} ({more} more issues)",
                    issue.code, issue.location, issue.message
                ))
            }
        }
        None => Error::Invalid(format!(
            "invalid photon source: {} issues reported but none stored",
            report.issue_count
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn edges_for(groups: usize) -> Vec<f64> {
        (0..=groups).map(|index| index as f64).collect()
    }

    fn valid_draft() -> PhotonSourceDraft {
        PhotonSourceDraft::new(
            [2, 3, 2],
            [0.0, 0.0, 0.0],
            [2.0, 3.0, 2.0],
            vec![1.0; 12],
            vec![0.0, 1.0, 2.0],
            0.0,
            vec![vec![1.0, 2.0]; 12],
        )
    }

    fn codes(report: &SourceReport) -> Vec<&'static str> {
        report.issues.iter().map(|issue| issue.code).collect()
    }

    // ---- normalization -------------------------------------------------

    #[test]
    fn normalization_hand_example_divides_by_volume() {
        let tally = vec![vec![4.0, 9.0]];
        let flux = normalize_flux(
            &tally,
            10.0,
            &[2.0],
            TallyNormalization::VolumeIntegratedPerSource,
        )
        .unwrap();
        assert_eq!(flux, vec![vec![20.0, 45.0]]);
    }

    #[test]
    fn normalization_uses_each_voxel_volume() {
        let tally = vec![vec![4.0, 9.0], vec![2.0, 3.0]];
        let flux = normalize_flux(
            &tally,
            10.0,
            &[2.0, 4.0],
            TallyNormalization::VolumeIntegratedPerSource,
        )
        .unwrap();
        assert_eq!(flux, vec![vec![20.0, 45.0], vec![5.0, 7.5]]);
    }

    #[test]
    fn normalization_zero_rate_preserves_shape_and_zeros() {
        let tally = vec![vec![4.0, 9.0, 0.0], vec![2.0, 3.0, 5.0]];
        let flux = normalize_flux(
            &tally,
            0.0,
            &[2.0, 4.0],
            TallyNormalization::VolumeIntegratedPerSource,
        )
        .unwrap();
        assert_eq!(flux, vec![vec![0.0; 3], vec![0.0; 3]]);
    }

    #[test]
    fn normalization_volume_averaged_branch_skips_division() {
        let tally = vec![vec![4.0, 9.0]];
        let flux = normalize_flux(
            &tally,
            10.0,
            &[2.0],
            TallyNormalization::VolumeAveragedPerSource,
        )
        .unwrap();
        assert_eq!(flux, vec![vec![40.0, 90.0]]);
    }

    #[test]
    fn normalization_rate_scaled_variants_skip_the_rate_multiplication() {
        let tally = vec![vec![4.0, 9.0]];
        let integrated = normalize_flux(
            &tally,
            10.0,
            &[2.0],
            TallyNormalization::VolumeIntegratedRateScaled,
        )
        .unwrap();
        assert_eq!(integrated, vec![vec![2.0, 4.5]]);
        let averaged = normalize_flux(
            &tally,
            10.0,
            &[2.0],
            TallyNormalization::VolumeAveragedRateScaled,
        )
        .unwrap();
        assert_eq!(averaged, vec![vec![4.0, 9.0]]);
    }

    #[test]
    fn normalization_rejects_physical_input() {
        let error = normalize_flux(&[vec![1.0]], 1.0, &[1.0], TallyNormalization::Physical)
            .unwrap_err()
            .to_string();
        assert!(error.contains("already normalized"));
        assert!(error.contains("do not normalize twice"));
    }

    #[test]
    fn normalization_rejects_bad_shapes_and_values() {
        let integrated = TallyNormalization::VolumeIntegratedPerSource;
        assert!(normalize_flux(&[], 1.0, &[], integrated).is_err());
        assert!(normalize_flux(&[vec![]], 1.0, &[1.0], integrated).is_err());
        assert!(
            normalize_flux(&[vec![1.0, 2.0], vec![3.0]], 1.0, &[1.0, 1.0], integrated).is_err()
        );
        assert!(normalize_flux(&[vec![1.0]], 1.0, &[], integrated).is_err());
        assert!(normalize_flux(&[vec![1.0]], 1.0, &[1.0, 2.0], integrated).is_err());
        for volume in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(normalize_flux(&[vec![1.0]], 1.0, &[volume], integrated).is_err());
        }
        for rate in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(normalize_flux(&[vec![1.0]], rate, &[1.0], integrated).is_err());
        }
        for value in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(normalize_flux(&[vec![value]], 1.0, &[1.0], integrated).is_err());
        }
    }

    #[test]
    fn normalization_detects_intermediate_and_result_overflow() {
        let integrated = TallyNormalization::VolumeIntegratedPerSource;
        assert!(normalize_flux(&[vec![f64::MAX]], 2.0, &[1.0], integrated).is_err());
        assert!(normalize_flux(&[vec![f64::MAX]], 1.0, &[f64::MIN_POSITIVE], integrated).is_err());
    }

    #[test]
    fn normalization_does_not_mutate_inputs() {
        let tally = vec![vec![4.0, 9.0], vec![2.0, 3.0]];
        let volumes = vec![2.0, 4.0];
        let tally_before = tally.clone();
        let volumes_before = volumes.clone();
        normalize_flux(
            &tally,
            10.0,
            &volumes,
            TallyNormalization::VolumeIntegratedPerSource,
        )
        .unwrap();
        assert_eq!(tally, tally_before);
        assert_eq!(volumes, volumes_before);
    }

    // ---- construction and validation ----------------------------------

    #[test]
    fn construction_accepts_a_dense_mesh_product() {
        let source = PhotonSource::from_draft(&valid_draft()).unwrap();
        assert_eq!(source.mesh_shape, [2, 3, 2]);
        assert_eq!(source.n_voxels(), 12);
        assert_eq!(source.n_groups(), 2);
        assert_eq!(source.checked_total_photons_per_s().unwrap(), 36.0);
        assert_eq!(source.voxel_totals().unwrap(), vec![3.0; 12]);
        assert_eq!(source.group_totals().unwrap(), vec![12.0, 24.0]);
    }

    #[test]
    fn construction_rejects_bad_mesh_and_geometry() {
        let good = valid_draft();

        let mut bad = good.clone();
        bad.mesh_shape = [0, 3, 2];
        assert!(PhotonSource::from_draft(&bad).is_err());

        let mut bad = good.clone();
        bad.voxel_volumes_cm3.pop();
        assert!(PhotonSource::from_draft(&bad).is_err());

        let mut bad = good.clone();
        bad.strengths_photons_per_s.pop();
        assert!(PhotonSource::from_draft(&bad).is_err());

        let mut bad = good.clone();
        bad.upper_right_cm[1] = bad.lower_left_cm[1];
        assert!(PhotonSource::from_draft(&bad).is_err());

        for volume in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let mut bad = good.clone();
            bad.voxel_volumes_cm3[2] = volume;
            assert!(PhotonSource::from_draft(&bad).is_err());
        }
        for edges in [
            vec![0.0],
            vec![0.0, 1.0, 1.0],
            vec![1.0, 0.0],
            vec![-1.0, 1.0],
            vec![0.0, f64::NAN],
        ] {
            let mut bad = good.clone();
            let groups = edges.len().saturating_sub(1);
            bad.edges_ev = edges;
            bad.strengths_photons_per_s = vec![vec![1.0; groups]; 12];
            assert!(PhotonSource::from_draft(&bad).is_err());
        }
        for cooling in [-1.0, f64::NAN, f64::INFINITY] {
            let mut bad = good.clone();
            bad.cooling_s = cooling;
            assert!(PhotonSource::from_draft(&bad).is_err());
        }

        let mut bad = good.clone();
        bad.strengths_photons_per_s[0].push(3.0);
        assert!(PhotonSource::from_draft(&bad).is_err());

        for value in [-1.0, f64::NAN, f64::INFINITY] {
            let mut bad = good.clone();
            bad.strengths_photons_per_s[3][1] = value;
            assert!(PhotonSource::from_draft(&bad).is_err());
        }
    }

    #[test]
    fn zero_source_voxel_is_allowed() {
        let source = PhotonSource::new(
            [1, 1, 2],
            [0.0; 3],
            [1.0, 1.0, 2.0],
            vec![1.0, 1.0],
            vec![0.0, 1.0, 2.0],
            0.0,
            vec![vec![0.0, 0.0], vec![3.0, 4.0]],
        )
        .unwrap();
        assert_eq!(source.voxel_totals().unwrap(), vec![0.0, 7.0]);
        assert_eq!(
            source.density_photons_per_s_cm3().unwrap(),
            vec![vec![0.0, 0.0], vec![3.0, 4.0]]
        );
    }

    #[test]
    fn all_zero_source_constructs_but_cannot_be_sampled() {
        let source = PhotonSource::new(
            [1, 1, 2],
            [0.0; 3],
            [1.0, 1.0, 2.0],
            vec![1.0, 1.0],
            vec![0.0, 1.0, 2.0],
            0.0,
            vec![vec![0.0, 0.0]; 2],
        )
        .unwrap();
        assert_eq!(source.checked_total_photons_per_s().unwrap(), 0.0);
        assert!(source.joint_probabilities().is_err());
    }

    #[test]
    fn joint_probabilities_are_joint_not_factored() {
        let source = PhotonSource::new(
            [1, 1, 2],
            [0.0; 3],
            [1.0, 1.0, 2.0],
            vec![2.0, 4.0],
            vec![0.0, 1.0, 2.0],
            0.0,
            vec![vec![100.0, 0.0], vec![0.0, 300.0]],
        )
        .unwrap();
        assert_eq!(
            source.joint_probabilities().unwrap(),
            vec![vec![0.25, 0.0], vec![0.0, 0.75]]
        );
    }

    #[test]
    fn energy_probabilities_are_per_voxel() {
        let source = PhotonSource::new(
            [1, 1, 2],
            [0.0; 3],
            [1.0, 1.0, 2.0],
            vec![1.0, 1.0],
            vec![0.0, 1.0, 2.0],
            0.0,
            vec![vec![100.0, 0.0], vec![0.0, 300.0]],
        )
        .unwrap();
        assert_eq!(
            source.energy_probabilities().unwrap(),
            vec![vec![1.0, 0.0], vec![0.0, 1.0]]
        );
        let zero_voxel = PhotonSource::new(
            [1, 1, 2],
            [0.0; 3],
            [1.0, 1.0, 2.0],
            vec![1.0, 1.0],
            vec![0.0, 1.0, 2.0],
            0.0,
            vec![vec![100.0, 0.0], vec![0.0, 0.0]],
        )
        .unwrap();
        assert_eq!(
            zero_voxel.energy_probabilities().unwrap(),
            vec![vec![1.0, 0.0], vec![0.0, 0.0]]
        );
    }

    #[test]
    fn density_divides_by_volume_exactly_once() {
        let source = PhotonSource::new(
            [1, 1, 2],
            [0.0; 3],
            [1.0, 1.0, 2.0],
            vec![2.0, 4.0],
            vec![0.0, 1.0, 2.0],
            0.0,
            vec![vec![100.0, 0.0], vec![0.0, 300.0]],
        )
        .unwrap();
        let density = source.density_photons_per_s_cm3().unwrap();
        assert_eq!(density, vec![vec![50.0, 0.0], vec![0.0, 75.0]]);
        for (voxel, row) in density.iter().enumerate() {
            for (group, value) in row.iter().enumerate() {
                assert_eq!(
                    value * source.voxel_volumes_cm3[voxel],
                    source.strengths_photons_per_s[voxel][group]
                );
            }
        }
    }

    #[test]
    fn irregular_edges_preserve_group_order() {
        let source = PhotonSource::new(
            [1, 1, 2],
            [0.0; 3],
            [1.0, 1.0, 2.0],
            vec![1.0, 1.0],
            vec![0.0, 0.5, 7.0, 20.0],
            0.0,
            vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]],
        )
        .unwrap();
        assert_eq!(source.n_groups(), 3);
        assert_eq!(source.group_totals().unwrap(), vec![5.0, 7.0, 9.0]);
        assert_eq!(source.checked_total_photons_per_s().unwrap(), 21.0);
    }

    #[test]
    fn checked_total_compensates_many_small_groups() {
        let small = f64::EPSILON / 4.0;
        let groups = 100_001usize;
        let mut strengths = vec![small; groups];
        strengths[0] = 1.0;
        let source = PhotonSource::new(
            [1, 1, 1],
            [0.0; 3],
            [1.0; 3],
            vec![1.0],
            edges_for(groups),
            0.0,
            vec![strengths],
        )
        .unwrap();
        let expected = 1.0 + 100_000.0 * small;
        assert_eq!(source.checked_total_photons_per_s().unwrap(), expected);
        // The ordinary floating-point sum loses the small terms.
        assert_eq!(source.strengths_photons_per_s[0].iter().sum::<f64>(), 1.0);
    }

    #[test]
    fn construction_rejects_non_finite_totals() {
        let draft = PhotonSourceDraft::new(
            [1, 1, 1],
            [0.0; 3],
            [1.0; 3],
            vec![1.0],
            vec![0.0, 1.0, 2.0],
            0.0,
            vec![vec![f64::MAX, f64::MAX]],
        );
        assert!(PhotonSource::from_draft(&draft).is_err());
        let report = SourceReport::from_draft(&draft);
        assert!(!report.is_valid());
        assert!(report.totals.is_none());
        assert!(codes(&report).contains(&"source-total-not-finite"));
    }

    // ---- sparse group export ------------------------------------------

    #[test]
    fn exported_groups_fill_omitted_bins_with_zeros() {
        let edges = edges_for(24);
        let exported = vec![vec![10.0, 20.0, 30.0], vec![40.0, 50.0, 60.0]];
        let source = PhotonSource::from_exported_groups(
            &edges,
            [1, 1, 2],
            [0.0; 3],
            [1.0, 1.0, 2.0],
            vec![1.0, 2.0],
            0.0,
            &[1, 5, 23],
            &exported,
        )
        .unwrap();
        assert_eq!(source.n_groups(), 24);
        assert_eq!(source.n_voxels(), 2);
        assert_eq!(source.strengths_photons_per_s[0][1], 10.0);
        assert_eq!(source.strengths_photons_per_s[0][5], 20.0);
        assert_eq!(source.strengths_photons_per_s[0][23], 30.0);
        assert_eq!(source.strengths_photons_per_s[1][1], 40.0);
        assert_eq!(source.strengths_photons_per_s[1][5], 50.0);
        assert_eq!(source.strengths_photons_per_s[1][23], 60.0);
        for group in 0..24 {
            if ![1, 5, 23].contains(&group) {
                assert_eq!(source.strengths_photons_per_s[0][group], 0.0);
                assert_eq!(source.strengths_photons_per_s[1][group], 0.0);
            }
        }
        assert_eq!(source.checked_total_photons_per_s().unwrap(), 210.0);
    }

    #[test]
    fn exported_groups_reject_ambiguous_indices() {
        let edges = edges_for(24);
        let exported = vec![vec![1.0, 2.0]];
        let call = |indices: &[usize], strengths: &[Vec<f64>]| {
            PhotonSource::from_exported_groups(
                &edges,
                [1, 1, 1],
                [0.0; 3],
                [1.0; 3],
                vec![1.0],
                0.0,
                indices,
                strengths,
            )
        };
        assert!(call(&[], &[]).is_err());
        assert!(call(&[1, 5, 5], &[vec![1.0, 2.0, 3.0]]).is_err());
        assert!(call(&[1, 24], &exported).is_err());
        assert!(call(&[5, 1], &exported).is_err());
        assert!(call(&[1, 5], &[vec![1.0]]).is_err());
        assert!(PhotonSource::from_exported_groups(
            &[0.0],
            [1, 1, 1],
            [0.0; 3],
            [1.0; 3],
            vec![1.0],
            0.0,
            &[0],
            &[vec![1.0]],
        )
        .is_err());
    }

    // ---- cooling-time selection ---------------------------------------

    #[test]
    fn cooling_selection_preserves_requested_order() {
        assert_eq!(
            select_cooling_indices(&[0.0, 10.0, 100.0], &[100.0, 0.0]).unwrap(),
            vec![2, 0]
        );
        assert_eq!(select_cooling_indices(&[0.0], &[0.0]).unwrap(), vec![0]);
        assert_eq!(
            select_cooling_indices(&[0.0, 10.0], &[]).unwrap(),
            Vec::<usize>::new()
        );
    }

    #[test]
    fn cooling_selection_is_exact_only() {
        let near = f64::MIN_POSITIVE;
        let error = select_cooling_indices(&[0.0], &[near])
            .unwrap_err()
            .to_string();
        assert!(error.contains(&near.to_string()));
        assert!(error.contains("available"));
    }

    #[test]
    fn cooling_selection_rejects_duplicates_and_invalid_times() {
        assert!(select_cooling_indices(&[0.0, 10.0], &[10.0, 10.0]).is_err());
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(select_cooling_indices(&[bad], &[0.0]).is_err());
            assert!(select_cooling_indices(&[0.0], &[bad]).is_err());
        }
        let error = select_cooling_indices(&[0.0, 10.0], &[2.5])
            .unwrap_err()
            .to_string();
        assert!(error.contains("2.5"));
        assert!(error.contains("0, 10"));
    }

    // ---- structured report --------------------------------------------

    #[test]
    fn report_valid_source_has_assumptions_and_totals() {
        let report = SourceReport::from_draft(&valid_draft());
        assert!(report.is_valid());
        assert!(!report.truncated);
        assert_eq!(report.issue_count, 2);
        assert_eq!(report.issues.len(), 2);
        assert_eq!(
            codes(&report),
            vec![
                "source-assumption-uniform-within-bin",
                "source-assumption-uniform-within-voxel"
            ]
        );
        assert!(report
            .issues
            .iter()
            .all(|issue| issue.severity == IssueSeverity::Assumption));
        let totals = report.totals.unwrap();
        assert_eq!(totals.total_photons_per_s, 36.0);
        assert_eq!(totals.voxel_totals_photons_per_s, vec![3.0; 12]);
        assert_eq!(totals.group_totals_photons_per_s, vec![12.0, 24.0]);
    }

    #[test]
    fn report_has_stable_codes_and_deterministic_order() {
        let draft = PhotonSourceDraft::new(
            [0, 2, 2],
            [0.0, 2.0, 0.0],
            [0.0, 1.0, 1.0],
            vec![1.0, 0.0, -2.0, f64::NAN],
            vec![1.0, 1.0, 0.0],
            -1.0,
            vec![vec![1.0], vec![f64::NAN]],
        );
        let report = SourceReport::from_draft(&draft);
        assert_eq!(
            codes(&report),
            vec![
                "source-mesh-degenerate",
                "source-bounds-not-ascending",
                "source-bounds-not-ascending",
                "source-volume-count-mismatch",
                "source-volume-not-positive",
                "source-volume-not-positive",
                "source-volume-not-finite",
                "source-edges-not-ascending",
                "source-edges-not-ascending",
                "source-cooling-negative",
                "source-strength-count-mismatch",
                "source-group-count-mismatch",
                "source-group-count-mismatch",
                "source-strength-not-finite",
            ]
        );
        assert_eq!(report.issue_count, 14);
        assert_eq!(report.issue_count, report.issues.len());
        assert!(!report.truncated);
        assert!(report.totals.is_none());
        assert!(!report.is_valid());

        let repeated = SourceReport::from_draft(&draft);
        assert_eq!(report.issues, repeated.issues);
        assert_eq!(report.issue_count, repeated.issue_count);
        assert_eq!(report.truncated, repeated.truncated);
    }

    #[test]
    fn report_flags_non_finite_bounds_and_edges() {
        let draft = PhotonSourceDraft::new(
            [1, 1, 1],
            [f64::NAN, 0.0, 0.0],
            [1.0, f64::INFINITY, 1.0],
            vec![1.0],
            vec![0.0, f64::NAN, 2.0],
            0.0,
            vec![vec![1.0, 2.0]],
        );
        let report = SourceReport::from_draft(&draft);
        let all = codes(&report);
        assert!(all.contains(&"source-bounds-non-finite"));
        assert!(all.contains(&"source-edge-not-finite"));
        assert!(!all.contains(&"source-edges-not-ascending"));
        assert!(!report.is_valid());

        let negative_start = PhotonSourceDraft::new(
            [1, 1, 1],
            [0.0; 3],
            [1.0; 3],
            vec![1.0],
            vec![-1.0, 1.0],
            0.0,
            vec![vec![1.0]],
        );
        let report = SourceReport::from_draft(&negative_start);
        assert!(codes(&report).contains(&"source-edge-negative"));
    }

    #[test]
    fn report_flags_mesh_overflow_and_short_edges() {
        let draft = PhotonSourceDraft::new(
            [usize::MAX, 2, 2],
            [0.0; 3],
            [1.0; 3],
            Vec::new(),
            vec![1.0],
            0.0,
            Vec::new(),
        );
        let report = SourceReport::from_draft(&draft);
        assert!(codes(&report).contains(&"source-mesh-overflow"));
        assert!(codes(&report).contains(&"source-edges-too-short"));
        assert!(!report.is_valid());
        assert!(report.totals.is_none());
    }

    #[test]
    fn report_notes_all_zero_sources() {
        let draft = PhotonSourceDraft::new(
            [1, 1, 1],
            [0.0; 3],
            [1.0; 3],
            vec![1.0],
            vec![0.0, 1.0],
            0.0,
            vec![vec![0.0]],
        );
        let report = SourceReport::from_draft(&draft);
        assert!(report.is_valid());
        assert!(report.issues.iter().any(|issue| {
            issue.code == "source-all-zero" && issue.severity == IssueSeverity::Warning
        }));
        assert_eq!(report.totals.unwrap().total_photons_per_s, 0.0);
        assert!(PhotonSource::from_draft(&draft).is_ok());
        assert!(PhotonSource::from_draft(&draft)
            .unwrap()
            .joint_probabilities()
            .is_err());
    }

    #[test]
    fn report_caps_stored_diagnostics_but_counts_every_issue() {
        let voxels = 1000;
        let draft = PhotonSourceDraft::new(
            [voxels, 1, 1],
            [0.0; 3],
            [1.0; 3],
            vec![-1.0; voxels],
            vec![0.0, 1.0],
            0.0,
            vec![vec![1.0]; voxels],
        );
        let report = SourceReport::from_draft(&draft);
        assert_eq!(report.issue_count, voxels);
        assert_eq!(report.issues.len(), SOURCE_ISSUE_LIMIT);
        assert!(report.truncated);
        assert!(!report.is_valid());
        assert!(report
            .issues
            .iter()
            .all(|issue| issue.code == "source-volume-not-positive"));
        assert_eq!(report.issues[0].location, "voxel 0");
        assert_eq!(
            report.issues[SOURCE_ISSUE_LIMIT - 1].location,
            format!("voxel {}", SOURCE_ISSUE_LIMIT - 1)
        );
        // Totals remain computable even though volumes are invalid.
        assert_eq!(report.totals.unwrap().total_photons_per_s, 1000.0);
    }

    #[test]
    fn report_and_constructor_share_one_interpretation() {
        let draft = PhotonSourceDraft::new(
            [1, 1, 1],
            [0.0; 3],
            [1.0; 3],
            vec![0.0],
            vec![0.0, 1.0],
            0.0,
            vec![vec![1.0]],
        );
        let report = SourceReport::from_draft(&draft);
        assert!(!report.is_valid());
        let from_draft = PhotonSource::from_draft(&draft).unwrap_err().to_string();
        let from_new = PhotonSource::new(
            [1, 1, 1],
            [0.0; 3],
            [1.0; 3],
            vec![0.0],
            vec![0.0, 1.0],
            0.0,
            vec![vec![1.0]],
        )
        .unwrap_err()
        .to_string();
        assert!(from_draft.contains("source-volume-not-positive"));
        assert_eq!(from_draft, from_new);
    }
}
