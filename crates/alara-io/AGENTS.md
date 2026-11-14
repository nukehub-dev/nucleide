# ALARA interoperability

## Purpose

Read ALARA inputs and outputs and provide explicit activation-workflow
transformations and radiological summaries over caller-supplied inputs.

## Ownership

- This NAD owns `crates/alara-io/`, including parser and regrouping contracts.
- Root NAD owns bindings, Python facades, capability inventories and release
  records. Bindings expose Rust results without implementing numerical logic.

## Local Contracts

- Preserve parser ordering, spelling, units and established file-format parity.
  Parsing and explicit cross-reference validation remain separate operations.
- Python deck summaries expose ordered `mat_loading` entries `{zone, mixture}`
  and explicit `volumes` entries `{zone, volume_cm3}`; absent blocks give `[]`.
  Preserve case and `void`, infer no volume, and resolve no external composition.
- `flux::regroup_flux_flat_energy` requires explicit finite, nonnegative,
  strictly ascending eV edges with identical endpoints and finite nonnegative
  group-integrated source flux. No inferred edges, tolerance merging,
  clipping, extrapolation or implicit conversion.
- Regrouping assumes constant differential flux per eV inside each source
  interval. Use a linear overlap sweep, check non-finite arithmetic, preserve
  identity-grid values, use compensated totals/target accumulation, and report
  zero-based split source intervals even
  when their flux is zero. Conserved total flux does not imply conserved
  reaction rates or recovered fine spectral structure. Document that
  representational underflow can prevent relative conservation.
- Evaluated material compositions, cross-section libraries and activation
  solves remain caller or external-solver responsibilities. Photon centroids
  are labels and cannot substitute for group boundaries.

## Work Guidance

- Keep numerical logic in Rust; add no external solver or runtime download
  dependency for regrouping. Committed numerical examples stay synthetic.
- Keep existing deck, flux-file and schedule APIs backward compatible.
- Update public docs, Python contracts and release records through their
  owners when Rust exports or numerical assumptions change.

## Verification

- Follow the root canonical checks and numeric/parser coverage requirement.
- Regrouping tests use independent analytic overlap and constant-density
  vectors, aligned coarsening, conservation and flat-model recovery, plus
  invalid grids, mismatched support, bad flux and overflow rejection.
- Flux conservation uses scale-aware relative tolerance `1e-12` and exact
  zero checks; identity tests preserve each input value's bits.

## Child NAD Index

None.
