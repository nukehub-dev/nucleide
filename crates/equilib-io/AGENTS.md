# `crates/equilib-io` AGENTS.md

## Purpose

Equilibrium data layer: read equilibrium data, never solve it. A
classic-netCDF `wout` reader (CDF-1/CDF-2 only, pure classic parser —
no HDF5 in Rust) plus the VMEC `&INDATA` input-text grammar, yielding
flux-surface Jacobians (J1–J3) for volume weighting.

## Ownership

Owns `crates/equilib-io/src/` (`classic.rs` probe + header walk,
`wout.rs` field extraction + Jacobian helpers, `wall.rs` wall-load
accumulation kernel, `indata.rs` namelist grammar, `error.rs`), the
Python surface (`nucleide.equilib`, `equilib_*` in `_internal`),
`tests/test_equilib.py`, and the
`validation/equilib_vs_desc.py` oracle.

## Local Contracts

- Gate pins (verbatim, do not re-probe): accept magic `43 44 46 01`
  (CDF-1) / `43 44 46 02` (CDF-2); reject the HDF5 signature
  `89 48 44 46 0D 0A 1A 0A` (netCDF-4) and `43 44 46 05` (CDF-5)
  loudly via `Error::WrongVariant`, pointing at facade-side
  conversion. `ncdump -k` must report `classic`/`64-bit offset`
  (`ClassicVariant::ncdump_kind` returns exactly those spellings).
- Reader minimum for Jacobians: `nfp`, `ns`, `xm`, `xn`, `mn_mode`,
  `mn_mode_nyq`, `rmnc`, `zmns`, `lmns`, `gmnc` (dims
  `radius`/`mn_mode`/`mn_mode_nyq`); later-use `bmnc`, `bsubumnc`,
  `bsubvmnc`, `bsubsmns`, `currumnc`, `currvmnc` ride through in
  `Wout::fields`, plus optional scalars `mpol`/`ntor`/`phiedge`/
  `volume_p`. Fourier tables must sit on `(radius, mode)` in that
  order (`Error::UnexpectedDims` otherwise).
- Nyquist mode maps resolve from `xm_nyq`/`xn_nyq` when present, else
  from `xm`/`xn` when their length already matches `mn_mode_nyq`;
  otherwise `jacobian`/`jacobian_grid`/`mean_jacobian` fail with
  `MissingVariable("xm_nyq")` — never a guessed truncation.
- Jacobian kernel (J1): `√g(θ,ζ) = Σ c·cos(m·θ − n·ζ)` over the
  resolved Nyquist maps; `ζ` follows the file's toroidal convention
  (one field period per the documented `wout` layout — the caller
  scales for full-torus work). J2 grids `ntheta` points over
  `[0, 2π)` by `nzeta` over `[0, 2π/nfp)`; J3 reads the `(0,0)`
  coefficient (loud when absent). The free `fourier_jacobian` kernel
  stays caller-array-clean for the wall-load consumer.
- INDATA grammar per the pinned sections (STELLOPT Input Data Format
  and v8.47 table, DESC VMEC Inputs, simsopt Running VMEC, pyQSC
  `to_vmec`): power-series profiles only, `NCURR = 1` means the
  `I'(s)` profile (`ncurr_is_iprime`), fixed-boundary stance
  (`is_free_boundary`, absent `LFREEB` reads fixed), scalar indices
  only (`Error::SlicedIndex` on `:`).
- Wall-load kernel (`wall.rs`, W1–W3): closed-form per-cell accumulation
  `q[j][k] = Σ_i S·J·Δs` of caller birth-rate densities onto the caller
  wall surface in the same flux coordinates (radial-major voxels,
  one-field-period `ζ`), reusing the crate `Error` vocabulary (no new
  variants) and the free `fourier_jacobian` kernel for caller Jacobian
  voxels. No transport, no shadowing, no FEM/thermal, no CAD, no HDF5;
  negative densities/Jacobians are loud. Gates: axisymmetric analytic
  limit, cosine hand vector, discrete conservation, malformed-input
  vectors (in-crate + `tests/test_equilib.py` + `equilib_vs_desc.py`
  E9–E10, same `equilib` report).
- OUT until re-probed: DESC-saved files and VMEC++-era outputs (probe
  every file regardless), any solving/optimization/coil design,
  HDF5/netCDF-4 in Rust, vendored solver outputs (synthetic `wout`
  bytes are built in-test from the published var lists only).

## Work Guidance

- New readers add a module plus Python/tests/docs in the same change
  (bindings thin, no logic).
- Keep the layering: this crate depends on no workspace crate;
  bindings depend on it, never the reverse.
- The in-test classic writer (`wout.rs` `mod tests`) is the only
  emitter: hand-packed headers + payloads, never solver bytes.

## Verification

- `cargo test -p nucleide-equilib-io` (magic accept/reject vectors,
  synthetic round-trips CDF-1/CDF-2, Jacobian hand vectors, INDATA
  grammar vectors).
- `pytest tests/test_equilib.py` after `maturin develop`.
- `validation/equilib_vs_desc.py` runs inside `run_all.sh`
  (auto-discovered; report name `equilib` is registered in
  `render_results.py` `SECTION_ORDER`; DESC/simsopt legs SKIP loudly
  outside the container).

## Child NAD Index

None.
