---
title: 'Nucleide: A Rust toolkit for nuclear-engineering data and workflow glue, with Python and WebAssembly interfaces'
tags:
  - Rust
  - Python
  - WebAssembly
  - nuclear engineering
  - Monte Carlo
  - depletion
  - isotope enrichment
  - variance reduction
authors:
  - name: Ahnaf Tahmid Chowdhury
    orcid: 0000-0003-1070-5576
    affiliation: 1
affiliations:
  - name: NukeHub
    index: 1
date: 29 August 2026
bibliography: paper.bib
---

<!-- markdownlint-disable MD025 -->

# Summary

Nucleide is a toolkit for the nuclear-engineering workflow *around* particle
transport codes: reading and writing legacy code files (MCNP, Serpent, FLUKA),
canonical nuclide identification across a dozen naming dialects, embedded
nuclear reference data (AME2020 masses, IUPAC 2013 abundances, ENDF/B-VIII.0
half-lives), material construction, one-step burnup-matrix depletion with the
Chebyshev Rational Approximation Method (CRAM) [@pusa2010cram; @pusa2016cram],
multicomponent isotope enrichment cascades (MARC/SWU) [@wood1999marc], and
variance-reduction utilities (MAGIC weight windows [@cooper2001magic] and
alias-table mesh source sampling [@walker1977alias; @vose1991alias]).

The core is written in Rust as a composable workspace of nine crates. A thin
PyO3 layer exposes a typed Python API (wheels for Linux, macOS, and Windows via
PyPI), and a `wasm-bindgen` build powers interactive tutorials that run
entirely in the browser. Correctness is anchored by byte-exact golden fixtures,
a strict CI pipeline (formatting, linting, unit and Python tests, coverage,
end-to-end browser tests), and a cross-code validation harness whose results
are committed to the repository.

# Statement of need

Analysts who work with Monte Carlo transport codes spend much of their time on
the surrounding file formats and data conversions rather than on transport
itself. PyNE [@pyne2014] pioneered this "toolkit around transport codes" role,
but its C++/Cython/Fortran build chain (CMake, MOAB, generated nuclear-data
HDF5) makes installation and embedding difficult, and the project sees little
active maintenance. OpenMC [@romano2015openmc; @romano2021depletion] absorbed
some of these capabilities natively (depletion, weight windows) but does not
provide legacy-code I/O or enrichment analytics, and its transport-centric API
is not designed as an embeddable utility library.

Nucleide fills this gap with a memory-safe, dependency-light Rust core that
installs from PyPI in seconds (`pip install nucleide`), has no CMake or Fortran
toolchain, and — uniquely among comparable tools — runs in the browser through
WebAssembly, enabling zero-install interactive teaching materials. It is a
complement to PyNE and OpenMC, not a competitor: it deliberately ports their
well-validated algorithms and validates against them (see below), while
omitting transport itself.

# Software design

The Rust workspace enforces strict layering: capability crates (`nuclei`,
`material`, `mcnp-io`, `serpent-io`, `fluka-io`, `vr-tools`, `enrichment`,
`depletion`, `linalg`) never depend on the bindings; `bindings/python` and
`bindings/wasm` are thin facades with no business logic; the pure-Python
package re-exports the compiled module behind `.pyi` stubs so the public API is
fully typed and `mypy --strict` clean. Parsers reproduce legacy output
byte-for-byte where the reference codes have formatting quirks, guarded by
golden-byte fixtures.

Two design choices improve on the reference implementations. The depletion
crate implements CRAM in incomplete-partial-fraction product form with a sparse
LU factorization whose symbolic pattern is computed once and reused across all
poles, and it validates inputs that crash or silently corrupt results in the
reference codes (invalid half-lives, duplicate reaction entries, malformed
tallies). The enrichment solver adds a golden-section polish to the classic
sign-tracking descent for the optimal mass separation factor $M^*$
[@wood1999marc; @zeng2014cascade].

# Validation and performance

The repository contains a runnable cross-code validation harness
(`validation/`, results committed in `validation/results.md`) comparing
Nucleide 0.1.0 against PyNE 0.7.5 (numerically identical to the 0.7.8 release
for the exercised modules) and OpenMC 0.16.0:

- **Depletion**: CRAM-48 on a realistic nickel activation chain agrees with
  OpenMC's CRAM-48 solver to a maximum relative difference of $8.3\times10^{-15}$;
  a three-nuclide chain matches the closed-form Bateman solution to
  $\sim10^{-15}$.
- **Enrichment**: the $M^*$-optimizing multicomponent solver agrees with PyNE's
  `multicomponent()` to $\sim10^{-4}$ or better in stage counts, $M^*$, and
  separative work for uranium and tungsten feeds; the small residuals stem from
  the AME2020-vs-AME2016 mass tables and the golden-section polish.
- **MAGIC weight windows**: output matches the reference formula exactly on a
  shared test tally (PyNE's MOAB-dependent path was unavailable, so a
  formula-equivalent reference was used).
- **Nuclear data**: natural abundances and half-lives match OpenMC exactly
  (both derive from IUPAC 2013 [@meija2016iupac] and ENDF/B-VIII.0
  [@brown2018endf]); masses match OpenMC's AME2020 [@huang2021ame2020;
  @wang2021ame2020] tables exactly and PyNE's AME2016 tables to
  $4\times10^{-6}$ u; all name-dialect conversions match `pyne.nucname`
  exactly.
- **Performance**: single-step CRAM-48 solves run in $\sim$147 µs from Python
  ($\sim$76 µs native) versus $\sim$3.0 ms for OpenMC's Python path; the
  default uranium enrichment solve runs in $\sim$109 µs versus $\sim$5.5 ms
  for PyNE; MAGIC weight-window generation runs in $\sim$0.8 µs versus
  $\sim$4.5 µs for an equivalent pure-Python implementation.

# Documentation

The documentation website (built with Astro, deployed to GitHub Pages) provides
tutorials, an API reference, and theory pages deriving the implemented
mathematics, plus six interactive browser tutorials powered by the WebAssembly
build that let users run depletion, enrichment, MAGIC, and file-parsing
examples with no installation.

# Availability

Nucleide is BSD-2-Clause-licensed and developed at
<https://github.com/nukehub-dev/nucleide>. Python wheels are published to PyPI;
the Rust crates can be published to crates.io from the same release workflow.

# Acknowledgements

Nucleide's algorithms are ports of work by the PyNE and OpenMC communities; the
author thanks both projects for their openly available code and documentation.

# References
