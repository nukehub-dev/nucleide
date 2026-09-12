# Nucleide Documentation

Welcome! Nucleide is a Rust toolkit for nuclear-engineering workflow glue,
exposed through a typed Python API: parse legacy code output (MCNP, Serpent,
FLUKA, ALARA, and friends), build materials, and run solvers for depletion,
point kinetics, enrichment cascades, spectroscopy, and variance reduction.
These docs are organized by audience, so you can jump straight to what you
need.

<svg viewBox="0 0 720 310" width="100%" role="img" font-family="sans-serif" aria-label="Nucleide capability map">
  <style>
    .cm-edge { stroke: var(--primary, #d97706); stroke-width: 2; fill: none; }
    .cm-head { fill: var(--primary, #d97706); }
    .cm-box { fill: none; stroke: currentColor; stroke-opacity: 0.35; stroke-width: 1.25; }
    .cm-bar { fill: var(--primary, #d97706); fill-opacity: 0.1;
      stroke: var(--primary, #d97706); stroke-opacity: 0.5; stroke-width: 1.25; }
    .cm-t { fill: currentColor; font-size: 14px; font-weight: 600; }
    .cm-i { fill: currentColor; fill-opacity: 0.75; font-size: 12px; }
    .cm-s { fill: currentColor; fill-opacity: 0.55; font-size: 11.5px; font-style: italic; }
  </style>
  <defs>
    <marker id="cm-arrow" markerWidth="9" markerHeight="9" refX="6.5" refY="3.5" orient="auto">
      <path d="M0,0 L7,3.5 L0,7 Z" class="cm-head"/>
    </marker>
  </defs>
  <rect x="110" y="20" width="500" height="44" rx="10" class="cm-bar"/>
  <text x="360" y="47" text-anchor="middle" class="cm-t">Python API (nucleide) &#xb7; WASM interactive tutorials</text>
  <line x1="185" y1="92" x2="185" y2="70" class="cm-edge" marker-end="url(#cm-arrow)"/>
  <line x1="475" y1="92" x2="475" y2="70" class="cm-edge" marker-end="url(#cm-arrow)"/>
  <rect x="60" y="96" width="260" height="130" rx="10" class="cm-box"/>
  <text x="190" y="124" text-anchor="middle" class="cm-t">Code I/O</text>
  <text x="190" y="148" text-anchor="middle" class="cm-i">MCNP &#xb7; Serpent &#xb7; FLUKA</text>
  <text x="190" y="168" text-anchor="middle" class="cm-i">ALARA &#xb7; CCCC</text>
  <text x="190" y="188" text-anchor="middle" class="cm-i">FISPACT-II &#xb7; ORIGEN</text>
  <text x="190" y="210" text-anchor="middle" class="cm-s">readers, writers, card emitters</text>
  <rect x="400" y="96" width="260" height="130" rx="10" class="cm-box"/>
  <text x="530" y="124" text-anchor="middle" class="cm-t">Solvers and analysis</text>
  <text x="530" y="148" text-anchor="middle" class="cm-i">Depletion (CRAM) &#xb7; Point kinetics</text>
  <text x="530" y="168" text-anchor="middle" class="cm-i">Enrichment cascades &#xb7; Spectroscopy</text>
  <text x="530" y="188" text-anchor="middle" class="cm-i">Variance reduction</text>
  <text x="530" y="210" text-anchor="middle" class="cm-s">plus the R2S workflow glue</text>
  <line x1="185" y1="252" x2="185" y2="232" class="cm-edge" marker-end="url(#cm-arrow)"/>
  <line x1="475" y1="252" x2="475" y2="232" class="cm-edge" marker-end="url(#cm-arrow)"/>
  <rect x="110" y="256" width="500" height="40" rx="10" class="cm-bar"/>
  <text x="360" y="281" text-anchor="middle" class="cm-t">Nuclear data (nuclei) &#xb7; Materials (material)</text>
</svg>

## Where to start

- **New here?** [Getting started](tutorials/getting-started.md) installs
  Nucleide from PyPI and runs your first snippet in minutes.
- **Parsing legacy code output?** The [tutorials](tutorials/index.md) walk
  through real files format by format.
- **Want the physics and math?** The [Theory](theory/index.mdx) pages explain
  the equations and assumptions behind each solver.
- **Contributing code?** Read [Local development](development/local-dev.md),
  then [Contributing](development/contributing.md).
- **Reviewing architecture?** Start with the
  [Architecture overview](architecture/overview.mdx).
- **Cutting a release?** See the [Roadmap](plan/roadmap.md) and the
  [changelog](../CHANGELOG.md).

The [project README](../README.md) covers the full feature list, project
status, and license.

## Full index

### Tutorials

| Document | Purpose |
| --- | --- |
| [Tutorials](tutorials/index.md) | Tutorial index and suggested reading order |
| [Getting started](tutorials/getting-started.md) | Install Nucleide from PyPI and run your first Python snippet |
| [Python tutorials](tutorials/python/index.md) | Python tutorial index and suggested reading order |
| [Parse MCNP output](tutorials/python/parse-mcnp-output.md) | Read xsdir, meshtal, MCTAL, WWINP, PTRAC, and SSW files, including the NumPy meshtal bridge |
| [Parse Serpent output](tutorials/python/parse-serpent-output.md) | Read Serpent `_res.m`, `_dep.m`, and `_det.m` output files |
| [Parse FLUKA output](tutorials/python/parse-fluka-output.md) | Read FLUKA USRBIN `.lis` tally files |
| [Build materials](tutorials/python/build-materials.md) | Build materials from formulae, mix compositions, export XML, and screen dose per gram |
| [Run depletion](tutorials/python/run-depletion.md) | Load a depletion chain and run single-step and multi-step CRAM/Bateman solves |
| [Enrichment cascade](tutorials/python/enrichment-cascade.md) | Set up and solve a multicomponent enrichment cascade |
| [Activation analysis](tutorials/python/activation-analysis.md) | Read ALARA, FISPACT-II, and ORIGEN files and assemble an R2S workflow, including ARMI snapshots |
| [Deterministic I/O](tutorials/python/deterministic-io.md) | Read ISOTXS and flux files and write PARTISN decks |
| [Emit code cards](tutorials/python/emit-cards.md) | Render one material to MCNP, Serpent, FLUKA, ALARA, and PARTISN cards, including ARMI blueprint keys |
| [Run kinetics](tutorials/python/run-kinetics.md) | Solve prescribed-reactivity point-kinetics transients |
| [Run spectroscopy](tutorials/python/run-spectroscopy.md) | Smooth spectra, count peaks, calibrate energy/efficiency, evaluate X-ray lines, read `.spe` files |
| [MCPL particle interchange](tutorials/python/mcpl-interchange.md) | Read and write MCPL particle lists |
| [UQ sampling](tutorials/python/uq-sampling.md) | Draw seeded MVN samples over caller-supplied covariance blocks and perturb decay data |
| [Interactive tutorials](tutorials/interactive/index.mdx) | Run Nucleide in the browser through the WASM build |
| [Interactive — nuclides](tutorials/interactive/nuclides.mdx) | Nuclide identifiers and nuclear data |
| [Interactive — materials](tutorials/interactive/materials.mdx) | Formulas, fractions, mixing, and XML export |
| [Interactive — enrichment](tutorials/interactive/enrichment.mdx) | Solve MARC cascades live |
| [Interactive — depletion](tutorials/interactive/depletion.mdx) | One-step CRAM depletion |
| [Interactive — MCNP I/O](tutorials/interactive/mcnp-io.mdx) | Parse MCNP file snippets |
| [Interactive — variance reduction](tutorials/interactive/variance-reduction.mdx) | MAGIC bounds and alias-table sampling |
| [Interactive — activation](tutorials/interactive/activation.mdx) | ALARA decks/outputs, FISPACT, R2S workflows |
| [Interactive — deterministic I/O](tutorials/interactive/deterministic.mdx) | ISOTXS fluxes and PARTISN decks |
| [Interactive — code-card emission](tutorials/interactive/emitter.mdx) | Five-dialect card emission and mass drift |
| [Interactive — point kinetics](tutorials/interactive/kinetics.mdx) | Step-reactivity transients and the prompt jump |
| [Interactive — spectroscopy](tutorials/interactive/spectroscopy.mdx) | Spectrum smoothing and peak counting |
| [Interactive — UQ sampling](tutorials/interactive/uq.mdx) | Seeded MVN sampling over caller-supplied covariance blocks |

### Theory

| Document | Purpose |
| --- | --- |
| [Theory](theory/index.mdx) | Theory index and suggested reading order |
| [Depletion](theory/depletion.mdx) | Burnup matrices, the Bateman equation, CRAM, and time-series integrators |
| [Enrichment cascades](theory/enrichment.mdx) | MARC cascades, separation factors, and SWU |
| [Variance reduction](theory/variance-reduction.mdx) | MAGIC weight windows and alias-table source sampling |
| [Nuclear data](theory/nuclear-data.mdx) | Nuclide IDs, name dialects, masses, half-lives, and screening data |
| [Point kinetics](theory/kinetics.mdx) | Prescribed-reactivity PKE system, inhour relation, prompt jump, and the stiff-aware solver |
| [Gamma-ray spectroscopy](theory/spectroscopy.mdx) | Smoothing, peak counting, energy/efficiency calibration, X-ray lines, and SPE readers |

### Reference

| Document | Purpose |
| --- | --- |
| [Reference](reference/index.md) | Reference index and quick links |
| [Crate overview](reference/crate-overview.mdx) | One-line responsibilities for every workspace crate |
| [Python API](reference/python-api.mdx) | Python facade overview and module map |
| [Fixtures](reference/fixtures.mdx) | Golden-byte test data index with per-file source links |

### Development

| Document | Purpose |
| --- | --- |
| [Local development](development/local-dev.md) | Toolchain, local build, test, and lint commands |
| [Contributing](development/contributing.md) | Branch workflow, commit style, and PR checklist |
| [Cross-code validation](development/validation.md) | Reproduce the PyNE/OpenMC validation harness and measured tables |

### Architecture

| Document | Purpose |
| --- | --- |
| [Architecture overview](architecture/overview.mdx) | High-level system overview and layer boundaries |
| [Crate responsibilities](architecture/crate-responsibilities.md) | Crate-level responsibilities and dependency rules |

### Plan

| Document | Purpose |
| --- | --- |
| [Roadmap](plan/roadmap.md) | Current status and upcoming priorities |
| [Decision log](plan/decision-log.md) | Architecture and process decisions with rationale |

See the project [AGENTS.md](../AGENTS.md) for ownership and contract details.
Rules for maintaining these docs live in
[Contributing](development/contributing.md#maintaining-the-docs).
