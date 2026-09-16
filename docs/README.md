# Nucleide Documentation

Welcome! Nucleide is a Rust toolkit for nuclear-engineering workflow glue,
exposed through a typed Python API: parse legacy code output (MCNP, Serpent,
FLUKA, ALARA, and friends), build materials, and run solvers for depletion,
point kinetics, enrichment cascades, spectroscopy, tritium transport, fusion
neutron sources, spectrum unfolding, and variance reduction — plus
damage/gas metrics (dpa, appm), clearance/waste-classification screening,
and uncertainty-quantification sampling.
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
  <line x1="190" y1="92" x2="190" y2="70" class="cm-edge" marker-end="url(#cm-arrow)"/>
  <line x1="530" y1="92" x2="530" y2="70" class="cm-edge" marker-end="url(#cm-arrow)"/>
  <rect x="40" y="96" width="300" height="130" rx="10" class="cm-box"/>
  <text x="190" y="124" text-anchor="middle" class="cm-t">Code I/O</text>
  <text x="190" y="148" text-anchor="middle" class="cm-i">MCNP &#xb7; Serpent &#xb7; FLUKA</text>
  <text x="190" y="168" text-anchor="middle" class="cm-i">ALARA &#xb7; CCCC &#xb7; MCPL</text>
  <text x="190" y="188" text-anchor="middle" class="cm-i">FISPACT-II &#xb7; ORIGEN</text>
  <text x="190" y="210" text-anchor="middle" class="cm-s">readers, writers, emitters, translation</text>
  <rect x="380" y="96" width="300" height="130" rx="10" class="cm-box"/>
  <text x="530" y="122" text-anchor="middle" class="cm-t">Solvers and analysis</text>
  <text x="530" y="144" text-anchor="middle" class="cm-i">Depletion (CRAM) &#xb7; Point kinetics</text>
  <text x="530" y="160" text-anchor="middle" class="cm-i">Enrichment &#xb7; Spectroscopy &#xb7; UQ sampling</text>
  <text x="530" y="176" text-anchor="middle" class="cm-i">Variance reduction &#xb7; Tritium transport</text>
  <text x="530" y="192" text-anchor="middle" class="cm-i">Fusion sources &#xb7; Spectrum unfolding</text>
  <text x="530" y="208" text-anchor="middle" class="cm-i">Damage metrics (dpa/appm) &#xb7; Clearance screening</text>
  <text x="530" y="221" text-anchor="middle" class="cm-s">plus the R2S workflow glue</text>
  <line x1="190" y1="252" x2="190" y2="232" class="cm-edge" marker-end="url(#cm-arrow)"/>
  <line x1="530" y1="252" x2="530" y2="232" class="cm-edge" marker-end="url(#cm-arrow)"/>
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
| [Enrichment cascade](tutorials/python/enrichment-cascade.md) | Set up and solve a multicomponent enrichment cascade, including standalone SWU and mass-ratio helpers |
| [Activation analysis](tutorials/python/activation-analysis.md) | Read ALARA, FISPACT-II, and ORIGEN files and assemble an R2S workflow, including ARMI snapshots |
| [Deterministic I/O](tutorials/python/deterministic-io.md) | Read ISOTXS and flux files and write PARTISN decks |
| [Emit code cards](tutorials/python/emit-cards.md) | Render one material to MCNP, Serpent, FLUKA, ALARA, and PARTISN cards, including ARMI blueprint keys |
| [Run kinetics](tutorials/python/run-kinetics.md) | Solve prescribed-reactivity point-kinetics transients |
| [Run spectroscopy](tutorials/python/run-spectroscopy.md) | Smooth spectra, count peaks, calibrate energy/efficiency, evaluate X-ray lines, read `.spe` files |
| [MCPL particle interchange](tutorials/python/mcpl-interchange.md) | Read and write MCPL particle lists; merge, subset, and summarize them; repair interrupted writes |
| [UQ sampling](tutorials/python/uq-sampling.md) | Draw seeded MVN, log-normal, and LHS samples over caller-supplied covariance blocks and perturb decay and fission-yield data |
| [VR and MAGIC](tutorials/python/vr-magic.md) | Derive weight-window lower bounds with MAGIC, emit them for OpenMC or Serpent, and sample birth voxels through alias tables |
| [Reaction names](tutorials/python/rxname.md) | Resolve reaction names/ids/MT numbers and walk the parent/daughter reaction graph |
| [Translate CSG](tutorials/python/translate-csg.md) | Translate scoped MCNP CSG decks to OpenMC, Serpent, PHITS, or GDML (Geant4) geometry, with a drift report listing every approximation the translator made |
| [Run tritium](tutorials/python/run-tritium.md) | Solve steady-state and transient 1D tritium permeation with traps and recombination boundaries |
| [Damage metrics](tutorials/python/damage-metrics.md) | Fold multigroup fluxes into dpa and gas production, and propagate uncertainty through the fold |
| [Fusion sources](tutorials/python/fusion-sources.md) | Build ring, point, and parametric tokamak neutron sources and emit MCNP/Serpent source cards |
| [Unfold a spectrum](tutorials/python/unfold-spectrum.md) | Recover a neutron spectrum from activation detector measurements with SAND-II, least-squares, or GRAVEL adjustment, including an IRDFF response pack |
| [Clearance screening](tutorials/python/clearance-screening.md) | Classify parsed activation inventories against EU or Spanish clearance tables |
| [Interactive tutorials](tutorials/interactive/index.mdx) | Run Nucleide in the browser through the WASM build |
| [Interactive — nuclides](tutorials/interactive/nuclides.mdx) | Nuclide identifiers and nuclear data |
| [Interactive — materials](tutorials/interactive/materials.mdx) | Formulas, fractions, mixing, and XML export |
| [Interactive — enrichment](tutorials/interactive/enrichment.mdx) | Solve MARC cascades live |
| [Interactive — depletion](tutorials/interactive/depletion.mdx) | One-step CRAM depletion |
| [Interactive — MCNP I/O](tutorials/interactive/mcnp-io.mdx) | Parse MCNP file snippets |
| [Interactive — Serpent I/O](tutorials/interactive/serpent-io.mdx) | Parse Serpent `_res.m`, `_dep.m`, and `_det.m` output |
| [Interactive — FLUKA I/O](tutorials/interactive/fluka-io.mdx) | Read USRBIN mesh tallies from `.lis` files |
| [Interactive — variance reduction](tutorials/interactive/variance-reduction.mdx) | MAGIC bounds and alias-table sampling |
| [Interactive — activation](tutorials/interactive/activation.mdx) | ALARA decks/outputs, FISPACT, R2S workflows |
| [Interactive — deterministic I/O](tutorials/interactive/deterministic.mdx) | ISOTXS fluxes and PARTISN decks |
| [Interactive — deck editor](tutorials/interactive/deck-editor.mdx) | Parse, validate, and edit full MCNP decks |
| [Interactive — code-card emission](tutorials/interactive/emitter.mdx) | Five-dialect card emission with a per-dialect mass-drift audit |
| [Interactive — point kinetics](tutorials/interactive/kinetics.mdx) | Step-reactivity transients and the prompt jump |
| [Interactive — tritium permeation](tutorials/interactive/tritium.mdx) | 1D breakthrough curve and time lag |
| [Interactive — spectroscopy](tutorials/interactive/spectroscopy.mdx) | Spectrum smoothing and peak counting |
| [Interactive — UQ sampling](tutorials/interactive/uq.mdx) | Seeded MVN sampling over caller-supplied covariance blocks (Python API and theory add log-normal + LHS) |
| [Interactive — MCPL particle lists](tutorials/interactive/mcpl-io.mdx) | Read and write MCPL particle lists |
| [Interactive — damage metrics](tutorials/interactive/damage.mdx) | Fold a multigroup flux into dpa and gas production |
| [Interactive — fusion sources](tutorials/interactive/fusion-sources.mdx) | Sample D-D/D-T ring sources and emit SDEF cards |
| [Interactive — spectrum unfolding](tutorials/interactive/unfold.mdx) | SAND-II, least-squares, and GRAVEL adjustment over detector measurements |
| [Interactive — clearance screening](tutorials/interactive/clearance.mdx) | Sum-of-fractions classification against EU or Spanish limits |

### Theory

| Document | Purpose |
| --- | --- |
| [Theory](theory/index.mdx) | Theory index and suggested reading order |
| [Depletion](theory/depletion.mdx) | Burnup matrices, the Bateman equation, CRAM, and time-series integrators |
| [Enrichment cascades](theory/enrichment.mdx) | MARC cascades, separation factors, and SWU |
| [Variance reduction](theory/variance-reduction.mdx) | MAGIC weight windows, OpenMC/Serpent weight-window emission, and alias-table source sampling |
| [Nuclear data](theory/nuclear-data.mdx) | Nuclide IDs, name dialects, masses, half-lives, and screening data |
| [Point kinetics](theory/kinetics.mdx) | Prescribed-reactivity PKE system, inhour relation, prompt jump, and the stiff-aware solver |
| [Tritium transport](theory/tritium.mdx) | 1D diffusion-trapping equations, McNabb–Foster traps, surface taxonomy, multi-layer series stacks with Sieverts or Henry interfaces, and the closed-form permeation checks |
| [Gamma-ray spectroscopy](theory/spectroscopy.mdx) | Smoothing, peak counting, energy/efficiency calibration, X-ray lines, and SPE readers |
| [UQ-lite sampling](theory/uq-sampling.mdx) | Seeded MVN, log-normal, and LHS sampling over caller-supplied covariance blocks, SANDY-compatible estimators, and the decay perturbation consumer |
| [Damage metrics](theory/damage-metrics.mdx) | NRT/arc-dpa displacement functions, the Lindhard partition, and spectral folding conventions |
| [Fusion neutron sources](theory/fusion-sources.mdx) | Ballabio fusion spectra, Bosch–Hale reactivity, the Miller-geometry parametric plasma map, and D-T/D-D fuel mixtures |
| [Neutron spectrum unfolding](theory/unfolding.mdx) | SAND-II, least-squares, and GRAVEL adjustment, convergence, and underdetermined systems |
| [Clearance screening](theory/clearance-screening.mdx) | Clearance index, the sum-of-fractions rule, and the EU and Spanish clearance tables |

### Reference

| Document | Purpose |
| --- | --- |
| [Reference](reference/index.md) | Reference index and quick links |
| [Crate overview](reference/crate-overview.mdx) | One-line responsibilities for every workspace crate |
| [Python API](reference/python-api.mdx) | Python facade overview and module map |
| [Sample data files](reference/fixtures.mdx) | Index of the sample files the test suite reads, with per-file source links |

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
