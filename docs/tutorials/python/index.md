---
title: Python tutorials
sidebar:
  order: 3
---

Hands-on Python guides for Nucleide. Each tutorial is short, self-contained,
and assumes you have already installed the project (see
[Getting started](../getting-started.md)).

## Suggested order

1. [Parse MCNP output](parse-mcnp-output.md) — read common MCNP output files.
2. [Parse Serpent output](parse-serpent-output.md) — read Serpent `_res.m`,
   `_dep.m`, and `_det.m` files.
3. [Parse FLUKA output](parse-fluka-output.md) — read FLUKA USRBIN `.lis`
   tally files.
4. [Build materials](build-materials.md) — build, mix, and serialize materials.
5. [Run depletion](run-depletion.md) — run a CRAM depletion solve.
6. [Enrichment cascade](enrichment-cascade.md) — solve a multicomponent
   enrichment cascade.
7. [Activation analysis](activation-analysis.md) — read ALARA, FISPACT-II,
   and ORIGEN files and assemble an R2S workflow.
8. [Deterministic I/O](deterministic-io.md) — read ISOTXS and flux files
   and write PARTISN decks.
9. [Emit code cards](emit-cards.md) — render one material to MCNP, Serpent,
   FLUKA, ALARA, and PARTISN cards with a mass-drift report.
10. [Run kinetics](run-kinetics.md) — solve prescribed-reactivity
    point-kinetics transients.
11. [Run spectroscopy](run-spectroscopy.md) — smooth spectra, count peaks,
    calibrate energy/efficiency, evaluate X-ray lines, and read `.spe` files.
12. [MCPL particle interchange](mcpl-interchange.md) — read and write MCPL
    particle lists.
13. [UQ sampling](uq-sampling.md) — draw seeded MVN samples over
    caller-supplied covariance blocks and perturb decay data.
14. [VR and MAGIC](vr-magic.md) — derive weight-window lower bounds with
    MAGIC and sample birth voxels through alias tables.
15. [Reaction names](rxname.md) — resolve reaction names/ids/MT numbers and
    walk the parent/daughter reaction graph.
