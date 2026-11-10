---
title: Fusion first wall
sidebar:
  order: 22
---

This guide ties the first-wall pieces together end to end: a tokamak
neutron source, the wall geometry, weight windows, particle exchange,
activation parsing, damage, clearance screening, tritium scoping, and a
spectrum check. Nucleide provides the glue between steps; running
transport or activation codes stays outside Nucleide, so each handoff
names the file or table the next code takes over.

**Screening-level example.** Every number below is synthetic hand
numbers, not evaluated nuclear data, and no code output is shipped with
the example. Use the workflow shape, then rerun it with your own
profiles, fluxes, response data, and limit tables. Sections H0–H7 name
each step so you can start from any section.

The stack is synthetic: tungsten armor 2 mm thick, then CuCrZr, then
EUROFER97, under a D-T wall load of about 1–2 MW/m². The wall load is
the one caller-override slot for the whole guide: change
`WALL_LOAD_MW_M2` in H0 and the fluxes in H5a, and every downstream
damage number, rescale with it.

## H0 — Tokamak source and wall-load flux

`nucleide.plasma_source` builds the birth source: a parametric
Miller-geometry plasma (the `kind="parametric"` spec, with
caller-supplied density and temperature profiles) emitting D-T neutrons
near 14 MeV. The wall load converts to a neutron flux through the D-T
neutron energy: each neutron carries 14.1 MeV, so

```text
flux [n/cm²/s] = (wall load [W/m²] / 1e4) / (14.1 MeV in joules)
```

```python
import nucleide.plasma_source as ps

# Caller override slot: D-T wall load [MW/m²]. Valid span is 1-2.
WALL_LOAD_MW_M2 = 1.5  # synthetic default
JOULE_PER_MEV = 1.602176634e-13

# 14.1 MeV per D-T neutron -> total neutron flux [n/cm²/s].
flux_n_cm2_s = WALL_LOAD_MW_M2 * 1e6 / (14.1 * JOULE_PER_MEV) / 1e4

# Synthetic ITER-like profiles: every geometry and profile value is a
# hand number for the example, not machine data.
spec = {
    "kind": "parametric",
    "reaction": "dt",
    "major_radius": 620.0,  # cm, synthetic
    "minor_radius": 200.0,  # cm, synthetic
    "elongation": 1.7,
    "triangularity": 0.3,
    "shafranov_factor": 10.0,  # cm, synthetic
    "mode": "H",
    "pedestal_radius": 150.0,  # cm, synthetic
    "ion_density_centre": 1.0e20,  # m^-3, synthetic
    "ion_density_peaking_factor": 1.0,
    "ion_density_pedestal": 5.0e19,  # m^-3, synthetic
    "ion_density_separatrix": 2.0e19,  # m^-3, synthetic
    "ion_temperature_centre": 20.0,  # keV, synthetic
    "ion_temperature_peaking_factor": 1.5,
    "ion_temperature_beta": 1.5,
    "ion_temperature_pedestal": 3.0,  # keV, synthetic
    "ion_temperature_separatrix": 0.5,  # keV, synthetic
}
parts = ps.particles(spec, 2048, seed=7)  # fixed seed, synthetic
energies = parts["energy"]
print(f"wall load {WALL_LOAD_MW_M2} MW/m2 -> {flux_n_cm2_s:.5g} n/cm2/s")
print(f"{len(energies)} sampled, mean {float(energies.mean()):.4f} MeV")

cards = ps.emit_source_cards(spec, bins=9)
print([row["quantity"] for row in cards["sdef"]["drift"]])
print(cards["spectrum"])
```

```text
wall load 1.5 MW/m2 -> 6.6399e+13 n/cm2/s
2048 sampled, mean 14.0744 MeV
['emission probability', 'spatial marginals', 'joint correlation']
{'nominal_mev': 14.021, 'mean_mev': 14.072874253353246, 'sigma_mev': 0.33772176292161377, 'mono': False}
```

## H1 — Wall geometry and material cards

`nucleide.mcnp` translates MCNP constructive-solid-geometry decks to
OpenMC, Serpent, PHITS, or GDML geometry. The deck below is a
centimetre slab stack: tungsten 0.2 cm, CuCrZr 0.5 cm, EUROFER97 1.3 cm
(all three thicknesses synthetic). Translation returns the emitted text
plus a drift report — a list of every approximation the translator
made, such as expanding each rectangular box into planes. Assert on
structure, never on emitted bytes; the drift report is the record to
review.

`nucleide.emit` renders one composition to MCNP, Serpent, FLUKA,
ALARA, and PARTISN cards. The first-wall alloys are caller nuclide
lists, so the executable card check below uses a synthetic two-nuclide
plumbing mix instead — it proves the emitter runs without claiming to
be a first-wall material. A pure-tungsten list then shows the known
FLUKA gap: the FLUKA name table carries no tungsten isotopes, so the
five-dialect call fails for that mix. Production tungsten cards use
the elemental FLUKA material name on the caller side.

```python
import nucleide.emit as emit
import nucleide.mcnp as mcnp

# Synthetic slab stack: 0.2 cm W + 0.5 cm CuCrZr + 1.3 cm EUROFER97.
# ZAIDs and densities are synthetic hand numbers.
fw_deck = """first wall slab stack synthetic
W CuCrZr EUROFER97 thicknesses synthetic, centimetre units
1 1 -19.3 -1
2 2 -8.9 -2
3 3 -7.8 -3
4 0 1

1 rpp 0 10 0 10 0 0.2
2 rpp 0 10 0 10 0.2 0.7
3 rpp 0 10 0 10 0.7 2.0

m1 74182 0.7 74184 0.3
m2 29100 0.2 29000 0.8
m3 26100 0.6 24000 0.4
"""
openmc_xml, openmc_drift = mcnp.parse_csg_to_openmc(fw_deck)
print(f"OpenMC {len(openmc_xml)} chars, drift {[d['action'] for d in openmc_drift]}")
assert "macrobody-expansion" in [d["action"] for d in openmc_drift]

serpent_text, _ = mcnp.parse_csg_to_serpent(fw_deck)
print(f"Serpent {len(serpent_text)} chars")
assert "surf 1 cuboid 0 10 0 10 0 0.2" in serpent_text.splitlines()

# Synthetic plumbing mix only: proves the emitter runs. Not a first-wall
# material -- production alloy lists are caller data.
proxy_cards = emit.emit_cards({"U235": 5.0, "U238": 95.0}, "fw_proxy", density=19.1)
print(proxy_cards["MCNP"][:50])
proxy_rows = emit.emit_drift_table({"U235": 5.0, "U238": 95.0}, "fw_proxy", density=19.1)
print([(row["code"], row["rel_drift"]) for row in proxy_rows])

try:
    emit.emit_cards({"W182": 93.0, "W184": 7.0}, "tungsten", density=19.3)
except ValueError as exc:
    print("ValueError:", exc)
```

```text
OpenMC 1208 chars, drift ['macrobody-expansion', 'macrobody-expansion', 'macrobody-expansion']
Serpent 416 chars
m1 92235.80c -0.05 92238.80c -0.95

[('MCNP', 0.0), ('Serpent', 0.0), ('FLUKA', 0.0), ('ALARA', 0.0), ('PARTISN', 0.0)]
ValueError: compound card requested with no components
```

## H2 — Weight windows from a mesh tally

`nucleide.vr` derives MCNP weight-window lower bounds from a mesh tally
with the MAGIC algorithm, then reformats the windows for OpenMC
(`settings.xml` fragment) or Serpent (WWINP file). The tally below is
the committed sample mesh tally; in a real run it is the H1 transport
calculation's tally on the translated geometry.

```python
import nucleide.vr as vr

meshtal = mcnp.read_meshtal("fixtures/mcnp/meshtal/mcnp_meshtal_single_meshtal.txt")
tally = meshtal.tallies[4]

magic_out = vr.magic(tally)
print(magic_out.groups_per_ve, len(magic_out.lower_bounds_ww))

grouped = vr.magic_with(tally, selection="per_group", tolerance=0.5, null_value=0.0)
openmc_windows = vr.emit_openmc_weight_windows(tally, grouped, mesh_id=1, window_id=1)
print(openmc_windows["notes"][0])

serpent_windows = vr.emit_serpent_wwin(tally, grouped, name="ww1", file="windows.wwd")
print(serpent_windows["card"])
assert serpent_windows["card"] == 'wwin ww1 wf "windows.wwd" 2'
```

```text
1 45
upper-bounds-synthesized: upper = lower * 5 (MAGIC produces lower bounds only; OpenMC requires both)
wwin ww1 wf "windows.wwd" 2
```

## H3 — Particle exchange through MCPL

`nucleide.mcpl` reads and writes MCPL particle lists, the code-neutral
interchange format for particle tracks. The crest here projects the H0
sampled births into a list: each sampled particle becomes one MCPL
record (energy in MeV, position in cm, neutron PDG code 2112). Files
live in a temporary directory only, so the example leaves nothing
behind.

```python
import os
import tempfile

import nucleide.mcpl as mcpl

tmp_dir = tempfile.mkdtemp(prefix="nucleide-fw-")
# 64 sampled births -> 64 MCPL records (synthetic stream, seed fixed in H0).
header = {
    "srcname": "fw_synth",
    "comments": [],
    "has_userflags": False,
    "has_polarisation": True,
    "double_prec": False,
    "universal_pdgcode": None,
    "universal_weight": None,
    "blobs": [],
}
mcpl_particles = []
for idx in range(64):
    mcpl_particles.append(
        {
            "ekin": float(parts["energy"][idx]),
            "polarisation": [0.0, 0.0, 0.0],
            "position": [float(parts["x"][idx]), float(parts["y"][idx]), float(parts["z"][idx])],
            "direction": [float(parts["u"][idx]), float(parts["v"][idx]), float(parts["w"][idx])],
            "time": 0.0,
            "weight": float(parts["weight"][idx]),
            "pdgcode": 2112,
            "userflags": 0,
        }
    )
first_path = os.path.join(tmp_dir, "fw1.mcpl")
mcpl.write_mcpl(first_path, header, mcpl_particles)
mcpl_file = mcpl.read_mcpl(first_path)
print(mcpl_file.nparticles, mcpl_file.srcname)

first_stats = mcpl.mcpl_stats(first_path)
print(first_stats["nparticles"], f"{first_stats['ekin_mean']:.4f}", first_stats["pdg_counts"])

second_path = os.path.join(tmp_dir, "fw2.mcpl")
mcpl.write_mcpl(second_path, header, mcpl_particles[:32])
merged_path = os.path.join(tmp_dir, "fw_merged.mcpl")
print(mcpl.merge_mcpl([first_path, second_path], merged_path))

subset_path = os.path.join(tmp_dir, "fw_subset.mcpl")
print(mcpl.extract_mcpl(merged_path, subset_path, {"start": 1, "stop": 5}))
print(mcpl.mcpl_stats(subset_path)["nparticles"])
```

```text
64 fw_synth
64 14.1204 [(2112, 64)]
96
4
4
```

## H4 — Activation workflow from committed samples

`nucleide.r2s` builds the two-step shutdown workflow description from an
ALARA deck, and `nucleide.alara` / `nucleide.fispact` parse ALARA and
FISPACT-II output without running either code. The committed sample
deck irradiates for one year and cools to 100 years — the multi-year
activation tail in miniature — and the FISPACT-II sample carries the
clearance-bearing inventory table (all of its numbers are synthetic
hand values). Your run supplies your own deck, flux files, and activation
output; the same calls parse them.

```python
from pathlib import Path

import nucleide.alara as alara
import nucleide.fispact as fispact
import nucleide.r2s as r2s

alara_deck_text = Path("fixtures/alara/decks/sample2").read_text()
workflow = r2s.r2s_from_deck(alara_deck_text)
r2s.r2s_validate(workflow, alara_deck_text)
print(workflow["steps"], workflow["cooling_s"])
print(r2s.r2s_expand(alara_deck_text))

alara_output = Path("fixtures/alara/output/sample2.out").read_text()
source = r2s.r2s_assemble(alara_output, "sample2", "inner_zone", 4)
print(source["zone"], f"{source['total']:.6g}")

alara_rows = alara.alara_parse_output(alara_output, "sample2")
print(len(alara_rows))

clearance_text = Path("fixtures/fispact/clearance.out").read_text()
clearance_rows = fispact.fispact_parse_clearance(clearance_text)
print(len(clearance_rows), sorted({r["interval"] for r in clearance_rows}))
```

```text
[{'zone': 'inner_zone', 'flux': 'flux_1'}, {'zone': 'outer_zone', 'flux': 'flux_1'}] [86400.0, 8640000.0, 31557600.0, 3155760000.0]
[{'duration_s': 31557600.0, 'flux': 'flux_1', 'is_cooling': False}]
inner_zone 7.87567e+12
2838
7 [1, 2]
```

## H5a — Damage and gas for one full-power year

`nucleide.damage` folds a caller-supplied multigroup flux against
caller-supplied response cross sections into displacement damage in dpa
and gas production in appm (atomic parts per million). The fluxes below
derive from H0: the first-wall group split and the backplate
attenuation are synthetic hand fractions, and the response cross
sections are synthetic round numbers from your own condensation
pipeline in a real run. The exposure is one full-power year
(`FPY_SECONDS`). `fold_uq` perturbs the folded `[flux, response]`
vector with seeded draws so each number carries a spread; the seed
pins the stream.

```python
import nucleide.damage as damage

# Three-group grid [MeV]. Group fractions and attenuation are synthetic.
bounds = [0.0, 0.1, 1.0, 20.0]
fw_flux = [flux_n_cm2_s * frac for frac in (0.15, 0.25, 0.60)]
bp_flux = [fw * att for fw, att in zip(fw_flux, (0.20, 0.15, 0.10))]
dpa_xs = [100.0, 200.0, 50.0]  # synthetic NRT responses [barns]
arc_xs = [80.0, 160.0, 40.0]  # synthetic arc-corrected responses [barns]
he_xs = [0.5, 1.0, 0.25]  # synthetic He responses [barns]
seconds = damage.FPY_SECONDS  # one full-power year

fw_nrt = damage.nrt_dpa(fw_flux, dpa_xs, bounds, seconds)
fw_arc = damage.arc_dpa(fw_flux, arc_xs, bounds, seconds)
bp_nrt = damage.nrt_dpa(bp_flux, dpa_xs, bounds, seconds)
fw_appm = damage.gas_appm(fw_flux, he_xs, bounds, seconds)
fw_ratio = damage.he_dpa_ratio(fw_flux, he_xs, dpa_xs, bounds, seconds)
print(f"first wall 1 fpy: NRT {fw_nrt:.5g} dpa, arc {fw_arc:.5g} dpa, He {fw_appm:.5g} appm")
print(f"backplate 1 fpy: NRT {bp_nrt:.5g} dpa (synthetic attenuation)")
print(f"He/dpa {fw_ratio:.5g} appm/dpa")
assert fw_ratio == fw_appm / fw_nrt

mean = [0.0] * 6
cov = [[0.0004 if i == j else 0.0 for j in range(6)] for i in range(6)]
uq_out = damage.fold_uq("nrt_dpa", fw_flux, dpa_xs, bounds, seconds, mean, cov, 4000, 20260915, 5.0)
print(f"fold UQ mean {uq_out['mean']:.5g} std {uq_out['std']:.5g} passed={uq_out['passed']}")
assert uq_out["passed"]
```

```text
first wall 1 fpy: NRT 0.19906 dpa, arc 0.15925 dpa, He 995.31 appm
backplate 1 fpy: NRT 0.028288 dpa (synthetic attenuation)
He/dpa 5000 appm/dpa
fold UQ mean 0.19898 std 0.0035589 passed=True
```

## H5b — Clearance screening over the parsed inventory

`nucleide.alara` compares an inventory against per-nuclide limits with
the sum-of-fractions rule: each nuclide contributes `activity / limit`,
and a sum up to 1 classifies as satisfied. The inventory below is the
parsed FISPACT-II step from H4 (total activities in Bq); the limits are
synthetic toy values chosen to land the example on a 2.0 sum, not a
governing table. Screening arithmetic is not a compliance decision —
production screening uses the governing table for the material and
pathway with a complete inventory.

```python
step1_inv = {r["nuclide"]: r["activity_bq"] for r in clearance_rows if r["interval"] == 1}
# Synthetic toy limits [Bq] matching the total-activity basis above.
toy_limits = {"V-55": 1.0, "Fe-56": 1.0, "Co-60": 4.0e6, "Rb-86m": 8.0e9}
print(alara.alara_clearance_index(step1_inv, toy_limits))
sof = alara.alara_sum_of_fractions(step1_inv, toy_limits)
print(sof["sum"], sof["class"], sof["max_nuclide"], sof["max_fraction"])
assert sof["sum"] == 2.0 and sof["class"] == "exceeded"
```

```text
2.0
2.0 exceeded Co-60 1.0
```

## H6 — Tritium scoping

`nucleide.tritium` solves 1D diffusion-trapping permeation in a slab.
The fluxes below scope the tungsten arm: a 2 mm single slab plus a
two-layer stack (2 mm tungsten over 10 mm EUROFER97), all thicknesses,
diffusivities, and solubilities synthetic hand numbers, with Dirichlet
face concentrations 1.0 and 0.0. The steady flux and the permeation
time lag give an order of magnitude only — this is a scoping check, not
a breeding or TBR claim.

```python
import nucleide.tritium as tritium

t_lag = tritium.time_lag(2e-3, 1e-9)  # L=2 mm, D=1e-9 m²/s, both synthetic
print(f"time lag {t_lag:.5g} s")

ss_slab = tritium.steady(
    2e-3,
    32,
    1e-9,
    {"kind": "dirichlet", "value": 1.0},
    {"kind": "dirichlet", "value": 0.0},
)
print(f"slab flux {ss_slab['flux_right']:.5g} mol/m2/s")

layers = [
    {"thickness": 2e-3, "cells": 16, "D": 1e-9, "solubility": 0.1},  # W, synthetic
    {"thickness": 10e-3, "cells": 32, "D": 5e-10, "solubility": 0.1},  # EUROFER97, synthetic
]
ss_stack = tritium.steady_layers(
    layers, {"kind": "dirichlet", "value": 1.0}, {"kind": "dirichlet", "value": 0.0}
)
print(f"stack flux {ss_stack['flux_right']:.5g} mol/m2/s (order-of-magnitude only)")
assert ss_stack["flux_right"] > 0.0
```

```text
time lag 666.67 s
slab flux 5e-07 mol/m2/s
stack flux 4.5455e-08 mol/m2/s (order-of-magnitude only)
```

## H7 — Foil check with SAND-II

`nucleide.unfold` recovers a neutron spectrum from activation-detector
measurements: each detector responds over the energy groups, and the
adjustment moves a starting guess until folding it through the
responses reproduces the measured rates. The check below folds a
hand-built synthetic spectrum through a hand-built synthetic detector
set, then unfolds it. In a real run the responses are caller data —
for example the runtime IRDFF-II pack — and the rates come from your
foil measurements. The reduced CuCrZr layer from H1 has no copper foil
in its synthetic activation chain; bring your own foil set for a real
measurement.

```python
import math

import nucleide.unfold as unfold


def synthetic_response(n_det, n_groups, width):
    rows = []
    for i in range(n_det):
        center = i * (n_groups - 1) / max(n_det - 1, 1)
        rows.append([math.exp(-(((j - center) / width) ** 2)) + 1e-3 for j in range(n_groups)])
    return rows


def log_midpoints(n_groups, lo, hi):
    step = math.log(hi / lo) / n_groups
    return [math.exp(math.log(lo) + (j + 0.5) * step) for j in range(n_groups)]


# Synthetic truth, responses, and guess: every value hand-built.
midpoints = log_midpoints(6, 1e-6, 10.0)  # MeV
truth = [en * math.exp(-en / 2.53e-5) + 1e-6 / en for en in midpoints]
response = synthetic_response(6, 6, 0.9)  # one detector per group
rates = unfold.forward_fold(response, truth)
unfolded = unfold.sandii(
    response, rates, [3.0 * p for p in truth], tolerance=1e-12, max_iterations=10_000
)
worst = max(abs(a - b) / b for a, b in zip(unfolded["spectrum"], truth))
print(unfolded["iterations"], f"worst error {worst:.3e}")
print(unfolded["rate_factors"])
assert worst < 1e-9
```

```text
2 worst error 3.284e-16
[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
```

## Caller inputs

Every stage takes caller data; nothing below is evaluated data:

- H0 profiles and `WALL_LOAD_MW_M2` (today: 1.5 MW/m² synthetic; set 1–2).
- H1 alloy compositions and the MCNP deck (today: synthetic thicknesses
  and ZAIDs; FLUKA tungsten cards are caller-side elemental names).
- H2 mesh tally from your transport run.
- H3 MCPL particle lists, written only to temporary files.
- H4 activation deck, flux files, and activation output.
- H5a group fluxes, response cross sections, and uncertainty blocks.
- H5b inventory and governing clearance limits (today: synthetic toys).
- H6 thicknesses, diffusivities, solubilities, and face laws.
- H7 detector responses and measured rates.

Change the H0 wall load and rerun the chain: the H5a fluxes and damage
numbers rescale linearly with it, so you can bracket the 1–2 MW/m²
span without re-deriving anything.

## See also

- [Fusion sources](fusion-sources.md) for the source sampler and card
  emission in depth.
- [Translate CSG](translate-csg.md) for the translation scope and drift
  reports.
- [VR and MAGIC](vr-magic.md) for weight windows and source sampling.
- [Damage metrics](damage-metrics.md) for the dpa and gas folds.
- [Clearance screening](clearance-screening.md) for the limit tables.
- [Run tritium](run-tritium.md) for the interface laws and time scans.
- [Unfold a spectrum](unfold-spectrum.md) for the adjustment methods.
- `notebooks/first-wall-workflow.ipynb` for the same chain as a runnable
  workflow.
