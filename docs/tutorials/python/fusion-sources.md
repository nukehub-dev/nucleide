---
title: Fusion sources
sidebar:
  order: 19
---

Nucleide's plasma-source capability builds tokamak fusion neutron sources over
the D-D (2.45 MeV) and D-T (14.1 MeV) reactions: a simple point or ring
source, and a parametric Miller-geometry plasma driven by caller-supplied
profiles. Ion-temperature broadening follows the closed-form fits of Ballabio
et al. (1998), and the parametric emission is weighted by the Bosch & Hale
(1992) thermonuclear reactivity. This tutorial covers the Python API; the
implementation lives in `crates/plasma-source`. For the physics — the
spectrum fits, the Miller map, and the profile family — see the
[Fusion neutron sources theory](../../theory/fusion-sources.mdx) page.

Units follow the transport-code card convention: lengths in centimetres,
energies in MeV, ion temperature in keV. Profile densities are caller data in
m⁻³ and reactivity is returned in m³/s.

## Spectrum moments and reactivity

`spectrum_moments(reaction, ion_temperature_kev)` returns the closed-form
Gaussian moments of a reaction at an ion temperature: the `nominal_mev` cold
line, the broadened `mean_mev`, and `sigma_mev` (0 for a monoenergetic
spectrum). `reactivity(reaction, ion_temperature_kev)` returns the
thermonuclear reactivity ⟨σv⟩ in m³/s — the weight that makes hot plasma
regions emit more neutrons:

```python
import nucleide.plasma_source as ps

m = ps.spectrum_moments("dt", 20.0)  # D-T at T_i = 20 keV
print(m["nominal_mev"], m["mean_mev"], m["sigma_mev"])
print(ps.spectrum_moments("dd", 0.0))  # T_i = 0: monoenergetic 2.4495 MeV line
print(ps.reactivity("dt", 10.0), ps.reactivity("dd", 10.0))  # <sigma v> [m^3/s]
```

```text
14.021 14.072874253353246 0.33772176292161377
{'reaction': 'dd', 'label': 'D-D', 'nominal_mev': 2.4495, 'mean_mev': 2.4495, 'sigma_mev': 0.0}
1.1361734886996056e-22 6.022478475561915e-25
```

## Sample a ring source

`particles(spec, n, seed)` draws `n` source particles from a spec dict and
returns per-field float64 NumPy arrays: positions `x`/`y`/`z` [cm], unit
direction cosines `u`/`v`/`w`, `energy` [MeV], and `weight`. The same seed
reproduces the same stream. `kind="ring"` takes `radius` and `height` [cm];
`kind="point"` takes a three-entry `position` [cm]; both take `reaction`
(`"dt"` or `"dd"`), `ion_temperature_kev` (0 for the monoenergetic nominal
line), and an optional `weight`:

```python
spec = {
    "kind": "ring",
    "radius": 300.0,  # cm
    "height": 25.0,  # cm
    "reaction": "dt",
    "ion_temperature_kev": 20.0,
}
out = ps.particles(spec, 8, seed=42)
print(out["x"])
print(out["z"])  # constant: every birth sits exactly on the ring plane
print(out["energy"])  # Gaussian around 14.07 MeV with sigma ~ 0.34 MeV
```

```text
[ 259.30701166   37.11268691   21.42195648 -130.13635523 -222.82374279
  129.39802418 -117.39284125   35.1029737 ]
[25. 25. 25. 25. 25. 25. 25. 25.]
[13.93259638 13.77455742 13.88885062 13.66562868 13.89902385 13.43674875
 13.75764472 14.30496707]
```

## Emit source cards

`emit_source_cards(spec, bins=21)` renders the same spec as MCNP `SDEF` and
Serpent `src` cards. It returns `sdef` and `serpent`, each with the card text
and a drift report — rows that say how much of the source each card accounts
for (`accounted`, `rel_drift`) and whether the row was machine-verified by
re-reading the emitted text (`reparsed`). The spectrum tail a tabulated card
cannot carry is reported, not silently dropped:

```python
import nucleide.mcnp as mcnp

cards = ps.emit_source_cards(spec, bins=7)  # fewer bins for a compact card
print(cards["spectrum"])
print(cards["sdef"]["card"])
print(cards["serpent"]["card"])

parsed = mcnp.parse_sdef(cards["sdef"]["card"])  # typed SDEF reader
print(parsed["rad"], parsed["pos"], parsed["erg"])
print(parsed["card"] == cards["sdef"]["card"])  # round-trip is byte-identical
for row in cards["sdef"]["drift"]:
    print(row["quantity"], row["accounted"], row["reparsed"])
```

```text
{'nominal_mev': 14.021, 'mean_mev': 14.072874253353246, 'sigma_mev': 0.33772176292161377, 'mono': False}
SDEF POS=0 0 25
     AXS=0 0 1
     RAD=D1
     ERG=D2
     WGT=1
     PAR=n
SI1 L 300 300
SP1 D 0 1
SI2 L 12.915 13.3009 13.6869 14.0729 14.4588 14.8448 15.2308
SP2 D 0.00210575 0.0411007 0.240616 0.432291 0.240616 0.0411007 0.00210575
src 1 pos 0 0 25
src 1 rad d1
src 1 erg d2
src 1 wgt 1
SI1 300 300
SP1 1
SI2 12.915 13.3009 13.6869 14.0729 14.4588 14.8448 15.2308
SP2 0.00210575 0.0411007 0.240616 0.432291 0.240616 0.0411007 0.00210575
D1 0 0 25 D2
True
emission probability 0.9999366279278736 True
spatial distribution 1.0 True
```

The ring itself is reproduced exactly (`SI1 L 300 300`); the Gaussian is
tabulated at bin centres over ±4σ and the report accounts for 0.99994 of the
probability mass. SDEF rows are re-read through `nucleide.mcnp.parse_sdef`
(`reparsed=True`); Serpent rows are checked analytically instead
(`reparsed=False`, since Nucleide has no Serpent source reader). Pass
`mcnp_version=6` in the spec for an MCNP 6 deck — for neutron sources the
card text is identical. To write particles to a file instead, project the
vectors with `nucleide.mcpl` (the crate returns data, never files).

## Parametric plasma

`kind="parametric"` models the plasma cross-section with Miller flux
surfaces and caller-supplied L/H/A-mode profiles of ion density [m⁻³] and ion
temperature [keV]. Birth positions weight the local neutron production
`f_fuel·n²·⟨σv⟩` times the flux-surface volume element, so the source follows
the profiles and the geometry. A flat profile is the peaking factor set to 0
with the pedestal and separatrix values equal to the centre value:

```python
parametric = {
    "kind": "parametric",
    "reaction": "dt",
    "major_radius": 620.0,  # cm
    "minor_radius": 200.0,  # cm
    "elongation": 1.7,
    "triangularity": 0.3,
    "shafranov_factor": 10.0,  # cm
    "mode": "L",
    "pedestal_radius": 100.0,  # required key; L mode does not use it
    "ion_density_centre": 1.0e20,  # m^-3
    "ion_density_peaking_factor": 0.0,  # 0 -> flat profile
    "ion_density_pedestal": 1.0e20,
    "ion_density_separatrix": 1.0e20,
    "ion_temperature_centre": 20.0,  # keV
    "ion_temperature_peaking_factor": 0.0,  # 0 -> flat profile
    "ion_temperature_beta": 1.0,
    "ion_temperature_pedestal": 20.0,
    "ion_temperature_separatrix": 20.0,
}
out = ps.particles(parametric, 100_000, seed=2024)
major = (out["x"] ** 2 + out["y"] ** 2) ** 0.5
print(major.min(), major.max())  # inside [R0 - a - esh, R0 + a + esh]
print(out["z"].min(), out["z"].max())  # inside +/- elongation * minor_radius
print(out["energy"].mean(), out["energy"].std())  # the T_i = 20 keV D-T line

cards = ps.emit_source_cards(parametric, bins=15)
print([row["quantity"] for row in cards["sdef"]["drift"]])
print(cards["spectrum"])
```

```text
420.0527008214557 819.8867815714732
-339.85606581391625 339.8906979122606
14.073671176722137 0.33789976086959445
['emission probability', 'spatial marginals', 'joint correlation']
{'nominal_mev': 14.021, 'mean_mev': 14.072874253353246, 'sigma_mev': 0.33772176292161377, 'mono': False}
```

With flat profiles every volume element emits equally, so the sampled bounds
fill the mapped flux-surface region. A product-form card cannot carry the
correlation between radius and height, so the parametric drift report adds a
"joint correlation" row quantifying that loss; `spectrum` reports the
magnetic-axis moments (the centre temperature), not the birth-weighted mean.

## Not yet supported

The parametric model covers equimolar D-T, pure D-D, and arbitrary D/T fuel
mixtures on the full torus. Toroidal sectors are not yet supported and raise
a clear error instead of guessing:

```python
mixture = dict(parametric, fuel={"D": 0.7, "T": 0.3})
print(ps.particles(mixture, 4, seed=0)["energy"][:2])

try:
    ps.particles(dict(parametric, rotation_angle=1.57), 4, seed=0)
except ValueError as e:
    print("sector:", e)
```

```text
[14.62170927 14.07200675]
sector: plasma-source: not yet supported: `rotation_angle` (sectors need a toroidal-angle distribution — outside the parametric model)
```

Malformed inputs raise a clear error at every entry point — for example a
zero ring radius (`ring radius must be > 0 cm`) or a two-entry point
position. Invalid `mcnp_version` values (anything other than 5 or 6) are
rejected the same way by `emit_source_cards`.

## See also

- [Fusion neutron sources theory](../../theory/fusion-sources.mdx) for the
  spectrum fits, the Miller map, and the profile equations.
- `tests/test_plasma_source.py` for replays of the moment checks, the card
  round trips, and the input-error cases.
