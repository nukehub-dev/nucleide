---
title: Run Spectroscopy
sidebar:
  order: 11
---

Nucleide spectroscopy smooths spectra, counts peaks, calibrates energy and
efficiency, estimates activity from a known gamma line, evaluates X-ray
lines, and reads `.spe` files with the
`nucleide-spectroscopy` crate. This tutorial covers the Python API; the
implementation lives in `crates/spectroscopy`. For the equations and the
pinned file-format quirks, see the
[Gamma-ray spectroscopy theory](../../theory/spectroscopy.mdx) page.

## Smooth a spectrum

```python
from nucleide.spectroscopy import five_point_smooth, rect_smooth

counts = [2.0, 5.0, 1.0, 6.0, 3.0, 8.0, 4.0]
print(rect_smooth(counts, 5))  # E1: window mean, edges copied
print(five_point_smooth(counts))  # E2: low-statistics weights
```

The rectangular width `m` must be odd and at least 3. Both passes copy the
edge channels that have no full neighbourhood.

## Count a peak

```python
from nucleide.spectroscopy import calc_bg, gross_count, net_counts

channels = [float(c) for c in range(len(counts))]
bg = calc_bg(counts, channels, 2, 5, 1)  # E3: m == 1 only
gc = gross_count(counts, channels, 2, 5)  # E4: half-open, excludes c2
nc = net_counts(counts, channels, 2, 5, 1)  # E5: gc - bg
print(bg, gc, nc)
```

Bounds index `counts` positionally: a nonzero `start_chan_num` never shifts
the slice.

## Calibrate energy and efficiency

```python
from nucleide.spectroscopy import detector_efficiency, energy_bins

ebin = energy_bins([0.0, 1.0, 2.0], [1.5, 2.0, 0.5])  # E6: a0 + a1*ch + a2*ch^2
eff = detector_efficiency(1.0, [-2.81861504261204, -0.727352820018942], 1)  # E7, MeV
print(ebin, eff)
```

Supply efficiency coefficients directly or fit them with `fit_efficiency`.
The `eff_fit` argument selects the `ln E` (1) or `1/E` (2) expansion.
The `.spe` readers return the FWHM coefficients stored in the header.

## Estimate activity from a known gamma line

Use the background-subtracted full-energy peak area, its counting standard
uncertainty, and the recorded live time. Supply the absolute full-energy
peak efficiency for the same measurement geometry and the probability of
that line per decay as fractions.

```python
from math import log

from nucleide.spectroscopy import gamma_line_activity

absolute_efficiency = detector_efficiency(1.0, [log(0.10)], 1)
estimate = gamma_line_activity(120.0, 12.0, 60.0, absolute_efficiency, 0.50)
assert abs(estimate.activity_bq - 40.0) < 1e-12
assert abs(estimate.counting_uncertainty_bq - 4.0) < 1e-12
```

The calculation retains a negative peak area as a signed estimate.
The reported uncertainty covers peak counting, with live time,
efficiency and emission probability held fixed. Apply any required
measurement corrections to the inputs.

## Evaluate X-ray lines

```python
from nucleide.spectroscopy import xray_lines

atomic = {
    "k_shell_fluor": 0.9,
    "l_shell_fluor": 0.4,
    "prob": 0.8,
    "kb_to_ka": 0.2,
    "ka2_to_ka1": 0.5,
    "ka1_en_kev": 10.0,
    "ka2_en_kev": 20.0,
    "kb_en_kev": 30.0,
    "l_en_kev": 40.0,
}
print(xray_lines(atomic, k_conv=2.0, l_conv=3.0))  # E8: [(energy_kev, intensity)] * 4
```

No atomic table is vendored — every constant above is explicit. `None` (or
NaN) marks a conversion absent. Upstream exposes no combined function for
this routine (only a material method), so this explicit entry point is the
documented Nucleide surface.

## Read `.spe` files

```python
from nucleide.spectroscopy import parse_dollar_spe, parse_spe, read_dollar_spe, read_spe

dollar = read_dollar_spe("fixtures/spectroscopy/dollar_min.spe")  # first line must be $SPEC_ID:
plain = read_spe("fixtures/spectroscopy/plain_min.spe")  # rejects the $SPEC_ID: magic
text_parsed = parse_dollar_spe(open("fixtures/spectroscopy/dollar_min.spe").read())
print(dollar["counts"][:5], dollar["ebin"][:5], dollar["dead_time"])
```

Each reader returns a dict with `spec_name`, `start_chan_num`,
`num_channels`, `channels`, `counts`, `ebin`, the time/detector header, both
calibration fits, and `file_name`. Note the pinned quirks: `$MEAS_TIM:` is
live-then-real, dollar channel labels stay positional under a nonzero
`start_chan_num`, and `$ROI:`/`$PRESETS:`/`$ENER_FIT:` are ignored.

## Emit an SDEF decay source

```python
from nucleide.spectroscopy import parse_lines_tsv, read_decay_lines, sdef_decay_source

lines = read_decay_lines("fixtures/spectroscopy/decay_lines_sample.tsv")  # energy_MeV intensity rows, # comments ok
bins, card = sdef_decay_source(lines, particle="Photon")
print(card)
```

`read_decay_lines` (or `parse_lines_tsv` on text) reads the runtime
decay-lines TSV interchange into `(energy_MeV, intensity)` pairs — no decay
data is vendored. `sdef_decay_source` merges duplicates, sorts, normalizes
to probabilities (E9), and renders the MCNP point-source card (inline
`ERG=<E>` for one line, `ERG=D1` + `SI1 L`/`SP1 D` for several). Card
*syntax* is verified; MCNP sampling *semantics* stay yours.

## See also

- [Gamma-ray spectroscopy theory](../../theory/spectroscopy.mdx) for E1–E10
  and the full quirk list.
- [Run kinetics](run-kinetics.md) for the companion counting-era toolkit.
- [Cross-code validation results](https://github.com/nukehub-dev/nucleide/blob/main/validation/results.md)
  for these routines, benchmarked against PyNE `spectanalysis`/`gammaspec`.
