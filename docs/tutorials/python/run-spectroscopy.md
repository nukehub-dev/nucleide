---
title: Run Spectroscopy
sidebar:
  order: 11
---

Nucleide spectroscopy smooths spectra, counts peaks, calibrates energy and
efficiency, evaluates X-ray lines, and reads `.spe` files with the
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

Efficiency coefficients are caller-supplied (there is no fitting routine);
`eff_fit` selects the `ln E` (1) or `1/E` (2) expansion. FWHM coefficients
are parsed from `.spe` headers but never evaluated.

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

dollar = read_dollar_spe("gv_format.spe")  # first line must be $SPEC_ID:
plain = read_spe("plain.spe")  # rejects the $SPEC_ID: magic
text_parsed = parse_dollar_spe(open("gv_format.spe").read())
print(dollar["counts"][:5], dollar["ebin"][:5], dollar["dead_time"])
```

Each reader returns a dict with `spec_name`, `start_chan_num`,
`num_channels`, `channels`, `counts`, `ebin`, the time/detector header, both
calibration fits, and `file_name`. Note the pinned quirks: `$MEAS_TIM:` is
live-then-real, dollar channel labels stay positional under a nonzero
`start_chan_num`, and `$ROI:`/`$PRESETS:`/`$ENER_FIT:` are ignored.

## See also

- [Gamma-ray spectroscopy theory](../../theory/spectroscopy.mdx) for E1–E8
  and the full quirk list.
- [Run kinetics](run-kinetics.md) for the companion counting-era toolkit.
- [Cross-code validation results](https://github.com/nukehub-dev/nucleide/blob/main/validation/results.md)
  for these routines, benchmarked against PyNE `spectanalysis`/`gammaspec`.
