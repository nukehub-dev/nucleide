# fixtures/spectroscopy/ — synthetic gamma-spectroscopy inputs and oracles

Hand-built synthetic data only — no laboratory or evaluated-library values.
Counts are small round numbers chosen for this repo; expected values are
closed-form consequences of the E1–E5 equations, computed once in float64
with the provenance recorded in each file.

- `smooth_oracle.json` — 7-channel counts plus the E1 rectangular (`m = 5`)
  and E2 five-point smoothed vectors, and the E3–E5 background/gross/net
  values at `c1 = 2`, `c2 = 5`, `m = 1`.
- `dollar_min.spe` — minimal dollar-format spectrum (magic `$SPEC_ID:`,
  live-real `$MEAS_TIM:`, `$DATA: 0 7` + 8 one-float lines, `keV`-suffixed
  triplets, ignored `$ROI:`/`$PRESETS:`/`$ENER_FIT:` sections).
- `plain_min.spe` — minimal plain-format spectrum (`key: value` headers,
  an ignored unknown key, double-space `Energy Fit:`, `SPECTRUM` channel
  lines) holding the same counts as `dollar_min.spe`, so the cross-format
  counts-equality check mirrors the upstream `counts[100]` pin.
