---
title: Unfold a spectrum
sidebar:
  order: 20
---

Nucleide unfolds neutron spectra from activation-detector measurements with
three adjustment methods sharing one convergence contract: the SAND-II
iterative adjustment (McElroy et al., AFWL-TR-67-41, 1967), a damped
least-squares adjustment, and the GRAVEL chi-square-weighted iteration. A
caller-supplied guess spectrum is adjusted until folding it through the
detector response matrix reproduces the measured rates. This tutorial covers
the SAND-II Python API; `staysl` and `gravel` take the same response/rates/
guess inputs (plus caller measurement sigmas) and are documented in the
[Python API reference](../../reference/python-api.mdx#nucleideunfold). The
implementation lives in `crates/unfold`. For the
fold/update equations (S1–S3) and the pinned correctness checks (U1–U18), see
the [Neutron spectrum unfolding theory](../../theory/unfolding.mdx) page.

Every input is caller data — the response matrix, the measured rates, the
guess, and the energy-group bounds. Response rows can be hand-built as below
or produced at runtime from the IRDFF-II response pack
(`nucleide.data.fetch_irdff` + `parse_irdff_g725`); no detector-response
library is vendored.

## The forward problem

Each detector (or detector reaction) responds to the neutron spectrum across
the energy groups: the response matrix holds **one row per detector**, with
one non-negative entry per energy group. Folding a spectrum through the
matrix sums each row against the spectrum, giving the calculated rate of that
detector — a plain matrix-vector product:

```python
from nucleide.unfold import forward_fold

response = [
    [1.0, 0.5, 0.0],  # detector 0 responds to groups 0 and 1
    [0.0, 0.5, 1.0],  # detector 1 responds to groups 1 and 2
]
spectrum = [2.0, 4.0, 1.0]
print(forward_fold(response, spectrum))
```

```text
[4.0, 3.0]
```

Unfolding inverts this map. `sandii(response, rates, guess)` adjusts the
guess — one strictly positive value per energy group — until its fold matches
the measured `rates` (one per detector). The update is multiplicative, so a
strictly positive guess stays strictly positive. The same fold is also the
natural way to build a practice data set: fold a known spectrum into
synthetic rates, then unfold and compare against the truth.

## Recover a determined system

When every group is seen by its own detector (a square, nearly-diagonal
response), recovery is fully determined. Here the response rows are
hand-built Gaussian bumps, the truth is a thermal Maxwellian plus a `1/E`
tail sampled at log-spaced group midpoints (MeV), and the guess is the truth
scaled by 3:

```python
import math

from nucleide.unfold import sandii

def synthetic_response(n_det, n_groups, width):
    rows = []
    for i in range(n_det):
        center = i * (n_groups - 1) / max(n_det - 1, 1)
        rows.append([math.exp(-(((j - center) / width) ** 2)) + 1e-3
                     for j in range(n_groups)])
    return rows

def log_midpoints(n_groups, lo, hi):
    step = math.log(hi / lo) / n_groups
    return [math.exp(math.log(lo) + (j + 0.5) * step) for j in range(n_groups)]

def max_rel_err(a, b):
    return max(abs(x - y) / y for x, y in zip(a, b))

midpoints = log_midpoints(6, 1e-6, 10.0)  # MeV
truth = [e * math.exp(-e / 2.53e-5) + 1e-6 / e for e in midpoints]
response = synthetic_response(6, 6, 0.9)  # one detector per group
rates = forward_fold(response, truth)

out = sandii(response, rates, [3.0 * p for p in truth],
             tolerance=1e-12, max_iterations=10_000)
print(out["iterations"], out["max_rel_change"])
print(out["rate_factors"])
print(f"worst relative error vs truth: {max_rel_err(out['spectrum'], truth):.3e}")
```

```text
2 0.0
[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
worst relative error vs truth: 3.284e-16
```

The first adjustment lands on the truth (up to round-off); the second is a
no-op that confirms convergence, so `iterations` is 2 and the final
`max_rel_change` is `0.0`. The returned dict carries the adjusted `spectrum`,
the refolded `rates`, the per-detector `rate_factors` (measured over folded
rate — all ones at the solution), plus `iterations`, `tolerance`, and
`max_rel_change`.

## More groups than detectors

The realistic case has far fewer detectors than energy groups, so infinitely
many spectra fit the rates. The adjustment moves only what the detectors can
see: any part of the guess lying in directions no detector responds to — the
response matrix's *null space* — survives iteration exactly. With one
detector responding to a single group:

```python
out = sandii([[1.0, 0.0, 0.0]], [4.0], [2.0, 7.0, 3.0])
print(out["spectrum"])  # group 0 moves; the unseen groups keep the guess
print(out["rate_factors"], out["iterations"])
```

```text
[4.0, 7.0, 3.0]
[1.0] 2
```

Group 0 moves to the measured rate; groups 1 and 2 keep the guess values
bit-for-bit. At a realistic scale — 6 detectors over 24 groups, with a
non-uniformly biased guess:

```python
response = synthetic_response(6, 24, 3.0)
midpoints = log_midpoints(24, 1e-6, 12.0)
truth = [e * math.exp(-e / 2.53e-5) + 1e-6 / e for e in midpoints]
rates = forward_fold(response, truth)
guess = [p * (1.3 + 0.7 * (j % 5) / 4.0) for j, p in enumerate(truth)]

out = sandii(response, rates, guess, tolerance=1e-9, max_iterations=50_000)
print(out["iterations"])
print([round(f, 9) for f in out["rate_factors"]])
print(f"{max_rel_err(guess, truth):.4f} -> {max_rel_err(out['spectrum'], truth):.4f}")
```

```text
33644
[1.000000001, 1.0, 1.0, 1.0, 1.0, 0.999999999]
1.0000 -> 0.3989
```

The rates are reproduced to `1e-9` and the guess error is cut by about 2.5x,
but no further: the remaining bias lies in directions the detectors cannot
see and keeps the guess. The guess is not just a starting point — in unseen
directions it is the answer, so a physically reasoned guess matters.

## Tune the stopping criteria

A run converges when the largest per-group relative change between successive
adjustments drops strictly below `tolerance`. The defaults are `tolerance=1e-3`
and `max_iterations=200` — library defaults, not values from the report. The
same 6-detector, 24-group system shows the trade:

```python
out = sandii(response, rates, guess)
print(out["iterations"], f"{out['max_rel_change']:.3e}", out["tolerance"])
out = sandii(response, rates, guess, tolerance=1e-9, max_iterations=50_000)
print(out["iterations"], f"{out['max_rel_change']:.3e}")
```

```text
31 9.251e-04 0.001
33644 9.997e-10
```

Loosening the tolerance stops after 31 adjustments; tightening it to `1e-9`
keeps iterating. `max_iterations` is the explicit cap on adjustments — raise
it when you tighten the tolerance, or the run ends in the error below.

## When the cap runs out

Exhausting the iteration cap without meeting the tolerance raises a clear
error — non-convergence is a hard fail, never a silent partial spectrum:

```python
try:
    sandii(response, rates, guess, tolerance=1e-12, max_iterations=1)
except ValueError as exc:
    print(exc)
```

```text
unfold: spectral adjustment did not converge within its iteration cap
```

Likewise, a set of rates no positive spectrum can produce (a detector with a
zero response row against a nonzero rate) raises a clear error naming that
detector. Malformed inputs — ragged or negative response entries, non-positive
guess values, mismatched lengths — are rejected the same way, each naming the
cause.

## See also

- [Neutron spectrum unfolding theory](../../theory/unfolding.mdx) for the
  S1–S3 equation set and the U1–U8 correctness checks.
- `tests/test_unfold.py` for replays of those checks and the input-error
  cases.
- [Python API reference](../../reference/python-api.mdx#nucleideunfold) for
  the full signatures.
