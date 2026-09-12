# mcpl — Synthetic SSW↔MCPL conversion pair authored for Nucleide (no license needed)

Hand-built fixtures for the neutron/gamma-only SSW↔MCPL conversion v1
(`crates/mcpl-io/src/ssw.rs`). No MCNP run, no vendored upstream bytes:
every value below is the recorded generator input, so any drift is a loud
regression, never a mystery.

`reference.w` (655 bytes): minimal plain-MCNP5-layout SSW file — one header
record (`mcnp` / `5` / `01012026`, dump 2), table 1 (`orignp1 = -1000`,
`nrss = 2`, `ncrd = 11`, `njsw = 2`, `niss = 2`, plus the trailing zero word
the 32-byte vendored MCNP5 records carry), table 2 (`niwr = 0`, `mipts = 3`,
`kjaq = 0` plus 17 zero words, matching all four vendored `niwr == 0`
records; `mipts` passes through the converter opaquely and carries no
per-track meaning here), two surface records (`id = 100`, `type = 1`,
params `[0.0]`; `id = 200`, `type = 1`, params `[10.0]`), a 30-integer
summary (`nrss` in the count slots by analogy with the 1-surface vendored
samples, the rest zeroed opaque padding), then two 11-double tracks:

| Track | nps | bitarray | wgt | erg (MeV) | tme (shakes) | x | y | z | u | v | cs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1 | 800000800.0 | 1.0 | 2.5 | 300000.0 | 1.0 | -2.0 | 0.5 | 0.0 | 0.0 | 1.0 |
| 2 | 2 | 1600001600.0 | 0.5 | 0.662 | 0.0 | 0.0 | 0.0 | 0.0 | 1.0 | 0.0 | 0.0 |

The `bitarray` high words are oracle-verified, not guessed: probing the
upstream 2.2.8 `ssw2mcpl` binary maps `8*cell + TYPE*8e8` with `TYPE = 1` to
PDG 2112 and `TYPE = 2` to PDG 22 (cells decode to 100/200 under the local
`cell()` rule; signs keep the derived `w` on the stored `cs`). Nucleide
itself never decodes these words — the pairing below stays the explicit
caller parameter, and no upstream table is transcribed anywhere in code.

The SSW format stores no per-track surface id or particle kind, so tests
pair each track with an explicit caller parameter: track 1 → `surf = 100`,
`kind = neutron`; track 2 → `surf = 200`, `kind = gamma`. The pairing is
consistent with the file's surface list, and the upstream binary converts
this file to field-identical particles (see `validation/mcpl_vs_refs.py`
tier 4).

`ssw2mcpl_expected.mcpl`: golden output of `ssw2mcpl` over that pairing with
default options (single-prec, `surf_to_userflags`, no deck blob), encoded by
the crate's own writer. Expected particles: `(ekin 2.5, time 3.0 ms,
pdg 2112, userflags 100)` and `(ekin 0.662, time 0.0 ms, pdg 22,
userflags 200)`; positions/directions/weights verbatim from the table above.
Regenerate only deliberately (same pairing + defaults), and say why in the
commit message.
