# endl — Synthetic EEDL-style tables authored for Nucleide (no license needed)

`synthetic_eedl.txt` is a hand-built ENDL file (EEDL/EPDL scope, `pyne.endl`
framing contract) for one nucleus id (`820000000`, natural Pb): a 2-column
integrated table (`rdesc=10`, `rprop=0`, `rmod=0` with a nonzero `x1` field to
pin the `x1=None` rule), a 3-column spectra table (`rdesc=82`, `rprop=21`),
and two subshell-indexed tables (`rdesc=81`, `rmod=1`, `x1=1/2` with
`yo=9/19`) to pin `get_rx` selector filtering. Table terminators are 71
spaces + `1`; every body field is a fixed 11-character ENDL number.
