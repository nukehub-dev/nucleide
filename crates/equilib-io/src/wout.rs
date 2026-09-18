//! Classic-netCDF `wout` reader plus flux-surface Jacobian helpers.
//!
//! Reads the documented classic-format `wout` variables (dims `radius` /
//! `mn_mode` / `mn_mode_nyq`; mode maps `xm`/`xn`; full-mesh geometry
//! `rmnc`/`zmns`; half-mesh stream function `lmns`; Nyquist-spectrum
//! Jacobian `gmnc`; later-use fields `bmnc`, `bsubumnc`, `bsubvmnc`,
//! `bsubsmns`, `currumnc`, `currvmnc`; scalars `nfp`, `ns`, `mpol`,
//! `ntor`, `phiedge`, `volume_p`) and evaluates the Jacobian Fourier sum
//! the wall-load mapping is built on. Reads data, never solves
//! equilibria.

use std::collections::BTreeMap;

use crate::classic::{ClassicFile, ClassicVariant};
use crate::error::{Error, Result};

/// Later-use `wout` fields carried through for the wall-load consumer
/// (`bmnc`, `bsubumnc`, `bsubvmnc`, `bsubsmns`, `currumnc`, `currvmnc`).
/// Only the variables present in the file appear here.
pub const EXTRA_FIELDS: [&str; 6] = [
    "bmnc", "bsubumnc", "bsubvmnc", "bsubsmns", "currumnc", "currvmnc",
];

/// Row-major 2-D field: `values[row * n_cols + col]`.
#[derive(Debug, Clone, PartialEq)]
pub struct Matrix {
    /// Row count (radius entries as stored in the file).
    pub n_rows: usize,
    /// Column count (mode entries).
    pub n_cols: usize,
    /// Row-major values.
    pub values: Vec<f64>,
}

impl Matrix {
    /// Borrow one row (radius entry) of coefficients.
    pub fn row(&self, row: usize) -> Result<&[f64]> {
        if row >= self.n_rows {
            return Err(Error::OutOfRange {
                what: "matrix row",
                index: row,
                len: self.n_rows,
            });
        }
        Ok(&self.values[row * self.n_cols..(row + 1) * self.n_cols])
    }

    /// Split into nested rows (outer over radius, inner over modes).
    pub fn to_nested(&self) -> Vec<Vec<f64>> {
        self.values
            .chunks(self.n_cols)
            .map(<[f64]>::to_vec)
            .collect()
    }
}

/// Parsed classic `wout` file: the reader-minimum variables plus optional
/// scalars and later-use fields, with the Nyquist mode maps resolved.
#[derive(Debug, Clone, PartialEq)]
pub struct Wout {
    /// CDF-1 or CDF-2, as probed.
    pub variant: ClassicVariant,
    /// Field periods (`nfp` scalar).
    pub nfp: i32,
    /// Flux-surface count (`ns` scalar, full mesh).
    pub ns: usize,
    /// Full-spectrum mode count (`mn_mode`).
    pub n_modes: usize,
    /// Nyquist-spectrum mode count (`mn_mode_nyq`).
    pub n_modes_nyq: usize,
    /// Poloidal resolution (`mpol` scalar, when present).
    pub mpol: Option<i32>,
    /// Toroidal resolution (`ntor` scalar, when present).
    pub ntor: Option<i32>,
    /// Toroidal flux at the edge (`phiedge` scalar, when present).
    pub phiedge: Option<f64>,
    /// Plasma volume (`volume_p` scalar, when present).
    pub volume_p: Option<f64>,
    /// Poloidal mode numbers (full spectrum, length `mn_mode`).
    pub xm: Vec<f64>,
    /// Toroidal mode numbers (full spectrum, length `mn_mode`).
    pub xn: Vec<f64>,
    /// Full-mesh cylindrical-R cosine components, `(radius, mn_mode)`.
    pub rmnc: Matrix,
    /// Full-mesh cylindrical-Z sine components, `(radius, mn_mode)`.
    pub zmns: Matrix,
    /// Half-mesh stream-function sine components, `(radius, mn_mode)`
    /// with radius entries as stored.
    pub lmns: Matrix,
    /// Half-mesh Jacobian cosine components, `(radius, mn_mode_nyq)`
    /// with radius entries as stored.
    pub gmnc: Matrix,
    /// Later-use fields present in the file, keyed by variable name.
    pub fields: BTreeMap<String, Matrix>,
    /// Resolved Nyquist mode maps: `xm_nyq`/`xn_nyq` when the file
    /// carries them, else `xm`/`xn` when their length already matches
    /// `mn_mode_nyq` (small synthetics), else absent — [`Wout::jacobian`]
    /// then fails loudly with `MissingVariable("xm_nyq")`.
    nyq: Option<(Vec<f64>, Vec<f64>)>,
}

impl Wout {
    /// Read a `wout` file from disk.
    pub fn read_file(path: &std::path::Path) -> Result<Wout> {
        let bytes = std::fs::read(path).map_err(|e| Error::Io(e.to_string()))?;
        Wout::from_bytes(&bytes)
    }

    /// Parse a `wout` file from memory.
    pub fn from_bytes(bytes: &[u8]) -> Result<Wout> {
        Wout::from_classic(&crate::classic::parse(bytes)?)
    }

    /// Build from an already-parsed classic file.
    pub fn from_classic(file: &ClassicFile) -> Result<Wout> {
        let radius = file.dim_id("radius")?;
        let mn_mode = file.dim_id("mn_mode")?;
        let mn_nyq = file.dim_id("mn_mode_nyq")?;
        let ns = file.var("ns")?.data.as_i64_scalar("ns")?;
        if ns < 0 {
            return Err(Error::BadValue("ns is negative".into()));
        }
        let ns = ns as usize;
        let n_modes = file.dims[mn_mode].len as usize;
        let n_modes_nyq = file.dims[mn_nyq].len as usize;
        if file.dims[radius].len as usize != ns {
            return Err(Error::BadValue("ns disagrees with radius length".into()));
        }
        let nfp = file.var("nfp")?.data.as_i64_scalar("nfp")?;
        let nfp = i32::try_from(nfp).map_err(|_| Error::BadValue("nfp out of range".into()))?;

        let xm = file.var("xm")?.data.as_f64_vec()?;
        let xn = file.var("xn")?.data.as_f64_vec()?;
        if xm.len() != n_modes {
            return Err(Error::BadShape {
                what: "xm",
                expected: n_modes,
                got: xm.len(),
            });
        }
        if xn.len() != n_modes {
            return Err(Error::BadShape {
                what: "xn",
                expected: n_modes,
                got: xn.len(),
            });
        }

        let rmnc = full_matrix(file, "rmnc", radius, mn_mode, ns, n_modes)?;
        let zmns = full_matrix(file, "zmns", radius, mn_mode, ns, n_modes)?;
        let lmns = wide_matrix(file, "lmns", radius, mn_mode, n_modes)?;
        let gmnc = wide_matrix(file, "gmnc", radius, mn_nyq, n_modes_nyq)?;

        let nyq = match (
            file.vars.iter().find(|v| v.name == "xm_nyq"),
            file.vars.iter().find(|v| v.name == "xn_nyq"),
        ) {
            (Some(xm_v), Some(xn_v)) => {
                let a = xm_v.data.as_f64_vec()?;
                let b = xn_v.data.as_f64_vec()?;
                if a.len() != n_modes_nyq {
                    return Err(Error::BadShape {
                        what: "xm_nyq",
                        expected: n_modes_nyq,
                        got: a.len(),
                    });
                }
                if b.len() != n_modes_nyq {
                    return Err(Error::BadShape {
                        what: "xn_nyq",
                        expected: n_modes_nyq,
                        got: b.len(),
                    });
                }
                Some((a, b))
            }
            (None, None) if xm.len() == n_modes_nyq => Some((xm.clone(), xn.clone())),
            _ => None,
        };

        let mut fields = BTreeMap::new();
        for name in EXTRA_FIELDS {
            if let Some(var) = file.vars.iter().find(|v| v.name == name) {
                fields.insert(name.to_string(), any_matrix(var, file)?);
            }
        }

        Ok(Wout {
            variant: file.variant,
            nfp,
            ns,
            n_modes,
            n_modes_nyq,
            mpol: opt_int(file, "mpol")?,
            ntor: opt_int(file, "ntor")?,
            phiedge: opt_float(file, "phiedge")?,
            volume_p: opt_float(file, "volume_p")?,
            xm,
            xn,
            rmnc,
            zmns,
            lmns,
            gmnc,
            fields,
            nyq,
        })
    }

    /// Borrow the resolved Nyquist mode maps (`(xm_nyq, xn_nyq)`).
    pub fn nyq_modes(&self) -> Result<(&[f64], &[f64])> {
        match &self.nyq {
            Some((m, n)) => Ok((m, n)),
            None => Err(Error::MissingVariable("xm_nyq".to_string())),
        }
    }

    /// Borrow one half-mesh Jacobian coefficient row.
    pub fn jacobian_row(&self, row: usize) -> Result<&[f64]> {
        self.gmnc.row(row)
    }

    /// Evaluate the Jacobian Fourier sum (J1) at one half-mesh row:
    /// `sqrt(g)(θ,ζ) = Σ_j c_j·cos(m_j·θ − n_j·ζ)` over the resolved
    /// Nyquist maps. Angles are in radians; `ζ` follows the file's
    /// toroidal mode convention (one field period per the documented
    /// `wout` layout — the caller scales for full-torus work).
    pub fn jacobian(&self, row: usize, theta: f64, zeta: f64) -> Result<f64> {
        let (m, n) = self.nyq_modes()?;
        fourier_jacobian(m, n, self.gmnc.row(row)?, theta, zeta)
    }

    /// Evaluate the Jacobian on a tensor grid (J2): `ntheta` poloidal
    /// points over `[0, 2π)` (endpoint excluded, periodic) by `nzeta`
    /// toroidal points over one field period `[0, 2π/nfp)`.
    pub fn jacobian_grid(&self, row: usize, ntheta: usize, nzeta: usize) -> Result<Matrix> {
        if ntheta == 0 {
            return Err(Error::BadValue("ntheta must be positive".into()));
        }
        if nzeta == 0 {
            return Err(Error::BadValue("nzeta must be positive".into()));
        }
        if self.nfp <= 0 {
            return Err(Error::BadValue("nfp must be positive".into()));
        }
        let coeffs = self.gmnc.row(row)?.to_vec();
        let (m, n) = self.nyq_modes()?;
        let mut values = Vec::with_capacity(ntheta * nzeta);
        for i in 0..ntheta {
            let theta = 2.0 * std::f64::consts::PI * i as f64 / ntheta as f64;
            for j in 0..nzeta {
                let zeta = 2.0 * std::f64::consts::PI * j as f64 / (nzeta as f64 * self.nfp as f64);
                values.push(fourier_jacobian(m, n, &coeffs, theta, zeta)?);
            }
        }
        Ok(Matrix {
            n_rows: ntheta,
            n_cols: nzeta,
            values,
        })
    }

    /// Surface average of the Jacobian (J3): the `(m,n) = (0,0)`
    /// coefficient of the half-mesh row. Fails loudly when the Nyquist
    /// spectrum carries no axisymmetric mode.
    pub fn mean_jacobian(&self, row: usize) -> Result<f64> {
        let (m, n) = self.nyq_modes()?;
        let coeffs = self.gmnc.row(row)?;
        for (j, (&mm, &nn)) in m.iter().zip(n.iter()).enumerate() {
            if mm == 0.0 && nn == 0.0 {
                let c = coeffs[j];
                if !c.is_finite() {
                    return Err(Error::BadValue("non-finite gmnc entry".into()));
                }
                return Ok(c);
            }
        }
        Err(Error::BadValue(
            "Nyquist spectrum has no (m,n)=(0,0) mode".into(),
        ))
    }

    /// Half-mesh row count (Jacobian rows), as stored.
    pub fn half_rows(&self) -> usize {
        self.gmnc.n_rows
    }
}

/// Evaluate `Σ_j c_j·cos(m_j·θ − n_j·ζ)` over caller-supplied mode maps.
/// Length mismatches are [`Error::BadShape`]; non-finite inputs are
/// [`Error::BadValue`]. This is the J1 kernel behind [`Wout::jacobian`],
/// kept free so the wall-load consumer can fold caller-side spectra.
pub fn fourier_jacobian(
    m: &[f64],
    n: &[f64],
    coeffs: &[f64],
    theta: f64,
    zeta: f64,
) -> Result<f64> {
    if m.len() != coeffs.len() {
        return Err(Error::BadShape {
            what: "xm",
            expected: coeffs.len(),
            got: m.len(),
        });
    }
    if n.len() != coeffs.len() {
        return Err(Error::BadShape {
            what: "xn",
            expected: coeffs.len(),
            got: n.len(),
        });
    }
    if !theta.is_finite() || !zeta.is_finite() {
        return Err(Error::BadValue("non-finite angle".into()));
    }
    let mut acc = 0.0;
    for ((&mm, &nn), &c) in m.iter().zip(n.iter()).zip(coeffs.iter()) {
        if !mm.is_finite() || !nn.is_finite() || !c.is_finite() {
            return Err(Error::BadValue("non-finite mode map or coefficient".into()));
        }
        acc += c * (mm * theta - nn * zeta).cos();
    }
    if !acc.is_finite() {
        return Err(Error::BadValue("non-finite Jacobian sum".into()));
    }
    Ok(acc)
}

fn opt_int(file: &ClassicFile, name: &str) -> Result<Option<i32>> {
    match file.vars.iter().find(|v| v.name == name) {
        None => Ok(None),
        Some(var) => {
            let v = var.data.as_i64_scalar("scalar")?;
            i32::try_from(v)
                .map(Some)
                .map_err(|_| Error::BadValue(format!("{name} out of range")))
        }
    }
}

fn opt_float(file: &ClassicFile, name: &str) -> Result<Option<f64>> {
    match file.vars.iter().find(|v| v.name == name) {
        None => Ok(None),
        Some(var) => var.data.as_f64_scalar("scalar").map(Some),
    }
}

/// Strict `(radius, mode)` matrix with an exact row count (full mesh).
fn full_matrix(
    file: &ClassicFile,
    name: &str,
    radius: usize,
    mode: usize,
    n_rows: usize,
    n_cols: usize,
) -> Result<Matrix> {
    let var = file.var(name)?;
    expect_dims(var, &[radius, mode])?;
    let values = var.data.as_f64_vec()?;
    if values.len() != n_rows * n_cols {
        return Err(Error::BadShape {
            what: "wout field",
            expected: n_rows * n_cols,
            got: values.len(),
        });
    }
    Ok(Matrix {
        n_rows,
        n_cols,
        values,
    })
}

/// `(radius, mode)` matrix with radius entries as stored (half mesh may
/// carry fewer rows than the full `ns`; the file is authoritative).
fn wide_matrix(
    file: &ClassicFile,
    name: &str,
    radius: usize,
    mode: usize,
    n_cols: usize,
) -> Result<Matrix> {
    let var = file.var(name)?;
    expect_dims(var, &[radius, mode])?;
    let values = var.data.as_f64_vec()?;
    if n_cols == 0 || values.len() % n_cols != 0 {
        return Err(Error::BadShape {
            what: "wout field",
            expected: n_cols,
            got: values.len(),
        });
    }
    Ok(Matrix {
        n_rows: values.len() / n_cols,
        n_cols,
        values,
    })
}

/// Later-use field with its stored shape (row-major as read).
fn any_matrix(var: &crate::classic::Var, file: &ClassicFile) -> Result<Matrix> {
    let lens: Vec<usize> = var
        .dim_ids
        .iter()
        .map(|id| {
            file.dims
                .get(*id as usize)
                .map(|d| d.len as usize)
                .ok_or(Error::BadHeader("var dim id out of range"))
        })
        .collect::<Result<_>>()?;
    let values = var.data.as_f64_vec()?;
    let (n_rows, n_cols) = match lens.as_slice() {
        [] => (1, 1),
        [n] => (1, *n),
        [r, c] => (*r, *c),
        _ => {
            let last = *lens.last().expect("non-empty");
            let rows = lens[..lens.len() - 1].iter().product();
            (rows, last)
        }
    };
    if values.len() != n_rows * n_cols {
        return Err(Error::BadShape {
            what: "wout field",
            expected: n_rows * n_cols,
            got: values.len(),
        });
    }
    Ok(Matrix {
        n_rows,
        n_cols,
        values,
    })
}

fn expect_dims(var: &crate::classic::Var, want: &[usize]) -> Result<()> {
    if var.dim_ids.len() != want.len()
        || var
            .dim_ids
            .iter()
            .zip(want.iter())
            .any(|(got, w)| *got as usize != *w)
    {
        return Err(Error::UnexpectedDims(var.name.clone()));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::classic::{parse, probe_variant, ClassicVariant};

    /// Minimal classic writer for synthetic `wout` gates: header layout
    /// follows the same field order the reader walks (never solver data).
    struct Synth {
        variant: ClassicVariant,
        dims: Vec<(String, u64)>,
        vars: Vec<SynthVar>,
    }

    struct SynthVar {
        name: &'static str,
        dims: Vec<u32>,
        code: u32,
        payload: Vec<u8>,
    }

    fn be32(buf: &mut Vec<u8>, v: u32) {
        buf.extend_from_slice(&v.to_be_bytes());
    }

    fn be64(buf: &mut Vec<u8>, v: u64) {
        buf.extend_from_slice(&v.to_be_bytes());
    }

    fn pad_name(buf: &mut Vec<u8>, variant: ClassicVariant, name: &str, len: u64) {
        be32(buf, name.len() as u32);
        buf.extend_from_slice(name.as_bytes());
        while buf.len() % 4 != 0 {
            buf.push(0);
        }
        if variant == ClassicVariant::Cdf2 {
            be64(buf, len);
        } else {
            be32(buf, len as u32);
        }
    }

    fn doubles(vals: &[f64]) -> Vec<u8> {
        let mut out = Vec::with_capacity(vals.len() * 8);
        for v in vals {
            out.extend_from_slice(&v.to_be_bytes());
        }
        out
    }

    fn ints(vals: &[i32]) -> Vec<u8> {
        let mut out = Vec::with_capacity(vals.len() * 4);
        for v in vals {
            out.extend_from_slice(&v.to_be_bytes());
        }
        out
    }

    impl Synth {
        fn bytes(&self) -> Vec<u8> {
            let mut head = vec![0x43, 0x44, 0x46];
            head.push(match self.variant {
                ClassicVariant::Cdf1 => 0x01,
                ClassicVariant::Cdf2 => 0x02,
            });
            if self.variant == ClassicVariant::Cdf2 {
                be64(&mut head, 0);
            } else {
                be32(&mut head, 0);
            }
            // dim list
            be32(&mut head, 10);
            be32(&mut head, self.dims.len() as u32);
            for (name, len) in &self.dims {
                pad_name(&mut head, self.variant, name, *len);
            }
            // global attrs: absent
            be32(&mut head, 0);
            be32(&mut head, 0);
            // var list
            be32(&mut head, 11);
            be32(&mut head, self.vars.len() as u32);
            // Data starts right after the header; lay payloads out back to
            // back at 4-byte alignment.
            let mut header_len = head.len();
            for var in &self.vars {
                header_len += 4 + var.name.len() + ((4 - var.name.len() % 4) % 4);
                header_len += 4 + 4 * var.dims.len() + 8 + 4;
                header_len += if self.variant == ClassicVariant::Cdf2 {
                    16
                } else {
                    8
                };
            }
            let mut begin = header_len;
            let mut entries = Vec::new();
            for var in &self.vars {
                entries.push((var, begin));
                begin += var.payload.len() + ((4 - var.payload.len() % 4) % 4);
            }
            for (var, off) in &entries {
                be32(&mut head, var.name.len() as u32);
                head.extend_from_slice(var.name.as_bytes());
                while head.len() % 4 != 0 {
                    head.push(0);
                }
                be32(&mut head, var.dims.len() as u32);
                for d in &var.dims {
                    be32(&mut head, *d);
                }
                be32(&mut head, 0);
                be32(&mut head, 0); // no per-var attrs
                be32(&mut head, var.code);
                if self.variant == ClassicVariant::Cdf2 {
                    be64(&mut head, var.payload.len() as u64);
                    be64(&mut head, *off as u64);
                } else {
                    be32(&mut head, var.payload.len() as u32);
                    be32(&mut head, *off as u32);
                }
            }
            assert_eq!(head.len(), header_len);
            let mut out = head;
            for (var, _) in &entries {
                out.extend_from_slice(&var.payload);
                while out.len() % 4 != 0 {
                    out.push(0);
                }
            }
            out
        }
    }

    /// Small synthetic `wout`: `ns = 3` rows on the shared `radius` dim
    /// (half-mesh fields ride the same dim, per the pinned dim listing),
    /// two full modes, two Nyquist modes. `with_nyq` adds explicit
    /// `xm_nyq`/`xn_nyq` maps.
    fn synth(variant: ClassicVariant, with_nyq: bool) -> Vec<u8> {
        let mut vars = vec![
            SynthVar {
                name: "nfp",
                dims: vec![],
                code: 4,
                payload: ints(&[3]),
            },
            SynthVar {
                name: "ns",
                dims: vec![],
                code: 4,
                payload: ints(&[3]),
            },
            SynthVar {
                name: "mpol",
                dims: vec![],
                code: 4,
                payload: ints(&[5]),
            },
            SynthVar {
                name: "ntor",
                dims: vec![],
                code: 4,
                payload: ints(&[4]),
            },
            SynthVar {
                name: "phiedge",
                dims: vec![],
                code: 6,
                payload: doubles(&[1.5]),
            },
            SynthVar {
                name: "volume_p",
                dims: vec![],
                code: 6,
                payload: doubles(&[30.0]),
            },
            SynthVar {
                name: "xm",
                dims: vec![1],
                code: 6,
                payload: doubles(&[0.0, 2.0]),
            },
            SynthVar {
                name: "xn",
                dims: vec![1],
                code: 6,
                payload: doubles(&[0.0, -1.0]),
            },
            SynthVar {
                name: "rmnc",
                dims: vec![0, 1],
                code: 6,
                payload: doubles(&[1.0, 0.1, 1.1, 0.1, 1.2, 0.1]),
            },
            SynthVar {
                name: "zmns",
                dims: vec![0, 1],
                code: 6,
                payload: doubles(&[0.0, 0.2, 0.0, 0.25, 0.0, 0.3]),
            },
            SynthVar {
                name: "lmns",
                dims: vec![0, 1],
                code: 6,
                payload: doubles(&[0.0, 0.01, 0.0, 0.02, 0.0, 0.03]),
            },
            SynthVar {
                name: "gmnc",
                dims: vec![0, 2],
                code: 6,
                payload: doubles(&[2.0, 0.0, 1.0, 0.5, 1.0, -1.0]),
            },
            SynthVar {
                name: "bmnc",
                dims: vec![0, 1],
                code: 6,
                payload: doubles(&[5.0, 0.0, 5.1, 0.0, 5.2, 0.0]),
            },
        ];
        if with_nyq {
            vars.push(SynthVar {
                name: "xm_nyq",
                dims: vec![2],
                code: 6,
                payload: doubles(&[0.0, 2.0]),
            });
            vars.push(SynthVar {
                name: "xn_nyq",
                dims: vec![2],
                code: 6,
                payload: doubles(&[0.0, -1.0]),
            });
        }
        Synth {
            variant,
            dims: vec![
                ("radius".to_string(), 3),
                ("mn_mode".to_string(), 2),
                ("mn_mode_nyq".to_string(), 2),
            ],
            vars,
        }
        .bytes()
    }

    #[test]
    fn synthetic_round_trip_cdf1_and_cdf2() {
        for variant in [ClassicVariant::Cdf1, ClassicVariant::Cdf2] {
            let bytes = synth(variant, true);
            assert_eq!(probe_variant(&bytes).unwrap(), variant);
            let w = Wout::from_bytes(&bytes).unwrap();
            assert_eq!(w.variant, variant);
            assert_eq!(w.nfp, 3);
            assert_eq!(w.ns, 3);
            assert_eq!((w.n_modes, w.n_modes_nyq), (2, 2));
            assert_eq!(w.mpol, Some(5));
            assert_eq!(w.ntor, Some(4));
            assert_eq!(w.phiedge, Some(1.5));
            assert_eq!(w.volume_p, Some(30.0));
            assert_eq!(w.xm, vec![0.0, 2.0]);
            assert_eq!(w.xn, vec![0.0, -1.0]);
            assert_eq!((w.rmnc.n_rows, w.rmnc.n_cols), (3, 2));
            assert_eq!((w.gmnc.n_rows, w.gmnc.n_cols), (3, 2));
            assert_eq!(w.half_rows(), 3);
            assert_eq!(w.jacobian_row(0).unwrap(), &[2.0, 0.0]);
            let (m, n) = w.nyq_modes().unwrap();
            assert_eq!((m, n), (&[0.0, 2.0][..], &[0.0, -1.0][..]));
            let b = w.fields.get("bmnc").expect("later-use field kept");
            assert_eq!((b.n_rows, b.n_cols), (3, 2));
        }
    }

    #[test]
    fn nyq_falls_back_to_xm_when_lengths_match() {
        let w = Wout::from_bytes(&synth(ClassicVariant::Cdf1, false)).unwrap();
        assert_eq!(w.jacobian(0, 0.7, -0.2).unwrap(), 2.0);
    }

    #[test]
    fn jacobian_hand_vectors() {
        let w = Wout::from_bytes(&synth(ClassicVariant::Cdf1, true)).unwrap();
        // Row 0 is pure (0,0): identically 2.0, mean 2.0.
        assert_eq!(w.jacobian(0, 0.0, 0.0).unwrap(), 2.0);
        assert_eq!(w.jacobian(0, 1.3, 2.1).unwrap(), 2.0);
        assert_eq!(w.mean_jacobian(0).unwrap(), 2.0);
        // Row 1 is 1.0 + 0.5*cos(2θ + ζ): hand checks.
        assert_eq!(w.jacobian(1, 0.0, 0.0).unwrap(), 1.5);
        assert!((w.jacobian(1, std::f64::consts::FRAC_PI_2, 0.0).unwrap() - 0.5).abs() < 1e-12);
        assert!((w.jacobian(1, 0.0, std::f64::consts::PI).unwrap() - 0.5).abs() < 1e-12);
        assert_eq!(w.mean_jacobian(1).unwrap(), 1.0);
        // Row 2 is 1.0 - cos(2θ + ζ): zero at the origin, 2.0 at θ=π/2.
        assert_eq!(w.jacobian(2, 0.0, 0.0).unwrap(), 0.0);
        assert!((w.jacobian(2, std::f64::consts::FRAC_PI_2, 0.0).unwrap() - 2.0).abs() < 1e-12);
        // Grid: 4 poloidal by 3 toroidal points over one field period.
        let grid = w.jacobian_grid(1, 4, 3).unwrap();
        assert_eq!((grid.n_rows, grid.n_cols), (4, 3));
        assert_eq!(grid.values[0], 1.5);
        assert!(
            (grid.values[1]
                - w.jacobian(1, 0.0, 2.0 * std::f64::consts::PI / 9.0)
                    .unwrap())
            .abs()
                < 1e-12
        );
        // Free kernel agrees and validates shapes.
        assert_eq!(
            fourier_jacobian(&[0.0, 2.0], &[0.0, -1.0], &[1.0, 0.5], 0.0, 0.0).unwrap(),
            1.5
        );
        assert!(fourier_jacobian(&[0.0], &[0.0, 1.0], &[1.0, 0.5], 0.0, 0.0).is_err());
        assert!(fourier_jacobian(&[0.0], &[0.0], &[f64::NAN], 0.0, 0.0).is_err());
        // Out-of-range rows and empty grids fail loudly.
        assert!(w.jacobian(7, 0.0, 0.0).is_err());
        assert!(w.jacobian_grid(0, 0, 3).is_err());
        assert!(w.jacobian_grid(0, 3, 0).is_err());
    }

    #[test]
    fn missing_and_misshapen_inputs_fail_loudly() {
        // No gmnc at all.
        let mut file = parse(&synth(ClassicVariant::Cdf1, true)).unwrap();
        file.vars.retain(|v| v.name != "gmnc");
        let err = Wout::from_classic(&file).unwrap_err();
        assert_eq!(err, Error::MissingVariable("gmnc".to_string()));
        // Transposed dims on rmnc.
        let mut file = parse(&synth(ClassicVariant::Cdf1, true)).unwrap();
        file.vars
            .iter_mut()
            .find(|v| v.name == "rmnc")
            .expect("rmnc")
            .dim_ids = vec![1, 0];
        assert_eq!(
            Wout::from_classic(&file).unwrap_err(),
            Error::UnexpectedDims("rmnc".to_string())
        );
        // xm length disagrees with mn_mode.
        let mut file = parse(&synth(ClassicVariant::Cdf1, true)).unwrap();
        file.vars
            .iter_mut()
            .find(|v| v.name == "xm")
            .expect("xm")
            .data = crate::classic::Decoded::Floats(vec![0.0]);
        assert_eq!(
            Wout::from_classic(&file).unwrap_err(),
            Error::BadShape {
                what: "xm",
                expected: 2,
                got: 1
            }
        );
        // Mean without an axisymmetric mode.
        let mut w = Wout::from_bytes(&synth(ClassicVariant::Cdf1, true)).unwrap();
        w.nyq = Some((vec![1.0, 2.0], vec![0.0, -1.0]));
        assert!(w.mean_jacobian(0).is_err());
    }

    #[test]
    fn record_variables_rejected() {
        let file = Synth {
            variant: ClassicVariant::Cdf1,
            dims: vec![
                ("radius".to_string(), 3),
                ("mn_mode".to_string(), 2),
                ("mn_mode_nyq".to_string(), 2),
                ("time".to_string(), 0),
            ],
            vars: vec![SynthVar {
                name: "trace",
                dims: vec![3],
                code: 6,
                payload: doubles(&[]),
            }],
        }
        .bytes();
        assert_eq!(
            parse(&file).unwrap_err(),
            Error::RecordVariable("trace".to_string())
        );
    }
}
