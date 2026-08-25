//! MCNP WWINP weight-window file parsing (Cartesian meshes only),
//! validated against the vendored `mcnp_wwinp_wwinp_{n,p,np}.txt` fixtures.
//!
//! Layout notes:
//! - Block 1: header counts; block 2: per-dimension coarse/fine mesh
//!   definition where MCNP duplicates boundary entries (the raw stream is
//!   `[rx_k, nf_k, rx_{k+1}]` triplets with the middle value repeated);
//!   block 3: energy upper bounds then weight-window lower bounds.
//! - Fine bounds between coarse points are linearly interpolated:
//!   `(b_hi - b_lo) * k / nf + b_lo`.
//! - Lower-bound storage order matches the file: `[particle][group][ve]`
//!   with volume elements ordered z slowest → x fastest.
//! - The writer reproduces Python's `{0:13.5E}` field formatting including
//!   two-digit exponents, enabling byte-stable round trips modulo trailing
//!   whitespace.

use std::fmt;
use std::path::Path;

/// Errors raised while reading or writing WWINP files.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    Io(String),
    /// Cylindrical meshes (`nr = 16`) are unsupported upstream too.
    Unsupported(&'static str),
    /// Structural problem with a block.
    BadStructure(String),
    /// Numeric field failed to parse.
    BadNumber {
        context: &'static str,
        text: String,
    },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(m) => write!(f, "io error: {m}"),
            Error::Unsupported(m) => write!(f, "unsupported: {m}"),
            Error::BadStructure(m) => write!(f, "malformed WWINP: {m}"),
            Error::BadNumber { context, text } => {
                write!(f, "cannot parse {context} from `{text}`")
            }
        }
    }
}

impl std::error::Error for Error {}

fn num(context: &'static str, tok: &str) -> Result<f64, Error> {
    tok.parse::<f64>().map_err(|_| Error::BadNumber {
        context,
        text: tok.to_string(),
    })
}

/// Parsed WWINP contents.
#[derive(Debug, Clone, PartialEq)]
pub struct Wwinp {
    /// 1 = neutrons only, 2 = photons or neutrons+photons.
    pub ni: u32,
    /// 10 = rectangular (only supported kind).
    pub nr: u32,
    /// Energy-group count per particle (len 1 or 2).
    pub ne: Vec<u32>,
    /// Fine mesh points per dimension.
    pub nf: [u32; 3],
    pub nft: u64,
    /// Mesh minimum corner.
    pub origin: [f64; 3],
    /// Coarse mesh points per dimension (excluding origin, MCNP style).
    pub nc: [u32; 3],
    pub nwg: u32,
    /// Raw header tail from line 1 (date/time), preserved for round trips.
    pub date_time: String,
    /// Interior coarse boundaries per dimension (MCNP convention: no origin).
    pub cm: Vec<Vec<f64>>,
    /// Fine-interval counts per dimension.
    pub fm: Vec<Vec<f64>>,
    /// Fully expanded spatial bounds per dimension.
    pub bounds: Vec<Vec<f64>>,
    /// Energy upper bounds per particle.
    pub e: Vec<Vec<f64>>,
    /// Weight-window lower bounds `[particle][group][ve]`, ve ordered
    /// z slowest → x fastest.
    pub ww: Vec<Vec<Vec<f64>>>,
}

impl Wwinp {
    /// Read and parse a WWINP file.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self, Error> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path)
            .map_err(|e| Error::Io(format!("{}: {}", path.display(), e)))?;
        Wwinp::parse(&text)
    }

    /// Parse WWINP text in memory.
    pub fn parse(text: &str) -> Result<Self, Error> {
        let mut tokens = TokenFeed::new(text.lines().peekable());

        // ---- Block 1 ----
        let l1 = tokens.line()?;
        let w1: Vec<&str> = l1.split_whitespace().collect();
        if w1.len() < 4 {
            return Err(Error::BadStructure("block-1 line 1 too short".into()));
        }
        let ni = num("ni", w1[2])? as u32;
        let nr = num("nr", w1[3])? as u32;
        let date_time = if w1.len() > 4 {
            w1[4..].join(" ")
        } else {
            String::new()
        };
        if nr != 10 {
            return Err(Error::Unsupported(
                "cylindrical WWINP (nr=16) not currently supported",
            ));
        }

        let ne: Vec<u32> = tokens
            .line()?
            .split_whitespace()
            .map(|t| num("ne", t).map(|v| v as u32))
            .collect::<Result<_, _>>()?;
        if ne.is_empty() || ne.len() > 2 {
            return Err(Error::BadStructure(format!("bad ne length {}", ne.len())));
        }

        let l3: Vec<&str> = tokens.line()?.split_whitespace().collect();
        if l3.len() < 6 {
            return Err(Error::BadStructure("block-1 line 3 too short".into()));
        }
        let mut nf = [0u32; 3];
        let mut origin = [0f64; 3];
        for i in 0..3 {
            nf[i] = num("nf", l3[i])? as u32;
            origin[i] = num("origin", l3[3 + i])?;
        }
        let nft = nf.iter().map(|v| *v as u64).product();

        let l4: Vec<&str> = tokens.line()?.split_whitespace().collect();
        if l4.len() < 4 {
            return Err(Error::BadStructure("block-1 line 4 too short".into()));
        }
        let mut nc = [0u32; 3];
        for i in 0..3 {
            nc[i] = num("nc", l4[i])? as u32;
        }
        let nwg = num("nwg", l4[3])? as u32;

        // ---- Block 2 ----
        let mut bounds = vec![Vec::<f64>::new(); 3];
        let mut cm = vec![Vec::<f64>::new(); 3];
        let mut fm = vec![Vec::<f64>::new(); 3];
        for i in 0..3 {
            let need = 3 * nc[i] as usize + 1;
            let mut raw: Vec<f64> = Vec::with_capacity(need);
            while raw.len() < need {
                for t in tokens.line()?.split_whitespace() {
                    raw.push(num("block-2", t)?);
                }
            }
            // Drop duplicated boundaries: keep j=0 and every j%3 != 0.
            let mut removed = vec![raw[0]];
            for (j, v) in raw.iter().enumerate().skip(1) {
                if j % 3 != 0 {
                    removed.push(*v);
                }
            }
            for (j, v) in removed.iter().enumerate() {
                if j % 2 == 0 {
                    bounds[i].push(*v);
                    if j != 0 {
                        cm[i].push(*v);
                    }
                } else {
                    let n_fine = *v;
                    fm[i].push(n_fine);
                    let lo = removed[j - 1];
                    let hi = removed[j + 1];
                    for k in 1..n_fine as usize {
                        bounds[i].push((hi - lo) * k as f64 / n_fine + lo);
                    }
                }
            }
        }

        // ---- Block 3 ----
        let nparticles = if ni == 1 { 1 } else { 2 };
        let mut e: Vec<Vec<f64>> = Vec::new();
        let mut ww: Vec<Vec<Vec<f64>>> = Vec::new();
        for p in 0..nparticles {
            let groups = ne[p.min(ne.len() - 1)] as usize;
            if groups == 0 {
                continue;
            }
            let mut energies = Vec::with_capacity(groups);
            while energies.len() < groups {
                for t in tokens.line()?.split_whitespace() {
                    energies.push(num("energy", t)?);
                }
            }
            let mut data = Vec::with_capacity(groups);
            for _ in 0..groups {
                let mut row = Vec::with_capacity(nft as usize);
                while row.len() < nft as usize {
                    for t in tokens.line()?.split_whitespace() {
                        row.push(num("ww bound", t)?);
                    }
                }
                data.push(row);
            }
            e.push(energies);
            ww.push(data);
        }

        Ok(Wwinp {
            ni,
            nr,
            ne,
            nf,
            nft,
            origin,
            nc,
            nwg,
            date_time,
            cm,
            fm,
            bounds,
            e,
            ww,
        })
    }

    /// Per-volume-element lower-bound vector (mesh-tag semantics).
    pub fn ww_column(&self, particle: usize, ve: usize) -> Vec<f64> {
        self.ww
            .get(particle)
            .map(|groups| groups.iter().map(|row| row[ve]).collect())
            .unwrap_or_default()
    }

    /// Write canonical WWINP text (Python-compatible numeric formatting).
    pub fn to_text(&self) -> Result<String, Error> {
        let mut out = String::new();

        // Block 1
        out += &format!(
            "{:>10}{:>10}{:>10}{:>10}{:>38}\n",
            1, 1, self.ni, self.nr, self.date_time
        );
        for g in &self.ne {
            out += &format!("{:>10}", g);
        }
        out += "\n";
        out += &fmt13(self.nf[0] as f64);
        out += &fmt13(self.nf[1] as f64);
        out += &fmt13(self.nf[2] as f64);
        out += &fmt13(self.origin[0]);
        out += &fmt13(self.origin[1]);
        out += &fmt13(self.origin[2]);
        out += "\n";
        out += &fmt13(self.nc[0] as f64);
        out += &fmt13(self.nc[1] as f64);
        out += &fmt13(self.nc[2] as f64);
        out += &fmt13(self.nwg as f64);
        out += "\n";

        // Block 2
        for i in 0..3 {
            let mut arr = vec![self.origin[i]];
            for j in 0..self.cm[i].len() {
                arr.push(self.fm[i][j]);
                arr.push(self.cm[i][j]);
                arr.push(1.0000);
            }
            out += &wrap6(&arr.iter().copied().map(fmt13).collect::<Vec<_>>());
        }

        // Block 3
        for p in 0..self.ww.len() {
            out += &wrap6(&self.e[p].iter().copied().map(fmt13).collect::<Vec<_>>());
            for group in &self.ww[p] {
                out += &wrap6(&group.iter().copied().map(fmt13).collect::<Vec<_>>());
            }
        }

        Ok(out)
    }

    /// Write to disk in canonical form.
    pub fn write_file(&self, path: impl AsRef<Path>) -> Result<(), Error> {
        std::fs::write(path, self.to_text()?).map_err(|e| Error::Io(e.to_string()))
    }
}

// ---- shared helpers ----

struct TokenFeed<'a, I: Iterator<Item = &'a str>> {
    lines: std::iter::Peekable<I>,
}

impl<'a, I: Iterator<Item = &'a str>> TokenFeed<'a, I> {
    fn new(lines: std::iter::Peekable<I>) -> Self {
        Self { lines }
    }

    fn line(&mut self) -> Result<&'a str, Error> {
        self.lines
            .next()
            .ok_or_else(|| Error::BadStructure("unexpected EOF".into()))
    }
}

/// Format like Python `{0:13.5E}`: width 13, upper-case E, two-digit exponent.
fn fmt13(v: f64) -> String {
    let s = format!("{v:.5E}");
    let (mant, exp) = s.split_once('E').expect("uppercase E always present");
    let (sign, digits) = match exp.strip_prefix('-') {
        Some(d) => ('-', d),
        None => ('+', exp),
    };
    let body = if digits.len() < 2 {
        format!("{mant}E{sign}0{digits}")
    } else {
        format!("{mant}E{sign}{digits}")
    };
    format!("{body:>13}")
}

/// Join preformatted 13-wide fields, six per line (WWINP wrapping rule).
fn wrap6(fields: &[String]) -> String {
    let mut out = String::new();
    let mut count = 0;
    for f in fields {
        out += f;
        count += 1;
        if count == 6 {
            out += "\n";
            count = 0;
        }
    }
    if count != 0 {
        out += "\n";
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str) -> String {
        format!(
            "{}/../../fixtures/mcnp/wwinp/{name}",
            env!("CARGO_MANIFEST_DIR")
        )
    }

    #[test]
    fn neutron_header_matches_oracle() {
        let w = Wwinp::from_file(fixture("mcnp_wwinp_wwinp_n.txt")).unwrap();
        assert_eq!(w.ni, 1);
        assert_eq!(w.nr, 10);
        assert_eq!(w.ne, vec![7]);
        assert_eq!(w.nf, [15, 8, 6]);
        assert_eq!(w.nft, 720);
        assert_eq!(w.origin, [-100.0, -100.0, -100.0]);
        assert_eq!(w.nc, [5, 3, 1]);
        assert_eq!(w.nwg, 1);
    }

    #[test]
    fn neutron_mesh_definition_matches_oracle() {
        let w = Wwinp::from_file(fixture("mcnp_wwinp_wwinp_n.txt")).unwrap();
        assert_eq!(w.cm[0], vec![-99.0, -97.0, 97.0, 99.0, 100.0]);
        assert_eq!(w.fm[0], vec![1.0, 1.0, 11.0, 1.0, 1.0]);
        assert_eq!(w.fm[1], vec![1.0, 3.0, 4.0]);
        assert_eq!(w.fm[2], vec![6.0]);
        // Interpolated bounds (oracle values).
        assert!((w.bounds[1][2] - (-13.333333333333336)).abs() < 1e-9);
        assert_eq!(w.bounds[0].len(), 16);
        assert_eq!(w.bounds[1].len(), 9);
        assert_eq!(w.bounds[2].len(), 7);
        // First x bounds exact.
        assert_eq!(&w.bounds[0][..3], &[-100.0, -99.0, -97.0]);
        assert_eq!(*w.bounds[0].last().unwrap(), 100.0);
    }

    #[test]
    fn neutron_energy_and_ww_shapes() {
        let w = Wwinp::from_file(fixture("mcnp_wwinp_wwinp_n.txt")).unwrap();
        let exp_e = [0.1, 0.14678, 0.21544, 0.31623, 0.46416, 0.68129, 1.0000];
        assert_eq!(w.e[0], exp_e.to_vec());
        assert_eq!(w.ww.len(), 1); // one particle
        assert_eq!(w.ww[0].len(), 7); // groups
        for g in &w.ww[0] {
            assert_eq!(g.len(), 720); // nft
        }
        // Per-ve column accessor.
        let col = w.ww_column(0, 0);
        assert_eq!(col.len(), 7);
        assert_eq!(col, w.ww[0].iter().map(|r| r[0]).collect::<Vec<_>>());
    }

    #[test]
    fn photon_fixture_parses() {
        // Photon-only WWINPs still declare two particle slots with zero
        // neutron groups.
        let w = Wwinp::from_file(fixture("mcnp_wwinp_wwinp_p.txt")).unwrap();
        assert_eq!(w.ni, 2);
        assert_eq!(w.ne, vec![0, 7]);
        assert_eq!(w.ww.len(), 1); // only photons present
        assert_eq!(w.e.len(), 1);
        assert_eq!(w.ww[0].len(), 7);
    }

    #[test]
    fn np_fixture_has_two_particles() {
        let w = Wwinp::from_file(fixture("mcnp_wwinp_wwinp_np.txt")).unwrap();
        assert_eq!(w.ni, 2);
        assert_eq!(w.ne.len(), 2);
        assert_eq!(w.ww.len(), 2);
        assert_eq!(w.e.len(), 2);
    }

    #[test]
    fn round_trip_preserves_everything() {
        let original = Wwinp::from_file(fixture("mcnp_wwinp_wwinp_n.txt")).unwrap();
        let text = original.to_text().unwrap();
        let reparsed = Wwinp::parse(&text).unwrap();
        assert_eq!(original, reparsed);

        let np = Wwinp::from_file(fixture("mcnp_wwinp_wwinp_np.txt")).unwrap();
        let reparsed_np = Wwinp::parse(&np.to_text().unwrap()).unwrap();
        assert_eq!(np, reparsed_np);
    }

    #[test]
    fn writer_uses_python_style_exponents() {
        let s = fmt13(-100.0);
        assert_eq!(s.trim(), "-1.00000E+02");
        assert_eq!(s.len(), 13);
        assert_eq!(fmt13(15.0).trim(), "1.50000E+01");
        assert_eq!(wrap6(&["x".to_string()]), "x\n");
        let padded = wrap6(&(0..6).map(|_| fmt13(1.0)).collect::<Vec<_>>());
        assert_eq!(padded.lines().count(), 1);
    }

    #[test]
    fn cylindrical_rejected() {
        let text = "         1         1         1        16                     06/23/13 16:49:26 \n         7\n";
        assert!(matches!(Wwinp::parse(text), Err(Error::Unsupported(_))));
    }

    #[test]
    fn truncated_block2_errors() {
        let text = "         1         1         1        10                     06/23/13 16:49:26 \n         7\n   15.000       8.0000       6.0000      -100.00      -100.00      -100.00    \n   5.0000       3.0000       1.0000       1.0000    \n";
        assert!(matches!(
            Wwinp::parse(text),
            Err(Error::BadStructure(m)) if m.contains("EOF")
        ));
    }
}
