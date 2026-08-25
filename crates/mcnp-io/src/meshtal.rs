//! MCNP meshtal file parsing (`Meshtal` → native structured tally data),
//! validated against the vendored single- and multiple-tally fixtures.
//! MOAB-backed meshing is out of scope: results land in [`MeshTallyData`]
//! with an explicit, documented cell-ordering convention (x slowest → z
//! fastest, matching the legacy fill order used by mesh tools).
//!
//! Single-energy-group tallies carry no `Total` block in MCNP output; their
//! totals mirror the sole energy group.

use std::collections::BTreeMap;
use std::fmt;
use std::path::Path;

/// Particle identified in a mesh tally header.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParticleKind {
    Neutron,
    Photon,
}

impl ParticleKind {
    fn parse(line: &str) -> Result<Self, Error> {
        if line.contains("neutron") {
            Ok(ParticleKind::Neutron)
        } else if line.contains("photon") {
            Ok(ParticleKind::Photon)
        } else {
            Err(Error::BadParticleLine(line.trim().to_string()))
        }
    }

    /// Short particle letter ('n' or 'p').
    pub fn letter(&self) -> char {
        match self {
            ParticleKind::Neutron => 'n',
            ParticleKind::Photon => 'p',
        }
    }
}

/// Errors raised while parsing meshtal text.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    Io(String),
    /// File did not start with a recognizable MCNP header.
    BadHeader(String),
    /// A `Mesh Tally Number` sub-block was malformed.
    BadTallyBlock(String),
    /// Unknown particle in a tally header.
    BadParticleLine(String),
    /// A numeric field failed to parse.
    BadNumber {
        context: &'static str,
        text: String,
    },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(m) => write!(f, "io error: {m}"),
            Error::BadHeader(m) => write!(f, "malformed meshtal header: {m}"),
            Error::BadTallyBlock(m) => write!(f, "malformed tally block: {m}"),
            Error::BadParticleLine(l) => write!(f, "unknown particle in `{l}`"),
            Error::BadNumber { context, text } => {
                write!(f, "cannot parse {context} from `{text}`")
            }
        }
    }
}

impl std::error::Error for Error {}

fn f64_of(context: &'static str, token: &str) -> Result<f64, Error> {
    token.parse::<f64>().map_err(|_| Error::BadNumber {
        context,
        text: token.to_string(),
    })
}

/// One parsed fmesh4 tally: bounds plus result/error arrays.
///
/// Cell ordering convention: flat volume-element index
/// `ve = (i * ny_cells + j) * nz_cells + k`, i.e. x slowest, z fastest —
/// identical to the row order MCNP writes and to the
/// `structured_iterate_hex("xyz")` assignment order.
#[derive(Debug, Clone, PartialEq)]
pub struct MeshTallyData {
    pub tally_number: u32,
    pub particle: ParticleKind,
    /// True when flux-to-dose conversion factors modified this tally.
    pub dose_response: bool,
    pub x_bounds: Vec<f64>,
    pub y_bounds: Vec<f64>,
    pub z_bounds: Vec<f64>,
    pub e_bounds: Vec<f64>,
    /// Column-name → index mapping from the table header
    /// (`"Rel Error"` normalizes to `"Rel_Error"`).
    pub column_idx: BTreeMap<String, usize>,
    /// `result[ve][energy_group]`.
    pub result: Vec<Vec<f64>>,
    /// `rel_error[ve][energy_group]`.
    pub rel_error: Vec<Vec<f64>>,
    /// Per-cell energy-integrated totals (mirrors the single group when
    /// the file has only one energy bin).
    pub total_result: Vec<f64>,
    pub total_rel_error: Vec<f64>,
}

impl MeshTallyData {
    /// `[nx, ny, nz]` cell counts.
    pub fn dims(&self) -> [usize; 3] {
        [
            self.x_bounds.len() - 1,
            self.y_bounds.len() - 1,
            self.z_bounds.len() - 1,
        ]
    }

    /// Total number of volume elements.
    pub fn num_ves(&self) -> usize {
        let d = self.dims();
        d[0] * d[1] * d[2]
    }

    /// Number of energy groups.
    pub fn num_e_groups(&self) -> usize {
        self.e_bounds.len() - 1
    }

    /// Flat ve index for logical cell `(i, j, k)` (x/y/z indices).
    pub fn ve_index(&self, i: usize, j: usize, k: usize) -> usize {
        let d = self.dims();
        (i * d[1] + j) * d[2] + k
    }

    /// All energy-group results/errors for one cell.
    pub fn cell(&self, i: usize, j: usize, k: usize) -> (&[f64], &[f64]) {
        let ve = self.ve_index(i, j, k);
        (&self.result[ve], &self.rel_error[ve])
    }

    /// Energy-integrated total (result, rel_error) for one cell.
    pub fn cell_total(&self, i: usize, j: usize, k: usize) -> (f64, f64) {
        let ve = self.ve_index(i, j, k);
        (self.total_result[ve], self.total_rel_error[ve])
    }
}

/// Parsed contents of a meshtal file: run metadata plus all tallies.
#[derive(Debug, Clone, PartialEq)]
pub struct Meshtal {
    /// MCNP version string (e.g. `"5.mpi"`).
    pub version: String,
    /// Version date (`ld=` value).
    pub ld: String,
    /// Title card from the input deck.
    pub title: String,
    /// Normalizing history count.
    pub histories: u64,
    /// Tallies keyed by `fmesh4` number.
    pub tallies: BTreeMap<u32, MeshTallyData>,
}

const DOSE_RESPONSE_LINE: &str = "This mesh tally is modified by a dose response function.";

impl Meshtal {
    /// Read and parse a meshtal file from disk.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self, Error> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path)
            .map_err(|e| Error::Io(format!("{}: {}", path.display(), e)))?;
        Meshtal::parse(&text)
    }

    /// Parse meshtal text in memory.
    pub fn parse(text: &str) -> Result<Self, Error> {
        let mut lines = text.lines();

        // Header: version / ld / title / histories.
        let line1 = lines
            .next()
            .ok_or_else(|| Error::BadHeader("empty file".into()))?;
        let t1: Vec<&str> = line1.split_whitespace().collect();
        if t1.len() < 4 {
            return Err(Error::BadHeader(format!(
                "expected MCNP banner, got `{line1}`"
            )));
        }
        let version = t1[2].to_string();
        let ld = t1[3].strip_prefix("ld=").unwrap_or(t1[3]).to_string();
        let title = lines
            .next()
            .ok_or_else(|| Error::BadHeader("missing title line".into()))?
            .trim()
            .to_string();
        let line3 = lines
            .next()
            .ok_or_else(|| Error::BadHeader("missing histories line".into()))?;
        let t3: Vec<&str> = line3.split_whitespace().collect();
        let histories = f64_of(
            "histories",
            t3.last()
                .ok_or_else(|| Error::BadHeader(format!("no histories value in `{line3}`")))?,
        )? as u64;

        let mut tallies = BTreeMap::new();
        let mut it = lines.peekable();
        while let Some(line) = it.next() {
            let w: Vec<&str> = line.split_whitespace().collect();
            if w.len() >= 4 && w[0] == "Mesh" && w[1] == "Tally" && w[2] == "Number" {
                let num_token = w[3];
                let tally_number = num_token.parse::<u32>().map_err(|_| Error::BadNumber {
                    context: "tally number",
                    text: num_token.to_string(),
                })?;
                let tally = parse_tally_block(tally_number, &mut it)?;
                tallies.insert(tally_number, tally);
            }
        }

        if tallies.is_empty() {
            return Err(Error::BadTallyBlock(
                "no `Mesh Tally Number` blocks found".into(),
            ));
        }

        Ok(Meshtal {
            version,
            ld,
            title,
            histories,
            tallies,
        })
    }
}

fn parse_tally_block(
    tally_number: u32,
    it: &mut std::iter::Peekable<std::str::Lines<'_>>,
) -> Result<MeshTallyData, Error> {
    let bad = |m: &str| Error::BadTallyBlock(format!("tally {tally_number}: {m}"));

    // Particle line + optional dose-response line.
    let p_line = it
        .next()
        .ok_or_else(|| bad("truncated before particle line"))?;
    let particle = ParticleKind::parse(p_line)?;
    let dr_line = it
        .next()
        .ok_or_else(|| bad("truncated before response line"))?;
    let dose_response = dr_line.trim() == DOSE_RESPONSE_LINE;

    // Advance to "Tally bin boundaries:".
    loop {
        match it.next() {
            Some(l) if l.trim() == "Tally bin boundaries:" => break,
            Some(_) => {}
            None => return Err(bad("no `Tally bin boundaries:` found")),
        }
    }

    let mut next_bounds = |prefix: &'static str| -> Result<Vec<f64>, Error> {
        let l = it
            .next()
            .ok_or_else(|| bad("truncated in bounds section"))?;
        let w: Vec<&str> = l.split_whitespace().collect();
        let skip = if prefix == "Energy" { 3 } else { 2 };
        let start = w
            .iter()
            .position(|t| *t == prefix)
            .map(|p| p + 1)
            .unwrap_or(skip);
        let start = start.max(skip);
        let vals: Vec<f64> = w[start..]
            .iter()
            .map(|t| f64_of("bounds", t))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(vals)
    };
    let x_bounds = next_bounds("direction:")?;
    let y_bounds = next_bounds("direction:")?;
    let z_bounds = next_bounds("direction:")?;
    let e_bounds = next_bounds("Energy")?;

    // Blank separator, then column header line.
    while matches!(it.peek(), Some(l) if l.trim().is_empty()) {
        it.next();
    }
    let header = it
        .next()
        .ok_or_else(|| bad("truncated before column headers"))?;
    let normalized = header.replace("Rel ", "Rel_").replace("Rslt * ", "Rslt_*_");
    let column_idx: BTreeMap<String, usize> = normalized
        .split_whitespace()
        .enumerate()
        .map(|(i, name)| (name.to_string(), i))
        .collect();
    let col = |name: &str| -> Result<usize, Error> {
        column_idx
            .get(name)
            .copied()
            .ok_or_else(|| bad(&format!("column `{name}` missing from header")))
    };
    let c_result = col("Result")?;
    let c_rel = col("Rel_Error")?;

    let dims_ok =
        !x_bounds.is_empty() && !y_bounds.is_empty() && !z_bounds.is_empty() && e_bounds.len() >= 2;
    if !dims_ok {
        return Err(bad("degenerate bounds"));
    }
    let n_cells = (x_bounds.len() - 1) * (y_bounds.len() - 1) * (z_bounds.len() - 1);
    let num_egs = e_bounds.len() - 1;

    // Data rows: one block per energy group, then totals when grouped.
    let mut read_block = |which: &'static str| -> Result<(Vec<f64>, Vec<f64>), Error> {
        let mut res = Vec::with_capacity(n_cells);
        let mut err = Vec::with_capacity(n_cells);
        for _ in 0..n_cells {
            let l = it.next().ok_or_else(|| bad("truncated in data rows"))?;
            let w: Vec<&str> = l.split_whitespace().collect();
            if c_result >= w.len() || c_rel >= w.len() {
                return Err(bad("data row too short"));
            }
            res.push(f64_of(which, w[c_result])?);
            err.push(f64_of("rel error", w[c_rel])?);
        }
        Ok((res, err))
    };

    let mut result = Vec::with_capacity(n_cells);
    let mut rel_error = Vec::with_capacity(n_cells);
    for _ in 0..num_egs {
        let (r, e) = read_block("result")?;
        result.push(r);
        rel_error.push(e);
    }
    // Transpose to [ve][group].
    let result_t: Vec<Vec<f64>> = (0..n_cells)
        .map(|ve| (0..num_egs).map(|eg| result[eg][ve]).collect())
        .collect();
    let rel_error_t: Vec<Vec<f64>> = (0..n_cells)
        .map(|ve| (0..num_egs).map(|eg| rel_error[eg][ve]).collect())
        .collect();

    let (total_result, total_rel_error) = if num_egs > 1 {
        read_block("total result")?
    } else {
        // Totals are the lone energy group across all cells.
        let tot_r = result_t.iter().map(|v| v[0]).collect();
        let tot_e = rel_error_t.iter().map(|v| v[0]).collect();
        (tot_r, tot_e)
    };

    Ok(MeshTallyData {
        tally_number,
        particle,
        dose_response,
        x_bounds,
        y_bounds,
        z_bounds,
        e_bounds,
        column_idx,
        result: result_t,
        rel_error: rel_error_t,
        total_result,
        total_rel_error,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str) -> String {
        format!(
            "{}/../../fixtures/mcnp/meshtal/{name}",
            env!("CARGO_MANIFEST_DIR")
        )
    }

    #[test]
    fn single_header_metadata() {
        let m = Meshtal::from_file(fixture("mcnp_meshtal_single_meshtal.txt")).unwrap();
        assert_eq!(m.version, "5.mpi");
        assert_eq!(m.ld, "09282010");
        assert_eq!(m.title, "Input file to general test meshtal file");
        assert_eq!(m.histories, 100000);
    }

    #[test]
    fn single_tally_bounds_match_oracle() {
        let m = Meshtal::from_file(fixture("mcnp_meshtal_single_meshtal.txt")).unwrap();
        let t = &m.tallies[&4];
        assert_eq!(t.tally_number, 4);
        assert_eq!(t.particle, ParticleKind::Neutron);
        assert!(t.dose_response);
        assert_eq!(t.x_bounds, vec![-200.00, -66.67, 66.67, 200.00]);
        assert_eq!(
            t.y_bounds,
            vec![-200.00, -120.00, -40.00, 40.00, 120.00, 200.00]
        );
        assert_eq!(t.z_bounds, vec![-200.00, -50.00, 100.00, 200.00]);
        assert_eq!(t.e_bounds, vec![0.00e0, 1.00e-1, 2.00e-1, 1.00e0]);
        assert_eq!(t.dims(), [3, 5, 3]);
        assert_eq!(t.num_ves(), 45);
        assert_eq!(t.num_e_groups(), 3);
    }

    #[test]
    fn single_cell_values_in_file_order() {
        let m = Meshtal::from_file(fixture("mcnp_meshtal_single_meshtal.txt")).unwrap();
        let t = &m.tallies[&4];
        // First data row of group 0.
        let (r, e) = t.cell(0, 0, 0);
        assert_eq!(r[0], 4.96471e-9);
        assert_eq!(e[0], 1.98750e-1);
        // First Total row corresponds to the same first cell.
        let (tr, te) = t.cell_total(0, 0, 0);
        assert_eq!(tr, 1.91370e-7);
        assert_eq!(te, 4.29395e-2);
        // Last data row overall (last cell of last group, from file tail).
        let last = t.num_ves() - 1;
        assert_eq!(t.result[last][2], 1.43674e-7);
        assert_eq!(t.rel_error[last][2], 5.51659e-2);
        // Last Total row (same cell) is the final line of the file.
        assert_eq!(t.total_result[last], 1.61927e-7);
        assert_eq!(t.total_rel_error[last], 5.35151e-2);
    }

    #[test]
    fn multiple_header_and_tally_set() {
        let m = Meshtal::from_file(fixture("mcnp_meshtal_multiple_meshtal.txt")).unwrap();
        assert_eq!(m.version, "5");
        let keys: Vec<u32> = m.tallies.keys().copied().collect();
        assert_eq!(keys, vec![4, 14, 24, 34]);
        assert_eq!(m.tallies[&4].particle, ParticleKind::Neutron);
        assert_eq!(m.tallies[&14].particle, ParticleKind::Neutron);
        assert_eq!(m.tallies[&24].particle, ParticleKind::Photon);
        assert_eq!(m.tallies[&34].particle, ParticleKind::Photon);
    }

    #[test]
    fn multiple_energy_group_counts_and_totals() {
        let m = Meshtal::from_file(fixture("mcnp_meshtal_multiple_meshtal.txt")).unwrap();
        // Tally 4: 6 groups with explicit Total block.
        let t4 = &m.tallies[&4];
        assert_eq!(t4.num_e_groups(), 6);
        assert_eq!(t4.dims(), [1, 4, 5]);
        assert_eq!(t4.result[0][0], 6.00211e3);
        assert_eq!(t4.rel_error[0][0], 1.27494e-1);
        let (tr, te) = t4.cell_total(0, 1, 0);
        assert_eq!(tr, 3.24329e5);
        assert_eq!(te, 2.32000e-2);

        // Tally 14: single group (bounds 0 .. 1E36) — no Total block;
        // totals must mirror the lone group.
        let t14 = &m.tallies[&14];
        assert_eq!(t14.num_e_groups(), 1);
        assert!(!t14.dose_response);
        assert_eq!(t14.e_bounds.last(), Some(&1.00e36));
        assert_eq!(
            t14.total_result,
            t14.result.iter().map(|v| v[0]).collect::<Vec<_>>()
        );

        // Tally 24/34 are photon tallies; 34 single-group.
        assert_eq!(m.tallies[&34].num_e_groups(), 1);
        assert_eq!(m.tallies[&34].particle.letter(), 'p');
    }

    #[test]
    fn column_order_map_parsed() {
        let m = Meshtal::from_file(fixture("mcnp_meshtal_single_meshtal.txt")).unwrap();
        let t = &m.tallies[&4];
        assert_eq!(t.column_idx["Result"], 4);
        assert_eq!(t.column_idx["Rel_Error"], 5);
        assert_eq!(t.column_idx["Energy"], 0);
    }

    #[test]
    fn ve_index_mapping_matches_xyz_iteration() {
        let m = Meshtal::from_file(fixture("mcnp_meshtal_single_meshtal.txt")).unwrap();
        let t = &m.tallies[&4];
        // z fastest: (0,0,1) is the second row of the file.
        assert_eq!(t.ve_index(0, 0, 1), 1);
        assert_eq!(t.ve_index(0, 1, 0), 3);
        assert_eq!(t.ve_index(1, 0, 0), 15);
        // Second data row's values land at (0,0,1).
        assert_eq!(t.result[t.ve_index(0, 0, 1)][0], 7.73879e-9);
    }

    #[test]
    fn truncated_data_rows_error() {
        let text = "\
mcnp   version 5 ld=010101  probid = 01/01/01
title
 Number of histories used for normalizing tallies = 10.00

 Mesh Tally Number         4
 This is a neutron mesh tally.

 Tally bin boundaries:
    X direction:   0.00   1.00
    Y direction:   0.00   1.00
    Z direction:   0.00   1.00
    Energy bin boundaries: 0.00E+00 1.00E+00

   Energy         X         Y         Z     Result     Rel Error
";
        assert!(matches!(
            Meshtal::parse(text),
            Err(Error::BadTallyBlock(m)) if m.contains("truncated")
        ));
    }

    #[test]
    fn unknown_particle_errors() {
        let text = "\
mcnp   version 5 ld=010101  probid = 01/01/01
title
 Number of histories used for normalizing tallies = 10.00

 Mesh Tally Number         4
 This is an electron mesh tally.
";
        assert!(matches!(
            Meshtal::parse(text),
            Err(Error::BadParticleLine(_))
        ));
    }

    #[test]
    fn no_tallies_errors() {
        let text = "\
mcnp   version 5 ld=010101  probid = 01/01/01
title
 Number of histories used for normalizing tallies = 10.00
";
        assert!(matches!(Meshtal::parse(text), Err(Error::BadTallyBlock(_))));
    }
}
