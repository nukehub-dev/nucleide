//! USRBIN tally parsing from FLUKA `.lis` output files,
//! validated against the vendored single/multiple/degenerate fixtures.
//!
//! # Format notes: text, not binary records
//!
//! A FLUKA `.lis` container mixes Fortran formatted and unformatted
//! records elsewhere, but the USRBIN section is fully text:
//! each tally starts at a page-break line whose first character is `'1'`,
//! followed by a quoted header line, three `X/Y/Z coordinate:` lines with
//! 11 whitespace-separated fields (min at index 3, max at 5, bin count at
//! 7, bin width at 10), and two equal-length blocks of Fortran `e11.4`
//! floats — track-length binned data first, then percentage errors. The
//! data blocks may be preceded by `accurate deposition along the tracks
//! requested` / `this is a track-length binning` banner lines; their exact
//! skip sequence is reproduced. Reading is therefore purely
//! line-oriented with no binary record unwrapping, matching upstream.
//!
//! # Data layout
//!
//! Volume-element ordering follows the file's `A(ix,iy,iz)` row order:
//! flat index `ve = (i * ny + j) * nz + k`, i.e. x slowest → z fastest,
//! matching the legacy zyx fill order. Meshing layers are out of scope;
//! results land in [`UsrbinTally`] vectors.

use std::fmt;
use std::path::Path;

/// Coordinate system declared by a tally header.
///
/// Only Cartesian is supported; anything else (`R-Z`, `R-Phi-Z`,
/// user-defined) raises.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoordSys {
    /// Plain x-y-z binning.
    Cartesian,
}

impl fmt::Display for CoordSys {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CoordSys::Cartesian => f.write_str("Cartesian"),
        }
    }
}

/// Errors raised while reading USRBIN output.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// Underlying file access failed.
    Io(String),
    /// No page-break-delimited USRBIN blocks were found.
    ///
    /// Unlike lenient readers that silently yield an empty tally collection,
    /// this is an error.
    NoTallies,
    /// Tally header line malformed (missing quoted name/particle fields).
    BadHeader(String),
    /// Non-Cartesian coordinate system in a tally header.
    NotCartesian(String),
    /// X/Y/Z coordinate line malformed or missing fields.
    BadDimensionLine(String),
    /// Bin-count product overflowed while sizing the data block.
    TooManyBins(String),
    /// File ended before a tally block was complete.
    Truncated {
        /// Which section was being read when EOF hit.
        context: &'static str,
    },
    /// A numeric field failed to parse.
    BadNumber {
        /// Which field was being parsed.
        context: &'static str,
        /// Offending token.
        text: String,
    },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(m) => write!(f, "io error: {m}"),
            Error::NoTallies => write!(f, "no USRBIN blocks found"),
            Error::BadHeader(l) => write!(f, "malformed USRBIN header `{l}`"),
            Error::NotCartesian(c) => {
                write!(
                    f,
                    "only cartesian coordinate system currently supported, got `{c}`"
                )
            }
            Error::BadDimensionLine(l) => write!(f, "malformed coordinate line `{l}`"),
            Error::TooManyBins(m) => write!(f, "bin counts overflow: {m}"),
            Error::Truncated { context } => {
                write!(f, "truncated file while reading {context}")
            }
            Error::BadNumber { context, text } => {
                write!(f, "cannot parse {context} from `{text}`")
            }
        }
    }
}

impl std::error::Error for Error {}

/// Binning description of one grid axis as written by the header:
/// `[min, max]` bounds, bin count along the axis, and the uniform width.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DimInfo {
    /// Lower bound of the first bin.
    pub min: f64,
    /// Upper bound of the last bin.
    pub max: f64,
    /// Number of bins along the axis.
    pub bins: usize,
    /// Width of each evenly spaced bin.
    pub width: f64,
}

impl DimInfo {
    /// Vertex bounds `min + i * width` for `i in 0..=bins`.
    ///
    /// Mesh bound generation; note the values come straight from
    /// the header's min/width rather than being rescaled to hit `max`.
    pub fn bounds(&self) -> Vec<f64> {
        (0..=self.bins)
            .map(|i| self.min + i as f64 * self.width)
            .collect()
    }
}

/// One USRBIN detector: header metadata, mesh bounds, and the paired
/// track-length / percentage-error data blocks.
///
/// Tally metadata (`coord_sys`, `name`, `particle`,
/// `x_bounds`, ...); the MOAB tag objects become plain [`Vec`]s plus the
/// [`UsrbinTally::part_data_tag`] / [`UsrbinTally::error_data_tag`] names.
#[derive(Debug, Clone, PartialEq)]
pub struct UsrbinTally {
    /// Declared coordinate system (only Cartesian is accepted).
    pub coord_sys: CoordSys,
    /// User-defined detector name from the quoted header field.
    pub name: String,
    /// Generalized-particle number code as a string (kept textual as written).
    pub particle: String,
    /// X-axis header information.
    pub x_info: DimInfo,
    /// Y-axis header information.
    pub y_info: DimInfo,
    /// Z-axis header information.
    pub z_info: DimInfo,
    /// Mesh vertex locations along x (`nx + 1` entries).
    pub x_bounds: Vec<f64>,
    /// Mesh vertex locations along y (`ny + 1` entries).
    pub y_bounds: Vec<f64>,
    /// Mesh vertex locations along z (`nz + 1` entries).
    pub z_bounds: Vec<f64>,
    /// Track-length binned tally data, x slowest → z fastest.
    pub part_data: Vec<f64>,
    /// Percentage error data, same ordering as [`UsrbinTally::part_data`].
    pub error_data: Vec<f64>,
}

impl UsrbinTally {
    /// `[nx, ny, nz]` cell counts.
    pub fn dims(&self) -> [usize; 3] {
        [self.x_info.bins, self.y_info.bins, self.z_info.bins]
    }

    /// Total number of volume elements (`nx * ny * nz`).
    pub fn num_ves(&self) -> usize {
        self.dims().iter().product()
    }

    /// Flat volume-element index for logical cell `(i, j, k)`; x slowest,
    /// z fastest, matching file order.
    pub fn ve_index(&self, i: usize, j: usize, k: usize) -> usize {
        let [_, ny, nz] = self.dims();
        (i * ny + j) * nz + k
    }

    /// Tag name for the data block (`"part_data_{particle}"`).
    pub fn part_data_tag(&self) -> String {
        format!("part_data_{}", self.particle)
    }

    /// Tag name for the error block (`"error_data_{particle}"`).
    pub fn error_data_tag(&self) -> String {
        format!("error_data_{}", self.particle)
    }
}

/// Read and parse all USRBIN tallies from a `.lis` file on disk.
pub fn read_usrbin_file(path: impl AsRef<Path>) -> Result<Vec<UsrbinTally>, Error> {
    let path = path.as_ref();
    let text = std::fs::read_to_string(path)
        .map_err(|e| Error::Io(format!("{}: {}", path.display(), e)))?;
    parse_usrbin(&text)
}

/// Parse all USRBIN tallies from `.lis` text in memory.
pub fn parse_usrbin(text: &str) -> Result<Vec<UsrbinTally>, Error> {
    let mut reader = Reader {
        lines: text.lines(),
    };
    let mut tallies = Vec::new();

    // A new tally begins at every page-break line whose first
    // character is '1'; anything else (including EOF) stops the scan.
    let mut current = reader.readline();
    while let Some(line) = current {
        if !line.starts_with('1') {
            break;
        }
        tallies.push(parse_tally(&mut reader)?);
        current = reader.readline();
    }

    if tallies.is_empty() {
        return Err(Error::NoTallies);
    }
    Ok(tallies)
}

/// Line source reproducing Fortran-file `readline()` semantics.
struct Reader<'a> {
    lines: std::str::Lines<'a>,
}

impl<'a> Reader<'a> {
    fn readline(&mut self) -> Option<&'a str> {
        self.lines.next()
    }

    fn expect(&mut self, context: &'static str) -> Result<&'a str, Error> {
        self.readline().ok_or(Error::Truncated { context })
    }
}

fn parse_tally(reader: &mut Reader<'_>) -> Result<UsrbinTally, Error> {
    // Header: `   Cartesian binning n.   1  "single_n  " , generalized
    // particle n.    8` — three double-quote-separated segments carry the
    // coordinate system, tally name, and trailing particle number.
    let header = reader.expect("tally header")?;
    let segs: Vec<&str> = header.split('"').collect();
    if segs.len() != 3 {
        return Err(Error::BadHeader(header.trim_end().to_string()));
    }
    let coord_token = segs[0]
        .split_whitespace()
        .next()
        .ok_or_else(|| Error::BadHeader(header.trim_end().to_string()))?;
    let coord_sys = match coord_token {
        "Cartesian" => CoordSys::Cartesian,
        other => return Err(Error::NotCartesian(other.to_string())),
    };
    let name = segs[1].trim().to_string();
    let particle = segs[2]
        .split_whitespace()
        .last()
        .ok_or_else(|| Error::BadHeader(header.trim_end().to_string()))?
        .to_string();

    let x_info = parse_dim(reader.expect("X coordinate line")?)?;
    let y_info = parse_dim(reader.expect("Y coordinate line")?)?;
    let z_info = parse_dim(reader.expect("Z coordinate line")?)?;

    // Banner: "Data follow in a matrix A(ix,iy,iz), format (...)".
    reader.expect("data-follows banner")?;

    // Preamble handling: two unconditional reads (the first
    // is discarded), then optional skips over the track-length banners.
    reader.expect("separator after banner")?;
    let mut cur = reader.expect("preamble or data line")?;
    if cur.contains("accurate deposition") {
        cur = reader.expect("post-deposition line")?;
    }
    if cur.contains("track-length binning") {
        cur = reader.expect("first data line")?;
    }

    let [nx, ny, nz] = [x_info.bins, y_info.bins, z_info.bins];
    let num_volume_element = nx
        .checked_mul(ny)
        .and_then(|v| v.checked_mul(nz))
        .ok_or_else(|| Error::TooManyBins(format!("{name}: {nx}x{ny}x{nz}")))?;

    let mut part_data = Vec::with_capacity(num_volume_element.min(4096));
    read_block(cur, &mut part_data, reader, num_volume_element)?;

    // Skip blank / "Percentage errors follow..." / blank lines.
    for _ in 0..3 {
        reader.expect("error-block separator")?;
    }
    let mut error_data = Vec::with_capacity(num_volume_element.min(4096));
    while error_data.len() < num_volume_element {
        let line = reader.expect("usrbin error datum")?;
        read_floats(line, &mut error_data)?;
    }

    Ok(UsrbinTally {
        coord_sys,
        name,
        particle,
        x_bounds: x_info.bounds(),
        y_bounds: y_info.bounds(),
        z_bounds: z_info.bounds(),
        x_info,
        y_info,
        z_info,
        part_data,
        error_data,
    })
}

fn read_block(
    first: &str,
    into: &mut Vec<f64>,
    reader: &mut Reader<'_>,
    expected_len: usize,
) -> Result<(), Error> {
    read_floats(first, into)?;
    while into.len() < expected_len {
        let line = reader.expect("usrbin datum")?;
        read_floats(line, into)?;
    }
    Ok(())
}

/// Parse one `X/Y/Z coordinate:` line: min at token 3, max at 5, bin
/// count at 7, width at 10.
fn parse_dim(line: &str) -> Result<DimInfo, Error> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 11 {
        return Err(Error::BadDimensionLine(line.trim_end().to_string()));
    }
    let bad = |text: &str| Error::BadNumber {
        context: "coordinate field",
        text: text.to_string(),
    };
    Ok(DimInfo {
        min: tokens[3].parse::<f64>().map_err(|_| bad(tokens[3]))?,
        max: tokens[5].parse::<f64>().map_err(|_| bad(tokens[5]))?,
        bins: tokens[7].parse::<usize>().map_err(|_| bad(tokens[7]))?,
        width: tokens[10].parse::<f64>().map_err(|_| bad(tokens[10]))?,
    })
}

fn read_floats(line: &str, into: &mut Vec<f64>) -> Result<(), Error> {
    for token in line.split_whitespace() {
        let value = token.parse::<f64>().map_err(|_| Error::BadNumber {
            context: "usrbin datum",
            text: token.to_string(),
        })?;
        into.push(value);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture_path(name: &str) -> String {
        format!("{}/../../fixtures/fluka/{name}", env!("CARGO_MANIFEST_DIR"))
    }

    fn load(name: &str) -> Vec<UsrbinTally> {
        read_usrbin_file(fixture_path(name)).unwrap()
    }

    const SINGLE_PART: [f64; 27] = [
        1.0984e-02, 4.1051e-03, 1.0636e-03, 2.1837e-02, 5.5610e-03, 1.9119e-03, 1.0971e-02,
        3.3943e-03, 1.2456e-03, 1.6615e-02, 2.9501e-03, 7.4597e-04, 1.0395e-01, 6.1186e-03,
        1.4997e-03, 1.7421e-02, 3.0824e-03, 7.3878e-04, 1.8097e-02, 5.2532e-03, 2.1572e-03,
        1.0465e-01, 6.2611e-03, 1.8829e-03, 1.7323e-02, 5.5092e-03, 2.1418e-03,
    ];
    const SINGLE_ERROR: [f64; 27] = [
        5.0179e+00, 1.6521e+01, 1.3973e+01, 4.2025e+00, 8.1766e+00, 1.1465e+01, 7.2005e+00,
        1.0479e+01, 1.5640e+01, 5.5994e+00, 1.3275e+01, 2.7617e+01, 7.3788e-01, 6.7200e+00,
        1.9092e+01, 7.3670e+00, 1.3018e+01, 2.8866e+01, 5.7221e+00, 1.5916e+01, 2.6001e+01,
        8.3490e-01, 1.6715e+01, 1.2759e+01, 5.0763e+00, 1.1420e+01, 1.0040e+01,
    ];
    const MULTI_P_PART: [f64; 27] = [
        7.5083e-04, 1.7570e-04, 3.3361e-05, 1.1232e-03, 3.4735e-04, 1.5816e-04, 6.2264e-04,
        2.3071e-04, 8.3469e-05, 1.6700e-03, 4.1785e-04, 7.6990e-05, 3.3842e-03, 9.2931e-04,
        2.4958e-04, 1.0121e-03, 2.7993e-04, 6.1043e-05, 7.7401e-04, 3.2480e-04, 9.3145e-06,
        1.4245e-03, 4.3352e-04, 1.7392e-04, 7.3166e-04, 2.4210e-04, 1.4804e-04,
    ];
    const MULTI_P_ERROR: [f64; 27] = [
        2.2149e+01, 7.4509e+01, 1.0000e+02, 2.4621e+01, 4.6383e+01, 3.3621e+01, 2.1616e+01,
        7.5885e+01, 1.0000e+02, 2.0067e+01, 3.3654e+01, 6.1265e+01, 1.8407e+01, 1.6239e+01,
        5.2119e+01, 1.5791e+01, 3.8452e+01, 1.0000e+02, 7.6577e+00, 3.5290e+01, 1.0000e+02,
        8.3702e+00, 5.3283e+01, 6.2602e+01, 1.1655e+01, 6.2289e+01, 6.7541e+01,
    ];

    #[test]
    fn single_header_and_bounds() {
        let tallies = load("fluka_usrbin_single.lis");
        assert_eq!(tallies.len(), 1);
        let t = &tallies[0];
        assert_eq!(t.coord_sys, CoordSys::Cartesian);
        assert_eq!(CoordSys::Cartesian.to_string(), "Cartesian");
        assert_eq!(t.name, "single_n");
        assert_eq!(t.particle, "8");
        assert_eq!(t.x_bounds, vec![-3.0, 0.0, 3.0, 6.0]);
        assert_eq!(t.y_bounds, vec![-3.0, -1.0, 1.0, 3.0]);
        assert_eq!(t.z_bounds, vec![-3.0, -2.0, -1.0, 0.0]);
        assert_eq!(
            t.x_info,
            DimInfo {
                min: -3.0,
                max: 6.0,
                bins: 3,
                width: 3.0
            }
        );
        assert_eq!(
            t.y_info,
            DimInfo {
                min: -3.0,
                max: 3.0,
                bins: 3,
                width: 2.0
            }
        );
        assert_eq!(
            t.z_info,
            DimInfo {
                min: -3.0,
                max: 0.0,
                bins: 3,
                width: 1.0
            }
        );
        assert_eq!(t.dims(), [3, 3, 3]);
        assert_eq!(t.num_ves(), 27);
        assert_eq!(t.part_data_tag(), "part_data_8");
        assert_eq!(t.error_data_tag(), "error_data_8");
    }

    #[test]
    fn single_part_data_matches_oracle() {
        let t = &load("fluka_usrbin_single.lis")[0];
        assert_eq!(t.part_data.len(), 27);
        assert_eq!(t.part_data, SINGLE_PART.to_vec());
        assert_eq!(t.part_data[t.ve_index(0, 0, 0)], SINGLE_PART[0]);
        assert_eq!(t.part_data[t.ve_index(2, 2, 2)], SINGLE_PART[26]);
        assert_eq!(t.part_data[t.ve_index(1, 1, 1)], SINGLE_PART[13]);
    }

    #[test]
    fn single_error_data_matches_oracle() {
        let t = &load("fluka_usrbin_single.lis")[0];
        assert_eq!(t.error_data.len(), 27);
        assert_eq!(t.error_data, SINGLE_ERROR.to_vec());
    }

    #[test]
    fn multiple_yields_two_tallies_with_oracle_values() {
        let tallies = load("fluka_usrbin_multiple.lis");
        assert_eq!(tallies.len(), 2);

        let p = &tallies[0];
        assert_eq!((p.name.as_str(), p.particle.as_str()), ("multi_p", "7"));
        assert_eq!(p.part_data_tag(), "part_data_7");
        assert_eq!(p.x_bounds, vec![-3.0, 0.0, 3.0, 6.0]);
        assert_eq!(p.y_bounds, vec![-3.0, -1.0, 1.0, 3.0]);
        assert_eq!(p.z_bounds, vec![-3.0, -2.0, -1.0, 0.0]);
        assert_eq!(p.part_data, MULTI_P_PART.to_vec());
        assert_eq!(p.error_data, MULTI_P_ERROR.to_vec());

        let n = &tallies[1];
        assert_eq!((n.name.as_str(), n.particle.as_str()), ("multi_n", "8"));
        assert_eq!(n.part_data_tag(), "part_data_8");
        assert_eq!(n.part_data, SINGLE_PART.to_vec());
        assert_eq!(n.error_data, SINGLE_ERROR.to_vec());
    }

    #[test]
    fn degenerate_tallies_match_oracle_bounds_and_data() {
        let tallies = load("fluka_usrbin_degenerate.lis");
        assert_eq!(tallies.len(), 3);

        let d1 = &tallies[0];
        assert_eq!(d1.name, "degen1");
        assert_eq!(d1.x_bounds, vec![-3.0, 0.0, 3.0, 6.0]);
        assert_eq!(d1.y_bounds, vec![-3.0, 0.0, 3.0]);
        assert_eq!(d1.z_bounds, vec![-3.0, 0.0]);
        assert_eq!(d1.num_ves(), 6);
        assert_eq!(
            d1.part_data,
            [3.5279e-02, 4.7334e-03, 1.4458e-03, 3.6242e-02, 4.6521e-03, 1.5292e-03].to_vec()
        );
        assert_eq!(
            d1.error_data,
            [1.2016e+00, 6.4313e+00, 7.7312e+00, 2.0235e+00, 9.4199e+00, 8.0514e+00].to_vec()
        );

        let d2 = &tallies[1];
        assert_eq!(d2.name, "degen2");
        assert_eq!(d2.x_bounds, vec![-3.0, 1.5, 6.0]);
        assert_eq!(d2.y_bounds, vec![-3.0, 3.0]);
        assert_eq!(d2.z_bounds, vec![-3.0, -2.0, -1.0, 0.0]);
        assert_eq!(
            d2.part_data,
            [1.1543e-02, 2.0295e-03, 3.2603e-02, 1.4229e-03, 3.3492e-02, 2.7923e-03].to_vec()
        );
        assert_eq!(
            d2.error_data,
            [2.7321e+00, 5.2342e+00, 7.4679e-01, 4.2862e+00, 1.3090e+00, 1.4151e+01].to_vec()
        );

        let d3 = &tallies[2];
        assert_eq!(d3.name, "degen3");
        assert_eq!(d3.x_bounds, vec![-3.0, 6.0]);
        assert_eq!(d3.y_bounds, vec![-3.0, -1.0, 1.0, 3.0]);
        assert_eq!(d3.z_bounds, vec![-3.0, -1.5, 0.0]);
        assert_eq!(
            d3.part_data,
            [5.8037e-03, 1.3260e-02, 5.6046e-03, 7.9677e-03, 4.3111e-02, 8.1349e-03].to_vec()
        );
        assert_eq!(
            d3.error_data,
            [6.1913e+00, 2.3684e+00, 4.6124e+00, 3.2523e+00, 1.3714e+00, 4.3161e+00].to_vec()
        );
    }

    #[test]
    fn ve_index_is_row_major_x_slowest() {
        let t = &load("fluka_usrbin_single.lis")[0];
        assert_eq!(t.ve_index(0, 0, 0), 0);
        assert_eq!(t.ve_index(0, 0, 1), 1);
        assert_eq!(t.ve_index(0, 1, 0), 3);
        assert_eq!(t.ve_index(1, 0, 0), 9);
        assert_eq!(t.ve_index(2, 2, 2), 26);
    }

    #[test]
    fn truncated_file_reports_context() {
        let full = std::fs::read_to_string(fixture_path("fluka_usrbin_single.lis")).unwrap();
        // Cut mid-data-block: the float loop runs off the end.
        let half = &full[..full.len() / 2];
        assert!(matches!(parse_usrbin(half), Err(Error::Truncated { .. })));
        // Cut right after part_data completes: the error-block separators
        // are missing instead.
        let cut = full.find("Percentage errors").unwrap();
        assert!(matches!(
            parse_usrbin(&full[..cut]),
            Err(Error::Truncated {
                context: "error-block separator"
            })
        ));
        // Cut immediately before the Y line (line start, not the text
        // occurrence): the X line parsed fine, EOF hits next.
        let early_cut = full.find("\n      Y coordinate").unwrap() + 1;
        assert!(matches!(
            parse_usrbin(&full[..early_cut]),
            Err(Error::Truncated {
                context: "Y coordinate line"
            })
        ));
    }

    #[test]
    fn empty_or_blockless_text_errors() {
        assert_eq!(parse_usrbin(""), Err(Error::NoTallies));
        assert_eq!(
            parse_usrbin("just some text\nwithout usrbin pages\n"),
            Err(Error::NoTallies)
        );
    }

    #[test]
    fn non_cartesian_header_rejected() {
        let text = concat!(
            "1\n",
            "   R-Z binning n.   1  \"rz_det    \" , generalized particle n.    8\n",
            "      R coordinate: from -3.0000E+00 to  6.0000E+00 cm,     3 bins ( 3.0000E+00 cm wide)\n",
        );
        match parse_usrbin(text) {
            Err(Error::NotCartesian(c)) => assert_eq!(c, "R-Z"),
            other => panic!("expected NotCartesian, got {other:?}"),
        }
    }

    #[test]
    fn malformed_dimension_line_errors() {
        let good_header =
            "   Cartesian binning n.   1  \"bad_dim    \" , generalized particle n.    8\n";
        let short_line = "      X coordinate: from -3.0000E+00\n";
        let text = format!("1\n{good_header}{short_line}");
        assert!(matches!(
            parse_usrbin(&text),
            Err(Error::BadDimensionLine(_))
        ));

        let bad_number = "      X coordinate: from -3.0000E+00 to  6.0000E+00 cm,     xx bins ( 3.0000E+00 cm wide)\n";
        let text = format!("1\n{good_header}{bad_number}");
        match parse_usrbin(&text) {
            Err(Error::BadNumber { text, .. }) => assert_eq!(text, "xx"),
            other => panic!("expected BadNumber, got {other:?}"),
        }
    }

    #[test]
    fn missing_quote_segments_error() {
        let text = "1\n   Cartesian binning without any quotes here\n";
        assert!(matches!(parse_usrbin(text), Err(Error::BadHeader(_))));
    }

    #[test]
    fn junk_after_last_tally_is_ignored() {
        let mut text = std::fs::read_to_string(fixture_path("fluka_usrbin_single.lis")).unwrap();
        text.push_str("\ntrailing garbage that is not a page break\n");
        let tallies = parse_usrbin(&text).unwrap();
        assert_eq!(tallies.len(), 1);
        assert_eq!(tallies[0].part_data, SINGLE_PART.to_vec());
    }
}
