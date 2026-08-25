//! MCNP input-deck generation from a structured mesh.
//!
//! Validated against the reference 3x2x1 voxel grid oracle
//! (`[0,1,2,3]x[0,1,2]x[0,1]`). The API takes explicit bounds plus a
//! per-voxel `(material name, density)` table; everything downstream
//! (numbering, ordering, text layout) follows the legacy convention:
//!
//! - Cells are enumerated x slowest -> z fastest (`(i * ny + j) * nz + k`,
//!   the same convention as [`crate::meshtal`]), numbered from 1. Each cell's
//!   material number equals its cell number.
//! - Surfaces are numbered one per axis boundary plane in x, y, z order; each
//!   is an axis-aligned `px`/`py`/`pz` plane.
//! - A final graveyard cell complements the mesh bounding box.
//! - Material cards carry the upstream comment block (`C name: ...`,
//!   `C density = %.5f`) followed by `m<voxel+1>`.
//!
//! Deviations from upstream (documented, driven by the reduced inputs):
//! - No nuclide fraction lines: compositions are not part of the input
//!   contract, so [`FracType`] is recorded but cannot alter output yet.
//! - A `None` voxel reproduces the auto-created default-material
//!   semantics: density falls back to `-1.0` on the cell card and the
//!   material card degenerates to a bare `m<N>` line (no name metadata).
//! - Densities equal to exactly `-1.0` suppress the `C density =` comment,
//!   matching the upstream sentinel check.

use std::fmt::Write as _;

/// Fraction type used for material definitions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FracType {
    /// Fractions by mass (MCNP mass fractions, printed negative upstream).
    Mass,
    /// Fractions by atom (printed positive upstream).
    Atom,
}

/// Options controlling deck generation.
#[derive(Debug, Clone)]
pub struct DeckOptions {
    /// MCNP title card placed on the first line.
    pub title_card: String,
    /// Mass or atom fractions for material definitions.
    pub frac_type: FracType,
}

impl Default for DeckOptions {
    fn default() -> Self {
        DeckOptions {
            title_card: "Generated from mesh".into(),
            frac_type: FracType::Mass,
        }
    }
}

/// Generate the geometry portion of an MCNP input file (title, cells,
/// surfaces, materials) from an axis-aligned structured mesh.
///
/// * `x_bounds` / `y_bounds` / `z_bounds`: mesh boundary coordinates
///   (`len == voxels + 1` per axis).
/// * `cell_materials`: one entry per voxel, ordered x slowest -> z fastest
///   (`idx = (i * ny + j) * nz + k`, same convention as
///   [`crate::meshtal`]). `None` marks an unassigned voxel (default-material
///   semantics, see module docs). Entries past the table end are treated as
///   unassigned rather than panicking.
pub fn mesh_to_geom(
    x_bounds: &[f64],
    y_bounds: &[f64],
    z_bounds: &[f64],
    cell_materials: &[Option<(String, f64)>],
    opts: &DeckOptions,
) -> String {
    let cell_cards = cell_cards(
        x_bounds.len(),
        y_bounds.len(),
        z_bounds.len(),
        cell_materials,
    );
    let surf_cards = surf_cards(x_bounds, y_bounds, z_bounds);
    let mat_cards = mat_cards(
        x_bounds.len(),
        y_bounds.len(),
        z_bounds.len(),
        cell_materials,
    );

    format!(
        "{}\n{}\n{}\n{}",
        opts.title_card, cell_cards, surf_cards, mat_cards
    )
}

fn voxel_at(cell_materials: &[Option<(String, f64)>], idx: usize) -> Option<&(String, f64)> {
    cell_materials.get(idx).and_then(Option::as_ref)
}

/// Voxel cells followed by the graveyard.
///
/// Surface-number extremes are established from division counts:
/// x surfaces are `1..=nx_planes`, y surfaces follow, then z. A cell spanning
/// voxel `(i, j, k)` (1-based planes) references `i/-i+1`,
/// `j+x_max/-(j+x_max+1)`, `k+y_max/-(k+y_max+1)`.
fn cell_cards(
    nx_planes: usize,
    ny_planes: usize,
    nz_planes: usize,
    cell_materials: &[Option<(String, f64)>],
) -> String {
    let mut cards = String::new();
    let mut count = 1usize;

    let x_min = 1usize;
    let x_max = nx_planes;
    let y_min = x_max + 1;
    let y_max = x_max + ny_planes;
    let z_min = y_max + 1;
    let z_max = y_max + nz_planes;

    let nx = nx_planes.saturating_sub(1);
    let ny = ny_planes.saturating_sub(1);
    let nz = nz_planes.saturating_sub(1);

    for i in 1..=nx {
        for j in 1..=ny {
            for k in 1..=nz {
                // Cell number, mat number, density. Unassigned voxels get
                // Default-material density sentinel (-1.0).
                let idx = (i - 1) * ny * nz + (j - 1) * nz + (k - 1);
                let density = voxel_at(cell_materials, idx).map_or(-1.0, |(_, d)| *d);
                write!(cards, "{count} {count} {} ", py_float(density)).unwrap();
                // x, y, and z surfaces.
                writeln!(
                    cards,
                    "{i} -{ip1} {ysurf} -{ysurf1} {zsurf} -{zsurf1}",
                    ip1 = i + 1,
                    ysurf = j + x_max,
                    ysurf1 = j + x_max + 1,
                    zsurf = k + y_max,
                    zsurf1 = k + y_max + 1,
                )
                .unwrap();
                count += 1;
            }
        }
    }

    // Append graveyard.
    writeln!(
        cards,
        "{count} 0 -{x_min}:{x_max}:-{y_min}:{y_max}:-{z_min}:{z_max}"
    )
    .unwrap();
    cards
}

/// `_mesh_to_surf_cards`: one plane surface per boundary, x then y then z.
fn surf_cards(x_bounds: &[f64], y_bounds: &[f64], z_bounds: &[f64]) -> String {
    let mut cards = String::new();
    let mut count = 1usize;
    for (dim, divs) in [("x", x_bounds), ("y", y_bounds), ("z", z_bounds)] {
        for div in divs {
            writeln!(cards, "{count} p{dim} {}", py_float(*div)).unwrap();
            count += 1;
        }
    }
    cards
}

/// `_mesh_to_mat_cards`: one material card block per voxel, numbered by
/// voxel index + 1 in iteration order.
fn mat_cards(
    nx_planes: usize,
    ny_planes: usize,
    nz_planes: usize,
    cell_materials: &[Option<(String, f64)>],
) -> String {
    let mut cards = String::new();
    let nx = nx_planes.saturating_sub(1);
    let ny = ny_planes.saturating_sub(1);
    let nz = nz_planes.saturating_sub(1);
    for idx in 0..nx * ny * nz {
        let mat_number = idx + 1;
        if let Some((name, density)) = voxel_at(cell_materials, idx) {
            writeln!(cards, "C name: {name}").unwrap();
            // Upstream prints the density unless it holds the -1.0 sentinel.
            if *density != -1.0 {
                writeln!(cards, "C density = {density:.5}").unwrap();
            }
        }
        // Default-constructed material: no name metadata, no density line.
        writeln!(cards, "m{mat_number}").unwrap();
        // Nuclide fraction lines require composition data outside the
        // simplified input contract; FracType will shape them once
        // compositions are plumbed through.
    }
    cards
}

/// Format a float like Python `str(float)`: integral values gain `.0`, and
/// magnitudes below 1e-4 or at/above 1e16 switch to scientific notation with
/// a signed two-digit exponent.
fn py_float(v: f64) -> String {
    if !v.is_finite() {
        return v.to_string();
    }
    let mag = v.abs();
    if mag != 0.0 && !(1e-4..1e16).contains(&mag) {
        // Python repr style scientific: shortest mantissa, e±NN exponent.
        let s = format!("{v:e}");
        let (mant, exp) = s.split_once('e').expect("lowercase e always present");
        let (sign, digits) = match exp.strip_prefix('-') {
            Some(d) => ('-', d),
            None => ('+', exp),
        };
        if digits.len() < 2 {
            format!("{mant}e{sign}0{digits}")
        } else {
            format!("{mant}e{sign}{digits}")
        }
    } else if mag == 0.0 || v == v.trunc() {
        format!("{v:.1}")
    } else {
        // Rust's Display shares Python repr's shortest-roundtrip algorithm.
        format!("{v}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The exact reference 3x2x1 grid and materials.
    fn oracle_grid() -> String {
        let mats: Vec<Option<(String, f64)>> = [
            ("0", 42.0),
            ("1", 43.0),
            ("2", 44.0),
            ("3", 45.0),
            ("4", 47.0),
            ("5", 5.0),
        ]
        .iter()
        .map(|(n, d)| Some((n.to_string(), *d)))
        .collect();
        mesh_to_geom(
            &[0.0, 1.0, 2.0, 3.0],
            &[0.0, 1.0, 2.0],
            &[0.0, 1.0],
            &mats,
            &DeckOptions::default(),
        )
    }

    #[test]
    fn oracle_cell_and_surface_sections_match_exactly() {
        let g = oracle_grid();
        let lines: Vec<&str> = g.lines().collect();
        assert_eq!(lines[0], "Generated from mesh");
        let expected_cells = [
            "1 1 42.0 1 -2 5 -6 8 -9",
            "2 2 43.0 1 -2 6 -7 8 -9",
            "3 3 44.0 2 -3 5 -6 8 -9",
            "4 4 45.0 2 -3 6 -7 8 -9",
            "5 5 47.0 3 -4 5 -6 8 -9",
            "6 6 5.0 3 -4 6 -7 8 -9",
            "7 0 -1:4:-5:7:-8:9",
        ];
        assert_eq!(&lines[1..8], &expected_cells[..]);
        assert_eq!(lines[8], "");
        let expected_surfs = [
            "1 px 0.0", "2 px 1.0", "3 px 2.0", "4 px 3.0", "5 py 0.0", "6 py 1.0", "7 py 2.0",
            "8 pz 0.0", "9 pz 1.0",
        ];
        assert_eq!(&lines[9..18], &expected_surfs[..]);
        assert_eq!(lines[18], "");
    }

    #[test]
    fn oracle_material_card_headers() {
        let g = oracle_grid();
        // Upstream appends nuclide lines we cannot emit without composition
        // data; the surrounding structure must still match exactly.
        let expected_blocks: &[&[&str]] = &[
            &["C name: 0", "C density = 42.00000", "m1"],
            &["C name: 1", "C density = 43.00000", "m2"],
            &["C name: 2", "C density = 44.00000", "m3"],
            &["C name: 3", "C density = 45.00000", "m4"],
            &["C name: 4", "C density = 47.00000", "m5"],
            &["C name: 5", "C density = 5.00000", "m6"],
        ];
        let mut lines = g.lines().skip(19);
        for block in expected_blocks {
            for expected in *block {
                assert_eq!(lines.next(), Some(*expected));
            }
        }
        assert_eq!(lines.next(), None);
    }

    #[test]
    fn section_card_counts_match_mesh_shape() {
        let g = oracle_grid();
        let sections: Vec<Vec<&str>> = g.split("\n\n").map(|s| s.lines().collect()).collect();
        assert_eq!(sections.len(), 3); // title+cells | surfaces | materials
        assert_eq!(sections[0].len(), 8); // title + 6 voxel cells + graveyard
        assert_eq!(sections[1].len(), 9); // one surface per boundary plane
        assert_eq!(sections[2].len(), 18); // 3 mat-card lines x 6 voxels
    }

    #[test]
    fn first_and_last_cards() {
        let g = oracle_grid();
        assert!(g.starts_with("Generated from mesh\n1 1 42.0 1 -2 5 -6 8 -9\n"));
        // Last card is m6, the final voxel's material number.
        assert!(g.ends_with("C name: 5\nC density = 5.00000\nm6\n"));
    }

    #[test]
    fn iteration_order_is_z_fastest() {
        // Distinct density per voxel; idx = (i*ny + j)*nz + k.
        let mats: Vec<Option<(String, f64)>> =
            (1..=8).map(|d| Some((format!("m{d}"), d as f64))).collect();
        let g = mesh_to_geom(
            &[0.0, 1.0, 2.0],
            &[0.0, 1.0, 2.0],
            &[0.0, 1.0, 2.0],
            &mats,
            &DeckOptions::default(),
        );
        let lines: Vec<&str> = g.lines().collect();
        // Voxel idx 1 differs from idx 0 in z alone -> cell 2 carries 2.0.
        assert!(lines[2].starts_with("2 2 2.0 "));
        // Voxel idx 3 is (x0, y1, z1) -> cell 4 carries 4.0.
        assert!(lines[4].starts_with("4 4 4.0 "));
        // Graveyard after 8 cells: cell 9 over planes 1..=9.
        assert_eq!(lines[9], "9 0 -1:3:-4:6:-7:9");
    }

    #[test]
    fn unassigned_voxels_follow_default_material_semantics() {
        // Missing voxels fall back to the bare default material: density -1.0 on
        // the cell card, bare m-line without comments.
        let mats = vec![
            None,
            Some(("water".to_string(), 0.997)),
            None,
            Some(("u235".to_string(), -1.0)), // sentinel density
        ];
        let g = mesh_to_geom(
            &[0.0, 1.0, 2.0],
            &[0.0, 1.0],
            &[0.0, 1.0, 2.0],
            &mats,
            &DeckOptions::default(),
        );
        let lines: Vec<&str> = g.lines().collect();
        // 2x1x2 voxels; x planes 1..=3, y 4..=5, z 6..=8.
        assert_eq!(lines[1], "1 1 -1.0 1 -2 4 -5 6 -7");
        assert_eq!(lines[2], "2 2 0.997 1 -2 4 -5 7 -8");
        assert_eq!(lines[3], "3 3 -1.0 2 -3 4 -5 6 -7");
        assert_eq!(lines[4], "4 4 -1.0 2 -3 4 -5 7 -8");
        assert_eq!(lines[5], "5 0 -1:3:-4:5:-6:8");
        // Mat blocks: bare m1, commented m2, bare m3, name-only m4 (the
        // -1.0 sentinel suppresses the density comment).
        assert_eq!(
            &lines[16..23],
            &[
                "m1",
                "C name: water",
                "C density = 0.99700",
                "m2",
                "m3",
                "C name: u235",
                "m4"
            ][..]
        );
    }

    #[test]
    fn short_material_table_treated_as_unassigned() {
        // Fewer entries than voxels must not panic; extras behave like None.
        let g = mesh_to_geom(
            &[0.0, 1.0, 2.0],
            &[0.0, 1.0],
            &[0.0, 1.0],
            &[Some(("a".into(), 1.25))],
            &DeckOptions::default(),
        );
        assert!(g.contains("1 1 1.25 1 -2 4 -5 6 -7"));
        assert!(g.contains("\n2 2 -1.0 "));
        assert!(g.ends_with("C name: a\nC density = 1.25000\nm1\nm2\n"));
    }

    #[test]
    fn density_formatting_matches_python_str() {
        assert_eq!(py_float(42.0), "42.0");
        assert_eq!(py_float(0.997), "0.997");
        assert_eq!(py_float(-1.0), "-1.0");
        assert_eq!(py_float(0.0), "0.0");
        assert_eq!(py_float(0.0001), "0.0001");
        assert_eq!(py_float(1e-05), "1e-05");
        assert_eq!(py_float(1e16), "1e+16");
        assert_eq!(py_float(9999999999999998.0), "9999999999999998.0");
        assert_eq!(py_float(5.0e-5), "5e-05");
        assert_eq!(py_float(1.5e20), "1.5e+20");
    }

    #[test]
    fn custom_options_are_honored() {
        let opts = DeckOptions {
            title_card: "My custom title".into(),
            frac_type: FracType::Atom,
        };
        let g = mesh_to_geom(&[0.0, 1.0], &[0.0, 1.0], &[0.0, 1.0], &[], &opts);
        assert!(g.starts_with("My custom title\n"));
        // frac_type is accepted but cannot alter output until compositions
        // exist; the deck itself stays structurally valid.
        assert!(g.contains("1 1 -1.0"));
    }
}
