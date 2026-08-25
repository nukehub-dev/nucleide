//! Parser for Serpent `*_dep.m` depletion files.
//!
//! Depletion files carry plain vectors (`BU`, `DAYS`, `ZAI`, `NAMES`),
//! per-material scalars/vectors (`*_VOLUME`, `*_BURNUP`) and large inventory
//! matrices (`*_MDENS`, `*_ADENS`, ...). Arrays containing `%` comments are
//! stored as [`Entry::Matrix`] with one row per commented line; uncommented
//! arrays stay flat [`Entry::Vector`]s.
//!
//! Two Serpent quirks are handled explicitly:
//!
//! * redundant index variables of the form `iU234 = 1;` (name starting with
//!   `i` followed by a letter) are dropped, while numeric
//!   indices such as `i952421 = 123;` are kept;
//! * the running totals Serpent writes as accumulations
//!   (`TOT_MASS = TOT_MASS + MAT_x_VOLUME.*MAT_x_MDENS;`,
//!   `TOT_ADENS = TOT_ADENS./TOT_VOLUME;`, `zeros(r, c)`) are evaluated, so
//!   `TOT_*` entries hold final values.

use crate::matlab::{self, Rhs};
use crate::{Error, Result, Table};

/// Parse Serpent dep-file contents into a [`Table`].
pub fn parse_dep(text: &str) -> Result<Table> {
    let scan = matlab::scan(text)?;
    let mut table = Table::new();
    for st in scan.statements {
        if st.indexed {
            return Err(Error::Syntax {
                line: st.line,
                message: format!("unexpected indexed assignment `{}` in dep file", st.name),
            });
        }
        if is_imaterial_index(&st.name, &st.rhs) {
            continue;
        }
        matlab::apply_simple(&mut table, &st)?;
    }
    Ok(table)
}

/// True for removable index variables: names beginning with `i`
/// followed by a letter (e.g. `iU234`, `iLOST`), assigned an integer literal.
fn is_imaterial_index(name: &str, rhs: &Rhs) -> bool {
    let mut chars = name.chars();
    if chars.next() != Some('i') {
        return false;
    }
    match chars.next() {
        Some(c) if c.is_ascii_alphabetic() => {}
        _ => return false,
    }
    if !name
        .chars()
        .skip(2)
        .all(|c| c.is_ascii_alphanumeric() || c == '_')
    {
        return false;
    }
    matches!(rhs, Rhs::Expr(expr) if expr.trim().parse::<u64>().is_ok())
}

#[cfg(test)]
mod tests {
    use super::parse_dep;
    use crate::testing;

    #[test]
    fn dep1_vectors_and_indices() {
        let dep = testing::load_dep("sample1_dep.m");
        testing::assert_close(&dep.get_vec_f64("BU").unwrap(), &[0.0, 84.0, 168.0]);
        testing::assert_close(&dep.get_vec_f64("DAYS").unwrap(), &[0.0, 2100.0, 4200.0]);
        let zai = dep.get_vec_f64("ZAI").unwrap();
        assert_eq!(zai.len(), 146);
        assert_eq!(zai[0], 10010.0);
        assert_eq!(zai[zai.len() - 2], 666.0);
        assert_eq!(zai[zai.len() - 1], 0.0);
        assert_eq!(zai[122], 952_421.0);
        assert_eq!(dep.get_f64("i952421").unwrap(), 123.0);
    }

    #[test]
    fn dep1_names_and_imaterial_filter() {
        let dep = testing::load_dep("sample1_dep.m");
        let names = dep.get_vec_str("NAMES").unwrap();
        assert_eq!(names.len(), 146);
        assert_eq!(names[0], "H-1     ");
        assert_eq!(names[93].trim(), "U-235");
        assert_eq!(names[110].trim(), "Pu-239");
        assert_eq!(names[names.len() - 1], "total   ");
        assert!(dep.get("iLOST").is_none());
        assert!(dep.get("iTOT").is_none());
        assert!(dep.get("i10010").is_some());
    }

    #[test]
    fn dep1_material_shapes_and_values() {
        let dep = testing::load_dep("sample1_dep.m");
        let h = dep.get_matrix("MAT_fuelp1r2_H").unwrap();
        assert_eq!((h.rows(), h.cols()), (146, 3));
        testing::assert_close(&h.row_f64(3).unwrap(), &[0.0, 5.56191e-11, 3.22483e-10]);
        let bu = dep.get_matrix("MAT_fuelp1r2_BURNUP").unwrap();
        assert_eq!((bu.rows(), bu.cols()), (1, 3));
        testing::assert_close(&bu.row_f64(0).unwrap(), &[0.0, 7.44008e+01, 1.47936e+02]);
        testing::assert_close(
            &[dep.get_f64("MAT_fuelp1r1_VOLUME").unwrap()],
            &[1.394189e+01],
        );
    }

    #[test]
    fn dep1_totals_are_accumulated() {
        let dep = testing::load_dep("sample1_dep.m");
        testing::assert_close(&[dep.get_f64("TOT_VOLUME").unwrap()], &[139.4189]);
        let mass = dep.get_matrix("TOT_MASS").unwrap();
        assert_eq!((mass.rows(), mass.cols()), (146, 3));
        testing::assert_close(
            &mass.row_f64(93).unwrap(),
            &[37.4748243877, 1.546190455725, 0.23779287584],
        );
        testing::assert_close(
            &mass.row_f64(110).unwrap(),
            &[0.0, 3.407477384773, 3.385095074567],
        );
        let adens = dep.get_matrix("TOT_ADENS").unwrap();
        assert_eq!((adens.rows(), adens.cols()), (146, 3));
        testing::assert_close(
            &adens.row_f64(4).unwrap(),
            &[0.045381700000000004, 0.04537348, 0.045360979999999995],
        );
    }

    #[test]
    fn dep2_basics() {
        let dep = testing::load_dep("sample2_dep.m");
        testing::assert_close(
            &dep.get_vec_f64("BU").unwrap(),
            &[0.0, 1.00000e-01, 1.00000e+00],
        );
        testing::assert_close(
            &dep.get_vec_f64("DAYS").unwrap(),
            &[0.0, 2.50000e+00, 2.50000e+01],
        );
        let zai = dep.get_vec_f64("ZAI").unwrap();
        assert_eq!(zai.len(), 29);
        assert_eq!(zai[27], 666.0);
        assert_eq!(zai[28], 0.0);
        assert_eq!(dep.get_f64("i621510").unwrap(), 22.0);
        assert!(dep.get("iSm151").is_none());
        assert!(dep.get("iTOT").is_none());
        assert_eq!(dep.get_vec_str("NAMES").unwrap()[1].trim(), "U235");
    }

    #[test]
    fn dep2_shapes_and_values() {
        let dep = testing::load_dep("sample2_dep.m");
        let vol = dep.get_vec_f64("MAT_fuel_VOLUME").unwrap();
        assert_eq!(vol.len(), 3);
        testing::assert_close(&vol, &[5.33267e-01; 3]);
        let burnup = dep.get_vec_f64("MAT_fuel_BURNUP").unwrap();
        testing::assert_close(&burnup, &[0.0, 1.00874e-01, 1.00857e+00]);
        let a = dep.get_matrix("MAT_fuel_A").unwrap();
        assert_eq!((a.rows(), a.cols()), (29, 3));
        testing::assert_close(&a.row_f64(6).unwrap(), &[0.0, 1.73176e+05, 4.96485e+06]);
        let mass = dep.get_matrix("TOT_MASS").unwrap();
        assert_eq!((mass.rows(), mass.cols()), (29, 3));
        testing::assert_close(
            &mass.row_f64(0).unwrap(),
            &[1.27486e-03, 1.27302e-03, 1.25753e-03],
        );
    }

    #[test]
    fn garbage_dep_is_rejected() {
        assert!(parse_dep("not matlab at all\n").is_err());
        assert!(parse_dep("A = zeros(2,\n").is_err());
        assert!(parse_dep("").unwrap().is_empty());
    }

    #[test]
    fn unsupported_expression_reports_error() {
        assert!(parse_dep("A + B\n").is_err());
        let err = parse_dep("X = Y + Z;\n").unwrap_err();
        assert_eq!(err.to_string(), "variable `Y` not found");
    }
}
