//! Parser for Serpent `*_res.m` results files.
//!
//! A res file contains one "block" per transport/burnup step, each introduced
//! by an `if (exist('idx', 'var')) ... end;` counter. Every variable is
//! written once per block with an `(idx, ...)` left-hand side, so parsed
//! entries accumulate across blocks:
//!
//! * `NAME(idx, 1) = v;` (scalar per block) becomes a [`Entry::Vector`] of
//!   length equal to the number of blocks;
//! * `NAME(idx, [1: n]) = [...];` becomes an `n`-column [`Entry::Matrix`]
//!   with one row per block;
//! * quoted strings become vectors of [`Value::Str`] regardless of their
//!   declared width.
//!
//! The number of blocks is stored under the key `"IDX"`.

use crate::matlab::{self, Rhs};
use crate::{Entry, Error, Matrix, Result, Table, Value};

/// Parse Serpent res-file contents into a [`Table`].
pub fn parse_res(text: &str) -> Result<Table> {
    let scan = matlab::scan(text)?;
    let mut table = Table::new();
    for st in scan.statements {
        apply_res(&mut table, st)?;
    }
    table.insert(
        "IDX".to_string(),
        Entry::Scalar(Value::Num(scan.blocks as f64)),
    );
    Ok(table)
}

fn apply_res(table: &mut Table, st: matlab::Statement) -> Result<()> {
    if !st.indexed {
        return matlab::apply_simple(table, &st);
    }
    match st.rhs {
        Rhs::Quoted(s) => push_scalar(table, st.name, Value::Str(s)),
        Rhs::Flat(tokens) => append_row(
            table,
            st.name,
            matlab::tokens_to_values(&tokens, st.line)?,
            st.line,
        ),
        Rhs::Rows(rows) => {
            let flat: Vec<String> = rows.into_iter().flatten().collect();
            append_row(
                table,
                st.name,
                matlab::tokens_to_values(&flat, st.line)?,
                st.line,
            )
        }
        Rhs::Expr(expr) => {
            let value = matlab::eval_expr(&expr, table, st.line)?;
            let num = value.as_f64().map_err(|_| Error::Syntax {
                line: st.line,
                message: format!("non-scalar value assigned to `{}`", st.name),
            })?;
            push_scalar(table, st.name, Value::Num(num))
        }
    }
}

fn push_scalar(table: &mut Table, name: String, value: Value) -> Result<()> {
    if !table.contains_key(&name) {
        table.insert(name.clone(), Entry::Vector(Vec::new()));
    }
    match table.get_mut(&name) {
        Some(Entry::Vector(vals)) => {
            vals.push(value);
            Ok(())
        }
        Some(other) => Err(Error::Type {
            context: name,
            expected: "vector",
            found: other.type_name(),
        }),
        None => unreachable!("just inserted"),
    }
}

fn append_row(table: &mut Table, name: String, values: Vec<Value>, line: usize) -> Result<()> {
    if let Some(existing) = table.get_mut(&name) {
        let matrix = match existing {
            Entry::Matrix(m) => m,
            other => {
                return Err(Error::Type {
                    context: name.clone(),
                    expected: "matrix",
                    found: other.type_name(),
                })
            }
        };
        if values.len() != matrix.cols() {
            return Err(Error::Syntax {
                line,
                message: format!(
                    "row of {} values does not fit `{name}` with {} columns",
                    values.len(),
                    matrix.cols()
                ),
            });
        }
        matrix.push_row(values);
        return Ok(());
    }
    let mut m = Matrix::from_parts(0, values.len(), Vec::new());
    m.push_row(values);
    table.insert(name, Entry::Matrix(m));
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::parse_res;
    use crate::{testing, Entry};

    #[test]
    fn res1_meta_and_block_count() {
        let res = testing::load_res("sample_res.m");
        assert_eq!(res.get_f64("IDX").unwrap(), 3.0);
        let version = res.get_vec_str("VERSION").unwrap();
        assert_eq!(version.len(), 3);
        assert_eq!(version[0], "Serpent 1.1.7");
        assert_eq!(
            res.get_vec_str("TITLE").unwrap()[1],
            "[CHAR] lwr2 Burnup Calculation"
        );
        assert_eq!(
            res.get_vec_str("DATE").unwrap()[2],
            "Sun May 22 22:33:50 2011"
        );
        assert_eq!(res.get_vec_f64("POP").unwrap(), vec![5000.0; 3]);
        testing::assert_close(
            &res.get_vec_f64("BURN_DAYS").unwrap(),
            &[0.0, 2100.0, 4200.0],
        );
        testing::assert_close(&res.get_vec_f64("BURNUP").unwrap(), &[0.0, 84.0, 168.0]);
    }

    #[test]
    fn res1_keff_values() {
        let res = testing::load_res("sample_res.m");
        let six_ff = res.get_matrix("SIX_FF_ETA").unwrap();
        assert_eq!(six_ff.rows(), 3);
        assert_eq!(six_ff.cols(), 2);
        testing::assert_close(&six_ff.row_f64(1).unwrap(), &[1.16446e0, 0.00186]);
        testing::assert_close(
            &res.get_matrix("IMP_KEFF").unwrap().row_f64(0).unwrap(),
            &[1.24207e0, 0.00053],
        );
        testing::assert_close(
            &res.get_matrix("COL_KEFF").unwrap().row_f64(2).unwrap(),
            &[6.59713e-01, 0.00156],
        );
    }

    #[test]
    fn res1_peakf10_last_block() {
        let res = testing::load_res("sample_res.m");
        let peakf = res.get_matrix("PEAKF10").unwrap();
        assert_eq!((peakf.rows(), peakf.cols()), (3, 4));
        testing::assert_close(
            &peakf.row_f64(peakf.rows() - 1).unwrap(),
            &[12.0, 11.0, 1.09824e0, 0.01768],
        );
        testing::assert_close(
            &peakf.row_f64(0).unwrap(),
            &[13.0, 12.0, 1.08448e0, 0.01643],
        );
        testing::assert_close(&peakf.row_f64(1).unwrap(), &[5.0, 14.0, 1.07986e0, 0.01784]);
    }

    #[test]
    fn res1_all_shapes_match_idx() {
        let res = testing::load_res("sample_res.m");
        for (key, entry) in &res {
            if key == "IDX" {
                continue;
            }
            match entry {
                Entry::Vector(v) => assert_eq!(v.len(), 3, "vector `{key}`"),
                Entry::Matrix(m) => assert_eq!(m.rows(), 3, "matrix `{key}`"),
                Entry::Scalar(_) => panic!("unexpected scalar `{key}`"),
            }
        }
        assert_eq!(res.get_matrix("GTRANSFP").unwrap().cols(), 200);
        assert_eq!(res.get_matrix("FLUX").unwrap().cols(), 22);
        assert_eq!(res.get_matrix("GC_BOUNDS").unwrap().cols(), 11);
        assert_eq!(
            res.get_vec_str("CPU_TYPE").unwrap()[0],
            "Intel(R) Core(TM) i7 CPU       M 640  @ 2.80GHz"
        );
    }

    #[test]
    fn res2_layout_and_values() {
        let res = testing::load_res("serp2_res.m");
        assert_eq!(res.get_f64("IDX").unwrap(), 1.0);
        assert_eq!(res.get_vec_str("VERSION").unwrap()[0], "Serpent 2.1.21");
        assert_eq!(res.get_vec_str("INPUT_FILE_NAME").unwrap()[0], "serpInput");
        testing::assert_close(
            &res.get_matrix("MEAN_POP_SIZE").unwrap().row_f64(0).unwrap(),
            &[1.00140e+02, 0.00359],
        );
        assert_eq!(res.get_vec_f64("POP").unwrap(), vec![100.0]);
        let ana = res.get_matrix("ANA_KEFF").unwrap();
        assert_eq!((ana.rows(), ana.cols()), (1, 6));
        testing::assert_close(&ana.row_f64(0).unwrap()[..2], &[9.19924e-01, 0.00477]);
        assert_eq!(res.get_matrix("MICRO_E").unwrap().cols(), 71);
        assert_eq!(res.get_vec_str("GC_UNIVERSE_NAME").unwrap(), vec!["0"]);
        for (key, entry) in &res {
            if key == "IDX" {
                continue;
            }
            match entry {
                Entry::Vector(v) => assert_eq!(v.len(), 1, "vector `{key}`"),
                Entry::Matrix(m) => assert_eq!(m.rows(), 1, "matrix `{key}`"),
                Entry::Scalar(_) => panic!("unexpected scalar `{key}`"),
            }
        }
    }

    #[test]
    fn garbage_res_is_rejected() {
        assert!(parse_res("this is not matlab\n").is_err());
        assert!(parse_res("X = [1 2 3;\n").is_err());
        assert_eq!(parse_res("").unwrap().get_f64("IDX").unwrap(), 0.0);
    }
}
