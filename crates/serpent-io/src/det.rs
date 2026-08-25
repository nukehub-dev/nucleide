//! Parser for Serpent `*_det.m` detector files.
//!
//! A det file contains scalar bin counts (`DET<name>_EBINS`, `_VALS`, ...)
//! and one flat value array per detector (`DET<name> = [...]`). After
//! scanning, each detector array is reshaped into an [`Entry::Matrix`]
//! following the reference layout rules:
//!
//! * Serpent 1 (detected via `<name>_VALS`/`<name>_EBINS` scalars): 13
//!   columns when a matching `<name>E` energy-bin array exists, with row
//!   count taken from the bin-count scalars; otherwise 3 columns;
//! * Serpent 2: 13 columns when both `E` and `T` bin arrays exist, 12 when
//!   only `E` exists, otherwise 3, with row count derived from the array
//!   length.

use crate::matlab::{self};
use crate::{Entry, Error, Matrix, Result, Table};

/// Parse Serpent det-file contents into a [`Table`].
pub fn parse_det(text: &str) -> Result<Table> {
    let scan = matlab::scan(text)?;
    let mut table = Table::new();
    for st in scan.statements {
        if st.indexed {
            return Err(Error::Syntax {
                line: st.line,
                message: format!("unexpected indexed assignment `{}` in det file", st.name),
            });
        }
        matlab::apply_simple(&mut table, &st)?;
    }
    reshape_detectors(&mut table)?;
    Ok(table)
}

fn reshape_detectors(table: &mut Table) -> Result<()> {
    let detector_names: Vec<String> = table
        .iter()
        .filter(|(key, entry)| key.starts_with("DET") && matches!(entry, Entry::Vector(_)))
        .map(|(key, _)| key.clone())
        .collect();
    if detector_names.is_empty() {
        return Ok(());
    }

    let serpent1 = table.keys().any(|key| {
        (key.ends_with("_VALS") && detector_names.contains(&key[..key.len() - 5].to_string()))
            || (key.ends_with("_EBINS")
                && detector_names.contains(&key[..key.len() - 6].to_string()))
    });

    for name in detector_names {
        let values = match table.get(&name) {
            Some(Entry::Vector(vals)) => vals.clone(),
            _ => continue,
        };
        let len = values.len();
        let (rows, cols) = detector_shape(table, &name, len, serpent1)?;
        if rows * cols != len {
            return Err(Error::Syntax {
                line: 0,
                message: format!(
                    "detector `{name}` holds {len} values which do not fit {rows}x{cols}"
                ),
            });
        }
        table.insert(name, Entry::Matrix(Matrix::from_parts(rows, cols, values)));
    }
    Ok(())
}

fn detector_shape(table: &Table, name: &str, len: usize, serpent1: bool) -> Result<(usize, usize)> {
    let has_bins = |suffix: &str| table.contains_key(format!("{name}{suffix}").as_str());
    if serpent1 {
        if has_bins("E") {
            let rows = table.get_f64(format!("{name}_VALS"))? as usize;
            Ok((rows, 13))
        } else {
            let base = &name[..name.len() - 1];
            let rows = table.get_f64(format!("{base}_EBINS"))? as usize;
            Ok((rows, 3))
        }
    } else if has_bins("T") {
        Ok((len / 13, 13))
    } else if has_bins("E") {
        Ok((len / 12, 12))
    } else {
        Ok((len / 3, 3))
    }
}

#[cfg(test)]
mod tests {
    use super::parse_det;
    use crate::{testing, Entry};

    #[test]
    fn det1_shapes() {
        let det = testing::load_det("sample_det.m");
        assert_eq!(det.get_f64("DETphi_EBINS").unwrap(), 63.0);
        assert_eq!(det.get_f64("DETphi_VALS").unwrap(), 63.0);
        assert_eq!(det.get_f64("DETphi_RBINS").unwrap(), 1.0);
        let phi = det.get_matrix("DETphi").unwrap();
        assert_eq!((phi.rows(), phi.cols()), (63, 13));
        let phi_e = det.get_matrix("DETphiE").unwrap();
        assert_eq!((phi_e.rows(), phi_e.cols()), (63, 3));
        for (key, entry) in &det {
            if key.contains('_') {
                assert!(
                    matches!(entry, Entry::Scalar(_)),
                    "`{key}` should stay scalar"
                );
            } else {
                assert!(
                    matches!(entry, Entry::Matrix(_)),
                    "`{key}` should be reshaped"
                );
            }
        }
    }

    #[test]
    fn det1_values() {
        let det = testing::load_det("sample_det.m");
        let phi = det.get_matrix("DETphi").unwrap();
        testing::assert_close(
            &phi.row_f64(6).unwrap(),
            &[
                7.0,
                7.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                2.92709e-02,
                0.00857,
                16768.0,
            ],
        );
        let phi_e = det.get_matrix("DETphiE").unwrap();
        testing::assert_close(
            &phi_e.row_f64(phi_e.rows() - 3).unwrap(),
            &[1.49182e+01, 1.69046e+01, 1.49182e+01],
        );
    }

    #[test]
    fn det2_shapes() {
        let det = testing::load_det("serp2_det.m");
        let expected = [
            ("DET1", (15, 13)),
            ("DET1E", (5, 3)),
            ("DET1T", (3, 3)),
            ("DET2", (240, 13)),
            ("DET2E", (5, 3)),
            ("DET2T", (3, 3)),
            ("DET2X", (4, 3)),
            ("DET2Y", (4, 3)),
            ("DET3", (3, 13)),
            ("DET3T", (3, 3)),
        ];
        for (key, (rows, cols)) in expected {
            let m = det
                .get_matrix(key)
                .expect("detector must exist after reshape");
            assert_eq!((m.rows(), m.cols()), (rows, cols), "shape of {key}");
        }
    }

    #[test]
    fn det2_values() {
        let det = testing::load_det("serp2_det.m");
        let d1 = det.get_matrix("DET1").unwrap();
        testing::assert_close(
            &d1.row_f64(4).unwrap(),
            &[
                5.0,
                1.0,
                5.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                5.11865e+05,
                0.00417,
            ],
        );
        let d1e = det.get_matrix("DET1E").unwrap();
        testing::assert_close(
            &d1e.row_f64(d1e.rows() - 3).unwrap(),
            &[5.25306e-05, 3.80731e-03, 1.92992e-03],
        );
    }

    #[test]
    fn garbage_det_is_rejected() {
        assert!(parse_det("nope nope\n").is_err());
        assert!(parse_det("DET1 = [1 2\n").is_err());
        assert!(parse_det("").unwrap().is_empty());
    }
}
