//! Internal scanner for Serpent's MATLAB-flavoured assignment syntax.
//!
//! Serpent output files are sequences of statements of the forms
//!
//! ```text
//! NAME = <rhs>;
//! NAME(idx, 1)     = <scalar>;
//! NAME(idx, [1: n]) = [ v0 v1 ... vn-1 ];
//! ```
//!
//! interspersed with `%` comments and (in `*_res.m` files) repeated
//! `if (exist('idx', 'var')); idx = idx + 1; else; idx = 1; end;` counter
//! blocks. The scanner turns that text into [`Statement`]s; the per-file-kind
//! modules decide how to accumulate them into a [`Table`].

use crate::{Entry, Error, Matrix, Result, Table, Value};

/// Right-hand side of an assignment, as written in the source file.
pub(crate) enum Rhs {
    /// Single-quoted MATLAB character array (outer quotes retained).
    Quoted(String),
    /// Bracketed array without interior comments, flattened to tokens.
    Flat(Vec<String>),
    /// Bracketed array containing `%` comments; one token row per line.
    Rows(Vec<Vec<String>>),
    /// Bare expression: number literal, `zeros(r, c)` or an accumulation
    /// such as `TOT_MASS + MAT_x_VOLUME.*MAT_x_MDENS`.
    Expr(String),
}

/// A single parsed assignment.
pub(crate) struct Statement {
    pub name: String,
    /// True when the left-hand side carried an `(idx, ...)` suffix.
    pub indexed: bool,
    /// One-based source line, for error reporting.
    pub line: usize,
    pub rhs: Rhs,
}

/// Result of scanning a whole file.
pub(crate) struct Scan {
    pub statements: Vec<Statement>,
    /// Number of `if (exist(...))` counter blocks encountered.
    pub blocks: usize,
}

pub(crate) fn scan(text: &str) -> Result<Scan> {
    let mut cur = Cursor::new(text);
    let mut statements = Vec::new();
    let mut blocks = 0usize;
    loop {
        cur.skip_ws();
        if cur.peek().is_none() {
            break;
        }
        if cur.peek() == Some('%') {
            cur.skip_to_newline();
            continue;
        }
        if cur.starts_with("if (exist(") {
            cur.skip_counter_block()?;
            blocks += 1;
            continue;
        }
        statements.push(cur.read_statement()?);
    }
    Ok(Scan { statements, blocks })
}

struct Cursor {
    chars: Vec<char>,
    pos: usize,
    line: usize,
}

impl Cursor {
    fn new(text: &str) -> Self {
        Self {
            chars: text.chars().collect(),
            pos: 0,
            line: 1,
        }
    }

    fn peek(&self) -> Option<char> {
        self.chars.get(self.pos).copied()
    }

    fn bump(&mut self) -> Option<char> {
        let c = self.peek()?;
        self.pos += 1;
        if c == '\n' {
            self.line += 1;
        }
        Some(c)
    }

    fn starts_with(&self, prefix: &str) -> bool {
        self.chars[self.pos..]
            .iter()
            .zip(prefix.chars())
            .all(|(a, b)| a == &b)
    }

    fn skip_ws(&mut self) {
        while matches!(self.peek(), Some(c) if c.is_ascii_whitespace()) {
            self.bump();
        }
    }

    fn skip_to_newline(&mut self) {
        while matches!(self.peek(), Some(c) if c != '\n') {
            self.pos += 1;
        }
    }

    fn expect(&mut self, want: char) -> Result<()> {
        match self.peek() {
            Some(c) if c == want => {
                self.bump();
                Ok(())
            }
            other => Err(Error::Syntax {
                line: self.line,
                message: format!(
                    "expected `{want}`, found `{}`",
                    other
                        .map(|c| c.to_string())
                        .unwrap_or_else(|| "end of file".into())
                ),
            }),
        }
    }

    /// Consume the `if (exist('idx','var')) ... end;` block verbatim.
    fn skip_counter_block(&mut self) -> Result<()> {
        while self.pos < self.chars.len() {
            if self.starts_with("end;") {
                self.pos += 4;
                return Ok(());
            }
            self.bump();
        }
        Err(Error::Syntax {
            line: self.line,
            message: "unterminated `if (exist(` block".into(),
        })
    }

    fn read_ident(&mut self) -> Result<String> {
        let mut name = String::new();
        match self.peek() {
            Some(c) if c.is_ascii_alphabetic() || c == '_' => name.push(self.bump().unwrap()),
            other => {
                return Err(Error::Syntax {
                    line: self.line,
                    message: format!(
                        "expected variable name, found `{}`",
                        other
                            .map(|c| c.to_string())
                            .unwrap_or_else(|| "end of file".into()),
                    ),
                })
            }
        }
        while matches!(self.peek(), Some(c) if c.is_ascii_alphanumeric() || c == '_') {
            name.push(self.bump().unwrap());
        }
        Ok(name)
    }

    /// Consume a balanced `( idx, [1: n] )` index spec.
    fn skip_index_group(&mut self) -> Result<()> {
        self.bump();
        let mut bracket_depth = 0usize;
        while let Some(c) = self.bump() {
            match c {
                '[' => bracket_depth += 1,
                ']' => bracket_depth = bracket_depth.saturating_sub(1),
                ')' if bracket_depth == 0 => return Ok(()),
                _ => {}
            }
        }
        Err(Error::Syntax {
            line: self.line,
            message: "unterminated `(` in assignment index".into(),
        })
    }

    fn read_statement(&mut self) -> Result<Statement> {
        let line = self.line;
        let name = self.read_ident()?;
        self.skip_ws();
        let mut indexed = false;
        if self.peek() == Some('(') {
            self.skip_index_group()?;
            indexed = true;
            self.skip_ws();
        }
        self.expect('=')?;
        self.bump();
        self.skip_ws();
        let rhs = match self.peek() {
            Some('\'') => {
                self.bump();
                Rhs::Quoted(self.take_quoted()?)
            }
            Some('[') => {
                self.bump();
                self.read_array()?
            }
            _ => Rhs::Expr(self.take_expr()),
        };
        // Optional statement terminator, possibly separated by spaces.
        loop {
            match self.peek() {
                Some(c) if c == ' ' || c == '\t' || c == '\r' => {
                    self.bump();
                }
                Some(';') => {
                    self.bump();
                    break;
                }
                _ => break,
            }
        }
        Ok(Statement {
            name,
            indexed,
            line,
            rhs,
        })
    }

    fn take_quoted(&mut self) -> Result<String> {
        let start = self.pos;
        while let Some(c) = self.bump() {
            if c == '\'' {
                let s: String = self.chars[start..self.pos - 1].iter().collect();
                if s.is_empty() {
                    return Err(Error::Syntax {
                        line: self.line,
                        message: "empty string literal".into(),
                    });
                }
                return Ok(s);
            }
        }
        Err(Error::Syntax {
            line: self.line,
            message: "unterminated string literal".into(),
        })
    }

    fn take_expr(&mut self) -> String {
        let start = self.pos;
        while let Some(c) = self.peek() {
            if c == ';' || c == '\n' {
                break;
            }
            self.pos += 1;
        }
        self.chars[start..self.pos]
            .iter()
            .collect::<String>()
            .trim()
            .to_string()
    }

    fn read_array(&mut self) -> Result<Rhs> {
        let start = self.pos;
        let mut depth = 1usize;
        while let Some(c) = self.bump() {
            match c {
                '[' => depth += 1,
                ']' => {
                    depth -= 1;
                    if depth == 0 {
                        let content: String = self.chars[start..self.pos - 1].iter().collect();
                        return finalize_array(&content, self.line);
                    }
                }
                _ => {}
            }
        }
        Err(Error::Syntax {
            line: self.line,
            message: "unterminated `[` array".into(),
        })
    }
}

fn finalize_array(content: &str, line: usize) -> Result<Rhs> {
    if content.contains('%') {
        let rows: Vec<Vec<String>> = content
            .lines()
            .map(|line| tokenize(line.split('%').next().unwrap_or("")))
            .filter(|toks| !toks.is_empty())
            .collect();
        if rows.is_empty() {
            return Err(Error::Syntax {
                line,
                message: "array contains only comments".into(),
            });
        }
        Ok(Rhs::Rows(rows))
    } else {
        let toks = tokenize(content);
        if toks.is_empty() {
            return Err(Error::Syntax {
                line,
                message: "empty array".into(),
            });
        }
        Ok(Rhs::Flat(toks))
    }
}

/// Split bracket contents into tokens; quoted spans may contain separators.
fn tokenize(s: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut in_quote = false;
    for c in s.chars() {
        if in_quote {
            cur.push(c);
            if c == '\'' {
                out.push(std::mem::take(&mut cur));
                in_quote = false;
            }
        } else {
            match c {
                '\'' => {
                    if !cur.is_empty() {
                        out.push(std::mem::take(&mut cur));
                    }
                    in_quote = true;
                    cur.push('\'');
                }
                c if c.is_ascii_whitespace() || c == ',' => {
                    if !cur.is_empty() {
                        out.push(std::mem::take(&mut cur));
                    }
                }
                _ => cur.push(c),
            }
        }
    }
    if !cur.is_empty() {
        out.push(cur);
    }
    out
}

pub(crate) fn parse_number(tok: &str, line: usize) -> Result<f64> {
    tok.trim().parse::<f64>().map_err(|_| Error::Syntax {
        line,
        message: format!("invalid number `{}`", tok.trim()),
    })
}

pub(crate) fn tokens_to_values(tokens: &[String], line: usize) -> Result<Vec<Value>> {
    tokens
        .iter()
        .map(|t| {
            if t.starts_with('\'') && t.ends_with('\'') && t.len() >= 2 {
                Ok(Value::Str(t[1..t.len() - 1].to_string()))
            } else {
                parse_number(t, line).map(Value::Num)
            }
        })
        .collect()
}

/// Apply a non-indexed statement to `table`.
pub(crate) fn apply_simple(table: &mut Table, st: &Statement) -> Result<()> {
    let entry = match &st.rhs {
        Rhs::Quoted(s) => Entry::Scalar(Value::Str(s.clone())),
        Rhs::Flat(tokens) => Entry::Vector(tokens_to_values(tokens, st.line)?),
        Rhs::Rows(rows) => Entry::Matrix(matrix_from_rows(rows, st.line)?),
        Rhs::Expr(expr) => eval_expr(expr, table, st.line)?,
    };
    table.insert(st.name.clone(), entry);
    Ok(())
}

pub(crate) fn matrix_from_rows(rows: &[Vec<String>], line: usize) -> Result<Matrix> {
    let values: Vec<Vec<Value>> = rows
        .iter()
        .map(|row| tokens_to_values(row, line))
        .collect::<Result<Vec<_>>>()?;
    Matrix::from_rows(values)
}

// --- minimal expression evaluation (dep-file accumulations) -----------------

enum Numeric {
    Scalar(f64),
    Vector(Vec<f64>),
    Tensor(usize, usize, Vec<f64>),
}

/// Evaluate the subset of MATLAB expressions Serpent writes in `*_dep.m`
/// files: numeric literals, `zeros(r, c)`, variable references and the
/// elementwise operators `+`, `.*`, `./`. Operators are expected to be
/// space-separated (`a + b.*c`) exactly as Serpent emits them.
pub(crate) fn eval_expr(expr: &str, table: &Table, line: usize) -> Result<Entry> {
    let expr = expr.trim();
    if expr.is_empty() {
        return Err(Error::Syntax {
            line,
            message: "empty expression".into(),
        });
    }
    let mut total: Option<Numeric> = None;
    for term in expr.split(" + ") {
        let value = eval_term(term.trim(), table, line)?;
        total = Some(match total {
            None => value,
            Some(acc) => add(acc, value, line)?,
        });
    }
    let total = total.ok_or_else(|| Error::Syntax {
        line,
        message: "empty expression".into(),
    })?;
    Ok(numeric_to_entry(total))
}

fn eval_term(term: &str, table: &Table, line: usize) -> Result<Numeric> {
    if let Some(i) = term.find(".*") {
        let lhs = eval_atom(&term[..i], table, line)?;
        let rhs = eval_atom(&term[i + 2..], table, line)?;
        mul(lhs, rhs, line)
    } else if let Some(i) = term.find("./") {
        let lhs = eval_atom(&term[..i], table, line)?;
        let rhs = eval_atom(&term[i + 2..], table, line)?;
        div(lhs, rhs, line)
    } else {
        eval_atom(term, table, line)
    }
}

fn eval_atom(atom: &str, table: &Table, line: usize) -> Result<Numeric> {
    let atom = atom.trim();
    if atom.is_empty() {
        return Err(Error::Syntax {
            line,
            message: "missing operand".into(),
        });
    }
    let first = atom.chars().next().unwrap();
    if first.is_ascii_digit() || first == '-' || first == '+' || first == '.' {
        return parse_number(atom, line).map(Numeric::Scalar);
    }
    if let Some(rest) = atom.strip_prefix("zeros(") {
        let dims = rest.strip_suffix(')').ok_or_else(|| Error::Syntax {
            line,
            message: format!("malformed `{atom}`"),
        })?;
        let mut parts = dims.split(',');
        let r: usize = parts
            .next()
            .unwrap_or("")
            .trim()
            .parse()
            .map_err(|_| Error::Syntax {
                line,
                message: format!("bad row count in `{atom}`"),
            })?;
        let c: usize = parts
            .next()
            .unwrap_or("")
            .trim()
            .parse()
            .map_err(|_| Error::Syntax {
                line,
                message: format!("bad column count in `{atom}`"),
            })?;
        return Ok(Numeric::Tensor(r, c, vec![0.0; r.saturating_mul(c)]));
    }
    match table.get(atom) {
        Some(entry) => entry_to_numeric(entry).ok_or_else(|| Error::Type {
            context: atom.to_string(),
            expected: "numeric value",
            found: entry.type_name(),
        }),
        None if !atom.chars().all(|c| c.is_ascii_alphanumeric() || c == '_') => {
            Err(Error::UnsupportedExpr(atom.to_string()))
        }
        None => Err(Error::Missing(atom.to_string())),
    }
}

fn entry_to_numeric(entry: &Entry) -> Option<Numeric> {
    match entry {
        Entry::Scalar(Value::Num(v)) => Some(Numeric::Scalar(*v)),
        Entry::Vector(vals) if vals.iter().all(|v| matches!(v, Value::Num(_))) => Some(
            Numeric::Vector(vals.iter().filter_map(|v| v.as_f64().ok()).collect()),
        ),
        Entry::Matrix(m) if m.data.iter().all(|v| matches!(v, Value::Num(_))) => {
            Some(Numeric::Tensor(
                m.rows,
                m.cols,
                m.data.iter().filter_map(|v| v.as_f64().ok()).collect(),
            ))
        }
        _ => None,
    }
}

fn numeric_to_entry(n: Numeric) -> Entry {
    match n {
        Numeric::Scalar(v) => Entry::Scalar(Value::Num(v)),
        Numeric::Vector(data) => Entry::Vector(data.into_iter().map(Value::Num).collect()),
        Numeric::Tensor(rows, cols, data) => Entry::Matrix(Matrix::from_parts(
            rows,
            cols,
            data.into_iter().map(Value::Num).collect(),
        )),
    }
}

macro_rules! zip_elements {
    ($line:expr, $a:expr, $b:expr, $len_a:expr, $len_b:expr, $op:expr) => {{
        if $len_a != $len_b {
            return Err(Error::Syntax {
                line: $line,
                message: format!("shape mismatch in expression ({} vs {})", $len_a, $len_b),
            });
        }
        $a.iter_mut()
            .zip($b.iter())
            .for_each(|(x, y)| *x = $op(*x, *y));
    }};
}

fn add(a: Numeric, b: Numeric, line: usize) -> Result<Numeric> {
    use Numeric::*;
    Ok(match (a, b) {
        (Scalar(x), Scalar(y)) => Scalar(x + y),
        (Scalar(x), Tensor(r, c, mut d)) | (Tensor(r, c, mut d), Scalar(x)) => {
            d.iter_mut().for_each(|v| *v += x);
            Tensor(r, c, d)
        }
        (Vector(mut x), Vector(y)) => {
            zip_elements!(line, x, y, x.len(), y.len(), |p: f64, q: f64| p + q);
            Vector(x)
        }
        (Vector(mut x), Scalar(y)) | (Scalar(y), Vector(mut x)) => {
            x.iter_mut().for_each(|v| *v += y);
            Vector(x)
        }
        (Tensor(ra, ca, mut x), Tensor(rb, cb, y)) => {
            if ra != rb || ca != cb {
                return Err(Error::Syntax {
                    line,
                    message: format!("matrix shape mismatch ({ra}x{ca} vs {rb}x{cb})"),
                });
            }
            zip_elements!(line, x, y, x.len(), y.len(), |p: f64, q: f64| p + q);
            Tensor(ra, ca, x)
        }
        (Tensor(..), Vector(_)) | (Vector(_), Tensor(..)) => {
            return Err(Error::Syntax {
                line,
                message: "cannot mix matrix and vector operands".into(),
            })
        }
    })
}

fn mul(a: Numeric, b: Numeric, line: usize) -> Result<Numeric> {
    use Numeric::*;
    Ok(match (a, b) {
        (Scalar(x), Scalar(y)) => Scalar(x * y),
        (Scalar(x), Tensor(r, c, mut d)) | (Tensor(r, c, mut d), Scalar(x)) => {
            d.iter_mut().for_each(|v| *v *= x);
            Tensor(r, c, d)
        }
        (Scalar(x), Vector(mut d)) | (Vector(mut d), Scalar(x)) => {
            d.iter_mut().for_each(|v| *v *= x);
            Vector(d)
        }
        (Vector(mut x), Vector(y)) => {
            zip_elements!(line, x, y, x.len(), y.len(), |p: f64, q: f64| p * q);
            Vector(x)
        }
        (Tensor(ra, ca, mut x), Tensor(rb, cb, y)) => {
            if ra != rb || ca != cb {
                return Err(Error::Syntax {
                    line,
                    message: format!("matrix shape mismatch ({ra}x{ca} vs {rb}x{cb})"),
                });
            }
            zip_elements!(line, x, y, x.len(), y.len(), |p: f64, q: f64| p * q);
            Tensor(ra, ca, x)
        }
        (Tensor(..), Vector(_)) | (Vector(_), Tensor(..)) => {
            return Err(Error::Syntax {
                line,
                message: "cannot mix matrix and vector operands".into(),
            })
        }
    })
}

fn div(a: Numeric, b: Numeric, line: usize) -> Result<Numeric> {
    use Numeric::*;
    Ok(match (a, b) {
        (Scalar(x), Scalar(y)) => Scalar(x / y),
        (Tensor(r, c, mut d), Scalar(x)) => {
            d.iter_mut().for_each(|v| *v /= x);
            Tensor(r, c, d)
        }
        (Vector(mut d), Scalar(x)) => {
            d.iter_mut().for_each(|v| *v /= x);
            Vector(d)
        }
        (Vector(mut x), Vector(y)) => {
            zip_elements!(line, x, y, x.len(), y.len(), |p: f64, q: f64| p / q);
            Vector(x)
        }
        (Tensor(ra, ca, mut x), Tensor(rb, cb, y)) => {
            if ra != rb || ca != cb {
                return Err(Error::Syntax {
                    line,
                    message: format!("matrix shape mismatch ({ra}x{ca} vs {rb}x{cb})"),
                });
            }
            zip_elements!(line, x, y, x.len(), y.len(), |p: f64, q: f64| p / q);
            Tensor(ra, ca, x)
        }
        (Scalar(_), Tensor(..)) | (Scalar(_), Vector(_)) => {
            return Err(Error::Syntax {
                line,
                message: "unsupported scalar / container operand".into(),
            })
        }
        (Tensor(..), Vector(_)) | (Vector(_), Tensor(..)) => {
            return Err(Error::Syntax {
                line,
                message: "cannot mix matrix and vector operands".into(),
            })
        }
    })
}
