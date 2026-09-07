//! MCNP cell-card parsing.
//!
//! A cell card is `num mat [dens] geom... [params...]`: the density is
//! absent for void cells (`mat == 0`). Geometry is a CSG expression over
//! signed surface numbers with `:` union (lowest precedence),
//! juxtaposition intersection, parentheses, `#n` cell complement, and an
//! optional `*` reflecting marker on half-spaces. Remaining tokens are
//! carried as cell parameters (`imp:n=1`, `u=3`, `fill=...`, ...).
//!
//! # Line conventions (shared with [`crate::inp`])
//!
//! - A card starts on a line with fewer than five leading spaces; lines
//!   with five or more leading spaces continue the previous card. Comment
//!   lines (first non-blank character `c`/`C` followed by blank or EOL)
//!   are skipped by the standalone parsers and preserved by
//!   [`crate::problem`].
//! - `$` truncates the line (trailing comments contribute no tokens).

use crate::inp::Error;

/// A signed surface half-space, with MCNP's `*` reflecting marker.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HalfSpace {
    /// Signed surface number (sign = sense).
    pub surf: i32,
    /// `*n` reflecting marker.
    pub reflecting: bool,
}

/// CSG geometry expression over surface half-spaces.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GeomExpr {
    /// One signed half-space (`-12`, `34`, `*-5`).
    HalfSpace(HalfSpace),
    /// Juxtaposition intersection (`1 -2 3`).
    Intersect(Vec<GeomExpr>),
    /// Colon union (`1 : -2`).
    Union(Box<GeomExpr>, Box<GeomExpr>),
    /// Cell complement (`#5`).
    Complement(Box<GeomExpr>),
}

/// One parsed MCNP cell card.
#[derive(Debug, Clone, PartialEq)]
pub struct CellCard {
    /// Cell number.
    pub num: u32,
    /// Material number (`0` = void).
    pub mat: u32,
    /// Density (atom or mass density as written); absent for void cells.
    pub dens: Option<f64>,
    /// CSG geometry expression.
    pub geom: GeomExpr,
    /// Trailing parameter tokens in card order (`imp:n=1`, `vol=...`).
    pub params: Vec<String>,
    /// Source lines forming this card (for format-preserving write-back).
    pub raw_lines: Vec<String>,
    /// Comment/blank lines preceding this card in its block, verbatim, so
    /// write-back preserves them in order (see [`crate::problem`]).
    pub prefix_lines: Vec<String>,
}

impl GeomExpr {
    /// Render in canonical MCNP form (single spaces, minimal parentheses).
    pub fn render(&self) -> String {
        match self {
            GeomExpr::HalfSpace(h) => {
                format!("{}{}", if h.reflecting { "*" } else { "" }, h.surf)
            }
            GeomExpr::Intersect(parts) => parts
                .iter()
                .map(|p| match p {
                    GeomExpr::Union(..) => format!("({})", p.render()),
                    _ => p.render(),
                })
                .collect::<Vec<_>>()
                .join(" "),
            GeomExpr::Union(a, b) => format!("{} : {}", a.render(), b.render()),
            GeomExpr::Complement(inner) => format!("#{}", inner.render()),
        }
    }
}

/// Join physical lines into logical card lines: continuations (five or more
/// leading spaces) append to the current card. Returns
/// `(1-based start line, joined text)`; blank and comment lines are
/// skipped (see [`is_comment`]).
pub(crate) fn join_continuations(lines: &[&str]) -> Vec<(usize, String)> {
    let mut out = Vec::new();
    let mut current: Option<(usize, String)> = None;
    for (idx, line) in lines.iter().enumerate() {
        let stripped = line.split('$').next().unwrap_or("").trim_end();
        if stripped.trim().is_empty() || is_comment(stripped) {
            continue;
        }
        if line.starts_with("     ") {
            if let Some((_, ref mut text)) = current {
                text.push(' ');
                text.push_str(stripped.trim_start());
                continue;
            }
        }
        if let Some(done) = current.take() {
            out.push(done);
        }
        current = Some((idx + 1, stripped.trim().to_string()));
    }
    if let Some(done) = current {
        out.push(done);
    }
    out
}

/// True for MCNP comment lines: first non-blank character `c`/`C`
/// followed by blank or end of line (so `cut:n` cards are not comments).
pub(crate) fn is_comment(line: &str) -> bool {
    let trimmed = line.trim_start();
    let mut chars = trimmed.chars();
    match chars.next() {
        Some('c') | Some('C') => chars.next().is_none_or(|c| c.is_whitespace()),
        _ => false,
    }
}

/// Parse every cell card in a cell block (lines up to the first blank).
pub fn parse_cells(text: &str) -> Result<Vec<CellCard>, Error> {
    let lines: Vec<&str> = text.lines().collect();
    let mut cards = Vec::new();
    for (lineno, logical) in join_continuations(&lines) {
        let mut card = parse_cell_line(&logical, lineno)?;
        card.raw_lines = vec![logical];
        cards.push(card);
    }
    Ok(cards)
}

/// Parse one logical cell-card line (continuations already joined).
pub fn parse_cell_line(logical: &str, lineno: usize) -> Result<CellCard, Error> {
    let mut tokens = logical.split_whitespace();
    let num: u32 = tokens
        .next()
        .map(str::parse)
        .transpose()
        .map_err(|_| Error::BadCard {
            line: lineno,
            message: format!("invalid cell number in `{logical}`"),
        })?
        .ok_or_else(|| Error::BadCard {
            line: lineno,
            message: "empty cell card".to_string(),
        })?;
    let mat: u32 = tokens
        .next()
        .map(str::parse)
        .transpose()
        .map_err(|_| Error::BadCard {
            line: lineno,
            message: format!("invalid cell material in `{logical}`"),
        })?
        .ok_or_else(|| Error::BadCard {
            line: lineno,
            message: format!("cell {num} is missing its material"),
        })?;
    // Void cells carry no density token.
    let mut rest: Vec<&str> = tokens.collect();
    let dens = if mat == 0 {
        None
    } else {
        let token = rest.first().ok_or_else(|| Error::BadCard {
            line: lineno,
            message: format!("cell {num} is missing its density"),
        })?;
        let dens: f64 = token.parse().map_err(|_| Error::BadCard {
            line: lineno,
            message: format!("invalid cell density `{token}`"),
        })?;
        rest.remove(0);
        Some(dens)
    };
    // Geometry tokens lead; the first `key=value` or bare keyword token
    // (anything starting an ASCII letter that is not part of geometry)
    // starts the parameter list. Geometry tokens never contain `=`.
    let split = rest
        .iter()
        .position(|t| t.contains('=') || PARAM_KEYWORDS.iter().any(|k| t.eq_ignore_ascii_case(k)));
    let (geom_tokens, param_tokens) = match split {
        Some(i) => (&rest[..i], &rest[i..]),
        None => (rest.as_slice(), &[][..]),
    };
    // Real decks glue operators to numbers (`-1:2:-3`, `(1:-2)`).
    let expanded: Vec<String> = geom_tokens
        .iter()
        .flat_map(|t| split_operators(t))
        .collect();
    let expanded_refs: Vec<&str> = expanded.iter().map(String::as_str).collect();
    let geom = parse_geom(&expanded_refs, lineno, num)?;
    Ok(CellCard {
        num,
        mat,
        dens,
        geom,
        params: param_tokens.iter().map(|s| (*s).to_string()).collect(),
        raw_lines: Vec::new(),
        prefix_lines: Vec::new(),
    })
}

/// Bare-keyword cell parameters (everything else uses `key=value`).
const PARAM_KEYWORDS: [&str; 4] = ["vol", "pwt", "ext", "fcl"];

/// Split attached CSG operators off a whitespace token: `-1:2` becomes
/// `-1`, `:`, `2`. Only `:`, `(`, `)` detach; `#`/`*`/signs stay glued.
fn split_operators(token: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut current = String::new();
    for ch in token.chars() {
        if ch == ':' || ch == '(' || ch == ')' {
            if !current.is_empty() {
                out.push(std::mem::take(&mut current));
            }
            out.push(ch.to_string());
        } else {
            current.push(ch);
        }
    }
    if !current.is_empty() {
        out.push(current);
    }
    out
}

/// Recursive-descent CSG parser: union (`:`) binds loosest, juxtaposition
/// intersection tighter, `#`/parens/numbers tightest.
fn parse_geom(tokens: &[&str], lineno: usize, cell: u32) -> Result<GeomExpr, Error> {
    let mut parser = GeomParser {
        tokens,
        pos: 0,
        lineno,
        cell,
    };
    let expr = parser.parse_union()?;
    if parser.pos != tokens.len() {
        return Err(Error::BadCard {
            line: lineno,
            message: format!(
                "trailing geometry token `{}` in cell {cell}",
                tokens[parser.pos]
            ),
        });
    }
    Ok(expr)
}

struct GeomParser<'a> {
    tokens: &'a [&'a str],
    pos: usize,
    lineno: usize,
    cell: u32,
}

impl<'a> GeomParser<'a> {
    fn parse_union(&mut self) -> Result<GeomExpr, Error> {
        let mut left = self.parse_intersect()?;
        while self.peek() == Some(":") {
            self.pos += 1;
            let right = self.parse_intersect()?;
            left = GeomExpr::Union(Box::new(left), Box::new(right));
        }
        Ok(left)
    }

    fn parse_intersect(&mut self) -> Result<GeomExpr, Error> {
        let mut parts = vec![self.parse_unary()?];
        while let Some(tok) = self.peek() {
            if tok == ":" || tok == ")" {
                break;
            }
            parts.push(self.parse_unary()?);
        }
        if parts.len() == 1 {
            Ok(parts.pop().expect("single parsed part"))
        } else {
            Ok(GeomExpr::Intersect(parts))
        }
    }

    fn parse_unary(&mut self) -> Result<GeomExpr, Error> {
        match self.peek() {
            None => Err(Error::BadCard {
                line: self.lineno,
                message: format!("cell {} is missing geometry", self.cell),
            }),
            Some("(") => {
                self.pos += 1;
                let expr = self.parse_union()?;
                match self.peek() {
                    Some(")") => {
                        self.pos += 1;
                        Ok(expr)
                    }
                    _ => Err(Error::BadCard {
                        line: self.lineno,
                        message: format!("unbalanced parenthesis in cell {}", self.cell),
                    }),
                }
            }
            Some(tok) if tok.starts_with('#') => {
                let num: i32 = tok[1..].parse().map_err(|_| Error::BadCard {
                    line: self.lineno,
                    message: format!("invalid cell complement `{tok}`"),
                })?;
                self.pos += 1;
                Ok(GeomExpr::Complement(Box::new(GeomExpr::HalfSpace(
                    HalfSpace {
                        surf: num,
                        reflecting: false,
                    },
                ))))
            }
            Some(tok) => {
                let (reflecting, digits) = match tok.strip_prefix('*') {
                    Some(rest) => (true, rest),
                    None => (false, tok),
                };
                let surf: i32 = digits.parse().map_err(|_| Error::BadCard {
                    line: self.lineno,
                    message: format!("invalid geometry token `{tok}` in cell {}", self.cell),
                })?;
                self.pos += 1;
                Ok(GeomExpr::HalfSpace(HalfSpace { surf, reflecting }))
            }
        }
    }

    fn peek(&self) -> Option<&'a str> {
        self.tokens.get(self.pos).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse_one(logical: &str) -> CellCard {
        parse_cell_line(logical, 1).unwrap()
    }

    #[test]
    fn simple_cell_with_density_and_params() {
        let card = parse_one("1 1 -19.1 -1 imp:n=1");
        assert_eq!(card.num, 1);
        assert_eq!(card.mat, 1);
        assert_eq!(card.dens, Some(-19.1));
        assert_eq!(
            card.geom,
            GeomExpr::HalfSpace(HalfSpace {
                surf: -1,
                reflecting: false
            })
        );
        assert_eq!(card.params, vec!["imp:n=1".to_string()]);
    }

    #[test]
    fn void_cell_has_no_density() {
        let card = parse_one("1 0 -1:2:-3 imp:n=0");
        assert_eq!(card.dens, None);
        assert_eq!(
            card.geom,
            GeomExpr::Union(
                Box::new(GeomExpr::Union(
                    Box::new(GeomExpr::HalfSpace(HalfSpace {
                        surf: -1,
                        reflecting: false
                    })),
                    Box::new(GeomExpr::HalfSpace(HalfSpace {
                        surf: 2,
                        reflecting: false
                    })),
                )),
                Box::new(GeomExpr::HalfSpace(HalfSpace {
                    surf: -3,
                    reflecting: false
                })),
            )
        );
    }

    #[test]
    fn intersect_binds_tighter_than_union() {
        // 1 -2 : 3  ==  (1 -2) : (3)
        let card = parse_one("5 0 1 -2 : 3");
        assert_eq!(
            card.geom,
            GeomExpr::Union(
                Box::new(GeomExpr::Intersect(vec![
                    GeomExpr::HalfSpace(HalfSpace {
                        surf: 1,
                        reflecting: false
                    }),
                    GeomExpr::HalfSpace(HalfSpace {
                        surf: -2,
                        reflecting: false
                    }),
                ])),
                Box::new(GeomExpr::HalfSpace(HalfSpace {
                    surf: 3,
                    reflecting: false
                })),
            )
        );
        assert_eq!(card.geom.render(), "1 -2 : 3");
    }

    #[test]
    fn parentheses_and_complement() {
        let card = parse_one("2 0 (1 : -2) 3");
        assert_eq!(card.geom.render(), "(1 : -2) 3");
        let card = parse_one("3 0 #5 -1");
        assert_eq!(
            card.geom,
            GeomExpr::Intersect(vec![
                GeomExpr::Complement(Box::new(GeomExpr::HalfSpace(HalfSpace {
                    surf: 5,
                    reflecting: false
                }))),
                GeomExpr::HalfSpace(HalfSpace {
                    surf: -1,
                    reflecting: false
                }),
            ])
        );
        assert_eq!(card.geom.render(), "#5 -1");
    }

    #[test]
    fn reflecting_marker_round_trips() {
        let card = parse_one("4 1 -2.7 *-10 20");
        assert_eq!(card.geom.render(), "*-10 20");
    }

    #[test]
    fn attached_colon_tokens_split() {
        // Real decks glue colons to numbers: -1:2:-3.
        let card = parse_one("1 0 -1:2:-3");
        assert_eq!(card.geom.render(), "-1 : 2 : -3");
    }

    #[test]
    fn errors_have_lines() {
        assert!(parse_cell_line("", 7).is_err());
        assert!(parse_cell_line("1", 7).is_err());
        assert!(parse_cell_line("1 1", 7).is_err());
        assert!(parse_cell_line("1 1 -2.7 (1 : 2", 7).is_err());
        assert!(parse_cell_line("1 1 -2.7 foo", 7).is_err());
    }
}
