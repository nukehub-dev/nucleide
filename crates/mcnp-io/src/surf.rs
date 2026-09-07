//! MCNP surface-card parsing.
//!
//! A surface card is `[*]num [transform] KIND coeffs...`: an optional `*`
//! reflecting marker, an optional transform number, a surface-type keyword,
//! and numeric coefficients. See [`SurfKind::arity`] for the enforced
//! coefficient counts.

use crate::cell::join_continuations;
use crate::inp::Error;

/// MCNP surface types: planes, quadrics, and macrobodies.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfKind {
    /// Plane normal to x/y/z (`PX`, `PY`, `PZ`).
    Px,
    Py,
    Pz,
    /// General plane (`P`).
    P,
    /// Sphere centered at origin / general / x/y/z-shifted (`S`, `SO`,
    /// `SX`, `SY`, `SZ`).
    S,
    So,
    Sx,
    Sy,
    Sz,
    /// Cylinder parallel to x/y/z (`CX`, `CY`, `CZ`).
    Cx,
    Cy,
    Cz,
    /// Cone parallel to x/y/z (`KX`, `KY`, `KZ`).
    Kx,
    Ky,
    Kz,
    /// Quadrics (`SQ`, `GQ`).
    Sq,
    Gq,
    /// Tori parallel to x/y/z (`TX`, `TY`, `TZ`).
    Tx,
    Ty,
    Tz,
    /// Planes `x/y/z = const` (`X`, `Y`, `Z`).
    X,
    Y,
    Z,
    /// Macrobodies: box, hexahedron, wedge, prism, ellipsoid, truncated
    /// cone, cylinder, torus (`BOX`, `RHP`, `HEX`, `WED`, `RPP`, `SPH`,
    /// `RCC`, `REC`, `TRC`, `ELL`, `ARB`).
    McBox,
    Rhp,
    Hex,
    Wed,
    Rpp,
    Sph,
    Rcc,
    Rec,
    Trc,
    Ell,
    Arb,
}

impl SurfKind {
    /// Parse a surface-type keyword (case-insensitive).
    pub fn from_keyword(keyword: &str) -> Option<SurfKind> {
        Some(match keyword.to_ascii_uppercase().as_str() {
            "PX" => SurfKind::Px,
            "PY" => SurfKind::Py,
            "PZ" => SurfKind::Pz,
            "P" => SurfKind::P,
            "S" => SurfKind::S,
            "SO" => SurfKind::So,
            "SX" => SurfKind::Sx,
            "SY" => SurfKind::Sy,
            "SZ" => SurfKind::Sz,
            "CX" => SurfKind::Cx,
            "CY" => SurfKind::Cy,
            "CZ" => SurfKind::Cz,
            "KX" => SurfKind::Kx,
            "KY" => SurfKind::Ky,
            "KZ" => SurfKind::Kz,
            "SQ" => SurfKind::Sq,
            "GQ" => SurfKind::Gq,
            "TX" => SurfKind::Tx,
            "TY" => SurfKind::Ty,
            "TZ" => SurfKind::Tz,
            "X" => SurfKind::X,
            "Y" => SurfKind::Y,
            "Z" => SurfKind::Z,
            "BOX" => SurfKind::McBox,
            "RHP" => SurfKind::Rhp,
            "HEX" => SurfKind::Hex,
            "WED" => SurfKind::Wed,
            "RPP" => SurfKind::Rpp,
            "SPH" => SurfKind::Sph,
            "RCC" => SurfKind::Rcc,
            "REC" => SurfKind::Rec,
            "TRC" => SurfKind::Trc,
            "ELL" => SurfKind::Ell,
            "ARB" => SurfKind::Arb,
            _ => return None,
        })
    }

    /// Canonical keyword spelling.
    pub fn keyword(self) -> &'static str {
        match self {
            SurfKind::Px => "PX",
            SurfKind::Py => "PY",
            SurfKind::Pz => "PZ",
            SurfKind::P => "P",
            SurfKind::S => "S",
            SurfKind::So => "SO",
            SurfKind::Sx => "SX",
            SurfKind::Sy => "SY",
            SurfKind::Sz => "SZ",
            SurfKind::Cx => "CX",
            SurfKind::Cy => "CY",
            SurfKind::Cz => "CZ",
            SurfKind::Kx => "KX",
            SurfKind::Ky => "KY",
            SurfKind::Kz => "KZ",
            SurfKind::Sq => "SQ",
            SurfKind::Gq => "GQ",
            SurfKind::Tx => "TX",
            SurfKind::Ty => "TY",
            SurfKind::Tz => "TZ",
            SurfKind::X => "X",
            SurfKind::Y => "Y",
            SurfKind::Z => "Z",
            SurfKind::McBox => "BOX",
            SurfKind::Rhp => "RHP",
            SurfKind::Hex => "HEX",
            SurfKind::Wed => "WED",
            SurfKind::Rpp => "RPP",
            SurfKind::Sph => "SPH",
            SurfKind::Rcc => "RCC",
            SurfKind::Rec => "REC",
            SurfKind::Trc => "TRC",
            SurfKind::Ell => "ELL",
            SurfKind::Arb => "ARB",
        }
    }

    /// Required coefficient count, or `None` when the kind accepts a
    /// variable tail (kept permissive for exotic macrobodies).
    pub fn arity(self) -> Option<usize> {
        Some(match self {
            SurfKind::Px
            | SurfKind::Py
            | SurfKind::Pz
            | SurfKind::X
            | SurfKind::Y
            | SurfKind::Z => 1,
            SurfKind::Cx | SurfKind::Cy | SurfKind::Cz => 1,
            SurfKind::So => 1,
            SurfKind::Sx | SurfKind::Sy | SurfKind::Sz => 2,
            SurfKind::S => 4,
            SurfKind::Kx | SurfKind::Ky | SurfKind::Kz => 5,
            SurfKind::Tx | SurfKind::Ty | SurfKind::Tz => 6,
            SurfKind::P => 9,
            SurfKind::Sq => 10,
            SurfKind::Gq => 16,
            SurfKind::Rpp => 6,
            SurfKind::Sph => 4,
            SurfKind::Rcc => 7,
            SurfKind::Trc => 8,
            SurfKind::Ell => 7,
            SurfKind::McBox => 12,
            SurfKind::Rec => 12,
            SurfKind::Wed => 12,
            SurfKind::Rhp | SurfKind::Hex => 15,
            SurfKind::Arb => return None,
        })
    }
}

/// One parsed MCNP surface card.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfCard {
    /// Surface number.
    pub num: u32,
    /// `*` reflecting marker.
    pub reflecting: bool,
    /// Transform number when the card carries one.
    pub transform: Option<u32>,
    /// Surface type.
    pub kind: SurfKind,
    /// Numeric coefficients.
    pub coeffs: Vec<f64>,
    /// Source lines forming this card (for format-preserving write-back).
    pub raw_lines: Vec<String>,
    /// Comment/blank lines preceding this card in its block, verbatim.
    pub prefix_lines: Vec<String>,
}

impl SurfCard {
    /// Render in canonical MCNP form.
    pub fn render(&self) -> String {
        let mut out = format!(
            "{}{}{} {}",
            if self.reflecting { "*" } else { "" },
            self.num,
            match self.transform {
                Some(t) => format!(" {t}"),
                None => String::new(),
            },
            self.kind.keyword()
        );
        for c in &self.coeffs {
            out.push_str(&format!(" {c}"));
        }
        out
    }
}

/// Parse every surface card in a surface block (lines up to the first
/// blank).
pub fn parse_surfs(text: &str) -> Result<Vec<SurfCard>, Error> {
    let lines: Vec<&str> = text.lines().collect();
    let mut cards = Vec::new();
    for (lineno, logical) in join_continuations(&lines) {
        let mut card = parse_surf_line(&logical, lineno)?;
        card.raw_lines = vec![logical];
        cards.push(card);
    }
    Ok(cards)
}

/// Parse one logical surface-card line (continuations already joined).
pub fn parse_surf_line(logical: &str, lineno: usize) -> Result<SurfCard, Error> {
    let mut tokens = logical.split_whitespace();
    let first = tokens.next().ok_or_else(|| Error::BadCard {
        line: lineno,
        message: "empty surface card".to_string(),
    })?;
    let (reflecting, num_text) = match first.strip_prefix('*') {
        Some(rest) => (true, rest),
        None => (false, first),
    };
    let num: u32 = num_text.parse().map_err(|_| Error::BadCard {
        line: lineno,
        message: format!("invalid surface number `{first}`"),
    })?;
    let rest: Vec<&str> = tokens.collect();
    if rest.is_empty() {
        return Err(Error::BadCard {
            line: lineno,
            message: format!("surface {num} is missing its type"),
        });
    }
    // Optional transform number between the id and the type keyword.
    let (transform, kind_token, coeff_tokens) = match SurfKind::from_keyword(rest[0]) {
        Some(_) => (None, rest[0], &rest[1..]),
        None => {
            let transform: u32 = rest[0].parse().map_err(|_| Error::BadCard {
                line: lineno,
                message: format!("invalid surface type `{}` on surface {num}", rest[0]),
            })?;
            let kind_token = *rest.get(1).ok_or_else(|| Error::BadCard {
                line: lineno,
                message: format!("surface {num} is missing its type"),
            })?;
            (Some(transform), kind_token, &rest[2..])
        }
    };
    let kind = SurfKind::from_keyword(kind_token).ok_or_else(|| Error::BadCard {
        line: lineno,
        message: format!("invalid surface type `{kind_token}` on surface {num}"),
    })?;
    let mut coeffs = Vec::with_capacity(coeff_tokens.len());
    for token in coeff_tokens {
        coeffs.push(token.parse::<f64>().map_err(|_| Error::BadCard {
            line: lineno,
            message: format!("invalid surface coefficient `{token}` on surface {num}"),
        })?);
    }
    if let Some(arity) = kind.arity() {
        if coeffs.len() != arity {
            return Err(Error::BadCard {
                line: lineno,
                message: format!(
                    "surface {num} type {} needs {arity} coefficients, found {}",
                    kind.keyword(),
                    coeffs.len()
                ),
            });
        }
    }
    Ok(SurfCard {
        num,
        reflecting,
        transform,
        kind,
        coeffs,
        raw_lines: Vec::new(),
        prefix_lines: Vec::new(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn planes_and_transforms() {
        let card = parse_surf_line("1 px -5.0", 1).unwrap();
        assert_eq!(card.num, 1);
        assert_eq!(card.kind, SurfKind::Px);
        assert_eq!(card.coeffs, vec![-5.0]);
        assert_eq!(card.transform, None);
        assert!(!card.reflecting);

        // `99 7 PX 180`: surface with a transform link (regression: this is
        // not a cell card).
        let card = parse_surf_line("99 7 PX 180", 1).unwrap();
        assert_eq!(card.num, 99);
        assert_eq!(card.transform, Some(7));
        assert_eq!(card.coeffs, vec![180.0]);
    }

    #[test]
    fn reflecting_and_macrobodies() {
        let card = parse_surf_line("*1 so 100", 1).unwrap();
        assert!(card.reflecting);
        assert_eq!(card.kind, SurfKind::So);
        assert_eq!(card.render(), "*1 SO 100");

        let card = parse_surf_line("10 rpp -5 5 -5 5 -5 5", 1).unwrap();
        assert_eq!(card.kind, SurfKind::Rpp);
        assert_eq!(card.coeffs.len(), 6);

        let card = parse_surf_line("11 sph 0 0 0 10", 1).unwrap();
        assert_eq!(card.render(), "11 SPH 0 0 0 10");
    }

    #[test]
    fn arity_is_checked() {
        assert!(parse_surf_line("1 px 1 2", 3).is_err());
        assert!(parse_surf_line("1 rpp 1 2 3", 3).is_err());
        assert!(parse_surf_line("1 zz 1", 3).is_err());
        assert!(parse_surf_line("1", 3).is_err());
        assert!(parse_surf_line("x px 1", 3).is_err());
    }
}
