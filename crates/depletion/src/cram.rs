//! IPF (incomplete partial factorization) CRAM solvers, orders 16 and 48.
//!
//! Chebyshev rational-approximation coefficients for the CRAM matrix exponent (Pusa & Li).
//! Algorithm per timestep: with `A' = dt*A`, factorize `A' - theta_k I` for
//! each pole, then iterate
//! `y += 2*Re(alpha_k * (A' - theta_k I)^-1 y)` for all k and scale by
//! `alpha0`. Symbolic LU analysis is reused across all poles.

use linalg::{ComplexCsc, ComplexLu, SymbolicLu, C64};

use crate::matrix::DepletionSystem;

/// CRAM approximation order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Order {
    /// Order 16 — cheaper, adequate for short chains.
    Order16,
    /// Order 48 — production accuracy.
    Order48,
}

// ---- order 16 coefficients ----
#[allow(clippy::excessive_precision)] // published CRAM constants
const C16_ALPHA: [C64; 8] = [
    C64 {
        re: 5.464930576870210e3,
        im: -3.797983575308356e4,
    },
    C64 {
        re: 9.045112476907548e1,
        im: -1.115537522430261e3,
    },
    C64 {
        re: 2.344818070467641e2,
        im: -4.228020157070496e2,
    },
    C64 {
        re: 9.453304067358312e1,
        im: -2.951294291446048e2,
    },
    C64 {
        re: 7.283792954673409e2,
        im: -1.205646080220011e5,
    },
    C64 {
        re: 3.648229059594851e1,
        im: -1.155509621409682e2,
    },
    C64 {
        re: 2.547321630156819e1,
        im: -2.639500283021502e1,
    },
    C64 {
        re: 2.394538338734709e1,
        im: -5.650522971778156e0,
    },
];

const C16_THETA: [C64; 8] = [
    C64 {
        re: 3.509103608414918,
        im: 8.436198985884374,
    },
    C64 {
        re: 5.948152268951177,
        im: 3.587457362018322,
    },
    C64 {
        re: -5.264971343442647,
        im: 16.22022147316793,
    },
    C64 {
        re: 1.419375897185666,
        im: 10.92536348449672,
    },
    C64 {
        re: 6.416177699099435,
        im: 1.194122393370139,
    },
    C64 {
        re: 4.993174737717997,
        im: 5.996881713603942,
    },
    C64 {
        re: -1.413928462488886,
        im: 13.49772569889275,
    },
    C64 {
        re: -10.84391707869699,
        im: 19.27744616718165,
    },
];

const C16_ALPHA0: f64 = 2.124853710495224e-16;

// ---- order 48 coefficients ----
#[rustfmt::skip]
#[allow(clippy::excessive_precision)] // published CRAM constants
const C48_THETA_RE: [f64; 24] = [
    -4.465731934165702e+01, -5.284616241568964e+00,
    -8.867715667624458e+00,3.493013124279215e+00,
   1.564102508858634e+01,1.742097597385893e+01,
    -2.834466755180654e+01,1.661569367939544e+01,
   8.011836167974721e+00, -2.056267541998229e+00,
   1.449208170441839e+01,1.853807176907916e+01,
   9.932562704505182e+00, -2.244223871767187e+01,
   8.590014121680897e-01, -1.286192925744479e+01,
   1.164596909542055e+01,1.806076684783089e+01,
   5.870672154659249e+00, -3.542938819659747e+01,
   1.901323489060250e+01,1.885508331552577e+01,
    -1.734689708174982e+01,1.316284237125190e+01,
];

#[rustfmt::skip]
#[allow(clippy::excessive_precision)] // published CRAM constants
const C48_THETA_IM: [f64; 24] = [
   6.233225190695437e+01,4.057499381311059e+01,
   4.325515754166724e+01,3.281615453173585e+01,
   1.558061616372237e+01,1.076629305714420e+01,
   5.492841024648724e+01,1.316994930024688e+01,
   2.780232111309410e+01,3.794824788914354e+01,
   1.799988210051809e+01,5.974332563100539e+00,
   2.532823409972962e+01,5.179633600312162e+01,
   3.536456194294350e+01,4.600304902833652e+01,
   2.287153304140217e+01,8.368200580099821e+00,
   3.029700159040121e+01,5.834381701800013e+01,
   1.194282058271408e+00,3.583428564427879e+00,
   4.883941101108207e+01,2.042951874827759e+01,
];

#[rustfmt::skip]
#[allow(clippy::excessive_precision)] // published CRAM constants
const C48_ALPHA_RE: [f64; 24] = [
   6.387380733878774e+02,1.909896179065730e+02,
   4.236195226571914e+02,4.645770595258726e+02,
   7.765163276752433e+02,1.907115136768522e+03,
   2.909892685603256e+03,1.944772206620450e+02,
   1.382799786972332e+05,5.628442079602433e+03,
   2.151681283794220e+02,1.324720240514420e+03,
   1.617548476343347e+04,1.112729040439685e+02,
   1.074624783191125e+02,8.835727765158191e+01,
   9.354078136054179e+01,9.418142823531573e+01,
   1.040012390717851e+02,6.861882624343235e+01,
   8.766654491283722e+01,1.056007619389650e+02,
   7.738987569039419e+01,1.041366366475571e+02,
];

#[rustfmt::skip]
#[allow(clippy::excessive_precision)] // published CRAM constants
const C48_ALPHA_IM: [f64; 24] = [
    -6.743912502859256e+02, -3.973203432721332e+02,
    -2.041233768918671e+03, -1.652917287299683e+03,
    -1.783617639907328e+04, -5.887068595142284e+04,
    -9.953255345514560e+03, -1.427131226068449e+03,
    -3.256885197214938e+06, -2.924284515884309e+04,
    -1.121774011188224e+03, -6.370088443140973e+04,
    -1.008798413156542e+06, -8.837109731680418e+01,
    -1.457246116408180e+02, -6.388286188419360e+01,
    -2.195424319460237e+02, -6.719055740098035e+02,
    -1.693747595553868e+02, -1.177598523430493e+01,
    -4.596464999363902e+03, -1.738294585524067e+03,
    -4.311715386228984e+01, -2.777743732451969e+02,
];

const C48_ALPHA0: f64 = 2.258038182743983e-47;

fn coefficients(order: Order) -> (&'static [C64], &'static [C64], f64) {
    match order {
        Order::Order16 => (&C16_ALPHA, &C16_THETA, C16_ALPHA0),
        Order::Order48 => (alpha48(), theta48(), C48_ALPHA0),
    }
}

/// Lazily-built order-48 alpha residues.
fn alpha48() -> &'static [C64; 24] {
    static CELL: std::sync::OnceLock<[C64; 24]> = std::sync::OnceLock::new();
    CELL.get_or_init(|| {
        let mut out = [C64 { re: 0.0, im: 0.0 }; 24];
        for i in 0..24 {
            out[i] = C64 {
                re: C48_ALPHA_RE[i],
                im: C48_ALPHA_IM[i],
            };
        }
        out
    })
}

/// Lazily-built order-48 poles.
fn theta48() -> &'static [C64; 24] {
    static CELL: std::sync::OnceLock<[C64; 24]> = std::sync::OnceLock::new();
    CELL.get_or_init(|| {
        let mut out = [C64 { re: 0.0, im: 0.0 }; 24];
        for i in 0..24 {
            out[i] = C64 {
                re: C48_THETA_RE[i],
                im: C48_THETA_IM[i],
            };
        }
        out
    })
}

/// Errors from CRAM evaluation.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// Backend/factorization failure.
    Linalg(String),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "CRAM error: {}",
            match self {
                Error::Linalg(m) => m,
            }
        )
    }
}
impl std::error::Error for Error {}

/// Solve one depletion step `n(dt) = exp(A*dt) n0` with IPF CRAM.
///
/// `n0` is the initial nuclide vector (atoms). Returns the vector after
/// `dt` seconds. The symbolic analysis is computed internally once; use
/// [`cram_with_symbolic`] to hoist it across many timesteps.
pub fn cram(sys: &DepletionSystem, order: Order, n0: &[f64], dt: f64) -> Result<Vec<f64>, Error> {
    let sym = SymbolicLu::try_new(&sys.pattern).map_err(|e| Error::Linalg(e.to_string()))?;
    cram_with_symbolic(sys, &sym, order, n0, dt)
}

/// [`cram`] with a caller-managed symbolic factorization.
pub fn cram_with_symbolic(
    sys: &DepletionSystem,
    sym: &SymbolicLu,
    order: Order,
    n0: &[f64],
    dt: f64,
) -> Result<Vec<f64>, Error> {
    if n0.len() != sys.pattern.nrows() {
        return Err(Error::Linalg(format!(
            "n0 length {} != chain size {}",
            n0.len(),
            sys.pattern.nrows()
        )));
    }

    let (alphas, thetas, alpha0) = coefficients(order);

    // Build and factorize A*dt - theta*I for every pole.
    let mut factors: Vec<ComplexLu> = Vec::with_capacity(thetas.len());
    for theta in thetas.iter().copied() {
        let values = sys.shifted_values(dt, theta);
        let mat = ComplexCsc::from_entries(&sys.pattern, &values)
            .map_err(|e| Error::Linalg(e.to_string()))?;
        factors.push(
            ComplexLu::try_new_with_symbolic(sym, &mat)
                .map_err(|e| Error::Linalg(e.to_string()))?,
        );
    }

    // IPF iteration on the working vector. NOTE: only the real part of
    // alpha*x is accumulated (matching the reference real-part semantics); keeping y real
    // is intrinsic to Pusa's incomplete partial factorization scheme.
    let mut y: Vec<C64> = n0.iter().map(|v| C64 { re: *v, im: 0.0 }).collect();
    for (alpha, lu) in alphas.iter().zip(&factors) {
        let x = lu.solve(&y);
        for (yi, xi) in y.iter_mut().zip(x) {
            yi.re += 2.0 * (alpha.re * xi.re - alpha.im * xi.im);
        }
    }
    for yi in &mut y {
        *yi = C64 {
            re: yi.re * alpha0,
            im: yi.im,
        };
    }

    Ok(y.into_iter().map(|v| v.re).collect())
}
