//! Sublet S1+S2+S4+S5+S6+S7 radiological totals: total activity, decay heat,
//! committed ingestion/inhalation hazards, the transport ratio, and the IAEA
//! clearance-index variant.
//!
//! Pure arithmetic over caller-supplied activation inventories, pinned to two
//! open sources (never the paywalled journal pages):
//!
//! - **S1 — total activity**: the inventory-table row `activity, Ai = Ni λi,
//!   Bq` of Table X ("Inventory table output definitions", §II L.1 "Step
//!   Output" of §II L "Inventory Run Output") of CCFE-PR(16)53 (Sublet,
//!   Eastwood, Morgan, Gilbert, Fleming, Arter — the author preprint of the
//!   FISPACT-II system paper), whose totals paragraph (§II L.1, pp. 18–19)
//!   states the activity, dose, and heating-power totals are sums over all
//!   nuclides of the Table X quantities (validation anchor: Table XVII total
//!   activity in Bq kg−1 for copper cooling times). The α/β/γ split follows
//!   the "Activity break-down, fission and hazards" section of the open
//!   `output_interpretation` reference (`fispact.github.io`): IRT 4 → alpha;
//!   IRT 1, 2, 11, 14, 16, 17, 19, 20 → beta; IRT 3 → gamma; IRT 12, 13 split
//!   α/β and IRT 15 splits α/γ by caller-supplied branch fractions. IRT
//!   classes are Table VI ("Decay Types (MF=8, MT=457)", IRT 1–27 with 8, 9
//!   unused and 10 unknown) of the same preprint.
//! - **S2 — decay heat**: the Table X rows `β-power, Ai Eβ,i C1, kW`,
//!   `α-power, Ai Eα,i C1, kW`, `γ-power, Ai Eγ,i C1, kW` (`C1` = eV→kJ
//!   conversion; supporting paragraph §II F.1.a "Decay heating", cf.
//!   §II L), with the per-class kW columns (cols 4/5/6) and the `TOTAL
//!   ALPHA/BETA/GAMMA HEAT PRODUCTION` + `TOTAL HEAT PRODUCTION` + `TOTAL
//!   HEAT EX TRITIUM` sums of the "Time line and nuclide inventory" and
//!   "Activity break-down" sections of `output_interpretation`. Radiation
//!   classes are Table VII ("Decay Radiation Types (MF=8, MT=457)", STYP
//!   0–9: 0 γ, 1 β−, 2 ec/β+, 3 not known, 4 α, 5 n, 6 SF, 7 p, 8 e−, 9 x)
//!   of the same preprint.
//!
//! - **S4 — ingestion hazard** and **S5 — inhalation hazard**: the Table X
//!   inventory row `activity, Ai = Ni λi, Bq` folded with the caller-supplied
//!   50-year committed dose coefficients `e^ing_i` / `e^inh_i` (Sv/Bq) into
//!   the wide-table `INGESTION DOSE(Sv)` / `INHALATION DOSE(Sv)` columns
//!   (cols 8/9 of the radiological block) and their `INGESTION HAZARD FOR
//!   ALL MATERIALS` / `INHALATION HAZARD FOR ALL MATERIALS` totals of the
//!   open `output_interpretation` reference (`fispact.github.io`) and the
//!   `HAZARDS` keyword page. Dose coefficients are ICRP-copyrighted and so
//!   stay caller inputs, never vendored.
//! - **S6 — transport ratio**: the Table X inventory row `activity, Ai =
//!   Ni λi, Bq` folded over the caller-supplied transport limits `A2,i`
//!   (TBq) through Eqs 78/79: the per-nuclide `Bq/A2` term `Ai/(A2,i·C2)`
//!   (`C2` = TBq→Bq conversion) sums to the `Total Bq/A2` ratio of the
//!   `ATWO` output block of the open `output_interpretation` reference
//!   (`fispact.github.io`), and the effective A2 `hA2 = ΣAi/ratio` (TBq) is
//!   defined by that ratio, exactly by construction. Transport limits are
//!   regulation-tabled values and so stay caller inputs, never vendored.
//! - **S7 — IAEA clearance index**: the Table X inventory row `activity, Ai =
//!   Ni λi, Bq` folded with the caller-supplied IAEA clearance levels `Li`
//!   (Bq/kg) and the total mass `Mtot` (kg) through Eq. 80 — the per-nuclide
//!   `Ai/(Mtot·Li)` term sums to the `CLEARANCE INDEX` of the `CLEAR` output
//!   block of the open `output_interpretation` reference (`fispact.github.io`).
//!   Same arithmetic as the EU-table [`crate::clearance`] variant — only the
//!   limit table differs (the IAEA levels are permission-gated, so they stay
//!   caller inputs, never vendored). The `<= 1` screening boundary classifies
//!   both sides, boundary included.
//!
//! The S1/S2/S4/S5 totals carry the `TOTAL ... EXCLUDING TRITIUM` companion
//! (total minus the tritium entry). Per-nuclide decay energies and dose
//! coefficients are caller-supplied, never vendored. Parts always sum to the
//! total by construction (`total == alpha + beta + gamma`, exactly in `f64`;
//! each hazard total is the exact running `ΣAi·e_i`).
//!
//! Explicitly out of scope: ICRP/IAEA coefficient tables (caller inputs,
//! never vendored), bremsstrahlung correction, DPA/KERMA/gas (`damage` owns),
//! fission count / burn-up.
//! (The S3 gamma dose-rate kernel lives in [`crate::dose`].)

use crate::clearance::ClearanceClass;
use nucleide_nuclei::NuclideId;

use crate::error::{Error, Result};

/// Conversion factor C1 from eV to kJ (Table X: `Ai E C1`, kW).
///
/// `Ai` (decays/s) × `E` (eV/decay) × C1 (kJ/eV) = kJ/s = kW. The exact SI
/// elementary charge `1.602176634e-19` J/eV scaled to kJ.
pub const EV_TO_KJ: f64 = 1.602_176_634e-22;

/// Cached tritium identity for the ex-tritium companion totals.
fn tritium_id() -> NuclideId {
    static TRITIUM: std::sync::OnceLock<NuclideId> = std::sync::OnceLock::new();
    *TRITIUM.get_or_init(|| NuclideId::from_name("H3").expect("nucleide-nuclei must resolve H3"))
}

/// One inventory row for the S1 total-activity sum.
///
/// `irt` is the Table VI decay-type identifier (1–27); `alpha_frac` carries
/// the caller's α branch fraction for the split IRTs only (12, 13 → α/β;
/// 15 → α/γ) and must be `None` for wholly-mapped IRTs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ActivityEntry {
    /// Nuclide carrying this activity (tritium feeds the ex-tritium total).
    pub nuclide: NuclideId,
    /// Per-nuclide activity `Ai = Ni λi` in Bq (finite, `>= 0`).
    pub activity_bq: f64,
    /// Table VI decay-type identifier.
    pub irt: u8,
    /// α branch fraction in `[0, 1]` for IRT 12, 13, 15; `None` otherwise.
    pub alpha_frac: Option<f64>,
}

/// S1 total activity with the α/β/γ IRT split (all Bq).
///
/// `total_bq == alpha_bq + beta_bq + gamma_bq` exactly; `ex_tritium_bq` is the
/// total minus the tritium entry (zero tritium ⇒ equals the total).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ActivityTotal {
    /// `TOTAL ACTIVITY FOR ALL MATERIALS`: ΣAi over all entries, Bq.
    pub total_bq: f64,
    /// `ALPHA BECQUERELS`: IRT 4 plus the α share of IRT 12, 13, 15, Bq.
    pub alpha_bq: f64,
    /// `BETA BECQUERELS`: IRT 1, 2, 11, 14, 16, 17, 19, 20 plus the β share
    /// of IRT 12, 13, Bq.
    pub beta_bq: f64,
    /// `GAMMA BECQUERELS`: IRT 3 plus the γ share of IRT 15, Bq.
    pub gamma_bq: f64,
    /// `TOTAL ACTIVITY EXCLUDING TRITIUM`, Bq.
    pub ex_tritium_bq: f64,
}

/// Whole-vs-split mapping of one Table VI IRT to the S1 classes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum IrtMapping {
    /// Whole activity counts as alpha (IRT 4).
    Alpha,
    /// Whole activity counts as beta (IRT 1, 2, 11, 14, 16, 17, 19, 20).
    Beta,
    /// Whole activity counts as gamma (IRT 3).
    Gamma,
    /// Split α/β by the caller's fraction (IRT 12, 13).
    AlphaBeta,
    /// Split α/γ by the caller's fraction (IRT 15).
    AlphaGamma,
}

/// Map a Table VI IRT to its S1 class, or `None` when the open
/// `output_interpretation` pin defines no mapping (8, 9 unused; 10 unknown;
/// 5–7, 18, 21–27 unlisted in the split paragraph).
fn irt_mapping(irt: u8) -> Option<IrtMapping> {
    match irt {
        4 => Some(IrtMapping::Alpha),
        1 | 2 | 11 | 14 | 16 | 17 | 19 | 20 => Some(IrtMapping::Beta),
        3 => Some(IrtMapping::Gamma),
        12 | 13 => Some(IrtMapping::AlphaBeta),
        15 => Some(IrtMapping::AlphaGamma),
        _ => None,
    }
}

/// S1 total activity `ΣAi` with the Source-B IRT→α/β/γ split (Bq).
///
/// Every entry's activity must be finite and `>= 0`; every IRT must carry a
/// pinned mapping ([`irt_mapping`]); split IRTs (12, 13, 15) require a finite
/// `alpha_frac` in `[0, 1]` while wholly-mapped IRTs reject one. Violations
/// are loud [`Error`]s, never silent drops. The split remainder is `Ai − α`
/// (not `Ai·(1−f)`) so each entry conserves exactly, and the reported total
/// is the class sum, so [`ActivityTotal`] parts always sum to the total.
pub fn total_activity(entries: &[ActivityEntry]) -> Result<ActivityTotal> {
    let tritium = tritium_id();
    let mut alpha = 0.0f64;
    let mut beta = 0.0f64;
    let mut gamma = 0.0f64;
    let mut tritium_bq = 0.0f64;
    for entry in entries {
        let name = entry.nuclide.to_name();
        if !entry.activity_bq.is_finite() || entry.activity_bq < 0.0 {
            return Err(Error::BadActivityValue {
                nuclide: name,
                msg: format!(
                    "activity must be finite and >= 0, got {}",
                    entry.activity_bq
                ),
            });
        }
        let mapping = irt_mapping(entry.irt).ok_or_else(|| Error::UnmappedIrt {
            nuclide: name.clone(),
            irt: entry.irt,
        })?;
        match mapping {
            IrtMapping::Alpha => {
                reject_split(entry)?;
                alpha += entry.activity_bq;
            }
            IrtMapping::Beta => {
                reject_split(entry)?;
                beta += entry.activity_bq;
            }
            IrtMapping::Gamma => {
                reject_split(entry)?;
                gamma += entry.activity_bq;
            }
            IrtMapping::AlphaBeta | IrtMapping::AlphaGamma => {
                let frac = entry.alpha_frac.ok_or_else(|| Error::BadIrtSplit {
                    nuclide: name.clone(),
                    irt: entry.irt,
                    msg: "split IRT requires alpha_frac".to_string(),
                })?;
                if !frac.is_finite() || frac < 0.0 || frac > 1.0 {
                    return Err(Error::BadIrtSplit {
                        nuclide: name.clone(),
                        irt: entry.irt,
                        msg: format!("alpha_frac must be finite and in [0, 1], got {frac}"),
                    });
                }
                let alpha_part = entry.activity_bq * frac;
                let rest = entry.activity_bq - alpha_part;
                alpha += alpha_part;
                if mapping == IrtMapping::AlphaBeta {
                    beta += rest;
                } else {
                    gamma += rest;
                }
            }
        }
        if entry.nuclide == tritium {
            tritium_bq += entry.activity_bq;
        }
    }
    let total = alpha + beta + gamma;
    if !total.is_finite() {
        return Err(Error::BadActivityValue {
            nuclide: "total".to_string(),
            msg: format!(
                "activity sum overflowed to {total} over {} entries",
                entries.len()
            ),
        });
    }
    Ok(ActivityTotal {
        total_bq: total,
        alpha_bq: alpha,
        beta_bq: beta,
        gamma_bq: gamma,
        ex_tritium_bq: total - tritium_bq,
    })
}

/// Reject a caller split fraction on a wholly-mapped IRT.
fn reject_split(entry: &ActivityEntry) -> Result<()> {
    if entry.alpha_frac.is_some() {
        return Err(Error::BadIrtSplit {
            nuclide: entry.nuclide.to_name(),
            irt: entry.irt,
            msg: "wholly-mapped IRT takes no alpha_frac (must be None)".to_string(),
        });
    }
    Ok(())
}

/// One inventory row for the S2 decay-heat sum.
///
/// Average decay energies per radiation class are caller-supplied in eV
/// (never vendored); classes with no emission carry `0.0`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DecayHeatEntry {
    /// Nuclide carrying this activity (tritium feeds the ex-tritium total).
    pub nuclide: NuclideId,
    /// Per-nuclide activity `Ai = Ni λi` in Bq (finite, `>= 0`).
    pub activity_bq: f64,
    /// Average α decay energy `Eα,i` in eV (finite, `>= 0`).
    pub e_alpha_ev: f64,
    /// Average β decay energy `Eβ,i` in eV (finite, `>= 0`).
    pub e_beta_ev: f64,
    /// Average γ decay energy `Eγ,i` in eV (finite, `>= 0`).
    pub e_gamma_ev: f64,
}

/// S2 decay heat with the α/β/γ split (all kW).
///
/// `total_kw == alpha_kw + beta_kw + gamma_kw` exactly; `ex_tritium_kw` is
/// the total minus the tritium entry.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DecayHeat {
    /// `TOTAL ALPHA HEAT PRODUCTION`: ΣAi·Eα,i·C1, kW.
    pub alpha_kw: f64,
    /// `TOTAL BETA HEAT PRODUCTION`: ΣAi·Eβ,i·C1, kW.
    pub beta_kw: f64,
    /// `TOTAL GAMMA HEAT PRODUCTION`: ΣAi·Eγ,i·C1, kW.
    pub gamma_kw: f64,
    /// `TOTAL HEAT PRODUCTION`, kW.
    pub total_kw: f64,
    /// `TOTAL HEAT EX TRITIUM`, kW.
    pub ex_tritium_kw: f64,
}

/// S2 decay heat `ΣAi·E·C1` per radiation class (kW).
///
/// Every activity and decay energy must be finite and `>= 0`; violations are
/// loud [`Error`]s. The reported total is the class sum, so [`DecayHeat`]
/// parts always sum to the total.
pub fn decay_heat(entries: &[DecayHeatEntry]) -> Result<DecayHeat> {
    let tritium = tritium_id();
    let mut alpha = 0.0f64;
    let mut beta = 0.0f64;
    let mut gamma = 0.0f64;
    let mut tritium_kw = 0.0f64;
    for entry in entries {
        let name = entry.nuclide.to_name();
        if !entry.activity_bq.is_finite() || entry.activity_bq < 0.0 {
            return Err(Error::BadHeatValue {
                nuclide: name,
                msg: format!(
                    "activity must be finite and >= 0, got {}",
                    entry.activity_bq
                ),
            });
        }
        for (label, energy) in [
            ("e_alpha_ev", entry.e_alpha_ev),
            ("e_beta_ev", entry.e_beta_ev),
            ("e_gamma_ev", entry.e_gamma_ev),
        ] {
            if !energy.is_finite() || energy < 0.0 {
                return Err(Error::BadHeatValue {
                    nuclide: name.clone(),
                    msg: format!("{label} must be finite and >= 0, got {energy}"),
                });
            }
        }
        let heat_alpha = entry.activity_bq * entry.e_alpha_ev * EV_TO_KJ;
        let heat_beta = entry.activity_bq * entry.e_beta_ev * EV_TO_KJ;
        let heat_gamma = entry.activity_bq * entry.e_gamma_ev * EV_TO_KJ;
        alpha += heat_alpha;
        beta += heat_beta;
        gamma += heat_gamma;
        if entry.nuclide == tritium {
            tritium_kw += heat_alpha + heat_beta + heat_gamma;
        }
    }
    let total = alpha + beta + gamma;
    if !total.is_finite() {
        return Err(Error::BadHeatValue {
            nuclide: "total".to_string(),
            msg: format!(
                "heat sum overflowed to {total} over {} entries",
                entries.len()
            ),
        });
    }
    Ok(DecayHeat {
        alpha_kw: alpha,
        beta_kw: beta,
        gamma_kw: gamma,
        total_kw: total,
        ex_tritium_kw: total - tritium_kw,
    })
}

/// One inventory row for the S4/S5 committed-dose hazard sums.
///
/// One shared `HAZARDS` shape covers both kernels: each kernel folds its own
/// caller-supplied 50-year committed dose-coefficient column (S4 ingestion
/// `e^ing_i`, S5 inhalation `e^inh_i`, both Sv/Bq). Coefficients are
/// ICRP-copyrighted and so stay caller inputs, never vendored. There is no
/// split key: every entry contributes its full `Ai·e_i`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HazardEntry {
    /// Nuclide carrying this activity (tritium feeds the ex-tritium total).
    pub nuclide: NuclideId,
    /// Per-nuclide activity `Ai = Ni λi` in Bq (finite, `>= 0`).
    pub activity_bq: f64,
    /// Committed dose coefficient `e_i` in Sv/Bq (finite, `>= 0`).
    pub coeff_sv_per_bq: f64,
}

/// One S4/S5 committed-dose hazard total (Sv).
///
/// `total_sv` is the `TOTAL ... HAZARD FOR ALL MATERIALS` sum
/// (`ΣAi·e_i`, exactly); `ex_tritium_sv` is the total minus the tritium
/// entry (the S1/S2 precedent).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HazardDose {
    /// `TOTAL INGESTION/INHALATION HAZARD FOR ALL MATERIALS`: `ΣAi·e_i`, Sv.
    pub total_sv: f64,
    /// Total minus the tritium entry, Sv.
    pub ex_tritium_sv: f64,
}

/// Shared S4/S5 fold: `ΣAi·e_i` plus the ex-tritium companion.
///
/// `what` names the hazard for loud errors (`"ingestion"` / `"inhalation"`).
/// Every activity and coefficient must be finite and `>= 0`, and the summed
/// total must stay finite; violations are loud [`Error::BadHazardValue`]s,
/// never silent drops.
fn fold_hazard(entries: &[HazardEntry], what: &str) -> Result<HazardDose> {
    let tritium = tritium_id();
    let mut total = 0.0f64;
    let mut tritium_dose = 0.0f64;
    for entry in entries {
        let name = entry.nuclide.to_name();
        if !entry.activity_bq.is_finite() || entry.activity_bq < 0.0 {
            return Err(Error::BadHazardValue {
                nuclide: name,
                msg: format!(
                    "{what}: activity must be finite and >= 0, got {}",
                    entry.activity_bq
                ),
            });
        }
        if !entry.coeff_sv_per_bq.is_finite() || entry.coeff_sv_per_bq < 0.0 {
            return Err(Error::BadHazardValue {
                nuclide: name,
                msg: format!(
                    "{what}: coefficient must be finite and >= 0, got {}",
                    entry.coeff_sv_per_bq
                ),
            });
        }
        let dose = entry.activity_bq * entry.coeff_sv_per_bq;
        total += dose;
        if entry.nuclide == tritium {
            tritium_dose += dose;
        }
    }
    if !total.is_finite() {
        return Err(Error::BadHazardValue {
            nuclide: "total".to_string(),
            msg: format!(
                "{what}: hazard sum overflowed to {total} over {} entries",
                entries.len()
            ),
        });
    }
    Ok(HazardDose {
        total_sv: total,
        ex_tritium_sv: total - tritium_dose,
    })
}

/// S4 ingestion hazard `ΣAi·e^ing_i` (Sv, 50-year committed).
///
/// The Table X inventory row `Ai` folded with the caller-supplied ingestion
/// coefficients into the open `output_interpretation` wide-table `INGESTION
/// DOSE(Sv)` column (col 8 of the radiological block) and its `INGESTION
/// HAZARD FOR ALL MATERIALS` total (see the `HAZARDS` keyword page).
pub fn ingestion_hazard(entries: &[HazardEntry]) -> Result<HazardDose> {
    fold_hazard(entries, "ingestion")
}

/// S5 inhalation hazard `ΣAi·e^inh_i` (Sv, 50-year committed).
///
/// The Table X inventory row `Ai` folded with the caller-supplied inhalation
/// coefficients into the open `output_interpretation` wide-table `INHALATION
/// DOSE(Sv)` column (col 9 of the radiological block) and its `INHALATION
/// HAZARD FOR ALL MATERIALS` total (see the `HAZARDS` keyword page).
pub fn inhalation_hazard(entries: &[HazardEntry]) -> Result<HazardDose> {
    fold_hazard(entries, "inhalation")
}

/// Conversion factor C2 from TBq to Bq (S6 transport denominator).
///
/// Each per-nuclide `Bq/A2` term divides `Ai` (Bq) by `A2,i·C2` (Bq), so the
/// ratio is dimensionless and the effective A2 below converts back to TBq.
pub const TBQ_TO_BQ: f64 = 1e12;

/// One inventory row for the S6 transport-ratio sum.
///
/// There is no split key: every entry contributes its full `Ai/(A2,i·C2)`
/// term (the S4/S5 `HazardEntry` precedent).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransportEntry {
    /// Nuclide carrying this activity.
    pub nuclide: NuclideId,
    /// Per-nuclide activity `Ai = Ni λi` in Bq (finite, `>= 0`).
    pub activity_bq: f64,
    /// Transport limit `A2,i` in TBq (finite, `> 0`; zero is a loud
    /// division error, never `inf`).
    pub a2_tbq: f64,
}

/// S6 transport ratio with the effective A2.
///
/// `ratio` is the dimensionless `Total Bq/A2` sum; `total_bq` is `ΣAi` (Bq);
/// `effective_a2_tbq` is `hA2 = ΣAi/ratio` in TBq, defined by the reported
/// ratio — exactly `total_bq / ratio / C2` whenever `ratio != 0.0` (and
/// `0.0` for an empty inventory, where the quotient is undefined).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransportRatio {
    /// `Total Bq/A2`: `ΣAi/(A2,i·C2)`, dimensionless.
    pub ratio: f64,
    /// `ΣAi` over all entries, Bq.
    pub total_bq: f64,
    /// Effective A2 `hA2 = ΣAi/ratio`, TBq (`0.0` when `ratio == 0.0`).
    pub effective_a2_tbq: f64,
}

/// S6 transport ratio `ΣAi/(A2,i·C2)` with the effective A2 (TBq).
///
/// The Table X inventory row `Ai` folded over the caller-supplied transport
/// limits through Eqs 78/79 into the `ATWO`-block `Total Bq/A2` ratio; the
/// effective A2 is defined by that ratio. Every activity must be finite and
/// `>= 0`, every A2 limit finite and `> 0`, and both sums must stay finite;
/// violations are loud [`Error::BadTransportValue`]s, never silent drops.
/// Transport limits are caller inputs (never vendored); duplicate nuclides
/// sum like any other entries.
pub fn transport_ratio(entries: &[TransportEntry]) -> Result<TransportRatio> {
    let mut ratio = 0.0f64;
    let mut total = 0.0f64;
    for entry in entries {
        let name = entry.nuclide.to_name();
        if !entry.activity_bq.is_finite() || entry.activity_bq < 0.0 {
            return Err(Error::BadTransportValue {
                nuclide: name,
                msg: format!(
                    "activity must be finite and >= 0, got {}",
                    entry.activity_bq
                ),
            });
        }
        if !entry.a2_tbq.is_finite() || entry.a2_tbq <= 0.0 {
            return Err(Error::BadTransportValue {
                nuclide: name,
                msg: format!("A2 limit must be finite and > 0, got {}", entry.a2_tbq),
            });
        }
        let denom = entry.a2_tbq * TBQ_TO_BQ;
        if !denom.is_finite() {
            return Err(Error::BadTransportValue {
                nuclide: "total".to_string(),
                msg: format!(
                    "A2·C2 product overflowed to {denom} for {}",
                    entry.nuclide.to_name()
                ),
            });
        }
        ratio += entry.activity_bq / denom;
        total += entry.activity_bq;
    }
    if !ratio.is_finite() || !total.is_finite() {
        return Err(Error::BadTransportValue {
            nuclide: "total".to_string(),
            msg: format!(
                "transport sum overflowed to ratio {ratio} / total {total} over {} entries",
                entries.len()
            ),
        });
    }
    Ok(TransportRatio {
        ratio,
        total_bq: total,
        effective_a2_tbq: if ratio == 0.0 {
            0.0
        } else {
            total / ratio / TBQ_TO_BQ
        },
    })
}

/// One inventory row for the S7 IAEA clearance-index sum.
///
/// The same arithmetic as the EU-table [`crate::clearance`] variant — only
/// the limit table differs (IAEA values, permission-gated, caller-supplied).
/// There is no split key: every entry contributes its full `Ai/(Mtot·Li)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IaeaEntry {
    /// Nuclide carrying this activity.
    pub nuclide: NuclideId,
    /// Per-nuclide activity `Ai = Ni λi` in Bq (finite, `>= 0`).
    pub activity_bq: f64,
    /// IAEA clearance level `Li` in Bq/kg (finite, `> 0`; zero is a loud
    /// division error, never `inf`).
    pub limit_bq_per_kg: f64,
}

/// S7 IAEA clearance index with its screening class.
///
/// `index` is the dimensionless `ΣAi/(Mtot·Li)` sum; `class` screens it at
/// the landed `<= 1` boundary (boundary included); `max_fraction` /
/// `max_nuclide` carry the dominant per-nuclide term (`None` for an empty
/// inventory, where `max_fraction` is `0.0`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IaeaClearance {
    /// `CLEARANCE INDEX`: `ΣAi/(Mtot·Li)`, dimensionless.
    pub index: f64,
    /// Screening classification of `index` (the landed boundary rule).
    pub class: ClearanceClass,
    /// Largest per-nuclide term `Ai/(Mtot·Li)`.
    pub max_fraction: f64,
    /// Nuclide carrying [`Self::max_fraction`]; `None` for an empty inventory.
    pub max_nuclide: Option<NuclideId>,
}

/// S7 IAEA clearance index `ΣAi/(Mtot·Li)` (dimensionless).
///
/// The Table X inventory row `Ai` folded with the caller-supplied IAEA levels
/// through Eq. 80 into the `CLEAR`-block `CLEARANCE INDEX`. `total_mass_kg`
/// is `Mtot` (finite, `> 0`); every activity must be finite and `>= 0` and
/// every level finite and `> 0`, and the summed index must stay finite;
/// violations are loud [`Error::BadClearanceValue`]s, never silent drops.
/// IAEA levels are caller inputs (never vendored); duplicate nuclides sum
/// like any other entries. Screening arithmetic only, never a compliance
/// decision.
pub fn iaea_clearance_index(total_mass_kg: f64, entries: &[IaeaEntry]) -> Result<IaeaClearance> {
    if !total_mass_kg.is_finite() || total_mass_kg <= 0.0 {
        return Err(Error::BadClearanceValue {
            nuclide: "mass".to_string(),
            msg: format!("total mass must be finite and > 0, got {total_mass_kg}"),
        });
    }
    let mut index = 0.0f64;
    let mut max_fraction = 0.0f64;
    let mut max_nuclide = None;
    for entry in entries {
        let name = entry.nuclide.to_name();
        if !entry.activity_bq.is_finite() || entry.activity_bq < 0.0 {
            return Err(Error::BadClearanceValue {
                nuclide: name,
                msg: format!(
                    "activity must be finite and >= 0, got {}",
                    entry.activity_bq
                ),
            });
        }
        if !entry.limit_bq_per_kg.is_finite() || entry.limit_bq_per_kg <= 0.0 {
            return Err(Error::BadClearanceValue {
                nuclide: name,
                msg: format!(
                    "IAEA level must be finite and > 0, got {}",
                    entry.limit_bq_per_kg
                ),
            });
        }
        let denom = total_mass_kg * entry.limit_bq_per_kg;
        if !denom.is_finite() {
            return Err(Error::BadClearanceValue {
                nuclide: "total".to_string(),
                msg: format!(
                    "mass·level product overflowed to {denom} for {}",
                    entry.nuclide.to_name()
                ),
            });
        }
        let fraction = entry.activity_bq / denom;
        index += fraction;
        if fraction > max_fraction {
            max_fraction = fraction;
            max_nuclide = Some(entry.nuclide);
        }
    }
    if !index.is_finite() {
        return Err(Error::BadClearanceValue {
            nuclide: "total".to_string(),
            msg: format!(
                "clearance sum overflowed to {index} over {} entries",
                entries.len()
            ),
        });
    }
    let class = if index <= 1.0 {
        ClearanceClass::Satisfied
    } else {
        ClearanceClass::Exceeded
    };
    Ok(IaeaClearance {
        index,
        class,
        max_fraction,
        max_nuclide,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn nuc(name: &str) -> NuclideId {
        NuclideId::from_name(name).unwrap()
    }

    fn act(nuclide: &str, activity_bq: f64, irt: u8, alpha_frac: Option<f64>) -> ActivityEntry {
        ActivityEntry {
            nuclide: nuc(nuclide),
            activity_bq,
            irt,
            alpha_frac,
        }
    }

    #[test]
    fn ev_to_kj_pins_the_si_elementary_charge() {
        assert_eq!(EV_TO_KJ, 1.602_176_634e-22);
    }

    #[test]
    fn hand_computed_activity_vectors_at_exact_equality() {
        // Single beta emitter at a round activity.
        let out = total_activity(&[act("Co60", 10.0, 1, None)]).unwrap();
        assert_eq!(out.total_bq, 10.0);
        assert_eq!(out.beta_bq, 10.0);
        assert_eq!((out.alpha_bq, out.gamma_bq), (0.0, 0.0));
        // Mixed whole-class vector: 10 (β) + 5 (α) + 4 (γ) == 19 exactly.
        let out = total_activity(&[
            act("Co60", 10.0, 1, None),
            act("Po210", 5.0, 4, None),
            act("Tc99_m1", 4.0, 3, None),
        ])
        .unwrap();
        assert_eq!(out.total_bq, 19.0);
        assert_eq!(out.alpha_bq, 5.0);
        assert_eq!(out.beta_bq, 10.0);
        assert_eq!(out.gamma_bq, 4.0);
        // Empty inventory: all totals zero.
        let out = total_activity(&[]).unwrap();
        assert_eq!(out.total_bq, 0.0);
        assert_eq!(out.ex_tritium_bq, 0.0);
    }

    #[test]
    fn split_irts_conserve_parts_to_total() {
        // IRT 12 with α fraction 1/4: 25 α + 75 β from 100 exactly.
        let out = total_activity(&[act("U235", 100.0, 12, Some(0.25))]).unwrap();
        assert_eq!(out.alpha_bq, 25.0);
        assert_eq!(out.beta_bq, 75.0);
        assert_eq!(out.total_bq, 100.0);
        // IRT 13 full-beta end (frac 0) and full-alpha end (frac 1).
        let out = total_activity(&[act("U235", 8.0, 13, Some(0.0))]).unwrap();
        assert_eq!((out.alpha_bq, out.beta_bq), (0.0, 8.0));
        let out = total_activity(&[act("U235", 8.0, 13, Some(1.0))]).unwrap();
        assert_eq!((out.alpha_bq, out.beta_bq), (8.0, 0.0));
        // IRT 15 half split: 5 α + 5 γ from 10 exactly.
        let out = total_activity(&[act("Pu240", 10.0, 15, Some(0.5))]).unwrap();
        assert_eq!(out.alpha_bq, 5.0);
        assert_eq!(out.gamma_bq, 5.0);
        assert_eq!(out.total_bq, 10.0);
        // Mixed whole + split vector conserves exactly.
        let out = total_activity(&[
            act("Co60", 10.0, 1, None),
            act("Po210", 5.0, 4, None),
            act("U235", 100.0, 12, Some(0.25)),
            act("Pu240", 10.0, 15, Some(0.5)),
        ])
        .unwrap();
        assert_eq!(out.alpha_bq, 5.0 + 25.0 + 5.0);
        assert_eq!(out.beta_bq, 10.0 + 75.0);
        assert_eq!(out.gamma_bq, 5.0);
        assert_eq!(out.total_bq, out.alpha_bq + out.beta_bq + out.gamma_bq);
    }

    #[test]
    fn ex_tritium_excludes_h3_only() {
        let out = total_activity(&[act("Co60", 10.0, 1, None), act("H3", 5.0, 1, None)]).unwrap();
        assert_eq!(out.total_bq, 15.0);
        assert_eq!(out.beta_bq, 15.0);
        assert_eq!(out.ex_tritium_bq, 10.0);
        // No tritium: ex-tritium equals the total exactly.
        let out = total_activity(&[act("Co60", 10.0, 1, None)]).unwrap();
        assert_eq!(out.ex_tritium_bq, out.total_bq);
    }

    #[test]
    fn all_beta_irts_map_to_beta() {
        // Every pinned whole-beta IRT routes to the beta total.
        let entries: Vec<ActivityEntry> = [1, 2, 11, 14, 16, 17, 19, 20]
            .iter()
            .map(|irt| act("Co60", 1.0, *irt, None))
            .collect();
        let out = total_activity(&entries).unwrap();
        assert_eq!(out.beta_bq, 8.0);
        assert_eq!(out.total_bq, 8.0);
    }

    #[test]
    fn malformed_activity_inputs_are_loud() {
        // Negative / non-finite activities.
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(
                matches!(
                    total_activity(&[act("Co60", bad, 1, None)]),
                    Err(Error::BadActivityValue { .. })
                ),
                "activity {bad} should fail"
            );
        }
        // Unmapped IRTs: unused (8, 9), unknown (10), unlisted (0, 5, 27, 28).
        for irt in [0, 5, 6, 7, 8, 9, 10, 18, 21, 27, 28, 255] {
            match total_activity(&[act("Co60", 1.0, irt, None)]) {
                Err(Error::UnmappedIrt { irt: got, .. }) => assert_eq!(got, irt),
                other => panic!("IRT {irt} should fail loud, got {other:?}"),
            }
        }
        // Split IRT without a fraction; whole IRT with one.
        assert!(matches!(
            total_activity(&[act("U235", 1.0, 12, None)]),
            Err(Error::BadIrtSplit { .. })
        ));
        assert!(matches!(
            total_activity(&[act("Pu240", 1.0, 15, None)]),
            Err(Error::BadIrtSplit { .. })
        ));
        assert!(matches!(
            total_activity(&[act("Co60", 1.0, 1, Some(0.5))]),
            Err(Error::BadIrtSplit { .. })
        ));
        // Out-of-range / non-finite fractions.
        for bad in [-0.5, 1.5, f64::NAN, f64::INFINITY] {
            assert!(
                matches!(
                    total_activity(&[act("U235", 1.0, 12, Some(bad))]),
                    Err(Error::BadIrtSplit { .. })
                ),
                "alpha_frac {bad} should fail"
            );
        }
    }

    #[test]
    fn overflowing_activity_sum_is_a_loud_error() {
        let out = total_activity(&[act("Co60", f64::MAX, 1, None), act("H3", f64::MAX, 1, None)]);
        match out {
            Err(Error::BadActivityValue { nuclide, msg }) => {
                assert_eq!(nuclide, "total");
                assert!(msg.contains("overflowed"), "msg was `{msg}`");
            }
            other => panic!("expected overflow error, got {other:?}"),
        }
    }

    fn heat(
        nuclide: &str,
        activity_bq: f64,
        e_alpha: f64,
        e_beta: f64,
        e_gamma: f64,
    ) -> DecayHeatEntry {
        DecayHeatEntry {
            nuclide: nuc(nuclide),
            activity_bq,
            e_alpha_ev: e_alpha,
            e_beta_ev: e_beta,
            e_gamma_ev: e_gamma,
        }
    }

    #[test]
    fn hand_computed_heat_vectors_at_exact_equality() {
        // 1e10 Bq × 1e6 eV × C1 == 1.602176634e-6 kW per class.
        let unit = 1e10 * 1e6 * EV_TO_KJ;
        let out = decay_heat(&[heat("Co60", 1e10, 0.0, 1e6, 1e6)]).unwrap();
        assert_eq!(out.beta_kw, unit);
        assert_eq!(out.gamma_kw, unit);
        assert_eq!(out.alpha_kw, 0.0);
        assert_eq!(out.total_kw, unit + unit);
        assert_eq!(out.total_kw, out.alpha_kw + out.beta_kw + out.gamma_kw);
        // Mixed vector: alpha emitter + beta/gamma emitter.
        let out = decay_heat(&[
            heat("Po210", 2e10, 5e6, 0.0, 0.0),
            heat("Co60", 1e10, 0.0, 1e6, 2e6),
        ])
        .unwrap();
        assert_eq!(out.alpha_kw, 2e10 * 5e6 * EV_TO_KJ);
        assert_eq!(out.beta_kw, 1e10 * 1e6 * EV_TO_KJ);
        assert_eq!(out.gamma_kw, 1e10 * 2e6 * EV_TO_KJ);
        assert_eq!(out.total_kw, out.alpha_kw + out.beta_kw + out.gamma_kw);
        // Empty inventory: all heats zero.
        let out = decay_heat(&[]).unwrap();
        assert_eq!(out.total_kw, 0.0);
        assert_eq!(out.ex_tritium_kw, 0.0);
    }

    #[test]
    fn ex_tritium_heat_excludes_h3_only() {
        // H3: 5 Bq × 6e3 eV β; Co60: 10 Bq × 1e6 eV β.
        let h3 = 5.0 * 6e3 * EV_TO_KJ;
        let co = 10.0 * 1e6 * EV_TO_KJ;
        let out = decay_heat(&[
            heat("Co60", 10.0, 0.0, 1e6, 0.0),
            heat("H3", 5.0, 0.0, 6e3, 0.0),
        ])
        .unwrap();
        assert_eq!(out.total_kw, co + h3);
        assert_eq!(out.ex_tritium_kw, co);
        assert_eq!(out.ex_tritium_kw, out.total_kw - h3);
    }

    #[test]
    fn malformed_heat_inputs_are_loud() {
        // Negative / non-finite activities and energies.
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(
                matches!(
                    decay_heat(&[heat("Co60", bad, 0.0, 1e6, 0.0)]),
                    Err(Error::BadHeatValue { .. })
                ),
                "activity {bad} should fail"
            );
            for entry in [
                heat("Co60", 1.0, bad, 0.0, 0.0),
                heat("Co60", 1.0, 0.0, bad, 0.0),
                heat("Co60", 1.0, 0.0, 0.0, bad),
            ] {
                assert!(
                    matches!(decay_heat(&[entry]), Err(Error::BadHeatValue { .. })),
                    "energy {bad} should fail"
                );
            }
        }
        // Zero energies are valid (no emission in that class).
        let out = decay_heat(&[heat("Fe55", 3.0, 0.0, 0.0, 0.0)]).unwrap();
        assert_eq!(out.total_kw, 0.0);
    }

    fn hazard(nuclide: &str, activity_bq: f64, coeff: f64) -> HazardEntry {
        HazardEntry {
            nuclide: nuc(nuclide),
            activity_bq,
            coeff_sv_per_bq: coeff,
        }
    }

    #[test]
    fn hand_computed_hazard_vectors_at_exact_equality() {
        // S4: 10 Bq × 3 Sv/Bq + 5 Bq × 2 Sv/Bq == 40 Sv exactly.
        let out = ingestion_hazard(&[hazard("Co60", 10.0, 3.0), hazard("H3", 5.0, 2.0)]).unwrap();
        assert_eq!(out.total_sv, 40.0);
        assert_eq!(out.ex_tritium_sv, 30.0);
        assert_eq!(out.ex_tritium_sv, out.total_sv - 5.0 * 2.0);
        // S5 over the same activities with its own column: 10×5 + 5×7 == 85.
        let out = inhalation_hazard(&[hazard("Co60", 10.0, 5.0), hazard("H3", 5.0, 7.0)]).unwrap();
        assert_eq!(out.total_sv, 85.0);
        assert_eq!(out.ex_tritium_sv, 50.0);
        // Single-entry and empty inventories.
        let out = ingestion_hazard(&[hazard("Co60", 10.0, 2.0)]).unwrap();
        assert_eq!(out.total_sv, 20.0);
        assert_eq!(out.ex_tritium_sv, out.total_sv);
        for kernel in [ingestion_hazard, inhalation_hazard] {
            let out = kernel(&[]).unwrap();
            assert_eq!(out.total_sv, 0.0);
            assert_eq!(out.ex_tritium_sv, 0.0);
        }
    }

    #[test]
    fn ex_tritium_hazard_excludes_h3_only() {
        // No tritium: ex-tritium equals the total exactly, for both kernels.
        for kernel in [ingestion_hazard, inhalation_hazard] {
            let out = kernel(&[hazard("Co60", 10.0, 3.0)]).unwrap();
            assert_eq!(out.ex_tritium_sv, out.total_sv);
        }
        // Zero coefficient contributes nothing but stays valid.
        let out = ingestion_hazard(&[hazard("Fe55", 3.0, 0.0)]).unwrap();
        assert_eq!(out.total_sv, 0.0);
    }

    #[test]
    fn malformed_hazard_inputs_are_loud() {
        // Negative / non-finite activities and coefficients, for both kernels.
        for kernel in [ingestion_hazard, inhalation_hazard] {
            for bad in [-1.0, f64::NAN, f64::INFINITY] {
                assert!(
                    matches!(
                        kernel(&[hazard("Co60", bad, 1.0)]),
                        Err(Error::BadHazardValue { .. })
                    ),
                    "activity {bad} should fail"
                );
                assert!(
                    matches!(
                        kernel(&[hazard("Co60", 1.0, bad)]),
                        Err(Error::BadHazardValue { .. })
                    ),
                    "coefficient {bad} should fail"
                );
            }
        }
    }

    #[test]
    fn overflowing_hazard_sum_is_a_loud_error() {
        let out = ingestion_hazard(&[hazard("Co60", f64::MAX, f64::MAX)]);
        match out {
            Err(Error::BadHazardValue { nuclide, msg }) => {
                assert_eq!(nuclide, "total");
                assert!(msg.contains("overflowed"), "msg was `{msg}`");
            }
            other => panic!("expected overflow error, got {other:?}"),
        }
    }

    fn transport(nuclide: &str, activity_bq: f64, a2_tbq: f64) -> TransportEntry {
        TransportEntry {
            nuclide: nuc(nuclide),
            activity_bq,
            a2_tbq,
        }
    }

    #[test]
    fn tbq_to_bq_pins_the_transport_conversion() {
        assert_eq!(TBQ_TO_BQ, 1e12);
    }

    #[test]
    fn hand_computed_transport_vectors_at_exact_equality() {
        // 1 TBq of activity against a 1 TBq limit: term == 1 exactly.
        let out = transport_ratio(&[transport("Co60", 1e12, 1.0)]).unwrap();
        assert_eq!(out.ratio, 1.0);
        assert_eq!(out.total_bq, 1e12);
        assert_eq!(out.effective_a2_tbq, 1.0);
        // Mixed vector: 1e12/1 + 2e12/1 == 3 exactly; effective == 1 TBq.
        let out =
            transport_ratio(&[transport("Co60", 1e12, 1.0), transport("H3", 2e12, 1.0)]).unwrap();
        assert_eq!(out.ratio, 3.0);
        assert_eq!(out.total_bq, 3e12);
        assert_eq!(out.effective_a2_tbq, 1.0);
        // Sub-unit term: 5e11 Bq against 2 TBq == 0.25 exactly.
        let out = transport_ratio(&[transport("Fe55", 5e11, 2.0)]).unwrap();
        assert_eq!(out.ratio, 0.25);
        assert_eq!(out.effective_a2_tbq, 2.0);
        // Effective is defined by the ratio: total / ratio / C2, exactly.
        let out = transport_ratio(&[transport("Co60", 3e12, 2.0), transport("Cs137", 1e12, 4.0)])
            .unwrap();
        assert_eq!(out.ratio, 1.5 + 0.25);
        assert_eq!(out.effective_a2_tbq, out.total_bq / out.ratio / TBQ_TO_BQ);
        // Empty inventory: ratio, total, and effective all zero.
        let out = transport_ratio(&[]).unwrap();
        assert_eq!(out.ratio, 0.0);
        assert_eq!(out.total_bq, 0.0);
        assert_eq!(out.effective_a2_tbq, 0.0);
        // Zero activity contributes nothing but stays valid.
        let out = transport_ratio(&[transport("Fe55", 0.0, 1.0)]).unwrap();
        assert_eq!(out.ratio, 0.0);
        assert_eq!(out.effective_a2_tbq, 0.0);
    }

    #[test]
    fn malformed_transport_inputs_are_loud() {
        // Negative / non-finite activities.
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(
                matches!(
                    transport_ratio(&[transport("Co60", bad, 1.0)]),
                    Err(Error::BadTransportValue { .. })
                ),
                "activity {bad} should fail"
            );
        }
        // Zero, negative, and non-finite A2 limits (zero divides loudly).
        for bad in [0.0, -2.0, f64::NAN, f64::INFINITY] {
            assert!(
                matches!(
                    transport_ratio(&[transport("Co60", 1.0, bad)]),
                    Err(Error::BadTransportValue { .. })
                ),
                "A2 limit {bad} should fail"
            );
        }
    }

    #[test]
    fn overflowing_transport_sum_is_a_loud_error() {
        let out = transport_ratio(&[transport("Co60", f64::MAX, 1e-308)]);
        match out {
            Err(Error::BadTransportValue { nuclide, msg }) => {
                assert_eq!(nuclide, "total");
                assert!(msg.contains("overflowed"), "msg was `{msg}`");
            }
            other => panic!("expected overflow error, got {other:?}"),
        }
        // An absurd A2 overflows the A2·C2 product itself (the term would
        // silently vanish to 0 against an infinite denominator).
        match transport_ratio(&[transport("Co60", 1.0, f64::MAX)]) {
            Err(Error::BadTransportValue { nuclide, msg }) => {
                assert_eq!(nuclide, "total");
                assert!(msg.contains("overflowed"), "msg was `{msg}`");
            }
            other => panic!("expected overflow error, got {other:?}"),
        }
    }

    fn iaea(nuclide: &str, activity_bq: f64, limit_bq_per_kg: f64) -> IaeaEntry {
        IaeaEntry {
            nuclide: nuc(nuclide),
            activity_bq,
            limit_bq_per_kg,
        }
    }

    #[test]
    fn hand_computed_iaea_vectors_at_exact_equality() {
        // Mass 2 kg: 10 Bq at 10 Bq/kg -> 10/(2·10) == 0.5 exactly.
        let out = iaea_clearance_index(2.0, &[iaea("Co60", 10.0, 10.0)]).unwrap();
        assert_eq!(out.index, 0.5);
        assert_eq!(out.class, ClearanceClass::Satisfied);
        assert_eq!(out.max_fraction, 0.5);
        assert_eq!(out.max_nuclide, Some(nuc("Co60")));
        // Two-entry vector at the boundary: 0.5 + 0.5 == 1 exactly.
        let out =
            iaea_clearance_index(2.0, &[iaea("Co60", 10.0, 10.0), iaea("H3", 5.0, 5.0)]).unwrap();
        assert_eq!(out.index, 1.0);
        assert_eq!(out.class, ClearanceClass::Satisfied);
        assert_eq!(out.max_nuclide, Some(nuc("Co60")));
        // Empty inventory: index 0, satisfied, no dominant nuclide.
        let out = iaea_clearance_index(2.0, &[]).unwrap();
        assert_eq!(out.index, 0.0);
        assert_eq!(out.class, ClearanceClass::Satisfied);
        assert_eq!(out.max_fraction, 0.0);
        assert_eq!(out.max_nuclide, None);
        // Zero activity contributes nothing but stays valid.
        let out = iaea_clearance_index(2.0, &[iaea("Fe55", 0.0, 1.0)]).unwrap();
        assert_eq!(out.index, 0.0);
        assert_eq!(out.class, ClearanceClass::Satisfied);
    }

    #[test]
    fn iaea_boundary_probes_both_sides() {
        // Exactly == 1 on the boundary: satisfied (boundary included).
        let out = iaea_clearance_index(2.0, &[iaea("Co60", 20.0, 10.0)]).unwrap();
        assert_eq!(out.index, 1.0);
        assert_eq!(out.class, ClearanceClass::Satisfied);
        // One ulp above: exceeded; dominant contributor tracked.
        let out = iaea_clearance_index(2.0, &[iaea("Co60", 20.000000000000004, 10.0)]).unwrap();
        assert!(out.index > 1.0);
        assert_eq!(out.class, ClearanceClass::Exceeded);
        assert_eq!(out.max_nuclide, Some(nuc("Co60")));
        // One ulp below: satisfied.
        let out = iaea_clearance_index(2.0, &[iaea("Co60", 19.999999999999996, 10.0)]).unwrap();
        assert!(out.index < 1.0);
        assert_eq!(out.class, ClearanceClass::Satisfied);
        // Accumulation both sides: 0.4 + 0.5 = 0.9 vs 0.6 + 0.5 = 1.1.
        let below = [iaea("Co60", 8.0, 10.0), iaea("H3", 5.0, 5.0)];
        assert_eq!(
            iaea_clearance_index(2.0, &below).unwrap().class,
            ClearanceClass::Satisfied
        );
        let above = [iaea("Co60", 12.0, 10.0), iaea("H3", 5.0, 5.0)];
        let out = iaea_clearance_index(2.0, &above).unwrap();
        assert!(out.index > 1.0);
        assert_eq!(out.class, ClearanceClass::Exceeded);
        assert_eq!(out.max_nuclide, Some(nuc("Co60")));
    }

    #[test]
    fn malformed_iaea_inputs_are_loud() {
        // Non-positive / non-finite total mass.
        for bad in [0.0, -2.0, f64::NAN, f64::INFINITY] {
            match iaea_clearance_index(bad, &[iaea("Co60", 1.0, 10.0)]) {
                Err(Error::BadClearanceValue { nuclide, .. }) => assert_eq!(nuclide, "mass"),
                other => panic!("mass {bad} should fail loud, got {other:?}"),
            }
        }
        // Negative / non-finite activities.
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(
                matches!(
                    iaea_clearance_index(2.0, &[iaea("Co60", bad, 10.0)]),
                    Err(Error::BadClearanceValue { .. })
                ),
                "activity {bad} should fail"
            );
        }
        // Zero, negative, and non-finite IAEA levels (zero divides loudly).
        for bad in [0.0, -2.0, f64::NAN, f64::INFINITY] {
            assert!(
                matches!(
                    iaea_clearance_index(2.0, &[iaea("Co60", 1.0, bad)]),
                    Err(Error::BadClearanceValue { .. })
                ),
                "IAEA level {bad} should fail"
            );
        }
    }

    #[test]
    fn overflowing_iaea_sum_is_a_loud_error() {
        let out = iaea_clearance_index(1e-308, &[iaea("Co60", f64::MAX, 1.0)]);
        match out {
            Err(Error::BadClearanceValue { nuclide, msg }) => {
                assert_eq!(nuclide, "total");
                assert!(msg.contains("overflowed"), "msg was `{msg}`");
            }
            other => panic!("expected overflow error, got {other:?}"),
        }
        // An absurd mass·level product overflows itself (the term would
        // silently vanish to 0 against an infinite denominator).
        let out = iaea_clearance_index(f64::MAX, &[iaea("Co60", 1.0, f64::MAX)]);
        match out {
            Err(Error::BadClearanceValue { nuclide, msg }) => {
                assert_eq!(nuclide, "total");
                assert!(msg.contains("overflowed"), "msg was `{msg}`");
            }
            other => panic!("expected overflow error, got {other:?}"),
        }
    }
}
