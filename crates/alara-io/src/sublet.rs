//! Sublet S1+S2 radiological totals: total activity and decay heat.
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
//! Both totals carry the `TOTAL ... EXCLUDING TRITIUM` companion (total minus
//! the tritium entry). Per-nuclide decay energies are caller-supplied, never
//! vendored. Parts always sum to the total by construction
//! (`total == alpha + beta + gamma`, exactly in `f64`).
//!
//! Explicitly out of scope: S3–S7 (later work, one per family), ICRP/IAEA
//! coefficient tables (caller inputs, never vendored), bremsstrahlung
//! correction, DPA/KERMA/gas (`damage` owns), fission count / burn-up.

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
}
