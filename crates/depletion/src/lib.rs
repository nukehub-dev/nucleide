//! Depletion: CRAM solver, chain files, results.
//!
//! Key decisions:
//! - chain files use the depletion-chain XML format (see [`chain`])
//! - backend is `faer` sparse LU with symbolic reuse, behind the `linalg` facade

pub mod chain;
pub mod cram;
pub mod integrate;
pub mod inventory;
pub mod matrix;

pub use chain::{Chain, ChainNuclide, DecayMode, Error, FissionYields, Reaction};
pub use cram::{cram, cram_with_symbolic, Order};
pub use integrate::{
    activity_vec, decay_constants, decay_energies_by_name, decay_energies_mev, decay_energy_mev,
    decay_heat_vec, integrate, Integrator, Step, TimeSeries, EV_PER_MEV, MEV_TO_JOULE,
};
pub use inventory::{
    branching_fraction, chain_edges, cumulative_decays, decay_mode, progeny, time_unit_from_str,
    DecayInventory, QuantityUnit, TimeUnit, AVOGADRO, CURIE_TO_BQ, SECONDS_PER_DAY,
    SECONDS_PER_YEAR,
};
pub use matrix::{DepletionSystem, ReactionRates};

use std::collections::BTreeMap;

/// Nuclide-vector result keyed by name, after a depletion step.
#[derive(Debug, Clone, PartialEq)]
pub struct DepletionResult {
    /// Final atom counts by nuclide name (all chain nuclides present).
    pub atoms: BTreeMap<String, f64>,
}

/// Convenience driver: build the system, solve one step, key results by name.
///
/// Returns [`DepletionResult`] so callers can diff/serialize without
/// touching indices.
pub fn deplete(
    sys: &DepletionSystem,
    order: Order,
    n0: &BTreeMap<String, f64>,
    dt: f64,
) -> Result<DepletionResult, Error> {
    let n = sys.chain.len();
    let mut n0_vec = vec![0.0; n];
    for (name, value) in n0 {
        let idx = sys
            .chain
            .index_of(name)
            .ok_or_else(|| Error::UnknownNuclide {
                name: name.clone(),
                context: "n0",
            })?;
        n0_vec[idx] = *value;
    }
    let final_vec =
        cram(sys, order, &n0_vec, dt).map_err(|e| Error::BadStructure(e.to_string()))?;
    let atoms: BTreeMap<String, f64> = sys
        .chain
        .nuclides
        .iter()
        .zip(final_vec)
        .map(|(nuc, v)| (nuc.name.clone(), v))
        .collect();
    Ok(DepletionResult { atoms })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chain::{Chain, ChainNuclide, DecayMode, Reaction};
    use crate::matrix::{DepletionSystem, ReactionRates};
    use std::collections::BTreeMap;

    fn simple_chain() -> Chain {
        // A -lam-> B -lam2-> C (stable), plus capture A->B
        let a = ChainNuclide {
            name: "A".into(),
            half_life: Some(std::f64::consts::LN_2 / 1e-6), // lambda = 1e-6
            decay_modes: vec![DecayMode {
                kind: "beta".into(),
                target: "B".into(),
                branching_ratio: 1.0,
            }],
            ..Default::default()
        };
        let b = ChainNuclide {
            name: "B".into(),
            half_life: Some(std::f64::consts::LN_2 / 1e-5), // lambda = 1e-5
            decay_modes: vec![DecayMode {
                kind: "beta".into(),
                target: "C".into(),
                branching_ratio: 1.0,
            }],
            reactions: vec![Reaction {
                kind: "(n,gamma)".into(),
                target: Some("C".into()),
                q: 0.0,
                branching_ratio: 1.0,
            }],
            ..Default::default()
        };
        let c = ChainNuclide {
            name: "C".into(),
            ..Default::default()
        };
        Chain::from_nuclides(vec![a, b, c]).unwrap()
    }

    #[test]
    fn analytic_two_step_decay_cram48() {
        let chain = Chain::from_xml(&xml_of(simple_chain())).unwrap();
        let sys = DepletionSystem::build(chain, &ReactionRates::new()).unwrap();
        let n0 = vec![1.0e15, 0.0, 0.0];
        let dt = 1.0e5;

        let n_final = crate::cram(&sys, Order::Order48, &n0, dt).unwrap();

        // Analytic solution for pure sequential decay:
        let l1 = 1.0e-6;
        let l2 = 1.0e-5;
        let na = n0[0] * (-l1 * dt).exp();
        let nb = n0[0] * l1 / (l2 - l1) * ((-l1 * dt).exp() - (-l2 * dt).exp());
        let nc = n0[0] - na - nb;

        assert!(
            (n_final[0] - na).abs() / na < 1e-8,
            "{} vs {}",
            n_final[0],
            na
        );
        assert!(
            (n_final[1] - nb).abs() / nb < 1e-7,
            "{} vs {}",
            n_final[1],
            nb
        );
        assert!(
            (n_final[2] - nc).abs() / nc < 1e-9,
            "{} vs {}",
            n_final[2],
            nc
        );
        // Conservation of atoms
        let total: f64 = n_final.iter().sum();
        assert!((total - n0[0]).abs() / n0[0] < 1e-8);
    }

    #[test]
    fn cram16_agrees_within_tolerance() {
        let chain = Chain::from_xml(&xml_of(simple_chain())).unwrap();
        let sys = DepletionSystem::build(chain, &ReactionRates::new()).unwrap();
        let n0 = vec![1.0e12, 0.0, 0.0];
        let dt = 5.0e4;
        let r48 = crate::cram(&sys, Order::Order48, &n0, dt).unwrap();
        let r16 = crate::cram(&sys, Order::Order16, &n0, dt).unwrap();
        for (a, b) in r48.iter().zip(&r16) {
            assert!((a - b).abs() / a.max(1e-30) < 1e-5, "{a} vs {b}");
        }
    }

    #[test]
    fn reaction_channel_and_capture_loss() {
        let mut rates = ReactionRates::new();
        // sigma*phi for B capture = 1e-5 (comparable to its decay rate)
        rates
            .entry(1usize)
            .or_default()
            .insert("(n,gamma)".to_string(), 1e-5);
        let chain = Chain::from_xml(&xml_of(simple_chain())).unwrap();
        let sys = DepletionSystem::build(chain, &rates).unwrap();
        let n0 = vec![0.0, 1.0e15, 0.0];
        let dt = 1.0e5;
        let n_final = crate::cram(&sys, Order::Order48, &n0, dt).unwrap();

        // Effective loss rate for B: lam2 + capture
        let l2p = 1e-5 + 1e-5;
        let nb = n0[1] * (-l2p * dt).exp();
        let nc = n0[1] - nb;
        assert!((n_final[1] - nb).abs() / nb < 1e-7);
        assert!((n_final[2] - nc).abs() / nc < 1e-7);
        assert!(n_final[0].abs() < 1e-3);
    }

    #[test]
    fn symbolic_reuse_matches_fresh() {
        let chain = Chain::from_xml(&xml_of(simple_chain())).unwrap();
        let sys = DepletionSystem::build(chain, &ReactionRates::new()).unwrap();
        let sym = nucleide_linalg::SymbolicLu::try_new(&sys.pattern).unwrap();
        let n0 = vec![1.0e14, 5e13, 1e10];
        let a = crate::cram(&sys, Order::Order48, &n0, 2.0e5).unwrap();
        let b = crate::cram_with_symbolic(&sys, &sym, Order::Order48, &n0, 2.0e5).unwrap();
        for (x, y) in a.iter().zip(&b) {
            assert_eq!(x, y); // bit-identical by construction
        }
    }

    #[test]
    fn vendored_chain_parses() {
        let path = format!(
            "{}/../../fixtures/depletion/chain_simple.xml",
            env!("CARGO_MANIFEST_DIR")
        );
        let chain = Chain::from_file(path).unwrap();
        assert_eq!(chain.len(), 9);
        let i135 = chain.index_of("I135").unwrap();
        assert_eq!(chain.nuclides[i135].decay_modes.len(), 1);
        assert_eq!(chain.nuclides[i135].decay_modes[0].target, "Xe135");
        // Gd157 capture goes to "Nothing" => target None
        let gd = chain.index_of("Gd157").unwrap();
        assert_eq!(chain.nuclides[gd].reactions[0].target, None);
        // U235 fission yields present
        let u235 = chain.index_of("U235").unwrap();
        assert!(!chain.nuclides[u235].neutron_fission_yields.is_empty());
        let fy = &chain.nuclides[u235].neutron_fission_yields[0];
        assert_eq!(fy.products["I135"], 0.0292737);
    }

    #[test]
    fn deplete_driver_keys_by_name() {
        let path = format!(
            "{}/../../fixtures/depletion/chain_simple.xml",
            env!("CARGO_MANIFEST_DIR")
        );
        let chain = Chain::from_file(path).unwrap();
        let sys = DepletionSystem::build(chain, &ReactionRates::new()).unwrap();
        let mut n0 = BTreeMap::new();
        n0.insert("I135".to_string(), 1e15);
        let res = crate::deplete(&sys, Order::Order48, &n0, 1e5).unwrap();
        // I135 decays to Xe135 with half life 23652 s -> mostly gone after 1e5 s
        assert!(res.atoms["I135"] < 5e14);
        assert!(res.atoms["Xe135"] > 1e13);
    }

    #[test]
    fn reaction_branching_ratio_scales_gain() {
        // Parent A captures with branching ratio 0.3 to B; loss uses full rate.
        let a = ChainNuclide {
            name: "A".into(),
            reactions: vec![Reaction {
                kind: "(n,gamma)".into(),
                target: Some("B".into()),
                q: 0.0,
                branching_ratio: 0.3,
            }],
            ..Default::default()
        };
        let b = ChainNuclide {
            name: "B".into(),
            ..Default::default()
        };
        let chain = Chain::from_nuclides(vec![a, b]).unwrap();
        let mut rates = ReactionRates::new();
        rates
            .entry(0usize)
            .or_default()
            .insert("(n,gamma)".to_string(), 1e-5);

        let sys = DepletionSystem::build(chain, &rates).unwrap();
        let dense = sys.matrix_for_dt(1.0).unwrap().to_dense();
        assert!((dense[0][0].re + 1e-5).abs() < 1e-18);
        assert!((dense[1][0].re - 3e-6).abs() < 1e-18);
    }

    #[test]
    fn duplicate_reaction_entries_share_loss() {
        // Two (n,gamma) entries: 0.4 to B and 0.6 to C. Loss is subtracted once.
        let a = ChainNuclide {
            name: "A".into(),
            reactions: vec![
                Reaction {
                    kind: "(n,gamma)".into(),
                    target: Some("B".into()),
                    q: 0.0,
                    branching_ratio: 0.4,
                },
                Reaction {
                    kind: "(n,gamma)".into(),
                    target: Some("C".into()),
                    q: 0.0,
                    branching_ratio: 0.6,
                },
            ],
            ..Default::default()
        };
        let b = ChainNuclide {
            name: "B".into(),
            ..Default::default()
        };
        let c = ChainNuclide {
            name: "C".into(),
            ..Default::default()
        };
        let chain = Chain::from_nuclides(vec![a, b, c]).unwrap();
        let mut rates = ReactionRates::new();
        rates
            .entry(0usize)
            .or_default()
            .insert("(n,gamma)".to_string(), 1e-5);

        let sys = DepletionSystem::build(chain, &rates).unwrap();
        let dense = sys.matrix_for_dt(1.0).unwrap().to_dense();
        assert!((dense[0][0].re + 1e-5).abs() < 1e-18);
        assert!((dense[1][0].re - 4e-6).abs() < 1e-18);
        assert!((dense[2][0].re - 6e-6).abs() < 1e-18);
    }

    #[test]
    fn alpha_decay_produces_he4_and_proton_decay_h1() {
        let a = ChainNuclide {
            name: "A".into(),
            half_life: Some(std::f64::consts::LN_2 / 1e-6),
            decay_modes: vec![
                DecayMode {
                    kind: "alpha".into(),
                    target: "B".into(),
                    branching_ratio: 0.5,
                },
                DecayMode {
                    kind: "p".into(),
                    target: "C".into(),
                    branching_ratio: 0.5,
                },
            ],
            ..Default::default()
        };
        let b = ChainNuclide {
            name: "B".into(),
            ..Default::default()
        };
        let c = ChainNuclide {
            name: "C".into(),
            ..Default::default()
        };
        let he4 = ChainNuclide {
            name: "He4".into(),
            ..Default::default()
        };
        let h1 = ChainNuclide {
            name: "H1".into(),
            ..Default::default()
        };
        let chain = Chain::from_nuclides(vec![a, b, c, he4, h1]).unwrap();
        let sys = DepletionSystem::build(chain, &ReactionRates::new()).unwrap();
        let dense = sys.matrix_for_dt(1.0).unwrap().to_dense();
        let i_a = 0usize;
        let i_b = 1;
        let i_c = 2;
        let i_he4 = 3;
        let i_h1 = 4;
        let lam = 1e-6;
        // Daughter gains
        assert!((dense[i_b][i_a].re - 0.5 * lam).abs() < 1e-18);
        assert!((dense[i_c][i_a].re - 0.5 * lam).abs() < 1e-18);
        // Light-particle gains
        assert!((dense[i_he4][i_a].re - 0.5 * lam).abs() < 1e-18);
        assert!((dense[i_h1][i_a].re - 0.5 * lam).abs() < 1e-18);
    }

    #[test]
    fn reaction_secondaries_produce_light_nuclides() {
        let a = ChainNuclide {
            name: "A".into(),
            reactions: vec![Reaction {
                kind: "(n,a)".into(),
                target: Some("B".into()),
                q: 0.0,
                branching_ratio: 1.0,
            }],
            ..Default::default()
        };
        let b = ChainNuclide {
            name: "B".into(),
            ..Default::default()
        };
        let he4 = ChainNuclide {
            name: "He4".into(),
            ..Default::default()
        };
        let chain = Chain::from_nuclides(vec![a, b, he4]).unwrap();
        let mut rates = ReactionRates::new();
        rates
            .entry(0usize)
            .or_default()
            .insert("(n,a)".to_string(), 2e-5);

        let sys = DepletionSystem::build(chain, &rates).unwrap();
        let dense = sys.matrix_for_dt(1.0).unwrap().to_dense();
        assert!((dense[0][0].re + 2e-5).abs() < 1e-18);
        assert!((dense[1][0].re - 2e-5).abs() < 1e-18);
        assert!((dense[2][0].re - 2e-5).abs() < 1e-18);
    }

    #[test]
    fn invalid_half_life_rejected() {
        let bad = ChainNuclide {
            name: "A".into(),
            half_life: Some(0.0),
            ..Default::default()
        };
        let chain = Chain::from_nuclides(vec![bad]).unwrap();
        assert!(matches!(
            DepletionSystem::build(chain, &ReactionRates::new()),
            Err(Error::InvalidHalfLife { name, value }) if name == "A" && value == 0.0
        ));

        let bad = ChainNuclide {
            name: "A".into(),
            half_life: Some(-1.0),
            ..Default::default()
        };
        let chain = Chain::from_nuclides(vec![bad]).unwrap();
        assert!(DepletionSystem::build(chain, &ReactionRates::new()).is_err());
    }

    #[test]
    fn decay_branching_ratios_renormalized() {
        // Sum is 0.6; the largest branch is adjusted to make the total 1.0.
        let a = ChainNuclide {
            name: "A".into(),
            half_life: Some(std::f64::consts::LN_2 / 1e-6),
            decay_modes: vec![
                DecayMode {
                    kind: "beta".into(),
                    target: "B".into(),
                    branching_ratio: 0.3,
                },
                DecayMode {
                    kind: "beta".into(),
                    target: "C".into(),
                    branching_ratio: 0.3,
                },
            ],
            ..Default::default()
        };
        let chain = Chain::from_nuclides(vec![a]).unwrap();
        assert!((chain.nuclides[0].decay_modes[0].branching_ratio - 0.7).abs() < 1e-12);
        assert!((chain.nuclides[0].decay_modes[1].branching_ratio - 0.3).abs() < 1e-12);
    }

    #[test]
    fn invalid_dt_rejected() {
        let chain = Chain::from_xml(&xml_of(simple_chain())).unwrap();
        let sys = DepletionSystem::build(chain, &ReactionRates::new()).unwrap();
        let n0 = vec![1.0e12, 0.0, 0.0];
        assert!(crate::cram(&sys, Order::Order16, &n0, 0.0).is_err());
        assert!(crate::cram(&sys, Order::Order16, &n0, -1.0).is_err());
        assert!(crate::cram(&sys, Order::Order16, &n0, f64::NAN).is_err());
        assert!(crate::cram(&sys, Order::Order16, &n0, f64::INFINITY).is_err());
    }

    #[test]
    fn fission_yield_parent_reference_resolves() {
        // CASL/VERA chain style: a nuclide borrows another nuclide's yields
        // via <neutron_fission_yields parent="U235"/>.
        let xml = r#"<depletion_chain>
  <nuclide name="U235" reactions="1">
    <reaction type="fission" Q="2.0e8"/>
    <neutron_fission_yields>
      <energies>0.0253 1.4e7</energies>
      <fission_yields energy="0.0253">
        <products>I135 Xe135</products>
        <data>0.03 0.06</data>
      </fission_yields>
      <fission_yields energy="1.4e7">
        <products>I135 Xe135</products>
        <data>0.04 0.05</data>
      </fission_yields>
    </neutron_fission_yields>
  </nuclide>
  <nuclide name="Cm247" reactions="1">
    <reaction type="fission" Q="2.0e8"/>
    <neutron_fission_yields parent="U235"/>
  </nuclide>
  <nuclide name="I135" reactions="0"/>
  <nuclide name="Xe135" reactions="0"/>
</depletion_chain>"#;
        let chain = Chain::from_xml(xml).unwrap();
        let cm = &chain.nuclides[chain.index_of("Cm247").unwrap()];
        assert_eq!(cm.neutron_fission_yields.len(), 2);
        assert_eq!(cm.neutron_fission_yields[0].products["I135"], 0.03);
        assert_eq!(cm.neutron_fission_yields[1].products["Xe135"], 0.05);
    }

    #[test]
    fn fission_yield_missing_parent_errors() {
        let xml = r#"<depletion_chain>
  <nuclide name="Cm247" reactions="1">
    <reaction type="fission" Q="2.0e8"/>
    <neutron_fission_yields parent="U235"/>
  </nuclide>
</depletion_chain>"#;
        match Chain::from_xml(xml) {
            Err(Error::BadStructure(m)) => assert!(m.contains("U235")),
            other => panic!("expected BadStructure, got {other:?}"),
        }
    }

    #[test]
    fn from_xml_keeps_branching_ratios_verbatim() {
        // OpenMC renormalizes decay branching ratios only when *generating*
        // a chain from ENDF; `Chain.from_xml` uses the written values
        // verbatim. The CASL/VERA chain relies on this (e.g. I128 beta- to
        // Xe128 with branching ratio 0.931, the remainder leaving the chain).
        let xml = r#"<depletion_chain>
  <nuclide name="I128" half_life="1499.4" reactions="0">
    <decay type="beta-" target="Xe128" branching_ratio="0.931"/>
  </nuclide>
  <nuclide name="Xe128" reactions="0"/>
</depletion_chain>"#;
        let chain = Chain::from_xml(xml).unwrap();
        let i128 = &chain.nuclides[chain.index_of("I128").unwrap()];
        assert_eq!(i128.decay_modes[0].branching_ratio, 0.931);
    }

    #[test]
    fn fission_uses_lowest_energy_yields() {
        // Yield blocks deliberately listed high-energy first: the matrix must
        // still use the lowest-energy set (OpenMC get_default_fission_yields).
        let xml = r#"<depletion_chain>
  <nuclide name="U235" reactions="1">
    <reaction type="fission" Q="2.0e8"/>
    <neutron_fission_yields>
      <energies>0.0253 1.4e7</energies>
      <fission_yields energy="1.4e7">
        <products>I135</products>
        <data>0.04</data>
      </fission_yields>
      <fission_yields energy="0.0253">
        <products>I135</products>
        <data>0.03</data>
      </fission_yields>
    </neutron_fission_yields>
  </nuclide>
  <nuclide name="I135" reactions="0"/>
</depletion_chain>"#;
        let chain = Chain::from_xml(xml).unwrap();
        let mut rates = ReactionRates::new();
        rates
            .entry(0usize)
            .or_default()
            .insert("fission".to_string(), 1e-5);
        let sys = DepletionSystem::build(chain, &rates).unwrap();
        let dense = sys.matrix_for_dt(1.0).unwrap().to_dense();
        assert!((dense[1][0].re - 3e-7).abs() < 1e-18);
    }

    #[test]
    fn coverage_chain_error_display_and_accessors() {
        // chain.rs Display arms (Io/Xml/UnknownNuclide/BadStructure/InvalidHalfLife).
        let errs = vec![
            (Error::Io("missing.xml: not found".into()), "io error"),
            (Error::Xml("bad float".into()), "XML error"),
            (
                Error::UnknownNuclide {
                    name: "X".into(),
                    context: "decay target",
                },
                "unknown nuclide `X`",
            ),
            (Error::BadStructure("nope".into()), "malformed chain"),
            (
                Error::InvalidHalfLife {
                    name: "A".into(),
                    value: -1.0,
                },
                "invalid half-life",
            ),
        ];
        for (e, needle) in errs {
            assert!(format!("{e}").contains(needle), "{e:?}");
            let _: &dyn std::error::Error = &e;
        }

        // len / is_empty / index_of.
        let chain = simple_chain();
        assert_eq!(chain.len(), 3);
        assert!(!chain.is_empty());
        assert_eq!(chain.index_of("A"), Some(0));
        assert_eq!(chain.index_of("Nope"), None);
        let empty = Chain::from_nuclides(vec![]).unwrap();
        assert!(empty.is_empty());
    }

    #[test]
    fn coverage_chain_from_xml_rejections() {
        // Wrong root.
        assert!(matches!(
            Chain::from_xml("<wrong/>"),
            Err(Error::BadStructure(_))
        ));
        // Unparseable XML.
        assert!(matches!(
            Chain::from_xml("<depletion_chain>"),
            Err(Error::Xml(_))
        ));
        // Empty chain.
        assert!(matches!(
            Chain::from_xml("<depletion_chain></depletion_chain>"),
            Err(Error::BadStructure(_))
        ));
        // Duplicate nuclide.
        let dup = r#"<depletion_chain>
  <nuclide name="A"/><nuclide name="A"/>
</depletion_chain>"#;
        assert!(matches!(Chain::from_xml(dup), Err(Error::BadStructure(_))));
        // Missing required attribute (nuclide name).
        let no_name = r#"<depletion_chain><nuclide/></depletion_chain>"#;
        assert!(matches!(Chain::from_xml(no_name), Err(Error::Xml(_))));
        // Bad float for half-life.
        let bad_float = r#"<depletion_chain><nuclide name="A" half_life="abc"/></depletion_chain>"#;
        assert!(matches!(Chain::from_xml(bad_float), Err(Error::Xml(_))));
        // Bad branching_ratio float on decay.
        let bad_br = r#"<depletion_chain><nuclide name="A">
  <decay type="beta" target="B" branching_ratio="abc"/>
</nuclide><nuclide name="B"/></depletion_chain>"#;
        assert!(matches!(Chain::from_xml(bad_br), Err(Error::Xml(_))));
        // Bad Q float on reaction.
        let bad_q = r#"<depletion_chain><nuclide name="A">
  <reaction type="(n,gamma)" target="B" Q="abc"/>
</nuclide><nuclide name="B"/></depletion_chain>"#;
        assert!(matches!(Chain::from_xml(bad_q), Err(Error::Xml(_))));
        // Fission yields products/data length mismatch.
        let mismatch = r#"<depletion_chain><nuclide name="U235">
  <reaction type="fission" Q="2.0e8"/>
  <neutron_fission_yields>
    <energies>0.0253</energies>
    <fission_yields energy="0.0253"><products>A B</products><data>0.1</data></fission_yields>
  </neutron_fission_yields>
</nuclide></depletion_chain>"#;
        assert!(matches!(
            Chain::from_xml(mismatch),
            Err(Error::BadStructure(_))
        ));
        // Fission yields missing <energies>.
        let no_energies = r#"<depletion_chain><nuclide name="U235">
  <neutron_fission_yields>
    <fission_yields energy="0.0253"><products>A</products><data>0.1</data></fission_yields>
  </neutron_fission_yields>
</nuclide></depletion_chain>"#;
        assert!(matches!(
            Chain::from_xml(no_energies),
            Err(Error::BadStructure(_))
        ));
        // Missing file maps to Io.
        assert!(matches!(
            Chain::from_file("/nonexistent/path/chain.xml"),
            Err(Error::Io(_))
        ));
    }

    #[test]
    fn coverage_chain_parent_cycle_errors() {
        // A borrows from B, B borrows from A: hop cap trips.
        let xml = r#"<depletion_chain>
  <nuclide name="A"><neutron_fission_yields parent="B"/></nuclide>
  <nuclide name="B"><neutron_fission_yields parent="A"/></nuclide>
</depletion_chain>"#;
        assert!(matches!(Chain::from_xml(xml), Err(Error::BadStructure(_))));
    }

    #[test]
    fn coverage_chain_normalize_largest_second() {
        // Largest branch is second: exercises the max_idx update arm.
        let a = ChainNuclide {
            name: "A".into(),
            half_life: Some(std::f64::consts::LN_2 / 1e-6),
            decay_modes: vec![
                DecayMode {
                    kind: "beta".into(),
                    target: "B".into(),
                    branching_ratio: 0.2,
                },
                DecayMode {
                    kind: "beta".into(),
                    target: "C".into(),
                    branching_ratio: 0.4,
                },
            ],
            ..Default::default()
        };
        let chain = Chain::from_nuclides(vec![a]).unwrap();
        let sum: f64 = chain.nuclides[0]
            .decay_modes
            .iter()
            .map(|m| m.branching_ratio)
            .sum();
        assert!((sum - 1.0).abs() < 1e-12);
        // Second (largest) branch absorbed the correction: 0.4 + 0.4 = 0.8.
        assert!((chain.nuclides[0].decay_modes[1].branching_ratio - 0.8).abs() < 1e-12);
    }

    #[test]
    fn coverage_matrix_decay_edge_cases() {
        // Zero-branch decay entry is skipped (paired with a 1.0 branch so
        // from_nuclides renormalization leaves the zero branch untouched).
        let a = ChainNuclide {
            name: "A".into(),
            half_life: Some(std::f64::consts::LN_2 / 1e-6),
            decay_modes: vec![
                DecayMode {
                    kind: "beta".into(),
                    target: "B".into(),
                    branching_ratio: 0.0,
                },
                DecayMode {
                    kind: "beta".into(),
                    target: "C".into(),
                    branching_ratio: 1.0,
                },
            ],
            ..Default::default()
        };
        let b = ChainNuclide {
            name: "B".into(),
            ..Default::default()
        };
        let c = ChainNuclide {
            name: "C".into(),
            ..Default::default()
        };
        let sys = DepletionSystem::build(
            Chain::from_nuclides(vec![a, b, c]).unwrap(),
            &ReactionRates::new(),
        )
        .unwrap();
        let dense = sys.matrix_for_dt(1.0).unwrap().to_dense();
        assert!(dense[1][0].re.abs() < 1e-18);
        assert!((dense[2][0].re - 1e-6).abs() < 1e-18);

        // Unknown decay target errors.
        let a = ChainNuclide {
            name: "A".into(),
            half_life: Some(std::f64::consts::LN_2 / 1e-6),
            decay_modes: vec![DecayMode {
                kind: "beta".into(),
                target: "Missing".into(),
                branching_ratio: 1.0,
            }],
            ..Default::default()
        };
        assert!(matches!(
            DepletionSystem::build(
                Chain::from_nuclides(vec![a]).unwrap(),
                &ReactionRates::new()
            ),
            Err(Error::UnknownNuclide { .. })
        ));

        // Spontaneous fission daughter is skipped (no gain), alpha without
        // He4 and proton without H1 exercise the missing-light-nuclide arms.
        let a = ChainNuclide {
            name: "A".into(),
            half_life: Some(std::f64::consts::LN_2 / 1e-6),
            decay_modes: vec![
                DecayMode {
                    kind: "sf".into(),
                    target: "B".into(),
                    branching_ratio: 0.5,
                },
                DecayMode {
                    kind: "alpha".into(),
                    target: "B".into(),
                    branching_ratio: 0.25,
                },
                DecayMode {
                    kind: "p".into(),
                    target: "B".into(),
                    branching_ratio: 0.25,
                },
            ],
            ..Default::default()
        };
        let b = ChainNuclide {
            name: "B".into(),
            ..Default::default()
        };
        let sys = DepletionSystem::build(
            Chain::from_nuclides(vec![a, b]).unwrap(),
            &ReactionRates::new(),
        )
        .unwrap();
        let dense = sys.matrix_for_dt(1.0).unwrap().to_dense();
        // Only the two non-sf daughters contribute: 0.5 * lambda total.
        assert!((dense[1][0].re - 0.5e-6).abs() < 1e-18);
    }

    #[test]
    fn coverage_matrix_fission_and_reaction_edges() {
        // Zero fission yield is skipped.
        let u = ChainNuclide {
            name: "U235".into(),
            reactions: vec![Reaction {
                kind: "fission".into(),
                target: None,
                q: 0.0,
                branching_ratio: 1.0,
            }],
            neutron_fission_yields: vec![FissionYields {
                energy: 0.0253,
                products: [("I135".to_string(), 0.0)].into_iter().collect(),
            }],
            ..Default::default()
        };
        let i = ChainNuclide {
            name: "I135".into(),
            ..Default::default()
        };
        let mut rates = ReactionRates::new();
        rates
            .entry(0usize)
            .or_default()
            .insert("fission".to_string(), 1e-5);
        let sys =
            DepletionSystem::build(Chain::from_nuclides(vec![u, i]).unwrap(), &rates).unwrap();
        let dense = sys.matrix_for_dt(1.0).unwrap().to_dense();
        assert!(dense[1][0].re.abs() < 1e-18);

        // Unknown fission-yield product errors.
        let u = ChainNuclide {
            name: "U235".into(),
            reactions: vec![Reaction {
                kind: "fission".into(),
                target: None,
                q: 0.0,
                branching_ratio: 1.0,
            }],
            neutron_fission_yields: vec![FissionYields {
                energy: 0.0253,
                products: [("Missing".to_string(), 0.5)].into_iter().collect(),
            }],
            ..Default::default()
        };
        assert!(matches!(
            DepletionSystem::build(Chain::from_nuclides(vec![u]).unwrap(), &rates),
            Err(Error::UnknownNuclide { .. })
        ));

        // Yield-less fission with explicit target behaves as transmutation.
        let u = ChainNuclide {
            name: "U235".into(),
            reactions: vec![Reaction {
                kind: "fission".into(),
                target: Some("I135".into()),
                q: 0.0,
                branching_ratio: 1.0,
            }],
            ..Default::default()
        };
        let i = ChainNuclide {
            name: "I135".into(),
            ..Default::default()
        };
        let sys =
            DepletionSystem::build(Chain::from_nuclides(vec![u, i]).unwrap(), &rates).unwrap();
        let dense = sys.matrix_for_dt(1.0).unwrap().to_dense();
        assert!((dense[1][0].re - 1e-5).abs() < 1e-18);

        // Yield-less fission with unknown target errors.
        let u = ChainNuclide {
            name: "U235".into(),
            reactions: vec![Reaction {
                kind: "fission".into(),
                target: Some("Missing".into()),
                q: 0.0,
                branching_ratio: 1.0,
            }],
            ..Default::default()
        };
        assert!(matches!(
            DepletionSystem::build(Chain::from_nuclides(vec![u]).unwrap(), &rates),
            Err(Error::UnknownNuclide { .. })
        ));

        // Unknown reaction target errors; unknown reaction kind yields no
        // secondaries (covers the `_ => &[]` arm) and missing light nuclide
        // covers the `index_of -> None` arm.
        let a = ChainNuclide {
            name: "A".into(),
            reactions: vec![Reaction {
                kind: "(n,gamma)".into(),
                target: Some("Missing".into()),
                q: 0.0,
                branching_ratio: 1.0,
            }],
            ..Default::default()
        };
        let mut rates2 = ReactionRates::new();
        rates2
            .entry(0usize)
            .or_default()
            .insert("(n,gamma)".to_string(), 1e-5);
        assert!(matches!(
            DepletionSystem::build(Chain::from_nuclides(vec![a]).unwrap(), &rates2),
            Err(Error::UnknownNuclide { .. })
        ));

        let a = ChainNuclide {
            name: "A".into(),
            reactions: vec![Reaction {
                kind: "(n,weird)".into(),
                target: Some("B".into()),
                q: 0.0,
                branching_ratio: 1.0,
            }],
            ..Default::default()
        };
        let b = ChainNuclide {
            name: "B".into(),
            ..Default::default()
        };
        let sys = DepletionSystem::build(Chain::from_nuclides(vec![a, b]).unwrap(), &rates2);
        // rates2 has no "(n,weird)" entry, so the channel contributes zero.
        assert!(sys.is_ok());

        // (n,a) without He4 in the chain: secondary lookup misses cleanly.
        let a = ChainNuclide {
            name: "A".into(),
            reactions: vec![Reaction {
                kind: "(n,a)".into(),
                target: Some("B".into()),
                q: 0.0,
                branching_ratio: 1.0,
            }],
            ..Default::default()
        };
        let b = ChainNuclide {
            name: "B".into(),
            ..Default::default()
        };
        let mut rates3 = ReactionRates::new();
        rates3
            .entry(0usize)
            .or_default()
            .insert("(n,a)".to_string(), 1e-5);
        let sys =
            DepletionSystem::build(Chain::from_nuclides(vec![a, b]).unwrap(), &rates3).unwrap();
        let dense = sys.matrix_for_dt(1.0).unwrap().to_dense();
        assert!((dense[1][0].re - 1e-5).abs() < 1e-18);
    }

    #[test]
    fn coverage_matrix_helpers_debug_shift() {
        let chain = simple_chain();
        let sys = DepletionSystem::build(chain, &ReactionRates::new()).unwrap();
        // shifted_values wrapper, diagonal entries, Debug impl.
        let vals = sys.shifted_values(1.0, nucleide_linalg::C64_ZERO);
        assert_eq!(vals.len(), sys.entries.len());
        assert_eq!(sys.diagonal_entries().len(), sys.chain.len());
        let dbg = format!("{sys:?}");
        assert!(dbg.contains("DepletionSystem"));
        // deplete() unknown-n0 errors.
        let mut n0 = BTreeMap::new();
        n0.insert("Missing".to_string(), 1.0);
        assert!(crate::deplete(&sys, crate::Order::Order16, &n0, 1.0).is_err());
    }

    #[test]
    fn coverage_matrix_remaining_arms() {
        // Non-nuclide top-level element is skipped (chain.rs `continue`).
        let xml = r#"<depletion_chain>
  <source>metadata ignored</source>
  <nuclide name="A"/>
</depletion_chain>"#;
        assert_eq!(Chain::from_xml(xml).unwrap().len(), 1);

        // Bad reaction branching_ratio float.
        let bad = r#"<depletion_chain><nuclide name="A">
  <reaction type="(n,gamma)" target="B" branching_ratio="abc"/>
</nuclide><nuclide name="B"/></depletion_chain>"#;
        assert!(matches!(Chain::from_xml(bad), Err(Error::Xml(_))));

        // Fission with no yields and no target: pure loss, no gain.
        let u = ChainNuclide {
            name: "U235".into(),
            reactions: vec![Reaction {
                kind: "fission".into(),
                target: None,
                q: 0.0,
                branching_ratio: 1.0,
            }],
            ..Default::default()
        };
        let mut rates = ReactionRates::new();
        rates
            .entry(0usize)
            .or_default()
            .insert("fission".to_string(), 1e-5);
        let sys = DepletionSystem::build(Chain::from_nuclides(vec![u]).unwrap(), &rates).unwrap();
        let dense = sys.matrix_for_dt(1.0).unwrap().to_dense();
        assert!((dense[0][0].re + 1e-5).abs() < 1e-18);

        // Target-None reaction with nonzero rate covers the `None` arm.
        let a = ChainNuclide {
            name: "A".into(),
            reactions: vec![Reaction {
                kind: "(n,gamma)".into(),
                target: None,
                q: 0.0,
                branching_ratio: 1.0,
            }],
            ..Default::default()
        };
        let mut rates = ReactionRates::new();
        rates
            .entry(0usize)
            .or_default()
            .insert("(n,gamma)".to_string(), 2e-5);
        let sys = DepletionSystem::build(Chain::from_nuclides(vec![a]).unwrap(), &rates).unwrap();
        let dense = sys.matrix_for_dt(1.0).unwrap().to_dense();
        assert!((dense[0][0].re + 2e-5).abs() < 1e-18);

        // Unknown reaction kind with nonzero rate covers `_ => &[]`.
        let a = ChainNuclide {
            name: "A".into(),
            reactions: vec![Reaction {
                kind: "(n,weird)".into(),
                target: Some("B".into()),
                q: 0.0,
                branching_ratio: 1.0,
            }],
            ..Default::default()
        };
        let b = ChainNuclide {
            name: "B".into(),
            ..Default::default()
        };
        let mut rates = ReactionRates::new();
        rates
            .entry(0usize)
            .or_default()
            .insert("(n,weird)".to_string(), 3e-5);
        let sys =
            DepletionSystem::build(Chain::from_nuclides(vec![a, b]).unwrap(), &rates).unwrap();
        let dense = sys.matrix_for_dt(1.0).unwrap().to_dense();
        assert!((dense[1][0].re - 3e-5).abs() < 1e-18);
    }

    // helper: serialize Chain back to minimal XML so build() gets owned data
    fn xml_of(chain: Chain) -> String {
        let mut s = String::from("<depletion_chain>\n");
        for n in &chain.nuclides {
            if let Some(t) = n.half_life {
                s += &format!("  <nuclide name=\"{}\" half_life=\"{t:e}\">\n", n.name);
            } else {
                s += &format!("  <nuclide name=\"{}\">\n", n.name);
            }
            for d in &n.decay_modes {
                s += &format!(
                    "    <decay type=\"{}\" target=\"{}\" branching_ratio=\"{}\"/>\n",
                    d.kind, d.target, d.branching_ratio
                );
            }
            for r in &n.reactions {
                match &r.target {
                    Some(t) => {
                        s += &format!(
                            "    <reaction type=\"{}\" target=\"{t}\" branching_ratio=\"{}\"/>\n",
                            r.kind, r.branching_ratio
                        )
                    }
                    None => {
                        s += &format!(
                            "    <reaction type=\"{}\" branching_ratio=\"{}\"/>\n",
                            r.kind, r.branching_ratio
                        )
                    }
                }
            }
            s += "  </nuclide>\n";
        }
        s += "</depletion_chain>";
        s
    }

    // silence unused warning for helper
    #[allow(dead_code)]
    fn _unused(_: Vec<usize>) {}
}
