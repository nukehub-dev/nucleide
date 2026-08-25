//! Depletion: CRAM solver, chain files, results.
//!
//! Key decisions:
//! - chain files use the depletion-chain XML format (see [`chain`])
//! - backend is `faer` sparse LU with symbolic reuse, behind the `linalg` facade

pub mod chain;
pub mod cram;
pub mod matrix;

pub use chain::{Chain, ChainNuclide, DecayMode, Error, FissionYields, Reaction};
pub use cram::{cram, cram_with_symbolic, Order};
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
        rates.insert((1usize, "(n,gamma)".to_string()), 1e-5);
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
        let sym = linalg::SymbolicLu::try_new(&sys.pattern).unwrap();
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
                        s += &format!("    <reaction type=\"{}\" target=\"{t}\"/>\n", r.kind)
                    }
                    None => s += &format!("    <reaction type=\"{}\"/>\n", r.kind),
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
