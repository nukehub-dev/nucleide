//! Transmutation matrix construction from a [`Chain`] plus one-group rates.

use std::collections::{BTreeMap, HashMap, HashSet};

use linalg::Pattern;

use crate::chain::{Chain, Error};

/// One-group reaction rate lookup: `nuclide_index → (reaction_name → rate [1/s])`.
///
/// The two-level map lets the matrix builder look up rates by `&str` without
/// allocating a temporary `String` for each reaction entry.
pub type ReactionRates = HashMap<usize, HashMap<String, f64>>;

/// Assembled depletion system: sparse pattern + unscaled column-stochastic
/// style entries `A[j, i]` (loss on the diagonal, production off-diagonal).
///
/// The pattern explicitly includes every diagonal position so that CRAM can
/// shift by `-theta` without changing sparsity.
#[derive(Clone)]
pub struct DepletionSystem {
    pub chain: Chain,
    pub pattern: Pattern,
    /// Entries parallel to the pattern's entry order.
    pub entries: Vec<Entry>,
    /// Unscaled A values in entry order [1/s].
    pub base_values: Vec<linalg::C64>,
    /// Entry index of each diagonal element (for fast theta shifts).
    diag_entry: Vec<usize>,
}

/// One matrix entry location.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Entry {
    pub row: usize,
    pub col: usize,
    pub is_diagonal: bool,
}

impl DepletionSystem {
    /// Build `A` from the chain. Fission production uses single-energy
    /// yields when present: for nuclide `i` with fission rate `F`,
    /// `A[product][i] += F * y_product`.
    ///
    /// Rates for plain reactions are looked up as `rates[nuclide_idx][kind]` in
    /// `rates`; missing entries contribute zero (pure-loss-only channels).
    pub fn build(chain: Chain, rates: &ReactionRates) -> Result<Self, Error> {
        let n = chain.len();
        // Accumulate dense-map then flatten deterministically.
        let mut acc: BTreeMap<(usize, usize), f64> = BTreeMap::new();
        let mut add = |r: usize, c: usize, v: f64| {
            *acc.entry((r, c)).or_insert(0.0) += v;
        };

        for i in 0..n {
            let nuc = &chain.nuclides[i];
            let lambda = nuc.decay_constant()?;

            if lambda > 0.0 {
                add(i, i, -lambda);
                for mode in &nuc.decay_modes {
                    let branch_val = lambda * mode.branching_ratio;
                    if branch_val == 0.0 {
                        continue;
                    }
                    // Gain from explicit decay daughter (skip spontaneous fission).
                    if !mode.kind.contains("sf") {
                        let j =
                            chain
                                .index_of(&mode.target)
                                .ok_or_else(|| Error::UnknownNuclide {
                                    name: mode.target.clone(),
                                    context: "decay target",
                                })?;
                        add(j, i, branch_val);
                    }
                    // Light-particle secondaries from alpha / proton decay,
                    // mirroring OpenMC chain.py:648-657.
                    if mode.kind.contains("alpha") {
                        if let Some(j) = chain.index_of("He4") {
                            let count = mode.kind.matches("alpha").count();
                            add(j, i, count as f64 * branch_val);
                        }
                    } else if mode.kind.contains('p') {
                        if let Some(j) = chain.index_of("H1") {
                            let count = mode.kind.matches('p').count();
                            add(j, i, count as f64 * branch_val);
                        }
                    }
                }
            }

            // Track reaction types already debited for loss on this nuclide.
            let mut seen_reactions: HashSet<&str> = HashSet::new();
            let nuc_rates = rates.get(&i);

            for reaction in &nuc.reactions {
                let sigma_phi = nuc_rates
                    .and_then(|m| m.get(reaction.kind.as_str()))
                    .copied()
                    .unwrap_or(0.0);
                if sigma_phi == 0.0 {
                    continue;
                }

                // Loss term: subtract once per reaction type, not once per
                // <reaction> entry (OpenMC chain.py:720-723).
                if seen_reactions.insert(reaction.kind.as_str()) {
                    add(i, i, -sigma_phi);
                }

                match reaction.kind.as_str() {
                    "fission" => {
                        // Single-energy yield set drives production.
                        if let Some(fy) = nuc.neutron_fission_yields.first() {
                            for (product, y) in &fy.products {
                                if *y == 0.0 {
                                    continue;
                                }
                                let j = chain.index_of(product).ok_or_else(|| {
                                    Error::UnknownNuclide {
                                        name: product.clone(),
                                        context: "fission yield product",
                                    }
                                })?;
                                add(j, i, sigma_phi * y);
                            }
                        } else if let Some(t) = &reaction.target {
                            // Yield-less fission with explicit target acts
                            // like a transmutation channel.
                            let j = chain.index_of(t).ok_or_else(|| Error::UnknownNuclide {
                                name: t.clone(),
                                context: "reaction target",
                            })?;
                            add(j, i, sigma_phi);
                        }
                    }
                    _ => {
                        let br = reaction.branching_ratio;
                        if let Some(t) = &reaction.target {
                            let j = chain.index_of(t).ok_or_else(|| Error::UnknownNuclide {
                                name: t.clone(),
                                context: "reaction target",
                            })?;
                            add(j, i, sigma_phi * br);
                        }
                        // Light-particle production for reactions like (n,a),
                        // (n,p), (n,d), etc., mirroring OpenMC chain.py:733-738.
                        for secondary in reaction_secondaries(&reaction.kind) {
                            if let Some(j) = chain.index_of(secondary) {
                                add(j, i, sigma_phi * br);
                            }
                        }
                    }
                }
            }
        }

        // Ensure every diagonal exists even for isolated stable nuclides.
        for i in 0..n {
            acc.entry((i, i)).or_insert(0.0);
        }

        let entries: Vec<Entry> = acc
            .keys()
            .map(|&(r, c)| Entry {
                row: r,
                col: c,
                is_diagonal: r == c,
            })
            .collect();
        let diag_entry: Vec<usize> = entries
            .iter()
            .enumerate()
            .filter(|(_, e)| e.is_diagonal)
            .map(|(k, _)| k)
            .collect();
        let base_values: Vec<f64> = acc.values().copied().collect();
        let pattern = Pattern::from_entries(
            n,
            &entries.iter().map(|e| (e.row, e.col)).collect::<Vec<_>>(),
        )
        .map_err(|e| Error::BadStructure(e.to_string()))?;

        Ok(Self {
            chain,
            pattern,
            entries,
            base_values: base_values.into_iter().map(linalg::C64::from).collect(),
            diag_entry,
        })
    }

    /// Values of `A*dt - theta*I` in entry order, written into `out`.
    pub fn shifted_values_into(&self, dt: f64, theta: linalg::C64, out: &mut [linalg::C64]) {
        assert_eq!(out.len(), self.entries.len());
        for (k, (e, v)) in self.entries.iter().zip(&self.base_values).enumerate() {
            let scaled = *v * dt;
            out[k] = if e.is_diagonal {
                scaled - theta
            } else {
                scaled
            };
        }
    }

    /// Values of `A*dt - theta*I` in entry order.
    pub fn shifted_values(&self, dt: f64, theta: linalg::C64) -> Vec<linalg::C64> {
        let mut out = vec![linalg::C64_ZERO; self.entries.len()];
        self.shifted_values_into(dt, theta, &mut out);
        out
    }

    /// Matrix for a given timestep (no shift) — useful for inspection/tests.
    pub fn matrix_for_dt(&self, dt: f64) -> Result<linalg::ComplexCsc, Error> {
        linalg::ComplexCsc::from_entries(
            &self.pattern,
            &self
                .entries
                .iter()
                .zip(&self.base_values)
                .map(|(_, v)| *v * dt)
                .collect::<Vec<_>>(),
        )
        .map_err(|e| Error::BadStructure(e.to_string()))
    }

    /// Entry indices of the diagonal elements.
    pub fn diagonal_entries(&self) -> &[usize] {
        &self.diag_entry
    }
}

impl std::fmt::Debug for DepletionSystem {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DepletionSystem")
            .field("nuclides", &self.chain.len())
            .field("entries", &self.entries.len())
            .finish()
    }
}

/// Light-particle secondaries emitted by a transmutation reaction. Mirrors a
/// minimal subset of OpenMC's `REACTIONS` table (chain.py:33-118) for the
/// reaction types commonly present in depletion chains.
fn reaction_secondaries(kind: &str) -> &'static [&'static str] {
    match kind {
        // Deuterium-producing.
        "(n,2nd)" | "(n,nd)" | "(n,4nd)" | "(n,5nd)" | "(n,6nd)" | "(n,d)" | "(n,dt)"
        | "(n,nda)" => &["H2"],
        // Tritium-producing.
        "(n,nt)" | "(n,2nt)" | "(n,3nt)" | "(n,4nt)" | "(n,5nt)" | "(n,6nt)" | "(n,t)"
        | "(n,ta)" | "(n,nta)" => &["H3"],
        // Helium-3-producing.
        "(n,n3He)" | "(n,3He)" | "(n,2n3He)" | "(n,3n3He)" | "(n,4n3He)" | "(n,p3He)"
        | "(n,3Hea)" => &["He3"],
        // Single alpha-producing.
        "(n,na)" | "(n,2na)" | "(n,3na)" | "(n,4na)" | "(n,5na)" | "(n,6na)" | "(n,7na)"
        | "(n,a)" | "(n,2a)" | "(n,3a)" => &["He4"],
        // Multiple alpha-producing.
        "(n,n3a)" => &["He4", "He4", "He4"],
        "(n,n2a)" | "(n,2n2a)" | "(n,t2a)" | "(n,d2a)" | "(n,4n2a)" | "(n,3n2a)" => &["He4", "He4"],
        // Proton-producing.
        "(n,np)" | "(n,2np)" | "(n,3np)" | "(n,p)" | "(n,4np)" | "(n,5np)" | "(n,6np)"
        | "(n,7np)" => &["H1"],
        // Two-proton-producing.
        "(n,n2p)" | "(n,2p)" | "(n,3n2p)" | "(n,3p)" | "(n,n3p)" | "(n,4n2p)" | "(n,5n2p)"
        | "(n,2n2p)" => &["H1", "H1"],
        // Mixed H1 + He4.
        "(n,npa)" | "(n,pa)" | "(n,2npa)" | "(n,3npa)" | "(n,4npa)" | "(n,3n2pa)" => &["H1", "He4"],
        // Mixed H1 + H2 / H3.
        "(n,npd)" | "(n,pd)" => &["H1", "H2"],
        "(n,npt)" | "(n,pt)" => &["H1", "H3"],
        // Mixed H2 + He4.
        "(n,nd2a)" | "(n,da)" => &["H2", "He4"],
        // Mixed H3 + He4.
        "(n,nt2a)" => &["H3", "He4"],
        // Mixed H2 + H3.
        "(n,ndt)" | "(n,d3He)" => &["H2", "H3"],
        // Mixed H1 + He3.
        "(n,np3He)" => &["H1", "He3"],
        // Mixed H2 + He3.
        "(n,nd3He)" => &["H2", "He3"],
        // Mixed H3 + He3.
        "(n,nt3He)" => &["H3", "He3"],
        // No light secondaries (neutrons and photons are not chain nuclides).
        "(n,2n)" | "(n,3n)" | "(n,4n)" | "(n,5n)" | "(n,6n)" | "(n,7n)" | "(n,8n)"
        | "(n,gamma)" | "fission" => &[],
        _ => &[],
    }
}
