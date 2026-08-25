//! Transmutation matrix construction from a [`Chain`] plus one-group rates.

use std::collections::BTreeMap;

use linalg::Pattern;

use crate::chain::{Chain, Error};

/// One-group reaction rate lookup: `(nuclide_index, reaction_name) → rate [1/s]`.
pub type ReactionRates = BTreeMap<(usize, String), f64>;

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
    /// Rates for plain reactions are looked up as `(nuclide_idx, kind)` in
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
            let lambda = nuc.decay_constant();

            if lambda > 0.0 {
                add(i, i, -lambda);
                for mode in &nuc.decay_modes {
                    let j = chain
                        .index_of(&mode.target)
                        .ok_or_else(|| Error::UnknownNuclide {
                            name: mode.target.clone(),
                            context: "decay target",
                        })?;
                    add(j, i, lambda * mode.branching_ratio);
                }
            }

            for reaction in &nuc.reactions {
                let sigma_phi = rates
                    .get(&(i, reaction.kind.clone()))
                    .copied()
                    .unwrap_or(0.0);
                if sigma_phi == 0.0 {
                    continue;
                }
                match reaction.kind.as_str() {
                    "fission" => {
                        add(i, i, -sigma_phi);
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
                        add(i, i, -sigma_phi);
                        if let Some(t) = &reaction.target {
                            let j = chain.index_of(t).ok_or_else(|| Error::UnknownNuclide {
                                name: t.clone(),
                                context: "reaction target",
                            })?;
                            add(j, i, sigma_phi);
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

    /// Values of `A*dt - theta*I` in entry order.
    pub fn shifted_values(&self, dt: f64, theta: linalg::C64) -> Vec<linalg::C64> {
        self.entries
            .iter()
            .zip(&self.base_values)
            .map(|(e, v)| {
                let scaled = *v * dt;
                if e.is_diagonal {
                    scaled - theta
                } else {
                    scaled
                }
            })
            .collect()
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
