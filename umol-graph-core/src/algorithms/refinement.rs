//! Weisfeiler–Leman and circular refinement.
//!
//! Each node is recolored by its own color plus the multiset of its incident
//! `(edge label, neighbor color)` pairs, iterated until the partition stops
//! refining (or a fixed number of rounds). The stable coloring is the equitable
//! partition; an order-independent digest of it is a fast graph invariant —
//! **sound but not complete**: 1-WL fails to separate some non-isomorphic graphs
//! (e.g. C6 vs two triangles), so a hash match must be confirmed against an exact
//! canonical form (`AutomorphismOutput`).
//!
//! The hashing is pluggable through [`RefinementHash`] (plain refinement) and
//! [`CircularRefinementHash`] (extended connectivity): a custom impl can reproduce
//! an external scheme exactly (the iterative-recoloring part — invariants, folding,
//! and dedup live in the caller). Concrete frozen schemes live downstream.
//!
//! The current selectors provide one-dimensional Weisfeiler--Leman refinement
//! and extended-connectivity circular refinement. See
//! [Weisfeiler and Leman (1968)](https://www.iti.zcu.cz/wl2018/pdf/wl_paper_translation.pdf)
//! and [Rogers and Hahn, *Extended-Connectivity Fingerprints*
//! (2010)](https://doi.org/10.1021/ci100050t).

use std::collections::{BTreeSet, HashMap};
use std::fmt::Debug;
use std::hash::Hash;

use crate::algorithms::traversal::TraversalAlgorithm;
use crate::graph::{EdgeId, Graph, NodeId};

/// A refinement algorithm and its configuration. Parameterized variants carry
/// their own parameters (unlike parameter-free unit-variant algorithm enums).
#[derive(Clone, Copy, Debug)]
pub enum RefinementAlgorithm<H> {
    /// 1-dimensional Weisfeiler–Leman (Weisfeiler & Leman 1968).
    WeisfeilerLehman { rounds: RefinementRounds, scheme: H },
}

/// Refinement rounds to run: `ToFixpoint` until stabilization; `Fixed(n)` for exactly `n` rounds.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RefinementRounds {
    ToFixpoint,
    Fixed(u32),
}

/// Hashing scheme for color refinement.
pub trait RefinementHash {
    /// Identifier width — the scheme's native type (`u32`/`u64`/`u128`).
    type Color: Copy + Ord + Eq + Hash + Debug;

    /// Initial color from a caller-supplied node label.
    fn seed(&self, label: u64) -> Self::Color;

    /// New color from the current color and the node's incident
    /// `(edge label, neighbor color)` pairs. The impl owns the aggregation
    /// (sorted vs commutative) and the byte layout, so it can match an external
    /// scheme exactly.
    fn refine(&self, current: Self::Color, neighbors: &[(u64, Self::Color)]) -> Self::Color;

    /// Order-independent digest of a coloring — the dedup key.
    fn graph_hash(&self, colors: &[Self::Color]) -> Self::Color;
}

/// Per-round colorings (`coloring_at(0)` is the seed) and pre-computed digest.
#[derive(Clone, Debug)]
pub struct Refinement<C> {
    colorings: Vec<Vec<C>>,
    digest: C,
}

impl Graph {
    /// Color-refine under `algorithm`, seeding node/edge colors from caller closures.
    pub fn refine<H: RefinementHash>(
        &self,
        node_label: impl Fn(NodeId) -> u64,
        edge_label: impl Fn(EdgeId) -> u64,
        algorithm: RefinementAlgorithm<H>,
    ) -> Refinement<H::Color> {
        match algorithm {
            RefinementAlgorithm::WeisfeilerLehman { rounds, scheme } => {
                self.refine_wl(node_label, edge_label, rounds, &scheme)
            }
        }
    }

    fn refine_wl<H: RefinementHash>(
        &self,
        node_label: impl Fn(NodeId) -> u64,
        edge_label: impl Fn(EdgeId) -> u64,
        rounds: RefinementRounds,
        scheme: &H,
    ) -> Refinement<H::Color> {
        let n = self.node_count();
        let seed: Vec<H::Color> = (0..n)
            .map(|i| scheme.seed(node_label(NodeId(i as u32))))
            .collect();
        let mut colorings = vec![seed];

        // 1-WL converges in fewer than n rounds.
        let max_rounds = match rounds {
            RefinementRounds::Fixed(r) => r as usize,
            RefinementRounds::ToFixpoint => n,
        };

        for _ in 0..max_rounds {
            let prev = colorings.last().expect("seed present").clone();
            let mut next = Vec::with_capacity(n);
            for i in 0..n {
                let neighbors: Vec<(u64, H::Color)> = self
                    .neighbors(NodeId(i as u32))
                    .iter()
                    .map(|nb| (edge_label(nb.edge), prev[nb.node.index()]))
                    .collect();
                next.push(scheme.refine(prev[i], &neighbors));
            }
            let converged = matches!(rounds, RefinementRounds::ToFixpoint)
                && cell_count(&next) == cell_count(&prev);
            colorings.push(next);
            if converged {
                break;
            }
        }
        let digest = scheme.graph_hash(colorings.last().expect("seed present"));
        Refinement { colorings, digest }
    }
}

impl<C: Copy + Ord + Eq + Hash> Refinement<C> {
    /// The stable (final) coloring.
    pub fn stable_coloring(&self) -> &[C] {
        self.colorings.last().expect("seed present")
    }

    /// Coloring after `round` rounds (`0` = seed).
    pub fn coloring_at(&self, round: usize) -> &[C] {
        &self.colorings[round]
    }

    /// Refinement rounds performed, excluding the seed.
    pub fn round_count(&self) -> usize {
        self.colorings.len() - 1
    }

    /// Distinct colors in the stable coloring — the equitable partition size.
    pub fn cell_count(&self) -> usize {
        cell_count(self.stable_coloring())
    }

    /// Order-independent graph digest — the dedup key.
    pub fn graph_hash(&self) -> C {
        self.digest
    }

    /// Identifier → occurrence count over every round.
    pub fn counts(&self) -> HashMap<C, u32> {
        let mut counts = HashMap::new();
        for coloring in &self.colorings {
            for &c in coloring {
                *counts.entry(c).or_insert(0) += 1;
            }
        }
        counts
    }

    /// Distinct identifiers over every round (the binary-fingerprint substrate).
    pub fn features(&self) -> Vec<C> {
        let mut features: Vec<C> = self.counts().into_keys().collect();
        features.sort_unstable();
        features
    }
}

fn cell_count<C: Ord + Copy>(coloring: &[C]) -> usize {
    let mut c = coloring.to_vec();
    c.sort_unstable();
    c.dedup();
    c.len()
}

/// Hash recipe for extended-connectivity refinement: `seed_hash` makes a node's
/// round-0 identifier from its invariant components; `combine` makes the next
/// identifier from the round, the previous identifier, and the sorted
/// `(edge label, neighbor identifier)` pairs. The impl owns the byte layout, so it
/// can match an external scheme exactly.
pub trait CircularRefinementHash {
    fn seed_hash(&self, components: &[u32]) -> u64;
    fn combine(&self, round: u32, current: u64, neighbors: &[(u32, u64)]) -> u64;
}

/// A circular refinement algorithm. One graph algorithm today (extended
/// connectivity); the hash recipe is the type parameter `H`, not an algorithmic
/// alternative.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CircularRefinementAlgorithm<H> {
    /// Extended-connectivity refinement (the Morgan-algorithm family underlying
    /// ECFP / Morgan fingerprints): iterate `radius` rounds, hashing each node from
    /// the round, its previous identifier, and its sorted `(edge label, neighbor's
    /// previous identifier)` pairs; then remove structurally-duplicate features.
    Ec { radius: u32, scheme: H },
}

impl Graph {
    /// Circular refinement; dispatches to the selected algorithm.
    pub fn circular_refine<H: CircularRefinementHash>(
        &self,
        node_components: impl Fn(NodeId) -> Vec<u32>,
        edge_label: impl Fn(EdgeId) -> u32,
        algorithm: CircularRefinementAlgorithm<H>,
    ) -> Vec<u64> {
        match algorithm {
            CircularRefinementAlgorithm::Ec { radius, scheme } => {
                self.circular_refine_ec(node_components, edge_label, radius, &scheme)
            }
        }
    }

    /// Extended-connectivity refinement plus structural duplicate removal — the
    /// graph algorithm behind ECFP / Morgan fingerprints. `node_components` gives
    /// each node's invariant component vector, `edge_label` each edge's label
    /// (the caller computes these); `scheme` is the hash.
    ///
    /// Rounds run `0..=radius`: round 0 hashes the components; round r hashes the
    /// round, the previous identifier, and the sorted `(edge label, neighbor's
    /// previous id)` pairs. Features whose covered bond set coincides are then
    /// reduced to the one with the smallest `(round, identifier)`. The result is the
    /// surviving feature identifiers (a multiset — equal ids may recur where
    /// distinct environments collide).
    fn circular_refine_ec<H: CircularRefinementHash>(
        &self,
        node_components: impl Fn(NodeId) -> Vec<u32>,
        edge_label: impl Fn(EdgeId) -> u32,
        radius: u32,
        scheme: &H,
    ) -> Vec<u64> {
        let node_count = self.node_count();

        let mut rounds: Vec<Vec<u64>> = Vec::with_capacity(radius as usize + 1);
        rounds.push(
            (0..node_count)
                .map(|i| scheme.seed_hash(&node_components(NodeId(i as u32))))
                .collect(),
        );
        for round in 1..=radius {
            let next = {
                let previous = rounds.last().expect("round 0 present");
                (0..node_count)
                    .map(|i| {
                        let mut neighbors: Vec<(u32, u64)> = self
                            .neighbors(NodeId(i as u32))
                            .iter()
                            .map(|nb| (edge_label(nb.edge), previous[nb.node.index()]))
                            .collect();
                        neighbors.sort_unstable();
                        scheme.combine(round, previous[i], &neighbors)
                    })
                    .collect()
            };
            rounds.push(next);
        }
        self.remove_duplicate_environments(&rounds, radius)
    }

    /// Rogers & Hahn duplicate-structure removal: round-0 identifiers are kept
    /// directly; for rounds ≥ 1, features whose covered bond set coincides collapse
    /// to the one with the smallest `(round, identifier)`. Bond sets come from a BFS.
    fn remove_duplicate_environments(&self, rounds: &[Vec<u64>], radius: u32) -> Vec<u64> {
        let mut identifiers: Vec<u64> = rounds[0].clone();
        let mut kept: HashMap<Vec<u32>, (u32, u64)> = HashMap::new();
        if radius >= 1 {
            for atom in 0..self.node_count() {
                let source = NodeId(atom as u32);
                // EC duplicate removal is defined over shortest-path radius
                // shells, so BFS is fixed by the operation rather than exposed
                // as an independent refinement choice.
                let neighborhood = self.neighborhood(source, radius - 1, TraversalAlgorithm::Bfs);
                let mut bond_set: BTreeSet<u32> = BTreeSet::new();
                let mut shell = 0;
                for round in 1..=radius {
                    while shell < neighborhood.len() && neighborhood[shell].1 == round - 1 {
                        for neighbor in self.neighbors(neighborhood[shell].0) {
                            bond_set.insert(neighbor.edge.index() as u32);
                        }
                        shell += 1;
                    }
                    let identifier = rounds[round as usize][source.index()];
                    let key: Vec<u32> = bond_set.iter().copied().collect();
                    kept.entry(key)
                        .and_modify(|best| {
                            if (round, identifier) < *best {
                                *best = (round, identifier);
                            }
                        })
                        .or_insert((round, identifier));
                }
            }
        }
        identifiers.extend(kept.values().map(|&(_, identifier)| identifier));
        identifiers
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    fn uniform(_: NodeId) -> u64 {
        0
    }
    fn no_edge_color(_: EdgeId) -> u64 {
        0
    }

    // A custom RefinementHash impl exercises the trait seam (and proves pluggability).
    #[derive(Clone, Copy, Debug)]
    struct CountingScheme;
    impl RefinementHash for CountingScheme {
        type Color = u64;
        fn seed(&self, label: u64) -> u64 {
            label
        }
        fn refine(&self, current: u64, neighbors: &[(u64, u64)]) -> u64 {
            current.wrapping_add(neighbors.len() as u64)
        }
        fn graph_hash(&self, colors: &[u64]) -> u64 {
            colors.iter().sum()
        }
    }

    #[rstest]
    fn test_refinement_custom_scheme() {
        let g = Graph::new(3, &[[0, 1], [1, 2]]);
        let algorithm = RefinementAlgorithm::WeisfeilerLehman {
            rounds: RefinementRounds::Fixed(1),
            scheme: CountingScheme,
        };
        let refinement = g.refine(uniform, no_edge_color, algorithm);
        // round 1 colors = degree: ends 1, middle 2 → sum = 4.
        assert_eq!(refinement.coloring_at(1), &[1, 2, 1]);
        assert_eq!(refinement.graph_hash(), 4);
    }

    // A custom CircularRefinementHash impl exercises that trait seam. Radius 0 returns
    // the per-node seed hash directly (no duplicate-environment removal).
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    struct SummingScheme;
    impl CircularRefinementHash for SummingScheme {
        fn seed_hash(&self, components: &[u32]) -> u64 {
            components.iter().map(|&c| u64::from(c)).sum()
        }
        fn combine(&self, _round: u32, current: u64, neighbors: &[(u32, u64)]) -> u64 {
            current + neighbors.iter().map(|&(_, color)| color).sum::<u64>()
        }
    }

    #[rstest]
    fn test_graph_circular_refine() {
        let graph = Graph::new(2, &[[0, 1]]);
        let ids = graph.circular_refine(
            |n: NodeId| vec![n.0 + 1],
            |_| 1,
            CircularRefinementAlgorithm::Ec {
                radius: 0,
                scheme: SummingScheme,
            },
        );
        assert_eq!(ids, vec![1, 2]);
    }
}
