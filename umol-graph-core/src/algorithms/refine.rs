//! Weisfeiler–Leman color refinement (1-WL): equitable partition + a content hash.
//!
//! Each node is recolored by its own color plus the multiset of its incident
//! `(edge label, neighbor color)` pairs, iterated until the partition stops
//! refining (or a fixed number of rounds). The stable coloring is the equitable
//! partition; an order-independent digest of it is a fast graph invariant —
//! **sound but not complete**: 1-WL fails to separate some non-isomorphic graphs
//! (e.g. C6 vs two triangles), so a hash match must be confirmed against an exact
//! canonical form (`Automorphism`).
//!
//! The hashing is pluggable through [`RefinementHash`]: the built-in
//! [`RefinementXxh3Scheme`] families give frozen, reproducible results, and a
//! custom impl can reproduce an external scheme exactly (the iterative-recoloring
//! part of it — invariants/folding/dedup live in the caller).

use std::collections::HashMap;
use std::fmt::Debug;
use std::hash::Hash;
use std::marker::PhantomData;

use xxhash_rust::xxh3::{xxh3_128_with_seed, xxh3_64_with_seed};

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
    Fixed(usize),
}

/// Hashing scheme for color refinement. The built-in [`RefinementXxh3Scheme`] families or
/// custom schemes must implement it.
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
            RefinementRounds::Fixed(r) => r,
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

/// Digest width of a [`RefinementXxh3Scheme`].
pub trait RefinementWidth {
    type Color: Copy + Ord + Eq + Hash + Debug;
    fn hash(seed: u64, bytes: &[u8]) -> Self::Color;
    fn push_le(buf: &mut Vec<u8>, color: Self::Color);
    fn zero() -> Self::Color;
    fn wrapping_add(a: Self::Color, b: Self::Color) -> Self::Color;
}

/// 64-bit digests (`xxh3_64`).
#[derive(Clone, Copy, Debug)]
pub enum RefinementWidth64 {}
impl RefinementWidth for RefinementWidth64 {
    type Color = u64;
    fn hash(seed: u64, bytes: &[u8]) -> u64 {
        xxh3_64_with_seed(bytes, seed)
    }
    fn push_le(buf: &mut Vec<u8>, color: u64) {
        buf.extend_from_slice(&color.to_le_bytes());
    }
    fn zero() -> u64 {
        0
    }
    fn wrapping_add(a: u64, b: u64) -> u64 {
        a.wrapping_add(b)
    }
}

/// 128-bit digests (`xxh3_128`) — the collision budget for large-scale dedup.
#[derive(Clone, Copy, Debug)]
pub enum RefinementWidth128 {}
impl RefinementWidth for RefinementWidth128 {
    type Color = u128;
    fn hash(seed: u64, bytes: &[u8]) -> u128 {
        xxh3_128_with_seed(bytes, seed)
    }
    fn push_le(buf: &mut Vec<u8>, color: u128) {
        buf.extend_from_slice(&color.to_le_bytes());
    }
    fn zero() -> u128 {
        0
    }
    fn wrapping_add(a: u128, b: u128) -> u128 {
        a.wrapping_add(b)
    }
}

/// Neighbor-multiset aggregation. `Sorted` is exact (canonical multiset);
/// `SumSketch` is commutative, O(degree), and incrementally updatable.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RefinementAggregation {
    Sorted,
    SumSketch,
}

/// Built-in reproducible scheme over xxh3: `seed` is the family selector,
/// `W` the width, `aggregation` the multiset combine.
#[derive(Clone, Copy, Debug)]
pub struct RefinementXxh3Scheme<W: RefinementWidth> {
    seed: u64,
    aggregation: RefinementAggregation,
    _width: PhantomData<W>,
}

// Placeholder family seeds — not stable scheme identities yet.
const ALBATROSS_SEED: u64 = 0xA1BA_7305_5EED_0001;
const BULLFINCH_SEED: u64 = 0xB011_F114_5EED_0002;

impl<W: RefinementWidth> RefinementXxh3Scheme<W> {
    pub fn new(seed: u64, aggregation: RefinementAggregation) -> Self {
        Self {
            seed,
            aggregation,
            _width: PhantomData,
        }
    }

    /// Default family: `Sorted` aggregation. Placeholder identity.
    pub fn default_scheme() -> Self {
        Self::new(ALBATROSS_SEED, RefinementAggregation::Sorted)
    }

    /// Placeholder family name — not a stable scheme identity yet.
    pub fn albatross() -> Self {
        Self::new(ALBATROSS_SEED, RefinementAggregation::Sorted)
    }

    /// Placeholder family name — not a stable scheme identity yet.
    pub fn bullfinch() -> Self {
        Self::new(BULLFINCH_SEED, RefinementAggregation::Sorted)
    }
}

impl<W: RefinementWidth> RefinementHash for RefinementXxh3Scheme<W> {
    type Color = W::Color;

    fn seed(&self, label: u64) -> W::Color {
        W::hash(self.seed, &label.to_le_bytes())
    }

    fn refine(&self, current: W::Color, neighbors: &[(u64, W::Color)]) -> W::Color {
        let element = |edge: u64, color: W::Color| {
            let mut buf = Vec::new();
            buf.extend_from_slice(&edge.to_le_bytes());
            W::push_le(&mut buf, color);
            W::hash(self.seed, &buf)
        };
        match self.aggregation {
            RefinementAggregation::Sorted => {
                let mut elements: Vec<W::Color> =
                    neighbors.iter().map(|&(edge, c)| element(edge, c)).collect();
                elements.sort_unstable();
                let mut buf = Vec::new();
                W::push_le(&mut buf, current);
                for e in elements {
                    W::push_le(&mut buf, e);
                }
                W::hash(self.seed, &buf)
            }
            RefinementAggregation::SumSketch => {
                let mut acc = W::zero();
                for &(edge, c) in neighbors {
                    acc = W::wrapping_add(acc, element(edge, c));
                }
                let mut buf = Vec::new();
                W::push_le(&mut buf, current);
                W::push_le(&mut buf, acc);
                W::hash(self.seed, &buf)
            }
        }
    }

    fn graph_hash(&self, colors: &[W::Color]) -> W::Color {
        let mut sorted = colors.to_vec();
        sorted.sort_unstable();
        let mut buf = Vec::new();
        buf.extend_from_slice(&(sorted.len() as u64).to_le_bytes());
        for c in sorted {
            W::push_le(&mut buf, c);
        }
        W::hash(self.seed, &buf)
    }
}

fn cell_count<C: Ord + Copy>(coloring: &[C]) -> usize {
    let mut c = coloring.to_vec();
    c.sort_unstable();
    c.dedup();
    c.len()
}

#[cfg(test)]
mod tests {
    use pretty_assertions::{assert_eq, assert_ne};
    use rstest::*;

    use super::*;

    fn uniform(_: NodeId) -> u64 {
        0
    }
    fn no_edge_color(_: EdgeId) -> u64 {
        0
    }

    fn wl_128(
        rounds: RefinementRounds,
    ) -> RefinementAlgorithm<RefinementXxh3Scheme<RefinementWidth128>> {
        RefinementAlgorithm::WeisfeilerLehman {
            rounds,
            scheme: RefinementXxh3Scheme::default_scheme(),
        }
    }

    fn hash_of(g: &Graph) -> u128 {
        g.refine(uniform, no_edge_color, wl_128(RefinementRounds::ToFixpoint))
            .graph_hash()
    }

    #[rstest]
    #[case::triangle(Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), 1)]
    #[case::path3(Graph::new(3, &[[0, 1], [1, 2]]), 2)]
    #[case::path4(Graph::new(4, &[[0, 1], [1, 2], [2, 3]]), 2)]
    #[case::star4(Graph::new(5, &[[0, 1], [0, 2], [0, 3], [0, 4]]), 2)]
    #[case::six_cycle(Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]), 1)]
    fn test_graph_refine(#[case] g: Graph, #[case] expected_cells: usize) {
        let refinement = g.refine(uniform, no_edge_color, wl_128(RefinementRounds::ToFixpoint));
        assert_eq!(refinement.cell_count(), expected_cells);
    }

    #[rstest]
    fn test_refinement_graph_hash_isomorphic() {
        let a = Graph::new(3, &[[0, 1], [1, 2]]);
        let b = Graph::new(3, &[[2, 1], [1, 0]]);
        assert_eq!(hash_of(&a), hash_of(&b));
    }

    #[rstest]
    fn test_refinement_graph_hash_distinguishes() {
        let path = Graph::new(3, &[[0, 1], [1, 2]]);
        let triangle = Graph::new(3, &[[0, 1], [1, 2], [0, 2]]);
        assert_ne!(hash_of(&path), hash_of(&triangle));
    }

    // 1-WL is sound but not complete: C6 and two triangles are non-isomorphic yet collide.
    #[rstest]
    fn test_refinement_graph_hash_wl_incompleteness() {
        let c6 = Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]);
        let two_c3 = Graph::new(6, &[[0, 1], [1, 2], [0, 2], [3, 4], [4, 5], [3, 5]]);
        assert_eq!(hash_of(&c6), hash_of(&two_c3));
    }

    #[rstest]
    fn test_refinement_graph_hash_colors_matter() {
        let g = Graph::new(2, &[[0, 1]]);
        let plain = hash_of(&g);
        let node_colored = g
            .refine(|nd| nd.index() as u64, no_edge_color, wl_128(RefinementRounds::ToFixpoint))
            .graph_hash();
        let edge_colored = g
            .refine(uniform, |_| 7, wl_128(RefinementRounds::ToFixpoint))
            .graph_hash();
        assert_ne!(plain, node_colored);
        assert_ne!(plain, edge_colored);
    }

    #[rstest]
    fn test_refinement_counts() {
        // path 0-1-2, radius 1: round 0 all-same (3); round 1 splits ends (2) from middle (1).
        let g = Graph::new(3, &[[0, 1], [1, 2]]);
        let refinement = g.refine(uniform, no_edge_color, wl_128(RefinementRounds::Fixed(1)));
        let mut values: Vec<u32> = refinement.counts().into_values().collect();
        values.sort_unstable();
        assert_eq!(values, vec![1, 2, 3]);
    }

    #[rstest]
    fn test_refinement_sum_sketch_matches_partition() {
        // SumSketch must still split path ends from the middle.
        let g = Graph::new(3, &[[0, 1], [1, 2]]);
        let algorithm = RefinementAlgorithm::WeisfeilerLehman {
            rounds: RefinementRounds::ToFixpoint,
            scheme: RefinementXxh3Scheme::<RefinementWidth128>::new(
                ALBATROSS_SEED,
                RefinementAggregation::SumSketch,
            ),
        };
        let refinement = g.refine(uniform, no_edge_color, algorithm);
        assert_eq!(refinement.cell_count(), 2);
    }

    // The 64- and 128-bit widths both distinguish path from triangle.
    fn distinguishes<W: RefinementWidth>() {
        let algorithm = || RefinementAlgorithm::WeisfeilerLehman {
            rounds: RefinementRounds::ToFixpoint,
            scheme: RefinementXxh3Scheme::<W>::default_scheme(),
        };
        let path = Graph::new(3, &[[0, 1], [1, 2]]);
        let triangle = Graph::new(3, &[[0, 1], [1, 2], [0, 2]]);
        let hp = path.refine(uniform, no_edge_color, algorithm()).graph_hash();
        let ht = triangle.refine(uniform, no_edge_color, algorithm()).graph_hash();
        assert_ne!(hp, ht);
    }

    #[rstest]
    fn test_refinement_width64() {
        distinguishes::<RefinementWidth64>();
    }

    #[rstest]
    fn test_refinement_width128() {
        distinguishes::<RefinementWidth128>();
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

    // Freeze guard: the albatross-128 digest of a fixed graph must not drift.
    #[rstest]
    fn test_refinement_frozen_albatross_128() {
        let g = Graph::new(3, &[[0, 1], [1, 2]]);
        let algorithm = RefinementAlgorithm::WeisfeilerLehman {
            rounds: RefinementRounds::ToFixpoint,
            scheme: RefinementXxh3Scheme::<RefinementWidth128>::albatross(),
        };
        let h = g.refine(uniform, no_edge_color, algorithm).graph_hash();
        assert_eq!(h, 313131582038434349855774725390837831516);
    }
}
