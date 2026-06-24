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

use std::collections::{BTreeSet, HashMap};
use std::fmt::Debug;
use std::hash::Hash;
use std::marker::PhantomData;

use xxhash_rust::xxh3::{xxh3_128_with_seed, xxh3_64_with_seed};

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

/// A circular refinement algorithm. One graph algorithm today (extended
/// connectivity); the hash recipe is a parameter ([`EcScheme`]), not an
/// algorithmic alternative.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CircularRefinementAlgorithm {
    /// Extended-connectivity refinement (the Morgan-algorithm family underlying
    /// ECFP / Morgan fingerprints): iterate `radius` rounds, hashing each node from
    /// the round, its previous identifier, and its sorted `(edge label, neighbor's
    /// previous identifier)` pairs; then remove structurally-duplicate features.
    Ec { radius: u32, scheme: EcScheme },
}

/// Hash recipe for extended-connectivity refinement. The graph algorithm is the
/// same for every variant; only the hashing differs. The paper/tool leaves the
/// hash to the implementer, so these are frozen choices.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EcScheme {
    /// Rogers & Hahn 2010: `xxh3_64` over the round array, seeded by `seed`.
    RogersHahn { seed: u64 },
    /// RDKit Morgan: the vendored 32-bit boost hash with incremental combine.
    Morgan,
}

impl Graph {
    /// Circular refinement; dispatches to the selected algorithm.
    pub fn circular_refine(
        &self,
        node_components: impl Fn(NodeId) -> Vec<u32>,
        edge_label: impl Fn(EdgeId) -> u32,
        algorithm: CircularRefinementAlgorithm,
    ) -> Vec<u64> {
        match algorithm {
            CircularRefinementAlgorithm::Ec { radius, scheme } => {
                self.circular_refine_ec(node_components, edge_label, radius, scheme)
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
    fn circular_refine_ec(
        &self,
        node_components: impl Fn(NodeId) -> Vec<u32>,
        edge_label: impl Fn(EdgeId) -> u32,
        radius: u32,
        scheme: EcScheme,
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

impl EcScheme {
    /// Hash a node's invariant components into its round-0 identifier.
    fn seed_hash(&self, components: &[u32]) -> u64 {
        match self {
            EcScheme::RogersHahn { seed } => {
                let mut buffer = Vec::with_capacity(components.len() * 4);
                for &component in components {
                    buffer.extend_from_slice(&component.to_le_bytes());
                }
                xxh3_64_with_seed(&buffer, *seed)
            }
            EcScheme::Morgan => u64::from(gboost_hash(components)),
        }
    }

    /// Combine a node's previous identifier and its sorted neighbor pairs into the
    /// next-round identifier.
    fn combine(&self, round: u32, current: u64, neighbors: &[(u32, u64)]) -> u64 {
        match self {
            EcScheme::RogersHahn { seed } => {
                let mut buffer = Vec::with_capacity(16 + neighbors.len() * 12);
                buffer.extend_from_slice(&u64::from(round).to_le_bytes());
                buffer.extend_from_slice(&current.to_le_bytes());
                for &(edge, color) in neighbors {
                    buffer.extend_from_slice(&edge.to_le_bytes());
                    buffer.extend_from_slice(&color.to_le_bytes());
                }
                xxh3_64_with_seed(&buffer, *seed)
            }
            EcScheme::Morgan => {
                // RDKit salts with the 0-based layer (round - 1) and hashes in 32 bits.
                let mut invariant = gboost_combine(round - 1, current as u32);
                for &(edge, color) in neighbors {
                    invariant = gboost_combine(invariant, gboost_hash(&[edge, color as u32]));
                }
                u64::from(invariant)
            }
        }
    }
}

/// RDKit's vendored 32-bit `boost::hash_combine` (frozen formula).
fn gboost_combine(seed: u32, value: u32) -> u32 {
    seed ^ value
        .wrapping_add(0x9e37_79b9)
        .wrapping_add(seed << 6)
        .wrapping_add(seed >> 2)
}

/// `boost::hash` over a sequence of 32-bit values: combine each from seed 0.
fn gboost_hash(values: &[u32]) -> u32 {
    let mut seed = 0;
    for &value in values {
        seed = gboost_combine(seed, value);
    }
    seed
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

    #[rstest]
    fn test_graph_circular_refine() {
        // Path 0-1-2, uniform labels: round 0 contributes 3 equal seeds; round 1 adds
        // three features (two symmetric ends + a distinct middle, each with a distinct
        // bond set, so none are deduped) — 6 identifiers, 3 distinct values.
        let graph = Graph::new(3, &[[0, 1], [1, 2]]);
        let ids = graph.circular_refine(
            |_| vec![1u32],
            |_| 1,
            CircularRefinementAlgorithm::Ec { radius: 1, scheme: EcScheme::Morgan },
        );
        assert_eq!(ids.len(), 6);
        let distinct: BTreeSet<u64> = ids.iter().copied().collect();
        assert_eq!(distinct.len(), 3);
    }

    // RDKit 2026.03.3 connectivity invariants (radius-0 Morgan ids) for these
    // component vectors; pins the gboost hash bit-exactly to RDKit.
    #[rstest]
    #[case(&[6, 4, 3, 0, 0], 2246728737)]
    #[case(&[6, 4, 2, 0, 0], 2245384272)]
    #[case(&[8, 2, 1, 0, 0], 864662311)]
    #[case(&[6, 4, 2, 0, 0, 1], 2968968094)]
    fn test_gboost_hash(#[case] components: &[u32], #[case] expected: u32) {
        assert_eq!(gboost_hash(components), expected);
    }

    #[rstest]
    fn test_graph_circular_refine_distinguishes_components() {
        // Radius 0: only the seed hash of the components; distinct components must give
        // distinct identifiers.
        let graph = Graph::new(2, &[[0, 1]]);
        let ids = graph.circular_refine(
            |n: NodeId| vec![6 + n.0],
            |_| 1,
            CircularRefinementAlgorithm::Ec { radius: 0, scheme: EcScheme::Morgan },
        );
        assert_eq!(ids.len(), 2);
        assert_ne!(ids[0], ids[1]);
    }
}
