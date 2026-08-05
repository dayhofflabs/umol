//! Frozen hashing schemes for reproducibility.

use std::fmt::Debug;
use std::hash::Hash;
use std::marker::PhantomData;

use umol_graph_core::{CircularRefinementHash, RefinementHash};
use xxhash_rust::xxh3::{xxh3_128_with_seed, xxh3_64_with_seed};

/// RDKit's vendored 32-bit `boost::hash_combine` (frozen formula).
pub(crate) fn gboost_combine(seed: u32, value: u32) -> u32 {
    seed ^ value
        .wrapping_add(0x9e37_79b9)
        .wrapping_add(seed << 6)
        .wrapping_add(seed >> 2)
}

/// `boost::hash` over a sequence of 32-bit values: combine each from seed 0.
pub(crate) fn gboost_hash(values: &[u32]) -> u32 {
    let mut seed = 0;
    for &value in values {
        seed = gboost_combine(seed, value);
    }
    seed
}

const WL_XXH3_SORTED_V1_SEED: u64 = 0xA1BA_7305_5EED_0001;
const ECFP_XXH3_64_V1_SEED: u64 = 0xECF0_5EED_0000_0001;
const XXH3_SCHEME_VERSION_1: u16 = 1;
const HASH_WIDTH_64: u16 = 64;

/// Circular-refinement recipe: RDKit Morgan via the 32-bit boost hash with
/// incremental combine. Bit-exact to RDKit 2026.03.x.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Morgan;

impl CircularRefinementHash for Morgan {
    fn seed_hash(&self, components: &[u32]) -> u64 {
        u64::from(gboost_hash(components))
    }

    fn combine(&self, round: u32, current: u64, neighbors: &[(u32, u64)]) -> u64 {
        // RDKit salts with the 0-based layer (round - 1) and hashes in 32 bits.
        let mut invariant = gboost_combine(round - 1, current as u32);
        for &(edge, color) in neighbors {
            invariant = gboost_combine(invariant, gboost_hash(&[edge, color as u32]));
        }
        u64::from(invariant)
    }
}

/// Circular-refinement recipe: Rogers & Hahn 2010 ECFP via `xxh3_64` over the round
/// array, seeded by `seed`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RogersHahn {
    pub seed: u64,
}

impl CircularRefinementHash for RogersHahn {
    fn seed_hash(&self, components: &[u32]) -> u64 {
        let mut buffer = Vec::with_capacity(components.len() * 4);
        for &component in components {
            buffer.extend_from_slice(&component.to_le_bytes());
        }
        xxh3_64_with_seed(&buffer, self.seed)
    }

    fn combine(&self, round: u32, current: u64, neighbors: &[(u32, u64)]) -> u64 {
        let mut buffer = Vec::with_capacity(16 + neighbors.len() * 12);
        buffer.extend_from_slice(&u64::from(round).to_le_bytes());
        buffer.extend_from_slice(&current.to_le_bytes());
        for &(edge, color) in neighbors {
            buffer.extend_from_slice(&edge.to_le_bytes());
            buffer.extend_from_slice(&color.to_le_bytes());
        }
        xxh3_64_with_seed(&buffer, self.seed)
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

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum WlHashScheme {
    #[default]
    Xxh3SortedWidth64V1,
}

impl WlHashScheme {
    pub const fn version(self) -> u16 {
        match self {
            Self::Xxh3SortedWidth64V1 => XXH3_SCHEME_VERSION_1,
        }
    }

    pub const fn identifier_width(self) -> u16 {
        match self {
            Self::Xxh3SortedWidth64V1 => HASH_WIDTH_64,
        }
    }

    pub(crate) const fn seed(self) -> u64 {
        match self {
            Self::Xxh3SortedWidth64V1 => WL_XXH3_SORTED_V1_SEED,
        }
    }

    pub(crate) const fn aggregation(self) -> RefinementAggregation {
        match self {
            Self::Xxh3SortedWidth64V1 => RefinementAggregation::Sorted,
        }
    }

    pub(crate) fn refinement_scheme(self) -> RefinementXxh3Scheme<RefinementWidth64> {
        RefinementXxh3Scheme::new(self.seed(), self.aggregation())
    }
}

#[allow(non_camel_case_types)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum EcfpHashScheme {
    #[default]
    Xxh3Width64V1,
}

impl EcfpHashScheme {
    pub const fn version(self) -> u16 {
        match self {
            Self::Xxh3Width64V1 => XXH3_SCHEME_VERSION_1,
        }
    }

    pub const fn identifier_width(self) -> u16 {
        match self {
            Self::Xxh3Width64V1 => HASH_WIDTH_64,
        }
    }

    pub(crate) const fn seed(self) -> u64 {
        match self {
            Self::Xxh3Width64V1 => ECFP_XXH3_64_V1_SEED,
        }
    }

    pub(crate) const fn recipe(self) -> RogersHahn {
        match self {
            Self::Xxh3Width64V1 => RogersHahn { seed: self.seed() },
        }
    }
}

/// Plain-refinement recipe over xxh3: `seed` is the family selector, `W` the width,
/// `aggregation` the multiset combine.
#[derive(Clone, Copy, Debug)]
pub struct RefinementXxh3Scheme<W: RefinementWidth> {
    seed: u64,
    aggregation: RefinementAggregation,
    _width: PhantomData<W>,
}

impl<W: RefinementWidth> RefinementXxh3Scheme<W> {
    pub fn new(seed: u64, aggregation: RefinementAggregation) -> Self {
        Self {
            seed,
            aggregation,
            _width: PhantomData,
        }
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
                let mut elements: Vec<W::Color> = neighbors
                    .iter()
                    .map(|&(edge, c)| element(edge, c))
                    .collect();
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

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use pretty_assertions::{assert_eq, assert_ne};
    use rstest::*;
    use umol_graph_core::{
        CircularRefinementAlgorithm, Graph, RefinementAlgorithm, RefinementRounds,
    };

    use super::*;

    fn wl_128(
        rounds: RefinementRounds,
    ) -> RefinementAlgorithm<RefinementXxh3Scheme<RefinementWidth128>> {
        RefinementAlgorithm::WeisfeilerLehman {
            rounds,
            scheme: RefinementXxh3Scheme::new(
                WL_XXH3_SORTED_V1_SEED,
                RefinementAggregation::Sorted,
            ),
        }
    }

    fn hash_of(g: &Graph) -> u128 {
        g.refine(|_| 0, |_| 0, wl_128(RefinementRounds::ToFixpoint))
            .graph_hash()
    }

    #[rstest]
    #[case::triangle(Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), 1)]
    #[case::path3(Graph::new(3, &[[0, 1], [1, 2]]), 2)]
    #[case::path4(Graph::new(4, &[[0, 1], [1, 2], [2, 3]]), 2)]
    #[case::star4(Graph::new(5, &[[0, 1], [0, 2], [0, 3], [0, 4]]), 2)]
    #[case::six_cycle(Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]), 1)]
    fn test_graph_refine(#[case] g: Graph, #[case] expected_cells: usize) {
        let refinement = g.refine(|_| 0, |_| 0, wl_128(RefinementRounds::ToFixpoint));
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
            .refine(
                |nd| nd.index() as u64,
                |_| 0,
                wl_128(RefinementRounds::ToFixpoint),
            )
            .graph_hash();
        let edge_colored = g
            .refine(|_| 0, |_| 7, wl_128(RefinementRounds::ToFixpoint))
            .graph_hash();
        assert_ne!(plain, node_colored);
        assert_ne!(plain, edge_colored);
    }

    #[rstest]
    #[case::xxh3_sorted_width64_v1(
        WlHashScheme::Xxh3SortedWidth64V1,
        XXH3_SCHEME_VERSION_1,
        HASH_WIDTH_64,
        WL_XXH3_SORTED_V1_SEED,
        RefinementAggregation::Sorted
    )]
    fn test_wl_hash_scheme(
        #[case] scheme: WlHashScheme,
        #[case] expected_version: u16,
        #[case] expected_width: u16,
        #[case] expected_seed: u64,
        #[case] expected_aggregation: RefinementAggregation,
    ) {
        assert_eq!(scheme.version(), expected_version);
        assert_eq!(scheme.identifier_width(), expected_width);
        assert_eq!(scheme.seed(), expected_seed);
        assert_eq!(scheme.aggregation(), expected_aggregation);
        assert_eq!(scheme.refinement_scheme().seed, expected_seed);
        assert_eq!(scheme.refinement_scheme().aggregation, expected_aggregation);
    }

    #[rstest]
    #[case::xxh3_width64_v1(
        EcfpHashScheme::Xxh3Width64V1,
        XXH3_SCHEME_VERSION_1,
        HASH_WIDTH_64,
        ECFP_XXH3_64_V1_SEED
    )]
    fn test_ecfp_hash_scheme(
        #[case] scheme: EcfpHashScheme,
        #[case] expected_version: u16,
        #[case] expected_width: u16,
        #[case] expected_seed: u64,
    ) {
        assert_eq!(scheme.version(), expected_version);
        assert_eq!(scheme.identifier_width(), expected_width);
        assert_eq!(scheme.seed(), expected_seed);
        assert_eq!(
            scheme.recipe(),
            RogersHahn {
                seed: expected_seed
            }
        );
    }

    #[rstest]
    fn test_refinement_counts() {
        // path 0-1-2, radius 1: round 0 all-same (3); round 1 splits ends (2) from middle (1).
        let g = Graph::new(3, &[[0, 1], [1, 2]]);
        let refinement = g.refine(|_| 0, |_| 0, wl_128(RefinementRounds::Fixed(1)));
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
                WL_XXH3_SORTED_V1_SEED,
                RefinementAggregation::SumSketch,
            ),
        };
        let refinement = g.refine(|_| 0, |_| 0, algorithm);
        assert_eq!(refinement.cell_count(), 2);
    }

    // The 64- and 128-bit widths both distinguish path from triangle.
    fn distinguishes<W: RefinementWidth>() {
        let algorithm = || RefinementAlgorithm::WeisfeilerLehman {
            rounds: RefinementRounds::ToFixpoint,
            scheme: RefinementXxh3Scheme::<W>::new(
                WL_XXH3_SORTED_V1_SEED,
                RefinementAggregation::Sorted,
            ),
        };
        let path = Graph::new(3, &[[0, 1], [1, 2]]);
        let triangle = Graph::new(3, &[[0, 1], [1, 2], [0, 2]]);
        let hp = path.refine(|_| 0, |_| 0, algorithm()).graph_hash();
        let ht = triangle.refine(|_| 0, |_| 0, algorithm()).graph_hash();
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

    #[rstest]
    fn test_refinement_xxh3_scheme_graph_hash() {
        let g = Graph::new(3, &[[0, 1], [1, 2]]);
        let algorithm = RefinementAlgorithm::WeisfeilerLehman {
            rounds: RefinementRounds::ToFixpoint,
            scheme: RefinementXxh3Scheme::<RefinementWidth128>::new(
                WL_XXH3_SORTED_V1_SEED,
                RefinementAggregation::Sorted,
            ),
        };
        let h = g.refine(|_| 0, |_| 0, algorithm).graph_hash();
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
            CircularRefinementAlgorithm::Ec {
                radius: 1,
                scheme: Morgan,
            },
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
}
