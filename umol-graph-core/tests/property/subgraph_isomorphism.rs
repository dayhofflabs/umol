//! Differential, soundness, completeness, symmetry, and delivery-agreement
//! properties for subgraph-isomorphism search.
//!
//! Validation spans four distinct operational domains and methods:
//!
//! 1. **Bounded random and planted (possibly non-induced) labeled graphs (`n <= 6`)**:
//!    compared directly against definition-level exhaustive monomorphism and
//!    automorphism-group orbit enumeration (`exhaustive`), establishing
//!    soundness, completeness, delivery agreement (`enumerate`, `visit`, `has`,
//!    `first`, `count`, early `ControlFlow::Break`), anchored `_at` restriction,
//!    predicate-hook equivalence, `TrailIr` == `Nauty` symmetry invariance, and
//!    exact canonical orbit representatives across all four `(pattern_sym, host_sym)`
//!    combinations.
//! 2. **Checked-in simple-graph corpus (`simple-through-8.g6`)**: exhaustively
//!    checks all 52 × 52 pairs through order 5 and all 18 × 156 pairs of
//!    order-1..=4 queries in order-6 targets across the four core engines, both
//!    automorphism algorithms, and all four symmetry modes against `exhaustive`.
//! 3. **Medium and wide generated regimes (`65..=128` and `129..=196`)**:
//!    sparse labeled graphs with repeated symmetric motifs and shuffled vertex
//!    IDs, checking cross-engine parity (`Vf2Rdkit`, `RayKirsch`, `BitsetRi`,
//!    `Glasgow`), `TrailIr` == `Nauty`, all four symmetry modes, pointwise vs.
//!    bitset-hook predicates, direct vs. prepared parity, anchored `_at`, and
//!    independent monomorphism validity on every emitted embedding.
//! 4. **Explicit boundary and lifecycle cases**: empty graphs, exact bitset
//!    word boundaries (`64`, `65`, `128`, `129`), `n1 > 16` queries, sparse wide
//!    targets above both eager (`1024`) and lazy (`4096`) adjacency cutoffs,
//!    proper subgroups via `try_from_generators`, and `SubgraphIsomorphismScratch`
//!    reuse across regimes and after early `Break`.

#[path = "subgraph_isomorphism/exhaustive.rs"]
mod exhaustive;

use std::array::from_fn;
use std::collections::HashMap;
use std::ops::ControlFlow;

use proptest::prelude::*;
use rstest::rstest;
use umol_graph_core::SubgraphIsomorphismAlgorithm::{BitsetRi, Glasgow, RayKirsch, Vf2Rdkit};
use umol_graph_core::{
    AutomorphismAlgorithm, AutomorphismScratch, EdgeId, Graph, NodeId, PreparedSubgraphQuery,
    PreparedSubgraphTarget, QuerySymmetryConstraints, SubgraphIsomorphismAlgorithm,
    SubgraphIsomorphismPredicate, SubgraphIsomorphismScratch, TargetSymmetry,
};

use self::exhaustive::{
    canonical_orbit_representatives, exhaustive_automorphisms, exhaustive_monomorphisms,
    is_valid_monomorphism, query_priority_order, subgroup_from_generators,
};
use super::corpus::simple_graphs;

const CORE_ALGORITHMS: [SubgraphIsomorphismAlgorithm; 4] = [Vf2Rdkit, RayKirsch, BitsetRi, Glasgow];

const AUTOMORPHISM_ALGORITHMS: [AutomorphismAlgorithm; 2] =
    [AutomorphismAlgorithm::TrailIr, AutomorphismAlgorithm::Nauty];

const SYMMETRY_MODES: [(bool, bool); 4] =
    [(false, false), (true, false), (false, true), (true, true)];

const NUM_NODE_LABELS: usize = 3;
const NUM_EDGE_LABELS: usize = 3;
const ALL_EDGE_MASK: u8 = (1u8 << NUM_EDGE_LABELS) - 1;

#[derive(Clone, Debug)]
struct LabeledGraph {
    graph: Graph,
    node_labels: Vec<u8>,
    edge_labels: Vec<u8>,
}

#[derive(Clone, Debug)]
struct MaskedQuery {
    graph: Graph,
    node_masks: Vec<u8>,
    edge_masks: Vec<u8>,
}

#[derive(Clone, Copy, Debug)]
enum EdgeHookMode {
    AllSome,
    Mixed,
    AllNone,
}

/// Pointwise-only mask predicate using default trait hooks.
struct PointwiseMaskPredicate<'a> {
    query_node_masks: &'a [u8],
    query_edge_masks: &'a [u8],
    target_node_labels: &'a [u8],
    target_edge_labels: &'a [u8],
}

impl SubgraphIsomorphismPredicate for PointwiseMaskPredicate<'_> {
    fn match_node(&mut self, query_node: NodeId, target_node: NodeId) -> bool {
        let label = self.target_node_labels[target_node.index()];
        (self.query_node_masks[query_node.index()] & (1u8 << label)) != 0
    }

    fn match_edge(&mut self, query_edge: EdgeId, target_edge: EdgeId) -> bool {
        let label = self.target_edge_labels[target_edge.index()];
        (self.query_edge_masks[query_edge.index()] & (1u8 << label)) != 0
    }
}

/// Optimized mask predicate exercising undetermined-edge hooks, query-node
/// equivalence classes, bulk bitset candidate fills, and edge-neighbor bitmasks.
struct BitsetMaskPredicate<'a> {
    query_node_masks: &'a [u8],
    query_edge_masks: &'a [u8],
    target_node_labels: &'a [u8],
    target_edge_labels: &'a [u8],
    all_edges_undetermined: bool,
    class_count: usize,
    node_class: Vec<usize>,
    node_masks_u64: [u64; NUM_NODE_LABELS],
    node_masks_u128: [u128; NUM_NODE_LABELS],
    node_masks_words: [Vec<u64>; NUM_NODE_LABELS],
    edge_adj_u64: Vec<[u64; NUM_EDGE_LABELS]>,
    edge_adj_u128: Vec<[u128; NUM_EDGE_LABELS]>,
    edge_hook_mode: EdgeHookMode,
}

impl<'a> BitsetMaskPredicate<'a> {
    fn new(query: &'a MaskedQuery, target: &'a LabeledGraph, edge_hook_mode: EdgeHookMode) -> Self {
        let mut class_masks: Vec<u8> = Vec::new();
        let mut node_class = Vec::with_capacity(query.node_masks.len());
        for &mask in &query.node_masks {
            let cls = match class_masks.iter().position(|&existing| existing == mask) {
                Some(idx) => idx,
                None => {
                    class_masks.push(mask);
                    class_masks.len() - 1
                }
            };
            node_class.push(cls);
        }

        let n2 = target.graph.node_count();
        let words = n2.div_ceil(64);
        let mut node_masks_u64 = [0u64; NUM_NODE_LABELS];
        let mut node_masks_u128 = [0u128; NUM_NODE_LABELS];
        let mut node_masks_words: [Vec<u64>; NUM_NODE_LABELS] = from_fn(|_| vec![0u64; words]);

        for (t, &label) in target.node_labels.iter().enumerate() {
            let lbl = label as usize;
            if t < 64 {
                node_masks_u64[lbl] |= 1u64 << t;
            }
            if t < 128 {
                node_masks_u128[lbl] |= 1u128 << t;
            }
            node_masks_words[lbl][t / 64] |= 1u64 << (t % 64);
        }

        let mut edge_adj_u64 = Vec::new();
        let mut edge_adj_u128 = Vec::new();
        if n2 <= 64 {
            edge_adj_u64 = vec![[0u64; NUM_EDGE_LABELS]; n2];
            for te in target.graph.edge_ids() {
                let [u, v] = target.graph.edge_endpoints(te);
                let lbl = target.edge_labels[te.index()] as usize;
                edge_adj_u64[u.index()][lbl] |= 1u64 << v.index();
                edge_adj_u64[v.index()][lbl] |= 1u64 << u.index();
            }
        }
        if n2 <= 128 {
            edge_adj_u128 = vec![[0u128; NUM_EDGE_LABELS]; n2];
            for te in target.graph.edge_ids() {
                let [u, v] = target.graph.edge_endpoints(te);
                let lbl = target.edge_labels[te.index()] as usize;
                edge_adj_u128[u.index()][lbl] |= 1u128 << v.index();
                edge_adj_u128[v.index()][lbl] |= 1u128 << u.index();
            }
        }

        let all_edges_undetermined = query.edge_masks.iter().all(|&m| m == ALL_EDGE_MASK);
        Self {
            query_node_masks: &query.node_masks,
            query_edge_masks: &query.edge_masks,
            target_node_labels: &target.node_labels,
            target_edge_labels: &target.edge_labels,
            all_edges_undetermined,
            class_count: class_masks.len(),
            node_class,
            node_masks_u64,
            node_masks_u128,
            node_masks_words,
            edge_adj_u64,
            edge_adj_u128,
            edge_hook_mode,
        }
    }
}

impl SubgraphIsomorphismPredicate for BitsetMaskPredicate<'_> {
    fn match_node(&mut self, query_node: NodeId, target_node: NodeId) -> bool {
        let label = self.target_node_labels[target_node.index()];
        (self.query_node_masks[query_node.index()] & (1u8 << label)) != 0
    }

    fn match_edge(&mut self, query_edge: EdgeId, target_edge: EdgeId) -> bool {
        let label = self.target_edge_labels[target_edge.index()];
        (self.query_edge_masks[query_edge.index()] & (1u8 << label)) != 0
    }

    fn all_edges_undetermined(&self) -> bool {
        self.all_edges_undetermined
    }

    fn edge_is_undetermined(&self, query_edge: EdgeId) -> bool {
        self.query_edge_masks[query_edge.index()] == ALL_EDGE_MASK
    }

    fn query_node_class_count(&self) -> Option<usize> {
        Some(self.class_count)
    }

    fn query_node_class(&self, query_node: NodeId) -> usize {
        self.node_class[query_node.index()]
    }

    fn node_candidates_u64(&mut self, query_node: NodeId, eligible: u64) -> u64 {
        let mask = self.query_node_masks[query_node.index()];
        let mut combined = 0u64;
        for (lbl, &bits) in self.node_masks_u64.iter().enumerate() {
            if (mask & (1u8 << lbl)) != 0 {
                combined |= bits;
            }
        }
        combined & eligible
    }

    fn node_candidates_u128(&mut self, query_node: NodeId, eligible: u128) -> u128 {
        let mask = self.query_node_masks[query_node.index()];
        let mut combined = 0u128;
        for (lbl, &bits) in self.node_masks_u128.iter().enumerate() {
            if (mask & (1u8 << lbl)) != 0 {
                combined |= bits;
            }
        }
        combined & eligible
    }

    fn fill_node_candidates_words(
        &mut self,
        query_node: NodeId,
        eligible: &[u64],
        out: &mut [u64],
    ) {
        let mask = self.query_node_masks[query_node.index()];
        let len = eligible.len().min(out.len());
        for word_idx in 0..len {
            let mut combined = 0u64;
            for (lbl, words) in self.node_masks_words.iter().enumerate() {
                if (mask & (1u8 << lbl)) != 0 {
                    combined |= words[word_idx];
                }
            }
            out[word_idx] = combined & eligible[word_idx];
        }
        out[len..].fill(0);
    }

    fn edge_neighbors_u64(&mut self, query_edge: EdgeId, target_node: NodeId) -> Option<u64> {
        let mask = self.query_edge_masks[query_edge.index()];
        match self.edge_hook_mode {
            EdgeHookMode::AllNone => return None,
            EdgeHookMode::Mixed if mask.count_ones() != 1 => return None,
            _ => {}
        }
        let mut combined = 0u64;
        for (lbl, &bits) in self.edge_adj_u64[target_node.index()].iter().enumerate() {
            if (mask & (1u8 << lbl)) != 0 {
                combined |= bits;
            }
        }
        Some(combined)
    }

    fn edge_neighbors_u128(&mut self, query_edge: EdgeId, target_node: NodeId) -> Option<u128> {
        let mask = self.query_edge_masks[query_edge.index()];
        match self.edge_hook_mode {
            EdgeHookMode::AllNone => return None,
            EdgeHookMode::Mixed if mask.count_ones() != 1 => return None,
            _ => {}
        }
        let mut combined = 0u128;
        for (lbl, &bits) in self.edge_adj_u128[target_node.index()].iter().enumerate() {
            if (mask & (1u8 << lbl)) != 0 {
                combined |= bits;
            }
        }
        Some(combined)
    }
}

/// Unlabeled predicate admitting every node and edge, reporting all edges
/// undetermined.
struct UnlabeledPredicate;

impl SubgraphIsomorphismPredicate for UnlabeledPredicate {
    fn match_node(&mut self, _query_node: NodeId, _target_node: NodeId) -> bool {
        true
    }

    fn match_edge(&mut self, _query_edge: EdgeId, _target_edge: EdgeId) -> bool {
        true
    }

    fn all_edges_undetermined(&self) -> bool {
        true
    }

    fn query_node_class_count(&self) -> Option<usize> {
        Some(1)
    }

    fn query_node_class(&self, _query_node: NodeId) -> usize {
        0
    }

    fn node_candidates_u64(&mut self, _query_node: NodeId, eligible: u64) -> u64 {
        eligible
    }

    fn node_candidates_u128(&mut self, _query_node: NodeId, eligible: u128) -> u128 {
        eligible
    }

    fn fill_node_candidates_words(
        &mut self,
        _query_node: NodeId,
        eligible: &[u64],
        out: &mut [u64],
    ) {
        let len = eligible.len().min(out.len());
        out[..len].copy_from_slice(&eligible[..len]);
        out[len..].fill(0);
    }
}

fn labeled_graph(
    min_n: usize,
    max_n: usize,
    node_labels: u8,
    edge_labels: u8,
) -> impl Strategy<Value = LabeledGraph> {
    (min_n..=max_n)
        .prop_flat_map(move |n| {
            let pairs: Vec<(u32, u32)> = (0..n as u32)
                .flat_map(|i| ((i + 1)..n as u32).map(move |j| (i, j)))
                .collect();
            let mask = prop::collection::vec(
                prop_oneof![Just(None::<u8>), (0..edge_labels).prop_map(Some)],
                pairs.len(),
            );
            let labels = prop::collection::vec(0..node_labels, n);
            (Just(n), Just(pairs), labels, mask)
        })
        .prop_map(|(n, pairs, node_labels, mask)| {
            let mut edges = Vec::new();
            let mut edge_labels = Vec::new();
            for (pair, maybe_label) in pairs.iter().zip(mask) {
                if let Some(label) = maybe_label {
                    edges.push([pair.0, pair.1]);
                    edge_labels.push(label);
                }
            }
            LabeledGraph {
                graph: Graph::new(n, &edges),
                node_labels,
                edge_labels,
            }
        })
}

fn masked_query(
    min_n: usize,
    max_n: usize,
    node_label_count: usize,
    edge_label_count: usize,
) -> impl Strategy<Value = MaskedQuery> {
    let max_node_mask = (1u8 << node_label_count) - 1;
    let max_edge_mask = (1u8 << edge_label_count) - 1;
    (min_n..=max_n)
        .prop_flat_map(move |n| {
            let pairs: Vec<(u32, u32)> = (0..n as u32)
                .flat_map(|i| ((i + 1)..n as u32).map(move |j| (i, j)))
                .collect();
            let edges = prop::collection::vec(
                prop_oneof![Just(None::<u8>), (1..=max_edge_mask).prop_map(Some)],
                pairs.len(),
            );
            let node_masks = prop::collection::vec(1..=max_node_mask, n);
            (Just(n), Just(pairs), node_masks, edges)
        })
        .prop_map(|(n, pairs, node_masks, mask)| {
            let mut edges = Vec::new();
            let mut edge_masks = Vec::new();
            for (pair, maybe_mask) in pairs.iter().zip(mask) {
                if let Some(edge_mask) = maybe_mask {
                    edges.push([pair.0, pair.1]);
                    edge_masks.push(edge_mask);
                }
            }
            MaskedQuery {
                graph: Graph::new(n, &edges),
                node_masks,
                edge_masks,
            }
        })
}

/// Extracts a possibly non-induced subgraph on `nodes` (with `keep_edge_bits`
/// choosing which induced edges survive) and relabels vertices by `order`.
fn planted_subgraph(
    target: &LabeledGraph,
    nodes: &[usize],
    keep_edge_bits: u64,
    extra_node_masks: &[u8],
    extra_edge_masks: &[u8],
) -> (LabeledGraph, MaskedQuery) {
    let remap: HashMap<usize, u32> = nodes
        .iter()
        .enumerate()
        .map(|(new_id, &old_id)| (old_id, new_id as u32))
        .collect();
    let mut edges = Vec::new();
    let mut edge_labels = Vec::new();
    let mut edge_masks = Vec::new();
    let mut induced_idx = 0usize;
    for eid in 0..target.graph.edge_count() {
        let [a, b] = target.graph.edge_endpoints(EdgeId(eid as u32));
        if let (Some(&na), Some(&nb)) = (remap.get(&a.index()), remap.get(&b.index())) {
            let keep = induced_idx >= 64 || ((keep_edge_bits >> induced_idx) & 1) == 1;
            if keep {
                let lbl = target.edge_labels[eid];
                let extra = extra_edge_masks.get(edges.len()).copied().unwrap_or(0);
                edges.push([na, nb]);
                edge_labels.push(lbl);
                edge_masks.push((1u8 << lbl) | (extra & ALL_EDGE_MASK));
            }
            induced_idx += 1;
        }
    }
    let node_labels: Vec<u8> = nodes.iter().map(|&old| target.node_labels[old]).collect();
    let all_node_mask = (1u8 << NUM_NODE_LABELS) - 1;
    let node_masks: Vec<u8> = node_labels
        .iter()
        .enumerate()
        .map(|(i, &lbl)| {
            let extra = extra_node_masks.get(i).copied().unwrap_or(0);
            (1u8 << lbl) | (extra & all_node_mask)
        })
        .collect();
    let graph = Graph::new(nodes.len(), &edges);
    (
        LabeledGraph {
            graph: graph.clone(),
            node_labels,
            edge_labels,
        },
        MaskedQuery {
            graph,
            node_masks,
            edge_masks,
        },
    )
}

/// Sparse labeled regime case with target size in `min_n..=max_n`, repeated
/// symmetric motif copies (yielding non-trivial `G_H`), shuffled vertex IDs,
/// and a planted or deliberately unsatisfied masked query.
#[derive(Clone, Debug)]
struct SparseRegimeCase {
    target: LabeledGraph,
    query: MaskedQuery,
    expect_match: bool,
    hook_mode: EdgeHookMode,
    anchor: (NodeId, NodeId),
}

fn sparse_regime_case(min_n: usize, max_n: usize) -> impl Strategy<Value = SparseRegimeCase> {
    (
        min_n..=max_n,
        0u8..5,
        2usize..=3,
        any::<bool>(),
        any::<bool>(),
        0u8..3,
        0u8..=ALL_EDGE_MASK,
        any::<u64>(),
    )
        .prop_map(
            |(n2, motif_kind, copies, make_unsat, two_components, hook_idx, extra_edge, seed)| {
                let mut rng_state = seed.wrapping_add(0x9e37_79b9_7f4a_7c15);
                let mut next_u64 = || {
                    rng_state ^= rng_state << 13;
                    rng_state ^= rng_state >> 7;
                    rng_state ^= rng_state << 17;
                    rng_state
                };

                // Motif topology + symmetric coloring. Motif 4 intentionally places a
                // rare label (2) on low-degree leaves while the background uses only
                // labels 0..=1 with degree-3 hubs so BitsetRi RI-DS and Glasgow MRV
                // branch on a non-priority vertex first while G_Q and G_H are active.
                let rare_leaf_motif = motif_kind == 4;
                let (motif_n, motif_edges, motif_node_labels, motif_edge_labels): (
                    usize,
                    Vec<[u32; 2]>,
                    Vec<u8>,
                    Vec<u8>,
                ) = match motif_kind {
                    0 => (
                        // 4-cycle with symmetric alternating node colors
                        4,
                        vec![[0, 1], [1, 2], [2, 3], [3, 0]],
                        vec![0, 1, 0, 1],
                        vec![0, 0, 0, 0],
                    ),
                    1 => (
                        // Diamond (K_4 - e) with symmetric apexes and waist
                        4,
                        vec![[0, 1], [1, 2], [2, 3], [3, 0], [0, 2]],
                        vec![0, 1, 0, 1],
                        vec![0, 0, 0, 0, 1],
                    ),
                    2 => (
                        // 5-node bowtie/branched symmetric motif
                        5,
                        vec![[0, 1], [0, 2], [1, 2], [0, 3], [0, 4]],
                        vec![0, 1, 1, 2, 2],
                        vec![0, 0, 1, 2, 2],
                    ),
                    3 => (
                        // 6-cycle with C_3 symmetry
                        6,
                        vec![[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]],
                        vec![0, 1, 0, 1, 0, 1],
                        vec![0, 1, 0, 1, 0, 1],
                    ),
                    _ => (
                        // Degree-3 common hub (label 0) + triangle + two symmetric
                        // degree-1 rare leaves (label 2).
                        5,
                        vec![[0, 1], [0, 2], [1, 2], [1, 3], [2, 4]],
                        vec![0, 0, 0, 2, 2],
                        vec![0, 0, 1, 0, 0],
                    ),
                };

                // Build raw target: `copies` identical motif components + asymmetric
                // sparse background path with a few chords.
                let mut raw_edges: Vec<[u32; 2]> = Vec::new();
                let mut raw_edge_labels: Vec<u8> = Vec::new();
                let mut raw_node_labels: Vec<u8> = Vec::with_capacity(n2);

                let motif_total = copies * motif_n;
                for copy in 0..copies {
                    let base = (copy * motif_n) as u32;
                    raw_node_labels.extend_from_slice(&motif_node_labels);
                    for (idx, &[u, v]) in motif_edges.iter().enumerate() {
                        raw_edges.push([base + u, base + v]);
                        raw_edge_labels.push(motif_edge_labels[idx]);
                    }
                }

                // Background nodes carry distinct color patterns along a sparse path
                // so `Aut(H)` stays fast even at n2 = 200 while acting non-trivially
                // on the repeated symmetric motifs.
                let bg_start = motif_total;
                let bg_label_mod = if rare_leaf_motif { 2 } else { NUM_NODE_LABELS };
                for i in bg_start..n2 {
                    let offset = i - bg_start;
                    raw_node_labels.push(((offset * 5 + (offset / 7)) % bg_label_mod) as u8);
                    if i > bg_start {
                        raw_edges.push([(i - 1) as u32, i as u32]);
                        raw_edge_labels.push(((offset + 1) % NUM_EDGE_LABELS) as u8);
                    }
                    if offset >= 3 && offset.is_multiple_of(5) {
                        raw_edges.push([(i - 3) as u32, i as u32]);
                        raw_edge_labels.push(((offset + 2) % NUM_EDGE_LABELS) as u8);
                    }
                }

                // Shuffle vertex IDs so motif nodes span all bitset words.
                let mut perm: Vec<usize> = (0..n2).collect();
                for i in (1..n2).rev() {
                    let j = (next_u64() as usize) % (i + 1);
                    perm.swap(i, j);
                }
                let mut inv_perm = vec![0usize; n2];
                for (new_id, &old_id) in perm.iter().enumerate() {
                    inv_perm[old_id] = new_id;
                }

                let mut target_node_labels = vec![0u8; n2];
                for (old_id, &lbl) in raw_node_labels.iter().enumerate() {
                    target_node_labels[inv_perm[old_id]] = lbl;
                }
                let target_edges: Vec<[u32; 2]> = raw_edges
                    .iter()
                    .map(|&[u, v]| [inv_perm[u as usize] as u32, inv_perm[v as usize] as u32])
                    .collect();
                let target = LabeledGraph {
                    graph: Graph::new(n2, &target_edges),
                    node_labels: target_node_labels,
                    edge_labels: raw_edge_labels,
                };

                // Build a masked query from the motif (either 1 component or 2
                // identical symmetric components; optionally drop one edge for
                // non-induced monomorphism, or add a chord/triangle for UNSAT).
                let mut q_edges: Vec<[u32; 2]> = Vec::new();
                let mut q_node_masks: Vec<u8> = Vec::new();
                let mut q_edge_masks: Vec<u8> = Vec::new();

                let q_copies = if two_components { 2 } else { 1 };
                for c in 0..q_copies {
                    let base = (c * motif_n) as u32;
                    for &lbl in &motif_node_labels {
                        q_node_masks.push(1u8 << lbl);
                    }
                    let drop_last = !make_unsat && !two_components && motif_edges.len() > 4;
                    let keep_len = if drop_last {
                        motif_edges.len() - 1
                    } else {
                        motif_edges.len()
                    };
                    for idx in 0..keep_len {
                        let [u, v] = motif_edges[idx];
                        q_edges.push([base + u, base + v]);
                        let base_mask = 1u8 << motif_edge_labels[idx];
                        // Optionally widen all edges equally to preserve motif symmetry.
                        let mask = if idx.is_multiple_of(2) {
                            base_mask | (extra_edge & ALL_EDGE_MASK)
                        } else {
                            base_mask
                        };
                        q_edge_masks.push(mask);
                    }
                }

                let expect_match = if make_unsat {
                    // Add an impossible clique of size 4 with all edge label 2 so
                    // no motif or sparse background chord contains it.
                    let base = q_node_masks.len() as u32;
                    q_node_masks.extend_from_slice(&[1u8, 1u8, 1u8, 1u8]);
                    for a in 0..4u32 {
                        for b in (a + 1)..4u32 {
                            q_edges.push([base + a, base + b]);
                            q_edge_masks.push(1u8 << 2);
                        }
                    }
                    false
                } else {
                    true
                };

                let query = MaskedQuery {
                    graph: Graph::new(q_node_masks.len(), &q_edges),
                    node_masks: q_node_masks,
                    edge_masks: q_edge_masks,
                };
                let hook_mode = match hook_idx {
                    0 => EdgeHookMode::AllSome,
                    1 => EdgeHookMode::Mixed,
                    _ => EdgeHookMode::AllNone,
                };
                // Anchor query node 0 to the shuffled image of motif copy 0 node 0.
                let anchor = (NodeId(0), NodeId(inv_perm[0] as u32));

                SparseRegimeCase {
                    target,
                    query,
                    expect_match,
                    hook_mode,
                    anchor,
                }
            },
        )
}

fn collect_direct_enumerate(
    query: &LabeledGraph,
    target: &LabeledGraph,
    alg: SubgraphIsomorphismAlgorithm,
) -> Vec<Vec<usize>> {
    let mut node_match =
        |q: NodeId, t: NodeId| query.node_labels[q.index()] == target.node_labels[t.index()];
    let mut edge_match =
        |qe: EdgeId, te: EdgeId| query.edge_labels[qe.index()] == target.edge_labels[te.index()];
    let mut found: Vec<Vec<usize>> = target
        .graph
        .enumerate_subgraph_isomorphisms(&query.graph, &mut node_match, &mut edge_match, alg)
        .iter()
        .map(|c| c.matched_pairs().iter().map(|&(_, t)| t.index()).collect())
        .collect();
    found.sort();
    found
}

fn collect_direct_visit_in_order(
    query: &LabeledGraph,
    target: &LabeledGraph,
    alg: SubgraphIsomorphismAlgorithm,
) -> Vec<Vec<usize>> {
    let mut node_match =
        |q: NodeId, t: NodeId| query.node_labels[q.index()] == target.node_labels[t.index()];
    let mut edge_match =
        |qe: EdgeId, te: EdgeId| query.edge_labels[qe.index()] == target.edge_labels[te.index()];
    let mut found = Vec::new();
    let result: ControlFlow<()> = target.graph.visit_subgraph_isomorphisms(
        &query.graph,
        &mut node_match,
        &mut edge_match,
        alg,
        |embedding| {
            found.push(embedding.iter().map(|n| n.index()).collect());
            ControlFlow::Continue(())
        },
    );
    assert_eq!(result, ControlFlow::Continue(()));
    found
}

fn verify_direct_modes_against_expected(
    query: &LabeledGraph,
    target: &LabeledGraph,
    alg: SubgraphIsomorphismAlgorithm,
    expected: &[Vec<usize>],
) {
    let enumerated = collect_direct_enumerate(query, target, alg);
    assert_eq!(
        enumerated, expected,
        "{alg:?} direct enumerate disagrees with expected"
    );

    let visited_in_order = collect_direct_visit_in_order(query, target, alg);
    let mut visited_sorted = visited_in_order.clone();
    visited_sorted.sort();
    assert_eq!(
        visited_sorted, expected,
        "{alg:?} direct visit disagrees with expected"
    );

    let mut nm =
        |q: NodeId, t: NodeId| query.node_labels[q.index()] == target.node_labels[t.index()];
    let mut em =
        |qe: EdgeId, te: EdgeId| query.edge_labels[qe.index()] == target.edge_labels[te.index()];
    assert_eq!(
        target
            .graph
            .has_subgraph_isomorphism(&query.graph, &mut nm, &mut em, alg),
        !expected.is_empty(),
        "{alg:?} direct has_subgraph_isomorphism disagrees"
    );
    assert_eq!(
        target
            .graph
            .count_subgraph_isomorphisms(&query.graph, &mut nm, &mut em, alg),
        expected.len() as u64,
        "{alg:?} direct count_subgraph_isomorphisms disagrees"
    );
    let first = target
        .graph
        .first_subgraph_isomorphism(&query.graph, &mut nm, &mut em, alg);
    match first {
        None => assert!(
            expected.is_empty(),
            "{alg:?} direct first_subgraph_isomorphism returned None on non-empty match set"
        ),
        Some(corr) => {
            let first_vec: Vec<usize> = corr
                .matched_pairs()
                .iter()
                .map(|&(_, t)| t.index())
                .collect();
            assert_eq!(
                first_vec, visited_in_order[0],
                "{alg:?} direct first_subgraph_isomorphism disagrees with first visited match"
            );
        }
    }

    if visited_in_order.len() >= 2 {
        let mut early = Vec::new();
        let flow = target.graph.visit_subgraph_isomorphisms(
            &query.graph,
            &mut nm,
            &mut em,
            alg,
            |embedding| {
                early.push(embedding.iter().map(|n| n.index()).collect::<Vec<_>>());
                if early.len() == 2 {
                    ControlFlow::Break(2usize)
                } else {
                    ControlFlow::Continue(())
                }
            },
        );
        assert_eq!(flow, ControlFlow::Break(2));
        assert_eq!(early, &visited_in_order[..2]);
    }
}

fn verify_prepared_modes_against_expected<P: SubgraphIsomorphismPredicate>(
    prepared_query: &PreparedSubgraphQuery,
    prepared_target: &PreparedSubgraphTarget,
    target_sym: Option<&TargetSymmetry>,
    predicate: &mut P,
    scratch: &mut SubgraphIsomorphismScratch,
    expected: &[Vec<usize>],
    context: &str,
) {
    let mut visited_in_order: Vec<Vec<usize>> = Vec::new();
    let flow: ControlFlow<()> =
        prepared_query.visit_matches(prepared_target, target_sym, predicate, scratch, |emb| {
            visited_in_order.push(emb.iter().map(|n| n.index()).collect());
            ControlFlow::Continue(())
        });
    assert_eq!(flow, ControlFlow::Continue(()), "{context}: visit_matches");

    let mut visited_sorted = visited_in_order.clone();
    visited_sorted.sort();
    assert_eq!(
        visited_sorted, expected,
        "{context}: visit_matches sorted set mismatch"
    );

    let count = prepared_query.count_matches(prepared_target, target_sym, predicate, scratch);
    assert_eq!(
        count,
        expected.len() as u64,
        "{context}: count_matches mismatch"
    );

    let has = prepared_query.has_match(prepared_target, target_sym, predicate, scratch);
    assert_eq!(has, !expected.is_empty(), "{context}: has_match mismatch");

    let first = prepared_query.first_match(prepared_target, target_sym, predicate, scratch);
    match first {
        None => assert!(
            expected.is_empty(),
            "{context}: first_match returned None on non-empty expected set"
        ),
        Some(corr) => {
            let first_vec: Vec<usize> = corr
                .matched_pairs()
                .iter()
                .map(|&(_, t)| t.index())
                .collect();
            assert_eq!(
                first_vec, visited_in_order[0],
                "{context}: first_match disagrees with first visited match"
            );
        }
    }

    if visited_in_order.len() >= 2 {
        let mut early: Vec<Vec<usize>> = Vec::new();
        let stop =
            prepared_query.visit_matches(prepared_target, target_sym, predicate, scratch, |emb| {
                early.push(emb.iter().map(|n| n.index()).collect());
                if early.len() == 2 {
                    ControlFlow::Break(2usize)
                } else {
                    ControlFlow::Continue(())
                }
            });
        assert_eq!(stop, ControlFlow::Break(2), "{context}: early break flow");
        assert_eq!(
            early,
            &visited_in_order[..2],
            "{context}: early break prefix"
        );
    }
}

fn target_orbits(sym: &TargetSymmetry) -> Vec<NodeId> {
    (0..sym.node_count())
        .map(|i| sym.orbit_of(NodeId(i as u32)))
        .collect()
}

/// Checks full cross-engine parity, `TrailIr` == `Nauty`, all four symmetry
/// modes, pointwise vs. bitset-hook predicates, direct parity, anchored `_at`,
/// and independent monomorphism validity on a `(query, target)` case.
fn verify_sparse_regime_case(case: &SparseRegimeCase) {
    let mut auto_scratch = AutomorphismScratch::new();
    let q_syms: Vec<QuerySymmetryConstraints> = AUTOMORPHISM_ALGORITHMS
        .iter()
        .map(|&auto_alg| {
            QuerySymmetryConstraints::from_graph_with_scratch(
                &case.query.graph,
                |q| case.query.node_masks[q.index()],
                |qe| case.query.edge_masks[qe.index()],
                auto_alg,
                &mut auto_scratch,
            )
        })
        .collect();
    assert_eq!(
        q_syms[0].constraints(),
        q_syms[1].constraints(),
        "TrailIr and Nauty produced different QuerySymmetryConstraints"
    );

    let t_syms: Vec<TargetSymmetry> = AUTOMORPHISM_ALGORITHMS
        .iter()
        .map(|&auto_alg| {
            TargetSymmetry::from_graph_with_scratch(
                &case.target.graph,
                |t| case.target.node_labels[t.index()],
                |te| case.target.edge_labels[te.index()],
                auto_alg,
                &mut auto_scratch,
            )
        })
        .collect();
    assert_eq!(
        target_orbits(&t_syms[0]),
        target_orbits(&t_syms[1]),
        "TrailIr and Nauty produced different TargetSymmetry orbits"
    );
    assert_eq!(
        t_syms[0].orbit_count(),
        t_syms[1].orbit_count(),
        "TrailIr and Nauty produced different TargetSymmetry orbit counts"
    );

    let prepared_target = PreparedSubgraphTarget::new(&case.target.graph);
    let mut scratch = SubgraphIsomorphismScratch::new();

    for &(pattern_sym, host_sym) in &SYMMETRY_MODES {
        let mut baseline: Option<Vec<Vec<usize>>> = None;

        for (auto_idx, &auto_alg) in AUTOMORPHISM_ALGORITHMS.iter().enumerate() {
            let target_sym = host_sym.then_some(&t_syms[auto_idx]);

            for &alg in &CORE_ALGORITHMS {
                let mut prepared_query = PreparedSubgraphQuery::new(&case.query.graph, alg);
                if pattern_sym {
                    prepared_query = prepared_query.with_symmetry(q_syms[auto_idx].clone());
                }

                let mut bitset_pred =
                    BitsetMaskPredicate::new(&case.query, &case.target, case.hook_mode);

                let expected = match &baseline {
                    Some(prev) => prev.as_slice(),
                    None => {
                        let mut first_set: Vec<Vec<usize>> = Vec::new();
                        let _: ControlFlow<()> = prepared_query.visit_matches(
                            &prepared_target,
                            target_sym,
                            &mut bitset_pred,
                            &mut scratch,
                            |emb| {
                                first_set.push(emb.iter().map(|n| n.index()).collect());
                                ControlFlow::Continue(())
                            },
                        );
                        first_set.sort();
                        let unique_len = {
                            let mut deduped = first_set.clone();
                            deduped.dedup();
                            deduped.len()
                        };
                        assert_eq!(
                            first_set.len(),
                            unique_len,
                            "duplicate embedding emitted under {alg:?} ({pattern_sym}, {host_sym})"
                        );
                        assert_eq!(
                            !first_set.is_empty(),
                            case.expect_match,
                            "unexpected satisfiability under {alg:?} ({pattern_sym}, {host_sym})"
                        );
                        let mut nm = |q: NodeId, t: NodeId| {
                            (case.query.node_masks[q.index()]
                                & (1u8 << case.target.node_labels[t.index()]))
                                != 0
                        };
                        let mut em = |qe: EdgeId, te: EdgeId| {
                            (case.query.edge_masks[qe.index()]
                                & (1u8 << case.target.edge_labels[te.index()]))
                                != 0
                        };
                        for emb in &first_set {
                            assert!(
                                is_valid_monomorphism(
                                    &case.query.graph,
                                    &case.target.graph,
                                    emb,
                                    &mut nm,
                                    &mut em,
                                ),
                                "emitted embedding {emb:?} is not a valid monomorphism"
                            );
                        }
                        baseline = Some(first_set);
                        baseline.as_deref().unwrap()
                    }
                };

                let ctx = format!(
                    "n2={} alg={alg:?} auto={auto_alg:?} sym=({pattern_sym},{host_sym}) hook={:?}",
                    case.target.graph.node_count(),
                    case.hook_mode
                );
                verify_prepared_modes_against_expected(
                    &prepared_query,
                    &prepared_target,
                    target_sym,
                    &mut bitset_pred,
                    &mut scratch,
                    expected,
                    &ctx,
                );

                // Verify pointwise default-hook predicate parity on TrailIr pass.
                if auto_idx == 0 {
                    let mut pointwise_pred = PointwiseMaskPredicate {
                        query_node_masks: &case.query.node_masks,
                        query_edge_masks: &case.query.edge_masks,
                        target_node_labels: &case.target.node_labels,
                        target_edge_labels: &case.target.edge_labels,
                    };
                    verify_prepared_modes_against_expected(
                        &prepared_query,
                        &prepared_target,
                        target_sym,
                        &mut pointwise_pred,
                        &mut scratch,
                        expected,
                        &format!("{ctx} pointwise"),
                    );
                }
            }
        }

        // For the unsymmetrized mode, verify direct `Graph` methods and anchored `_at`.
        if !pattern_sym && !host_sym {
            let expected_all = baseline.as_deref().unwrap();
            let expected_anchored: Vec<Vec<usize>> = if case.query.graph.node_count() > 0 {
                expected_all
                    .iter()
                    .filter(|emb| emb[case.anchor.0.index()] == case.anchor.1.index())
                    .cloned()
                    .collect()
            } else {
                Vec::new()
            };

            for &alg in &CORE_ALGORITHMS {
                let mut nm = |q: NodeId, t: NodeId| {
                    (case.query.node_masks[q.index()] & (1u8 << case.target.node_labels[t.index()]))
                        != 0
                };
                let mut em = |qe: EdgeId, te: EdgeId| {
                    (case.query.edge_masks[qe.index()]
                        & (1u8 << case.target.edge_labels[te.index()]))
                        != 0
                };

                let mut direct_visit: Vec<Vec<usize>> = Vec::new();
                let _: ControlFlow<()> = case.target.graph.visit_subgraph_isomorphisms(
                    &case.query.graph,
                    &mut nm,
                    &mut em,
                    alg,
                    |emb| {
                        direct_visit.push(emb.iter().map(|n| n.index()).collect());
                        ControlFlow::Continue(())
                    },
                );
                let mut direct_sorted = direct_visit.clone();
                direct_sorted.sort();
                assert_eq!(direct_sorted, expected_all, "{alg:?} direct visit mismatch");
                assert_eq!(
                    case.target.graph.has_subgraph_isomorphism(
                        &case.query.graph,
                        &mut nm,
                        &mut em,
                        alg
                    ),
                    !expected_all.is_empty(),
                    "{alg:?} direct has mismatch"
                );
                assert_eq!(
                    case.target.graph.count_subgraph_isomorphisms(
                        &case.query.graph,
                        &mut nm,
                        &mut em,
                        alg
                    ),
                    expected_all.len() as u64,
                    "{alg:?} direct count mismatch"
                );
                let direct_first = case.target.graph.first_subgraph_isomorphism(
                    &case.query.graph,
                    &mut nm,
                    &mut em,
                    alg,
                );
                match direct_first {
                    None => assert!(expected_all.is_empty()),
                    Some(c) => {
                        let got: Vec<usize> =
                            c.matched_pairs().iter().map(|&(_, t)| t.index()).collect();
                        assert_eq!(got, direct_visit[0]);
                    }
                }

                if case.query.graph.node_count() > 0 && case.target.graph.node_count() > 0 {
                    let mut anc_visit: Vec<Vec<usize>> = Vec::new();
                    let _: ControlFlow<()> = case.target.graph.visit_subgraph_isomorphisms_at(
                        &case.query.graph,
                        case.anchor,
                        &mut nm,
                        &mut em,
                        alg,
                        |emb| {
                            anc_visit.push(emb.iter().map(|n| n.index()).collect());
                            ControlFlow::Continue(())
                        },
                    );
                    let mut anc_sorted = anc_visit.clone();
                    anc_sorted.sort();
                    assert_eq!(
                        anc_sorted, expected_anchored,
                        "{alg:?} anchored visit mismatch"
                    );
                    assert_eq!(
                        case.target.graph.has_subgraph_isomorphism_at(
                            &case.query.graph,
                            case.anchor,
                            &mut nm,
                            &mut em,
                            alg
                        ),
                        !expected_anchored.is_empty(),
                        "{alg:?} anchored has mismatch"
                    );
                    assert_eq!(
                        case.target.graph.count_subgraph_isomorphisms_at(
                            &case.query.graph,
                            case.anchor,
                            &mut nm,
                            &mut em,
                            alg
                        ),
                        expected_anchored.len() as u64,
                        "{alg:?} anchored count mismatch"
                    );
                    let anc_first = case.target.graph.first_subgraph_isomorphism_at(
                        &case.query.graph,
                        case.anchor,
                        &mut nm,
                        &mut em,
                        alg,
                    );
                    match anc_first {
                        None => assert!(expected_anchored.is_empty()),
                        Some(c) => {
                            let got: Vec<usize> =
                                c.matched_pairs().iter().map(|&(_, t)| t.index()).collect();
                            assert_eq!(got, anc_visit[0]);
                        }
                    }
                }
            }
        }
    }
}

/// Builds a deterministic boundary target of exact order `n2` containing two
/// symmetric copies of a labeled `q_len`-node cycle/ladder motif (with the
/// second copy ending at vertex `n2 - 1` so the top bit of the final word is
/// exercised) plus an asymmetric sparse background path.
fn boundary_case(q_len: usize, n2: usize) -> SparseRegimeCase {
    if q_len == 0 {
        let edges: Vec<[u32; 2]> = (1..n2 as u32).map(|i| [i - 1, i]).collect();
        let target = LabeledGraph {
            graph: Graph::new(n2, &edges),
            node_labels: vec![0; n2],
            edge_labels: vec![0; edges.len()],
        };
        let query = MaskedQuery {
            graph: Graph::new(0, &[]),
            node_masks: Vec::new(),
            edge_masks: Vec::new(),
        };
        return SparseRegimeCase {
            target,
            query,
            expect_match: true,
            hook_mode: EdgeHookMode::Mixed,
            anchor: (NodeId(0), NodeId(0)),
        };
    }

    assert!(
        n2 >= 2 * q_len,
        "boundary target must fit two copies of the motif"
    );

    // Motif: an even-length cycle of length `q_len` with alternating node
    // labels `[0, 1, 0, 1, ...]` and alternating edge labels `[0, 1, 0, 1, ...]`,
    // giving non-trivial cyclic symmetry on the query AND transposition symmetry
    // between the two target copies.
    let mut motif_edges: Vec<[u32; 2]> = Vec::with_capacity(q_len);
    let mut motif_node_labels: Vec<u8> = Vec::with_capacity(q_len);
    let mut motif_edge_labels: Vec<u8> = Vec::with_capacity(q_len);
    for i in 0..q_len {
        motif_node_labels.push((i % 2) as u8);
        let next = (i + 1) % q_len;
        motif_edges.push([i as u32, next as u32]);
        motif_edge_labels.push((i % 2) as u8);
    }

    // Place copy 0 at `0..q_len` and copy 1 at `(n2 - q_len)..n2` so both low
    // bit 0 and high bit `n2 - 1` participate in real matches and symmetries.
    let copy1_start = n2 - q_len;
    let mut target_node_labels = vec![2u8; n2];
    let mut target_edges = Vec::new();
    let mut target_edge_labels = Vec::new();

    for base in [0usize, copy1_start] {
        for (i, &lbl) in motif_node_labels.iter().enumerate() {
            target_node_labels[base + i] = lbl;
        }
        for (idx, &[u, v]) in motif_edges.iter().enumerate() {
            target_edges.push([(base + u as usize) as u32, (base + v as usize) as u32]);
            target_edge_labels.push(motif_edge_labels[idx]);
        }
    }

    // Asymmetric background path on `q_len..copy1_start` with label 2 (disjoint
    // from motif labels 0 and 1 so `Aut(H)` acts cleanly on the two motif copies
    // even at n2 = 4120).
    for i in (q_len + 1)..copy1_start {
        target_edges.push([(i - 1) as u32, i as u32]);
        target_edge_labels.push(2);
    }

    let query = MaskedQuery {
        graph: Graph::new(q_len, &motif_edges),
        node_masks: motif_node_labels.iter().map(|&l| 1u8 << l).collect(),
        edge_masks: motif_edge_labels
            .iter()
            .enumerate()
            .map(|(i, &l)| {
                if i.is_multiple_of(2) {
                    1u8 << l
                } else {
                    ALL_EDGE_MASK
                }
            })
            .collect(),
    };

    SparseRegimeCase {
        target: LabeledGraph {
            graph: Graph::new(n2, &target_edges),
            node_labels: target_node_labels,
            edge_labels: target_edge_labels,
        },
        query,
        expect_match: true,
        hook_mode: EdgeHookMode::Mixed,
        anchor: (NodeId(0), NodeId(copy1_start as u32)),
    }
}

#[rstest]
fn test_prepared_subgraph_query_corpus_through_order_5() {
    let graphs: Vec<Graph> = simple_graphs()
        .take_while(|graph| graph.node_count() <= 5)
        .collect();
    assert_eq!(graphs.len(), 52);

    let mut auto_scratch = AutomorphismScratch::new();
    let trivial_groups: Vec<Vec<Vec<usize>>> = graphs
        .iter()
        .map(|g| vec![(0..g.node_count()).collect()])
        .collect();
    let exact_groups: Vec<Vec<Vec<usize>>> = graphs
        .iter()
        .map(|g| exhaustive_automorphisms(g, |_| 0u8, |_| 0u8))
        .collect();
    let priority_orders: Vec<Vec<usize>> = graphs.iter().map(query_priority_order).collect();

    let q_syms_trail: Vec<QuerySymmetryConstraints> = graphs
        .iter()
        .map(|g| {
            QuerySymmetryConstraints::from_graph_with_scratch(
                g,
                |_| 0u8,
                |_| 0u8,
                AutomorphismAlgorithm::TrailIr,
                &mut auto_scratch,
            )
        })
        .collect();
    let q_syms_nauty: Vec<QuerySymmetryConstraints> = graphs
        .iter()
        .map(|g| {
            QuerySymmetryConstraints::from_graph_with_scratch(
                g,
                |_| 0u8,
                |_| 0u8,
                AutomorphismAlgorithm::Nauty,
                &mut auto_scratch,
            )
        })
        .collect();
    let t_syms_trail: Vec<TargetSymmetry> = graphs
        .iter()
        .map(|g| {
            let group = g.automorphisms_with_scratch(
                |_| 0u8,
                |_| 0u8,
                AutomorphismAlgorithm::TrailIr,
                &mut auto_scratch,
            );
            TargetSymmetry::from_automorphism_group(&group)
        })
        .collect();
    let t_syms_nauty: Vec<TargetSymmetry> = graphs
        .iter()
        .map(|g| {
            TargetSymmetry::from_graph_with_scratch(
                g,
                |_| 0u8,
                |_| 0u8,
                AutomorphismAlgorithm::Nauty,
                &mut auto_scratch,
            )
        })
        .collect();

    for i in 0..graphs.len() {
        assert_eq!(q_syms_trail[i].constraints(), q_syms_nauty[i].constraints());
        assert_eq!(
            target_orbits(&t_syms_trail[i]),
            target_orbits(&t_syms_nauty[i])
        );
        assert_eq!(t_syms_trail[i].orbit_count(), t_syms_nauty[i].orbit_count());
    }

    let prepared_targets: Vec<PreparedSubgraphTarget> =
        graphs.iter().map(PreparedSubgraphTarget::new).collect();
    let prepared_queries_plain: Vec<Vec<PreparedSubgraphQuery>> = graphs
        .iter()
        .map(|g| {
            CORE_ALGORITHMS
                .iter()
                .map(|&alg| PreparedSubgraphQuery::new(g, alg))
                .collect()
        })
        .collect();
    let prepared_queries_sym: Vec<Vec<PreparedSubgraphQuery>> = graphs
        .iter()
        .enumerate()
        .map(|(idx, g)| {
            CORE_ALGORITHMS
                .iter()
                .map(|&alg| {
                    PreparedSubgraphQuery::new(g, alg).with_symmetry(q_syms_trail[idx].clone())
                })
                .collect()
        })
        .collect();

    let mut scratch = SubgraphIsomorphismScratch::new();
    let mut pred = UnlabeledPredicate;

    for qi in 0..graphs.len() {
        for ti in 0..graphs.len() {
            let all_embeddings = exhaustive_monomorphisms(
                &graphs[qi],
                &graphs[ti],
                &mut |_, _| true,
                &mut |_, _| true,
            );

            for &(pattern_sym, host_sym) in &SYMMETRY_MODES {
                let q_group = if pattern_sym {
                    &exact_groups[qi]
                } else {
                    &trivial_groups[qi]
                };
                let t_group = if host_sym {
                    &exact_groups[ti]
                } else {
                    &trivial_groups[ti]
                };
                let expected = canonical_orbit_representatives(
                    &all_embeddings,
                    q_group,
                    t_group,
                    &priority_orders[qi],
                );
                let target_sym = host_sym.then_some(&t_syms_trail[ti]);
                let queries_for_q = if pattern_sym {
                    &prepared_queries_sym[qi]
                } else {
                    &prepared_queries_plain[qi]
                };

                for (alg_idx, &alg) in CORE_ALGORITHMS.iter().enumerate() {
                    let ctx = format!(
                        "corpus<=5 q={qi} t={ti} alg={alg:?} sym=({pattern_sym},{host_sym})"
                    );
                    verify_prepared_modes_against_expected(
                        &queries_for_q[alg_idx],
                        &prepared_targets[ti],
                        target_sym,
                        &mut pred,
                        &mut scratch,
                        &expected,
                        &ctx,
                    );
                }
            }
        }
    }
}

#[rstest]
fn test_prepared_subgraph_query_corpus_order_6_targets() {
    let queries: Vec<Graph> = simple_graphs()
        .take_while(|graph| graph.node_count() <= 4)
        .collect();
    let targets: Vec<Graph> = simple_graphs()
        .skip_while(|graph| graph.node_count() < 6)
        .take_while(|graph| graph.node_count() == 6)
        .collect();
    assert_eq!(queries.len(), 18);
    assert_eq!(targets.len(), 156);

    let mut auto_scratch = AutomorphismScratch::new();
    let q_trivial: Vec<Vec<Vec<usize>>> = queries
        .iter()
        .map(|g| vec![(0..g.node_count()).collect()])
        .collect();
    let q_exact: Vec<Vec<Vec<usize>>> = queries
        .iter()
        .map(|g| exhaustive_automorphisms(g, |_| 0u8, |_| 0u8))
        .collect();
    let q_priority: Vec<Vec<usize>> = queries.iter().map(query_priority_order).collect();
    let q_syms: Vec<QuerySymmetryConstraints> = queries
        .iter()
        .map(|g| {
            let trail = QuerySymmetryConstraints::from_graph_with_scratch(
                g,
                |_| 0u8,
                |_| 0u8,
                AutomorphismAlgorithm::TrailIr,
                &mut auto_scratch,
            );
            let nauty = QuerySymmetryConstraints::from_graph_with_scratch(
                g,
                |_| 0u8,
                |_| 0u8,
                AutomorphismAlgorithm::Nauty,
                &mut auto_scratch,
            );
            assert_eq!(trail.constraints(), nauty.constraints());
            trail
        })
        .collect();

    let prepared_queries_plain: Vec<Vec<PreparedSubgraphQuery>> = queries
        .iter()
        .map(|g| {
            CORE_ALGORITHMS
                .iter()
                .map(|&alg| PreparedSubgraphQuery::new(g, alg))
                .collect()
        })
        .collect();
    let prepared_queries_sym: Vec<Vec<PreparedSubgraphQuery>> = queries
        .iter()
        .enumerate()
        .map(|(qi, g)| {
            CORE_ALGORITHMS
                .iter()
                .map(|&alg| PreparedSubgraphQuery::new(g, alg).with_symmetry(q_syms[qi].clone()))
                .collect()
        })
        .collect();

    let mut scratch = SubgraphIsomorphismScratch::new();
    let mut pred = UnlabeledPredicate;
    let t_trivial = vec![vec![0usize, 1, 2, 3, 4, 5]];

    for (ti, target) in targets.iter().enumerate() {
        let t_exact = exhaustive_automorphisms(target, |_| 0u8, |_| 0u8);
        let t_sym_trail = TargetSymmetry::from_graph_with_scratch(
            target,
            |_| 0u8,
            |_| 0u8,
            AutomorphismAlgorithm::TrailIr,
            &mut auto_scratch,
        );
        let t_sym_nauty = TargetSymmetry::from_graph_with_scratch(
            target,
            |_| 0u8,
            |_| 0u8,
            AutomorphismAlgorithm::Nauty,
            &mut auto_scratch,
        );
        assert_eq!(target_orbits(&t_sym_trail), target_orbits(&t_sym_nauty));
        assert_eq!(t_sym_trail.orbit_count(), t_sym_nauty.orbit_count());
        let prepared_target = PreparedSubgraphTarget::new(target);

        for (qi, query) in queries.iter().enumerate() {
            let all_embeddings =
                exhaustive_monomorphisms(query, target, &mut |_, _| true, &mut |_, _| true);

            for &(pattern_sym, host_sym) in &SYMMETRY_MODES {
                let q_group = if pattern_sym {
                    &q_exact[qi]
                } else {
                    &q_trivial[qi]
                };
                let t_group = if host_sym { &t_exact } else { &t_trivial };
                let expected = canonical_orbit_representatives(
                    &all_embeddings,
                    q_group,
                    t_group,
                    &q_priority[qi],
                );
                let target_sym = host_sym.then_some(&t_sym_trail);
                let queries_for_q = if pattern_sym {
                    &prepared_queries_sym[qi]
                } else {
                    &prepared_queries_plain[qi]
                };

                for (alg_idx, &alg) in CORE_ALGORITHMS.iter().enumerate() {
                    let ctx =
                        format!("corpus6 q={qi} t={ti} alg={alg:?} sym=({pattern_sym},{host_sym})");
                    verify_prepared_modes_against_expected(
                        &queries_for_q[alg_idx],
                        &prepared_target,
                        target_sym,
                        &mut pred,
                        &mut scratch,
                        &expected,
                        &ctx,
                    );
                }
            }
        }
    }
}

#[rstest]
#[case::empty_in_empty(0, 0)]
#[case::empty_in_u64(0, 8)]
#[case::u64_exact_64(6, 64)]
#[case::u64_large_query_18(18, 64)]
#[case::u128_boundary_65(6, 65)]
#[case::u128_exact_128(6, 128)]
#[case::u128_large_query_18(18, 128)]
#[case::wide_boundary_129(6, 129)]
#[case::wide_large_query_18(18, 160)]
#[case::wide_above_eager_limit_1050(6, 1050)]
#[case::wide_above_lazy_limit_4120(6, 4120)]
fn test_subgraph_isomorphism_width_boundaries(#[case] q_len: usize, #[case] n2: usize) {
    let case = boundary_case(q_len, n2);
    verify_sparse_regime_case(&case);
}

#[rstest]
#[case::u64_high_degree(64)]
#[case::u128_high_degree(128)]
#[case::wide_high_degree(140)]
fn test_subgraph_isomorphism_high_degree_across_widths(#[case] n2: usize) {
    // Query: wheel W_17 (1 hub of degree 17 + 17-cycle rim = 18 nodes).
    // Target: two symmetric copies of W_17 (one ending at vertex `n2 - 1`) plus
    // a decoy W_16 (hub degree 16, rejected only by exact `deg(q) <= deg(t)`).
    let wheel_rim = 17usize;
    let wheel_n = wheel_rim + 1;
    let decoy_rim = 16usize;
    let decoy_n = decoy_rim + 1;

    let mut q_edges = Vec::new();
    let mut q_edge_masks = Vec::new();
    for i in 0..wheel_rim {
        let u = (i + 1) as u32;
        let v = (((i + 1) % wheel_rim) + 1) as u32;
        q_edges.push([0, u]);
        q_edge_masks.push(1u8);
        q_edges.push([u, v]);
        q_edge_masks.push(ALL_EDGE_MASK);
    }
    let mut q_node_masks = vec![1u8 << 1; wheel_n];
    q_node_masks[0] = 1u8;

    let mut t_edges = Vec::new();
    let mut t_edge_labels = Vec::new();
    let mut t_node_labels = vec![2u8; n2];
    let copy1_base = n2 - wheel_n;
    for base in [0usize, copy1_base] {
        t_node_labels[base] = 0;
        for i in 0..wheel_rim {
            t_node_labels[base + 1 + i] = 1;
            let u = (base + 1 + i) as u32;
            let v = (base + 1 + ((i + 1) % wheel_rim)) as u32;
            t_edges.push([base as u32, u]);
            t_edge_labels.push(0);
            t_edges.push([u, v]);
            t_edge_labels.push(1);
        }
    }
    // Decoy W_16 sits directly after copy 0.
    let decoy_base = wheel_n;
    t_node_labels[decoy_base] = 0;
    for i in 0..decoy_rim {
        t_node_labels[decoy_base + 1 + i] = 1;
        let u = (decoy_base + 1 + i) as u32;
        let v = (decoy_base + 1 + ((i + 1) % decoy_rim)) as u32;
        t_edges.push([decoy_base as u32, u]);
        t_edge_labels.push(0);
        t_edges.push([u, v]);
        t_edge_labels.push(1);
    }
    // Asymmetric tail on remaining unused vertices so `Aut(H)` stays small.
    let tail_start = decoy_base + decoy_n;
    for i in (tail_start + 1)..copy1_base {
        t_edges.push([(i - 1) as u32, i as u32]);
        t_edge_labels.push(2);
    }

    let case = SparseRegimeCase {
        target: LabeledGraph {
            graph: Graph::new(n2, &t_edges),
            node_labels: t_node_labels,
            edge_labels: t_edge_labels,
        },
        query: MaskedQuery {
            graph: Graph::new(wheel_n, &q_edges),
            node_masks: q_node_masks,
            edge_masks: q_edge_masks,
        },
        expect_match: true,
        hook_mode: EdgeHookMode::Mixed,
        anchor: (NodeId(0), NodeId(copy1_base as u32)),
    };
    verify_sparse_regime_case(&case);
}

#[rstest]
#[case::c4_rotation_in_two_c4(
    4,
    &[[0, 1], [1, 2], [2, 3], [3, 0]],
    vec![vec![1, 2, 3, 0]],
    8,
    &[[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4]],
    vec![vec![1, 2, 3, 0, 5, 6, 7, 4], vec![4, 5, 6, 7, 0, 1, 2, 3]],
)]
#[case::a4_even_subgroup_of_k4_in_k5(
    4,
    &[[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]],
    vec![vec![1, 2, 0, 3], vec![1, 0, 3, 2]],
    5,
    &[[0, 1], [0, 2], [0, 3], [0, 4], [1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]],
    vec![vec![1, 2, 3, 4, 0]],
)]
#[case::disconnected_pair_proper_subgroup(
    4,
    &[[0, 1], [2, 3]],
    vec![vec![2, 3, 0, 1]],
    6,
    &[[0, 1], [1, 2], [2, 0], [3, 4], [4, 5], [5, 3]],
    vec![vec![1, 2, 0, 4, 5, 3]],
)]
fn test_prepared_subgraph_query_proper_subgroups(
    #[case] q_n: usize,
    #[case] q_edges: &[[u32; 2]],
    #[case] q_gen_raw: Vec<Vec<u32>>,
    #[case] t_n: usize,
    #[case] t_edges: &[[u32; 2]],
    #[case] t_gen_raw: Vec<Vec<u32>>,
) {
    let query = Graph::new(q_n, q_edges);
    let target = Graph::new(t_n, t_edges);
    let q_gens: Vec<Vec<NodeId>> = q_gen_raw
        .iter()
        .map(|g| g.iter().map(|&v| NodeId(v)).collect())
        .collect();
    let t_gens: Vec<Vec<NodeId>> = t_gen_raw
        .iter()
        .map(|g| g.iter().map(|&v| NodeId(v)).collect())
        .collect();

    let q_sym = QuerySymmetryConstraints::try_from_generators(&query, &q_gens)
        .expect("valid proper query subgroup");
    let t_sym = TargetSymmetry::try_from_generators(&target, &t_gens)
        .expect("valid proper target subgroup");

    let q_exact = subgroup_from_generators(q_n, &q_gens);
    let t_exact = subgroup_from_generators(t_n, &t_gens);
    let q_trivial = vec![(0..q_n).collect()];
    let t_trivial = vec![(0..t_n).collect()];
    let priority = query_priority_order(&query);

    let all_embeddings =
        exhaustive_monomorphisms(&query, &target, &mut |_, _| true, &mut |_, _| true);
    let prepared_target = PreparedSubgraphTarget::new(&target);
    let mut scratch = SubgraphIsomorphismScratch::new();
    let mut pred = UnlabeledPredicate;

    for &(pattern_sym, host_sym) in &SYMMETRY_MODES {
        let expected = canonical_orbit_representatives(
            &all_embeddings,
            if pattern_sym { &q_exact } else { &q_trivial },
            if host_sym { &t_exact } else { &t_trivial },
            &priority,
        );
        let target_sym = host_sym.then_some(&t_sym);
        for &alg in &CORE_ALGORITHMS {
            let mut prepared_query = PreparedSubgraphQuery::new(&query, alg);
            if pattern_sym {
                prepared_query = prepared_query.with_symmetry(q_sym.clone());
            }
            let ctx = format!("proper_subgroup alg={alg:?} sym=({pattern_sym},{host_sym})");
            verify_prepared_modes_against_expected(
                &prepared_query,
                &prepared_target,
                target_sym,
                &mut pred,
                &mut scratch,
                &expected,
                &ctx,
            );
        }
    }
}

#[rstest]
fn test_prepared_subgraph_query_scratch_reuse_across_regimes() {
    let cases = [
        boundary_case(6, 1050),
        boundary_case(6, 64),
        boundary_case(6, 128),
        boundary_case(6, 140),
        boundary_case(0, 8),
        boundary_case(6, 64),
    ];

    let mut shared_scratch = SubgraphIsomorphismScratch::new();
    for (idx, case) in cases.iter().enumerate() {
        let prepared_target = PreparedSubgraphTarget::new(&case.target.graph);
        let q_sym = QuerySymmetryConstraints::from_graph(
            &case.query.graph,
            |q| case.query.node_masks[q.index()],
            |qe| case.query.edge_masks[qe.index()],
            AutomorphismAlgorithm::TrailIr,
        );
        let t_sym = TargetSymmetry::from_graph(
            &case.target.graph,
            |t| case.target.node_labels[t.index()],
            |te| case.target.edge_labels[te.index()],
            AutomorphismAlgorithm::TrailIr,
        );

        for &alg in &CORE_ALGORITHMS {
            let prepared_query =
                PreparedSubgraphQuery::new(&case.query.graph, alg).with_symmetry(q_sym.clone());

            // Force early termination on the shared scratch before checking full results.
            let mut pred_early =
                BitsetMaskPredicate::new(&case.query, &case.target, case.hook_mode);
            let _: ControlFlow<()> = prepared_query.visit_matches(
                &prepared_target,
                Some(&t_sym),
                &mut pred_early,
                &mut shared_scratch,
                |_| ControlFlow::Break(()),
            );

            let mut pred_fresh =
                BitsetMaskPredicate::new(&case.query, &case.target, case.hook_mode);
            let mut fresh_scratch = SubgraphIsomorphismScratch::new();
            let mut expected: Vec<Vec<usize>> = Vec::new();
            let _: ControlFlow<()> = prepared_query.visit_matches(
                &prepared_target,
                Some(&t_sym),
                &mut pred_fresh,
                &mut fresh_scratch,
                |emb| {
                    expected.push(emb.iter().map(|n| n.index()).collect());
                    ControlFlow::Continue(())
                },
            );
            expected.sort();

            let mut pred_shared =
                BitsetMaskPredicate::new(&case.query, &case.target, case.hook_mode);
            verify_prepared_modes_against_expected(
                &prepared_query,
                &prepared_target,
                Some(&t_sym),
                &mut pred_shared,
                &mut shared_scratch,
                &expected,
                &format!("scratch_reuse step={idx} alg={alg:?}"),
            );
        }
    }
}

proptest! {
    #[test]
    fn test_graph_subgraph_isomorphism_exact_and_mode_splits(
        query in labeled_graph(0, 4, 2, 2),
        target in labeled_graph(0, 6, 2, 2),
    ) {
        let mut nm =
            |q: NodeId, t: NodeId| query.node_labels[q.index()] == target.node_labels[t.index()];
        let mut em =
            |qe: EdgeId, te: EdgeId| query.edge_labels[qe.index()] == target.edge_labels[te.index()];
        let expected = exhaustive_monomorphisms(&query.graph, &target.graph, &mut nm, &mut em);

        for alg in CORE_ALGORITHMS {
            verify_direct_modes_against_expected(&query, &target, alg, &expected);
        }
    }

    #[test]
    fn test_graph_subgraph_isomorphism_planted_non_induced(
        (target, nodes, keep_bits) in labeled_graph(1, 6, 2, 2).prop_flat_map(|t| {
            let n = t.graph.node_count();
            let cap = n.min(4);
            (
                Just(t),
                prop::sample::subsequence((0..n).collect::<Vec<_>>(), 1..=cap),
                any::<u64>(),
            )
        }),
    ) {
        let (query, _) = planted_subgraph(&target, &nodes, keep_bits, &[], &[]);
        let mut nm =
            |q: NodeId, t: NodeId| query.node_labels[q.index()] == target.node_labels[t.index()];
        let mut em =
            |qe: EdgeId, te: EdgeId| query.edge_labels[qe.index()] == target.edge_labels[te.index()];
        let expected = exhaustive_monomorphisms(&query.graph, &target.graph, &mut nm, &mut em);
        prop_assert!(!expected.is_empty(), "planted non-induced query must embed in its source");

        for alg in CORE_ALGORITHMS {
            verify_direct_modes_against_expected(&query, &target, alg, &expected);
        }
    }

    #[test]
    fn test_graph_subgraph_isomorphism_anchored(
        query in labeled_graph(1, 3, 2, 2),
        target in labeled_graph(1, 5, 2, 2),
    ) {
        let mut nm =
            |q: NodeId, t: NodeId| query.node_labels[q.index()] == target.node_labels[t.index()];
        let mut em =
            |qe: EdgeId, te: EdgeId| query.edge_labels[qe.index()] == target.edge_labels[te.index()];
        let all_expected = exhaustive_monomorphisms(&query.graph, &target.graph, &mut nm, &mut em);

        for alg in CORE_ALGORITHMS {
            for qa in 0..query.graph.node_count() {
                for ta in 0..target.graph.node_count() {
                    let anchor = (NodeId(qa as u32), NodeId(ta as u32));
                    let expected: Vec<Vec<usize>> = all_expected
                        .iter()
                        .filter(|emb| emb[qa] == ta)
                        .cloned()
                        .collect();

                    let mut via_visit_in_order: Vec<Vec<usize>> = Vec::new();
                    let _: ControlFlow<()> = target.graph.visit_subgraph_isomorphisms_at(
                        &query.graph,
                        anchor,
                        &mut nm,
                        &mut em,
                        alg,
                        |embedding| {
                            via_visit_in_order.push(embedding.iter().map(|n| n.index()).collect());
                            ControlFlow::Continue(())
                        },
                    );
                    let mut via_visit_sorted = via_visit_in_order.clone();
                    via_visit_sorted.sort();
                    prop_assert_eq!(
                        &via_visit_sorted,
                        &expected,
                        "{:?} anchored at {:?}: visit mismatch",
                        alg,
                        anchor
                    );

                    let mut via_enumerate: Vec<Vec<usize>> = target
                        .graph
                        .enumerate_subgraph_isomorphisms_at(
                            &query.graph,
                            anchor,
                            &mut nm,
                            &mut em,
                            alg,
                        )
                        .iter()
                        .map(|c| c.matched_pairs().iter().map(|&(_, t)| t.index()).collect())
                        .collect();
                    via_enumerate.sort();
                    prop_assert_eq!(
                        &via_enumerate,
                        &expected,
                        "{:?} anchored at {:?}: enumerate mismatch",
                        alg,
                        anchor
                    );

                    prop_assert_eq!(
                        target.graph.has_subgraph_isomorphism_at(
                            &query.graph,
                            anchor,
                            &mut nm,
                            &mut em,
                            alg
                        ),
                        !expected.is_empty(),
                        "{:?} anchored at {:?}: has mismatch",
                        alg,
                        anchor
                    );
                    prop_assert_eq!(
                        target.graph.count_subgraph_isomorphisms_at(
                            &query.graph,
                            anchor,
                            &mut nm,
                            &mut em,
                            alg
                        ),
                        expected.len() as u64,
                        "{:?} anchored at {:?}: count mismatch",
                        alg,
                        anchor
                    );
                    let first = target.graph.first_subgraph_isomorphism_at(
                        &query.graph,
                        anchor,
                        &mut nm,
                        &mut em,
                        alg,
                    );
                    match first {
                        None => prop_assert!(expected.is_empty()),
                        Some(corr) => {
                            let first_vec: Vec<usize> = corr
                                .matched_pairs()
                                .iter()
                                .map(|&(_, t)| t.index())
                                .collect();
                            prop_assert_eq!(first_vec, via_visit_in_order[0].clone());
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn test_prepared_subgraph_query_symmetry_and_predicates_u64(
        query in masked_query(0, 4, NUM_NODE_LABELS, NUM_EDGE_LABELS),
        target in labeled_graph(0, 6, NUM_NODE_LABELS as u8, NUM_EDGE_LABELS as u8),
        hook_idx in 0u8..3,
    ) {
        let hook_mode = match hook_idx {
            0 => EdgeHookMode::AllSome,
            1 => EdgeHookMode::Mixed,
            _ => EdgeHookMode::AllNone,
        };

        let mut nm = |q: NodeId, t: NodeId| {
            (query.node_masks[q.index()] & (1u8 << target.node_labels[t.index()])) != 0
        };
        let mut em = |qe: EdgeId, te: EdgeId| {
            (query.edge_masks[qe.index()] & (1u8 << target.edge_labels[te.index()])) != 0
        };
        let all_embeddings = exhaustive_monomorphisms(&query.graph, &target.graph, &mut nm, &mut em);

        let q_exact = exhaustive_automorphisms(
            &query.graph,
            |q| query.node_masks[q.index()],
            |qe| query.edge_masks[qe.index()],
        );
        let t_exact = exhaustive_automorphisms(
            &target.graph,
            |t| target.node_labels[t.index()],
            |te| target.edge_labels[te.index()],
        );
        let q_trivial = vec![(0..query.graph.node_count()).collect()];
        let t_trivial = vec![(0..target.graph.node_count()).collect()];
        let priority = query_priority_order(&query.graph);

        let mut auto_scratch = AutomorphismScratch::new();
        let q_syms: Vec<QuerySymmetryConstraints> = AUTOMORPHISM_ALGORITHMS
            .iter()
            .map(|&auto_alg| {
                QuerySymmetryConstraints::from_graph_with_scratch(
                    &query.graph,
                    |q| query.node_masks[q.index()],
                    |qe| query.edge_masks[qe.index()],
                    auto_alg,
                    &mut auto_scratch,
                )
            })
            .collect();
        prop_assert_eq!(q_syms[0].constraints(), q_syms[1].constraints());

        let t_syms: Vec<TargetSymmetry> = AUTOMORPHISM_ALGORITHMS
            .iter()
            .map(|&auto_alg| {
                TargetSymmetry::from_graph_with_scratch(
                    &target.graph,
                    |t| target.node_labels[t.index()],
                    |te| target.edge_labels[te.index()],
                    auto_alg,
                    &mut auto_scratch,
                )
            })
            .collect();
        prop_assert_eq!(target_orbits(&t_syms[0]), target_orbits(&t_syms[1]));
        prop_assert_eq!(t_syms[0].orbit_count(), t_syms[1].orbit_count());

        let prepared_target = PreparedSubgraphTarget::new(&target.graph);
        let mut scratch = SubgraphIsomorphismScratch::new();

        for &(pattern_sym, host_sym) in &SYMMETRY_MODES {
            let expected = canonical_orbit_representatives(
                &all_embeddings,
                if pattern_sym { &q_exact } else { &q_trivial },
                if host_sym { &t_exact } else { &t_trivial },
                &priority,
            );
            for (auto_idx, &auto_alg) in AUTOMORPHISM_ALGORITHMS.iter().enumerate() {
                let target_sym = host_sym.then_some(&t_syms[auto_idx]);
                for &alg in &CORE_ALGORITHMS {
                    let mut prepared_query = PreparedSubgraphQuery::new(&query.graph, alg);
                    if pattern_sym {
                        prepared_query = prepared_query.with_symmetry(q_syms[auto_idx].clone());
                    }
                    let ctx = format!(
                        "u64 alg={alg:?} auto={auto_alg:?} sym=({pattern_sym},{host_sym}) hook={hook_mode:?}"
                    );
                    let mut bitset_pred = BitsetMaskPredicate::new(&query, &target, hook_mode);
                    verify_prepared_modes_against_expected(
                        &prepared_query,
                        &prepared_target,
                        target_sym,
                        &mut bitset_pred,
                        &mut scratch,
                        &expected,
                        &ctx,
                    );
                    if auto_idx == 0 {
                        let mut pointwise_pred = PointwiseMaskPredicate {
                            query_node_masks: &query.node_masks,
                            query_edge_masks: &query.edge_masks,
                            target_node_labels: &target.node_labels,
                            target_edge_labels: &target.edge_labels,
                        };
                        verify_prepared_modes_against_expected(
                            &prepared_query,
                            &prepared_target,
                            target_sym,
                            &mut pointwise_pred,
                            &mut scratch,
                            &expected,
                            &format!("{ctx} pointwise"),
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn test_prepared_subgraph_query_planted_symmetry_u64(
        (target, nodes, keep_bits, extra_node_masks, extra_edge_masks, hook_idx) in
            labeled_graph(1, 6, 2, 2).prop_flat_map(|t| {
                let n = t.graph.node_count();
                let cap = n.min(4);
                (
                    Just(t),
                    prop::sample::subsequence((0..n).collect::<Vec<_>>(), 1..=cap),
                    any::<u64>(),
                    prop::collection::vec(0u8..=ALL_EDGE_MASK, cap),
                    prop::collection::vec(0u8..=ALL_EDGE_MASK, 6),
                    0u8..3,
                )
            }),
    ) {
        let (_, query) = planted_subgraph(
            &target,
            &nodes,
            keep_bits,
            &extra_node_masks,
            &extra_edge_masks,
        );
        let hook_mode = match hook_idx {
            0 => EdgeHookMode::AllSome,
            1 => EdgeHookMode::Mixed,
            _ => EdgeHookMode::AllNone,
        };

        let mut nm = |q: NodeId, t: NodeId| {
            (query.node_masks[q.index()] & (1u8 << target.node_labels[t.index()])) != 0
        };
        let mut em = |qe: EdgeId, te: EdgeId| {
            (query.edge_masks[qe.index()] & (1u8 << target.edge_labels[te.index()])) != 0
        };
        let all_embeddings = exhaustive_monomorphisms(&query.graph, &target.graph, &mut nm, &mut em);
        prop_assert!(!all_embeddings.is_empty(), "planted masked query must embed");

        let q_exact = exhaustive_automorphisms(
            &query.graph,
            |q| query.node_masks[q.index()],
            |qe| query.edge_masks[qe.index()],
        );
        let t_exact = exhaustive_automorphisms(
            &target.graph,
            |t| target.node_labels[t.index()],
            |te| target.edge_labels[te.index()],
        );
        let q_trivial = vec![(0..query.graph.node_count()).collect()];
        let t_trivial = vec![(0..target.graph.node_count()).collect()];
        let priority = query_priority_order(&query.graph);

        let mut auto_scratch = AutomorphismScratch::new();
        let q_sym = QuerySymmetryConstraints::from_graph_with_scratch(
            &query.graph,
            |q| query.node_masks[q.index()],
            |qe| query.edge_masks[qe.index()],
            AutomorphismAlgorithm::TrailIr,
            &mut auto_scratch,
        );
        let q_sym_nauty = QuerySymmetryConstraints::from_graph_with_scratch(
            &query.graph,
            |q| query.node_masks[q.index()],
            |qe| query.edge_masks[qe.index()],
            AutomorphismAlgorithm::Nauty,
            &mut auto_scratch,
        );
        prop_assert_eq!(q_sym.constraints(), q_sym_nauty.constraints());

        let t_sym = TargetSymmetry::from_graph_with_scratch(
            &target.graph,
            |t| target.node_labels[t.index()],
            |te| target.edge_labels[te.index()],
            AutomorphismAlgorithm::TrailIr,
            &mut auto_scratch,
        );
        let t_sym_nauty = TargetSymmetry::from_graph_with_scratch(
            &target.graph,
            |t| target.node_labels[t.index()],
            |te| target.edge_labels[te.index()],
            AutomorphismAlgorithm::Nauty,
            &mut auto_scratch,
        );
        prop_assert_eq!(target_orbits(&t_sym), target_orbits(&t_sym_nauty));

        let prepared_target = PreparedSubgraphTarget::new(&target.graph);
        let mut scratch = SubgraphIsomorphismScratch::new();

        for &(pattern_sym, host_sym) in &SYMMETRY_MODES {
            let expected = canonical_orbit_representatives(
                &all_embeddings,
                if pattern_sym { &q_exact } else { &q_trivial },
                if host_sym { &t_exact } else { &t_trivial },
                &priority,
            );
            prop_assert!(!expected.is_empty(), "non-empty embedding set must have canonical reps");
            let target_sym = host_sym.then_some(&t_sym);
            for &alg in &CORE_ALGORITHMS {
                let mut prepared_query = PreparedSubgraphQuery::new(&query.graph, alg);
                if pattern_sym {
                    prepared_query = prepared_query.with_symmetry(q_sym.clone());
                }
                let ctx = format!(
                    "planted_u64 alg={alg:?} sym=({pattern_sym},{host_sym}) hook={hook_mode:?}"
                );
                let mut bitset_pred = BitsetMaskPredicate::new(&query, &target, hook_mode);
                verify_prepared_modes_against_expected(
                    &prepared_query,
                    &prepared_target,
                    target_sym,
                    &mut bitset_pred,
                    &mut scratch,
                    &expected,
                    &ctx,
                );
                let mut pointwise_pred = PointwiseMaskPredicate {
                    query_node_masks: &query.node_masks,
                    query_edge_masks: &query.edge_masks,
                    target_node_labels: &target.node_labels,
                    target_edge_labels: &target.edge_labels,
                };
                verify_prepared_modes_against_expected(
                    &prepared_query,
                    &prepared_target,
                    target_sym,
                    &mut pointwise_pred,
                    &mut scratch,
                    &expected,
                    &format!("{ctx} pointwise"),
                );
            }
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(32))]

    #[test]
    fn test_prepared_subgraph_query_u128_regime(case in sparse_regime_case(65, 128)) {
        verify_sparse_regime_case(&case);
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn test_prepared_subgraph_query_wide_regime(case in sparse_regime_case(129, 196)) {
        verify_sparse_regime_case(&case);
    }
}
