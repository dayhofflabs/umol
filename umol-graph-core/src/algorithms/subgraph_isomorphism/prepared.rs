//! Reusable scratch storage and prepared query/target carriers for subgraph-isomorphism search.
//!
//! Separating query preparation ([`PreparedSubgraphQuery`]), target indexing
//! ([`PreparedSubgraphTarget`]), and transient search buffers ([`SubgraphIsomorphismScratch`])
//! amortizes plan compilation, symmetry analysis, degree histograms, and bit-parallel adjacency
//! across 1-pattern × many-targets and many-patterns × 1-target loops without caching semantic
//! state inside scratch buffers.

use std::cmp::Reverse;

use super::symmetry::{query_priority_order, QuerySymmetryConstraints, TargetStabilizerChain};
use super::SubgraphIsomorphismAlgorithm;
use crate::graph::Graph;
use crate::id::{EdgeId, NodeId};

/// Maximum target node count for eagerly materializing wide (`> 128`) bitset adjacency rows in
/// [`PreparedSubgraphTarget`]. Larger sparse targets build only touched rows lazily in scratch.
const WIDE_EAGER_ADJ_LIMIT: usize = 1024;

/// Opaque, stateless reusable buffer carrier for subgraph-isomorphism search.
///
/// Safe to reuse across arbitrary queries, targets, algorithms, and predicates.
///
/// # Semantic properties
///
/// Prior calls on a scratch instance never affect the output of subsequent calls; matching with
/// a reused scratch produces results identical to a fresh scratch.
#[derive(Clone, Debug, Default)]
#[allow(dead_code, reason = "buffers consumed by matching engines")]
pub struct SubgraphIsomorphismScratch {
    pub(crate) order: Vec<usize>,
    pub(crate) pos: Vec<usize>,
    pub(crate) parents_offsets: Vec<usize>,
    pub(crate) parents_data: Vec<(usize, EdgeId)>,
    pub(crate) mapping: Vec<usize>,
    pub(crate) used: Vec<bool>,
    pub(crate) out_nodes: Vec<NodeId>,
    pub(crate) placed: Vec<bool>,
    pub(crate) ordered_neighbors: Vec<usize>,
    pub(crate) cand_counts: Vec<usize>,
    pub(crate) cand_u64: Vec<u64>,
    pub(crate) cand_u128: Vec<u128>,
    pub(crate) cand_wide: Vec<u64>,
    pub(crate) class_u64: Vec<u64>,
    pub(crate) class_u128: Vec<u128>,
    pub(crate) class_wide: Vec<u64>,
    pub(crate) adj_wide: Vec<u64>,
    pub(crate) adj_built_wide: Vec<u64>,
    pub(crate) trail_u64: Vec<u64>,
    pub(crate) trail_u128: Vec<u128>,
    pub(crate) trail_wide: Vec<u64>,
    pub(crate) target_chain: TargetStabilizerChain,
}

impl SubgraphIsomorphismScratch {
    /// Allocate empty reusable scratch storage.
    pub fn new() -> Self {
        Self::default()
    }
}

/// Width-specialized degree-threshold and adjacency bitsets for a prepared target graph.
#[derive(Clone, Debug)]
#[allow(dead_code, reason = "consumed by bitset matching engines")]
pub(crate) enum TargetBitsets {
    U64 {
        deg_ge: Vec<u64>,
        adj: Vec<u64>,
    },
    U128 {
        deg_ge: Vec<u128>,
        adj: Vec<u128>,
    },
    Wide {
        words: usize,
        deg_ge: Vec<u64>,
        adj: Option<Vec<u64>>,
    },
}

/// Immutable, algorithm-independent prepared target graph for repeated subgraph-isomorphism
/// matching.
///
/// Holds an `Arc`-shared snapshot of the source [`Graph`] alongside its cumulative degree
/// histogram and width-specialized degree/adjacency bitsets, so mutating the source graph after
/// construction cannot affect the prepared target.
///
/// # Semantic properties
///
/// Matching any query against `PreparedSubgraphTarget::new(target)` produces the exact same
/// embedding set, existence, and count as matching directly against `target`.
#[derive(Clone, Debug)]
pub struct PreparedSubgraphTarget {
    graph: Graph,
    degree_hist: Vec<u32>,
    bitsets: TargetBitsets,
}

impl PreparedSubgraphTarget {
    /// Prepare `target` for repeated subgraph-isomorphism matching.
    pub fn new(target: &Graph) -> Self {
        let graph = target.clone();
        let degree_hist = build_degree_hist(&graph);
        let bitsets = build_target_bitsets(&graph);
        Self {
            graph,
            degree_hist,
            bitsets,
        }
    }

    /// Snapshotted target [`Graph`].
    pub fn graph(&self) -> &Graph {
        &self.graph
    }

    /// Number of nodes in the target graph.
    pub fn node_count(&self) -> usize {
        self.graph.node_count()
    }

    /// Number of edges in the target graph.
    pub fn edge_count(&self) -> usize {
        self.graph.edge_count()
    }

    /// Maximum vertex degree in the target graph (`0` for empty or edgeless graphs).
    pub fn max_degree(&self) -> usize {
        self.graph.max_degree()
    }

    #[allow(dead_code, reason = "consumed by matching engines")]
    pub(crate) fn degree_hist(&self) -> &[u32] {
        &self.degree_hist
    }

    #[allow(dead_code, reason = "consumed by matching engines")]
    pub(crate) fn bitsets(&self) -> &TargetBitsets {
        &self.bitsets
    }

    /// Necessary degree-sequence check: for every degree `d`, the target must contain at least as
    /// many vertices of degree `>= d` as the query.
    #[allow(dead_code, reason = "consumed by matching engines")]
    pub(crate) fn satisfies_degree_sequence(&self, query_hist: &[u32]) -> bool {
        if query_hist.len() > self.degree_hist.len() {
            return false;
        }
        query_hist
            .iter()
            .zip(self.degree_hist.iter())
            .all(|(&q, &t)| q <= t)
    }
}

/// Static vertex-order plan with flat CSR parent-edge lists (used by `Vf2Rdkit`).
#[derive(Clone, Debug)]
#[allow(dead_code, reason = "consumed by Vf2Rdkit matching engine")]
pub(crate) struct StaticOrderPlan {
    pub(crate) order: Vec<usize>,
    pub(crate) parents_offsets: Vec<usize>,
    pub(crate) parents_data: Vec<(usize, EdgeId)>,
}

/// Compiled bond-oriented instruction for `RayKirsch`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[allow(dead_code, reason = "consumed by RayKirsch matching engine")]
pub(crate) enum RkStep {
    /// Start a new connected component at `node` (with query degree `degree`).
    Seed { node: usize, degree: usize },
    /// Extend from already-mapped `from` along `edge` to unmapped `to` (with query degree
    /// `degree_to`).
    Extend {
        from: usize,
        to: usize,
        degree_to: usize,
        edge: EdgeId,
    },
    /// Verify a ring-closing edge between two already-mapped query nodes.
    Close {
        first: usize,
        second: usize,
        edge: EdgeId,
    },
}

/// Algorithm-specific precompiled query plan stored inside [`PreparedSubgraphQuery`].
#[derive(Clone, Debug)]
#[allow(dead_code, reason = "consumed by matching engines")]
pub(crate) enum PreparedQueryPlan {
    Vf2Rdkit(StaticOrderPlan),
    RayKirsch(Vec<RkStep>),
    /// Placeholder for algorithms compiled dynamically per target (`BitsetRi`, `Glasgow`) or
    /// legacy variants scheduled for removal.
    Dynamic,
}

/// Immutable, algorithm-tagged prepared query graph for repeated subgraph-isomorphism matching.
///
/// Holds an `Arc`-shared snapshot of the query [`Graph`], its selected
/// [`SubgraphIsomorphismAlgorithm`], degree metadata, algorithm-specific plan, and optional
/// [`QuerySymmetryConstraints`]. Mutating the source graph after construction cannot affect the
/// prepared query.
///
/// # Semantic properties
///
/// Matching a `PreparedSubgraphQuery` against a [`PreparedSubgraphTarget`] produces the exact
/// same embedding set, existence, and count as matching the underlying graphs directly with the
/// same algorithm, predicate, and symmetry constraints.
#[derive(Clone, Debug)]
pub struct PreparedSubgraphQuery {
    graph: Graph,
    algorithm: SubgraphIsomorphismAlgorithm,
    degree_hist: Vec<u32>,
    degrees: Vec<u32>,
    priority_order: Vec<u32>,
    priority_rank: Vec<u32>,
    symmetry: QuerySymmetryConstraints,
    /// Flat CSR of strict lower-bound symmetry predecessors `lo` (`phi(lo) < phi(q)`) per `q`.
    sym_gt_offsets: Vec<u32>,
    sym_gt_nodes: Vec<u32>,
    /// Flat CSR of strict upper-bound symmetry successors `hi` (`phi(q) < phi(hi)`) per `q`.
    sym_lt_offsets: Vec<u32>,
    sym_lt_nodes: Vec<u32>,
    plan: PreparedQueryPlan,
}

impl PreparedSubgraphQuery {
    /// Prepare `query` for subgraph-isomorphism search with `algorithm` and trivial query
    /// symmetry (`{1}`).
    pub fn new(query: &Graph, algorithm: SubgraphIsomorphismAlgorithm) -> Self {
        let graph = query.clone();
        let node_count = graph.node_count();
        let degree_hist = build_degree_hist(&graph);
        let degrees: Vec<u32> = graph
            .node_ids()
            .map(|node| graph.degree(node) as u32)
            .collect();
        let priority_order = query_priority_order(&graph);
        let mut priority_rank = vec![0u32; node_count];
        for (rank, &node) in priority_order.iter().enumerate() {
            priority_rank[node as usize] = rank as u32;
        }
        let plan = match algorithm {
            SubgraphIsomorphismAlgorithm::Vf2Rdkit => {
                PreparedQueryPlan::Vf2Rdkit(build_vf2rdkit_plan(&graph, None))
            }
            SubgraphIsomorphismAlgorithm::RayKirsch => {
                PreparedQueryPlan::RayKirsch(build_ray_kirsch_plan(&graph, None))
            }
            _ => PreparedQueryPlan::Dynamic,
        };
        let symmetry = QuerySymmetryConstraints::none(node_count);
        let empty_offsets = vec![0u32; node_count + 1];

        Self {
            graph,
            algorithm,
            degree_hist,
            degrees,
            priority_order,
            priority_rank,
            symmetry,
            sym_gt_offsets: empty_offsets.clone(),
            sym_gt_nodes: Vec::new(),
            sym_lt_offsets: empty_offsets,
            sym_lt_nodes: Vec::new(),
            plan,
        }
    }

    /// Attach query symmetry-breaking constraints to this prepared query.
    ///
    /// # Panics
    ///
    /// Panics if `symmetry` is incompatible with this query's node count or degree-priority order.
    /// Use [`Self::try_with_symmetry`] for independently constructed inputs.
    pub fn with_symmetry(self, symmetry: QuerySymmetryConstraints) -> Self {
        self.try_with_symmetry(symmetry)
            .expect("query symmetry constraints must match the prepared query graph")
    }

    /// Checked form of [`Self::with_symmetry`].
    ///
    /// Returns `None` if `symmetry.node_count()` differs from `self.node_count()` or if non-trivial
    /// constraints were constructed for a different vertex degree-priority order.
    pub fn try_with_symmetry(mut self, symmetry: QuerySymmetryConstraints) -> Option<Self> {
        let node_count = self.graph.node_count();
        if symmetry.node_count() != node_count {
            return None;
        }
        if !symmetry.is_trivial() && symmetry.priority_order() != self.priority_order.as_slice() {
            return None;
        }
        let (sym_gt_offsets, sym_gt_nodes, sym_lt_offsets, sym_lt_nodes) =
            build_symmetry_csr(node_count, symmetry.constraints());
        self.symmetry = symmetry;
        self.sym_gt_offsets = sym_gt_offsets;
        self.sym_gt_nodes = sym_gt_nodes;
        self.sym_lt_offsets = sym_lt_offsets;
        self.sym_lt_nodes = sym_lt_nodes;
        Some(self)
    }

    /// Snapshotted query [`Graph`].
    pub fn graph(&self) -> &Graph {
        &self.graph
    }

    /// Selected [`SubgraphIsomorphismAlgorithm`].
    pub fn algorithm(&self) -> SubgraphIsomorphismAlgorithm {
        self.algorithm
    }

    /// Number of nodes in the query graph.
    pub fn node_count(&self) -> usize {
        self.graph.node_count()
    }

    /// Number of edges in the query graph.
    pub fn edge_count(&self) -> usize {
        self.graph.edge_count()
    }

    /// Maximum vertex degree in the query graph (`0` for empty or edgeless graphs).
    pub fn max_degree(&self) -> usize {
        self.graph.max_degree()
    }

    /// Attached query symmetry-breaking constraints.
    pub fn symmetry(&self) -> &QuerySymmetryConstraints {
        &self.symmetry
    }

    #[allow(dead_code, reason = "consumed by matching engines")]
    pub(crate) fn degree_hist(&self) -> &[u32] {
        &self.degree_hist
    }

    #[allow(dead_code, reason = "consumed by matching engines")]
    pub(crate) fn degrees(&self) -> &[u32] {
        &self.degrees
    }

    #[allow(dead_code, reason = "consumed by matching engines")]
    pub(crate) fn priority_order(&self) -> &[u32] {
        &self.priority_order
    }

    #[allow(dead_code, reason = "consumed by matching engines")]
    pub(crate) fn priority_rank(&self) -> &[u32] {
        &self.priority_rank
    }

    /// Query nodes `lo` constrained to satisfy `phi(lo) < phi(node)`.
    #[allow(dead_code, reason = "consumed by matching engines")]
    pub(crate) fn sym_lower_bounds(&self, node: usize) -> &[u32] {
        let start = self.sym_gt_offsets[node] as usize;
        let end = self.sym_gt_offsets[node + 1] as usize;
        &self.sym_gt_nodes[start..end]
    }

    /// Query nodes `hi` constrained to satisfy `phi(node) < phi(hi)`.
    #[allow(dead_code, reason = "consumed by matching engines")]
    pub(crate) fn sym_upper_bounds(&self, node: usize) -> &[u32] {
        let start = self.sym_lt_offsets[node] as usize;
        let end = self.sym_lt_offsets[node + 1] as usize;
        &self.sym_lt_nodes[start..end]
    }

    #[allow(dead_code, reason = "consumed by matching engines")]
    pub(crate) fn plan(&self) -> &PreparedQueryPlan {
        &self.plan
    }
}

fn build_degree_hist(graph: &Graph) -> Vec<u32> {
    let max_deg = graph.max_degree();
    let mut hist = vec![0u32; max_deg + 1];
    for node in graph.node_ids() {
        hist[graph.degree(node)] += 1;
    }
    for d in (0..max_deg).rev() {
        hist[d] += hist[d + 1];
    }
    hist
}

fn build_target_bitsets(graph: &Graph) -> TargetBitsets {
    let n = graph.node_count();
    let max_deg = graph.max_degree();
    if n <= 64 {
        let mut deg_ge = vec![0u64; max_deg + 1];
        let mut adj = vec![0u64; n];
        for (v, adj_row) in adj.iter_mut().enumerate() {
            let node = NodeId(v as u32);
            deg_ge[graph.degree(node)] |= 1u64 << v;
            let mut mask = 0u64;
            for nb in graph.neighbors(node) {
                mask |= 1u64 << nb.node.index();
            }
            *adj_row = mask;
        }
        for d in (0..max_deg).rev() {
            deg_ge[d] |= deg_ge[d + 1];
        }
        TargetBitsets::U64 { deg_ge, adj }
    } else if n <= 128 {
        let mut deg_ge = vec![0u128; max_deg + 1];
        let mut adj = vec![0u128; n];
        for (v, adj_row) in adj.iter_mut().enumerate() {
            let node = NodeId(v as u32);
            deg_ge[graph.degree(node)] |= 1u128 << v;
            let mut mask = 0u128;
            for nb in graph.neighbors(node) {
                mask |= 1u128 << nb.node.index();
            }
            *adj_row = mask;
        }
        for d in (0..max_deg).rev() {
            deg_ge[d] |= deg_ge[d + 1];
        }
        TargetBitsets::U128 { deg_ge, adj }
    } else {
        let words = n.div_ceil(64);
        let mut deg_ge = vec![0u64; (max_deg + 1) * words];
        for v in 0..n {
            let deg = graph.degree(NodeId(v as u32));
            deg_ge[deg * words + v / 64] |= 1u64 << (v % 64);
        }
        for d in (0..max_deg).rev() {
            let (lo, hi) = deg_ge.split_at_mut((d + 1) * words);
            let dst = &mut lo[d * words..(d + 1) * words];
            let src = &hi[..words];
            for (a, &b) in dst.iter_mut().zip(src.iter()) {
                *a |= b;
            }
        }
        let adj = (n <= WIDE_EAGER_ADJ_LIMIT).then(|| {
            let mut rows = vec![0u64; n * words];
            for v in 0..n {
                let row = &mut rows[v * words..(v + 1) * words];
                for nb in graph.neighbors(NodeId(v as u32)) {
                    let u = nb.node.index();
                    row[u / 64] |= 1u64 << (u % 64);
                }
            }
            rows
        });
        TargetBitsets::Wide { words, deg_ge, adj }
    }
}

fn build_symmetry_csr(
    node_count: usize,
    constraints: &[(NodeId, NodeId)],
) -> (Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>) {
    let mut gt_deg = vec![0u32; node_count];
    let mut lt_deg = vec![0u32; node_count];
    for &(lo, hi) in constraints {
        gt_deg[hi.index()] += 1;
        lt_deg[lo.index()] += 1;
    }

    let mut gt_offsets = Vec::with_capacity(node_count + 1);
    let mut lt_offsets = Vec::with_capacity(node_count + 1);
    gt_offsets.push(0);
    lt_offsets.push(0);
    for i in 0..node_count {
        gt_offsets.push(gt_offsets[i] + gt_deg[i]);
        lt_offsets.push(lt_offsets[i] + lt_deg[i]);
    }

    let mut gt_nodes = vec![0u32; constraints.len()];
    let mut lt_nodes = vec![0u32; constraints.len()];
    let mut gt_cursor = gt_offsets[..node_count].to_vec();
    let mut lt_cursor = lt_offsets[..node_count].to_vec();
    for &(lo, hi) in constraints {
        let hi_idx = hi.index();
        gt_nodes[gt_cursor[hi_idx] as usize] = lo.0;
        gt_cursor[hi_idx] += 1;

        let lo_idx = lo.index();
        lt_nodes[lt_cursor[lo_idx] as usize] = hi.0;
        lt_cursor[lo_idx] += 1;
    }

    (gt_offsets, gt_nodes, lt_offsets, lt_nodes)
}

/// Build the static BFS/index-0 vertex ordering and flat CSR parent table for `Vf2Rdkit`.
pub(crate) fn build_vf2rdkit_plan(query: &Graph, anchor_q: Option<usize>) -> StaticOrderPlan {
    let n1 = query.node_count();
    let mut order = Vec::with_capacity(n1);
    let mut placed = vec![false; n1];
    if let Some(qa) = anchor_q {
        order.push(qa);
        placed[qa] = true;
    }
    while order.len() < n1 {
        let mut fallback = None;
        let mut chosen = None;
        for q in 0..n1 {
            if placed[q] {
                continue;
            }
            if fallback.is_none() {
                fallback = Some(q);
            }
            if query
                .neighbors(NodeId(q as u32))
                .iter()
                .any(|nb| placed[nb.node.index()])
            {
                chosen = Some(q);
                break;
            }
        }
        let next = chosen.or(fallback).expect("unplaced node exists");
        order.push(next);
        placed[next] = true;
    }

    let mut pos = vec![0usize; n1];
    for (d, &q) in order.iter().enumerate() {
        pos[q] = d;
    }
    let mut parents_offsets = Vec::with_capacity(n1 + 1);
    let mut parents_data = Vec::with_capacity(query.edge_count());
    parents_offsets.push(0);
    for (d, &q) in order.iter().enumerate() {
        for nb in query.neighbors(NodeId(q as u32)) {
            let p = pos[nb.node.index()];
            if p < d {
                parents_data.push((p, nb.edge));
            }
        }
        parents_offsets.push(parents_data.len());
    }

    StaticOrderPlan {
        order,
        parents_offsets,
        parents_data,
    }
}

/// Build the RI-guided `Seed`/`Extend`/`Close` step program for `RayKirsch`.
pub(crate) fn build_ray_kirsch_plan(query: &Graph, anchor_q: Option<usize>) -> Vec<RkStep> {
    let n1 = query.node_count();
    if n1 == 0 {
        return Vec::new();
    }

    let mut order = Vec::with_capacity(n1);
    let mut placed = vec![false; n1];
    let mut ordered_neighbors = vec![0usize; n1];

    let start = anchor_q.unwrap_or_else(|| {
        (0..n1)
            .min_by_key(|&u| (Reverse(query.degree(NodeId(u as u32))), u))
            .expect("non-empty query")
    });
    order.push(start);
    placed[start] = true;
    for nb in query.neighbors(NodeId(start as u32)) {
        ordered_neighbors[nb.node.index()] += 1;
    }

    while order.len() < n1 {
        let mut best = None;
        for u in 0..n1 {
            if placed[u] {
                continue;
            }
            let (mut vm, mut vn, mut vo) = (0usize, 0usize, 0usize);
            for nb in query.neighbors(NodeId(u as u32)) {
                let x = nb.node.index();
                if placed[x] {
                    vm += 1;
                } else if ordered_neighbors[x] > 0 {
                    vn += 1;
                } else {
                    vo += 1;
                }
            }
            let key = (vm, vn, vo, Reverse(u));
            if best.is_none_or(|(best_key, _)| key > best_key) {
                best = Some((key, u));
            }
        }
        let next = best.expect("unplaced node exists").1;
        order.push(next);
        placed[next] = true;
        for nb in query.neighbors(NodeId(next as u32)) {
            ordered_neighbors[nb.node.index()] += 1;
        }
    }

    let mut pos = vec![0usize; n1];
    for (d, &u) in order.iter().enumerate() {
        pos[u] = d;
    }

    let mut steps = Vec::with_capacity(n1 + query.edge_count());
    for (d, &u) in order.iter().enumerate() {
        let mut prev_edges = Vec::new();
        for nb in query.neighbors(NodeId(u as u32)) {
            let p = pos[nb.node.index()];
            if p < d {
                prev_edges.push((nb.node.index(), nb.edge));
            }
        }
        let degree = query.degree(NodeId(u as u32));
        if prev_edges.is_empty() {
            steps.push(RkStep::Seed { node: u, degree });
        } else {
            let (from, edge) = prev_edges[0];
            steps.push(RkStep::Extend {
                from,
                to: u,
                degree_to: degree,
                edge,
            });
            for &(prev, close_edge) in &prev_edges[1..] {
                steps.push(RkStep::Close {
                    first: prev,
                    second: u,
                    edge: close_edge,
                });
            }
        }
    }
    steps
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;

    use super::*;
    use crate::algorithms::automorphism::AutomorphismAlgorithm;

    #[rstest]
    fn test_subgraph_isomorphism_scratch_new() {
        let scratch = SubgraphIsomorphismScratch::new();
        assert!(scratch.order.is_empty());
        assert!(scratch.mapping.is_empty());
        let cloned = scratch.clone();
        assert!(cloned.out_nodes.is_empty());
    }

    #[rstest]
    fn test_prepared_subgraph_target_new() {
        // Empty target.
        let empty = PreparedSubgraphTarget::new(&Graph::new(0, &[]));
        assert_eq!(empty.node_count(), 0);
        assert_eq!(empty.edge_count(), 0);
        assert_eq!(empty.max_degree(), 0);
        assert_eq!(empty.degree_hist(), &[0]);

        // Small target (<= 64): path 0-1-2 plus isolated node 3.
        let mut source = Graph::new(4, &[[0, 1], [1, 2]]);
        let prepared = PreparedSubgraphTarget::new(&source);
        // Mutating the source graph after preparation must not affect the snapshot.
        source.add_edge(NodeId(0), NodeId(3));
        assert_eq!(prepared.node_count(), 4);
        assert_eq!(prepared.edge_count(), 2);
        assert_eq!(prepared.max_degree(), 2);
        assert_eq!(prepared.degree_hist(), &[4, 3, 1]);
        match prepared.bitsets() {
            TargetBitsets::U64 { deg_ge, adj } => {
                assert_eq!(deg_ge, &[0b1111, 0b0111, 0b0010]);
                assert_eq!(adj, &[0b0010, 0b0101, 0b0010, 0b0000]);
            }
            other => panic!("expected U64 bitsets, got {other:?}"),
        }

        // Medium target (65..=128): path of 70 nodes.
        let edges_70: Vec<[u32; 2]> = (0..69).map(|i| [i, i + 1]).collect();
        let med = PreparedSubgraphTarget::new(&Graph::new(70, &edges_70));
        assert_eq!(med.degree_hist(), &[70, 70, 68]);
        match med.bitsets() {
            TargetBitsets::U128 { deg_ge, adj } => {
                assert_eq!(deg_ge[0], (1u128 << 70) - 1);
                assert_eq!(deg_ge[2], ((1u128 << 69) - 1) ^ 1u128);
                assert_eq!(adj[0], 1u128 << 1);
                assert_eq!(adj[69], 1u128 << 68);
                assert_eq!(adj[64], (1u128 << 63) | (1u128 << 65));
            }
            other => panic!("expected U128 bitsets, got {other:?}"),
        }

        // Wide target with eager adj (129..=1024): 144 nodes, single edge (0, 130).
        let wide = PreparedSubgraphTarget::new(&Graph::new(144, &[[0, 130]]));
        assert_eq!(wide.degree_hist(), &[144, 2]);
        match wide.bitsets() {
            TargetBitsets::Wide { words, deg_ge, adj } => {
                assert_eq!(*words, 3);
                assert_eq!(deg_ge[3], 1u64);
                assert_eq!(deg_ge[4], 0u64);
                assert_eq!(deg_ge[5], 1u64 << 2);
                let adj_rows = adj.as_ref().expect("n=144 populates eager adj");
                assert_eq!(adj_rows[2], 1u64 << 2);
                assert_eq!(adj_rows[130 * 3], 1u64);
            }
            other => panic!("expected Wide bitsets, got {other:?}"),
        }

        // Very large sparse target (> WIDE_EAGER_ADJ_LIMIT): skips eager adj matrix.
        let huge = PreparedSubgraphTarget::new(&Graph::new(1050, &[[0, 1049]]));
        match huge.bitsets() {
            TargetBitsets::Wide { words, deg_ge, adj } => {
                assert_eq!(*words, 17);
                assert!(adj.is_none());
                assert_eq!(deg_ge[17], 1u64);
                assert_eq!(deg_ge[17 + 16], 1u64 << (1049 % 64));
            }
            other => panic!("expected Wide bitsets, got {other:?}"),
        }
    }

    #[rstest]
    #[case::empty_in_empty(Graph::new(0, &[]), Graph::new(0, &[]), true)]
    #[case::path_in_cycle(Graph::new(3, &[[0, 1], [1, 2]]), Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]), true)]
    #[case::triangle_not_in_path(Graph::new(3, &[[0, 1], [1, 2], [0, 2]]), Graph::new(3, &[[0, 1], [1, 2]]), false)]
    #[case::star_exceeds_max_degree(Graph::new(4, &[[0, 1], [0, 2], [0, 3]]), Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]), false)]
    fn test_prepared_subgraph_target_satisfies_degree_sequence(
        #[case] query: Graph,
        #[case] target: Graph,
        #[case] expected: bool,
    ) {
        let pq = PreparedSubgraphQuery::new(&query, SubgraphIsomorphismAlgorithm::Vf2Rdkit);
        let pt = PreparedSubgraphTarget::new(&target);
        assert_eq!(pt.satisfies_degree_sequence(pq.degree_hist()), expected);
    }

    #[rstest]
    fn test_prepared_subgraph_query_plans_and_symmetry() {
        let mut source = Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]);
        let sym = QuerySymmetryConstraints::from_graph(
            &source,
            |_| 0u8,
            |_| 0u8,
            AutomorphismAlgorithm::TrailIr,
        );
        let pq_rk = PreparedSubgraphQuery::new(&source, SubgraphIsomorphismAlgorithm::RayKirsch)
            .with_symmetry(sym.clone());
        let pq_vf2 = PreparedSubgraphQuery::new(&source, SubgraphIsomorphismAlgorithm::Vf2Rdkit)
            .with_symmetry(sym);

        // Mutating source does not affect prepared snapshots.
        source.add_node();
        assert_eq!(pq_rk.node_count(), 4);
        assert_eq!(pq_rk.edge_count(), 4);
        assert_eq!(pq_rk.max_degree(), 2);
        assert_eq!(pq_rk.algorithm(), SubgraphIsomorphismAlgorithm::RayKirsch);
        assert_eq!(pq_vf2.algorithm(), SubgraphIsomorphismAlgorithm::Vf2Rdkit);

        // Symmetry CSR bounds for C_4 (D_4 constraints: (0,1), (0,2), (0,3), (1,3)):
        assert_eq!(pq_rk.sym_lower_bounds(0), &[]);
        assert_eq!(pq_rk.sym_upper_bounds(0), &[1, 2, 3]);
        assert_eq!(pq_rk.sym_lower_bounds(1), &[0]);
        assert_eq!(pq_rk.sym_upper_bounds(1), &[3]);
        assert_eq!(pq_rk.sym_lower_bounds(3), &[0, 1]);
        assert_eq!(pq_rk.sym_upper_bounds(3), &[]);

        // RayKirsch plan covers every vertex once (Seed/Extend) and every edge once (Extend/Close).
        match pq_rk.plan() {
            PreparedQueryPlan::RayKirsch(steps) => {
                assert_eq!(steps.len(), 5);
                assert_eq!(steps[0], RkStep::Seed { node: 0, degree: 2 });
                let closes = steps
                    .iter()
                    .filter(|s| matches!(s, RkStep::Close { .. }))
                    .count();
                assert_eq!(closes, 1);
            }
            other => panic!("expected RayKirsch plan, got {other:?}"),
        }

        // Vf2Rdkit plan covers all 4 vertices and all 4 parent edges.
        match pq_vf2.plan() {
            PreparedQueryPlan::Vf2Rdkit(plan) => {
                assert_eq!(plan.order.len(), 4);
                assert_eq!(plan.parents_data.len(), 4);
            }
            other => panic!("expected Vf2Rdkit plan, got {other:?}"),
        }
        // Empty and disconnected queries produce well-formed plans without panicking.
        assert!(build_ray_kirsch_plan(&Graph::new(0, &[]), None).is_empty());
        let disc = Graph::new(5, &[[0, 1], [2, 3]]);
        let disc_steps = build_ray_kirsch_plan(&disc, Some(4));
        assert_eq!(disc_steps[0], RkStep::Seed { node: 4, degree: 0 });
        assert_eq!(
            disc_steps
                .iter()
                .filter(|s| matches!(s, RkStep::Seed { .. }))
                .count(),
            3
        );
        let disc_vf2 = build_vf2rdkit_plan(&disc, Some(4));
        assert_eq!(disc_vf2.order[0], 4);
        assert_eq!(disc_vf2.order.len(), 5);
        assert_eq!(disc_vf2.parents_data.len(), 2);
    }

    #[rstest]
    fn test_prepared_subgraph_query_try_with_symmetry_error() {
        let path4 = Graph::new(4, &[[0, 1], [1, 2], [2, 3]]);
        let cycle4 = Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]);
        let triangle = Graph::new(3, &[[0, 1], [1, 2], [0, 2]]);

        let pq = PreparedSubgraphQuery::new(&path4, SubgraphIsomorphismAlgorithm::Vf2Rdkit);

        // Wrong node count rejected.
        let tri_sym = QuerySymmetryConstraints::from_graph(
            &triangle,
            |_| 0u8,
            |_| 0u8,
            AutomorphismAlgorithm::TrailIr,
        );
        assert!(pq.clone().try_with_symmetry(tri_sym).is_none());

        // Same node count but different degree-priority order rejected.
        let cycle_sym = QuerySymmetryConstraints::from_graph(
            &cycle4,
            |_| 0u8,
            |_| 0u8,
            AutomorphismAlgorithm::TrailIr,
        );
        assert!(pq.clone().try_with_symmetry(cycle_sym).is_none());

        // Trivial symmetry with matching node count accepted.
        assert!(pq
            .try_with_symmetry(QuerySymmetryConstraints::none(4))
            .is_some());
    }

    #[rstest]
    #[should_panic(expected = "query symmetry constraints must match the prepared query graph")]
    fn test_prepared_subgraph_query_with_symmetry_error() {
        let path4 = Graph::new(4, &[[0, 1], [1, 2], [2, 3]]);
        let pq = PreparedSubgraphQuery::new(&path4, SubgraphIsomorphismAlgorithm::Vf2Rdkit);
        let _ = pq.with_symmetry(QuerySymmetryConstraints::none(3));
    }
}
