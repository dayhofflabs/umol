//! Bit-parallel RI-DS subgraph-isomorphism engine (`BitsetRi`).
//!
//! Combines RI-DS rarity-and-connectivity query ordering (Bonnici & Giugno 2017) with
//! width-specialized `u64`/`u128`/wide bitsets, descending-degree candidate initialization
//! (with optional query-node class deduplication), lazy or precomputed target adjacency, direct
//! sparse neighbor scanning for selective single-parent steps, 2-instruction bit-parallel
//! monomorphism lookahead on cyclic queries, push-down query/target symmetry breaking, and
//! leaf `POPCNT` counting when leaf edges are bitmask-resolved.

use std::cmp::Reverse;
use std::mem;
use std::ops::{BitAnd, BitOr, ControlFlow, Not, Shl};

use super::predicate::SubgraphIsomorphismPredicate;
use super::prepared::{
    PreparedSubgraphQuery, PreparedSubgraphTarget, SubgraphIsomorphismScratch, TargetBitsets,
};
use super::symmetry::{
    is_canonical_double_orbit_with_built_chain, QuerySymmetryConstraints, TargetStabilizerChain,
    TargetSymmetry,
};
use crate::graph::Graph;
use crate::id::{EdgeId, NodeId};

/// Maximum target degree for direct sparse CSR neighbor scans on selective single-parent steps.
const SPARSE_SINGLE_PARENT_DEG_LIMIT: usize = 8;

/// Maximum target node count for lazily allocating a dense `n2 * words` wide adjacency matrix in
/// scratch; larger sparse targets use `Graph::find_edge` directly to avoid quadratic memory.
const WIDE_LAZY_ADJ_LIMIT: usize = 4096;

/// Stack-allocated parent-edge buffer threshold for direct `u64`/`u128` calls.
const STACK_PARENTS_LIMIT: usize = 128;

trait WordMask:
    Copy
    + Eq
    + BitAnd<Output = Self>
    + BitOr<Output = Self>
    + Not<Output = Self>
    + Shl<usize, Output = Self>
{
    const ZERO: Self;

    fn bit(index: usize) -> Self;
    fn all_below(count: usize) -> Self;
    fn range_mask(min_v: usize, max_v: usize) -> Self;
    fn is_zero(self) -> bool;
    fn count_ones(self) -> u32;
    fn trailing_zeros(self) -> u32;
    fn clear_lowest(self) -> Self;
    fn test(self, index: usize) -> bool;
    fn orbit_reps(target_sym: &TargetSymmetry) -> Self;
    fn node_candidates<P: SubgraphIsomorphismPredicate>(
        pred: &mut P,
        query_node: NodeId,
        eligible: Self,
    ) -> Self;
    fn edge_neighbors<P: SubgraphIsomorphismPredicate>(
        pred: &mut P,
        query_edge: EdgeId,
        target_node: NodeId,
    ) -> Option<Self>;
}

impl WordMask for u64 {
    const ZERO: Self = 0;

    #[inline(always)]
    fn bit(index: usize) -> Self {
        1u64 << index
    }

    #[inline(always)]
    fn all_below(count: usize) -> Self {
        if count == 0 {
            0
        } else if count >= 64 {
            u64::MAX
        } else {
            (1u64 << count) - 1
        }
    }

    #[inline(always)]
    fn range_mask(min_v: usize, max_v: usize) -> Self {
        if min_v >= max_v || min_v >= 64 {
            0
        } else {
            Self::all_below(max_v.min(64)) & (u64::MAX << min_v)
        }
    }

    #[inline(always)]
    fn is_zero(self) -> bool {
        self == 0
    }

    #[inline(always)]
    fn count_ones(self) -> u32 {
        u64::count_ones(self)
    }

    #[inline(always)]
    fn trailing_zeros(self) -> u32 {
        u64::trailing_zeros(self)
    }

    #[inline(always)]
    fn clear_lowest(self) -> Self {
        self & self.wrapping_sub(1)
    }

    #[inline(always)]
    fn test(self, index: usize) -> bool {
        ((self >> index) & 1) != 0
    }

    #[inline(always)]
    fn orbit_reps(target_sym: &TargetSymmetry) -> Self {
        target_sym.orbit_reps_u64()
    }

    #[inline(always)]
    fn node_candidates<P: SubgraphIsomorphismPredicate>(
        pred: &mut P,
        query_node: NodeId,
        eligible: Self,
    ) -> Self {
        pred.node_candidates_u64(query_node, eligible)
    }

    #[inline(always)]
    fn edge_neighbors<P: SubgraphIsomorphismPredicate>(
        pred: &mut P,
        query_edge: EdgeId,
        target_node: NodeId,
    ) -> Option<Self> {
        pred.edge_neighbors_u64(query_edge, target_node)
    }
}

impl WordMask for u128 {
    const ZERO: Self = 0;

    #[inline(always)]
    fn bit(index: usize) -> Self {
        1u128 << index
    }

    #[inline(always)]
    fn all_below(count: usize) -> Self {
        if count == 0 {
            0
        } else if count >= 128 {
            u128::MAX
        } else {
            (1u128 << count) - 1
        }
    }

    #[inline(always)]
    fn range_mask(min_v: usize, max_v: usize) -> Self {
        if min_v >= max_v || min_v >= 128 {
            0
        } else {
            Self::all_below(max_v.min(128)) & (u128::MAX << min_v)
        }
    }

    #[inline(always)]
    fn is_zero(self) -> bool {
        self == 0
    }

    #[inline(always)]
    fn count_ones(self) -> u32 {
        u128::count_ones(self)
    }

    #[inline(always)]
    fn trailing_zeros(self) -> u32 {
        u128::trailing_zeros(self)
    }

    #[inline(always)]
    fn clear_lowest(self) -> Self {
        self & self.wrapping_sub(1)
    }

    #[inline(always)]
    fn test(self, index: usize) -> bool {
        ((self >> index) & 1) != 0
    }

    #[inline(always)]
    fn orbit_reps(target_sym: &TargetSymmetry) -> Self {
        target_sym.orbit_reps_u128()
    }

    #[inline(always)]
    fn node_candidates<P: SubgraphIsomorphismPredicate>(
        pred: &mut P,
        query_node: NodeId,
        eligible: Self,
    ) -> Self {
        pred.node_candidates_u128(query_node, eligible)
    }

    #[inline(always)]
    fn edge_neighbors<P: SubgraphIsomorphismPredicate>(
        pred: &mut P,
        query_edge: EdgeId,
        target_node: NodeId,
    ) -> Option<Self> {
        pred.edge_neighbors_u128(query_edge, target_node)
    }
}

trait WordAdj<W: WordMask> {
    fn row(&mut self, u: usize) -> W;
}

struct PrecomputedWordAdj<'a, W>(&'a [W]);

impl<W: WordMask> WordAdj<W> for PrecomputedWordAdj<'_, W> {
    #[inline(always)]
    fn row(&mut self, u: usize) -> W {
        self.0[u]
    }
}

struct LazyWordAdj<'a, W> {
    target: &'a Graph,
    adj: &'a mut [W],
    built: W,
}

impl<W: WordMask> WordAdj<W> for LazyWordAdj<'_, W> {
    #[inline(always)]
    fn row(&mut self, u: usize) -> W {
        let bit = W::bit(u);
        if (self.built & bit).is_zero() {
            let mut mask = W::ZERO;
            for nb in self.target.neighbors(NodeId(u as u32)) {
                mask = mask | W::bit(nb.node.index());
            }
            self.adj[u] = mask;
            self.built = self.built | bit;
            mask
        } else {
            self.adj[u]
        }
    }
}

enum WideAdj<'a> {
    Precomputed(&'a [u64]),
    Lazy {
        adj: &'a mut [u64],
        built: &'a mut [u64],
    },
    None,
}

impl WideAdj<'_> {
    #[inline(always)]
    fn ensure_row(&mut self, target: &Graph, u: usize, words: usize) {
        if let Self::Lazy { adj, built } = self {
            let bw = u >> 6;
            let bm = 1u64 << (u & 63);
            if (built[bw] & bm) == 0 {
                let base = u * words;
                adj[base..base + words].fill(0);
                for nb in target.neighbors(NodeId(u as u32)) {
                    let v = nb.node.index();
                    adj[base + (v >> 6)] |= 1u64 << (v & 63);
                }
                built[bw] |= bm;
            }
        }
    }

    #[inline(always)]
    fn has_edge_bit(&self, u: usize, vw: usize, vb: u64, words: usize) -> Option<bool> {
        match self {
            Self::Precomputed(adj) => Some((adj[u * words + vw] & vb) != 0),
            Self::Lazy { adj, .. } => Some((adj[u * words + vw] & vb) != 0),
            Self::None => None,
        }
    }
}

/// Active push-down symmetry state for a single `BitsetRi` search.
///
/// Only constructed when at least one of query symmetry (`G_Q`) or target symmetry (`G_H`) is
/// non-trivial, so asymmetric searches incur a single predicted `None` branch per depth.
struct SearchSymmetry<'a> {
    has_q_sym: bool,
    has_h_sym: bool,
    sym_gt_offsets: Vec<usize>,
    sym_gt_depths: Vec<usize>,
    sym_lt_offsets: Vec<usize>,
    sym_lt_depths: Vec<usize>,
    priority_order: &'a [u32],
    priority_pos: Vec<usize>,
    prefix_len_after: Vec<usize>,
    query_sym: Option<&'a QuerySymmetryConstraints>,
    target_sym: Option<&'a TargetSymmetry>,
    target_chain: &'a mut TargetStabilizerChain,
    by_query: Vec<usize>,
}

impl<'a> SearchSymmetry<'a> {
    fn build(
        prepared_query: Option<&'a PreparedSubgraphQuery>,
        target_sym: Option<&'a TargetSymmetry>,
        order: &[usize],
        pos: &[usize],
        target_chain: &'a mut TargetStabilizerChain,
    ) -> Option<Self> {
        let query_sym = prepared_query
            .map(PreparedSubgraphQuery::symmetry)
            .filter(|sym| !sym.is_trivial());
        let active_target_sym = target_sym.filter(|sym| !sym.is_trivial());
        let has_q_sym = query_sym.is_some();
        let has_h_sym = active_target_sym.is_some();
        if !has_q_sym && !has_h_sym {
            return None;
        }

        let n1 = order.len();
        let mut sym_gt_offsets = Vec::new();
        let mut sym_gt_depths = Vec::new();
        let mut sym_lt_offsets = Vec::new();
        let mut sym_lt_depths = Vec::new();

        if let Some(pq) = prepared_query {
            if has_q_sym {
                sym_gt_offsets.reserve(n1 + 1);
                sym_lt_offsets.reserve(n1 + 1);
                sym_gt_offsets.push(0);
                sym_lt_offsets.push(0);
                for (d, &q) in order.iter().enumerate() {
                    for &lo in pq.sym_lower_bounds(q) {
                        let p = pos[lo as usize];
                        if p < d {
                            sym_gt_depths.push(p);
                        }
                    }
                    for &hi in pq.sym_upper_bounds(q) {
                        let p = pos[hi as usize];
                        if p < d {
                            sym_lt_depths.push(p);
                        }
                    }
                    sym_gt_offsets.push(sym_gt_depths.len());
                    sym_lt_offsets.push(sym_lt_depths.len());
                }
            }
        }

        let priority_order = prepared_query
            .map(PreparedSubgraphQuery::priority_order)
            .unwrap_or(&[]);
        let mut priority_pos = Vec::new();
        let mut prefix_len_after = Vec::new();
        if let Some(ts) = active_target_sym {
            debug_assert_eq!(priority_order.len(), n1);
            target_chain.reset(ts);
            priority_pos.reserve(n1);
            for &q in priority_order {
                priority_pos.push(pos[q as usize]);
            }
            prefix_len_after.reserve(n1);
            let mut k = 0;
            for d in 0..n1 {
                while k < n1 && priority_pos[k] <= d {
                    k += 1;
                }
                prefix_len_after.push(k);
            }
        }

        let by_query = if has_q_sym && has_h_sym {
            vec![0usize; n1]
        } else {
            Vec::new()
        };

        Some(Self {
            has_q_sym,
            has_h_sym,
            sym_gt_offsets,
            sym_gt_depths,
            sym_lt_offsets,
            sym_lt_depths,
            priority_order,
            priority_pos,
            prefix_len_after,
            query_sym,
            target_sym: active_target_sym,
            target_chain,
            by_query,
        })
    }

    #[inline(always)]
    fn depth_bounds(&self, depth: usize, mapping: &[usize], n2: usize) -> (usize, usize) {
        if !self.has_q_sym {
            return (0, n2);
        }
        let mut min_v = 0usize;
        for &p in &self.sym_gt_depths[self.sym_gt_offsets[depth]..self.sym_gt_offsets[depth + 1]] {
            min_v = min_v.max(mapping[p] + 1);
        }
        let mut max_v = n2;
        for &p in &self.sym_lt_depths[self.sym_lt_offsets[depth]..self.sym_lt_offsets[depth + 1]] {
            max_v = max_v.min(mapping[p]);
        }
        (min_v, max_v)
    }

    #[inline(always)]
    fn can_leaf_popcnt(&self) -> bool {
        !(self.has_q_sym && self.has_h_sym)
            && (!self.has_h_sym || self.target_chain.is_current_trivial())
    }

    /// Fast pre-edge check when `order[depth]` is the next unstabilized priority vertex `p_{prev_k}`.
    #[inline(always)]
    fn rejects_candidate_before_edge_check(&self, depth: usize, v: usize) -> bool {
        if !self.has_h_sym || self.target_chain.is_current_trivial() {
            return false;
        }
        let prev_k = if depth == 0 {
            0
        } else {
            self.prefix_len_after[depth - 1]
        };
        if self.prefix_len_after[depth] > prev_k {
            !self.target_chain.is_orbit_min_at_current(v)
        } else {
            false
        }
    }

    /// Advance `target_chain` across all priority vertices completed by assigning `mapping[depth]`.
    ///
    /// Returns `Some(saved_depth)` on success (`usize::MAX` if no level was pushed), or `None` if
    /// pruned (leaving `target_chain` restored).
    #[inline]
    fn try_advance_host_chain(&mut self, depth: usize, mapping: &[usize]) -> Option<usize> {
        if !self.has_h_sym || self.target_chain.is_current_trivial() {
            return Some(usize::MAX);
        }
        let prev_k = if depth == 0 {
            0
        } else {
            self.prefix_len_after[depth - 1]
        };
        let new_k = self.prefix_len_after[depth];
        if new_k == prev_k {
            return Some(usize::MAX);
        }

        let n1 = self.priority_pos.len();
        let saved_depth = self.target_chain.depth();
        let mut pushed = false;
        for r in prev_k..new_k {
            if self.target_chain.is_current_trivial() {
                break;
            }
            let t_r = mapping[self.priority_pos[r]];
            if !self.has_q_sym && r + 1 == n1 {
                if !self.target_chain.is_orbit_min_at_current(t_r) {
                    if pushed {
                        self.target_chain.truncate(saved_depth);
                    }
                    return None;
                }
            } else if self.target_chain.try_push_orbit_min(t_r) {
                pushed = true;
            } else {
                if pushed {
                    self.target_chain.truncate(saved_depth);
                }
                return None;
            }
        }
        Some(if pushed { saved_depth } else { usize::MAX })
    }

    #[inline(always)]
    fn restore_host_chain(&mut self, saved_depth: usize) {
        if saved_depth != usize::MAX {
            self.target_chain.truncate(saved_depth);
        }
    }

    #[inline]
    fn accepts_complete_embedding(&mut self, order: &[usize], mapping: &[usize]) -> bool {
        if !(self.has_q_sym && self.has_h_sym) {
            return true;
        }
        for (d, &q) in order.iter().enumerate() {
            self.by_query[q] = mapping[d];
        }
        is_canonical_double_orbit_with_built_chain(
            &self.by_query,
            self.priority_order,
            self.query_sym
                .expect("query symmetry present when has_q_sym is true"),
            self.target_sym
                .expect("target symmetry present when has_h_sym is true"),
            self.target_chain,
        )
    }
}

/// Visit every subgraph isomorphism of `query` in `target` directly (without prepared carriers).
pub(crate) fn visit_direct<B, P: SubgraphIsomorphismPredicate>(
    target: &Graph,
    query: &Graph,
    anchor: Option<(NodeId, NodeId)>,
    pred: &mut P,
    visitor: &mut impl FnMut(&[NodeId]) -> ControlFlow<B>,
) -> ControlFlow<B> {
    let mut count = 0u64;
    solve_direct::<false, B, P>(target, query, anchor, pred, &mut count, visitor)
}

/// Count every subgraph isomorphism of `query` in `target` directly (without prepared carriers).
#[allow(dead_code, reason = "wired into Graph::count_subgraph_isomorphisms")]
pub(crate) fn count_direct<P: SubgraphIsomorphismPredicate>(
    target: &Graph,
    query: &Graph,
    anchor: Option<(NodeId, NodeId)>,
    pred: &mut P,
) -> u64 {
    let mut count = 0u64;
    let mut noop = |_: &[NodeId]| ControlFlow::<()>::Continue(());
    let _ = solve_direct::<true, (), P>(target, query, anchor, pred, &mut count, &mut noop);
    count
}

/// Visit every subgraph isomorphism of a prepared query in a prepared target.
#[allow(
    dead_code,
    reason = "wired into PreparedSubgraphQuery matching methods"
)]
pub(crate) fn visit_prepared<B, P: SubgraphIsomorphismPredicate>(
    target: &PreparedSubgraphTarget,
    query: &PreparedSubgraphQuery,
    target_sym: Option<&TargetSymmetry>,
    pred: &mut P,
    scratch: &mut SubgraphIsomorphismScratch,
    visitor: &mut impl FnMut(&[NodeId]) -> ControlFlow<B>,
) -> ControlFlow<B> {
    let mut count = 0u64;
    solve_prepared::<false, B, P>(
        target, query, target_sym, pred, scratch, &mut count, visitor,
    )
}

/// Count every subgraph isomorphism of a prepared query in a prepared target.
#[allow(
    dead_code,
    reason = "wired into PreparedSubgraphQuery matching methods"
)]
pub(crate) fn count_prepared<P: SubgraphIsomorphismPredicate>(
    target: &PreparedSubgraphTarget,
    query: &PreparedSubgraphQuery,
    target_sym: Option<&TargetSymmetry>,
    pred: &mut P,
    scratch: &mut SubgraphIsomorphismScratch,
) -> u64 {
    let mut count = 0u64;
    let mut noop = |_: &[NodeId]| ControlFlow::<()>::Continue(());
    let _ = solve_prepared::<true, (), P>(
        target, query, target_sym, pred, scratch, &mut count, &mut noop,
    );
    count
}

fn solve_direct<const COUNT_ONLY: bool, B, P: SubgraphIsomorphismPredicate>(
    target: &Graph,
    query: &Graph,
    anchor: Option<(NodeId, NodeId)>,
    pred: &mut P,
    count: &mut u64,
    visitor: &mut impl FnMut(&[NodeId]) -> ControlFlow<B>,
) -> ControlFlow<B> {
    let n1 = query.node_count();
    let n2 = target.node_count();
    if n1 == 0 {
        if anchor.is_some() {
            return ControlFlow::Continue(());
        }
        if COUNT_ONLY {
            *count += 1;
            return ControlFlow::Continue(());
        }
        return visitor(&[]);
    }
    if n1 > n2
        || query.edge_count() > target.edge_count()
        || query.max_degree() > target.max_degree()
    {
        return ControlFlow::Continue(());
    }
    if let Some((qa, ta)) = anchor {
        if qa.index() >= n1
            || ta.index() >= n2
            || query.degree(qa) > target.degree(ta)
            || !pred.match_node(qa, ta)
        {
            return ControlFlow::Continue(());
        }
    }
    if !degree_sequence_ok(query, target) {
        return ControlFlow::Continue(());
    }

    if n2 <= 64 {
        solve_word_direct::<u64, 64, COUNT_ONLY, B, P>(target, query, anchor, pred, count, visitor)
    } else if n2 <= 128 {
        solve_word_direct::<u128, 128, COUNT_ONLY, B, P>(
            target, query, anchor, pred, count, visitor,
        )
    } else {
        let mut scratch = SubgraphIsomorphismScratch::new();
        solve_wide::<COUNT_ONLY, B, P>(
            target,
            query,
            None,
            None,
            None,
            anchor,
            None,
            pred,
            &mut scratch,
            count,
            visitor,
        )
    }
}

fn solve_prepared<const COUNT_ONLY: bool, B, P: SubgraphIsomorphismPredicate>(
    prepared_target: &PreparedSubgraphTarget,
    prepared_query: &PreparedSubgraphQuery,
    target_sym: Option<&TargetSymmetry>,
    pred: &mut P,
    scratch: &mut SubgraphIsomorphismScratch,
    count: &mut u64,
    visitor: &mut impl FnMut(&[NodeId]) -> ControlFlow<B>,
) -> ControlFlow<B> {
    let target = prepared_target.graph();
    let query = prepared_query.graph();
    if let Some(ts) = target_sym {
        assert_eq!(
            ts.node_count(),
            target.node_count(),
            "target symmetry node count must match the prepared target graph"
        );
    }

    let n1 = query.node_count();
    let n2 = target.node_count();
    if n1 == 0 {
        if COUNT_ONLY {
            *count += 1;
            return ControlFlow::Continue(());
        }
        return visitor(&[]);
    }
    if n1 > n2
        || query.edge_count() > target.edge_count()
        || !prepared_target.satisfies_degree_sequence(prepared_query.degree_hist())
    {
        return ControlFlow::Continue(());
    }

    match prepared_target.bitsets() {
        TargetBitsets::U64 { deg_ge, adj } => solve_word_prepared::<u64, 64, COUNT_ONLY, B, P>(
            target,
            prepared_query,
            deg_ge,
            adj,
            target_sym,
            pred,
            scratch,
            count,
            visitor,
        ),
        TargetBitsets::U128 { deg_ge, adj } => solve_word_prepared::<u128, 128, COUNT_ONLY, B, P>(
            target,
            prepared_query,
            deg_ge,
            adj,
            target_sym,
            pred,
            scratch,
            count,
            visitor,
        ),
        TargetBitsets::Wide { words, deg_ge, adj } => solve_wide::<COUNT_ONLY, B, P>(
            target,
            query,
            Some(prepared_query),
            Some((deg_ge.as_slice(), *words)),
            adj.as_deref(),
            None,
            target_sym,
            pred,
            scratch,
            count,
            visitor,
        ),
    }
}

/// O(n1 + n2) necessary degree-sequence domination filter for direct calls.
#[inline]
fn degree_sequence_ok(query: &Graph, target: &Graph) -> bool {
    let mut q_hist = [0u32; 16];
    let mut t_hist = [0u32; 16];
    for q in query.node_ids() {
        q_hist[query.degree(q).min(15)] += 1;
    }
    for t in target.node_ids() {
        t_hist[target.degree(t).min(15)] += 1;
    }
    let mut q_ge = 0u32;
    let mut t_ge = 0u32;
    for d in (1..16).rev() {
        q_ge += q_hist[d];
        t_ge += t_hist[d];
        if t_ge < q_ge {
            return false;
        }
    }
    true
}

fn compute_lookahead_counts(
    query: &Graph,
    order: &[usize],
    pos: &[usize],
    q_term: &mut [u32],
    q_rem: &mut [u32],
) -> bool {
    let n1 = order.len();
    let mut in_frontier = [false; 128];
    let mut any_term = false;
    for d in 0..n1 {
        let q = order[d];
        let mut term = 0u32;
        let mut rem = 0u32;
        for nb in query.neighbors(NodeId(q as u32)) {
            let u = nb.node.index();
            if pos[u] > d {
                rem += 1;
                if in_frontier[u] {
                    term += 1;
                }
            }
        }
        q_term[d] = term;
        q_rem[d] = rem;
        if term > 0 {
            any_term = true;
        }
        in_frontier[q] = false;
        for nb in query.neighbors(NodeId(q as u32)) {
            let u = nb.node.index();
            if pos[u] > d {
                in_frontier[u] = true;
            }
        }
    }
    any_term
}

/// Build RI-DS rarity+connectivity ordering on stack slices for `n1 <= 128`.
fn build_word_ri_order(
    query: &Graph,
    q_deg: &[usize],
    cand_cnt: &[usize],
    anchor: Option<(NodeId, NodeId)>,
    order: &mut [usize],
    pos: &mut [usize],
) {
    let n1 = query.node_count();
    let mut placed = [false; 128];
    let mut conn = [0usize; 128];

    let start = if let Some((qa, _)) = anchor {
        qa.index()
    } else {
        let mut best = 0;
        for q in 1..n1 {
            if (cand_cnt[q], Reverse(q_deg[q]), q) < (cand_cnt[best], Reverse(q_deg[best]), best) {
                best = q;
            }
        }
        best
    };

    order[0] = start;
    pos[start] = 0;
    placed[start] = true;
    for nb in query.neighbors(NodeId(start as u32)) {
        conn[nb.node.index()] += 1;
    }

    for (step, slot) in order.iter_mut().enumerate().take(n1).skip(1) {
        let mut best = None;
        for q in 0..n1 {
            if placed[q] {
                continue;
            }
            let c = conn[q];
            let mut frontier_deg = 0usize;
            for nb in query.neighbors(NodeId(q as u32)) {
                let u = nb.node.index();
                if !placed[u] && conn[u] > 0 {
                    frontier_deg += 1;
                }
            }
            let key = (c, Reverse(cand_cnt[q]), frontier_deg, q_deg[q], Reverse(q));
            if best.is_none_or(|(best_key, _)| key > best_key) {
                best = Some((key, q));
            }
        }
        let next = best.expect("unplaced query vertex exists").1;
        *slot = next;
        pos[next] = step;
        placed[next] = true;
        for nb in query.neighbors(NodeId(next as u32)) {
            conn[nb.node.index()] += 1;
        }
    }
}

/// Populate per-node candidate bitmasks for `n2 <= 128` using descending-degree fail-fast and
/// optional query-node class deduplication. Returns `false` if any query node has an empty domain.
#[allow(clippy::too_many_arguments)]
fn init_word_candidates_by_class<W: WordMask, const N: usize, P: SubgraphIsomorphismPredicate>(
    n1: usize,
    q_deg: &[usize; N],
    q_by_deg: &[usize; N],
    anchor: Option<(NodeId, NodeId)>,
    anchor_clear: W,
    root_orbit_filter: Option<(usize, W)>,
    mut eligible_for_deg: impl FnMut(usize) -> W,
    pred: &mut P,
    cand: &mut [W; N],
    cand_cnt: &mut [usize; N],
    class_min_deg: &mut [usize],
    class_count: &mut [usize],
    class_computed: &mut [bool],
    class_masks: &mut [W],
) -> bool {
    for &q in &q_by_deg[..n1] {
        if anchor.is_some_and(|(qa, _)| q == qa.index()) {
            continue;
        }
        let c = pred.query_node_class(NodeId(q as u32));
        debug_assert!(c < class_min_deg.len());
        class_min_deg[c] = class_min_deg[c].min(q_deg[q]);
        class_count[c] += 1;
    }

    for &q in &q_by_deg[..n1] {
        if let Some((qa, ta)) = anchor {
            if q == qa.index() {
                cand[q] = W::bit(ta.index());
                cand_cnt[q] = 1;
                continue;
            }
        }
        let c = pred.query_node_class(NodeId(q as u32));
        if !class_computed[c] {
            let mut eligible = eligible_for_deg(class_min_deg[c]) & !anchor_clear;
            if class_count[c] == 1 {
                if let Some((p0, orbit_reps)) = root_orbit_filter {
                    if q == p0 {
                        eligible = eligible & orbit_reps;
                    }
                }
            }
            class_masks[c] = if eligible.is_zero() {
                W::ZERO
            } else {
                W::node_candidates(pred, NodeId(q as u32), eligible)
            };
            class_computed[c] = true;
        }
        let mut mask = class_masks[c] & eligible_for_deg(q_deg[q]);
        if let Some((p0, orbit_reps)) = root_orbit_filter {
            if q == p0 {
                mask = mask & orbit_reps;
            }
        }
        if mask.is_zero() {
            return false;
        }
        cand[q] = mask;
        cand_cnt[q] = mask.count_ones() as usize;
    }
    true
}

#[allow(clippy::too_many_arguments)]
fn init_word_candidates<W: WordMask, const N: usize, P: SubgraphIsomorphismPredicate>(
    n1: usize,
    q_deg: &[usize; N],
    q_by_deg: &[usize; N],
    anchor: Option<(NodeId, NodeId)>,
    root_orbit_filter: Option<(usize, W)>,
    mut eligible_for_deg: impl FnMut(usize) -> W,
    pred: &mut P,
    cand: &mut [W; N],
    cand_cnt: &mut [usize; N],
) -> bool {
    let anchor_clear = anchor.map_or(W::ZERO, |(_, ta)| W::bit(ta.index()));

    if let Some(num_classes) = pred.query_node_class_count() {
        if num_classes <= N {
            let mut class_min_deg = [usize::MAX; N];
            let mut class_count = [0usize; N];
            let mut class_computed = [false; N];
            let mut class_masks = [W::ZERO; N];
            return init_word_candidates_by_class::<W, N, P>(
                n1,
                q_deg,
                q_by_deg,
                anchor,
                anchor_clear,
                root_orbit_filter,
                eligible_for_deg,
                pred,
                cand,
                cand_cnt,
                &mut class_min_deg[..num_classes],
                &mut class_count[..num_classes],
                &mut class_computed[..num_classes],
                &mut class_masks[..num_classes],
            );
        }
        let mut class_min_deg = vec![usize::MAX; num_classes];
        let mut class_count = vec![0usize; num_classes];
        let mut class_computed = vec![false; num_classes];
        let mut class_masks = vec![W::ZERO; num_classes];
        return init_word_candidates_by_class::<W, N, P>(
            n1,
            q_deg,
            q_by_deg,
            anchor,
            anchor_clear,
            root_orbit_filter,
            eligible_for_deg,
            pred,
            cand,
            cand_cnt,
            &mut class_min_deg,
            &mut class_count,
            &mut class_computed,
            &mut class_masks,
        );
    }

    for &q in &q_by_deg[..n1] {
        if let Some((qa, ta)) = anchor {
            if q == qa.index() {
                cand[q] = W::bit(ta.index());
                cand_cnt[q] = 1;
                continue;
            }
        }
        let mut eligible = eligible_for_deg(q_deg[q]) & !anchor_clear;
        if let Some((p0, orbit_reps)) = root_orbit_filter {
            if q == p0 {
                eligible = eligible & orbit_reps;
            }
        }
        if eligible.is_zero() {
            return false;
        }
        let mask = W::node_candidates(pred, NodeId(q as u32), eligible);
        if mask.is_zero() {
            return false;
        }
        cand[q] = mask;
        cand_cnt[q] = mask.count_ones() as usize;
    }
    true
}

fn solve_word_direct<W: WordMask, const N: usize, const COUNT_ONLY: bool, B, P>(
    target: &Graph,
    query: &Graph,
    anchor: Option<(NodeId, NodeId)>,
    pred: &mut P,
    count: &mut u64,
    visitor: &mut impl FnMut(&[NodeId]) -> ControlFlow<B>,
) -> ControlFlow<B>
where
    P: SubgraphIsomorphismPredicate,
{
    let n1 = query.node_count();
    let n2 = target.node_count();

    let mut deg_ge = [W::ZERO; 16];
    let mut t_deg = [0u16; N];
    for (t, deg_slot) in t_deg.iter_mut().enumerate().take(n2) {
        let d = target.degree(NodeId(t as u32));
        *deg_slot = d as u16;
        let bit = W::bit(t);
        for slot in &mut deg_ge[..=d.min(15)] {
            *slot = *slot | bit;
        }
    }

    let mut q_by_deg = [0usize; N];
    let mut q_deg = [0usize; N];
    for q in 0..n1 {
        q_by_deg[q] = q;
        q_deg[q] = query.degree(NodeId(q as u32));
    }
    q_by_deg[..n1].sort_unstable_by_key(|&q| Reverse(q_deg[q]));

    let mut cand = [W::ZERO; N];
    let mut cand_cnt = [0usize; N];
    let eligible_for_deg = |qd: usize| -> W {
        let base = deg_ge[qd.min(15)];
        if qd < 15 {
            base
        } else {
            let mut exact = W::ZERO;
            let mut rem = base;
            while !rem.is_zero() {
                let t = rem.trailing_zeros() as usize;
                rem = rem.clear_lowest();
                if (t_deg[t] as usize) >= qd {
                    exact = exact | W::bit(t);
                }
            }
            exact
        }
    };
    if !init_word_candidates::<W, N, P>(
        n1,
        &q_deg,
        &q_by_deg,
        anchor,
        None,
        eligible_for_deg,
        pred,
        &mut cand,
        &mut cand_cnt,
    ) {
        return ControlFlow::Continue(());
    }

    let mut order = [0usize; N];
    let mut pos = [0usize; N];
    build_word_ri_order(
        query,
        &q_deg[..n1],
        &cand_cnt[..n1],
        anchor,
        &mut order[..n1],
        &mut pos[..n1],
    );

    let mut parents_offsets = [0usize; 129];
    let mut stack_parents = [(0usize, EdgeId(0)); STACK_PARENTS_LIMIT];
    let mut heap_parents = Vec::new();
    let use_stack_parents = query.edge_count() <= STACK_PARENTS_LIMIT;
    let mut parent_len = 0usize;
    parents_offsets[0] = 0;
    for d in 0..n1 {
        let q = order[d];
        for nb in query.neighbors(NodeId(q as u32)) {
            let p = pos[nb.node.index()];
            if p < d {
                if use_stack_parents {
                    stack_parents[parent_len] = (p, nb.edge);
                    parent_len += 1;
                } else {
                    heap_parents.push((p, nb.edge));
                }
            }
        }
        parents_offsets[d + 1] = if use_stack_parents {
            parent_len
        } else {
            heap_parents.len()
        };
    }
    let parents_data: &[(usize, EdgeId)] = if use_stack_parents {
        &stack_parents[..parent_len]
    } else {
        &heap_parents
    };

    let mut cand_ord = [W::ZERO; N];
    for d in 0..n1 {
        cand_ord[d] = cand[order[d]];
    }

    let mut q_term = [0u32; 128];
    let mut q_rem = [0u32; 128];
    let use_lookahead =
        compute_lookahead_counts(query, &order[..n1], &pos[..n1], &mut q_term, &mut q_rem);

    let mut adj_buf = [W::ZERO; N];
    let mut adj = LazyWordAdj {
        target,
        adj: &mut adj_buf[..n2],
        built: W::ZERO,
    };
    let mut mapping = [0usize; N];
    let mut out_nodes = [NodeId(0); N];

    let mut ctx = WordDfsCtx {
        target,
        n2,
        order: &order[..n1],
        cand_ord: &cand_ord[..n1],
        parents_offsets: &parents_offsets[..=n1],
        parents_data,
        q_term: &q_term[..n1],
        q_rem: &q_rem[..n1],
        all_edges_undetermined: pred.all_edges_undetermined(),
        adj: &mut adj,
        mapping: &mut mapping[..n1],
        out_nodes: &mut out_nodes[..n1],
        sym: None,
        pred,
        count,
        visitor,
    };
    if use_lookahead {
        word_dfs::<W, _, true, COUNT_ONLY, B, P>(&mut ctx, 0, W::ZERO, W::ZERO)
    } else {
        word_dfs::<W, _, false, COUNT_ONLY, B, P>(&mut ctx, 0, W::ZERO, W::ZERO)
    }
}

#[allow(clippy::too_many_arguments)]
fn solve_word_prepared<W: WordMask, const N: usize, const COUNT_ONLY: bool, B, P>(
    target: &Graph,
    prepared_query: &PreparedSubgraphQuery,
    deg_ge: &[W],
    adj_rows: &[W],
    target_sym: Option<&TargetSymmetry>,
    pred: &mut P,
    scratch: &mut SubgraphIsomorphismScratch,
    count: &mut u64,
    visitor: &mut impl FnMut(&[NodeId]) -> ControlFlow<B>,
) -> ControlFlow<B>
where
    P: SubgraphIsomorphismPredicate,
{
    let query = prepared_query.graph();
    let n1 = query.node_count();
    let n2 = target.node_count();

    let mut q_by_deg = [0usize; N];
    let mut q_deg = [0usize; N];
    for (q, &deg) in prepared_query.degrees().iter().enumerate() {
        q_by_deg[q] = q;
        q_deg[q] = deg as usize;
    }
    q_by_deg[..n1].sort_unstable_by_key(|&q| Reverse(q_deg[q]));

    let root_orbit_filter = target_sym.filter(|sym| !sym.is_trivial()).map(|sym| {
        (
            prepared_query.priority_order()[0] as usize,
            W::orbit_reps(sym),
        )
    });

    let mut cand = [W::ZERO; N];
    let mut cand_cnt = [0usize; N];
    if !init_word_candidates::<W, N, P>(
        n1,
        &q_deg,
        &q_by_deg,
        None,
        root_orbit_filter,
        |qd| deg_ge[qd],
        pred,
        &mut cand,
        &mut cand_cnt,
    ) {
        return ControlFlow::Continue(());
    }

    let mut order = [0usize; N];
    let mut pos = [0usize; N];
    build_word_ri_order(
        query,
        &q_deg[..n1],
        &cand_cnt[..n1],
        None,
        &mut order[..n1],
        &mut pos[..n1],
    );

    scratch.parents_offsets.clear();
    scratch.parents_data.clear();
    scratch.parents_offsets.push(0);
    for (d, &q) in order.iter().enumerate().take(n1) {
        for nb in query.neighbors(NodeId(q as u32)) {
            let p = pos[nb.node.index()];
            if p < d {
                scratch.parents_data.push((p, nb.edge));
            }
        }
        scratch.parents_offsets.push(scratch.parents_data.len());
    }

    let mut cand_ord = [W::ZERO; N];
    for d in 0..n1 {
        cand_ord[d] = cand[order[d]];
    }

    let mut q_term = [0u32; 128];
    let mut q_rem = [0u32; 128];
    let use_lookahead =
        compute_lookahead_counts(query, &order[..n1], &pos[..n1], &mut q_term, &mut q_rem);

    let parents_offsets = mem::take(&mut scratch.parents_offsets);
    let parents_data = mem::take(&mut scratch.parents_data);
    scratch.out_nodes.resize(n1, NodeId(0));

    let sym = SearchSymmetry::build(
        Some(prepared_query),
        target_sym,
        &order[..n1],
        &pos[..n1],
        &mut scratch.target_chain,
    );

    let mut adj = PrecomputedWordAdj(adj_rows);
    let mut mapping = [0usize; N];
    let res = {
        let mut ctx = WordDfsCtx {
            target,
            n2,
            order: &order[..n1],
            cand_ord: &cand_ord[..n1],
            parents_offsets: &parents_offsets,
            parents_data: &parents_data,
            q_term: &q_term[..n1],
            q_rem: &q_rem[..n1],
            all_edges_undetermined: pred.all_edges_undetermined(),
            adj: &mut adj,
            mapping: &mut mapping[..n1],
            out_nodes: &mut scratch.out_nodes,
            sym,
            pred,
            count,
            visitor,
        };
        if use_lookahead {
            word_dfs::<W, _, true, COUNT_ONLY, B, P>(&mut ctx, 0, W::ZERO, W::ZERO)
        } else {
            word_dfs::<W, _, false, COUNT_ONLY, B, P>(&mut ctx, 0, W::ZERO, W::ZERO)
        }
    };

    scratch.parents_offsets = parents_offsets;
    scratch.parents_data = parents_data;
    res
}

struct WordDfsCtx<'a, W, A, P, F> {
    target: &'a Graph,
    n2: usize,
    order: &'a [usize],
    cand_ord: &'a [W],
    parents_offsets: &'a [usize],
    parents_data: &'a [(usize, EdgeId)],
    q_term: &'a [u32],
    q_rem: &'a [u32],
    all_edges_undetermined: bool,
    adj: &'a mut A,
    mapping: &'a mut [usize],
    out_nodes: &'a mut [NodeId],
    sym: Option<SearchSymmetry<'a>>,
    pred: &'a mut P,
    count: &'a mut u64,
    visitor: &'a mut F,
}

#[inline(always)]
fn try_descend_word<W, A, const LOOKAHEAD: bool, const COUNT_ONLY: bool, B, P, F>(
    ctx: &mut WordDfsCtx<'_, W, A, P, F>,
    depth: usize,
    v: usize,
    used: W,
    term_t: W,
) -> ControlFlow<B>
where
    W: WordMask,
    A: WordAdj<W>,
    P: SubgraphIsomorphismPredicate,
    F: FnMut(&[NodeId]) -> ControlFlow<B>,
{
    let is_last = depth + 1 == ctx.order.len();
    let next_used = used | W::bit(v);
    let next_term = if LOOKAHEAD && !is_last {
        let av = ctx.adj.row(v);
        let qt = ctx.q_term[depth];
        if qt > 0 && (av & term_t).count_ones() < qt {
            return ControlFlow::Continue(());
        }
        let qr = ctx.q_rem[depth];
        if qr > 0 && (av & !next_used).count_ones() < qr {
            return ControlFlow::Continue(());
        }
        (term_t | av) & !next_used
    } else {
        W::ZERO
    };

    ctx.mapping[depth] = v;
    let saved_depth = if let Some(sym) = ctx.sym.as_mut() {
        let Some(saved) = sym.try_advance_host_chain(depth, ctx.mapping) else {
            return ControlFlow::Continue(());
        };
        saved
    } else {
        usize::MAX
    };

    let res = if is_last {
        let accepted = ctx
            .sym
            .as_mut()
            .is_none_or(|sym| sym.accepts_complete_embedding(ctx.order, ctx.mapping));
        if accepted {
            if COUNT_ONLY {
                *ctx.count += 1;
                ControlFlow::Continue(())
            } else {
                for (d, &q) in ctx.order.iter().enumerate() {
                    ctx.out_nodes[q] = NodeId(ctx.mapping[d] as u32);
                }
                (ctx.visitor)(ctx.out_nodes)
            }
        } else {
            ControlFlow::Continue(())
        }
    } else {
        word_dfs::<W, A, LOOKAHEAD, COUNT_ONLY, B, P>(ctx, depth + 1, next_used, next_term)
    };

    if let Some(sym) = ctx.sym.as_mut() {
        sym.restore_host_chain(saved_depth);
    }
    res
}

fn word_dfs<W, A, const LOOKAHEAD: bool, const COUNT_ONLY: bool, B, P>(
    ctx: &mut WordDfsCtx<'_, W, A, P, impl FnMut(&[NodeId]) -> ControlFlow<B>>,
    depth: usize,
    used: W,
    term_t: W,
) -> ControlFlow<B>
where
    W: WordMask,
    A: WordAdj<W>,
    P: SubgraphIsomorphismPredicate,
{
    let mut avail = ctx.cand_ord[depth] & !used;
    if let Some(sym) = ctx.sym.as_ref() {
        if sym.has_q_sym {
            let (min_v, max_v) = sym.depth_bounds(depth, ctx.mapping, ctx.n2);
            avail = avail & W::range_mask(min_v, max_v);
        }
    }
    if avail.is_zero() {
        return ControlFlow::Continue(());
    }

    let is_last = depth + 1 == ctx.order.len();
    let allow_leaf_popcnt =
        COUNT_ONLY && is_last && ctx.sym.as_ref().is_none_or(SearchSymmetry::can_leaf_popcnt);

    let p_start = ctx.parents_offsets[depth];
    let p_end = ctx.parents_offsets[depth + 1];
    let parent_count = p_end - p_start;

    if parent_count == 0 {
        if allow_leaf_popcnt {
            *ctx.count += avail.count_ones() as u64;
            return ControlFlow::Continue(());
        }
        let mut cands = avail;
        while !cands.is_zero() {
            let v = cands.trailing_zeros() as usize;
            cands = cands.clear_lowest();
            if ctx
                .sym
                .as_ref()
                .is_some_and(|sym| sym.rejects_candidate_before_edge_check(depth, v))
            {
                continue;
            }
            let res = try_descend_word::<W, A, LOOKAHEAD, COUNT_ONLY, B, P, _>(
                ctx, depth, v, used, term_t,
            );
            if let ControlFlow::Break(val) = res {
                return ControlFlow::Break(val);
            }
        }
        return ControlFlow::Continue(());
    }

    if parent_count == 1 {
        let (p, pe) = ctx.parents_data[p_start];
        let mp = ctx.mapping[p];
        let edge_undetermined = ctx.all_edges_undetermined || ctx.pred.edge_is_undetermined(pe);
        let exact_mask = if edge_undetermined {
            Some(ctx.adj.row(mp))
        } else {
            W::edge_neighbors(ctx.pred, pe, NodeId(mp as u32))
        };

        if let Some(pmask) = exact_mask {
            let mut cands = avail & pmask;
            if cands.is_zero() {
                return ControlFlow::Continue(());
            }
            if allow_leaf_popcnt {
                *ctx.count += cands.count_ones() as u64;
                return ControlFlow::Continue(());
            }
            while !cands.is_zero() {
                let v = cands.trailing_zeros() as usize;
                cands = cands.clear_lowest();
                if ctx
                    .sym
                    .as_ref()
                    .is_some_and(|sym| sym.rejects_candidate_before_edge_check(depth, v))
                {
                    continue;
                }
                let res = try_descend_word::<W, A, LOOKAHEAD, COUNT_ONLY, B, P, _>(
                    ctx, depth, v, used, term_t,
                );
                if let ControlFlow::Break(val) = res {
                    return ControlFlow::Break(val);
                }
            }
            return ControlFlow::Continue(());
        }

        let nbs = ctx.target.neighbors(NodeId(mp as u32));
        if nbs.len() <= SPARSE_SINGLE_PARENT_DEG_LIMIT {
            for nb in nbs {
                let v = nb.node.index();
                if !avail.test(v) {
                    continue;
                }
                if ctx
                    .sym
                    .as_ref()
                    .is_some_and(|sym| sym.rejects_candidate_before_edge_check(depth, v))
                {
                    continue;
                }
                if !ctx.pred.match_edge(pe, nb.edge) {
                    continue;
                }
                let res = try_descend_word::<W, A, LOOKAHEAD, COUNT_ONLY, B, P, _>(
                    ctx, depth, v, used, term_t,
                );
                if let ControlFlow::Break(val) = res {
                    return ControlFlow::Break(val);
                }
            }
            return ControlFlow::Continue(());
        }

        let mut cands = avail & ctx.adj.row(mp);
        while !cands.is_zero() {
            let v = cands.trailing_zeros() as usize;
            cands = cands.clear_lowest();
            if ctx
                .sym
                .as_ref()
                .is_some_and(|sym| sym.rejects_candidate_before_edge_check(depth, v))
            {
                continue;
            }
            let te = ctx
                .target
                .find_edge(NodeId(mp as u32), NodeId(v as u32))
                .expect("candidate from adjacency bitset is a neighbor");
            if !ctx.pred.match_edge(pe, te) {
                continue;
            }
            let res = try_descend_word::<W, A, LOOKAHEAD, COUNT_ONLY, B, P, _>(
                ctx, depth, v, used, term_t,
            );
            if let ControlFlow::Break(val) = res {
                return ControlFlow::Break(val);
            }
        }
        return ControlFlow::Continue(());
    }

    // Multi-parent step: intersect all parent adjacency / exact-edge bitmasks first.
    let mut cands = avail;
    let mut unverified_mask = 0u128;
    if ctx.all_edges_undetermined {
        for idx in p_start..p_end {
            let (p, _) = ctx.parents_data[idx];
            cands = cands & ctx.adj.row(ctx.mapping[p]);
            if cands.is_zero() {
                return ControlFlow::Continue(());
            }
        }
    } else {
        for offset in 0..parent_count {
            let (p, pe) = ctx.parents_data[p_start + offset];
            let mp = ctx.mapping[p];
            let pmask = if ctx.pred.edge_is_undetermined(pe) {
                ctx.adj.row(mp)
            } else if let Some(emask) = W::edge_neighbors(ctx.pred, pe, NodeId(mp as u32)) {
                emask
            } else {
                unverified_mask |= 1u128 << offset;
                ctx.adj.row(mp)
            };
            cands = cands & pmask;
            if cands.is_zero() {
                return ControlFlow::Continue(());
            }
        }
    }

    if unverified_mask == 0 && allow_leaf_popcnt {
        *ctx.count += cands.count_ones() as u64;
        return ControlFlow::Continue(());
    }

    while !cands.is_zero() {
        let v = cands.trailing_zeros() as usize;
        cands = cands.clear_lowest();
        if ctx
            .sym
            .as_ref()
            .is_some_and(|sym| sym.rejects_candidate_before_edge_check(depth, v))
        {
            continue;
        }
        if unverified_mask != 0 {
            let v_node = NodeId(v as u32);
            let mut rem = unverified_mask;
            let mut edges_ok = true;
            while rem != 0 {
                let offset = rem.trailing_zeros() as usize;
                rem &= rem - 1;
                let (p, pe) = ctx.parents_data[p_start + offset];
                let te = ctx
                    .target
                    .find_edge(NodeId(ctx.mapping[p] as u32), v_node)
                    .expect("candidate from adjacency bitset is a neighbor");
                if !ctx.pred.match_edge(pe, te) {
                    edges_ok = false;
                    break;
                }
            }
            if !edges_ok {
                continue;
            }
        }
        let res =
            try_descend_word::<W, A, LOOKAHEAD, COUNT_ONLY, B, P, _>(ctx, depth, v, used, term_t);
        if let ControlFlow::Break(val) = res {
            return ControlFlow::Break(val);
        }
    }
    ControlFlow::Continue(())
}

// Wide (`n2 > 128`) BitsetRi path.

#[allow(clippy::too_many_arguments)]
fn solve_wide<const COUNT_ONLY: bool, B, P: SubgraphIsomorphismPredicate>(
    target: &Graph,
    query: &Graph,
    prepared_query: Option<&PreparedSubgraphQuery>,
    prepared_deg_ge: Option<(&[u64], usize)>,
    prepared_adj: Option<&[u64]>,
    anchor: Option<(NodeId, NodeId)>,
    target_sym: Option<&TargetSymmetry>,
    pred: &mut P,
    scratch: &mut SubgraphIsomorphismScratch,
    count: &mut u64,
    visitor: &mut impl FnMut(&[NodeId]) -> ControlFlow<B>,
) -> ControlFlow<B> {
    let n1 = query.node_count();
    let n2 = target.node_count();
    let words = prepared_deg_ge.map_or_else(|| n2.div_ceil(64), |(_, w)| w);

    let q_deg: Vec<usize> = match prepared_query {
        Some(pq) => pq.degrees().iter().map(|&d| d as usize).collect(),
        None => query.node_ids().map(|q| query.degree(q)).collect(),
    };
    let mut q_by_deg: Vec<usize> = (0..n1).collect();
    q_by_deg.sort_unstable_by_key(|&q| Reverse(q_deg[q]));

    let root_orbit_filter = if anchor.is_none() {
        target_sym.filter(|sym| !sym.is_trivial()).map(|sym| {
            let p0 = prepared_query
                .expect("target symmetry is paired with prepared query")
                .priority_order()[0] as usize;
            (p0, sym.orbit_reps_words())
        })
    } else {
        None
    };

    scratch.cand_wide.clear();
    scratch.cand_wide.resize(n1 * words, 0);
    scratch.cand_counts.clear();
    scratch.cand_counts.resize(n1, 0);

    let t_deg: Vec<usize> = if prepared_deg_ge.is_none() {
        target.node_ids().map(|t| target.degree(t)).collect()
    } else {
        Vec::new()
    };
    let mut eligible_buf = if prepared_deg_ge.is_none() || anchor.is_some() {
        vec![0u64; words]
    } else {
        Vec::new()
    };

    let fill_eligible =
        |qd: usize, clear_anchor: Option<usize>, orbit_mask: Option<&[u64]>, buf: &mut [u64]| {
            if let Some((deg_ge, _)) = prepared_deg_ge {
                buf.copy_from_slice(&deg_ge[qd * words..(qd + 1) * words]);
            } else {
                buf.fill(0);
                for (t, &td) in t_deg.iter().enumerate() {
                    if td >= qd {
                        buf[t >> 6] |= 1u64 << (t & 63);
                    }
                }
            }
            if let Some(ta) = clear_anchor {
                buf[ta >> 6] &= !(1u64 << (ta & 63));
            }
            if let Some(orbits) = orbit_mask {
                for (dst, &m) in buf.iter_mut().zip(orbits.iter()) {
                    *dst &= m;
                }
            }
        };

    let anchor_ta = anchor.map(|(_, ta)| ta.index());
    if let Some(num_classes) = pred.query_node_class_count() {
        let mut class_min_deg = vec![usize::MAX; num_classes];
        let mut class_count = vec![0usize; num_classes];
        for &q in &q_by_deg {
            if anchor.is_some_and(|(qa, _)| q == qa.index()) {
                continue;
            }
            let c = pred.query_node_class(NodeId(q as u32));
            debug_assert!(c < num_classes);
            class_min_deg[c] = class_min_deg[c].min(q_deg[q]);
            class_count[c] += 1;
        }
        let mut class_computed = vec![false; num_classes];
        scratch.class_wide.clear();
        scratch.class_wide.resize(num_classes * words, 0);

        for &q in &q_by_deg {
            let base = q * words;
            if let Some((qa, ta)) = anchor {
                if q == qa.index() {
                    let t = ta.index();
                    scratch.cand_wide[base + (t >> 6)] = 1u64 << (t & 63);
                    scratch.cand_counts[q] = 1;
                    continue;
                }
            }
            let c = pred.query_node_class(NodeId(q as u32));
            let cbase = c * words;
            if !class_computed[c] {
                let solo_orbit = root_orbit_filter
                    .filter(|&(p0, _)| class_count[c] == 1 && q == p0)
                    .map(|(_, mask)| mask);
                if let (Some((deg_ge, _)), None, None) = (prepared_deg_ge, anchor_ta, solo_orbit) {
                    let elig = &deg_ge[class_min_deg[c] * words..(class_min_deg[c] + 1) * words];
                    pred.fill_node_candidates_words(
                        NodeId(q as u32),
                        elig,
                        &mut scratch.class_wide[cbase..cbase + words],
                    );
                } else {
                    if eligible_buf.len() != words {
                        eligible_buf.resize(words, 0);
                    }
                    fill_eligible(class_min_deg[c], anchor_ta, solo_orbit, &mut eligible_buf);
                    pred.fill_node_candidates_words(
                        NodeId(q as u32),
                        &eligible_buf,
                        &mut scratch.class_wide[cbase..cbase + words],
                    );
                }
                class_computed[c] = true;
            }

            let qd = q_deg[q];
            let orbit_slice = root_orbit_filter
                .filter(|&(p0, _)| q == p0)
                .map(|(_, mask)| mask);
            let mut cnt = 0usize;
            if let Some((deg_ge, _)) = prepared_deg_ge {
                let deg_slice = &deg_ge[qd * words..(qd + 1) * words];
                for w in 0..words {
                    let mut word = scratch.class_wide[cbase + w] & deg_slice[w];
                    if let Some(orbits) = orbit_slice {
                        word &= orbits[w];
                    }
                    scratch.cand_wide[base + w] = word;
                    cnt += word.count_ones() as usize;
                }
            } else {
                fill_eligible(qd, anchor_ta, orbit_slice, &mut eligible_buf);
                for (w, &elig_word) in eligible_buf.iter().enumerate().take(words) {
                    let word = scratch.class_wide[cbase + w] & elig_word;
                    scratch.cand_wide[base + w] = word;
                    cnt += word.count_ones() as usize;
                }
            }
            if cnt == 0 {
                return ControlFlow::Continue(());
            }
            scratch.cand_counts[q] = cnt;
        }
    } else {
        for &q in &q_by_deg {
            let base = q * words;
            if let Some((qa, ta)) = anchor {
                if q == qa.index() {
                    let t = ta.index();
                    scratch.cand_wide[base + (t >> 6)] = 1u64 << (t & 63);
                    scratch.cand_counts[q] = 1;
                    continue;
                }
            }
            let qd = q_deg[q];
            let orbit_slice = root_orbit_filter
                .filter(|&(p0, _)| q == p0)
                .map(|(_, mask)| mask);
            if let (Some((deg_ge, _)), None, None) = (prepared_deg_ge, anchor_ta, orbit_slice) {
                let elig = &deg_ge[qd * words..(qd + 1) * words];
                pred.fill_node_candidates_words(
                    NodeId(q as u32),
                    elig,
                    &mut scratch.cand_wide[base..base + words],
                );
            } else {
                if eligible_buf.len() != words {
                    eligible_buf.resize(words, 0);
                }
                fill_eligible(qd, anchor_ta, orbit_slice, &mut eligible_buf);
                pred.fill_node_candidates_words(
                    NodeId(q as u32),
                    &eligible_buf,
                    &mut scratch.cand_wide[base..base + words],
                );
            }
            let mut cnt = 0usize;
            for &word in &scratch.cand_wide[base..base + words] {
                cnt += word.count_ones() as usize;
            }
            if cnt == 0 {
                return ControlFlow::Continue(());
            }
            scratch.cand_counts[q] = cnt;
        }
    }

    build_wide_ri_order(query, &q_deg, anchor, scratch);

    scratch.parents_offsets.clear();
    scratch.parents_data.clear();
    scratch.parents_offsets.push(0);
    for d in 0..n1 {
        let q = scratch.order[d];
        for nb in query.neighbors(NodeId(q as u32)) {
            let p = scratch.pos[nb.node.index()];
            if p < d {
                scratch.parents_data.push((p, nb.edge));
            }
        }
        scratch.parents_offsets.push(scratch.parents_data.len());
    }

    let order = mem::take(&mut scratch.order);
    let pos = mem::take(&mut scratch.pos);
    let parents_offsets = mem::take(&mut scratch.parents_offsets);
    let parents_data = mem::take(&mut scratch.parents_data);
    let cand_wide = mem::take(&mut scratch.cand_wide);

    let has_multi_parent = (0..n1).any(|d| parents_offsets[d + 1] - parents_offsets[d] >= 2);
    let use_lazy_adj = prepared_adj.is_none() && has_multi_parent && n2 <= WIDE_LAZY_ADJ_LIMIT;
    if use_lazy_adj {
        scratch.adj_wide.resize(n2 * words, 0);
        scratch.adj_built_wide.clear();
        scratch.adj_built_wide.resize(words, 0);
    }
    let mut adj_wide = mem::take(&mut scratch.adj_wide);
    let mut adj_built_wide = mem::take(&mut scratch.adj_built_wide);

    scratch.mapping.clear();
    scratch.mapping.resize(n1, 0);
    scratch.used.clear();
    scratch.used.resize(n2, false);
    scratch.out_nodes.resize(n1, NodeId(0));

    let adj = if let Some(pre) = prepared_adj {
        WideAdj::Precomputed(pre)
    } else if use_lazy_adj {
        WideAdj::Lazy {
            adj: &mut adj_wide,
            built: &mut adj_built_wide,
        }
    } else {
        WideAdj::None
    };

    let sym = SearchSymmetry::build(
        prepared_query,
        target_sym,
        &order,
        &pos,
        &mut scratch.target_chain,
    );

    let res = {
        let mut ctx = WideDfsCtx {
            target,
            n2,
            words,
            order: &order,
            cand_wide: &cand_wide,
            parents_offsets: &parents_offsets,
            parents_data: &parents_data,
            all_edges_undetermined: pred.all_edges_undetermined(),
            adj,
            mapping: &mut scratch.mapping,
            used: &mut scratch.used,
            out_nodes: &mut scratch.out_nodes,
            sym,
            pred,
            count,
            visitor,
        };
        wide_dfs::<COUNT_ONLY, B, P>(&mut ctx, 0)
    };

    scratch.order = order;
    scratch.pos = pos;
    scratch.parents_offsets = parents_offsets;
    scratch.parents_data = parents_data;
    scratch.cand_wide = cand_wide;
    scratch.adj_wide = adj_wide;
    scratch.adj_built_wide = adj_built_wide;
    res
}

fn build_wide_ri_order(
    query: &Graph,
    q_deg: &[usize],
    anchor: Option<(NodeId, NodeId)>,
    scratch: &mut SubgraphIsomorphismScratch,
) {
    let n1 = query.node_count();
    scratch.order.clear();
    scratch.pos.clear();
    scratch.pos.resize(n1, 0);
    scratch.placed.clear();
    scratch.placed.resize(n1, false);
    scratch.ordered_neighbors.clear();
    scratch.ordered_neighbors.resize(n1, 0);

    let start = if let Some((qa, _)) = anchor {
        qa.index()
    } else {
        (0..n1)
            .min_by_key(|&q| (scratch.cand_counts[q], Reverse(q_deg[q]), q))
            .expect("non-empty query")
    };
    scratch.order.push(start);
    scratch.pos[start] = 0;
    scratch.placed[start] = true;
    for nb in query.neighbors(NodeId(start as u32)) {
        scratch.ordered_neighbors[nb.node.index()] += 1;
    }

    for step in 1..n1 {
        let mut best = None;
        for (q, &qd) in q_deg.iter().enumerate().take(n1) {
            if scratch.placed[q] {
                continue;
            }
            let c = scratch.ordered_neighbors[q];
            let mut frontier_deg = 0usize;
            for nb in query.neighbors(NodeId(q as u32)) {
                let u = nb.node.index();
                if !scratch.placed[u] && scratch.ordered_neighbors[u] > 0 {
                    frontier_deg += 1;
                }
            }
            let key = (
                c,
                Reverse(scratch.cand_counts[q]),
                frontier_deg,
                qd,
                Reverse(q),
            );
            if best.is_none_or(|(best_key, _)| key > best_key) {
                best = Some((key, q));
            }
        }
        let next = best.expect("unplaced query vertex exists").1;
        scratch.order.push(next);
        scratch.pos[next] = step;
        scratch.placed[next] = true;
        for nb in query.neighbors(NodeId(next as u32)) {
            scratch.ordered_neighbors[nb.node.index()] += 1;
        }
    }
}

struct WideDfsCtx<'a, P, F> {
    target: &'a Graph,
    n2: usize,
    words: usize,
    order: &'a [usize],
    cand_wide: &'a [u64],
    parents_offsets: &'a [usize],
    parents_data: &'a [(usize, EdgeId)],
    all_edges_undetermined: bool,
    adj: WideAdj<'a>,
    mapping: &'a mut [usize],
    used: &'a mut [bool],
    out_nodes: &'a mut [NodeId],
    sym: Option<SearchSymmetry<'a>>,
    pred: &'a mut P,
    count: &'a mut u64,
    visitor: &'a mut F,
}

#[inline(always)]
fn try_descend_wide<const COUNT_ONLY: bool, B, P, F>(
    ctx: &mut WideDfsCtx<'_, P, F>,
    depth: usize,
    v: usize,
) -> ControlFlow<B>
where
    P: SubgraphIsomorphismPredicate,
    F: FnMut(&[NodeId]) -> ControlFlow<B>,
{
    let is_last = depth + 1 == ctx.order.len();
    ctx.mapping[depth] = v;
    let saved_depth = if let Some(sym) = ctx.sym.as_mut() {
        let Some(saved) = sym.try_advance_host_chain(depth, ctx.mapping) else {
            return ControlFlow::Continue(());
        };
        saved
    } else {
        usize::MAX
    };

    let res = if is_last {
        let accepted = ctx
            .sym
            .as_mut()
            .is_none_or(|sym| sym.accepts_complete_embedding(ctx.order, ctx.mapping));
        if accepted {
            if COUNT_ONLY {
                *ctx.count += 1;
                ControlFlow::Continue(())
            } else {
                for (d, &q) in ctx.order.iter().enumerate() {
                    ctx.out_nodes[q] = NodeId(ctx.mapping[d] as u32);
                }
                (ctx.visitor)(ctx.out_nodes)
            }
        } else {
            ControlFlow::Continue(())
        }
    } else {
        ctx.used[v] = true;
        let inner = wide_dfs::<COUNT_ONLY, B, P>(ctx, depth + 1);
        ctx.used[v] = false;
        inner
    };

    if let Some(sym) = ctx.sym.as_mut() {
        sym.restore_host_chain(saved_depth);
    }
    res
}

fn wide_dfs<const COUNT_ONLY: bool, B, P: SubgraphIsomorphismPredicate>(
    ctx: &mut WideDfsCtx<'_, P, impl FnMut(&[NodeId]) -> ControlFlow<B>>,
    depth: usize,
) -> ControlFlow<B> {
    let (min_v, max_v) = match ctx.sym.as_ref() {
        Some(sym) if sym.has_q_sym => sym.depth_bounds(depth, ctx.mapping, ctx.n2),
        _ => (0, ctx.n2),
    };
    if min_v >= max_v {
        return ControlFlow::Continue(());
    }

    let q = ctx.order[depth];
    let cbase = q * ctx.words;
    let p_start = ctx.parents_offsets[depth];
    let p_end = ctx.parents_offsets[depth + 1];
    let parent_count = p_end - p_start;
    let is_last = depth + 1 == ctx.order.len();
    let allow_leaf_shortcut =
        COUNT_ONLY && is_last && ctx.sym.as_ref().is_none_or(SearchSymmetry::can_leaf_popcnt);

    if parent_count == 0 {
        let w_start = min_v >> 6;
        let w_end = (max_v - 1) >> 6;
        for w in w_start..=w_end {
            let mut bits = ctx.cand_wide[cbase + w];
            if w == w_start && (min_v & 63) != 0 {
                bits &= u64::MAX << (min_v & 63);
            }
            if w == w_end && (max_v & 63) != 0 {
                bits &= (1u64 << (max_v & 63)) - 1;
            }
            while bits != 0 {
                let tz = bits.trailing_zeros() as usize;
                bits &= bits - 1;
                let v = (w << 6) + tz;
                if ctx.used[v] {
                    continue;
                }
                if allow_leaf_shortcut {
                    *ctx.count += 1;
                    continue;
                }
                if ctx
                    .sym
                    .as_ref()
                    .is_some_and(|sym| sym.rejects_candidate_before_edge_check(depth, v))
                {
                    continue;
                }
                let res = try_descend_wide::<COUNT_ONLY, B, P, _>(ctx, depth, v);
                if let ControlFlow::Break(val) = res {
                    return ControlFlow::Break(val);
                }
            }
        }
        return ControlFlow::Continue(());
    }

    if parent_count == 1 {
        let (p, pe) = ctx.parents_data[p_start];
        let mp = ctx.mapping[p];
        let edge_undetermined = ctx.all_edges_undetermined || ctx.pred.edge_is_undetermined(pe);
        for nb in ctx.target.neighbors(NodeId(mp as u32)) {
            let v = nb.node.index();
            if v < min_v
                || v >= max_v
                || ctx.used[v]
                || ((ctx.cand_wide[cbase + (v >> 6)] >> (v & 63)) & 1) == 0
            {
                continue;
            }
            if ctx
                .sym
                .as_ref()
                .is_some_and(|sym| sym.rejects_candidate_before_edge_check(depth, v))
            {
                continue;
            }
            if !edge_undetermined && !ctx.pred.match_edge(pe, nb.edge) {
                continue;
            }
            if allow_leaf_shortcut {
                *ctx.count += 1;
                continue;
            }
            let res = try_descend_wide::<COUNT_ONLY, B, P, _>(ctx, depth, v);
            if let ControlFlow::Break(val) = res {
                return ControlFlow::Break(val);
            }
        }
        return ControlFlow::Continue(());
    }

    // Multi-parent step: iterate the minimum-degree parent's neighbors and check other parents.
    let mut min_offset = 0usize;
    let mut min_deg = ctx
        .target
        .degree(NodeId(ctx.mapping[ctx.parents_data[p_start].0] as u32));
    for offset in 1..parent_count {
        let d = ctx.target.degree(NodeId(
            ctx.mapping[ctx.parents_data[p_start + offset].0] as u32,
        ));
        if d < min_deg {
            min_deg = d;
            min_offset = offset;
        }
    }
    for offset in 0..parent_count {
        if offset != min_offset {
            let pk = ctx.parents_data[p_start + offset].0;
            ctx.adj.ensure_row(ctx.target, ctx.mapping[pk], ctx.words);
        }
    }

    let (p0, pe0) = ctx.parents_data[p_start + min_offset];
    let mp0 = ctx.mapping[p0];
    let pe0_undetermined = ctx.all_edges_undetermined || ctx.pred.edge_is_undetermined(pe0);
    for nb in ctx.target.neighbors(NodeId(mp0 as u32)) {
        let v = nb.node.index();
        if v < min_v
            || v >= max_v
            || ctx.used[v]
            || ((ctx.cand_wide[cbase + (v >> 6)] >> (v & 63)) & 1) == 0
        {
            continue;
        }
        let vw = v >> 6;
        let vb = 1u64 << (v & 63);
        let mut adj_ok = true;
        for offset in 0..parent_count {
            if offset == min_offset {
                continue;
            }
            let pk = ctx.parents_data[p_start + offset].0;
            if let Some(has_bit) = ctx.adj.has_edge_bit(ctx.mapping[pk], vw, vb, ctx.words) {
                if !has_bit {
                    adj_ok = false;
                    break;
                }
            }
        }
        if !adj_ok {
            continue;
        }
        if ctx
            .sym
            .as_ref()
            .is_some_and(|sym| sym.rejects_candidate_before_edge_check(depth, v))
        {
            continue;
        }
        if !pe0_undetermined && !ctx.pred.match_edge(pe0, nb.edge) {
            continue;
        }
        let v_node = nb.node;
        let mut edges_ok = true;
        for offset in 0..parent_count {
            if offset == min_offset {
                continue;
            }
            let (pk, pek) = ctx.parents_data[p_start + offset];
            let mpk = ctx.mapping[pk];
            let pek_undetermined = ctx.all_edges_undetermined || ctx.pred.edge_is_undetermined(pek);
            if pek_undetermined && !matches!(ctx.adj, WideAdj::None) {
                continue;
            }
            let Some(te) = ctx.target.find_edge(NodeId(mpk as u32), v_node) else {
                edges_ok = false;
                break;
            };
            if !pek_undetermined && !ctx.pred.match_edge(pek, te) {
                edges_ok = false;
                break;
            }
        }
        if !edges_ok {
            continue;
        }
        if allow_leaf_shortcut {
            *ctx.count += 1;
            continue;
        }
        let res = try_descend_wide::<COUNT_ONLY, B, P, _>(ctx, depth, v);
        if let ControlFlow::Break(val) = res {
            return ControlFlow::Break(val);
        }
    }
    ControlFlow::Continue(())
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::ops::ControlFlow;

    use pretty_assertions::assert_eq;
    use rstest::rstest;

    use super::*;
    use crate::algorithms::automorphism::AutomorphismAlgorithm;
    use crate::algorithms::subgraph_isomorphism::{
        ClosureSubgraphIsomorphismPredicate, SubgraphIsomorphismAlgorithm,
    };
    use crate::union_find::UnionFind;

    const SIMPLE_GRAPHS_GRAPH6: &str = include_str!("../../../tests/data/simple-through-8.g6");

    fn parse_graph6(source: &str) -> Graph {
        let bytes = source.as_bytes();
        let node_count = (bytes[0] - 63) as usize;
        let mut edges = Vec::new();
        let mut bit = 0;
        for second in 1..node_count {
            for first in 0..second {
                let value = bytes[1 + bit / 6] - 63;
                if value & (1u8 << (5 - bit % 6)) != 0 {
                    edges.push([first as u32, second as u32]);
                }
                bit += 1;
            }
        }
        Graph::new(node_count, &edges)
    }

    struct UndeterminedPredicate;

    impl SubgraphIsomorphismPredicate for UndeterminedPredicate {
        fn match_node(&mut self, _query_node: NodeId, _target_node: NodeId) -> bool {
            true
        }

        fn match_edge(&mut self, _query_edge: EdgeId, _target_edge: EdgeId) -> bool {
            true
        }

        fn all_edges_undetermined(&self) -> bool {
            true
        }
    }

    fn collect_vf2rdkit(
        target: &Graph,
        query: &Graph,
        anchor: Option<(NodeId, NodeId)>,
        node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
    ) -> Vec<Vec<usize>> {
        let mut out = Vec::new();
        let mut emit = |emb: &[NodeId]| {
            out.push(emb.iter().map(|n| n.index()).collect());
            ControlFlow::<()>::Continue(())
        };
        let _ = match anchor {
            None => target.visit_subgraph_isomorphisms(
                query,
                node_match,
                edge_match,
                SubgraphIsomorphismAlgorithm::Vf2Rdkit,
                &mut emit,
            ),
            Some(a) => target.visit_subgraph_isomorphisms_at(
                query,
                a,
                node_match,
                edge_match,
                SubgraphIsomorphismAlgorithm::Vf2Rdkit,
                &mut emit,
            ),
        };
        out.sort_unstable();
        out
    }

    fn canonical_orbit_reps(
        embeddings: &[Vec<usize>],
        priority_order: &[u32],
        query_gens: &[Vec<NodeId>],
        target_gens: &[Vec<NodeId>],
    ) -> Vec<Vec<usize>> {
        let index: HashMap<&[usize], usize> = embeddings
            .iter()
            .enumerate()
            .map(|(i, e)| (e.as_slice(), i))
            .collect();
        let mut uf = UnionFind::new(embeddings.len());
        let mut buf = Vec::new();
        for (i, emb) in embeddings.iter().enumerate() {
            for g in query_gens {
                buf.clear();
                buf.extend(g.iter().map(|&q| emb[q.index()]));
                let j = index[buf.as_slice()];
                uf.union(i, j);
            }
            for h in target_gens {
                buf.clear();
                buf.extend(emb.iter().map(|&t| h[t].index()));
                let j = index[buf.as_slice()];
                uf.union(i, j);
            }
        }
        let mut groups: HashMap<usize, Vec<&Vec<usize>>> = HashMap::new();
        for (i, emb) in embeddings.iter().enumerate() {
            groups.entry(uf.find(i)).or_default().push(emb);
        }
        let mut reps: Vec<Vec<usize>> = groups
            .into_values()
            .map(|members| {
                members
                    .into_iter()
                    .min_by(|a, b| {
                        priority_order
                            .iter()
                            .map(|&q| a[q as usize])
                            .cmp(priority_order.iter().map(|&q| b[q as usize]))
                    })
                    .expect("non-empty orbit")
                    .clone()
            })
            .collect();
        reps.sort_unstable();
        reps
    }

    struct HookedParityPredicate {
        target: Graph,
        q_classes: Vec<usize>,
        t_classes: Vec<usize>,
        undetermined_edges: Vec<bool>,
        masked_edges: Vec<bool>,
        all_undetermined: bool,
        node_mask_calls: usize,
        match_edge_calls: usize,
    }

    impl SubgraphIsomorphismPredicate for HookedParityPredicate {
        fn match_node(&mut self, query_node: NodeId, target_node: NodeId) -> bool {
            self.q_classes[query_node.index()] == self.t_classes[target_node.index()]
        }

        fn match_edge(&mut self, query_edge: EdgeId, target_edge: EdgeId) -> bool {
            self.match_edge_calls += 1;
            if self.all_undetermined || self.undetermined_edges[query_edge.index()] {
                true
            } else {
                query_edge.index() % 2 == target_edge.index() % 2
            }
        }

        fn all_edges_undetermined(&self) -> bool {
            self.all_undetermined
        }

        fn edge_is_undetermined(&self, query_edge: EdgeId) -> bool {
            self.all_undetermined || self.undetermined_edges[query_edge.index()]
        }

        fn query_node_class_count(&self) -> Option<usize> {
            Some(self.q_classes.iter().copied().max().map_or(0, |m| m + 1))
        }

        fn query_node_class(&self, query_node: NodeId) -> usize {
            self.q_classes[query_node.index()]
        }

        fn node_candidates_u64(&mut self, query_node: NodeId, eligible: u64) -> u64 {
            self.node_mask_calls += 1;
            let qc = self.q_classes[query_node.index()];
            let mut mask = 0u64;
            let mut rem = eligible;
            while rem != 0 {
                let t = rem.trailing_zeros() as usize;
                rem &= rem - 1;
                if self.t_classes[t] == qc {
                    mask |= 1u64 << t;
                }
            }
            mask
        }

        fn node_candidates_u128(&mut self, query_node: NodeId, eligible: u128) -> u128 {
            self.node_mask_calls += 1;
            let qc = self.q_classes[query_node.index()];
            let mut mask = 0u128;
            let mut rem = eligible;
            while rem != 0 {
                let t = rem.trailing_zeros() as usize;
                rem &= rem - 1;
                if self.t_classes[t] == qc {
                    mask |= 1u128 << t;
                }
            }
            mask
        }

        fn fill_node_candidates_words(
            &mut self,
            query_node: NodeId,
            eligible: &[u64],
            out: &mut [u64],
        ) {
            self.node_mask_calls += 1;
            let qc = self.q_classes[query_node.index()];
            let len = eligible.len().min(out.len());
            for w in 0..len {
                let mut mask = 0u64;
                let mut rem = eligible[w];
                while rem != 0 {
                    let b = rem.trailing_zeros() as usize;
                    rem &= rem - 1;
                    let t = (w << 6) + b;
                    if self.t_classes[t] == qc {
                        mask |= 1u64 << b;
                    }
                }
                out[w] = mask;
            }
            out[len..].fill(0);
        }

        fn edge_neighbors_u64(&mut self, query_edge: EdgeId, target_node: NodeId) -> Option<u64> {
            if !self.masked_edges[query_edge.index()] {
                return None;
            }
            let parity = query_edge.index() % 2;
            let mut mask = 0u64;
            for nb in self.target.neighbors(target_node) {
                if nb.edge.index() % 2 == parity {
                    mask |= 1u64 << nb.node.index();
                }
            }
            Some(mask)
        }

        fn edge_neighbors_u128(&mut self, query_edge: EdgeId, target_node: NodeId) -> Option<u128> {
            if !self.masked_edges[query_edge.index()] {
                return None;
            }
            let parity = query_edge.index() % 2;
            let mut mask = 0u128;
            for nb in self.target.neighbors(target_node) {
                if nb.edge.index() % 2 == parity {
                    mask |= 1u128 << nb.node.index();
                }
            }
            Some(mask)
        }
    }

    #[rstest]
    #[case::empty_query(Graph::new(4, &[[0, 1], [1, 2]]), Graph::new(0, &[]), false)]
    #[case::singleton_in_edgeless(Graph::new(4, &[]), Graph::new(1, &[]), false)]
    #[case::cubane_square(
        Graph::new(
            8,
            &[
                [0, 1], [1, 2], [2, 3], [3, 0],
                [4, 5], [5, 6], [6, 7], [7, 4],
                [0, 4], [1, 5], [2, 6], [3, 7],
            ],
        ),
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
        false
    )]
    #[case::naphthalene_hexagon(
        Graph::new(
            10,
            &[
                [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0],
                [3, 6], [6, 7], [7, 8], [8, 9], [9, 4],
            ],
        ),
        Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]),
        false
    )]
    #[case::naphthalene_hexagon_labeled(
        Graph::new(
            10,
            &[
                [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0],
                [3, 6], [6, 7], [7, 8], [8, 9], [9, 4],
            ],
        ),
        Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]),
        true
    )]
    #[case::rook_k5_unsat(
        Graph::new(
            16,
            &(0..4u32)
                .flat_map(|r| (0..4u32).flat_map(move |c1| ((c1 + 1)..4u32).map(move |c2| [r * 4 + c1, r * 4 + c2])))
                .chain((0..4u32).flat_map(|c| (0..4u32).flat_map(move |r1| ((r1 + 1)..4u32).map(move |r2| [r1 * 4 + c, r2 * 4 + c]))))
                .collect::<Vec<_>>(),
        ),
        Graph::new(
            5,
            &[[0, 1], [0, 2], [0, 3], [0, 4], [1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]],
        ),
        false
    )]
    #[case::disconnected_two_edges(
        Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]),
        Graph::new(4, &[[0, 1], [2, 3]]),
        false
    )]
    #[case::high_degree_wheel_16_with_deg15_decoy(
        Graph::new(
            33,
            &(0..16u32)
                .flat_map(|i| [[0, i + 1], [i + 1, (i + 1) % 16 + 1]])
                .chain((0..15u32).flat_map(|i| [[17, 18 + i], [18 + i, 18 + (i + 1) % 15]]))
                .collect::<Vec<_>>(),
        ),
        Graph::new(
            17,
            &(0..16u32)
                .flat_map(|i| [[0, i + 1], [i + 1, (i + 1) % 16 + 1]])
                .collect::<Vec<_>>(),
        ),
        false
    )]
    #[case::u128_cycle_path(
        Graph::new(70, &(0..70u32).map(|i| [i, (i + 1) % 70]).collect::<Vec<_>>()),
        Graph::new(4, &[[0, 1], [1, 2], [2, 3]]),
        false
    )]
    #[case::u128_cycle_labeled(
        Graph::new(70, &(0..70u32).map(|i| [i, (i + 1) % 70]).collect::<Vec<_>>()),
        Graph::new(4, &[[0, 1], [1, 2], [2, 3]]),
        true
    )]
    #[case::wide_eager_adj_144(
        Graph::new(144, &[[0, 1], [1, 2], [2, 0], [100, 130], [130, 140], [140, 100]]),
        Graph::new(3, &[[0, 1], [1, 2], [2, 0]]),
        false
    )]
    #[case::wide_lazy_adj_1050(
        Graph::new(1050, &[[0, 500], [500, 1049], [1049, 0], [10, 20], [20, 30]]),
        Graph::new(3, &[[0, 1], [1, 2], [2, 0]]),
        false
    )]
    fn test_visit_direct(#[case] target: Graph, #[case] query: Graph, #[case] use_labels: bool) {
        let mut nm_ref = |q: NodeId, t: NodeId| !use_labels || (q.index() % 2 == t.index() % 2);
        let mut em_ref = |qe: EdgeId, te: EdgeId| !use_labels || (qe.index() % 2 == te.index() % 2);
        let expected = collect_vf2rdkit(&target, &query, None, &mut nm_ref, &mut em_ref);

        let mut direct_out = Vec::new();
        let mut pred_direct = ClosureSubgraphIsomorphismPredicate::new(
            |q: NodeId, t: NodeId| !use_labels || (q.index() % 2 == t.index() % 2),
            |qe: EdgeId, te: EdgeId| !use_labels || (qe.index() % 2 == te.index() % 2),
        );
        let flow: ControlFlow<()> =
            visit_direct(&target, &query, None, &mut pred_direct, &mut |emb| {
                direct_out.push(emb.iter().map(|n| n.index()).collect::<Vec<_>>());
                ControlFlow::Continue(())
            });
        assert_eq!(flow, ControlFlow::Continue(()));
        direct_out.sort_unstable();
        assert_eq!(direct_out, expected);
        assert_eq!(
            count_direct(&target, &query, None, &mut pred_direct),
            expected.len() as u64
        );

        let pt = PreparedSubgraphTarget::new(&target);
        let pq = PreparedSubgraphQuery::new(&query, SubgraphIsomorphismAlgorithm::BitsetRi);
        let mut scratch = SubgraphIsomorphismScratch::new();
        let mut prep_out = Vec::new();
        let flow: ControlFlow<()> =
            visit_prepared(&pt, &pq, None, &mut pred_direct, &mut scratch, &mut |emb| {
                prep_out.push(emb.iter().map(|n| n.index()).collect::<Vec<_>>());
                ControlFlow::Continue(())
            });
        assert_eq!(flow, ControlFlow::Continue(()));
        prep_out.sort_unstable();
        assert_eq!(prep_out, expected);
        assert_eq!(
            count_prepared(&pt, &pq, None, &mut pred_direct, &mut scratch),
            expected.len() as u64
        );
    }

    #[rstest]
    fn test_visit_direct_wide_sparse() {
        let target = Graph::new(4200, &[[1, 2000], [2000, 4199], [4199, 1], [10, 11]]);
        let query = Graph::new(3, &[[0, 1], [1, 2], [2, 0]]);
        let mut pred = ClosureSubgraphIsomorphismPredicate::new(|_, _| true, |_, _| true);
        let mut out = Vec::new();
        let _: ControlFlow<()> = visit_direct(&target, &query, None, &mut pred, &mut |emb| {
            out.push(emb.iter().map(|n| n.index()).collect::<Vec<_>>());
            ControlFlow::Continue(())
        });
        out.sort_unstable();
        assert_eq!(
            out,
            vec![
                vec![1, 2000, 4199],
                vec![1, 4199, 2000],
                vec![2000, 1, 4199],
                vec![2000, 4199, 1],
                vec![4199, 1, 2000],
                vec![4199, 2000, 1],
            ]
        );
        assert_eq!(count_direct(&target, &query, None, &mut pred), 6);
    }

    #[rstest]
    #[case::u64_cubane(
        Graph::new(
            8,
            &[
                [0, 1], [1, 2], [2, 3], [3, 0],
                [4, 5], [5, 6], [6, 7], [7, 4],
                [0, 4], [1, 5], [2, 6], [3, 7],
            ],
        ),
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
        (NodeId(0), NodeId(3))
    )]
    #[case::u64_singleton(Graph::new(4, &[[0, 1]]), Graph::new(1, &[]), (NodeId(0), NodeId(2)))]
    #[case::u64_empty_query(Graph::new(3, &[[0, 1]]), Graph::new(0, &[]), (NodeId(0), NodeId(0)))]
    #[case::u64_oob_anchor(Graph::new(3, &[[0, 1]]), Graph::new(2, &[[0, 1]]), (NodeId(2), NodeId(0)))]
    #[case::u64_degree_reject(Graph::new(4, &[[0, 1], [1, 2]]), Graph::new(3, &[[0, 1], [1, 2]]), (NodeId(1), NodeId(0)))]
    #[case::u128_cycle(
        Graph::new(70, &(0..70u32).map(|i| [i, (i + 1) % 70]).collect::<Vec<_>>()),
        Graph::new(3, &[[0, 1], [1, 2]]),
        (NodeId(1), NodeId(65))
    )]
    #[case::wide_cycle(
        Graph::new(144, &(0..144u32).map(|i| [i, (i + 1) % 144]).collect::<Vec<_>>()),
        Graph::new(3, &[[0, 1], [1, 2]]),
        (NodeId(1), NodeId(130))
    )]
    fn test_visit_direct_anchored(
        #[case] target: Graph,
        #[case] query: Graph,
        #[case] anchor: (NodeId, NodeId),
    ) {
        let mut nm = |_: NodeId, _: NodeId| true;
        let mut em = |_: EdgeId, _: EdgeId| true;
        let expected = collect_vf2rdkit(&target, &query, Some(anchor), &mut nm, &mut em);

        let mut pred = ClosureSubgraphIsomorphismPredicate::new(|_, _| true, |_, _| true);
        let mut out = Vec::new();
        let _: ControlFlow<()> =
            visit_direct(&target, &query, Some(anchor), &mut pred, &mut |emb| {
                out.push(emb.iter().map(|n| n.index()).collect::<Vec<_>>());
                ControlFlow::Continue(())
            });
        out.sort_unstable();
        assert_eq!(out, expected);
        assert_eq!(
            count_direct(&target, &query, Some(anchor), &mut pred),
            expected.len() as u64
        );
    }

    #[rstest]
    #[case::u64_all_undetermined(
        Graph::new(
            8,
            &[
                [0, 1], [1, 2], [2, 3], [3, 0],
                [4, 5], [5, 6], [6, 7], [7, 4],
                [0, 4], [1, 5], [2, 6], [3, 7],
            ],
        ),
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
        true,
        false,
        true
    )]
    #[case::u64_masked_edges(
        Graph::new(
            10,
            &[
                [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0],
                [3, 6], [6, 7], [7, 8], [8, 9], [9, 4],
            ],
        ),
        Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]),
        false,
        true,
        true
    )]
    #[case::u64_mixed_edges(
        Graph::new(
            10,
            &[
                [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0],
                [3, 6], [6, 7], [7, 8], [8, 9], [9, 4],
            ],
        ),
        Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]),
        false,
        false,
        false
    )]
    #[case::u128_masked_edges(
        Graph::new(70, &(0..70u32).map(|i| [i, (i + 1) % 70]).collect::<Vec<_>>()),
        Graph::new(4, &[[0, 1], [1, 2], [2, 3]]),
        false,
        true,
        true
    )]
    #[case::wide_all_undetermined(
        Graph::new(
            144,
            &[[0, 1], [1, 2], [2, 3], [3, 0], [130, 131], [131, 132], [132, 133], [133, 130]],
        ),
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
        true,
        false,
        true
    )]
    fn test_visit_prepared_predicate_hooks(
        #[case] target: Graph,
        #[case] query: Graph,
        #[case] all_undetermined: bool,
        #[case] all_masked: bool,
        #[case] expect_zero_match_edge: bool,
    ) {
        let n1 = query.node_count();
        let n2 = target.node_count();
        let m1 = query.edge_count();
        let q_classes: Vec<usize> = (0..n1).map(|q| q % 2).collect();
        let t_classes: Vec<usize> = (0..n2).map(|t| t % 2).collect();
        let undetermined_edges: Vec<bool> = (0..m1).map(|e| all_undetermined || e == 0).collect();
        let masked_edges: Vec<bool> = (0..m1)
            .map(|e| !undetermined_edges[e] && (all_masked || e % 2 == 1))
            .collect();

        let mut nm_ref = |q: NodeId, t: NodeId| q_classes[q.index()] == t_classes[t.index()];
        let mut em_ref = |qe: EdgeId, te: EdgeId| {
            all_undetermined || undetermined_edges[qe.index()] || qe.index() % 2 == te.index() % 2
        };
        let expected = collect_vf2rdkit(&target, &query, None, &mut nm_ref, &mut em_ref);

        let pt = PreparedSubgraphTarget::new(&target);
        let pq = PreparedSubgraphQuery::new(&query, SubgraphIsomorphismAlgorithm::BitsetRi);
        let mut scratch = SubgraphIsomorphismScratch::new();

        let mut pred = HookedParityPredicate {
            target: target.clone(),
            q_classes,
            t_classes,
            undetermined_edges,
            masked_edges,
            all_undetermined,
            node_mask_calls: 0,
            match_edge_calls: 0,
        };
        let mut out = Vec::new();
        let _: ControlFlow<()> =
            visit_prepared(&pt, &pq, None, &mut pred, &mut scratch, &mut |emb| {
                out.push(emb.iter().map(|n| n.index()).collect::<Vec<_>>());
                ControlFlow::Continue(())
            });
        out.sort_unstable();
        assert_eq!(out, expected);
        // Two distinct query node classes should be probed at most twice total.
        assert!(pred.node_mask_calls <= 2);
        if expect_zero_match_edge {
            assert_eq!(pred.match_edge_calls, 0);
        }

        pred.match_edge_calls = 0;
        let cnt = count_prepared(&pt, &pq, None, &mut pred, &mut scratch);
        assert_eq!(cnt, expected.len() as u64);
        if expect_zero_match_edge {
            assert_eq!(pred.match_edge_calls, 0);
        }
    }

    #[rstest]
    #[case::cross_inversion_counterexample(
        Graph::new(4, &[[0, 2], [1, 3]]),
        Graph::new(2, &[[0, 1]])
    )]
    #[case::c4_in_c4(
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]])
    )]
    #[case::path3_in_c6(
        Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]),
        Graph::new(3, &[[0, 1], [1, 2]])
    )]
    #[case::square_in_cubane(
        Graph::new(
            8,
            &[
                [0, 1], [1, 2], [2, 3], [3, 0],
                [4, 5], [5, 6], [6, 7], [7, 4],
                [0, 4], [1, 5], [2, 6], [3, 7],
            ],
        ),
        Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]])
    )]
    #[case::hexagon_in_naphthalene(
        Graph::new(
            10,
            &[
                [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0],
                [3, 6], [6, 7], [7, 8], [8, 9], [9, 4],
            ],
        ),
        Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]])
    )]
    #[case::k3_in_k5(
        Graph::new(
            5,
            &[[0, 1], [0, 2], [0, 3], [0, 4], [1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]],
        ),
        Graph::new(3, &[[0, 1], [1, 2], [2, 0]])
    )]
    #[case::disconnected_two_edges_in_c6(
        Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]),
        Graph::new(4, &[[0, 1], [2, 3]])
    )]
    #[case::u128_c70_path3(
        Graph::new(70, &(0..70u32).map(|i| [i, (i + 1) % 70]).collect::<Vec<_>>()),
        Graph::new(3, &[[0, 1], [1, 2]])
    )]
    #[case::wide_two_c72_path3(
        Graph::new(
            144,
            &(0..72u32)
                .map(|i| [i, (i + 1) % 72])
                .chain((0..72u32).map(|i| [72 + i, 72 + (i + 1) % 72]))
                .collect::<Vec<_>>()
        ),
        Graph::new(3, &[[0, 1], [1, 2]])
    )]
    fn test_visit_prepared_symmetry(#[case] target: Graph, #[case] query: Graph) {
        let mut nm = |_: NodeId, _: NodeId| true;
        let mut em = |_: EdgeId, _: EdgeId| true;
        let all = collect_vf2rdkit(&target, &query, None, &mut nm, &mut em);

        let pt = PreparedSubgraphTarget::new(&target);
        let pq_base = PreparedSubgraphQuery::new(&query, SubgraphIsomorphismAlgorithm::BitsetRi);
        let priority = pq_base.priority_order().to_vec();
        let q_group = query.automorphisms(|_| 0u8, |_| 0u8, AutomorphismAlgorithm::TrailIr);
        let t_group = target.automorphisms(|_| 0u8, |_| 0u8, AutomorphismAlgorithm::TrailIr);

        let mut scratch = SubgraphIsomorphismScratch::new();
        let mut closure_pred = ClosureSubgraphIsomorphismPredicate::new(|_, _| true, |_, _| true);
        let mut undetermined_pred = UndeterminedPredicate;

        for aut_alg in [AutomorphismAlgorithm::TrailIr, AutomorphismAlgorithm::Nauty] {
            let q_sym = QuerySymmetryConstraints::from_graph(&query, |_| 0u8, |_| 0u8, aut_alg);
            let t_sym = TargetSymmetry::from_graph(&target, |_| 0u8, |_| 0u8, aut_alg);

            for (use_q, use_h) in [(false, false), (true, false), (false, true), (true, true)] {
                let pq = if use_q {
                    pq_base.clone().with_symmetry(q_sym.clone())
                } else {
                    pq_base.clone()
                };
                let ts = use_h.then_some(&t_sym);
                let expected = canonical_orbit_reps(
                    &all,
                    &priority,
                    if use_q { q_group.generators() } else { &[] },
                    if use_h { t_group.generators() } else { &[] },
                );

                let mut got = Vec::new();
                let _: ControlFlow<()> =
                    visit_prepared(&pt, &pq, ts, &mut closure_pred, &mut scratch, &mut |emb| {
                        got.push(emb.iter().map(|n| n.index()).collect::<Vec<_>>());
                        ControlFlow::Continue(())
                    });
                got.sort_unstable();
                assert_eq!(got, expected, "alg={aut_alg:?} use_q={use_q} use_h={use_h}");
                assert_eq!(
                    count_prepared(&pt, &pq, ts, &mut closure_pred, &mut scratch),
                    expected.len() as u64,
                    "count closure alg={aut_alg:?} use_q={use_q} use_h={use_h}"
                );
                assert_eq!(
                    count_prepared(&pt, &pq, ts, &mut undetermined_pred, &mut scratch),
                    expected.len() as u64,
                    "count undetermined alg={aut_alg:?} use_q={use_q} use_h={use_h}"
                );
            }
        }
    }

    #[rstest]
    fn test_visit_prepared_symmetry_colored() {
        // Target: K_2,4-like symmetric graph where nodes {0,1} have color 1 (rare) and {2..7}
        // have color 0. Query: 4-vertex path 0-1-2-3 where node 0 (deg 1) has rare color 1 and
        // node 1 (deg 2, priority_order[0]) has color 0, forcing RI-DS to map node 0 BEFORE
        // priority_order[0] and exercising multi-rank stabilizer prefix jumps.
        let target = Graph::new(
            8,
            &[
                [0, 2],
                [0, 3],
                [1, 4],
                [1, 5],
                [2, 3],
                [2, 6],
                [3, 7],
                [4, 5],
                [4, 6],
                [5, 7],
                [6, 7],
            ],
        );
        let query = Graph::new(4, &[[0, 1], [1, 2], [2, 3]]);
        let q_color = |q: NodeId| u8::from(q.index() == 0);
        let t_color = |t: NodeId| u8::from(t.index() < 2);

        let mut nm = |q: NodeId, t: NodeId| q_color(q) == t_color(t);
        let mut em = |_: EdgeId, _: EdgeId| true;
        let all = collect_vf2rdkit(&target, &query, None, &mut nm, &mut em);
        assert!(!all.is_empty());

        let q_group = query.automorphisms(q_color, |_| 0u8, AutomorphismAlgorithm::TrailIr);
        let t_group = target.automorphisms(t_color, |_| 0u8, AutomorphismAlgorithm::TrailIr);
        let q_sym = QuerySymmetryConstraints::from_graph(
            &query,
            q_color,
            |_| 0u8,
            AutomorphismAlgorithm::TrailIr,
        );
        let t_sym =
            TargetSymmetry::from_graph(&target, t_color, |_| 0u8, AutomorphismAlgorithm::TrailIr);
        assert!(!t_sym.is_trivial());

        let pt = PreparedSubgraphTarget::new(&target);
        let pq_base = PreparedSubgraphQuery::new(&query, SubgraphIsomorphismAlgorithm::BitsetRi);
        let priority = pq_base.priority_order().to_vec();
        assert_ne!(priority[0], 0);

        let mut scratch = SubgraphIsomorphismScratch::new();
        let mut pred = ClosureSubgraphIsomorphismPredicate::new(
            |q: NodeId, t: NodeId| q_color(q) == t_color(t),
            |_, _| true,
        );

        for (use_q, use_h) in [(false, true), (true, true)] {
            let pq = if use_q {
                pq_base.clone().with_symmetry(q_sym.clone())
            } else {
                pq_base.clone()
            };
            let expected = canonical_orbit_reps(
                &all,
                &priority,
                if use_q { q_group.generators() } else { &[] },
                if use_h { t_group.generators() } else { &[] },
            );
            let mut got = Vec::new();
            let _: ControlFlow<()> = visit_prepared(
                &pt,
                &pq,
                Some(&t_sym),
                &mut pred,
                &mut scratch,
                &mut |emb| {
                    got.push(emb.iter().map(|n| n.index()).collect::<Vec<_>>());
                    ControlFlow::Continue(())
                },
            );
            got.sort_unstable();
            assert_eq!(got, expected);
            assert_eq!(
                count_prepared(&pt, &pq, Some(&t_sym), &mut pred, &mut scratch),
                expected.len() as u64
            );
        }
    }

    #[rstest]
    fn test_visit_prepared_symmetry_small_graphs() {
        let graphs: Vec<Graph> = SIMPLE_GRAPHS_GRAPH6
            .lines()
            .filter(|line| !line.is_empty())
            .map(parse_graph6)
            .filter(|g| g.node_count() <= 5)
            .collect();
        assert_eq!(graphs.len(), 52);

        let prepared_targets: Vec<PreparedSubgraphTarget> =
            graphs.iter().map(PreparedSubgraphTarget::new).collect();
        let target_groups: Vec<_> = graphs
            .iter()
            .map(|g| g.automorphisms(|_| 0u8, |_| 0u8, AutomorphismAlgorithm::TrailIr))
            .collect();
        let target_syms: Vec<TargetSymmetry> = target_groups
            .iter()
            .map(TargetSymmetry::from_automorphism_group)
            .collect();
        let query_syms: Vec<QuerySymmetryConstraints> = graphs
            .iter()
            .map(|g| {
                QuerySymmetryConstraints::from_graph(
                    g,
                    |_| 0u8,
                    |_| 0u8,
                    AutomorphismAlgorithm::TrailIr,
                )
            })
            .collect();
        let prepared_queries_plain: Vec<PreparedSubgraphQuery> = graphs
            .iter()
            .map(|g| PreparedSubgraphQuery::new(g, SubgraphIsomorphismAlgorithm::BitsetRi))
            .collect();
        let prepared_queries_sym: Vec<PreparedSubgraphQuery> = prepared_queries_plain
            .iter()
            .zip(query_syms.iter())
            .map(|(pq, sym)| pq.clone().with_symmetry(sym.clone()))
            .collect();

        let mut scratch = SubgraphIsomorphismScratch::new();
        let mut pred = ClosureSubgraphIsomorphismPredicate::new(|_, _| true, |_, _| true);

        for (ti, target) in graphs.iter().enumerate() {
            for (qi, query) in graphs.iter().enumerate() {
                if query.node_count() > target.node_count() {
                    continue;
                }
                let mut nm = |_: NodeId, _: NodeId| true;
                let mut em = |_: EdgeId, _: EdgeId| true;
                let all = collect_vf2rdkit(target, query, None, &mut nm, &mut em);
                let priority = prepared_queries_plain[qi].priority_order();

                for (use_q, use_h) in [(false, false), (true, false), (false, true), (true, true)] {
                    let pq = if use_q {
                        &prepared_queries_sym[qi]
                    } else {
                        &prepared_queries_plain[qi]
                    };
                    let ts = use_h.then_some(&target_syms[ti]);
                    let expected = canonical_orbit_reps(
                        &all,
                        priority,
                        if use_q {
                            target_groups[qi].generators()
                        } else {
                            &[]
                        },
                        if use_h {
                            target_groups[ti].generators()
                        } else {
                            &[]
                        },
                    );
                    let mut got = Vec::new();
                    let _: ControlFlow<()> = visit_prepared(
                        &prepared_targets[ti],
                        pq,
                        ts,
                        &mut pred,
                        &mut scratch,
                        &mut |emb| {
                            got.push(emb.iter().map(|n| n.index()).collect::<Vec<_>>());
                            ControlFlow::Continue(())
                        },
                    );
                    got.sort_unstable();
                    assert_eq!(got, expected, "ti={ti} qi={qi} q={use_q} h={use_h}");
                    assert_eq!(
                        count_prepared(&prepared_targets[ti], pq, ts, &mut pred, &mut scratch),
                        expected.len() as u64,
                        "count ti={ti} qi={qi} q={use_q} h={use_h}"
                    );
                }
            }
        }
    }

    #[rstest]
    fn test_visit_prepared_scratch_reuse() {
        let mut scratch = SubgraphIsomorphismScratch::new();
        let mut pred = ClosureSubgraphIsomorphismPredicate::new(|_, _| true, |_, _| true);

        let pairs = [
            (
                Graph::new(
                    8,
                    &[
                        [0, 1],
                        [1, 2],
                        [2, 3],
                        [3, 0],
                        [4, 5],
                        [5, 6],
                        [6, 7],
                        [7, 4],
                        [0, 4],
                        [1, 5],
                        [2, 6],
                        [3, 7],
                    ],
                ),
                Graph::new(4, &[[0, 1], [1, 2], [2, 3], [3, 0]]),
                true,
            ),
            (
                Graph::new(
                    70,
                    &(0..70u32).map(|i| [i, (i + 1) % 70]).collect::<Vec<_>>(),
                ),
                Graph::new(3, &[[0, 1], [1, 2]]),
                false,
            ),
            (
                Graph::new(
                    144,
                    &[[0, 1], [1, 2], [2, 0], [130, 131], [131, 132], [132, 130]],
                ),
                Graph::new(3, &[[0, 1], [1, 2], [2, 0]]),
                true,
            ),
            (
                Graph::new(
                    10,
                    &[
                        [0, 1],
                        [1, 2],
                        [2, 3],
                        [3, 4],
                        [4, 5],
                        [5, 0],
                        [3, 6],
                        [6, 7],
                        [7, 8],
                        [8, 9],
                        [9, 4],
                    ],
                ),
                Graph::new(6, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]),
                false,
            ),
        ];

        for (target, query, use_sym) in pairs {
            let pt = PreparedSubgraphTarget::new(&target);
            let mut pq = PreparedSubgraphQuery::new(&query, SubgraphIsomorphismAlgorithm::BitsetRi);
            let t_sym = if use_sym {
                pq = pq.with_symmetry(QuerySymmetryConstraints::from_graph(
                    &query,
                    |_| 0u8,
                    |_| 0u8,
                    AutomorphismAlgorithm::TrailIr,
                ));
                Some(TargetSymmetry::from_graph(
                    &target,
                    |_| 0u8,
                    |_| 0u8,
                    AutomorphismAlgorithm::TrailIr,
                ))
            } else {
                None
            };

            // Break early on first match to verify scratch/chain state does not leak.
            let early: ControlFlow<Vec<NodeId>> = visit_prepared(
                &pt,
                &pq,
                t_sym.as_ref(),
                &mut pred,
                &mut scratch,
                &mut |emb| ControlFlow::Break(emb.to_vec()),
            );
            assert!(matches!(early, ControlFlow::Break(_)));

            let mut reused_out = Vec::new();
            let _: ControlFlow<()> = visit_prepared(
                &pt,
                &pq,
                t_sym.as_ref(),
                &mut pred,
                &mut scratch,
                &mut |emb| {
                    reused_out.push(emb.to_vec());
                    ControlFlow::Continue(())
                },
            );
            let mut fresh = SubgraphIsomorphismScratch::new();
            let mut fresh_out = Vec::new();
            let _: ControlFlow<()> = visit_prepared(
                &pt,
                &pq,
                t_sym.as_ref(),
                &mut pred,
                &mut fresh,
                &mut |emb| {
                    fresh_out.push(emb.to_vec());
                    ControlFlow::Continue(())
                },
            );
            assert_eq!(reused_out, fresh_out);
            assert_eq!(
                count_prepared(&pt, &pq, t_sym.as_ref(), &mut pred, &mut scratch),
                fresh_out.len() as u64
            );
        }
    }

    #[rstest]
    #[should_panic(expected = "target symmetry node count must match")]
    fn test_count_prepared_error() {
        let pt = PreparedSubgraphTarget::new(&Graph::new(3, &[[0, 1], [1, 2]]));
        let pq =
            PreparedSubgraphQuery::new(&Graph::new(0, &[]), SubgraphIsomorphismAlgorithm::BitsetRi);
        let wrong_sym = TargetSymmetry::none(4);
        let mut scratch = SubgraphIsomorphismScratch::new();
        let mut pred = ClosureSubgraphIsomorphismPredicate::new(|_, _| true, |_, _| true);
        let _ = count_prepared(&pt, &pq, Some(&wrong_sym), &mut pred, &mut scratch);
    }
}
