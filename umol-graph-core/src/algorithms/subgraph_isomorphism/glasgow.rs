//! Glasgow-style constraint-programming subgraph-isomorphism engine (`Glasgow`).
//!
//! Combines compact flat-trail (`u64`/`u128`) and sparse undo-trail (`Wide`) bit-parallel
//! forward checking, dynamic MRV variable selection, unit propagation, Glasgow bit-parallel
//! Hall/AllDifferent wipeout pruning (McCreesh, Prosser, and Trimble 2020), push-down query
//! and target symmetry reduction, anchored search, and leaf `POPCNT` counting.

use std::cmp::Reverse;
use std::mem;
use std::ops::ControlFlow;

use super::bitset::{
    degree_sequence_ok, init_wide_candidates, init_word_candidates, LazyWordAdj,
    PrecomputedWordAdj, WideAdj, WordAdj, WordMask, WIDE_LAZY_ADJ_LIMIT,
};
use super::predicate::SubgraphIsomorphismPredicate;
use super::prepared::{
    PreparedSubgraphQuery, PreparedSubgraphTarget, SubgraphIsomorphismScratch, TargetBitsets,
};
use super::symmetry::{
    is_canonical_double_orbit_with_built_chain, QuerySymmetryConstraints, TargetStabilizerChain,
    TargetSymmetry,
};
use crate::graph::Graph;
use crate::id::NodeId;

/// Maximum query node count using a stack-allocated flat trail `(n1 + 1) * n1` in direct calls.
const STACK_TRAIL_MAX_Q: usize = 16;
const STACK_TRAIL_LEN: usize = (STACK_TRAIL_MAX_Q + 1) * STACK_TRAIL_MAX_Q;

/// Active dynamic push-down symmetry state for a single `Glasgow` search.
struct DynamicSymmetry<'a> {
    has_q_sym: bool,
    has_h_sym: bool,
    prepared_query: &'a PreparedSubgraphQuery,
    query_sym: Option<&'a QuerySymmetryConstraints>,
    target_sym: Option<&'a TargetSymmetry>,
    target_chain: &'a mut TargetStabilizerChain,
}

impl<'a> DynamicSymmetry<'a> {
    fn build(
        prepared_query: Option<&'a PreparedSubgraphQuery>,
        target_sym: Option<&'a TargetSymmetry>,
        target_chain: &'a mut TargetStabilizerChain,
    ) -> Option<Self> {
        let pq = prepared_query?;
        let query_sym = (!pq.symmetry().is_trivial()).then_some(pq.symmetry());
        let active_target_sym = target_sym.filter(|sym| !sym.is_trivial());
        let has_q_sym = query_sym.is_some();
        let has_h_sym = active_target_sym.is_some();
        if !has_q_sym && !has_h_sym {
            return None;
        }
        if let Some(ts) = active_target_sym {
            target_chain.reset(ts);
        }
        Some(Self {
            has_q_sym,
            has_h_sym,
            prepared_query: pq,
            query_sym,
            target_sym: active_target_sym,
            target_chain,
        })
    }

    #[inline(always)]
    fn can_leaf_popcnt(&self) -> bool {
        !(self.has_q_sym && self.has_h_sym)
            && (!self.has_h_sym || self.target_chain.is_current_trivial())
    }

    #[inline(always)]
    fn rejects_before_edge_check(&self, prefix_len: usize, best_q: usize, v: usize) -> bool {
        if !self.has_h_sym || self.target_chain.is_current_trivial() {
            return false;
        }
        let priority = self.prepared_query.priority_order();
        if prefix_len < priority.len() && priority[prefix_len] as usize == best_q {
            !self.target_chain.is_orbit_min_at_current(v)
        } else {
            false
        }
    }

    /// Advance `target_chain` across all priority vertices now mapped in `mapping` starting from
    /// `prefix_len`. Returns `Some((next_prefix_len, saved_depth))` (`saved_depth == usize::MAX`
    /// when no level was pushed), or `None` if pruned (with `target_chain` restored).
    #[inline]
    fn try_advance_host_chain(
        &mut self,
        prefix_len: usize,
        mapping: &[usize],
    ) -> Option<(usize, usize)> {
        if !self.has_h_sym {
            return Some((prefix_len, usize::MAX));
        }
        let priority = self.prepared_query.priority_order();
        let n1 = priority.len();
        let mut next_prefix_len = prefix_len;
        while next_prefix_len < n1 && mapping[priority[next_prefix_len] as usize] != usize::MAX {
            next_prefix_len += 1;
        }
        if next_prefix_len == prefix_len || self.target_chain.is_current_trivial() {
            return Some((next_prefix_len, usize::MAX));
        }

        let saved_depth = self.target_chain.depth();
        let mut pushed = false;
        for &p_r in &priority[prefix_len..next_prefix_len] {
            if self.target_chain.is_current_trivial() {
                break;
            }
            let r_is_last = next_prefix_len == n1 && p_r == priority[n1 - 1];
            let t_r = mapping[p_r as usize];
            if !self.has_q_sym && r_is_last {
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
        Some((
            next_prefix_len,
            if pushed { saved_depth } else { usize::MAX },
        ))
    }

    #[inline(always)]
    fn restore_host_chain(&mut self, saved_depth: usize) {
        if saved_depth != usize::MAX {
            self.target_chain.truncate(saved_depth);
        }
    }

    #[inline]
    fn accepts_complete_embedding(&self, mapping: &[usize]) -> bool {
        if !(self.has_q_sym && self.has_h_sym) {
            return true;
        }
        is_canonical_double_orbit_with_built_chain(
            mapping,
            self.prepared_query.priority_order(),
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

    let mut stack_trail = [W::ZERO; STACK_TRAIL_LEN];
    let mut heap_trail = Vec::new();
    let trail_slice: &mut [W] = if n1 <= STACK_TRAIL_MAX_Q {
        &mut stack_trail[..(n1 + 1) * n1]
    } else {
        heap_trail.resize((n1 + 1) * n1, W::ZERO);
        &mut heap_trail
    };

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
        &mut trail_slice[..n1],
        None,
    ) {
        return ControlFlow::Continue(());
    }

    let mut adj_buf = [W::ZERO; N];
    let mut adj = LazyWordAdj {
        target,
        adj: &mut adj_buf[..n2],
        built: W::ZERO,
    };
    let mut mapping = [usize::MAX; N];
    let mut out_nodes = [NodeId(0); N];

    let mut ctx = WordGlasgowCtx {
        target,
        query,
        n2,
        q_deg: &q_deg[..n1],
        all_edges_undetermined: pred.all_edges_undetermined(),
        adj: &mut adj,
        trail: trail_slice,
        mapping: &mut mapping[..n1],
        out_nodes: &mut out_nodes[..n1],
        sym: None,
        pred,
        count,
        visitor,
    };

    let mut init_unit_done = W::ZERO;
    let mut init_verified_edges = 0u128;
    if n1 > 1
        && !propagate_word::<W, N, _, _, _>(
            &mut ctx,
            0,
            W::ZERO,
            &mut init_unit_done,
            &mut init_verified_edges,
        )
    {
        return ControlFlow::Continue(());
    }

    word_dfs::<W, N, _, COUNT_ONLY, B, P>(
        &mut ctx,
        0,
        W::ZERO,
        init_unit_done,
        init_verified_edges,
        0,
    )
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

    let mut trail = W::take_trail(scratch);
    trail.clear();
    trail.resize((n1 + 1) * n1, W::ZERO);

    let ok = init_word_candidates::<W, N, P>(
        n1,
        &q_deg,
        &q_by_deg,
        None,
        root_orbit_filter,
        |qd| deg_ge[qd],
        pred,
        &mut trail[..n1],
        None,
    );
    if !ok {
        W::restore_trail(scratch, trail);
        return ControlFlow::Continue(());
    }

    scratch.out_nodes.resize(n1, NodeId(0));
    let sym = DynamicSymmetry::build(Some(prepared_query), target_sym, &mut scratch.target_chain);

    let mut adj = PrecomputedWordAdj(adj_rows);
    let mut mapping = [usize::MAX; N];
    let res = {
        let mut ctx = WordGlasgowCtx {
            target,
            query,
            n2,
            q_deg: &q_deg[..n1],
            all_edges_undetermined: pred.all_edges_undetermined(),
            adj: &mut adj,
            trail: &mut trail,
            mapping: &mut mapping[..n1],
            out_nodes: &mut scratch.out_nodes,
            sym,
            pred,
            count,
            visitor,
        };

        let mut init_unit_done = W::ZERO;
        let mut init_verified_edges = 0u128;
        if n1 > 1
            && !propagate_word::<W, N, _, _, _>(
                &mut ctx,
                0,
                W::ZERO,
                &mut init_unit_done,
                &mut init_verified_edges,
            )
        {
            ControlFlow::Continue(())
        } else {
            word_dfs::<W, N, _, COUNT_ONLY, B, P>(
                &mut ctx,
                0,
                W::ZERO,
                init_unit_done,
                init_verified_edges,
                0,
            )
        }
    };

    W::restore_trail(scratch, trail);
    res
}

struct WordGlasgowCtx<'a, W, A, P, F> {
    target: &'a Graph,
    query: &'a Graph,
    n2: usize,
    q_deg: &'a [usize],
    all_edges_undetermined: bool,
    adj: &'a mut A,
    trail: &'a mut [W],
    mapping: &'a mut [usize],
    out_nodes: &'a mut [NodeId],
    sym: Option<DynamicSymmetry<'a>>,
    pred: &'a mut P,
    count: &'a mut u64,
    visitor: &'a mut F,
}

/// Run unit propagation (AllDifferent + adjacency/edge FC + Grochow--Kellis bounds) and Glasgow
/// bit-parallel Hall/AllDifferent propagation on `ctx.trail[base..base + n1]`.
///
/// Returns `false` on domain wipeout or Hall violation.
fn propagate_word<W, const N: usize, A, P, F>(
    ctx: &mut WordGlasgowCtx<'_, W, A, P, F>,
    base: usize,
    assigned_q: W,
    unit_done: &mut W,
    verified_edges: &mut u128,
) -> bool
where
    W: WordMask,
    A: WordAdj<W>,
    P: SubgraphIsomorphismPredicate,
{
    let n1 = ctx.q_deg.len();
    let unassigned_mask = W::all_below(n1) & !assigned_q;
    let mut unit_queue = W::ZERO;

    let mut scan = unassigned_mask & !*unit_done;
    while !scan.is_zero() {
        let q = scan.trailing_zeros() as usize;
        scan = scan.clear_lowest();
        let d = ctx.trail[base + q];
        if d.is_zero() {
            return false;
        }
        if d.count_ones() == 1 {
            unit_queue = unit_queue | W::bit(q);
        }
    }

    let pq_q_sym = ctx
        .sym
        .as_ref()
        .filter(|s| s.has_q_sym)
        .map(|s| s.prepared_query);

    loop {
        // Phase 1: Drain forced singletons via AllDifferent + adjacency/edge FC + G_Q bounds.
        while !unit_queue.is_zero() {
            let u = unit_queue.trailing_zeros() as usize;
            unit_queue = unit_queue.clear_lowest();
            if unit_done.test(u) {
                continue;
            }
            let du = ctx.trail[base + u];
            if du.is_zero() {
                return false;
            }
            *unit_done = *unit_done | W::bit(u);
            let w = du.trailing_zeros() as usize;
            let not_w = !W::bit(w);

            let mut rem_q = unassigned_mask & !W::bit(u);
            while !rem_q.is_zero() {
                let q = rem_q.trailing_zeros() as usize;
                rem_q = rem_q.clear_lowest();
                let old_d = ctx.trail[base + q];
                let new_d = old_d & not_w;
                if new_d != old_d {
                    if new_d.is_zero() {
                        return false;
                    }
                    ctx.trail[base + q] = new_d;
                    if new_d.count_ones() == 1 && !unit_done.test(q) {
                        unit_queue = unit_queue | W::bit(q);
                    }
                }
            }

            let mut aw_cached: Option<W> = None;
            let w_node = NodeId(w as u32);
            for nb in ctx.query.neighbors(NodeId(u as u32)) {
                let x = nb.node.index();
                if assigned_q.test(x) {
                    continue;
                }
                let pmask = if ctx.all_edges_undetermined || ctx.pred.edge_is_undetermined(nb.edge)
                {
                    if nb.edge.index() < 128 {
                        *verified_edges |= 1u128 << nb.edge.index();
                    }
                    *aw_cached.get_or_insert_with(|| ctx.adj.row(w))
                } else if let Some(emask) = W::edge_neighbors(ctx.pred, nb.edge, w_node) {
                    if nb.edge.index() < 128 {
                        *verified_edges |= 1u128 << nb.edge.index();
                    }
                    emask
                } else {
                    *aw_cached.get_or_insert_with(|| ctx.adj.row(w))
                };
                let old_d = ctx.trail[base + x];
                let new_d = old_d & pmask;
                if new_d != old_d {
                    if new_d.is_zero() {
                        return false;
                    }
                    ctx.trail[base + x] = new_d;
                    if new_d.count_ones() == 1 && !unit_done.test(x) {
                        unit_queue = unit_queue | W::bit(x);
                    }
                }
            }

            if let Some(pq) = pq_q_sym {
                let gt_mask = W::range_mask(w + 1, ctx.n2);
                for &hi in pq.sym_upper_bounds(u) {
                    let x = hi as usize;
                    if assigned_q.test(x) {
                        continue;
                    }
                    let old_d = ctx.trail[base + x];
                    let new_d = old_d & gt_mask;
                    if new_d != old_d {
                        if new_d.is_zero() {
                            return false;
                        }
                        ctx.trail[base + x] = new_d;
                        if new_d.count_ones() == 1 && !unit_done.test(x) {
                            unit_queue = unit_queue | W::bit(x);
                        }
                    }
                }
                let lt_mask = W::range_mask(0, w);
                for &lo in pq.sym_lower_bounds(u) {
                    let x = lo as usize;
                    if assigned_q.test(x) {
                        continue;
                    }
                    let old_d = ctx.trail[base + x];
                    let new_d = old_d & lt_mask;
                    if new_d != old_d {
                        if new_d.is_zero() {
                            return false;
                        }
                        ctx.trail[base + x] = new_d;
                        if new_d.count_ones() == 1 && !unit_done.test(x) {
                            unit_queue = unit_queue | W::bit(x);
                        }
                    }
                }
            }
        }

        // Phase 2: Bit-parallel Hall / AllDifferent sweep over remaining non-unit unassigned vars.
        let non_unit_mask = unassigned_mask & !*unit_done;
        let k_nu = non_unit_mask.count_ones();
        if k_nu <= 1 {
            return true;
        }

        let mut union_all = W::ZERO;
        let mut small = [(0usize, 0u32); N];
        let mut small_len = 0usize;
        let mut large_mask = W::ZERO;
        let mut rem = non_unit_mask;
        while !rem.is_zero() {
            let q = rem.trailing_zeros() as usize;
            rem = rem.clear_lowest();
            let d = ctx.trail[base + q];
            union_all = union_all | d;
            let cnt = d.count_ones();
            if cnt < k_nu {
                small[small_len] = (q, cnt);
                small_len += 1;
            } else {
                large_mask = large_mask | W::bit(q);
            }
        }
        if union_all.count_ones() < k_nu {
            return false;
        }
        if small_len < 2 {
            return true;
        }

        // Insertion sort small[..small_len] by ascending domain size.
        for i in 1..small_len {
            let key = small[i];
            let mut j = i;
            while j > 0 && small[j - 1].1 > key.1 {
                small[j] = small[j - 1];
                j -= 1;
            }
            small[j] = key;
        }

        let mut left = W::ZERO;
        let mut hall = W::ZERO;
        let mut processed = 0u32;
        for &(q, _) in &small[..small_len] {
            let old_d = ctx.trail[base + q];
            let d = old_d & !hall;
            if d.is_zero() {
                return false;
            }
            if d != old_d {
                ctx.trail[base + q] = d;
                if d.count_ones() == 1 && !unit_done.test(q) {
                    unit_queue = unit_queue | W::bit(q);
                }
            }
            left = left | d;
            processed += 1;
            let left_sz = left.count_ones();
            if left_sz < processed {
                return false;
            }
            if left_sz == processed {
                hall = hall | left;
            }
        }

        if !hall.is_zero() {
            let not_hall = !hall;
            let mut rem_large = large_mask;
            while !rem_large.is_zero() {
                let q = rem_large.trailing_zeros() as usize;
                rem_large = rem_large.clear_lowest();
                let old_d = ctx.trail[base + q];
                let d = old_d & not_hall;
                if d.is_zero() {
                    return false;
                }
                if d != old_d {
                    ctx.trail[base + q] = d;
                    if d.count_ones() == 1 && !unit_done.test(q) {
                        unit_queue = unit_queue | W::bit(q);
                    }
                }
            }
        }

        if unit_queue.is_zero() {
            return true;
        }
    }
}

fn word_dfs<W, const N: usize, A, const COUNT_ONLY: bool, B, P>(
    ctx: &mut WordGlasgowCtx<'_, W, A, P, impl FnMut(&[NodeId]) -> ControlFlow<B>>,
    depth: usize,
    assigned_q: W,
    unit_done: W,
    verified_edges: u128,
    prefix_len: usize,
) -> ControlFlow<B>
where
    W: WordMask,
    A: WordAdj<W>,
    P: SubgraphIsomorphismPredicate,
{
    let n1 = ctx.q_deg.len();
    let base = depth * n1;

    // Dynamic MRV: pick unassigned query vertex with smallest live domain, tie-break by max
    // degree then smallest node id.
    let mut best_q = usize::MAX;
    let mut best_cnt = u32::MAX;
    let mut best_deg = 0usize;
    let mut rem_q = W::all_below(n1) & !assigned_q;
    while !rem_q.is_zero() {
        let q = rem_q.trailing_zeros() as usize;
        rem_q = rem_q.clear_lowest();
        let cnt = ctx.trail[base + q].count_ones();
        if cnt == 0 {
            return ControlFlow::Continue(());
        }
        let deg = ctx.q_deg[q];
        if cnt < best_cnt || (cnt == best_cnt && deg > best_deg) {
            best_q = q;
            best_cnt = cnt;
            best_deg = deg;
            if cnt == 1 {
                break;
            }
        }
    }

    let mut cands = ctx.trail[base + best_q];
    let q_nbs = ctx.query.neighbors(NodeId(best_q as u32));
    let is_last = depth + 1 == n1;

    let all_incident_verified = ctx.all_edges_undetermined
        || q_nbs.iter().all(|nb| {
            !assigned_q.test(nb.node.index())
                || (nb.edge.index() < 128 && ((verified_edges >> nb.edge.index()) & 1) != 0)
                || ctx.pred.edge_is_undetermined(nb.edge)
        });

    if COUNT_ONLY
        && is_last
        && all_incident_verified
        && ctx
            .sym
            .as_ref()
            .is_none_or(DynamicSymmetry::can_leaf_popcnt)
    {
        *ctx.count += cands.count_ones() as u64;
        return ControlFlow::Continue(());
    }

    let next_assigned_q = assigned_q | W::bit(best_q);
    let next_unassigned = W::all_below(n1) & !next_assigned_q;
    let next_base = (depth + 1) * n1;

    while !cands.is_zero() {
        let v = cands.trailing_zeros() as usize;
        cands = cands.clear_lowest();

        if ctx
            .sym
            .as_ref()
            .is_some_and(|sym| sym.rejects_before_edge_check(prefix_len, best_q, v))
        {
            continue;
        }

        let v_node = NodeId(v as u32);
        if !all_incident_verified {
            let mut edges_ok = true;
            for nb in q_nbs {
                let u = nb.node.index();
                if !assigned_q.test(u) {
                    continue;
                }
                if (nb.edge.index() < 128 && ((verified_edges >> nb.edge.index()) & 1) != 0)
                    || ctx.pred.edge_is_undetermined(nb.edge)
                {
                    continue;
                }
                let te = ctx
                    .target
                    .find_edge(v_node, NodeId(ctx.mapping[u] as u32))
                    .expect("forward checking guarantees adjacency");
                if !ctx.pred.match_edge(nb.edge, te) {
                    edges_ok = false;
                    break;
                }
            }
            if !edges_ok {
                continue;
            }
        }

        ctx.mapping[best_q] = v;
        let (next_prefix_len, saved_depth) = if let Some(sym) = ctx.sym.as_mut() {
            let Some(pair) = sym.try_advance_host_chain(prefix_len, ctx.mapping) else {
                ctx.mapping[best_q] = usize::MAX;
                continue;
            };
            pair
        } else {
            (prefix_len, usize::MAX)
        };

        if is_last {
            let accepted = ctx
                .sym
                .as_ref()
                .is_none_or(|sym| sym.accepts_complete_embedding(ctx.mapping));
            let res = if accepted {
                if COUNT_ONLY {
                    *ctx.count += 1;
                    ControlFlow::Continue(())
                } else {
                    for (slot, &target_idx) in ctx.out_nodes.iter_mut().zip(ctx.mapping.iter()) {
                        *slot = NodeId(target_idx as u32);
                    }
                    (ctx.visitor)(ctx.out_nodes)
                }
            } else {
                ControlFlow::Continue(())
            };
            if let Some(sym) = ctx.sym.as_mut() {
                sym.restore_host_chain(saved_depth);
            }
            ctx.mapping[best_q] = usize::MAX;
            if let ControlFlow::Break(val) = res {
                return ControlFlow::Break(val);
            }
            continue;
        }

        // Copy compact n1-word domain frame while clearing `v` (AllDifferent).
        let not_v = !W::bit(v);
        let mut wipeout = false;
        let mut rem = next_unassigned;
        while !rem.is_zero() {
            let q = rem.trailing_zeros() as usize;
            rem = rem.clear_lowest();
            let d = ctx.trail[base + q] & not_v;
            if d.is_zero() {
                wipeout = true;
                break;
            }
            ctx.trail[next_base + q] = d;
        }

        let mut next_verified_edges = verified_edges;
        if !wipeout {
            // Adjacency / exact-edge forward checking to unassigned neighbors of `best_q`.
            let mut av_cached: Option<W> = None;
            for nb in q_nbs {
                let u = nb.node.index();
                if next_assigned_q.test(u) {
                    continue;
                }
                let pmask = if ctx.all_edges_undetermined || ctx.pred.edge_is_undetermined(nb.edge)
                {
                    if nb.edge.index() < 128 {
                        next_verified_edges |= 1u128 << nb.edge.index();
                    }
                    *av_cached.get_or_insert_with(|| ctx.adj.row(v))
                } else if let Some(emask) = W::edge_neighbors(ctx.pred, nb.edge, v_node) {
                    if nb.edge.index() < 128 {
                        next_verified_edges |= 1u128 << nb.edge.index();
                    }
                    emask
                } else {
                    *av_cached.get_or_insert_with(|| ctx.adj.row(v))
                };
                let d = ctx.trail[next_base + u] & pmask;
                if d.is_zero() {
                    wipeout = true;
                    break;
                }
                ctx.trail[next_base + u] = d;
            }
        }

        if !wipeout {
            // Push-down Grochow--Kellis G_Q inequalities and newly stabilized G_H orbit reps.
            if let Some(sym) = ctx.sym.as_ref() {
                if sym.has_q_sym {
                    let gt_mask = W::range_mask(v + 1, ctx.n2);
                    for &hi in sym.prepared_query.sym_upper_bounds(best_q) {
                        let u = hi as usize;
                        if next_assigned_q.test(u) {
                            continue;
                        }
                        let d = ctx.trail[next_base + u] & gt_mask;
                        if d.is_zero() {
                            wipeout = true;
                            break;
                        }
                        ctx.trail[next_base + u] = d;
                    }
                    if !wipeout {
                        let lt_mask = W::range_mask(0, v);
                        for &lo in sym.prepared_query.sym_lower_bounds(best_q) {
                            let u = lo as usize;
                            if next_assigned_q.test(u) {
                                continue;
                            }
                            let d = ctx.trail[next_base + u] & lt_mask;
                            if d.is_zero() {
                                wipeout = true;
                                break;
                            }
                            ctx.trail[next_base + u] = d;
                        }
                    }
                }
                if !wipeout
                    && sym.has_h_sym
                    && next_prefix_len > prefix_len
                    && next_prefix_len < n1
                    && !sym.target_chain.is_current_trivial()
                {
                    let p_next = sym.prepared_query.priority_order()[next_prefix_len] as usize;
                    let mut filtered = W::ZERO;
                    let mut dom_bits = ctx.trail[next_base + p_next];
                    while !dom_bits.is_zero() {
                        let cand_t = dom_bits.trailing_zeros() as usize;
                        dom_bits = dom_bits.clear_lowest();
                        if sym.target_chain.is_orbit_min_at_current(cand_t) {
                            filtered = filtered | W::bit(cand_t);
                        }
                    }
                    if filtered.is_zero() {
                        wipeout = true;
                    } else {
                        ctx.trail[next_base + p_next] = filtered;
                    }
                }
            }
        }

        let mut next_unit_done = unit_done & !W::bit(best_q);
        if !wipeout {
            wipeout = !propagate_word::<W, N, A, P, _>(
                ctx,
                next_base,
                next_assigned_q,
                &mut next_unit_done,
                &mut next_verified_edges,
            );
        }

        let res = if wipeout {
            ControlFlow::Continue(())
        } else {
            word_dfs::<W, N, A, COUNT_ONLY, B, P>(
                ctx,
                depth + 1,
                next_assigned_q,
                next_unit_done,
                next_verified_edges,
                next_prefix_len,
            )
        };

        if let Some(sym) = ctx.sym.as_mut() {
            sym.restore_host_chain(saved_depth);
        }
        ctx.mapping[best_q] = usize::MAX;
        if let ControlFlow::Break(val) = res {
            return ControlFlow::Break(val);
        }
    }

    ControlFlow::Continue(())
}

// Wide (`n2 > 128`) Glasgow engine using O(n1 * words) domain bitsets + sparse word undo trail.

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
    if !init_wide_candidates(
        target,
        &q_deg,
        prepared_query,
        prepared_deg_ge,
        anchor,
        target_sym,
        pred,
        scratch,
    ) {
        return ControlFlow::Continue(());
    }

    let mut dom_wide = mem::take(&mut scratch.cand_wide);
    let mut dom_cnt: Vec<u32> = scratch.cand_counts.iter().map(|&c| c as u32).collect();

    let use_lazy_adj = prepared_adj.is_none() && n2 <= WIDE_LAZY_ADJ_LIMIT;
    if use_lazy_adj {
        scratch.adj_wide.resize(n2 * words, 0);
        scratch.adj_built_wide.clear();
        scratch.adj_built_wide.resize(words, 0);
    }
    let mut adj_wide = mem::take(&mut scratch.adj_wide);
    let mut adj_built_wide = mem::take(&mut scratch.adj_built_wide);
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

    scratch.mapping.clear();
    scratch.mapping.resize(n1, usize::MAX);
    scratch.out_nodes.resize(n1, NodeId(0));

    let sym = DynamicSymmetry::build(prepared_query, target_sym, &mut scratch.target_chain);

    let mut assigned = vec![false; n1];
    let mut unit_done = vec![false; n1];
    let mut undo_words = Vec::new();
    let mut undo_cnts = Vec::new();
    let mut undo_unit_done = Vec::new();
    let mut adj_fallback = vec![0u64; words];
    let mut hall_left = vec![0u64; words];
    let mut hall_set = vec![0u64; words];
    let mut hall_union = vec![0u64; words];

    let res = {
        let mut ctx = WideGlasgowCtx {
            target,
            query,
            n2,
            words,
            q_deg: &q_deg,
            all_edges_undetermined: pred.all_edges_undetermined(),
            adj,
            dom_wide: &mut dom_wide,
            dom_cnt: &mut dom_cnt,
            assigned: &mut assigned,
            unit_done: &mut unit_done,
            undo_words: &mut undo_words,
            undo_cnts: &mut undo_cnts,
            undo_unit_done: &mut undo_unit_done,
            adj_fallback: &mut adj_fallback,
            hall_left: &mut hall_left,
            hall_set: &mut hall_set,
            hall_union: &mut hall_union,
            mapping: &mut scratch.mapping,
            out_nodes: &mut scratch.out_nodes,
            sym,
            pred,
            count,
            visitor,
        };

        if n1 > 1 && !propagate_wide(&mut ctx) {
            ControlFlow::Continue(())
        } else {
            wide_dfs::<COUNT_ONLY, B, P>(&mut ctx, 0, 0)
        }
    };

    scratch.cand_wide = dom_wide;
    scratch.adj_wide = adj_wide;
    scratch.adj_built_wide = adj_built_wide;
    res
}

struct WideGlasgowCtx<'a, P, F> {
    target: &'a Graph,
    query: &'a Graph,
    n2: usize,
    words: usize,
    q_deg: &'a [usize],
    all_edges_undetermined: bool,
    adj: WideAdj<'a>,
    dom_wide: &'a mut [u64],
    dom_cnt: &'a mut [u32],
    assigned: &'a mut [bool],
    unit_done: &'a mut [bool],
    undo_words: &'a mut Vec<(usize, u64)>,
    undo_cnts: &'a mut Vec<(usize, u32)>,
    undo_unit_done: &'a mut Vec<usize>,
    adj_fallback: &'a mut [u64],
    hall_left: &'a mut [u64],
    hall_set: &'a mut [u64],
    hall_union: &'a mut [u64],
    mapping: &'a mut [usize],
    out_nodes: &'a mut [NodeId],
    sym: Option<DynamicSymmetry<'a>>,
    pred: &'a mut P,
    count: &'a mut u64,
    visitor: &'a mut F,
}

impl<P, F> WideGlasgowCtx<'_, P, F> {
    #[inline(always)]
    fn restore_marks(&mut self, word_mark: usize, cnt_mark: usize, unit_mark: usize) {
        while self.undo_words.len() > word_mark {
            let (idx, old) = self.undo_words.pop().expect("len > word_mark");
            self.dom_wide[idx] = old;
        }
        while self.undo_cnts.len() > cnt_mark {
            let (q, old) = self.undo_cnts.pop().expect("len > cnt_mark");
            self.dom_cnt[q] = old;
        }
        while self.undo_unit_done.len() > unit_mark {
            let q = self.undo_unit_done.pop().expect("len > unit_mark");
            self.unit_done[q] = false;
        }
    }

    #[inline(always)]
    fn clear_bit(&mut self, q: usize, v: usize) -> bool {
        let idx = q * self.words + (v >> 6);
        let bit = 1u64 << (v & 63);
        let old = self.dom_wide[idx];
        if (old & bit) == 0 {
            return true;
        }
        self.undo_words.push((idx, old));
        self.dom_wide[idx] = old & !bit;
        let old_cnt = self.dom_cnt[q];
        self.undo_cnts.push((q, old_cnt));
        let new_cnt = old_cnt - 1;
        self.dom_cnt[q] = new_cnt;
        new_cnt > 0
    }

    #[inline]
    fn intersect_row_with_adj(&mut self, q: usize, target_u: usize) -> bool {
        let words = self.words;
        let base = q * words;
        let adj_slice = self
            .adj
            .row_into(self.target, target_u, words, self.adj_fallback);
        let old_cnt = self.dom_cnt[q];
        let mut removed = 0u32;
        for (w, &adj_word) in adj_slice.iter().enumerate() {
            let idx = base + w;
            let old = self.dom_wide[idx];
            let new = old & adj_word;
            if new != old {
                self.undo_words.push((idx, old));
                self.dom_wide[idx] = new;
                removed += (old ^ new).count_ones();
            }
        }
        if removed > 0 {
            self.undo_cnts.push((q, old_cnt));
            self.dom_cnt[q] = old_cnt - removed;
        }
        self.dom_cnt[q] > 0
    }

    #[inline]
    fn narrow_ge(&mut self, q: usize, min_v: usize) -> bool {
        if min_v == 0 {
            return true;
        }
        let words = self.words;
        let base = q * words;
        let old_cnt = self.dom_cnt[q];
        let mut removed = 0u32;
        if min_v >= self.n2 {
            for w in 0..words {
                let idx = base + w;
                let old = self.dom_wide[idx];
                if old != 0 {
                    self.undo_words.push((idx, old));
                    self.dom_wide[idx] = 0;
                    removed += old.count_ones();
                }
            }
        } else {
            let first_w = min_v >> 6;
            let rem = min_v & 63;
            for w in 0..first_w {
                let idx = base + w;
                let old = self.dom_wide[idx];
                if old != 0 {
                    self.undo_words.push((idx, old));
                    self.dom_wide[idx] = 0;
                    removed += old.count_ones();
                }
            }
            if rem > 0 {
                let idx = base + first_w;
                let old = self.dom_wide[idx];
                let new = old & (u64::MAX << rem);
                if new != old {
                    self.undo_words.push((idx, old));
                    self.dom_wide[idx] = new;
                    removed += (old ^ new).count_ones();
                }
            }
        }
        if removed > 0 {
            self.undo_cnts.push((q, old_cnt));
            self.dom_cnt[q] = old_cnt - removed;
        }
        self.dom_cnt[q] > 0
    }

    #[inline]
    fn narrow_lt(&mut self, q: usize, max_v: usize) -> bool {
        if max_v >= self.n2 {
            return true;
        }
        let words = self.words;
        let base = q * words;
        let old_cnt = self.dom_cnt[q];
        let mut removed = 0u32;
        let last_w = max_v >> 6;
        let rem = max_v & 63;
        let zero_start = if rem > 0 {
            let idx = base + last_w;
            let old = self.dom_wide[idx];
            let new = old & ((1u64 << rem) - 1);
            if new != old {
                self.undo_words.push((idx, old));
                self.dom_wide[idx] = new;
                removed += (old ^ new).count_ones();
            }
            last_w + 1
        } else {
            last_w
        };
        for w in zero_start..words {
            let idx = base + w;
            let old = self.dom_wide[idx];
            if old != 0 {
                self.undo_words.push((idx, old));
                self.dom_wide[idx] = 0;
                removed += old.count_ones();
            }
        }
        if removed > 0 {
            self.undo_cnts.push((q, old_cnt));
            self.dom_cnt[q] = old_cnt - removed;
        }
        self.dom_cnt[q] > 0
    }

    #[inline]
    fn first_set_bit(&self, q: usize) -> usize {
        let base = q * self.words;
        for w in 0..self.words {
            let bits = self.dom_wide[base + w];
            if bits != 0 {
                return (w << 6) + bits.trailing_zeros() as usize;
            }
        }
        unreachable!("non-empty wide domain has a set bit")
    }
}

fn propagate_wide<P: SubgraphIsomorphismPredicate, F>(ctx: &mut WideGlasgowCtx<'_, P, F>) -> bool {
    let n1 = ctx.q_deg.len();
    let words = ctx.words;
    let mut unit_queue = Vec::new();
    for q in 0..n1 {
        if !ctx.assigned[q] && !ctx.unit_done[q] {
            let cnt = ctx.dom_cnt[q];
            if cnt == 0 {
                return false;
            }
            if cnt == 1 {
                unit_queue.push(q);
            }
        }
    }

    let pq_q_sym = ctx
        .sym
        .as_ref()
        .filter(|s| s.has_q_sym)
        .map(|s| s.prepared_query);

    loop {
        while let Some(u) = unit_queue.pop() {
            if ctx.unit_done[u] {
                continue;
            }
            if ctx.dom_cnt[u] == 0 {
                return false;
            }
            ctx.unit_done[u] = true;
            ctx.undo_unit_done.push(u);
            let w = ctx.first_set_bit(u);

            for q in 0..n1 {
                if q == u || ctx.assigned[q] {
                    continue;
                }
                let prev_cnt = ctx.dom_cnt[q];
                if !ctx.clear_bit(q, w) {
                    return false;
                }
                if prev_cnt > 1 && ctx.dom_cnt[q] == 1 && !ctx.unit_done[q] {
                    unit_queue.push(q);
                }
            }

            for nb in ctx.query.neighbors(NodeId(u as u32)) {
                let x = nb.node.index();
                if ctx.assigned[x] {
                    continue;
                }
                let prev_cnt = ctx.dom_cnt[x];
                if !ctx.intersect_row_with_adj(x, w) {
                    return false;
                }
                if prev_cnt > 1 && ctx.dom_cnt[x] == 1 && !ctx.unit_done[x] {
                    unit_queue.push(x);
                }
            }

            if let Some(pq) = pq_q_sym {
                for &hi in pq.sym_upper_bounds(u) {
                    let x = hi as usize;
                    if ctx.assigned[x] {
                        continue;
                    }
                    let prev_cnt = ctx.dom_cnt[x];
                    if !ctx.narrow_ge(x, w + 1) {
                        return false;
                    }
                    if prev_cnt > 1 && ctx.dom_cnt[x] == 1 && !ctx.unit_done[x] {
                        unit_queue.push(x);
                    }
                }
                for &lo in pq.sym_lower_bounds(u) {
                    let x = lo as usize;
                    if ctx.assigned[x] {
                        continue;
                    }
                    let prev_cnt = ctx.dom_cnt[x];
                    if !ctx.narrow_lt(x, w) {
                        return false;
                    }
                    if prev_cnt > 1 && ctx.dom_cnt[x] == 1 && !ctx.unit_done[x] {
                        unit_queue.push(x);
                    }
                }
            }
        }

        // Phase 2: Multi-word Glasgow Hall / AllDifferent sweep over non-unit unassigned vars.
        let mut non_unit_count = 0u32;
        for q in 0..n1 {
            if !ctx.assigned[q] && !ctx.unit_done[q] {
                non_unit_count += 1;
            }
        }
        if non_unit_count <= 1 {
            return true;
        }

        ctx.hall_union.fill(0);
        let mut small = Vec::new();
        let mut large = Vec::new();
        for q in 0..n1 {
            if ctx.assigned[q] || ctx.unit_done[q] {
                continue;
            }
            let base = q * words;
            for (dst, &src) in ctx
                .hall_union
                .iter_mut()
                .zip(&ctx.dom_wide[base..base + words])
            {
                *dst |= src;
            }
            let cnt = ctx.dom_cnt[q];
            if cnt < non_unit_count {
                small.push((q, cnt));
            } else {
                large.push(q);
            }
        }
        let union_sz: u32 = ctx.hall_union.iter().map(|w| w.count_ones()).sum();
        if union_sz < non_unit_count {
            return false;
        }
        if small.len() < 2 {
            return true;
        }

        small.sort_unstable_by_key(|&(_, cnt)| cnt);
        ctx.hall_left.fill(0);
        ctx.hall_set.fill(0);
        let mut hall_nonempty = false;
        let mut processed = 0u32;

        for &(q, _) in &small {
            let base = q * words;
            if hall_nonempty {
                let old_cnt = ctx.dom_cnt[q];
                let mut removed = 0u32;
                for (w, &hword) in ctx.hall_set.iter().enumerate() {
                    let idx = base + w;
                    let old = ctx.dom_wide[idx];
                    let new = old & !hword;
                    if new != old {
                        ctx.undo_words.push((idx, old));
                        ctx.dom_wide[idx] = new;
                        removed += (old ^ new).count_ones();
                    }
                }
                if removed > 0 {
                    ctx.undo_cnts.push((q, old_cnt));
                    let new_cnt = old_cnt - removed;
                    ctx.dom_cnt[q] = new_cnt;
                    if new_cnt == 0 {
                        return false;
                    }
                    if new_cnt == 1 && !ctx.unit_done[q] {
                        unit_queue.push(q);
                    }
                }
            }

            let mut left_sz = 0u32;
            for (dst, &src) in ctx
                .hall_left
                .iter_mut()
                .zip(&ctx.dom_wide[base..base + words])
            {
                *dst |= src;
                left_sz += dst.count_ones();
            }
            processed += 1;
            if left_sz < processed {
                return false;
            }
            if left_sz == processed {
                for (h, &l) in ctx.hall_set.iter_mut().zip(ctx.hall_left.iter()) {
                    *h |= l;
                }
                hall_nonempty = true;
            }
        }

        if hall_nonempty {
            for &q in &large {
                let base = q * words;
                let old_cnt = ctx.dom_cnt[q];
                let mut removed = 0u32;
                for (w, &hword) in ctx.hall_set.iter().enumerate() {
                    let idx = base + w;
                    let old = ctx.dom_wide[idx];
                    let new = old & !hword;
                    if new != old {
                        ctx.undo_words.push((idx, old));
                        ctx.dom_wide[idx] = new;
                        removed += (old ^ new).count_ones();
                    }
                }
                if removed > 0 {
                    ctx.undo_cnts.push((q, old_cnt));
                    let new_cnt = old_cnt - removed;
                    ctx.dom_cnt[q] = new_cnt;
                    if new_cnt == 0 {
                        return false;
                    }
                    if new_cnt == 1 && !ctx.unit_done[q] {
                        unit_queue.push(q);
                    }
                }
            }
        }

        if unit_queue.is_empty() {
            return true;
        }
    }
}

fn wide_dfs<const COUNT_ONLY: bool, B, P: SubgraphIsomorphismPredicate>(
    ctx: &mut WideGlasgowCtx<'_, P, impl FnMut(&[NodeId]) -> ControlFlow<B>>,
    depth: usize,
    prefix_len: usize,
) -> ControlFlow<B> {
    let n1 = ctx.q_deg.len();
    let words = ctx.words;

    let mut best_q = usize::MAX;
    let mut best_cnt = u32::MAX;
    let mut best_deg = 0usize;
    for (q, &deg) in ctx.q_deg.iter().enumerate() {
        if ctx.assigned[q] {
            continue;
        }
        let cnt = ctx.dom_cnt[q];
        if cnt == 0 {
            return ControlFlow::Continue(());
        }
        if cnt < best_cnt || (cnt == best_cnt && deg > best_deg) {
            best_q = q;
            best_cnt = cnt;
            best_deg = deg;
            if cnt == 1 {
                break;
            }
        }
    }

    let q_nbs = ctx.query.neighbors(NodeId(best_q as u32));
    let is_last = depth + 1 == n1;
    let all_incident_undetermined = ctx.all_edges_undetermined
        || q_nbs
            .iter()
            .all(|nb| !ctx.assigned[nb.node.index()] || ctx.pred.edge_is_undetermined(nb.edge));

    if COUNT_ONLY
        && is_last
        && all_incident_undetermined
        && ctx
            .sym
            .as_ref()
            .is_none_or(DynamicSymmetry::can_leaf_popcnt)
    {
        *ctx.count += ctx.dom_cnt[best_q] as u64;
        return ControlFlow::Continue(());
    }

    // Snapshot candidate bits of `best_q` so propagation can mutate `dom_wide` safely.
    let base = best_q * words;
    let cands: Vec<u64> = ctx.dom_wide[base..base + words].to_vec();

    for (w, &word_bits) in cands.iter().enumerate() {
        let mut bits = word_bits;
        while bits != 0 {
            let tz = bits.trailing_zeros() as usize;
            bits &= bits - 1;
            let v = (w << 6) + tz;

            if ctx
                .sym
                .as_ref()
                .is_some_and(|sym| sym.rejects_before_edge_check(prefix_len, best_q, v))
            {
                continue;
            }

            let v_node = NodeId(v as u32);
            if !all_incident_undetermined {
                let mut edges_ok = true;
                for nb in q_nbs {
                    let u = nb.node.index();
                    if !ctx.assigned[u] || ctx.pred.edge_is_undetermined(nb.edge) {
                        continue;
                    }
                    let te = ctx
                        .target
                        .find_edge(v_node, NodeId(ctx.mapping[u] as u32))
                        .expect("forward checking guarantees adjacency");
                    if !ctx.pred.match_edge(nb.edge, te) {
                        edges_ok = false;
                        break;
                    }
                }
                if !edges_ok {
                    continue;
                }
            }

            ctx.mapping[best_q] = v;
            let (next_prefix_len, saved_depth) = if let Some(sym) = ctx.sym.as_mut() {
                let Some(pair) = sym.try_advance_host_chain(prefix_len, ctx.mapping) else {
                    ctx.mapping[best_q] = usize::MAX;
                    continue;
                };
                pair
            } else {
                (prefix_len, usize::MAX)
            };

            if is_last {
                let accepted = ctx
                    .sym
                    .as_ref()
                    .is_none_or(|sym| sym.accepts_complete_embedding(ctx.mapping));
                let res = if accepted {
                    if COUNT_ONLY {
                        *ctx.count += 1;
                        ControlFlow::Continue(())
                    } else {
                        for (slot, &target_idx) in ctx.out_nodes.iter_mut().zip(ctx.mapping.iter())
                        {
                            *slot = NodeId(target_idx as u32);
                        }
                        (ctx.visitor)(ctx.out_nodes)
                    }
                } else {
                    ControlFlow::Continue(())
                };
                if let Some(sym) = ctx.sym.as_mut() {
                    sym.restore_host_chain(saved_depth);
                }
                ctx.mapping[best_q] = usize::MAX;
                if let ControlFlow::Break(val) = res {
                    return ControlFlow::Break(val);
                }
                continue;
            }

            let word_mark = ctx.undo_words.len();
            let cnt_mark = ctx.undo_cnts.len();
            let unit_mark = ctx.undo_unit_done.len();
            ctx.assigned[best_q] = true;

            let mut wipeout = false;
            for q in 0..n1 {
                if !ctx.assigned[q] && !ctx.clear_bit(q, v) {
                    wipeout = true;
                    break;
                }
            }

            if !wipeout {
                for nb in q_nbs {
                    let u = nb.node.index();
                    if !ctx.assigned[u] && !ctx.intersect_row_with_adj(u, v) {
                        wipeout = true;
                        break;
                    }
                }
            }

            if !wipeout {
                let pq_q_sym = ctx
                    .sym
                    .as_ref()
                    .filter(|s| s.has_q_sym)
                    .map(|s| s.prepared_query);
                if let Some(pq) = pq_q_sym {
                    for &hi in pq.sym_upper_bounds(best_q) {
                        let u = hi as usize;
                        if !ctx.assigned[u] && !ctx.narrow_ge(u, v + 1) {
                            wipeout = true;
                            break;
                        }
                    }
                    if !wipeout {
                        for &lo in pq.sym_lower_bounds(best_q) {
                            let u = lo as usize;
                            if !ctx.assigned[u] && !ctx.narrow_lt(u, v) {
                                wipeout = true;
                                break;
                            }
                        }
                    }
                }
            }

            if !wipeout && next_prefix_len > prefix_len && next_prefix_len < n1 {
                let filter_p_next = ctx
                    .sym
                    .as_ref()
                    .filter(|s| s.has_h_sym && !s.target_chain.is_current_trivial())
                    .map(|s| s.prepared_query.priority_order()[next_prefix_len] as usize);
                if let Some(p_next) = filter_p_next {
                    let p_base = p_next * words;
                    let old_cnt = ctx.dom_cnt[p_next];
                    let mut removed = 0u32;
                    for word_idx in 0..words {
                        let idx = p_base + word_idx;
                        let old = ctx.dom_wide[idx];
                        if old == 0 {
                            continue;
                        }
                        let mut filtered = 0u64;
                        let mut rem_bits = old;
                        while rem_bits != 0 {
                            let b = rem_bits.trailing_zeros() as usize;
                            rem_bits &= rem_bits - 1;
                            let cand_t = (word_idx << 6) + b;
                            let is_min = ctx
                                .sym
                                .as_ref()
                                .expect("sym present")
                                .target_chain
                                .is_orbit_min_at_current(cand_t);
                            if is_min {
                                filtered |= 1u64 << b;
                            }
                        }
                        if filtered != old {
                            ctx.undo_words.push((idx, old));
                            ctx.dom_wide[idx] = filtered;
                            removed += (old ^ filtered).count_ones();
                        }
                    }
                    if removed > 0 {
                        ctx.undo_cnts.push((p_next, old_cnt));
                        ctx.dom_cnt[p_next] = old_cnt - removed;
                        if ctx.dom_cnt[p_next] == 0 {
                            wipeout = true;
                        }
                    }
                }
            }

            if !wipeout {
                wipeout = !propagate_wide(ctx);
            }

            let res = if wipeout {
                ControlFlow::Continue(())
            } else {
                wide_dfs::<COUNT_ONLY, B, P>(ctx, depth + 1, next_prefix_len)
            };

            ctx.assigned[best_q] = false;
            ctx.restore_marks(word_mark, cnt_mark, unit_mark);
            if let Some(sym) = ctx.sym.as_mut() {
                sym.restore_host_chain(saved_depth);
            }
            ctx.mapping[best_q] = usize::MAX;
            if let ControlFlow::Break(val) = res {
                return ControlFlow::Break(val);
            }
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
    use crate::id::EdgeId;
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
        let pq = PreparedSubgraphQuery::new(&query, SubgraphIsomorphismAlgorithm::Glasgow);
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
    #[case::u64_hall_wipeout(10, 5, vec![0, 0, 0, 1, 2], vec![0, 0, 1, 2, 1, 2, 1, 2, 1, 2], 0)]
    #[case::u64_hall_forcing_cascade(
        6,
        4,
        vec![0, 0, 1, 2],
        vec![0, 0, 1, 2, 1, 2],
        2
    )]
    #[case::u128_hall_wipeout(
        72,
        5,
        vec![0, 0, 0, 1, 2],
        (0..72).map(|t| if t < 2 { 0 } else if t % 2 == 0 { 1 } else { 2 }).collect(),
        0
    )]
    #[case::wide_hall_wipeout(
        144,
        5,
        vec![0, 0, 0, 1, 2],
        (0..144).map(|t| if t < 2 { 0 } else if t % 2 == 0 { 1 } else { 2 }).collect(),
        0
    )]
    #[case::wide_hall_forcing_cascade(
        144,
        4,
        vec![0, 0, 1, 2],
        (0..144).map(|t| if t < 2 { 0 } else if t == 2 { 1 } else if t == 3 { 2 } else { 3 }).collect(),
        2
    )]
    fn test_visit_prepared_hall_and_unit_propagation(
        #[case] target_nodes: usize,
        #[case] query_nodes: usize,
        #[case] q_colors: Vec<usize>,
        #[case] t_colors: Vec<usize>,
        #[case] expected_count: u64,
    ) {
        // Connect first 4 target nodes as a path 0-1-2-3 so forcing cascade has a unique path
        // embedding up to swapping {0, 1}, and Hall wipeout cases have 3 query vertices sharing
        // only 2 color-0 target vertices.
        let target = Graph::new(target_nodes, &[[0, 1], [1, 2], [2, 3], [0, 2]]);
        let query = if query_nodes == 4 {
            Graph::new(4, &[[0, 1], [1, 2], [2, 3]])
        } else {
            Graph::new(query_nodes, &[])
        };
        let pt = PreparedSubgraphTarget::new(&target);
        let pq = PreparedSubgraphQuery::new(&query, SubgraphIsomorphismAlgorithm::Glasgow);
        let mut scratch = SubgraphIsomorphismScratch::new();
        let mut pred = ClosureSubgraphIsomorphismPredicate::new(
            |q: NodeId, t: NodeId| q_colors[q.index()] == t_colors[t.index()],
            |_, _| true,
        );
        let mut out = Vec::new();
        let _: ControlFlow<()> =
            visit_prepared(&pt, &pq, None, &mut pred, &mut scratch, &mut |emb| {
                out.push(emb.to_vec());
                ControlFlow::Continue(())
            });
        assert_eq!(out.len() as u64, expected_count);
        assert_eq!(
            count_prepared(&pt, &pq, None, &mut pred, &mut scratch),
            expected_count
        );
    }

    #[rstest]
    #[case::u64_unsat(8, vec![vec![0, 1], vec![0, 1], vec![0, 1, 2], vec![0, 1, 2]], 0)]
    #[case::u64_sat(8, vec![vec![0, 1], vec![0, 1], vec![0, 1, 2], vec![0, 2, 3]], 2)]
    #[case::u128_unsat(70, vec![vec![62, 65], vec![62, 65], vec![62, 65, 68], vec![62, 65, 68]], 0)]
    #[case::u128_sat(70, vec![vec![62, 65], vec![62, 65], vec![62, 65, 68], vec![62, 68, 69]], 2)]
    #[case::wide_unsat(144, vec![vec![10, 130], vec![10, 130], vec![10, 130, 140], vec![10, 130, 140]], 0)]
    #[case::wide_sat(144, vec![vec![10, 130], vec![10, 130], vec![10, 130, 140], vec![10, 140, 141]], 2)]
    fn test_visit_prepared_hall_to_unit_cascade(
        #[case] target_nodes: usize,
        #[case] allowed: Vec<Vec<usize>>,
        #[case] expected_count: u64,
    ) {
        let target = Graph::new(target_nodes, &[]);
        let query = Graph::new(allowed.len(), &[]);
        let pt = PreparedSubgraphTarget::new(&target);
        let pq = PreparedSubgraphQuery::new(&query, SubgraphIsomorphismAlgorithm::Glasgow);
        let mut scratch = SubgraphIsomorphismScratch::new();
        let mut pred = ClosureSubgraphIsomorphismPredicate::new(
            |q: NodeId, t: NodeId| allowed[q.index()].contains(&t.index()),
            |_, _| true,
        );
        let mut out = Vec::new();
        let _: ControlFlow<()> =
            visit_prepared(&pt, &pq, None, &mut pred, &mut scratch, &mut |emb| {
                out.push(emb.to_vec());
                ControlFlow::Continue(())
            });
        assert_eq!(out.len() as u64, expected_count);
        assert_eq!(
            count_prepared(&pt, &pq, None, &mut pred, &mut scratch),
            expected_count
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
        let pq = PreparedSubgraphQuery::new(&query, SubgraphIsomorphismAlgorithm::Glasgow);
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
        let pq_base = PreparedSubgraphQuery::new(&query, SubgraphIsomorphismAlgorithm::Glasgow);
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
        let pq_base = PreparedSubgraphQuery::new(&query, SubgraphIsomorphismAlgorithm::Glasgow);
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
            .map(|g| PreparedSubgraphQuery::new(g, SubgraphIsomorphismAlgorithm::Glasgow))
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
            let mut pq = PreparedSubgraphQuery::new(&query, SubgraphIsomorphismAlgorithm::Glasgow);
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
            PreparedSubgraphQuery::new(&Graph::new(0, &[]), SubgraphIsomorphismAlgorithm::Glasgow);
        let wrong_sym = TargetSymmetry::none(4);
        let mut scratch = SubgraphIsomorphismScratch::new();
        let mut pred = ClosureSubgraphIsomorphismPredicate::new(|_, _| true, |_, _| true);
        let _ = count_prepared(&pt, &pq, Some(&wrong_sym), &mut pred, &mut scratch);
    }
}
