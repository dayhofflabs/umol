//! Upgraded Ray--Kirsch bond-oriented backtracking engine (`RayKirsch`).
//!
//! Compiles the query into a flat RI-guided [`RkStep`] (`Seed`/`Extend`/`Close`) program anchored
//! to `(Reverse(degree(q)), q)` (or the caller-supplied anchor vertex), checks `deg(q) <= deg(t)`
//! and degree-sequence domination, verifies extension bonds directly from the CSR neighbor slice,
//! and pushes down query (`G_Q`) and target (`G_H`) symmetry constraints.
//!
//! See [Ray and Kirsch (1957)](https://doi.org/10.1126/science.126.3278.814) and Open Babel's
//! [`parsmart.cpp`](https://github.com/openbabel/openbabel/blob/master/src/parsmart.cpp).

use std::ops::ControlFlow;

use super::predicate::SubgraphIsomorphismPredicate;
use super::prepared::{
    build_ray_kirsch_plan, PreparedQueryPlan, PreparedSubgraphQuery, PreparedSubgraphTarget,
    RkStep, SubgraphIsomorphismScratch,
};
use super::symmetry::{
    is_canonical_double_orbit_with_built_chain, QuerySymmetryConstraints, TargetStabilizerChain,
    TargetSymmetry,
};
use crate::graph::Graph;
use crate::id::NodeId;

/// Step-indexed push-down symmetry state for `RayKirsch`.
struct StepSymmetry<'a> {
    has_q_sym: bool,
    has_h_sym: bool,
    sym_gt_offsets: Vec<usize>,
    sym_gt_nodes: Vec<usize>,
    sym_lt_offsets: Vec<usize>,
    sym_lt_nodes: Vec<usize>,
    priority_order: &'a [u32],
    prefix_len_after: Vec<usize>,
    query_sym: Option<&'a QuerySymmetryConstraints>,
    target_sym: Option<&'a TargetSymmetry>,
    target_chain: &'a mut TargetStabilizerChain,
}

impl<'a> StepSymmetry<'a> {
    fn build(
        prepared_query: &'a PreparedSubgraphQuery,
        target_sym: Option<&'a TargetSymmetry>,
        steps: &[RkStep],
        target_chain: &'a mut TargetStabilizerChain,
    ) -> Option<Self> {
        let query_sym =
            (!prepared_query.symmetry().is_trivial()).then_some(prepared_query.symmetry());
        let active_target_sym = target_sym.filter(|sym| !sym.is_trivial());
        let has_q_sym = query_sym.is_some();
        let has_h_sym = active_target_sym.is_some();
        if !has_q_sym && !has_h_sym {
            return None;
        }

        let n1 = prepared_query.node_count();
        let step_count = steps.len();
        let mut node_step = vec![usize::MAX; n1];
        for (s, step) in steps.iter().enumerate() {
            match *step {
                RkStep::Seed { node, .. } | RkStep::Extend { to: node, .. } => {
                    node_step[node] = s;
                }
                RkStep::Close { .. } => {}
            }
        }

        let mut sym_gt_offsets = Vec::new();
        let mut sym_gt_nodes = Vec::new();
        let mut sym_lt_offsets = Vec::new();
        let mut sym_lt_nodes = Vec::new();
        if has_q_sym {
            sym_gt_offsets.reserve(step_count + 1);
            sym_lt_offsets.reserve(step_count + 1);
            sym_gt_offsets.push(0);
            sym_lt_offsets.push(0);
            for (s, step) in steps.iter().enumerate() {
                if let RkStep::Seed { node: q, .. } | RkStep::Extend { to: q, .. } = *step {
                    for &lo in prepared_query.sym_lower_bounds(q) {
                        let lo_idx = lo as usize;
                        if node_step[lo_idx] < s {
                            sym_gt_nodes.push(lo_idx);
                        }
                    }
                    for &hi in prepared_query.sym_upper_bounds(q) {
                        let hi_idx = hi as usize;
                        if node_step[hi_idx] < s {
                            sym_lt_nodes.push(hi_idx);
                        }
                    }
                }
                sym_gt_offsets.push(sym_gt_nodes.len());
                sym_lt_offsets.push(sym_lt_nodes.len());
            }
        }

        let priority_order = prepared_query.priority_order();
        let mut prefix_len_after = Vec::new();
        if let Some(ts) = active_target_sym {
            target_chain.reset(ts);
            prefix_len_after.reserve(step_count);
            let mut k = 0;
            for s in 0..step_count {
                while k < n1 && node_step[priority_order[k] as usize] <= s {
                    k += 1;
                }
                prefix_len_after.push(k);
            }
        }

        Some(Self {
            has_q_sym,
            has_h_sym,
            sym_gt_offsets,
            sym_gt_nodes,
            sym_lt_offsets,
            sym_lt_nodes,
            priority_order,
            prefix_len_after,
            query_sym,
            target_sym: active_target_sym,
            target_chain,
        })
    }

    #[inline(always)]
    fn step_bounds(&self, idx: usize, mapping: &[usize], n2: usize) -> (usize, usize) {
        if !self.has_q_sym {
            return (0, n2);
        }
        let mut min_v = 0usize;
        for &lo in &self.sym_gt_nodes[self.sym_gt_offsets[idx]..self.sym_gt_offsets[idx + 1]] {
            min_v = min_v.max(mapping[lo] + 1);
        }
        let mut max_v = n2;
        for &hi in &self.sym_lt_nodes[self.sym_lt_offsets[idx]..self.sym_lt_offsets[idx + 1]] {
            max_v = max_v.min(mapping[hi]);
        }
        (min_v, max_v)
    }

    #[inline(always)]
    fn rejects_before_predicate(&self, idx: usize, v: usize) -> bool {
        let Some(target_sym) = self.target_sym else {
            return false;
        };
        if self.target_chain.is_current_trivial() {
            return false;
        }
        let prev_k = if idx == 0 {
            0
        } else {
            self.prefix_len_after[idx - 1]
        };
        if self.prefix_len_after[idx] > prev_k {
            !self.target_chain.is_orbit_min_at_current(target_sym, v)
        } else {
            false
        }
    }

    #[inline]
    fn try_advance_host_chain(&mut self, idx: usize, mapping: &[usize]) -> Option<usize> {
        let Some(target_sym) = self.target_sym else {
            return Some(usize::MAX);
        };
        if self.target_chain.is_current_trivial() {
            return Some(usize::MAX);
        }
        let prev_k = if idx == 0 {
            0
        } else {
            self.prefix_len_after[idx - 1]
        };
        let new_k = self.prefix_len_after[idx];
        if new_k == prev_k {
            return Some(usize::MAX);
        }

        let n1 = self.priority_order.len();
        let saved_depth = self.target_chain.depth();
        let mut pushed = false;
        for r in prev_k..new_k {
            if self.target_chain.is_current_trivial() {
                break;
            }
            let t_r = mapping[self.priority_order[r] as usize];
            if !self.has_q_sym && r + 1 == n1 {
                if !self.target_chain.is_orbit_min_at_current(target_sym, t_r) {
                    if pushed {
                        self.target_chain.truncate(saved_depth);
                    }
                    return None;
                }
            } else if self.target_chain.try_push_orbit_min_with_tree(
                target_sym,
                t_r,
                self.has_q_sym,
            ) {
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
    fn accepts_complete_embedding(&self, mapping: &[usize]) -> bool {
        if !(self.has_q_sym && self.has_h_sym) {
            return true;
        }
        is_canonical_double_orbit_with_built_chain(
            mapping,
            self.priority_order,
            self.query_sym
                .expect("query symmetry present when has_q_sym is true"),
            self.target_sym
                .expect("target symmetry present when has_h_sym is true"),
            self.target_chain,
        )
    }
}

/// Visit every subgraph isomorphism of `query` in `target` using upgraded `RayKirsch`.
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

/// Count every subgraph isomorphism of `query` in `target` using upgraded `RayKirsch`.
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

/// Visit every subgraph isomorphism of a prepared query in a prepared target using `RayKirsch`.
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

/// Count every subgraph isomorphism of a prepared query in a prepared target using `RayKirsch`.
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

    let steps = build_ray_kirsch_plan(query, anchor.map(|(qa, _)| qa.index()));
    let mut mapping = vec![usize::MAX; n1];
    let mut used = vec![false; n2];
    let mut out_nodes = vec![NodeId(0); n1];
    let mut ctx = RayKirschCtx {
        target,
        steps: &steps,
        forced_root: anchor.map(|(_, ta)| ta.index()),
        all_edges_undetermined: pred.all_edges_undetermined(),
        mapping: &mut mapping,
        used: &mut used,
        out_nodes: &mut out_nodes,
        sym: None,
        pred,
        count,
        visitor,
    };
    dfs::<COUNT_ONLY, B, P>(&mut ctx, 0)
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

    let PreparedQueryPlan::RayKirsch(steps) = prepared_query.plan() else {
        unreachable!("RayKirsch engine invoked with non-RayKirsch prepared query");
    };

    scratch.mapping.clear();
    scratch.mapping.resize(n1, usize::MAX);
    scratch.used.clear();
    scratch.used.resize(n2, false);
    scratch.out_nodes.resize(n1, NodeId(0));

    let sym = StepSymmetry::build(prepared_query, target_sym, steps, &mut scratch.target_chain);
    let mut ctx = RayKirschCtx {
        target,
        steps,
        forced_root: None,
        all_edges_undetermined: pred.all_edges_undetermined(),
        mapping: &mut scratch.mapping,
        used: &mut scratch.used,
        out_nodes: &mut scratch.out_nodes,
        sym,
        pred,
        count,
        visitor,
    };
    dfs::<COUNT_ONLY, B, P>(&mut ctx, 0)
}

struct RayKirschCtx<'a, P, F> {
    target: &'a Graph,
    steps: &'a [RkStep],
    forced_root: Option<usize>,
    all_edges_undetermined: bool,
    mapping: &'a mut [usize],
    used: &'a mut [bool],
    out_nodes: &'a mut [NodeId],
    sym: Option<StepSymmetry<'a>>,
    pred: &'a mut P,
    count: &'a mut u64,
    visitor: &'a mut F,
}

#[inline(always)]
fn try_descend<const COUNT_ONLY: bool, B, P, F>(
    ctx: &mut RayKirschCtx<'_, P, F>,
    idx: usize,
    q: usize,
    v: usize,
) -> ControlFlow<B>
where
    P: SubgraphIsomorphismPredicate,
    F: FnMut(&[NodeId]) -> ControlFlow<B>,
{
    ctx.mapping[q] = v;
    let saved_depth = if let Some(sym) = ctx.sym.as_mut() {
        let Some(saved) = sym.try_advance_host_chain(idx, ctx.mapping) else {
            ctx.mapping[q] = usize::MAX;
            return ControlFlow::Continue(());
        };
        saved
    } else {
        usize::MAX
    };

    ctx.used[v] = true;
    let res = dfs::<COUNT_ONLY, B, P>(ctx, idx + 1);
    ctx.used[v] = false;
    if let Some(sym) = ctx.sym.as_mut() {
        sym.restore_host_chain(saved_depth);
    }
    ctx.mapping[q] = usize::MAX;
    res
}

fn dfs<const COUNT_ONLY: bool, B, P: SubgraphIsomorphismPredicate>(
    ctx: &mut RayKirschCtx<'_, P, impl FnMut(&[NodeId]) -> ControlFlow<B>>,
    idx: usize,
) -> ControlFlow<B> {
    if idx == ctx.steps.len() {
        let accepted = ctx
            .sym
            .as_ref()
            .is_none_or(|sym| sym.accepts_complete_embedding(ctx.mapping));
        if !accepted {
            return ControlFlow::Continue(());
        }
        if COUNT_ONLY {
            *ctx.count += 1;
            return ControlFlow::Continue(());
        }
        for (slot, &target_idx) in ctx.out_nodes.iter_mut().zip(ctx.mapping.iter()) {
            *slot = NodeId(target_idx as u32);
        }
        return (ctx.visitor)(ctx.out_nodes);
    }

    match ctx.steps[idx] {
        RkStep::Seed { node: q, degree } => {
            if idx == 0 {
                if let Some(ta) = ctx.forced_root {
                    return try_descend::<COUNT_ONLY, B, P, _>(ctx, 0, q, ta);
                }
            }
            let n2 = ctx.target.node_count();
            let (min_v, max_v) = match ctx.sym.as_ref() {
                Some(sym) if sym.has_q_sym => sym.step_bounds(idx, ctx.mapping, n2),
                _ => (0, n2),
            };
            if min_v >= max_v {
                return ControlFlow::Continue(());
            }
            let q_node = NodeId(q as u32);
            for v in min_v..max_v {
                if ctx.used[v] {
                    continue;
                }
                let v_node = NodeId(v as u32);
                if degree > ctx.target.degree(v_node) {
                    continue;
                }
                if ctx
                    .sym
                    .as_ref()
                    .is_some_and(|sym| sym.rejects_before_predicate(idx, v))
                {
                    continue;
                }
                if !ctx.pred.match_node(q_node, v_node) {
                    continue;
                }
                let res = try_descend::<COUNT_ONLY, B, P, _>(ctx, idx, q, v);
                if let ControlFlow::Break(val) = res {
                    return ControlFlow::Break(val);
                }
            }
            ControlFlow::Continue(())
        }
        RkStep::Extend {
            from,
            to: q,
            degree_to,
            edge,
        } => {
            let n2 = ctx.target.node_count();
            let (min_v, max_v) = match ctx.sym.as_ref() {
                Some(sym) if sym.has_q_sym => sym.step_bounds(idx, ctx.mapping, n2),
                _ => (0, n2),
            };
            if min_v >= max_v {
                return ControlFlow::Continue(());
            }
            let from_node = NodeId(ctx.mapping[from] as u32);
            let q_node = NodeId(q as u32);
            let edge_undetermined =
                ctx.all_edges_undetermined || ctx.pred.edge_is_undetermined(edge);
            for nb in ctx.target.neighbors(from_node) {
                let v = nb.node.index();
                if v < min_v || v >= max_v || ctx.used[v] || degree_to > ctx.target.degree(nb.node)
                {
                    continue;
                }
                if ctx
                    .sym
                    .as_ref()
                    .is_some_and(|sym| sym.rejects_before_predicate(idx, v))
                {
                    continue;
                }
                if !ctx.pred.match_node(q_node, nb.node) {
                    continue;
                }
                if !edge_undetermined && !ctx.pred.match_edge(edge, nb.edge) {
                    continue;
                }
                let res = try_descend::<COUNT_ONLY, B, P, _>(ctx, idx, q, v);
                if let ControlFlow::Break(val) = res {
                    return ControlFlow::Break(val);
                }
            }
            ControlFlow::Continue(())
        }
        RkStep::Close {
            first,
            second,
            edge,
        } => {
            let Some(te) = ctx.target.find_edge(
                NodeId(ctx.mapping[first] as u32),
                NodeId(ctx.mapping[second] as u32),
            ) else {
                return ControlFlow::Continue(());
            };
            if !ctx.all_edges_undetermined
                && !ctx.pred.edge_is_undetermined(edge)
                && !ctx.pred.match_edge(edge, te)
            {
                return ControlFlow::Continue(());
            }
            dfs::<COUNT_ONLY, B, P>(ctx, idx + 1)
        }
    }
}
