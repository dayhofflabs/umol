//! Cleaned RDKit-compatible VF2 engine (`Vf2Rdkit`).
//!
//! Uses a precompiled static BFS query order with flat CSR parent-edge lists, checks the primary
//! parent edge directly from the target CSR neighbor slice (avoiding a `find_edge` lookup for the
//! first parent), enforces `deg(q) <= deg(t)` and degree-sequence domination, and pushes down
//! query (`G_Q`) and target (`G_H`) symmetry constraints.

use std::ops::ControlFlow;

use super::predicate::SubgraphIsomorphismPredicate;
use super::prepared::{
    build_vf2rdkit_plan, PreparedQueryPlan, PreparedSubgraphQuery, PreparedSubgraphTarget,
    StaticOrderPlan, SubgraphIsomorphismScratch,
};
use super::symmetry::{
    is_canonical_double_orbit_with_built_chain, QuerySymmetryConstraints, TargetStabilizerChain,
    TargetSymmetry,
};
use crate::graph::Graph;
use crate::id::NodeId;

/// Active static-order push-down symmetry state for `Vf2Rdkit`.
struct StaticSymmetry<'a> {
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

impl<'a> StaticSymmetry<'a> {
    fn build(
        prepared_query: &'a PreparedSubgraphQuery,
        target_sym: Option<&'a TargetSymmetry>,
        order: &[usize],
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

        let n1 = order.len();
        let mut pos = vec![0usize; n1];
        for (d, &q) in order.iter().enumerate() {
            pos[q] = d;
        }

        let mut sym_gt_offsets = Vec::new();
        let mut sym_gt_depths = Vec::new();
        let mut sym_lt_offsets = Vec::new();
        let mut sym_lt_depths = Vec::new();
        if has_q_sym {
            sym_gt_offsets.reserve(n1 + 1);
            sym_lt_offsets.reserve(n1 + 1);
            sym_gt_offsets.push(0);
            sym_lt_offsets.push(0);
            for (d, &q) in order.iter().enumerate() {
                for &lo in prepared_query.sym_lower_bounds(q) {
                    let p = pos[lo as usize];
                    if p < d {
                        sym_gt_depths.push(p);
                    }
                }
                for &hi in prepared_query.sym_upper_bounds(q) {
                    let p = pos[hi as usize];
                    if p < d {
                        sym_lt_depths.push(p);
                    }
                }
                sym_gt_offsets.push(sym_gt_depths.len());
                sym_lt_offsets.push(sym_lt_depths.len());
            }
        }

        let priority_order = prepared_query.priority_order();
        let mut priority_pos = Vec::new();
        let mut prefix_len_after = Vec::new();
        if let Some(ts) = active_target_sym {
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
    fn rejects_before_predicate(&self, depth: usize, v: usize) -> bool {
        let Some(target_sym) = self.target_sym else {
            return false;
        };
        if self.target_chain.is_current_trivial() {
            return false;
        }
        let prev_k = if depth == 0 {
            0
        } else {
            self.prefix_len_after[depth - 1]
        };
        if self.prefix_len_after[depth] > prev_k {
            !self.target_chain.is_orbit_min_at_current(target_sym, v)
        } else {
            false
        }
    }

    #[inline]
    fn try_advance_host_chain(&mut self, depth: usize, mapping: &[usize]) -> Option<usize> {
        let Some(target_sym) = self.target_sym else {
            return Some(usize::MAX);
        };
        if self.target_chain.is_current_trivial() {
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

/// Visit every subgraph isomorphism of `query` in `target` using cleaned `Vf2Rdkit`.
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

/// Count every subgraph isomorphism of `query` in `target` using cleaned `Vf2Rdkit`.
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

/// Visit every subgraph isomorphism of a prepared query in a prepared target using `Vf2Rdkit`.
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

/// Count every subgraph isomorphism of a prepared query in a prepared target using `Vf2Rdkit`.
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

    let plan = build_vf2rdkit_plan(query, anchor.map(|(qa, _)| qa.index()));
    let mut mapping = vec![0usize; n1];
    let mut used = vec![false; n2];
    let mut out_nodes = vec![NodeId(0); n1];
    let mut ctx = Vf2RdkitCtx {
        target,
        query,
        plan: &plan,
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

    let PreparedQueryPlan::Vf2Rdkit(plan) = prepared_query.plan() else {
        unreachable!("Vf2Rdkit engine invoked with non-Vf2Rdkit prepared query");
    };

    scratch.mapping.clear();
    scratch.mapping.resize(n1, 0);
    scratch.used.clear();
    scratch.used.resize(n2, false);
    scratch.out_nodes.resize(n1, NodeId(0));

    let sym = StaticSymmetry::build(
        prepared_query,
        target_sym,
        &plan.order,
        &mut scratch.target_chain,
    );
    let mut ctx = Vf2RdkitCtx {
        target,
        query,
        plan,
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

struct Vf2RdkitCtx<'a, P, F> {
    target: &'a Graph,
    query: &'a Graph,
    plan: &'a StaticOrderPlan,
    forced_root: Option<usize>,
    all_edges_undetermined: bool,
    mapping: &'a mut [usize],
    used: &'a mut [bool],
    out_nodes: &'a mut [NodeId],
    sym: Option<StaticSymmetry<'a>>,
    pred: &'a mut P,
    count: &'a mut u64,
    visitor: &'a mut F,
}

#[inline(always)]
fn try_descend<const COUNT_ONLY: bool, B, P, F>(
    ctx: &mut Vf2RdkitCtx<'_, P, F>,
    depth: usize,
    v: usize,
) -> ControlFlow<B>
where
    P: SubgraphIsomorphismPredicate,
    F: FnMut(&[NodeId]) -> ControlFlow<B>,
{
    let is_last = depth + 1 == ctx.plan.order.len();
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
            .is_none_or(|sym| sym.accepts_complete_embedding(&ctx.plan.order, ctx.mapping));
        if accepted {
            if COUNT_ONLY {
                *ctx.count += 1;
                ControlFlow::Continue(())
            } else {
                for (d, &q) in ctx.plan.order.iter().enumerate() {
                    ctx.out_nodes[q] = NodeId(ctx.mapping[d] as u32);
                }
                (ctx.visitor)(ctx.out_nodes)
            }
        } else {
            ControlFlow::Continue(())
        }
    } else {
        ctx.used[v] = true;
        let inner = dfs::<COUNT_ONLY, B, P>(ctx, depth + 1);
        ctx.used[v] = false;
        inner
    };

    if let Some(sym) = ctx.sym.as_mut() {
        sym.restore_host_chain(saved_depth);
    }
    res
}

fn dfs<const COUNT_ONLY: bool, B, P: SubgraphIsomorphismPredicate>(
    ctx: &mut Vf2RdkitCtx<'_, P, impl FnMut(&[NodeId]) -> ControlFlow<B>>,
    depth: usize,
) -> ControlFlow<B> {
    if depth == 0 {
        if let Some(ta) = ctx.forced_root {
            return try_descend::<COUNT_ONLY, B, P, _>(ctx, 0, ta);
        }
    }

    let n2 = ctx.target.node_count();
    let (min_v, max_v) = match ctx.sym.as_ref() {
        Some(sym) if sym.has_q_sym => sym.depth_bounds(depth, ctx.mapping, n2),
        _ => (0, n2),
    };
    if min_v >= max_v {
        return ControlFlow::Continue(());
    }

    let q = ctx.plan.order[depth];
    let q_node = NodeId(q as u32);
    let q_deg = ctx.query.degree(q_node);
    let p_start = ctx.plan.parents_offsets[depth];
    let p_end = ctx.plan.parents_offsets[depth + 1];

    if p_start == p_end {
        for v in min_v..max_v {
            if ctx.used[v] {
                continue;
            }
            let v_node = NodeId(v as u32);
            if q_deg > ctx.target.degree(v_node) {
                continue;
            }
            if ctx
                .sym
                .as_ref()
                .is_some_and(|sym| sym.rejects_before_predicate(depth, v))
            {
                continue;
            }
            if !ctx.pred.match_node(q_node, v_node) {
                continue;
            }
            let res = try_descend::<COUNT_ONLY, B, P, _>(ctx, depth, v);
            if let ControlFlow::Break(val) = res {
                return ControlFlow::Break(val);
            }
        }
        return ControlFlow::Continue(());
    }

    let (p0, pe0) = ctx.plan.parents_data[p_start];
    let mp0 = ctx.mapping[p0];
    let pe0_undetermined = ctx.all_edges_undetermined || ctx.pred.edge_is_undetermined(pe0);

    for nb in ctx.target.neighbors(NodeId(mp0 as u32)) {
        let v = nb.node.index();
        if v < min_v || v >= max_v || ctx.used[v] || q_deg > ctx.target.degree(nb.node) {
            continue;
        }
        if ctx
            .sym
            .as_ref()
            .is_some_and(|sym| sym.rejects_before_predicate(depth, v))
        {
            continue;
        }
        if !ctx.pred.match_node(q_node, nb.node) {
            continue;
        }
        if !pe0_undetermined && !ctx.pred.match_edge(pe0, nb.edge) {
            continue;
        }
        let mut edges_ok = true;
        for &(pk, pek) in &ctx.plan.parents_data[(p_start + 1)..p_end] {
            let Some(te) = ctx
                .target
                .find_edge(NodeId(ctx.mapping[pk] as u32), nb.node)
            else {
                edges_ok = false;
                break;
            };
            if !ctx.all_edges_undetermined
                && !ctx.pred.edge_is_undetermined(pek)
                && !ctx.pred.match_edge(pek, te)
            {
                edges_ok = false;
                break;
            }
        }
        if !edges_ok {
            continue;
        }
        let res = try_descend::<COUNT_ONLY, B, P, _>(ctx, depth, v);
        if let ControlFlow::Break(val) = res {
            return ControlFlow::Break(val);
        }
    }

    ControlFlow::Continue(())
}
