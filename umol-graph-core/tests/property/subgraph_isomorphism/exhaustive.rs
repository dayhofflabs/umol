//! Independent definition-level subgraph-isomorphism and orbit calculations for
//! bounded test graphs.
//!
//! All calculations in this module check the mathematical definitions directly
//! without calling any production subgraph-isomorphism engine, automorphism
//! solver, or Schreier--Sims stabilizer code.

use std::cmp::Reverse;
use std::collections::{HashSet, VecDeque};

use umol_graph_core::{EdgeId, Graph, NodeId};

/// Deterministic query vertex priority `(Reverse(degree(q)), q)` matching the
/// documented canonical orbit ordering.
pub(super) fn query_priority_order(query: &Graph) -> Vec<usize> {
    let mut order: Vec<usize> = (0..query.node_count()).collect();
    order.sort_unstable_by_key(|&q| (Reverse(query.degree(NodeId(q as u32))), q));
    order
}

/// Checks directly whether `embedding` is a valid labeled monomorphism from
/// `query` into `target`.
pub(super) fn is_valid_monomorphism(
    query: &Graph,
    target: &Graph,
    embedding: &[usize],
    node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
    edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
) -> bool {
    let n1 = query.node_count();
    let n2 = target.node_count();
    if embedding.len() != n1 {
        return false;
    }
    if embedding.iter().any(|&t| t >= n2) {
        return false;
    }

    let mut sorted = embedding.to_vec();
    sorted.sort_unstable();
    if sorted.windows(2).any(|pair| pair[0] == pair[1]) {
        return false;
    }

    for (q, &t) in embedding.iter().enumerate() {
        if !node_match(NodeId(q as u32), NodeId(t as u32)) {
            return false;
        }
    }

    for qe in query.edge_ids() {
        let [u, v] = query.edge_endpoints(qe);
        let tu = NodeId(embedding[u.index()] as u32);
        let tv = NodeId(embedding[v.index()] as u32);
        match target.find_edge(tu, tv) {
            Some(te) if edge_match(qe, te) => {}
            _ => return false,
        }
    }

    true
}

/// Exhaustively enumerates every labeled monomorphism of `query` into `target`
/// by scanning all injective mappings on `0..query.node_count()` and checking
/// every query node and edge against the definition.
pub(super) fn exhaustive_monomorphisms(
    query: &Graph,
    target: &Graph,
    node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
    edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
) -> Vec<Vec<usize>> {
    let n1 = query.node_count();
    let n2 = target.node_count();
    if n1 == 0 {
        return vec![Vec::new()];
    }
    if n1 > n2 {
        return Vec::new();
    }

    let mut compat = vec![vec![false; n2]; n1];
    for (q, row) in compat.iter_mut().enumerate() {
        for (t, cell) in row.iter_mut().enumerate() {
            *cell = node_match(NodeId(q as u32), NodeId(t as u32));
        }
    }

    let query_edges: Vec<(EdgeId, usize, usize)> = query
        .edge_ids()
        .map(|qe| {
            let [u, v] = query.edge_endpoints(qe);
            (qe, u.index(), v.index())
        })
        .collect();

    let mut mapping = vec![0usize; n1];
    let mut used = vec![false; n2];
    let mut results = Vec::new();

    struct MonoContext<'a> {
        compat: &'a [Vec<bool>],
        query_edges: &'a [(EdgeId, usize, usize)],
        target: &'a Graph,
    }

    fn backtrack(
        q: usize,
        ctx: &MonoContext<'_>,
        edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
        mapping: &mut [usize],
        used: &mut [bool],
        results: &mut Vec<Vec<usize>>,
    ) {
        if q == ctx.compat.len() {
            for &(qe, u, v) in ctx.query_edges {
                let tu = NodeId(mapping[u] as u32);
                let tv = NodeId(mapping[v] as u32);
                match ctx.target.find_edge(tu, tv) {
                    Some(te) if edge_match(qe, te) => {}
                    _ => return,
                }
            }
            results.push(mapping.to_vec());
            return;
        }
        for t in 0..used.len() {
            if !used[t] && ctx.compat[q][t] {
                used[t] = true;
                mapping[q] = t;
                backtrack(q + 1, ctx, edge_match, mapping, used, results);
                used[t] = false;
            }
        }
    }

    let ctx = MonoContext {
        compat: &compat,
        query_edges: &query_edges,
        target,
    };
    backtrack(0, &ctx, edge_match, &mut mapping, &mut used, &mut results);
    results.sort();
    results
}

/// Exhaustively enumerates the full colored automorphism group of a small
/// `graph` by testing bijections directly against adjacency and color equality.
pub(super) fn exhaustive_automorphisms<NC: Eq, EC: Eq>(
    graph: &Graph,
    node_color: impl Fn(NodeId) -> NC,
    edge_color: impl Fn(EdgeId) -> EC,
) -> Vec<Vec<usize>> {
    let n = graph.node_count();
    if n == 0 {
        return vec![Vec::new()];
    }

    let node_colors: Vec<NC> = graph.node_ids().map(&node_color).collect();
    let degrees: Vec<usize> = graph.node_ids().map(|v| graph.degree(v)).collect();
    let edge_colors: Vec<EC> = graph.edge_ids().map(&edge_color).collect();

    let mut perm = vec![0usize; n];
    let mut used = vec![false; n];
    let mut group = Vec::new();

    struct AutoContext<'a, NC, EC> {
        graph: &'a Graph,
        node_colors: &'a [NC],
        degrees: &'a [usize],
        edge_colors: &'a [EC],
    }

    fn backtrack<NC: Eq, EC: Eq>(
        v: usize,
        ctx: &AutoContext<'_, NC, EC>,
        perm: &mut [usize],
        used: &mut [bool],
        group: &mut Vec<Vec<usize>>,
    ) {
        if v == used.len() {
            group.push(perm.to_vec());
            return;
        }
        for img in 0..used.len() {
            if used[img]
                || ctx.degrees[v] != ctx.degrees[img]
                || ctx.node_colors[v] != ctx.node_colors[img]
            {
                continue;
            }
            let mut ok = true;
            for (u, &img_u) in perm.iter().enumerate().take(v) {
                let src_edge = ctx.graph.find_edge(NodeId(u as u32), NodeId(v as u32));
                let dst_edge = ctx
                    .graph
                    .find_edge(NodeId(img_u as u32), NodeId(img as u32));
                match (src_edge, dst_edge) {
                    (None, None) => {}
                    (Some(se), Some(de))
                        if ctx.edge_colors[se.index()] == ctx.edge_colors[de.index()] => {}
                    _ => {
                        ok = false;
                        break;
                    }
                }
            }
            if ok {
                used[img] = true;
                perm[v] = img;
                backtrack(v + 1, ctx, perm, used, group);
                used[img] = false;
            }
        }
    }

    let ctx = AutoContext {
        graph,
        node_colors: &node_colors,
        degrees: &degrees,
        edge_colors: &edge_colors,
    };
    backtrack(0, &ctx, &mut perm, &mut used, &mut group);
    group.sort();
    group
}

/// Computes the exact finite subgroup generated by `generators` via breadth-first
/// closure under permutation composition.
pub(super) fn subgroup_from_generators(
    node_count: usize,
    generators: &[Vec<NodeId>],
) -> Vec<Vec<usize>> {
    let gens: Vec<Vec<usize>> = generators
        .iter()
        .map(|gen| gen.iter().map(|node| node.index()).collect())
        .collect();
    let identity: Vec<usize> = (0..node_count).collect();
    let mut seen: HashSet<Vec<usize>> = HashSet::new();
    let mut queue: VecDeque<Vec<usize>> = VecDeque::new();
    seen.insert(identity.clone());
    queue.push_back(identity);

    while let Some(current) = queue.pop_front() {
        for gen in &gens {
            let next: Vec<usize> = current.iter().map(|&img| gen[img]).collect();
            if seen.insert(next.clone()) {
                queue.push_back(next);
            }
        }
    }

    let mut elements: Vec<Vec<usize>> = seen.into_iter().collect();
    elements.sort();
    elements
}

/// Selects the unique lexicographically minimum embedding along `priority_order`
/// from each `(query_group × target_group)` orbit in `embeddings`.
pub(super) fn canonical_orbit_representatives(
    embeddings: &[Vec<usize>],
    query_group: &[Vec<usize>],
    target_group: &[Vec<usize>],
    priority_order: &[usize],
) -> Vec<Vec<usize>> {
    assert!(!query_group.is_empty(), "query group must contain identity");
    assert!(
        !target_group.is_empty(),
        "target group must contain identity"
    );

    if query_group.len() == 1 && target_group.len() == 1 {
        let mut all = embeddings.to_vec();
        all.sort();
        return all;
    }

    let mut reps = Vec::new();
    for phi in embeddings {
        let mut is_min = true;
        'outer: for g in query_group {
            for h in target_group {
                for &p in priority_order {
                    let lhs = phi[p];
                    let rhs = h[phi[g[p]]];
                    if lhs != rhs {
                        if rhs < lhs {
                            is_min = false;
                            break 'outer;
                        }
                        break;
                    }
                }
            }
        }
        if is_min {
            reps.push(phi.clone());
        }
    }
    reps.sort();
    reps
}
