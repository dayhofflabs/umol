//! Property-based cross-validation of the subgraph-isomorphism algorithms: on
//! random labeled graphs, every named algorithm must return the same match set as
//! VF2 (the reference). For ArcMatch this is a complete correctness oracle — its
//! reduction is provably safe (Bonnici 2024 Thm 2-3), so the only possible defect
//! is returning *fewer* matches, which a set mismatch against VF2 catches.

use std::collections::HashMap;

use proptest::prelude::*;
use umol_graph_core::SubgraphIsomorphismAlgorithm::{
    ArcMatch, RayKirsch, Ri, Ullmann, Vf2, Vf2Rdkit,
};
use umol_graph_core::{
    EdgeId, Graph, NodeId, SubgraphIsomorphismAlgorithm, ARCMATCH_DEFAULT_PATH_LENGTH,
};

/// A graph with a label per node and per edge. Edge `k` is the `k`-th edge passed
/// to `Graph::new`, so `edge_labels[k]` is the label of `EdgeId(k)`.
#[derive(Clone, Debug)]
struct LabeledGraph {
    graph: Graph,
    node_labels: Vec<u32>,
    edge_labels: Vec<u32>,
}

/// Random labeled graph of up to `max_n` nodes: each undirected pair is present
/// with even odds and, when present, carries an edge label from `0..edge_labels`;
/// each node is labeled from `0..node_labels`.
fn labeled_graph(
    max_n: usize,
    node_labels: u32,
    edge_labels: u32,
) -> impl Strategy<Value = LabeledGraph> {
    (1..=max_n)
        .prop_flat_map(move |n| {
            let pairs: Vec<(u32, u32)> = (0..n as u32)
                .flat_map(|i| ((i + 1)..n as u32).map(move |j| (i, j)))
                .collect();
            let mask = prop::collection::vec(
                prop_oneof![Just(None::<u32>), (0..edge_labels).prop_map(Some)],
                pairs.len(),
            );
            let labels = prop::collection::vec(0..node_labels, n);
            (Just(n), Just(pairs), labels, mask)
        })
        .prop_map(|(n, pairs, node_labels, mask)| {
            let mut edges = Vec::new();
            let mut edge_labels = Vec::new();
            for (pair, slot) in pairs.iter().zip(mask) {
                if let Some(label) = slot {
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

/// The induced subgraph of `target` on `nodes` (distinct, ascending), relabeled to
/// `0..nodes.len()` and inheriting node/edge labels — embeds in `target` by
/// construction, guaranteeing at least one match.
fn induced(target: &LabeledGraph, nodes: &[usize]) -> LabeledGraph {
    let remap: HashMap<usize, u32> = nodes
        .iter()
        .enumerate()
        .map(|(new, &old)| (old, new as u32))
        .collect();
    let mut edges = Vec::new();
    let mut edge_labels = Vec::new();
    for eid in 0..target.graph.edge_count() {
        let [a, b] = target.graph.edge_endpoints(EdgeId(eid as u32));
        if let (Some(&na), Some(&nb)) = (remap.get(&a.index()), remap.get(&b.index())) {
            edges.push([na, nb]);
            edge_labels.push(target.edge_labels[eid]);
        }
    }
    let node_labels = nodes.iter().map(|&old| target.node_labels[old]).collect();
    LabeledGraph {
        graph: Graph::new(nodes.len(), &edges),
        node_labels,
        edge_labels,
    }
}

/// Sorted match set under `alg`, with label-equality node and edge matching.
fn matches(
    query: &LabeledGraph,
    target: &LabeledGraph,
    alg: SubgraphIsomorphismAlgorithm,
) -> Vec<Vec<usize>> {
    let mut node_match =
        |q: NodeId, t: NodeId| query.node_labels[q.index()] == target.node_labels[t.index()];
    let mut edge_match =
        |qe: EdgeId, te: EdgeId| query.edge_labels[qe.index()] == target.edge_labels[te.index()];
    let mut found =
        target
            .graph
            .subgraph_isomorphisms(&query.graph, &mut node_match, &mut edge_match, alg);
    found.sort();
    found
}

proptest! {
    #[test]
    fn test_subgraph_isomorphisms_cross_validation(
        query in labeled_graph(4, 2, 2),
        target in labeled_graph(6, 2, 2),
    ) {
        let reference = matches(&query, &target, Vf2);
        for alg in [Ullmann, Ri, ArcMatch { path_length: ARCMATCH_DEFAULT_PATH_LENGTH }, Vf2Rdkit, RayKirsch] {
            prop_assert_eq!(
                matches(&query, &target, alg),
                reference.clone(),
                "{:?} disagrees with Vf2",
                alg
            );
        }
    }

    #[test]
    fn test_subgraph_isomorphisms_cross_validation_planted(
        (target, nodes) in labeled_graph(6, 2, 2).prop_flat_map(|t| {
            let cap = t.graph.node_count().min(4);
            let n = t.graph.node_count();
            (Just(t), prop::sample::subsequence((0..n).collect::<Vec<_>>(), 1..=cap))
        }),
    ) {
        let query = induced(&target, &nodes);
        let reference = matches(&query, &target, Vf2);
        prop_assert!(!reference.is_empty(), "planted query must embed in its source");
        for alg in [Ullmann, Ri, ArcMatch { path_length: ARCMATCH_DEFAULT_PATH_LENGTH }, Vf2Rdkit, RayKirsch] {
            prop_assert_eq!(
                matches(&query, &target, alg),
                reference.clone(),
                "{:?} disagrees with Vf2",
                alg
            );
        }
    }
}
