#[path = "matching/fixture.rs"]
#[allow(dead_code)]
mod fixture;

#[cfg(feature = "proptest")]
use proptest::prelude::*;
#[cfg(feature = "proptest")]
use proptest::test_runner::{Config, TestRunner};
use rstest::rstest;
use umol_graph_core::{
    Correspondence, EdgeId, Graph, GraphCorrespondence, MatchingEnumerationAlgorithm, NodeId,
};
use MatchingEnumerationAlgorithm::BranchAndBound;

#[rstest]
#[case::bipartite(fixture::BENZENE, vec![5, 4, 3, 2, 1, 0])]
#[case::non_bipartite(fixture::AZULENE, vec![3, 8, 1, 6, 9, 0, 4, 7, 2, 5])]
#[case::disconnected(
    fixture::DISCONNECTED_CYCLES,
    (0..24_u32).rev().collect(),
)]
fn test_matching_enumeration_relabeling(#[case] source: &str, #[case] permutation: Vec<u32>) {
    let original = fixture::parse(source).graph();
    let relabeled_edges: Vec<_> = original
        .edge_ids()
        .map(|edge| {
            let [first, second] = original.edge_endpoints(edge);
            [permutation[first.index()], permutation[second.index()]]
        })
        .collect();
    let relabeled = Graph::new(original.node_count(), &relabeled_edges);
    let correspondence = GraphCorrespondence::induced(
        &original,
        &relabeled,
        Correspondence::from_images(
            &permutation.iter().copied().map(NodeId).collect::<Vec<_>>(),
            relabeled.node_count(),
        ),
    )
    .expect("relabeling a simple graph induces a unique graph correspondence");
    let canonical_original = |matchings: Vec<umol_graph_core::Matching>| {
        let mut canonical: Vec<_> = matchings
            .into_iter()
            .map(|matching| matching.edges().to_vec())
            .collect();
        canonical.sort_unstable();
        canonical
    };
    let canonical_relabeled = |matchings: Vec<umol_graph_core::Matching>| {
        let mut canonical: Vec<_> = matchings
            .into_iter()
            .map(|matching| {
                let mut edges: Vec<_> = matching
                    .edges()
                    .iter()
                    .map(|&edge| correspondence.edges().left_of(edge).unwrap())
                    .collect();
                edges.sort_unstable();
                edges
            })
            .collect();
        canonical.sort_unstable();
        canonical
    };

    assert_eq!(
        canonical_relabeled(relabeled.enumerate_perfect_matchings(BranchAndBound)),
        canonical_original(original.enumerate_perfect_matchings(BranchAndBound)),
    );
    assert_eq!(
        canonical_relabeled(relabeled.enumerate_maximum_matchings(BranchAndBound)),
        canonical_original(original.enumerate_maximum_matchings(BranchAndBound)),
    );
}

#[cfg(feature = "proptest")]
#[rstest]
fn test_matching_enumeration_relabeling_property() {
    const PROPERTY_CASES: u32 = 32;

    let original = fixture::parse(fixture::AZULENE).graph();
    let strategy = prop::collection::vec(any::<u64>(), original.node_count());
    let mut runner = TestRunner::new(Config {
        cases: PROPERTY_CASES,
        ..Config::default()
    });

    runner
        .run(&strategy, |keys| {
            let mut order: Vec<_> = (0..original.node_count()).collect();
            order.sort_unstable_by_key(|&node| (keys[node], node));
            let mut permutation = vec![0_u32; original.node_count()];
            for (relabeled, original_node) in order.into_iter().enumerate() {
                permutation[original_node] = relabeled as u32;
            }
            let relabeled_edges: Vec<_> = original
                .edge_ids()
                .map(|edge| {
                    let [first, second] = original.edge_endpoints(edge);
                    [permutation[first.index()], permutation[second.index()]]
                })
                .collect();
            let relabeled = Graph::new(original.node_count(), &relabeled_edges);
            let correspondence = GraphCorrespondence::induced(
                &original,
                &relabeled,
                Correspondence::from_images(
                    &permutation.iter().copied().map(NodeId).collect::<Vec<_>>(),
                    relabeled.node_count(),
                ),
            )
            .expect("relabeling a simple graph induces a unique graph correspondence");

            for (original_matchings, relabeled_matchings) in [
                (
                    original.enumerate_perfect_matchings(BranchAndBound),
                    relabeled.enumerate_perfect_matchings(BranchAndBound),
                ),
                (
                    original.enumerate_maximum_matchings(BranchAndBound),
                    relabeled.enumerate_maximum_matchings(BranchAndBound),
                ),
            ] {
                let mut expected: Vec<_> = original_matchings
                    .into_iter()
                    .map(|matching| matching.edges().to_vec())
                    .collect();
                expected.sort_unstable();
                let mut actual: Vec<_> = relabeled_matchings
                    .into_iter()
                    .map(|matching| {
                        let mut edges: Vec<_> = matching
                            .edges()
                            .iter()
                            .map(|&edge| correspondence.edges().left_of(edge).unwrap())
                            .collect();
                        edges.sort_unstable();
                        edges
                    })
                    .collect();
                actual.sort_unstable();
                prop_assert_eq!(actual, expected);
            }
            Ok(())
        })
        .unwrap();
}

#[rstest]
#[case::benzene(fixture::BENZENE, vec![NodeId(0), NodeId(1)], 1)]
#[case::ladder(fixture::LADDER, vec![NodeId(0), NodeId(4)], 3)]
fn test_matching_enumeration_holes(
    #[case] source: &str,
    #[case] holes: Vec<NodeId>,
    #[case] expected_count: usize,
) {
    let graph = fixture::parse(source).graph();
    let retained: Vec<_> = graph
        .node_ids()
        .filter(|node| !holes.contains(node))
        .collect();
    let subgraph = graph.induced_subgraph(&retained);
    let residual = graph.extract(&subgraph);
    let matchings = residual.enumerate_perfect_matchings(BranchAndBound);

    assert_eq!(matchings.len(), expected_count);
    for matching in matchings {
        let mut covered = vec![false; graph.node_count()];
        for &residual_edge in matching.edges() {
            let original_edge = subgraph.edges().right_of(residual_edge).unwrap();
            let [first, second] = graph.edge_endpoints(original_edge);
            covered[first.index()] = true;
            covered[second.index()] = true;
        }
        let exposed: Vec<_> = graph
            .node_ids()
            .filter(|node| !covered[node.index()])
            .collect();
        assert_eq!(exposed, holes);
        assert!(retained.iter().all(|node| covered[node.index()]));
    }
}

#[rstest]
fn test_matching_enumeration_mobile_hole() {
    let graph = Graph::new(5, &[[0, 1], [1, 2], [2, 3], [3, 4], [4, 0]]);
    let mut deleted_vertex_matchings = Vec::new();

    for hole in graph.node_ids() {
        let retained: Vec<_> = graph.node_ids().filter(|&node| node != hole).collect();
        let subgraph = graph.induced_subgraph(&retained);
        let residual = graph.extract(&subgraph);
        for matching in residual.enumerate_perfect_matchings(BranchAndBound) {
            let mut covered = vec![false; graph.node_count()];
            let mut edges: Vec<EdgeId> = matching
                .edges()
                .iter()
                .map(|&edge| subgraph.edges().right_of(edge).unwrap())
                .collect();
            edges.sort_unstable();
            for &edge in &edges {
                let [first, second] = graph.edge_endpoints(edge);
                covered[first.index()] = true;
                covered[second.index()] = true;
            }
            assert_eq!(
                graph
                    .node_ids()
                    .filter(|node| !covered[node.index()])
                    .collect::<Vec<_>>(),
                vec![hole],
            );
            deleted_vertex_matchings.push(edges);
        }
    }
    deleted_vertex_matchings.sort_unstable();
    let mut maximum_matchings: Vec<_> = graph
        .enumerate_maximum_matchings(BranchAndBound)
        .into_iter()
        .map(|matching| matching.edges().to_vec())
        .collect();
    maximum_matchings.sort_unstable();

    assert_eq!(deleted_vertex_matchings, maximum_matchings);
    assert_eq!(deleted_vertex_matchings.len(), 5);
}
