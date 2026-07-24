//! Captured external simple-cycle results.

use std::fs;

use rstest::rstest;
use umol_graph_core::SimpleCycleEnumerationAlgorithm::ReadTarjan;
use umol_graph_core::{EdgeId, Graph, NodeId};

use crate::corpus::parse_graph6;

const SIMPLE_CYCLES: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/data/simple-cycles/simple-through-8.tsv"
);
const MULTIGRAPH_CYCLES: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/data/simple-cycles/multigraph-through-4-edges-5.tsv"
);

fn parse_node_cycles(source: &str) -> Vec<Vec<NodeId>> {
    if source == "-" {
        return Vec::new();
    }
    source
        .split(';')
        .map(|cycle| {
            cycle
                .split(',')
                .map(|node| NodeId(node.parse().expect("captured node must be an integer")))
                .collect()
        })
        .collect()
}

fn parse_edge_cycles(source: &str) -> Vec<Vec<EdgeId>> {
    if source == "-" {
        return Vec::new();
    }
    source
        .split(';')
        .map(|cycle| {
            cycle
                .split(',')
                .map(|edge| EdgeId(edge.parse().expect("captured edge must be an integer")))
                .collect()
        })
        .collect()
}

fn parse_edges(source: &str) -> Vec<[u32; 2]> {
    if source.is_empty() {
        return Vec::new();
    }
    source
        .split(';')
        .map(|edge| {
            let (first, second) = edge
                .split_once(',')
                .expect("captured edge must contain two endpoints");
            [
                first.parse().expect("captured endpoint must be an integer"),
                second
                    .parse()
                    .expect("captured endpoint must be an integer"),
            ]
        })
        .collect()
}

#[rstest]
fn test_graph_enumerate_simple_cycles_simple_corpus() {
    let results =
        fs::read_to_string(SIMPLE_CYCLES).expect("captured simple cycles must be readable");
    let mut graph_count = 0;
    for line in results.lines() {
        let (source, expected) = line
            .split_once('\t')
            .expect("captured simple-cycle row must contain expected cycles");
        let graph = parse_graph6(source);
        let mut actual: Vec<Vec<NodeId>> = graph
            .enumerate_simple_cycles(usize::MAX, ReadTarjan)
            .into_iter()
            .map(|cycle| cycle.nodes().to_vec())
            .collect();
        actual.sort_by(|left, right| left.len().cmp(&right.len()).then_with(|| left.cmp(right)));

        assert_eq!(actual, parse_node_cycles(expected), "{source}");
        graph_count += 1;
    }
    assert_eq!(graph_count, 13_598);
}

#[rstest]
fn test_graph_enumerate_simple_cycles_multigraph_corpus() {
    let results =
        fs::read_to_string(MULTIGRAPH_CYCLES).expect("captured multigraph cycles must be readable");
    let mut graph_count = 0;
    for line in results.lines() {
        let mut fields = line.splitn(3, '\t');
        let node_count: usize = fields
            .next()
            .expect("captured multigraph row must contain a node count")
            .parse()
            .expect("captured node count must be an integer");
        let edges = parse_edges(
            fields
                .next()
                .expect("captured multigraph row must contain edges"),
        );
        let expected = fields
            .next()
            .expect("captured multigraph row must contain expected cycles");
        let graph = Graph::new(node_count, &edges);
        let mut actual: Vec<Vec<EdgeId>> = graph
            .enumerate_simple_cycles(usize::MAX, ReadTarjan)
            .into_iter()
            .map(|cycle| {
                let mut edges = cycle.edges().to_vec();
                edges.sort_unstable();
                edges
            })
            .collect();
        actual.sort_by(|left, right| left.len().cmp(&right.len()).then_with(|| left.cmp(right)));

        assert_eq!(actual, parse_edge_cycles(expected), "{graph:?}");
        graph_count += 1;
    }
    assert_eq!(graph_count, 3_453);
}
