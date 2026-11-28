//! Captured external simple-cycle results.

use std::fs;

use rstest::rstest;
use umol_graph_core::NodeId;
use umol_graph_core::SimpleCycleEnumerationAlgorithm::ReadTarjan;

use crate::corpus::parse_graph6;

const SIMPLE_CYCLES: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/data/simple-cycles/simple-through-8.tsv"
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
