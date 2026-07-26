//! Captured external cycle-family results.

use std::cmp::Ordering;
use std::fs;
use std::ops::ControlFlow;

use rstest::rstest;
use umol_graph_core::{
    MinimumCycleBasisAlgorithm, RelevantCycleEnumerationAlgorithm, UniqueRingFamilyAlgorithm,
};

use crate::corpus::parse_graph6;

const CYCLE_FAMILIES: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/data/cycles/cycle-families-through-8.tsv"
);

#[derive(Debug, PartialEq, Eq)]
struct CapturedFamily {
    weight: usize,
    count: usize,
    nodes: Vec<u32>,
    edges: Vec<u32>,
    cycles: Vec<Vec<u32>>,
}

fn parse_numbers(source: &str) -> Vec<u32> {
    if source.is_empty() {
        return Vec::new();
    }
    source
        .split(',')
        .map(|value| {
            value
                .parse()
                .expect("captured identifier must be an integer")
        })
        .collect()
}

fn parse_cycles(source: &str, separator: char) -> Vec<Vec<u32>> {
    if source == "-" {
        return Vec::new();
    }
    source.split(separator).map(parse_numbers).collect()
}

fn parse_families(source: &str) -> Vec<CapturedFamily> {
    if source == "-" {
        return Vec::new();
    }
    source
        .split('/')
        .map(|family| {
            let mut fields = family.splitn(5, ':');
            CapturedFamily {
                weight: fields
                    .next()
                    .expect("captured family must contain a weight")
                    .parse()
                    .expect("captured family weight must be an integer"),
                count: fields
                    .next()
                    .expect("captured family must contain a cycle count")
                    .parse()
                    .expect("captured family count must be an integer"),
                nodes: parse_numbers(fields.next().expect("captured family must contain nodes")),
                edges: parse_numbers(fields.next().expect("captured family must contain edges")),
                cycles: parse_cycles(
                    fields.next().expect("captured family must contain cycles"),
                    '.',
                ),
            }
        })
        .collect()
}

fn compare_cycles(left: &[u32], right: &[u32]) -> Ordering {
    left.len().cmp(&right.len()).then_with(|| left.cmp(right))
}

fn compare_families(left: &CapturedFamily, right: &CapturedFamily) -> Ordering {
    left.weight
        .cmp(&right.weight)
        .then_with(|| left.cycles.cmp(&right.cycles))
        .then_with(|| left.nodes.cmp(&right.nodes))
        .then_with(|| left.edges.cmp(&right.edges))
}

#[rstest]
fn test_graph_cycle_families_corpus() {
    let results =
        fs::read_to_string(CYCLE_FAMILIES).expect("captured cycle families must be readable");
    let mut graph_count = 0;
    for line in results.lines() {
        let mut fields = line.splitn(6, '\t');
        let source = fields
            .next()
            .expect("captured cycle-family row must contain graph6");
        let expected_dimension = fields
            .next()
            .expect("captured cycle-family row must contain an MCB dimension")
            .parse()
            .expect("captured MCB dimension must be an integer");
        let expected_total_length = fields
            .next()
            .expect("captured cycle-family row must contain an MCB total length")
            .parse()
            .expect("captured MCB total length must be an integer");
        let expected_relevant = parse_cycles(
            fields
                .next()
                .expect("captured cycle-family row must contain relevant cycles"),
            ';',
        );
        let mut expected_families = parse_families(
            fields
                .next()
                .expect("captured cycle-family row must contain URFs"),
        );
        let validation = fields
            .next()
            .expect("captured cycle-family row must contain validation status");
        let graph = parse_graph6(source);

        let basis = graph.minimum_cycle_basis(MinimumCycleBasisAlgorithm::Horton);
        assert_eq!(basis.dimension(), expected_dimension, "{source}");
        assert_eq!(basis.total_length(), expected_total_length, "{source}");

        let mut relevant = graph
            .enumerate_relevant_cycles(usize::MAX, RelevantCycleEnumerationAlgorithm::Vismara)
            .into_iter()
            .map(|cycle| {
                let mut edges = cycle
                    .edges()
                    .iter()
                    .map(|edge| edge.index() as u32)
                    .collect::<Vec<_>>();
                edges.sort_unstable();
                edges
            })
            .collect::<Vec<_>>();
        relevant.sort_by(|left, right| compare_cycles(left, right));
        assert_eq!(relevant, expected_relevant, "{source}");

        let decomposition = graph.unique_ring_families(UniqueRingFamilyAlgorithm::Kolodzik);
        let mut families = decomposition
            .ids()
            .map(|id| {
                let family = decomposition
                    .get(id)
                    .expect("a returned family id must be valid");
                let mut cycles = Vec::new();
                let flow = decomposition.visit_relevant_cycles(id, |cycle| {
                    let mut edges = cycle
                        .edges()
                        .iter()
                        .map(|edge| edge.index() as u32)
                        .collect::<Vec<_>>();
                    edges.sort_unstable();
                    cycles.push(edges);
                    ControlFlow::<()>::Continue(())
                });
                assert_eq!(flow, ControlFlow::Continue(()), "{source}");
                cycles.sort_by(|left, right| compare_cycles(left, right));

                CapturedFamily {
                    weight: family.weight(),
                    count: family
                        .relevant_cycle_count()
                        .0
                        .to_string()
                        .parse()
                        .expect("captured-range family count must fit usize"),
                    nodes: family
                        .nodes()
                        .iter()
                        .map(|node| node.index() as u32)
                        .collect(),
                    edges: family
                        .edges()
                        .iter()
                        .map(|edge| edge.index() as u32)
                        .collect(),
                    cycles,
                }
            })
            .collect::<Vec<_>>();
        families.sort_by(compare_families);
        expected_families.sort_by(compare_families);
        assert_eq!(families, expected_families, "{source}");

        let expected_validation = if expected_relevant.is_empty() {
            "-"
        } else if graph.node_count() <= 5 {
            "1"
        } else {
            "0"
        };
        assert_eq!(validation, expected_validation, "{source}");
        graph_count += 1;
    }
    assert_eq!(graph_count, 13_598);
}
