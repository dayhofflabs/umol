//! Cycle properties use definition-level exhaustive enumeration for bounded
//! random simple graphs and the simple corpus through order six. `captured` and
//! `families` independently compare the production operations with captured
//! results over the larger simple corpus. `literature` preserves
//! source-stated results for selected chemical and graph-theoretic examples;
//! its intentional overlap with the exhaustive corpus records provenance, while
//! the para-bridged macrocycle extends coverage to a larger chemical topology.

#[path = "cycles/captured.rs"]
mod captured;
#[path = "cycles/exhaustive.rs"]
mod exhaustive;
#[path = "cycles/families.rs"]
mod families;
#[path = "cycles/literature.rs"]
mod literature;

use std::collections::HashSet;
use std::ops::ControlFlow;

use num_bigint::BigUint;
use proptest::prelude::*;
use rstest::rstest;
use umol_graph_core::{
    EdgeId, Graph, MinimumCycleBasisAlgorithm, RelevantCycleEnumerationAlgorithm,
    SimpleCycleEnumerationAlgorithm, UniqueRingFamilyAlgorithm, UniqueRingFamilyId,
};

use self::exhaustive::{
    are_linearly_independent, cycle_space_rank, enumerate_cycles, minimum_cycle_bases,
    relevant_cycles, unique_ring_families,
};
use super::corpus::simple_graphs;
use super::strategy::simple_graph;

const SIMPLE_CYCLE_ALGORITHM: SimpleCycleEnumerationAlgorithm =
    SimpleCycleEnumerationAlgorithm::ReadTarjan;

#[rstest]
fn test_graph_enumerate_simple_cycles_corpus() {
    let mut graph_count = 0;
    for graph in simple_graphs().take_while(|graph| graph.node_count() <= 6) {
        let cycles = graph.enumerate_simple_cycles(usize::MAX, SIMPLE_CYCLE_ALGORITHM);

        let mut actual: Vec<Vec<EdgeId>> = cycles
            .into_iter()
            .map(|cycle| {
                let mut edges = cycle.edges().to_vec();
                edges.sort_unstable();
                edges
            })
            .collect();
        actual.sort_by(|left, right| left.len().cmp(&right.len()).then_with(|| left.cmp(right)));

        assert_eq!(actual, enumerate_cycles(&graph), "{graph:?}");
        graph_count += 1;
    }
    assert_eq!(graph_count, 208);
}

proptest! {
    #[test]
    fn test_graph_enumerate_simple_cycles(
        graph in simple_graph(5, 10),
        max_cycle_size in 0usize..=5,
    ) {
        let cycles = graph.enumerate_simple_cycles(max_cycle_size, SIMPLE_CYCLE_ALGORITHM);
        let mut actual: Vec<Vec<EdgeId>> = cycles
            .iter()
            .map(|cycle| {
                let mut edges = cycle.edges().to_vec();
                edges.sort_unstable();
                edges
            })
            .collect();
        actual.sort_by(|left, right| left.len().cmp(&right.len()).then_with(|| left.cmp(right)));

        let expected: Vec<Vec<EdgeId>> = enumerate_cycles(&graph)
            .into_iter()
            .filter(|cycle| cycle.len() <= max_cycle_size)
            .collect();
        prop_assert_eq!(actual, expected);

        let mut unique = HashSet::new();
        for cycle in &cycles {
            prop_assert!(unique.insert(cycle));
        }
    }

    #[test]
    fn test_graph_visit_simple_cycles(
        graph in simple_graph(5, 10),
        max_cycle_size in 0usize..=5,
    ) {
        let expected = graph.enumerate_simple_cycles(max_cycle_size, SIMPLE_CYCLE_ALGORITHM);

        let mut visited = Vec::new();
        let result = graph.visit_simple_cycles(max_cycle_size, SIMPLE_CYCLE_ALGORITHM, |cycle| {
            visited.push(cycle);
            ControlFlow::<()>::Continue(())
        });
        prop_assert_eq!(result, ControlFlow::Continue(()));
        prop_assert_eq!(&visited, &expected);
    }

    #[test]
    fn test_graph_enumerate_simple_cycles_relabeling(
        graph in simple_graph(5, 10),
        max_cycle_size in 0usize..=5,
    ) {
        let node_count = graph.node_count();
        let relabeled_edges: Vec<[u32; 2]> = graph
            .edge_ids()
            .map(|edge| {
                let [first, second] = graph.edge_endpoints(edge);
                [
                    (node_count - 1 - first.index()) as u32,
                    (node_count - 1 - second.index()) as u32,
                ]
            })
            .collect();
        let relabeled = Graph::new(node_count, &relabeled_edges);

        let edge_sets = |source: &Graph| {
            let mut cycles: Vec<Vec<EdgeId>> = source
                .enumerate_simple_cycles(max_cycle_size, SIMPLE_CYCLE_ALGORITHM)
                .into_iter()
                .map(|cycle| {
                    let mut edges = cycle.edges().to_vec();
                    edges.sort_unstable();
                    edges
                })
                .collect();
            cycles.sort_by(|left, right| left.len().cmp(&right.len()).then_with(|| left.cmp(right)));
            cycles
        };

        prop_assert_eq!(edge_sets(&relabeled), edge_sets(&graph));
    }

    #[test]
    fn test_graph_minimum_cycle_basis(graph in simple_graph(5, 10)) {
        let result = graph.minimum_cycle_basis(MinimumCycleBasisAlgorithm::Horton);
        let cycles = result
            .iter()
            .map(|cycle| cycle.edges().to_vec())
            .collect::<Vec<_>>();
        let expected_rank = cycle_space_rank(&graph);
        let expected_length = minimum_cycle_bases(&graph)
            .first()
            .expect("the exhaustive cycle set spans the cycle space")
            .iter()
            .map(Vec::len)
            .sum::<usize>();

        prop_assert_eq!(result.dimension(), expected_rank);
        prop_assert!(are_linearly_independent(&cycles, graph.edge_count()));
        prop_assert_eq!(cycles.len(), expected_rank);
        prop_assert_eq!(result.total_length(), expected_length);
    }

    #[test]
    fn test_graph_visit_relevant_cycles(
        graph in simple_graph(5, 10),
        max_cycle_size in 0usize..=5,
    ) {
        let expected = graph.enumerate_relevant_cycles(
            max_cycle_size,
            RelevantCycleEnumerationAlgorithm::Vismara,
        );

        let mut visited = Vec::new();
        let result = graph.visit_relevant_cycles(
            max_cycle_size,
            RelevantCycleEnumerationAlgorithm::Vismara,
            |cycle| {
                visited.push(cycle);
                ControlFlow::<()>::Continue(())
            },
        );
        prop_assert_eq!(result, ControlFlow::Continue(()));
        prop_assert_eq!(&visited, &expected);
    }

    #[test]
    fn test_graph_enumerate_relevant_cycles(
        graph in simple_graph(5, 10),
        max_cycle_size in 0usize..=5,
    ) {
        let mut actual = graph
            .enumerate_relevant_cycles(
                max_cycle_size,
                RelevantCycleEnumerationAlgorithm::Vismara,
            )
            .into_iter()
            .map(|cycle| {
                let mut edges = cycle.edges().to_vec();
                edges.sort_unstable();
                edges
            })
            .collect::<Vec<_>>();
        actual.sort_by(|left, right| {
            left.len()
                .cmp(&right.len())
                .then_with(|| left.cmp(right))
        });

        let expected = relevant_cycles(&graph)
            .into_iter()
            .filter(|cycle| cycle.len() <= max_cycle_size)
            .collect::<Vec<_>>();
        prop_assert_eq!(actual, expected);
    }

    #[test]
    fn test_graph_unique_ring_families(
        graph in simple_graph(5, 10),
        prefix in any::<usize>(),
    ) {
        let result = graph.unique_ring_families(UniqueRingFamilyAlgorithm::Kolodzik);
        let count = result.count();
        let prefix = prefix.min(count);
        let mut ids = result.ids();
        prop_assert_eq!(ids.len(), count);
        prop_assert_eq!(ids.size_hint(), (count, Some(count)));
        for index in 0..prefix {
            prop_assert_eq!(ids.next(), Some(UniqueRingFamilyId(index as u32)));
            let remaining = count - index - 1;
            prop_assert_eq!(ids.len(), remaining);
            prop_assert_eq!(ids.size_hint(), (remaining, Some(remaining)));
        }
        prop_assert_eq!(
            ids.collect::<Vec<_>>(),
            (prefix..count)
                .map(|index| UniqueRingFamilyId(index as u32))
                .collect::<Vec<_>>(),
        );

        let mut actual = Vec::new();
        for id in result.ids() {
            let family = result.get(id).expect("a returned family id must be valid");
            let mut cycles = Vec::new();
            let flow = result.visit_relevant_cycles(id, |cycle| {
                let mut edges = cycle.edges().to_vec();
                edges.sort_unstable();
                cycles.push(edges);
                ControlFlow::<()>::Continue(())
            });
            cycles.sort();

            let mut nodes = cycles
                .iter()
                .flatten()
                .flat_map(|&edge| graph.edge_endpoints(edge))
                .collect::<Vec<_>>();
            nodes.sort_unstable();
            nodes.dedup();
            let mut edges = cycles.iter().flatten().copied().collect::<Vec<_>>();
            edges.sort_unstable();
            edges.dedup();

            prop_assert_eq!(flow, ControlFlow::Continue(()));
            prop_assert_eq!(family.nodes(), nodes);
            prop_assert_eq!(family.edges(), edges);
            prop_assert_eq!(family.weight(), cycles[0].len());
            prop_assert_eq!(
                &family.relevant_cycle_count().0,
                &BigUint::from(cycles.len())
            );
            actual.push(cycles);
        }
        actual.sort();

        let mut expected = unique_ring_families(&graph);
        expected.sort();
        prop_assert_eq!(actual, expected);
    }

    #[test]
    fn test_unique_ring_families_families_containing_node(graph in simple_graph(5, 10)) {
        let families = graph.unique_ring_families(UniqueRingFamilyAlgorithm::Kolodzik);
        for node in graph.node_ids() {
            let expected = families
                .ids()
                .filter(|&id| {
                    families
                        .get(id)
                        .expect("a returned family id must be valid")
                        .nodes()
                        .contains(&node)
                })
                .collect::<Vec<_>>();
            prop_assert_eq!(families.families_containing_node(node), expected);
        }
    }

    #[test]
    fn test_unique_ring_families_families_containing_edge(graph in simple_graph(5, 10)) {
        let families = graph.unique_ring_families(UniqueRingFamilyAlgorithm::Kolodzik);
        for edge in graph.edge_ids() {
            let expected = families
                .ids()
                .filter(|&id| {
                    families
                        .get(id)
                        .expect("a returned family id must be valid")
                        .edges()
                        .contains(&edge)
                })
                .collect::<Vec<_>>();
            prop_assert_eq!(families.families_containing_edge(edge), expected);
        }
    }

    #[test]
    fn test_graph_unique_ring_families_relabeling(graph in simple_graph(5, 10)) {
        let node_count = graph.node_count();
        let relabeled_edges = graph
            .edge_ids()
            .map(|edge| {
                let [first, second] = graph.edge_endpoints(edge);
                [
                    (node_count - 1 - first.index()) as u32,
                    (node_count - 1 - second.index()) as u32,
                ]
            })
            .collect::<Vec<_>>();
        let relabeled = Graph::new(node_count, &relabeled_edges);

        let decomposition = |source: &Graph| {
            let result = source.unique_ring_families(UniqueRingFamilyAlgorithm::Kolodzik);
            let mut families = result
                .ids()
                .map(|id| {
                    let family = result.get(id).expect("a returned family id must be valid");
                    let mut cycles = Vec::new();
                    let _: ControlFlow<()> = result.visit_relevant_cycles(id, |cycle| {
                        let mut edges = cycle.edges().to_vec();
                        edges.sort_unstable();
                        cycles.push(edges);
                        ControlFlow::Continue(())
                    });
                    cycles.sort();
                    (
                        family.edges().to_vec(),
                        family.weight(),
                        family.relevant_cycle_count().clone(),
                        cycles,
                    )
                })
                .collect::<Vec<_>>();
            families.sort();
            families
        };
        prop_assert_eq!(decomposition(&graph), decomposition(&relabeled));
    }

    #[test]
    fn test_unique_ring_families(graph in simple_graph(5, 10)) {
        let relevant = relevant_cycles(&graph);
        let families = unique_ring_families(&graph);

        for family in &families {
            prop_assert!(!family.is_empty());
            let weight = family[0].len();
            prop_assert!(family.iter().all(|cycle| cycle.len() == weight));
        }

        let mut partition: Vec<_> = families.into_iter().flatten().collect();
        partition.sort_by(|left, right| {
            left.len()
                .cmp(&right.len())
                .then_with(|| left.cmp(right))
        });
        prop_assert_eq!(partition, relevant);
    }
}
