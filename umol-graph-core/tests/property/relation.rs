//! Exact-size iterator properties and participant-transport laws.
//!
//! Transport laws use permutations of eight node and edge ids; unit cases cover partial mappings.

use proptest::prelude::*;
use proptest::test_runner::TestCaseResult;
use umol_graph_core::{
    Correspondence, EdgeId, FixedFixedBirelationSet, FixedRelationSet, FixedVarBirelationSet,
    GraphCorrespondence, GraphRemapping, NodeId, RelationId, VarRelationSet, VarVarBirelationSet,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct TestData(usize);

fn participant_mapping_strategy() -> impl Strategy<Value = (GraphCorrespondence, GraphRemapping)> {
    (
        Just((0..8).map(NodeId).collect::<Vec<_>>()).prop_shuffle(),
        Just((0..8).map(EdgeId).collect::<Vec<_>>()).prop_shuffle(),
    )
        .prop_map(|(nodes, edges)| {
            (
                GraphCorrespondence::new(
                    Correspondence::from_images(&nodes, 8),
                    Correspondence::from_images(&edges, 8),
                ),
                GraphRemapping::new(nodes, edges),
            )
        })
}

fn assert_relation_ids(
    mut iterator: impl ExactSizeIterator<Item = RelationId>,
    count: usize,
    prefix: usize,
) -> TestCaseResult {
    prop_assert_eq!(iterator.len(), count);
    prop_assert_eq!(iterator.size_hint(), (count, Some(count)));
    for index in 0..prefix {
        prop_assert_eq!(iterator.next(), Some(RelationId(index as u32)));
        let remaining = count - index - 1;
        prop_assert_eq!(iterator.len(), remaining);
        prop_assert_eq!(iterator.size_hint(), (remaining, Some(remaining)));
    }
    prop_assert_eq!(
        iterator.collect::<Vec<_>>(),
        (prefix..count)
            .map(|index| RelationId(index as u32))
            .collect::<Vec<_>>(),
    );
    Ok(())
}

fn assert_data_iter_mut<'a>(
    mut iterator: impl ExactSizeIterator<Item = &'a mut TestData>,
    count: usize,
    prefix: usize,
) -> TestCaseResult {
    prop_assert_eq!(iterator.len(), count);
    prop_assert_eq!(iterator.size_hint(), (count, Some(count)));
    for index in 0..prefix {
        let item = iterator.next().expect("expected generated data item");
        prop_assert_eq!(*item, TestData(index));
        item.0 += count;
        let remaining = count - index - 1;
        prop_assert_eq!(iterator.len(), remaining);
        prop_assert_eq!(iterator.size_hint(), (remaining, Some(remaining)));
    }
    for index in prefix..count {
        let item = iterator.next().expect("expected remaining data item");
        prop_assert_eq!(*item, TestData(index));
        item.0 += count;
        let remaining = count - index - 1;
        prop_assert_eq!(iterator.len(), remaining);
        prop_assert_eq!(iterator.size_hint(), (remaining, Some(remaining)));
    }
    prop_assert_eq!(iterator.next(), None);
    prop_assert_eq!(iterator.len(), 0);
    Ok(())
}

proptest! {

    #[test]
    fn test_fixed_relation_set_map_composition(
        rows in prop::collection::vec((prop::array::uniform2(0u32..8), prop::array::uniform2(0u32..8), any::<[u32; 2]>()), 0..8),
        (first, remapping) in participant_mapping_strategy(),
        (second, _) in participant_mapping_strategy(),
    ) {
        let input: FixedRelationSet<NodeId, [u32; 2], 2> = FixedRelationSet::new(rows.into_iter().map(|(_edges, nodes, data)| (nodes.map(NodeId), data)).collect());
        let mapped = input.map(&first);
        prop_assert_eq!(&mapped, &input.remap(&remapping));
        let composed = first.compose(&second).unwrap();
        prop_assert_eq!(mapped.map(&second), input.map(&composed));
        let inverse = GraphCorrespondence::new(first.nodes().reverse(), first.edges().reverse());
        prop_assert_eq!(mapped.map(&inverse), input);
    }

    #[test]
    fn test_var_relation_set_map_composition(
        rows in prop::collection::vec((prop::array::uniform2(0u32..8), prop::array::uniform2(0u32..8), any::<[u32; 2]>()), 0..8),
        (first, remapping) in participant_mapping_strategy(),
        (second, _) in participant_mapping_strategy(),
    ) {
        let input: VarRelationSet<NodeId, [u32; 2]> = VarRelationSet::new(rows.into_iter().map(|(_edges, nodes, data)| (nodes.map(NodeId).to_vec(), data)).collect());
        let mapped = input.map(&first);
        prop_assert_eq!(&mapped, &input.remap(&remapping));
        let composed = first.compose(&second).unwrap();
        prop_assert_eq!(mapped.map(&second), input.map(&composed));
        let inverse = GraphCorrespondence::new(first.nodes().reverse(), first.edges().reverse());
        prop_assert_eq!(mapped.map(&inverse), input);
    }

    #[test]
    fn test_fixed_fixed_birelation_set_map_composition(
        rows in prop::collection::vec((prop::array::uniform2(0u32..8), prop::array::uniform2(0u32..8), any::<[u32; 2]>()), 0..8),
        (first, remapping) in participant_mapping_strategy(),
        (second, _) in participant_mapping_strategy(),
    ) {
        let input: FixedFixedBirelationSet<EdgeId, 2, NodeId, 2, [u32; 2]> = FixedFixedBirelationSet::new(rows.into_iter().map(|(edges, nodes, data)| (edges.map(EdgeId), nodes.map(NodeId), data)).collect());
        let mapped = input.map(&first);
        prop_assert_eq!(&mapped, &input.remap(&remapping));
        let composed = first.compose(&second).unwrap();
        prop_assert_eq!(mapped.map(&second), input.map(&composed));
        let inverse = GraphCorrespondence::new(first.nodes().reverse(), first.edges().reverse());
        prop_assert_eq!(mapped.map(&inverse), input);
    }

    #[test]
    fn test_fixed_var_birelation_set_map_composition(
        rows in prop::collection::vec((prop::array::uniform2(0u32..8), prop::array::uniform2(0u32..8), any::<[u32; 2]>()), 0..8),
        (first, remapping) in participant_mapping_strategy(),
        (second, _) in participant_mapping_strategy(),
    ) {
        let input: FixedVarBirelationSet<EdgeId, 2, NodeId, [u32; 2]> = FixedVarBirelationSet::new(rows.into_iter().map(|(edges, nodes, data)| (edges.map(EdgeId), nodes.map(NodeId).to_vec(), data)).collect());
        let mapped = input.map(&first);
        prop_assert_eq!(&mapped, &input.remap(&remapping));
        let composed = first.compose(&second).unwrap();
        prop_assert_eq!(mapped.map(&second), input.map(&composed));
        let inverse = GraphCorrespondence::new(first.nodes().reverse(), first.edges().reverse());
        prop_assert_eq!(mapped.map(&inverse), input);
    }

    #[test]
    fn test_var_var_birelation_set_map_composition(
        rows in prop::collection::vec((prop::array::uniform2(0u32..8), prop::array::uniform2(0u32..8), any::<[u32; 2]>()), 0..8),
        (first, remapping) in participant_mapping_strategy(),
        (second, _) in participant_mapping_strategy(),
    ) {
        let input: VarVarBirelationSet<EdgeId, NodeId, [u32; 2]> = VarVarBirelationSet::new(rows.into_iter().map(|(edges, nodes, data)| (edges.map(EdgeId).to_vec(), nodes.map(NodeId).to_vec(), data)).collect());
        let mapped = input.map(&first);
        prop_assert_eq!(&mapped, &input.remap(&remapping));
        let composed = first.compose(&second).unwrap();
        prop_assert_eq!(mapped.map(&second), input.map(&composed));
        let inverse = GraphCorrespondence::new(first.nodes().reverse(), first.edges().reverse());
        prop_assert_eq!(mapped.map(&inverse), input);
    }
    #[test]
    fn test_fixed_relation_set_iterators_exact_size(
        count in 0usize..16,
        prefix in any::<usize>(),
    ) {
        let entries = (0..count)
            .map(|index| {
                (
                    [NodeId((2 * index) as u32), NodeId((2 * index + 1) as u32)],
                    TestData(index),
                )
            })
            .collect::<Vec<_>>();
        let mut relations = FixedRelationSet::<NodeId, TestData, 2>::new(entries);
        let prefix = prefix.min(count);

        assert_relation_ids(relations.ids(), count, prefix)?;
        assert_data_iter_mut(
            relations.iter_mut().map(|(_, _, data)| data),
            count,
            prefix,
        )?;
        for index in 0..count {
            prop_assert_eq!(
                relations.data(RelationId(index as u32)),
                &TestData(index + count),
            );
        }
    }

    #[test]
    fn test_var_relation_set_iterators_exact_size(
        count in 0usize..16,
        prefix in any::<usize>(),
    ) {
        let entries = (0..count)
            .map(|index| {
                (
                    vec![NodeId((2 * index) as u32), NodeId((2 * index + 1) as u32)],
                    TestData(index),
                )
            })
            .collect::<Vec<_>>();
        let mut relations = VarRelationSet::<NodeId, TestData>::new(entries);
        let prefix = prefix.min(count);

        assert_relation_ids(relations.ids(), count, prefix)?;
        assert_data_iter_mut(
            relations.iter_mut().map(|(_, _, data)| data),
            count,
            prefix,
        )?;
        for index in 0..count {
            prop_assert_eq!(
                relations.data(RelationId(index as u32)),
                &TestData(index + count),
            );
        }
    }

    #[test]
    fn test_fixed_fixed_birelation_set_iterators_exact_size(
        count in 0usize..16,
        prefix in any::<usize>(),
    ) {
        let entries = (0..count)
            .map(|index| {
                (
                    [NodeId(index as u32)],
                    [EdgeId(index as u32)],
                    TestData(index),
                )
            })
            .collect::<Vec<_>>();
        let mut relations = FixedFixedBirelationSet::<NodeId, 1, EdgeId, 1, TestData, >::new(entries);
        let prefix = prefix.min(count);

        assert_relation_ids(relations.ids(), count, prefix)?;
        assert_data_iter_mut(
            relations.iter_mut().map(|(_, _, _, data)| data),
            count,
            prefix,
        )?;
        for index in 0..count {
            prop_assert_eq!(
                relations.data(RelationId(index as u32)),
                &TestData(index + count),
            );
        }
    }

    #[test]
    fn test_fixed_var_birelation_set_iterators_exact_size(
        count in 0usize..16,
        prefix in any::<usize>(),
    ) {
        let entries = (0..count)
            .map(|index| {
                (
                    [NodeId(index as u32)],
                    vec![EdgeId(index as u32)],
                    TestData(index),
                )
            })
            .collect::<Vec<_>>();
        let mut relations = FixedVarBirelationSet::<NodeId, 1, EdgeId, TestData, >::new(entries);
        let prefix = prefix.min(count);

        assert_relation_ids(relations.ids(), count, prefix)?;
        assert_data_iter_mut(
            relations.iter_mut().map(|(_, _, _, data)| data),
            count,
            prefix,
        )?;
        for index in 0..count {
            prop_assert_eq!(
                relations.data(RelationId(index as u32)),
                &TestData(index + count),
            );
        }
    }

    #[test]
    fn test_var_var_birelation_set_iterators_exact_size(
        count in 0usize..16,
        prefix in any::<usize>(),
    ) {
        let entries = (0..count)
            .map(|index| {
                (
                    vec![NodeId(index as u32)],
                    vec![EdgeId(index as u32)],
                    TestData(index),
                )
            })
            .collect::<Vec<_>>();
        let mut relations = VarVarBirelationSet::<NodeId, EdgeId, TestData, >::new(entries);
        let prefix = prefix.min(count);

        assert_relation_ids(relations.ids(), count, prefix)?;
        assert_data_iter_mut(
            relations.iter_mut().map(|(_, _, _, data)| data),
            count,
            prefix,
        )?;
        for index in 0..count {
            prop_assert_eq!(
                relations.data(RelationId(index as u32)),
                &TestData(index + count),
            );
        }
    }
}
