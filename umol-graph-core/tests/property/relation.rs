//! Exact-size properties of stored relation iterators.

use proptest::prelude::*;
use proptest::test_runner::TestCaseResult;
use umol_graph_core::{
    BiRelationData, EdgeId, FixedFixedBirelationSet, FixedRelationSet, FixedVarBirelationSet,
    NodeId, ParticipantPosition, RelationData, RelationId, Unordered, VarRelationSet,
    VarVarBirelationSet,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct TestData(usize);

impl RelationData for TestData {
    fn on_permutation(&mut self, _: &[ParticipantPosition]) {}
}

impl BiRelationData for TestData {
    fn on_permutation(&mut self, _: &[ParticipantPosition], _: &[ParticipantPosition]) {}
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
        let mut relations = FixedRelationSet::<NodeId, Unordered, TestData, 2>::new(entries);
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
        let mut relations = VarRelationSet::<NodeId, Unordered, TestData>::new(entries);
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
        let mut relations = FixedFixedBirelationSet::<
            NodeId,
            Unordered,
            1,
            EdgeId,
            Unordered,
            1,
            TestData,
        >::new(entries);
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
        let mut relations = FixedVarBirelationSet::<
            NodeId,
            Unordered,
            1,
            EdgeId,
            Unordered,
            TestData,
        >::new(entries);
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
        let mut relations = VarVarBirelationSet::<
            NodeId,
            Unordered,
            EdgeId,
            Unordered,
            TestData,
        >::new(entries);
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
