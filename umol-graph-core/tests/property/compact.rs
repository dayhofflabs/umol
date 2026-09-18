//! Finite compaction laws against survivor enumeration, independent of shift arithmetic.
//!
//! Generated columns exercise id/value compaction and correspondence conversion through size 64.
//! Exhaustive removal subsets through size 12 check both inverse lookups and their result-domain
//! boundary against the complete surviving source-id sequence.

use proptest::prelude::*;
use rstest::rstest;
use umol_graph_core::{Compaction, Correspondence, NodeId};

proptest! {
    #[test]
    fn test_compaction_compact(
        entries in prop::collection::vec((any::<bool>(), any::<i64>()), 0..65),
    ) {
        let removed = entries.iter().enumerate()
            .filter(|(_, (remove, _))| *remove)
            .map(|(idx, _)| NodeId::from(idx)).collect::<Vec<_>>();
        let survivors = entries.iter().enumerate()
            .filter(|(_, (remove, _))| !*remove)
            .map(|(idx, (_, value))| (NodeId::from(idx), *value)).collect::<Vec<_>>();
        let pairs = survivors.iter().enumerate()
            .map(|(idx, &(source, _))| (source, NodeId::from(idx))).collect::<Vec<_>>();
        let expected_values = survivors.iter().map(|&(_, value)| value).collect::<Vec<_>>();
        let values = entries.iter().map(|&(_, value)| value).collect::<Vec<_>>();
        let compaction = Compaction::new(entries.len(), removed).unwrap();

        prop_assert_eq!(compaction.source_count(), entries.len());
        prop_assert_eq!(compaction.result_count(), survivors.len());
        prop_assert_eq!(compaction.compact_vec(&values), expected_values.clone());
        prop_assert_eq!(compaction.try_compact_vec(&values), Some(expected_values));
        for (idx, &(remove, _)) in entries.iter().enumerate() {
            let source = NodeId::from(idx);
            let expected = survivors.iter().position(|&(id, _)| id == source).map(NodeId::from);
            prop_assert_eq!(compaction.compact(source), expected);
            prop_assert_eq!(expected.is_none(), remove);
        }
        for &(source, result) in &pairs {
            prop_assert_eq!(compaction.uncompact(result), source);
            prop_assert_eq!(compaction.try_uncompact(result), Some(source));
        }
        prop_assert_eq!(compaction.compact(NodeId::from(entries.len())), None);
        prop_assert_eq!(compaction.try_uncompact(NodeId::from(survivors.len())), None);
        prop_assert_eq!(Correspondence::from(&compaction),
            Correspondence::new(pairs, entries.len(), survivors.len()).unwrap());
    }
}

#[rstest]
fn test_compaction_uncompact() {
    for count in 0..=12 {
        for mask in 0..(1usize << count) {
            let removed = (0..count)
                .filter(|id| mask & (1 << id) != 0)
                .map(NodeId::from)
                .collect();
            let expected: Vec<_> = (0..count)
                .filter(|id| mask & (1 << id) == 0)
                .map(NodeId::from)
                .collect();
            let compaction = Compaction::new(count, removed).unwrap();
            let actual: Vec<_> = (0..expected.len())
                .map(|id| compaction.uncompact(NodeId::from(id)))
                .collect();
            assert_eq!(actual, expected);
            for (id, source) in expected
                .into_iter()
                .map(Some)
                .chain([None, None])
                .enumerate()
            {
                assert_eq!(compaction.try_uncompact(NodeId::from(id)), source);
            }
        }
    }
}
