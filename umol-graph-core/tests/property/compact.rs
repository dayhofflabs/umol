//! Finite compaction laws against survivor enumeration, independent of shift arithmetic.

use proptest::prelude::*;
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
