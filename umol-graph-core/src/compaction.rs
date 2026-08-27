//! Dense renumbering of an id space by removal.
//!
//! `Compaction<Id>` is the operation over one id space; `GraphCompaction` pairs one over nodes with
//! one over edges. A removal leaves surviving ids closed up in place, so a stale reference is either
//! shifted to its new position or reported as removed.

use std::ops::{Add, Sub};

use crate::graph::{EdgeId, NodeId};

/// Dense renumbering of one id space by removal: a surviving id maps to its position in the
/// closed-up post-removal table, a removed id has no image.
///
/// # Semantic properties
///
/// `compact` is strictly monotonic on survivors and `uncompact` inverts it, so
/// `uncompact(compact(id).unwrap()) == id` for every surviving `id`. An empty removal list is the
/// identity on both.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Compaction<Id> {
    removed: Vec<Id>,
}

impl<Id> Compaction<Id>
where
    Id: Copy + Ord + Add<usize, Output = Id> + Sub<usize, Output = Id>,
{
    /// The compaction removing `removed`; duplicates and input order do not matter.
    pub fn new(mut removed: Vec<Id>) -> Self {
        removed.sort_unstable();
        removed.dedup();
        Self { removed }
    }

    /// The post-removal image of `old`, or `None` when `old` was removed.
    pub fn compact(&self, old: Id) -> Option<Id> {
        if self.removed.binary_search(&old).is_ok() {
            return None;
        }
        Some(old - self.removed.partition_point(|&r| r < old))
    }

    /// The pre-removal id of a surviving post-removal id, re-adding every removed id at or below it.
    pub fn uncompact(&self, post: Id) -> Id {
        let mut old = post;
        loop {
            let next = post + self.removed.partition_point(|&r| r <= old);
            if next == old {
                return old;
            }
            old = next;
        }
    }

    /// The removed ids, ascending and deduplicated.
    pub fn removed(&self) -> &[Id] {
        &self.removed
    }
}

impl<Id> Default for Compaction<Id> {
    fn default() -> Self {
        Self {
            removed: Vec::new(),
        }
    }
}

/// Dense renumbering of a graph's node and edge spaces by one removal.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct GraphCompaction {
    nodes: Compaction<NodeId>,
    edges: Compaction<EdgeId>,
}

impl GraphCompaction {
    pub fn new(removed_nodes: Vec<NodeId>, removed_edges: Vec<EdgeId>) -> Self {
        Self {
            nodes: Compaction::new(removed_nodes),
            edges: Compaction::new(removed_edges),
        }
    }

    pub fn nodes(&self) -> &Compaction<NodeId> {
        &self.nodes
    }

    pub fn edges(&self) -> &Compaction<EdgeId> {
        &self.edges
    }

    pub fn compact_node(&self, old: NodeId) -> Option<NodeId> {
        self.nodes.compact(old)
    }

    pub fn compact_edge(&self, old: EdgeId) -> Option<EdgeId> {
        self.edges.compact(old)
    }

    pub fn uncompact_node(&self, post: NodeId) -> NodeId {
        self.nodes.uncompact(post)
    }

    pub fn uncompact_edge(&self, post: EdgeId) -> EdgeId {
        self.edges.uncompact(post)
    }
}

/// Compact a node-indexed data column to the post-removal layout (drop removed, keep order).
pub fn compact_node_vec<T: Clone>(compaction: &GraphCompaction, data: &[T]) -> Vec<T> {
    data.iter()
        .enumerate()
        .filter(|(i, _)| compaction.compact_node(NodeId(*i as u32)).is_some())
        .map(|(_, v)| v.clone())
        .collect()
}

/// Compact an edge-indexed data column to the post-removal layout (drop removed, keep order).
pub fn compact_edge_vec<T: Clone>(compaction: &GraphCompaction, data: &[T]) -> Vec<T> {
    data.iter()
        .enumerate()
        .filter(|(i, _)| compaction.compact_edge(EdgeId(*i as u32)).is_some())
        .map(|(_, v)| v.clone())
        .collect()
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::empty(NodeId(0), vec![], Some(NodeId(0)))]
    #[case::before_removed(NodeId(0), vec![NodeId(2)], Some(NodeId(0)))]
    #[case::removed(NodeId(2), vec![NodeId(2)], None)]
    #[case::after_removed(NodeId(3), vec![NodeId(2)], Some(NodeId(2)))]
    #[case::multi_removed(NodeId(5), vec![NodeId(1), NodeId(3)], Some(NodeId(3)))]
    fn test_compaction_node(
        #[case] old: NodeId,
        #[case] removed: Vec<NodeId>,
        #[case] expected: Option<NodeId>,
    ) {
        let compaction = GraphCompaction::new(removed, vec![]);
        assert_eq!(compaction.compact_node(old), expected);
    }

    #[rstest]
    #[case::empty(NodeId(0), vec![], NodeId(0))]
    #[case::before_gap(NodeId(0), vec![NodeId(2)], NodeId(0))]
    #[case::at_gap(NodeId(2), vec![NodeId(2)], NodeId(3))]
    #[case::after_gap(NodeId(3), vec![NodeId(2)], NodeId(4))]
    #[case::multi_removed(NodeId(3), vec![NodeId(1), NodeId(3)], NodeId(5))]
    fn test_uncompaction_node(
        #[case] post: NodeId,
        #[case] removed: Vec<NodeId>,
        #[case] expected: NodeId,
    ) {
        let compaction = GraphCompaction::new(removed, vec![]);
        assert_eq!(compaction.uncompact_node(post), expected);
    }

    #[rstest]
    #[case::empty(vec![], NodeId(3), Some(NodeId(3)))]
    #[case::front(vec![NodeId(0)], NodeId(3), Some(NodeId(2)))]
    #[case::middle(vec![NodeId(2)], NodeId(3), Some(NodeId(2)))]
    #[case::end(vec![NodeId(4)], NodeId(3), Some(NodeId(3)))]
    #[case::removed(vec![NodeId(3)], NodeId(3), None)]
    #[case::beyond_every_removal(vec![NodeId(0), NodeId(1)], NodeId(9), Some(NodeId(7)))]
    #[case::unsorted_input(vec![NodeId(3), NodeId(0)], NodeId(4), Some(NodeId(2)))]
    #[case::duplicate_input(vec![NodeId(0), NodeId(0)], NodeId(2), Some(NodeId(1)))]
    fn test_compaction_compact(
        #[case] removed: Vec<NodeId>,
        #[case] old: NodeId,
        #[case] expected: Option<NodeId>,
    ) {
        assert_eq!(Compaction::new(removed).compact(old), expected);
    }

    #[rstest]
    #[case::empty(vec![], NodeId(3), NodeId(3))]
    #[case::front(vec![NodeId(0)], NodeId(2), NodeId(3))]
    #[case::middle(vec![NodeId(2)], NodeId(2), NodeId(3))]
    #[case::end(vec![NodeId(4)], NodeId(3), NodeId(3))]
    #[case::consecutive(vec![NodeId(1), NodeId(2)], NodeId(1), NodeId(3))]
    fn test_compaction_uncompact(
        #[case] removed: Vec<NodeId>,
        #[case] post: NodeId,
        #[case] expected: NodeId,
    ) {
        assert_eq!(Compaction::new(removed).uncompact(post), expected);
    }

    #[rstest]
    #[case::empty(vec![])]
    #[case::front(vec![NodeId(0)])]
    #[case::middle(vec![NodeId(3)])]
    #[case::scattered(vec![NodeId(1), NodeId(4), NodeId(5)])]
    fn test_compaction_compact_roundtrip(#[case] removed: Vec<NodeId>) {
        let compaction = Compaction::new(removed.clone());
        for old in (0..8).map(NodeId) {
            if let Some(post) = compaction.compact(old) {
                assert_eq!(compaction.uncompact(post), old);
            } else {
                assert!(removed.contains(&old));
            }
        }
    }

    #[rstest]
    #[case::unsorted(vec![NodeId(2), NodeId(0)], vec![NodeId(0), NodeId(2)])]
    #[case::duplicated(vec![NodeId(1), NodeId(1)], vec![NodeId(1)])]
    #[case::empty(vec![], vec![])]
    fn test_compaction_removed(#[case] removed: Vec<NodeId>, #[case] expected: Vec<NodeId>) {
        assert_eq!(Compaction::new(removed).removed(), expected);
    }
}
