//! Count-bearing, order-preserving removal from dense id spaces.

use std::error::Error;
use std::fmt::{self, Debug, Display, Formatter};
use std::ops::{Add, Sub};

use crate::graph::{EdgeId, NodeId};

/// Failure to construct a compaction over a finite source domain.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CompactionError<Id> {
    RemovedIdOutOfRange { id: Id, source_count: usize },
}

impl<Id: Debug> Display for CompactionError<Id> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::RemovedIdOutOfRange { id, source_count } => write!(
                f,
                "removed id {id:?} is out of range for {source_count} entries"
            ),
        }
    }
}

impl<Id: Debug> Error for CompactionError<Id> {}

/// Order-preserving removal between two finite dense id spaces.
///
/// Stores the source count and sorted, distinct removed ids. The result count is the source
/// count minus the number removed. Compatibility with a supplied object's counts is contextual.
///
/// # Semantic properties
///
/// On survivors, `compact` is strictly monotonic and `uncompact` is its inverse.
/// Identity preserves every id in its declared domain.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Compaction<Id> {
    source_count: usize,
    removed: Vec<Id>,
}

impl<Id> Compaction<Id>
where
    Id: Copy + Ord + Into<usize> + Add<usize, Output = Id> + Sub<usize, Output = Id>,
{
    /// Construct a compaction; input order and repeated removed ids do not matter.
    ///
    /// # Errors
    ///
    /// Returns an error for the first supplied removed id outside `0..source_count`.
    pub fn new(source_count: usize, mut removed: Vec<Id>) -> Result<Self, CompactionError<Id>> {
        for &id in &removed {
            if id.into() >= source_count {
                return Err(CompactionError::RemovedIdOutOfRange { id, source_count });
            }
        }
        removed.sort_unstable();
        removed.dedup();
        Ok(Self {
            source_count,
            removed,
        })
    }

    /// Identity on the declared source domain.
    pub fn identity(source_count: usize) -> Self {
        Self {
            source_count,
            removed: Vec::new(),
        }
    }

    /// Number of ids before removal.
    pub fn source_count(&self) -> usize {
        self.source_count
    }

    /// Number of surviving ids.
    pub fn result_count(&self) -> usize {
        self.source_count - self.removed.len()
    }

    /// The survivor's result id, or `None` for a removed or out-of-source-range id.
    pub fn compact(&self, old: Id) -> Option<Id> {
        if old.into() >= self.source_count || self.removed.binary_search(&old).is_ok() {
            return None;
        }
        Some(old - self.removed.partition_point(|&r| r < old))
    }

    /// Recover a survivor's source id.
    ///
    /// # Panics
    ///
    /// Panics when `post` is outside the result domain.
    pub fn uncompact(&self, post: Id) -> Id {
        self.try_uncompact(post)
            .expect("id outside compaction result domain")
    }

    /// Checked form of [`Self::uncompact`]; returns `None` outside the result domain.
    pub fn try_uncompact(&self, post: Id) -> Option<Id> {
        if post.into() >= self.result_count() {
            return None;
        }
        let mut old = post;
        loop {
            let next = post + self.removed.partition_point(|&r| r <= old);
            if next == old {
                return Some(old);
            }
            old = next;
        }
    }

    /// Compact a source-id-indexed data column, preserving survivor order.
    ///
    /// # Panics
    ///
    /// Panics when the column length differs from `source_count`.
    pub fn compact_vec<T: Clone>(&self, data: &[T]) -> Vec<T> {
        self.try_compact_vec(data)
            .expect("value count differs from compaction source count")
    }

    /// Checked form of [`Self::compact_vec`]; returns `None` on a source-count mismatch.
    pub fn try_compact_vec<T: Clone>(&self, data: &[T]) -> Option<Vec<T>> {
        if data.len() != self.source_count {
            return None;
        }
        Some(
            data.iter()
                .enumerate()
                .filter(|(index, _)| {
                    self.removed
                        .binary_search_by_key(index, |&id| id.into())
                        .is_err()
                })
                .map(|(_, value)| value.clone())
                .collect(),
        )
    }

    /// Removed ids, ascending and deduplicated.
    pub fn removed(&self) -> &[Id] {
        &self.removed
    }
}

/// Independent compactions of a graph's node and edge spaces.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GraphCompaction {
    nodes: Compaction<NodeId>,
    edges: Compaction<EdgeId>,
}

impl GraphCompaction {
    /// Assemble validated node and edge compactions.
    pub fn new(nodes: Compaction<NodeId>, edges: Compaction<EdgeId>) -> Self {
        Self { nodes, edges }
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
    /// Recover a node's source id; panics outside the result node domain.
    pub fn uncompact_node(&self, post: NodeId) -> NodeId {
        self.nodes.uncompact(post)
    }
    /// Checked form of [`Self::uncompact_node`].
    pub fn try_uncompact_node(&self, post: NodeId) -> Option<NodeId> {
        self.nodes.try_uncompact(post)
    }
    /// Recover an edge's source id; panics outside the result edge domain.
    pub fn uncompact_edge(&self, post: EdgeId) -> EdgeId {
        self.edges.uncompact(post)
    }
    /// Checked form of [`Self::uncompact_edge`].
    pub fn try_uncompact_edge(&self, post: EdgeId) -> Option<EdgeId> {
        self.edges.try_uncompact(post)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case(CompactionError::RemovedIdOutOfRange { id: NodeId(2), source_count: 2 }, "removed id NodeId(2) is out of range for 2 entries")]
    fn test_compaction_error_display(
        #[case] error: CompactionError<NodeId>,
        #[case] expected: &str,
    ) {
        assert_eq!(error.to_string(), expected);
    }

    #[rstest]
    #[case::empty(0, vec![], vec![])]
    #[case::full(3, vec![NodeId(2), NodeId(0), NodeId(1)], vec![NodeId(0), NodeId(1), NodeId(2)])]
    #[case::duplicates(4, vec![NodeId(2), NodeId(0), NodeId(2)], vec![NodeId(0), NodeId(2)])]
    fn test_compaction_new(
        #[case] source_count: usize,
        #[case] removed: Vec<NodeId>,
        #[case] expected: Vec<NodeId>,
    ) {
        assert_eq!(
            Compaction::new(source_count, removed),
            Ok(Compaction {
                source_count,
                removed: expected
            })
        );
    }

    #[rstest]
    #[case::empty_domain(0, vec![NodeId(0)], NodeId(0))]
    #[case::boundary(3, vec![NodeId(3)], NodeId(3))]
    #[case::outside_column(4, vec![NodeId(1), NodeId(8)], NodeId(8))]
    #[case::first_invalid(3, vec![NodeId(5), NodeId(3)], NodeId(5))]
    fn test_compaction_new_error(
        #[case] source_count: usize,
        #[case] removed: Vec<NodeId>,
        #[case] id: NodeId,
    ) {
        assert_eq!(
            Compaction::new(source_count, removed),
            Err(CompactionError::RemovedIdOutOfRange { id, source_count })
        );
    }

    #[rstest]
    #[case(0)]
    #[case(3)]
    fn test_compaction_identity(#[case] source_count: usize) {
        assert_eq!(
            Compaction::<NodeId>::identity(source_count),
            Compaction {
                source_count,
                removed: vec![]
            }
        );
    }

    #[rstest]
    #[case::empty(0, vec![], 0)]
    #[case::identity(3, vec![], 3)]
    #[case::partial(4, vec![NodeId(1), NodeId(1)], 3)]
    #[case::full(2, vec![NodeId(0), NodeId(1)], 0)]
    fn test_compaction_result_count(
        #[case] count: usize,
        #[case] removed: Vec<NodeId>,
        #[case] expected: usize,
    ) {
        let compaction = Compaction::new(count, removed).unwrap();
        assert_eq!(compaction.source_count(), count);
        assert_eq!(compaction.result_count(), expected);
    }

    #[rstest]
    #[case::empty(vec![])]
    #[case::front(vec![NodeId(0)])]
    #[case::middle(vec![NodeId(3)])]
    #[case::scattered(vec![NodeId(1), NodeId(4), NodeId(5)])]
    fn test_compaction_compact_roundtrip(#[case] removed: Vec<NodeId>) {
        let compaction = Compaction::new(8, removed.clone()).unwrap();
        for old in (0..8).map(NodeId) {
            if let Some(post) = compaction.compact(old) {
                assert_eq!(compaction.uncompact(post), old);
            } else {
                assert!(removed.contains(&old));
            }
        }
    }

    #[rstest]
    #[case::empty(0, vec![], NodeId(0), None)]
    #[case::boundary(3, vec![], NodeId(3), None)]
    #[case::fully_removed(2, vec![NodeId(0), NodeId(1)], NodeId(1), None)]
    fn test_compaction_compact_domain(
        #[case] count: usize,
        #[case] removed: Vec<NodeId>,
        #[case] id: NodeId,
        #[case] expected: Option<NodeId>,
    ) {
        assert_eq!(
            Compaction::new(count, removed).unwrap().compact(id),
            expected
        );
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
        assert_eq!(Compaction::new(10, removed).unwrap().compact(old), expected);
    }

    #[rstest]
    #[case::short(vec![10])]
    #[case::long(vec![10, 20, 30])]
    #[should_panic(expected = "value count differs from compaction source count")]
    fn test_compaction_compact_vec_error(#[case] data: Vec<i32>) {
        Compaction::<NodeId>::identity(2).compact_vec(&data);
    }

    #[rstest]
    #[case::empty(vec![])]
    fn test_compaction_compact_vec_identity(#[case] removed: Vec<NodeId>) {
        let data = vec![10, 20, 30, 40];
        assert_eq!(
            Compaction::new(4, removed).unwrap().compact_vec(&data),
            data
        );
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
        assert_eq!(
            Compaction::new(10, removed).unwrap().uncompact(post),
            expected
        );
    }

    #[rstest]
    #[case::boundary(NodeId(2))]
    #[case::beyond(NodeId(3))]
    #[should_panic(expected = "id outside compaction result domain")]
    fn test_compaction_uncompact_error(#[case] id: NodeId) {
        Compaction::new(3, vec![NodeId(1)]).unwrap().uncompact(id);
    }

    #[rstest]
    #[case::empty(0, vec![], NodeId(0), None)]
    #[case::survivor(3, vec![NodeId(1)], NodeId(1), Some(NodeId(2)))]
    #[case::boundary(3, vec![NodeId(1)], NodeId(2), None)]
    #[case::full(2, vec![NodeId(0), NodeId(1)], NodeId(0), None)]
    fn test_compaction_try_uncompact(
        #[case] count: usize,
        #[case] removed: Vec<NodeId>,
        #[case] id: NodeId,
        #[case] expected: Option<NodeId>,
    ) {
        assert_eq!(
            Compaction::new(count, removed).unwrap().try_uncompact(id),
            expected
        );
    }

    #[rstest]
    #[case::first(vec![NodeId(0)], vec![20, 30, 40])]
    #[case::middle(vec![NodeId(1)], vec![10, 30, 40])]
    #[case::scattered(vec![NodeId(0), NodeId(2)], vec![20, 40])]
    fn test_compaction_compact_vec(#[case] removed: Vec<NodeId>, #[case] expected: Vec<i32>) {
        assert_eq!(
            Compaction::new(4, removed)
                .unwrap()
                .compact_vec(&[10, 20, 30, 40]),
            expected
        );
    }

    #[rstest]
    #[case::empty(0, vec![], vec![], Some(vec![]))]
    #[case::full(2, vec![NodeId(0), NodeId(1)], vec![10, 20], Some(vec![]))]
    #[case::survivors(3, vec![NodeId(1)], vec![10, 20, 30], Some(vec![10, 30]))]
    #[case::short(3, vec![NodeId(1)], vec![10, 20], None)]
    #[case::long(1, vec![], vec![10, 20], None)]
    fn test_compaction_try_compact_vec(
        #[case] count: usize,
        #[case] removed: Vec<NodeId>,
        #[case] data: Vec<i32>,
        #[case] expected: Option<Vec<i32>>,
    ) {
        assert_eq!(
            Compaction::new(count, removed)
                .unwrap()
                .try_compact_vec(&data),
            expected
        );
    }

    #[rstest]
    #[case::unsorted(vec![NodeId(2), NodeId(0)], vec![NodeId(0), NodeId(2)])]
    #[case::duplicated(vec![NodeId(1), NodeId(1)], vec![NodeId(1)])]
    #[case::empty(vec![], vec![])]
    fn test_compaction_removed(#[case] removed: Vec<NodeId>, #[case] expected: Vec<NodeId>) {
        assert_eq!(Compaction::new(10, removed).unwrap().removed(), expected);
    }
    #[rstest]
    #[case::empty(NodeId(0), vec![], Some(NodeId(0)))]
    #[case::before_removed(NodeId(0), vec![NodeId(2)], Some(NodeId(0)))]
    #[case::removed(NodeId(2), vec![NodeId(2)], None)]
    #[case::after_removed(NodeId(3), vec![NodeId(2)], Some(NodeId(2)))]
    #[case::multi_removed(NodeId(5), vec![NodeId(1), NodeId(3)], Some(NodeId(3)))]
    fn test_graph_compaction_compact_node(
        #[case] old: NodeId,
        #[case] removed: Vec<NodeId>,
        #[case] expected: Option<NodeId>,
    ) {
        let compaction = GraphCompaction::new(
            Compaction::new(6, removed).unwrap(),
            Compaction::identity(0),
        );
        assert_eq!(compaction.compact_node(old), expected);
    }

    #[rstest]
    #[case::empty(NodeId(0), vec![], NodeId(0))]
    #[case::before_gap(NodeId(0), vec![NodeId(2)], NodeId(0))]
    #[case::at_gap(NodeId(2), vec![NodeId(2)], NodeId(3))]
    #[case::after_gap(NodeId(3), vec![NodeId(2)], NodeId(4))]
    #[case::multi_removed(NodeId(3), vec![NodeId(1), NodeId(3)], NodeId(5))]
    fn test_graph_compaction_uncompact_node(
        #[case] post: NodeId,
        #[case] removed: Vec<NodeId>,
        #[case] expected: NodeId,
    ) {
        let compaction = GraphCompaction::new(
            Compaction::new(6, removed).unwrap(),
            Compaction::identity(0),
        );
        assert_eq!(compaction.uncompact_node(post), expected);
    }
}
