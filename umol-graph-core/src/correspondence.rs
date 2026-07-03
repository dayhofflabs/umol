//! A partial bijection between the node id spaces of two graphs.
//!
//! Framed as a matching in the bipartite graph over the two node sets: a node is **mated** (paired
//! with a partner on the other side — the shared interface) or **exposed** (unpaired, present on
//! only one side).

use crate::graph::{EdgeId, Graph, NodeId};

/// A partial bijection between two graphs' node id spaces: the **mated** `(left, right)` pairs;
/// every unmated node is **exposed** on its side. Only the mated pairs are stored — exposed nodes
/// and the induced edge correspondence are derived on demand, so the carrier stays cheap to produce
/// on the hot enumeration path.
///
/// Invariant: `mates` is sorted by left node (no duplicate lefts — it is a bijection), so `right_of`
/// is a binary search and `left_exposed` a single merge. Producers already emit this order (a
/// subgraph-isomorphism match is query-index order; a maximum-common-subgraph is sorted by the first
/// graph's node), so `new` only confirms it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Correspondence {
    mates: Vec<(NodeId, NodeId)>,
    left_count: usize,
    right_count: usize,
}

impl Correspondence {
    /// A correspondence from its mated `(left, right)` pairs over two graphs with the given node
    /// counts. Every node of either graph not appearing in `mates` is exposed on its side. The
    /// pairs are sorted by left node to establish the lookup invariant (cheap when already sorted).
    pub fn new(mut mates: Vec<(NodeId, NodeId)>, left_count: usize, right_count: usize) -> Self {
        mates.sort_unstable_by_key(|&(left, _)| left);
        Self {
            mates,
            left_count,
            right_count,
        }
    }

    /// The mated `(left, right)` pairs (sorted by left) — the shared interface.
    pub fn mates(&self) -> &[(NodeId, NodeId)] {
        &self.mates
    }

    /// The number of mated pairs (the interface size).
    pub fn node_count(&self) -> usize {
        self.mates.len()
    }

    /// The right partner of a left node, if mated. Binary search (mates sorted by left).
    pub fn right_of(&self, left: NodeId) -> Option<NodeId> {
        self.mates
            .binary_search_by_key(&left, |&(l, _)| l)
            .ok()
            .map(|index| self.mates[index].1)
    }

    /// The left partner of a right node, if mated. Linear — the un-indexed reverse direction.
    pub fn left_of(&self, right: NodeId) -> Option<NodeId> {
        self.mates
            .iter()
            .find(|&&(_, r)| r == right)
            .map(|&(l, _)| l)
    }

    /// Left nodes with no partner — the ones deleted when read as a transformation. Merge over the
    /// sorted left column.
    pub fn left_exposed(&self) -> Vec<NodeId> {
        exposed(self.left_count, self.mates.iter().map(|&(left, _)| left))
    }

    /// Right nodes with no partner — the ones created when read as a transformation. Sorts the right
    /// column (unindexed) once, then merges.
    pub fn right_exposed(&self) -> Vec<NodeId> {
        let mut rights: Vec<NodeId> = self.mates.iter().map(|&(_, right)| right).collect();
        rights.sort_unstable();
        exposed(self.right_count, rights.into_iter())
    }

    /// The induced edge correspondence: `(left_edge, right_edge)` pairs whose endpoints are mated
    /// to an edge on the other side.
    pub fn edge_mates(&self, left: &Graph, right: &Graph) -> Vec<(EdgeId, EdgeId)> {
        left.edge_ids()
            .filter_map(|left_edge| {
                let [u, v] = left.edge_endpoints(left_edge);
                let right_edge = right.find_edge(self.right_of(u)?, self.right_of(v)?)?;
                Some((left_edge, right_edge))
            })
            .collect()
    }

    /// The number of edges shared under the node mating (the maximum-common-edge-subgraph objective).
    pub fn shared_edge_count(&self, left: &Graph, right: &Graph) -> usize {
        self.edge_mates(left, right).len()
    }
}

/// The nodes `0..count` absent from `sorted_mated` (which must be ascending, no duplicates) — a
/// single merge pass, no per-node search.
fn exposed(count: usize, sorted_mated: impl Iterator<Item = NodeId>) -> Vec<NodeId> {
    let mut mated = sorted_mated.peekable();
    (0..count as u32)
        .map(NodeId)
        .filter(|&node| {
            if mated.peek() == Some(&node) {
                mated.next();
                false
            } else {
                true
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    fn n(i: u32) -> NodeId {
        NodeId(i)
    }

    fn e(i: u32) -> EdgeId {
        EdgeId(i)
    }

    #[fixture]
    fn paths() -> (Graph, Graph, Correspondence) {
        // left path 0-1-2 mapped onto the interior of right path 0-1-2-3; right node 0 exposed.
        (
            Graph::new(3, &[[0, 1], [1, 2]]),
            Graph::new(4, &[[0, 1], [1, 2], [2, 3]]),
            Correspondence::new(vec![(n(0), n(1)), (n(1), n(2)), (n(2), n(3))], 3, 4),
        )
    }

    #[rstest]
    fn test_correspondence_mates() {
        let c = Correspondence::new(vec![(n(0), n(2)), (n(1), n(3))], 3, 4);
        assert_eq!(c.mates(), &[(n(0), n(2)), (n(1), n(3))]);
        assert_eq!(c.node_count(), 2);
    }

    #[rstest]
    fn test_correspondence_new_sorts() {
        let c = Correspondence::new(vec![(n(2), n(0)), (n(0), n(3)), (n(1), n(1))], 3, 4);
        assert_eq!(c.mates(), &[(n(0), n(3)), (n(1), n(1)), (n(2), n(0))]);
        assert_eq!(c.right_of(n(2)), Some(n(0)));
    }

    #[rstest]
    #[case::mated_first(n(0), Some(n(2)))]
    #[case::mated_second(n(1), Some(n(3)))]
    #[case::unmated(n(2), None)]
    fn test_correspondence_right_of(#[case] left: NodeId, #[case] expected: Option<NodeId>) {
        let c = Correspondence::new(vec![(n(0), n(2)), (n(1), n(3))], 3, 4);
        assert_eq!(c.right_of(left), expected);
    }

    #[rstest]
    #[case::mated_first(n(2), Some(n(0)))]
    #[case::mated_second(n(3), Some(n(1)))]
    #[case::unmated(n(0), None)]
    fn test_correspondence_left_of(#[case] right: NodeId, #[case] expected: Option<NodeId>) {
        let c = Correspondence::new(vec![(n(0), n(2)), (n(1), n(3))], 3, 4);
        assert_eq!(c.left_of(right), expected);
    }

    #[rstest]
    fn test_correspondence_left_exposed() {
        let c = Correspondence::new(vec![(n(0), n(2)), (n(1), n(3))], 3, 4);
        assert_eq!(c.left_exposed(), vec![n(2)]);
    }

    #[rstest]
    fn test_correspondence_right_exposed() {
        let c = Correspondence::new(vec![(n(0), n(2)), (n(1), n(3))], 3, 4);
        assert_eq!(c.right_exposed(), vec![n(0), n(1)]);
    }

    #[rstest]
    fn test_correspondence_edge_mates(paths: (Graph, Graph, Correspondence)) {
        let (left, right, c) = paths;
        assert_eq!(
            c.edge_mates(&left, &right),
            vec![(e(0), e(1)), (e(1), e(2))]
        );
    }

    #[rstest]
    fn test_correspondence_shared_edge_count(paths: (Graph, Graph, Correspondence)) {
        let (left, right, c) = paths;
        assert_eq!(c.shared_edge_count(&left, &right), 2);
    }
}
