//! A partial bijection between two id spaces.
//!
//! Framed as a matching in the bipartite graph over the two id sets: an id is **matched** (paired
//! with a partner on the other side — the shared interface) or **unmatched** (present on only one
//! side). Generic in the id type, so the same carrier serves node correspondences (atoms) and every
//! entity family (bonds, overlays) one layer up; `Correspondence<NodeId>` additionally exposes the
//! induced edge correspondence over the two graphs.

use crate::graph::{EdgeId, Graph, NodeId, Remapping};

/// A partial bijection between two `Id` spaces: the matched `(left, right)` pairs; every unmatched
/// id is reported on its side. Only the matched pairs are stored — unmatched ids are derived on
/// demand, so the carrier stays cheap to produce on the hot enumeration path.
///
/// Invariant: `matched_pairs` is sorted by left id (no duplicate lefts — it is a bijection), so
/// `right_of` is a binary search and `left_unmatched` a single merge. Producers already emit this
/// order (a subgraph-isomorphism match is query-index order; a maximum-common-subgraph is sorted by
/// the first graph's node), so `new` only confirms it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Correspondence<Id> {
    matched_pairs: Vec<(Id, Id)>,
    left_count: usize,
    right_count: usize,
}

impl<Id: Copy + Ord + From<usize>> Correspondence<Id> {
    /// A correspondence from its matched `(left, right)` pairs over two id spaces of the given sizes.
    /// Every id of either side not appearing in `matched_pairs` is unmatched on that side. Pairs are
    /// sorted by left id to establish the lookup invariant (cheap when already sorted).
    pub fn new(mut matched_pairs: Vec<(Id, Id)>, left_count: usize, right_count: usize) -> Self {
        matched_pairs.sort_unstable_by_key(|&(left, _)| left);
        Self {
            matched_pairs,
            left_count,
            right_count,
        }
    }

    /// A correspondence whose left space is dense (`0..images.len()`), pairing each left id `i` with
    /// `images[i]` — the shape a subgraph match / embedding induces (query index → host id). Already
    /// sorted by left, so no sort. Injectivity is a property of `images`, not enforced (as in `new`).
    pub fn from_images(images: &[Id], right_count: usize) -> Self {
        Self {
            matched_pairs: images
                .iter()
                .enumerate()
                .map(|(left, &right)| (Id::from(left), right))
                .collect(),
            left_count: images.len(),
            right_count,
        }
    }

    /// The matched `(left, right)` pairs (sorted by left) — the shared interface.
    pub fn matched_pairs(&self) -> &[(Id, Id)] {
        &self.matched_pairs
    }

    /// The number of matched pairs (the interface size).
    pub fn matched_pair_count(&self) -> usize {
        self.matched_pairs.len()
    }

    /// The size of the left id space.
    pub fn left_count(&self) -> usize {
        self.left_count
    }

    /// The size of the right id space.
    pub fn right_count(&self) -> usize {
        self.right_count
    }

    /// Whether every id on both sides is matched — a total bijection with no unmatched ids. A diff
    /// through such a correspondence adds and removes nothing.
    pub fn is_total(&self) -> bool {
        self.matched_pairs.len() == self.left_count && self.matched_pairs.len() == self.right_count
    }

    /// The right partner of a left id, if matched. Binary search (pairs sorted by left).
    pub fn right_of(&self, left: Id) -> Option<Id> {
        self.matched_pairs
            .binary_search_by_key(&left, |&(l, _)| l)
            .ok()
            .map(|index| self.matched_pairs[index].1)
    }

    /// The left partner of a right id, if matched. Linear — the un-indexed reverse direction.
    pub fn left_of(&self, right: Id) -> Option<Id> {
        self.matched_pairs
            .iter()
            .find(|&&(_, r)| r == right)
            .map(|&(l, _)| l)
    }

    /// Left ids with no partner — the ones deleted when read as a transformation. Merge over the
    /// sorted left column.
    pub fn left_unmatched(&self) -> Vec<Id> {
        unmatched(
            self.left_count,
            self.matched_pairs.iter().map(|&(left, _)| left),
        )
    }

    /// Right ids with no partner — the ones created when read as a transformation. Sorts the right
    /// column (unindexed) once, then merges.
    pub fn right_unmatched(&self) -> Vec<Id> {
        let mut rights: Vec<Id> = self.matched_pairs.iter().map(|&(_, right)| right).collect();
        rights.sort_unstable();
        unmatched(self.right_count, rights.into_iter())
    }

    /// Relational composition: `self` (left↔middle) followed by `other` (middle↔right), yielding a
    /// left↔right correspondence. A left id matched to a middle id that `other` leaves unmatched
    /// becomes unmatched. `self`'s right space and `other`'s left space must be the same.
    pub fn compose(&self, other: &Correspondence<Id>) -> Correspondence<Id> {
        let matched_pairs = self
            .matched_pairs
            .iter()
            .filter_map(|&(left, middle)| other.right_of(middle).map(|right| (left, right)))
            .collect();
        Correspondence::new(matched_pairs, self.left_count, other.right_count)
    }

    /// The inverse correspondence (right↔left): each matched pair swapped and the two id-space sizes
    /// exchanged. A left-unmatched id becomes right-unmatched and vice versa.
    pub fn reverse(&self) -> Correspondence<Id> {
        let matched_pairs = self
            .matched_pairs
            .iter()
            .map(|&(left, right)| (right, left))
            .collect();
        Correspondence::new(matched_pairs, self.right_count, self.left_count)
    }
}

impl Correspondence<NodeId> {
    /// The induced edge correspondence: `(left_edge, right_edge)` pairs whose endpoints are matched
    /// to an edge on the other side.
    pub fn edge_matched_pairs(&self, left: &Graph, right: &Graph) -> Vec<(EdgeId, EdgeId)> {
        left.edge_ids()
            .filter_map(|left_edge| {
                let [u, v] = left.edge_endpoints(left_edge);
                let right_edge = right.find_edge(self.right_of(u)?, self.right_of(v)?)?;
                Some((left_edge, right_edge))
            })
            .collect()
    }

    /// The number of edges shared under the node matching (the maximum-common-edge-subgraph
    /// objective).
    pub fn shared_edge_count(&self, left: &Graph, right: &Graph) -> usize {
        self.edge_matched_pairs(left, right).len()
    }
}

/// A subgraph↔host correspondence over a `Graph`: its node and edge families. The graph-core base
/// that the molecule-level `MoleculeCorrespondence` (atoms + bonds + overlays) extends — produced by
/// induced subgraphs, subiso matches, and common-subgraph search. The objective of each is a family
/// size: `nodes().matched_pair_count()` (induced / MCIS), `edges().matched_pair_count()` (MCES).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GraphCorrespondence {
    nodes: Correspondence<NodeId>,
    edges: Correspondence<EdgeId>,
}

impl GraphCorrespondence {
    pub fn new(nodes: Correspondence<NodeId>, edges: Correspondence<EdgeId>) -> Self {
        Self { nodes, edges }
    }

    /// The graph correspondence a node correspondence induces over `left` / `right`: its edge family
    /// is the induced edge correspondence (an edge matched when both endpoints are matched). Exact
    /// for a subiso match or a common *induced* subgraph — the two cases where every structurally
    /// matched edge is admissible; not for an edge-subgraph result under a nontrivial edge predicate
    /// (there the producer supplies the edge family directly).
    pub fn induced(left: &Graph, right: &Graph, nodes: Correspondence<NodeId>) -> Self {
        let edges = Correspondence::new(
            nodes.edge_matched_pairs(left, right),
            left.edge_count(),
            right.edge_count(),
        );
        Self { nodes, edges }
    }

    pub fn nodes(&self) -> &Correspondence<NodeId> {
        &self.nodes
    }

    pub fn edges(&self) -> &Correspondence<EdgeId> {
        &self.edges
    }

    /// This correspondence as a [`Remapping`] — a dense old→new relabel of both id spaces. Requires
    /// it be **total on the left** (every left id matched), as a pushout's coprojection is: left id `i`
    /// maps to its partner.
    pub fn to_remapping(&self) -> Remapping {
        Remapping::new(dense_images(&self.nodes), dense_images(&self.edges))
    }
}

/// The image column of a total-on-left correspondence: `out[i]` is the partner of left id `i` (the
/// matched pairs are sorted by left, and totality fills `0..left_count`).
fn dense_images<Id: Copy + Ord + From<usize>>(correspondence: &Correspondence<Id>) -> Vec<Id> {
    debug_assert_eq!(
        correspondence.matched_pair_count(),
        correspondence.left_count(),
        "to_remapping requires a total-on-left correspondence",
    );
    correspondence
        .matched_pairs()
        .iter()
        .map(|&(_, right)| right)
        .collect()
}

/// The ids `0..count` absent from `sorted_matched` (which must be ascending, no duplicates) — a
/// single merge pass, no per-id search.
fn unmatched<Id: Copy + Ord + From<usize>>(
    count: usize,
    sorted_matched: impl Iterator<Item = Id>,
) -> Vec<Id> {
    let mut matched = sorted_matched.peekable();
    (0..count)
        .map(Id::from)
        .filter(|node| {
            if matched.peek() == Some(node) {
                matched.next();
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
    fn paths() -> (Graph, Graph, Correspondence<NodeId>) {
        // left path 0-1-2 mapped onto the interior of right path 0-1-2-3; right node 0 unmatched.
        (
            Graph::new(3, &[[0, 1], [1, 2]]),
            Graph::new(4, &[[0, 1], [1, 2], [2, 3]]),
            Correspondence::new(vec![(n(0), n(1)), (n(1), n(2)), (n(2), n(3))], 3, 4),
        )
    }

    #[rstest]
    fn test_correspondence_matched_pairs() {
        let c = Correspondence::new(vec![(n(0), n(2)), (n(1), n(3))], 3, 4);
        assert_eq!(c.matched_pairs(), &[(n(0), n(2)), (n(1), n(3))]);
        assert_eq!(c.matched_pair_count(), 2);
    }

    #[rstest]
    fn test_correspondence_new_sorts() {
        let c = Correspondence::new(vec![(n(2), n(0)), (n(0), n(3)), (n(1), n(1))], 3, 4);
        assert_eq!(
            c.matched_pairs(),
            &[(n(0), n(3)), (n(1), n(1)), (n(2), n(0))]
        );
        assert_eq!(c.right_of(n(2)), Some(n(0)));
    }

    #[rstest]
    fn test_correspondence_from_images() {
        // dense left 0,1,2 → images 3,1,0; every left matched, host id 2 right-unmatched.
        let c = Correspondence::from_images(&[n(3), n(1), n(0)], 4);
        assert_eq!(
            c.matched_pairs(),
            &[(n(0), n(3)), (n(1), n(1)), (n(2), n(0))]
        );
        assert_eq!(c.matched_pair_count(), 3);
        assert_eq!(c.left_unmatched(), Vec::<NodeId>::new());
        assert_eq!(c.right_unmatched(), vec![n(2)]);
    }

    #[rstest]
    fn test_graph_correspondence() {
        let c = GraphCorrespondence::new(
            Correspondence::from_images(&[n(1), n(0)], 3),
            Correspondence::from_images(&[e(2)], 4),
        );
        assert_eq!(c.nodes().matched_pairs(), &[(n(0), n(1)), (n(1), n(0))]);
        assert_eq!(c.edges().matched_pairs(), &[(e(0), e(2))]);
    }

    #[rstest]
    fn test_graph_correspondence_to_remapping() {
        // total-on-left node map 0→2, 1→0, 2→1 and edge map 0→1, 1→0.
        let c = GraphCorrespondence::new(
            Correspondence::from_images(&[n(2), n(0), n(1)], 3),
            Correspondence::from_images(&[e(1), e(0)], 2),
        );
        let remapping = c.to_remapping();
        assert_eq!(remapping.map_node(n(0)), n(2));
        assert_eq!(remapping.map_node(n(1)), n(0));
        assert_eq!(remapping.map_node(n(2)), n(1));
        assert_eq!(remapping.map_edge(e(0)), e(1));
        assert_eq!(remapping.map_edge(e(1)), e(0));
    }

    #[rstest]
    #[case::matched_first(n(0), Some(n(2)))]
    #[case::matched_second(n(1), Some(n(3)))]
    #[case::unmatched(n(2), None)]
    fn test_correspondence_right_of(#[case] left: NodeId, #[case] expected: Option<NodeId>) {
        let c = Correspondence::new(vec![(n(0), n(2)), (n(1), n(3))], 3, 4);
        assert_eq!(c.right_of(left), expected);
    }

    #[rstest]
    #[case::matched_first(n(2), Some(n(0)))]
    #[case::matched_second(n(3), Some(n(1)))]
    #[case::unmatched(n(0), None)]
    fn test_correspondence_left_of(#[case] right: NodeId, #[case] expected: Option<NodeId>) {
        let c = Correspondence::new(vec![(n(0), n(2)), (n(1), n(3))], 3, 4);
        assert_eq!(c.left_of(right), expected);
    }

    #[rstest]
    fn test_correspondence_left_unmatched() {
        let c = Correspondence::new(vec![(n(0), n(2)), (n(1), n(3))], 3, 4);
        assert_eq!(c.left_unmatched(), vec![n(2)]);
    }

    #[rstest]
    fn test_correspondence_right_unmatched() {
        let c = Correspondence::new(vec![(n(0), n(2)), (n(1), n(3))], 3, 4);
        assert_eq!(c.right_unmatched(), vec![n(0), n(1)]);
    }

    #[rstest]
    fn test_correspondence_edge_matched_pairs(paths: (Graph, Graph, Correspondence<NodeId>)) {
        let (left, right, c) = paths;
        assert_eq!(
            c.edge_matched_pairs(&left, &right),
            vec![(e(0), e(1)), (e(1), e(2))]
        );
    }

    #[rstest]
    fn test_correspondence_shared_edge_count(paths: (Graph, Graph, Correspondence<NodeId>)) {
        let (left, right, c) = paths;
        assert_eq!(c.shared_edge_count(&left, &right), 2);
    }

    #[rstest]
    #[case::total(vec![(n(0), n(0)), (n(1), n(1))], 2, 2, true)]
    #[case::left_unmatched(vec![(n(0), n(0))], 2, 1, false)]
    #[case::right_unmatched(vec![(n(0), n(0))], 1, 2, false)]
    fn test_correspondence_is_total(
        #[case] matched_pairs: Vec<(NodeId, NodeId)>,
        #[case] left_count: usize,
        #[case] right_count: usize,
        #[case] expected: bool,
    ) {
        assert_eq!(
            Correspondence::new(matched_pairs, left_count, right_count).is_total(),
            expected
        );
    }

    #[rstest]
    fn test_correspondence_compose() {
        // A⇌B then B⇌C; A-node 2 maps to B-node 12, which B⇌C leaves unmatched, so 2 drops out.
        let ab = Correspondence::new(vec![(n(0), n(10)), (n(1), n(11)), (n(2), n(12))], 3, 13);
        let bc = Correspondence::new(vec![(n(10), n(100)), (n(11), n(101))], 13, 102);
        let ac = ab.compose(&bc);
        assert_eq!(ac.matched_pairs(), &[(n(0), n(100)), (n(1), n(101))]);
        assert_eq!(ac.left_unmatched(), vec![n(2)]);
    }

    #[rstest]
    fn test_correspondence_reverse() {
        // pairs and counts swap; the left-unmatched id 2 becomes right-unmatched.
        let c = Correspondence::new(vec![(n(0), n(3)), (n(1), n(1))], 3, 4);
        let reversed = c.reverse();
        assert_eq!(reversed.matched_pairs(), &[(n(1), n(1)), (n(3), n(0))]);
        assert_eq!(reversed.left_unmatched(), vec![n(0), n(2)]);
        assert_eq!(reversed.right_unmatched(), vec![n(2)]);
    }
}
