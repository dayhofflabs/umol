//! Construction and algebraic properties of correspondences.

use std::collections::BTreeSet;
use std::fmt::Debug;

use proptest::prelude::*;
use umol_graph_core::{Correspondence, EdgeId, Graph, GraphCorrespondence, NodeId};

use super::strategy;

fn correspondence_strategy<Id>() -> impl Strategy<Value = Correspondence<Id>>
where
    Id: Copy + Debug + Ord + From<usize> + 'static,
{
    (0usize..8, 0usize..8).prop_flat_map(|(left_count, right_count)| {
        correspondence_with_counts_strategy(left_count, right_count)
    })
}

fn correspondence_with_counts_strategy<Id>(
    left_count: usize,
    right_count: usize,
) -> impl Strategy<Value = Correspondence<Id>>
where
    Id: Copy + Debug + Ord + From<usize> + 'static,
{
    let pair_count = left_count.min(right_count);
    (
        Just((0..left_count).map(Id::from).collect::<Vec<_>>()).prop_shuffle(),
        Just((0..right_count).map(Id::from).collect::<Vec<_>>()).prop_shuffle(),
        0usize..=pair_count,
    )
        .prop_map(move |(mut left, mut right, count)| {
            left.truncate(count);
            right.truncate(count);
            Correspondence::new(
                left.into_iter().zip(right).collect(),
                left_count,
                right_count,
            )
            .expect("correspondence producer preserves partial-bijection invariants")
        })
}

fn graph_context_strategy() -> impl Strategy<Value = (Graph, Graph, Correspondence<NodeId>)> {
    (strategy::multigraph(6, 10), strategy::multigraph(6, 10)).prop_flat_map(|(left, right)| {
        let correspondence =
            correspondence_with_counts_strategy(left.node_count(), right.node_count());
        (Just(left), Just(right), correspondence)
    })
}

fn correspondence_chain_strategy(
    length: usize,
) -> impl Strategy<Value = Vec<Correspondence<NodeId>>> {
    prop::collection::vec(0usize..8, length + 1).prop_flat_map(|counts| {
        counts
            .windows(2)
            .map(|pair| correspondence_with_counts_strategy::<NodeId>(pair[0], pair[1]).boxed())
            .collect::<Vec<_>>()
    })
}

fn graph_correspondence_strategy() -> impl Strategy<Value = GraphCorrespondence> {
    (
        correspondence_strategy::<NodeId>(),
        correspondence_strategy::<EdgeId>(),
    )
        .prop_map(|(nodes, edges)| GraphCorrespondence::new(nodes, edges))
}

fn reference_edge_matched_pairs(
    left: &Graph,
    right: &Graph,
    correspondence: &Correspondence<NodeId>,
) -> Option<Vec<(EdgeId, EdgeId)>> {
    let mut matched_pairs = Vec::new();
    for left_edge in left.edge_ids() {
        let [left_first, left_second] = left.edge_endpoints(left_edge);
        let (Some(right_first), Some(right_second)) = (
            correspondence.right_of(left_first),
            correspondence.right_of(left_second),
        ) else {
            continue;
        };
        let mut mapped_endpoints = [right_first, right_second];
        mapped_endpoints.sort_unstable();
        for right_edge in right.edge_ids() {
            if right.edge_endpoints(right_edge) == mapped_endpoints {
                matched_pairs.push((left_edge, right_edge));
            }
        }
    }

    let mut left_edges = BTreeSet::new();
    let mut right_edges = BTreeSet::new();
    if matched_pairs
        .iter()
        .any(|&(left, right)| !left_edges.insert(left) || !right_edges.insert(right))
    {
        None
    } else {
        Some(matched_pairs)
    }
}

fn correspondence_images_strategy() -> impl Strategy<Value = (Vec<NodeId>, usize)> {
    (0usize..8).prop_flat_map(|right_count| {
        (
            Just((0..right_count).map(NodeId::from).collect::<Vec<_>>()).prop_shuffle(),
            0usize..=right_count,
        )
            .prop_map(move |(mut images, left_count)| {
                images.truncate(left_count);
                (images, right_count)
            })
    })
}

proptest! {
    #[test]
    fn test_correspondence_from_images(
        (images, right_count) in correspondence_images_strategy(),
    ) {
        let correspondence = Correspondence::from_images(&images, right_count);
        let matched_pairs = correspondence.matched_pairs();
        let expected_pairs = images
            .iter()
            .enumerate()
            .map(|(left, &right)| (NodeId::from(left), right))
            .collect::<Vec<_>>();
        let left_matched = matched_pairs
            .iter()
            .map(|&(left, _)| left)
            .collect::<BTreeSet<_>>();
        let right_matched = matched_pairs
            .iter()
            .map(|&(_, right)| right)
            .collect::<BTreeSet<_>>();
        let mut left_partition = left_matched
            .iter()
            .copied()
            .chain(correspondence.left_unmatched())
            .collect::<Vec<_>>();
        let mut right_partition = right_matched
            .iter()
            .copied()
            .chain(correspondence.right_unmatched())
            .collect::<Vec<_>>();
        left_partition.sort_unstable();
        right_partition.sort_unstable();

        prop_assert_eq!(matched_pairs, expected_pairs);
        prop_assert!(matched_pairs.windows(2).all(|pair| pair[0].0 < pair[1].0));
        prop_assert_eq!(left_matched.len(), matched_pairs.len());
        prop_assert_eq!(right_matched.len(), matched_pairs.len());
        prop_assert!(matched_pairs
            .iter()
            .all(|&(left, right)| left < NodeId::from(images.len())
                && right < NodeId::from(right_count)));
        prop_assert_eq!(
            left_partition,
            (0..images.len()).map(NodeId::from).collect::<Vec<_>>(),
        );
        prop_assert_eq!(
            right_partition,
            (0..right_count).map(NodeId::from).collect::<Vec<_>>(),
        );
    }

    #[test]
    fn test_correspondence_is_total(
        correspondence in correspondence_strategy::<NodeId>(),
    ) {
        let total_on_left = correspondence.matched_pair_count() == correspondence.left_count();
        let total_on_right = correspondence.matched_pair_count() == correspondence.right_count();
        let reverse = correspondence.reverse();

        prop_assert_eq!(correspondence.is_total_on_left(), total_on_left);
        prop_assert_eq!(correspondence.is_total_on_right(), total_on_right);
        prop_assert_eq!(correspondence.is_total(), total_on_left && total_on_right);
        prop_assert_eq!(reverse.is_total_on_left(), total_on_right);
        prop_assert_eq!(reverse.is_total_on_right(), total_on_left);
    }

    #[test]
    fn test_correspondence_compose_associativity(
        chain in correspondence_chain_strategy(3),
    ) {
        let [first, second, third] = chain.as_slice() else { unreachable!() };
        prop_assert_eq!(
            first.compose(second).unwrap().compose(third),
            first.compose(&second.compose(third).unwrap()),
        );
    }

    #[test]
    fn test_correspondence_compose_identity(
        correspondence in correspondence_strategy::<NodeId>(),
    ) {
        let left_images = (0..correspondence.left_count())
            .map(NodeId::from)
            .collect::<Vec<_>>();
        let right_images = (0..correspondence.right_count())
            .map(NodeId::from)
            .collect::<Vec<_>>();
        let left_identity =
            Correspondence::from_images(&left_images, correspondence.left_count());
        let right_identity =
            Correspondence::from_images(&right_images, correspondence.right_count());

        prop_assert_eq!(left_identity.compose(&correspondence), Ok(correspondence.clone()));
        prop_assert_eq!(correspondence.compose(&right_identity), Ok(correspondence));
    }

    #[test]
    fn test_correspondence_compose_all(
        chain in correspondence_chain_strategy(3),
    ) {
        let expected = chain[0].compose(&chain[1]).unwrap().compose(&chain[2]).unwrap();

        prop_assert_eq!(
            Correspondence::compose_all(chain),
            Ok(Some(expected)),
        );
    }

    #[test]
    fn test_correspondence_compose_all_concatenation(
        correspondences in (0usize..8).prop_flat_map(correspondence_chain_strategy),
        split in any::<usize>(),
    ) {
        let split = split.min(correspondences.len());
        let left = Correspondence::compose_all(correspondences[..split].iter().cloned()).unwrap();
        let right = Correspondence::compose_all(correspondences[split..].iter().cloned()).unwrap();
        let expected = match (left, right) {
            (Some(left), Some(right)) => Some(left.compose(&right).unwrap()),
            (Some(correspondence), None) | (None, Some(correspondence)) => Some(correspondence),
            (None, None) => None,
        };

        prop_assert_eq!(
            Correspondence::compose_all(correspondences),
            Ok(expected),
        );
    }

    #[test]
    fn test_correspondence_reverse_involution(
        correspondence in correspondence_strategy::<NodeId>(),
    ) {
        prop_assert_eq!(correspondence.reverse().reverse(), correspondence);
    }

    #[test]
    fn test_correspondence_edge_matched_pairs(
        (left, right, correspondence) in graph_context_strategy(),
    ) {
        let expected = reference_edge_matched_pairs(&left, &right, &correspondence);

        prop_assert_eq!(
            correspondence.edge_matched_pairs(&left, &right),
            expected.clone(),
        );
        prop_assert_eq!(
            correspondence.shared_edge_count(&left, &right),
            expected.as_ref().map(Vec::len),
        );

        let expected = expected.map(|matched_pairs| {
            GraphCorrespondence::new(
                correspondence.clone(),
                Correspondence::new(matched_pairs, left.edge_count(), right.edge_count())
                    .expect("reference relation is a partial bijection"),
            )
        });
        prop_assert_eq!(
            GraphCorrespondence::induce(&left, &right, correspondence),
            expected,
        );
    }

    #[test]
    fn test_graph_correspondence_is_total(
        correspondence in graph_correspondence_strategy(),
    ) {
        let total_on_left = correspondence.nodes().is_total_on_left()
            && correspondence.edges().is_total_on_left();
        let total_on_right = correspondence.nodes().is_total_on_right()
            && correspondence.edges().is_total_on_right();

        prop_assert_eq!(correspondence.is_total_on_left(), total_on_left);
        prop_assert_eq!(correspondence.is_total_on_right(), total_on_right);
        prop_assert_eq!(correspondence.is_total(), total_on_left && total_on_right);
    }

}
