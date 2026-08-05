//! Algebraic properties of correspondence composition.

use proptest::prelude::*;
use umol_graph_core::{Correspondence, NodeId};

fn correspondence_strategy() -> impl Strategy<Value = Correspondence<NodeId>> {
    (0usize..8, 0usize..8).prop_flat_map(|(left_count, right_count)| {
        let pair_count = left_count.min(right_count);
        (
            Just((0..left_count).map(NodeId::from).collect::<Vec<_>>()).prop_shuffle(),
            Just((0..right_count).map(NodeId::from).collect::<Vec<_>>()).prop_shuffle(),
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
    })
}

proptest! {
    #[test]
    fn test_correspondence_compose_associativity(
        first in correspondence_strategy(),
        second in correspondence_strategy(),
        third in correspondence_strategy(),
    ) {
        prop_assert_eq!(
            first.compose(&second).compose(&third),
            first.compose(&second.compose(&third)),
        );
    }

    #[test]
    fn test_correspondence_compose_identity(
        correspondence in correspondence_strategy(),
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

        prop_assert_eq!(left_identity.compose(&correspondence), correspondence.clone());
        prop_assert_eq!(correspondence.compose(&right_identity), correspondence);
    }

    #[test]
    fn test_correspondence_compose_all(
        first in correspondence_strategy(),
        second in correspondence_strategy(),
        third in correspondence_strategy(),
    ) {
        let expected = first.compose(&second).compose(&third);

        prop_assert_eq!(
            Correspondence::compose_all([first, second, third]),
            Some(expected),
        );
    }

    #[test]
    fn test_correspondence_compose_all_concatenation(
        correspondences in prop::collection::vec(correspondence_strategy(), 0..8),
        split in any::<usize>(),
    ) {
        let split = split.min(correspondences.len());
        let left = Correspondence::compose_all(correspondences[..split].iter().cloned());
        let right = Correspondence::compose_all(correspondences[split..].iter().cloned());
        let expected = match (left, right) {
            (Some(left), Some(right)) => Some(left.compose(&right)),
            (Some(correspondence), None) | (None, Some(correspondence)) => Some(correspondence),
            (None, None) => None,
        };

        prop_assert_eq!(
            Correspondence::compose_all(correspondences),
            expected,
        );
    }
}
