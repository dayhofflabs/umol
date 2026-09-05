//! Value transport under dense permutations.
//!
//! The vector-remapping law is checked for shuffled permutations of zero through 64 entries.
//! Sorting image/value pairs supplies a reference independent of indexed placement.

use proptest::prelude::*;
use umol_graph_core::{NodeId, Remapping};

proptest! {
    #[test]
    fn test_remapping_remap_vec(
        (values, images) in prop::collection::vec(any::<i64>(), 0..65)
            .prop_flat_map(|values| {
                let images = (0..values.len()).map(NodeId::from).collect::<Vec<_>>();
                (Just(values), Just(images).prop_shuffle())
            }),
    ) {
        let values = values.into_iter().enumerate().collect::<Vec<_>>();
        let mut paired = images.iter().copied().zip(values.iter().copied()).collect::<Vec<_>>();
        paired.sort_by_key(|&(image, _)| image);
        let expected = paired.into_iter().map(|(_, value)| value).collect::<Vec<_>>();
        let remapping = Remapping::new(images).unwrap();

        prop_assert_eq!(remapping.remap_vec(values.clone()), expected.clone());
        prop_assert_eq!(remapping.try_remap_vec(values), Some(expected));
    }
}
