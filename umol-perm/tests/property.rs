//! Property tests for the permutation / coset algebra: structural round-trips
//! and group / coset-space laws.

use proptest::prelude::*;
use umol_perm::{
    space, ClassKey, Orientation, OrientedPermutation, OrientedPermutationGroup, Permutation,
};

fn factorial(n: usize) -> usize {
    (1..=n).product::<usize>().max(1)
}

/// A uniform random permutation of `degree` positions, via its Lehmer rank.
fn perm_of(degree: usize) -> impl Strategy<Value = Permutation> {
    (0..factorial(degree)).prop_map(move |rank| Permutation::unrank(degree, rank))
}

fn permutation() -> impl Strategy<Value = Permutation> {
    (2usize..=6).prop_flat_map(perm_of)
}

fn perm_pair() -> impl Strategy<Value = (Permutation, Permutation)> {
    (2usize..=6).prop_flat_map(|d| (perm_of(d), perm_of(d)))
}

fn perm_triple() -> impl Strategy<Value = (Permutation, Permutation, Permutation)> {
    (2usize..=6).prop_flat_map(|d| (perm_of(d), perm_of(d), perm_of(d)))
}

fn oriented_of(degree: usize) -> impl Strategy<Value = OrientedPermutation> {
    (
        perm_of(degree),
        prop_oneof![Just(Orientation::Proper), Just(Orientation::Improper)],
    )
        .prop_map(|(permutation, orientation)| OrientedPermutation::new(permutation, orientation))
}

/// A generated oriented group, paired with its degree. Degree is held low so the
/// closure check (over all element pairs) stays cheap.
fn oriented_group() -> impl Strategy<Value = (usize, OrientedPermutationGroup)> {
    (2usize..=4).prop_flat_map(|degree| {
        prop::collection::vec(oriented_of(degree), 0..=3)
            .prop_map(move |gens| (degree, OrientedPermutationGroup::generate(degree, &gens)))
    })
}

fn class_key() -> impl Strategy<Value = ClassKey> {
    prop_oneof![
        Just(ClassKey::Tetrahedral),
        Just(ClassKey::CisTrans),
        Just(ClassKey::Axial),
        Just(ClassKey::SquarePlanar),
        Just(ClassKey::TrigonalBipyramidal),
        Just(ClassKey::Octahedral),
    ]
}

/// A class key paired with a valid coset index in its space.
fn coset_index() -> impl Strategy<Value = (ClassKey, u32)> {
    class_key().prop_flat_map(|key| (Just(key), 0..space(key).count() as u32))
}

proptest! {
    #[test]
    fn test_permutation_cycle_round_trip(p in permutation()) {
        prop_assert_eq!(Permutation::from_cycles(p.degree(), &p.cycles()), p);
    }

    #[test]
    fn test_permutation_rank_round_trip(p in permutation()) {
        prop_assert_eq!(Permutation::unrank(p.degree(), p.rank()), p);
    }

    #[test]
    fn test_permutation_inverse_involution(p in permutation()) {
        prop_assert_eq!(p.inverse().inverse(), p);
    }

    #[test]
    fn test_permutation_compose_inverse_identity(p in permutation()) {
        let identity = Permutation::identity(p.degree());
        prop_assert_eq!(p.compose(p.inverse()), identity);
        prop_assert_eq!(p.inverse().compose(p), identity);
    }

    #[test]
    fn test_permutation_compose_associative((a, b, c) in perm_triple()) {
        prop_assert_eq!(a.compose(b).compose(c), a.compose(b.compose(c)));
    }

    #[test]
    fn test_permutation_sign_homomorphism((a, b) in perm_pair()) {
        prop_assert_eq!(a.compose(b).sign(), a.sign() * b.sign());
    }

    #[test]
    fn test_permutation_sign_inverse(p in permutation()) {
        prop_assert_eq!(p.inverse().sign(), p.sign());
    }

    #[test]
    fn test_permutation_identity_sign(degree in 2usize..=6) {
        prop_assert_eq!(Permutation::identity(degree).sign(), 1);
    }

    #[test]
    fn test_oriented_permutation_group_contains_identity((degree, group) in oriented_group()) {
        prop_assert!(group.contains(OrientedPermutation::identity(degree)));
    }

    #[test]
    fn test_oriented_permutation_group_closed((_degree, group) in oriented_group()) {
        let elements = group.elements();
        prop_assert_eq!(group.order(), elements.len());
        for &a in &elements {
            prop_assert!(group.contains(a.inverse()));
            for &b in &elements {
                prop_assert!(group.contains(a.compose(b)));
            }
        }
    }

    #[test]
    fn test_coset_space_index_round_trip((key, i) in coset_index()) {
        let s = space(key);
        prop_assert_eq!(s.index(s.unindex(i)), i);
    }

    #[test]
    fn test_coset_space_enantiomer_involution((key, i) in coset_index()) {
        let s = space(key);
        prop_assert_eq!(s.enantiomer(s.enantiomer(i)), i);
    }

    #[test]
    fn test_coset_space_reindex_identity((key, i) in coset_index()) {
        let s = space(key);
        prop_assert_eq!(s.reindex(i, Permutation::identity(s.degree())), i);
    }

    #[test]
    fn test_coset_space_observable_coset_no_fluxional((key, i) in coset_index()) {
        let s = space(key);
        prop_assert_eq!(s.observable_coset(i, &[]), i);
    }

    #[test]
    fn test_coset_space_merge_under_empty(key in class_key()) {
        let s = space(key);
        prop_assert_eq!(s.merge_under(&[]), (0..s.count() as u32).collect::<Vec<u32>>());
    }
}
