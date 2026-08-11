//! Property tests for the permutation / coset algebra: structural round-trips
//! and group / coset-space laws.

use std::ptr;

use proptest::prelude::*;
use umol_perm::{
    ClassKey, Coset, Orientation, OrientedPermutation, OrientedPermutationGroup, Permutation,
    PermutationGroup,
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

fn permutation_image() -> impl Strategy<Value = Vec<usize>> {
    (0usize..=6).prop_flat_map(|degree| Just((0..degree).collect::<Vec<_>>()).prop_shuffle())
}

fn repeated_orderings() -> impl Strategy<Value = (Vec<u8>, Vec<u8>)> {
    prop::collection::vec(0u8..3, 0..=6).prop_flat_map(|source| {
        Just(source.clone())
            .prop_shuffle()
            .prop_map(move |target| (source.clone(), target))
    })
}

fn perm_pair() -> impl Strategy<Value = (Permutation, Permutation)> {
    (2usize..=6).prop_flat_map(|d| (perm_of(d), perm_of(d)))
}

fn perm_triple() -> impl Strategy<Value = (Permutation, Permutation, Permutation)> {
    (2usize..=6).prop_flat_map(|d| (perm_of(d), perm_of(d), perm_of(d)))
}

fn perm_pair_with_items() -> impl Strategy<Value = (Permutation, Permutation, Vec<u16>)> {
    (2usize..=6).prop_flat_map(|degree| {
        (
            perm_of(degree),
            perm_of(degree),
            prop::collection::vec(any::<u16>(), degree),
        )
    })
}

fn permutation_group() -> impl Strategy<Value = (usize, Vec<Permutation>, PermutationGroup)> {
    (2usize..=4).prop_flat_map(|degree| {
        prop::collection::vec(perm_of(degree), 0..=3).prop_map(move |generators| {
            let group = PermutationGroup::generate(degree, &generators);
            (degree, generators, group)
        })
    })
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
fn oriented_group(
) -> impl Strategy<Value = (usize, Vec<OrientedPermutation>, OrientedPermutationGroup)> {
    (2usize..=4).prop_flat_map(|degree| {
        prop::collection::vec(oriented_of(degree), 0..=3).prop_map(move |generators| {
            let group = OrientedPermutationGroup::generate(degree, &generators);
            (degree, generators, group)
        })
    })
}

fn class_key_text() -> impl Strategy<Value = ClassKey> {
    prop_oneof![
        (0u8..=6).prop_map(ClassKey::Symmetric),
        (0u8..=6).prop_map(ClassKey::Alternating),
        (0u8..=6).prop_map(ClassKey::Cyclic),
        (0u8..=6).prop_map(ClassKey::Dihedral),
        Just(ClassKey::Tetrahedral),
        Just(ClassKey::CisTrans),
        Just(ClassKey::Axial),
        Just(ClassKey::SquarePlanar),
        Just(ClassKey::TrigonalBipyramidal),
        Just(ClassKey::Octahedral),
    ]
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
    class_key().prop_flat_map(|key| (Just(key), 0..key.space().count() as u32))
}

fn coset_indices() -> impl Strategy<Value = (ClassKey, u32, u32)> {
    class_key().prop_flat_map(|key| {
        let count = key.space().count() as u32;
        (Just(key), 0..count, 0..count)
    })
}

fn coset_generators() -> impl Strategy<Value = (ClassKey, Vec<Permutation>)> {
    class_key().prop_flat_map(|key| {
        let space = key.space();
        let group_order = space.group().order();
        let parent_order = space.count() * group_order;
        prop::collection::vec(0..parent_order, 0..=4).prop_map(move |indices| {
            let generators = indices
                .into_iter()
                .map(|index| {
                    let group_element = space.group().elements()[index % group_order];
                    let coset = space
                        .unindex((index / group_order) as u32)
                        .expect("generated coset index is in range");
                    group_element.compose(coset)
                })
                .collect();
            (key, generators)
        })
    })
}

proptest! {
    #[test]
    fn test_permutation_cycle_round_trip(p in permutation()) {
        prop_assert_eq!(Permutation::from_cycles(p.degree(), &p.cycles()), Ok(p));
    }

    #[test]
    fn test_permutation_rank_round_trip(p in permutation()) {
        prop_assert_eq!(Permutation::unrank(p.degree(), p.rank()), p);
    }

    #[test]
    fn test_permutation_image_round_trip(image in permutation_image()) {
        let permutation = Permutation::try_from(image.as_slice()).unwrap();
        let recovered = (0..permutation.degree())
            .map(|point| permutation.apply(point))
            .collect::<Vec<_>>();
        prop_assert_eq!(permutation.degree(), image.len());
        prop_assert_eq!(recovered, image);
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
    fn test_permutation_between(p in permutation()) {
        let source = (0..p.degree()).collect::<Vec<_>>();
        let target = p.act(&source);
        prop_assert_eq!(Permutation::between(&source, &target), Some(p));
        prop_assert_eq!(Permutation::between(&target, &source), Some(p.inverse()));
    }

    #[test]
    fn test_permutation_between_all((source, target) in repeated_orderings()) {
        let mut actual = Permutation::between_all(&source, &target);
        let mut expected = (0..factorial(source.len()))
            .map(|rank| Permutation::unrank(source.len(), rank))
            .filter(|permutation| permutation.act(&source) == target)
            .collect::<Vec<_>>();
        actual.sort_unstable();
        expected.sort_unstable();

        prop_assert_eq!(&actual, &expected);
        prop_assert_eq!(
            Permutation::between(&source, &target),
            (actual.len() == 1).then(|| actual[0])
        );
    }

    #[test]
    fn test_permutation_act_composition((a, b, items) in perm_pair_with_items()) {
        prop_assert_eq!(a.compose(b).act(&items), b.act(&a.act(&items)));
    }

    #[test]
    fn test_permutation_group_generate((degree, generators, group) in permutation_group()) {
        prop_assert!(group.contains(&Permutation::identity(degree)));
        for generator in generators {
            prop_assert!(group.contains(&generator));
        }
        for &a in group.elements() {
            prop_assert!(group.contains(&a.inverse()));
            for &b in group.elements() {
                prop_assert!(group.contains(&a.compose(b)));
            }
        }
    }

    #[test]
    fn test_oriented_permutation_group_generate((degree, generators, group) in oriented_group()) {
        prop_assert!(group.contains(OrientedPermutation::identity(degree)));
        for generator in generators {
            prop_assert!(group.contains(generator));
        }
    }

    #[test]
    fn test_oriented_permutation_group_closed((_degree, _generators, group) in oriented_group()) {
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
        let s = key.space();
        prop_assert_eq!(s.index(s.unindex(i).unwrap()), Some(i));
    }

    #[test]
    fn test_coset_space_enantiomer_involution((key, i) in coset_index()) {
        let s = key.space();
        prop_assert_eq!(s.enantiomer(s.enantiomer(i).unwrap()), Some(i));
    }

    #[test]
    fn test_coset_space_reindex_identity((key, i) in coset_index()) {
        let s = key.space();
        prop_assert_eq!(s.reindex(i, Permutation::identity(s.degree())), Some(i));
    }

    #[test]
    fn test_coset_space_observable_coset((key, index, generator_index) in coset_indices()) {
        let s = key.space();
        let generator = s.unindex(generator_index).unwrap();
        let expected = s.orbit_reps(&[generator]).unwrap()[index as usize];
        prop_assert_eq!(s.observable_coset(index, &[generator]), Some(expected));
    }

    #[test]
    fn test_coset_space_observable_coset_index(key in class_key()) {
        let s = key.space();
        prop_assert_eq!(s.observable_coset(s.count() as u32, &[]), None);
    }

    #[test]
    fn test_coset_space_observable_coset_generator((key, index) in coset_index()) {
        let s = key.space();
        let generator = Permutation::identity(s.degree() - 1);
        prop_assert_eq!(s.observable_coset(index, &[generator]), None);
    }

    #[test]
    fn test_coset_space_orbit_reps_identity(key in class_key()) {
        let s = key.space();
        let expected = (0..s.count() as u32).collect::<Vec<u32>>();
        prop_assert_eq!(s.orbit_reps(&[]), Some(expected.clone()));
        prop_assert_eq!(s.orbit_reps(&[Permutation::identity(s.degree())]), Some(expected));
    }

    #[test]
    fn test_coset_space_orbit_reps((key, generators) in coset_generators()) {
        let space = key.space();
        let mut expected = Vec::with_capacity(space.count());
        for start in 0..space.count() as u32 {
            let mut reached = vec![start];
            let mut frontier = vec![start];
            while let Some(index) = frontier.pop() {
                for &generator in &generators {
                    let moved = space.reindex(index, generator).unwrap();
                    if !reached.contains(&moved) {
                        reached.push(moved);
                        frontier.push(moved);
                    }
                }
            }
            expected.push(*reached.iter().min().unwrap());
        }
        prop_assert_eq!(space.orbit_reps(&generators), Some(expected));
    }

    #[test]
    fn test_class_key_display_from_str_roundtrip(key in class_key_text()) {
        prop_assert_eq!(key.to_string().parse::<ClassKey>(), Ok(key));
    }

    #[test]
    fn test_class_key_space_interning(key in class_key()) {
        prop_assert!(ptr::eq(key.space(), key.space()));
    }

    #[test]
    fn test_coset_space((key, index) in coset_index()) {
        let coset = Coset::new(key, index);
        prop_assert!(ptr::eq(coset.space(), key.space()));
    }
}
