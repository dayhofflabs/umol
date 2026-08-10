//! Molecule comparison properties.
//!
//! The identity-frame laws for `equiv`, the correspondence laws for
//! `equiv_under`, and agreement with `==` on normalized graph-IR values deliberately use
//! overlapping molecule domains. They establish distinct relations: semantic
//! equivalence in a shared frame, semantic equivalence under an explicit frame
//! mapping, and a normal-form oracle, respectively.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_core::Correspondence;
use umol_graph_ir::ir::MoleculeCorrespondence;

use crate::strategies::*;

fn identity_correspondence(molecule: &Molecule) -> MoleculeCorrespondence {
    fn identity<Id>(count: usize) -> Correspondence<Id>
    where
        Id: Copy + Ord + From<usize>,
    {
        let images: Vec<Id> = (0..count).map(Id::from).collect();
        Correspondence::from_images(&images, count)
    }

    MoleculeCorrespondence::new(
        identity::<AtomId>(molecule.atoms().count()),
        identity::<BondId>(molecule.bonds().count()),
        identity::<DativeBondId>(molecule.dative_bonds().count()),
        identity::<AromaticSystemId>(molecule.aromatic_systems().count()),
        identity::<MulticenterBondId>(molecule.multicenter_bonds().count()),
        identity::<NoncovalentBondId>(molecule.noncovalent_bonds().count()),
        identity::<StereoAtomId>(molecule.stereo_atoms().count()),
        identity::<StereoBondId>(molecule.stereo_bonds().count()),
    )
}

fn atom_only_correspondence(images: &[AtomId], count: usize) -> MoleculeCorrespondence {
    fn empty<Id>() -> Correspondence<Id>
    where
        Id: Copy + Ord + From<usize>,
    {
        Correspondence::from_images(&[], 0)
    }
    MoleculeCorrespondence::new(
        Correspondence::from_images(images, count),
        empty(),
        empty(),
        empty(),
        empty(),
        empty(),
        empty(),
        empty(),
    )
}

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]
    #[test]
    fn test_molecule_equiv_reflexive(molecule in molecule_with_constraints_strategy()) {
        prop_assert!(molecule.equiv(&molecule));
    }

    #[test]
    fn test_molecule_equiv_symmetric(
        left in molecule_with_constraints_strategy(),
        right in molecule_with_constraints_strategy(),
    ) {
        prop_assert_eq!(left.equiv(&right), right.equiv(&left));
    }

    #[test]
    fn test_molecule_equiv_under_transitive(
        atoms in prop::collection::vec(atom_form_strategy(), 0..=5),
    ) {
        let count = atoms.len();
        let first_order = (0..count).collect::<Vec<_>>();
        let second_order = (0..count).rev().collect::<Vec<_>>();
        let mut third_order = first_order.clone();
        if count > 0 {
            third_order.rotate_left(1);
        }
        let molecule = |order: &[usize]| {
            Molecule::from_entries(MoleculeEntries {
                atoms: order.iter().map(|&index| atoms[index].clone()).collect(),
                ..Default::default()
            })
        };
        let correspondence = |left: &[usize], right: &[usize]| {
            let images = left
                .iter()
                .map(|original| {
                    AtomId::from(
                        right
                            .iter()
                            .position(|candidate| candidate == original)
                            .expect("orders contain the same indices"),
                    )
                })
                .collect::<Vec<_>>();
            atom_only_correspondence(&images, count)
        };

        let first = molecule(&first_order);
        let second = molecule(&second_order);
        let third = molecule(&third_order);
        let first_second = correspondence(&first_order, &second_order);
        let second_third = correspondence(&second_order, &third_order);

        prop_assert!(first.equiv_under(&second, &first_second));
        prop_assert!(second.equiv_under(&third, &second_third));
        prop_assert!(first.equiv_under(&third, &first_second.compose(&second_third)));
    }

    #[test]
    fn test_molecule_equiv_agrees_with_equality_for_normalized_molecules(
        left in molecule_strategy(),
        right in molecule_strategy(),
    ) {
        prop_assert_eq!(left.equiv(&right), left == right);
    }

    #[test]
    fn test_molecule_equiv_under_identity_reduces_to_equiv(
        molecule in molecule_with_constraints_strategy(),
    ) {
        let correspondence = identity_correspondence(&molecule);
        let mut other = molecule.clone();
        if other.atoms().count() > 0 {
            other.atom_mut(AtomId(0)).attributes.charge = NumForm::Lit(99);
        }

        prop_assert_eq!(
            molecule.equiv_under(&other, &correspondence),
            molecule.equiv(&other),
        );
    }

    #[test]
    fn test_molecule_equiv_under_symmetric_under_reverse(
        atoms in prop::collection::vec(atom_form_strategy(), 0..=5),
        change_mapped_atom in any::<bool>(),
    ) {
        let count = atoms.len();
        let left = Molecule::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
            ..Default::default()
        });
        let mut right = Molecule::from_entries(MoleculeEntries {
            atoms: atoms.into_iter().rev().collect(),
            ..Default::default()
        });
        if change_mapped_atom && count > 0 {
            right.atom_mut(AtomId((count - 1) as u32)).attributes.charge = NumForm::Lit(99);
        }
        let images: Vec<AtomId> = (0..count).rev().map(AtomId::from).collect();
        let correspondence = MoleculeCorrespondence::new(
            Correspondence::from_images(&images, count),
            Correspondence::from_images(&[], 0),
            Correspondence::from_images(&[], 0),
            Correspondence::from_images(&[], 0),
            Correspondence::from_images(&[], 0),
            Correspondence::from_images(&[], 0),
            Correspondence::from_images(&[], 0),
            Correspondence::from_images(&[], 0),
        );

        let forward = left.equiv_under(&right, &correspondence);
        let reverse = right.equiv_under(&left, &correspondence.reverse());
        prop_assert_eq!(forward, reverse);
        prop_assert_eq!(forward, !change_mapped_atom || count == 0);
    }
}
