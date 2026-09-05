//! Molecule-correspondence algebra properties.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};

use crate::strategies::*;

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_molecule_correspondence_reverse(
        (molecule, atoms) in molecule_with_atom_subset_strategy(),
    ) {
        let correspondence = molecule.induced_subgraph(&atoms);
        prop_assert_eq!(correspondence.reverse().reverse(), correspondence);
    }

    #[test]
    fn test_molecule_correspondence_induce(
        (molecule, atoms) in molecule_with_atom_subset_strategy(),
    ) {
        let correspondence = molecule.induced_subgraph(&atoms);
        let extracted = molecule.extract(&correspondence);
        let induced = MoleculeCorrespondence::induce(
            &extracted,
            &molecule,
            correspondence.atoms().clone(),
        ).expect("extraction preserves unique entity incidence");
        prop_assert_eq!(induced, correspondence);
    }

    #[test]
    fn test_molecule_correspondence_is_total(
        (molecule, atoms) in molecule_with_atom_subset_strategy(),
    ) {
        let correspondence = molecule.induced_subgraph(&atoms);
        let reverse = correspondence.reverse();

        prop_assert_eq!(
            correspondence.is_total(),
            correspondence.is_total_on_left() && correspondence.is_total_on_right(),
        );
        prop_assert_eq!(correspondence.is_total_on_left(), reverse.is_total_on_right());
        prop_assert_eq!(correspondence.is_total_on_right(), reverse.is_total_on_left());
    }

    #[test]
    fn test_molecule_correspondence_compose_identity(
        (molecule, atoms) in molecule_with_atom_subset_strategy(),
    ) {
        let correspondence = molecule.induced_subgraph(&atoms);
        let identity = correspondence.compose(&correspondence.reverse()).unwrap();

        prop_assert!(identity.is_total());
        for index in 0..identity.atoms().left_count() {
            let id = AtomId(index as u32);
            prop_assert_eq!(identity.atoms().right_of(id), Some(id));
        }
        for index in 0..identity.bonds().left_count() {
            let id = BondId(index as u32);
            prop_assert_eq!(identity.bonds().right_of(id), Some(id));
        }
        for index in 0..identity.dative_bonds().left_count() {
            let id = DativeBondId(index as u32);
            prop_assert_eq!(identity.dative_bonds().right_of(id), Some(id));
        }
        for index in 0..identity.aromatic_systems().left_count() {
            let id = AromaticSystemId(index as u32);
            prop_assert_eq!(identity.aromatic_systems().right_of(id), Some(id));
        }
        for index in 0..identity.multicenter_bonds().left_count() {
            let id = MulticenterBondId(index as u32);
            prop_assert_eq!(identity.multicenter_bonds().right_of(id), Some(id));
        }
        for index in 0..identity.noncovalent_bonds().left_count() {
            let id = NoncovalentBondId(index as u32);
            prop_assert_eq!(identity.noncovalent_bonds().right_of(id), Some(id));
        }
        for index in 0..identity.stereo_atoms().left_count() {
            let id = StereoAtomId(index as u32);
            prop_assert_eq!(identity.stereo_atoms().right_of(id), Some(id));
        }
        for index in 0..identity.stereo_bonds().left_count() {
            let id = StereoBondId(index as u32);
            prop_assert_eq!(identity.stereo_bonds().right_of(id), Some(id));
        }
    }

    #[test]
    fn test_molecule_correspondence_compose_associativity(
        (molecule, atoms) in molecule_with_atom_subset_strategy(),
    ) {
        let correspondence = molecule.induced_subgraph(&atoms);
        let reverse = correspondence.reverse();

        prop_assert_eq!(
            correspondence.compose(&reverse).unwrap().compose(&correspondence),
            correspondence.compose(&reverse.compose(&correspondence).unwrap()),
        );
    }

    #[test]
    fn test_molecule_correspondence_compose_all(
        (molecule, atoms) in molecule_with_atom_subset_strategy(),
    ) {
        let correspondence = molecule.induced_subgraph(&atoms);
        let reverse = correspondence.reverse();
        let expected = correspondence.compose(&reverse).unwrap().compose(&correspondence).unwrap();

        prop_assert_eq!(
            MoleculeCorrespondence::compose_all([
                correspondence.clone(),
                reverse,
                correspondence,
            ]),
            Ok(Some(expected)),
        );
    }
}
