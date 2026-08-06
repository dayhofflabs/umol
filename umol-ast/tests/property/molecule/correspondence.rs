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
        (ast, atoms) in molecule_ast_with_atom_subset_strategy(),
    ) {
        let correspondence = ast.induced_subgraph(&atoms);
        prop_assert_eq!(correspondence.reverse().reverse(), correspondence);
    }

    #[test]
    fn test_molecule_correspondence_induce(
        (ast, atoms) in molecule_ast_with_atom_subset_strategy(),
    ) {
        let correspondence = ast.induced_subgraph(&atoms);
        let extracted = ast.extract(&correspondence);
        let induced = MoleculeCorrespondence::induce(
            &extracted,
            &ast,
            correspondence.atoms().clone(),
        ).expect("extraction preserves unique entity incidence");
        prop_assert_eq!(induced, correspondence);
    }

    #[test]
    fn test_molecule_correspondence_is_total(
        (ast, atoms) in molecule_ast_with_atom_subset_strategy(),
    ) {
        let correspondence = ast.induced_subgraph(&atoms);
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
        (ast, atoms) in molecule_ast_with_atom_subset_strategy(),
    ) {
        let correspondence = ast.induced_subgraph(&atoms);
        let identity = correspondence.compose(&correspondence.reverse());

        prop_assert!(identity.is_total());
        let remapping = identity
            .to_remapping()
            .expect("the composed identity is total on the left");
        for index in 0..identity.atoms().left_count() {
            let id = AtomId(index as u32);
            prop_assert_eq!(remapping.map_atom(id), id);
        }
        for index in 0..identity.bonds().left_count() {
            let id = BondId(index as u32);
            prop_assert_eq!(remapping.map_bond(id), id);
        }
        for index in 0..identity.dative_bonds().left_count() {
            let id = DativeBondId(index as u32);
            prop_assert_eq!(remapping.map_dative(id), id);
        }
        for index in 0..identity.aromatic_systems().left_count() {
            let id = AromaticSystemId(index as u32);
            prop_assert_eq!(remapping.map_aromatic(id), id);
        }
        for index in 0..identity.multicenter_bonds().left_count() {
            let id = MulticenterBondId(index as u32);
            prop_assert_eq!(remapping.map_multicenter(id), id);
        }
        for index in 0..identity.noncovalent_bonds().left_count() {
            let id = NoncovalentBondId(index as u32);
            prop_assert_eq!(remapping.map_noncovalent(id), id);
        }
        for index in 0..identity.stereo_atoms().left_count() {
            let id = StereoAtomId(index as u32);
            prop_assert_eq!(remapping.map_stereo_atom(id), id);
        }
        for index in 0..identity.stereo_bonds().left_count() {
            let id = StereoBondId(index as u32);
            prop_assert_eq!(remapping.map_stereo_bond(id), id);
        }
    }

    #[test]
    fn test_molecule_correspondence_compose_associativity(
        (ast, atoms) in molecule_ast_with_atom_subset_strategy(),
    ) {
        let correspondence = ast.induced_subgraph(&atoms);
        let reverse = correspondence.reverse();

        prop_assert_eq!(
            correspondence.compose(&reverse).compose(&correspondence),
            correspondence.compose(&reverse.compose(&correspondence)),
        );
    }

    #[test]
    fn test_molecule_correspondence_compose_all(
        (ast, atoms) in molecule_ast_with_atom_subset_strategy(),
    ) {
        let correspondence = ast.induced_subgraph(&atoms);
        let reverse = correspondence.reverse();
        let expected = correspondence.compose(&reverse).compose(&correspondence);

        prop_assert_eq!(
            MoleculeCorrespondence::compose_all([
                correspondence.clone(),
                reverse,
                correspondence,
            ]),
            Some(expected),
        );
    }
}
