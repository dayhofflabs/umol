//! Molecule extraction, combination, and splitting properties.

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
    fn test_molecule_ast_extract(
        (ast, atoms) in molecule_ast_with_atom_subset_strategy(),
    ) {
        let correspondence = ast.induced_subgraph(&atoms);
        let extracted = ast.extract(&correspondence);
        let reinduced = MoleculeCorrespondence::induce(
            &extracted,
            &ast,
            correspondence.atoms().clone(),
        );

        prop_assert_eq!(&reinduced, &correspondence);
        prop_assert_eq!(ast.extract(&reinduced), extracted);
    }

    #[test]
    fn test_molecule_ast_combine_all(
        molecules in prop::collection::vec(
            molecule_ast_structurally_unambiguous_strategy(),
            0..5,
        ),
    ) {
        let (combined, correspondences) = MoleculeAst::combine_all(&molecules);
        prop_assert_eq!(correspondences.len(), molecules.len());
        for (molecule, correspondence) in molecules.iter().zip(&correspondences) {
            prop_assert_eq!(&combined.extract(correspondence), molecule);
        }
    }

    #[test]
    fn test_molecule_ast_combine(
        left in molecule_ast_structurally_unambiguous_strategy(),
        right in molecule_ast_structurally_unambiguous_strategy(),
    ) {
        let (combined, correspondence) = left.combine(&right);
        prop_assert_eq!(combined.extract(&correspondence), right);
    }

    #[test]
    fn test_molecule_ast_combine_from(
        left in molecule_ast_structurally_unambiguous_strategy(),
        right in molecule_ast_structurally_unambiguous_strategy(),
    ) {
        let (expected, expected_correspondence) = left.combine(&right);
        let mut combined = left;
        let correspondence = combined.combine_from(&right);

        prop_assert_eq!(combined, expected);
        prop_assert_eq!(correspondence, expected_correspondence);
    }

    #[test]
    fn test_molecule_ast_split(ast in molecule_ast_structurally_unambiguous_strategy()) {
        let components = ast.split();
        let mut covered_atoms = Vec::new();

        for (component, correspondence) in &components {
            prop_assert_eq!(&ast.extract(correspondence), component);
            covered_atoms.extend(
                correspondence
                    .atoms()
                    .matched_pairs()
                    .iter()
                    .map(|&(_, host)| AtomId::from(host)),
            );
        }

        covered_atoms.sort_unstable();
        prop_assert_eq!(
            covered_atoms,
            (0..ast.atoms().count()).map(AtomId::from).collect::<Vec<_>>(),
        );
    }
}
