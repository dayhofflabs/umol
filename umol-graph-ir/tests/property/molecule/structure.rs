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
    fn test_molecule_ast_try_from_entries(
        entries in molecule_entries_with_constraints_strategy(),
    ) {
        let expected = MoleculeAst::from_entries(entries.clone());

        prop_assert_eq!(MoleculeAst::try_from_entries(entries), Ok(expected));
    }

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
        ).expect("extraction preserves unique entity incidence");

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
                    .map(|&(_, host)| host),
            );
        }

        covered_atoms.sort_unstable();
        prop_assert_eq!(
            covered_atoms,
            (0..ast.atoms().count()).map(AtomId::from).collect::<Vec<_>>(),
        );
    }

    #[test]
    fn test_molecule_ast_constraint_atoms_unpaired_electron_coupling(
        (atom_count, subset_mask) in (1usize..=8).prop_flat_map(|atom_count| (
            Just(atom_count),
            prop::collection::vec(any::<bool>(), atom_count),
        )),
    ) {
        let ast = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomAst::from_element(Element::C); atom_count],
            ..Default::default()
        });
        let all_atoms = (0..atom_count).map(AtomId::from).collect::<Vec<_>>();
        let subset = subset_mask
            .into_iter()
            .enumerate()
            .filter_map(|(index, include)| include.then_some(AtomId::from(index)))
            .collect::<Vec<_>>();
        let unpaired_electrons = UnpairedElectronsForm::from((0_u8, 1_u8));

        prop_assert_eq!(
            ast.constraint_atoms(&Constraint::Molecule(
                MoleculeConstraint::UnpairedElectronCoupling {
                    atoms: None,
                    unpaired_electrons: unpaired_electrons.clone(),
                },
            )),
            all_atoms,
        );
        prop_assert_eq!(
            ast.constraint_atoms(&Constraint::Molecule(
                MoleculeConstraint::UnpairedElectronCoupling {
                    atoms: Some(subset.clone()),
                    unpaired_electrons,
                },
            )),
            subset,
        );
    }
}
