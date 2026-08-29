//! Molecule extraction, combination, and splitting properties.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_perm::MAX_DEGREE;

use crate::strategies::*;

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_molecule_try_from_entries(
        entries in molecule_entries_with_constraints_strategy(),
    ) {
        let expected = Molecule::from_entries(entries.clone());

        prop_assert_eq!(Molecule::try_from_entries(entries), Ok(expected));
    }

    #[test]
    fn test_molecule_stereo_frame_integrity(molecule in molecule_with_constraints_strategy()) {
        for stereo_atom in molecule.stereo_atoms().iter() {
            let frame = stereo_atom.ligand_frame();
            prop_assert!(frame.len() <= MAX_DEGREE);
            for (position, ligand) in frame.iter().enumerate() {
                prop_assert!(!frame[..position].contains(ligand));
            }
        }
        for stereo_bond in molecule.stereo_bonds().iter() {
            let frame = stereo_bond.ligand_frame();
            prop_assert!(frame.len() <= MAX_DEGREE);
            for (position, ligand) in frame.iter().enumerate() {
                prop_assert!(!frame[..position].contains(ligand));
            }
        }
    }

    #[test]
    fn test_molecule_extract(
        (molecule, atoms) in molecule_with_atom_subset_strategy(),
    ) {
        let correspondence = molecule.induced_subgraph(&atoms);
        let extracted = molecule.extract(&correspondence);
        let reinduced = MoleculeCorrespondence::induce(
            &extracted,
            &molecule,
            correspondence.atoms().clone(),
        ).expect("extraction preserves unique entity incidence");

        prop_assert_eq!(&reinduced, &correspondence);
        prop_assert_eq!(molecule.extract(&reinduced), extracted);
    }

    #[test]
    fn test_molecule_combine_all(
        molecules in prop::collection::vec(
            molecule_structurally_unambiguous_strategy(),
            0..5,
        ),
    ) {
        let (combined, correspondences) = Molecule::combine_all(&molecules);
        prop_assert_eq!(correspondences.len(), molecules.len());
        for (molecule, correspondence) in molecules.iter().zip(&correspondences) {
            prop_assert_eq!(&combined.extract(correspondence), molecule);
        }
    }

    #[test]
    fn test_molecule_combine(
        left in molecule_structurally_unambiguous_strategy(),
        right in molecule_structurally_unambiguous_strategy(),
    ) {
        let (combined, correspondence) = left.combine(&right);
        prop_assert_eq!(combined.extract(&correspondence), right);
    }

    #[test]
    fn test_molecule_combine_from(
        left in molecule_structurally_unambiguous_strategy(),
        right in molecule_structurally_unambiguous_strategy(),
    ) {
        let (expected, expected_correspondence) = left.combine(&right);
        let mut combined = left;
        let correspondence = combined.combine_from(&right);

        prop_assert_eq!(combined, expected);
        prop_assert_eq!(correspondence, expected_correspondence);
    }

    #[test]
    fn test_molecule_split(molecule in molecule_structurally_unambiguous_strategy()) {
        let components = molecule.split();
        let mut covered_atoms = Vec::new();

        for (component, correspondence) in &components {
            prop_assert_eq!(&molecule.extract(correspondence), component);
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
            (0..molecule.atoms().count()).map(AtomId::from).collect::<Vec<_>>(),
        );
    }

    #[test]
    fn test_molecule_constraint_atoms_unpaired_electron_coupling(
        (atom_count, subset_mask) in (1usize..=8).prop_flat_map(|atom_count| (
            Just(atom_count),
            prop::collection::vec(any::<bool>(), atom_count),
        )),
    ) {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); atom_count],
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
            molecule.constraint_atoms(&Constraint::Molecule(
                MoleculeConstraint::UnpairedElectronCoupling {
                    atoms: None,
                    unpaired_electrons: unpaired_electrons.clone(),
                },
            )),
            all_atoms,
        );
        prop_assert_eq!(
            molecule.constraint_atoms(&Constraint::Molecule(
                MoleculeConstraint::UnpairedElectronCoupling {
                    atoms: Some(subset.clone()),
                    unpaired_electrons,
                },
            )),
            subset,
        );
    }
}
