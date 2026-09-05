//! Molecule extraction, combination, and splitting properties.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_core::Correspondence;
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
        let (tracked, compaction) = molecule.tracked_extract(&correspondence);
        prop_assert_eq!(&tracked, &extracted);
        prop_assert_eq!(MoleculeCorrespondence::from(&compaction), correspondence.reverse());
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
        let combined = Molecule::combine_all(&molecules);
        prop_assert_eq!(
            combined.atoms().count(),
            molecules.iter().map(|input| input.atoms().count()).sum::<usize>(),
        );
        prop_assert_eq!(
            combined.bonds().count(),
            molecules.iter().map(|input| input.bonds().count()).sum::<usize>(),
        );
        prop_assert_eq!(
            combined.dative_bonds().count(),
            molecules.iter().map(|input| input.dative_bonds().count()).sum::<usize>(),
        );
        prop_assert_eq!(
            combined.aromatic_systems().count(),
            molecules.iter().map(|input| input.aromatic_systems().count()).sum::<usize>(),
        );
        prop_assert_eq!(
            combined.multicenter_bonds().count(),
            molecules.iter().map(|input| input.multicenter_bonds().count()).sum::<usize>(),
        );
        prop_assert_eq!(
            combined.noncovalent_bonds().count(),
            molecules.iter().map(|input| input.noncovalent_bonds().count()).sum::<usize>(),
        );
        prop_assert_eq!(
            combined.stereo_atoms().count(),
            molecules.iter().map(|input| input.stereo_atoms().count()).sum::<usize>(),
        );
        prop_assert_eq!(
            combined.stereo_bonds().count(),
            molecules.iter().map(|input| input.stereo_bonds().count()).sum::<usize>(),
        );
        for (source_idx, molecule) in molecules.iter().enumerate() {
            let correspondence = MoleculeCorrespondence::new(
                Correspondence::from_images(
                    &(0..molecule.atoms().count())
                        .map(|idx| {
                            let offset = molecules[..source_idx]
                                .iter()
                                .map(|input| input.atoms().count())
                                .sum::<usize>();
                            AtomId::from(idx + offset)
                        })
                        .collect::<Vec<_>>(),
                    combined.atoms().count(),
                ),
                Correspondence::from_images(
                    &(0..molecule.bonds().count())
                        .map(|idx| {
                            let offset = molecules[..source_idx]
                                .iter()
                                .map(|input| input.bonds().count())
                                .sum::<usize>();
                            BondId::from(idx + offset)
                        })
                        .collect::<Vec<_>>(),
                    combined.bonds().count(),
                ),
                Correspondence::from_images(
                    &(0..molecule.dative_bonds().count())
                        .map(|idx| {
                            let offset = molecules[..source_idx]
                                .iter()
                                .map(|input| input.dative_bonds().count())
                                .sum::<usize>();
                            DativeBondId::from(idx + offset)
                        })
                        .collect::<Vec<_>>(),
                    combined.dative_bonds().count(),
                ),
                Correspondence::from_images(
                    &(0..molecule.aromatic_systems().count())
                        .map(|idx| {
                            let offset = molecules[..source_idx]
                                .iter()
                                .map(|input| input.aromatic_systems().count())
                                .sum::<usize>();
                            AromaticSystemId::from(idx + offset)
                        })
                        .collect::<Vec<_>>(),
                    combined.aromatic_systems().count(),
                ),
                Correspondence::from_images(
                    &(0..molecule.multicenter_bonds().count())
                        .map(|idx| {
                            let offset = molecules[..source_idx]
                                .iter()
                                .map(|input| input.multicenter_bonds().count())
                                .sum::<usize>();
                            MulticenterBondId::from(idx + offset)
                        })
                        .collect::<Vec<_>>(),
                    combined.multicenter_bonds().count(),
                ),
                Correspondence::from_images(
                    &(0..molecule.noncovalent_bonds().count())
                        .map(|idx| {
                            let offset = molecules[..source_idx]
                                .iter()
                                .map(|input| input.noncovalent_bonds().count())
                                .sum::<usize>();
                            NoncovalentBondId::from(idx + offset)
                        })
                        .collect::<Vec<_>>(),
                    combined.noncovalent_bonds().count(),
                ),
                Correspondence::from_images(
                    &(0..molecule.stereo_atoms().count())
                        .map(|idx| {
                            let offset = molecules[..source_idx]
                                .iter()
                                .map(|input| input.stereo_atoms().count())
                                .sum::<usize>();
                            StereoAtomId::from(idx + offset)
                        })
                        .collect::<Vec<_>>(),
                    combined.stereo_atoms().count(),
                ),
                Correspondence::from_images(
                    &(0..molecule.stereo_bonds().count())
                        .map(|idx| {
                            let offset = molecules[..source_idx]
                                .iter()
                                .map(|input| input.stereo_bonds().count())
                                .sum::<usize>();
                            StereoBondId::from(idx + offset)
                        })
                        .collect::<Vec<_>>(),
                    combined.stereo_bonds().count(),
                ),
            );
            prop_assert_eq!(&combined.extract(&correspondence), molecule);
        }
    }

    #[test]
    fn test_molecule_combine(
        left in molecule_structurally_unambiguous_strategy(),
        right in molecule_structurally_unambiguous_strategy(),
    ) {
        let combined = left.combine(&right);
        let correspondence = MoleculeCorrespondence::new(
            Correspondence::from_images(
                &(0..right.atoms().count())
                    .map(|idx| AtomId::from(idx + left.atoms().count()))
                    .collect::<Vec<_>>(),
                combined.atoms().count(),
            ),
            Correspondence::from_images(
                &(0..right.bonds().count())
                    .map(|idx| BondId::from(idx + left.bonds().count()))
                    .collect::<Vec<_>>(),
                combined.bonds().count(),
            ),
            Correspondence::from_images(
                &(0..right.dative_bonds().count())
                    .map(|idx| DativeBondId::from(idx + left.dative_bonds().count()))
                    .collect::<Vec<_>>(),
                combined.dative_bonds().count(),
            ),
            Correspondence::from_images(
                &(0..right.aromatic_systems().count())
                    .map(|idx| AromaticSystemId::from(idx + left.aromatic_systems().count()))
                    .collect::<Vec<_>>(),
                combined.aromatic_systems().count(),
            ),
            Correspondence::from_images(
                &(0..right.multicenter_bonds().count())
                    .map(|idx| MulticenterBondId::from(idx + left.multicenter_bonds().count()))
                    .collect::<Vec<_>>(),
                combined.multicenter_bonds().count(),
            ),
            Correspondence::from_images(
                &(0..right.noncovalent_bonds().count())
                    .map(|idx| NoncovalentBondId::from(idx + left.noncovalent_bonds().count()))
                    .collect::<Vec<_>>(),
                combined.noncovalent_bonds().count(),
            ),
            Correspondence::from_images(
                &(0..right.stereo_atoms().count())
                    .map(|idx| StereoAtomId::from(idx + left.stereo_atoms().count()))
                    .collect::<Vec<_>>(),
                combined.stereo_atoms().count(),
            ),
            Correspondence::from_images(
                &(0..right.stereo_bonds().count())
                    .map(|idx| StereoBondId::from(idx + left.stereo_bonds().count()))
                    .collect::<Vec<_>>(),
                combined.stereo_bonds().count(),
            ),
        );
        prop_assert_eq!(combined.extract(&correspondence), right);
    }

    #[test]
    fn test_molecule_combine_from(
        left in molecule_structurally_unambiguous_strategy(),
        right in molecule_structurally_unambiguous_strategy(),
    ) {
        let expected = left.combine(&right);
        let mut combined = left;
        combined.combine_from(&right);

        prop_assert_eq!(combined, expected);
    }

    #[test]
    fn test_molecule_split(molecule in molecule_structurally_unambiguous_strategy()) {
        let components = molecule.tracked_split();
        prop_assert_eq!(
            molecule.split(),
            components.iter().map(|(component, _)| component.clone()).collect::<Vec<_>>(),
        );
        let mut covered_atoms = Vec::new();

        for (component, correspondence) in &components {
            prop_assert_eq!(&molecule.extract(&correspondence.reverse()), component);
            covered_atoms.extend(
                correspondence
                    .atoms()
                    .matched_pairs()
                    .iter()
                    .map(|&(source, _)| source),
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
