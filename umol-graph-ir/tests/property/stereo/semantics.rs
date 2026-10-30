//! Stereo graph IR semantic properties.

use std::{array, iter};

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_chem::element::Element;
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{
    AtomForm, AtomId, BondForm, ConstitutionColoring, GraphSymmetryConfig, Lattice,
    LigandPermutation, Molecule, MoleculeEntries, StereoAtomForm, StereoAtomId,
    StereoConfigurationForm, StereoCoset, StereoKind, StereoLigand, StereoLigandKind,
    StereoLigandPair, StereoLigandPosition,
};

use crate::strategies::{element_strategy, ligand_symmetry_strategy, permutation_strategy};

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]
    /// `StereoLigandPair::new` normalizes to `first <= second` and is symmetric.
    #[test]
    fn test_ligand_pair_normalization(a in 0u32..6, b in 0u32..6) {
        let pair = StereoLigandPair::new(StereoLigandPosition(a), StereoLigandPosition(b));
        prop_assert!(pair.first().0 <= pair.second().0);
        prop_assert_eq!(pair, StereoLigandPair::new(StereoLigandPosition(b), StereoLigandPosition(a)));
    }

    /// Concrete (non-lattice) literal `matches` is exactly equality.
    #[test]
    fn test_permutation_matches_is_equality(
        (a, b) in (2usize..=6).prop_flat_map(|d| (permutation_strategy(d), permutation_strategy(d))),
    ) {
        let (x, y) = (LigandPermutation(a), LigandPermutation(b));
        prop_assert!(x.matches(&x));
        prop_assert_eq!(x.matches(&y), x == y);
    }

    #[test]
    fn test_ligand_symmetry_form_matches_is_equality(
        (x, y) in (2usize..=6)
            .prop_flat_map(|d| (ligand_symmetry_strategy(d), ligand_symmetry_strategy(d))),
    ) {
        prop_assert!(x.matches(&x));
        prop_assert_eq!(x.matches(&y), x == y);
    }

    #[test]
    fn test_stereo_symmetry_is_stereogenic(
        elements in prop::collection::vec(element_strategy(), 4),
        coset in prop_oneof![
            Just(StereoCoset::Undetermined),
            (0u32..2).prop_map(StereoCoset::Lit),
            Just(StereoCoset::lit_set([0, 1])),
        ],
    ) {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: iter::once(AtomForm::from_element(Element::C))
                .chain(elements.into_iter().map(AtomForm::from_element))
                .collect(),
            bonds: (1..=4)
                .map(|ligand| (AtomId(0), AtomId(ligand), BondForm::from_order(1)))
                .collect(),
            stereo_atoms: vec![(
                AtomId(0),
                (1..=4)
                    .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, coset),
            )],
            ..Default::default()
        });
        let graph = molecule.graph_symmetry(&GraphSymmetryConfig {
            coloring: ConstitutionColoring::full(),
            iterate_to_fixpoint: true,
            max_iterations: 16,
            automorphism_algorithm: AutomorphismAlgorithm::Nauty,
        });
        let symmetry = molecule.stereo_atom_symmetry(&graph, StereoAtomId(0));
        let expected = symmetry.group().elements().iter().all(|operation| {
            let image = array::from_fn(|position| operation.permutation().apply(position));
            (0..2).all(|index| {
                super::reference::binary_coset_action(StereoKind::Tetrahedral, index, &image) == Some(index)
            })
        });

        prop_assert_eq!(symmetry.is_stereogenic(), expected);
    }

    #[test]
    fn test_molecule_try_from_entries_preserves_stereo_coset(
        elements in prop::collection::vec(element_strategy(), 4),
        coset in 2u32..=32,
    ) {
        let molecule = Molecule::try_from_entries(MoleculeEntries {
            atoms: iter::once(AtomForm::from_element(Element::C))
                .chain(elements.into_iter().map(AtomForm::from_element))
                .collect(),
            bonds: (1..=4)
                .map(|ligand| (AtomId(0), AtomId(ligand), BondForm::from_order(1)))
                .collect(),
            stereo_atoms: vec![(
                AtomId(0),
                (1..=4)
                    .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, coset),
            )],
            ..Default::default()
        }).unwrap();

        prop_assert_eq!(
            &molecule.stereo_atom(StereoAtomId(0)).attributes().configuration,
            &StereoConfigurationForm::kinded(StereoKind::Tetrahedral, coset),
        );
    }

}
