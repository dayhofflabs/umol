//! Stereo graph IR semantic properties.

use std::collections::BTreeSet;
use std::iter;

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{ConstitutionColoring, GraphSymmetryConfig};

use crate::strategies::*;

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
    fn test_stereo_symmetry_stereogenicity_agrees_with_coset_action(
        elements in prop::collection::vec(element_strategy(), 4),
        coset in 0u32..2,
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
            symmetry
                .kind()
                .class_key()
                .space()
                .reindex(coset, operation.permutation())
                == Some(coset)
        });

        prop_assert_eq!(symmetry.is_stereogenic(), expected);
    }

    #[test]
    fn test_molecule_try_from_entries_rejects_stereo_coset_out_of_range(
        elements in prop::collection::vec(element_strategy(), 4),
        coset in 2u32..=32,
    ) {
        let result = Molecule::try_from_entries(MoleculeEntries {
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

        prop_assert_eq!(
            result,
            Err(MoleculeIntegrityError::StereoCosetOutOfRange {
                entity: Entity::StereoAtom(StereoAtomId(0)),
                kind: StereoKind::Tetrahedral,
                coset,
                count: 2,
            }),
        );
    }

    /// A ligand frame that repeats a ligand is a valid arrangement record — distinctness is
    /// neither required nor asserted, and the site being non-stereogenic does not make the stored
    /// coset vacuous. Frame selection must therefore accept it, and the value it selects must not
    /// depend on which presentation of the arrangement it started from.
    #[test]
    fn test_stereo_atom_form_select_frame_repeated_ligands(
        (form, frame, action) in repeated_ligand_stereo_atom_strategy(),
    ) {
        prop_assert!(
            frame.len() > frame.iter().collect::<BTreeSet<_>>().len(),
            "the generator must repeat a ligand",
        );

        let selected = form
            .select_frame(&frame)
            .expect("a repeated ligand frame is a valid arrangement record");
        let canonical = form.clone().reframe_by(selected);
        prop_assert!(canonical.is_some());

        let restated = form.reframe_by(action).expect("a parent-group action is admissible");
        let restated_frame = action.act(&frame);
        let restated_selected = restated
            .select_frame(&restated_frame)
            .expect("a repeated ligand frame is a valid arrangement record");
        prop_assert_eq!(restated.clone().reframe_by(restated_selected), canonical);

        let converged = restated
            .reframe_by(restated_selected)
            .expect("selection produced an admissible action");
        let converged_frame = restated_selected.act(&restated_frame);
        let again = converged
            .select_frame(&converged_frame)
            .expect("a repeated ligand frame is a valid arrangement record");
        prop_assert_eq!(converged.clone().reframe_by(again), Some(converged));
    }

    #[test]
    fn test_stereo_atom_form_reframe_to(
        args in stereo_atom_form_strategy().prop_flat_map(|form| {
            let kind = form.configuration.kind().expect("strategy generates a kinded form");
            (Just(form), stereo_frame_permutation_strategy(kind))
        }),
    ) {
        let (form, permutation) = args;
        let kind = form.configuration.kind().expect("strategy generates a kinded form");
        let before: Vec<StereoLigand> = (0..kind.degree() as u32)
            .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
            .collect();
        let after = permutation.act(&before);
        let transformed = form.clone().reframe_to(&before, &after);

        prop_assert_eq!(transformed, form.clone().reframe_by(permutation));
        prop_assert_eq!(
            form.clone()
                .reframe_by(permutation)
                .and_then(|form| form.reframe_by(permutation.inverse())),
            Some(form),
        );
    }

    #[test]
    fn test_stereo_bond_form_reframe_to(
        args in stereo_bond_form_strategy().prop_flat_map(|form| (
            Just(form),
            stereo_frame_permutation_strategy(StereoKind::CisTrans),
        )),
    ) {
        let (form, permutation) = args;
        let before: Vec<StereoLigand> = (0..StereoKind::CisTrans.degree() as u32)
            .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
            .collect();
        let after = permutation.act(&before);
        let transformed = form.clone().reframe_to(&before, &after);

        prop_assert_eq!(transformed, form.clone().reframe_by(permutation));
        prop_assert_eq!(
            form.clone()
                .reframe_by(permutation)
                .and_then(|form| form.reframe_by(permutation.inverse())),
            Some(form),
        );
    }
}
