//! Molecule meet-pushout properties.
//!
//! Besides the operation laws, disjoint generated inputs exercise composition from each input
//! through the pushout to every split component. The composed witnesses must induce all entity
//! correspondences from their atom components.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_core::{Correspondence, EdgeId, GraphCorrespondence, NodeId};
use umol_graph_ir::ir::MoleculePushoutCorrespondence;

use crate::strategies::*;

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]
    /// Disjoint gluing agrees with concatenation and maps every input entity to its append position.
    #[test]
    fn test_molecule_tracked_meet_pushout_disjoint(
        left in molecule_structurally_unambiguous_strategy(),
        right in molecule_structurally_unambiguous_strategy(),
    ) {
        let overlap = GraphCorrespondence::new(
            Correspondence::new(vec![], left.atoms().count(), right.atoms().count()).unwrap(),
            Correspondence::new(vec![], left.bonds().count(), right.bonds().count()).unwrap(),
        );
        let expected = left.combine(&right);
        let expected_left = MoleculeCorrespondence::new(
            Correspondence::from_images(
                &(0..left.atoms().count())
                    .map(AtomId::from)
                    .collect::<Vec<_>>(),
                expected.atoms().count(),
            ),
            Correspondence::from_images(
                &(0..left.bonds().count())
                    .map(BondId::from)
                    .collect::<Vec<_>>(),
                expected.bonds().count(),
            ),
            Correspondence::from_images(
                &(0..left.dative_bonds().count())
                    .map(DativeBondId::from)
                    .collect::<Vec<_>>(),
                expected.dative_bonds().count(),
            ),
            Correspondence::from_images(
                &(0..left.aromatic_systems().count())
                    .map(AromaticSystemId::from)
                    .collect::<Vec<_>>(),
                expected.aromatic_systems().count(),
            ),
            Correspondence::from_images(
                &(0..left.multicenter_bonds().count())
                    .map(MulticenterBondId::from)
                    .collect::<Vec<_>>(),
                expected.multicenter_bonds().count(),
            ),
            Correspondence::from_images(
                &(0..left.noncovalent_bonds().count())
                    .map(NoncovalentBondId::from)
                    .collect::<Vec<_>>(),
                expected.noncovalent_bonds().count(),
            ),
            Correspondence::from_images(
                &(0..left.stereo_atoms().count())
                    .map(StereoAtomId::from)
                    .collect::<Vec<_>>(),
                expected.stereo_atoms().count(),
            ),
            Correspondence::from_images(
                &(0..left.stereo_bonds().count())
                    .map(StereoBondId::from)
                    .collect::<Vec<_>>(),
                expected.stereo_bonds().count(),
            ),
        );
        let expected_right = MoleculeCorrespondence::new(
            Correspondence::from_images(
                &(0..right.atoms().count())
                    .map(|idx| AtomId::from(idx + left.atoms().count()))
                    .collect::<Vec<_>>(),
                expected.atoms().count(),
            ),
            Correspondence::from_images(
                &(0..right.bonds().count())
                    .map(|idx| BondId::from(idx + left.bonds().count()))
                    .collect::<Vec<_>>(),
                expected.bonds().count(),
            ),
            Correspondence::from_images(
                &(0..right.dative_bonds().count())
                    .map(|idx| DativeBondId::from(idx + left.dative_bonds().count()))
                    .collect::<Vec<_>>(),
                expected.dative_bonds().count(),
            ),
            Correspondence::from_images(
                &(0..right.aromatic_systems().count())
                    .map(|idx| AromaticSystemId::from(idx + left.aromatic_systems().count()))
                    .collect::<Vec<_>>(),
                expected.aromatic_systems().count(),
            ),
            Correspondence::from_images(
                &(0..right.multicenter_bonds().count())
                    .map(|idx| MulticenterBondId::from(idx + left.multicenter_bonds().count()))
                    .collect::<Vec<_>>(),
                expected.multicenter_bonds().count(),
            ),
            Correspondence::from_images(
                &(0..right.noncovalent_bonds().count())
                    .map(|idx| NoncovalentBondId::from(idx + left.noncovalent_bonds().count()))
                    .collect::<Vec<_>>(),
                expected.noncovalent_bonds().count(),
            ),
            Correspondence::from_images(
                &(0..right.stereo_atoms().count())
                    .map(|idx| StereoAtomId::from(idx + left.stereo_atoms().count()))
                    .collect::<Vec<_>>(),
                expected.stereo_atoms().count(),
            ),
            Correspondence::from_images(
                &(0..right.stereo_bonds().count())
                    .map(|idx| StereoBondId::from(idx + left.stereo_bonds().count()))
                    .collect::<Vec<_>>(),
                expected.stereo_bonds().count(),
            ),
        );

        prop_assert_eq!(left.meet_pushout(&right, &overlap), Some(expected.clone()));
        prop_assert_eq!(
            left.tracked_meet_pushout(&right, &overlap),
            Some((expected, MoleculePushoutCorrespondence {
                left: expected_left,
                right: expected_right,
            })),
        );
    }

    #[test]
    fn test_molecule_tracked_meet_pushout_split_composition(
        left in molecule_structurally_unambiguous_strategy(),
        right in molecule_structurally_unambiguous_strategy(),
    ) {
        prop_assume!(left.atoms().count() + right.atoms().count() > 0);
        let overlap = GraphCorrespondence::new(
            Correspondence::new(vec![], left.atoms().count(), right.atoms().count()).unwrap(),
            Correspondence::new(vec![], left.bonds().count(), right.bonds().count()).unwrap(),
        );
        let (pushout, witness) = left
            .tracked_meet_pushout(&right, &overlap)
            .expect("disjoint molecule gluing is admissible");

        for (component, split) in pushout.tracked_split() {
            let left_to_component = witness.left.compose(&split).unwrap();
            let right_to_component = witness.right.compose(&split).unwrap();

            prop_assert!(left_to_component.is_compatible(&left, &component));
            prop_assert!(right_to_component.is_compatible(&right, &component));
            prop_assert_eq!(
                &left_to_component,
                &MoleculeCorrespondence::induce(
                    &left,
                    &component,
                    left_to_component.atoms().clone(),
                )
                .expect("disjoint input incidence remains unique in each split component"),
            );
            prop_assert_eq!(
                &right_to_component,
                &MoleculeCorrespondence::induce(
                    &right,
                    &component,
                    right_to_component.atoms().clone(),
                )
                .expect("disjoint input incidence remains unique in each split component"),
            );
        }
    }

    #[test]
    fn test_molecule_meet_pushout_stereo_atom_frame(
        coset in 0..StereoKind::Tetrahedral.count() as u32,
        permutation in stereo_frame_permutation_strategy(StereoKind::Tetrahedral),
    ) {
        let atoms = vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
            AtomForm::from_element(Element::I),
            AtomForm::from_element(Element::N),
        ];
        let bonds: Vec<(AtomId, AtomId, BondForm)> = (1..=4)
            .map(|ligand| (AtomId(0), AtomId(ligand), BondForm::from_order(1)))
            .collect();
        let left_frame: Vec<StereoLigand> = (1..=4)
            .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
            .collect();
        let left_form = StereoAtomForm::new(StereoKind::Tetrahedral, coset);
        let left = Molecule::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_atoms: vec![(AtomId(0), left_frame.clone(), left_form.clone())],
            ..Default::default()
        });
        let right = Molecule::from_entries(MoleculeEntries {
            atoms,
            bonds,
            stereo_atoms: vec![(
                AtomId(0),
                permutation.act(&left_frame),
                left_form.apply(permutation).expect("the permutation is a parent-group action of the form's kind"),
            )],
            ..Default::default()
        });
        let overlap = GraphCorrespondence::new(
            Correspondence::from_images(&(0..6u32).map(NodeId).collect::<Vec<_>>(), 6),
            Correspondence::from_images(&(0..4u32).map(EdgeId).collect::<Vec<_>>(), 4),
        );

        prop_assert_eq!(
            left.meet_pushout(&right, &overlap),
            Some(left.clone()),
        );
        prop_assert_eq!(
            left.tracked_meet_pushout(&right, &overlap).map(|(object, _)| object),
            Some(left),
        );
    }

    #[test]
    fn test_molecule_meet_pushout_stereo_atom_ligand(
        coset in 0..StereoKind::Tetrahedral.count() as u32,
        permutation in stereo_frame_permutation_strategy(StereoKind::Tetrahedral),
    ) {
        let atoms = vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
            AtomForm::from_element(Element::I),
            AtomForm::from_element(Element::N),
        ];
        let bonds: Vec<(AtomId, AtomId, BondForm)> = (1..=4)
            .map(|ligand| (AtomId(0), AtomId(ligand), BondForm::from_order(1)))
            .chain([(AtomId(0), AtomId(5), BondForm::from_order(1))])
            .collect();
        let left_frame: Vec<StereoLigand> = (1..=4)
            .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
            .collect();
        let left_form = StereoAtomForm::new(StereoKind::Tetrahedral, coset);
        let left = Molecule::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_atoms: vec![(AtomId(0), left_frame.clone(), left_form.clone())],
            ..Default::default()
        });
        let mut right_frame = permutation.act(&left_frame);
        right_frame[0] = StereoLigand::new(AtomId(5), StereoLigandKind::Atom);
        let right = Molecule::from_entries(MoleculeEntries {
            atoms,
            bonds,
            stereo_atoms: vec![(AtomId(0), right_frame, left_form.apply(permutation).expect("the permutation is a parent-group action of the form's kind"))],
            ..Default::default()
        });
        let overlap = GraphCorrespondence::new(
            Correspondence::from_images(&(0..6u32).map(NodeId).collect::<Vec<_>>(), 6),
            Correspondence::from_images(&(0..5u32).map(EdgeId).collect::<Vec<_>>(), 5),
        );

        prop_assert_eq!(left.meet_pushout(&right, &overlap), None);
        prop_assert_eq!(left.tracked_meet_pushout(&right, &overlap), None);
    }

    #[test]
    fn test_molecule_meet_pushout_stereo_bond_frame(
        coset in 0..StereoKind::CisTrans.count() as u32,
        permutation in stereo_frame_permutation_strategy(StereoKind::CisTrans),
    ) {
        let atoms = vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
            AtomForm::from_element(Element::I),
        ];
        let bonds = vec![
            (AtomId(0), AtomId(1), BondForm::from_order(2)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(1), AtomId(4), BondForm::from_order(1)),
            (AtomId(1), AtomId(5), BondForm::from_order(1)),
        ];
        let left_frame: Vec<StereoLigand> = (2..=5)
            .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
            .collect();
        let left_form = StereoBondForm::new(StereoKind::CisTrans, coset);
        let left = Molecule::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_bonds: vec![(BondId(0), left_frame.clone(), left_form.clone())],
            ..Default::default()
        });
        let right = Molecule::from_entries(MoleculeEntries {
            atoms,
            bonds,
            stereo_bonds: vec![(
                BondId(0),
                permutation.act(&left_frame),
                left_form.apply(permutation).expect("the permutation is a parent-group action of the form's kind"),
            )],
            ..Default::default()
        });
        let overlap = GraphCorrespondence::new(
            Correspondence::from_images(&(0..6u32).map(NodeId).collect::<Vec<_>>(), 6),
            Correspondence::from_images(&(0..5u32).map(EdgeId).collect::<Vec<_>>(), 5),
        );

        prop_assert_eq!(
            left.meet_pushout(&right, &overlap),
            Some(left.clone()),
        );
        prop_assert_eq!(
            left.tracked_meet_pushout(&right, &overlap).map(|(object, _)| object),
            Some(left),
        );
    }

    #[test]
    fn test_molecule_meet_pushout_stereo_bond_ligand(
        coset in 0..StereoKind::CisTrans.count() as u32,
        permutation in stereo_frame_permutation_strategy(StereoKind::CisTrans),
    ) {
        let atoms = vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
            AtomForm::from_element(Element::I),
            AtomForm::from_element(Element::N),
        ];
        let bonds = vec![
            (AtomId(0), AtomId(1), BondForm::from_order(2)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(1), AtomId(4), BondForm::from_order(1)),
            (AtomId(1), AtomId(5), BondForm::from_order(1)),
            (AtomId(0), AtomId(6), BondForm::from_order(1)),
            (AtomId(1), AtomId(6), BondForm::from_order(1)),
        ];
        let left_frame: Vec<StereoLigand> = (2..=5)
            .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
            .collect();
        let left_form = StereoBondForm::new(StereoKind::CisTrans, coset);
        let left = Molecule::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_bonds: vec![(BondId(0), left_frame.clone(), left_form.clone())],
            ..Default::default()
        });
        let mut right_frame = permutation.act(&left_frame);
        right_frame[0] = StereoLigand::new(AtomId(6), StereoLigandKind::Atom);
        let right = Molecule::from_entries(MoleculeEntries {
            atoms,
            bonds,
            stereo_bonds: vec![(BondId(0), right_frame, left_form.apply(permutation).expect("the permutation is a parent-group action of the form's kind"))],
            ..Default::default()
        });
        let overlap = GraphCorrespondence::new(
            Correspondence::from_images(&(0..7u32).map(NodeId).collect::<Vec<_>>(), 7),
            Correspondence::from_images(&(0..7u32).map(EdgeId).collect::<Vec<_>>(), 7),
        );

        prop_assert_eq!(left.meet_pushout(&right, &overlap), None);
        prop_assert_eq!(left.tracked_meet_pushout(&right, &overlap), None);
    }
}
