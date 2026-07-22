//! Molecule meet-pushout properties.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_core::{Correspondence, EdgeId, GraphCorrespondence, NodeId};

use crate::strategies::*;

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]
    #[test]
    fn test_molecule_ast_meet_pushout_reframes_stereo_atom(
        coset in 0..StereoKind::Tetrahedral.count() as u32,
        permutation in stereo_frame_permutation_strategy(StereoKind::Tetrahedral),
    ) {
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::F),
            AtomAst::from_element(Element::Cl),
            AtomAst::from_element(Element::Br),
            AtomAst::from_element(Element::I),
            AtomAst::from_element(Element::N),
        ];
        let bonds: Vec<(AtomId, AtomId, BondAst)> = (1..=4)
            .map(|ligand| (AtomId(0), AtomId(ligand), BondAst::from_order(1)))
            .collect();
        let left_frame: Vec<StereoLigand> = (1..=4)
            .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
            .collect();
        let left_ast = StereoAtomAst::new(StereoKind::Tetrahedral, coset);
        let left = MoleculeAst::from_parts(MoleculeParts {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_atoms: vec![(AtomId(0), left_frame.clone(), left_ast.clone())],
            ..Default::default()
        });
        let right = MoleculeAst::from_parts(MoleculeParts {
            atoms,
            bonds,
            stereo_atoms: vec![(
                AtomId(0),
                permutation.act(&left_frame),
                left_ast.apply(permutation),
            )],
            ..Default::default()
        });
        let overlap = GraphCorrespondence::new(
            Correspondence::from_images(&(0..6u32).map(NodeId).collect::<Vec<_>>(), 6),
            Correspondence::from_images(&(0..4u32).map(EdgeId).collect::<Vec<_>>(), 4),
        );

        prop_assert_eq!(
            left.meet_pushout(&right, &overlap).map(|pushout| pushout.object),
            Some(left),
        );
    }

    #[test]
    fn test_molecule_ast_meet_pushout_rejects_changed_stereo_atom_ligand(
        coset in 0..StereoKind::Tetrahedral.count() as u32,
        permutation in stereo_frame_permutation_strategy(StereoKind::Tetrahedral),
    ) {
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::F),
            AtomAst::from_element(Element::Cl),
            AtomAst::from_element(Element::Br),
            AtomAst::from_element(Element::I),
            AtomAst::from_element(Element::N),
        ];
        let bonds: Vec<(AtomId, AtomId, BondAst)> = (1..=4)
            .map(|ligand| (AtomId(0), AtomId(ligand), BondAst::from_order(1)))
            .collect();
        let left_frame: Vec<StereoLigand> = (1..=4)
            .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
            .collect();
        let left_ast = StereoAtomAst::new(StereoKind::Tetrahedral, coset);
        let left = MoleculeAst::from_parts(MoleculeParts {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_atoms: vec![(AtomId(0), left_frame.clone(), left_ast.clone())],
            ..Default::default()
        });
        let mut right_frame = permutation.act(&left_frame);
        right_frame[0] = StereoLigand::new(AtomId(5), StereoLigandKind::Atom);
        let right = MoleculeAst::from_parts(MoleculeParts {
            atoms,
            bonds,
            stereo_atoms: vec![(AtomId(0), right_frame, left_ast.apply(permutation))],
            ..Default::default()
        });
        let overlap = GraphCorrespondence::new(
            Correspondence::from_images(&(0..6u32).map(NodeId).collect::<Vec<_>>(), 6),
            Correspondence::from_images(&(0..4u32).map(EdgeId).collect::<Vec<_>>(), 4),
        );

        prop_assert!(left.meet_pushout(&right, &overlap).is_none());
    }

    #[test]
    fn test_molecule_ast_meet_pushout_reframes_stereo_bond(
        coset in 0..StereoKind::CisTrans.count() as u32,
        permutation in stereo_frame_permutation_strategy(StereoKind::CisTrans),
    ) {
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::F),
            AtomAst::from_element(Element::Cl),
            AtomAst::from_element(Element::Br),
            AtomAst::from_element(Element::I),
        ];
        let bonds = vec![
            (AtomId(0), AtomId(1), BondAst::from_order(2)),
            (AtomId(0), AtomId(2), BondAst::from_order(1)),
            (AtomId(0), AtomId(3), BondAst::from_order(1)),
            (AtomId(1), AtomId(4), BondAst::from_order(1)),
            (AtomId(1), AtomId(5), BondAst::from_order(1)),
        ];
        let left_frame: Vec<StereoLigand> = (2..=5)
            .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
            .collect();
        let left_ast = StereoBondAst::new(StereoKind::CisTrans, coset);
        let left = MoleculeAst::from_parts(MoleculeParts {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_bonds: vec![(BondId(0), left_frame.clone(), left_ast.clone())],
            ..Default::default()
        });
        let right = MoleculeAst::from_parts(MoleculeParts {
            atoms,
            bonds,
            stereo_bonds: vec![(
                BondId(0),
                permutation.act(&left_frame),
                left_ast.apply(permutation),
            )],
            ..Default::default()
        });
        let overlap = GraphCorrespondence::new(
            Correspondence::from_images(&(0..6u32).map(NodeId).collect::<Vec<_>>(), 6),
            Correspondence::from_images(&(0..5u32).map(EdgeId).collect::<Vec<_>>(), 5),
        );

        prop_assert_eq!(
            left.meet_pushout(&right, &overlap).map(|pushout| pushout.object),
            Some(left),
        );
    }

    #[test]
    fn test_molecule_ast_meet_pushout_rejects_changed_stereo_bond_ligand(
        coset in 0..StereoKind::CisTrans.count() as u32,
        permutation in stereo_frame_permutation_strategy(StereoKind::CisTrans),
    ) {
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::F),
            AtomAst::from_element(Element::Cl),
            AtomAst::from_element(Element::Br),
            AtomAst::from_element(Element::I),
            AtomAst::from_element(Element::N),
        ];
        let bonds = vec![
            (AtomId(0), AtomId(1), BondAst::from_order(2)),
            (AtomId(0), AtomId(2), BondAst::from_order(1)),
            (AtomId(0), AtomId(3), BondAst::from_order(1)),
            (AtomId(1), AtomId(4), BondAst::from_order(1)),
            (AtomId(1), AtomId(5), BondAst::from_order(1)),
        ];
        let left_frame: Vec<StereoLigand> = (2..=5)
            .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
            .collect();
        let left_ast = StereoBondAst::new(StereoKind::CisTrans, coset);
        let left = MoleculeAst::from_parts(MoleculeParts {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_bonds: vec![(BondId(0), left_frame.clone(), left_ast.clone())],
            ..Default::default()
        });
        let mut right_frame = permutation.act(&left_frame);
        right_frame[0] = StereoLigand::new(AtomId(6), StereoLigandKind::Atom);
        let right = MoleculeAst::from_parts(MoleculeParts {
            atoms,
            bonds,
            stereo_bonds: vec![(BondId(0), right_frame, left_ast.apply(permutation))],
            ..Default::default()
        });
        let overlap = GraphCorrespondence::new(
            Correspondence::from_images(&(0..7u32).map(NodeId).collect::<Vec<_>>(), 7),
            Correspondence::from_images(&(0..5u32).map(EdgeId).collect::<Vec<_>>(), 5),
        );

        prop_assert!(left.meet_pushout(&right, &overlap).is_none());
    }
}
