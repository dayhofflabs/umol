//! Stereo AST semantic properties.

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
    fn test_ligand_pair_ast_normalization(a in 0u32..6, b in 0u32..6) {
        let pair = StereoLigandPair::new(StereoLigandPosition(a), StereoLigandPosition(b));
        prop_assert!(pair.first().0 <= pair.second().0);
        prop_assert_eq!(pair, StereoLigandPair::new(StereoLigandPosition(b), StereoLigandPosition(a)));
    }

    /// Concrete (non-lattice) literal `matches` is exactly equality.
    #[test]
    fn test_permutation_ast_matches_is_equality(
        (a, b) in (2usize..=6).prop_flat_map(|d| (permutation_strategy(d), permutation_strategy(d))),
    ) {
        let (x, y) = (LigandPermutation(a), LigandPermutation(b));
        prop_assert!(x.matches(&x));
        prop_assert_eq!(x.matches(&y), x == y);
    }

    #[test]
    fn test_ligand_symmetry_ast_matches_is_equality(
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
        let molecule = MoleculeAst::from_entries(MoleculeEntries {
            atoms: iter::once(AtomAst::from_element(Element::C))
                .chain(elements.into_iter().map(AtomAst::from_element))
                .collect(),
            bonds: (1..=4)
                .map(|ligand| (AtomId(0), AtomId(ligand), BondAst::from_order(1)))
                .collect(),
            stereo_atoms: vec![(
                AtomId(0),
                (1..=4)
                    .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomAst::new(StereoKind::Tetrahedral, coset),
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
    fn test_stereo_symmetry_malformed_coset_is_not_stereogenic(
        elements in prop::collection::vec(element_strategy(), 4),
        coset in 2u32..=32,
    ) {
        let molecule = MoleculeAst::from_entries(MoleculeEntries {
            atoms: iter::once(AtomAst::from_element(Element::C))
                .chain(elements.into_iter().map(AtomAst::from_element))
                .collect(),
            bonds: (1..=4)
                .map(|ligand| (AtomId(0), AtomId(ligand), BondAst::from_order(1)))
                .collect(),
            stereo_atoms: vec![(
                AtomId(0),
                (1..=4)
                    .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomAst::new(StereoKind::Tetrahedral, coset),
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

        prop_assert!(!symmetry.is_stereogenic());
    }

    #[test]
    fn test_stereo_atom_ast_transform_frame(
        args in stereo_atom_kind_strategy().prop_flat_map(|kind| {
            (
                Just(kind),
                0..kind.count() as u32,
                stereo_frame_permutation_strategy(kind),
            )
        }),
    ) {
        let (kind, coset, permutation) = args;
        let before: Vec<StereoLigand> = (0..kind.degree() as u32)
            .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
            .collect();
        let after = permutation.act(&before);
        let ast = StereoAtomAst::new(kind, coset);
        let transformed = ast.transform_frame(&before, &after);

        prop_assert_eq!(transformed.as_ref(), Some(&ast.apply(permutation)));
        prop_assert_eq!(
            transformed.and_then(|ast| ast.transform_frame(&after, &before)),
            Some(ast),
        );
    }

    #[test]
    fn test_stereo_bond_ast_transform_frame(
        coset in 0..StereoKind::CisTrans.count() as u32,
        permutation in stereo_frame_permutation_strategy(StereoKind::CisTrans),
    ) {
        let before: Vec<StereoLigand> = (0..StereoKind::CisTrans.degree() as u32)
            .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
            .collect();
        let after = permutation.act(&before);
        let ast = StereoBondAst::new(StereoKind::CisTrans, coset);
        let transformed = ast.transform_frame(&before, &after);

        prop_assert_eq!(transformed.as_ref(), Some(&ast.apply(permutation)));
        prop_assert_eq!(
            transformed.and_then(|ast| ast.transform_frame(&after, &before)),
            Some(ast),
        );
    }
}
