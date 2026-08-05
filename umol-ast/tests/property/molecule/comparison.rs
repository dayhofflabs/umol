//! Molecule comparison properties.
//!
//! The identity-frame laws for `equiv`, the correspondence laws for
//! `equiv_under`, and agreement with `==` on canonical ASTs deliberately use
//! overlapping molecule domains. They establish distinct relations: semantic
//! equivalence in a shared frame, semantic equivalence under an explicit frame
//! mapping, and a canonical-representation oracle, respectively.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_ast::ast::MoleculeCorrespondence;
use umol_graph_core::Correspondence;

use crate::strategies::*;

fn identity_correspondence(ast: &MoleculeAst) -> MoleculeCorrespondence {
    fn identity<Id>(count: usize) -> Correspondence<Id>
    where
        Id: Copy + Ord + From<usize>,
    {
        let images: Vec<Id> = (0..count).map(Id::from).collect();
        Correspondence::from_images(&images, count)
    }

    MoleculeCorrespondence::new(
        identity::<AtomId>(ast.atoms().count()),
        identity::<BondId>(ast.bonds().count()),
        identity::<DativeBondId>(ast.dative_bonds().count()),
        identity::<AromaticSystemId>(ast.aromatic_systems().count()),
        identity::<MulticenterBondId>(ast.multicenter_bonds().count()),
        identity::<NoncovalentBondId>(ast.noncovalent_bonds().count()),
        identity::<StereoAtomId>(ast.stereo_atoms().count()),
        identity::<StereoBondId>(ast.stereo_bonds().count()),
    )
}

fn atom_only_correspondence(images: &[AtomId], count: usize) -> MoleculeCorrespondence {
    fn empty<Id>() -> Correspondence<Id>
    where
        Id: Copy + Ord + From<usize>,
    {
        Correspondence::from_images(&[], 0)
    }
    MoleculeCorrespondence::new(
        Correspondence::from_images(images, count),
        empty(),
        empty(),
        empty(),
        empty(),
        empty(),
        empty(),
        empty(),
    )
}

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]
    #[test]
    fn test_molecule_ast_equiv_reflexive(ast in molecule_ast_with_constraints_strategy()) {
        prop_assert!(ast.equiv(&ast));
    }

    #[test]
    fn test_molecule_ast_equiv_symmetric(
        left in molecule_ast_with_constraints_strategy(),
        right in molecule_ast_with_constraints_strategy(),
    ) {
        prop_assert_eq!(left.equiv(&right), right.equiv(&left));
    }

    #[test]
    fn test_molecule_ast_equiv_under_transitive(
        atoms in prop::collection::vec(atom_ast_strategy(), 0..=5),
    ) {
        let count = atoms.len();
        let first_order = (0..count).collect::<Vec<_>>();
        let second_order = (0..count).rev().collect::<Vec<_>>();
        let mut third_order = first_order.clone();
        if count > 0 {
            third_order.rotate_left(1);
        }
        let molecule = |order: &[usize]| {
            MoleculeAst::from_parts(MoleculeParts {
                atoms: order.iter().map(|&index| atoms[index].clone()).collect(),
                ..Default::default()
            })
        };
        let correspondence = |left: &[usize], right: &[usize]| {
            let images = left
                .iter()
                .map(|original| {
                    AtomId::from(
                        right
                            .iter()
                            .position(|candidate| candidate == original)
                            .expect("orders contain the same indices"),
                    )
                })
                .collect::<Vec<_>>();
            atom_only_correspondence(&images, count)
        };

        let first = molecule(&first_order);
        let second = molecule(&second_order);
        let third = molecule(&third_order);
        let first_second = correspondence(&first_order, &second_order);
        let second_third = correspondence(&second_order, &third_order);

        prop_assert!(first.equiv_under(&second, &first_second));
        prop_assert!(second.equiv_under(&third, &second_third));
        prop_assert!(first.equiv_under(&third, &first_second.compose(&second_third)));
    }

    #[test]
    fn test_molecule_ast_equiv_agrees_with_equality_for_canonical_asts(
        left in molecule_ast_strategy(),
        right in molecule_ast_strategy(),
    ) {
        prop_assert_eq!(left.equiv(&right), left == right);
    }

    #[test]
    fn test_molecule_ast_equiv_under_identity_reduces_to_equiv(
        ast in molecule_ast_with_constraints_strategy(),
    ) {
        let correspondence = identity_correspondence(&ast);
        let mut other = ast.clone();
        if other.atoms().count() > 0 {
            other.atom_mut(AtomId(0)).ast.charge = ValueAst::Lit(99);
        }

        prop_assert_eq!(
            ast.equiv_under(&other, &correspondence),
            ast.equiv(&other),
        );
    }

    #[test]
    fn test_molecule_ast_equiv_under_symmetric_under_reverse(
        atoms in prop::collection::vec(atom_ast_strategy(), 0..=5),
        change_mapped_atom in any::<bool>(),
    ) {
        let count = atoms.len();
        let left = MoleculeAst::from_parts(MoleculeParts {
            atoms: atoms.clone(),
            ..Default::default()
        });
        let mut right = MoleculeAst::from_parts(MoleculeParts {
            atoms: atoms.into_iter().rev().collect(),
            ..Default::default()
        });
        if change_mapped_atom && count > 0 {
            right.atom_mut(AtomId((count - 1) as u32)).ast.charge = ValueAst::Lit(99);
        }
        let images: Vec<AtomId> = (0..count).rev().map(AtomId::from).collect();
        let correspondence = MoleculeCorrespondence::new(
            Correspondence::from_images(&images, count),
            Correspondence::from_images(&[], 0),
            Correspondence::from_images(&[], 0),
            Correspondence::from_images(&[], 0),
            Correspondence::from_images(&[], 0),
            Correspondence::from_images(&[], 0),
            Correspondence::from_images(&[], 0),
            Correspondence::from_images(&[], 0),
        );

        let forward = left.equiv_under(&right, &correspondence);
        let reverse = right.equiv_under(&left, &correspondence.reverse());
        prop_assert_eq!(forward, reverse);
        prop_assert_eq!(forward, !change_mapped_atom || count == 0);
    }
}
