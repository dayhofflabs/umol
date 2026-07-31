use proptest::prelude::*;
use umol_ast::ast::StereoKind;

use crate::strategies::*;

// Lattice-law sweep: every `impl Lattice` type satisfies commutativity,
// associativity, absorption, idempotence, and `matches`↔`meet` consistency,
// checked by `assert_lattice_laws` over a generated value triple.

proptest! {
    #[test]
    fn test_value_ast_lattice_laws(
        a in any_value_ast_strategy(),
        b in any_value_ast_strategy(),
        c in any_value_ast_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_boolean_ast_lattice_laws(
        a in prop_oneof![Just(BooleanAst::Undetermined), any::<bool>().prop_map(BooleanAst::Lit)],
        b in prop_oneof![Just(BooleanAst::Undetermined), any::<bool>().prop_map(BooleanAst::Lit)],
        c in prop_oneof![Just(BooleanAst::Undetermined), any::<bool>().prop_map(BooleanAst::Lit)],
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    /// The universal (input-canonicality-independent) laws on raw, non-canonical
    /// inputs — `canonical()` fold path, `equiv`, `matches`↔meet, and meet/join
    /// canonicality — which the canonicalized strategies above never reach. One
    /// per `canonical()`-overriding leaf with a fold path.
    #[test]
    fn test_value_ast_lattice_laws_raw(
        a in raw_value_ast_strategy(),
        b in raw_value_ast_strategy(),
        c in raw_value_ast_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_element_ast_lattice_laws_raw(
        a in raw_element_ast_strategy(),
        b in raw_element_ast_strategy(),
        c in raw_element_ast_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_isotope_mass_ast_lattice_laws_raw(
        a in raw_isotope_strategy(),
        b in raw_isotope_strategy(),
        c in raw_isotope_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_aromatic_valence_ast_lattice_laws_raw(
        a in raw_aromatic_valence_ast_strategy(),
        b in raw_aromatic_valence_ast_strategy(),
        c in raw_aromatic_valence_ast_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_multicenter_valence_ast_lattice_laws_raw(
        a in raw_multicenter_valence_ast_strategy(),
        b in raw_multicenter_valence_ast_strategy(),
        c in raw_multicenter_valence_ast_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_tetrahedral_stereo_ast_lattice_laws_raw(
        a in raw_tetrahedral_stereo_strategy(),
        b in raw_tetrahedral_stereo_strategy(),
        c in raw_tetrahedral_stereo_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_cis_trans_stereo_ast_lattice_laws_raw(
        a in raw_cis_trans_stereo_strategy(),
        b in raw_cis_trans_stereo_strategy(),
        c in raw_cis_trans_stereo_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_stereo_configuration_ast_lattice_laws_raw(
        a in raw_stereo_configuration_strategy(),
        b in raw_stereo_configuration_strategy(),
        c in raw_stereo_configuration_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_isotope_mass_ast_as_lit_laws(
        a in raw_isotope_strategy(),
        b in raw_isotope_strategy(),
    ) {
        prop_assert_eq!(a.is_ground(), a.as_lit().is_some());
        prop_assert_eq!(b.is_ground(), b.as_lit().is_some());
        let a = a.canonicalize().unwrap();
        let b = b.canonicalize().unwrap();
        if a.is_ground() && b.is_ground() && a != b {
            prop_assert_ne!(a.as_lit(), b.as_lit());
        }
    }

    #[test]
    fn test_aromatic_valence_ast_as_lit_laws(
        a in raw_aromatic_valence_ast_strategy(),
        b in raw_aromatic_valence_ast_strategy(),
    ) {
        prop_assert_eq!(a.is_ground(), a.as_lit().is_some());
        prop_assert_eq!(b.is_ground(), b.as_lit().is_some());
        let a = a.canonicalize().unwrap();
        let b = b.canonicalize().unwrap();
        if a.is_ground() && b.is_ground() && a != b {
            prop_assert_ne!(a.as_lit(), b.as_lit());
        }
    }

    #[test]
    fn test_aromatic_valence_ast_aromatic_covalence(
        valence in any::<i64>(),
    ) {
        let ast = AromaticValenceAst::from(AromaticValence::Aromatic(valence));
        let expected = if valence == 1 { 1 } else { 0 };
        prop_assert_eq!(aromatic_covalence(valence), expected);
        prop_assert_eq!(ast.aromatic_covalence(), ValueAst::Lit(expected));
    }

    #[test]
    fn test_multicenter_valence_ast_as_lit_laws(
        a in raw_multicenter_valence_ast_strategy(),
        b in raw_multicenter_valence_ast_strategy(),
    ) {
        prop_assert_eq!(a.is_ground(), a.as_lit().is_some());
        prop_assert_eq!(b.is_ground(), b.as_lit().is_some());
        let a = a.canonicalize().unwrap();
        let b = b.canonicalize().unwrap();
        if a.is_ground() && b.is_ground() && a != b {
            prop_assert_ne!(a.as_lit(), b.as_lit());
        }
    }

    #[test]
    fn test_tetrahedral_stereo_ast_as_lit_laws(
        a in raw_tetrahedral_stereo_strategy(),
        b in raw_tetrahedral_stereo_strategy(),
    ) {
        prop_assert_eq!(a.is_ground(), a.as_lit().is_some());
        prop_assert_eq!(b.is_ground(), b.as_lit().is_some());
        let a = a.canonicalize().unwrap();
        let b = b.canonicalize().unwrap();
        if a.is_ground() && b.is_ground() && a != b {
            prop_assert_ne!(a.as_lit(), b.as_lit());
        }
    }

    #[test]
    fn test_cis_trans_stereo_ast_as_lit_laws(
        a in raw_cis_trans_stereo_strategy(),
        b in raw_cis_trans_stereo_strategy(),
    ) {
        prop_assert_eq!(a.is_ground(), a.as_lit().is_some());
        prop_assert_eq!(b.is_ground(), b.as_lit().is_some());
        let a = a.canonicalize().unwrap();
        let b = b.canonicalize().unwrap();
        if a.is_ground() && b.is_ground() && a != b {
            prop_assert_ne!(a.as_lit(), b.as_lit());
        }
    }

    #[test]
    fn test_topicity_relation_ast_lattice_laws_raw(
        a in raw_topicity_relation_strategy(),
        b in raw_topicity_relation_strategy(),
        c in raw_topicity_relation_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_stereogenicity_ast_lattice_laws_raw(
        a in raw_stereogenicity_relation_strategy(),
        b in raw_stereogenicity_relation_strategy(),
        c in raw_stereogenicity_relation_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_element_ast_lattice_laws(
        a in element_ast_strategy(),
        b in element_ast_strategy(),
        c in element_ast_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_isotope_mass_ast_lattice_laws(
        a in isotope_strategy(),
        b in isotope_strategy(),
        c in isotope_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_aromatic_valence_ast_lattice_laws(
        a in aromatic_valence_ast_strategy(),
        b in aromatic_valence_ast_strategy(),
        c in aromatic_valence_ast_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_multicenter_valence_ast_lattice_laws(
        a in multicenter_valence_ast_strategy(),
        b in multicenter_valence_ast_strategy(),
        c in multicenter_valence_ast_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_noncovalent_bond_kind_ast_lattice_laws(
        a in noncovalent_bond_kind_ast_strategy(),
        b in noncovalent_bond_kind_ast_strategy(),
        c in noncovalent_bond_kind_ast_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_electron_counts_ast_lattice_laws(
        a in electron_counts_ast_strategy(),
        b in electron_counts_ast_strategy(),
        c in electron_counts_ast_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_stereo_configuration_ast_lattice_laws(
        a in stereo_configuration_lattice_strategy(),
        b in stereo_configuration_lattice_strategy(),
        c in stereo_configuration_lattice_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_tetrahedral_stereo_ast_lattice_laws(
        a in tetrahedral_stereo_lattice_strategy(),
        b in tetrahedral_stereo_lattice_strategy(),
        c in tetrahedral_stereo_lattice_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_cis_trans_stereo_ast_lattice_laws(
        a in cis_trans_stereo_lattice_strategy(),
        b in cis_trans_stereo_lattice_strategy(),
        c in cis_trans_stereo_lattice_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_stereogenicity_ast_lattice_laws(
        ra in stereogenicity_relation_lattice_strategy(),
        rb in stereogenicity_relation_lattice_strategy(),
        rc in stereogenicity_relation_lattice_strategy(),
    ) {
        let a = ra;
        let b = rb;
        let c = rc;
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_dative_bond_ast_lattice_laws(
        a in dative_bond_strategy(),
        b in dative_bond_strategy(),
        c in dative_bond_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_multicenter_bond_ast_lattice_laws(
        a in multicenter_bond_ast_strategy(),
        b in multicenter_bond_ast_strategy(),
        c in multicenter_bond_ast_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_aromatic_system_ast_lattice_laws(
        a in aromatic_system_ast_strategy(),
        b in aromatic_system_ast_strategy(),
        c in aromatic_system_ast_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_noncovalent_bond_ast_lattice_laws(
        a in noncovalent_bond_ast_strategy(),
        b in noncovalent_bond_ast_strategy(),
        c in noncovalent_bond_ast_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_atom_constraints_lattice_laws(
        a in atom_constraints_strategy(),
        b in atom_constraints_strategy(),
        c in atom_constraints_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_atom_constraint_ast_lattice_laws(
        a in atom_constraint_strategy().prop_map(|value| value.canonicalize().unwrap()),
        b in atom_constraint_strategy().prop_map(|value| value.canonicalize().unwrap()),
        c in atom_constraint_strategy().prop_map(|value| value.canonicalize().unwrap()),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_bond_constraints_lattice_laws(
        a in bond_constraints_strategy(),
        b in bond_constraints_strategy(),
        c in bond_constraints_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_bond_constraint_ast_lattice_laws(
        a in bond_constraint_strategy().prop_map(|value| value.canonicalize().unwrap()),
        b in bond_constraint_strategy().prop_map(|value| value.canonicalize().unwrap()),
        c in bond_constraint_strategy().prop_map(|value| value.canonicalize().unwrap()),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_aromatic_system_constraints_lattice_laws(
        a in optional_aromatic_electron_count(),
        b in optional_aromatic_electron_count(),
        c in optional_aromatic_electron_count(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_aromatic_system_constraint_ast_lattice_laws(
        a in constraint_value_strategy(0..=8).prop_map(|value| AromaticSystemConstraintAst::ElectronCount(value).canonicalize().unwrap()),
        b in constraint_value_strategy(0..=8).prop_map(|value| AromaticSystemConstraintAst::ElectronCount(value).canonicalize().unwrap()),
        c in constraint_value_strategy(0..=8).prop_map(|value| AromaticSystemConstraintAst::ElectronCount(value).canonicalize().unwrap()),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_dative_bond_constraints_lattice_laws(
        a in dative_bond_constraints_strategy(),
        b in dative_bond_constraints_strategy(),
        c in dative_bond_constraints_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_dative_bond_constraint_ast_lattice_laws(
        a in dative_bond_constraint_strategy().prop_map(|value| value.canonicalize().unwrap()),
        b in dative_bond_constraint_strategy().prop_map(|value| value.canonicalize().unwrap()),
        c in dative_bond_constraint_strategy().prop_map(|value| value.canonicalize().unwrap()),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_multicenter_bond_constraints_lattice_laws(
        a in optional_multicenter_electron_count(),
        b in optional_multicenter_electron_count(),
        c in optional_multicenter_electron_count(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_multicenter_bond_constraint_ast_lattice_laws(
        a in constraint_value_strategy(0..=8).prop_map(|value| MulticenterBondConstraintAst::ElectronCount(value).canonicalize().unwrap()),
        b in constraint_value_strategy(0..=8).prop_map(|value| MulticenterBondConstraintAst::ElectronCount(value).canonicalize().unwrap()),
        c in constraint_value_strategy(0..=8).prop_map(|value| MulticenterBondConstraintAst::ElectronCount(value).canonicalize().unwrap()),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_noncovalent_bond_constraints_lattice_laws(
        a in noncovalent_bond_constraints_strategy(),
        b in noncovalent_bond_constraints_strategy(),
        c in noncovalent_bond_constraints_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_noncovalent_bond_constraint_ast_lattice_laws(
        a in noncovalent_bond_constraint_strategy(),
        b in noncovalent_bond_constraint_strategy(),
        c in noncovalent_bond_constraint_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_stereo_atom_constraints_lattice_laws(
        a in stereo_atom_constraints_strategy(StereoKind::Tetrahedral),
        b in stereo_atom_constraints_strategy(StereoKind::Tetrahedral),
        c in stereo_atom_constraints_strategy(StereoKind::Tetrahedral),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_stereo_atom_constraint_ast_lattice_laws(
        a in stereo_atom_constraint_strategy(StereoKind::Tetrahedral),
        b in stereo_atom_constraint_strategy(StereoKind::Tetrahedral),
        c in stereo_atom_constraint_strategy(StereoKind::Tetrahedral),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_stereo_bond_constraints_lattice_laws(
        a in stereo_bond_constraints_strategy(StereoKind::CisTrans),
        b in stereo_bond_constraints_strategy(StereoKind::CisTrans),
        c in stereo_bond_constraints_strategy(StereoKind::CisTrans),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_stereo_bond_constraint_ast_lattice_laws(
        a in stereo_bond_constraint_strategy(StereoKind::CisTrans),
        b in stereo_bond_constraint_strategy(StereoKind::CisTrans),
        c in stereo_bond_constraint_strategy(StereoKind::CisTrans),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_topicity_relation_ast_lattice_laws(
        a in topicity_relation_lattice_strategy(),
        b in topicity_relation_lattice_strategy(),
        c in topicity_relation_lattice_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    // Keyed value semilattices: the strategies vary the sub-key (scope /
    // permutation / pair) so a triple spans fibers, exercising the cross-fiber
    // `meet` → `None` / `join` → `Err(NoJoin)` path the containers never reach.
    #[test]
    fn test_ring_membership_ast_lattice_laws(
        a in ring_membership_lattice_strategy(),
        b in ring_membership_lattice_strategy(),
        c in ring_membership_lattice_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_ligand_symmetry_ast_lattice_laws(
        a in ligand_symmetry_strategy(4),
        b in ligand_symmetry_strategy(4),
        c in ligand_symmetry_strategy(4),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_fluxionality_ast_lattice_laws(
        a in fluxionality_strategy(4),
        b in fluxionality_strategy(4),
        c in fluxionality_strategy(4),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_topicity_ast_lattice_laws(
        a in topicity_strategy(4),
        b in topicity_strategy(4),
        c in topicity_strategy(4),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_atom_ast_lattice_laws(
        a in atom_ast_strategy(),
        b in atom_ast_strategy(),
        c in atom_ast_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

    #[test]
    fn test_bond_ast_lattice_laws(
        a in bond_ast_strategy(),
        b in bond_ast_strategy(),
        c in bond_ast_strategy(),
    ) {
        assert_lattice_laws(&a, &b, &c)?;
        assert_canonical_lattice_laws(&a, &b, &c)?;
    }

}
