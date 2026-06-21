use proptest::prelude::*;
use umol_ast::ast::StereoKind;

use crate::strategies::*;

// Lattice-law sweep: every `impl Lattice` type satisfies commutativity,
// associativity, absorption, idempotence, and `matches`↔`meet` consistency,
// checked by `assert_lattice_laws` over a generated value triple.
// `NoncovalentBondConstraints` is omitted: its inner enum is uninhabited, so the
// collection has the single empty value and the laws would be vacuous.

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
    fn test_atom_constraints_lattice_laws(
        a in atom_constraints_strategy(),
        b in atom_constraints_strategy(),
        c in atom_constraints_strategy(),
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
    fn test_aromatic_system_constraints_lattice_laws(
        a in optional_aromatic_electron_count(),
        b in optional_aromatic_electron_count(),
        c in optional_aromatic_electron_count(),
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
    fn test_multicenter_bond_constraints_lattice_laws(
        a in optional_multicenter_electron_count(),
        b in optional_multicenter_electron_count(),
        c in optional_multicenter_electron_count(),
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
    fn test_stereo_bond_constraints_lattice_laws(
        a in stereo_bond_constraints_strategy(StereoKind::CisTrans),
        b in stereo_bond_constraints_strategy(StereoKind::CisTrans),
        c in stereo_bond_constraints_strategy(StereoKind::CisTrans),
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
