use proptest::prelude::*;
use crate::strategies::*;

proptest! {
    #[test]
    fn test_stereo_atom_dsl_display_from_str_roundtrip(
        stereo in stereo_atom_ast_strategy(),
    ) {
        let dsl = StereoAtomDsl(stereo);
        let rendered = dsl.to_string();
        let parsed: StereoAtomDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_stereo_bond_dsl_display_from_str_roundtrip(
        stereo in stereo_bond_ast_strategy(),
    ) {
        let dsl = StereoBondDsl(stereo);
        let rendered = dsl.to_string();
        let parsed: StereoBondDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    /// Canonical `Th0`/`Th1` render to the `:ccw`/`:cw` EDN keyword shorthand
    /// and parse back to the same AST.
    #[test]
    fn test_stereo_atom_dsl_keyword_to_edn_from_edn_roundtrip(
        coset in prop_oneof![Just(0u32), Just(1u32)],
    ) {
        let dsl = StereoAtomDsl(StereoAtomAst::new(
            StereoKind::Tetrahedral,
            StereoCosetAst::Lit(coset),
        ));
        let edn = dsl.to_edn();
        prop_assert!(
            matches!(&edn, Edn::Keyword(_)),
            "expected keyword render for canonical stereo atom, got {edn:?}",
        );
        let parsed = StereoAtomDsl::from_edn(&edn).map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    /// Canonical `Ct0`/`Ct1` render to the `:z`/`:e` EDN keyword shorthand.
    #[test]
    fn test_stereo_bond_dsl_keyword_to_edn_from_edn_roundtrip(
        coset in prop_oneof![Just(0u32), Just(1u32)],
    ) {
        let dsl = StereoBondDsl(StereoBondAst::new(
            StereoKind::CisTrans,
            StereoCosetAst::Lit(coset),
        ));
        let edn = dsl.to_edn();
        prop_assert!(
            matches!(&edn, Edn::Keyword(_)),
            "expected keyword render for canonical stereo bond, got {edn:?}",
        );
        let parsed = StereoBondDsl::from_edn(&edn).map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    /// Molecule-scope (EDN-shaped) stereo atom constraint: `to_edn` → `from_edn`
    /// roundtrips for every kind and constraint variant.
    #[test]
    fn test_stereo_atom_constraint_dsl_to_edn_from_edn_roundtrip(
        args in stereo_atom_kind_strategy().prop_flat_map(|kind| {
            stereo_atom_constraint_strategy(kind).prop_map(move |c| (kind, c))
        }),
    ) {
        let (kind, constraint) = args;
        let dsl = StereoAtomConstraintDsl(kind, constraint);
        let edn = dsl.to_edn();
        let parsed = StereoAtomConstraintDsl::from_edn(&edn).map_err(|e| {
            TestCaseError::fail(format!("from_edn failed for {edn:?}: {e}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_stereo_bond_constraint_dsl_to_edn_from_edn_roundtrip(
        constraint in stereo_bond_constraint_strategy(StereoKind::CisTrans),
    ) {
        let dsl = StereoBondConstraintDsl(StereoKind::CisTrans, constraint);
        let edn = dsl.to_edn();
        let parsed = StereoBondConstraintDsl::from_edn(&edn).map_err(|e| {
            TestCaseError::fail(format!("from_edn failed for {edn:?}: {e}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    /// `StereoLigandPair::new` normalizes to `first <= second` and is symmetric.
    #[test]
    fn test_ligand_pair_ast_normalization(a in 0u8..6, b in 0u8..6) {
        let pair = StereoLigandPair::new(StereoLigandId(a), StereoLigandId(b));
        prop_assert!(pair.first().0 <= pair.second().0);
        prop_assert_eq!(pair, StereoLigandPair::new(StereoLigandId(b), StereoLigandId(a)));
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
}
