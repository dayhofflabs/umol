use proptest::prelude::*;
use rstest::rstest;

use crate::strategies::*;

proptest! {
    // The entity DSL types carry a compact string form (Display) parsed by their
    // own `FromStr`; the invariant is `parse(display(x)) == x` for any generated AST.
    #[test]
    fn test_atom_dsl_display_from_str_roundtrip(atom in atom_ast_strategy()) {
        let dsl = AtomDsl(atom);
        let rendered = dsl.to_string();
        let parsed: AtomDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_bond_dsl_display_from_str_roundtrip(bond in bond_ast_strategy()) {
        let dsl = BondDsl(bond);
        let rendered = dsl.to_string();
        let parsed: BondDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    /// `BondDsl::ToEdn` ↔ `FromEdn` round-trips for any generated bond
    /// shape. Non-canonical bonds render as bond strings; canonical
    /// shapes (order-only, no charge / spin / non-aromatic constraints,
    /// or order-1 with the `Aromatic` flag) render as keyword shorthands.
    #[test]
    fn test_bond_dsl_to_edn_from_edn_roundtrip(bond in bond_ast_strategy()) {
        let dsl = BondDsl(bond);
        let edn = dsl.to_edn();
        let parsed = BondDsl::from_edn(&edn).map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    /// Canonical-shape bonds render as keyword shorthands, and the
    /// keyword form parses back to the same AST.
    #[test]
    fn test_bond_dsl_keyword_to_edn_from_edn_roundtrip(
        bond in canonical_keyword_bond_strategy(),
    ) {
        let dsl = BondDsl(bond);
        let edn = dsl.to_edn();
        prop_assert!(
            matches!(&edn, Edn::Keyword(_)),
            "expected keyword render for canonical bond, got {edn:?}",
        );
        let parsed = BondDsl::from_edn(&edn).map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_aromatic_system_dsl_display_from_str_roundtrip(
        system in aromatic_system_ast_strategy(),
    ) {
        let dsl = AromaticSystemDsl(system);
        let rendered = dsl.to_string();
        let parsed: AromaticSystemDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_multicenter_bond_dsl_display_from_str_roundtrip(
        bond in multicenter_bond_ast_strategy(),
    ) {
        let dsl = MulticenterBondDsl(bond);
        let rendered = dsl.to_string();
        let parsed: MulticenterBondDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_dative_bond_dsl_display_from_str_roundtrip(
        bond in dative_bond_strategy(),
    ) {
        let dsl = DativeBondDsl(bond);
        let rendered = dsl.to_string();
        let parsed: DativeBondDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_noncovalent_bond_dsl_display_from_str_roundtrip(
        bond in noncovalent_bond_ast_strategy(),
    ) {
        let dsl = NoncovalentBondDsl(bond);
        let rendered = dsl.to_string();
        let parsed: NoncovalentBondDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    /// Atom-update DSL round-trips through its compact string form, including
    /// omitted and explicitly undetermined fields.
    #[test]
    fn test_atom_update_dsl_display_from_str_roundtrip(update in atom_update_strategy()) {
        let dsl = AtomUpdateDsl(update);
        let rendered = dsl.to_string();
        let parsed: AtomUpdateDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_bond_update_dsl_display_from_str_roundtrip(update in bond_update_strategy()) {
        let dsl = BondUpdateDsl(update);
        let rendered = dsl.to_string();
        let parsed: BondUpdateDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    /// Atom-update DSL round-trips through its EDN string leaf.
    #[test]
    fn test_atom_update_dsl_to_edn_from_edn_roundtrip(update in atom_update_strategy()) {
        let dsl = AtomUpdateDsl(update);
        let edn = dsl.to_edn();
        let parsed = AtomUpdateDsl::from_edn(&edn).map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_bond_update_dsl_to_edn_from_edn_roundtrip(update in bond_update_strategy()) {
        let dsl = BondUpdateDsl(update);
        let edn = dsl.to_edn();
        let parsed = BondUpdateDsl::from_edn(&edn).map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

}

/// Vacuous-payload `AtomConstraintAst` variants render to nothing in the
/// canonical entity-string form. The proptest generator excludes these from
/// roundtrip strategies; this asserts the elision invariant directly so a
/// regression in `fmt_value_field_required` / `fmt_ring_count` / the
/// AromaticValence / MulticenterValence formatters can't slip through.
#[rstest]
#[case::valence(AtomConstraintAst::Valence(ValueAst::Undetermined))]
#[case::total_valence(AtomConstraintAst::TotalValence(ValueAst::Undetermined))]
#[case::donated_pairs(AtomConstraintAst::DonatedPairs(ValueAst::Undetermined))]
#[case::accepted_pairs(AtomConstraintAst::AcceptedPairs(ValueAst::Undetermined))]
#[case::degree(AtomConstraintAst::Degree(ValueAst::Undetermined))]
#[case::total_degree(AtomConstraintAst::TotalDegree(ValueAst::Undetermined))]
#[case::ring_degree(AtomConstraintAst::RingDegree(ValueAst::Undetermined))]
#[case::ring_valence(AtomConstraintAst::RingValence(ValueAst::Undetermined))]
#[case::total_hydrogens(AtomConstraintAst::TotalHydrogens(ValueAst::Undetermined))]
#[case::ring_membership_all(AtomConstraintAst::ring_membership(
    RingScope::All,
    ValueAst::Undetermined
))]
#[case::ring_membership_size(AtomConstraintAst::ring_membership(
    RingScope::All,
    ValueAst::Undetermined
))]
#[case::aromatic_valence_undetermined(AtomConstraintAst::AromaticValence(
    AromaticValenceAst::Undetermined
))]
#[case::multicenter_valence_undetermined(AtomConstraintAst::MulticenterValence(
    MulticenterValenceAst::Undetermined
))]
fn test_atom_dsl_vacuous_constraint_renders_empty(#[case] vacuous: AtomConstraintAst) {
    let mut atom = AtomAst::default();
    atom.constraints.set(vacuous);
    let with_vacuous = AtomDsl(atom).to_string();
    let bare = AtomDsl(AtomAst::default()).to_string();
    assert_eq!(with_vacuous, bare);
}
