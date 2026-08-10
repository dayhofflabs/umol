use proptest::prelude::*;
use proptest::test_runner::TestCaseResult;
use rstest::rstest;

use crate::strategies::*;

fn assert_exact_values<T>(
    mut iterator: impl ExactSizeIterator<Item = T>,
    expected: &[T],
    prefix: usize,
) -> TestCaseResult
where
    T: Debug + PartialEq,
{
    let prefix = prefix.min(expected.len());
    prop_assert_eq!(iterator.len(), expected.len());
    prop_assert_eq!(iterator.size_hint(), (expected.len(), Some(expected.len())),);
    for expected_item in &expected[..prefix] {
        let previous = iterator.len();
        let actual = iterator.next();
        prop_assert_eq!(actual.as_ref(), Some(expected_item));
        let remaining = iterator.len();
        prop_assert_eq!(remaining, previous - 1);
        prop_assert_eq!(iterator.size_hint(), (remaining, Some(remaining)));
    }
    for expected_item in &expected[prefix..] {
        let previous = iterator.len();
        let actual = iterator.next();
        prop_assert_eq!(actual.as_ref(), Some(expected_item));
        let remaining = iterator.len();
        prop_assert_eq!(remaining, previous - 1);
        prop_assert_eq!(iterator.size_hint(), (remaining, Some(remaining)));
    }
    prop_assert_eq!(iterator.next(), None);
    prop_assert_eq!(iterator.len(), 0);
    Ok(())
}

proptest! {
    // The entity DSL types carry a compact string form (Display) parsed by their
    // own `FromStr`; the invariant is `parse(display(x)) == x` for any generated form.
    #[test]
    fn test_atom_dsl_display_from_str_roundtrip(atom in atom_form_strategy()) {
        let dsl = AtomDsl(atom);
        let rendered = dsl.to_string();
        let parsed: AtomDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_bond_dsl_display_from_str_roundtrip(bond in bond_form_strategy()) {
        let dsl = BondDsl(bond);
        let rendered = dsl.to_string();
        let parsed: BondDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    /// `BondDsl::ToEdn` ↔ `FromEdn` round-trips for any generated bond
    /// shape. Non-canonical bonds render as bond strings; canonical
    /// shapes (order-only, no charge / unpaired-electron fields / non-aromatic constraints,
    /// or order-1 with the `Aromatic` flag) render as keyword shorthands.
    #[test]
    fn test_bond_dsl_to_edn_from_edn_roundtrip(bond in bond_form_strategy()) {
        let dsl = BondDsl(bond);
        let edn = dsl.to_edn();
        let parsed = BondDsl::from_edn(&edn).map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    /// Canonical-shape bonds render as keyword shorthands, and the
    /// keyword form parses back to the same form.
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
        system in aromatic_system_form_strategy(),
    ) {
        let dsl = AromaticSystemDsl(system);
        let rendered = dsl.to_string();
        let parsed: AromaticSystemDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_aromatic_system_update_dsl_display_from_str_roundtrip(
        update in aromatic_system_update_strategy(),
    ) {
        let dsl = AromaticSystemUpdateDsl(update);
        let rendered = dsl.to_string();
        let parsed: AromaticSystemUpdateDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_aromatic_system_update_display_from_str_roundtrip(
        update in aromatic_system_update_strategy(),
    ) {
        let rendered = update.to_string();
        let parsed: AromaticSystemUpdate = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(update, parsed);
    }

    #[test]
    fn test_aromatic_system_update_dsl_to_edn_from_edn_roundtrip(
        update in aromatic_system_update_strategy(),
    ) {
        let dsl = AromaticSystemUpdateDsl(update);
        let edn = dsl.to_edn();
        let parsed = AromaticSystemUpdateDsl::from_edn(&edn)
            .map_err(|e| TestCaseError::fail(format!("parse failed: {e}")))?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_multicenter_bond_dsl_display_from_str_roundtrip(
        bond in multicenter_bond_form_strategy(),
    ) {
        let dsl = MulticenterBondDsl(bond);
        let rendered = dsl.to_string();
        let parsed: MulticenterBondDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_multicenter_bond_update_dsl_display_from_str_roundtrip(
        update in multicenter_bond_update_strategy(),
    ) {
        let dsl = MulticenterBondUpdateDsl(update);
        let rendered = dsl.to_string();
        let parsed: MulticenterBondUpdateDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_multicenter_bond_update_display_from_str_roundtrip(
        update in multicenter_bond_update_strategy(),
    ) {
        let rendered = update.to_string();
        let parsed: MulticenterBondUpdate = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(update, parsed);
    }

    #[test]
    fn test_multicenter_bond_update_dsl_to_edn_from_edn_roundtrip(
        update in multicenter_bond_update_strategy(),
    ) {
        let dsl = MulticenterBondUpdateDsl(update);
        let edn = dsl.to_edn();
        let parsed = MulticenterBondUpdateDsl::from_edn(&edn)
            .map_err(|e| TestCaseError::fail(format!("parse failed: {e}")))?;
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
    fn test_dative_bond_update_dsl_display_from_str_roundtrip(
        update in dative_bond_update_strategy(),
    ) {
        let dsl = DativeBondUpdateDsl(update);
        let rendered = dsl.to_string();
        let parsed: DativeBondUpdateDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_dative_bond_update_display_from_str_roundtrip(
        update in dative_bond_update_strategy(),
    ) {
        let rendered = update.to_string();
        let parsed: DativeBondUpdate = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(update, parsed);
    }

    #[test]
    fn test_dative_bond_update_dsl_to_edn_from_edn_roundtrip(
        update in dative_bond_update_strategy(),
    ) {
        let dsl = DativeBondUpdateDsl(update);
        let edn = dsl.to_edn();
        let parsed = DativeBondUpdateDsl::from_edn(&edn)
            .map_err(|e| TestCaseError::fail(format!("parse failed: {e}")))?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_noncovalent_bond_dsl_display_from_str_roundtrip(
        bond in noncovalent_bond_form_strategy(),
    ) {
        let dsl = NoncovalentBondDsl(bond);
        let rendered = dsl.to_string();
        let parsed: NoncovalentBondDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_noncovalent_bond_update_dsl_display_from_str_roundtrip(
        update in noncovalent_bond_update_strategy(),
    ) {
        let dsl = NoncovalentBondUpdateDsl(update);
        let rendered = dsl.to_string();
        let parsed: NoncovalentBondUpdateDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_noncovalent_bond_update_display_from_str_roundtrip(
        update in noncovalent_bond_update_strategy(),
    ) {
        let rendered = update.to_string();
        let parsed: NoncovalentBondUpdate = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(update, parsed);
    }

    #[test]
    fn test_noncovalent_bond_update_dsl_to_edn_from_edn_roundtrip(
        update in noncovalent_bond_update_strategy(),
    ) {
        let dsl = NoncovalentBondUpdateDsl(update);
        let edn = dsl.to_edn();
        let parsed = NoncovalentBondUpdateDsl::from_edn(&edn)
            .map_err(|e| TestCaseError::fail(format!("parse failed: {e}")))?;
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
    fn test_atom_update_display_from_str_roundtrip(update in atom_update_strategy()) {
        let rendered = update.to_string();
        let parsed: AtomUpdate = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(update, parsed);
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

    #[test]
    fn test_bond_update_display_from_str_roundtrip(update in bond_update_strategy()) {
        let rendered = update.to_string();
        let parsed: BondUpdate = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(update, parsed);
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

proptest! {
    #[test]
    fn test_atom_constraints_form_iterators_exact_size(
        mut constraints in atom_constraints_strategy(),
        prefix in any::<usize>(),
    ) {
        let expected = constraints.iter().cloned().collect::<Vec<_>>();
        assert_exact_values(constraints.iter().cloned(), &expected, prefix)?;
        assert_exact_values(constraints.take(), &expected, prefix)?;
        prop_assert_eq!(constraints, AtomConstraintsForm::new());
    }

    #[test]
    fn test_bond_constraints_form_take_exact_size(
        mut constraints in bond_constraints_strategy(),
        prefix in any::<usize>(),
    ) {
        let expected = constraints.iter().cloned().collect::<Vec<_>>();
        assert_exact_values(constraints.take(), &expected, prefix)?;
        prop_assert_eq!(constraints, BondConstraintsForm::new());
    }

    #[test]
    fn test_dative_bond_constraints_form_take_exact_size(
        mut constraints in dative_bond_constraints_strategy(),
        prefix in any::<usize>(),
    ) {
        let expected = constraints.iter().cloned().collect::<Vec<_>>();
        assert_exact_values(constraints.take(), &expected, prefix)?;
        prop_assert_eq!(constraints, DativeBondConstraintsForm::new());
    }

    #[test]
    fn test_aromatic_system_constraints_form_take_exact_size(
        mut constraints in aromatic_system_form_strategy().prop_map(|form| form.constraints),
        prefix in any::<usize>(),
    ) {
        let expected = constraints.iter().cloned().collect::<Vec<_>>();
        assert_exact_values(constraints.take(), &expected, prefix)?;
        prop_assert_eq!(constraints, AromaticSystemConstraintsForm::new());
    }

    #[test]
    fn test_multicenter_bond_constraints_form_take_exact_size(
        mut constraints in multicenter_bond_form_strategy().prop_map(|form| form.constraints),
        prefix in any::<usize>(),
    ) {
        let expected = constraints.iter().cloned().collect::<Vec<_>>();
        assert_exact_values(constraints.take(), &expected, prefix)?;
        prop_assert_eq!(constraints, MulticenterBondConstraintsForm::new());
    }

    #[test]
    fn test_noncovalent_bond_constraints_form_take_exact_size(
        mut constraints in noncovalent_bond_constraints_strategy(),
        prefix in any::<usize>(),
    ) {
        let expected = constraints.iter().cloned().collect::<Vec<_>>();
        assert_exact_values(constraints.take(), &expected, prefix)?;
        prop_assert_eq!(constraints, NoncovalentBondConstraintsForm::new());
    }

    #[test]
    fn test_stereo_atom_constraints_form_take_exact_size(
        mut constraints in stereo_atom_constraints_strategy(StereoKind::Tetrahedral),
        prefix in any::<usize>(),
    ) {
        let expected = constraints.iter().cloned().collect::<Vec<_>>();
        assert_exact_values(constraints.take(), &expected, prefix)?;
        prop_assert_eq!(constraints, StereoAtomConstraintsForm::new());
    }

    #[test]
    fn test_stereo_bond_constraints_form_take_exact_size(
        mut constraints in stereo_bond_constraints_strategy(StereoKind::CisTrans),
        prefix in any::<usize>(),
    ) {
        let expected = constraints.iter().cloned().collect::<Vec<_>>();
        assert_exact_values(constraints.take(), &expected, prefix)?;
        prop_assert_eq!(constraints, StereoBondConstraintsForm::new());
    }
}

/// Vacuous-payload `AtomConstraintForm` variants render to nothing in the
/// canonical entity-string form. The proptest generator excludes these from
/// roundtrip strategies; this asserts the elision invariant directly so a
/// regression in `fmt_num_field_required` / `fmt_ring_count` / the
/// AromaticValence / MulticenterValence formatters can't slip through.
#[rstest]
#[case::valence(AtomConstraintForm::Valence(NumForm::Undetermined))]
#[case::total_valence(AtomConstraintForm::TotalValence(NumForm::Undetermined))]
#[case::donated_pairs(AtomConstraintForm::DonatedPairs(NumForm::Undetermined))]
#[case::accepted_pairs(AtomConstraintForm::AcceptedPairs(NumForm::Undetermined))]
#[case::degree(AtomConstraintForm::Degree(NumForm::Undetermined))]
#[case::total_degree(AtomConstraintForm::TotalDegree(NumForm::Undetermined))]
#[case::ring_degree(AtomConstraintForm::RingDegree(NumForm::Undetermined))]
#[case::ring_valence(AtomConstraintForm::RingValence(NumForm::Undetermined))]
#[case::total_hydrogens(AtomConstraintForm::TotalHydrogens(NumForm::Undetermined))]
#[case::ring_membership_all(AtomConstraintForm::ring_membership(
    RingScope::All,
    NumForm::Undetermined
))]
#[case::ring_membership_size(AtomConstraintForm::ring_membership(
    RingScope::All,
    NumForm::Undetermined
))]
#[case::aromatic_valence_undetermined(AtomConstraintForm::AromaticValence(
    AromaticValenceForm::Undetermined
))]
#[case::multicenter_valence_undetermined(AtomConstraintForm::MulticenterValence(
    MulticenterValenceForm::Undetermined
))]
fn test_atom_dsl_vacuous_constraint_renders_empty(#[case] vacuous: AtomConstraintForm) {
    let mut atom = AtomForm::default();
    atom.constraints.set(vacuous);
    let with_vacuous = AtomDsl(atom).to_string();
    let bare = AtomDsl(AtomForm::default()).to_string();
    assert_eq!(with_vacuous, bare);
}
