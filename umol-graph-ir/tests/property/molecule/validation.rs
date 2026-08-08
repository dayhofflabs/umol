//! Constraint validation preserves the three-valued logical laws independently
//! of combinator ordering and storage location.

use proptest::prelude::*;
use proptest::test_runner::FileFailurePersistence;
use umol_chem::element::Element;
use umol_graph_core::{
    ConnectedComponentsAlgorithm, RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm,
};
use umol_graph_ir::ir::{
    AtomAst, AtomConstraintAst, AtomId, Constraint, ConstraintValidateConfig, ConstraintValidator,
    Constraints, MoleculeAst, MoleculeEntries, SubstructureMatchAlgorithm,
};

use super::REGRESSION_FILE;

const CONFIG: ConstraintValidateConfig = ConstraintValidateConfig {
    relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
    connected_components_algorithm: ConnectedComponentsAlgorithm::Bfs,
    substructure_match_algorithm: SubstructureMatchAlgorithm::GraphAndOverlays,
    subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::Vf2Rdkit,
};

fn molecule_with(constraint: Constraint) -> MoleculeAst {
    MoleculeAst::from_entries(MoleculeEntries {
        atoms: vec![AtomAst::from_element(Element::C)],
        constraints: Constraints::from(constraint),
        ..MoleculeEntries::default()
    })
}

proptest! {
    #![proptest_config(ProptestConfig {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(REGRESSION_FILE))),
        ..ProptestConfig::default()
    })]

    #[test]
    fn test_constraint_and_permutation_invariant(values in prop::collection::vec(0i64..=3, 0..=8)) {
        let constraints: Vec<_> = values
            .iter()
            .map(|&value| Constraint::Atom(AtomId(0), AtomConstraintAst::valence(value)))
            .collect();
        let mut reversed = constraints.clone();
        reversed.reverse();
        let validator = ConstraintValidator::new(CONFIG);
        let forward = validator.validate(&molecule_with(Constraint::And(constraints))).unwrap();
        let reverse = validator.validate(&molecule_with(Constraint::And(reversed))).unwrap();

        prop_assert_eq!(forward.is_determined(), reverse.is_determined());
        prop_assert_eq!(forward.is_underdetermined(), reverse.is_underdetermined());
        prop_assert_eq!(forward.is_contradictory(), reverse.is_contradictory());
    }

    #[test]
    fn test_constraint_or_permutation_invariant(values in prop::collection::vec(0i64..=3, 0..=8)) {
        let constraints: Vec<_> = values
            .iter()
            .map(|&value| Constraint::Atom(AtomId(0), AtomConstraintAst::valence(value)))
            .collect();
        let mut reversed = constraints.clone();
        reversed.reverse();
        let validator = ConstraintValidator::new(CONFIG);
        let forward = validator.validate(&molecule_with(Constraint::Or(constraints))).unwrap();
        let reverse = validator.validate(&molecule_with(Constraint::Or(reversed))).unwrap();

        prop_assert_eq!(forward.is_determined(), reverse.is_determined());
        prop_assert_eq!(forward.is_underdetermined(), reverse.is_underdetermined());
        prop_assert_eq!(forward.is_contradictory(), reverse.is_contradictory());
    }

    #[test]
    fn test_constraint_double_negation(value in 0i64..=3) {
        let constraint = Constraint::Atom(AtomId(0), AtomConstraintAst::valence(value));
        let double_negation = Constraint::Not(Box::new(Constraint::Not(Box::new(
            constraint.clone(),
        ))));
        let validator = ConstraintValidator::new(CONFIG);
        let direct = validator.validate(&molecule_with(constraint)).unwrap();
        let negated = validator.validate(&molecule_with(double_negation)).unwrap();

        prop_assert_eq!(direct.is_determined(), negated.is_determined());
        prop_assert_eq!(direct.is_underdetermined(), negated.is_underdetermined());
        prop_assert_eq!(direct.is_contradictory(), negated.is_contradictory());
    }

    #[test]
    fn test_constraint_inline_top_level_leaf_agreement(value in 0i64..=3) {
        let constraint = AtomConstraintAst::valence(value);
        let inline = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomAst::from_element(Element::C).with_constraint(constraint.clone())],
            ..MoleculeEntries::default()
        });
        let top_level = molecule_with(Constraint::Atom(AtomId(0), constraint));
        let validator = ConstraintValidator::new(CONFIG);

        prop_assert_eq!(validator.validate(&inline), validator.validate(&top_level));
    }
}

#[test]
fn test_constraint_vacuous_conjunction() {
    assert!(ConstraintValidator::new(CONFIG)
        .validate(&molecule_with(Constraint::And(Vec::new())))
        .unwrap()
        .is_determined());
}
