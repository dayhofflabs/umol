use rstest::*;

use super::super::ApplicationState;
use super::removal_entries;
use crate::ir::constraint::{AtomConstraintForm, BondConstraintForm, Constraint};
use crate::ir::edit::{CascadedConstraints, Edit, RemovedConstraint, Undo};
use crate::ir::entity::EntityKind;
use crate::ir::id::{AtomId, BondId};
use crate::ir::molecule::MoleculeEntries;
use crate::ir::num::NumForm;
use crate::ir::{Molecule, TransactionError};

#[rstest]
fn test_molecule_apply_edit_add_molecule_constraint(
    mut removal_entries: MoleculeEntries,
    #[values(false, true)] journaled: bool,
) {
    let repeated = Constraint::Atom(AtomId(0), AtomConstraintForm::degree(1));
    let middle = Constraint::Atom(AtomId(1), AtomConstraintForm::degree(2));
    removal_entries.constraints = vec![repeated.clone(), middle, repeated.clone()].into();
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::AddMoleculeConstraint {
        constraint: Box::new(repeated.clone().into()),
    };
    let undo = if journaled {
        let undo = molecule.apply_edit_with_undo(edit, &mut state).unwrap();
        assert_eq!(
            undo,
            Some(Undo::RemoveAddedMoleculeConstraint { position: 3 })
        );
        undo
    } else {
        molecule.apply_edit(edit, &mut state).unwrap();
        None
    };
    removal_entries.constraints.push(repeated);
    assert_eq!(molecule, Molecule::from_entries(removal_entries));
    if let Some(undo) = undo {
        molecule.apply_undo(undo);
        assert_eq!(molecule, initial);
    }
}

#[rstest]
#[case::repeated(Constraint::Atom(AtomId(0), AtomConstraintForm::degree(1)))]
#[case::equivalent(Constraint::Atom(
    AtomId(0),
    AtomConstraintForm::degree(NumForm::lit_set([1]))
))]
fn test_molecule_apply_edit_remove_molecule_constraint(
    mut removal_entries: MoleculeEntries,
    #[case] copy: Constraint,
    #[values(false, true)] journaled: bool,
) {
    let removed = Constraint::Atom(AtomId(0), AtomConstraintForm::degree(1));
    let middle = Constraint::Atom(AtomId(1), AtomConstraintForm::degree(2));
    removal_entries.constraints = vec![removed.clone(), middle.clone(), copy.clone()].into();
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::RemoveMoleculeConstraint {
        constraint: Box::new(removed.clone().into()),
    };
    let undo = if journaled {
        let undo = molecule.apply_edit_with_undo(edit, &mut state).unwrap();
        assert_eq!(
            undo,
            Some(Undo::RestoreMoleculeConstraints(CascadedConstraints {
                removed: vec![
                    RemovedConstraint {
                        position: 0,
                        constraint: removed,
                    },
                    RemovedConstraint {
                        position: 2,
                        constraint: copy,
                    },
                ],
                modified: Vec::new(),
            }))
        );
        undo
    } else {
        molecule.apply_edit(edit, &mut state).unwrap();
        None
    };
    removal_entries.constraints = vec![middle].into();
    assert_eq!(molecule, Molecule::from_entries(removal_entries));
    if let Some(undo) = undo {
        molecule.apply_undo(undo);
        assert_eq!(molecule, initial);
    }
}

#[rstest]
#[case::missing(Edit::RemoveMoleculeConstraint {
    constraint: Box::new(Constraint::Atom(AtomId(0), AtomConstraintForm::degree(2)).into()),
}, TransactionError::MissingEntry)]
#[case::add_nested(Edit::AddMoleculeConstraint {
    constraint: Box::new(
        Constraint::And(vec![
            Constraint::Atom(AtomId(0), AtomConstraintForm::degree(1)),
            Constraint::Or(vec![Constraint::Bond(
                BondId(2),
                BondConstraintForm::aromatic(true),
            )]),
        ])
        .into(),
    ),
}, TransactionError::HandleOutOfRange { kind: EntityKind::Bond, index: 2, count: 2 })]
#[case::remove_nested(Edit::RemoveMoleculeConstraint {
    constraint: Box::new(
        Constraint::And(vec![
            Constraint::Atom(AtomId(0), AtomConstraintForm::degree(1)),
            Constraint::Or(vec![Constraint::Bond(
                BondId(2),
                BondConstraintForm::aromatic(true),
            )]),
        ])
        .into(),
    ),
}, TransactionError::HandleOutOfRange { kind: EntityKind::Bond, index: 2, count: 2 })]
fn test_molecule_apply_edit_molecule_constraint_error(
    mut removal_entries: MoleculeEntries,
    #[case] edit: Edit,
    #[case] expected: TransactionError,
    #[values(false, true)] journaled: bool,
) {
    removal_entries.constraints =
        vec![Constraint::Atom(AtomId(0), AtomConstraintForm::degree(1))].into();
    let initial = Molecule::from_entries(removal_entries);
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let result = if journaled {
        molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
    } else {
        molecule.apply_edit(edit, &mut state)
    };
    assert_eq!(result, Err(expected));
    assert_eq!(molecule, initial);
}
