use rstest::*;
use umol_graph_core::Compaction;

use super::super::{ApplicationState, HandleTable};
use crate::ir::atom::AtomForm;
use crate::ir::constraint::{AtomConstraintForm, AtomConstraintKey};
use crate::ir::edit::{AtomHandle, Edit, Edits};
use crate::ir::entity::EntityKind;
use crate::ir::id::AtomId;
use crate::ir::molecule::MoleculeEntries;
use crate::ir::num::NumForm;
use crate::ir::{Molecule, RollbackMethod, TransactionError};

#[rstest]
#[case::first(3, 0, Ok(AtomId(0)))]
#[case::last(3, 2, Ok(AtomId(2)))]
#[case::outside(3, 3, Err(TransactionError::HandleOutOfRange {
    kind: EntityKind::Atom, index: 3, count: 3,
}))]
#[case::empty(0, 0, Err(TransactionError::HandleOutOfRange {
    kind: EntityKind::Atom, index: 0, count: 0,
}))]
fn test_handle_table_initial(
    #[case] count: usize,
    #[case] index: usize,
    #[case] expected: Result<AtomId, TransactionError>,
) {
    let table = HandleTable::new(count);

    assert_eq!(table.initial(EntityKind::Atom, index), expected);
    assert_eq!(table.initial, None);
}

#[rstest]
#[case::first(0, Ok(AtomId(3)))]
#[case::second(1, Ok(AtomId(4)))]
#[case::outside(2, Err(TransactionError::HandleOutOfRange {
    kind: EntityKind::Atom, index: 2, count: 2,
}))]
fn test_handle_table_created(
    #[case] index: usize,
    #[case] expected: Result<AtomId, TransactionError>,
) {
    let mut table = HandleTable::new(3);
    table.push(AtomId(3));
    table.push(AtomId(4));

    assert_eq!(table.created(EntityKind::Atom, index), expected);
    assert_eq!(table.initial, None);
}

#[rstest]
#[case::identity(vec![], vec![Some(AtomId(0)), Some(AtomId(1)), Some(AtomId(2))], vec![Some(AtomId(3)), Some(AtomId(4))])]
#[case::initial_and_created(vec![AtomId(1), AtomId(3)], vec![Some(AtomId(0)), None, Some(AtomId(1))], vec![None, Some(AtomId(2))])]
#[case::all(vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4)], vec![None, None, None], vec![None, None])]
fn test_handle_table_compact(
    #[case] removed: Vec<AtomId>,
    #[case] initial: Vec<Option<AtomId>>,
    #[case] created: Vec<Option<AtomId>>,
) {
    let mut table = HandleTable::new(3);
    table.push(AtomId(3));
    table.push(AtomId(4));
    let compaction = Compaction::new(5, removed).unwrap();

    table.compact(|id| compaction.compact(id));

    assert_eq!(table.initial_count, 3);
    assert_eq!(table.initial, Some(initial));
    assert_eq!(table.created, created);
}

#[rstest]
fn test_handle_table_compact_sequence() {
    let mut table = HandleTable::new(3);
    table.push(AtomId(3));
    table.push(AtomId(4));
    let first = Compaction::new(5, vec![AtomId(1), AtomId(3)]).unwrap();
    table.compact(|id| first.compact(id));
    table.push(AtomId(3));
    let second = Compaction::new(4, vec![AtomId(0), AtomId(2)]).unwrap();
    table.compact(|id| second.compact(id));

    assert_eq!(table.initial_count, 3);
    assert_eq!(table.initial, Some(vec![None, None, Some(AtomId(0))]));
    assert_eq!(table.created, vec![None, None, Some(AtomId(1))]);
    assert_eq!(table.initial(EntityKind::Atom, 2), Ok(AtomId(0)));
    assert_eq!(table.created(EntityKind::Atom, 2), Ok(AtomId(1)));
    assert_eq!(
        table.initial(EntityKind::Atom, 1),
        Err(TransactionError::HandleRemoved {
            kind: EntityKind::Atom,
            index: 1
        }),
    );
    assert_eq!(
        table.created(EntityKind::Atom, 1),
        Err(TransactionError::HandleRemoved {
            kind: EntityKind::Atom,
            index: 1
        }),
    );
    assert_eq!(
        table.initial(EntityKind::Atom, 3),
        Err(TransactionError::HandleOutOfRange {
            kind: EntityKind::Atom,
            index: 3,
            count: 3
        }),
    );
    assert_eq!(
        table.created(EntityKind::Atom, 3),
        Err(TransactionError::HandleOutOfRange {
            kind: EntityKind::Atom,
            index: 3,
            count: 3
        }),
    );
}

#[rstest]
fn test_molecule_apply_atom_batches_removed(
    #[values(false, true)] created: bool,
    #[values(false, true)] constraint: bool,
    #[values(false, true)] journaled: bool,
) {
    let source = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::default(); 2],
        ..Default::default()
    });
    let mut prefix = Edits::new();
    let removed = if created {
        prefix.add_atom(AtomForm::default())
    } else {
        AtomHandle::Id(AtomId(1))
    };
    prefix.remove_atom(removed.clone());
    let batch = if constraint {
        Edit::ModifyAtomConstraints {
            key: AtomConstraintKey::Valence,
            changes: vec![
                (
                    AtomHandle::Id(AtomId(0)),
                    None,
                    Some(AtomConstraintForm::valence(1)),
                ),
                (removed, None, None),
            ],
        }
    } else {
        Edit::ModifyAtomCharges {
            changes: vec![
                (
                    AtomHandle::Id(AtomId(0)),
                    NumForm::Undetermined,
                    NumForm::Lit(1),
                ),
                (removed, NumForm::Undetermined, NumForm::Lit(1)),
            ],
        }
    };
    let error = TransactionError::HandleRemoved {
        kind: EntityKind::Atom,
        index: if created { 0 } else { 1 },
    };
    let mut molecule = source.clone();
    let mut state = ApplicationState::new(&molecule);
    for edit in prefix.clone() {
        molecule.apply_edit(edit, &mut state).unwrap();
    }
    let result = if journaled {
        molecule
            .apply_edit_with_undo(batch.clone(), &mut state)
            .map(|_| ())
    } else {
        molecule.apply_edit(batch.clone(), &mut state)
    };
    assert_eq!(result, Err(error.clone()));
    let mut atoms = vec![AtomForm::default(); if created { 2 } else { 1 }];
    if !journaled {
        if constraint {
            atoms[0].constraints = vec![AtomConstraintForm::valence(1)].into();
        } else {
            atoms[0].charge = NumForm::Lit(1);
        }
    }
    assert_eq!(
        molecule,
        Molecule::from_entries(MoleculeEntries {
            atoms,
            ..Default::default()
        })
    );
    prefix.push(batch);
    let mut molecule = source.clone();
    assert_eq!(
        molecule.transact([prefix], RollbackMethod::Journal),
        Err(error.into())
    );
    assert_eq!(molecule, source);
}
