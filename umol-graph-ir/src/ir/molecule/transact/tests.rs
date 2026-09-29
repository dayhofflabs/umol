//! Transaction batch execution, undo entries, and rollback.
use std::collections::HashMap;

use rstest::*;
use umol_chem::element::Element;

use super::{Transaction, TransactionError};
use crate::ir::aromatic::AromaticSystemForm;
use crate::ir::atom::{AtomForm, ElementForm};
use crate::ir::bond::BondForm;
use crate::ir::constraint::{
    AromaticSystemConstraintForm, AtomConstraintForm, BondConstraintForm, Constraint,
    DativeBondConstraintForm, MoleculeConstraint, MulticenterBondConstraintForm,
    NoncovalentBondConstraintForm, RelationalConstraint, RingScope, StereoAtomConstraintForm,
    StereoBondConstraintForm, StereogenicityForm,
};
use crate::ir::dative::DativeBondForm;
use crate::ir::edit::{
    AddBond, AromaticSystemFieldChange, AromaticSystemHandle, AtomFieldChange, AtomHandle,
    BondFieldChange, BondHandle, ConstraintEdit, DativeBondFieldChange, DativeBondHandle, Edit,
    Edits, EntityHandle, MulticenterBondFieldChange, MulticenterBondHandle,
    NoncovalentBondFieldChange, NoncovalentBondHandle, StereoAtomFieldChange, StereoAtomHandle,
    StereoBondFieldChange, StereoBondHandle, Undo,
};
use crate::ir::entity::{Entity, EntityKind};
use crate::ir::error::MoleculeApplyError;
use crate::ir::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ir::ligand::{StereoLigand, StereoLigandKind};
use crate::ir::molecule::editor::MoleculeEditor;
use crate::ir::molecule::Molecule;
use crate::ir::multicenter::MulticenterBondForm;
use crate::ir::noncovalent::{NoncovalentBondForm, NoncovalentBondKind, NoncovalentBondKindForm};
use crate::ir::num::NumForm;
use crate::ir::stereo::{
    CisTransStereoForm, StereoAtomForm, StereoBondForm, StereoConfigurationForm, StereoCoset,
    StereoKind,
};
use crate::ir::traits::Normalize;
use crate::ir::BooleanForm;

#[fixture]
fn empty() -> MoleculeEditor {
    Molecule::default().edit()
}

#[fixture]
fn one_atom() -> MoleculeEditor {
    let mut b = Molecule::default().edit();
    b.add_atom(AtomForm::from_element(Element::C));
    b
}

#[fixture]
fn diatomic() -> MoleculeEditor {
    let mut b = Molecule::default().edit();
    b.add_atom(AtomForm::from_element(Element::C));
    b.add_atom(AtomForm::from_element(Element::C));
    b.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));
    b
}

#[rstest]
fn test_transaction_apply_add_atom(empty: MoleculeEditor) {
    let mut edits = Edits::new();
    edits.add_atom(AtomForm::from_element(Element::C));
    let mut molecule = empty.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction.apply(edits).unwrap();
            assert!(matches!(
                transaction.journal.as_slice(),
                [Undo::RemoveAddedTopology { atoms, bonds }]
                    if atoms.iter().map(|a| a.id).collect::<Vec<_>>() == vec![AtomId(0)]
                        && bonds.is_empty()
            ));
            let built = transaction.probe().unwrap();
            assert_eq!(built.atoms().count(), 1);
            assert_eq!(
                built.atom(AtomId(0)).attributes().element,
                ElementForm::Lit(Element::C)
            );

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_add_atoms(empty: MoleculeEditor) {
    let mut edits = Edits::new();
    edits.add_atoms([
        AtomForm::from_element(Element::C),
        AtomForm::from_element(Element::N),
    ]);
    let mut molecule = empty.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction.apply(edits).unwrap();
            assert!(matches!(
                transaction.journal.as_slice(),
                [Undo::RemoveAddedTopology { atoms, bonds }]
                    if atoms.iter().map(|a| a.id).collect::<Vec<_>>() == vec![AtomId(0), AtomId(1)]
                        && bonds.is_empty()
            ));
            let built = transaction.probe().unwrap();
            assert_eq!(built.atoms().count(), 2);
            assert_eq!(
                built.atom(AtomId(0)).attributes().element,
                ElementForm::Lit(Element::C)
            );
            assert_eq!(
                built.atom(AtomId(1)).attributes().element,
                ElementForm::Lit(Element::N)
            );

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_add_bond_via_handle(empty: MoleculeEditor) {
    let mut edits = Edits::new();
    let atoms = edits.add_atoms([
        AtomForm::from_element(Element::C),
        AtomForm::from_element(Element::C),
    ]);
    edits.add_bond(atoms[0].clone(), atoms[1].clone(), BondForm::from_order(1));
    let mut molecule = empty.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction.apply(edits).unwrap();
            assert!(matches!(
                transaction.journal.as_slice(),
                [
                    Undo::RemoveAddedTopology { atoms, bonds },
                    Undo::RemoveAddedTopology { atoms: bond_atoms, bonds: added_bonds },
                ] if atoms.iter().map(|a| a.id).collect::<Vec<_>>() == vec![AtomId(0), AtomId(1)]
                    && bonds.is_empty()
                    && bond_atoms.is_empty()
                    && added_bonds.iter().map(|b| b.id).collect::<Vec<_>>() == vec![BondId(0)]
            ));

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_rollback(one_atom: MoleculeEditor) {
    let before = one_atom.clone().build();
    // Mid-batch failure (out-of-range id on edit 2) rolls back the
    // already-applied AddAtom on edit 1.
    let mut edits = Edits::new();
    edits.add_atom(AtomForm::from_element(Element::N));
    edits.remove_atom(AtomHandle::Id(AtomId(99)));
    let mut molecule = one_atom.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            let err = transaction.apply(edits).unwrap_err();
            assert_eq!(
                err,
                MoleculeApplyError::Transaction(TransactionError::HandleOutOfRange {
                    kind: EntityKind::Atom,
                    index: 99,
                    count: 1,
                })
            );
            assert_eq!(transaction.probe().unwrap(), &before);

            Ok(())
        },
    );
    assert_eq!(result, Err(TransactionError::Aborted.into()));
}

#[rstest]
fn test_transaction_apply_set_atom_field(one_atom: MoleculeEditor) {
    let mut molecule = one_atom.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::ModifyAtomField {
                    id: AtomHandle::Id(AtomId(0)),
                    change: AtomFieldChange::Charge {
                        old: NumForm::default(),
                        new: NumForm::Lit(1),
                    },
                }]))
                .unwrap();
            assert_eq!(
                transaction.journal.as_slice(),
                &[Undo::ModifyAtomField {
                    id: AtomId(0),
                    change: AtomFieldChange::Charge {
                        old: NumForm::Lit(1),
                        new: NumForm::default(),
                    },
                }],
            );
            assert_eq!(
                transaction
                    .probe()
                    .unwrap()
                    .clone()
                    .atom(AtomId(0))
                    .attributes()
                    .charge,
                NumForm::Lit(1)
            );

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_set_atom_field_error(one_atom: MoleculeEditor) {
    let mut molecule = one_atom.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            let err = transaction
                .apply(Edits::from_iter([Edit::ModifyAtomField {
                    id: AtomHandle::Id(AtomId(0)),
                    change: AtomFieldChange::Charge {
                        old: NumForm::Lit(99),
                        new: NumForm::Lit(1),
                    },
                }]))
                .unwrap_err();
            assert_eq!(
                err,
                MoleculeApplyError::Transaction(TransactionError::OldStateMismatch)
            );

            Ok(())
        },
    );
    assert_eq!(result, Err(TransactionError::Aborted.into()));
}

#[rstest]
#[case::created_out_of_range(
    Edits::from_iter([Edit::AddBonds {
        bonds: vec![AddBond {
            endpoints: [AtomHandle::New(5), AtomHandle::New(6)],
            attributes: BondForm::default(),
        }],
    }]),
    MoleculeApplyError::Transaction(TransactionError::HandleOutOfRange {
        kind: EntityKind::Atom,
        index: 5,
        count: 0,
    }),
)]
#[case::initial_out_of_range(
    Edits::from_iter([Edit::RemoveTopology {
        atoms: vec![AtomHandle::Id(AtomId(0))],
        bonds: Vec::new(),
    }]),
    MoleculeApplyError::Transaction(TransactionError::HandleOutOfRange {
        kind: EntityKind::Atom,
        index: 0,
        count: 0,
    }),
)]
fn test_transaction_apply_handle_error(
    empty: MoleculeEditor,
    #[case] edits: Edits,
    #[case] expected: MoleculeApplyError,
) {
    let before = empty.clone().build();
    let mut molecule = empty.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            assert_eq!(transaction.apply(edits).unwrap_err(), expected);
            assert_eq!(transaction.probe().unwrap(), &before);

            Ok(())
        },
    );
    assert_eq!(result, Err(TransactionError::Aborted.into()));
}

#[rstest]
#[case::initial(
    1,
    Edits::from_iter([
        Edit::RemoveTopology {
            atoms: vec![AtomHandle::Id(AtomId(0))],
            bonds: Vec::new(),
        },
        Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(1),
            },
        },
    ]),
)]
#[case::created(
    0,
    Edits::from_iter([
        Edit::AddAtoms {
            atoms: vec![AtomForm::from_element(Element::C)],
        },
        Edit::RemoveTopology {
            atoms: vec![AtomHandle::New(0)],
            bonds: Vec::new(),
        },
        Edit::ModifyAtomField {
            id: AtomHandle::New(0),
            change: AtomFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(1),
            },
        },
    ]),
)]
fn test_transaction_apply_handle_removed_error(
    #[case] initial_atom_count: usize,
    #[case] edits: Edits,
) {
    let mut editor = Molecule::default().edit();
    for _ in 0..initial_atom_count {
        editor.add_atom(AtomForm::from_element(Element::C));
    }
    let before = editor.clone().build();
    let mut molecule = editor.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            assert_eq!(
                transaction.apply(edits).unwrap_err(),
                MoleculeApplyError::Transaction(TransactionError::HandleRemoved {
                    kind: EntityKind::Atom,
                    index: 0,
                })
            );
            assert_eq!(transaction.probe().unwrap(), &before);

            Ok(())
        },
    );
    assert_eq!(result, Err(TransactionError::Aborted.into()));
}

#[rstest]
fn test_transaction_apply_handles_initial() {
    let mut editor = Molecule::default().edit();
    editor.add_atom(AtomForm::from_element(Element::C));
    editor.add_atom(AtomForm::from_element(Element::N));
    editor.add_atom(AtomForm::from_element(Element::O));
    let edits = Edits::from_iter([
        Edit::RemoveTopology {
            atoms: vec![AtomHandle::Id(AtomId(0))],
            bonds: Vec::new(),
        },
        Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(1)),
            change: AtomFieldChange::Element {
                old: ElementForm::Lit(Element::N),
                new: ElementForm::Lit(Element::F),
            },
        },
        Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(2)),
            change: AtomFieldChange::Element {
                old: ElementForm::Lit(Element::O),
                new: ElementForm::Lit(Element::Cl),
            },
        },
    ]);
    let mut molecule = editor.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction.apply(edits).unwrap();

            assert_eq!(
                (0..transaction.molecule.atoms().count())
                    .map(|index| transaction
                        .molecule
                        .atom(AtomId(index as u32))
                        .attributes()
                        .element
                        .clone())
                    .collect::<Vec<_>>(),
                vec![ElementForm::Lit(Element::F), ElementForm::Lit(Element::Cl)]
            );

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_handles_created() {
    let editor = Molecule::default().edit();
    let mut edits = Edits::new();
    let atoms = edits.add_atoms([
        AtomForm::from_element(Element::C),
        AtomForm::from_element(Element::N),
    ]);
    edits.remove_atom(atoms[0].clone());
    edits.push(Edit::ModifyAtomField {
        id: atoms[1].clone(),
        change: AtomFieldChange::Charge {
            old: NumForm::default(),
            new: NumForm::Lit(1),
        },
    });
    let mut molecule = editor.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction.apply(edits).unwrap();

            assert_eq!(transaction.molecule.atoms().count(), 1);
            assert_eq!(
                (
                    transaction
                        .molecule
                        .atom(AtomId(0))
                        .attributes()
                        .element
                        .clone(),
                    transaction
                        .molecule
                        .atom(AtomId(0))
                        .attributes()
                        .charge
                        .clone(),
                ),
                (ElementForm::Lit(Element::N), NumForm::Lit(1))
            );

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_handles_reuse() {
    let editor = Molecule::default().edit();
    let mut edits = Edits::new();
    let removed = edits.add_atom(AtomForm::from_element(Element::C));
    edits.remove_atom(removed);
    let surviving = edits.add_atom(AtomForm::from_element(Element::N));
    edits.push(Edit::ModifyAtomField {
        id: surviving,
        change: AtomFieldChange::Charge {
            old: NumForm::default(),
            new: NumForm::Lit(-1),
        },
    });
    let mut molecule = editor.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction.apply(edits).unwrap();

            assert_eq!(transaction.molecule.atoms().count(), 1);
            assert_eq!(
                (
                    transaction
                        .molecule
                        .atom(AtomId(0))
                        .attributes()
                        .element
                        .clone(),
                    transaction
                        .molecule
                        .atom(AtomId(0))
                        .attributes()
                        .charge
                        .clone(),
                ),
                (ElementForm::Lit(Element::N), NumForm::Lit(-1))
            );

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_handles_per_kind() {
    let editor = Molecule::default().edit();
    let before = editor.clone().build();
    let mut edits = Edits::new();
    let atoms = edits.add_atoms([
        AtomForm::from_element(Element::C),
        AtomForm::from_element(Element::N),
        AtomForm::from_element(Element::O),
        AtomForm::from_element(Element::F),
    ]);
    let bonds = edits.add_bonds([
        AddBond {
            endpoints: [atoms[0].clone(), atoms[1].clone()],
            attributes: BondForm::from_order(1),
        },
        AddBond {
            endpoints: [atoms[1].clone(), atoms[2].clone()],
            attributes: BondForm::from_order(1),
        },
        AddBond {
            endpoints: [atoms[2].clone(), atoms[3].clone()],
            attributes: BondForm::from_order(1),
        },
    ]);
    let dative = edits.add_dative_bond(
        vec![atoms[0].clone()],
        atoms[1].clone(),
        DativeBondForm::from_order(1),
    );
    let aromatic = edits.add_aromatic_system(
        vec![atoms[0].clone(), atoms[1].clone()],
        AromaticSystemForm::default(),
    );
    let multicenter = edits.add_multicenter_bond(
        vec![atoms[0].clone(), atoms[1].clone()],
        MulticenterBondForm::default(),
    );
    let noncovalent = edits.add_noncovalent_bond(
        [atoms[0].clone(), atoms[1].clone()],
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
    );
    let stereo_atom = edits.add_stereo_atom(
        atoms[1].clone(),
        vec![
            (atoms[0].clone(), StereoLigandKind::Atom),
            (atoms[2].clone(), StereoLigandKind::Atom),
            (atoms[3].clone(), StereoLigandKind::Atom),
            (atoms[1].clone(), StereoLigandKind::ImplicitHydrogen),
        ],
        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
    );
    let stereo_bond = edits.add_stereo_bond(
        bonds[1].clone(),
        vec![
            (atoms[0].clone(), StereoLigandKind::Atom),
            (atoms[1].clone(), StereoLigandKind::ImplicitHydrogen),
            (atoms[3].clone(), StereoLigandKind::Atom),
            (atoms[2].clone(), StereoLigandKind::ImplicitHydrogen),
        ],
        StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
    );
    edits.push(Edit::ModifyAtomField {
        id: atoms[0].clone(),
        change: AtomFieldChange::Charge {
            old: NumForm::default(),
            new: NumForm::Lit(1),
        },
    });
    edits.push(Edit::ModifyBondField {
        id: bonds[0].clone(),
        change: BondFieldChange::Order {
            old: NumForm::Lit(1),
            new: NumForm::Lit(2),
        },
    });
    edits.push(Edit::ModifyDativeBondField {
        id: dative,
        change: DativeBondFieldChange::Order {
            old: NumForm::Lit(1),
            new: NumForm::Lit(2),
        },
    });
    edits.push(Edit::ModifyAromaticSystemField {
        id: aromatic,
        change: AromaticSystemFieldChange::Charge {
            old: NumForm::default(),
            new: NumForm::Lit(1),
        },
    });
    edits.push(Edit::ModifyMulticenterBondField {
        id: multicenter,
        change: MulticenterBondFieldChange::Charge {
            old: NumForm::default(),
            new: NumForm::Lit(-1),
        },
    });
    edits.push(Edit::ModifyNoncovalentBondField {
        id: noncovalent,
        change: NoncovalentBondFieldChange::Kind {
            old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
            new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
        },
    });
    edits.push(Edit::ModifyStereoAtomField {
        id: stereo_atom,
        change: StereoAtomFieldChange::Configuration {
            old: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            new: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        },
    });
    edits.push(Edit::ModifyStereoBondField {
        id: stereo_bond,
        change: StereoBondFieldChange::Configuration {
            old: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(1)),
            new: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0)),
        },
    });
    let mut molecule = editor.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction.apply(edits).unwrap();

            assert_eq!(
                transaction.molecule.atom(AtomId(0)).attributes().charge,
                NumForm::Lit(1)
            );
            assert_eq!(
                transaction.molecule.bond(BondId(0)).attributes().order,
                NumForm::Lit(2)
            );
            assert_eq!(
                transaction
                    .molecule
                    .dative_bond(DativeBondId(0))
                    .attributes()
                    .order,
                NumForm::Lit(2)
            );
            assert_eq!(
                transaction
                    .molecule
                    .aromatic_system(AromaticSystemId(0))
                    .attributes()
                    .charge,
                NumForm::Lit(1)
            );
            assert_eq!(
                transaction
                    .molecule
                    .multicenter_bond(MulticenterBondId(0))
                    .attributes()
                    .charge,
                NumForm::Lit(-1)
            );
            assert_eq!(
                transaction
                    .molecule
                    .noncovalent_bond(NoncovalentBondId(0))
                    .attributes()
                    .kind,
                NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic)
            );
            assert_eq!(
                transaction
                    .molecule
                    .stereo_atom(StereoAtomId(0))
                    .attributes()
                    .configuration,
                StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0))
            );
            assert_eq!(
                transaction
                    .molecule
                    .stereo_bond(StereoBondId(0))
                    .attributes()
                    .configuration,
                StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0))
            );

            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[rstest]
#[case::first(0)]
#[case::middle(1)]
#[case::last(2)]
fn test_transaction_apply_add_bonds_error(
    diatomic: MoleculeEditor,
    #[case] invalid_position: usize,
) {
    let before = diatomic.clone().build();
    let mut bonds = vec![
        AddBond {
            endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
            attributes: BondForm::from_order(1),
        },
        AddBond {
            endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
            attributes: BondForm::from_order(2),
        },
        AddBond {
            endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
            attributes: BondForm::from_order(3),
        },
    ];
    bonds[invalid_position].endpoints[1] = AtomHandle::Id(AtomId(9));
    let mut molecule = diatomic.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            assert_eq!(
                transaction
                    .apply(Edits::from_iter([Edit::AddBonds { bonds }]))
                    .unwrap_err(),
                MoleculeApplyError::Transaction(TransactionError::HandleOutOfRange {
                    kind: EntityKind::Atom,
                    index: 9,
                    count: 2,
                })
            );
            assert_eq!(transaction.probe().unwrap(), &before);

            Ok(())
        },
    );
    assert_eq!(result, Err(TransactionError::Aborted.into()));
}

#[rstest]
fn test_transaction_apply_remove_topology(diatomic: MoleculeEditor) {
    let mut molecule = diatomic.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::RemoveTopology {
                    atoms: Vec::new(),
                    bonds: vec![BondHandle::Id(BondId(0))],
                }]))
                .unwrap();

            assert_eq!(transaction.molecule.bonds().count(), 0);
            let [Undo::RestoreRemovedTopology {
                atoms,
                bonds,
                overlays,
                compaction,
                ..
            }] = transaction.journal.as_slice()
            else {
                panic!("RemoveTopology should produce one topology-restore undo")
            };
            assert!(atoms.is_empty());
            assert_eq!(
                bonds.iter().map(|b| b.id).collect::<Vec<_>>(),
                vec![BondId(0)]
            );
            assert!(overlays.dative_bonds.is_empty());
            assert_eq!(compaction.compact_bond(BondId(0)), None);

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
#[case::atom(Edit::ModifyAtomConstraint { id: AtomHandle::Id(AtomId(0)), old: None, new: None })]
#[case::bond(Edit::ModifyBondConstraint { id: BondHandle::Id(BondId(0)), old: None, new: None })]
#[case::dative_bond(Edit::ModifyDativeBondConstraint { id: DativeBondHandle::Id(DativeBondId(0)), old: None, new: None })]
#[case::aromatic_system(Edit::ModifyAromaticSystemConstraint { id: AromaticSystemHandle::Id(AromaticSystemId(0)), old: None, new: None })]
#[case::multicenter_bond(Edit::ModifyMulticenterBondConstraint { id: MulticenterBondHandle::Id(MulticenterBondId(0)), old: None, new: None })]
#[case::noncovalent_bond(Edit::ModifyNoncovalentBondConstraint { id: NoncovalentBondHandle::Id(NoncovalentBondId(0)), old: None, new: None })]
#[case::stereo_atom(Edit::ModifyStereoAtomConstraint { id: StereoAtomHandle::Id(StereoAtomId(0)), kind: None, old: None, new: None })]
#[case::stereo_bond(Edit::ModifyStereoBondConstraint { id: StereoBondHandle::Id(StereoBondId(0)), kind: None, old: None, new: None })]
fn test_transaction_apply_constraint_identity(
    batched_overlays: MoleculeEditor,
    #[case] edit: Edit,
) {
    let expected = batched_overlays.snapshot().unwrap();
    let mut molecule = batched_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction.apply(Edits::from_iter([edit])).unwrap();

            assert_eq!(transaction.journal.as_slice(), &[]);
            assert_eq!(transaction.probe().cloned(), Ok(expected));

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_add_atom_constraint(one_atom: MoleculeEditor) {
    let mut molecule = one_atom.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([
                    Edit::ModifyAtomConstraint {
                        id: AtomHandle::Id(AtomId(0)),
                        old: None,
                        new: Some(AtomConstraintForm::ring_membership(RingScope::Size(5), 1)),
                    },
                    Edit::ModifyAtomConstraint {
                        id: AtomHandle::Id(AtomId(0)),
                        old: None,
                        new: Some(AtomConstraintForm::ring_membership(RingScope::Size(6), 1)),
                    },
                ]))
                .unwrap();
            let next = transaction.probe().unwrap().clone();
            let cs: Vec<_> = next
                .atom(AtomId(0))
                .attributes()
                .constraints
                .iter()
                .cloned()
                .collect();
            assert_eq!(
                cs,
                vec![
                    AtomConstraintForm::ring_membership(RingScope::Size(5), 1),
                    AtomConstraintForm::ring_membership(RingScope::Size(6), 1),
                ]
            );

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
#[case::singleton_set(NumForm::Lit(1), NumForm::lit_set([1]))]
fn test_transaction_apply_modify_atom_field_canonical(
    mut one_atom: MoleculeEditor,
    #[case] current: NumForm,
    #[case] old: NumForm,
) {
    // The modify's recorded `old` is equivalent to — but structurally distinct from — the
    // stored charge, so the old-state check passes (structural `!=` would raise `OldStateMismatch`).
    one_atom.atom_mut(AtomId(0)).attributes_mut().charge = current;
    let mut molecule = one_atom.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::ModifyAtomField {
                    id: AtomHandle::Id(AtomId(0)),
                    change: AtomFieldChange::Charge {
                        old,
                        new: NumForm::Lit(2),
                    },
                }]))
                .unwrap();
            assert_eq!(
                transaction
                    .molecule
                    .atom_mut(AtomId(0))
                    .attributes_mut()
                    .charge,
                NumForm::Lit(2)
            );

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_modify_atom_constraint_absent_error(one_atom: MoleculeEditor) {
    let mut molecule = one_atom.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            let err = transaction
                .apply(Edits::from_iter([Edit::ModifyAtomConstraint {
                    id: AtomHandle::Id(AtomId(0)),
                    old: Some(AtomConstraintForm::ring_membership(RingScope::Size(5), 1)),
                    new: None,
                }]))
                .unwrap_err();
            assert_eq!(
                err,
                MoleculeApplyError::Transaction(TransactionError::OldStateMismatch)
            );

            Ok(())
        },
    );
    assert_eq!(result, Err(TransactionError::Aborted.into()));
}

#[rstest]
#[case::introduce(None, Some(AtomConstraintForm::valence(4)), Some(NumForm::Lit(4)))]
#[case::replace(
    Some(AtomConstraintForm::valence(3)),
    Some(AtomConstraintForm::valence(4)),
    Some(NumForm::Lit(4))
)]
#[case::remove(Some(AtomConstraintForm::valence(3)), None, None)]
fn test_transaction_apply_set_atom_constraint(
    one_atom: MoleculeEditor,
    #[case] old: Option<AtomConstraintForm>,
    #[case] new: Option<AtomConstraintForm>,
    #[case] expected: Option<NumForm>,
) {
    let mut molecule = one_atom.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            if let Some(c) = old.clone() {
                transaction
                    .molecule
                    .atom_mut(AtomId(0))
                    .attributes_mut()
                    .constraints
                    .set(c);
            }
            transaction
                .apply(Edits::from_iter([Edit::ModifyAtomConstraint {
                    id: AtomHandle::Id(AtomId(0)),
                    old,
                    new,
                }]))
                .unwrap();
            assert_eq!(
                transaction
                    .molecule
                    .atom_mut(AtomId(0))
                    .attributes_mut()
                    .constraints
                    .valence(),
                expected.as_ref()
            );

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_set_bond_constraint(diatomic: MoleculeEditor) {
    let mut molecule = diatomic.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(0)),
                    old: None,
                    new: Some(BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
                }]))
                .unwrap();
            assert!(transaction
                .molecule
                .bond_mut(BondId(0))
                .attributes_mut()
                .constraints
                .iter()
                .any(|c| *c == BondConstraintForm::Aromatic(BooleanForm::Lit(true))));

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_add_molecule_constraint(empty: MoleculeEditor) {
    let c = Constraint::Molecule(MoleculeConstraint::Connected { atoms: None });
    let mut molecule = empty.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::AddMoleculeConstraint {
                    constraint: c.clone().into(),
                }]))
                .unwrap();
            assert_eq!(transaction.molecule.constraints_mut().as_slice(), &[c]);

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_remove_molecule_constraint(mut empty: MoleculeEditor) {
    let c = Constraint::Molecule(MoleculeConstraint::Connected { atoms: None });
    empty.constraints_mut().push(c.clone());
    let mut molecule = empty.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::RemoveMoleculeConstraint {
                    constraint: c.clone().into(),
                }]))
                .unwrap();
            assert!(transaction.molecule.constraints_mut().as_slice().is_empty());

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_remove_molecule_constraint_absent_error(mut empty: MoleculeEditor) {
    let c = Constraint::Molecule(MoleculeConstraint::Connected { atoms: None });
    empty.constraints_mut().push(c.clone());
    let mut molecule = empty.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            let err = transaction
                .apply(Edits::from_iter([Edit::RemoveMoleculeConstraint {
                    constraint: Constraint::Molecule(MoleculeConstraint::ChargeSum {
                        atoms: None,
                        sum: NumForm::Lit(0),
                    })
                    .into(),
                }]))
                .unwrap_err();
            assert_eq!(
                err,
                MoleculeApplyError::Transaction(TransactionError::MissingEntry)
            );
            assert_eq!(transaction.molecule.constraints_mut().as_slice(), &[c]);

            Ok(())
        },
    );
    assert_eq!(result, Err(TransactionError::Aborted.into()));
}

#[rstest]
fn test_transaction_apply_molecule_constraint_initial(batched_overlays: MoleculeEditor) {
    let constraint = Constraint::And(vec![
        Constraint::Atom(AtomId(1), AtomConstraintForm::valence(3_i64)),
        Constraint::Bond(BondId(1), BondConstraintForm::aromatic(true)),
        Constraint::DativeBond(DativeBondId(1), DativeBondConstraintForm::aromatic(true)),
        Constraint::AromaticSystem(
            AromaticSystemId(1),
            AromaticSystemConstraintForm::electron_count(6_i64),
        ),
        Constraint::MulticenterBond(
            MulticenterBondId(1),
            MulticenterBondConstraintForm::electron_count(2_i64),
        ),
        Constraint::NoncovalentBond(
            NoncovalentBondId(1),
            NoncovalentBondConstraintForm::intramolecular(true),
        ),
        Constraint::StereoAtom(
            StereoAtomId(1),
            StereoKind::Tetrahedral,
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
        ),
        Constraint::StereoBond(
            StereoBondId(1),
            StereoKind::CisTrans,
            StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
        ),
    ]);
    let before = batched_overlays.clone().build();
    let mut edits = Edits::new();
    edits.add_molecule_constraint(constraint.clone().into());
    let mut molecule = batched_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction.apply(edits).unwrap();
            assert_eq!(transaction.molecule.constraints().as_slice(), &[constraint]);

            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[rstest]
fn test_transaction_apply_molecule_constraint_created(empty: MoleculeEditor) {
    let before = empty.clone().build();
    let mut edits = Edits::new();
    let atoms = edits.add_atoms([
        AtomForm::from_element(Element::C),
        AtomForm::from_element(Element::N),
    ]);
    let bond = edits.add_bond(atoms[0].clone(), atoms[1].clone(), BondForm::from_order(1));
    let dative = edits.add_dative_bond(
        vec![atoms[0].clone()],
        atoms[1].clone(),
        DativeBondForm::from_order(1),
    );
    let aromatic = edits.add_aromatic_system(atoms.clone(), AromaticSystemForm::default());
    let multicenter = edits.add_multicenter_bond(atoms.clone(), MulticenterBondForm::default());
    let noncovalent = edits.add_noncovalent_bond(
        [atoms[0].clone(), atoms[1].clone()],
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
    );
    let stereo_atom = edits.add_stereo_atom(
        atoms[0].clone(),
        vec![
            (atoms[1].clone(), StereoLigandKind::Atom),
            (atoms[0].clone(), StereoLigandKind::ImplicitHydrogen),
            (atoms[0].clone(), StereoLigandKind::LonePair),
            (atoms[1].clone(), StereoLigandKind::ImplicitHydrogen),
        ],
        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
    );
    let stereo_bond = edits.add_stereo_bond(
        bond.clone(),
        vec![
            (atoms[0].clone(), StereoLigandKind::Atom),
            (atoms[0].clone(), StereoLigandKind::ImplicitHydrogen),
            (atoms[1].clone(), StereoLigandKind::Atom),
            (atoms[1].clone(), StereoLigandKind::ImplicitHydrogen),
        ],
        StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
    );
    let source = Constraint::And(vec![
        Constraint::Atom(AtomId(7), AtomConstraintForm::valence(3_i64)),
        Constraint::Bond(BondId(7), BondConstraintForm::aromatic(true)),
        Constraint::DativeBond(DativeBondId(7), DativeBondConstraintForm::aromatic(true)),
        Constraint::AromaticSystem(
            AromaticSystemId(7),
            AromaticSystemConstraintForm::electron_count(6_i64),
        ),
        Constraint::MulticenterBond(
            MulticenterBondId(7),
            MulticenterBondConstraintForm::electron_count(2_i64),
        ),
        Constraint::NoncovalentBond(
            NoncovalentBondId(7),
            NoncovalentBondConstraintForm::intramolecular(true),
        ),
        Constraint::StereoAtom(
            StereoAtomId(7),
            StereoKind::Tetrahedral,
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
        ),
        Constraint::StereoBond(
            StereoBondId(7),
            StereoKind::CisTrans,
            StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
        ),
        Constraint::Relational(RelationalConstraint::DativeBondParallels {
            dative: DativeBondId(7),
            parallel: BondId(7),
        }),
    ]);
    let mappings = HashMap::from([
        (
            Entity::Atom(AtomId(7)),
            EntityHandle::Atom(atoms[0].clone()),
        ),
        (Entity::Bond(BondId(7)), EntityHandle::Bond(bond)),
        (
            Entity::DativeBond(DativeBondId(7)),
            EntityHandle::DativeBond(dative),
        ),
        (
            Entity::AromaticSystem(AromaticSystemId(7)),
            EntityHandle::AromaticSystem(aromatic),
        ),
        (
            Entity::MulticenterBond(MulticenterBondId(7)),
            EntityHandle::MulticenterBond(multicenter),
        ),
        (
            Entity::NoncovalentBond(NoncovalentBondId(7)),
            EntityHandle::NoncovalentBond(noncovalent),
        ),
        (
            Entity::StereoAtom(StereoAtomId(7)),
            EntityHandle::StereoAtom(stereo_atom),
        ),
        (
            Entity::StereoBond(StereoBondId(7)),
            EntityHandle::StereoBond(stereo_bond),
        ),
    ]);
    edits.add_molecule_constraint(
        ConstraintEdit::new(source, |entity| mappings.get(&entity).cloned()).unwrap(),
    );
    let expected = Constraint::And(vec![
        Constraint::Atom(AtomId(0), AtomConstraintForm::valence(3_i64)),
        Constraint::Bond(BondId(0), BondConstraintForm::aromatic(true)),
        Constraint::DativeBond(DativeBondId(0), DativeBondConstraintForm::aromatic(true)),
        Constraint::AromaticSystem(
            AromaticSystemId(0),
            AromaticSystemConstraintForm::electron_count(6_i64),
        ),
        Constraint::MulticenterBond(
            MulticenterBondId(0),
            MulticenterBondConstraintForm::electron_count(2_i64),
        ),
        Constraint::NoncovalentBond(
            NoncovalentBondId(0),
            NoncovalentBondConstraintForm::intramolecular(true),
        ),
        Constraint::StereoAtom(
            StereoAtomId(0),
            StereoKind::Tetrahedral,
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
        ),
        Constraint::StereoBond(
            StereoBondId(0),
            StereoKind::CisTrans,
            StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
        ),
        Constraint::Relational(RelationalConstraint::DativeBondParallels {
            dative: DativeBondId(0),
            parallel: BondId(0),
        }),
    ]);
    let mut molecule = empty.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction.apply(edits).unwrap();
            assert_eq!(transaction.molecule.constraints().as_slice(), &[expected]);

            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[rstest]
fn test_transaction_apply_molecule_constraint_compaction(mut batched_overlays: MoleculeEditor) {
    let removed = Constraint::AromaticSystem(
        AromaticSystemId(1),
        AromaticSystemConstraintForm::electron_count(6_i64),
    );
    let added = Constraint::AromaticSystem(
        AromaticSystemId(1),
        AromaticSystemConstraintForm::electron_count(4_i64),
    );
    batched_overlays.constraints_mut().push(removed.clone());
    let before = batched_overlays.clone().build();
    let mut edits = Edits::from_iter([Edit::RemoveAromaticSystems {
        removes: vec![(
            AromaticSystemHandle::Id(AromaticSystemId(0)),
            vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
            AromaticSystemForm::default(),
        )],
    }]);
    edits.add_molecule_constraint(added.into());
    edits.remove_molecule_constraint(removed.into());
    let mut molecule = batched_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction.apply(edits).unwrap();
            assert_eq!(
                transaction.molecule.constraints().as_slice(),
                &[Constraint::AromaticSystem(
                    AromaticSystemId(0),
                    AromaticSystemConstraintForm::electron_count(4_i64),
                )],
            );

            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[rstest]
#[case::forward(
    0,
    Edits::from_iter([Edit::AddMoleculeConstraint {
        constraint: ConstraintEdit::new(
            Constraint::Atom(AtomId(7), AtomConstraintForm::valence(3_i64)),
            |_| Some(EntityHandle::Atom(AtomHandle::New(0))),
        ).unwrap(),
    }]),
    MoleculeApplyError::Transaction(TransactionError::HandleOutOfRange { kind: EntityKind::Atom, index: 0, count: 0 }),
)]
#[case::removed(
    1,
    Edits::from_iter([
        Edit::RemoveTopology {
            atoms: vec![AtomHandle::Id(AtomId(0))],
            bonds: Vec::new(),
        },
        Edit::AddMoleculeConstraint {
            constraint: ConstraintEdit::new(
                Constraint::Atom(AtomId(7), AtomConstraintForm::valence(3_i64)),
                |_| Some(EntityHandle::Atom(AtomHandle::Id(AtomId(0)))),
            ).unwrap(),
        },
    ]),
    MoleculeApplyError::Transaction(TransactionError::HandleRemoved { kind: EntityKind::Atom, index: 0 }),
)]
fn test_transaction_apply_molecule_constraint_error(
    #[case] initial_atom_count: usize,
    #[case] edits: Edits,
    #[case] expected: MoleculeApplyError,
) {
    let mut editor = Molecule::default().edit();
    for _ in 0..initial_atom_count {
        editor.add_atom(AtomForm::from_element(Element::C));
    }
    let before = editor.clone().build();
    let mut molecule = editor.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            assert_eq!(transaction.apply(edits), Err(expected));
            assert_eq!(transaction.probe().unwrap(), &before);

            Ok(())
        },
    );
    assert_eq!(result, Err(TransactionError::Aborted.into()));
}

#[rstest]
fn test_transaction_apply_remove_topology_atom_error(one_atom: MoleculeEditor) {
    let before = one_atom.clone().build();
    let mut edits = Edits::new();
    edits.remove_atom(AtomHandle::Id(AtomId(9)));
    let mut molecule = one_atom.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            let err = transaction.apply(edits).unwrap_err();
            assert_eq!(
                err,
                MoleculeApplyError::Transaction(TransactionError::HandleOutOfRange {
                    kind: EntityKind::Atom,
                    index: 9,
                    count: 1,
                })
            );
            assert_eq!(transaction.probe().unwrap(), &before);

            Ok(())
        },
    );
    assert_eq!(result, Err(TransactionError::Aborted.into()));
}

#[rstest]
fn test_transaction_apply_remove_topology_bond_error(diatomic: MoleculeEditor) {
    let before = diatomic.clone().build();
    let mut edits = Edits::new();
    edits.remove_bond(BondHandle::Id(BondId(9)));
    let mut molecule = diatomic.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            let err = transaction.apply(edits).unwrap_err();
            assert_eq!(
                err,
                MoleculeApplyError::Transaction(TransactionError::HandleOutOfRange {
                    kind: EntityKind::Bond,
                    index: 9,
                    count: 1,
                })
            );
            assert_eq!(transaction.probe().unwrap(), &before);

            Ok(())
        },
    );
    assert_eq!(result, Err(TransactionError::Aborted.into()));
}

#[rstest]
fn test_transaction_apply_add_dative_bond_acceptor_error(one_atom: MoleculeEditor) {
    let mut molecule = one_atom.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            let err = transaction
                .apply(Edits::from_iter([Edit::AddDativeBond {
                    donors: vec![],
                    acceptor: AtomHandle::Id(AtomId(9)),
                    attributes: DativeBondForm::from_order(1),
                }]))
                .unwrap_err();
            assert_eq!(
                err,
                MoleculeApplyError::Transaction(TransactionError::HandleOutOfRange {
                    kind: EntityKind::Atom,
                    index: 9,
                    count: 1,
                })
            );

            Ok(())
        },
    );
    assert_eq!(result, Err(TransactionError::Aborted.into()));
}

#[rstest]
fn test_transaction_apply_set_bond_field(diatomic: MoleculeEditor) {
    let mut molecule = diatomic.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::ModifyBondField {
                    id: BondHandle::Id(BondId(0)),
                    change: BondFieldChange::Order {
                        old: NumForm::Lit(1),
                        new: NumForm::Lit(2),
                    },
                }]))
                .unwrap();
            assert_eq!(
                transaction.molecule.bond(BondId(0)).attributes().order,
                NumForm::Lit(2)
            );

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_set_bond_field_error(diatomic: MoleculeEditor) {
    let mut molecule = diatomic.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            let err = transaction
                .apply(Edits::from_iter([Edit::ModifyBondField {
                    id: BondHandle::Id(BondId(0)),
                    change: BondFieldChange::Order {
                        old: NumForm::Lit(99),
                        new: NumForm::Lit(2),
                    },
                }]))
                .unwrap_err();
            assert_eq!(
                err,
                MoleculeApplyError::Transaction(TransactionError::OldStateMismatch)
            );

            Ok(())
        },
    );
    assert_eq!(result, Err(TransactionError::Aborted.into()));
}

#[fixture]
fn stereo_atom_skeleton() -> MoleculeEditor {
    let mut b = Molecule::default().edit();
    for el in [Element::C, Element::F, Element::Cl, Element::Br, Element::I] {
        b.add_atom(AtomForm::from_element(el));
    }
    for t in 1u32..=4 {
        b.add_bond(AtomId(0), AtomId(t), BondForm::from_order(1));
    }
    b
}

#[fixture]
fn tetrahedral_ligands() -> Vec<StereoLigand> {
    (1u32..=4)
        .map(|t| StereoLigand::new(AtomId(t), StereoLigandKind::Atom))
        .collect()
}

#[rstest]
fn test_transaction_apply_add_stereo_atom(stereo_atom_skeleton: MoleculeEditor) {
    let before = stereo_atom_skeleton.clone().build();
    let mut molecule = stereo_atom_skeleton.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::AddStereoAtom {
                    site: AtomHandle::Id(AtomId(0)),
                    ligands: (1u32..=4)
                        .map(|t| (AtomHandle::Id(AtomId(t)), StereoLigandKind::Atom))
                        .collect(),
                    attributes: StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                }]))
                .unwrap();
            assert_eq!(transaction.molecule.stereo_atoms().count(), 1);
            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[rstest]
fn test_transaction_apply_remove_stereo_atom(mut stereo_atom_skeleton: MoleculeEditor) {
    stereo_atom_skeleton.add_stereo_atom(
        AtomId(0),
        &tetrahedral_ligands(),
        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
    );
    let before = stereo_atom_skeleton.clone().build();
    let mut molecule = stereo_atom_skeleton.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::RemoveStereoAtoms {
                    removes: vec![(
                        StereoAtomHandle::Id(StereoAtomId(0)),
                        AtomHandle::Id(AtomId(0)),
                        (1u32..=4)
                            .map(|t| (AtomHandle::Id(AtomId(t)), StereoLigandKind::Atom))
                            .collect(),
                        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                    )],
                }]))
                .unwrap();
            assert_eq!(transaction.molecule.stereo_atoms().count(), 0);
            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[rstest]
fn test_transaction_apply_remove_stereo_atom_error(mut stereo_atom_skeleton: MoleculeEditor) {
    stereo_atom_skeleton.add_stereo_atom(
        AtomId(0),
        &tetrahedral_ligands(),
        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
    );
    let mut molecule = stereo_atom_skeleton.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            let err = transaction
                .apply(Edits::from_iter([Edit::RemoveStereoAtoms {
                    removes: vec![(
                        StereoAtomHandle::Id(StereoAtomId(0)),
                        AtomHandle::Id(AtomId(0)),
                        (1u32..=4)
                            .map(|t| (AtomHandle::Id(AtomId(t)), StereoLigandKind::Atom))
                            .collect(),
                        // Wrong recorded coset (Th0 vs the stored Th1).
                        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                    )],
                }]))
                .unwrap_err();
            assert_eq!(
                err,
                MoleculeApplyError::Transaction(TransactionError::OldStateMismatch)
            );

            Ok(())
        },
    );
    assert_eq!(result, Err(TransactionError::Aborted.into()));
}

#[rstest]
fn test_transaction_apply_topology_removal_restores_stereo_atom(
    mut stereo_atom_skeleton: MoleculeEditor,
) {
    stereo_atom_skeleton.add_stereo_atom(
        AtomId(0),
        &tetrahedral_ligands(),
        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
    );
    let before = stereo_atom_skeleton.clone().build();
    let mut molecule = stereo_atom_skeleton.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            // Removing a ligand atom cascades the stereo element away.
            transaction
                .apply(Edits::from_iter([Edit::RemoveTopology {
                    atoms: vec![AtomHandle::Id(AtomId(1))],
                    bonds: Vec::new(),
                }]))
                .unwrap();
            assert_eq!(transaction.molecule.stereo_atoms().count(), 0);
            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[fixture]
fn stereo_bond_skeleton() -> MoleculeEditor {
    let mut b = Molecule::default().edit();
    for _ in 0..4 {
        b.add_atom(AtomForm::from_element(Element::C));
    }
    b.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));
    b.add_bond(AtomId(1), AtomId(2), BondForm::from_order(2));
    b.add_bond(AtomId(2), AtomId(3), BondForm::from_order(1));
    b
}

#[rstest]
fn test_transaction_apply_add_stereo_bond(stereo_bond_skeleton: MoleculeEditor) {
    let before = stereo_bond_skeleton.clone().build();
    let mut molecule = stereo_bond_skeleton.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::AddStereoBond {
                    site: BondHandle::Id(BondId(1)),
                    ligands: vec![
                        (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                        (
                            AtomHandle::Id(AtomId(1)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                        (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                        (
                            AtomHandle::Id(AtomId(2)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                    ],
                    attributes: StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                }]))
                .unwrap();
            assert_eq!(transaction.molecule.stereo_bonds().count(), 1);
            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[rstest]
fn test_transaction_apply_remove_stereo_bond(mut stereo_bond_skeleton: MoleculeEditor) {
    stereo_bond_skeleton.add_stereo_bond(
        BondId(1),
        &[
            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
        ],
        StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
    );
    let before = stereo_bond_skeleton.clone().build();
    let mut molecule = stereo_bond_skeleton.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::RemoveStereoBonds {
                    removes: vec![(
                        StereoBondHandle::Id(StereoBondId(0)),
                        BondHandle::Id(BondId(1)),
                        vec![
                            (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                            (
                                AtomHandle::Id(AtomId(1)),
                                StereoLigandKind::ImplicitHydrogen,
                            ),
                            (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                            (
                                AtomHandle::Id(AtomId(2)),
                                StereoLigandKind::ImplicitHydrogen,
                            ),
                        ],
                        StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                    )],
                }]))
                .unwrap();
            assert_eq!(transaction.molecule.stereo_bonds().count(), 0);
            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[rstest]
fn test_transaction_apply_remove_stereo_bond_error(mut stereo_bond_skeleton: MoleculeEditor) {
    stereo_bond_skeleton.add_stereo_bond(
        BondId(1),
        &[
            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
        ],
        StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
    );
    let mut molecule = stereo_bond_skeleton.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            let error = transaction
                .apply(Edits::from_iter([Edit::RemoveStereoBonds {
                    removes: vec![(
                        StereoBondHandle::Id(StereoBondId(0)),
                        BondHandle::Id(BondId(1)),
                        vec![
                            (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                            (
                                AtomHandle::Id(AtomId(1)),
                                StereoLigandKind::ImplicitHydrogen,
                            ),
                            (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                            (
                                AtomHandle::Id(AtomId(2)),
                                StereoLigandKind::ImplicitHydrogen,
                            ),
                        ],
                        // Wrong recorded coset (Ct0 vs the stored Ct1).
                        StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
                    )],
                }]))
                .unwrap_err();

            assert_eq!(
                error,
                MoleculeApplyError::Transaction(TransactionError::OldStateMismatch)
            );

            Ok(())
        },
    );
    assert_eq!(result, Err(TransactionError::Aborted.into()));
}

#[rstest]
fn test_transaction_apply_set_stereo_atom_field(mut stereo_atom_skeleton: MoleculeEditor) {
    stereo_atom_skeleton.add_stereo_atom(
        AtomId(0),
        &tetrahedral_ligands(),
        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
    );
    let before = stereo_atom_skeleton.clone().build();
    let mut molecule = stereo_atom_skeleton.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::ModifyStereoAtomField {
                    id: StereoAtomHandle::Id(StereoAtomId(0)),
                    change: StereoAtomFieldChange::Configuration {
                        old: StereoConfigurationForm::kinded(
                            StereoKind::Tetrahedral,
                            StereoCoset::Lit(1),
                        ),
                        new: StereoConfigurationForm::kinded(
                            StereoKind::Tetrahedral,
                            StereoCoset::Lit(0),
                        ),
                    },
                }]))
                .unwrap();
            assert_eq!(
                transaction
                    .molecule
                    .stereo_atom(StereoAtomId(0))
                    .attributes()
                    .configuration,
                StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0),),
            );
            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[rstest]
fn test_transaction_apply_set_stereo_atom_field_error(mut stereo_atom_skeleton: MoleculeEditor) {
    stereo_atom_skeleton.add_stereo_atom(
        AtomId(0),
        &tetrahedral_ligands(),
        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
    );
    let mut molecule = stereo_atom_skeleton.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            let err = transaction
                .apply(Edits::from_iter([Edit::ModifyStereoAtomField {
                    id: StereoAtomHandle::Id(StereoAtomId(0)),
                    change: StereoAtomFieldChange::Configuration {
                        // Wrong recorded coset (Th0 vs the stored Th1).
                        old: StereoConfigurationForm::kinded(
                            StereoKind::Tetrahedral,
                            StereoCoset::Lit(0),
                        ),
                        new: StereoConfigurationForm::kinded(
                            StereoKind::Tetrahedral,
                            StereoCoset::Lit(1),
                        ),
                    },
                }]))
                .unwrap_err();
            assert_eq!(
                err,
                MoleculeApplyError::Transaction(TransactionError::OldStateMismatch)
            );

            Ok(())
        },
    );
    assert_eq!(result, Err(TransactionError::Aborted.into()));
}

#[rstest]
fn test_transaction_apply_set_stereo_bond_field(mut stereo_bond_skeleton: MoleculeEditor) {
    stereo_bond_skeleton.add_stereo_bond(
        BondId(1),
        &[
            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
        ],
        StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
    );
    let before = stereo_bond_skeleton.clone().build();
    let mut molecule = stereo_bond_skeleton.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::ModifyStereoBondField {
                    id: StereoBondHandle::Id(StereoBondId(0)),
                    change: StereoBondFieldChange::Configuration {
                        old: StereoConfigurationForm::kinded(
                            StereoKind::CisTrans,
                            StereoCoset::Lit(1),
                        ),
                        new: StereoConfigurationForm::kinded(
                            StereoKind::CisTrans,
                            StereoCoset::Lit(0),
                        ),
                    },
                }]))
                .unwrap();
            assert_eq!(
                transaction
                    .molecule
                    .stereo_bond(StereoBondId(0))
                    .attributes()
                    .configuration,
                StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0)),
            );
            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[rstest]
fn test_transaction_apply_set_stereo_bond_field_error(mut stereo_bond_skeleton: MoleculeEditor) {
    stereo_bond_skeleton.add_stereo_bond(
        BondId(1),
        &[
            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
        ],
        StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
    );
    let mut molecule = stereo_bond_skeleton.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            let err = transaction
                .apply(Edits::from_iter([Edit::ModifyStereoBondField {
                    id: StereoBondHandle::Id(StereoBondId(0)),
                    change: StereoBondFieldChange::Configuration {
                        // Wrong recorded coset (vs the stored 1).
                        old: StereoConfigurationForm::kinded(
                            StereoKind::CisTrans,
                            StereoCoset::Lit(0),
                        ),
                        new: StereoConfigurationForm::kinded(
                            StereoKind::CisTrans,
                            StereoCoset::Lit(1),
                        ),
                    },
                }]))
                .unwrap_err();
            assert_eq!(
                err,
                MoleculeApplyError::Transaction(TransactionError::OldStateMismatch)
            );

            Ok(())
        },
    );
    assert_eq!(result, Err(TransactionError::Aborted.into()));
}

#[fixture]
fn diatomic_with_overlays() -> MoleculeEditor {
    let mut b = Molecule::default().edit();
    b.add_atom(AtomForm::from_element(Element::C));
    b.add_atom(AtomForm::from_element(Element::N));
    b.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));
    b.add_dative_bond(&[AtomId(0)], AtomId(1), DativeBondForm::from_order(1));
    b.add_aromatic_system(&[AtomId(0), AtomId(1)], AromaticSystemForm::default());
    b.add_multicenter_bond(&[AtomId(0), AtomId(1)], MulticenterBondForm::default());
    b.add_noncovalent_bond(
        [AtomId(0), AtomId(1)],
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
    );
    b
}

#[fixture]
fn batched_overlays() -> MoleculeEditor {
    let mut editor = Molecule::default().edit();
    for _ in 0..6 {
        editor.add_atom(AtomForm::from_element(Element::C));
    }
    for index in 0..3_u32 {
        let first = AtomId(index * 2);
        let second = AtomId(index * 2 + 1);
        let bond = editor.add_bond(first, second, BondForm::from_order(1));
        editor.add_dative_bond(&[first], second, DativeBondForm::from_order(1));
        editor.add_aromatic_system(&[first, second], AromaticSystemForm::default());
        editor.add_multicenter_bond(&[first, second], MulticenterBondForm::default());
        editor.add_noncovalent_bond(
            [first, second],
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        );
        editor.add_stereo_atom(
            first,
            &[
                StereoLigand::new(second, StereoLigandKind::Atom),
                StereoLigand::new(first, StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(first, StereoLigandKind::LonePair),
            ],
            StereoAtomForm::default(),
        );
        editor.add_stereo_bond(
            bond,
            &[
                StereoLigand::new(first, StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(first, StereoLigandKind::LonePair),
                StereoLigand::new(second, StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(second, StereoLigandKind::LonePair),
            ],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        );
    }
    editor
}

#[rstest]
fn test_transaction_apply_replace_dative_bond_donors(batched_overlays: MoleculeEditor) {
    let before = batched_overlays.clone().build();
    let donors = vec![AtomId(2), AtomId(4)];
    let mut molecule = batched_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::ReplaceDativeBondDonors {
                    id: DativeBondHandle::Id(DativeBondId(0)),
                    old: vec![AtomHandle::Id(AtomId(0))],
                    new: donors.iter().copied().map(AtomHandle::Id).collect(),
                }]))
                .unwrap();
            assert_eq!(
                transaction.journal.as_slice(),
                &[Undo::RestoreDativeBondDonors {
                    id: DativeBondId(0),
                    donors: vec![AtomId(0)],
                }]
            );
            assert_eq!(
                transaction
                    .molecule
                    .dative_bond(DativeBondId(0))
                    .donor_ids()
                    .collect::<Vec<_>>(),
                donors
            );
            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[rstest]
fn test_transaction_apply_replace_dative_bond_acceptor(batched_overlays: MoleculeEditor) {
    let before = batched_overlays.clone().build();
    let mut molecule = batched_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::ReplaceDativeBondAcceptor {
                    id: DativeBondHandle::Id(DativeBondId(0)),
                    old: AtomHandle::Id(AtomId(1)),
                    new: AtomHandle::Id(AtomId(3)),
                }]))
                .unwrap();
            assert_eq!(
                transaction.journal.as_slice(),
                &[Undo::RestoreDativeBondAcceptor {
                    id: DativeBondId(0),
                    acceptor: AtomId(1),
                }]
            );
            assert_eq!(
                transaction
                    .molecule
                    .dative_bond(DativeBondId(0))
                    .acceptor_id(),
                AtomId(3)
            );
            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[rstest]
fn test_transaction_apply_replace_aromatic_system_atoms(batched_overlays: MoleculeEditor) {
    let before = batched_overlays.clone().build();
    let mut molecule = batched_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::ReplaceAromaticSystemAtoms {
                    id: AromaticSystemHandle::Id(AromaticSystemId(0)),
                    old: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    new: vec![AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(0))],
                }]))
                .unwrap();
            assert_eq!(
                transaction.journal.as_slice(),
                &[Undo::RestoreAromaticSystemAtoms {
                    id: AromaticSystemId(0),
                    atoms: vec![AtomId(0), AtomId(1)],
                }]
            );
            assert_eq!(
                transaction
                    .molecule
                    .aromatic_system(AromaticSystemId(0))
                    .atom_ids()
                    .collect::<Vec<_>>(),
                vec![AtomId(1), AtomId(0)]
            );
            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[rstest]
fn test_transaction_apply_replace_multicenter_bond_atoms(batched_overlays: MoleculeEditor) {
    let before = batched_overlays.clone().build();
    let mut molecule = batched_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::ReplaceMulticenterBondAtoms {
                    id: MulticenterBondHandle::Id(MulticenterBondId(0)),
                    old: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    new: vec![AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(0))],
                }]))
                .unwrap();
            assert_eq!(
                transaction.journal.as_slice(),
                &[Undo::RestoreMulticenterBondAtoms {
                    id: MulticenterBondId(0),
                    atoms: vec![AtomId(0), AtomId(1)],
                }]
            );
            assert_eq!(
                transaction
                    .molecule
                    .multicenter_bond(MulticenterBondId(0))
                    .atom_ids()
                    .collect::<Vec<_>>(),
                vec![AtomId(1), AtomId(0)]
            );
            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[rstest]
fn test_transaction_apply_replace_noncovalent_bond_atoms(batched_overlays: MoleculeEditor) {
    let before = batched_overlays.clone().build();
    let mut molecule = batched_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::ReplaceNoncovalentBondAtoms {
                    id: NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                    old: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    new: [AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(0))],
                }]))
                .unwrap();
            assert_eq!(
                transaction.journal.as_slice(),
                &[Undo::RestoreNoncovalentBondAtoms {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                }]
            );
            assert_eq!(
                transaction
                    .molecule
                    .noncovalent_bond(NoncovalentBondId(0))
                    .atom_ids(),
                [AtomId(1), AtomId(0)]
            );
            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[rstest]
fn test_transaction_apply_replace_stereo_atom_site(batched_overlays: MoleculeEditor) {
    let before = batched_overlays.clone().build();
    let mut molecule = batched_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::ReplaceStereoAtomSite {
                    id: StereoAtomHandle::Id(StereoAtomId(0)),
                    old: AtomHandle::Id(AtomId(0)),
                    new: AtomHandle::Id(AtomId(2)),
                }]))
                .unwrap();
            assert_eq!(
                transaction.journal.as_slice(),
                &[Undo::RestoreStereoAtomSite {
                    id: StereoAtomId(0),
                    site: AtomId(0),
                }]
            );
            assert_eq!(
                transaction.molecule.stereo_atom(StereoAtomId(0)).site_id(),
                AtomId(2)
            );
            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[rstest]
fn test_transaction_apply_replace_stereo_atom_ligands(batched_overlays: MoleculeEditor) {
    let before = batched_overlays.clone().build();
    let old = vec![
        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
    ];
    let new = vec![old[2], old[0], old[1]];
    let mut molecule = batched_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::ReplaceStereoAtomLigands {
                    id: StereoAtomHandle::Id(StereoAtomId(0)),
                    old: old
                        .iter()
                        .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                        .collect(),
                    new: new
                        .iter()
                        .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                        .collect(),
                }]))
                .unwrap();
            assert_eq!(
                transaction.journal.as_slice(),
                &[Undo::RestoreStereoAtomLigands {
                    id: StereoAtomId(0),
                    ligands: old,
                }]
            );
            assert_eq!(
                transaction
                    .molecule
                    .stereo_atom(StereoAtomId(0))
                    .ligand_ids(),
                new
            );
            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[rstest]
fn test_transaction_apply_replace_stereo_bond_site(batched_overlays: MoleculeEditor) {
    let before = batched_overlays.clone().build();
    let mut molecule = batched_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::ReplaceStereoBondSite {
                    id: StereoBondHandle::Id(StereoBondId(0)),
                    old: BondHandle::Id(BondId(0)),
                    new: BondHandle::Id(BondId(1)),
                }]))
                .unwrap();
            assert_eq!(
                transaction.journal.as_slice(),
                &[Undo::RestoreStereoBondSite {
                    id: StereoBondId(0),
                    site: BondId(0),
                }]
            );
            assert_eq!(
                transaction.molecule.stereo_bond(StereoBondId(0)).site_id(),
                BondId(1)
            );
            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[rstest]
fn test_transaction_apply_replace_stereo_bond_ligands(batched_overlays: MoleculeEditor) {
    let before = batched_overlays.clone().build();
    let old = vec![
        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
        StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
    ];
    let new = vec![old[1], old[0], old[3], old[2]];
    let mut molecule = batched_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::ReplaceStereoBondLigands {
                    id: StereoBondHandle::Id(StereoBondId(0)),
                    old: old
                        .iter()
                        .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                        .collect(),
                    new: new
                        .iter()
                        .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                        .collect(),
                }]))
                .unwrap();
            assert_eq!(
                transaction.journal.as_slice(),
                &[Undo::RestoreStereoBondLigands {
                    id: StereoBondId(0),
                    ligands: old,
                }]
            );
            assert_eq!(
                transaction
                    .molecule
                    .stereo_bond(StereoBondId(0))
                    .ligand_ids(),
                new
            );
            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[rstest]
#[case::dative_donors(Edit::ReplaceDativeBondDonors {
        id: DativeBondHandle::Id(DativeBondId(0)),
        old: vec![],
        new: vec![AtomHandle::Id(AtomId(2))],
    })]
#[case::dative_acceptor(Edit::ReplaceDativeBondAcceptor {
        id: DativeBondHandle::Id(DativeBondId(0)),
        old: AtomHandle::Id(AtomId(2)),
        new: AtomHandle::Id(AtomId(3)),
    })]
#[case::aromatic_atoms(Edit::ReplaceAromaticSystemAtoms {
        id: AromaticSystemHandle::Id(AromaticSystemId(0)),
        old: vec![AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(0))],
        new: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(2))],
    })]
#[case::multicenter_atoms(Edit::ReplaceMulticenterBondAtoms {
        id: MulticenterBondHandle::Id(MulticenterBondId(0)),
        old: vec![AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(0))],
        new: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(2))],
    })]
#[case::noncovalent_atoms(Edit::ReplaceNoncovalentBondAtoms {
        id: NoncovalentBondHandle::Id(NoncovalentBondId(0)),
        old: [AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(0))],
        new: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(2))],
    })]
#[case::stereo_atom_site(Edit::ReplaceStereoAtomSite {
        id: StereoAtomHandle::Id(StereoAtomId(0)),
        old: AtomHandle::Id(AtomId(2)),
        new: AtomHandle::Id(AtomId(3)),
    })]
#[case::stereo_atom_ligands(Edit::ReplaceStereoAtomLigands {
        id: StereoAtomHandle::Id(StereoAtomId(0)),
        old: vec![],
        new: vec![(AtomHandle::Id(AtomId(2)), StereoLigandKind::Atom)],
    })]
#[case::stereo_bond_site(Edit::ReplaceStereoBondSite {
        id: StereoBondHandle::Id(StereoBondId(0)),
        old: BondHandle::Id(BondId(1)),
        new: BondHandle::Id(BondId(2)),
    })]
#[case::stereo_bond_ligands(Edit::ReplaceStereoBondLigands {
        id: StereoBondHandle::Id(StereoBondId(0)),
        old: vec![],
        new: vec![(AtomHandle::Id(AtomId(2)), StereoLigandKind::Atom)],
    })]
fn test_transaction_apply_replace_old_state_error(
    batched_overlays: MoleculeEditor,
    #[case] edit: Edit,
) {
    let before = batched_overlays.clone().build();
    let mut molecule = batched_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            assert_eq!(
                transaction.apply(Edits::from_iter([edit])),
                Err(MoleculeApplyError::Transaction(
                    TransactionError::OldStateMismatch
                )),
            );
            assert_eq!(transaction.probe().unwrap(), &before);

            Ok(())
        },
    );
    assert_eq!(result, Err(TransactionError::Aborted.into()));
}

#[rstest]
#[case::dative_first(EntityKind::DativeBond, 0)]
#[case::dative_middle(EntityKind::DativeBond, 1)]
#[case::dative_last(EntityKind::DativeBond, 2)]
#[case::aromatic_first(EntityKind::AromaticSystem, 0)]
#[case::aromatic_middle(EntityKind::AromaticSystem, 1)]
#[case::aromatic_last(EntityKind::AromaticSystem, 2)]
#[case::multicenter_first(EntityKind::MulticenterBond, 0)]
#[case::multicenter_middle(EntityKind::MulticenterBond, 1)]
#[case::multicenter_last(EntityKind::MulticenterBond, 2)]
#[case::noncovalent_first(EntityKind::NoncovalentBond, 0)]
#[case::noncovalent_middle(EntityKind::NoncovalentBond, 1)]
#[case::noncovalent_last(EntityKind::NoncovalentBond, 2)]
#[case::stereo_atom_first(EntityKind::StereoAtom, 0)]
#[case::stereo_atom_middle(EntityKind::StereoAtom, 1)]
#[case::stereo_atom_last(EntityKind::StereoAtom, 2)]
#[case::stereo_bond_first(EntityKind::StereoBond, 0)]
#[case::stereo_bond_middle(EntityKind::StereoBond, 1)]
#[case::stereo_bond_last(EntityKind::StereoBond, 2)]
fn test_transaction_apply_remove_overlays_error(
    batched_overlays: MoleculeEditor,
    #[case] kind: EntityKind,
    #[case] invalid_position: usize,
) {
    let before = batched_overlays.clone().build();
    let edit = match kind {
        EntityKind::DativeBond => Edit::RemoveDativeBonds {
            removes: (0..3_u32)
                .map(|index| {
                    (
                        DativeBondHandle::Id(DativeBondId(if index as usize == invalid_position {
                            9
                        } else {
                            index
                        })),
                        vec![AtomHandle::Id(AtomId(index * 2))],
                        AtomHandle::Id(AtomId(index * 2 + 1)),
                        DativeBondForm::from_order(1),
                    )
                })
                .collect(),
        },
        EntityKind::AromaticSystem => Edit::RemoveAromaticSystems {
            removes: (0..3_u32)
                .map(|index| {
                    (
                        AromaticSystemHandle::Id(AromaticSystemId(
                            if index as usize == invalid_position {
                                9
                            } else {
                                index
                            },
                        )),
                        vec![
                            AtomHandle::Id(AtomId(index * 2)),
                            AtomHandle::Id(AtomId(index * 2 + 1)),
                        ],
                        AromaticSystemForm::default(),
                    )
                })
                .collect(),
        },
        EntityKind::MulticenterBond => Edit::RemoveMulticenterBonds {
            removes: (0..3_u32)
                .map(|index| {
                    (
                        MulticenterBondHandle::Id(MulticenterBondId(
                            if index as usize == invalid_position {
                                9
                            } else {
                                index
                            },
                        )),
                        vec![
                            AtomHandle::Id(AtomId(index * 2)),
                            AtomHandle::Id(AtomId(index * 2 + 1)),
                        ],
                        MulticenterBondForm::default(),
                    )
                })
                .collect(),
        },
        EntityKind::NoncovalentBond => Edit::RemoveNoncovalentBonds {
            removes: (0..3_u32)
                .map(|index| {
                    (
                        NoncovalentBondHandle::Id(NoncovalentBondId(
                            if index as usize == invalid_position {
                                9
                            } else {
                                index
                            },
                        )),
                        [
                            AtomHandle::Id(AtomId(index * 2)),
                            AtomHandle::Id(AtomId(index * 2 + 1)),
                        ],
                        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                    )
                })
                .collect(),
        },
        EntityKind::StereoAtom => Edit::RemoveStereoAtoms {
            removes: (0..3_u32)
                .map(|index| {
                    (
                        StereoAtomHandle::Id(StereoAtomId(if index as usize == invalid_position {
                            9
                        } else {
                            index
                        })),
                        AtomHandle::Id(AtomId(index * 2)),
                        vec![
                            (
                                AtomHandle::Id(AtomId(index * 2 + 1)),
                                StereoLigandKind::Atom,
                            ),
                            (
                                AtomHandle::Id(AtomId(index * 2)),
                                StereoLigandKind::ImplicitHydrogen,
                            ),
                            (
                                AtomHandle::Id(AtomId(index * 2)),
                                StereoLigandKind::LonePair,
                            ),
                        ],
                        StereoAtomForm::default(),
                    )
                })
                .collect(),
        },
        EntityKind::StereoBond => Edit::RemoveStereoBonds {
            removes: (0..3_u32)
                .map(|index| {
                    (
                        StereoBondHandle::Id(StereoBondId(if index as usize == invalid_position {
                            9
                        } else {
                            index
                        })),
                        BondHandle::Id(BondId(index)),
                        vec![
                            (
                                AtomHandle::Id(AtomId(index * 2)),
                                StereoLigandKind::ImplicitHydrogen,
                            ),
                            (
                                AtomHandle::Id(AtomId(index * 2)),
                                StereoLigandKind::LonePair,
                            ),
                            (
                                AtomHandle::Id(AtomId(index * 2 + 1)),
                                StereoLigandKind::ImplicitHydrogen,
                            ),
                            (
                                AtomHandle::Id(AtomId(index * 2 + 1)),
                                StereoLigandKind::LonePair,
                            ),
                        ],
                        StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                    )
                })
                .collect(),
        },
        EntityKind::Atom | EntityKind::Bond => unreachable!(),
    };
    let mut molecule = batched_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            assert_eq!(
                transaction.apply(Edits::from_iter([edit])).unwrap_err(),
                MoleculeApplyError::Transaction(TransactionError::HandleOutOfRange {
                    kind,
                    index: 9,
                    count: 3,
                })
            );
            assert_eq!(transaction.probe().unwrap(), &before);

            Ok(())
        },
    );
    assert_eq!(result, Err(TransactionError::Aborted.into()));
}

#[rstest]
#[case::atom(EntityKind::Atom)]
#[case::bond(EntityKind::Bond)]
#[case::dative_bond(EntityKind::DativeBond)]
#[case::aromatic_system(EntityKind::AromaticSystem)]
#[case::multicenter_bond(EntityKind::MulticenterBond)]
#[case::noncovalent_bond(EntityKind::NoncovalentBond)]
#[case::stereo_atom(EntityKind::StereoAtom)]
#[case::stereo_bond(EntityKind::StereoBond)]
fn test_transaction_apply_duplicate_removal_error(
    batched_overlays: MoleculeEditor,
    #[case] kind: EntityKind,
) {
    let before = batched_overlays.clone().build();
    let edit = match kind {
        EntityKind::Atom => Edit::RemoveTopology {
            atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(0))],
            bonds: Vec::new(),
        },
        EntityKind::Bond => Edit::RemoveTopology {
            atoms: Vec::new(),
            bonds: vec![BondHandle::Id(BondId(0)), BondHandle::Id(BondId(0))],
        },
        EntityKind::DativeBond => Edit::RemoveDativeBonds {
            removes: vec![
                (
                    DativeBondHandle::Id(DativeBondId(0)),
                    vec![AtomHandle::Id(AtomId(0))],
                    AtomHandle::Id(AtomId(1)),
                    DativeBondForm::from_order(1),
                ),
                (
                    DativeBondHandle::Id(DativeBondId(0)),
                    vec![AtomHandle::Id(AtomId(0))],
                    AtomHandle::Id(AtomId(1)),
                    DativeBondForm::from_order(1),
                ),
            ],
        },
        EntityKind::AromaticSystem => Edit::RemoveAromaticSystems {
            removes: vec![
                (
                    AromaticSystemHandle::Id(AromaticSystemId(0)),
                    vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    AromaticSystemForm::default(),
                ),
                (
                    AromaticSystemHandle::Id(AromaticSystemId(0)),
                    vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    AromaticSystemForm::default(),
                ),
            ],
        },
        EntityKind::MulticenterBond => Edit::RemoveMulticenterBonds {
            removes: vec![
                (
                    MulticenterBondHandle::Id(MulticenterBondId(0)),
                    vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    MulticenterBondForm::default(),
                ),
                (
                    MulticenterBondHandle::Id(MulticenterBondId(0)),
                    vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    MulticenterBondForm::default(),
                ),
            ],
        },
        EntityKind::NoncovalentBond => Edit::RemoveNoncovalentBonds {
            removes: vec![
                (
                    NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                    [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                ),
                (
                    NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                    [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                ),
            ],
        },
        EntityKind::StereoAtom => Edit::RemoveStereoAtoms {
            removes: vec![
                (
                    StereoAtomHandle::Id(StereoAtomId(0)),
                    AtomHandle::Id(AtomId(0)),
                    vec![
                        (AtomHandle::Id(AtomId(1)), StereoLigandKind::Atom),
                        (
                            AtomHandle::Id(AtomId(0)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                        (AtomHandle::Id(AtomId(0)), StereoLigandKind::LonePair),
                    ],
                    StereoAtomForm::default(),
                ),
                (
                    StereoAtomHandle::Id(StereoAtomId(0)),
                    AtomHandle::Id(AtomId(0)),
                    vec![
                        (AtomHandle::Id(AtomId(1)), StereoLigandKind::Atom),
                        (
                            AtomHandle::Id(AtomId(0)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                        (AtomHandle::Id(AtomId(0)), StereoLigandKind::LonePair),
                    ],
                    StereoAtomForm::default(),
                ),
            ],
        },
        EntityKind::StereoBond => Edit::RemoveStereoBonds {
            removes: vec![
                (
                    StereoBondHandle::Id(StereoBondId(0)),
                    BondHandle::Id(BondId(0)),
                    vec![
                        (
                            AtomHandle::Id(AtomId(0)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                        (AtomHandle::Id(AtomId(0)), StereoLigandKind::LonePair),
                        (
                            AtomHandle::Id(AtomId(1)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                        (AtomHandle::Id(AtomId(1)), StereoLigandKind::LonePair),
                    ],
                    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                ),
                (
                    StereoBondHandle::Id(StereoBondId(0)),
                    BondHandle::Id(BondId(0)),
                    vec![
                        (
                            AtomHandle::Id(AtomId(0)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                        (AtomHandle::Id(AtomId(0)), StereoLigandKind::LonePair),
                        (
                            AtomHandle::Id(AtomId(1)),
                            StereoLigandKind::ImplicitHydrogen,
                        ),
                        (AtomHandle::Id(AtomId(1)), StereoLigandKind::LonePair),
                    ],
                    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                ),
            ],
        },
    };
    let mut molecule = batched_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            assert_eq!(
                transaction.apply(Edits::from_iter([edit])).unwrap_err(),
                MoleculeApplyError::Transaction(TransactionError::DuplicateRemoval { kind })
            );
            assert_eq!(transaction.probe().unwrap(), &before);

            Ok(())
        },
    );
    assert_eq!(result, Err(TransactionError::Aborted.into()));
}

#[rstest]
#[case::dative_bond(EntityKind::DativeBond)]
#[case::aromatic_system(EntityKind::AromaticSystem)]
#[case::multicenter_bond(EntityKind::MulticenterBond)]
#[case::noncovalent_bond(EntityKind::NoncovalentBond)]
#[case::stereo_atom(EntityKind::StereoAtom)]
#[case::stereo_bond(EntityKind::StereoBond)]
fn test_transaction_apply_handle_removed_error_cascade(
    batched_overlays: MoleculeEditor,
    #[case] kind: EntityKind,
) {
    let before = batched_overlays.clone().build();
    let mut edits = Edits::from_iter([Edit::RemoveTopology {
        atoms: vec![AtomHandle::Id(AtomId(0))],
        bonds: Vec::new(),
    }]);
    edits.push(match kind {
        EntityKind::DativeBond => Edit::RemoveDativeBonds {
            removes: vec![(
                DativeBondHandle::Id(DativeBondId(0)),
                vec![AtomHandle::Id(AtomId(0))],
                AtomHandle::Id(AtomId(1)),
                DativeBondForm::from_order(1),
            )],
        },
        EntityKind::AromaticSystem => Edit::RemoveAromaticSystems {
            removes: vec![(
                AromaticSystemHandle::Id(AromaticSystemId(0)),
                vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                AromaticSystemForm::default(),
            )],
        },
        EntityKind::MulticenterBond => Edit::RemoveMulticenterBonds {
            removes: vec![(
                MulticenterBondHandle::Id(MulticenterBondId(0)),
                vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                MulticenterBondForm::default(),
            )],
        },
        EntityKind::NoncovalentBond => Edit::RemoveNoncovalentBonds {
            removes: vec![(
                NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
        },
        EntityKind::StereoAtom => Edit::RemoveStereoAtoms {
            removes: vec![(
                StereoAtomHandle::Id(StereoAtomId(0)),
                AtomHandle::Id(AtomId(0)),
                vec![(AtomHandle::Id(AtomId(1)), StereoLigandKind::Atom)],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
        },
        EntityKind::StereoBond => Edit::RemoveStereoBonds {
            removes: vec![(
                StereoBondHandle::Id(StereoBondId(0)),
                BondHandle::Id(BondId(0)),
                vec![
                    (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                    (AtomHandle::Id(AtomId(1)), StereoLigandKind::Atom),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
        },
        EntityKind::Atom | EntityKind::Bond => unreachable!(),
    });
    let mut molecule = batched_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            assert_eq!(
                transaction.apply(edits).unwrap_err(),
                MoleculeApplyError::Transaction(TransactionError::HandleRemoved { kind, index: 0 })
            );
            assert_eq!(transaction.probe().unwrap(), &before);

            Ok(())
        },
    );
    assert_eq!(result, Err(TransactionError::Aborted.into()));
}

#[rstest]
fn test_transaction_apply_set_dative_bond_field(diatomic_with_overlays: MoleculeEditor) {
    let mut molecule = diatomic_with_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::ModifyDativeBondField {
                    id: DativeBondHandle::Id(DativeBondId(0)),
                    change: DativeBondFieldChange::Order {
                        old: NumForm::Lit(1),
                        new: NumForm::Lit(2),
                    },
                }]))
                .unwrap();
            assert_eq!(
                transaction
                    .molecule
                    .dative_bond(DativeBondId(0))
                    .attributes()
                    .order,
                NumForm::Lit(2),
            );

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_set_aromatic_system_field(diatomic_with_overlays: MoleculeEditor) {
    let mut molecule = diatomic_with_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::ModifyAromaticSystemField {
                    id: AromaticSystemHandle::Id(AromaticSystemId(0)),
                    change: AromaticSystemFieldChange::Charge {
                        old: NumForm::default(),
                        new: NumForm::Lit(1),
                    },
                }]))
                .unwrap();
            assert_eq!(
                transaction
                    .molecule
                    .aromatic_system(AromaticSystemId(0))
                    .attributes()
                    .charge,
                NumForm::Lit(1),
            );

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_set_multicenter_bond_field(diatomic_with_overlays: MoleculeEditor) {
    let mut molecule = diatomic_with_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::ModifyMulticenterBondField {
                    id: MulticenterBondHandle::Id(MulticenterBondId(0)),
                    change: MulticenterBondFieldChange::Charge {
                        old: NumForm::default(),
                        new: NumForm::Lit(-1),
                    },
                }]))
                .unwrap();
            assert_eq!(
                transaction
                    .molecule
                    .multicenter_bond(MulticenterBondId(0))
                    .attributes()
                    .charge,
                NumForm::Lit(-1),
            );

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_set_noncovalent_bond_field(diatomic_with_overlays: MoleculeEditor) {
    let mut molecule = diatomic_with_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::ModifyNoncovalentBondField {
                    id: NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                    change: NoncovalentBondFieldChange::Kind {
                        old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
                        new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
                    },
                }]))
                .unwrap();
            assert_eq!(
                transaction
                    .molecule
                    .noncovalent_bond(NoncovalentBondId(0))
                    .attributes()
                    .kind,
                NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
            );

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_add_dative_bond(diatomic: MoleculeEditor) {
    let mut molecule = diatomic.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::AddDativeBond {
                    donors: vec![AtomHandle::Id(AtomId(0))],
                    acceptor: AtomHandle::Id(AtomId(1)),
                    attributes: DativeBondForm::from_order(1),
                }]))
                .unwrap();
            assert!(matches!(
                transaction.journal.as_slice(),
                [Undo::RemoveAddedDativeBond(added)] if added.id == DativeBondId(0)
            ));
            assert_eq!(transaction.molecule.dative_bonds().count(), 1);

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_add_aromatic_system(diatomic: MoleculeEditor) {
    let mut molecule = diatomic.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::AddAromaticSystem {
                    atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    attributes: AromaticSystemForm::default(),
                }]))
                .unwrap();
            assert!(matches!(
                transaction.journal.as_slice(),
                [Undo::RemoveAddedAromaticSystem(added)] if added.id == AromaticSystemId(0)
            ));
            assert_eq!(transaction.molecule.aromatic_systems().count(), 1);

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_add_multicenter_bond(diatomic: MoleculeEditor) {
    let mut molecule = diatomic.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::AddMulticenterBond {
                    atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    attributes: MulticenterBondForm::default(),
                }]))
                .unwrap();
            assert!(matches!(
                transaction.journal.as_slice(),
                [Undo::RemoveAddedMulticenterBond(added)] if added.id == MulticenterBondId(0)
            ));
            assert_eq!(transaction.molecule.multicenter_bonds().count(), 1);

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_add_noncovalent_bond(diatomic: MoleculeEditor) {
    let mut molecule = diatomic.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::AddNoncovalentBond {
                    atoms: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                }]))
                .unwrap();
            assert!(matches!(
                transaction.journal.as_slice(),
                [Undo::RemoveAddedNoncovalentBond(added)] if added.id == NoncovalentBondId(0)
            ));
            assert_eq!(transaction.molecule.noncovalent_bonds().count(), 1);

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_remove_dative_bond(diatomic_with_overlays: MoleculeEditor) {
    let before = diatomic_with_overlays.clone().build();
    let mut molecule = diatomic_with_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::RemoveDativeBonds {
                    removes: vec![(
                        DativeBondHandle::Id(DativeBondId(0)),
                        vec![AtomHandle::Id(AtomId(0))],
                        AtomHandle::Id(AtomId(1)),
                        DativeBondForm {
                            order: NumForm::Lit(1),
                            constraints: Default::default(),
                        },
                    )],
                }]))
                .unwrap();
            assert_eq!(transaction.molecule.dative_bonds().count(), 0);
            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[rstest]
fn test_transaction_apply_remove_dative_bond_roles_error(diatomic_with_overlays: MoleculeEditor) {
    let mut molecule = diatomic_with_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            let err = transaction
                .apply(Edits::from_iter([Edit::RemoveDativeBonds {
                    removes: vec![(
                        DativeBondHandle::Id(DativeBondId(0)),
                        vec![AtomHandle::Id(AtomId(1))],
                        AtomHandle::Id(AtomId(0)),
                        DativeBondForm {
                            order: NumForm::Lit(1),
                            constraints: Default::default(),
                        },
                    )],
                }]))
                .unwrap_err();
            assert_eq!(
                err,
                MoleculeApplyError::Transaction(TransactionError::OldStateMismatch)
            );
            assert_eq!(transaction.molecule.dative_bonds().count(), 1);

            Ok(())
        },
    );
    assert_eq!(result, Err(TransactionError::Aborted.into()));
}

#[rstest]
fn test_transaction_apply_remove_aromatic_system(diatomic_with_overlays: MoleculeEditor) {
    let before = diatomic_with_overlays.clone().build();
    let mut molecule = diatomic_with_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::RemoveAromaticSystems {
                    removes: vec![(
                        AromaticSystemHandle::Id(AromaticSystemId(0)),
                        vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        AromaticSystemForm::default(),
                    )],
                }]))
                .unwrap();
            assert_eq!(transaction.molecule.aromatic_systems().count(), 0);
            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

// Batch removal of non-contiguous same-kind ids (0 and 2) in one edit: ids resolve against the
// pre-removal state and compact once, so the survivor (former id 1) remaps to id 0. A single-id
// sequence would stale id 2 after removing id 0.
#[rstest]
fn test_transaction_apply_remove_aromatic_systems() {
    let mut b = Molecule::default().edit();
    for _ in 0..6 {
        b.add_atom(AtomForm::from_element(Element::C));
    }
    b.add_aromatic_system(&[AtomId(0), AtomId(1)], AromaticSystemForm::default());
    b.add_aromatic_system(&[AtomId(2), AtomId(3)], AromaticSystemForm::default());
    b.add_aromatic_system(&[AtomId(4), AtomId(5)], AromaticSystemForm::default());
    let mut molecule = b.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::RemoveAromaticSystems {
                    removes: vec![
                        (
                            AromaticSystemHandle::Id(AromaticSystemId(0)),
                            vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                            AromaticSystemForm::default(),
                        ),
                        (
                            AromaticSystemHandle::Id(AromaticSystemId(2)),
                            vec![AtomHandle::Id(AtomId(4)), AtomHandle::Id(AtomId(5))],
                            AromaticSystemForm::default(),
                        ),
                    ],
                }]))
                .unwrap();
            assert_eq!(transaction.molecule.aromatic_systems().count(), 1);
            assert_eq!(
                transaction
                    .molecule
                    .aromatic_system(AromaticSystemId(0))
                    .atom_ids()
                    .collect::<Vec<_>>(),
                vec![AtomId(2), AtomId(3)],
            );

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

// Rolling back an aromatic-system removal restores a molecule constraint the removal dropped
// (`dropped`) or remapped (`remapped`) — the overlay-remove undo captures the constraint cascade.
#[rstest]
#[case::dropped(AromaticSystemId(0), 0)]
#[case::remapped(AromaticSystemId(1), 1)]
fn test_transaction_apply_remove_aromatic_system_rollback(
    #[case] constrained: AromaticSystemId,
    #[case] forward_constraint_count: usize,
) {
    let mut b = Molecule::default().edit();
    for _ in 0..6 {
        b.add_atom(AtomForm::from_element(Element::C));
    }
    b.add_aromatic_system(
        &[AtomId(0), AtomId(1), AtomId(2)],
        AromaticSystemForm::default(),
    );
    b.add_aromatic_system(
        &[AtomId(3), AtomId(4), AtomId(5)],
        AromaticSystemForm::default(),
    );
    let mut molecule = b.build();
    molecule
        .transact([Edits::from_iter([Edit::AddMoleculeConstraint {
            constraint: Constraint::AromaticSystem(
                constrained,
                AromaticSystemConstraintForm::ElectronCount(NumForm::Lit(6)),
            )
            .into(),
        }])])
        .unwrap();
    let before = molecule.clone();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::RemoveAromaticSystems {
                    removes: vec![(
                        AromaticSystemHandle::Id(AromaticSystemId(0)),
                        vec![
                            AtomHandle::Id(AtomId(0)),
                            AtomHandle::Id(AtomId(1)),
                            AtomHandle::Id(AtomId(2)),
                        ],
                        AromaticSystemForm::default(),
                    )],
                }]))
                .unwrap();
            assert_eq!(transaction.molecule.aromatic_systems().count(), 1);
            assert_eq!(
                transaction.molecule.constraints().iter().count(),
                forward_constraint_count
            );

            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[rstest]
fn test_transaction_apply_remove_multicenter_bond(diatomic_with_overlays: MoleculeEditor) {
    let before = diatomic_with_overlays.clone().build();
    let mut molecule = diatomic_with_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::RemoveMulticenterBonds {
                    removes: vec![(
                        MulticenterBondHandle::Id(MulticenterBondId(0)),
                        vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        MulticenterBondForm::default(),
                    )],
                }]))
                .unwrap();
            assert_eq!(transaction.molecule.multicenter_bonds().count(), 0);
            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[rstest]
fn test_transaction_apply_remove_noncovalent_bond(diatomic_with_overlays: MoleculeEditor) {
    let before = diatomic_with_overlays.clone().build();
    let mut molecule = diatomic_with_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::RemoveNoncovalentBonds {
                    removes: vec![(
                        NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                        [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                    )],
                }]))
                .unwrap();
            assert_eq!(transaction.molecule.noncovalent_bonds().count(), 0);
            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[rstest]
fn test_transaction_apply_remove_noncovalent_bond_form_mismatch_error(
    diatomic_with_overlays: MoleculeEditor,
) {
    let mut molecule = diatomic_with_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            let err = transaction
                .apply(Edits::from_iter([Edit::RemoveNoncovalentBonds {
                    removes: vec![(
                        NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                        [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        NoncovalentBondForm::from_kind(NoncovalentBondKind::Ionic), // wrong
                    )],
                }]))
                .unwrap_err();
            assert_eq!(
                err,
                MoleculeApplyError::Transaction(TransactionError::OldStateMismatch)
            );

            Ok(())
        },
    );
    assert_eq!(result, Err(TransactionError::Aborted.into()));
}

#[rstest]
fn test_transaction_apply_set_bond_constraint_value_bearing(diatomic: MoleculeEditor) {
    let mut molecule = diatomic.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(0)),
                    old: None,
                    new: Some(BondConstraintForm::cis_trans_stereo(
                        CisTransStereoForm::NotStereo,
                    )),
                }]))
                .unwrap();
            assert_eq!(
                transaction
                    .molecule
                    .bond(BondId(0))
                    .attributes()
                    .constraints
                    .iter()
                    .cloned()
                    .collect::<Vec<_>>(),
                vec![BondConstraintForm::cis_trans_stereo(
                    CisTransStereoForm::NotStereo
                )],
            );

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_add_bond_constraint(diatomic: MoleculeEditor) {
    let mut molecule = diatomic.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(0)),
                    old: None,
                    new: Some(BondConstraintForm::ring_membership(RingScope::Size(5), 1)),
                }]))
                .unwrap();
            assert!(transaction
                .molecule
                .bond(BondId(0))
                .attributes()
                .constraints
                .iter()
                .any(|c| *c == BondConstraintForm::ring_membership(RingScope::Size(5), 1)));

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_modify_bond_constraint_absent_error(diatomic: MoleculeEditor) {
    let mut molecule = diatomic.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            let err = transaction
                .apply(Edits::from_iter([Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(0)),
                    old: Some(BondConstraintForm::ring_membership(RingScope::Size(5), 1)),
                    new: None,
                }]))
                .unwrap_err();
            assert_eq!(
                err,
                MoleculeApplyError::Transaction(TransactionError::OldStateMismatch)
            );

            Ok(())
        },
    );
    assert_eq!(result, Err(TransactionError::Aborted.into()));
}

#[rstest]
fn test_transaction_apply_set_dative_bond_constraint(diatomic_with_overlays: MoleculeEditor) {
    let mut molecule = diatomic_with_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::ModifyDativeBondConstraint {
                    id: DativeBondHandle::Id(DativeBondId(0)),
                    old: None,
                    new: Some(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
                }]))
                .unwrap();
            assert!(transaction
                .molecule
                .dative_bond(DativeBondId(0))
                .attributes()
                .constraints
                .iter()
                .any(|c| *c == DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))));

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_set_aromatic_system_constraint(diatomic_with_overlays: MoleculeEditor) {
    let mut molecule = diatomic_with_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::ModifyAromaticSystemConstraint {
                    id: AromaticSystemHandle::Id(AromaticSystemId(0)),
                    old: None,
                    new: Some(AromaticSystemConstraintForm::ElectronCount(NumForm::Lit(6))),
                }]))
                .unwrap();
            assert_eq!(
                transaction
                    .molecule
                    .aromatic_system(AromaticSystemId(0))
                    .attributes()
                    .constraints
                    .iter()
                    .cloned()
                    .collect::<Vec<_>>(),
                vec![AromaticSystemConstraintForm::ElectronCount(NumForm::Lit(6))],
            );

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[rstest]
fn test_transaction_apply_set_multicenter_bond_constraint(diatomic_with_overlays: MoleculeEditor) {
    let mut molecule = diatomic_with_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::ModifyMulticenterBondConstraint {
                    id: MulticenterBondHandle::Id(MulticenterBondId(0)),
                    old: None,
                    new: Some(MulticenterBondConstraintForm::ElectronCount(NumForm::Lit(
                        2,
                    ))),
                }]))
                .unwrap();
            assert_eq!(
                transaction
                    .molecule
                    .multicenter_bond(MulticenterBondId(0))
                    .attributes()
                    .constraints
                    .iter()
                    .cloned()
                    .collect::<Vec<_>>(),
                vec![MulticenterBondConstraintForm::ElectronCount(NumForm::Lit(
                    2
                ))],
            );

            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
}

#[fixture]
fn triatomic_with_overlays() -> MoleculeEditor {
    let mut b = Molecule::default().edit();
    b.add_atom(AtomForm::from_element(Element::C));
    b.add_atom(AtomForm::from_element(Element::N));
    b.add_atom(AtomForm::from_element(Element::O));
    b.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));
    b.add_bond(AtomId(1), AtomId(2), BondForm::from_order(1));
    b.add_dative_bond(&[AtomId(0)], AtomId(1), DativeBondForm::from_order(1));
    b.add_aromatic_system(
        &[AtomId(0), AtomId(1), AtomId(2)],
        AromaticSystemForm::default(),
    );
    b.add_multicenter_bond(
        &[AtomId(0), AtomId(1), AtomId(2)],
        MulticenterBondForm::default(),
    );
    b.add_noncovalent_bond(
        [AtomId(0), AtomId(2)],
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
    );
    b
}

#[rstest]
fn test_transaction_rollback_batches(diatomic: MoleculeEditor) {
    let mut molecule = diatomic.build();
    let before = molecule.clone();
    Transaction::run(&mut molecule, |mut transaction| {
        transaction.apply(Edits::from_iter([Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(1),
            },
        }]))?;
        transaction.apply(Edits::from_iter([Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        }]))?;
        assert_eq!(
            transaction.probe()?.atom(AtomId(0)).attributes().charge,
            NumForm::Lit(2)
        );
        transaction.rollback();
        Ok::<_, MoleculeApplyError>(())
    })
    .unwrap();
    assert_eq!(molecule, before);
}

#[rstest]
fn test_transaction_apply_batches_error(diatomic: MoleculeEditor) {
    let mut molecule = diatomic.build();
    let before = molecule.clone();
    let result = Transaction::run(&mut molecule, |mut transaction| {
        transaction.apply(Edits::from_iter([Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(1),
            },
        }]))?;
        transaction.apply(Edits::from_iter([Edit::ModifyAtomConstraint {
            id: AtomHandle::Id(AtomId(0)),
            old: None,
            new: Some(AtomConstraintForm::degree(1)),
        }]))?;
        transaction.apply(Edits::from_iter([Edit::AddDativeBond {
            donors: vec![AtomHandle::Id(AtomId(0))],
            acceptor: AtomHandle::Id(AtomId(1)),
            attributes: DativeBondForm::from_order(1),
        }]))?;
        let mut rejected = Edits::from_iter([Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        }]);
        rejected.remove_atom(AtomHandle::Id(AtomId(99)));
        transaction.apply(rejected)
    });
    assert_eq!(
        result,
        Err(MoleculeApplyError::Transaction(
            TransactionError::HandleOutOfRange {
                kind: EntityKind::Atom,
                index: 99,
                count: 2,
            }
        ))
    );
    assert_eq!(molecule, before);
}

enum RollbackCase {
    RemoveTopology,
    RemoveBond,
    AddTopology,
    Field,
    AddOverlay,
    RemoveOverlay,
    Constraint,
    CascadedConstraints,
}

#[rstest]
#[case::remove_topology(RollbackCase::RemoveTopology)]
#[case::remove_bond(RollbackCase::RemoveBond)]
#[case::add_topology(RollbackCase::AddTopology)]
#[case::field(RollbackCase::Field)]
#[case::add_overlay(RollbackCase::AddOverlay)]
#[case::remove_overlay(RollbackCase::RemoveOverlay)]
#[case::constraint(RollbackCase::Constraint)]
#[case::cascade(RollbackCase::CascadedConstraints)]
fn test_transaction_rollback_entries(#[case] case: RollbackCase) {
    let editor = match case {
        RollbackCase::AddTopology => Molecule::default().edit(),
        RollbackCase::Field | RollbackCase::Constraint => {
            let mut b = Molecule::default().edit();
            b.add_atom(AtomForm::from_element(Element::C));
            b
        }
        RollbackCase::AddOverlay => {
            let mut b = Molecule::default().edit();
            b.add_atom(AtomForm::from_element(Element::C));
            b.add_atom(AtomForm::from_element(Element::C));
            b.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));
            b
        }
        RollbackCase::RemoveBond => {
            let mut b = Molecule::default().edit();
            b.add_atom(AtomForm::from_element(Element::C));
            b.add_atom(AtomForm::from_element(Element::N));
            b.add_atom(AtomForm::from_element(Element::O));
            b.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));
            b.add_bond(AtomId(1), AtomId(2), BondForm::from_order(2));
            b
        }
        RollbackCase::CascadedConstraints => {
            let mut b = Molecule::default().edit();
            b.add_atom(AtomForm::from_element(Element::C));
            b.add_atom(AtomForm::from_element(Element::N));
            b.constraints_mut()
                .push(Constraint::Atom(AtomId(1), AtomConstraintForm::degree(3)));
            b
        }
        RollbackCase::RemoveTopology | RollbackCase::RemoveOverlay => triatomic_with_overlays(),
    };
    let before = editor.clone().build();
    let edits = match case {
        RollbackCase::RemoveTopology => Edits::from_iter([Edit::RemoveTopology {
            atoms: vec![AtomHandle::Id(AtomId(1))],
            bonds: vec![],
        }]),
        RollbackCase::RemoveBond => Edits::from_iter([Edit::RemoveTopology {
            atoms: Vec::new(),
            bonds: vec![BondHandle::Id(BondId(0))],
        }]),
        RollbackCase::AddTopology => Edits::from_iter([
            Edit::AddAtoms {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::O),
                ],
            },
            Edit::AddBonds {
                bonds: vec![AddBond {
                    endpoints: [AtomHandle::New(0), AtomHandle::New(1)],
                    attributes: BondForm::from_order(2),
                }],
            },
        ]),
        RollbackCase::Field => Edits::from_iter([Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(1),
            },
        }]),
        RollbackCase::AddOverlay => Edits::from_iter([Edit::AddDativeBond {
            donors: vec![AtomHandle::Id(AtomId(0))],
            acceptor: AtomHandle::Id(AtomId(1)),
            attributes: DativeBondForm::from_order(1),
        }]),
        RollbackCase::RemoveOverlay => Edits::from_iter([Edit::RemoveDativeBonds {
            removes: vec![(
                DativeBondHandle::Id(DativeBondId(0)),
                vec![AtomHandle::Id(AtomId(0))],
                AtomHandle::Id(AtomId(1)),
                DativeBondForm {
                    order: NumForm::Lit(1),
                    constraints: Default::default(),
                },
            )],
        }]),
        RollbackCase::Constraint => Edits::from_iter([Edit::ModifyAtomConstraint {
            id: AtomHandle::Id(AtomId(0)),
            old: None,
            new: Some(AtomConstraintForm::ring_membership(RingScope::Size(5), 1)),
        }]),
        RollbackCase::CascadedConstraints => Edits::from_iter([Edit::RemoveTopology {
            atoms: vec![AtomHandle::Id(AtomId(1))],
            bonds: Vec::new(),
        }]),
    };
    let mut molecule = editor.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction.apply(edits).unwrap();
            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}

#[rstest]
#[case::first(vec![AtomHandle::Id(AtomId(0))], 2)]
#[case::middle(vec![AtomHandle::Id(AtomId(2))], 2)]
#[case::separated(vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(4))], 1)]
fn test_transaction_rollback_topology(
    mut batched_overlays: MoleculeEditor,
    #[case] atoms: Vec<AtomHandle>,
    #[case] remaining: usize,
) {
    batched_overlays.constraints_mut().extend(vec![
        Constraint::Atom(AtomId(0), AtomConstraintForm::degree(1)),
        Constraint::Atom(AtomId(5), AtomConstraintForm::degree(1)),
        Constraint::Atom(AtomId(0), AtomConstraintForm::degree(1)),
    ]);
    let before = batched_overlays.clone().build();
    let mut molecule = batched_overlays.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::RemoveTopology {
                    atoms,
                    bonds: vec![],
                }]))
                .unwrap();
            assert_eq!(
                [
                    transaction.molecule.dative_bonds().count(),
                    transaction.molecule.aromatic_systems().count(),
                    transaction.molecule.multicenter_bonds().count(),
                    transaction.molecule.noncovalent_bonds().count(),
                    transaction.molecule.stereo_atoms().count(),
                    transaction.molecule.stereo_bonds().count(),
                ],
                [remaining; 6],
            );
            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    let restored = molecule;
    assert!(restored.normalized_eq(&before));
    assert_eq!(restored, before);
}

#[rstest]
fn test_transaction_rollback_empty(one_atom: MoleculeEditor) {
    let mut molecule = one_atom.build();
    let before = molecule.clone();
    Transaction::run(&mut molecule, |transaction| {
        transaction.rollback();
        Ok::<_, MoleculeApplyError>(())
    })
    .unwrap();
    assert_eq!(molecule, before);
}

#[rstest]
fn test_transaction_rollback_molecule_constraint_order(mut one_atom: MoleculeEditor) {
    let repeated = Constraint::Atom(AtomId(0), AtomConstraintForm::degree(1));
    let middle = Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4));
    one_atom.constraints_mut().push(repeated.clone());
    one_atom.constraints_mut().push(middle.clone());
    one_atom.constraints_mut().push(repeated.clone());
    let before = one_atom.clone().build();
    let mut molecule = one_atom.build();
    let result = Transaction::run(
        &mut molecule,
        |mut transaction| -> Result<(), MoleculeApplyError> {
            transaction
                .apply(Edits::from_iter([Edit::RemoveMoleculeConstraint {
                    constraint: repeated.clone().into(),
                }]))
                .unwrap();
            assert_eq!(
                transaction.molecule.constraints().as_slice(),
                &[repeated, middle]
            );

            transaction.rollback();
            Ok(())
        },
    );
    assert_eq!(result, Ok(()));
    assert_eq!(molecule, before);
}
