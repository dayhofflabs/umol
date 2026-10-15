use rstest::*;
use umol_chem::element::Element;

use super::super::ApplicationState;
use super::removal_entries;
use crate::ir::aromatic::AromaticSystemForm;
use crate::ir::atom::{AtomForm, ElementForm, IsotopeMassForm};
use crate::ir::bond::BondForm;
use crate::ir::dative::DativeBondForm;
use crate::ir::edit::{
    AromaticSystemFieldChange, AromaticSystemHandle, AtomFieldChange, AtomHandle, BondFieldChange,
    BondHandle, DativeBondFieldChange, DativeBondHandle, Edit, Edits, MulticenterBondFieldChange,
    MulticenterBondHandle, NoncovalentBondFieldChange, NoncovalentBondHandle,
    StereoAtomFieldChange, StereoAtomHandle, StereoBondFieldChange, StereoBondHandle, Undo,
};
use crate::ir::electrons::ElectronCountsForm;
use crate::ir::entity::EntityKind;
use crate::ir::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ir::molecule::MoleculeEntries;
use crate::ir::multicenter::MulticenterBondForm;
use crate::ir::noncovalent::{NoncovalentBondForm, NoncovalentBondKind, NoncovalentBondKindForm};
use crate::ir::num::NumForm;
use crate::ir::spin::UnpairedElectronsForm;
use crate::ir::stereo::{
    StereoAtomForm, StereoBondForm, StereoConfigurationForm, StereoCoset, StereoKind,
};
use crate::ir::traits::Normalize;
use crate::ir::{Molecule, MoleculeApplyError, RollbackMethod, Transaction, TransactionError};

#[rstest]
#[case::element(
    AtomFieldChange::Element { old: ElementForm::Lit(Element::C), new: ElementForm::Lit(Element::N) },
    AtomForm::from_element(Element::N),
)]
#[case::element_equivalent(
    AtomFieldChange::Element { old: ElementForm::lit_set([Element::C]), new: ElementForm::Lit(Element::N) },
    AtomForm::from_element(Element::N),
)]
#[case::isotope_mass(
    AtomFieldChange::IsotopeMass { old: IsotopeMassForm::Undetermined, new: IsotopeMassForm::Lit(13) },
    AtomForm { isotope_mass: IsotopeMassForm::Lit(13), ..AtomForm::from_element(Element::C) },
)]
#[case::charge(
    AtomFieldChange::Charge { old: NumForm::Undetermined, new: NumForm::Lit(-1) },
    AtomForm { charge: NumForm::Lit(-1), ..AtomForm::from_element(Element::C) },
)]
#[case::implicit_hydrogens(
    AtomFieldChange::ImplicitHydrogens { old: NumForm::Undetermined, new: NumForm::Lit(2) },
    AtomForm { implicit_hydrogens: NumForm::Lit(2), ..AtomForm::from_element(Element::C) },
)]
#[case::lone_pairs(
    AtomFieldChange::LonePairs { old: NumForm::Undetermined, new: NumForm::Lit(1) },
    AtomForm { lone_pairs: NumForm::Lit(1), ..AtomForm::from_element(Element::C) },
)]
#[case::unpaired_electrons(
    AtomFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::default(), new: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) } },
    AtomForm { unpaired_electrons: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) }, ..AtomForm::from_element(Element::C) },
)]
fn test_molecule_apply_edit_modify_atom_field(
    mut removal_entries: MoleculeEntries,
    #[case] change: AtomFieldChange,
    #[case] expected: AtomForm,
    #[values(false, true)] journaled: bool,
) {
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyAtomField {
        id: AtomHandle::Id(AtomId(0)),
        change,
    };
    let undo = if journaled {
        molecule
            .apply_edit_with_undo(edit.clone(), &mut state)
            .unwrap()
    } else {
        molecule.apply_edit(edit.clone(), &mut state).unwrap();
        None
    };
    removal_entries.atoms[0] = expected;
    let expected = Molecule::from_entries(removal_entries);
    assert_eq!(molecule, expected);

    let repeated = if journaled {
        molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
    } else {
        molecule.apply_edit(edit, &mut state)
    };
    assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
    assert_eq!(molecule, expected);

    if journaled {
        molecule.apply_undo(undo.unwrap());
        assert!(molecule.normalized_eq(&initial));
    }
}

#[rstest]
#[case::elements(
    Edit::ModifyAtomElements { changes: vec![
        (AtomHandle::Id(AtomId(0)), ElementForm::Undetermined, ElementForm::Lit(Element::C)),
        (AtomHandle::New(0), ElementForm::Undetermined, ElementForm::Lit(Element::C)),
        (AtomHandle::Id(AtomId(0)), ElementForm::Lit(Element::C), ElementForm::Lit(Element::N)),
] },
    vec![
        Edit::ModifyAtomField { id: AtomHandle::Id(AtomId(0)), change: AtomFieldChange::Element { old: ElementForm::Undetermined, new: ElementForm::Lit(Element::C) } },
        Edit::ModifyAtomField { id: AtomHandle::New(0), change: AtomFieldChange::Element { old: ElementForm::Undetermined, new: ElementForm::Lit(Element::C) } },
        Edit::ModifyAtomField { id: AtomHandle::Id(AtomId(0)), change: AtomFieldChange::Element { old: ElementForm::Lit(Element::C), new: ElementForm::Lit(Element::N) } },
    ],
    AtomForm { element: ElementForm::Lit(Element::N), ..Default::default() },
    AtomForm { element: ElementForm::Lit(Element::C), ..Default::default() },
)]
#[case::isotope_masses(
    Edit::ModifyAtomIsotopeMasses { changes: vec![
        (AtomHandle::Id(AtomId(0)), IsotopeMassForm::Undetermined, IsotopeMassForm::Lit(12)),
        (AtomHandle::New(0), IsotopeMassForm::Undetermined, IsotopeMassForm::Lit(12)),
        (AtomHandle::Id(AtomId(0)), IsotopeMassForm::Lit(12), IsotopeMassForm::Lit(13)),
] },
    vec![
        Edit::ModifyAtomField { id: AtomHandle::Id(AtomId(0)), change: AtomFieldChange::IsotopeMass { old: IsotopeMassForm::Undetermined, new: IsotopeMassForm::Lit(12) } },
        Edit::ModifyAtomField { id: AtomHandle::New(0), change: AtomFieldChange::IsotopeMass { old: IsotopeMassForm::Undetermined, new: IsotopeMassForm::Lit(12) } },
        Edit::ModifyAtomField { id: AtomHandle::Id(AtomId(0)), change: AtomFieldChange::IsotopeMass { old: IsotopeMassForm::Lit(12), new: IsotopeMassForm::Lit(13) } },
    ],
    AtomForm { isotope_mass: IsotopeMassForm::Lit(13), ..Default::default() },
    AtomForm { isotope_mass: IsotopeMassForm::Lit(12), ..Default::default() },
)]
#[case::charges(
    Edit::ModifyAtomCharges { changes: vec![
        (AtomHandle::Id(AtomId(0)), NumForm::Undetermined, NumForm::Lit(1)),
        (AtomHandle::New(0), NumForm::Undetermined, NumForm::Lit(1)),
        (AtomHandle::Id(AtomId(0)), NumForm::Lit(1), NumForm::Lit(-1)),
] },
    vec![
        Edit::ModifyAtomField { id: AtomHandle::Id(AtomId(0)), change: AtomFieldChange::Charge { old: NumForm::Undetermined, new: NumForm::Lit(1) } },
        Edit::ModifyAtomField { id: AtomHandle::New(0), change: AtomFieldChange::Charge { old: NumForm::Undetermined, new: NumForm::Lit(1) } },
        Edit::ModifyAtomField { id: AtomHandle::Id(AtomId(0)), change: AtomFieldChange::Charge { old: NumForm::Lit(1), new: NumForm::Lit(-1) } },
    ],
    AtomForm { charge: NumForm::Lit(-1), ..Default::default() },
    AtomForm { charge: NumForm::Lit(1), ..Default::default() },
)]
#[case::implicit_hydrogens(
    Edit::ModifyAtomImplicitHydrogens { changes: vec![
        (AtomHandle::Id(AtomId(0)), NumForm::Undetermined, NumForm::Lit(1)),
        (AtomHandle::New(0), NumForm::Undetermined, NumForm::Lit(1)),
        (AtomHandle::Id(AtomId(0)), NumForm::Lit(1), NumForm::Lit(2)),
] },
    vec![
        Edit::ModifyAtomField { id: AtomHandle::Id(AtomId(0)), change: AtomFieldChange::ImplicitHydrogens { old: NumForm::Undetermined, new: NumForm::Lit(1) } },
        Edit::ModifyAtomField { id: AtomHandle::New(0), change: AtomFieldChange::ImplicitHydrogens { old: NumForm::Undetermined, new: NumForm::Lit(1) } },
        Edit::ModifyAtomField { id: AtomHandle::Id(AtomId(0)), change: AtomFieldChange::ImplicitHydrogens { old: NumForm::Lit(1), new: NumForm::Lit(2) } },
    ],
    AtomForm { implicit_hydrogens: NumForm::Lit(2), ..Default::default() },
    AtomForm { implicit_hydrogens: NumForm::Lit(1), ..Default::default() },
)]
#[case::lone_pairs(
    Edit::ModifyAtomLonePairs { changes: vec![
        (AtomHandle::Id(AtomId(0)), NumForm::Undetermined, NumForm::Lit(1)),
        (AtomHandle::New(0), NumForm::Undetermined, NumForm::Lit(1)),
        (AtomHandle::Id(AtomId(0)), NumForm::Lit(1), NumForm::Lit(2)),
] },
    vec![
        Edit::ModifyAtomField { id: AtomHandle::Id(AtomId(0)), change: AtomFieldChange::LonePairs { old: NumForm::Undetermined, new: NumForm::Lit(1) } },
        Edit::ModifyAtomField { id: AtomHandle::New(0), change: AtomFieldChange::LonePairs { old: NumForm::Undetermined, new: NumForm::Lit(1) } },
        Edit::ModifyAtomField { id: AtomHandle::Id(AtomId(0)), change: AtomFieldChange::LonePairs { old: NumForm::Lit(1), new: NumForm::Lit(2) } },
    ],
    AtomForm { lone_pairs: NumForm::Lit(2), ..Default::default() },
    AtomForm { lone_pairs: NumForm::Lit(1), ..Default::default() },
)]
#[case::unpaired_electrons(
    Edit::ModifyAtomUnpairedElectrons { changes: vec![
        (AtomHandle::Id(AtomId(0)), UnpairedElectronsForm::default(), UnpairedElectronsForm::closed_shell()),
        (AtomHandle::New(0), UnpairedElectronsForm::default(), UnpairedElectronsForm::closed_shell()),
        (AtomHandle::Id(AtomId(0)), UnpairedElectronsForm::closed_shell(), UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) }),
] },
    vec![
        Edit::ModifyAtomField { id: AtomHandle::Id(AtomId(0)), change: AtomFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::default(), new: UnpairedElectronsForm::closed_shell() } },
        Edit::ModifyAtomField { id: AtomHandle::New(0), change: AtomFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::default(), new: UnpairedElectronsForm::closed_shell() } },
        Edit::ModifyAtomField { id: AtomHandle::Id(AtomId(0)), change: AtomFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::closed_shell(), new: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) } } },
    ],
    AtomForm { unpaired_electrons: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) }, ..Default::default() },
    AtomForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() },
)]
fn test_molecule_apply_atom_batches(
    #[case] batch: Edit,
    #[case] singles: Vec<Edit>,
    #[case] changed: AtomForm,
    #[case] created: AtomForm,
) {
    let source = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::default(); 2],
        ..Default::default()
    });
    let addition = Edit::AddAtoms {
        atoms: vec![AtomForm::default()],
    };
    let batches = Edits::from_iter([addition.clone(), batch]);
    let singles = Edits::from_iter([addition].into_iter().chain(singles));
    let expected = Molecule::from_entries(MoleculeEntries {
        atoms: vec![changed, AtomForm::default(), created],
        ..Default::default()
    });
    let actual = source.clone().tracked_apply(batches.clone()).unwrap();
    assert_eq!(actual.0, expected);
    assert_eq!(actual, source.clone().tracked_apply(singles).unwrap());
    for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
        let mut molecule = source.clone();
        Transaction::run(
            &mut molecule,
            recovery,
            |mut transaction| -> Result<(), MoleculeApplyError> {
                transaction.apply(batches.clone())?;
                assert_eq!(transaction.probe()?, &expected);
                transaction.rollback();
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(molecule, source);
    }
}

#[rstest]
#[case(Edit::ModifyAtomElements { changes: vec![] })]
#[case(Edit::ModifyAtomIsotopeMasses { changes: vec![] })]
#[case(Edit::ModifyAtomCharges { changes: vec![] })]
#[case(Edit::ModifyAtomImplicitHydrogens { changes: vec![] })]
#[case(Edit::ModifyAtomLonePairs { changes: vec![] })]
#[case(Edit::ModifyAtomUnpairedElectrons { changes: vec![] })]
fn test_molecule_apply_atom_batches_identity(#[case] edit: Edit) {
    let source = Molecule::new();
    assert_eq!(
        source
            .clone()
            .apply(Edits::from_iter([edit.clone()]))
            .unwrap(),
        source
    );
    let mut molecule = source.clone();
    let mut state = ApplicationState::new(&molecule);
    let undo = molecule
        .apply_edit_with_undo(edit, &mut state)
        .unwrap()
        .unwrap();
    molecule.apply_undo(undo);
    assert_eq!(molecule, source);
}

#[rstest]
fn test_molecule_apply_atom_batches_normalized() {
    let accepted = NumForm::LitSet(Box::new([1].into()));
    let source = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm {
            charge: NumForm::Lit(1),
            ..Default::default()
        }],
        ..Default::default()
    });
    let mut molecule = source;
    let mut state = ApplicationState::new(&molecule);
    let undo = molecule
        .apply_edit_with_undo(
            Edit::ModifyAtomCharges {
                changes: vec![
                    (AtomHandle::Id(AtomId(0)), accepted.clone(), NumForm::Lit(2)),
                    (AtomHandle::Id(AtomId(0)), NumForm::Lit(2), NumForm::Lit(3)),
                ],
            },
            &mut state,
        )
        .unwrap()
        .unwrap();
    assert_eq!(
        undo,
        Undo::ModifyAtomCharges {
            changes: vec![
                (AtomId(0), NumForm::Lit(2), accepted.clone()),
                (AtomId(0), NumForm::Lit(3), NumForm::Lit(2)),
            ]
        }
    );
    molecule.apply_undo(undo);
    assert_eq!(
        molecule,
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm {
                charge: accepted,
                ..Default::default()
            }],
            ..Default::default()
        })
    );
}

#[rstest]
#[case::elements(Edit::ModifyAtomElements { changes: (0..3).map(|id| (AtomHandle::Id(AtomId(id)), ElementForm::Undetermined, ElementForm::Lit(Element::C))).collect() }, AtomForm { element: ElementForm::Lit(Element::C), ..Default::default() })]
#[case::isotope_masses(Edit::ModifyAtomIsotopeMasses { changes: (0..3).map(|id| (AtomHandle::Id(AtomId(id)), IsotopeMassForm::Undetermined, IsotopeMassForm::Lit(12))).collect() }, AtomForm { isotope_mass: IsotopeMassForm::Lit(12), ..Default::default() })]
#[case::charges(Edit::ModifyAtomCharges { changes: (0..3).map(|id| (AtomHandle::Id(AtomId(id)), NumForm::Undetermined, NumForm::Lit(1))).collect() }, AtomForm { charge: NumForm::Lit(1), ..Default::default() })]
#[case::implicit_hydrogens(Edit::ModifyAtomImplicitHydrogens { changes: (0..3).map(|id| (AtomHandle::Id(AtomId(id)), NumForm::Undetermined, NumForm::Lit(1))).collect() }, AtomForm { implicit_hydrogens: NumForm::Lit(1), ..Default::default() })]
#[case::lone_pairs(Edit::ModifyAtomLonePairs { changes: (0..3).map(|id| (AtomHandle::Id(AtomId(id)), NumForm::Undetermined, NumForm::Lit(1))).collect() }, AtomForm { lone_pairs: NumForm::Lit(1), ..Default::default() })]
#[case::unpaired_electrons(Edit::ModifyAtomUnpairedElectrons { changes: (0..3).map(|id| (AtomHandle::Id(AtomId(id)), UnpairedElectronsForm::default(), UnpairedElectronsForm::closed_shell())).collect() }, AtomForm { unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default() })]
fn test_molecule_apply_atom_batches_error(
    #[case] mut batch: Edit,
    #[case] changed: AtomForm,
    #[values(0, 1, 2)] failure: usize,
    #[values(false, true)] missing: bool,
    #[values(false, true)] journaled: bool,
) {
    match &mut batch {
        Edit::ModifyAtomElements { changes } => {
            if missing {
                changes[failure].0 = AtomHandle::Id(AtomId(4));
            } else {
                changes[failure].1 = changes[failure].2.clone();
            }
        }
        Edit::ModifyAtomIsotopeMasses { changes } => {
            if missing {
                changes[failure].0 = AtomHandle::Id(AtomId(4));
            } else {
                changes[failure].1 = changes[failure].2.clone();
            }
        }
        Edit::ModifyAtomCharges { changes } => {
            if missing {
                changes[failure].0 = AtomHandle::Id(AtomId(4));
            } else {
                changes[failure].1 = changes[failure].2.clone();
            }
        }
        Edit::ModifyAtomImplicitHydrogens { changes } => {
            if missing {
                changes[failure].0 = AtomHandle::Id(AtomId(4));
            } else {
                changes[failure].1 = changes[failure].2.clone();
            }
        }
        Edit::ModifyAtomLonePairs { changes } => {
            if missing {
                changes[failure].0 = AtomHandle::Id(AtomId(4));
            } else {
                changes[failure].1 = changes[failure].2.clone();
            }
        }
        Edit::ModifyAtomUnpairedElectrons { changes } => {
            if missing {
                changes[failure].0 = AtomHandle::Id(AtomId(4));
            } else {
                changes[failure].1 = changes[failure].2.clone();
            }
        }
        _ => unreachable!(),
    }
    let source = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::default(); 4],
        ..Default::default()
    });
    let prefix = Edit::ModifyAtomField {
        id: AtomHandle::Id(AtomId(3)),
        change: AtomFieldChange::Charge {
            old: NumForm::Undetermined,
            new: NumForm::Lit(1),
        },
    };
    let error = if missing {
        TransactionError::HandleOutOfRange {
            kind: EntityKind::Atom,
            index: 4,
            count: 4,
        }
    } else {
        TransactionError::OldStateMismatch
    };
    let mut molecule = source.clone();
    let mut state = ApplicationState::new(&molecule);
    molecule.apply_edit(prefix.clone(), &mut state).unwrap();
    let result = if journaled {
        molecule
            .apply_edit_with_undo(batch.clone(), &mut state)
            .map(|_| ())
    } else {
        molecule.apply_edit(batch.clone(), &mut state)
    };
    assert_eq!(result, Err(error.clone()));
    let mut atoms = vec![AtomForm::default(); 4];
    atoms[3].charge = NumForm::Lit(1);
    if !journaled {
        atoms[..failure].fill(changed);
    }
    assert_eq!(
        molecule,
        Molecule::from_entries(MoleculeEntries {
            atoms,
            ..Default::default()
        })
    );
    for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
        let mut molecule = source.clone();
        assert_eq!(
            molecule.transact(
                [Edits::from_iter([prefix.clone(), batch.clone()])],
                recovery
            ),
            Err(error.clone().into())
        );
        assert_eq!(molecule, source);
    }
}

#[rstest]
#[case::order(
    BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
    BondForm::from_order(2),
)]
#[case::order_equivalent(
    BondFieldChange::Order { old: NumForm::lit_set([1]), new: NumForm::Lit(2) },
    BondForm::from_order(2),
)]
#[case::charge(
    BondFieldChange::Charge { old: NumForm::Undetermined, new: NumForm::Lit(1) },
    BondForm { charge: NumForm::Lit(1), ..BondForm::from_order(1) },
)]
#[case::unpaired_electrons(
    BondFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::default(), new: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) } },
    BondForm { unpaired_electrons: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) }, ..BondForm::from_order(1) },
)]
fn test_molecule_apply_edit_modify_bond_field(
    mut removal_entries: MoleculeEntries,
    #[case] change: BondFieldChange,
    #[case] expected: BondForm,
    #[values(false, true)] journaled: bool,
) {
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyBondField {
        id: BondHandle::Id(BondId(0)),
        change,
    };
    let undo = if journaled {
        molecule
            .apply_edit_with_undo(edit.clone(), &mut state)
            .unwrap()
    } else {
        molecule.apply_edit(edit.clone(), &mut state).unwrap();
        None
    };
    removal_entries.bonds[0].2 = expected;
    let expected = Molecule::from_entries(removal_entries);
    assert_eq!(molecule, expected);

    let repeated = if journaled {
        molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
    } else {
        molecule.apply_edit(edit, &mut state)
    };
    assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
    assert_eq!(molecule, expected);

    if journaled {
        molecule.apply_undo(undo.unwrap());
        assert!(molecule.normalized_eq(&initial));
    }
}

#[rstest]
#[case::order(
    DativeBondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
    DativeBondForm::from_order(2),
)]
#[case::order_equivalent(
    DativeBondFieldChange::Order { old: NumForm::lit_set([1]), new: NumForm::Lit(2) },
    DativeBondForm::from_order(2),
)]
fn test_molecule_apply_edit_modify_dative_bond_field(
    mut removal_entries: MoleculeEntries,
    #[case] change: DativeBondFieldChange,
    #[case] expected: DativeBondForm,
    #[values(false, true)] journaled: bool,
) {
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyDativeBondField {
        id: DativeBondHandle::Id(DativeBondId(0)),
        change,
    };
    let undo = if journaled {
        molecule
            .apply_edit_with_undo(edit.clone(), &mut state)
            .unwrap()
    } else {
        molecule.apply_edit(edit.clone(), &mut state).unwrap();
        None
    };
    removal_entries.dative[0].2 = expected;
    let expected = Molecule::from_entries(removal_entries);
    assert_eq!(molecule, expected);

    let repeated = if journaled {
        molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
    } else {
        molecule.apply_edit(edit, &mut state)
    };
    assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
    assert_eq!(molecule, expected);

    if journaled {
        molecule.apply_undo(undo.unwrap());
        assert!(molecule.normalized_eq(&initial));
    }
}

#[rstest]
#[case::electrons(
    AromaticSystemFieldChange::Electrons { old: ElectronCountsForm::Lit(vec![1, 2]), new: ElectronCountsForm::Lit(vec![2, 2]) },
    AromaticSystemForm::from_electrons(vec![2, 2]),
)]
#[case::electron_count_length(
    AromaticSystemFieldChange::Electrons { old: ElectronCountsForm::Lit(vec![1, 2]), new: ElectronCountsForm::Lit(vec![2]) },
    AromaticSystemForm::from_electrons(vec![2]),
)]
#[case::charge(
    AromaticSystemFieldChange::Charge { old: NumForm::Undetermined, new: NumForm::Lit(1) },
    AromaticSystemForm { charge: NumForm::Lit(1), ..AromaticSystemForm::from_electrons(vec![1, 2]) },
)]
#[case::unpaired_electrons(
    AromaticSystemFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::default(), new: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) } },
    AromaticSystemForm { unpaired_electrons: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) }, ..AromaticSystemForm::from_electrons(vec![1, 2]) },
)]
fn test_molecule_apply_edit_modify_aromatic_system_field(
    mut removal_entries: MoleculeEntries,
    #[case] change: AromaticSystemFieldChange,
    #[case] expected: AromaticSystemForm,
    #[values(false, true)] journaled: bool,
) {
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyAromaticSystemField {
        id: AromaticSystemHandle::Id(AromaticSystemId(0)),
        change,
    };
    let undo = if journaled {
        molecule
            .apply_edit_with_undo(edit.clone(), &mut state)
            .unwrap()
    } else {
        molecule.apply_edit(edit.clone(), &mut state).unwrap();
        None
    };
    removal_entries.aromatic[0].1 = expected;
    let expected = Molecule::from_entries(removal_entries);
    assert_eq!(molecule, expected);

    let repeated = if journaled {
        molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
    } else {
        molecule.apply_edit(edit, &mut state)
    };
    assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
    assert_eq!(molecule, expected);

    if journaled {
        molecule.apply_undo(undo.unwrap());
        assert!(molecule.normalized_eq(&initial));
    }
}

#[rstest]
#[case::electrons(
    MulticenterBondFieldChange::Electrons { old: ElectronCountsForm::Lit(vec![1, 2]), new: ElectronCountsForm::Lit(vec![2, 2]) },
    MulticenterBondForm::from_electrons(vec![2, 2]),
)]
#[case::electron_count_length(
    MulticenterBondFieldChange::Electrons { old: ElectronCountsForm::Lit(vec![1, 2]), new: ElectronCountsForm::Lit(vec![2, 1, 2]) },
    MulticenterBondForm::from_electrons(vec![2, 1, 2]),
)]
#[case::charge(
    MulticenterBondFieldChange::Charge { old: NumForm::Undetermined, new: NumForm::Lit(-1) },
    MulticenterBondForm { charge: NumForm::Lit(-1), ..MulticenterBondForm::from_electrons(vec![1, 2]) },
)]
#[case::unpaired_electrons(
    MulticenterBondFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::default(), new: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) } },
    MulticenterBondForm { unpaired_electrons: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) }, ..MulticenterBondForm::from_electrons(vec![1, 2]) },
)]
fn test_molecule_apply_edit_modify_multicenter_bond_field(
    mut removal_entries: MoleculeEntries,
    #[case] change: MulticenterBondFieldChange,
    #[case] expected: MulticenterBondForm,
    #[values(false, true)] journaled: bool,
) {
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyMulticenterBondField {
        id: MulticenterBondHandle::Id(MulticenterBondId(0)),
        change,
    };
    let undo = if journaled {
        molecule
            .apply_edit_with_undo(edit.clone(), &mut state)
            .unwrap()
    } else {
        molecule.apply_edit(edit.clone(), &mut state).unwrap();
        None
    };
    removal_entries.multicenter[0].1 = expected;
    let expected = Molecule::from_entries(removal_entries);
    assert_eq!(molecule, expected);

    let repeated = if journaled {
        molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
    } else {
        molecule.apply_edit(edit, &mut state)
    };
    assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
    assert_eq!(molecule, expected);

    if journaled {
        molecule.apply_undo(undo.unwrap());
        assert!(molecule.normalized_eq(&initial));
    }
}

#[rstest]
#[case::kind(
    NoncovalentBondFieldChange::Kind { old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic) },
    NoncovalentBondForm::from_kind(NoncovalentBondKind::Ionic),
)]
fn test_molecule_apply_edit_modify_noncovalent_bond_field(
    mut removal_entries: MoleculeEntries,
    #[case] change: NoncovalentBondFieldChange,
    #[case] expected: NoncovalentBondForm,
    #[values(false, true)] journaled: bool,
) {
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyNoncovalentBondField {
        id: NoncovalentBondHandle::Id(NoncovalentBondId(0)),
        change,
    };
    let undo = if journaled {
        molecule
            .apply_edit_with_undo(edit.clone(), &mut state)
            .unwrap()
    } else {
        molecule.apply_edit(edit.clone(), &mut state).unwrap();
        None
    };
    removal_entries.noncovalent[0].1 = expected;
    let expected = Molecule::from_entries(removal_entries);
    assert_eq!(molecule, expected);

    let repeated = if journaled {
        molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
    } else {
        molecule.apply_edit(edit, &mut state)
    };
    assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
    assert_eq!(molecule, expected);

    if journaled {
        molecule.apply_undo(undo.unwrap());
        assert!(molecule.normalized_eq(&initial));
    }
}

#[rstest]
#[case::configuration(
    StereoAtomFieldChange::Configuration { old: StereoConfigurationForm::Undetermined, new: StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)) },
    StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
)]
#[case::coset_range(
    StereoAtomFieldChange::Configuration { old: StereoConfigurationForm::Undetermined, new: StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(2)) },
    StereoAtomForm::new(StereoKind::Tetrahedral, 2_u32),
)]
fn test_molecule_apply_edit_modify_stereo_atom_field(
    mut removal_entries: MoleculeEntries,
    #[case] change: StereoAtomFieldChange,
    #[case] expected: StereoAtomForm,
    #[values(false, true)] journaled: bool,
) {
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyStereoAtomField {
        id: StereoAtomHandle::Id(StereoAtomId(0)),
        change,
    };
    let undo = if journaled {
        molecule
            .apply_edit_with_undo(edit.clone(), &mut state)
            .unwrap()
    } else {
        molecule.apply_edit(edit.clone(), &mut state).unwrap();
        None
    };
    removal_entries.stereo_atoms[0].2 = expected;
    let expected = Molecule::from_entries(removal_entries);
    assert_eq!(molecule, expected);

    let repeated = if journaled {
        molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
    } else {
        molecule.apply_edit(edit, &mut state)
    };
    assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
    assert_eq!(molecule, expected);

    if journaled {
        molecule.apply_undo(undo.unwrap());
        assert!(molecule.normalized_eq(&initial));
    }
}

#[rstest]
#[case::configuration(
    StereoBondFieldChange::Configuration { old: StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(1)), new: StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(0)) },
    StereoBondForm::new(StereoKind::CisTrans, 0_u32),
)]
#[case::configuration_equivalent(
    StereoBondFieldChange::Configuration { old: StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::lit_set([1])), new: StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(0)) },
    StereoBondForm::new(StereoKind::CisTrans, 0_u32),
)]
#[case::coset_range(
    StereoBondFieldChange::Configuration { old: StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(1)), new: StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(2)) },
    StereoBondForm::new(StereoKind::CisTrans, 2_u32),
)]
fn test_molecule_apply_edit_modify_stereo_bond_field(
    mut removal_entries: MoleculeEntries,
    #[case] change: StereoBondFieldChange,
    #[case] expected: StereoBondForm,
    #[values(false, true)] journaled: bool,
) {
    let initial = Molecule::from_entries(removal_entries.clone());
    let mut molecule = initial.clone();
    let mut state = ApplicationState::new(&molecule);
    let edit = Edit::ModifyStereoBondField {
        id: StereoBondHandle::Id(StereoBondId(0)),
        change,
    };
    let undo = if journaled {
        molecule
            .apply_edit_with_undo(edit.clone(), &mut state)
            .unwrap()
    } else {
        molecule.apply_edit(edit.clone(), &mut state).unwrap();
        None
    };
    removal_entries.stereo_bonds[0].2 = expected;
    let expected = Molecule::from_entries(removal_entries);
    assert_eq!(molecule, expected);

    let repeated = if journaled {
        molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
    } else {
        molecule.apply_edit(edit, &mut state)
    };
    assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
    assert_eq!(molecule, expected);

    if journaled {
        molecule.apply_undo(undo.unwrap());
        assert!(molecule.normalized_eq(&initial));
    }
}
