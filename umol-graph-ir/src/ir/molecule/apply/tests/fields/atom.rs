use super::*;

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
