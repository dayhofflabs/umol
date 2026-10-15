use proptest::prelude::*;
use rstest::{fixture, rstest};
use umol_graph_ir::ir::{
    AddAromaticSystem, AddMulticenterBond, AddStereoAtom, AddStereoBond, AromaticSystemId,
    AtomHandle, AtomId, BondFieldChange, BondHandle, BondId, DativeBondId, Edits, EntityKind,
    Molecule, MoleculeApplyError, MulticenterBondId, NoncovalentBondId, RollbackMethod,
    StereoAtomId, StereoBondId, TrackedTransaction, Transaction, TransactionError,
};

use crate::strategies::{molecule_entries_strategy, *};

proptest! {

    /// Field batches agree with ordered singles and the independently updated atom table;
    /// rollback restores the original raw forms after repeated writes to arbitrary targets.
    #[test]
    fn test_molecule_apply_atom_batches_elements(
        atoms in prop::collection::vec(atom_form_strategy(), 3),
        updates in prop::collection::vec((0usize..3, atom_form_strategy()), 0..16),
    ) {
        let source = Molecule::from_entries(MoleculeEntries { atoms: atoms.clone(), ..Default::default() });
        let mut expected_atoms = atoms;
        let mut rows = Vec::new();
        let mut singles = Edits::new();
        for (index, form) in updates {
            let old = expected_atoms[index].element.clone();
            let new = form.element;
            expected_atoms[index].element = new.clone();
            let id = AtomHandle::Id(AtomId::from(index));
            rows.push((id.clone(), old.clone(), new.clone()));
            singles.push(Edit::ModifyAtomField { id, change: AtomFieldChange::Element { old, new } });
        }
        let expected = Molecule::from_entries(MoleculeEntries { atoms: expected_atoms, ..Default::default() });
        let mut batch = Edits::new();
        batch.modify_atom_elements(rows);
        let actual = source.clone().tracked_apply(batch.clone()).unwrap();
        prop_assert_eq!(&actual.0, &expected);
        prop_assert_eq!(actual, source.clone().tracked_apply(singles).unwrap());
        for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
            let mut restored = source.clone();
            Transaction::run(&mut restored, recovery, |mut transaction| -> Result<(), MoleculeApplyError> {
                transaction.apply(batch.clone())?;
                assert_eq!(transaction.probe()?, &expected);
                transaction.rollback();
                Ok(())
            }).unwrap();
            prop_assert_eq!(&restored, &source);
        }
    }

    /// Field batches agree with ordered singles and the independently updated atom table;
    /// rollback restores the original raw forms after repeated writes to arbitrary targets.
    #[test]
    fn test_molecule_apply_atom_batches_isotope_masses(
        atoms in prop::collection::vec(atom_form_strategy(), 3),
        updates in prop::collection::vec((0usize..3, atom_form_strategy()), 0..16),
    ) {
        let source = Molecule::from_entries(MoleculeEntries { atoms: atoms.clone(), ..Default::default() });
        let mut expected_atoms = atoms;
        let mut rows = Vec::new();
        let mut singles = Edits::new();
        for (index, form) in updates {
            let old = expected_atoms[index].isotope_mass.clone();
            let new = form.isotope_mass;
            expected_atoms[index].isotope_mass = new.clone();
            let id = AtomHandle::Id(AtomId::from(index));
            rows.push((id.clone(), old.clone(), new.clone()));
            singles.push(Edit::ModifyAtomField { id, change: AtomFieldChange::IsotopeMass { old, new } });
        }
        let expected = Molecule::from_entries(MoleculeEntries { atoms: expected_atoms, ..Default::default() });
        let mut batch = Edits::new();
        batch.modify_atom_isotope_masses(rows);
        let actual = source.clone().tracked_apply(batch.clone()).unwrap();
        prop_assert_eq!(&actual.0, &expected);
        prop_assert_eq!(actual, source.clone().tracked_apply(singles).unwrap());
        for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
            let mut restored = source.clone();
            Transaction::run(&mut restored, recovery, |mut transaction| -> Result<(), MoleculeApplyError> {
                transaction.apply(batch.clone())?;
                assert_eq!(transaction.probe()?, &expected);
                transaction.rollback();
                Ok(())
            }).unwrap();
            prop_assert_eq!(&restored, &source);
        }
    }

    /// Field batches agree with ordered singles and the independently updated atom table;
    /// rollback restores the original raw forms after repeated writes to arbitrary targets.
    #[test]
    fn test_molecule_apply_atom_batches_charges(
        atoms in prop::collection::vec(atom_form_strategy(), 3),
        updates in prop::collection::vec((0usize..3, atom_form_strategy()), 0..16),
    ) {
        let source = Molecule::from_entries(MoleculeEntries { atoms: atoms.clone(), ..Default::default() });
        let mut expected_atoms = atoms;
        let mut rows = Vec::new();
        let mut singles = Edits::new();
        for (index, form) in updates {
            let old = expected_atoms[index].charge.clone();
            let new = form.charge;
            expected_atoms[index].charge = new.clone();
            let id = AtomHandle::Id(AtomId::from(index));
            rows.push((id.clone(), old.clone(), new.clone()));
            singles.push(Edit::ModifyAtomField { id, change: AtomFieldChange::Charge { old, new } });
        }
        let expected = Molecule::from_entries(MoleculeEntries { atoms: expected_atoms, ..Default::default() });
        let mut batch = Edits::new();
        batch.modify_atom_charges(rows);
        let actual = source.clone().tracked_apply(batch.clone()).unwrap();
        prop_assert_eq!(&actual.0, &expected);
        prop_assert_eq!(actual, source.clone().tracked_apply(singles).unwrap());
        for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
            let mut restored = source.clone();
            Transaction::run(&mut restored, recovery, |mut transaction| -> Result<(), MoleculeApplyError> {
                transaction.apply(batch.clone())?;
                assert_eq!(transaction.probe()?, &expected);
                transaction.rollback();
                Ok(())
            }).unwrap();
            prop_assert_eq!(&restored, &source);
        }
    }

    /// Field batches agree with ordered singles and the independently updated atom table;
    /// rollback restores the original raw forms after repeated writes to arbitrary targets.
    #[test]
    fn test_molecule_apply_atom_batches_implicit_hydrogens(
        atoms in prop::collection::vec(atom_form_strategy(), 3),
        updates in prop::collection::vec((0usize..3, atom_form_strategy()), 0..16),
    ) {
        let source = Molecule::from_entries(MoleculeEntries { atoms: atoms.clone(), ..Default::default() });
        let mut expected_atoms = atoms;
        let mut rows = Vec::new();
        let mut singles = Edits::new();
        for (index, form) in updates {
            let old = expected_atoms[index].implicit_hydrogens.clone();
            let new = form.implicit_hydrogens;
            expected_atoms[index].implicit_hydrogens = new.clone();
            let id = AtomHandle::Id(AtomId::from(index));
            rows.push((id.clone(), old.clone(), new.clone()));
            singles.push(Edit::ModifyAtomField { id, change: AtomFieldChange::ImplicitHydrogens { old, new } });
        }
        let expected = Molecule::from_entries(MoleculeEntries { atoms: expected_atoms, ..Default::default() });
        let mut batch = Edits::new();
        batch.modify_atom_implicit_hydrogens(rows);
        let actual = source.clone().tracked_apply(batch.clone()).unwrap();
        prop_assert_eq!(&actual.0, &expected);
        prop_assert_eq!(actual, source.clone().tracked_apply(singles).unwrap());
        for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
            let mut restored = source.clone();
            Transaction::run(&mut restored, recovery, |mut transaction| -> Result<(), MoleculeApplyError> {
                transaction.apply(batch.clone())?;
                assert_eq!(transaction.probe()?, &expected);
                transaction.rollback();
                Ok(())
            }).unwrap();
            prop_assert_eq!(&restored, &source);
        }
    }

    /// Field batches agree with ordered singles and the independently updated atom table;
    /// rollback restores the original raw forms after repeated writes to arbitrary targets.
    #[test]
    fn test_molecule_apply_atom_batches_lone_pairs(
        atoms in prop::collection::vec(atom_form_strategy(), 3),
        updates in prop::collection::vec((0usize..3, atom_form_strategy()), 0..16),
    ) {
        let source = Molecule::from_entries(MoleculeEntries { atoms: atoms.clone(), ..Default::default() });
        let mut expected_atoms = atoms;
        let mut rows = Vec::new();
        let mut singles = Edits::new();
        for (index, form) in updates {
            let old = expected_atoms[index].lone_pairs.clone();
            let new = form.lone_pairs;
            expected_atoms[index].lone_pairs = new.clone();
            let id = AtomHandle::Id(AtomId::from(index));
            rows.push((id.clone(), old.clone(), new.clone()));
            singles.push(Edit::ModifyAtomField { id, change: AtomFieldChange::LonePairs { old, new } });
        }
        let expected = Molecule::from_entries(MoleculeEntries { atoms: expected_atoms, ..Default::default() });
        let mut batch = Edits::new();
        batch.modify_atom_lone_pairs(rows);
        let actual = source.clone().tracked_apply(batch.clone()).unwrap();
        prop_assert_eq!(&actual.0, &expected);
        prop_assert_eq!(actual, source.clone().tracked_apply(singles).unwrap());
        for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
            let mut restored = source.clone();
            Transaction::run(&mut restored, recovery, |mut transaction| -> Result<(), MoleculeApplyError> {
                transaction.apply(batch.clone())?;
                assert_eq!(transaction.probe()?, &expected);
                transaction.rollback();
                Ok(())
            }).unwrap();
            prop_assert_eq!(&restored, &source);
        }
    }

    /// Field batches agree with ordered singles and the independently updated atom table;
    /// rollback restores the original raw forms after repeated writes to arbitrary targets.
    #[test]
    fn test_molecule_apply_atom_batches_unpaired_electrons(
        atoms in prop::collection::vec(atom_form_strategy(), 3),
        updates in prop::collection::vec((0usize..3, atom_form_strategy()), 0..16),
    ) {
        let source = Molecule::from_entries(MoleculeEntries { atoms: atoms.clone(), ..Default::default() });
        let mut expected_atoms = atoms;
        let mut rows = Vec::new();
        let mut singles = Edits::new();
        for (index, form) in updates {
            let old = expected_atoms[index].unpaired_electrons.clone();
            let new = form.unpaired_electrons;
            expected_atoms[index].unpaired_electrons = new.clone();
            let id = AtomHandle::Id(AtomId::from(index));
            rows.push((id.clone(), old.clone(), new.clone()));
            singles.push(Edit::ModifyAtomField { id, change: AtomFieldChange::UnpairedElectrons { old, new } });
        }
        let expected = Molecule::from_entries(MoleculeEntries { atoms: expected_atoms, ..Default::default() });
        let mut batch = Edits::new();
        batch.modify_atom_unpaired_electrons(rows);
        let actual = source.clone().tracked_apply(batch.clone()).unwrap();
        prop_assert_eq!(&actual.0, &expected);
        prop_assert_eq!(actual, source.clone().tracked_apply(singles).unwrap());
        for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
            let mut restored = source.clone();
            Transaction::run(&mut restored, recovery, |mut transaction| -> Result<(), MoleculeApplyError> {
                transaction.apply(batch.clone())?;
                assert_eq!(transaction.probe()?, &expected);
                transaction.rollback();
                Ok(())
            }).unwrap();
            prop_assert_eq!(&restored, &source);
        }
    }

    /// Rejection inside a field batch restores its prefix and the preceding journal entries;
    /// destructive editor application retains exactly the independently computed prefix.
    #[test]
    fn test_molecule_apply_atom_batches_failure(
        values in prop::collection::vec((-3i64..=3, 0usize..3), 0..16),
        split in 0usize..16,
    ) {
        let source = Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::default(); 3], ..Default::default() });
        let mut prefix = Edits::new();
        prefix.modify_atom_lone_pairs(vec![(AtomHandle::Id(AtomId(0)), NumForm::Undetermined, NumForm::Lit(1))]);
        let split = split.min(values.len());
        let mut expected_atoms = vec![AtomForm::default(); 3];
        expected_atoms[0].lone_pairs = NumForm::Lit(1);
        let mut rows = Vec::new();
        for (new, index) in &values[..split] {
            let old = expected_atoms[*index].charge.clone();
            expected_atoms[*index].charge = NumForm::Lit(*new);
            rows.push((AtomHandle::Id(AtomId::from(*index)), old, NumForm::Lit(*new)));
        }
        rows.push((AtomHandle::Id(AtomId(3)), NumForm::Undetermined, NumForm::Lit(1)));
        for (new, index) in &values[split..] {
            rows.push((AtomHandle::Id(AtomId::from(*index)), NumForm::Undetermined, NumForm::Lit(*new)));
        }
        let mut batch = Edits::new();
        batch.modify_atom_charges(rows);
        let error = MoleculeApplyError::Transaction(TransactionError::HandleOutOfRange { kind: EntityKind::Atom, index: 3, count: 3 });
        for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
            let mut molecule = source.clone();
            prop_assert_eq!(molecule.transact([prefix.clone(), batch.clone()], recovery), Err(error.clone()));
            prop_assert_eq!(molecule, source.clone());
        }
        let mut editor = source.edit();
        editor.apply(prefix).unwrap();
        prop_assert_eq!(editor.apply(batch), Err(error));
        prop_assert_eq!(editor.finish().unwrap(), Molecule::from_entries(MoleculeEntries { atoms: expected_atoms, ..Default::default() }));
    }

    /// Optional constraint batches agree with ordered single-key edits and preserve the raw
    /// source on rollback, including repeated add, replace, remove, and absent-to-absent rows.
    #[test]
    fn test_molecule_apply_atom_constraint_batches(
        updates in prop::collection::vec((0usize..3, prop::option::of(-2i64..=6)), 0..20),
    ) {
        let source = Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm::default(); 3], ..Default::default() });
        let mut values = vec![None; 3];
        let mut rows = Vec::new();
        let mut singles = Edits::new();
        for (index, value) in updates {
            let old = values[index].clone();
            let new = value.map(AtomConstraintForm::valence);
            values[index] = new.clone();
            let id = AtomHandle::Id(AtomId::from(index));
            rows.push((id.clone(), old.clone(), new.clone()));
            singles.push(Edit::ModifyAtomConstraint { id, old, new });
        }
        let expected = Molecule::from_entries(MoleculeEntries {
            atoms: values.into_iter().map(|value| AtomForm { constraints: value.into_iter().collect(), ..Default::default() }).collect(),
            ..Default::default()
        });
        let mut batch = Edits::new();
        batch.modify_atom_constraints(AtomConstraintKey::Valence, rows);
        let actual = source.clone().tracked_apply(batch.clone()).unwrap();
        prop_assert_eq!(&actual.0, &expected);
        prop_assert_eq!(actual, source.clone().tracked_apply(singles).unwrap());
        for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
            let mut restored = source.clone();
            Transaction::run(&mut restored, recovery, |mut transaction| -> Result<(), MoleculeApplyError> {
                transaction.apply(batch.clone())?;
                transaction.rollback();
                Ok(())
            }).unwrap();
            prop_assert_eq!(&restored, &source);
        }
    }

    /// Bond field batches agree with ordered singles and an independently replaced bond
    /// table; rollback restores raw attributes and topology after repeated writes.
    #[test]
    fn test_molecule_apply_bond_batches_orders(
        bonds in prop::collection::vec(bond_form_strategy(), 3),
        updates in prop::collection::vec((0usize..3, bond_form_strategy()), 0..16),
    ) {
        let entries = MoleculeEntries {
            atoms: vec![AtomForm::default(); 4],
            bonds: bonds.iter().enumerate().map(|(i, b)| (AtomId::from(i), AtomId::from(i + 1), b.clone())).collect(),
            ..Default::default()
        };
        let source = Molecule::from_entries(entries.clone());
        let mut expected = entries;
        let mut rows = Vec::new();
        let mut singles = Edits::new();
        for (index, form) in updates {
            let old = expected.bonds[index].2.order.clone();
            let new = form.order;
            expected.bonds[index].2.order = new.clone();
            let id = BondHandle::Id(BondId::from(index));
            rows.push((id.clone(), old.clone(), new.clone()));
            singles.push(Edit::ModifyBondField { id, change: BondFieldChange::Order { old, new } });
        }
        let expected = Molecule::from_entries(expected);
        let mut batch = Edits::new();
        batch.modify_bond_orders(rows);
        let actual = source.clone().tracked_apply(batch.clone()).unwrap();
        prop_assert_eq!(&actual.0, &expected);
        prop_assert_eq!(actual, source.clone().tracked_apply(singles).unwrap());
        for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
            let mut molecule = source.clone();
            Transaction::run(&mut molecule, recovery, |mut transaction| -> Result<(), MoleculeApplyError> {
                transaction.apply(batch.clone())?;
                assert_eq!(transaction.probe()?, &expected);
                transaction.rollback();
                Ok(())
            }).unwrap();
            prop_assert_eq!(&molecule, &source);
        }
    }

    /// Bond field batches agree with ordered singles and an independently replaced bond
    /// table; rollback restores raw attributes and topology after repeated writes.
    #[test]
    fn test_molecule_apply_bond_batches_charges(
        bonds in prop::collection::vec(bond_form_strategy(), 3),
        updates in prop::collection::vec((0usize..3, bond_form_strategy()), 0..16),
    ) {
        let entries = MoleculeEntries {
            atoms: vec![AtomForm::default(); 4],
            bonds: bonds.iter().enumerate().map(|(i, b)| (AtomId::from(i), AtomId::from(i + 1), b.clone())).collect(),
            ..Default::default()
        };
        let source = Molecule::from_entries(entries.clone());
        let mut expected = entries;
        let mut rows = Vec::new();
        let mut singles = Edits::new();
        for (index, form) in updates {
            let old = expected.bonds[index].2.charge.clone();
            let new = form.charge;
            expected.bonds[index].2.charge = new.clone();
            let id = BondHandle::Id(BondId::from(index));
            rows.push((id.clone(), old.clone(), new.clone()));
            singles.push(Edit::ModifyBondField { id, change: BondFieldChange::Charge { old, new } });
        }
        let expected = Molecule::from_entries(expected);
        let mut batch = Edits::new();
        batch.modify_bond_charges(rows);
        let actual = source.clone().tracked_apply(batch.clone()).unwrap();
        prop_assert_eq!(&actual.0, &expected);
        prop_assert_eq!(actual, source.clone().tracked_apply(singles).unwrap());
        for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
            let mut molecule = source.clone();
            Transaction::run(&mut molecule, recovery, |mut transaction| -> Result<(), MoleculeApplyError> {
                transaction.apply(batch.clone())?;
                assert_eq!(transaction.probe()?, &expected);
                transaction.rollback();
                Ok(())
            }).unwrap();
            prop_assert_eq!(&molecule, &source);
        }
    }

    /// Bond field batches agree with ordered singles and an independently replaced bond
    /// table; rollback restores raw attributes and topology after repeated writes.
    #[test]
    fn test_molecule_apply_bond_batches_unpaired_electrons(
        bonds in prop::collection::vec(bond_form_strategy(), 3),
        updates in prop::collection::vec((0usize..3, bond_form_strategy()), 0..16),
    ) {
        let entries = MoleculeEntries {
            atoms: vec![AtomForm::default(); 4],
            bonds: bonds.iter().enumerate().map(|(i, b)| (AtomId::from(i), AtomId::from(i + 1), b.clone())).collect(),
            ..Default::default()
        };
        let source = Molecule::from_entries(entries.clone());
        let mut expected = entries;
        let mut rows = Vec::new();
        let mut singles = Edits::new();
        for (index, form) in updates {
            let old = expected.bonds[index].2.unpaired_electrons.clone();
            let new = form.unpaired_electrons;
            expected.bonds[index].2.unpaired_electrons = new.clone();
            let id = BondHandle::Id(BondId::from(index));
            rows.push((id.clone(), old.clone(), new.clone()));
            singles.push(Edit::ModifyBondField { id, change: BondFieldChange::UnpairedElectrons { old, new } });
        }
        let expected = Molecule::from_entries(expected);
        let mut batch = Edits::new();
        batch.modify_bond_unpaired_electrons(rows);
        let actual = source.clone().tracked_apply(batch.clone()).unwrap();
        prop_assert_eq!(&actual.0, &expected);
        prop_assert_eq!(actual, source.clone().tracked_apply(singles).unwrap());
        for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
            let mut molecule = source.clone();
            Transaction::run(&mut molecule, recovery, |mut transaction| -> Result<(), MoleculeApplyError> {
                transaction.apply(batch.clone())?;
                assert_eq!(transaction.probe()?, &expected);
                transaction.rollback();
                Ok(())
            }).unwrap();
            prop_assert_eq!(&molecule, &source);
        }
    }

    /// Optional ring constraints agree with ordered single-key edits over generated scopes
    /// and non-literal counts; rollback restores absent and stored constraints exactly.
    #[test]
    fn test_molecule_apply_bond_constraint_batches(
        scope in prop_oneof![Just(RingScope::All), (3u8..=8).prop_map(RingScope::Size)],
        updates in prop::collection::vec((0usize..3, prop::option::of(any_num_form_strategy())), 0..20),
    ) {
        let source = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::default(); 4],
            bonds: (0..3).map(|i| (AtomId(i), AtomId(i + 1), BondForm::default())).collect(),
            ..Default::default()
        });
        let mut values = vec![None; 3];
        let mut rows = Vec::new();
        let mut singles = Edits::new();
        for (index, value) in updates {
            let old = values[index].clone();
            let new = value.map(|v| BondConstraintForm::ring_membership(scope, v));
            values[index] = new.clone();
            let id = BondHandle::Id(BondId::from(index));
            rows.push((id.clone(), old.clone(), new.clone()));
            singles.push(Edit::ModifyBondConstraint { id, old, new });
        }
        let expected = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::default(); 4],
            bonds: values.into_iter().enumerate().map(|(i, value)| (AtomId::from(i), AtomId::from(i + 1), BondForm { constraints: value.into_iter().collect(), ..Default::default() })).collect(),
            ..Default::default()
        });
        let mut batch = Edits::new();
        batch.modify_bond_constraints(BondConstraintKey::RingMembership(scope), rows);
        let actual = source.clone().tracked_apply(batch.clone()).unwrap();
        prop_assert_eq!(&actual.0, &expected);
        prop_assert_eq!(actual, source.clone().tracked_apply(singles).unwrap());
        for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
            let mut restored = source.clone();
            Transaction::run(&mut restored, recovery, |mut transaction| -> Result<(), MoleculeApplyError> {
                transaction.apply(batch.clone())?;
                assert_eq!(transaction.probe()?, &expected);
                transaction.rollback();
                Ok(())
            }).unwrap();
            prop_assert_eq!(&restored, &source);
        }
    }

    /// A failing bond row restores the preceding atom edit and completed bond prefix in
    /// either transaction mode; a destructive editor retains precisely that prefix.
    #[test]
    fn test_molecule_apply_bond_batches_failure(
        values in prop::collection::vec((0usize..3, 1i64..=3), 1..16),
        split in 0usize..16,
    ) {
        let split = split.min(values.len());
        let entries = MoleculeEntries { atoms: vec![AtomForm::default(); 4], bonds: (0..3).map(|i| (AtomId(i), AtomId(i+1), BondForm::default())).collect(), ..Default::default() };
        let source = Molecule::from_entries(entries.clone());
        let mut prefix = Edits::new();
        prefix.push(Edit::ModifyAtomField { id: AtomHandle::Id(AtomId(0)), change: AtomFieldChange::Charge { old: NumForm::Undetermined, new: NumForm::Lit(1) } });
        let mut expected = entries;
        expected.atoms[0].charge = NumForm::Lit(1);
        let mut rows = Vec::new();
        for (index, value) in &values[..split] {
            let old = expected.bonds[*index].2.order.clone();
            expected.bonds[*index].2.order = NumForm::Lit(*value);
            rows.push((BondHandle::Id(BondId::from(*index)), old, NumForm::Lit(*value)));
        }
        rows.push((BondHandle::Id(BondId(3)), NumForm::Undetermined, NumForm::Lit(1)));
        for (index, value) in &values[split..] { rows.push((BondHandle::Id(BondId::from(*index)), NumForm::Undetermined, NumForm::Lit(*value))); }
        let mut batch = Edits::new();
        batch.modify_bond_orders(rows);
        let error = MoleculeApplyError::Transaction(TransactionError::HandleOutOfRange { kind: EntityKind::Bond, index: 3, count: 3 });
        for recovery in [RollbackMethod::Journal, RollbackMethod::Clone] {
            let mut molecule = source.clone();
            prop_assert_eq!(molecule.transact([prefix.clone(), batch.clone()], recovery), Err(error.clone()));
            prop_assert_eq!(molecule, source.clone());
        }
        let mut editor = source.edit();
        editor.apply(prefix).unwrap();
        prop_assert_eq!(editor.apply(batch), Err(error));
        prop_assert_eq!(editor.finish().unwrap(), Molecule::from_entries(expected));
    }

    /// Consecutive bulk additions reconstruct the supplied complete entries for any
    /// batch split. The constructor supplies the independent publication result.
    #[test]
    fn test_molecule_editor_finish_bulk_additions(
        entries in molecule_entries_strategy(),
        batch_size in 1usize..=4,
    ) {
        let expected = Molecule::from_entries(entries.clone());
        let mut editor = Molecule::new().edit();

        for (index, batch) in entries.atoms.chunks(batch_size).enumerate() {
            let ids = editor.add_atoms(batch.to_vec());
            let start = index * batch_size;
            let expected_ids: Vec<_> = (start..start + batch.len()).map(AtomId::from).collect();
            prop_assert_eq!(ids.collect::<Vec<_>>(), expected_ids);
        }

        for (index, batch) in entries.bonds.chunks(batch_size).enumerate() {
            let ids = editor.add_bonds(
                batch.iter()
                    .map(|(first, second, attributes)| ([*first, *second], attributes.clone()))
                    .collect(),
            );
            let start = index * batch_size;
            let expected_ids: Vec<_> = (start..start + batch.len()).map(BondId::from).collect();
            prop_assert_eq!(ids.collect::<Vec<_>>(), expected_ids);
        }

        for (index, batch) in entries.dative.chunks(batch_size).enumerate() {
            let ids = editor.add_dative_bonds(
                batch.iter()
                    .map(|(donors, acceptor, attributes)| (donors.as_slice(), *acceptor, attributes.clone()))
                    .collect(),
            );
            let start = index * batch_size;
            let expected_ids: Vec<_> = (start..start + batch.len()).map(DativeBondId::from).collect();
            prop_assert_eq!(ids.collect::<Vec<_>>(), expected_ids);
        }

        for (index, batch) in entries.aromatic.chunks(batch_size).enumerate() {
            let ids = editor.add_aromatic_systems(
                batch.iter()
                    .map(|(atoms, attributes)| (atoms.as_slice(), attributes.clone()))
                    .collect(),
            );
            let start = index * batch_size;
            let expected_ids: Vec<_> = (start..start + batch.len()).map(AromaticSystemId::from).collect();
            prop_assert_eq!(ids.collect::<Vec<_>>(), expected_ids);
        }

        for (index, batch) in entries.multicenter.chunks(batch_size).enumerate() {
            let ids = editor.add_multicenter_bonds(
                batch.iter()
                    .map(|(atoms, attributes)| (atoms.as_slice(), attributes.clone()))
                    .collect(),
            );
            let start = index * batch_size;
            let expected_ids: Vec<_> = (start..start + batch.len()).map(MulticenterBondId::from).collect();
            prop_assert_eq!(ids.collect::<Vec<_>>(), expected_ids);
        }

        for (index, batch) in entries.noncovalent.chunks(batch_size).enumerate() {
            let ids = editor.add_noncovalent_bonds(batch.to_vec());
            let start = index * batch_size;
            let expected_ids: Vec<_> = (start..start + batch.len()).map(NoncovalentBondId::from).collect();
            prop_assert_eq!(ids.collect::<Vec<_>>(), expected_ids);
        }

        for (index, batch) in entries.stereo_atoms.chunks(batch_size).enumerate() {
            let ids = editor.add_stereo_atoms(
                batch.iter()
                    .map(|(site, ligands, attributes)| (*site, ligands.as_slice(), attributes.clone()))
                    .collect(),
            );
            let start = index * batch_size;
            let expected_ids: Vec<_> = (start..start + batch.len()).map(StereoAtomId::from).collect();
            prop_assert_eq!(ids.collect::<Vec<_>>(), expected_ids);
        }

        for (index, batch) in entries.stereo_bonds.chunks(batch_size).enumerate() {
            let ids = editor.add_stereo_bonds(
                batch.iter()
                    .map(|(site, ligands, attributes)| (*site, ligands.as_slice(), attributes.clone()))
                    .collect(),
            );
            let start = index * batch_size;
            let expected_ids: Vec<_> = (start..start + batch.len()).map(StereoBondId::from).collect();
            prop_assert_eq!(ids.collect::<Vec<_>>(), expected_ids);
        }
        *editor.constraints_mut() = entries.constraints;
        prop_assert_eq!(editor.finish(), Ok(expected));
    }

    /// The standalone edit surface preserves every ordered raw edit, including repeated entries,
    /// while rebuilding the per-kind creation ordinals through `Edits` construction.
    #[test]
    fn test_edits_dsl_roundtrip(edits in edits_dsl_strategy()) {
        let defaults = MoleculeDefaults::new();
        let dsl = EditsDsl::from_ir(&edits, &defaults);
        let rendered = dsl.to_edn();
        let parsed = EditsDsl::from_edn(&rendered)
            .map_err(|error| TestCaseError::fail(format!("edit parse failed: {error}")))?;
        let rebuilt = parsed.into_ir(&defaults);

        prop_assert_eq!(rebuilt, edits);
    }

    /// Creation ordinals are reconstructed solely from the ordered entries: uninterrupted
    /// construction, raw pushes, and `FromIterator` agree for arbitrary interleavings, with one
    /// independent namespace per entity kind.
    #[test]
    fn test_edits_creation_ordinals(
        kinds in prop::collection::vec(
            prop::sample::select(vec![
                EntityKind::Atom,
                EntityKind::Bond,
                EntityKind::DativeBond,
                EntityKind::AromaticSystem,
                EntityKind::MulticenterBond,
                EntityKind::NoncovalentBond,
                EntityKind::StereoAtom,
                EntityKind::StereoBond,
            ]),
            0..64,
        ),
    ) {
        let mut direct = Edits::new();
        let mut entries = Vec::with_capacity(kinds.len());
        let mut counts = [0usize; 8];
        for kind in kinds {
            match kind {
                EntityKind::Atom => {
                    counts[0] += 1;
                    direct.add_atom(AtomForm::default());
                    entries.push(Edit::AddAtoms { atoms: vec![AtomForm::default()] });
                }
                EntityKind::Bond => {
                    counts[1] += 1;
                    direct.add_bond(
                        AtomHandle::Id(AtomId(0)),
                        AtomHandle::Id(AtomId(1)),
                        BondForm::default(),
                    );
                    entries.push(Edit::AddBonds {
                        bonds: vec![AddBond {
                            endpoints: [
                                AtomHandle::Id(AtomId(0)),
                                AtomHandle::Id(AtomId(1)),
                            ],
                            attributes: BondForm::default(),
                        }],
                    });
                }
                EntityKind::DativeBond => {
                    counts[2] += 1;
                    direct.add_dative_bond(Vec::new(), AtomHandle::Id(AtomId(0)), DativeBondForm::default());
                    entries.push(Edit::AddDativeBond {
                        donors: Vec::new(),
                        acceptor: AtomHandle::Id(AtomId(0)),
                        attributes: DativeBondForm::default(),
                    });
                }
                EntityKind::AromaticSystem => {
                    counts[3] += 1;
                    direct.add_aromatic_system(Vec::new(), AromaticSystemForm::default());
                    entries.push(Edit::AddAromaticSystem(Box::new(AddAromaticSystem {
                        atoms: Vec::new(),
                        attributes: AromaticSystemForm::default(),
                    })));
                }
                EntityKind::MulticenterBond => {
                    counts[4] += 1;
                    direct.add_multicenter_bond(Vec::new(), MulticenterBondForm::default());
                    entries.push(Edit::AddMulticenterBond(Box::new(AddMulticenterBond {
                        atoms: Vec::new(),
                        attributes: MulticenterBondForm::default(),
                    })));
                }
                EntityKind::NoncovalentBond => {
                    counts[5] += 1;
                    direct.add_noncovalent_bond(
                        [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        NoncovalentBondForm::default(),
                    );
                    entries.push(Edit::AddNoncovalentBond {
                        atoms: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        attributes: NoncovalentBondForm::default(),
                    });
                }
                EntityKind::StereoAtom => {
                    counts[6] += 1;
                    direct.add_stereo_atom(
                        AtomHandle::Id(AtomId(0)),
                        Vec::new(),
                        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                    );
                    entries.push(Edit::AddStereoAtom(Box::new(AddStereoAtom {
                        site: AtomHandle::Id(AtomId(0)),
                        ligands: Vec::new(),
                        attributes: StereoAtomForm::new(
                            StereoKind::Tetrahedral,
                            StereoCoset::Lit(0),
                        ),
                    })));
                }
                EntityKind::StereoBond => {
                    counts[7] += 1;
                    direct.add_stereo_bond(
                        BondHandle::Id(BondId(0)),
                        Vec::new(),
                        StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
                    );
                    entries.push(Edit::AddStereoBond(Box::new(AddStereoBond {
                        site: BondHandle::Id(BondId(0)),
                        ligands: Vec::new(),
                        attributes: StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
                    })));
                }
            }
        }

        let mut pushed = Edits::new();
        for entry in entries.clone() {
            pushed.push(entry);
        }
        let mut collected: Edits = entries.clone().into_iter().collect();
        prop_assert_eq!(direct.as_slice(), entries.as_slice());
        prop_assert_eq!(pushed.as_slice(), entries.as_slice());
        prop_assert_eq!(collected.as_slice(), entries.as_slice());

        let expected = (
            AtomHandle::New(counts[0]),
            BondHandle::New(counts[1]),
            DativeBondHandle::New(counts[2]),
            AromaticSystemHandle::New(counts[3]),
            MulticenterBondHandle::New(counts[4]),
            NoncovalentBondHandle::New(counts[5]),
            StereoAtomHandle::New(counts[6]),
            StereoBondHandle::New(counts[7]),
        );
        let direct_next = (
            direct.add_atom(AtomForm::default()),
            direct.add_bond(
                AtomHandle::Id(AtomId(0)),
                AtomHandle::Id(AtomId(1)),
                BondForm::default(),
            ),
            direct.add_dative_bond(Vec::new(), AtomHandle::Id(AtomId(0)), DativeBondForm::default()),
            direct.add_aromatic_system(Vec::new(), AromaticSystemForm::default()),
            direct.add_multicenter_bond(Vec::new(), MulticenterBondForm::default()),
            direct.add_noncovalent_bond(
                [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                NoncovalentBondForm::default(),
            ),
            direct.add_stereo_atom(
                AtomHandle::Id(AtomId(0)),
                Vec::new(),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            direct.add_stereo_bond(
                BondHandle::Id(BondId(0)),
                Vec::new(),
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
            ),
        );
        let pushed_next = (
            pushed.add_atom(AtomForm::default()),
            pushed.add_bond(
                AtomHandle::Id(AtomId(0)),
                AtomHandle::Id(AtomId(1)),
                BondForm::default(),
            ),
            pushed.add_dative_bond(Vec::new(), AtomHandle::Id(AtomId(0)), DativeBondForm::default()),
            pushed.add_aromatic_system(Vec::new(), AromaticSystemForm::default()),
            pushed.add_multicenter_bond(Vec::new(), MulticenterBondForm::default()),
            pushed.add_noncovalent_bond(
                [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                NoncovalentBondForm::default(),
            ),
            pushed.add_stereo_atom(
                AtomHandle::Id(AtomId(0)),
                Vec::new(),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            pushed.add_stereo_bond(
                BondHandle::Id(BondId(0)),
                Vec::new(),
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
            ),
        );
        let collected_next = (
            collected.add_atom(AtomForm::default()),
            collected.add_bond(
                AtomHandle::Id(AtomId(0)),
                AtomHandle::Id(AtomId(1)),
                BondForm::default(),
            ),
            collected.add_dative_bond(Vec::new(), AtomHandle::Id(AtomId(0)), DativeBondForm::default()),
            collected.add_aromatic_system(Vec::new(), AromaticSystemForm::default()),
            collected.add_multicenter_bond(Vec::new(), MulticenterBondForm::default()),
            collected.add_noncovalent_bond(
                [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                NoncovalentBondForm::default(),
            ),
            collected.add_stereo_atom(
                AtomHandle::Id(AtomId(0)),
                Vec::new(),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            collected.add_stereo_bond(
                BondHandle::Id(BondId(0)),
                Vec::new(),
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
            ),
        );
        prop_assert_eq!(&direct_next, &expected);
        prop_assert_eq!(&pushed_next, &expected);
        prop_assert_eq!(&collected_next, &expected);
    }

    /// Stable initial-host and created handles retain their labeled entity across arbitrary removal
    /// subsets. A later creation uses the next creation ordinal even when compaction reuses a
    /// concrete id; the expected molecule is obtained by filtering the original label sequences.
    #[test]
    fn test_molecule_transact_handle_identity(
        trace in stable_atom_handle_trace_strategy(false),
    ) {
        let expected = trace.expected();
        let mut molecule = trace.base();
        molecule.transact([trace.edits()], RollbackMethod::Journal).unwrap();
        prop_assert_eq!(molecule, expected);
    }

    /// Removing either an initial-host or created target leaves a tombstone in its own handle table.
    /// Failure after additional compaction and creation rolls the complete transaction back.
    #[test]
    fn test_molecule_transact_handle_error(
        trace in stable_atom_handle_trace_strategy(true),
    ) {
        let mut molecule = trace.base();
        let before = molecule.clone();
        prop_assert_eq!(
            molecule.transact([trace.edits()], RollbackMethod::Journal).unwrap_err(),
            MoleculeApplyError::Transaction(trace.expected_removed_error())
        );
        prop_assert_eq!(molecule, before);
    }

    /// A shuffled subset of entity creations can use `New(0)` independently for every kind; later
    /// field edits resolve each handle in its same-kind namespace regardless of creation order.
    #[test]
    fn test_molecule_editor_apply_handle_namespaces(
        kinds in transaction_entity_kind_order_strategy(),
    ) {
        let base = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
                AtomForm::from_element(Element::F),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        });
        let mut edits = Edits::new();
        for kind in &kinds {
            match kind {
                EntityKind::Atom => {
                    edits.add_atom(AtomForm::from_element(Element::P));
                }
                EntityKind::Bond => {
                    edits.add_bond(
                        AtomHandle::Id(AtomId(2)),
                        AtomHandle::Id(AtomId(3)),
                        BondForm::from_order(1),
                    );
                }
                EntityKind::DativeBond => {
                    edits.add_dative_bond(
                        vec![AtomHandle::Id(AtomId(0))],
                        AtomHandle::Id(AtomId(1)),
                        DativeBondForm::from_order(1),
                    );
                }
                EntityKind::AromaticSystem => {
                    edits.add_aromatic_system(
                        vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        AromaticSystemForm::default(),
                    );
                }
                EntityKind::MulticenterBond => {
                    edits.add_multicenter_bond(
                        vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        MulticenterBondForm::default(),
                    );
                }
                EntityKind::NoncovalentBond => {
                    edits.add_noncovalent_bond(
                        [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                    );
                }
                EntityKind::StereoAtom => {
                    edits.add_stereo_atom(
                        AtomHandle::Id(AtomId(0)),
                        (0..4_u32)
                            .map(|index| {
                                (AtomHandle::Id(AtomId(index)), StereoLigandKind::Atom)
                            })
                            .collect(),
                        StereoAtomForm::new(
                            StereoKind::Tetrahedral,
                            StereoCoset::Lit(1),
                        ),
                    );
                }
                EntityKind::StereoBond => {
                    edits.add_stereo_bond(
                        BondHandle::Id(BondId(0)),
                        (0..4_u32)
                            .map(|index| {
                                (AtomHandle::Id(AtomId(index)), StereoLigandKind::Atom)
                            })
                            .collect(),
                        StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                    );
                }
            }
        }
        for kind in &kinds {
            edits.push(match kind {
                EntityKind::Atom => Edit::ModifyAtomField {
                    id: AtomHandle::New(0),
                    change: AtomFieldChange::Charge {
                        old: NumForm::default(),
                        new: NumForm::Lit(1),
                    },
                },
                EntityKind::Bond => Edit::ModifyBondField {
                    id: BondHandle::New(0),
                    change: BondFieldChange::Order {
                        old: NumForm::Lit(1),
                        new: NumForm::Lit(2),
                    },
                },
                EntityKind::DativeBond => Edit::ModifyDativeBondField {
                    id: DativeBondHandle::New(0),
                    change: DativeBondFieldChange::Order {
                        old: NumForm::Lit(1),
                        new: NumForm::Lit(2),
                    },
                },
                EntityKind::AromaticSystem => Edit::ModifyAromaticSystemField {
                    id: AromaticSystemHandle::New(0),
                    change: AromaticSystemFieldChange::Charge {
                        old: NumForm::default(),
                        new: NumForm::Lit(1),
                    },
                },
                EntityKind::MulticenterBond => Edit::ModifyMulticenterBondField {
                    id: MulticenterBondHandle::New(0),
                    change: MulticenterBondFieldChange::Charge {
                        old: NumForm::default(),
                        new: NumForm::Lit(-1),
                    },
                },
                EntityKind::NoncovalentBond => Edit::ModifyNoncovalentBondField {
                    id: NoncovalentBondHandle::New(0),
                    change: NoncovalentBondFieldChange::Kind {
                        old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
                        new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
                    },
                },
                EntityKind::StereoAtom => Edit::ModifyStereoAtomField {
                    id: StereoAtomHandle::New(0),
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
                },
                EntityKind::StereoBond => Edit::ModifyStereoBondField {
                    id: StereoBondHandle::New(0),
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
                },
            });
        }

        let mut editor = base.edit();
        editor.apply(edits).unwrap();

        if kinds.contains(&EntityKind::Atom) {
            prop_assert_eq!(&editor.atom(AtomId(4)).attributes().charge, &NumForm::Lit(1));
        }
        if kinds.contains(&EntityKind::Bond) {
            prop_assert_eq!(&editor.bond(BondId(1)).attributes().order, &NumForm::Lit(2));
        }
        if kinds.contains(&EntityKind::DativeBond) {
            prop_assert_eq!(
                &editor.dative_bond(DativeBondId(0)).attributes().order,
                &NumForm::Lit(2),
            );
        }
        if kinds.contains(&EntityKind::AromaticSystem) {
            prop_assert_eq!(
                &editor.aromatic_system(AromaticSystemId(0)).attributes().charge,
                &NumForm::Lit(1),
            );
        }
        if kinds.contains(&EntityKind::MulticenterBond) {
            prop_assert_eq!(
                &editor.multicenter_bond(MulticenterBondId(0)).attributes().charge,
                &NumForm::Lit(-1),
            );
        }
        if kinds.contains(&EntityKind::NoncovalentBond) {
            prop_assert_eq!(
                &editor.noncovalent_bond(NoncovalentBondId(0)).attributes().kind,
                &NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
            );
        }
        if kinds.contains(&EntityKind::StereoAtom) {
            prop_assert_eq!(
                &editor.stereo_atom(StereoAtomId(0)).attributes().configuration,
                &StereoConfigurationForm::kinded(
                    StereoKind::Tetrahedral,
                    StereoCoset::Lit(0),
                ),
            );
        }
        if kinds.contains(&EntityKind::StereoBond) {
            prop_assert_eq!(
                &editor.stereo_bond(StereoBondId(0)).attributes().configuration,
                &StereoConfigurationForm::kinded(
                    StereoKind::CisTrans,
                    StereoCoset::Lit(0),
                ),
            );
        }
    }

    /// Every batched operation resolves all entries before mutation: an invalid handle at any
    /// generated position leaves the molecule equal to its transaction-entry state.
    #[test]
    fn test_molecule_transact_batch_error(
        batch in invalid_transaction_batch_strategy(),
    ) {
        let before = batch.base();
        let mut molecule = before.clone();
        let mut tracked = before.clone();
        prop_assert_eq!(
            before.clone().edit().tracked_apply(batch.edits()).err(),
            Some(MoleculeApplyError::Transaction(batch.expected_error()))
        );
        prop_assert_eq!(
            tracked.tracked_transact([batch.edits()], RollbackMethod::Journal).unwrap_err(),
            MoleculeApplyError::Transaction(batch.expected_error())
        );
        prop_assert_eq!(tracked, before.clone());
        prop_assert_eq!(
            molecule.transact([batch.edits()], RollbackMethod::Journal).unwrap_err(),
            MoleculeApplyError::Transaction(batch.expected_error())
        );
        prop_assert_eq!(molecule, before);
    }

    /// `lift_constraints` followed by `inline_constraints` is idempotent:
    /// running the pair twice yields the same `Molecule` as running it
    /// once. Duplicate (entity, kind) entries across the inline + molecule
    /// scopes combine by meet on the first pass, so the second pass is a
    /// fixpoint; a meet to `⊥` fails the first pass and never the second.
    #[test]
    fn test_lift_inline_idempotent(molecule in molecule_with_constraints_strategy()) {
        let Ok(once) = molecule.lift_constraints().inline_constraints() else {
            return Ok(());
        };
        let twice = once.clone().lift_constraints().inline_constraints().unwrap();

        prop_assert_eq!(once, twice);
    }

    /// `lift_constraints` drains every entity's inline `constraints` store.
    #[test]
    fn test_lift_drains_entity_stores(molecule in molecule_with_constraints_strategy()) {
        let a = molecule.lift_constraints();
        for view in a.atoms().iter() {
            prop_assert!(view.attributes().constraints.is_empty());
        }
        for view in a.bonds().iter() {
            prop_assert!(view.attributes().constraints.is_empty());
        }
        for view in a.dative_bonds().iter() {
            prop_assert!(view.attributes().constraints.is_empty());
        }
        for view in a.aromatic_systems().iter() {
            prop_assert!(view.attributes().constraints.is_empty());
        }
        for view in a.multicenter_bonds().iter() {
            prop_assert!(view.attributes().constraints.is_empty());
        }
        for view in a.noncovalent_bonds().iter() {
            prop_assert!(view.attributes().constraints.is_empty());
        }
    }

    #[test]
    fn test_transaction_rollback(
        (base, edits) in transaction_edits_strategy(),
        rollback_method in prop::sample::select(vec![RollbackMethod::Clone, RollbackMethod::Journal]),
    ) {
        let mut molecule = base.clone();
        Transaction::run(
            &mut molecule,
            rollback_method,
            |mut transaction| -> Result<(), MoleculeApplyError> {
                transaction.apply(edits)?;
                transaction.rollback();
                Ok(())
            },
        )
        .unwrap();
        prop_assert_eq!(molecule, base);
    }

    /// A valid prefix followed by a rejected edit reports the primary failure and restores the
    /// initial state without replacing the primary error.
    #[test]
    fn test_molecule_transact_error(
        (base, prefix) in transaction_edits_strategy(),
    ) {
        let atom_count = base.atoms().count();
        let edits = prefix
            .into_iter()
            .chain([Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(atom_count as u32)),
                change: AtomFieldChange::Charge {
                    old: NumForm::default(),
                    new: NumForm::Lit(1),
                },
            }])
            .collect();
        let expected = TransactionError::HandleOutOfRange {
            kind: EntityKind::Atom,
            index: atom_count,
            count: atom_count,
        };
        let mut molecule = base.clone();

        let error = molecule.transact([edits], RollbackMethod::Journal).unwrap_err();
        prop_assert_eq!(error, MoleculeApplyError::Transaction(expected));
        prop_assert_eq!(molecule, base);
    }

    /// Removing the first entity of any kind drops constraints on that entity, compacts every
    /// surviving reference, and preserves the exact order and multiplicity of duplicate entries.
    #[test]
    fn test_transaction_apply_constraint_compaction(
        case in constraint_compaction_case_strategy(),
        rollback_method in prop::sample::select(vec![RollbackMethod::Clone, RollbackMethod::Journal]),
    ) {
        let base = case.base();
        let mut editor = base.clone().edit();
        editor.apply(case.edits()).unwrap();
        let constraints = editor.constraints().iter().cloned().collect::<Vec<_>>();
        prop_assert_eq!(constraints.as_slice(), case.expected());

        let mut molecule = base.clone();
        Transaction::run(
            &mut molecule,
            rollback_method,
            |mut transaction| -> Result<(), MoleculeApplyError> {
                transaction.apply(case.edits())?;
                transaction.rollback();
                Ok(())
            },
        )
        .unwrap();
        prop_assert_eq!(molecule, base);
    }

    #[test]
    fn test_transaction_rollback_unpaired_electrons(
        atom_components in partial_unpaired_electrons_update_strategy(),
        bond_components in partial_unpaired_electrons_update_strategy(),
        aromatic_components in partial_unpaired_electrons_update_strategy(),
        multicenter_components in partial_unpaired_electrons_update_strategy(),
        rollback_method in prop::sample::select(vec![RollbackMethod::Clone, RollbackMethod::Journal]),
    ) {
        let atom = AtomForm::from_element(Element::C).with_unpaired_electrons((2_u8, 3_u8));
        let bond = BondForm::from_order(1).with_unpaired_electrons((2_u8, 3_u8));
        let aromatic =
            AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((2_u8, 3_u8));
        let multicenter =
            MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_unpaired_electrons((2_u8, 3_u8));
        let atom_update = AtomUpdate {
            unpaired_electrons: atom_components,
            ..Default::default()
        };
        let bond_update = BondUpdate {
            unpaired_electrons: bond_components,
            ..Default::default()
        };
        let aromatic_update = AromaticSystemUpdate {
            unpaired_electrons: aromatic_components,
            ..Default::default()
        };
        let multicenter_update = MulticenterBondUpdate {
            unpaired_electrons: multicenter_components,
            ..Default::default()
        };
        let base = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                atom.clone(),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![(AtomId(0), AtomId(1), bond.clone())],
            aromatic: vec![(vec![AtomId(0), AtomId(1), AtomId(2)], aromatic.clone())],
            multicenter: vec![(vec![AtomId(0), AtomId(1), AtomId(2)], multicenter.clone())],
            ..Default::default()
        });
        let expected = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                atom.update(&atom_update),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![(AtomId(0), AtomId(1), bond.update(&bond_update))],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                aromatic.update(&aromatic_update),
            )],
            multicenter: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                multicenter.update(&multicenter_update),
            )],
            ..Default::default()
        });
        let mut edits = Edits::new();
        edits.update_atom(AtomHandle::Id(AtomId(0)), &atom, &atom_update);
        edits.update_bond(BondHandle::Id(BondId(0)), &bond, &bond_update);
        edits.update_aromatic_system(
            AromaticSystemHandle::Id(AromaticSystemId(0)),
            &aromatic,
            &aromatic_update,
        );
        edits.update_multicenter_bond(
            MulticenterBondHandle::Id(MulticenterBondId(0)),
            &multicenter,
            &multicenter_update,
        );

        let mut molecule = base.clone();
        let applied = Transaction::run(
            &mut molecule,
            rollback_method,
            |mut transaction| -> Result<_, MoleculeApplyError> {
                transaction.apply(edits)?;
                let applied = transaction.probe()?.clone();
                transaction.rollback();
                Ok(applied)
            },
        )
        .unwrap();
        prop_assert_eq!(applied, expected);
        prop_assert_eq!(molecule, base);
    }

    #[test]
    fn test_transaction_rollback_batches(
        (base, first_edits, second_edits) in consecutive_transaction_strategy(),
        rollback_method in prop::sample::select(vec![RollbackMethod::Clone, RollbackMethod::Journal]),
    ) {
        let mut molecule = base.clone();
        Transaction::run(
            &mut molecule,
            rollback_method,
            |mut transaction| -> Result<(), MoleculeApplyError> {
                transaction.apply(first_edits)?;
                transaction.apply(second_edits)?;
                transaction.rollback();
                Ok(())
            },
        )
        .unwrap();
        prop_assert_eq!(molecule, base);
    }

    /// Scoped application agrees with journaled execution on generated valid edit batches.
    #[test]
    fn test_molecule_apply_into((base, edits) in transaction_edits_strategy()) {
        let mut expected = base.clone();
        expected.transact([edits.clone()], RollbackMethod::Journal).unwrap();
        let mut molecule = base;
        prop_assert_eq!(molecule.apply_into(edits), Ok(()));
        prop_assert_eq!(molecule, expected);
    }

    /// One callback containing consecutive batches and repeated apply_into calls have the same
    /// result as a transaction containing those batches, including their separate handle spaces.
    #[test]
    fn test_molecule_edit_into(
        (base, first, second) in consecutive_transaction_strategy(),
    ) {
        let mut expected = base.clone();
        expected.transact([first.clone(), second.clone()], RollbackMethod::Journal).unwrap();
        let mut molecule = base.clone();
        let result = molecule.edit_into(|editor| -> Result<_, MoleculeApplyError> {
            editor.apply(first.clone())?;
            editor.apply(second.clone())?;
            Ok(editor.atom_count())
        });
        prop_assert_eq!(result, Ok(expected.atoms().count()));
        prop_assert_eq!(&molecule, &expected);
        let mut repeated = base;
        repeated.apply_into(first).unwrap();
        repeated.apply_into(second).unwrap();
        prop_assert_eq!(repeated, expected);
    }

    /// A failing edit after a valid generated prefix preserves the complete entry representation.
    /// This checks recovery, separately from the successful journal/application agreement above.
    #[test]
    fn test_molecule_apply_into_error((base, mut edits) in transaction_edits_strategy()) {
        let count = base.atoms().count();
        edits.remove_atom(AtomHandle::Id(AtomId(count as u32)));
        for scoped in [false, true] {
            let mut molecule = base.clone();
            let result = if scoped {
                molecule.edit_into(|editor| editor.apply(edits.clone()))
            } else {
                molecule.apply_into(edits.clone())
            };
            prop_assert_eq!(result, Err(MoleculeApplyError::Transaction(TransactionError::HandleOutOfRange {
                kind: EntityKind::Atom,
                index: count,
                count,
            })));
            prop_assert_eq!(molecule, base.clone());
        }
    }

    #[test]
    fn test_molecule_editor_apply((base, edits) in transaction_edits_strategy()) {
        let mut checked = base.clone();
        checked.transact([edits.clone()], RollbackMethod::Journal).unwrap();
        let mut applied = base.edit();
        applied.apply(edits).unwrap();
        prop_assert_eq!(applied.finish().unwrap(), checked);
    }

    /// A rejected edit leaves the same completed prefix as successful transaction execution.
    /// A later batch can still be applied and published from that state.
    #[test]
    fn test_molecule_editor_apply_partial((base, mut edits) in transaction_edits_strategy()) {
        let mut expected = base.clone();
        expected.transact([edits.clone()], RollbackMethod::Journal).unwrap();
        let count = base.atoms().count();
        edits.remove_atom(AtomHandle::Id(AtomId(count as u32)));
        for tracked in [false, true] {
            let mut editor = base.clone().edit();
            let error = if tracked {
                editor.tracked_apply(edits.clone()).unwrap_err()
            } else {
                editor.apply(edits.clone()).unwrap_err()
            };
            prop_assert_eq!(error, MoleculeApplyError::Transaction(TransactionError::HandleOutOfRange {
                kind: EntityKind::Atom,
                index: count,
                count,
            }));
            prop_assert_eq!(editor.probe().unwrap(), &expected);
            editor.apply(Edits::new()).unwrap();
            prop_assert_eq!(editor.finish().unwrap(), expected.clone());
        }
    }

    /// A rejected edit after a valid prefix preserves the primary error and restores the entry
    /// state identically through both scoped handles.
    #[test]
    fn test_tracked_transaction_run_error(
        (base, prefix) in transaction_edits_strategy(),
        rollback_method in prop::sample::select(vec![RollbackMethod::Clone, RollbackMethod::Journal]),
    ) {
        let atom_count = base.atoms().count();
        let edits: Edits = prefix
            .into_iter()
            .chain([Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(atom_count as u32)),
                change: AtomFieldChange::Charge {
                    old: NumForm::default(),
                    new: NumForm::Lit(1),
                },
            }])
            .collect();
        let expected = MoleculeApplyError::Transaction(TransactionError::HandleOutOfRange {
            kind: EntityKind::Atom,
            index: atom_count,
            count: atom_count,
        });
        let mut ordinary = base.clone();
        let plain = Transaction::run(&mut ordinary, rollback_method, |mut transaction| {
            transaction.apply(edits.clone())?;
            transaction.commit()
        });
        let mut tracked = base.clone();
        let result = TrackedTransaction::run(&mut tracked, rollback_method, |mut transaction| {
            transaction.apply(edits)?;
            transaction.commit()
        });
        prop_assert_eq!(result.map(|_| ()), Err(expected.clone()));
        prop_assert_eq!(plain, Err(expected));
        prop_assert_eq!(tracked, base.clone());
        prop_assert_eq!(ordinary, base);
    }

    /// Discarding a tracked commit's correspondence gives the ordinary outcome; the mapping
    /// and resulting molecule also agree with consuming application of the same edits.
    #[test]
    fn test_tracked_transaction_commit(
        (base, edits) in transaction_edits_strategy(),
        rollback_method in prop::sample::select(vec![RollbackMethod::Clone, RollbackMethod::Journal]),
    ) {
        let expected = base.clone().tracked_apply(edits.clone()).unwrap();
        let mut ordinary = base.clone();
        Transaction::run(&mut ordinary, rollback_method, |mut transaction| {
            transaction.apply(edits.clone())?;
            transaction.commit()
        }).unwrap();
        let mut tracked = base;
        let correspondence = TrackedTransaction::run(&mut tracked, rollback_method, |mut transaction| {
            transaction.apply(edits)?;
            transaction.commit()
        }).unwrap();
        prop_assert_eq!(&tracked, &ordinary);
        prop_assert_eq!((tracked, correspondence), expected);
    }

    /// A transaction spanning consecutive batches has the composed batch correspondence and the
    /// same result as separate commits. Rolling both batches back restores the initial molecule.
    #[test]
    fn test_tracked_transaction_commit_composition(
        (base, first_edits, second_edits) in consecutive_transaction_strategy(),
        rollback_method in prop::sample::select(vec![RollbackMethod::Clone, RollbackMethod::Journal]),
    ) {
        let mut separate = base.clone();
        let first = separate.tracked_transact([first_edits.clone()], RollbackMethod::Journal).unwrap();
        let first_result = separate.clone();
        let mut applied = separate.clone().edit();
        let applied_witness = applied.tracked_apply(second_edits.clone()).unwrap();
        let second = separate.tracked_transact([second_edits.clone()], RollbackMethod::Journal).unwrap();
        prop_assert_eq!(applied_witness, second.clone());
        prop_assert_eq!(applied.finish().unwrap(), separate.clone());
        prop_assert_eq!(
            base.clone().tracked_apply(first_edits.clone()).unwrap(),
            (first_result, first.clone())
        );

        let mut together = base.clone();
        let composed = TrackedTransaction::run(&mut together, rollback_method, |mut transaction| {
            transaction.apply(first_edits.clone())?;
            transaction.apply(second_edits.clone())?;
            transaction.commit()
        }).unwrap();
        prop_assert_eq!(composed, first.compose(&second).unwrap());
        prop_assert_eq!(together, separate.clone());
        let mut plain = base.clone();
        plain
            .transact([first_edits.clone(), second_edits.clone()], RollbackMethod::Journal)
            .unwrap();
        prop_assert_eq!(plain, separate);

        let mut restored = base.clone();
        TrackedTransaction::run(
            &mut restored,
            rollback_method,
            |mut transaction| -> Result<(), MoleculeApplyError> {
                transaction.apply(first_edits)?;
                transaction.apply(second_edits)?;
                transaction.rollback();
                Ok(())
            },
        )
        .unwrap();
        prop_assert_eq!(restored, base);
    }

    /// Descending atom removals and subsequent additions name the same entities under any
    /// batch split. Accumulation and composition of separate commits agree, including cascades.
    #[test]
    fn test_tracked_transaction_commit_partition(
        (base, removed) in molecule_with_atom_subset_strategy(),
        added in prop::collection::vec(atom_form_strategy(), 0..=6),
        batch_size in 1usize..=4,
        rollback_method in prop::sample::select(vec![RollbackMethod::Clone, RollbackMethod::Journal]),
    ) {
        let mut edits = Edits::new();
        for id in removed.into_iter().rev() {
            edits.remove_atom(AtomHandle::Id(id));
        }
        for attributes in added {
            edits.add_atom(attributes);
        }
        let raw: Vec<_> = edits.iter().cloned().collect();
        let mut single = base.clone();
        let expected = single.tracked_transact([edits], rollback_method).unwrap();
        let mut separate = base.clone();
        let mut composed = separate.tracked_transact([], RollbackMethod::Journal).unwrap();
        for chunk in raw.chunks(batch_size) {
            let next = separate
                .tracked_transact([Edits::from_iter(chunk.iter().cloned())], RollbackMethod::Journal)
                .unwrap();
            composed = composed.compose(&next).unwrap();
        }
        let mut partitioned = base;
        let actual = TrackedTransaction::run(&mut partitioned, rollback_method, |mut transaction| {
            for chunk in raw.chunks(batch_size) {
                transaction.apply(Edits::from_iter(chunk.iter().cloned()))?;
            }
            transaction.commit()
        }).unwrap();
        prop_assert_eq!(&actual, &expected);
        prop_assert_eq!(composed, expected);
        prop_assert_eq!(&partitioned, &separate);
        prop_assert_eq!(partitioned, single);
    }

    /// Explicit cancellation restores the complete entry representation after a generated batch.
    #[test]
    fn test_tracked_transaction_rollback(
        (base, edits) in transaction_edits_strategy(),
        rollback_method in prop::sample::select(vec![RollbackMethod::Clone, RollbackMethod::Journal]),
    ) {
        let mut molecule = base.clone();
        TrackedTransaction::run(
            &mut molecule,
            rollback_method,
            |mut transaction| -> Result<(), MoleculeApplyError> {
                transaction.apply(edits)?;
                transaction.rollback();
                Ok(())
            },
        ).unwrap();
        prop_assert_eq!(molecule, base);
    }

    /// `inline_constraints` removes every TOP-LEVEL inline-capable narrow
    /// leaf from the molecule list; a meet to `⊥` fails. Combinator-nested entries, relational leaves,
    /// molecule-scope leaves are preserved.
    #[test]
    fn test_inline_removes_top_level_leaves(molecule in molecule_with_constraints_strategy()) {
        let Ok(a) = molecule.inline_constraints() else {
            return Ok(());
        };
        for c in a.constraints().iter() {
            prop_assert!(
                !matches!(
                    c,
                    Constraint::Atom(..)
                        | Constraint::Bond(..)
                        | Constraint::DativeBond(..)
                        | Constraint::AromaticSystem(..)
                        | Constraint::MulticenterBond(..)
                        | Constraint::NoncovalentBond(..)
                ),
                "inline-capable narrow leaf survived inline_constraints: {c:?}",
            );
        }
    }

    /// `inline_constraints` deposits each top-level narrow leaf into the
    /// targeted entity's inline `constraints` store, indexed by the leaf's
    /// `kind()`. Colliding same-`(id, kind)` entries combine by meet, so on
    /// success the kind is still present after the call; a meet to `⊥`
    /// fails and leaves the molecule unchanged.
    #[test]
    fn test_inline_deposits_leaves_into_entities(
        molecule in molecule_with_constraints_strategy(),
    ) {
        let mut atom_keys: HashSet<(AtomId, AtomConstraintKey)> = HashSet::new();
        let mut bond_keys: HashSet<(BondId, BondConstraintKey)> = HashSet::new();
        let mut dative_keys: HashSet<(DativeBondId, DativeBondConstraintKey)> = HashSet::new();
        let mut aromatic_keys: HashSet<(AromaticSystemId, AromaticSystemConstraintKey)> =
            HashSet::new();
        let mut multicenter_keys: HashSet<(MulticenterBondId, MulticenterBondConstraintKey)> =
            HashSet::new();
        for c in molecule.constraints().iter() {
            match c {
                Constraint::Atom(id, inner) => {
                    atom_keys.insert((*id, inner.key()));
                }
                Constraint::Bond(id, inner) => {
                    bond_keys.insert((*id, inner.key()));
                }
                Constraint::DativeBond(id, inner) => {
                    dative_keys.insert((*id, inner.key()));
                }
                Constraint::AromaticSystem(id, inner) => {
                    aromatic_keys.insert((*id, inner.key()));
                }
                Constraint::MulticenterBond(id, inner) => {
                    multicenter_keys.insert((*id, inner.key()));
                }
                _ => {}
            }
        }

        let Ok(a) = molecule.inline_constraints() else {
            return Ok(());
        };

        for (id, key) in atom_keys {
            prop_assert!(
                a.atom(id).attributes().constraints.contains(key),
                "atom {id:?} missing key {key:?} after inline",
            );
        }
        for (id, key) in bond_keys {
            prop_assert!(
                a.bond(id).attributes().constraints.contains(key),
                "bond {id:?} missing key {key:?} after inline",
            );
        }
        for (id, key) in dative_keys {
            prop_assert!(
                a.dative_bond(id).attributes().constraints.contains(key),
                "dative bond {id:?} missing key {key:?} after inline",
            );
        }
        for (id, key) in aromatic_keys {
            prop_assert!(
                a.aromatic_system(id).attributes().constraints.contains(key),
                "aromatic system {id:?} missing key {key:?} after inline",
            );
        }
        for (id, key) in multicenter_keys {
            prop_assert!(
                a.multicenter_bond(id).attributes().constraints.contains(key),
                "multicenter bond {id:?} missing key {key:?} after inline",
            );
        }
    }
}

#[fixture]
fn source() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::default(); 2],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::default())],
        ..Default::default()
    })
}

#[rstest]
fn test_molecule_apply_error(source: Molecule) {
    let edits: Edits = [
        Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(1),
            },
        },
        Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(2)),
            change: AtomFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(1),
            },
        },
    ]
    .into_iter()
    .collect();
    let expected = TransactionError::HandleOutOfRange {
        kind: EntityKind::Atom,
        index: 2,
        count: 2,
    };

    assert_eq!(
        source.clone().apply(edits.clone()),
        Err(MoleculeApplyError::Transaction(expected.clone()))
    );
    let mut molecule = source.clone();
    assert_eq!(
        molecule.transact([edits], RollbackMethod::Journal),
        Err(expected.into())
    );
    assert_eq!(molecule, source);
}

#[rstest]
fn test_molecule_editor_finish_error(source: Molecule) {
    let mut editor = source.edit();
    editor.add_bond(AtomId(0), AtomId(1), BondForm::default());

    assert_eq!(
        editor.finish(),
        Err(MoleculeIntegrityError::ParallelBonds {
            atoms: [AtomId(0), AtomId(1)],
        })
    );
}
