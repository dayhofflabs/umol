//! Source-relative lowering agrees with a direct raw update under normalized_eq.
//! Sources and results must normalize; raw kindless payloads remain in the domain.
//! Direct application and journaled application are checked against the same result.
//! Ordered batches also preserve raw rows and agree with independent entry updates;
//! rollback restores the source and DSL export omits only kindless constraint batches.

use proptest::prelude::*;
use umol_edn::{FromEdn, ToEdn};
use umol_graph_ir::dsl::{EditsDsl, MoleculeDefaults};
use umol_graph_ir::ir::{
    AtomForm, AtomId, BondForm, BondId, Edit, Edits, FromIr, IntoIr, Molecule, MoleculeApplyError,
    MoleculeEntries, Normalize, RollbackMethod, StereoAtomConstraintForm, StereoAtomFieldChange,
    StereoAtomForm, StereoAtomHandle, StereoAtomId, StereoBondHandle, StereoBondId, StereoLigand,
    StereoLigandKind, Transaction,
};

use crate::strategies::{
    raw_stereo_atom_form_strategy, raw_stereo_atom_update_strategy, raw_stereo_bond_form_strategy,
    raw_stereo_bond_update_strategy, raw_stereo_configuration_strategy,
    stereo_atom_constraint_strategy, stereo_atom_kind_strategy,
};

proptest! {

    /// The independent entry table, singles, and batches agree exactly; rollback
    /// restores raw configurations and the original site/ligand order.
    #[test]
    fn test_molecule_apply_stereo_atom_configurations(
        values in prop::collection::vec(raw_stereo_configuration_strategy(), 3),
        updates in prop::collection::vec((0usize..3, raw_stereo_configuration_strategy()), 0..16),
    ) {
        let mut entries = MoleculeEntries { atoms: vec![AtomForm::default(); 3], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1)), (AtomId(1), AtomId(2), BondForm::from_order(1)), (AtomId(2), AtomId(0), BondForm::from_order(1))], stereo_atoms: values.into_iter().enumerate().map(|(i, configuration)| (AtomId::from(i), vec![StereoLigand::new(AtomId::from((i + 2) % 3), StereoLigandKind::Atom), StereoLigand::new(AtomId::from(i), StereoLigandKind::LonePair)], StereoAtomForm { configuration, ..Default::default() })).collect(), ..Default::default() };
        let source = Molecule::from_entries(entries.clone());
        let mut rows = Vec::new();
        let mut singles = Edits::new();
        for (index, new) in updates {
            let old = entries.stereo_atoms[index].2.configuration.clone();
            entries.stereo_atoms[index].2.configuration = new.clone();
            let id = StereoAtomHandle::Id(StereoAtomId::from(index));
            rows.push((id.clone(), old.clone(), new.clone()));
            singles.push(Edit::ModifyStereoAtomField { id, change: StereoAtomFieldChange::Configuration { old, new } });
        }
        let expected = Molecule::from_entries(entries);
        let mut batch = Edits::new();
        batch.modify_stereo_atom_configurations(rows);
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
        let defaults = MoleculeDefaults::new();
        let dsl = EditsDsl::from_ir(&batch, &defaults);
        prop_assert_eq!(EditsDsl::from_edn(&dsl.to_edn()).unwrap().into_ir(&defaults), batch);
    }

    /// All keys retain their parameters through ordered mutation and restoration.
    /// A known-kind batch roundtrips exactly; a kindless batch is wholly omitted.
    #[test]
    fn test_molecule_apply_stereo_atom_constraints(
        (kind, value) in stereo_atom_kind_strategy().prop_flat_map(|kind| stereo_atom_constraint_strategy(kind).prop_map(move |value| (kind, value))),
        kindless in any::<bool>(),
        updates in prop::collection::vec((0usize..3, 0u8..3), 0..20),
    ) {
        let kind = if kindless { None } else { Some(kind) };
        let key = value.key();
        let mut entries = MoleculeEntries { atoms: vec![AtomForm::default(); 3], bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1)), (AtomId(1), AtomId(2), BondForm::from_order(1)), (AtomId(2), AtomId(0), BondForm::from_order(1))], stereo_atoms: (0..3).map(|i| (AtomId(i), vec![StereoLigand::new(AtomId((i + 2) % 3), StereoLigandKind::Atom)], StereoAtomForm::default())).collect(), ..Default::default() };
        let source = Molecule::from_entries(entries.clone());
        let mut rows = Vec::new();
        let mut singles = Edits::new();
        for (index, choice) in updates {
            let old = entries.stereo_atoms[index].2.constraints.get(key).cloned();
            let new = match choice { 0 => None, 1 => Some(value.clone()), _ => Some(StereoAtomConstraintForm::undetermined(key)) };
            entries.stereo_atoms[index].2.constraints = new.clone().into_iter().collect();
            let id = StereoAtomHandle::Id(StereoAtomId::from(index));
            rows.push((id.clone(), old.clone(), new.clone()));
            singles.push(Edit::ModifyStereoAtomConstraint { id, kind, old, new });
        }
        let expected = Molecule::from_entries(entries);
        let mut batch = Edits::new();
        batch.modify_stereo_atom_constraints(kind, key, rows);
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
        let defaults = MoleculeDefaults::new();
        let dsl = EditsDsl::from_ir(&batch, &defaults);
        let expected_export = if kindless { Edits::new() } else { batch };
        prop_assert_eq!(EditsDsl::from_edn(&dsl.to_edn()).unwrap().into_ir(&defaults), expected_export);
    }

    #[test]
    fn test_edits_update_stereo_atom(
        current in raw_stereo_atom_form_strategy(),
        update in raw_stereo_atom_update_strategy(),
    ) {
        let expected = current.update(&update);
        if current.clone().normalize().is_ok() && expected.clone().normalize().is_ok() {
            let source = Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::default(); 6],
                bonds: [(0, 1), (0, 2), (0, 3), (0, 4), (1, 4), (1, 5)]
                    .into_iter().map(|(a, b)| (AtomId(a), AtomId(b), BondForm::from_order(1))).collect(),
                stereo_atoms: vec![(AtomId(0), (1..=4).map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).collect(), current.clone())],
                ..Default::default()
            });
            let mut edits = Edits::new();
            edits.update_stereo_atom(StereoAtomHandle::Id(StereoAtomId(0)), &current, &update);
            for edit in edits.iter() {
                if let Edit::ModifyStereoAtomConstraint { kind, old, new, .. } = edit {
                    prop_assert_eq!(*kind, expected.configuration.kind());
                    let key = old.as_ref().or(new.as_ref()).unwrap().key();
                    prop_assert_eq!(old.as_ref(), current.constraints.get(key));
                }
            }
            let direct = source.clone().apply(edits.clone()).unwrap();
            prop_assert!(direct.stereo_atom(StereoAtomId(0)).attributes().normalized_eq(&expected));
            let mut journaled = source.clone();
            Transaction::run(&mut journaled, RollbackMethod::Journal, |mut transaction| {
                transaction.apply(edits.clone())?;
                transaction.commit()
            }).unwrap();
            prop_assert!(journaled.stereo_atom(StereoAtomId(0)).attributes().normalized_eq(&expected));
            prop_assert_eq!(journaled, direct);
        }
    }

    #[test]
    fn test_edits_update_stereo_bond(
        current in raw_stereo_bond_form_strategy(),
        update in raw_stereo_bond_update_strategy(),
    ) {
        let expected = current.update(&update);
        if current.clone().normalize().is_ok() && expected.clone().normalize().is_ok() {
            let source = Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::default(); 6],
                bonds: [(0, 1), (0, 2), (0, 3), (0, 4), (1, 4), (1, 5)]
                    .into_iter().map(|(a, b)| (AtomId(a), AtomId(b), BondForm::from_order(1))).collect(),
                stereo_bonds: vec![(BondId(0), (2..=5).map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).collect(), current.clone())],
                ..Default::default()
            });
            let mut edits = Edits::new();
            edits.update_stereo_bond(StereoBondHandle::Id(StereoBondId(0)), &current, &update);
            for edit in edits.iter() {
                if let Edit::ModifyStereoBondConstraint { kind, old, new, .. } = edit {
                    prop_assert_eq!(*kind, expected.configuration.kind());
                    let key = old.as_ref().or(new.as_ref()).unwrap().key();
                    prop_assert_eq!(old.as_ref(), current.constraints.get(key));
                }
            }
            let direct = source.clone().apply(edits.clone()).unwrap();
            prop_assert!(direct.stereo_bond(StereoBondId(0)).attributes().normalized_eq(&expected));
            let mut journaled = source.clone();
            Transaction::run(&mut journaled, RollbackMethod::Journal, |mut transaction| {
                transaction.apply(edits.clone())?;
                transaction.commit()
            }).unwrap();
            prop_assert!(journaled.stereo_bond(StereoBondId(0)).attributes().normalized_eq(&expected));
            prop_assert_eq!(journaled, direct);
        }
    }

}
