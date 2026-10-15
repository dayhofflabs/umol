//! Source-relative lowering agrees with a direct raw update under normalized_eq.
//! Sources and results must normalize; raw kindless payloads remain in the domain.
//! Direct application and journaled application are checked against the same result.

use proptest::prelude::*;
use umol_graph_ir::ir::{
    AtomForm, AtomId, BondForm, BondId, Edit, Edits, Molecule, MoleculeEntries, Normalize,
    RollbackMethod, StereoAtomHandle, StereoAtomId, StereoBondHandle, StereoBondId, StereoLigand,
    StereoLigandKind, Transaction,
};

use crate::strategies::{
    raw_stereo_atom_form_strategy, raw_stereo_atom_update_strategy, raw_stereo_bond_form_strategy,
    raw_stereo_bond_update_strategy,
};

proptest! {
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
