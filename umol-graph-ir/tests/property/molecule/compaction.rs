//! Molecule compaction and mixed witness-composition properties.
//!
//! The mixed sequence starts with a dense renumbering, removes an atom with a compaction, and
//! adds an atom through a batch correspondence. Composition must remain compatible
//! with every intermediate molecule and with the final source-to-result pair.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_chem::element::Element;
use umol_graph_core::{Compaction, GraphCompaction};
use umol_graph_ir::ir::{
    AromaticSystemId, AtomForm, AtomHandle, AtomId, BondHandle, BondId, DativeBondId, Edits,
    MoleculeCompaction, MoleculeCorrespondence, MulticenterBondId, NoncovalentBondId,
    RollbackMethod, StereoAtomId, StereoBondId,
};

use crate::strategies::{molecule_dense_renumbering_strategy, molecule_with_removals_strategy};

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_molecule_compaction_undo_compaction(
        (molecule, atoms, bonds) in molecule_with_removals_strategy(),
    ) {
        let counts = (
            molecule.atoms().count(),
            molecule.bonds().count(),
            molecule.dative_bonds().count(),
            molecule.aromatic_systems().count(),
            molecule.multicenter_bonds().count(),
            molecule.noncovalent_bonds().count(),
            molecule.stereo_atoms().count(),
            molecule.stereo_bonds().count(),
        );
        let mut editor = molecule.clone().edit();
        let mut edits = Edits::new();
        edits.remove_topology_cascading(
            atoms.iter().copied().map(AtomHandle::Id).collect(),
            bonds.iter().copied().map(BondHandle::Id).collect(),
        );
        let correspondence = editor.tracked_apply(edits).unwrap();
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(
                Compaction::new(counts.0, correspondence.atoms().left_unmatched().into_iter().map(Into::into).collect()).unwrap(),
                Compaction::new(counts.1, correspondence.bonds().left_unmatched().into_iter().map(Into::into).collect()).unwrap(),
            ),
            Compaction::new(counts.2, correspondence.dative_bonds().left_unmatched()).unwrap(),
            Compaction::new(counts.3, correspondence.aromatic_systems().left_unmatched()).unwrap(),
            Compaction::new(counts.4, correspondence.multicenter_bonds().left_unmatched()).unwrap(),
            Compaction::new(counts.5, correspondence.noncovalent_bonds().left_unmatched()).unwrap(),
            Compaction::new(counts.6, correspondence.stereo_atoms().left_unmatched()).unwrap(),
            Compaction::new(counts.7, correspondence.stereo_bonds().left_unmatched()).unwrap(),
        );
        prop_assert_eq!(correspondence, MoleculeCorrespondence::from(&compaction));
        let mut plain = molecule.edit();
        plain.remove_topology_cascading(&atoms, &bonds);
        let publication = editor.finish();
        let expected = plain.finish();
        prop_assert_eq!(publication, expected);
        let undo = compaction.undo_compaction();

        for index in 0..counts.0 {
            let original = AtomId(index as u32);
            if let Some(compacted) = compaction.compact_atom(original) {
                prop_assert_eq!(undo.uncompact_atom(compacted), original);
            }
        }
        for index in 0..counts.1 {
            let original = BondId(index as u32);
            if let Some(compacted) = compaction.compact_bond(original) {
                prop_assert_eq!(undo.uncompact_bond(compacted), original);
            }
        }
        for index in 0..counts.2 {
            let original = DativeBondId(index as u32);
            if let Some(compacted) = compaction.compact_dative_bond(original) {
                prop_assert_eq!(undo.uncompact_dative_bond(compacted), original);
            }
        }
        for index in 0..counts.3 {
            let original = AromaticSystemId(index as u32);
            if let Some(compacted) = compaction.compact_aromatic_system(original) {
                prop_assert_eq!(undo.uncompact_aromatic_system(compacted), original);
            }
        }
        for index in 0..counts.4 {
            let original = MulticenterBondId(index as u32);
            if let Some(compacted) = compaction.compact_multicenter_bond(original) {
                prop_assert_eq!(undo.uncompact_multicenter_bond(compacted), original);
            }
        }
        for index in 0..counts.5 {
            let original = NoncovalentBondId(index as u32);
            if let Some(compacted) = compaction.compact_noncovalent_bond(original) {
                prop_assert_eq!(undo.uncompact_noncovalent_bond(compacted), original);
            }
        }
        for index in 0..counts.6 {
            let original = StereoAtomId(index as u32);
            if let Some(compacted) = compaction.compact_stereo_atom(original) {
                prop_assert_eq!(undo.uncompact_stereo_atom(compacted), original);
            }
        }
        for index in 0..counts.7 {
            let original = StereoBondId(index as u32);
            if let Some(compacted) = compaction.compact_stereo_bond(original) {
                prop_assert_eq!(undo.uncompact_stereo_bond(compacted), original);
            }
        }
    }
    /// Selecting dependent rows from the source and removing them explicitly gives the same
    /// result as the stored-reference cascade, including participant frames and constraints.
    #[test]
    fn test_molecule_editor_remove_topology_cascading(
        (source, atoms, bonds) in molecule_with_removals_strategy(),
    ) {
        let removed_bonds = source.bonds().iter().filter(|view| {
            bonds.contains(&view.id()) || view.atom_ids().iter().any(|id| atoms.contains(id))
        }).map(|view| view.id()).collect::<Vec<_>>();
        let dative = source.dative_bonds().iter()
            .filter(|view| view.atom_ids().any(|id| atoms.contains(&id)))
            .map(|view| view.id()).collect::<Vec<_>>();
        let aromatic = source.aromatic_systems().iter()
            .filter(|view| view.atom_ids().any(|id| atoms.contains(&id)))
            .map(|view| view.id()).collect::<Vec<_>>();
        let multicenter = source.multicenter_bonds().iter()
            .filter(|view| view.atom_ids().any(|id| atoms.contains(&id)))
            .map(|view| view.id()).collect::<Vec<_>>();
        let noncovalent = source.noncovalent_bonds().iter()
            .filter(|view| view.atom_ids().iter().any(|id| atoms.contains(id)))
            .map(|view| view.id()).collect::<Vec<_>>();
        let stereo_atoms = source.stereo_atoms().iter()
            .filter(|view| view.atom_ids().any(|id| atoms.contains(&id)))
            .map(|view| view.id()).collect::<Vec<_>>();
        let stereo_bonds = source.stereo_bonds().iter()
            .filter(|view| removed_bonds.contains(&view.site_id()) || view.ligand_ids().iter().any(|ligand| atoms.contains(&ligand.atom_id)))
            .map(|view| view.id()).collect::<Vec<_>>();

        let mut exact = source.clone().edit();
        exact.remove_dative_bonds(&dative);
        exact.remove_aromatic_systems(&aromatic);
        exact.remove_multicenter_bonds(&multicenter);
        exact.remove_noncovalent_bonds(&noncovalent);
        exact.remove_stereo_atoms(&stereo_atoms);
        exact.remove_stereo_bonds(&stereo_bonds);
        let exact_compaction = exact.tracked_remove_topology(&atoms, &removed_bonds).unwrap();
        let mut cascade = source.edit();
        let compaction = cascade.tracked_remove_topology_cascading(&atoms, &bonds);
        prop_assert_eq!(exact_compaction.graph(), compaction.graph());
        prop_assert_eq!(exact.finish(), cascade.finish());
    }

    #[test]
    fn test_molecule_editor_tracked_apply_composition(
        (molecule, atoms, bonds) in molecule_with_removals_strategy(),
    ) {
        let mut first = Edits::new();
        first.remove_topology_cascading(
            atoms.into_iter().map(AtomHandle::Id).collect(),
            bonds.into_iter().map(BondHandle::Id).collect(),
        );
        let mut addition = Edits::new();
        addition.add_atom(AtomForm::from_element(Element::F));
        let mut last = Edits::new();
        last.remove_topology_cascading(vec![AtomHandle::Id(AtomId(0))], vec![]);
        let mut expected = molecule.clone();
        let complete = expected.tracked_transact([first.clone(), addition.clone(), last.clone()], RollbackMethod::Journal).unwrap();

        let mut editor = molecule.edit();
        let first = editor.tracked_apply(first).unwrap();
        let addition = editor.tracked_apply(addition).unwrap();
        let last = editor.tracked_apply(last).unwrap();
        let composed = first.compose(&addition).unwrap().compose(&last).unwrap();
        prop_assert_eq!(composed, complete);
        prop_assert_eq!(editor.finish(), Ok(expected));
    }

    #[test]
    fn test_mixed_witness_composition(
        (source, remapping) in molecule_dense_renumbering_strategy(),
    ) {
        let remapped = source.clone().remap(&remapping);
        let removed = (remapped.atoms().count() > 0)
            .then_some(AtomId(0))
            .into_iter()
            .collect::<Vec<_>>();
        let remapping = MoleculeCorrespondence::from(&remapping);

        let mut removal_editor = remapped.clone().edit();
        let mut edits = Edits::new();
        edits.remove_topology_cascading(removed.into_iter().map(AtomHandle::Id).collect(), vec![]);
        let compaction = removal_editor.tracked_apply(edits).unwrap();
        let compacted = removal_editor.finish().unwrap();

        let mut edits = Edits::new();
        edits.add_atom(AtomForm::from_element(Element::F));
        let (result, addition) = compacted.clone().tracked_apply(edits).unwrap();

        prop_assert!(remapping.is_compatible(&source, &remapped));
        prop_assert!(compaction.is_compatible(&remapped, &compacted));
        prop_assert!(addition.is_compatible(&compacted, &result));

        let composed = MoleculeCorrespondence::compose_all([
            remapping,
            compaction,
            addition,
        ])
        .unwrap()
        .expect("the sequence contains three correspondences");
        prop_assert!(composed.is_compatible(&source, &result));
        prop_assert_eq!(
            &composed,
            &MoleculeCorrespondence::induce(&source, &result, composed.atoms().clone())
                .expect("the composed atom pairs uniquely induce the surviving entities"),
        );
    }

}
