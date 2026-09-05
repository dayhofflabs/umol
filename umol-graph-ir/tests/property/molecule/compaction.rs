//! Molecule-editor compaction and mixed witness-composition properties.
//!
//! The mixed sequence starts with a dense renumbering, removes an atom with a compaction, and
//! adds an atom through the editor's general correspondence. Composition must remain compatible
//! with every intermediate molecule and with the final source-to-result pair.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};

use crate::strategies::*;

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_id_compaction_undo(
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
        let mut editor = molecule.edit();
        let compaction = editor.tracked_remove(&atoms, &bonds);
        let mut plain = molecule.edit();
        plain.remove(&atoms, &bonds);
        let publication = editor.try_tracked_build();
        let expected = plain.try_build().map(|molecule| (molecule, MoleculeCorrespondence::from(&compaction)));
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
    #[test]
    fn test_molecule_editor_tracked_build_composition(
        (molecule, atoms, bonds) in molecule_with_removals_strategy(),
    ) {
        let mut editor = molecule.edit();
        let first = editor.tracked_remove(&atoms, &bonds);
        editor.add_atom(AtomForm::from_element(Element::F));
        let second = editor.tracked_remove(&[AtomId(0)], &[]);
        let expected = MoleculeCorrespondence::from(&first)
            .extend_right(EntityKind::Atom, 1)
            .compose(&MoleculeCorrespondence::from(&second)).unwrap();
        let plain = editor.clone().try_build();
        prop_assert_eq!(editor.try_tracked_build(), plain.map(|molecule| (molecule, expected)));
    }

    #[test]
    fn test_mixed_witness_composition(
        (source, remapping) in molecule_dense_renumbering_strategy(),
    ) {
        let remapped = source.remap(&remapping);
        let removed = (remapped.atoms().count() > 0)
            .then_some(AtomId(0))
            .into_iter()
            .collect::<Vec<_>>();
        let remapping = MoleculeCorrespondence::from(&remapping);

        let mut removal_editor = remapped.edit();
        let compaction = removal_editor.tracked_remove(&removed, &[]);
        let compacted = removal_editor.build();
        let compaction = MoleculeCorrespondence::from(&compaction);

        let mut addition_editor = compacted.edit();
        addition_editor.add_atom(AtomForm::from_element(Element::F));
        let (result, addition) = addition_editor.tracked_build();

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
