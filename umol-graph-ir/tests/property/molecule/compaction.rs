//! Molecule-editor compaction properties.

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
        prop_assert_eq!(plain.try_build(), editor.try_build());
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
}
