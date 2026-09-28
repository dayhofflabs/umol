//! Context-free normalization of the molecule-scope constraint tree:
//! recurse, flatten same-combinator children, drop empty combinator children,
//! sort + dedup, and reduce trivial wrappers — a singleton `And`/`Or` is its
//! element. The raw operational domain includes singleton and duplicate-child
//! combinators.
//!
//! Constraint storage preserves list order and multiplicity through restoration.
//! Compaction roundtrips use bounded entity ids; independent list edits check
//! saved positions and values. Manipulated histories require only panic freedom.

use proptest::prelude::*;
use umol_graph_core::{Compaction, EdgeId, GraphCompaction, NodeId};
use umol_graph_ir::ir::{
    AromaticSystemId, CascadedConstraints, Constraint, Constraints, DativeBondId,
    ModifiedConstraint, MoleculeCompaction, MulticenterBondId, NoncovalentBondId, Normalize,
    RemovedConstraint, StereoAtomId, StereoBondId,
};

use crate::strategies::raw_constraint_tree_strategy;

/// Whether any `And`/`Or` node with exactly one child occurs at any depth.
fn has_singleton_combinator(constraint: &Constraint) -> bool {
    match constraint {
        Constraint::And(children) | Constraint::Or(children) => {
            children.len() == 1 || children.iter().any(has_singleton_combinator)
        }
        Constraint::Not(inner) => has_singleton_combinator(inner),
        _ => false,
    }
}

proptest! {
    /// `normalize` is idempotent on the constraint tree: normalizing the
    /// normal form is a no-op.
    #[test]
    fn test_constraint_normalize_idempotent(c in raw_constraint_tree_strategy()) {
        if let Ok(normalized) = c.normalize() {
            prop_assert_eq!(normalized.clone().normalize(), Ok(normalized));
        }
    }

    /// The normal form carries no trivial wrapper: no `And`/`Or` with exactly
    /// one child remains at any depth.
    #[test]
    fn test_constraint_normalize_normal_form(c in raw_constraint_tree_strategy()) {
        if let Ok(normalized) = c.normalize() {
            prop_assert!(!has_singleton_combinator(&normalized));
        }
    }

    /// Tracked and plain compaction agree; matching restoration recovers the raw list.
    #[test]
    fn test_constraints_tracked_compact(
        entries in prop::collection::vec(raw_constraint_tree_strategy(), 0..12),
        masks in prop::array::uniform8(0u8..16),
    ) {
        let ids = |mask: u8| (0..4).filter(move |id| mask & (1 << id) != 0);
        let compaction = MoleculeCompaction::new(
            GraphCompaction::new(
                Compaction::new(4, ids(masks[0]).map(NodeId).collect()).unwrap(),
                Compaction::new(4, ids(masks[1]).map(EdgeId).collect()).unwrap(),
            ),
            Compaction::new(4, ids(masks[2]).map(DativeBondId).collect()).unwrap(),
            Compaction::new(4, ids(masks[3]).map(AromaticSystemId).collect()).unwrap(),
            Compaction::new(4, ids(masks[4]).map(MulticenterBondId).collect()).unwrap(),
            Compaction::new(4, ids(masks[5]).map(NoncovalentBondId).collect()).unwrap(),
            Compaction::new(4, ids(masks[6]).map(StereoAtomId).collect()).unwrap(),
            Compaction::new(4, ids(masks[7]).map(StereoBondId).collect()).unwrap(),
        );
        let original = Constraints::from(entries);
        let mut plain = original.clone();
        plain.compact(&compaction);
        let mut restored = original.clone();
        let changes = restored.tracked_compact(&compaction);
        prop_assert_eq!(&restored, &plain);
        restored.restore(&changes);
        prop_assert_eq!(restored, original);
    }

    /// Independently removing and rewriting list entries records sufficient undo data.
    #[test]
    fn test_constraints_restore(
        entries in prop::collection::vec((raw_constraint_tree_strategy(), 0u8..3), 0..12),
        reverse in any::<bool>(),
    ) {
        let original = Constraints::from_iter(entries.iter().map(|(entry, _)| entry.clone()));
        let mut current = Constraints::new();
        let mut changes = CascadedConstraints::default();
        for (position, (constraint, action)) in entries.into_iter().enumerate() {
            match action {
                0 => current.push(constraint),
                1 => changes.removed.push(RemovedConstraint { position, constraint }),
                _ => {
                    let new = Constraint::Not(Box::new(constraint.clone()));
                    current.push(new.clone());
                    changes.modified.push(ModifiedConstraint { position, old: constraint, new });
                }
            }
        }
        if reverse {
            changes.removed.reverse();
            changes.modified.reverse();
        }
        current.restore(&changes);
        prop_assert_eq!(current, original);
    }

    /// Independent entries and saved positions impose no result requirement on restoration.
    #[test]
    fn test_constraints_restore_history(
        entries in prop::collection::vec(raw_constraint_tree_strategy(), 0..8),
        removed in prop::collection::vec((0usize..24, raw_constraint_tree_strategy()), 0..8),
        modified in prop::collection::vec(
            (0usize..24, raw_constraint_tree_strategy(), raw_constraint_tree_strategy()), 0..8,
        ),
    ) {
        let mut current = Constraints::from(entries);
        current.restore(&CascadedConstraints {
            removed: removed.into_iter().map(|(position, constraint)| {
                RemovedConstraint { position, constraint }
            }).collect(),
            modified: modified.into_iter().map(|(position, old, new)| {
                ModifiedConstraint { position, old, new }
            }).collect(),
        });
    }
}
