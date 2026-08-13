//! Context-free normalization of the molecule-scope constraint tree:
//! recurse, flatten same-combinator children, drop empty combinator children,
//! sort + dedup, and reduce trivial wrappers — a singleton `And`/`Or` is its
//! element. The raw operational domain includes singleton and duplicate-child
//! combinators.

use proptest::prelude::*;
use umol_graph_ir::ir::{Constraint, Normalize};

use crate::strategies::*;

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
}
