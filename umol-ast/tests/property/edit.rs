use proptest::prelude::*;

use crate::strategies::*;

proptest! {
    /// `lift_constraints` followed by `inline_constraints` is idempotent:
    /// running the pair twice yields the same `MoleculeAst` as running it
    /// once. This holds even if the original AST has duplicate (entity, kind)
    /// entries across the inline + molecule scopes — the first pass collapses
    /// them via the entity store's last-wins policy and the second pass is
    /// a fixpoint.
    #[test]
    fn test_lift_inline_idempotent(ast in molecule_ast_with_constraints_strategy()) {
        let mut once = ast.clone();
        once.lift_constraints();
        once.inline_constraints();

        let mut twice = once.clone();
        twice.lift_constraints();
        twice.inline_constraints();

        prop_assert_eq!(once, twice);
    }

    /// `lift_constraints` drains every entity's inline `constraints` store.
    #[test]
    fn test_lift_drains_entity_stores(ast in molecule_ast_with_constraints_strategy()) {
        let mut a = ast;
        a.lift_constraints();
        for view in a.atoms().iter() {
            prop_assert!(view.ast.constraints.is_empty());
        }
        for view in a.bonds().iter() {
            prop_assert!(view.ast.constraints.is_empty());
        }
        for view in a.dative_bonds().iter() {
            prop_assert!(view.ast.constraints.is_empty());
        }
        for view in a.aromatic_systems().iter() {
            prop_assert!(view.ast.constraints.is_empty());
        }
        for view in a.multicenter_bonds().iter() {
            prop_assert!(view.ast.constraints.is_empty());
        }
        for view in a.noncovalent_bonds().iter() {
            prop_assert!(view.ast.constraints.is_empty());
        }
    }

    #[test]
    fn test_molecule_builder_transact_rollback(
        (base, edits) in transaction_edits_strategy(),
    ) {
        let mut builder = base.edit();
        let before = builder.clone().build();
        let tx = builder
            .transact(edits)
            .map_err(|e| TestCaseError::fail(format!("transact failed: {e}")))?;

        tx.rollback(&mut builder)
            .map_err(|e| TestCaseError::fail(format!("rollback failed: {e}")))?;

        prop_assert_eq!(builder.build(), before);
    }

    #[test]
    fn test_molecule_builder_transact_unchecked((base, edits) in transaction_edits_strategy()) {
        let mut checked = base.edit();
        checked
            .transact(edits.clone())
            .map_err(|e| TestCaseError::fail(format!("checked transact failed: {e}")))?;

        let mut unchecked = base.edit();
        unchecked.transact_unchecked(edits);

        prop_assert_eq!(unchecked.build(), checked.build());
    }

    /// `inline_constraints` removes every TOP-LEVEL inline-capable narrow
    /// leaf from the molecule list. Combinator-nested entries, relational
    /// leaves, molecule-scope leaves are preserved.
    #[test]
    fn test_inline_removes_top_level_leaves(ast in molecule_ast_with_constraints_strategy()) {
        let mut a = ast;
        a.inline_constraints();
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
    /// `kind()`. Last-wins per kind: if the same `(id, kind)` appeared
    /// multiple times, or if the entity already had an inline same-kind
    /// entry, the kind is still present after the call.
    #[test]
    fn test_inline_deposits_leaves_into_entities(
        ast in molecule_ast_with_constraints_strategy(),
    ) {
        let mut atom_kinds: HashSet<(AtomId, AtomConstraintKind)> = HashSet::new();
        let mut bond_kinds: HashSet<(BondId, BondConstraintKind)> = HashSet::new();
        let mut dative_kinds: HashSet<(DativeBondId, DativeBondConstraintKind)> = HashSet::new();
        let mut aromatic_kinds: HashSet<(AromaticSystemId, AromaticSystemConstraintKind)> =
            HashSet::new();
        let mut multicenter_kinds: HashSet<(MulticenterBondId, MulticenterBondConstraintKind)> =
            HashSet::new();
        for c in ast.constraints().iter() {
            match c {
                Constraint::Atom(id, inner) => {
                    atom_kinds.insert((*id, inner.kind()));
                }
                Constraint::Bond(id, inner) => {
                    bond_kinds.insert((*id, inner.kind()));
                }
                Constraint::DativeBond(id, inner) => {
                    dative_kinds.insert((*id, inner.kind()));
                }
                Constraint::AromaticSystem(id, inner) => {
                    aromatic_kinds.insert((*id, inner.kind()));
                }
                Constraint::MulticenterBond(id, inner) => {
                    multicenter_kinds.insert((*id, inner.kind()));
                }
                _ => {}
            }
        }

        let mut a = ast;
        a.inline_constraints();

        for (id, kind) in atom_kinds {
            prop_assert!(
                a[id].constraints.iter().any(|c| c.kind() == kind),
                "atom {id:?} missing kind {kind:?} after inline",
            );
        }
        for (id, kind) in bond_kinds {
            prop_assert!(
                a[id].constraints.iter().any(|c| c.kind() == kind),
                "bond {id:?} missing kind {kind:?} after inline",
            );
        }
        for (id, kind) in dative_kinds {
            prop_assert!(
                a[id].constraints.iter().any(|c| c.kind() == kind),
                "dative bond {id:?} missing kind {kind:?} after inline",
            );
        }
        for (id, kind) in aromatic_kinds {
            prop_assert!(
                a[id].constraints.iter().any(|c| c.kind() == kind),
                "aromatic system {id:?} missing kind {kind:?} after inline",
            );
        }
        for (id, kind) in multicenter_kinds {
            prop_assert!(
                a[id].constraints.iter().any(|c| c.kind() == kind),
                "multicenter bond {id:?} missing kind {kind:?} after inline",
            );
        }
    }
}
