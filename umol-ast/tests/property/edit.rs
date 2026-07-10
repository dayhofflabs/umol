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
    fn test_molecule_editor_transact_rollback(
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
    fn test_molecule_editor_transact_unchecked((base, edits) in transaction_edits_strategy()) {
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
        let mut atom_keys: HashSet<(AtomId, AtomConstraintKey)> = HashSet::new();
        let mut bond_keys: HashSet<(BondId, BondConstraintKey)> = HashSet::new();
        let mut dative_keys: HashSet<(DativeBondId, DativeBondConstraintKey)> = HashSet::new();
        let mut aromatic_keys: HashSet<(AromaticSystemId, AromaticSystemConstraintKey)> =
            HashSet::new();
        let mut multicenter_keys: HashSet<(MulticenterBondId, MulticenterBondConstraintKey)> =
            HashSet::new();
        for c in ast.constraints().iter() {
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

        let mut a = ast;
        a.inline_constraints();

        for (id, key) in atom_keys {
            prop_assert!(
                a[id].constraints.contains(key),
                "atom {id:?} missing key {key:?} after inline",
            );
        }
        for (id, key) in bond_keys {
            prop_assert!(
                a[id].constraints.contains(key),
                "bond {id:?} missing key {key:?} after inline",
            );
        }
        for (id, key) in dative_keys {
            prop_assert!(
                a[id].constraints.contains(key),
                "dative bond {id:?} missing key {key:?} after inline",
            );
        }
        for (id, key) in aromatic_keys {
            prop_assert!(
                a[id].constraints.contains(key),
                "aromatic system {id:?} missing key {key:?} after inline",
            );
        }
        for (id, key) in multicenter_keys {
            prop_assert!(
                a[id].constraints.contains(key),
                "multicenter bond {id:?} missing key {key:?} after inline",
            );
        }
    }
}
