use proptest::prelude::*;
use umol_ast::ast::Transaction;

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
    fn test_molecule_editor_transact_rollback_unpaired_electrons(
        atom_components in partial_unpaired_electrons_update_strategy(),
        bond_components in partial_unpaired_electrons_update_strategy(),
        aromatic_components in partial_unpaired_electrons_update_strategy(),
        multicenter_components in partial_unpaired_electrons_update_strategy(),
    ) {
        let atom = AtomAst::from_element(Element::C).with_unpaired_electrons((2_u8, 3_u8));
        let bond = BondAst::from_order(1).with_unpaired_electrons((2_u8, 3_u8));
        let aromatic = AromaticSystemAst::from_electrons(vec![1, 1, 1])
            .with_unpaired_electrons((2_u8, 3_u8));
        let multicenter = MulticenterBondAst::from_electrons(vec![1, 1, 1])
            .with_unpaired_electrons((2_u8, 3_u8));
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
        let base = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![atom.clone(), AtomAst::from_element(Element::N), AtomAst::from_element(Element::O)],
            bonds: vec![(AtomId(0), AtomId(1), bond.clone())],
            aromatic: vec![(vec![AtomId(0), AtomId(1), AtomId(2)], aromatic.clone())],
            multicenter: vec![(vec![AtomId(0), AtomId(1), AtomId(2)], multicenter.clone())],
            ..Default::default()
        });
        let expected = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![atom.update(&atom_update), AtomAst::from_element(Element::N), AtomAst::from_element(Element::O)],
            bonds: vec![(AtomId(0), AtomId(1), bond.update(&bond_update))],
            aromatic: vec![(vec![AtomId(0), AtomId(1), AtomId(2)], aromatic.update(&aromatic_update))],
            multicenter: vec![(vec![AtomId(0), AtomId(1), AtomId(2)], multicenter.update(&multicenter_update))],
            ..Default::default()
        });
        let mut edits = Edit::for_atom_update(AtomHandle::Id(AtomId(0)), &atom, &atom_update);
        edits.extend(Edit::for_bond_update(
            BondHandle::Id(BondId(0)),
            &bond,
            &bond_update,
        ));
        edits.extend(Edit::for_aromatic_system_update(
            AromaticSystemHandle::Id(AromaticSystemId(0)),
            &aromatic,
            &aromatic_update,
        ));
        edits.extend(Edit::for_multicenter_bond_update(
            MulticenterBondHandle::Id(MulticenterBondId(0)),
            &multicenter,
            &multicenter_update,
        ));

        let mut editor = base.clone().edit();
        let transaction = editor
            .transact(edits)
            .map_err(|error| TestCaseError::fail(format!("transact failed: {error}")))?;
        prop_assert_eq!(editor.clone().build(), expected);

        transaction
            .rollback(&mut editor)
            .map_err(|error| TestCaseError::fail(format!("rollback failed: {error}")))?;
        prop_assert_eq!(editor.build(), base);
    }

    #[test]
    fn test_transaction_append_materialization(
        (base, edits) in overlay_transaction_strategy(),
    ) {
        prop_assume!(!edits.is_empty());

        let mut single = base.edit();
        single
            .transact(edits.clone())
            .map_err(|e| TestCaseError::fail(format!("single transact failed: {e}")))?;
        let expected = single.build();

        let mut staged = base.edit();
        let mut combined = Transaction::default();
        for edit in edits {
            let transaction = staged
                .transact(vec![edit])
                .map_err(|e| TestCaseError::fail(format!("staged transact failed: {e}")))?;
            combined.append(transaction);
            let state = staged.build();
            staged = state.edit();
        }
        let materialized = staged.build();
        prop_assert_eq!(&materialized, &expected);

        let mut staged = materialized.edit();
        combined
            .rollback(&mut staged)
            .map_err(|e| TestCaseError::fail(format!("combined rollback failed: {e}")))?;
        prop_assert_eq!(staged.build(), base);
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
                a.atom(id).ast.constraints.contains(key),
                "atom {id:?} missing key {key:?} after inline",
            );
        }
        for (id, key) in bond_keys {
            prop_assert!(
                a.bond(id).ast.constraints.contains(key),
                "bond {id:?} missing key {key:?} after inline",
            );
        }
        for (id, key) in dative_keys {
            prop_assert!(
                a.dative_bond(id).ast.constraints.contains(key),
                "dative bond {id:?} missing key {key:?} after inline",
            );
        }
        for (id, key) in aromatic_keys {
            prop_assert!(
                a.aromatic_system(id).ast.constraints.contains(key),
                "aromatic system {id:?} missing key {key:?} after inline",
            );
        }
        for (id, key) in multicenter_keys {
            prop_assert!(
                a.multicenter_bond(id).ast.constraints.contains(key),
                "multicenter bond {id:?} missing key {key:?} after inline",
            );
        }
    }
}
