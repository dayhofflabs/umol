use proptest::prelude::*;
use umol_ast::ast::Transaction;

use crate::strategies::*;

proptest! {
    /// Creation ordinals are reconstructed solely from the ordered entries: uninterrupted
    /// construction, raw pushes, and `FromIterator` agree for arbitrary interleavings, with one
    /// independent namespace per entity kind.
    #[test]
    fn test_edits_creation_ordinals(
        kinds in prop::collection::vec(
            prop::sample::select(vec![
                EntityKind::Atom,
                EntityKind::Bond,
                EntityKind::DativeBond,
                EntityKind::AromaticSystem,
                EntityKind::MulticenterBond,
                EntityKind::NoncovalentBond,
                EntityKind::StereoAtom,
                EntityKind::StereoBond,
            ]),
            0..64,
        ),
    ) {
        let mut direct = Edits::new();
        let mut entries = Vec::with_capacity(kinds.len());
        let mut counts = [0usize; 8];
        for kind in kinds {
            match kind {
                EntityKind::Atom => {
                    counts[0] += 1;
                    direct.add_atom(AtomAst::default());
                    entries.push(Edit::AddAtoms { atoms: vec![AtomAst::default()] });
                }
                EntityKind::Bond => {
                    counts[1] += 1;
                    direct.add_bond(
                        AtomHandle::Id(AtomId(0)),
                        AtomHandle::Id(AtomId(1)),
                        BondAst::default(),
                    );
                    entries.push(Edit::AddBonds {
                        bonds: vec![AddBond {
                            endpoints: [
                                AtomHandle::Id(AtomId(0)),
                                AtomHandle::Id(AtomId(1)),
                            ],
                            ast: BondAst::default(),
                        }],
                    });
                }
                EntityKind::DativeBond => {
                    counts[2] += 1;
                    direct.add_dative_bond(Vec::new(), DativeBondAst::default());
                    entries.push(Edit::AddDativeBond {
                        atoms: Vec::new(),
                        ast: DativeBondAst::default(),
                    });
                }
                EntityKind::AromaticSystem => {
                    counts[3] += 1;
                    direct.add_aromatic_system(Vec::new(), AromaticSystemAst::default());
                    entries.push(Edit::AddAromaticSystem {
                        atoms: Vec::new(),
                        ast: AromaticSystemAst::default(),
                    });
                }
                EntityKind::MulticenterBond => {
                    counts[4] += 1;
                    direct.add_multicenter_bond(Vec::new(), MulticenterBondAst::default());
                    entries.push(Edit::AddMulticenterBond {
                        atoms: Vec::new(),
                        ast: MulticenterBondAst::default(),
                    });
                }
                EntityKind::NoncovalentBond => {
                    counts[5] += 1;
                    direct.add_noncovalent_bond(
                        [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        NoncovalentBondAst::default(),
                    );
                    entries.push(Edit::AddNoncovalentBond {
                        atoms: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        ast: NoncovalentBondAst::default(),
                    });
                }
                EntityKind::StereoAtom => {
                    counts[6] += 1;
                    direct.add_stereo_atom(
                        AtomHandle::Id(AtomId(0)),
                        Vec::new(),
                        StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                    );
                    entries.push(Edit::AddStereoAtom {
                        site: AtomHandle::Id(AtomId(0)),
                        ligands: Vec::new(),
                        ast: StereoAtomAst::new(
                            StereoKind::Tetrahedral,
                            StereoCoset::Lit(0),
                        ),
                    });
                }
                EntityKind::StereoBond => {
                    counts[7] += 1;
                    direct.add_stereo_bond(
                        BondHandle::Id(BondId(0)),
                        Vec::new(),
                        StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
                    );
                    entries.push(Edit::AddStereoBond {
                        site: BondHandle::Id(BondId(0)),
                        ligands: Vec::new(),
                        ast: StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
                    });
                }
            }
        }

        let mut pushed = Edits::new();
        for entry in entries.clone() {
            pushed.push(entry);
        }
        let mut collected: Edits = entries.clone().into_iter().collect();
        prop_assert_eq!(direct.as_slice(), entries.as_slice());
        prop_assert_eq!(pushed.as_slice(), entries.as_slice());
        prop_assert_eq!(collected.as_slice(), entries.as_slice());

        let expected = (
            AtomHandle::New(counts[0]),
            BondHandle::New(counts[1]),
            DativeBondHandle::New(counts[2]),
            AromaticSystemHandle::New(counts[3]),
            MulticenterBondHandle::New(counts[4]),
            NoncovalentBondHandle::New(counts[5]),
            StereoAtomHandle::New(counts[6]),
            StereoBondHandle::New(counts[7]),
        );
        let direct_next = (
            direct.add_atom(AtomAst::default()),
            direct.add_bond(
                AtomHandle::Id(AtomId(0)),
                AtomHandle::Id(AtomId(1)),
                BondAst::default(),
            ),
            direct.add_dative_bond(Vec::new(), DativeBondAst::default()),
            direct.add_aromatic_system(Vec::new(), AromaticSystemAst::default()),
            direct.add_multicenter_bond(Vec::new(), MulticenterBondAst::default()),
            direct.add_noncovalent_bond(
                [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                NoncovalentBondAst::default(),
            ),
            direct.add_stereo_atom(
                AtomHandle::Id(AtomId(0)),
                Vec::new(),
                StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            direct.add_stereo_bond(
                BondHandle::Id(BondId(0)),
                Vec::new(),
                StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
            ),
        );
        let pushed_next = (
            pushed.add_atom(AtomAst::default()),
            pushed.add_bond(
                AtomHandle::Id(AtomId(0)),
                AtomHandle::Id(AtomId(1)),
                BondAst::default(),
            ),
            pushed.add_dative_bond(Vec::new(), DativeBondAst::default()),
            pushed.add_aromatic_system(Vec::new(), AromaticSystemAst::default()),
            pushed.add_multicenter_bond(Vec::new(), MulticenterBondAst::default()),
            pushed.add_noncovalent_bond(
                [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                NoncovalentBondAst::default(),
            ),
            pushed.add_stereo_atom(
                AtomHandle::Id(AtomId(0)),
                Vec::new(),
                StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            pushed.add_stereo_bond(
                BondHandle::Id(BondId(0)),
                Vec::new(),
                StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
            ),
        );
        let collected_next = (
            collected.add_atom(AtomAst::default()),
            collected.add_bond(
                AtomHandle::Id(AtomId(0)),
                AtomHandle::Id(AtomId(1)),
                BondAst::default(),
            ),
            collected.add_dative_bond(Vec::new(), DativeBondAst::default()),
            collected.add_aromatic_system(Vec::new(), AromaticSystemAst::default()),
            collected.add_multicenter_bond(Vec::new(), MulticenterBondAst::default()),
            collected.add_noncovalent_bond(
                [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                NoncovalentBondAst::default(),
            ),
            collected.add_stereo_atom(
                AtomHandle::Id(AtomId(0)),
                Vec::new(),
                StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            collected.add_stereo_bond(
                BondHandle::Id(BondId(0)),
                Vec::new(),
                StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
            ),
        );
        prop_assert_eq!(&direct_next, &expected);
        prop_assert_eq!(&pushed_next, &expected);
        prop_assert_eq!(&collected_next, &expected);
    }

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
        let mut edits = Edits::new();
        edits.update_atom(AtomHandle::Id(AtomId(0)), &atom, &atom_update);
        edits.update_bond(
            BondHandle::Id(BondId(0)),
            &bond,
            &bond_update,
        );
        edits.update_aromatic_system(
            AromaticSystemHandle::Id(AromaticSystemId(0)),
            &aromatic,
            &aromatic_update,
        );
        edits.update_multicenter_bond(
            MulticenterBondHandle::Id(MulticenterBondId(0)),
            &multicenter,
            &multicenter_update,
        );

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
                .transact(Edits::from_iter([edit]))
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
