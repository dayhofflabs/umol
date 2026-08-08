//! Property tests for reaction application, including host-relative updates and stereo-frame transport.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_core::{RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm};
use umol_graph_ir::ir::{SubstructureMatchAlgorithm, SubstructureMatchConfig};

use crate::strategies::*;

const MATCH_ALGORITHM: SubstructureMatchAlgorithm = SubstructureMatchAlgorithm::GraphAndOverlays;
const SUBISO_ALGORITHM: SubgraphIsomorphismAlgorithm = SubgraphIsomorphismAlgorithm::Vf2;
const RELEVANT_CYCLE_ALGORITHM: RelevantCycleEnumerationAlgorithm =
    RelevantCycleEnumerationAlgorithm::Vismara;
const MATCH_CONFIG: SubstructureMatchConfig = SubstructureMatchConfig {
    match_algorithm: MATCH_ALGORITHM,
    subgraph_isomorphism_algorithm: SUBISO_ALGORITHM,
    relevant_cycle_algorithm: RELEVANT_CYCLE_ALGORITHM,
};

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    // These eight host-refinement properties deliberately have parallel shapes.
    // Each exercises the lowering path for a distinct entity delta; keeping them
    // separate makes a violation of host-relative old-value semantics local to
    // the affected entity family.

    /// A pattern-relative atom update lowers against the matched host atom, including independent
    /// unpaired-electron components and keyed constraint set / replace / remove operations.
    #[test]
    fn test_reaction_ast_apply_atom_update(
        host_atom in atom_ast_strategy(),
        update in atom_update_strategy(),
    ) {
        let pattern_atom = AtomAst::default();
        let effective_update = pattern_atom.difference_to(&pattern_atom.update(&update));
        let expected_atom = host_atom.update(&effective_update).canonicalize().unwrap();
        let atom_deltas = AtomDelta::for_update(AtomId(0), &pattern_atom, &effective_update);
        let reaction = ReactionAst::new(
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![AtomAst::default()],
                ..Default::default()
            }),
            Deltas::from_iter(atom_deltas.into_iter().map(Delta::Atom)),
        );
        let host = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![host_atom],
            ..Default::default()
        });
        let expected = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![expected_atom],
            ..Default::default()
        });
        let products: Vec<MoleculeAst> = reaction
            .apply(
                &host,
                MATCH_CONFIG,
            )
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();

        prop_assert_eq!(products.len(), 1);
        prop_assert!(products[0].equiv(&expected));
    }

    /// A pattern-relative localized-bond update lowers against the matched host bond.
    #[test]
    fn test_reaction_ast_apply_bond_update(
        host_bond in bond_ast_strategy(),
        update in bond_update_strategy(),
    ) {
        let pattern_bond = BondAst::default();
        let effective_update = pattern_bond.difference_to(&pattern_bond.update(&update));
        let expected_bond = host_bond.update(&effective_update).canonicalize().unwrap();
        let bond_deltas = BondDelta::for_update(BondId(0), &pattern_bond, &effective_update);
        let reaction = ReactionAst::new(
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
                bonds: vec![(AtomId(0), AtomId(1), BondAst::default())],
                ..Default::default()
            }),
            Deltas::from_iter(bond_deltas.into_iter().map(Delta::Bond)),
        );
        let host = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            bonds: vec![(AtomId(0), AtomId(1), host_bond)],
            ..Default::default()
        });
        let expected = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            bonds: vec![(AtomId(0), AtomId(1), expected_bond)],
            ..Default::default()
        });
        let products: Vec<MoleculeAst> = reaction
            .apply(
                &host,
                MATCH_CONFIG,
            )
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();

        prop_assert_eq!(products.len(), 1);
        prop_assert!(products[0].equiv(&expected));
    }

    /// A pattern-relative dative-bond update lowers against the matched host relation.
    #[test]
    fn test_reaction_ast_apply_dative_bond_update(
        host_bond in dative_bond_strategy(),
        update in dative_bond_update_strategy(),
    ) {
        let pattern_bond = DativeBondAst::default();
        let effective_update = pattern_bond.difference_to(&pattern_bond.update(&update));
        let expected_bond = host_bond.update(&effective_update).canonicalize().unwrap();
        let dative_deltas = DativeBondDelta::for_update(
            DativeBondId(0),
            &pattern_bond,
            &effective_update,
        );
        let reaction = ReactionAst::new(
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
                dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondAst::default())],
                ..Default::default()
            }),
            Deltas::from_iter(dative_deltas.into_iter().map(Delta::DativeBond)),
        );
        let host = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            dative: vec![(vec![AtomId(0)], AtomId(1), host_bond)],
            ..Default::default()
        });
        let expected = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            dative: vec![(vec![AtomId(0)], AtomId(1), expected_bond)],
            ..Default::default()
        });
        let products: Vec<MoleculeAst> = reaction
            .apply(
                &host,
                MATCH_CONFIG,
            )
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();

        prop_assert_eq!(products.len(), 1);
        prop_assert!(products[0].equiv(&expected));
    }

    /// A pattern-relative aromatic-system update lowers against the matched host relation.
    #[test]
    fn test_reaction_ast_apply_aromatic_system_update(
        mut host_system in aromatic_system_ast_for(3),
        update in aromatic_system_update_for(3),
    ) {
        host_system.unpaired_electrons = UnpairedElectronsAst::from((2_u8, 3_u8));
        let pattern_system = AromaticSystemAst::default();
        let effective_update = pattern_system.difference_to(&pattern_system.update(&update));
        let expected_system = host_system.update(&effective_update).canonicalize().unwrap();
        let aromatic_deltas = AromaticSystemDelta::for_update(
            AromaticSystemId(0),
            &pattern_system,
            &effective_update,
        );
        let reaction = ReactionAst::new(
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::N),
                    AtomAst::from_element(Element::O),
                ],
                aromatic: vec![(vec![AtomId(0), AtomId(1), AtomId(2)], AromaticSystemAst::default())],
                ..Default::default()
            }),
            Deltas::from_iter(aromatic_deltas.into_iter().map(Delta::AromaticSystem)),
        );
        let host = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::O),
            ],
            aromatic: vec![(vec![AtomId(0), AtomId(1), AtomId(2)], host_system)],
            ..Default::default()
        });
        let expected = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::O),
            ],
            aromatic: vec![(vec![AtomId(0), AtomId(1), AtomId(2)], expected_system)],
            ..Default::default()
        });
        let products: Vec<MoleculeAst> = reaction
            .apply(
                &host,
                MATCH_CONFIG,
            )
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();

        prop_assert_eq!(products.len(), 1);
        prop_assert!(products[0].equiv(&expected));
    }

    /// A pattern-relative multicenter-bond update lowers against the matched host relation.
    #[test]
    fn test_reaction_ast_apply_multicenter_bond_update(
        mut host_bond in multicenter_bond_ast_for(3),
        update in multicenter_bond_update_for(3),
    ) {
        host_bond.unpaired_electrons = UnpairedElectronsAst::from((2_u8, 3_u8));
        let pattern_bond = MulticenterBondAst::default();
        let effective_update = pattern_bond.difference_to(&pattern_bond.update(&update));
        let expected_bond = host_bond.update(&effective_update).canonicalize().unwrap();
        let multicenter_deltas = MulticenterBondDelta::for_update(
            MulticenterBondId(0),
            &pattern_bond,
            &effective_update,
        );
        let reaction = ReactionAst::new(
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::N),
                    AtomAst::from_element(Element::O),
                ],
                multicenter: vec![(vec![AtomId(0), AtomId(1), AtomId(2)], MulticenterBondAst::default())],
                ..Default::default()
            }),
            Deltas::from_iter(multicenter_deltas.into_iter().map(Delta::MulticenterBond)),
        );
        let host = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::O),
            ],
            multicenter: vec![(vec![AtomId(0), AtomId(1), AtomId(2)], host_bond)],
            ..Default::default()
        });
        let expected = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::O),
            ],
            multicenter: vec![(vec![AtomId(0), AtomId(1), AtomId(2)], expected_bond)],
            ..Default::default()
        });
        let products: Vec<MoleculeAst> = reaction
            .apply(
                &host,
                MATCH_CONFIG,
            )
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();

        prop_assert_eq!(products.len(), 1);
        prop_assert!(products[0].equiv(&expected));
    }

    /// A pattern-relative noncovalent-bond update lowers against the matched host relation.
    #[test]
    fn test_reaction_ast_apply_noncovalent_bond_update(
        host_bond in noncovalent_bond_ast_strategy(),
        update in noncovalent_bond_update_strategy(),
    ) {
        let pattern_bond = NoncovalentBondAst::default();
        let effective_update = pattern_bond.difference_to(&pattern_bond.update(&update));
        let expected_bond = host_bond.update(&effective_update).canonicalize().unwrap();
        let noncovalent_deltas = NoncovalentBondDelta::for_update(
            NoncovalentBondId(0),
            &pattern_bond,
            &effective_update,
        );
        let reaction = ReactionAst::new(
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
                noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondAst::default())],
                ..Default::default()
            }),
            Deltas::from_iter(noncovalent_deltas.into_iter().map(Delta::NoncovalentBond)),
        );
        let host = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            noncovalent: vec![(AtomId(0), AtomId(1), host_bond)],
            ..Default::default()
        });
        let expected = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            noncovalent: vec![(AtomId(0), AtomId(1), expected_bond)],
            ..Default::default()
        });
        let products: Vec<MoleculeAst> = reaction
            .apply(
                &host,
                MATCH_CONFIG,
            )
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();

        prop_assert_eq!(products.len(), 1);
        prop_assert!(products[0].equiv(&expected));
    }

    /// A pattern-relative stereo-atom update lowers against the matched host configuration and
    /// keyed constraints.
    #[test]
    fn test_reaction_ast_apply_stereo_atom_update(
        host_coset in stereo_coset_for_kind(StereoKind::Tetrahedral),
        host_constraints in stereo_atom_constraints_strategy(StereoKind::Tetrahedral),
        update in stereo_atom_application_update_strategy(),
    ) {
        let ligands = vec![
            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        ];
        let pattern_atom = StereoAtomAst::new(
            StereoKind::Tetrahedral,
            StereoCoset::Undetermined,
        );
        let host_atom = StereoAtomAst::new(StereoKind::Tetrahedral, host_coset)
            .with_constraints(host_constraints);
        let effective_update = pattern_atom.difference_to(&pattern_atom.update(&update));
        let expected_atom = host_atom.update(&effective_update).canonicalize().unwrap();
        let stereo_atom_deltas =
            StereoAtomDelta::for_update(StereoAtomId(0), &pattern_atom, &effective_update);
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::N),
            AtomAst::from_element(Element::O),
            AtomAst::from_element(Element::F),
        ];
        let reaction = ReactionAst::new(
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: atoms.clone(),
                stereo_atoms: vec![(AtomId(0), ligands.clone(), pattern_atom.clone())],
                ..Default::default()
            }),
            Deltas::from_iter(stereo_atom_deltas.into_iter().map(Delta::StereoAtom)),
        );
        let host = MoleculeAst::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
            stereo_atoms: vec![(AtomId(0), ligands.clone(), host_atom)],
            ..Default::default()
        });
        let expected = MoleculeAst::from_entries(MoleculeEntries {
            atoms,
            stereo_atoms: vec![(AtomId(0), ligands, expected_atom)],
            ..Default::default()
        });
        let products: Vec<MoleculeAst> = reaction
            .apply(
                &host,
                MATCH_CONFIG,
            )
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();

        prop_assert_eq!(products.len(), 1);
        prop_assert!(products[0].equiv(&expected));
    }

    /// A pattern-relative stereo-bond update lowers against the matched host configuration and
    /// keyed constraints.
    #[test]
    fn test_reaction_ast_apply_stereo_bond_update(
        host_coset in stereo_coset_for_kind(StereoKind::CisTrans),
        host_constraints in stereo_bond_constraints_strategy(StereoKind::CisTrans),
        update in stereo_bond_application_update_strategy(),
    ) {
        let ligands = vec![
            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        ];
        let pattern_bond = StereoBondAst::new(
            StereoKind::CisTrans,
            StereoCoset::Undetermined,
        );
        let host_bond = StereoBondAst::new(StereoKind::CisTrans, host_coset)
            .with_constraints(host_constraints);
        let effective_update = pattern_bond.difference_to(&pattern_bond.update(&update));
        let expected_bond = host_bond.update(&effective_update).canonicalize().unwrap();
        let stereo_bond_deltas =
            StereoBondDelta::for_update(StereoBondId(0), &pattern_bond, &effective_update);
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::N),
            AtomAst::from_element(Element::O),
            AtomAst::from_element(Element::F),
        ];
        let bonds = vec![(AtomId(0), AtomId(1), BondAst::from_order(2))];
        let reaction = ReactionAst::new(
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: atoms.clone(),
                bonds: bonds.clone(),
                stereo_bonds: vec![(BondId(0), ligands.clone(), pattern_bond.clone())],
                ..Default::default()
            }),
            Deltas::from_iter(stereo_bond_deltas.into_iter().map(Delta::StereoBond)),
        );
        let host = MoleculeAst::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_bonds: vec![(BondId(0), ligands.clone(), host_bond)],
            ..Default::default()
        });
        let expected = MoleculeAst::from_entries(MoleculeEntries {
            atoms,
            bonds,
            stereo_bonds: vec![(BondId(0), ligands, expected_bond)],
            ..Default::default()
        });
        let products: Vec<MoleculeAst> = reaction
            .apply(
                &host,
                MATCH_CONFIG,
            )
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();

        prop_assert_eq!(products.len(), 1);
        prop_assert!(products[0].equiv(&expected));
    }

    /// Applying a reaction at the identity occurrence of its own `lhs` reproduces the span's
    /// `right()` — the `transact`-apply path agrees with the span projection.
    #[test]
    fn test_reaction_ast_apply_reproduces_right(reaction in reaction_strategy()) {
        if let Ok(span) = reaction.to_reaction_span() {
            let right = span.rhs();
            prop_assert!(reaction
                .apply(
                    &reaction.lhs,
                    MATCH_CONFIG,
                )
                .unwrap()
                .any(|derivation| derivation.unwrap().rhs() == &right));
        }
    }

    /// Isolation probe: a plain overlay reaction's `apply` at its own `lhs` reproduces its
    /// `right()`. If this fails, the discrepancy is in apply-vs-span for overlays, not compose.
    #[test]
    fn test_reaction_ast_apply_reproduces_right_overlay(reaction in overlay_reaction_strategy()) {
        if let Ok(span) = reaction.to_reaction_span() {
            let right = span.rhs();
            prop_assert!(reaction
                .apply(
                    &reaction.lhs,
                    MATCH_CONFIG,
                )
                .unwrap()
                .any(|derivation| derivation.unwrap().rhs() == &right));
        }
    }

    #[test]
    fn test_reaction_ast_apply_reframes_stereo_atom_modification(
        old in 0..StereoKind::Tetrahedral.count() as u32,
        permutation in stereo_frame_permutation_strategy(StereoKind::Tetrahedral),
    ) {
        let new = 1 - old;
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::F),
            AtomAst::from_element(Element::Cl),
            AtomAst::from_element(Element::Br),
            AtomAst::from_element(Element::I),
        ];
        let bonds: Vec<(AtomId, AtomId, BondAst)> = (1..=4)
            .map(|ligand| (AtomId(0), AtomId(ligand), BondAst::from_order(1)))
            .collect();
        let rule_frame: Vec<StereoLigand> = (1..=4)
            .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
            .collect();
        let old_ast = StereoAtomAst::new(StereoKind::Tetrahedral, old);
        let new_ast = StereoAtomAst::new(StereoKind::Tetrahedral, new);
        let lhs = MoleculeAst::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_atoms: vec![(AtomId(0), rule_frame.clone(), old_ast.clone())],
            ..Default::default()
        });
        let reaction = ReactionAst::new(
            lhs,
            Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::ModifyField {
                id: StereoAtomId(0),
                change: StereoAtomFieldChange::Configuration {
                    old: old_ast.configuration,
                    new: new_ast.configuration,
                },
            })]),
        );
        let host_frame = permutation.act(&rule_frame);
        let host = MoleculeAst::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_atoms: vec![(
                AtomId(0),
                host_frame.clone(),
                StereoAtomAst::new(StereoKind::Tetrahedral, old).apply(permutation),
            )],
            ..Default::default()
        });
        let expected = MoleculeAst::from_entries(MoleculeEntries {
            atoms,
            bonds,
            stereo_atoms: vec![(
                AtomId(0),
                host_frame,
                StereoAtomAst::new(StereoKind::Tetrahedral, new).apply(permutation),
            )],
            ..Default::default()
        });
        let products: Vec<MoleculeAst> = reaction
            .apply(
                &host,
                MATCH_CONFIG,
            )
            .map_err(|error| TestCaseError::fail(format!("application precondition: {error:?}")))?
            .map(|result| result.map(|derivation| derivation.rhs().clone()))
            .collect::<Result<_, _>>()
            .map_err(|error| TestCaseError::fail(format!("application failed: {error:?}")))?;

        prop_assert_eq!(products, vec![expected]);
    }

    #[test]
    fn test_reaction_ast_apply_reframes_stereo_atom_removal(
        coset in 0..StereoKind::Tetrahedral.count() as u32,
        permutation in stereo_frame_permutation_strategy(StereoKind::Tetrahedral),
    ) {
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::F),
            AtomAst::from_element(Element::Cl),
            AtomAst::from_element(Element::Br),
            AtomAst::from_element(Element::I),
        ];
        let bonds: Vec<(AtomId, AtomId, BondAst)> = (1..=4)
            .map(|ligand| (AtomId(0), AtomId(ligand), BondAst::from_order(1)))
            .collect();
        let rule_frame: Vec<StereoLigand> = (1..=4)
            .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
            .collect();
        let rule_ast = StereoAtomAst::new(StereoKind::Tetrahedral, coset);
        let lhs = MoleculeAst::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_atoms: vec![(AtomId(0), rule_frame.clone(), rule_ast.clone())],
            ..Default::default()
        });
        let reaction = ReactionAst::new(
            lhs,
            Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::Remove {
                id: StereoAtomId(0),
                site: AtomId(0),
                ligands: rule_frame.clone(),
                ast: rule_ast,
            })]),
        );
        let host = MoleculeAst::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_atoms: vec![(
                AtomId(0),
                permutation.act(&rule_frame),
                StereoAtomAst::new(StereoKind::Tetrahedral, coset).apply(permutation),
            )],
            ..Default::default()
        });
        let expected = MoleculeAst::from_entries(MoleculeEntries {
            atoms,
            bonds,
            ..Default::default()
        });
        let products: Vec<MoleculeAst> = reaction
            .apply(
                &host,
                MATCH_CONFIG,
            )
            .map_err(|error| TestCaseError::fail(format!("application precondition: {error:?}")))?
            .map(|result| result.map(|derivation| derivation.rhs().clone()))
            .collect::<Result<_, _>>()
            .map_err(|error| TestCaseError::fail(format!("application failed: {error:?}")))?;

        prop_assert_eq!(products, vec![expected]);
    }

    #[test]
    fn test_reaction_ast_apply_reframes_stereo_bond_modification(
        old in 0..StereoKind::CisTrans.count() as u32,
        permutation in stereo_frame_permutation_strategy(StereoKind::CisTrans),
    ) {
        let new = 1 - old;
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::F),
            AtomAst::from_element(Element::Cl),
            AtomAst::from_element(Element::Br),
            AtomAst::from_element(Element::I),
        ];
        let bonds = vec![
            (AtomId(0), AtomId(1), BondAst::from_order(2)),
            (AtomId(0), AtomId(2), BondAst::from_order(1)),
            (AtomId(0), AtomId(3), BondAst::from_order(1)),
            (AtomId(1), AtomId(4), BondAst::from_order(1)),
            (AtomId(1), AtomId(5), BondAst::from_order(1)),
        ];
        let rule_frame: Vec<StereoLigand> = (2..=5)
            .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
            .collect();
        let old_ast = StereoBondAst::new(StereoKind::CisTrans, old);
        let new_ast = StereoBondAst::new(StereoKind::CisTrans, new);
        let lhs = MoleculeAst::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_bonds: vec![(BondId(0), rule_frame.clone(), old_ast.clone())],
            ..Default::default()
        });
        let reaction = ReactionAst::new(
            lhs,
            Deltas::from_iter([Delta::StereoBond(StereoBondDelta::ModifyField {
                id: StereoBondId(0),
                change: StereoBondFieldChange::Configuration {
                    old: old_ast.configuration,
                    new: new_ast.configuration,
                },
            })]),
        );
        let host_frame = permutation.act(&rule_frame);
        let host = MoleculeAst::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_bonds: vec![(
                BondId(0),
                host_frame.clone(),
                StereoBondAst::new(StereoKind::CisTrans, old).apply(permutation),
            )],
            ..Default::default()
        });
        let expected = MoleculeAst::from_entries(MoleculeEntries {
            atoms,
            bonds,
            stereo_bonds: vec![(
                BondId(0),
                host_frame,
                StereoBondAst::new(StereoKind::CisTrans, new).apply(permutation),
            )],
            ..Default::default()
        });
        let products: Vec<MoleculeAst> = reaction
            .apply(
                &host,
                MATCH_CONFIG,
            )
            .map_err(|error| TestCaseError::fail(format!("application precondition: {error:?}")))?
            .map(|result| result.map(|derivation| derivation.rhs().clone()))
            .collect::<Result<_, _>>()
            .map_err(|error| TestCaseError::fail(format!("application failed: {error:?}")))?;

        prop_assert_eq!(products, vec![expected]);
    }

    #[test]
    fn test_reaction_ast_apply_reframes_stereo_bond_removal(
        coset in 0..StereoKind::CisTrans.count() as u32,
        permutation in stereo_frame_permutation_strategy(StereoKind::CisTrans),
    ) {
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::F),
            AtomAst::from_element(Element::Cl),
            AtomAst::from_element(Element::Br),
            AtomAst::from_element(Element::I),
        ];
        let bonds = vec![
            (AtomId(0), AtomId(1), BondAst::from_order(2)),
            (AtomId(0), AtomId(2), BondAst::from_order(1)),
            (AtomId(0), AtomId(3), BondAst::from_order(1)),
            (AtomId(1), AtomId(4), BondAst::from_order(1)),
            (AtomId(1), AtomId(5), BondAst::from_order(1)),
        ];
        let rule_frame: Vec<StereoLigand> = (2..=5)
            .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
            .collect();
        let rule_ast = StereoBondAst::new(StereoKind::CisTrans, coset);
        let lhs = MoleculeAst::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_bonds: vec![(BondId(0), rule_frame.clone(), rule_ast.clone())],
            ..Default::default()
        });
        let reaction = ReactionAst::new(
            lhs,
            Deltas::from_iter([Delta::StereoBond(StereoBondDelta::Remove {
                id: StereoBondId(0),
                site: BondId(0),
                ligands: rule_frame.clone(),
                ast: rule_ast,
            })]),
        );
        let host = MoleculeAst::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
            bonds: bonds.clone(),
            stereo_bonds: vec![(
                BondId(0),
                permutation.act(&rule_frame),
                StereoBondAst::new(StereoKind::CisTrans, coset).apply(permutation),
            )],
            ..Default::default()
        });
        let expected = MoleculeAst::from_entries(MoleculeEntries {
            atoms,
            bonds,
            ..Default::default()
        });
        let products: Vec<MoleculeAst> = reaction
            .apply(
                &host,
                MATCH_CONFIG,
            )
            .map_err(|error| TestCaseError::fail(format!("application precondition: {error:?}")))?
            .map(|result| result.map(|derivation| derivation.rhs().clone()))
            .collect::<Result<_, _>>()
            .map_err(|error| TestCaseError::fail(format!("application failed: {error:?}")))?;

        prop_assert_eq!(products, vec![expected]);
    }
}
