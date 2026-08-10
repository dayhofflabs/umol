//! Property tests for malformed reaction application inputs.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_core::{RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm};
use umol_graph_ir::ir::{
    ApplyError, ApplyPreconditionError, AromaticSystemDelta, AtomDelta, BondDelta, ConstraintDelta,
    DativeBondDelta, Entity, MulticenterBondDelta, NoncovalentBondDelta, StereoAtomDelta,
    StereoBondDelta, SubstructureMatchAlgorithm, SubstructureMatchConfig, TransactionError,
};

use crate::strategies::*;

fn unavailable_entity_strategy() -> impl Strategy<Value = (Reaction, ApplyPreconditionError)> {
    prop_oneof![
        (0u32..64).prop_map(|id| {
            (
                Reaction::new(
                    Molecule::default(),
                    Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                        id: AtomId(id),
                        attributes: AtomForm::default(),
                    })]),
                ),
                ApplyPreconditionError::InvalidReactionReference {
                    entity: Entity::Atom(AtomId(id)),
                },
            )
        }),
        (0u32..64).prop_map(|id| {
            (
                Reaction::new(
                    Molecule::default(),
                    Deltas::from_iter([Delta::Bond(BondDelta::Remove {
                        id: BondId(id),
                        atoms: [AtomId(0), AtomId(1)],
                        attributes: BondForm::default(),
                    })]),
                ),
                ApplyPreconditionError::InvalidReactionReference {
                    entity: Entity::Bond(BondId(id)),
                },
            )
        }),
        (0u32..64).prop_map(|id| {
            (
                Reaction::new(
                    Molecule::default(),
                    Deltas::from_iter([Delta::DativeBond(DativeBondDelta::Remove {
                        id: DativeBondId(id),
                        donors: vec![AtomId(0)],
                        acceptor: AtomId(1),
                        attributes: DativeBondForm::default(),
                    })]),
                ),
                ApplyPreconditionError::InvalidReactionReference {
                    entity: Entity::DativeBond(DativeBondId(id)),
                },
            )
        }),
        (0u32..64).prop_map(|id| {
            (
                Reaction::new(
                    Molecule::default(),
                    Deltas::from_iter([Delta::AromaticSystem(AromaticSystemDelta::Remove {
                        id: AromaticSystemId(id),
                        atoms: vec![AtomId(0), AtomId(1)],
                        attributes: AromaticSystemForm::default(),
                    })]),
                ),
                ApplyPreconditionError::InvalidReactionReference {
                    entity: Entity::AromaticSystem(AromaticSystemId(id)),
                },
            )
        }),
        (0u32..64).prop_map(|id| {
            (
                Reaction::new(
                    Molecule::default(),
                    Deltas::from_iter([Delta::MulticenterBond(MulticenterBondDelta::Remove {
                        id: MulticenterBondId(id),
                        atoms: vec![AtomId(0), AtomId(1)],
                        attributes: MulticenterBondForm::default(),
                    })]),
                ),
                ApplyPreconditionError::InvalidReactionReference {
                    entity: Entity::MulticenterBond(MulticenterBondId(id)),
                },
            )
        }),
        (0u32..64).prop_map(|id| {
            (
                Reaction::new(
                    Molecule::default(),
                    Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                        id: NoncovalentBondId(id),
                        atoms: [AtomId(0), AtomId(1)],
                        attributes: NoncovalentBondForm::default(),
                    })]),
                ),
                ApplyPreconditionError::InvalidReactionReference {
                    entity: Entity::NoncovalentBond(NoncovalentBondId(id)),
                },
            )
        }),
        (0u32..64).prop_map(|id| {
            (
                Reaction::new(
                    Molecule::default(),
                    Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::Remove {
                        id: StereoAtomId(id),
                        site: AtomId(0),
                        ligands: vec![],
                        attributes: StereoAtomForm::default(),
                    })]),
                ),
                ApplyPreconditionError::InvalidReactionReference {
                    entity: Entity::StereoAtom(StereoAtomId(id)),
                },
            )
        }),
        (0u32..64).prop_map(|id| {
            (
                Reaction::new(
                    Molecule::default(),
                    Deltas::from_iter([Delta::StereoBond(StereoBondDelta::Remove {
                        id: StereoBondId(id),
                        site: BondId(0),
                        ligands: vec![],
                        attributes: StereoBondForm::default(),
                    })]),
                ),
                ApplyPreconditionError::InvalidReactionReference {
                    entity: Entity::StereoBond(StereoBondId(id)),
                },
            )
        }),
    ]
}

fn unavailable_participant_strategy() -> impl Strategy<Value = (Reaction, ApplyPreconditionError)> {
    (1u32..64).prop_flat_map(|missing| {
        let lhs = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            ..Default::default()
        });
        prop_oneof![
            Just((
                Reaction::new(
                    lhs.clone(),
                    Deltas::from_iter([Delta::Bond(BondDelta::Add {
                        id: BondId(0),
                        atoms: [AtomId(0), AtomId(missing)],
                        attributes: BondForm::default(),
                    })]),
                ),
                ApplyPreconditionError::InvalidReactionReference {
                    entity: Entity::Atom(AtomId(missing)),
                },
            )),
            Just((
                Reaction::new(
                    lhs.clone(),
                    Deltas::from_iter([Delta::DativeBond(DativeBondDelta::Add {
                        id: DativeBondId(0),
                        donors: vec![AtomId(0)],
                        acceptor: AtomId(missing),
                        attributes: DativeBondForm::default(),
                    })]),
                ),
                ApplyPreconditionError::InvalidReactionReference {
                    entity: Entity::Atom(AtomId(missing)),
                },
            )),
            Just((
                Reaction::new(
                    lhs.clone(),
                    Deltas::from_iter([Delta::AromaticSystem(AromaticSystemDelta::Add {
                        id: AromaticSystemId(0),
                        atoms: vec![AtomId(0), AtomId(missing)],
                        attributes: AromaticSystemForm::default(),
                    })]),
                ),
                ApplyPreconditionError::InvalidReactionReference {
                    entity: Entity::Atom(AtomId(missing)),
                },
            )),
            Just((
                Reaction::new(
                    lhs.clone(),
                    Deltas::from_iter([Delta::MulticenterBond(MulticenterBondDelta::Add {
                        id: MulticenterBondId(0),
                        atoms: vec![AtomId(0), AtomId(missing)],
                        attributes: MulticenterBondForm::default(),
                    })]),
                ),
                ApplyPreconditionError::InvalidReactionReference {
                    entity: Entity::Atom(AtomId(missing)),
                },
            )),
            Just((
                Reaction::new(
                    lhs.clone(),
                    Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::Add {
                        id: NoncovalentBondId(0),
                        atoms: [AtomId(0), AtomId(missing)],
                        attributes: NoncovalentBondForm::default(),
                    })]),
                ),
                ApplyPreconditionError::InvalidReactionReference {
                    entity: Entity::Atom(AtomId(missing)),
                },
            )),
            Just((
                Reaction::new(
                    lhs.clone(),
                    Deltas::from_iter([Delta::Constraint(ConstraintDelta::Add(
                        Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling {
                            atoms: Some(vec![AtomId(0), AtomId(missing)]),
                            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
                        }),
                    ))]),
                ),
                ApplyPreconditionError::InvalidReactionReference {
                    entity: Entity::Atom(AtomId(missing)),
                },
            )),
            Just((
                Reaction::new(
                    lhs.clone(),
                    Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::Add {
                        id: StereoAtomId(0),
                        site: AtomId(missing),
                        ligands: vec![],
                        attributes: StereoAtomForm::default(),
                    })]),
                ),
                ApplyPreconditionError::InvalidReactionReference {
                    entity: Entity::Atom(AtomId(missing)),
                },
            )),
            Just((
                Reaction::new(
                    lhs.clone(),
                    Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::Add {
                        id: StereoAtomId(0),
                        site: AtomId(0),
                        ligands: vec![StereoLigand::new(AtomId(missing), StereoLigandKind::Atom,)],
                        attributes: StereoAtomForm::default(),
                    })]),
                ),
                ApplyPreconditionError::InvalidReactionReference {
                    entity: Entity::Atom(AtomId(missing)),
                },
            )),
            Just((
                Reaction::new(
                    lhs,
                    Deltas::from_iter([Delta::StereoBond(StereoBondDelta::Add {
                        id: StereoBondId(0),
                        site: BondId(missing),
                        ligands: vec![],
                        attributes: StereoBondForm::default(),
                    })]),
                ),
                ApplyPreconditionError::InvalidReactionReference {
                    entity: Entity::Bond(BondId(missing)),
                },
            )),
        ]
    })
}

fn incompatible_incidence_strategy() -> impl Strategy<Value = (Reaction, ApplyPreconditionError)> {
    prop_oneof![
        Just((
            Reaction::new(
                Molecule::from_entries(MoleculeEntries {
                    atoms: vec![AtomForm::from_element(Element::C); 3],
                    dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondForm::default(),)],
                    ..Default::default()
                }),
                Deltas::from_iter([Delta::DativeBond(DativeBondDelta::Remove {
                    id: DativeBondId(0),
                    donors: vec![AtomId(0)],
                    acceptor: AtomId(2),
                    attributes: DativeBondForm::default(),
                })]),
            ),
            ApplyPreconditionError::ReactionIncidenceMismatch {
                entity: Entity::DativeBond(DativeBondId(0)),
            },
        )),
        Just((
            Reaction::new(
                Molecule::from_entries(MoleculeEntries {
                    atoms: vec![AtomForm::from_element(Element::C); 3],
                    aromatic: vec![(vec![AtomId(0), AtomId(1)], AromaticSystemForm::default(),)],
                    ..Default::default()
                }),
                Deltas::from_iter([Delta::AromaticSystem(AromaticSystemDelta::Remove {
                    id: AromaticSystemId(0),
                    atoms: vec![AtomId(0), AtomId(2)],
                    attributes: AromaticSystemForm::default(),
                })]),
            ),
            ApplyPreconditionError::ReactionIncidenceMismatch {
                entity: Entity::AromaticSystem(AromaticSystemId(0)),
            },
        )),
        Just((
            Reaction::new(
                Molecule::from_entries(MoleculeEntries {
                    atoms: vec![AtomForm::from_element(Element::C); 3],
                    multicenter: vec![
                        (vec![AtomId(0), AtomId(1)], MulticenterBondForm::default(),)
                    ],
                    ..Default::default()
                }),
                Deltas::from_iter([Delta::MulticenterBond(MulticenterBondDelta::Remove {
                    id: MulticenterBondId(0),
                    atoms: vec![AtomId(0), AtomId(2)],
                    attributes: MulticenterBondForm::default(),
                })]),
            ),
            ApplyPreconditionError::ReactionIncidenceMismatch {
                entity: Entity::MulticenterBond(MulticenterBondId(0)),
            },
        )),
        Just((
            Reaction::new(
                Molecule::from_entries(MoleculeEntries {
                    atoms: vec![AtomForm::from_element(Element::C); 3],
                    noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondForm::default(),)],
                    ..Default::default()
                }),
                Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(2)],
                    attributes: NoncovalentBondForm::default(),
                })]),
            ),
            ApplyPreconditionError::ReactionIncidenceMismatch {
                entity: Entity::NoncovalentBond(NoncovalentBondId(0)),
            },
        )),
        Just((
            Reaction::new(
                Molecule::from_entries(MoleculeEntries {
                    atoms: vec![AtomForm::from_element(Element::C); 5],
                    stereo_atoms: vec![(
                        AtomId(0),
                        (1..=4)
                            .map(|atom| { StereoLigand::new(AtomId(atom), StereoLigandKind::Atom) })
                            .collect(),
                        StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
                    )],
                    ..Default::default()
                }),
                Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::Remove {
                    id: StereoAtomId(0),
                    site: AtomId(0),
                    ligands: [2, 1, 3, 4]
                        .into_iter()
                        .map(|atom| { StereoLigand::new(AtomId(atom), StereoLigandKind::Atom) })
                        .collect(),
                    attributes: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
                })]),
            ),
            ApplyPreconditionError::ReactionIncidenceMismatch {
                entity: Entity::StereoAtom(StereoAtomId(0)),
            },
        )),
        Just((
            Reaction::new(
                Molecule::from_entries(MoleculeEntries {
                    atoms: vec![AtomForm::from_element(Element::C); 6],
                    bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))],
                    stereo_bonds: vec![(
                        BondId(0),
                        (2..=5)
                            .map(|atom| { StereoLigand::new(AtomId(atom), StereoLigandKind::Atom) })
                            .collect(),
                        StereoBondForm::new(StereoKind::CisTrans, 0u32),
                    )],
                    ..Default::default()
                }),
                Deltas::from_iter([Delta::StereoBond(StereoBondDelta::Remove {
                    id: StereoBondId(0),
                    site: BondId(0),
                    ligands: [3, 2, 4, 5]
                        .into_iter()
                        .map(|atom| { StereoLigand::new(AtomId(atom), StereoLigandKind::Atom) })
                        .collect(),
                    attributes: StereoBondForm::new(StereoKind::CisTrans, 0u32),
                })]),
            ),
            ApplyPreconditionError::ReactionIncidenceMismatch {
                entity: Entity::StereoBond(StereoBondId(0)),
            },
        )),
    ]
}

fn malformed_update_strategy() -> impl Strategy<Value = Reaction> {
    let discontinuous_field = (-16i64..=16, 1i64..=4, 1i64..=4, 1i64..=4).prop_map(
        |(old, first_step, gap, second_step)| {
            let first = old + first_step;
            let discontinuous_old = first + gap;
            let new = discontinuous_old + second_step;
            Reaction::new(
                Molecule::from_entries(MoleculeEntries {
                    atoms: vec![AtomForm::default().with_charge(old)],
                    ..Default::default()
                }),
                Deltas::from_iter([
                    Delta::Atom(AtomDelta::ModifyField {
                        id: AtomId(0),
                        change: AtomFieldChange::Charge {
                            old: NumForm::Lit(old),
                            new: NumForm::Lit(first),
                        },
                    }),
                    Delta::Atom(AtomDelta::ModifyField {
                        id: AtomId(0),
                        change: AtomFieldChange::Charge {
                            old: NumForm::Lit(discontinuous_old),
                            new: NumForm::Lit(new),
                        },
                    }),
                ]),
            )
        },
    );
    let stereo_atom = (stereo_atom_kind_strategy(), 0u32..16).prop_map(|(kind, offset)| {
        let ligands: Vec<StereoLigand> = (1..=kind.degree() as u32)
            .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
            .collect();
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); kind.degree() + 1],
                stereo_atoms: vec![(AtomId(0), ligands, StereoAtomForm::new(kind, 0u32))],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::ModifyField {
                id: StereoAtomId(0),
                change: StereoAtomFieldChange::Configuration {
                    old: StereoConfigurationForm::kinded(kind, 0u32),
                    new: StereoConfigurationForm::kinded(kind, kind.count() as u32 + offset),
                },
            })]),
        )
    });
    let stereo_bond = (0u32..16).prop_map(|offset| {
        let kind = StereoKind::CisTrans;
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); 6],
                bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))],
                stereo_bonds: vec![(
                    BondId(0),
                    (2..=5)
                        .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
                        .collect(),
                    StereoBondForm::new(kind, 0u32),
                )],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::StereoBond(StereoBondDelta::ModifyField {
                id: StereoBondId(0),
                change: StereoBondFieldChange::Configuration {
                    old: StereoConfigurationForm::kinded(kind, 0u32),
                    new: StereoConfigurationForm::kinded(kind, kind.count() as u32 + offset),
                },
            })]),
        )
    });

    prop_oneof![discontinuous_field, stereo_atom, stereo_bond]
}

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_reaction_validate_application_entity_reference_error(
        (reaction, expected) in unavailable_entity_strategy(),
    ) {
        prop_assert_eq!(reaction.validate_application(&Molecule::default()), Err(expected));
    }

    #[test]
    fn test_reaction_validate_application_participant_reference_error(
        (reaction, expected) in unavailable_participant_strategy(),
    ) {
        prop_assert_eq!(reaction.validate_application(&Molecule::default()), Err(expected));
    }

    #[test]
    fn test_reaction_validate_application_incidence_error(
        (reaction, expected) in incompatible_incidence_strategy(),
    ) {
        prop_assert_eq!(reaction.validate_application(&Molecule::default()), Err(expected));
    }

    #[test]
    fn test_reaction_validate_application_update_error(
        reaction in malformed_update_strategy(),
    ) {
        prop_assert_eq!(
            reaction.validate_application(&Molecule::default()),
            Err(ApplyPreconditionError::InconsistentReaction),
        );
    }

    #[test]
    fn test_reaction_apply_error(host_atom_count in 1usize..=8) {
        let constraint = Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(0)]),
            sum: NumForm::Lit(0),
        });
        let reaction = Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C)],
                constraints: Constraints::from(constraint.clone()),
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Constraint(ConstraintDelta::Remove(constraint))]),
        );
        let host = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); host_atom_count],
            ..Default::default()
        });
        let mut applications = reaction
            .apply(
                &host,
                SubstructureMatchConfig {
                    match_algorithm: SubstructureMatchAlgorithm::GraphAndOverlays,
                    subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::Vf2,
                    relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                },
            )
            .map_err(|error| TestCaseError::fail(format!("application precondition: {error:?}")))?;

        prop_assert_eq!(
            applications.next(),
            Some(Err(ApplyError::Transaction(TransactionError::MissingEntry))),
        );
        prop_assert_eq!(applications.next(), None);
        prop_assert_eq!(applications.next(), None);
    }
}
