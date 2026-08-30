//! Property tests for malformed reaction application inputs.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_core::{RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm};
use umol_graph_ir::ir::{
    ApplyError, ApplyPreconditionError, AromaticSystemDelta, AtomDelta, BondDelta, ConstraintDelta,
    Contradiction, DativeBondDelta, Entity, MulticenterBondDelta, NoncovalentBondDelta, React,
    ReactionIntegrityError, StereoAtomDelta, StereoBondDelta, SubstructureMatchAlgorithm,
    SubstructureMatchConfig, TransactionError,
};

use crate::strategies::*;

const MATCH_CONFIG: SubstructureMatchConfig = SubstructureMatchConfig {
    match_algorithm: SubstructureMatchAlgorithm::GraphAndOverlays,
    subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::Vf2,
    relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
};

fn unavailable_entity_strategy() -> impl Strategy<
    Value = (
        Result<Reaction, ReactionIntegrityError>,
        ReactionIntegrityError,
    ),
> {
    prop_oneof![
        (0u32..64).prop_map(|id| {
            (
                Reaction::try_new(
                    Molecule::default(),
                    Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                        id: AtomId(id),
                        attributes: AtomForm::default(),
                    })]),
                ),
                ReactionIntegrityError::InvalidReference {
                    entity: Entity::Atom(AtomId(id)),
                },
            )
        }),
        (0u32..64).prop_map(|id| {
            (
                Reaction::try_new(
                    Molecule::default(),
                    Deltas::from_iter([Delta::Bond(BondDelta::Remove {
                        id: BondId(id),
                        atoms: [AtomId(0), AtomId(1)],
                        attributes: BondForm::default(),
                    })]),
                ),
                ReactionIntegrityError::InvalidReference {
                    entity: Entity::Bond(BondId(id)),
                },
            )
        }),
        (0u32..64).prop_map(|id| {
            (
                Reaction::try_new(
                    Molecule::default(),
                    Deltas::from_iter([Delta::DativeBond(DativeBondDelta::Remove {
                        id: DativeBondId(id),
                        donors: vec![AtomId(0)],
                        acceptor: AtomId(1),
                        attributes: DativeBondForm::default(),
                    })]),
                ),
                ReactionIntegrityError::InvalidReference {
                    entity: Entity::DativeBond(DativeBondId(id)),
                },
            )
        }),
        (0u32..64).prop_map(|id| {
            (
                Reaction::try_new(
                    Molecule::default(),
                    Deltas::from_iter([Delta::AromaticSystem(AromaticSystemDelta::Remove {
                        id: AromaticSystemId(id),
                        atoms: vec![AtomId(0), AtomId(1)],
                        attributes: AromaticSystemForm::default(),
                    })]),
                ),
                ReactionIntegrityError::InvalidReference {
                    entity: Entity::AromaticSystem(AromaticSystemId(id)),
                },
            )
        }),
        (0u32..64).prop_map(|id| {
            (
                Reaction::try_new(
                    Molecule::default(),
                    Deltas::from_iter([Delta::MulticenterBond(MulticenterBondDelta::Remove {
                        id: MulticenterBondId(id),
                        atoms: vec![AtomId(0), AtomId(1)],
                        attributes: MulticenterBondForm::default(),
                    })]),
                ),
                ReactionIntegrityError::InvalidReference {
                    entity: Entity::MulticenterBond(MulticenterBondId(id)),
                },
            )
        }),
        (0u32..64).prop_map(|id| {
            (
                Reaction::try_new(
                    Molecule::default(),
                    Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                        id: NoncovalentBondId(id),
                        atoms: [AtomId(0), AtomId(1)],
                        attributes: NoncovalentBondForm::default(),
                    })]),
                ),
                ReactionIntegrityError::InvalidReference {
                    entity: Entity::NoncovalentBond(NoncovalentBondId(id)),
                },
            )
        }),
        (0u32..64).prop_map(|id| {
            (
                Reaction::try_new(
                    Molecule::default(),
                    Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::Remove {
                        id: StereoAtomId(id),
                        site: AtomId(0),
                        ligands: vec![],
                        attributes: StereoAtomForm::default(),
                    })]),
                ),
                ReactionIntegrityError::InvalidReference {
                    entity: Entity::StereoAtom(StereoAtomId(id)),
                },
            )
        }),
        (0u32..64).prop_map(|id| {
            (
                Reaction::try_new(
                    Molecule::default(),
                    Deltas::from_iter([Delta::StereoBond(StereoBondDelta::Remove {
                        id: StereoBondId(id),
                        site: BondId(0),
                        ligands: vec![],
                        attributes: StereoBondForm::default(),
                    })]),
                ),
                ReactionIntegrityError::InvalidReference {
                    entity: Entity::StereoBond(StereoBondId(id)),
                },
            )
        }),
    ]
}

fn unavailable_participant_strategy() -> impl Strategy<
    Value = (
        Result<Reaction, ReactionIntegrityError>,
        ReactionIntegrityError,
    ),
> {
    (1u32..64).prop_flat_map(|missing| {
        let lhs = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            ..Default::default()
        });
        prop_oneof![
            Just((
                Reaction::try_new(
                    lhs.clone(),
                    Deltas::from_iter([Delta::Bond(BondDelta::Add {
                        id: BondId(0),
                        atoms: [AtomId(0), AtomId(missing)],
                        attributes: BondForm::default(),
                    })]),
                ),
                ReactionIntegrityError::InvalidReference {
                    entity: Entity::Atom(AtomId(missing)),
                },
            )),
            Just((
                Reaction::try_new(
                    lhs.clone(),
                    Deltas::from_iter([Delta::DativeBond(DativeBondDelta::Add {
                        id: DativeBondId(0),
                        donors: vec![AtomId(0)],
                        acceptor: AtomId(missing),
                        attributes: DativeBondForm::default(),
                    })]),
                ),
                ReactionIntegrityError::InvalidReference {
                    entity: Entity::Atom(AtomId(missing)),
                },
            )),
            Just((
                Reaction::try_new(
                    lhs.clone(),
                    Deltas::from_iter([Delta::AromaticSystem(AromaticSystemDelta::Add {
                        id: AromaticSystemId(0),
                        atoms: vec![AtomId(0), AtomId(missing)],
                        attributes: AromaticSystemForm::default(),
                    })]),
                ),
                ReactionIntegrityError::InvalidReference {
                    entity: Entity::Atom(AtomId(missing)),
                },
            )),
            Just((
                Reaction::try_new(
                    lhs.clone(),
                    Deltas::from_iter([Delta::MulticenterBond(MulticenterBondDelta::Add {
                        id: MulticenterBondId(0),
                        atoms: vec![AtomId(0), AtomId(missing)],
                        attributes: MulticenterBondForm::default(),
                    })]),
                ),
                ReactionIntegrityError::InvalidReference {
                    entity: Entity::Atom(AtomId(missing)),
                },
            )),
            Just((
                Reaction::try_new(
                    lhs.clone(),
                    Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::Add {
                        id: NoncovalentBondId(0),
                        atoms: [AtomId(0), AtomId(missing)],
                        attributes: NoncovalentBondForm::default(),
                    })]),
                ),
                ReactionIntegrityError::InvalidReference {
                    entity: Entity::Atom(AtomId(missing)),
                },
            )),
            Just((
                Reaction::try_new(
                    lhs.clone(),
                    Deltas::from_iter([Delta::Constraint(ConstraintDelta::Add(
                        Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling {
                            atoms: Some(vec![AtomId(0), AtomId(missing)]),
                            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
                        }),
                    ))]),
                ),
                ReactionIntegrityError::InvalidReference {
                    entity: Entity::Atom(AtomId(missing)),
                },
            )),
            Just((
                Reaction::try_new(
                    lhs.clone(),
                    Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::Add {
                        id: StereoAtomId(0),
                        site: AtomId(missing),
                        ligands: vec![],
                        attributes: StereoAtomForm::default(),
                    })]),
                ),
                ReactionIntegrityError::InvalidReference {
                    entity: Entity::Atom(AtomId(missing)),
                },
            )),
            Just((
                Reaction::try_new(
                    lhs.clone(),
                    Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::Add {
                        id: StereoAtomId(0),
                        site: AtomId(0),
                        ligands: vec![StereoLigand::new(AtomId(missing), StereoLigandKind::Atom,)],
                        attributes: StereoAtomForm::default(),
                    })]),
                ),
                ReactionIntegrityError::InvalidReference {
                    entity: Entity::Atom(AtomId(missing)),
                },
            )),
            Just((
                Reaction::try_new(
                    lhs,
                    Deltas::from_iter([Delta::StereoBond(StereoBondDelta::Add {
                        id: StereoBondId(0),
                        site: BondId(missing),
                        ligands: vec![],
                        attributes: StereoBondForm::default(),
                    })]),
                ),
                ReactionIntegrityError::InvalidReference {
                    entity: Entity::Bond(BondId(missing)),
                },
            )),
        ]
    })
}

fn incompatible_incidence_strategy() -> impl Strategy<
    Value = (
        Result<Reaction, ReactionIntegrityError>,
        ReactionIntegrityError,
    ),
> {
    prop_oneof![
        Just((
            Reaction::try_new(
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
            ReactionIntegrityError::IncidenceMismatch {
                entity: Entity::DativeBond(DativeBondId(0)),
            },
        )),
        Just((
            Reaction::try_new(
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
            ReactionIntegrityError::IncidenceMismatch {
                entity: Entity::AromaticSystem(AromaticSystemId(0)),
            },
        )),
        Just((
            Reaction::try_new(
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
            ReactionIntegrityError::IncidenceMismatch {
                entity: Entity::MulticenterBond(MulticenterBondId(0)),
            },
        )),
        Just((
            Reaction::try_new(
                Molecule::from_entries(MoleculeEntries {
                    atoms: vec![AtomForm::from_element(Element::C); 3],
                    noncovalent: vec![([AtomId(0), AtomId(1)], NoncovalentBondForm::default(),)],
                    ..Default::default()
                }),
                Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                    id: NoncovalentBondId(0),
                    atoms: [AtomId(0), AtomId(2)],
                    attributes: NoncovalentBondForm::default(),
                })]),
            ),
            ReactionIntegrityError::IncidenceMismatch {
                entity: Entity::NoncovalentBond(NoncovalentBondId(0)),
            },
        )),
        Just((
            Reaction::try_new(
                Molecule::from_entries(MoleculeEntries {
                    atoms: vec![AtomForm::from_element(Element::C); 6],
                    bonds: vec![
                        (AtomId(0), AtomId(1), BondForm::from_order(2)),
                        (AtomId(0), AtomId(2), BondForm::from_order(1)),
                        (AtomId(0), AtomId(3), BondForm::from_order(1)),
                        (AtomId(1), AtomId(4), BondForm::from_order(1)),
                        (AtomId(1), AtomId(5), BondForm::from_order(1)),
                    ],
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
                    ligands: [2, 4, 3, 5]
                        .into_iter()
                        .map(|atom| { StereoLigand::new(AtomId(atom), StereoLigandKind::Atom) })
                        .collect(),
                    attributes: StereoBondForm::new(StereoKind::CisTrans, 0u32),
                })]),
            ),
            ReactionIntegrityError::IncidenceMismatch {
                entity: Entity::StereoBond(StereoBondId(0)),
            },
        )),
    ]
}

fn malformed_update_strategy() -> impl Strategy<Value = Reaction> {
    let stereo_atom = (stereo_atom_kind_strategy(), 0u32..16).prop_map(|(kind, offset)| {
        let ligands: Vec<StereoLigand> = (1..=kind.degree() as u32)
            .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
            .collect();
        Reaction::try_new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); kind.degree() + 1],
                bonds: (1..=kind.degree() as u32)
                    .map(|atom| (AtomId(0), AtomId(atom), BondForm::from_order(1)))
                    .collect(),
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
        .expect("generated stereo-atom update is representation-valid")
    });
    let stereo_bond = (0u32..16).prop_map(|offset| {
        let kind = StereoKind::CisTrans;
        Reaction::try_new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); 6],
                bonds: vec![
                    (AtomId(0), AtomId(1), BondForm::from_order(2)),
                    (AtomId(0), AtomId(2), BondForm::from_order(1)),
                    (AtomId(0), AtomId(3), BondForm::from_order(1)),
                    (AtomId(1), AtomId(4), BondForm::from_order(1)),
                    (AtomId(1), AtomId(5), BondForm::from_order(1)),
                ],
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
        .expect("generated stereo-bond update is representation-valid")
    });

    prop_oneof![
        discontinuous_atom_update_reaction_strategy(),
        stereo_atom,
        stereo_bond
    ]
}

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_reaction_try_new_entity_reference_error(
        (actual, expected) in unavailable_entity_strategy(),
    ) {
        prop_assert_eq!(actual, Err(expected));
    }

    #[test]
    fn test_reaction_try_new_participant_reference_error(
        (actual, expected) in unavailable_participant_strategy(),
    ) {
        prop_assert_eq!(actual, Err(expected));
    }

    #[test]
    fn test_reaction_try_new_incidence_error(
        (actual, expected) in incompatible_incidence_strategy(),
    ) {
        prop_assert_eq!(actual, Err(expected));
    }

    #[test]
    fn test_reaction_check_preconditions_update_error(
        reaction in malformed_update_strategy(),
    ) {
        prop_assert_eq!(
            reaction.check_preconditions(),
            Err(ApplyPreconditionError::InconsistentReaction),
        );
    }

    /// One discontinuous field-update chain is rejected at the reaction/span materialization
    /// boundary with the exact normalization error.
    #[test]
    fn test_reaction_to_reaction_span_error(
        reaction in discontinuous_atom_update_reaction_strategy(),
    ) {
        prop_assert_eq!(reaction.to_reaction_span(), Err(Contradiction));
    }

    #[test]
    #[ignore = "re-enable when matching evaluates molecule-scope pattern constraints"]
    fn test_reaction_apply_error(host_atom_count in 1usize..=8) {
        let constraint = Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(0)]),
            sum: NumForm::Lit(0),
        });
        let reaction = Reaction::try_new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C)],
                constraints: Constraints::from(constraint.clone()),
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Constraint(ConstraintDelta::Remove(constraint))]),
        )
        .expect("generated constraint removal is representation-valid");
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

        let mut products = host
            .react(&reaction, MATCH_CONFIG)
            .map_err(|error| TestCaseError::fail(format!("application precondition: {error:?}")))?;

        prop_assert_eq!(
            products.next(),
            Some(Err(ApplyError::Transaction(TransactionError::MissingEntry))),
        );
        prop_assert_eq!(products.next(), None);
        prop_assert_eq!(products.next(), None);
    }
}
