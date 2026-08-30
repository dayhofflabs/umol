//! Reaction frame-transport and reframe properties.
//!
//! Comprehensive reactions exercise the action laws and the normalization/reframe pipeline over
//! the general well-formed domain. A separate generated domain pairs an owner-framed removal with
//! the same removal in a distinct compatible local frame for each of the six overlay kinds. The
//! latter checks the action laws in the coordinates where reaction transport must conjugate the
//! owner action, and checks that reframe converges with the owner-framed representation.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_chem::element::Element;
use umol_graph_ir::ir::{
    AromaticSystemDelta, AromaticSystemForm, AromaticSystemId, AtomForm, AtomId, BondForm, BondId,
    Contradiction, DativeBondDelta, DativeBondForm, DativeBondId, Delta, Deltas, FrameTransport,
    Molecule, MoleculeEntries, MulticenterBondDelta, MulticenterBondForm, MulticenterBondId,
    NoncovalentBondDelta, NoncovalentBondForm, NoncovalentBondId, NoncovalentBondKind, Normalize,
    Reaction, Reframe, StereoAtomDelta, StereoAtomForm, StereoAtomId, StereoBondDelta,
    StereoBondForm, StereoBondId, StereoKind, StereoLigand, StereoLigandKind,
};
use umol_perm::{DynPermutation, Permutation};

use crate::strategies::{
    comprehensive_reaction_strategy, intrinsic_contradiction_scenario_strategy,
    standardization_scenario_strategy,
};

const STEREO_BOND_FRAMES: [[u32; 4]; 8] = [
    [0, 1, 2, 3],
    [1, 0, 2, 3],
    [0, 1, 3, 2],
    [1, 0, 3, 2],
    [2, 3, 0, 1],
    [3, 2, 0, 1],
    [2, 3, 1, 0],
    [3, 2, 1, 0],
];

#[derive(Debug)]
struct RemovalFrameScenario {
    owner: Reaction,
    local: Reaction,
}

fn distinct_atom_frames_strategy(atoms: Vec<AtomId>) -> BoxedStrategy<(Vec<AtomId>, Vec<AtomId>)> {
    let degree = atoms.len();
    (Just(atoms).prop_shuffle(), 0..degree, 1..degree)
        .prop_map(|(owner, first, offset)| {
            let mut local = owner.clone();
            let second = (first + offset) % local.len();
            local.swap(first, second);
            (owner, local)
        })
        .boxed()
}

fn distinct_stereo_bond_frames_strategy() -> impl Strategy<Value = (Vec<AtomId>, Vec<AtomId>)> {
    (0..STEREO_BOND_FRAMES.len(), 1..STEREO_BOND_FRAMES.len()).prop_map(|(owner_index, offset)| {
        let local_index = (owner_index + offset) % STEREO_BOND_FRAMES.len();
        let owner = STEREO_BOND_FRAMES[owner_index]
            .map(AtomId)
            .into_iter()
            .collect();
        let local = STEREO_BOND_FRAMES[local_index]
            .map(AtomId)
            .into_iter()
            .collect();
        (owner, local)
    })
}

fn removal_frame_reaction_strategy() -> BoxedStrategy<RemovalFrameScenario> {
    let dative = distinct_atom_frames_strategy(vec![AtomId(0), AtomId(1), AtomId(2)]).prop_map(
        |(owner_donors, local_donors)| {
            let attributes = DativeBondForm::from_order(2);
            let lhs = Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); 4],
                dative: vec![(owner_donors.clone(), AtomId(3), attributes.clone())],
                ..Default::default()
            });
            RemovalFrameScenario {
                owner: Reaction::new(
                    lhs.clone(),
                    Deltas::from_iter([Delta::DativeBond(DativeBondDelta::Remove {
                        id: DativeBondId(0),
                        donors: owner_donors,
                        acceptor: AtomId(3),
                        attributes: attributes.clone(),
                    })]),
                ),
                local: Reaction::new(
                    lhs,
                    Deltas::from_iter([Delta::DativeBond(DativeBondDelta::Remove {
                        id: DativeBondId(0),
                        donors: local_donors,
                        acceptor: AtomId(3),
                        attributes,
                    })]),
                ),
            }
        },
    );

    let aromatic = distinct_atom_frames_strategy(vec![AtomId(0), AtomId(1), AtomId(2)]).prop_map(
        |(owner_atoms, local_atoms)| {
            let reference_atoms = vec![AtomId(0), AtomId(1), AtomId(2)];
            let reference_attributes = AromaticSystemForm::from_electrons(vec![1, 2, 3]);
            let owner_attributes = reference_attributes
                .clone()
                .reframe_by(
                    &DynPermutation::between(&reference_atoms, &owner_atoms)
                        .expect("the owner frame contains the reference atoms"),
                )
                .expect("the owner action has the form's degree");
            let local_attributes = reference_attributes
                .reframe_by(
                    &DynPermutation::between(&reference_atoms, &local_atoms)
                        .expect("the removal frame contains the reference atoms"),
                )
                .expect("the removal action has the form's degree");
            let lhs = Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); 3],
                aromatic: vec![(owner_atoms.clone(), owner_attributes.clone())],
                ..Default::default()
            });
            RemovalFrameScenario {
                owner: Reaction::new(
                    lhs.clone(),
                    Deltas::from_iter([Delta::AromaticSystem(AromaticSystemDelta::Remove {
                        id: AromaticSystemId(0),
                        atoms: owner_atoms,
                        attributes: owner_attributes,
                    })]),
                ),
                local: Reaction::new(
                    lhs,
                    Deltas::from_iter([Delta::AromaticSystem(AromaticSystemDelta::Remove {
                        id: AromaticSystemId(0),
                        atoms: local_atoms,
                        attributes: local_attributes,
                    })]),
                ),
            }
        },
    );

    let multicenter = distinct_atom_frames_strategy(vec![AtomId(0), AtomId(1), AtomId(2)])
        .prop_map(|(owner_atoms, local_atoms)| {
            let reference_atoms = vec![AtomId(0), AtomId(1), AtomId(2)];
            let reference_attributes = MulticenterBondForm::from_electrons(vec![3, 2, 1]);
            let owner_attributes = reference_attributes
                .clone()
                .reframe_by(
                    &DynPermutation::between(&reference_atoms, &owner_atoms)
                        .expect("the owner frame contains the reference atoms"),
                )
                .expect("the owner action has the form's degree");
            let local_attributes = reference_attributes
                .reframe_by(
                    &DynPermutation::between(&reference_atoms, &local_atoms)
                        .expect("the removal frame contains the reference atoms"),
                )
                .expect("the removal action has the form's degree");
            let lhs = Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); 3],
                multicenter: vec![(owner_atoms.clone(), owner_attributes.clone())],
                ..Default::default()
            });
            RemovalFrameScenario {
                owner: Reaction::new(
                    lhs.clone(),
                    Deltas::from_iter([Delta::MulticenterBond(MulticenterBondDelta::Remove {
                        id: MulticenterBondId(0),
                        atoms: owner_atoms,
                        attributes: owner_attributes,
                    })]),
                ),
                local: Reaction::new(
                    lhs,
                    Deltas::from_iter([Delta::MulticenterBond(MulticenterBondDelta::Remove {
                        id: MulticenterBondId(0),
                        atoms: local_atoms,
                        attributes: local_attributes,
                    })]),
                ),
            }
        });

    let noncovalent = distinct_atom_frames_strategy(vec![AtomId(0), AtomId(1)]).prop_map(
        |(owner_atoms, local_atoms)| {
            let attributes = NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond);
            let lhs = Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); 2],
                noncovalent: vec![(owner_atoms[0], owner_atoms[1], attributes.clone())],
                ..Default::default()
            });
            RemovalFrameScenario {
                owner: Reaction::new(
                    lhs.clone(),
                    Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                        id: NoncovalentBondId(0),
                        atoms: [owner_atoms[0], owner_atoms[1]],
                        attributes: attributes.clone(),
                    })]),
                ),
                local: Reaction::new(
                    lhs,
                    Deltas::from_iter([Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                        id: NoncovalentBondId(0),
                        atoms: [local_atoms[0], local_atoms[1]],
                        attributes,
                    })]),
                ),
            }
        },
    );

    let stereo_atom =
        distinct_atom_frames_strategy(vec![AtomId(1), AtomId(2), AtomId(3), AtomId(4)]).prop_map(
            |(owner_atoms, local_atoms)| {
                let reference_atoms = vec![AtomId(1), AtomId(2), AtomId(3), AtomId(4)];
                let owner_ligands = owner_atoms
                    .iter()
                    .copied()
                    .map(|atom| StereoLigand::new(atom, StereoLigandKind::Atom))
                    .collect::<Vec<_>>();
                let local_ligands = local_atoms
                    .iter()
                    .copied()
                    .map(|atom| StereoLigand::new(atom, StereoLigandKind::Atom))
                    .collect::<Vec<_>>();
                let reference_attributes = StereoAtomForm::new(StereoKind::Tetrahedral, 0u32);
                let owner_attributes = reference_attributes
                    .clone()
                    .reframe_by(
                        &Permutation::between(&reference_atoms, &owner_atoms)
                            .expect("the owner frame contains the reference atoms"),
                    )
                    .expect("the owner action is tetrahedral");
                let local_attributes = reference_attributes
                    .reframe_by(
                        &Permutation::between(&reference_atoms, &local_atoms)
                            .expect("the removal frame contains the reference atoms"),
                    )
                    .expect("the removal action is tetrahedral");
                let lhs = Molecule::from_entries(MoleculeEntries {
                    atoms: vec![AtomForm::from_element(Element::C); 5],
                    bonds: (1..=4)
                        .map(|ligand| (AtomId(0), AtomId(ligand), BondForm::from_order(1)))
                        .collect(),
                    stereo_atoms: vec![(
                        AtomId(0),
                        owner_ligands.clone(),
                        owner_attributes.clone(),
                    )],
                    ..Default::default()
                });
                RemovalFrameScenario {
                    owner: Reaction::new(
                        lhs.clone(),
                        Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::Remove {
                            id: StereoAtomId(0),
                            site: AtomId(0),
                            ligands: owner_ligands,
                            attributes: owner_attributes,
                        })]),
                    ),
                    local: Reaction::new(
                        lhs,
                        Deltas::from_iter([Delta::StereoAtom(StereoAtomDelta::Remove {
                            id: StereoAtomId(0),
                            site: AtomId(0),
                            ligands: local_ligands,
                            attributes: local_attributes,
                        })]),
                    ),
                }
            },
        );

    let stereo_bond =
        distinct_stereo_bond_frames_strategy().prop_map(|(owner_atoms, local_atoms)| {
            let reference_atoms = vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)];
            let owner_ligands = owner_atoms
                .iter()
                .copied()
                .map(|atom| StereoLigand::new(atom, StereoLigandKind::Atom))
                .collect::<Vec<_>>();
            let local_ligands = local_atoms
                .iter()
                .copied()
                .map(|atom| StereoLigand::new(atom, StereoLigandKind::Atom))
                .collect::<Vec<_>>();
            let reference_attributes = StereoBondForm::new(StereoKind::CisTrans, 0u32);
            let owner_attributes = reference_attributes
                .clone()
                .reframe_by(
                    &Permutation::between(&reference_atoms, &owner_atoms)
                        .expect("the owner frame contains the reference atoms"),
                )
                .expect("the owner action preserves the stereo-bond endpoint blocks");
            let local_attributes = reference_attributes
                .reframe_by(
                    &Permutation::between(&reference_atoms, &local_atoms)
                        .expect("the removal frame contains the reference atoms"),
                )
                .expect("the removal action preserves the stereo-bond endpoint blocks");
            let lhs = Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); 6],
                bonds: vec![
                    (AtomId(4), AtomId(5), BondForm::from_order(2)),
                    (AtomId(4), AtomId(0), BondForm::from_order(1)),
                    (AtomId(4), AtomId(1), BondForm::from_order(1)),
                    (AtomId(5), AtomId(2), BondForm::from_order(1)),
                    (AtomId(5), AtomId(3), BondForm::from_order(1)),
                ],
                stereo_bonds: vec![(BondId(0), owner_ligands.clone(), owner_attributes.clone())],
                ..Default::default()
            });
            RemovalFrameScenario {
                owner: Reaction::new(
                    lhs.clone(),
                    Deltas::from_iter([Delta::StereoBond(StereoBondDelta::Remove {
                        id: StereoBondId(0),
                        site: BondId(0),
                        ligands: owner_ligands,
                        attributes: owner_attributes,
                    })]),
                ),
                local: Reaction::new(
                    lhs,
                    Deltas::from_iter([Delta::StereoBond(StereoBondDelta::Remove {
                        id: StereoBondId(0),
                        site: BondId(0),
                        ligands: local_ligands,
                        attributes: local_attributes,
                    })]),
                ),
            }
        });

    prop_oneof![
        dative,
        aromatic,
        multicenter,
        noncovalent,
        stereo_atom,
        stereo_bond,
    ]
    .boxed()
}

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_reaction_reframe_by(reaction in comprehensive_reaction_strategy()) {
        let action = reaction.representative_action();
        let identity = action.identity();
        let inverse = action.inverse();
        let composite = action
            .compose(&inverse)
            .expect("an action and its inverse have the same domain");

        prop_assert_eq!(
            reaction.clone().reframe_by(&identity),
            Some(reaction.clone()),
        );
        prop_assert_eq!(
            reaction
                .clone()
                .reframe_by(&action)
                .and_then(|transported| transported.reframe_by(&inverse)),
            reaction.clone().reframe_by(&composite),
        );
        prop_assert_eq!(reaction.clone().reframe_by(&composite), Some(reaction));
    }

    #[test]
    fn test_reaction_reframe_by_removal(scenario in removal_frame_reaction_strategy()) {
        let reaction = scenario.local;
        let action = reaction.representative_action();
        let identity = action.identity();
        let inverse = action.inverse();
        let composite = action
            .compose(&inverse)
            .expect("an action and its inverse have the same domain");

        prop_assert_eq!(
            reaction.clone().reframe_by(&identity),
            Some(reaction.clone()),
        );
        prop_assert_eq!(
            reaction
                .clone()
                .reframe_by(&action)
                .and_then(|transported| transported.reframe_by(&inverse)),
            reaction.clone().reframe_by(&composite),
        );
        prop_assert_eq!(reaction.clone().reframe_by(&composite), Some(reaction));
    }

    #[test]
    fn test_reaction_representative_action_erased_entity(
        scenario in standardization_scenario_strategy(),
    ) {
        let normalized = scenario.reaction.clone().normalize().map_err(|_| {
            TestCaseError::fail("generated reaction is intrinsically contradictory")
        })?;
        let (reframed, action) = scenario.reaction.reframe_with_action().map_err(|_| {
            TestCaseError::fail("generated reaction is intrinsically contradictory")
        })?;
        let transported = normalized
            .clone()
            .reframe_by(&action)
            .ok_or_else(|| TestCaseError::fail("input-domain action did not cover its source"))?
            .normalize()
            .map_err(|_| TestCaseError::fail("transported reaction became contradictory"))?;

        prop_assert!(action.multicenter_bonds().contains(MulticenterBondId(7)));
        prop_assert!(!normalized
            .representative_action()
            .multicenter_bonds()
            .contains(MulticenterBondId(7)));
        prop_assert_eq!(transported, reframed);
    }

    #[test]
    fn test_reaction_representative_action_contradiction(
        scenario in intrinsic_contradiction_scenario_strategy(),
    ) {
        for reaction in scenario.reactions {
            let action = reaction.representative_action();

            prop_assert_eq!(
                action.compose(&action.identity()),
                Some(action.clone()),
            );
            prop_assert_eq!(reaction.clone().normalize(), Err(Contradiction));
            prop_assert_eq!(reaction.reframe(), Err(Contradiction));
        }
    }

    #[test]
    fn test_reaction_reframe_with_action(reaction in comprehensive_reaction_strategy()) {
        let fused = reaction.clone().reframe().map_err(|_| {
            TestCaseError::fail("generated reaction is intrinsically contradictory")
        })?;
        let (witnessed, action) = reaction.clone().reframe_with_action().map_err(|_| {
            TestCaseError::fail("generated reaction is intrinsically contradictory")
        })?;
        let transported = reaction
            .normalize()
            .map_err(|_| {
                TestCaseError::fail("generated reaction is intrinsically contradictory")
            })?
            .reframe_by(&action)
            .ok_or_else(|| TestCaseError::fail("representative action did not cover its source"))?
            .normalize()
            .map_err(|_| {
                TestCaseError::fail("transported reaction is intrinsically contradictory")
            })?;
        let selected_action = witnessed.representative_action();

        prop_assert_eq!(fused, witnessed.clone());
        prop_assert_eq!(transported, witnessed);
        prop_assert_eq!(selected_action.clone(), selected_action.identity());
    }

    #[test]
    fn test_reaction_reframe_removal(scenario in removal_frame_reaction_strategy()) {
        prop_assert_eq!(scenario.local.reframe(), scenario.owner.reframe());
    }

    #[test]
    fn test_reaction_framed_eq(reaction in comprehensive_reaction_strategy()) {
        let normalized = reaction.clone().normalize().map_err(|_| {
            TestCaseError::fail("generated reaction is intrinsically contradictory")
        })?;
        let reframed = reaction.clone().reframe().map_err(|_| {
            TestCaseError::fail("generated reaction is intrinsically contradictory")
        })?;

        prop_assert!(reaction.normalized_eq(&normalized));
        prop_assert!(reaction.framed_eq(&normalized));
        prop_assert!(reaction.framed_eq(&reframed));
    }
}
