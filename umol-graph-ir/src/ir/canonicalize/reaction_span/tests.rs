//! Unit tests for reaction-span and reaction canonicalization.

use rstest::{fixture, rstest};
use umol_chem::element::Element;

use super::super::tests::{canonicalize_context, initial_color_molecule};
use super::super::*;
use super::*;
use crate::ir::{
    AromaticSystemForm, AromaticSystemId, AtomFieldChange, AtomForm, AtomId, BondId,
    ConstraintDelta, DativeBondForm, DativeBondId, Deltas, Entity, IncidenceLevel,
    MoleculeCorrespondence, MoleculeEntries, MulticenterBondForm, MulticenterBondId,
    NoncovalentBondForm, NoncovalentBondId, ReactionSpanEntries, StereoAtomForm, StereoAtomId,
    StereoBondForm, StereoBondId, StereoKind, StereoLigand, StereoLigandKind, TopicityForm,
};

#[rstest]
#[case::topology(
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C).with_charge(0_i64),
                AtomForm::from_element(Element::O),
            ],
            ..Default::default()
        }),
        [Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(0),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(0),
                new: NumForm::Lit(1),
            },
        })]
        .into_iter()
        .collect(),
    ),
    DescriptionLevel::Topology,
)]
#[case::constitution(
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::O),
                AtomForm::from_element(Element::H),
            ],
            ..Default::default()
        }),
        [Delta::NoncovalentBond(NoncovalentBondDelta::Add {
            id: NoncovalentBondId(0),
            atoms: [AtomId(0), AtomId(1)],
            attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        })]
        .into_iter()
        .collect(),
    ),
    DescriptionLevel::Constitution,
)]
#[case::structure(
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::H),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
                AtomForm::from_element(Element::Br),
            ],
            bonds: (1..=4)
                .map(|id| (AtomId(0), AtomId(id), BondForm::from_order(1)))
                .collect(),
            ..Default::default()
        }),
        [Delta::StereoAtom(StereoAtomDelta::Add {
            id: StereoAtomId(0),
            site: AtomId(0),
            ligands: vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            attributes: StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
        })]
        .into_iter()
        .collect(),
    ),
    DescriptionLevel::Structure,
)]
#[case::full(
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            ..Default::default()
        }),
        [Delta::Constraint(ConstraintDelta::Add(
            Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }),
        ))]
        .into_iter()
        .collect(),
    ),
    DescriptionLevel::Full,
)]
fn test_reaction_canonicalize_level(
    #[case] reaction: Reaction,
    #[case] expected: DescriptionLevel,
) {
    assert_eq!(reaction_canonicalize_level(&reaction), expected);
}

#[rstest]
#[case::topology(
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C).with_charge(0_i64),
                AtomForm::from_element(Element::O),
            ],
            ..Default::default()
        }),
        [Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(0),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(0),
                new: NumForm::Lit(1),
            },
        })]
        .into_iter()
        .collect(),
    ),
    DescriptionLevel::Topology,
)]
#[case::constitution(
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::O),
                AtomForm::from_element(Element::H),
            ],
            ..Default::default()
        }),
        [Delta::NoncovalentBond(NoncovalentBondDelta::Add {
            id: NoncovalentBondId(0),
            atoms: [AtomId(0), AtomId(1)],
            attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        })]
        .into_iter()
        .collect(),
    ),
    DescriptionLevel::Constitution,
)]
#[case::structure(
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::H),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
                AtomForm::from_element(Element::Br),
            ],
            bonds: (1..=4)
                .map(|id| (AtomId(0), AtomId(id), BondForm::from_order(1)))
                .collect(),
            ..Default::default()
        }),
        [Delta::StereoAtom(StereoAtomDelta::Add {
            id: StereoAtomId(0),
            site: AtomId(0),
            ligands: vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            attributes: StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
        })]
        .into_iter()
        .collect(),
    ),
    DescriptionLevel::Structure,
)]
#[case::full(
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            ..Default::default()
        }),
        [Delta::Constraint(ConstraintDelta::Add(
            Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }),
        ))]
        .into_iter()
        .collect(),
    ),
    DescriptionLevel::Full,
)]
fn test_reaction_span_canonicalize_level(
    #[case] reaction: Reaction,
    #[case] expected: DescriptionLevel,
) {
    let span = reaction
        .to_reaction_span()
        .expect("case reaction materializes");

    assert_eq!(reaction_span_canonicalize_level(&span), expected);
}

#[rstest]
#[case::dative(Entity::DativeBond(DativeBondId(0)))]
#[case::aromatic(Entity::AromaticSystem(AromaticSystemId(0)))]
#[case::multicenter(Entity::MulticenterBond(MulticenterBondId(0)))]
#[case::noncovalent(Entity::NoncovalentBond(NoncovalentBondId(0)))]
#[case::stereo_atom(Entity::StereoAtom(StereoAtomId(0)))]
#[case::stereo_bond(Entity::StereoBond(StereoBondId(0)))]
fn test_remapping_from_order(initial_color_molecule: Molecule, #[case] excluded: Entity) {
    let span = Reaction::new(initial_color_molecule, Deltas::new())
        .to_reaction_span()
        .expect("identity reaction materializes");
    let incidence_graph = span.incidence_graph(IncidenceLevel::Topology);
    let order = incidence_graph.graph().node_ids().collect::<Vec<_>>();
    let remapping = remapping_from_order(&span, &incidence_graph, &order);

    assert_eq!(
        MoleculeCorrespondence::from(&remapping).right_of(excluded),
        Some(excluded)
    );
}

#[rstest]
#[case::identity(
    [0, 1, 2, 3, 4].map(NodeId).into(),
    molecule_remapping(&[
        vec![0, 1, 2],
        vec![0, 1],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    ])
)]
#[case::path_reversal(
    [2, 1, 0, 4, 3].map(NodeId).into(),
    molecule_remapping(&[
        vec![2, 1, 0],
        vec![1, 0],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    ])
)]
fn test_remapping_from_generator(
    #[case] generator: Vec<NodeId>,
    #[case] expected: MoleculeRemapping,
) {
    let span = ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 3],
        bonds: vec![
            (
                AtomId(0),
                AtomId(1),
                EntitySpan::Unchanged(BondForm::from_order(1)),
            ),
            (
                AtomId(1),
                AtomId(2),
                EntitySpan::Unchanged(BondForm::from_order(1)),
            ),
        ],
        ..Default::default()
    });
    let incidence_graph = span.incidence_graph(IncidenceLevel::Full);

    assert_eq!(
        remapping_from_generator(&span, &incidence_graph, &generator),
        expected
    );
}

#[rstest]
#[case::feature_free(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 2],
        bonds: vec![(AtomId(0), AtomId(1), EntitySpan::Unchanged(BondForm::from_order(1)))],
        ..Default::default()
    }),
    vec![NodeId(1), NodeId(0), NodeId(2)],
    true,
)]
#[case::distinct_ligands(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        bonds: (1..=4)
            .map(|ligand| {
                (AtomId(0), AtomId(ligand), EntitySpan::Unchanged(BondForm::from_order(1)))
            })
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            (1..=4)
                .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                .collect(),
            EntitySpan::Unchanged(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32)),
        )],
        ..Default::default()
    }),
    (0..10).map(NodeId::from).collect(),
    true,
)]
#[case::prochiral_explicit_hydrogens(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::H, Element::H, Element::F, Element::Cl]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        bonds: (1..=4)
            .map(|ligand| {
                (AtomId(0), AtomId(ligand), EntitySpan::Unchanged(BondForm::from_order(1)))
            })
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            (1..=4)
                .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                .collect(),
            EntitySpan::Unchanged(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32)),
        )],
        ..Default::default()
    }),
    vec![
        NodeId(0), NodeId(2), NodeId(1), NodeId(3), NodeId(4),
        NodeId(6), NodeId(5), NodeId(7), NodeId(8), NodeId(9),
    ],
    false,
)]
#[case::stereo_site_exchange(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [
            Element::C, Element::F, Element::Cl, Element::Br, Element::I,
            Element::C, Element::F, Element::Cl, Element::Br, Element::I,
        ]
        .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
        .into(),
        bonds: [(0, 1), (0, 2), (0, 3), (0, 4), (5, 6), (5, 7), (5, 8), (5, 9)]
            .map(|(site, ligand)| {
                (AtomId(site), AtomId(ligand), EntitySpan::Unchanged(BondForm::from_order(1)))
            })
            .into(),
        stereo_atoms: vec![
            (
                AtomId(0),
                [1, 2, 3, 4].map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).into(),
                EntitySpan::Unchanged(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32)),
            ),
            (
                AtomId(5),
                [6, 7, 8, 9].map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).into(),
                EntitySpan::Unchanged(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32)),
            ),
        ],
        ..Default::default()
    }),
    vec![
        NodeId(5), NodeId(6), NodeId(7), NodeId(8), NodeId(9),
        NodeId(0), NodeId(1), NodeId(2), NodeId(3), NodeId(4),
        NodeId(14), NodeId(15), NodeId(16), NodeId(17),
        NodeId(10), NodeId(11), NodeId(12), NodeId(13),
        NodeId(19), NodeId(18),
    ],
    true,
)]
#[case::stereo_bond_endpoint_exchange(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::C, Element::F, Element::Cl, Element::F, Element::Cl]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Unchanged(BondForm::from_order(2))),
            (AtomId(0), AtomId(2), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(0), AtomId(3), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(1), AtomId(4), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(1), AtomId(5), EntitySpan::Unchanged(BondForm::from_order(1))),
        ],
        stereo_bonds: vec![(
            BondId(0),
            [2, 3, 4, 5]
                .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                .into(),
            EntitySpan::Unchanged(StereoBondForm::new(StereoKind::CisTrans, 0u32)),
        )],
        ..Default::default()
    }),
    vec![
        NodeId(1), NodeId(0), NodeId(4), NodeId(5), NodeId(2), NodeId(3),
        NodeId(6), NodeId(9), NodeId(10), NodeId(7), NodeId(8), NodeId(11),
    ],
    true,
)]
#[case::lhs_configuration(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::H, Element::H, Element::F, Element::Cl]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        bonds: (1..=4)
            .map(|ligand| {
                (AtomId(0), AtomId(ligand), EntitySpan::Unchanged(BondForm::from_order(1)))
            })
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            (1..=4)
                .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                .collect(),
            EntitySpan::Modified {
                lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
                rhs: StereoAtomForm::default(),
            },
        )],
        ..Default::default()
    }),
    vec![
        NodeId(0), NodeId(2), NodeId(1), NodeId(3), NodeId(4),
        NodeId(6), NodeId(5), NodeId(7), NodeId(8), NodeId(9),
    ],
    false,
)]
#[case::rhs_configuration(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::H, Element::H, Element::F, Element::Cl]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        bonds: (1..=4)
            .map(|ligand| {
                (AtomId(0), AtomId(ligand), EntitySpan::Unchanged(BondForm::from_order(1)))
            })
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            (1..=4)
                .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                .collect(),
            EntitySpan::Modified {
                lhs: StereoAtomForm::default(),
                rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            },
        )],
        ..Default::default()
    }),
    vec![
        NodeId(0), NodeId(2), NodeId(1), NodeId(3), NodeId(4),
        NodeId(6), NodeId(5), NodeId(7), NodeId(8), NodeId(9),
    ],
    false,
)]
#[case::modified_unchanged_exchange(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [
            Element::C, Element::F, Element::Cl, Element::Br, Element::I,
            Element::C, Element::F, Element::Cl, Element::Br, Element::I,
        ]
        .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
        .into(),
        bonds: [(0, 1), (0, 2), (0, 3), (0, 4), (5, 6), (5, 7), (5, 8), (5, 9)]
            .map(|(site, ligand)| {
                (AtomId(site), AtomId(ligand), EntitySpan::Unchanged(BondForm::from_order(1)))
            })
            .into(),
        stereo_atoms: vec![
            (
                AtomId(0),
                [1, 2, 3, 4].map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).into(),
                EntitySpan::Unchanged(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32)),
            ),
            (
                AtomId(5),
                [6, 7, 8, 9].map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).into(),
                EntitySpan::Modified {
                    lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
                    rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
                },
            ),
        ],
        ..Default::default()
    }),
    vec![
        NodeId(5), NodeId(6), NodeId(7), NodeId(8), NodeId(9),
        NodeId(0), NodeId(1), NodeId(2), NodeId(3), NodeId(4),
        NodeId(14), NodeId(15), NodeId(16), NodeId(17),
        NodeId(10), NodeId(11), NodeId(12), NodeId(13),
        NodeId(19), NodeId(18),
    ],
    false,
)]
#[case::modified_stereo_bond_endpoint_exchange(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::C, Element::F, Element::Cl, Element::F, Element::Cl]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Unchanged(BondForm::from_order(2))),
            (AtomId(0), AtomId(2), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(0), AtomId(3), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(1), AtomId(4), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(1), AtomId(5), EntitySpan::Unchanged(BondForm::from_order(1))),
        ],
        stereo_bonds: vec![(
            BondId(0),
            [2, 3, 4, 5]
                .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                .into(),
            EntitySpan::Modified {
                lhs: StereoBondForm::new(StereoKind::CisTrans, 0u32),
                rhs: StereoBondForm::new(StereoKind::CisTrans, 1u32),
            },
        )],
        ..Default::default()
    }),
    vec![
        NodeId(1), NodeId(0), NodeId(4), NodeId(5), NodeId(2), NodeId(3),
        NodeId(6), NodeId(9), NodeId(10), NodeId(7), NodeId(8), NodeId(11),
    ],
    true,
)]
#[case::rhs_topicity_one_site(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [
            Element::C, Element::F, Element::Cl, Element::Br, Element::I,
            Element::C, Element::F, Element::Cl, Element::Br, Element::I,
        ]
        .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
        .into(),
        bonds: [(0, 1), (0, 2), (0, 3), (0, 4), (5, 6), (5, 7), (5, 8), (5, 9)]
            .map(|(site, ligand)| {
                (AtomId(site), AtomId(ligand), EntitySpan::Unchanged(BondForm::from_order(1)))
            })
            .into(),
        stereo_atoms: vec![
            (
                AtomId(0),
                [1, 2, 3, 4].map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).into(),
                EntitySpan::Unchanged(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32)),
            ),
            (
                AtomId(5),
                [6, 7, 8, 9].map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).into(),
                EntitySpan::Modified {
                    lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
                    rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32).with_constraint(
                        StereoAtomConstraintForm::Topicity(TopicityForm {
                            pair: StereoLigandPair::new(0usize.into(), 1usize.into()),
                            relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
                        }),
                    ),
                },
            ),
        ],
        ..Default::default()
    }),
    vec![
        NodeId(5), NodeId(6), NodeId(7), NodeId(8), NodeId(9),
        NodeId(0), NodeId(1), NodeId(2), NodeId(3), NodeId(4),
        NodeId(14), NodeId(15), NodeId(16), NodeId(17),
        NodeId(10), NodeId(11), NodeId(12), NodeId(13),
        NodeId(19), NodeId(18),
    ],
    false,
)]
#[case::rhs_topicity_both_sites(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [
            Element::C, Element::F, Element::Cl, Element::Br, Element::I,
            Element::C, Element::F, Element::Cl, Element::Br, Element::I,
        ]
        .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
        .into(),
        bonds: [(0, 1), (0, 2), (0, 3), (0, 4), (5, 6), (5, 7), (5, 8), (5, 9)]
            .map(|(site, ligand)| {
                (AtomId(site), AtomId(ligand), EntitySpan::Unchanged(BondForm::from_order(1)))
            })
            .into(),
        stereo_atoms: vec![
            (
                AtomId(0),
                [1, 2, 3, 4].map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).into(),
                EntitySpan::Modified {
                    lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
                    rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32).with_constraint(
                        StereoAtomConstraintForm::Topicity(TopicityForm {
                            pair: StereoLigandPair::new(0usize.into(), 1usize.into()),
                            relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
                        }),
                    ),
                },
            ),
            (
                AtomId(5),
                [6, 7, 8, 9].map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).into(),
                EntitySpan::Modified {
                    lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
                    rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32).with_constraint(
                        StereoAtomConstraintForm::Topicity(TopicityForm {
                            pair: StereoLigandPair::new(0usize.into(), 1usize.into()),
                            relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
                        }),
                    ),
                },
            ),
        ],
        ..Default::default()
    }),
    vec![
        NodeId(5), NodeId(6), NodeId(7), NodeId(8), NodeId(9),
        NodeId(0), NodeId(1), NodeId(2), NodeId(3), NodeId(4),
        NodeId(14), NodeId(15), NodeId(16), NodeId(17),
        NodeId(10), NodeId(11), NodeId(12), NodeId(13),
        NodeId(19), NodeId(18),
    ],
    true,
)]
fn test_retain_stereo_preserving_generators(
    #[case] span: ReactionSpan,
    #[case] generator: Vec<NodeId>,
    #[case] expected_retained: bool,
) {
    let incidence_graph = span.incidence_graph(IncidenceLevel::Full);
    let mut generators = vec![generator.clone()];
    let expected = if expected_retained {
        vec![generator]
    } else {
        Vec::new()
    };

    retain_stereo_preserving_generators(&span, &incidence_graph, &mut generators);

    assert_eq!(generators, expected);
}

#[rstest]
#[case::added_over_one_atom(
    vec![ConstraintSpan::Added(Constraint::Molecule(MoleculeConstraint::ChargeSum {
        atoms: Some(vec![AtomId(0)]),
        sum: NumForm::Lit(0),
    }))],
    Vec::new()
)]
#[case::added_over_both_atoms(
    vec![ConstraintSpan::Added(Constraint::Molecule(MoleculeConstraint::ChargeSum {
        atoms: Some(vec![AtomId(0), AtomId(1)]),
        sum: NumForm::Lit(0),
    }))],
    vec![[1, 0].map(NodeId).into()]
)]
#[case::removed_and_added_over_one_atom_each(
    vec![
        ConstraintSpan::Removed(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(0)]),
            sum: NumForm::Lit(0),
        })),
        ConstraintSpan::Added(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(1)]),
            sum: NumForm::Lit(0),
        })),
    ],
    Vec::new()
)]
fn test_retain_preserving_generators(
    #[case] constraints: Vec<ConstraintSpan>,
    #[case] expected: Vec<Vec<NodeId>>,
) {
    let span = ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 2],
        constraints,
        ..Default::default()
    });
    let incidence_graph = span.incidence_graph(IncidenceLevel::Full);
    let mut generators: Vec<Vec<NodeId>> = vec![[1, 0].map(NodeId).into()];

    retain_preserving_generators(&span, &incidence_graph, &mut generators);

    assert_eq!(generators, expected);
}

#[fixture]
fn initial_color_span() -> ReactionSpan {
    let carbon = AtomForm::from_element(Element::C);
    let nitrogen = AtomForm::from_element(Element::N);
    let ligand = StereoLigand::new;
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(carbon.clone()),
            EntitySpan::Unchanged(carbon.clone()),
            EntitySpan::Added(carbon.clone()),
            EntitySpan::Removed(carbon.clone()),
            EntitySpan::Modified {
                lhs: carbon.clone(),
                rhs: nitrogen.clone(),
            },
            EntitySpan::Modified {
                lhs: carbon.clone(),
                rhs: nitrogen.clone(),
            },
            EntitySpan::Modified {
                lhs: nitrogen,
                rhs: carbon.clone(),
            },
            EntitySpan::Unchanged(carbon.clone()),
            EntitySpan::Unchanged(carbon.clone()),
            EntitySpan::Unchanged(carbon.clone()),
            EntitySpan::Unchanged(carbon.clone()),
            EntitySpan::Unchanged(carbon.clone()),
            EntitySpan::Unchanged(carbon),
        ],
        bonds: vec![
            (
                AtomId(0),
                AtomId(1),
                EntitySpan::Unchanged(BondForm::from_order(1)),
            ),
            (
                AtomId(7),
                AtomId(8),
                EntitySpan::Unchanged(BondForm::from_order(1)),
            ),
            (
                AtomId(1),
                AtomId(2),
                EntitySpan::Added(BondForm::from_order(1)),
            ),
            (
                AtomId(1),
                AtomId(3),
                EntitySpan::Removed(BondForm::from_order(1)),
            ),
            (
                AtomId(4),
                AtomId(5),
                EntitySpan::Modified {
                    lhs: BondForm::from_order(1),
                    rhs: BondForm::from_order(2),
                },
            ),
            (
                AtomId(5),
                AtomId(6),
                EntitySpan::Modified {
                    lhs: BondForm::from_order(1),
                    rhs: BondForm::from_order(2),
                },
            ),
            (
                AtomId(9),
                AtomId(10),
                EntitySpan::Modified {
                    lhs: BondForm::from_order(2),
                    rhs: BondForm::from_order(1),
                },
            ),
            (
                AtomId(7),
                AtomId(11),
                EntitySpan::Unchanged(BondForm::from_order(1)),
            ),
            (
                AtomId(8),
                AtomId(12),
                EntitySpan::Unchanged(BondForm::from_order(1)),
            ),
            (
                AtomId(0),
                AtomId(9),
                EntitySpan::Unchanged(BondForm::from_order(1)),
            ),
        ],
        dative: vec![
            (
                vec![AtomId(0)],
                AtomId(7),
                EntitySpan::Unchanged(DativeBondForm::from_order(1)),
            ),
            (
                vec![AtomId(1)],
                AtomId(8),
                EntitySpan::Unchanged(DativeBondForm::from_order(1)),
            ),
            (
                vec![AtomId(9)],
                AtomId(10),
                EntitySpan::Added(DativeBondForm::from_order(1)),
            ),
        ],
        aromatic: vec![
            (
                vec![AtomId(0), AtomId(1)],
                EntitySpan::Unchanged(AromaticSystemForm::from_electrons(vec![1, 2])),
            ),
            (
                vec![AtomId(7), AtomId(8)],
                EntitySpan::Modified {
                    lhs: AromaticSystemForm::from_electrons(vec![1, 2]),
                    rhs: AromaticSystemForm::from_electrons(vec![2, 1]),
                },
            ),
            (
                vec![AtomId(9), AtomId(10)],
                EntitySpan::Modified {
                    lhs: AromaticSystemForm::from_electrons(vec![1, 2]),
                    rhs: AromaticSystemForm::from_electrons(vec![1, 2]).with_charge(1_i64),
                },
            ),
            (
                vec![AtomId(11), AtomId(12)],
                EntitySpan::Added(AromaticSystemForm::from_electrons(vec![1, 2])),
            ),
        ],
        multicenter: vec![
            (
                vec![AtomId(0), AtomId(7)],
                EntitySpan::Unchanged(MulticenterBondForm::from_electrons(vec![1, 2])),
            ),
            (
                vec![AtomId(1), AtomId(8)],
                EntitySpan::Modified {
                    lhs: MulticenterBondForm::from_electrons(vec![1, 2]),
                    rhs: MulticenterBondForm::from_electrons(vec![2, 1]),
                },
            ),
            (
                vec![AtomId(9), AtomId(11)],
                EntitySpan::Removed(MulticenterBondForm::from_electrons(vec![1, 2])),
            ),
        ],
        noncovalent: vec![
            (
                [AtomId(0), AtomId(8)],
                EntitySpan::Unchanged(NoncovalentBondForm::from_kind(
                    NoncovalentBondKind::HydrogenBond,
                )),
            ),
            (
                [AtomId(1), AtomId(7)],
                EntitySpan::Unchanged(NoncovalentBondForm::from_kind(
                    NoncovalentBondKind::HydrogenBond,
                )),
            ),
            (
                [AtomId(9), AtomId(12)],
                EntitySpan::Modified {
                    lhs: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                    rhs: NoncovalentBondForm::from_kind(NoncovalentBondKind::Ionic),
                },
            ),
        ],
        stereo_atoms: vec![
            (
                AtomId(7),
                vec![
                    ligand(AtomId(8), StereoLigandKind::Atom),
                    ligand(AtomId(11), StereoLigandKind::Atom),
                    ligand(AtomId(7), StereoLigandKind::ImplicitHydrogen),
                    ligand(AtomId(7), StereoLigandKind::LonePair),
                ],
                EntitySpan::Unchanged(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32)),
            ),
            (
                AtomId(8),
                vec![
                    ligand(AtomId(7), StereoLigandKind::Atom),
                    ligand(AtomId(12), StereoLigandKind::Atom),
                    ligand(AtomId(8), StereoLigandKind::ImplicitHydrogen),
                    ligand(AtomId(8), StereoLigandKind::LonePair),
                ],
                EntitySpan::Modified {
                    lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
                    rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
                },
            ),
            (
                AtomId(0),
                vec![
                    ligand(AtomId(1), StereoLigandKind::Atom),
                    ligand(AtomId(9), StereoLigandKind::Atom),
                    ligand(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    ligand(AtomId(0), StereoLigandKind::LonePair),
                ],
                EntitySpan::Added(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32)),
            ),
        ],
        stereo_bonds: vec![
            (
                BondId(1),
                vec![
                    ligand(AtomId(11), StereoLigandKind::Atom),
                    ligand(AtomId(7), StereoLigandKind::ImplicitHydrogen),
                    ligand(AtomId(12), StereoLigandKind::Atom),
                    ligand(AtomId(8), StereoLigandKind::ImplicitHydrogen),
                ],
                EntitySpan::Unchanged(StereoBondForm::new(StereoKind::CisTrans, 0u32)),
            ),
            (
                BondId(0),
                vec![
                    ligand(AtomId(9), StereoLigandKind::Atom),
                    ligand(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    ligand(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    ligand(AtomId(1), StereoLigandKind::LonePair),
                ],
                EntitySpan::Modified {
                    lhs: StereoBondForm::new(StereoKind::CisTrans, 0u32),
                    rhs: StereoBondForm::new(StereoKind::CisTrans, 1u32),
                },
            ),
        ],
        ..Default::default()
    })
}

#[rstest]
#[case::atom_unchanged(Entity::Atom(AtomId(0)), Entity::Atom(AtomId(1)), true)]
#[case::atom_added(Entity::Atom(AtomId(0)), Entity::Atom(AtomId(2)), false)]
#[case::atom_removed(Entity::Atom(AtomId(0)), Entity::Atom(AtomId(3)), false)]
#[case::atom_added_removed(Entity::Atom(AtomId(2)), Entity::Atom(AtomId(3)), false)]
#[case::atom_modified(Entity::Atom(AtomId(4)), Entity::Atom(AtomId(5)), true)]
#[case::atom_modified_reversed(Entity::Atom(AtomId(4)), Entity::Atom(AtomId(6)), false)]
#[case::bond_unchanged(Entity::Bond(BondId(0)), Entity::Bond(BondId(1)), true)]
#[case::bond_added(Entity::Bond(BondId(0)), Entity::Bond(BondId(2)), false)]
#[case::bond_removed(Entity::Bond(BondId(0)), Entity::Bond(BondId(3)), false)]
#[case::bond_modified(Entity::Bond(BondId(4)), Entity::Bond(BondId(5)), true)]
#[case::bond_modified_reversed(Entity::Bond(BondId(4)), Entity::Bond(BondId(6)), false)]
#[case::dative_unchanged(
    Entity::DativeBond(DativeBondId(0)),
    Entity::DativeBond(DativeBondId(1)),
    true
)]
#[case::dative_added(
    Entity::DativeBond(DativeBondId(0)),
    Entity::DativeBond(DativeBondId(2)),
    false
)]
#[case::aromatic_modified_electrons(
    Entity::AromaticSystem(AromaticSystemId(0)),
    Entity::AromaticSystem(AromaticSystemId(1)),
    true
)]
#[case::aromatic_modified_charge(
    Entity::AromaticSystem(AromaticSystemId(0)),
    Entity::AromaticSystem(AromaticSystemId(2)),
    false
)]
#[case::aromatic_added(
    Entity::AromaticSystem(AromaticSystemId(0)),
    Entity::AromaticSystem(AromaticSystemId(3)),
    false
)]
#[case::multicenter_modified_electrons(
    Entity::MulticenterBond(MulticenterBondId(0)),
    Entity::MulticenterBond(MulticenterBondId(1)),
    true
)]
#[case::multicenter_removed(
    Entity::MulticenterBond(MulticenterBondId(0)),
    Entity::MulticenterBond(MulticenterBondId(2)),
    false
)]
#[case::noncovalent_unchanged(
    Entity::NoncovalentBond(NoncovalentBondId(0)),
    Entity::NoncovalentBond(NoncovalentBondId(1)),
    true
)]
#[case::noncovalent_modified(
    Entity::NoncovalentBond(NoncovalentBondId(0)),
    Entity::NoncovalentBond(NoncovalentBondId(2)),
    false
)]
#[case::stereo_atom_modified_configuration(
    Entity::StereoAtom(StereoAtomId(0)),
    Entity::StereoAtom(StereoAtomId(1)),
    true
)]
#[case::stereo_atom_added(
    Entity::StereoAtom(StereoAtomId(0)),
    Entity::StereoAtom(StereoAtomId(2)),
    false
)]
#[case::stereo_bond_modified_configuration(
    Entity::StereoBond(StereoBondId(0)),
    Entity::StereoBond(StereoBondId(1)),
    true
)]
#[case::entity_kind(Entity::Atom(AtomId(0)), Entity::Bond(BondId(0)), false)]
fn test_initial_colors(
    initial_color_span: ReactionSpan,
    #[case] lhs: Entity,
    #[case] rhs: Entity,
    #[case] expected: bool,
) {
    let span = initial_color_span
        .normalize()
        .expect("fixed span normalizes");
    let incidence_graph = span.incidence_graph(IncidenceLevel::Full);
    let colors = initial_colors(&span, &incidence_graph, DescriptionLevel::Structure)
        .expect("fixed span has colors");

    assert_eq!(
        colors.entities[incidence_graph.node_of(lhs).index()]
            == colors.entities[incidence_graph.node_of(rhs).index()],
        expected
    );
}

#[rstest]
#[case::aromatic_unchanged(
    [Entity::AromaticSystem(AromaticSystemId(0)), Entity::Atom(AtomId(0))],
    [Entity::AromaticSystem(AromaticSystemId(2)), Entity::Atom(AtomId(9))],
    true
)]
#[case::aromatic_position(
    [Entity::AromaticSystem(AromaticSystemId(0)), Entity::Atom(AtomId(0))],
    [Entity::AromaticSystem(AromaticSystemId(0)), Entity::Atom(AtomId(1))],
    false
)]
#[case::aromatic_modified(
    [Entity::AromaticSystem(AromaticSystemId(0)), Entity::Atom(AtomId(0))],
    [Entity::AromaticSystem(AromaticSystemId(1)), Entity::Atom(AtomId(7))],
    false
)]
#[case::aromatic_modified_reversed(
    [Entity::AromaticSystem(AromaticSystemId(1)), Entity::Atom(AtomId(7))],
    [Entity::AromaticSystem(AromaticSystemId(1)), Entity::Atom(AtomId(8))],
    false
)]
#[case::aromatic_added(
    [Entity::AromaticSystem(AromaticSystemId(0)), Entity::Atom(AtomId(0))],
    [Entity::AromaticSystem(AromaticSystemId(3)), Entity::Atom(AtomId(11))],
    false
)]
#[case::multicenter_removed(
    [Entity::MulticenterBond(MulticenterBondId(0)), Entity::Atom(AtomId(0))],
    [Entity::MulticenterBond(MulticenterBondId(2)), Entity::Atom(AtomId(9))],
    false
)]
#[case::role(
    [Entity::AromaticSystem(AromaticSystemId(1)), Entity::Atom(AtomId(7))],
    [Entity::MulticenterBond(MulticenterBondId(1)), Entity::Atom(AtomId(1))],
    false
)]
#[case::bond_endpoint(
    [Entity::Bond(BondId(0)), Entity::Atom(AtomId(0))],
    [Entity::Bond(BondId(2)), Entity::Atom(AtomId(1))],
    true
)]
#[case::dative_role(
    [Entity::DativeBond(DativeBondId(0)), Entity::Atom(AtomId(0))],
    [Entity::DativeBond(DativeBondId(0)), Entity::Atom(AtomId(7))],
    false
)]
fn test_initial_colors_incidence(
    initial_color_span: ReactionSpan,
    #[case] lhs: [Entity; 2],
    #[case] rhs: [Entity; 2],
    #[case] expected: bool,
) {
    let span = initial_color_span
        .normalize()
        .expect("fixed span normalizes");
    let incidence_graph = span.incidence_graph(IncidenceLevel::Full);
    let colors = initial_colors(&span, &incidence_graph, DescriptionLevel::Structure)
        .expect("fixed span has colors");
    let color = |[owner, atom]: [Entity; 2]| {
        let nodes = [owner, atom].map(|entity| incidence_graph.node_of(entity));
        let (edge, _) = incidence_graph
            .incidences()
            .find(|&(edge, _)| {
                let [first, second] = incidence_graph.graph().edge_endpoints(edge);
                nodes == [first, second] || nodes == [second, first]
            })
            .expect("the case names an incidence of the fixed span");
        colors.incidences[edge.index()]
    };

    assert_eq!(color(lhs) == color(rhs), expected);
}

#[rstest]
#[case::aromatic_electron_count_length(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 3],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            EntitySpan::Modified {
                lhs: AromaticSystemForm::from_electrons(vec![1, 1, 1]),
                rhs: AromaticSystemForm::from_electrons(vec![2, 1]),
            },
        )],
        ..Default::default()
    }),
    DescriptionLevel::Constitution
)]
#[case::multicenter_electron_count_length(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::B, Element::H, Element::B]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        multicenter: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            EntitySpan::Added(MulticenterBondForm::from_electrons(vec![1, 0, 1, 0])),
        )],
        ..Default::default()
    }),
    DescriptionLevel::Constitution
)]
#[case::stereo_atom_coset_range(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        bonds: (1..=4)
            .map(|ligand| {
                (
                    AtomId(0),
                    AtomId(ligand),
                    EntitySpan::Unchanged(BondForm::from_order(1)),
                )
            })
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            (1..=4)
                .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                .collect(),
            EntitySpan::Modified {
                lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
                rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 2u32),
            },
        )],
        ..Default::default()
    }),
    DescriptionLevel::Structure
)]
fn test_initial_colors_error(#[case] span: ReactionSpan, #[case] level: DescriptionLevel) {
    let incidence_graph = span.incidence_graph(match level {
        DescriptionLevel::Topology => IncidenceLevel::Topology,
        DescriptionLevel::Constitution => IncidenceLevel::Constitution,
        DescriptionLevel::Structure | DescriptionLevel::Full => IncidenceLevel::Full,
    });

    assert_eq!(
        initial_colors(&span, &incidence_graph, level),
        Err(Contradiction)
    );
}

#[rstest]
#[case::unchanged_bond(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 2],
        bonds: vec![(AtomId(0), AtomId(1), EntitySpan::Unchanged(BondForm::from_order(1)))],
        ..Default::default()
    }),
    Graph::new(2, &[[0, 1]]),
    vec![AutomorphismColor::Entity(0); 2],
    vec![
        SubdivisionNodeSource::Node(NodeId(0)),
        SubdivisionNodeSource::Node(NodeId(1)),
    ],
    vec![vec![NodeId(0), NodeId(1)]],
)]
#[case::added_atom(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Added(AtomForm::from_element(Element::C)),
        ],
        bonds: vec![(AtomId(0), AtomId(1), EntitySpan::Added(BondForm::from_order(1)))],
        ..Default::default()
    }),
    Graph::new(2, &[[0, 1]]),
    vec![AutomorphismColor::Entity(0), AutomorphismColor::Entity(1)],
    vec![
        SubdivisionNodeSource::Node(NodeId(0)),
        SubdivisionNodeSource::Node(NodeId(1)),
    ],
    vec![vec![NodeId(0)], vec![NodeId(1)]],
)]
#[case::modified_bond(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 3],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(1), AtomId(2), EntitySpan::Unchanged(BondForm::from_order(1))),
            (
                AtomId(0),
                AtomId(2),
                EntitySpan::Modified {
                    lhs: BondForm::from_order(1),
                    rhs: BondForm::from_order(2),
                },
            ),
        ],
        ..Default::default()
    }),
    Graph::new(4, &[[0, 1], [1, 2], [0, 3], [3, 2]]),
    vec![
        AutomorphismColor::Entity(0),
        AutomorphismColor::Entity(0),
        AutomorphismColor::Entity(0),
        AutomorphismColor::Entity(2),
    ],
    vec![
        SubdivisionNodeSource::Node(NodeId(0)),
        SubdivisionNodeSource::Node(NodeId(1)),
        SubdivisionNodeSource::Node(NodeId(2)),
        SubdivisionNodeSource::Node(NodeId(3)),
    ],
    vec![vec![NodeId(0), NodeId(1), NodeId(2)], vec![NodeId(3)]],
)]
#[case::removed_and_added_bond(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 2],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Added(BondForm::from_order(2))),
            (AtomId(0), AtomId(1), EntitySpan::Removed(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    Graph::new(3, &[[0, 2], [2, 1], [0, 1]]),
    vec![
        AutomorphismColor::Entity(0),
        AutomorphismColor::Entity(0),
        AutomorphismColor::Entity(2),
    ],
    vec![
        SubdivisionNodeSource::Node(NodeId(0)),
        SubdivisionNodeSource::Node(NodeId(1)),
        SubdivisionNodeSource::Node(NodeId(2)),
    ],
    vec![vec![NodeId(0), NodeId(1)], vec![NodeId(2)]],
)]
fn test_topology_search_input(
    #[case] span: ReactionSpan,
    #[case] expected_graph: Graph,
    #[case] expected_colors: Vec<AutomorphismColor>,
    #[case] expected_node_sources: Vec<SubdivisionNodeSource>,
    #[case] expected_partition: Vec<Vec<NodeId>>,
) {
    let span = span.normalize().unwrap();
    let incidence_graph = span.incidence_graph(IncidenceLevel::Topology);
    let colors = initial_colors(&span, &incidence_graph, DescriptionLevel::Topology).unwrap();
    let input = topology_search_input(&span, &incidence_graph, &colors.entities);

    assert_eq!(input.adapter.graph, expected_graph);
    assert_eq!(input.adapter.colors, expected_colors);
    assert_eq!(input.adapter.node_sources, expected_node_sources);
    assert_eq!(input.adapter.source_node_count, span.atoms().len());
    assert_eq!(input.initial_partition.cells, expected_partition);
}

#[rstest]
#[case::bond_participants(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 3],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(1), AtomId(2), EntitySpan::Unchanged(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    vec![
        Entity::Atom(AtomId(2)),
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::Bond(BondId(1)),
        Entity::Bond(BondId(0)),
    ],
)]
#[case::added_bond(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 3],
        bonds: vec![
            (AtomId(1), AtomId(2), EntitySpan::Added(BondForm::from_order(1))),
            (AtomId(0), AtomId(1), EntitySpan::Unchanged(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    vec![
        Entity::Atom(AtomId(2)),
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::Bond(BondId(1)),
        Entity::Bond(BondId(0)),
    ],
)]
#[case::removed_and_added_bond(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 2],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Added(BondForm::from_order(2))),
            (AtomId(0), AtomId(1), EntitySpan::Removed(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    vec![
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::Bond(BondId(1)),
        Entity::Bond(BondId(0)),
    ],
)]
#[case::added_dative_bond(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 4],
        dative: vec![
            (vec![AtomId(1)], AtomId(3), EntitySpan::Added(DativeBondForm::from_order(1))),
            (vec![AtomId(0)], AtomId(2), EntitySpan::Unchanged(DativeBondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    vec![
        Entity::Atom(AtomId(3)),
        Entity::Atom(AtomId(2)),
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::DativeBond(DativeBondId(1)),
        Entity::DativeBond(DativeBondId(0)),
    ],
)]
#[case::added_aromatic_system(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 4],
        aromatic: vec![
            (vec![AtomId(2), AtomId(3)], EntitySpan::Added(AromaticSystemForm::default())),
            (vec![AtomId(0), AtomId(1)], EntitySpan::Unchanged(AromaticSystemForm::default())),
        ],
        ..Default::default()
    }),
    vec![
        Entity::Atom(AtomId(3)),
        Entity::Atom(AtomId(2)),
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::AromaticSystem(AromaticSystemId(1)),
        Entity::AromaticSystem(AromaticSystemId(0)),
    ],
)]
#[case::added_multicenter_bond(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 4],
        multicenter: vec![
            (vec![AtomId(2), AtomId(3)], EntitySpan::Added(MulticenterBondForm::default())),
            (vec![AtomId(0), AtomId(1)], EntitySpan::Unchanged(MulticenterBondForm::default())),
        ],
        ..Default::default()
    }),
    vec![
        Entity::Atom(AtomId(3)),
        Entity::Atom(AtomId(2)),
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::MulticenterBond(MulticenterBondId(1)),
        Entity::MulticenterBond(MulticenterBondId(0)),
    ],
)]
#[case::added_noncovalent_bond(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 4],
        noncovalent: vec![
            ([AtomId(2), AtomId(3)], EntitySpan::Added(NoncovalentBondForm::default())),
            ([AtomId(0), AtomId(1)], EntitySpan::Unchanged(NoncovalentBondForm::default())),
        ],
        ..Default::default()
    }),
    vec![
        Entity::Atom(AtomId(3)),
        Entity::Atom(AtomId(2)),
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::NoncovalentBond(NoncovalentBondId(1)),
        Entity::NoncovalentBond(NoncovalentBondId(0)),
    ],
)]
#[case::added_stereo_atom(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 2],
        stereo_atoms: vec![
            (AtomId(1), Vec::new(), EntitySpan::Added(StereoAtomForm::default())),
            (AtomId(0), Vec::new(), EntitySpan::Unchanged(StereoAtomForm::default())),
        ],
        ..Default::default()
    }),
    vec![
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::StereoAtom(StereoAtomId(1)),
        Entity::StereoAtom(StereoAtomId(0)),
    ],
)]
#[case::added_stereo_bond(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 4],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Unchanged(BondForm::from_order(2))),
            (AtomId(2), AtomId(3), EntitySpan::Unchanged(BondForm::from_order(2))),
        ],
        stereo_bonds: vec![
            (
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(3), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(3), StereoLigandKind::LonePair),
                ],
                EntitySpan::Added(StereoBondForm::default()),
            ),
            (
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                ],
                EntitySpan::Unchanged(StereoBondForm::default()),
            ),
        ],
        ..Default::default()
    }),
    vec![
        Entity::Atom(AtomId(3)),
        Entity::Atom(AtomId(2)),
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::Bond(BondId(1)),
        Entity::Bond(BondId(0)),
        Entity::StereoBond(StereoBondId(1)),
        Entity::StereoBond(StereoBondId(0)),
    ],
)]
fn test_entity_order(#[case] span: ReactionSpan, #[case] expected: Vec<Entity>) {
    let incidence_graph = span.incidence_graph(IncidenceLevel::Full);
    let mut order = incidence_graph.graph().node_ids().collect::<Vec<_>>();
    order.reverse();

    assert_eq!(
        entity_order(&span, &incidence_graph, DescriptionLevel::Full, &order)
            .into_iter()
            .map(|node| incidence_graph.entity(node))
            .collect::<Vec<_>>(),
        expected,
    );
}

#[rstest]
#[case::unchanged_atoms(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 3],
        ..Default::default()
    }),
    DescriptionLevel::Topology,
    [6, 1],
)]
#[case::added_atoms(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Added(AtomForm::from_element(Element::C)); 3],
        ..Default::default()
    }),
    DescriptionLevel::Topology,
    [6, 1],
)]
#[case::removed_atoms(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Removed(AtomForm::from_element(Element::C)); 3],
        ..Default::default()
    }),
    DescriptionLevel::Topology,
    [6, 1],
)]
#[case::modified_atom(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Modified {
                lhs: AtomForm::from_element(Element::C),
                rhs: AtomForm::from_element(Element::N),
            },
        ],
        ..Default::default()
    }),
    DescriptionLevel::Topology,
    [2, 1],
)]
#[case::removed_and_added_bond(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 2],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Added(BondForm::from_order(2))),
            (AtomId(0), AtomId(1), EntitySpan::Removed(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    DescriptionLevel::Topology,
    [2, 1],
)]
#[case::unchanged_dative_bond(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::N, Element::N, Element::B]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        dative: vec![(
            vec![AtomId(0), AtomId(1)],
            AtomId(2),
            EntitySpan::Unchanged(DativeBondForm::from_order(1)),
        )],
        ..Default::default()
    }),
    DescriptionLevel::Constitution,
    [2, 1],
)]
#[case::modified_dative_bond(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::N, Element::N, Element::B]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        dative: vec![(
            vec![AtomId(0), AtomId(1)],
            AtomId(2),
            EntitySpan::Modified {
                lhs: DativeBondForm::from_order(1),
                rhs: DativeBondForm::from_order(2),
            },
        )],
        ..Default::default()
    }),
    DescriptionLevel::Constitution,
    [2, 1],
)]
#[case::unchanged_aromatic_system(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 3],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            EntitySpan::Unchanged(AromaticSystemForm::default()),
        )],
        ..Default::default()
    }),
    DescriptionLevel::Constitution,
    [6, 1],
)]
#[case::modified_aromatic_system(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 3],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            EntitySpan::Modified {
                lhs: AromaticSystemForm::from_electrons(vec![1, 1, 1]),
                rhs: AromaticSystemForm::from_electrons(vec![1, 1, 2]),
            },
        )],
        ..Default::default()
    }),
    DescriptionLevel::Constitution,
    [2, 1],
)]
#[case::aromatic_electron_counts(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 5],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(1), AtomId(2), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(2), AtomId(3), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(3), AtomId(4), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(4), AtomId(0), EntitySpan::Unchanged(BondForm::from_order(1))),
        ],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4)],
            EntitySpan::Unchanged(AromaticSystemForm::from_electrons(vec![2, 1, 1, 1, 1])),
        )],
        ..Default::default()
    }),
    DescriptionLevel::Constitution,
    [2, 1],
)]
#[case::equal_stereocentres(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .repeat(3)
            .into_iter()
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .collect(),
        bonds: [0, 5, 10]
            .into_iter()
            .flat_map(|site| {
                (1..=4).map(move |ligand| {
                    (
                        AtomId(site),
                        AtomId(site + ligand),
                        EntitySpan::Unchanged(BondForm::from_order(1)),
                    )
                })
            })
            .collect(),
        stereo_atoms: vec![
            (
                AtomId(0),
                [1, 2, 3, 4].map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).into(),
                EntitySpan::Unchanged(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32)),
            ),
            (
                AtomId(5),
                [6, 7, 8, 9].map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).into(),
                EntitySpan::Unchanged(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32)),
            ),
            (
                AtomId(10),
                [11, 12, 13, 14]
                    .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                    .into(),
                EntitySpan::Unchanged(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32)),
            ),
        ],
        ..Default::default()
    }),
    DescriptionLevel::Structure,
    [6, 1],
)]
#[case::lhs_configuration(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .repeat(3)
            .into_iter()
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .collect(),
        bonds: [0, 5, 10]
            .into_iter()
            .flat_map(|site| {
                (1..=4).map(move |ligand| {
                    (
                        AtomId(site),
                        AtomId(site + ligand),
                        EntitySpan::Unchanged(BondForm::from_order(1)),
                    )
                })
            })
            .collect(),
        stereo_atoms: vec![
            (
                AtomId(0),
                [1, 2, 3, 4].map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).into(),
                EntitySpan::Unchanged(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32)),
            ),
            (
                AtomId(5),
                [6, 7, 8, 9].map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).into(),
                EntitySpan::Unchanged(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32)),
            ),
            (
                AtomId(10),
                [11, 12, 13, 14]
                    .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                    .into(),
                EntitySpan::Modified {
                    lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
                    rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
                },
            ),
        ],
        ..Default::default()
    }),
    DescriptionLevel::Structure,
    [6, 3],
)]
#[case::rhs_configuration(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .repeat(3)
            .into_iter()
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .collect(),
        bonds: [0, 5, 10]
            .into_iter()
            .flat_map(|site| {
                (1..=4).map(move |ligand| {
                    (
                        AtomId(site),
                        AtomId(site + ligand),
                        EntitySpan::Unchanged(BondForm::from_order(1)),
                    )
                })
            })
            .collect(),
        stereo_atoms: vec![
            (
                AtomId(0),
                [1, 2, 3, 4].map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).into(),
                EntitySpan::Unchanged(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32)),
            ),
            (
                AtomId(5),
                [6, 7, 8, 9].map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).into(),
                EntitySpan::Unchanged(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32)),
            ),
            (
                AtomId(10),
                [11, 12, 13, 14]
                    .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                    .into(),
                EntitySpan::Modified {
                    lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
                    rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
                },
            ),
        ],
        ..Default::default()
    }),
    DescriptionLevel::Structure,
    [6, 3],
)]
#[case::prochiral_rhs_configuration(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::H, Element::H, Element::F, Element::Cl]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        bonds: (1..=4)
            .map(|ligand| {
                (AtomId(0), AtomId(ligand), EntitySpan::Unchanged(BondForm::from_order(1)))
            })
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            (1..=4)
                .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                .collect(),
            EntitySpan::Modified {
                lhs: StereoAtomForm::default(),
                rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            },
        )],
        ..Default::default()
    }),
    DescriptionLevel::Structure,
    [2, 2],
)]
#[case::symmetric_constraint(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 3],
        constraints: vec![ConstraintSpan::Added(Constraint::Molecule(
            MoleculeConstraint::ChargeSum {
                atoms: None,
                sum: NumForm::Lit(0),
            },
        ))],
        ..Default::default()
    }),
    DescriptionLevel::Full,
    [6, 1],
)]
#[case::asymmetric_constraint(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 3],
        constraints: vec![ConstraintSpan::Added(Constraint::Molecule(
            MoleculeConstraint::ChargeSum {
                atoms: Some(vec![AtomId(0)]),
                sum: NumForm::Lit(0),
            },
        ))],
        ..Default::default()
    }),
    DescriptionLevel::Full,
    [6, 3],
)]
fn test_canonical_candidate(
    canonicalize_context: CanonicalizeContext,
    #[case] span: ReactionSpan,
    #[case] level: DescriptionLevel,
    #[case] expected_visited_leaves: [usize; 2],
) {
    let span = span.normalize().expect("fixed span normalizes");
    let incidence_graph = span.incidence_graph(match level {
        DescriptionLevel::Topology => IncidenceLevel::Topology,
        DescriptionLevel::Constitution => IncidenceLevel::Constitution,
        DescriptionLevel::Structure | DescriptionLevel::Full => IncidenceLevel::Full,
    });
    let [unpruned, pruned] = [false, true].map(|automorphism_pruning| {
        canonical_candidate(
            &span,
            level,
            &canonicalize_context,
            CanonicalSearchOptions {
                automorphism_pruning,
                branch_order: if automorphism_pruning {
                    BranchOrdering::BackendCanonical
                } else {
                    BranchOrdering::NodeId
                },
            },
        )
        .expect("fixed span has colors")
    });
    let canonical = |search: &CanonicalSearchResult<Option<ColorKey>>| {
        span.clone()
            .remap(&remapping_from_order(
                &span,
                &incidence_graph,
                &search.candidate.entity_order,
            ))
            .reframe()
            .expect("selected leaf preserves compatible frames")
    };

    assert_eq!(canonical(&pruned), canonical(&unpruned));
    assert_eq!(
        [unpruned.stats.visited_leaves, pruned.stats.visited_leaves],
        expected_visited_leaves,
    );
    for search in [&unpruned, &pruned] {
        assert_eq!(
            search.candidate.key.is_some(),
            search.stats.visited_leaves > 1
        );
    }
}

#[rstest]
#[case::topology(DescriptionLevel::Topology)]
#[case::constitution(DescriptionLevel::Constitution)]
#[case::structure(DescriptionLevel::Structure)]
#[case::full(DescriptionLevel::Full)]
fn test_canonicalize_checked_reaction_span_modified(
    canonicalize_context: CanonicalizeContext,
    #[case] level: DescriptionLevel,
) {
    let lhs_atom = AtomForm::from_element(Element::C).with_charge(NumForm::Lit(1));
    let rhs_atom = AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([1_i64]));
    let source = ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Modified {
            lhs: lhs_atom.clone(),
            rhs: rhs_atom,
        }],
        ..Default::default()
    });
    let expected = ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(lhs_atom)],
        ..Default::default()
    });

    assert_eq!(
        source
            .clone()
            .tracked_canonicalize_at(level, &canonicalize_context)
            .map(|(canonical, _)| canonical),
        Ok(expected),
    );
}
