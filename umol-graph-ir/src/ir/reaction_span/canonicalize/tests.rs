//! Unit tests for reaction-span and reaction canonicalization.

use rstest::{fixture, rstest};
use umol_chem::element::Element;
use umol_edn::FromEdn;

use super::*;
use crate::ir::canonicalize::tests::{canonicalize_context, explicitly_dense_equivalent};
use crate::ir::{
    ArithExpr, AromaticSystemForm, AromaticSystemId, AtomConstraintForm, AtomDelta,
    AtomFieldChange, AtomForm, AtomId, BondForm, BondId, Constraint, ConstraintDelta,
    DativeBondForm, DativeBondId, Delta, Deltas, ElementForm, Entity, IncidenceLevel, Molecule,
    MoleculeConstraint, MoleculeCorrespondence, MoleculeEntries, MulticenterBondForm,
    MulticenterBondId, NoncovalentBondDelta, NoncovalentBondForm, NoncovalentBondId,
    NoncovalentBondKind, NumForm, ReactionSpanEntries, RelationalConstraint, StereoAtomDelta,
    StereoAtomForm, StereoAtomId, StereoBondForm, StereoBondId, StereoCoset, StereoKind,
    StereoLigand, StereoLigandKind,
};
use crate::mol_dsl;

#[fixture]
fn initial_color_molecule() -> Molecule {
    let normalized_three = NumForm::ArithExpr(Box::new(ArithExpr::Sum(vec![
        ArithExpr::Lit(1),
        ArithExpr::Lit(2),
    ])));
    let normalized_one = NumForm::ArithExpr(Box::new(ArithExpr::Sum(vec![
        ArithExpr::Lit(0),
        ArithExpr::Lit(1),
    ])));
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_charge(3_i64),
            AtomForm::from_element(Element::C)
                .with_charge(normalized_three)
                .with_constraint(AtomConstraintForm::Valence(NumForm::Lit(4))),
            AtomForm::from_element(Element::O).with_charge(3_i64),
            AtomForm::from_element(Element::C).with_charge(4_i64),
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::H),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::new(normalized_one.clone())),
            (AtomId(2), AtomId(3), BondForm::from_order(2)),
            (AtomId(0), AtomId(5), BondForm::from_order(1)),
        ],
        dative: vec![
            (vec![AtomId(0)], AtomId(4), DativeBondForm::from_order(1)),
            (
                vec![AtomId(1)],
                AtomId(4),
                DativeBondForm::new(normalized_one),
            ),
            (vec![AtomId(2)], AtomId(4), DativeBondForm::from_order(2)),
        ],
        aromatic: vec![
            (
                vec![AtomId(0), AtomId(1)],
                AromaticSystemForm::from_electrons(vec![1, 2]),
            ),
            (
                vec![AtomId(2), AtomId(3)],
                AromaticSystemForm::from_electrons(vec![2, 1]),
            ),
            (
                vec![AtomId(4), AtomId(5)],
                AromaticSystemForm::from_electrons(vec![1, 2]).with_charge(1_i64),
            ),
        ],
        multicenter: vec![
            (
                vec![AtomId(0), AtomId(2)],
                MulticenterBondForm::from_electrons(vec![1, 2]),
            ),
            (
                vec![AtomId(1), AtomId(3)],
                MulticenterBondForm::from_electrons(vec![2, 1]),
            ),
        ],
        noncovalent: vec![
            (
                [AtomId(0), AtomId(5)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            ),
            (
                [AtomId(1), AtomId(5)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            ),
            (
                [AtomId(2), AtomId(5)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::Ionic),
            ),
        ],
        stereo_atoms: vec![
            (
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            (
                AtomId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            ),
            (
                AtomId(2),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
                ],
                StereoAtomForm::new(StereoKind::SquarePlanar, StereoCoset::Lit(0)),
            ),
        ],
        stereo_bonds: vec![
            (
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
            ),
            (
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            ),
        ],
        ..Default::default()
    })
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
#[case::unchanged(
    EntitySpan::Unchanged(AtomForm::new(ElementForm::lit_set([Element::C]))),
    EntitySpan::Unchanged(AtomForm::from_element(Element::C))
)]
#[case::removed(
    EntitySpan::Removed(AtomForm::new(ElementForm::lit_set([Element::C]))),
    EntitySpan::Removed(AtomForm::from_element(Element::C))
)]
#[case::added(
    EntitySpan::Added(AtomForm::new(ElementForm::lit_set([Element::C]))),
    EntitySpan::Added(AtomForm::from_element(Element::C))
)]
#[case::modified_raw(
    EntitySpan::Modified {
        lhs: AtomForm::from_element(Element::C),
        rhs: AtomForm::new(ElementForm::lit_set([Element::C])),
    },
    EntitySpan::Unchanged(AtomForm::from_element(Element::C))
)]
#[case::modified_constraint(
    EntitySpan::Modified {
        lhs: AtomForm::from_element(Element::C),
        rhs: AtomForm::from_element(Element::C)
            .with_constraint(AtomConstraintForm::Valence(NumForm::lit_set([4]))),
    },
    EntitySpan::Modified {
        lhs: AtomForm::from_element(Element::C),
        rhs: AtomForm::from_element(Element::C)
            .with_constraint(AtomConstraintForm::Valence(NumForm::Lit(4))),
    }
)]
fn test_reaction_span_canonicalize_colors(
    canonicalize_context: CanonicalizeContext,
    #[case] atom: EntitySpan<AtomForm>,
    #[case] expected_atom: EntitySpan<AtomForm>,
) {
    let span = ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![atom],
        ..Default::default()
    });
    let expected = ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![expected_atom],
        ..Default::default()
    });
    assert_eq!(
        span.clone().canonicalize(&canonicalize_context),
        Ok(expected.clone())
    );
    assert!(span.canonical_eq(&expected, &canonicalize_context));
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
fn test_same_entity_color(
    initial_color_span: ReactionSpan,
    #[case] entity: Entity,
    #[case] image: Entity,
    #[case] expected: bool,
) {
    let span = initial_color_span
        .normalize()
        .expect("fixed span normalizes");

    assert_eq!(same_entity_color(&span, entity, &span, image), expected);
}

#[rstest]
#[case::renumbered(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Removed(AtomForm::from_element(Element::F)),
            EntitySpan::Added(AtomForm::from_element(Element::O)),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Removed(BondForm::from_order(1))),
            (AtomId(0), AtomId(2), EntitySpan::Added(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Added(AtomForm::from_element(Element::O)),
            EntitySpan::Removed(AtomForm::from_element(Element::F)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
        ],
        bonds: vec![
            (AtomId(0), AtomId(2), EntitySpan::Added(BondForm::from_order(1))),
            (AtomId(1), AtomId(2), EntitySpan::Removed(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    true
)]
#[case::added_and_removed_exchanged(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Removed(AtomForm::from_element(Element::F)),
            EntitySpan::Added(AtomForm::from_element(Element::O)),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Removed(BondForm::from_order(1))),
            (AtomId(0), AtomId(2), EntitySpan::Added(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Added(AtomForm::from_element(Element::F)),
            EntitySpan::Removed(AtomForm::from_element(Element::O)),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Added(BondForm::from_order(1))),
            (AtomId(0), AtomId(2), EntitySpan::Removed(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    false
)]
#[case::different_element(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Removed(AtomForm::from_element(Element::F)),
            EntitySpan::Added(AtomForm::from_element(Element::O)),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Removed(BondForm::from_order(1))),
            (AtomId(0), AtomId(2), EntitySpan::Added(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Added(AtomForm::from_element(Element::S)),
            EntitySpan::Removed(AtomForm::from_element(Element::F)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
        ],
        bonds: vec![
            (AtomId(0), AtomId(2), EntitySpan::Added(BondForm::from_order(1))),
            (AtomId(1), AtomId(2), EntitySpan::Removed(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    false
)]
#[case::modified_and_unchanged(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Modified {
                lhs: AtomForm::from_element(Element::C).with_charge(0_i64),
                rhs: AtomForm::from_element(Element::C).with_charge(1_i64),
            },
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C).with_charge(0_i64)),
        ],
        ..Default::default()
    }),
    false
)]
#[case::electron_count_change_reversed(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::O)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
        ],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4)],
            EntitySpan::Modified {
                lhs: AromaticSystemForm::from_electrons(vec![2, 1, 1, 1, 1]),
                rhs: AromaticSystemForm::from_electrons(vec![1, 1, 1, 1, 1]),
            },
        )],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::O)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
        ],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4)],
            EntitySpan::Modified {
                lhs: AromaticSystemForm::from_electrons(vec![1, 1, 1, 1, 1]),
                rhs: AromaticSystemForm::from_electrons(vec![2, 1, 1, 1, 1]),
            },
        )],
        ..Default::default()
    }),
    false
)]
fn test_same_colors(
    #[case] source: ReactionSpan,
    #[case] target: ReactionSpan,
    #[case] expected: bool,
) {
    let source = source.normalize().unwrap();
    let target = target.normalize().unwrap();
    let source_incidence = source.incidence_graph(IncidenceLevel::Constitution);
    let target_incidence = target.incidence_graph(IncidenceLevel::Constitution);
    let source_colors = initial_colors(
        &source,
        &source_incidence,
        DescriptionLevel::Constitution,
        true,
    )
    .unwrap();
    let target_colors = initial_colors(
        &target,
        &target_incidence,
        DescriptionLevel::Constitution,
        true,
    )
    .unwrap();

    assert_eq!(
        same_colors(
            &source,
            &source_incidence,
            &source_colors,
            &target,
            &target_incidence,
            &target_colors,
        ),
        expected,
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
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 3],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Added(BondForm::from_order(2))),
            (AtomId(1), AtomId(2), EntitySpan::Removed(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    Graph::new(4, &[[0, 3], [3, 1], [1, 2]]),
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
fn test_topology_search_input(
    #[case] span: ReactionSpan,
    #[case] expected_graph: Graph,
    #[case] expected_colors: Vec<AutomorphismColor>,
    #[case] expected_node_sources: Vec<SubdivisionNodeSource>,
    #[case] expected_partition: Vec<Vec<NodeId>>,
) {
    let span = span.normalize().unwrap();
    let incidence_graph = span.incidence_graph(IncidenceLevel::Topology);
    let colors = initial_colors(&span, &incidence_graph, DescriptionLevel::Topology, true).unwrap();
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
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 3],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Added(BondForm::from_order(2))),
            (AtomId(1), AtomId(2), EntitySpan::Removed(BondForm::from_order(1))),
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
    1,
)]
#[case::added_atoms(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Added(AtomForm::from_element(Element::C)); 3],
        ..Default::default()
    }),
    DescriptionLevel::Topology,
    1,
)]
#[case::removed_atoms(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Removed(AtomForm::from_element(Element::C)); 3],
        ..Default::default()
    }),
    DescriptionLevel::Topology,
    1,
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
    1,
)]
#[case::removed_and_added_bond(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 3],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Added(BondForm::from_order(2))),
            (AtomId(1), AtomId(2), EntitySpan::Removed(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    DescriptionLevel::Topology,
    1,
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
    1,
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
    1,
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
    1,
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
    1,
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
    1,
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
    1,
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
    3,
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
    3,
)]
#[case::prostereogenic_rhs_configuration(
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
    2,
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
    1,
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
    3,
)]
fn test_canonical_candidate(
    canonicalize_context: CanonicalizeContext,
    #[case] span: ReactionSpan,
    #[case] level: DescriptionLevel,
    #[case] expected_visited_leaves: usize,
) {
    let span = span.normalize().expect("fixed span normalizes");
    let search =
        canonical_candidate(&span, level, &canonicalize_context).expect("fixed span has colors");

    assert_eq!(search.stats.visited_leaves, expected_visited_leaves);
    assert_eq!(
        search.candidate.key.is_some(),
        search.stats.visited_leaves > 1
    );
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

#[rstest]
fn test_reaction_span_canonicalize_contradiction(canonicalize_context: CanonicalizeContext) {
    let span = ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(
            AtomForm::from_element(Element::C).with_charge(NumForm::lit_set(Vec::<i64>::new())),
        )],
        ..Default::default()
    });

    assert_eq!(
        span.canonicalize(&canonicalize_context),
        Err(ReactionSpanCanonicalizeError::Contradiction(Contradiction)),
    );
}

#[rstest]
#[case::atom_charge_unchanged(
    r#"{:atoms ["C#c+" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C#c+" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::atom_charge_unchanged_absent(
    r#"{:atoms ["C#c+" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::atom_charge_modified(
    r#"{:atoms [{:modify ["C" "C#c+"]} "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" {:modify ["C" "C#c+"]} "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::atom_charge_modified_reversed(
    r#"{:atoms [{:modify ["C" "C#c+"]} "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms [{:modify ["C#c+" "C"]} "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::atom_charge_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:add "C#c+"}]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:add [0 6 "1"]}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:add "C#c+"}]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:add [3 6 "1"]}]}"#,
    true
)]
#[case::atom_charge_added_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:add "C#c+"}]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:add [0 6 "1"]}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:add "C"}]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:add [0 6 "1"]}]}"#,
    false
)]
#[case::atom_charge_removed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:remove "C#c+"}]
          :bonds
          [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:remove [0 6 "1"]}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:remove "C#c+"}]
          :bonds
          [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:remove [3 6 "1"]}]}"#,
    true
)]
#[case::atom_charge_removed_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:remove "C#c+"}]
          :bonds
          [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:remove [0 6 "1"]}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:add "C#c+"}]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:add [0 6 "1"]}]}"#,
    false
)]
#[case::atom_constraint_unchanged(
    r#"{:atoms ["C#v4" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C#v4" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::atom_constraint_unchanged_absent(
    r#"{:atoms ["C#v4" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::atom_constraint_modified(
    r#"{:atoms [{:modify ["C" "C#v4"]} "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" {:modify ["C" "C#v4"]} "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::atom_constraint_modified_reversed(
    r#"{:atoms [{:modify ["C" "C#v4"]} "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms [{:modify ["C#v4" "C"]} "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::atom_constraint_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:add "C#v4"}]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:add [0 6 "1"]}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:add "C#v4"}]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:add [3 6 "1"]}]}"#,
    true
)]
#[case::atom_constraint_added_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:add "C#v4"}]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:add [0 6 "1"]}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:add "C"}]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:add [0 6 "1"]}]}"#,
    false
)]
#[case::atom_constraint_removed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:remove "C#v4"}]
          :bonds
          [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:remove [0 6 "1"]}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:remove "C#v4"}]
          :bonds
          [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:remove [3 6 "1"]}]}"#,
    true
)]
#[case::atom_constraint_removed_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:remove "C#v4"}]
          :bonds
          [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:remove [0 6 "1"]}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C" {:add "C#v4"}]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"] {:add [0 6 "1"]}]}"#,
    false
)]
#[case::bond_order_unchanged(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "2"] [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::bond_order_unchanged_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::bond_order_modified(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:modify [0 1 ["1" "2"]]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] {:modify [3 4 ["1" "2"]]} [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::bond_order_modified_reversed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:modify [0 1 ["1" "2"]]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:modify [0 1 ["2" "1"]]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::bond_order_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "2"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] {:add [3 4 "2"]} [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::bond_order_added_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "2"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "1"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::bond_order_removed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:remove [0 1 "2"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] {:remove [3 4 "2"]} [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::bond_order_removed_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:remove [0 1 "2"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "2"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::bond_constraint_unchanged(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1#a"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1#a"] [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::bond_constraint_unchanged_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1#a"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::bond_constraint_modified(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:modify [0 1 ["1" "1#a"]]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] {:modify [3 4 ["1" "1#a"]]} [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::bond_constraint_modified_reversed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:modify [0 1 ["1" "1#a"]]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:modify [0 1 ["1#a" "1"]]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::bond_constraint_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "1#a"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] {:add [3 4 "1#a"]} [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::bond_constraint_added_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "1#a"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "1"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::bond_constraint_removed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:remove [0 1 "1#a"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] {:remove [3 4 "1#a"]} [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::bond_constraint_removed_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:remove [0 1 "1#a"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "1#a"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::atom_added(
    r#"{:atoms [{:add "C"} "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "1"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] {:add [0 5 "1"]}]}"#,
    r#"{:atoms ["C" "C" "C" {:add "C"} "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] {:add [2 3 "1"]} {:add [3 4 "1"]} [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::atom_added_removed(
    r#"{:atoms [{:add "C"} "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "1"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] {:add [0 5 "1"]}]}"#,
    r#"{:atoms [{:remove "C"} "C" "C" "C" "C" "C"]
          :bonds
          [{:remove [0 1 "1"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] {:remove [0 5 "1"]}]}"#,
    false
)]
#[case::bond_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "1"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] {:add [3 4 "1"]} [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::bond_added_pair(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "1"]} {:add [1 2 "1"]} [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "1"]} [1 2 "1"] [2 3 "1"] {:add [3 4 "1"]} [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::bond_removed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:remove [0 1 "1"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] {:remove [3 4 "1"]} [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::bond_removed_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:remove [0 1 "1"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [{:add [0 1 "1"]} [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::dative_bond_unchanged(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:donors [0] :acceptor 3 :attrs "1"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:donors [1] :acceptor 4 :attrs "1"}]}"#,
    true
)]
#[case::dative_bond_unchanged_distance(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:donors [0] :acceptor 3 :attrs "1"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:donors [0] :acceptor 2 :attrs "1"}]}"#,
    false
)]
#[case::dative_bond_modified(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:modify {:donors [0] :acceptor 3 :attrs ["1" "2"]}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:modify {:donors [1] :acceptor 4 :attrs ["1" "2"]}}]}"#,
    true
)]
#[case::dative_bond_modified_reversed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:modify {:donors [0] :acceptor 3 :attrs ["1" "2"]}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:modify {:donors [0] :acceptor 3 :attrs ["2" "1"]}}]}"#,
    false
)]
#[case::dative_bond_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:add {:donors [0] :acceptor 3 :attrs "1"}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:add {:donors [1] :acceptor 4 :attrs "1"}}]}"#,
    true
)]
#[case::dative_bond_added_distance(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:add {:donors [0] :acceptor 3 :attrs "1"}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:add {:donors [0] :acceptor 2 :attrs "1"}}]}"#,
    false
)]
#[case::dative_bond_removed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:remove {:donors [0] :acceptor 3 :attrs "1"}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:remove {:donors [1] :acceptor 4 :attrs "1"}}]}"#,
    true
)]
#[case::dative_bond_removed_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:remove {:donors [0] :acceptor 3 :attrs "1"}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:add {:donors [0] :acceptor 3 :attrs "1"}}]}"#,
    false
)]
#[case::aromatic_electron_counts_unchanged(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[2,1,1,1,1,1]"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,2,1,1]"}]}"#,
    true
)]
#[case::aromatic_electron_counts_unchanged_pair(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[2,2,1,1,1,1]"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[2,1,1,2,1,1]"}]}"#,
    false
)]
#[case::aromatic_electron_counts_modified(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems
          [{:modify {:atoms [0 1 2 3 4 5] :attrs ["[1,1,1,1,1,1]" "[2,1,1,1,1,1]"]}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems
          [{:modify {:atoms [0 1 2 3 4 5] :attrs ["[1,1,1,1,1,1]" "[1,1,1,2,1,1]"]}}]}"#,
    true
)]
#[case::aromatic_electron_counts_modified_shift(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems
          [{:modify {:atoms [0 1 2 3 4 5] :attrs ["[2,1,1,1,1,1]" "[1,2,1,1,1,1]"]}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems
          [{:modify {:atoms [0 1 2 3 4 5] :attrs ["[2,1,1,1,1,1]" "[1,1,1,2,1,1]"]}}]}"#,
    false
)]
#[case::aromatic_electron_counts_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2 3 4 5] :attrs "[2,1,1,1,1,1]"}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2 3 4 5] :attrs "[1,1,1,2,1,1]"}}]}"#,
    true
)]
#[case::aromatic_electron_counts_added_pair(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2 3 4 5] :attrs "[2,2,1,1,1,1]"}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2 3 4 5] :attrs "[2,1,1,2,1,1]"}}]}"#,
    false
)]
#[case::aromatic_electron_counts_removed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:remove {:atoms [0 1 2 3 4 5] :attrs "[2,1,1,1,1,1]"}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:remove {:atoms [0 1 2 3 4 5] :attrs "[1,1,1,2,1,1]"}}]}"#,
    true
)]
#[case::aromatic_electron_counts_removed_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:remove {:atoms [0 1 2 3 4 5] :attrs "[2,1,1,1,1,1]"}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2 3 4 5] :attrs "[2,1,1,1,1,1]"}}]}"#,
    false
)]
#[case::aromatic_system_charge_unchanged(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*#c-"} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*"} {:atoms [3 4 5] :attrs "*#c-"}]}"#,
    true
)]
#[case::aromatic_system_charge_unchanged_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*#c-"} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*"} {:atoms [3 4 5] :attrs "*"}]}"#,
    false
)]
#[case::aromatic_system_charge_modified(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:modify {:atoms [0 1 2] :attrs ["*" "*#c-"]}} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:atoms [0 1 2] :attrs "*"} {:modify {:atoms [3 4 5] :attrs ["*" "*#c-"]}}]}"#,
    true
)]
#[case::aromatic_system_charge_modified_reversed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:modify {:atoms [0 1 2] :attrs ["*" "*#c-"]}} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:modify {:atoms [0 1 2] :attrs ["*#c-" "*"]}} {:atoms [3 4 5] :attrs "*"}]}"#,
    false
)]
#[case::aromatic_system_charge_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2] :attrs "*#c-"}} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*"} {:add {:atoms [3 4 5] :attrs "*#c-"}}]}"#,
    true
)]
#[case::aromatic_system_charge_added_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2] :attrs "*#c-"}} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2] :attrs "*"}} {:atoms [3 4 5] :attrs "*"}]}"#,
    false
)]
#[case::aromatic_system_charge_removed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:remove {:atoms [0 1 2] :attrs "*#c-"}} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:atoms [0 1 2] :attrs "*"} {:remove {:atoms [3 4 5] :attrs "*#c-"}}]}"#,
    true
)]
#[case::aromatic_system_charge_removed_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:remove {:atoms [0 1 2] :attrs "*#c-"}} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2] :attrs "*#c-"}} {:atoms [3 4 5] :attrs "*"}]}"#,
    false
)]
#[case::aromatic_system_constraint_unchanged(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*#e2"} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*"} {:atoms [3 4 5] :attrs "*#e2"}]}"#,
    true
)]
#[case::aromatic_system_constraint_unchanged_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*#e2"} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*"} {:atoms [3 4 5] :attrs "*"}]}"#,
    false
)]
#[case::aromatic_system_constraint_modified(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:modify {:atoms [0 1 2] :attrs ["*" "*#e2"]}} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:atoms [0 1 2] :attrs "*"} {:modify {:atoms [3 4 5] :attrs ["*" "*#e2"]}}]}"#,
    true
)]
#[case::aromatic_system_constraint_modified_reversed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:modify {:atoms [0 1 2] :attrs ["*" "*#e2"]}} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:modify {:atoms [0 1 2] :attrs ["*#e2" "*"]}} {:atoms [3 4 5] :attrs "*"}]}"#,
    false
)]
#[case::aromatic_system_constraint_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2] :attrs "*#e2"}} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*"} {:add {:atoms [3 4 5] :attrs "*#e2"}}]}"#,
    true
)]
#[case::aromatic_system_constraint_added_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2] :attrs "*#e2"}} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2] :attrs "*"}} {:atoms [3 4 5] :attrs "*"}]}"#,
    false
)]
#[case::aromatic_system_constraint_removed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:remove {:atoms [0 1 2] :attrs "*#e2"}} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:atoms [0 1 2] :attrs "*"} {:remove {:atoms [3 4 5] :attrs "*#e2"}}]}"#,
    true
)]
#[case::aromatic_system_constraint_removed_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems
          [{:remove {:atoms [0 1 2] :attrs "*#e2"}} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:add {:atoms [0 1 2] :attrs "*#e2"}} {:atoms [3 4 5] :attrs "*"}]}"#,
    false
)]
#[case::multicenter_electron_counts_unchanged(
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:atoms [0 1 2] :attrs "[2,0,0]"}]}"#,
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:atoms [0 1 2] :attrs "[0,0,2]"}]}"#,
    true
)]
#[case::multicenter_electron_counts_unchanged_distinct(
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:atoms [0 1 2] :attrs "[2,0,0]"}]}"#,
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:atoms [0 1 2] :attrs "[1,1,0]"}]}"#,
    false
)]
#[case::multicenter_electron_counts_modified(
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:modify {:atoms [0 1 2] :attrs ["[2,0,0]" "[0,2,0]"]}}]}"#,
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:modify {:atoms [0 1 2] :attrs ["[0,0,2]" "[2,0,0]"]}}]}"#,
    true
)]
#[case::multicenter_electron_counts_modified_distinct(
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:modify {:atoms [0 1 2] :attrs ["[2,0,0]" "[0,2,0]"]}}]}"#,
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:modify {:atoms [0 1 2] :attrs ["[2,0,0]" "[1,1,0]"]}}]}"#,
    false
)]
#[case::multicenter_electron_counts_added(
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:add {:atoms [0 1 2] :attrs "[2,0,0]"}}]}"#,
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:add {:atoms [0 1 2] :attrs "[0,0,2]"}}]}"#,
    true
)]
#[case::multicenter_electron_counts_added_distinct(
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:add {:atoms [0 1 2] :attrs "[2,0,0]"}}]}"#,
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:add {:atoms [0 1 2] :attrs "[1,1,0]"}}]}"#,
    false
)]
#[case::multicenter_electron_counts_removed(
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:remove {:atoms [0 1 2] :attrs "[2,0,0]"}}]}"#,
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:remove {:atoms [0 1 2] :attrs "[0,0,2]"}}]}"#,
    true
)]
#[case::multicenter_electron_counts_removed_added(
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:remove {:atoms [0 1 2] :attrs "[2,0,0]"}}]}"#,
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:add {:atoms [0 1 2] :attrs "[2,0,0]"}}]}"#,
    false
)]
#[case::noncovalent_bond_unchanged(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:atoms [0 3] :attrs "Hbd"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:atoms [1 4] :attrs "Hbd"}]}"#,
    true
)]
#[case::noncovalent_bond_unchanged_distance(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:atoms [0 3] :attrs "Hbd"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:atoms [0 2] :attrs "Hbd"}]}"#,
    false
)]
#[case::noncovalent_bond_modified(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:modify {:atoms [0 3] :attrs ["Hbd" "Xbd"]}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:modify {:atoms [1 4] :attrs ["Hbd" "Xbd"]}}]}"#,
    true
)]
#[case::noncovalent_bond_modified_reversed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:modify {:atoms [0 3] :attrs ["Hbd" "Xbd"]}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:modify {:atoms [0 3] :attrs ["Xbd" "Hbd"]}}]}"#,
    false
)]
#[case::noncovalent_bond_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:add {:atoms [0 3] :attrs "Hbd"}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:add {:atoms [1 4] :attrs "Hbd"}}]}"#,
    true
)]
#[case::noncovalent_bond_added_distance(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:add {:atoms [0 3] :attrs "Hbd"}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:add {:atoms [0 2] :attrs "Hbd"}}]}"#,
    false
)]
#[case::noncovalent_bond_removed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:remove {:atoms [0 3] :attrs "Hbd"}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:remove {:atoms [1 4] :attrs "Hbd"}}]}"#,
    true
)]
#[case::noncovalent_bond_removed_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:remove {:atoms [0 3] :attrs "Hbd"}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:add {:atoms [0 3] :attrs "Hbd"}}]}"#,
    false
)]
#[case::stereo_configuration_unchanged(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th1"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    true
)]
#[case::stereo_configuration_unchanged_distinct(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th1"}]}"#,
    false
)]
#[case::stereo_configuration_modified(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:modify {:site 0 :ligands [1 2 3 4] :attrs ["Th0" "Th1"]}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:modify {:site 4 :ligands [5 6 7 0] :attrs ["Th0" "Th1"]}}]}"#,
    true
)]
#[case::stereo_configuration_modified_reversed(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:modify {:site 0 :ligands [1 2 3 4] :attrs ["Th0" "Th1"]}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:modify {:site 0 :ligands [1 2 3 4] :attrs ["Th1" "Th0"]}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    false
)]
#[case::stereo_configuration_added(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4] :attrs "Th0"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th1"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}
                         {:add {:site 4 :ligands [5 6 7 0] :attrs "Th0"}}]}"#,
    true
)]
#[case::stereo_configuration_added_distinct(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4] :attrs "Th0"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th1"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4] :attrs "Th1"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th1"}]}"#,
    false
)]
#[case::stereo_configuration_removed(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:remove {:site 0 :ligands [1 2 3 4] :attrs "Th0"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th1"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}
                         {:remove {:site 4 :ligands [5 6 7 0] :attrs "Th0"}}]}"#,
    true
)]
#[case::stereo_configuration_removed_added(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:remove {:site 0 :ligands [1 2 3 4] :attrs "Th0"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th1"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4] :attrs "Th0"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th1"}]}"#,
    false
)]
#[case::stereo_stereogenicity_unchanged(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0#g/"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0#g/"}]}"#,
    true
)]
#[case::stereo_stereogenicity_unchanged_absent(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0#g/"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    false
)]
#[case::stereo_stereogenicity_modified(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:modify {:site 0 :ligands [1 2 3 4] :attrs ["Th0" "Th0#g/"]}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:modify {:site 4 :ligands [5 6 7 0] :attrs ["Th0" "Th0#g/"]}}]}"#,
    true
)]
#[case::stereo_stereogenicity_modified_reversed(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:modify {:site 0 :ligands [1 2 3 4] :attrs ["Th0" "Th0#g/"]}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:modify {:site 0 :ligands [1 2 3 4] :attrs ["Th0#g/" "Th0"]}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    false
)]
#[case::stereo_stereogenicity_added(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4] :attrs "Th0#g/"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:add {:site 4 :ligands [5 6 7 0] :attrs "Th0#g/"}}]}"#,
    true
)]
#[case::stereo_stereogenicity_added_absent(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4] :attrs "Th0#g/"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4] :attrs "Th0"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    false
)]
#[case::stereo_stereogenicity_removed(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:remove {:site 0 :ligands [1 2 3 4] :attrs "Th0#g/"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:remove {:site 4 :ligands [5 6 7 0] :attrs "Th0#g/"}}]}"#,
    true
)]
#[case::stereo_stereogenicity_removed_added(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:remove {:site 0 :ligands [1 2 3 4] :attrs "Th0#g/"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4] :attrs "Th0#g/"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    false
)]
#[case::stereo_topicity_unchanged(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0#o(0,1)'"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0#o(0,1)'"}]}"#,
    true
)]
#[case::stereo_topicity_unchanged_pair(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0#o(0,1)'"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0#o(0,2)'"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    false
)]
#[case::stereo_topicity_modified(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:modify {:site 0 :ligands [1 2 3 4] :attrs ["Th0" "Th0#o(0,1)'"]}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:modify {:site 4 :ligands [5 6 7 0] :attrs ["Th0" "Th0#o(0,1)'"]}}]}"#,
    true
)]
#[case::stereo_topicity_modified_pair(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:modify {:site 0 :ligands [1 2 3 4] :attrs ["Th0" "Th0#o(0,1)'"]}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:modify {:site 0 :ligands [1 2 3 4] :attrs ["Th0" "Th0#o(0,2)'"]}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    false
)]
#[case::stereo_topicity_added(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4] :attrs "Th0#o(0,1)'"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:add {:site 4 :ligands [5 6 7 0] :attrs "Th0#o(0,1)'"}}]}"#,
    true
)]
#[case::stereo_topicity_added_pair(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4] :attrs "Th0#o(0,1)'"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4] :attrs "Th0#o(0,2)'"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    false
)]
#[case::stereo_topicity_removed(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:remove {:site 0 :ligands [1 2 3 4] :attrs "Th0#o(0,1)'"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:remove {:site 4 :ligands [5 6 7 0] :attrs "Th0#o(0,1)'"}}]}"#,
    true
)]
#[case::stereo_topicity_removed_added(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:remove {:site 0 :ligands [1 2 3 4] :attrs "Th0#o(0,1)'"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:add {:site 0 :ligands [1 2 3 4] :attrs "Th0#o(0,1)'"}}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    false
)]
#[case::constraint_unchanged(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [0] :sum 0}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [3] :sum 0}}]}"#,
    true
)]
#[case::constraint_unchanged_pair(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [0 1] :sum 0}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [0 3] :sum 0}}]}"#,
    false
)]
#[case::constraint_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:add {:charge-sum {:atoms [0] :sum 0}}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:add {:charge-sum {:atoms [3] :sum 0}}}]}"#,
    true
)]
#[case::constraint_added_pair(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:add {:charge-sum {:atoms [0 1] :sum 0}}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:add {:charge-sum {:atoms [0 3] :sum 0}}}]}"#,
    false
)]
#[case::constraint_removed(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:remove {:charge-sum {:atoms [0] :sum 0}}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:remove {:charge-sum {:atoms [3] :sum 0}}}]}"#,
    true
)]
#[case::constraint_removed_added(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:remove {:charge-sum {:atoms [0] :sum 0}}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:add {:charge-sum {:atoms [0] :sum 0}}}]}"#,
    false
)]
fn test_reaction_span_canonicalize_placement(
    canonicalize_context: CanonicalizeContext,
    #[case] left: &str,
    #[case] right: &str,
    #[case] expected: bool,
) {
    let [left, right] =
        [left, right].map(|source| ReactionSpan::from_edn_str(source).expect("case span parses"));

    assert_eq!(
        left.clone().canonicalize(&canonicalize_context)
            == right.clone().canonicalize(&canonicalize_context),
        expected
    );
    assert_eq!(left.canonical_eq(&right, &canonicalize_context), expected);
    if expected {
        assert_eq!(
            left.canonical_hash(&canonicalize_context),
            right.canonical_hash(&canonicalize_context),
        );
    }
}

fn reaction_canonicalization_fixture() -> Reaction {
    let atom = AtomForm::from_element(Element::C).with_charge(0_i64);
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![atom.clone(), atom],
            ..Default::default()
        }),
        [Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(1),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(0),
                new: NumForm::Lit(1),
            },
        })]
        .into_iter()
        .collect(),
    )
}

#[rustfmt::skip]
#[rstest]
#[case::molecule_subset(ReactionSpan::from_entries(ReactionSpanEntries {
    atoms: vec![
        EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
        EntitySpan::Unchanged(AtomForm::from_element(Element::N)),
        EntitySpan::Unchanged(AtomForm::from_element(Element::O)),
    ],
    constraints: vec![ConstraintSpan::Unchanged(Constraint::Molecule(MoleculeConstraint::ChargeSum {
        atoms: Some(vec![AtomId(0), AtomId(1)]),
        sum: NumForm::Lit(0),
    }))],
    ..Default::default()
}))]
#[case::relational_atom_list(ReactionSpan::from_entries(ReactionSpanEntries {
    atoms: vec![
        EntitySpan::Unchanged(AtomForm::from_element(Element::Fe)),
        EntitySpan::Unchanged(AtomForm::from_element(Element::N)),
        EntitySpan::Unchanged(AtomForm::from_element(Element::O)),
    ],
    dative: vec![(vec![AtomId(1), AtomId(2)], AtomId(0), EntitySpan::Unchanged(DativeBondForm::from_order(1)))],
    constraints: vec![ConstraintSpan::Unchanged(Constraint::Relational(RelationalConstraint::DativeBondDonors {
        bond: DativeBondId(0),
        atoms: vec![AtomId(1), AtomId(2)],
    }))],
    ..Default::default()
}))]
fn test_reaction_span_canonical_eq_renumbering(
    canonicalize_context: CanonicalizeContext,
    #[case] span: ReactionSpan,
) {
    let images = span.entity_counts().map(|count| (0..count).rev().collect::<Vec<_>>());
    let renumbered = span.clone().remap(&molecule_remapping(&images));

    assert!(span.canonical_eq(&renumbered, &canonicalize_context));
    assert!(span.to_reaction().canonical_eq(&renumbered.to_reaction(), &canonicalize_context));
}

#[rstest]
#[case::one_fails(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([]))),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
        ],
        ..Default::default()
    }),
    false
)]
#[case::both_fail_identical(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([]))),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([]))),
        ],
        ..Default::default()
    }),
    true
)]
#[case::both_fail_equal_normal_forms(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Modified { lhs: AtomForm::from_element(Element::C), rhs: AtomForm::from_element(Element::C) },
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
        ],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            EntitySpan::Unchanged(AromaticSystemForm::from_electrons(vec![2, 1])),
        )],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
        ],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            EntitySpan::Unchanged(AromaticSystemForm::from_electrons(vec![2, 1])),
        )],
        ..Default::default()
    }),
    true
)]
#[case::both_fail_different_frames(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
        ],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            EntitySpan::Unchanged(AromaticSystemForm::from_electrons(vec![2, 1])),
        )],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
        ],
        aromatic: vec![(
            vec![AtomId(1), AtomId(0), AtomId(2)],
            EntitySpan::Unchanged(AromaticSystemForm::from_electrons(vec![2, 1])),
        )],
        ..Default::default()
    }),
    false
)]
#[case::both_fail_different(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([]))),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::N).with_charge(NumForm::lit_set([]))),
        ],
        ..Default::default()
    }),
    false
)]
#[case::both_fail_different_counts(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([]))),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::N).with_charge(NumForm::lit_set([]))),
            EntitySpan::Unchanged(AtomForm::from_element(Element::O)),
        ],
        ..Default::default()
    }),
    false
)]
fn test_reaction_span_canonical_eq_contradiction(
    canonicalize_context: CanonicalizeContext,
    #[case] left: ReactionSpan,
    #[case] right: ReactionSpan,
    #[case] expected: bool,
) {
    assert_eq!(left.canonical_eq(&right, &canonicalize_context), expected);
}

#[rstest]
#[case::hydrolysis_renumbered(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Removed(AtomForm::from_element(Element::F)),
            EntitySpan::Added(AtomForm::from_element(Element::O)),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Removed(BondForm::from_order(1))),
            (AtomId(0), AtomId(2), EntitySpan::Added(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Added(AtomForm::from_element(Element::O)),
            EntitySpan::Removed(AtomForm::from_element(Element::F)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
        ],
        bonds: vec![
            (AtomId(0), AtomId(2), EntitySpan::Added(BondForm::from_order(1))),
            (AtomId(1), AtomId(2), EntitySpan::Removed(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    true
)]
#[case::hydrolysis_thiolysis(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Removed(AtomForm::from_element(Element::F)),
            EntitySpan::Added(AtomForm::from_element(Element::O)),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Removed(BondForm::from_order(1))),
            (AtomId(0), AtomId(2), EntitySpan::Added(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Added(AtomForm::from_element(Element::S)),
            EntitySpan::Removed(AtomForm::from_element(Element::F)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
        ],
        bonds: vec![
            (AtomId(0), AtomId(2), EntitySpan::Added(BondForm::from_order(1))),
            (AtomId(1), AtomId(2), EntitySpan::Removed(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    false
)]
#[case::hydrolysis_reverse_reaction(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Removed(AtomForm::from_element(Element::F)),
            EntitySpan::Added(AtomForm::from_element(Element::O)),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Removed(BondForm::from_order(1))),
            (AtomId(0), AtomId(2), EntitySpan::Added(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Added(AtomForm::from_element(Element::F)),
            EntitySpan::Removed(AtomForm::from_element(Element::O)),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Added(BondForm::from_order(1))),
            (AtomId(0), AtomId(2), EntitySpan::Removed(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    false
)]
#[case::hexane_2_methylpentane(
    Reaction::new(
        mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"]]}"#),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    Reaction::new(
        mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [1 5 "1"]]}"#),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    false
)]
#[case::inversion_renumbered(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(0), AtomId(2), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(0), AtomId(3), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(0), AtomId(4), EntitySpan::Unchanged(BondForm::from_order(1))),
        ],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            EntitySpan::Modified {
                lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
                rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            },
        )],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::I, Element::Br, Element::Cl, Element::F, Element::C]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        bonds: vec![
            (AtomId(4), AtomId(3), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(4), AtomId(2), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(4), AtomId(1), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(4), AtomId(0), EntitySpan::Unchanged(BondForm::from_order(1))),
        ],
        stereo_atoms: vec![(
            AtomId(4),
            vec![
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
            ],
            EntitySpan::Modified {
                lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
                rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            },
        )],
        ..Default::default()
    }),
    true
)]
#[case::inversion_reverse_reaction(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(0), AtomId(2), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(0), AtomId(3), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(0), AtomId(4), EntitySpan::Unchanged(BondForm::from_order(1))),
        ],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            EntitySpan::Modified {
                lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
                rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            },
        )],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(0), AtomId(2), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(0), AtomId(3), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(0), AtomId(4), EntitySpan::Unchanged(BondForm::from_order(1))),
        ],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            EntitySpan::Modified {
                lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
                rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
            },
        )],
        ..Default::default()
    }),
    false
)]
#[case::decalin_renumbered(
    Reaction::new(
        mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C" "C"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                        [5 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 4 "1"]]}"#
        ),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    Reaction::new(
        mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C" "C"]
                :bonds [[9 8 "1"] [8 7 "1"] [7 6 "1"] [6 5 "1"] [5 4 "1"] [4 9 "1"]
                        [4 3 "1"] [3 2 "1"] [2 1 "1"] [1 0 "1"] [0 5 "1"]]}"#
        ),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    true
)]
#[case::cuneane_renumbered(
    Reaction::new(
        mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 6 "1"]
                        [6 7 "1"] [7 0 "1"] [0 4 "1"] [1 3 "1"] [2 6 "1"] [5 7 "1"]]}"#
        ),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    Reaction::new(
        mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C"]
                :bonds [[7 6 "1"] [6 5 "1"] [5 4 "1"] [4 3 "1"] [3 2 "1"] [2 1 "1"]
                        [1 0 "1"] [0 7 "1"] [7 3 "1"] [6 4 "1"] [5 1 "1"] [2 0 "1"]]}"#
        ),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    true
)]
#[case::cubane_cuneane(
    Reaction::new(
        mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 0 "1"] [4 5 "1"] [5 6 "1"]
                        [6 7 "1"] [7 4 "1"] [0 4 "1"] [1 5 "1"] [2 6 "1"] [3 7 "1"]]}"#
        ),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    Reaction::new(
        mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 6 "1"]
                        [6 7 "1"] [7 0 "1"] [0 4 "1"] [1 3 "1"] [2 6 "1"] [5 7 "1"]]}"#
        ),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    false
)]
#[case::meso_dichlorobutane_mirror_image(
    Reaction::new(
        mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C" "Cl" "Cl"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [1 4 "1"] [2 5 "1"]]
                :stereo-atoms [{:site 1 :ligands [0 4 2 [:h 1]] :attrs "Th1"}
                               {:site 2 :ligands [3 5 1 [:h 2]] :attrs "Th0"}]}"#
        ),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    Reaction::new(
        mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C" "Cl" "Cl"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [1 4 "1"] [2 5 "1"]]
                :stereo-atoms [{:site 1 :ligands [0 4 2 [:h 1]] :attrs "Th0"}
                               {:site 2 :ligands [3 5 1 [:h 2]] :attrs "Th1"}]}"#
        ),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    true
)]
#[case::dichlorobutane_enantiomers(
    Reaction::new(
        mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C" "Cl" "Cl"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [1 4 "1"] [2 5 "1"]]
                :stereo-atoms [{:site 1 :ligands [0 4 2 [:h 1]] :attrs "Th1"}
                               {:site 2 :ligands [3 5 1 [:h 2]] :attrs "Th1"}]}"#
        ),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    Reaction::new(
        mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C" "Cl" "Cl"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [1 4 "1"] [2 5 "1"]]
                :stereo-atoms [{:site 1 :ligands [0 4 2 [:h 1]] :attrs "Th0"}
                               {:site 2 :ligands [3 5 1 [:h 2]] :attrs "Th0"}]}"#
        ),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    false
)]
#[case::charge_sum_rotated(
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 6],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
                (AtomId(3), AtomId(4), BondForm::from_order(1)),
                (AtomId(4), AtomId(5), BondForm::from_order(1)),
                (AtomId(5), AtomId(0), BondForm::from_order(1)),
            ],
            constraints: Constraint::Molecule(MoleculeConstraint::ChargeSum {
                atoms: Some(vec![AtomId(0), AtomId(1)]),
                sum: NumForm::Lit(0),
            })
            .into(),
            ..Default::default()
        }),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 6],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
                (AtomId(3), AtomId(4), BondForm::from_order(1)),
                (AtomId(4), AtomId(5), BondForm::from_order(1)),
                (AtomId(5), AtomId(0), BondForm::from_order(1)),
            ],
            constraints: Constraint::Molecule(MoleculeConstraint::ChargeSum {
                atoms: Some(vec![AtomId(3), AtomId(4)]),
                sum: NumForm::Lit(0),
            })
            .into(),
            ..Default::default()
        }),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    true
)]
#[case::charge_sum_other_atoms(
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 6],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
                (AtomId(3), AtomId(4), BondForm::from_order(1)),
                (AtomId(4), AtomId(5), BondForm::from_order(1)),
                (AtomId(5), AtomId(0), BondForm::from_order(1)),
            ],
            constraints: Constraint::Molecule(MoleculeConstraint::ChargeSum {
                atoms: Some(vec![AtomId(0), AtomId(1)]),
                sum: NumForm::Lit(0),
            })
            .into(),
            ..Default::default()
        }),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 6],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
                (AtomId(3), AtomId(4), BondForm::from_order(1)),
                (AtomId(4), AtomId(5), BondForm::from_order(1)),
                (AtomId(5), AtomId(0), BondForm::from_order(1)),
            ],
            constraints: Constraint::Molecule(MoleculeConstraint::ChargeSum {
                atoms: Some(vec![AtomId(0), AtomId(2)]),
                sum: NumForm::Lit(0),
            })
            .into(),
            ..Default::default()
        }),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    false
)]
fn test_reaction_span_canonical_eq_compounds(
    canonicalize_context: CanonicalizeContext,
    #[case] left: ReactionSpan,
    #[case] right: ReactionSpan,
    #[case] expected: bool,
) {
    assert_eq!(left.canonical_eq(&right, &canonicalize_context), expected);
    assert_eq!(right.canonical_eq(&left, &canonicalize_context), expected);
}

#[rstest]
#[case::hydrolysis_renumbered(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Removed(AtomForm::from_element(Element::F)),
            EntitySpan::Added(AtomForm::from_element(Element::O)),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Removed(BondForm::from_order(1))),
            (AtomId(0), AtomId(2), EntitySpan::Added(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Added(AtomForm::from_element(Element::O)),
            EntitySpan::Removed(AtomForm::from_element(Element::F)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
        ],
        bonds: vec![
            (AtomId(0), AtomId(2), EntitySpan::Added(BondForm::from_order(1))),
            (AtomId(1), AtomId(2), EntitySpan::Removed(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    true
)]
#[case::hydrolysis_thiolysis(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Removed(AtomForm::from_element(Element::F)),
            EntitySpan::Added(AtomForm::from_element(Element::O)),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Removed(BondForm::from_order(1))),
            (AtomId(0), AtomId(2), EntitySpan::Added(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Added(AtomForm::from_element(Element::S)),
            EntitySpan::Removed(AtomForm::from_element(Element::F)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
        ],
        bonds: vec![
            (AtomId(0), AtomId(2), EntitySpan::Added(BondForm::from_order(1))),
            (AtomId(1), AtomId(2), EntitySpan::Removed(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    false
)]
#[case::hydrolysis_reverse_reaction(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Removed(AtomForm::from_element(Element::F)),
            EntitySpan::Added(AtomForm::from_element(Element::O)),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Removed(BondForm::from_order(1))),
            (AtomId(0), AtomId(2), EntitySpan::Added(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Added(AtomForm::from_element(Element::F)),
            EntitySpan::Removed(AtomForm::from_element(Element::O)),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Added(BondForm::from_order(1))),
            (AtomId(0), AtomId(2), EntitySpan::Removed(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    false
)]
#[case::inversion_renumbered(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(0), AtomId(2), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(0), AtomId(3), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(0), AtomId(4), EntitySpan::Unchanged(BondForm::from_order(1))),
        ],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            EntitySpan::Modified {
                lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
                rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            },
        )],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::I, Element::Br, Element::Cl, Element::F, Element::C]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        bonds: vec![
            (AtomId(4), AtomId(3), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(4), AtomId(2), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(4), AtomId(1), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(4), AtomId(0), EntitySpan::Unchanged(BondForm::from_order(1))),
        ],
        stereo_atoms: vec![(
            AtomId(4),
            vec![
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
            ],
            EntitySpan::Modified {
                lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
                rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            },
        )],
        ..Default::default()
    }),
    true
)]
#[case::inversion_reverse_reaction(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(0), AtomId(2), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(0), AtomId(3), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(0), AtomId(4), EntitySpan::Unchanged(BondForm::from_order(1))),
        ],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            EntitySpan::Modified {
                lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
                rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            },
        )],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(0), AtomId(2), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(0), AtomId(3), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(0), AtomId(4), EntitySpan::Unchanged(BondForm::from_order(1))),
        ],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            EntitySpan::Modified {
                lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
                rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
            },
        )],
        ..Default::default()
    }),
    false
)]
fn test_reaction_span_canonical_eq_remapping_witness(
    canonicalize_context: CanonicalizeContext,
    #[case] left: ReactionSpan,
    #[case] right: ReactionSpan,
    #[case] expected: bool,
) {
    assert!(left.clone().canonicalize(&canonicalize_context).is_ok());
    assert!(right.clone().canonicalize(&canonicalize_context).is_ok());
    assert_eq!(
        explicitly_dense_equivalent(left.entity_counts(), |remapping| {
            left.framed_eq_under(&right, remapping)
        }),
        expected
    );
    assert_eq!(left.canonical_eq(&right, &canonicalize_context), expected);
}

#[rstest]
#[case::contradiction(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C).with_charge(0_i64)],
                ..Default::default()
            }),
            [Delta::Atom(AtomDelta::ModifyField {
                id: AtomId(0),
                change: AtomFieldChange::Charge {
                    old: NumForm::Lit(1),
                    new: NumForm::Lit(2),
                },
            })].into_iter().collect(),
        ),
        ReactionCanonicalizeError::Contradiction(Contradiction),
    )]
fn test_reaction_canonicalize_error(
    canonicalize_context: CanonicalizeContext,
    #[case] reaction: Reaction,
    #[case] expected: ReactionCanonicalizeError,
) {
    assert_eq!(reaction.canonicalize(&canonicalize_context), Err(expected));
}

#[rstest]
#[case::stereo_atom_configuration(ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        bonds: (1..=4)
            .map(|ligand| (AtomId(0), AtomId(ligand), EntitySpan::Unchanged(BondForm::from_order(1))))
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            (1..=4)
                .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                .collect(),
            EntitySpan::Unchanged(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::lit_set([]))),
        )],
        ..Default::default()
}))]
fn test_reaction_span_canonical_hash_error(
    canonicalize_context: CanonicalizeContext,
    #[case] span: ReactionSpan,
) {
    assert_eq!(
        span.canonical_hash(&canonicalize_context),
        Err(ReactionSpanCanonicalizeError::Contradiction(Contradiction)),
    );
}

#[rstest]
fn test_reaction_canonical_eq(canonicalize_context: CanonicalizeContext) {
    let source = reaction_canonicalization_fixture();
    let canonical = source
        .clone()
        .canonicalize(&canonicalize_context)
        .expect("fixed reaction canonicalizes");

    assert!(source.canonical_eq(&canonical, &canonicalize_context));
}

#[rstest]
#[case::hydrolysis_renumbered(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Removed(AtomForm::from_element(Element::F)),
            EntitySpan::Added(AtomForm::from_element(Element::O)),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Removed(BondForm::from_order(1))),
            (AtomId(0), AtomId(2), EntitySpan::Added(BondForm::from_order(1))),
        ],
        ..Default::default()
    })
    .to_reaction(),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Added(AtomForm::from_element(Element::O)),
            EntitySpan::Removed(AtomForm::from_element(Element::F)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
        ],
        bonds: vec![
            (AtomId(0), AtomId(2), EntitySpan::Added(BondForm::from_order(1))),
            (AtomId(1), AtomId(2), EntitySpan::Removed(BondForm::from_order(1))),
        ],
        ..Default::default()
    })
    .to_reaction(),
    true
)]
#[case::hydrolysis_thiolysis(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Removed(AtomForm::from_element(Element::F)),
            EntitySpan::Added(AtomForm::from_element(Element::O)),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Removed(BondForm::from_order(1))),
            (AtomId(0), AtomId(2), EntitySpan::Added(BondForm::from_order(1))),
        ],
        ..Default::default()
    })
    .to_reaction(),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Added(AtomForm::from_element(Element::S)),
            EntitySpan::Removed(AtomForm::from_element(Element::F)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
        ],
        bonds: vec![
            (AtomId(0), AtomId(2), EntitySpan::Added(BondForm::from_order(1))),
            (AtomId(1), AtomId(2), EntitySpan::Removed(BondForm::from_order(1))),
        ],
        ..Default::default()
    })
    .to_reaction(),
    false
)]
#[case::hydrolysis_reverse_reaction(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Removed(AtomForm::from_element(Element::F)),
            EntitySpan::Added(AtomForm::from_element(Element::O)),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Removed(BondForm::from_order(1))),
            (AtomId(0), AtomId(2), EntitySpan::Added(BondForm::from_order(1))),
        ],
        ..Default::default()
    })
    .to_reaction(),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Added(AtomForm::from_element(Element::F)),
            EntitySpan::Removed(AtomForm::from_element(Element::O)),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Added(BondForm::from_order(1))),
            (AtomId(0), AtomId(2), EntitySpan::Removed(BondForm::from_order(1))),
        ],
        ..Default::default()
    })
    .to_reaction(),
    false
)]
#[case::charge_set_to_the_same_value(
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C).with_charge(0_i64)],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(0),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(0),
                new: NumForm::lit_set([0]),
            },
        })]),
    ),
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C).with_charge(0_i64)],
            ..Default::default()
        }),
        Deltas::new(),
    ),
    true
)]
#[case::charge_changed(
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C).with_charge(0_i64)],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(0),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(0),
                new: NumForm::Lit(1),
            },
        })]),
    ),
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C).with_charge(0_i64)],
            ..Default::default()
        }),
        Deltas::new(),
    ),
    false
)]
fn test_reaction_canonical_eq_compounds(
    canonicalize_context: CanonicalizeContext,
    #[case] left: Reaction,
    #[case] right: Reaction,
    #[case] expected: bool,
) {
    assert_eq!(left.canonical_eq(&right, &canonicalize_context), expected);
    assert_eq!(right.canonical_eq(&left, &canonicalize_context), expected);
}

#[rstest]
#[case::one_fails(
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C).with_charge(0_i64),
            ],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(0),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        })]),
    ),
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C).with_charge(0_i64),
            ],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(0),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(0),
                new: NumForm::Lit(2),
            },
        })]),
    ),
    false
)]
#[case::both_fail_identical(
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C).with_charge(0_i64),
            ],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(0),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        })]),
    ),
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C).with_charge(0_i64),
            ],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(0),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        })]),
    ),
    true
)]
#[case::both_fail_equal_reframed_forms(
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::O).with_charge(0_i64),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
                (AtomId(3), AtomId(4), BondForm::from_order(1)),
                (AtomId(0), AtomId(4), BondForm::from_order(1)),
            ],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4)],
                AromaticSystemForm::from_electrons(vec![2, 1, 1, 1, 1]),
            )],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(0),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        })]),
    ),
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::O).with_charge(0_i64),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
                (AtomId(3), AtomId(4), BondForm::from_order(1)),
                (AtomId(0), AtomId(4), BondForm::from_order(1)),
            ],
            aromatic: vec![(
                vec![AtomId(1), AtomId(2), AtomId(3), AtomId(4), AtomId(0)],
                AromaticSystemForm::from_electrons(vec![1, 1, 1, 1, 2]),
            )],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(0),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        })]),
    ),
    true
)]
#[case::both_fail_different(
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C).with_charge(0_i64),
            ],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(0),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        })]),
    ),
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C).with_charge(0_i64),
            ],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(0),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(1),
                new: NumForm::Lit(3),
            },
        })]),
    ),
    false
)]
#[case::both_fail_equal_normal_forms(
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([1])),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemForm::from_electrons(vec![2, 1]),
            )],
            ..Default::default()
        }),
        Deltas::new(),
    ),
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C).with_charge(1_i64),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemForm::from_electrons(vec![2, 1]),
            )],
            ..Default::default()
        }),
        Deltas::new(),
    ),
    true
)]
#[case::both_fail_different_frames(
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            aromatic: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemForm::from_electrons(vec![2, 1]),
            )],
            ..Default::default()
        }),
        Deltas::new(),
    ),
    Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
            ],
            aromatic: vec![(
                vec![AtomId(1), AtomId(0), AtomId(2)],
                AromaticSystemForm::from_electrons(vec![2, 1]),
            )],
            ..Default::default()
        }),
        Deltas::new(),
    ),
    false
)]
fn test_reaction_canonical_eq_contradiction(
    canonicalize_context: CanonicalizeContext,
    #[case] left: Reaction,
    #[case] right: Reaction,
    #[case] expected: bool,
) {
    assert_eq!(left.canonical_eq(&right, &canonicalize_context), expected);
}
