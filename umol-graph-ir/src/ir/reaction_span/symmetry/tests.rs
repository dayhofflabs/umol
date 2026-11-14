//! Unit tests for the colors and the remapping checks of a reaction span.

use std::cmp::Ordering;

use rstest::{fixture, rstest};
use umol_chem::element::Element;

use super::*;
use crate::ir::molecule::symmetry::tests::constrained_color_molecule;
use crate::ir::stereo::{StereoAtomSpans, StereoBondSpans};
use crate::ir::{
    AromaticSystemForm, AromaticSystemId, AtomConstraintForm, AtomForm, BondForm, Constraint,
    ConstraintSpan, DativeBondForm, DativeBondId, Deltas, ElementForm, IncidenceLevel, Molecule,
    MoleculeConstraint, MulticenterBondForm, MulticenterBondId, NoncovalentBondForm,
    NoncovalentBondId, NoncovalentBondKind, NumForm, Reaction, ReactionSpanEntries,
    StereoAtomConstraintForm, StereoAtomForm, StereoAtomId, StereoBondConstraintForm,
    StereoBondForm, StereoBondId, StereoCoset, StereoKind, StereoLigandKind, StereoLigandPair,
    StereoLigandPosition, StereoTerm, StereogenicityForm, Topicity, TopicityForm,
    TopicityRelationForm,
};
use crate::mol_dsl;

#[fixture]
fn constrained_color_span(constrained_color_molecule: Molecule) -> ReactionSpan {
    Reaction::new(constrained_color_molecule, Deltas::default())
        .to_reaction_span()
        .unwrap()
}

#[rstest]
#[case::atom(Entity::Atom(AtomId(0)), Entity::Atom(AtomId(1)), [true, false])]
#[case::bond(Entity::Bond(BondId(0)), Entity::Bond(BondId(1)), [true, false])]
#[case::dative(Entity::DativeBond(DativeBondId(0)), Entity::DativeBond(DativeBondId(1)), [true, false])]
#[case::aromatic(Entity::AromaticSystem(AromaticSystemId(0)), Entity::AromaticSystem(AromaticSystemId(1)), [true, false])]
#[case::multicenter(Entity::MulticenterBond(MulticenterBondId(0)), Entity::MulticenterBond(MulticenterBondId(1)), [true, false])]
#[case::noncovalent(Entity::NoncovalentBond(NoncovalentBondId(0)), Entity::NoncovalentBond(NoncovalentBondId(1)), [true, false])]
#[case::stereo_atom(Entity::StereoAtom(StereoAtomId(0)), Entity::StereoAtom(StereoAtomId(1)), [true, false])]
#[case::stereo_bond(Entity::StereoBond(StereoBondId(0)), Entity::StereoBond(StereoBondId(1)), [true, false])]
fn test_initial_colors_constraints(
    constrained_color_span: ReactionSpan,
    #[case] lhs: Entity,
    #[case] rhs: Entity,
    #[case] expected: [bool; 2],
    #[values(false, true)] include_constraints: bool,
) {
    let span = constrained_color_span;
    let original = span.clone();
    let graph = span.incidence_graph(IncidenceLevel::Full);
    let colors = initial_colors(
        &span,
        &graph,
        DescriptionLevel::Structure,
        include_constraints,
    )
    .unwrap();
    assert_eq!(
        colors.entities[graph.node_of(lhs).index()] == colors.entities[graph.node_of(rhs).index()],
        expected[usize::from(include_constraints)]
    );
    assert_eq!(span, original);
}

#[rstest]
#[case::unchanged(
    [EntitySpan::Unchanged(AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(4))),
     EntitySpan::Unchanged(AtomForm::from_element(Element::C))], [true, false]
)]
#[case::removed(
    [EntitySpan::Removed(AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(4))),
     EntitySpan::Removed(AtomForm::from_element(Element::C))], [true, false]
)]
#[case::added(
    [EntitySpan::Added(AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(4))),
     EntitySpan::Added(AtomForm::from_element(Element::C))], [true, false]
)]
#[case::modified(
    [EntitySpan::Modified {
        lhs: AtomForm::from_element(Element::C),
        rhs: AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(4)),
     }, EntitySpan::Unchanged(AtomForm::from_element(Element::C))], [true, false]
)]
#[case::raw_modified_constraint(
    [EntitySpan::Modified {
        lhs: AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(4)),
        rhs: AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::Valence(NumForm::lit_set([4]))),
     }, EntitySpan::Unchanged(AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(4)))], [true, false]
)]
#[case::raw_modified_field(
    [EntitySpan::Modified {
        lhs: AtomForm::from_element(Element::C),
        rhs: AtomForm::new(ElementForm::lit_set([Element::C])),
     }, EntitySpan::Unchanged(AtomForm::from_element(Element::C))], [false, false]
)]
fn test_initial_colors_span(
    #[case] atoms: [EntitySpan<AtomForm>; 2],
    #[case] expected: [bool; 2],
    #[values(false, true)] include_constraints: bool,
) {
    let span = ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: atoms.into(),
        constraints: vec![ConstraintSpan::Unchanged(Constraint::Molecule(
            MoleculeConstraint::ChargeSum {
                atoms: None,
                sum: NumForm::lit_set([]),
            },
        ))],
        ..Default::default()
    });
    let original = span.clone();
    let graph = span.incidence_graph(IncidenceLevel::Full);
    let colors = initial_colors(
        &span,
        &graph,
        DescriptionLevel::Structure,
        include_constraints,
    )
    .unwrap();
    assert_eq!(
        colors.entities[graph.node_of(Entity::Atom(AtomId(0))).index()]
            == colors.entities[graph.node_of(Entity::Atom(AtomId(1))).index()],
        expected[usize::from(include_constraints)]
    );
    assert_eq!(span, original);
}

#[rstest]
#[case::unchanged(EntitySpan::Unchanged(StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32).with_constraint(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(7)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) }))), [Ok(true), Err(Contradiction)])]
#[case::removed(EntitySpan::Removed(StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32).with_constraint(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(7)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) }))), [Ok(false), Err(Contradiction)])]
#[case::added(EntitySpan::Added(StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32).with_constraint(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(7)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) }))), [Ok(false), Err(Contradiction)])]
#[case::modified_lhs(EntitySpan::Modified { lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32).with_constraint(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(7)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) })), rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32) }, [Ok(true), Err(Contradiction)])]
#[case::modified_rhs(EntitySpan::Modified { lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32), rhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32).with_constraint(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(7)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) })) }, [Ok(true), Err(Contradiction)])]
#[case::contradiction(EntitySpan::Unchanged(StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32).with_constraint(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::lit_set([])))), [Ok(true), Ok(false)])]
#[case::raw_configurations(EntitySpan::Modified { lhs: StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::term(StereoTerm::var_in("x", []))), rhs: StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::lit_set([0])) }, [Ok(true), Ok(true)])]
#[case::configuration_frame(EntitySpan::Modified { lhs: StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32), rhs: StereoAtomForm::new(StereoKind::Octahedral, 0_u32) }, [Err(Contradiction), Err(Contradiction)])]
fn test_initial_colors_stereo_atom(
    initial_color_span: ReactionSpan,
    #[case] attributes: EntitySpan<StereoAtomForm>,
    #[case] expected: [Result<bool, Contradiction>; 2],
    #[values(false, true)] include_constraints: bool,
) {
    let mut span = initial_color_span;
    let entries = span
        .stereo_atoms()
        .ids()
        .map(|id| {
            (
                span.stereo_atoms().site(id),
                span.stereo_atoms().ligands(id).to_vec(),
                if id == StereoAtomId(0) {
                    attributes.clone()
                } else {
                    span.stereo_atoms().attributes(id).clone()
                },
            )
        })
        .collect();
    span.stereo_atoms = StereoAtomSpans::new(entries);
    let original = span.clone();
    let graph = span.incidence_graph(IncidenceLevel::Full);
    let actual = initial_colors(
        &span,
        &graph,
        DescriptionLevel::Structure,
        include_constraints,
    )
    .map(|colors| {
        colors.entities[graph.node_of(Entity::StereoAtom(StereoAtomId(0))).index()]
            == colors.entities[graph.node_of(Entity::StereoAtom(StereoAtomId(1))).index()]
    });
    assert_eq!(actual, expected[usize::from(include_constraints)]);
    assert_eq!(span, original);
}

#[rstest]
#[case::unchanged(EntitySpan::Unchanged(StereoBondForm::new(StereoKind::CisTrans, 0_u32).with_constraint(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(7)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) }))), [Ok(true), Err(Contradiction)])]
#[case::removed(EntitySpan::Removed(StereoBondForm::new(StereoKind::CisTrans, 0_u32).with_constraint(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(7)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) }))), [Ok(false), Err(Contradiction)])]
#[case::added(EntitySpan::Added(StereoBondForm::new(StereoKind::CisTrans, 0_u32).with_constraint(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(7)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) }))), [Ok(false), Err(Contradiction)])]
#[case::modified_lhs(EntitySpan::Modified { lhs: StereoBondForm::new(StereoKind::CisTrans, 0_u32).with_constraint(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(7)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) })), rhs: StereoBondForm::new(StereoKind::CisTrans, 0_u32) }, [Ok(true), Err(Contradiction)])]
#[case::modified_rhs(EntitySpan::Modified { lhs: StereoBondForm::new(StereoKind::CisTrans, 0_u32), rhs: StereoBondForm::new(StereoKind::CisTrans, 0_u32).with_constraint(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(7)), relation: TopicityRelationForm::Lit(Topicity::Enantiotopic) })) }, [Ok(true), Err(Contradiction)])]
#[case::contradiction(EntitySpan::Unchanged(StereoBondForm::new(StereoKind::CisTrans, 0_u32).with_constraint(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::lit_set([])))), [Ok(true), Ok(false)])]
#[case::raw_configurations(EntitySpan::Modified { lhs: StereoBondForm::new(StereoKind::CisTrans, StereoCoset::term(StereoTerm::var_in("x", []))), rhs: StereoBondForm::new(StereoKind::CisTrans, StereoCoset::lit_set([0])) }, [Ok(true), Ok(true)])]
#[case::configuration_frame(EntitySpan::Modified { lhs: StereoBondForm::new(StereoKind::CisTrans, 0_u32), rhs: StereoBondForm::new(StereoKind::Tetrahedral, 0_u32) }, [Err(Contradiction), Err(Contradiction)])]
fn test_initial_colors_stereo_bond(
    initial_color_span: ReactionSpan,
    #[case] attributes: EntitySpan<StereoBondForm>,
    #[case] expected: [Result<bool, Contradiction>; 2],
    #[values(false, true)] include_constraints: bool,
) {
    let mut span = initial_color_span;
    let entries = span
        .stereo_bonds()
        .ids()
        .map(|id| {
            (
                span.stereo_bonds().site(id),
                span.stereo_bonds().ligands(id).to_vec(),
                if id == StereoBondId(0) {
                    attributes.clone()
                } else {
                    span.stereo_bonds().attributes(id).clone()
                },
            )
        })
        .collect();
    span.stereo_bonds = StereoBondSpans::new(entries);
    let original = span.clone();
    let graph = span.incidence_graph(IncidenceLevel::Full);
    let actual = initial_colors(
        &span,
        &graph,
        DescriptionLevel::Structure,
        include_constraints,
    )
    .map(|colors| {
        colors.entities[graph.node_of(Entity::StereoBond(StereoBondId(0))).index()]
            == colors.entities[graph.node_of(Entity::StereoBond(StereoBondId(1))).index()]
    });
    assert_eq!(actual, expected[usize::from(include_constraints)]);
    assert_eq!(span, original);
}

#[rstest]
#[case::rotation(
    Reaction::new(
        mol_dsl!(
            r#"{:atoms ["C" "F" "F" "F" "Cl"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
                :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
        ),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    Reaction::new(
        mol_dsl!(
            r#"{:atoms ["C" "F" "F" "F" "Cl"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
                :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
        ),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    vec![0, 2, 3, 1, 4, 6, 7, 5, 8, 9],
    true
)]
#[case::transposition(
    Reaction::new(
        mol_dsl!(
            r#"{:atoms ["C" "F" "F" "F" "Cl"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
                :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
        ),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    Reaction::new(
        mol_dsl!(
            r#"{:atoms ["C" "F" "F" "F" "Cl"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
                :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
        ),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    vec![0, 2, 1, 3, 4, 6, 5, 7, 8, 9],
    false
)]
#[case::site_moved(
    Reaction::new(
        mol_dsl!(
            r#"{:atoms ["C" "F" "Cl" "Br" "I"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
                :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
        ),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    Reaction::new(
        mol_dsl!(
            r#"{:atoms ["C" "F" "Cl" "Br" "I"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
                :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
        ),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    vec![1, 0, 2, 3, 4, 5, 6, 7, 8, 9],
    false
)]
#[case::renumbered(
    Reaction::new(
        mol_dsl!(
            r#"{:atoms ["C" "F" "Cl" "Br" "I"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
                :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
        ),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    Reaction::new(
        mol_dsl!(
            r#"{:atoms ["I" "Br" "Cl" "F" "C"]
                :bonds [[4 0 "1"] [4 1 "1"] [4 2 "1"] [4 3 "1"]]
                :stereo-atoms [{:site 4 :ligands [3 2 1 0] :attrs "Th1"}]}"#
        ),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    vec![4, 3, 2, 1, 0, 8, 7, 6, 5, 9],
    true
)]
#[case::renumbered_inverted(
    Reaction::new(
        mol_dsl!(
            r#"{:atoms ["C" "F" "Cl" "Br" "I"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
                :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
        ),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    Reaction::new(
        mol_dsl!(
            r#"{:atoms ["I" "Br" "Cl" "F" "C"]
                :bonds [[4 0 "1"] [4 1 "1"] [4 2 "1"] [4 3 "1"]]
                :stereo-atoms [{:site 4 :ligands [3 2 1 0] :attrs "Th0"}]}"#
        ),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    vec![4, 3, 2, 1, 0, 8, 7, 6, 5, 9],
    false
)]
#[case::stereo_bond(
    Reaction::new(
        mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C"]
                :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
                :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#
        ),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    Reaction::new(
        mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C"]
                :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
                :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#
        ),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    vec![0, 1, 2, 3, 4, 5, 6, 7],
    true
)]
#[case::stereo_bond_inverted(
    Reaction::new(
        mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C"]
                :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
                :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#
        ),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    Reaction::new(
        mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C"]
                :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
                :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct0"}]}"#
        ),
        Deltas::default(),
    )
    .to_reaction_span()
    .unwrap(),
    vec![0, 1, 2, 3, 4, 5, 6, 7],
    false
)]
#[case::modified_configuration(
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
    vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    true
)]
#[case::modified_configuration_reversed(
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
    vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    false
)]
fn test_remapping_preserves_stereo(
    #[case] source: ReactionSpan,
    #[case] target: ReactionSpan,
    #[case] images: Vec<u32>,
    #[case] expected: bool,
) {
    let incidence_graph = source.incidence_graph(IncidenceLevel::Full);
    let images = images.into_iter().map(NodeId).collect::<Vec<_>>();

    assert_eq!(
        remapping_preserves_stereo(&source, &target, &incidence_graph, &images),
        expected,
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
#[case::prostereogenic_explicit_hydrogens(
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
#[case::bond_endpoint(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
        ],
        bonds: vec![(
            AtomId(0),
            AtomId(1),
            EntitySpan::Unchanged(BondForm::from_order(1)),
        )],
        ..Default::default()
    }),
    Entity::Bond(BondId(0)),
    AtomId(0),
    (Incidence::BondEndpoint, None)
)]
#[case::aromatic_atom_unchanged(
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
            EntitySpan::Unchanged(AromaticSystemForm::from_electrons(vec![2, 1, 1, 1, 1])),
        )],
        ..Default::default()
    }),
    Entity::AromaticSystem(AromaticSystemId(0)),
    AtomId(0),
    (Incidence::AromaticAtom, Some(EntitySpan::Unchanged(Some(2))))
)]
#[case::aromatic_atom_modified(
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
    Entity::AromaticSystem(AromaticSystemId(0)),
    AtomId(0),
    (
        Incidence::AromaticAtom,
        Some(EntitySpan::Modified {
            lhs: Some(2),
            rhs: Some(1),
        })
    )
)]
#[case::aromatic_atom_equal_sides(
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
    Entity::AromaticSystem(AromaticSystemId(0)),
    AtomId(1),
    (Incidence::AromaticAtom, Some(EntitySpan::Unchanged(Some(1))))
)]
#[case::aromatic_atom_added(
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
            EntitySpan::Added(AromaticSystemForm::from_electrons(vec![2, 1, 1, 1, 1])),
        )],
        ..Default::default()
    }),
    Entity::AromaticSystem(AromaticSystemId(0)),
    AtomId(0),
    (Incidence::AromaticAtom, Some(EntitySpan::Added(Some(2))))
)]
#[case::aromatic_atom_removed(
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
            EntitySpan::Removed(AromaticSystemForm::from_electrons(vec![2, 1, 1, 1, 1])),
        )],
        ..Default::default()
    }),
    Entity::AromaticSystem(AromaticSystemId(0)),
    AtomId(0),
    (Incidence::AromaticAtom, Some(EntitySpan::Removed(Some(2))))
)]
fn test_incidence_color(
    #[case] span: ReactionSpan,
    #[case] owner: Entity,
    #[case] atom: AtomId,
    #[case] expected: (Incidence, Option<EntitySpan<Option<i64>>>),
) {
    let incidence_graph = span.incidence_graph(IncidenceLevel::Constitution);
    let edge = incidence_graph
        .graph()
        .find_edge(
            incidence_graph.node_of(owner),
            incidence_graph.node_of(Entity::Atom(atom)),
        )
        .unwrap();
    let (incidence, electron_count) = incidence_color(&span, &incidence_graph, edge);

    assert_eq!((incidence.clone(), electron_count), expected);
}

#[fixture]
pub(crate) fn initial_color_span() -> ReactionSpan {
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
    let colors = initial_colors(&span, &incidence_graph, DescriptionLevel::Structure, true)
        .expect("fixed span has colors");

    assert_eq!(
        colors.entities[incidence_graph.node_of(lhs).index()]
            == colors.entities[incidence_graph.node_of(rhs).index()],
        expected
    );
}

#[rstest]
#[case::removed(Entity::Atom(AtomId(0)), Entity::Atom(AtomId(3)), Ordering::Less)]
#[case::modified(Entity::Atom(AtomId(3)), Entity::Atom(AtomId(4)), Ordering::Less)]
#[case::added(Entity::Atom(AtomId(4)), Entity::Atom(AtomId(2)), Ordering::Less)]
#[case::modified_sides(Entity::Atom(AtomId(4)), Entity::Atom(AtomId(6)), Ordering::Less)]
#[case::aromatic_frame(
    Entity::AromaticSystem(AromaticSystemId(0)),
    Entity::AromaticSystem(AromaticSystemId(1)),
    Ordering::Equal
)]
#[case::aromatic_field(
    Entity::AromaticSystem(AromaticSystemId(1)),
    Entity::AromaticSystem(AromaticSystemId(2)),
    Ordering::Less
)]
#[case::stereo_atom_frame(
    Entity::StereoAtom(StereoAtomId(0)),
    Entity::StereoAtom(StereoAtomId(1)),
    Ordering::Equal
)]
#[case::stereo_bond_frame(
    Entity::StereoBond(StereoBondId(0)),
    Entity::StereoBond(StereoBondId(1)),
    Ordering::Equal
)]
fn test_initial_colors_order(
    initial_color_span: ReactionSpan,
    #[case] lhs: Entity,
    #[case] rhs: Entity,
    #[case] expected: Ordering,
    #[values(false, true)] include_constraints: bool,
) {
    let span = initial_color_span;
    let incidence_graph = span.incidence_graph(IncidenceLevel::Full);
    let colors = initial_colors(
        &span,
        &incidence_graph,
        DescriptionLevel::Structure,
        include_constraints,
    )
    .unwrap();
    assert_eq!(
        colors.entities[incidence_graph.node_of(lhs).index()]
            .cmp(&colors.entities[incidence_graph.node_of(rhs).index()]),
        expected,
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
    #[values(false, true)] include_constraints: bool,
) {
    let span = initial_color_span
        .normalize()
        .expect("fixed span normalizes");
    let incidence_graph = span.incidence_graph(IncidenceLevel::Full);
    let colors = initial_colors(
        &span,
        &incidence_graph,
        DescriptionLevel::Structure,
        include_constraints,
    )
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
fn test_initial_colors_error(
    #[case] span: ReactionSpan,
    #[case] level: DescriptionLevel,
    #[values(false, true)] include_constraints: bool,
) {
    let incidence_graph = span.incidence_graph(match level {
        DescriptionLevel::Topology => IncidenceLevel::Topology,
        DescriptionLevel::Constitution => IncidenceLevel::Constitution,
        DescriptionLevel::Structure | DescriptionLevel::Full => IncidenceLevel::Full,
    });

    assert_eq!(
        initial_colors(&span, &incidence_graph, level, include_constraints),
        Err(Contradiction)
    );
}
