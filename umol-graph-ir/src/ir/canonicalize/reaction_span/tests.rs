//! Unit tests for reaction-span and reaction canonicalization.

use rstest::rstest;
use umol_chem::element::Element;

use super::super::tests::canonicalize_context;
use super::super::*;
use super::*;
use crate::ir::{
    AromaticSystemForm, AromaticSystemId, AtomFieldChange, AtomForm, AtomId, ConstraintDelta,
    Entity, IncidenceLevel, MoleculeEntries, MulticenterBondForm, MulticenterBondId,
    NoncovalentBondForm, ReactionSpanEntries,
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
fn test_reaction_canonicalize_description_level(
    canonicalize_context: CanonicalizeContext,
    #[case] source: Reaction,
    #[case] expected: DescriptionLevel,
) {
    let forced_full = source
        .clone()
        .tracked_canonicalize_at(DescriptionLevel::Full, &canonicalize_context)
        .map(|(canonical, _)| canonical)
        .expect("representative reaction canonicalizes at Full");
    let canonical = source
        .clone()
        .canonicalize(&canonicalize_context)
        .expect("representative reaction canonicalizes by its effective level");
    let (with_remapping, remapping) = source
        .clone()
        .tracked_canonicalize(&canonicalize_context)
        .expect("representative reaction canonicalizes with a remapping");
    let transported = source
        .to_reaction_span()
        .expect("representative reaction materializes")
        .remap(&remapping)
        .reframe()
        .expect("canonical remapping preserves compatible frames")
        .to_reaction();

    assert_eq!(reaction_canonicalize_level(&source), expected);
    assert_eq!(canonical, forced_full);
    assert_eq!(with_remapping, canonical);
    assert_eq!(transported, canonical);
    assert_eq!(
        source.canonical_hash(&canonicalize_context),
        canonical.canonical_hash(&canonicalize_context),
    );
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
fn test_reaction_span_canonicalize_description_level(
    canonicalize_context: CanonicalizeContext,
    #[case] reaction: Reaction,
    #[case] expected: DescriptionLevel,
) {
    let source = reaction
        .to_reaction_span()
        .expect("representative reaction materializes");
    let forced_full = source
        .clone()
        .tracked_canonicalize_at(DescriptionLevel::Full, &canonicalize_context)
        .map(|(canonical, _)| canonical)
        .expect("representative span canonicalizes at Full");
    let canonical = source
        .clone()
        .canonicalize(&canonicalize_context)
        .expect("representative span canonicalizes by its effective level");
    let (with_remapping, remapping) = source
        .clone()
        .tracked_canonicalize(&canonicalize_context)
        .expect("representative span canonicalizes with a remapping");
    let images = reaction_span_counts(&source).map(|count| (0..count).rev().collect::<Vec<_>>());
    let remapped = source.clone().remap(&molecule_remapping(&images));

    assert_eq!(reaction_span_canonicalize_level(&source), expected);
    assert_eq!(reaction_span_canonicalize_level(&remapped), expected);
    assert_eq!(canonical, forced_full);
    assert_eq!(with_remapping, canonical);
    assert_eq!(
        source.clone().remap(&remapping).reframe(),
        Ok(canonical.clone())
    );
    assert_eq!(
        source.canonical_hash(&canonicalize_context),
        canonical.canonical_hash(&canonicalize_context),
    );
}

#[rstest]
#[case::unchanged(
    EntitySpan::Unchanged(AromaticSystemForm::from_electrons(vec![1, 2])),
    EntitySpan::Unchanged(MulticenterBondForm::from_electrons(vec![1, 2])),
    [(SpanTagPosition::UNCHANGED, vec![NumForm::Lit(1)]),
     (SpanTagPosition::UNCHANGED, vec![NumForm::Lit(2)])],
)]
#[case::added(
    EntitySpan::Added(AromaticSystemForm::from_electrons(vec![1, 2])),
    EntitySpan::Added(MulticenterBondForm::from_electrons(vec![1, 2])),
    [(SpanTagPosition::ADDED, vec![NumForm::Lit(1)]),
     (SpanTagPosition::ADDED, vec![NumForm::Lit(2)])],
)]
#[case::removed(
    EntitySpan::Removed(AromaticSystemForm::from_electrons(vec![1, 2])),
    EntitySpan::Removed(MulticenterBondForm::from_electrons(vec![1, 2])),
    [(SpanTagPosition::REMOVED, vec![NumForm::Lit(1)]),
     (SpanTagPosition::REMOVED, vec![NumForm::Lit(2)])],
)]
#[case::modified(
    EntitySpan::Modified {
        lhs: AromaticSystemForm::from_electrons(vec![1, 2]),
        rhs: AromaticSystemForm::from_electrons(vec![1, 0]),
    },
    EntitySpan::Modified {
        lhs: MulticenterBondForm::from_electrons(vec![1, 2]),
        rhs: MulticenterBondForm::from_electrons(vec![1, 0]),
    },
    [(SpanTagPosition::UNCHANGED, vec![NumForm::Lit(1)]),
     (SpanTagPosition::MODIFIED, vec![NumForm::Lit(2), NumForm::Lit(0)])],
)]
#[case::undetermined(
    EntitySpan::Unchanged(AromaticSystemForm::default()),
    EntitySpan::Unchanged(MulticenterBondForm::default()),
    [(SpanTagPosition::UNCHANGED, vec![NumForm::Undetermined]),
     (SpanTagPosition::UNCHANGED, vec![NumForm::Undetermined])],
)]
fn test_initial_color_keys_incidence(
    #[case] aromatic: EntitySpan<AromaticSystemForm>,
    #[case] multicenter: EntitySpan<MulticenterBondForm>,
    #[case] expected: [(SpanTagPosition, Vec<NumForm>); 2],
) {
    let span = ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![EntitySpan::Unchanged(AtomForm::from_element(Element::C)); 2],
        aromatic: vec![(vec![AtomId(1), AtomId(0)], aromatic)],
        multicenter: vec![(vec![AtomId(0), AtomId(1)], multicenter)],
        ..Default::default()
    });
    let incidence = span.incidence_graph(IncidenceLevel::Constitution);
    let (_, keys) = initial_color_keys(&span, &incidence).unwrap();
    assert_eq!(
        incidence
            .incidences()
            .map(|(_, role)| role.clone())
            .collect::<Vec<_>>(),
        vec![
            Incidence::AromaticAtom,
            Incidence::AromaticAtom,
            Incidence::MulticenterAtom,
            Incidence::MulticenterAtom
        ],
    );
    let expected = expected.map(|(position, values)| {
        CanonicalKeyValue::Span(SpanKey {
            position,
            values: values.iter().map(num_form_key).collect(),
        })
    });
    assert_eq!(
        incidence
            .incidences()
            .map(|(edge, _)| {
                let [atom, owner] = incidence.graph().edge_endpoints(edge);
                (
                    incidence.entity(atom),
                    incidence.entity(owner),
                    keys[edge.index()].clone(),
                )
            })
            .collect::<Vec<_>>(),
        vec![
            (
                Entity::Atom(AtomId(1)),
                Entity::AromaticSystem(AromaticSystemId(0)),
                InitialColorKey::Incidence(variant(3, [expected[0].clone()]))
            ),
            (
                Entity::Atom(AtomId(0)),
                Entity::AromaticSystem(AromaticSystemId(0)),
                InitialColorKey::Incidence(variant(3, [expected[1].clone()]))
            ),
            (
                Entity::Atom(AtomId(0)),
                Entity::MulticenterBond(MulticenterBondId(0)),
                InitialColorKey::Incidence(variant(4, [expected[0].clone()]))
            ),
            (
                Entity::Atom(AtomId(1)),
                Entity::MulticenterBond(MulticenterBondId(0)),
                InitialColorKey::Incidence(variant(4, [expected[1].clone()]))
            ),
        ],
    );
}

#[rstest]
#[case::topology(DescriptionLevel::Topology)]
#[case::constitution(DescriptionLevel::Constitution)]
#[case::structure(DescriptionLevel::Structure)]
#[case::full(DescriptionLevel::Full)]
fn test_canonicalize_checked_reaction_span_by(
    canonicalize_context: CanonicalizeContext,
    #[case] level: DescriptionLevel,
) {
    let source = ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Added(AtomForm::from_element(Element::O)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Removed(AtomForm::from_element(Element::N)),
        ],
        bonds: vec![(
            AtomId(1),
            AtomId(2),
            EntitySpan::Removed(BondForm::from_order(1)),
        )],
        ..Default::default()
    });
    let expected = ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Removed(AtomForm::from_element(Element::N)),
            EntitySpan::Added(AtomForm::from_element(Element::O)),
        ],
        bonds: vec![(
            AtomId(0),
            AtomId(1),
            EntitySpan::Removed(BondForm::from_order(1)),
        )],
        ..Default::default()
    });

    assert_eq!(
        source
            .clone()
            .tracked_canonicalize_at(level, &canonicalize_context)
            .map(|(canonical, _)| canonical),
        Ok(expected.clone()),
    );
    assert_eq!(
        source
            .clone()
            .tracked_canonicalize_at(level, &canonicalize_context)
            .map(|(canonical, _)| canonical),
        Ok(expected),
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
        Ok(expected.clone()),
    );
    assert_eq!(
        source.canonical_key_at(level, &canonicalize_context),
        expected.canonical_key_at(level, &canonicalize_context),
    );
}
