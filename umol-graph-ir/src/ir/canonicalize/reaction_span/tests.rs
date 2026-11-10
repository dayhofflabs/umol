//! Unit tests for reaction-span and reaction canonicalization.

use std::iter;

use rstest::rstest;
use umol_chem::element::Element;

use super::super::tests::{canonicalize_context, initial_color_molecule};
use super::super::*;
use super::*;
use crate::ir::{
    AromaticSystemForm, AromaticSystemId, AtomFieldChange, AtomForm, AtomId, BondId,
    ConstraintDelta, DativeBondForm, DativeBondId, Deltas, Entity, IncidenceLevel,
    MoleculeCorrespondence, MoleculeEntries, MulticenterBondForm, MulticenterBondId,
    NoncovalentBondForm, NoncovalentBondId, ReactionSpanEntries, StereoAtomForm, StereoAtomId,
    StereoBondForm, StereoBondId, StereoKind, StereoLigand, StereoLigandKind,
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
#[case::unchanged(
    EntitySpan::Unchanged(StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32)),
    vec![1, 0, 2, 3],
    SpanTagPosition::UNCHANGED,
    vec![0, 1, 2, 3],
    vec![StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1u32)],
)]
#[case::added(
    EntitySpan::Added(StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32)),
    vec![1, 0, 2, 3],
    SpanTagPosition::ADDED,
    vec![0, 1, 2, 3],
    vec![StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1u32)],
)]
#[case::removed(
    EntitySpan::Removed(StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32)),
    vec![1, 0, 2, 3],
    SpanTagPosition::REMOVED,
    vec![0, 1, 2, 3],
    vec![StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1u32)],
)]
#[case::modified(
    EntitySpan::Modified {
        lhs: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32),
        rhs: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1u32),
    },
    vec![1, 0, 2, 3],
    SpanTagPosition::MODIFIED,
    vec![0, 1, 2, 3],
    vec![
        StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1u32),
        StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32),
    ],
)]
#[case::modified_equal_sides(
    EntitySpan::Modified {
        lhs: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32),
        rhs: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32),
    },
    vec![1, 0, 2, 3],
    SpanTagPosition::UNCHANGED,
    vec![0, 1, 2, 3],
    vec![StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1u32)],
)]
#[case::modified_equal_ligand_classes(
    EntitySpan::Modified {
        lhs: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1u32),
        rhs: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32),
    },
    vec![0, 0, 1, 2],
    SpanTagPosition::MODIFIED,
    vec![0, 0, 1, 2],
    vec![
        StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32),
        StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1u32),
    ],
)]
#[case::modified_lhs_undetermined(
    EntitySpan::Modified {
        lhs: StereoConfigurationForm::Undetermined,
        rhs: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32),
    },
    vec![1, 0, 2, 3],
    SpanTagPosition::MODIFIED,
    vec![0, 1, 2, 3],
    vec![
        StereoConfigurationForm::Undetermined,
        StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1u32),
    ],
)]
#[case::modified_rhs_undetermined(
    EntitySpan::Modified {
        lhs: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32),
        rhs: StereoConfigurationForm::Undetermined,
    },
    vec![1, 0, 2, 3],
    SpanTagPosition::MODIFIED,
    vec![0, 1, 2, 3],
    vec![
        StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1u32),
        StereoConfigurationForm::Undetermined,
    ],
)]
#[case::undetermined(
    EntitySpan::Unchanged(StereoConfigurationForm::Undetermined),
    vec![2, 0, 1],
    SpanTagPosition::UNCHANGED,
    vec![0, 1, 2],
    vec![StereoConfigurationForm::Undetermined],
)]
#[case::modified_undetermined_sides(
    EntitySpan::Modified {
        lhs: StereoConfigurationForm::Undetermined,
        rhs: StereoConfigurationForm::Undetermined,
    },
    vec![2, 0, 1],
    SpanTagPosition::UNCHANGED,
    vec![0, 1, 2],
    vec![StereoConfigurationForm::Undetermined],
)]
fn test_stereo_refinement_descriptor(
    #[case] span: EntitySpan<StereoConfigurationForm>,
    #[case] ligand_classes: Vec<u32>,
    #[case] expected_position: SpanTagPosition,
    #[case] expected_ligand_classes: Vec<u32>,
    #[case] expected_configurations: Vec<StereoConfigurationForm>,
) {
    let ligands = ligand_classes
        .into_iter()
        .map(|class| (class, StereoLigandKind::Atom))
        .collect::<Vec<_>>();
    let expected = CanonicalKeyValue::Span(SpanKey {
        position: expected_position,
        values: expected_configurations
            .iter()
            .map(|configuration| {
                product([
                    CanonicalKeyValue::Unsigned(7),
                    sequence(
                        expected_ligand_classes
                            .iter()
                            .map(|&class| stereo_ligand_key(class, StereoLigandKind::Atom)),
                    ),
                    stereo_configuration_form_key(configuration),
                ])
            })
            .collect(),
    });

    assert_eq!(
        stereo_refinement_descriptor(7, &ligands, &span, |configuration| configuration),
        Ok(expected),
    );
}

#[rstest]
#[case::unchanged(|kind| EntitySpan::Unchanged(StereoConfigurationForm::kinded(kind, 0u32)))]
#[case::added(|kind| EntitySpan::Added(StereoConfigurationForm::kinded(kind, 0u32)))]
#[case::removed(|kind| EntitySpan::Removed(StereoConfigurationForm::kinded(kind, 0u32)))]
#[case::modified(|kind| EntitySpan::Modified {
    lhs: StereoConfigurationForm::kinded(kind, 0u32),
    rhs: StereoConfigurationForm::kinded(kind, 1u32),
})]
#[case::modified_lhs_undetermined(|kind| EntitySpan::Modified {
    lhs: StereoConfigurationForm::Undetermined,
    rhs: StereoConfigurationForm::kinded(kind, 0u32),
})]
#[case::modified_rhs_undetermined(|kind| EntitySpan::Modified {
    lhs: StereoConfigurationForm::kinded(kind, 0u32),
    rhs: StereoConfigurationForm::Undetermined,
})]
fn test_stereo_refinement_descriptor_frame_invariant(
    #[values(
        StereoKind::Tetrahedral,
        StereoKind::CisTrans,
        StereoKind::Axial,
        StereoKind::SquarePlanar,
        StereoKind::TrigonalBipyramidal,
        StereoKind::Octahedral
    )]
    kind: StereoKind,
    #[case] span: fn(StereoKind) -> EntitySpan<StereoConfigurationForm>,
) {
    let degree = kind.degree();
    let frame = if matches!(kind, StereoKind::CisTrans | StereoKind::Axial) {
        Permutation::from_image(&[2, 3, 0, 1])
    } else {
        Permutation::from_image(&(1..degree).chain(iter::once(0)).collect::<Vec<_>>())
    };
    let ligands = (0..degree)
        .map(|class| (class as u32, StereoLigandKind::Atom))
        .collect::<Vec<_>>();
    let span = span(kind);
    let transported = span
        .clone()
        .try_map(|configuration| configuration.apply(frame))
        .expect("the permutation is a parent-group action of the form's kind");

    assert_eq!(
        stereo_refinement_descriptor(7, &frame.act(&ligands), &transported, |configuration| {
            configuration
        })
        .expect("the frame has the degree of the kind"),
        stereo_refinement_descriptor(7, &ligands, &span, |configuration| configuration)
            .expect("the frame has the degree of the kind"),
    );
}

#[rstest]
#[case::unchanged_short(
    EntitySpan::Unchanged(StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32)),
    3
)]
#[case::unchanged_long(
    EntitySpan::Unchanged(StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32)),
    5
)]
#[case::added_short(
    EntitySpan::Added(StereoConfigurationForm::kinded(StereoKind::CisTrans, 0u32)),
    3
)]
#[case::removed_long(
    EntitySpan::Removed(StereoConfigurationForm::kinded(StereoKind::CisTrans, 0u32)),
    5
)]
#[case::modified_short(
    EntitySpan::Modified {
        lhs: StereoConfigurationForm::kinded(StereoKind::Octahedral, 0u32),
        rhs: StereoConfigurationForm::kinded(StereoKind::Octahedral, 1u32),
    },
    5
)]
#[case::modified_lhs_undetermined_long(
    EntitySpan::Modified {
        lhs: StereoConfigurationForm::Undetermined,
        rhs: StereoConfigurationForm::kinded(StereoKind::Octahedral, 0u32),
    },
    7
)]
fn test_stereo_refinement_descriptor_error(
    #[case] span: EntitySpan<StereoConfigurationForm>,
    #[case] degree: usize,
) {
    let ligands = (0..degree)
        .map(|class| (class as u32, StereoLigandKind::Atom))
        .collect::<Vec<_>>();

    assert_eq!(
        stereo_refinement_descriptor(0, &ligands, &span, |configuration| configuration),
        Err(Contradiction),
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
    [6, 3],
)]
fn test_canonical_candidate(
    canonicalize_context: CanonicalizeContext,
    #[case] span: ReactionSpan,
    #[case] level: DescriptionLevel,
    #[case] expected_visited_leaves: [usize; 2],
) {
    let unpruned = canonical_candidate(
        &span,
        level,
        &canonicalize_context,
        CanonicalSearchOptions {
            automorphism_pruning: false,
            branch_order: BranchOrdering::NodeId,
        },
    )
    .expect("fixed span has initial colors");
    let pruned = canonical_candidate(
        &span,
        level,
        &canonicalize_context,
        CanonicalSearchOptions {
            automorphism_pruning: true,
            branch_order: BranchOrdering::BackendCanonical,
        },
    )
    .expect("fixed span has initial colors");

    assert_eq!(pruned.candidate.key, unpruned.candidate.key);
    assert_eq!(
        [unpruned.stats.visited_leaves, pruned.stats.visited_leaves],
        expected_visited_leaves,
    );
}

#[rstest]
#[case::atom(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Added(AtomForm::from_element(Element::O)),
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Removed(AtomForm::from_element(Element::N)),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: vec![
            EntitySpan::Unchanged(AtomForm::from_element(Element::C)),
            EntitySpan::Removed(AtomForm::from_element(Element::N)),
            EntitySpan::Added(AtomForm::from_element(Element::O)),
        ],
        ..Default::default()
    }),
    &[
        DescriptionLevel::Topology,
        DescriptionLevel::Constitution,
        DescriptionLevel::Structure,
        DescriptionLevel::Full,
    ],
)]
#[case::bond(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::N, Element::O]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        bonds: vec![
            (AtomId(0), AtomId(1), EntitySpan::Added(BondForm::from_order(1))),
            (AtomId(1), AtomId(2), EntitySpan::Unchanged(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::N, Element::O]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        bonds: vec![
            (AtomId(1), AtomId(2), EntitySpan::Unchanged(BondForm::from_order(1))),
            (AtomId(0), AtomId(1), EntitySpan::Added(BondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    &[
        DescriptionLevel::Topology,
        DescriptionLevel::Constitution,
        DescriptionLevel::Structure,
        DescriptionLevel::Full,
    ],
)]
#[case::dative_bond(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::B, Element::N, Element::Al]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        dative: vec![
            (vec![AtomId(1)], AtomId(0), EntitySpan::Added(DativeBondForm::from_order(1))),
            (vec![AtomId(1)], AtomId(2), EntitySpan::Unchanged(DativeBondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::B, Element::N, Element::Al]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        dative: vec![
            (vec![AtomId(1)], AtomId(2), EntitySpan::Unchanged(DativeBondForm::from_order(1))),
            (vec![AtomId(1)], AtomId(0), EntitySpan::Added(DativeBondForm::from_order(1))),
        ],
        ..Default::default()
    }),
    &[DescriptionLevel::Constitution, DescriptionLevel::Structure, DescriptionLevel::Full],
)]
#[case::aromatic_system(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::N, Element::O, Element::S]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        aromatic: vec![
            (vec![AtomId(0), AtomId(1)], EntitySpan::Added(AromaticSystemForm::default())),
            (vec![AtomId(2), AtomId(3)], EntitySpan::Unchanged(AromaticSystemForm::default())),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::N, Element::O, Element::S]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        aromatic: vec![
            (vec![AtomId(2), AtomId(3)], EntitySpan::Unchanged(AromaticSystemForm::default())),
            (vec![AtomId(0), AtomId(1)], EntitySpan::Added(AromaticSystemForm::default())),
        ],
        ..Default::default()
    }),
    &[DescriptionLevel::Constitution, DescriptionLevel::Structure, DescriptionLevel::Full],
)]
#[case::multicenter_bond(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::N, Element::O, Element::S]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        multicenter: vec![
            (vec![AtomId(0), AtomId(1)], EntitySpan::Added(MulticenterBondForm::default())),
            (vec![AtomId(2), AtomId(3)], EntitySpan::Unchanged(MulticenterBondForm::default())),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::N, Element::O, Element::S]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        multicenter: vec![
            (vec![AtomId(2), AtomId(3)], EntitySpan::Unchanged(MulticenterBondForm::default())),
            (vec![AtomId(0), AtomId(1)], EntitySpan::Added(MulticenterBondForm::default())),
        ],
        ..Default::default()
    }),
    &[DescriptionLevel::Constitution, DescriptionLevel::Structure, DescriptionLevel::Full],
)]
#[case::noncovalent_bond(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::N, Element::O, Element::S]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        noncovalent: vec![
            ([AtomId(0), AtomId(1)], EntitySpan::Added(NoncovalentBondForm::default())),
            ([AtomId(2), AtomId(3)], EntitySpan::Unchanged(NoncovalentBondForm::default())),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [Element::C, Element::N, Element::O, Element::S]
            .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
            .into(),
        noncovalent: vec![
            ([AtomId(2), AtomId(3)], EntitySpan::Unchanged(NoncovalentBondForm::default())),
            ([AtomId(0), AtomId(1)], EntitySpan::Added(NoncovalentBondForm::default())),
        ],
        ..Default::default()
    }),
    &[DescriptionLevel::Constitution, DescriptionLevel::Structure, DescriptionLevel::Full],
)]
#[case::stereo_atom(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [
            Element::C, Element::N, Element::O, Element::F,
            Element::Si, Element::P, Element::S, Element::Cl,
        ]
        .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
        .into(),
        bonds: [(0, 1), (0, 2), (0, 3), (4, 5), (4, 6), (4, 7)]
            .map(|(first, second)| {
                (AtomId(first), AtomId(second), EntitySpan::Unchanged(BondForm::from_order(1)))
            })
            .into(),
        stereo_atoms: vec![
            (
                AtomId(4),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(6), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                ],
                EntitySpan::Added(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32)),
            ),
            (
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                ],
                EntitySpan::Unchanged(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32)),
            ),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [
            Element::C, Element::N, Element::O, Element::F,
            Element::Si, Element::P, Element::S, Element::Cl,
        ]
        .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
        .into(),
        bonds: [(0, 1), (0, 2), (0, 3), (4, 5), (4, 6), (4, 7)]
            .map(|(first, second)| {
                (AtomId(first), AtomId(second), EntitySpan::Unchanged(BondForm::from_order(1)))
            })
            .into(),
        stereo_atoms: vec![
            (
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                ],
                EntitySpan::Unchanged(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32)),
            ),
            (
                AtomId(4),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(6), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                ],
                EntitySpan::Added(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32)),
            ),
        ],
        ..Default::default()
    }),
    &[DescriptionLevel::Structure, DescriptionLevel::Full],
)]
#[case::stereo_bond(
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [
            Element::C, Element::N, Element::O, Element::F,
            Element::Si, Element::P, Element::S, Element::Cl,
        ]
        .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
        .into(),
        bonds: [(0, 2, 1), (1, 3, 1), (4, 6, 1), (5, 7, 1), (0, 1, 2), (4, 5, 2)]
            .map(|(first, second, order)| {
                (AtomId(first), AtomId(second), EntitySpan::Unchanged(BondForm::from_order(order)))
            })
            .into(),
        stereo_bonds: vec![
            (
                BondId(5),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(6), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                ],
                EntitySpan::Added(StereoBondForm::new(StereoKind::CisTrans, 0u32)),
            ),
            (
                BondId(4),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                ],
                EntitySpan::Unchanged(StereoBondForm::new(StereoKind::CisTrans, 0u32)),
            ),
        ],
        ..Default::default()
    }),
    ReactionSpan::from_entries(ReactionSpanEntries {
        atoms: [
            Element::C, Element::N, Element::O, Element::F,
            Element::Si, Element::P, Element::S, Element::Cl,
        ]
        .map(|element| EntitySpan::Unchanged(AtomForm::from_element(element)))
        .into(),
        bonds: [(0, 2, 1), (1, 3, 1), (4, 6, 1), (5, 7, 1), (0, 1, 2), (4, 5, 2)]
            .map(|(first, second, order)| {
                (AtomId(first), AtomId(second), EntitySpan::Unchanged(BondForm::from_order(order)))
            })
            .into(),
        stereo_bonds: vec![
            (
                BondId(4),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                ],
                EntitySpan::Unchanged(StereoBondForm::new(StereoKind::CisTrans, 0u32)),
            ),
            (
                BondId(5),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(6), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                ],
                EntitySpan::Added(StereoBondForm::new(StereoKind::CisTrans, 0u32)),
            ),
        ],
        ..Default::default()
    }),
    &[DescriptionLevel::Structure, DescriptionLevel::Full],
)]
fn test_reaction_span_tracked_canonicalize_at(
    canonicalize_context: CanonicalizeContext,
    #[case] source: ReactionSpan,
    #[case] expected: ReactionSpan,
    #[case] levels: &[DescriptionLevel],
) {
    for &level in levels {
        assert_eq!(
            source
                .clone()
                .tracked_canonicalize_at(level, &canonicalize_context)
                .map(|(canonical, _)| canonical),
            Ok(expected.clone()),
            "{level:?}",
        );
    }
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
