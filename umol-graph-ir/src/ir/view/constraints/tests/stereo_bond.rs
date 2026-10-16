use pretty_assertions::assert_eq;

use super::*;

#[rstest]
#[case::present(
    Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
        Stereogenicity::Stereogenic,
    ))),
    Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
        Stereogenicity::Stereogenic,
    )))
)]
#[case::absent(
    None,
    Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined))
)]
fn test_stereo_bond_constraints_view_asserted_complete(
    #[case] set: Option<StereoBondConstraintForm>,
    #[case] expected: Option<StereoBondConstraintForm>,
) {
    let mut form = StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1));
    if let Some(constraint) = set {
        form.constraints.set(constraint);
    }
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 4],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(2)),
            (AtomId(2), AtomId(3), BondForm::from_order(1)),
        ],
        stereo_bonds: vec![(
            BondId(1),
            vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
            ],
            form,
        )],
        ..Default::default()
    });
    assert_eq!(
        molecule
            .stereo_bond(StereoBondId(0))
            .constraints()
            .asserted_complete(StereoBondConstraintKey::Stereogenicity),
        expected
    );
}

#[rstest]
#[case::symmetry(
    StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm {
        permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper },
        invariant: BooleanForm::Lit(true),
    }),
    Some(StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm {
        permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper },
        invariant: BooleanForm::Lit(true),
    })),
)]
#[case::fluxionality(
    StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::identity(4)), active: BooleanForm::Lit(false) }),
    Some(StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::identity(4)), active: BooleanForm::Lit(false) })),
)]
#[case::topicity(
    StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::lit_set([Topicity::Homotopic]) }),
    Some(StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
)]
#[case::stereogenicity(
    StereoBondConstraintForm::Stereogenicity(StereogenicityForm::lit_set([Stereogenicity::Stereogenic])),
    Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
)]
#[case::topicity_bottom(
    StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::lit_set([]) }),
    None,
)]
#[case::stereogenicity_bottom(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::lit_set([])), None)]
fn test_stereo_bond_constraints_view_asserted_constraints(
    #[case] raw: StereoBondConstraintForm,
    #[case] normalized: Option<StereoBondConstraintForm>,
    #[values(false, true)] kindless: bool,
) {
    let form = StereoBondForm {
        configuration: if kindless {
            StereoConfigurationForm::Undetermined
        } else {
            StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Undetermined)
        },
        constraints: raw.clone().into(),
    };
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 6],
        bonds: [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)]
            .into_iter()
            .map(|(a, b)| (AtomId(a), AtomId(b), BondForm::from_order(1)))
            .collect(),
        stereo_bonds: vec![(
            BondId(0),
            [2, 3, 4, 5]
                .into_iter()
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                .collect(),
            form,
        )],
        ..Default::default()
    });
    let entity = molecule.stereo_bond(StereoBondId(0));
    let view = entity.constraints();
    let key = raw.key();
    let top = Some(StereoBondConstraintForm::undetermined(key));
    let expected = if kindless {
        top.clone()
    } else {
        normalized.clone()
    };
    assert_eq!(view.asserted(key), expected);
    assert_eq!(view.asserted_complete(key), expected);
    assert_eq!(view.evidenced(key), expected);
    assert_eq!(view.evidenced_complete(key), expected);
    assert_eq!(view.derived(key), top);
    assert_eq!(view.derived_complete(key), top);
    assert_eq!(entity.attributes().constraints.asserted(key), normalized);
    assert_eq!(
        entity.attributes().constraints.asserted_complete(key),
        normalized
    );
    assert_eq!(view.iter().collect::<Vec<_>>(), vec![&raw]);
    let pattern = StereoBondConstraintsForm::from(raw);
    assert_eq!(view.satisfies(&pattern), !kindless && normalized.is_some());
    assert_eq!(view.is_compatible(&pattern), normalized.is_some());
}

#[rstest]
fn test_stereo_bond_constraints_view_derived() {
    // No projection is defined for stereo constraint kinds.
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 4],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(2)),
            (AtomId(2), AtomId(3), BondForm::from_order(1)),
        ],
        stereo_bonds: vec![(
            BondId(1),
            vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        )],
        ..Default::default()
    });
    let view = molecule.stereo_bond(StereoBondId(0)).constraints();
    assert_eq!(
        view.derived(StereoBondConstraintKey::Stereogenicity),
        Some(StereoBondConstraintForm::Stereogenicity(
            StereogenicityForm::Undetermined
        ))
    );
    assert_eq!(
        view.derived_complete(StereoBondConstraintKey::Stereogenicity),
        Some(StereoBondConstraintForm::Stereogenicity(
            StereogenicityForm::Undetermined
        ))
    );
    assert!(view.is_empty());
}

#[rustfmt::skip]
#[rstest]
#[case::stored(
    vec![
        StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }),
        StereoBondConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(false) }),
        StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
        StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
    ],
    BooleanForm::Lit(true),
    BooleanForm::Lit(false),
    TopicityRelationForm::Lit(Topicity::Homotopic),
    StereogenicityForm::Lit(Stereogenicity::Stereogenic),
)]
#[case::absent(
    vec![],
    BooleanForm::Undetermined,
    BooleanForm::Undetermined,
    TopicityRelationForm::Undetermined,
    StereogenicityForm::Undetermined,
)]
fn test_stereo_bond_constraints_view_container_methods(
    #[case] constraints: Vec<StereoBondConstraintForm>,
    #[case] ligand_symmetry: BooleanForm,
    #[case] fluxionality: BooleanForm,
    #[case] topicity: TopicityRelationForm,
    #[case] stereogenicity: StereogenicityForm,
) {
    let mut form = StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1));
    form.constraints = StereoBondConstraintsForm::from(constraints.clone());
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 4],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(2)),
            (AtomId(2), AtomId(3), BondForm::from_order(1)),
        ],
        stereo_bonds: vec![(
            BondId(1),
            vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
            ],
            form,
        )],
        ..Default::default()
    });
    let view = molecule.stereo_bond(StereoBondId(0)).constraints();
    let permutation = LigandPermutation(Permutation::from_image(&[1, 0, 2, 3]));
    assert_eq!(
        view.ligand_symmetry(OrientedLigandPermutation {
            permutation,
            orientation: Orientation::Proper,
        }),
        &ligand_symmetry
    );
    assert_eq!(view.fluxionality(permutation), &fluxionality);
    assert_eq!(
        view.topicity(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1))),
        &topicity
    );
    assert_eq!(view.stereogenicity(), &stereogenicity);
    assert_eq!(view.len(), constraints.len());
    assert_eq!(view.is_empty(), constraints.is_empty());
    assert_eq!(view.iter().collect::<Vec<_>>(), constraints.iter().collect::<Vec<_>>());
}
