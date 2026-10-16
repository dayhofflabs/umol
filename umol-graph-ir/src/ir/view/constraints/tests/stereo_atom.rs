use pretty_assertions::assert_eq;

use super::*;

#[rstest]
#[case::present(
    Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
        Stereogenicity::Stereogenic,
    ))),
    Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
        Stereogenicity::Stereogenic,
    )))
)]
#[case::absent(
    None,
    Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined))
)]
fn test_stereo_atom_constraints_view_asserted_complete(
    #[case] set: Option<StereoAtomConstraintForm>,
    #[case] expected: Option<StereoAtomConstraintForm>,
) {
    let mut form = StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1));
    if let Some(constraint) = set {
        form.constraints.set(constraint);
    }
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 5],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(0), AtomId(4), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            form,
        )],
        ..Default::default()
    });
    assert_eq!(
        molecule
            .stereo_atom(StereoAtomId(0))
            .constraints()
            .asserted_complete(StereoAtomConstraintKey::Stereogenicity),
        expected
    );
}

#[rstest]
fn test_stereo_atom_constraints_view_asserted() {
    let mut form = StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1));
    form.constraints
        .set(StereoAtomConstraintForm::Stereogenicity(
            StereogenicityForm::Lit(Stereogenicity::Stereogenic),
        ));
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 5],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(0), AtomId(4), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            form,
        )],
        ..Default::default()
    });
    let view = molecule.stereo_atom(StereoAtomId(0)).constraints();
    assert_eq!(
        view.asserted(StereoAtomConstraintKey::Stereogenicity),
        Some(StereoAtomConstraintForm::Stereogenicity(
            StereogenicityForm::Lit(Stereogenicity::Stereogenic),
        ))
    );
    // No projection is defined for stereo constraint kinds.
    assert_eq!(
        view.derived(StereoAtomConstraintKey::Stereogenicity),
        Some(StereoAtomConstraintForm::Stereogenicity(
            StereogenicityForm::Undetermined
        ))
    );
    assert_eq!(
        view.derived_complete(StereoAtomConstraintKey::Stereogenicity),
        Some(StereoAtomConstraintForm::Stereogenicity(
            StereogenicityForm::Undetermined
        ))
    );
}

#[rstest]
#[case::symmetry(
    StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm {
        permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper },
        invariant: BooleanForm::Lit(true),
    }),
    Some(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm {
        permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper },
        invariant: BooleanForm::Lit(true),
    })),
)]
#[case::fluxionality(
    StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::identity(4)), active: BooleanForm::Lit(false) }),
    Some(StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::identity(4)), active: BooleanForm::Lit(false) })),
)]
#[case::topicity(
    StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::lit_set([Topicity::Homotopic]) }),
    Some(StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) })),
)]
#[case::stereogenicity(
    StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::lit_set([Stereogenicity::Stereogenic])),
    Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
)]
#[case::topicity_bottom(
    StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::lit_set([]) }),
    None,
)]
#[case::stereogenicity_bottom(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::lit_set([])), None)]
fn test_stereo_atom_constraints_view_asserted_constraints(
    #[case] raw: StereoAtomConstraintForm,
    #[case] normalized: Option<StereoAtomConstraintForm>,
    #[values(false, true)] kindless: bool,
) {
    let form = StereoAtomForm {
        configuration: if kindless {
            StereoConfigurationForm::Undetermined
        } else {
            StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined)
        },
        constraints: raw.clone().into(),
    };
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 6],
        bonds: [(0, 1), (0, 2), (0, 3), (0, 4)]
            .into_iter()
            .map(|(a, b)| (AtomId(a), AtomId(b), BondForm::from_order(1)))
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            [1, 2, 3, 4]
                .into_iter()
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                .collect(),
            form,
        )],
        ..Default::default()
    });
    let entity = molecule.stereo_atom(StereoAtomId(0));
    let view = entity.constraints();
    let key = raw.key();
    let top = Some(StereoAtomConstraintForm::undetermined(key));
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
    let pattern = StereoAtomConstraintsForm::from(raw);
    assert_eq!(view.satisfies(&pattern), !kindless && normalized.is_some());
    assert_eq!(view.is_compatible(&pattern), normalized.is_some());
}

#[rustfmt::skip]
#[rstest]
#[case::stored(
    vec![
        StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }),
        StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), active: BooleanForm::Lit(false) }),
        StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
        StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
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
fn test_stereo_atom_constraints_view_container_methods(
    #[case] constraints: Vec<StereoAtomConstraintForm>,
    #[case] ligand_symmetry: BooleanForm,
    #[case] fluxionality: BooleanForm,
    #[case] topicity: TopicityRelationForm,
    #[case] stereogenicity: StereogenicityForm,
) {
    let mut form = StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1));
    form.constraints = StereoAtomConstraintsForm::from(constraints.clone());
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 5],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(0), AtomId(4), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            form,
        )],
        ..Default::default()
    });
    let view = molecule.stereo_atom(StereoAtomId(0)).constraints();
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
