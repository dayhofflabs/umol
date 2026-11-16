use pretty_assertions::assert_eq;

use super::*;
use crate::ir::GraphSymmetryConfig;
use crate::mol_dsl_concrete;

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
    assert_eq!(entity.attributes().constraints.asserted(key), normalized);
    assert_eq!(
        entity.attributes().constraints.asserted_complete(key),
        normalized
    );
    assert_eq!(view.iter().collect::<Vec<_>>(), vec![&raw]);
    if kindless {
        assert_eq!(view.derived(key), top);
        assert_eq!(view.derived_complete(key), top);
        assert_eq!(view.evidenced(key), expected);
        assert_eq!(view.evidenced_complete(key), expected);
        let pattern = StereoAtomConstraintsForm::from(raw);
        assert!(!view.satisfies(&pattern));
        assert_eq!(view.is_compatible(&pattern), normalized.is_some());
    }
}

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

#[rustfmt::skip]
#[rstest]
#[case::stereogenicity(
    StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
    StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
    StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Nonstereogenic)),
)]
#[case::topicity(
    StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Distinct) }),
    StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined }),
    StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
)]
#[case::invariant(
    StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }),
    StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined }),
    StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(false) }),
)]
#[case::moved(
    StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(false) }),
    StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined }),
    StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }),
)]
fn test_stereo_atom_constraints_view_derived(
    #[case] expected: StereoAtomConstraintForm,
    #[case] top: StereoAtomConstraintForm,
    #[case] contrary: StereoAtomConstraintForm,
    #[values("absent", "matching", "conflicting")] assertion: &str,
    #[values(false, true)] include_constraints: bool,
    #[values(false, true)] attach_on_constraints: bool,
) {
    let mut molecule = mol_dsl_concrete!(
        r#"{:atoms ["C" "H" "F" "Cl" "Br"]
           :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
           :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#
    );
    let asserted = match assertion {
        "absent" => top,
        "matching" => expected.clone(),
        "conflicting" => contrary,
        _ => unreachable!(),
    };
    if assertion != "absent" {
        molecule
            .stereo_atom_mut(StereoAtomId(0))
            .attributes_mut()
            .constraints
            .set(asserted.clone());
    }
    let original = molecule.clone();
    let symmetry = molecule
        .graph_symmetry(&GraphSymmetryConfig {
            include_constraints,
            ..Default::default()
        })
        .into_graph_symmetry();
    let stereo = molecule.stereo_atom(StereoAtomId(0));
    let view = if attach_on_constraints {
        stereo.constraints().with_symmetry(&symmetry)
    } else {
        stereo.with_symmetry(&symmetry).constraints()
    };
    let key = expected.key();
    assert_eq!(view.asserted(key), Some(asserted.clone()));
    assert_eq!(view.asserted_complete(key), Some(asserted));
    assert_eq!(view.derived(key), Some(expected.clone()));
    assert_eq!(view.derived_complete(key), Some(expected.clone()));
    let evidence = if assertion == "conflicting" {
        None
    } else {
        Some(expected.clone())
    };
    assert_eq!(view.evidenced(key), evidence);
    assert_eq!(view.evidenced_complete(key), evidence);
    let pattern = StereoAtomConstraintsForm::from(expected);
    assert_eq!(view.satisfies(&pattern), assertion != "conflicting");
    assert_eq!(view.is_compatible(&pattern), assertion != "conflicting");
    assert!(view.satisfies(&StereoAtomConstraintsForm::default()));
    assert!(view.is_compatible(&StereoAtomConstraintsForm::default()));
    assert_eq!(molecule, original);
}

#[rstest]
#[case::stereogenicity(StereoAtomConstraintKey::Stereogenicity)]
#[case::topicity(StereoAtomConstraintKey::Topicity(StereoLigandPair::new(
    StereoLigandPosition(0),
    StereoLigandPosition(1)
)))]
#[case::ligand_symmetry(StereoAtomConstraintKey::LigandSymmetry(OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }))]
#[should_panic(expected = "stereo view has symmetry data")]
fn test_stereo_atom_constraints_view_derived_error(
    #[case] key: StereoAtomConstraintKey,
    #[values(
        "derived",
        "derived_complete",
        "evidenced",
        "evidenced_complete",
        "satisfies",
        "is_compatible"
    )]
    method: &str,
) {
    let molecule = mol_dsl_concrete!(
        r#"{:atoms ["C" "H" "F" "Cl" "Br"]
           :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
           :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#
    );
    let view = molecule.stereo_atom(StereoAtomId(0)).constraints();
    let pattern = StereoAtomConstraintsForm::from(StereoAtomConstraintForm::undetermined(key));
    match method {
        "derived" => {
            view.derived(key);
        }
        "derived_complete" => {
            view.derived_complete(key);
        }
        "evidenced" => {
            view.evidenced(key);
        }
        "evidenced_complete" => {
            view.evidenced_complete(key);
        }
        "satisfies" => {
            view.satisfies(&pattern);
        }
        "is_compatible" => {
            view.is_compatible(&pattern);
        }
        _ => unreachable!(),
    }
}

#[rstest]
#[case::stereogenicity(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::lit_set([])))]
#[case::topicity(StereoAtomConstraintForm::Topicity(TopicityForm {
    pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)),
    relation: TopicityRelationForm::lit_set([]),
}))]
fn test_stereo_atom_constraints_view_evidenced(#[case] assertion: StereoAtomConstraintForm) {
    let mut molecule = mol_dsl_concrete!(
        r#"{:atoms ["C" "H" "F" "Cl" "Br"]
           :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
           :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#
    );
    let key = assertion.key();
    molecule
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .constraints
        .set(assertion);
    let view = molecule.stereo_atom(StereoAtomId(0)).constraints();
    assert_eq!(view.asserted(key), None);
    assert_eq!(view.evidenced(key), None);
    assert_eq!(view.evidenced_complete(key), None);
    let pattern = StereoAtomConstraintsForm::from(StereoAtomConstraintForm::undetermined(key));
    assert!(!view.satisfies(&pattern));
    assert!(!view.is_compatible(&pattern));
    assert!(view.satisfies(&StereoAtomConstraintsForm::default()));
    assert!(view.is_compatible(&StereoAtomConstraintsForm::default()));
}

#[rustfmt::skip]
#[rstest]
#[case::stored(
    vec![
        StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }),
        StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
        StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
    ],
    BooleanForm::Lit(true),
    TopicityRelationForm::Lit(Topicity::Homotopic),
    StereogenicityForm::Lit(Stereogenicity::Stereogenic),
)]
#[case::absent(
    vec![],
    BooleanForm::Undetermined,
    TopicityRelationForm::Undetermined,
    StereogenicityForm::Undetermined,
)]
fn test_stereo_atom_constraints_view_container_methods(
    #[case] constraints: Vec<StereoAtomConstraintForm>,
    #[case] ligand_symmetry: BooleanForm,
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
    assert_eq!(
        view.topicity(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1))),
        &topicity
    );
    assert_eq!(view.stereogenicity(), &stereogenicity);
    assert_eq!(view.len(), constraints.len());
    assert_eq!(view.is_empty(), constraints.is_empty());
    assert_eq!(view.iter().collect::<Vec<_>>(), constraints.iter().collect::<Vec<_>>());
}
