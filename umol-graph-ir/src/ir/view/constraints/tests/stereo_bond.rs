use pretty_assertions::assert_eq;

use super::*;
use crate::ir::GraphSymmetryConfig;
use crate::mol_dsl_concrete;

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
        let pattern = StereoBondConstraintsForm::from(raw);
        assert!(!view.satisfies(&pattern));
        assert_eq!(view.is_compatible(&pattern), normalized.is_some());
    }
}

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

#[rustfmt::skip]
#[rstest]
#[case::stereogenicity(
    StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
    StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
    StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Nonstereogenic)),
)]
#[case::topicity(
    StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
    StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2)), relation: TopicityRelationForm::Undetermined }),
    StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2)), relation: TopicityRelationForm::Lit(Topicity::Distinct) }),
)]
#[case::invariant(
    StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[2, 3, 0, 1])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }),
    StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[2, 3, 0, 1])), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined }),
    StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[2, 3, 0, 1])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(false) }),
)]
#[case::moved(
    StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(false) }),
    StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Undetermined }),
    StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }),
)]
fn test_stereo_bond_constraints_view_derived(
    #[case] expected: StereoBondConstraintForm,
    #[case] top: StereoBondConstraintForm,
    #[case] contrary: StereoBondConstraintForm,
    #[values("absent", "matching", "conflicting")] assertion: &str,
    #[values(false, true)] include_constraints: bool,
    #[values(false, true)] attach_on_constraints: bool,
) {
    let mut molecule = mol_dsl_concrete!(
        r#"{:atoms ["C" "C" "F" "Cl" "F" "Cl"]
           :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
           :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct0"}]}"#
    );
    let asserted = match assertion {
        "absent" => top,
        "matching" => expected.clone(),
        "conflicting" => contrary,
        _ => unreachable!(),
    };
    if assertion != "absent" {
        molecule
            .stereo_bond_mut(StereoBondId(0))
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
    let stereo = molecule.stereo_bond(StereoBondId(0));
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
    let pattern = StereoBondConstraintsForm::from(expected);
    assert_eq!(view.satisfies(&pattern), assertion != "conflicting");
    assert_eq!(view.is_compatible(&pattern), assertion != "conflicting");
    assert!(view.satisfies(&StereoBondConstraintsForm::default()));
    assert!(view.is_compatible(&StereoBondConstraintsForm::default()));
    assert_eq!(molecule, original);
}

#[rstest]
#[case::stereogenicity(StereoBondConstraintKey::Stereogenicity)]
#[case::topicity(StereoBondConstraintKey::Topicity(StereoLigandPair::new(
    StereoLigandPosition(0),
    StereoLigandPosition(2)
)))]
#[case::ligand_symmetry(StereoBondConstraintKey::LigandSymmetry(OrientedLigandPermutation { permutation: LigandPermutation(Permutation::identity(4)), orientation: Orientation::Proper }))]
#[should_panic(expected = "stereo view has symmetry data")]
fn test_stereo_bond_constraints_view_derived_error(
    #[case] key: StereoBondConstraintKey,
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
        r#"{:atoms ["C" "C" "F" "Cl" "F" "Cl"]
           :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
           :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct0"}]}"#
    );
    let view = molecule.stereo_bond(StereoBondId(0)).constraints();
    let pattern = StereoBondConstraintsForm::from(StereoBondConstraintForm::undetermined(key));
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
#[case::stereogenicity(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::lit_set([])))]
#[case::topicity(StereoBondConstraintForm::Topicity(TopicityForm {
    pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)),
    relation: TopicityRelationForm::lit_set([]),
}))]
fn test_stereo_bond_constraints_view_evidenced(#[case] assertion: StereoBondConstraintForm) {
    let mut molecule = mol_dsl_concrete!(
        r#"{:atoms ["C" "C" "F" "Cl" "F" "Cl"]
           :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
           :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct0"}]}"#
    );
    let key = assertion.key();
    molecule
        .stereo_bond_mut(StereoBondId(0))
        .attributes_mut()
        .constraints
        .set(assertion);
    let view = molecule.stereo_bond(StereoBondId(0)).constraints();
    assert_eq!(view.asserted(key), None);
    assert_eq!(view.evidenced(key), None);
    assert_eq!(view.evidenced_complete(key), None);
    let pattern = StereoBondConstraintsForm::from(StereoBondConstraintForm::undetermined(key));
    assert!(!view.satisfies(&pattern));
    assert!(!view.is_compatible(&pattern));
    assert!(view.satisfies(&StereoBondConstraintsForm::default()));
    assert!(view.is_compatible(&StereoBondConstraintsForm::default()));
}

#[rustfmt::skip]
#[rstest]
#[case::stored(
    vec![
        StereoBondConstraintForm::LigandSymmetry(LigandSymmetryForm { permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])), orientation: Orientation::Proper }, invariant: BooleanForm::Lit(true) }),
        StereoBondConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }),
        StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
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
fn test_stereo_bond_constraints_view_container_methods(
    #[case] constraints: Vec<StereoBondConstraintForm>,
    #[case] ligand_symmetry: BooleanForm,
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
    assert_eq!(
        view.topicity(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1))),
        &topicity
    );
    assert_eq!(view.stereogenicity(), &stereogenicity);
    assert_eq!(view.len(), constraints.len());
    assert_eq!(view.is_empty(), constraints.is_empty());
    assert_eq!(view.iter().collect::<Vec<_>>(), constraints.iter().collect::<Vec<_>>());
}
