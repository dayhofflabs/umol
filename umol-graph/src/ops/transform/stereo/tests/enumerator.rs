use std::collections::BTreeSet;

use rstest::{fixture, rstest};
use umol_chem::element::Element;
use umol_graph_ir::ir::{
    AtomConstraintForm, AtomForm, AtomId, BondConstraintForm, BondForm, BondId, CisTransStereoForm,
    Constraint, DativeBondConstraintForm, DativeBondId, Molecule, MoleculeEntries, NumForm,
    RingScope, StereoAtomConstraintForm, StereoAtomForm, StereoAtomId, StereoBondId,
    StereoConfigurationForm, StereoCoset, StereoKind, StereoLigand, StereoLigandKind, StereoTerm,
    Stereogenicity, StereogenicityForm, TetrahedralStereoForm,
};
use umol_graph_ir::mol_dsl_concrete;

use crate::ops::model::StereoModel;
use crate::ops::stereo::{StereoInconsistency, StereoPerceptionError, StereoSelection, StereoSite};
use crate::ops::transform::tests::stereo::reference::absolute_assignments;
use crate::ops::transform::{
    EnumerateStereoConfig, EnumerateStereoError, EnumerateStereoIter, StereoEnumerator, Transformer,
};
use crate::ops::validate::constraint::ConstraintInvariantsError;

#[fixture]
fn mixed_stereo() -> Molecule {
    mol_dsl_concrete!(
        r#"{:atoms ["C#D4" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2#R!"] [5 7 "1"] [6 8 "1"]]
        :stereo-atoms [{:site 0 :ligands [2 1 3 4] :attrs "Th*"}]
        :stereo-bonds [{:site 4 :ligands [7 [:h 5] 8 [:h 6]] :attrs "Ct*"}]}"#
    )
}

#[rstest]
#[case::tetrahedral(
    mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]}"#),
    vec![StereoSite::Atom(AtomId(0))],
    vec![mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#)],
)]
#[case::cis_trans(
    mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"] :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
        :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct*"}]}"#),
    vec![StereoSite::Bond(BondId(0))],
    vec![mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"] :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
        :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct0"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"] :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
        :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct1"}]}"#)],
)]
fn test_stereo_enumerator_new(
    #[case] molecule: Molecule,
    #[case] order: Vec<StereoSite>,
    #[case] expected: Vec<Molecule>,
) {
    let transformer = StereoEnumerator::new(&StereoModel::default(), order.clone());
    let configured = StereoEnumerator::with_config(
        &StereoModel::default(),
        EnumerateStereoConfig::default(),
        order,
    );
    let original = molecule.clone();
    assert_eq!(
        transformer.transform_iter(&molecule).collect::<Vec<_>>(),
        expected
    );
    assert_eq!(
        configured.transform_iter(&molecule).collect::<Vec<_>>(),
        expected
    );
    assert_eq!(
        Transformer::transform_iter(&transformer, &molecule).collect::<Vec<_>>(),
        expected
    );
    assert_eq!(
        transformer.transform(molecule.clone()),
        Ok(expected[0].clone())
    );
    let mut in_place = molecule.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected[0]);
    assert_eq!(molecule, original);
    for result in expected {
        assert_eq!(result.clone().edit().finish(), Ok(result));
    }
}

#[rstest]
#[case::all(StereoSelection::default(), vec![(Some(0), Some(0)), (Some(0), Some(1)), (Some(1), Some(0)), (Some(1), Some(1))])]
#[case::atom(StereoSelection { sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(0))])), ..Default::default() }, vec![(Some(0), None), (Some(1), None)])]
#[case::bond(StereoSelection { sites: Some(BTreeSet::from([StereoSite::Bond(BondId(4))])), ..Default::default() }, vec![(None, Some(0)), (None, Some(1))])]
#[case::tetrahedral(StereoSelection { kinds: Some(BTreeSet::from([StereoKind::Tetrahedral])), ..Default::default() }, vec![(Some(0), None), (Some(1), None)])]
#[case::cis_trans(StereoSelection { kinds: Some(BTreeSet::from([StereoKind::CisTrans])), ..Default::default() }, vec![(None, Some(0)), (None, Some(1))])]
#[case::other_kind(StereoSelection { kinds: Some(BTreeSet::from([StereoKind::SquarePlanar])), ..Default::default() }, vec![(None, None)])]
#[case::stereogenic(StereoSelection { stereogenicities: Some(BTreeSet::from([Stereogenicity::Stereogenic])), ..Default::default() }, vec![(Some(0), Some(0)), (Some(0), Some(1)), (Some(1), Some(0)), (Some(1), Some(1))])]
#[case::symmetric(StereoSelection { stereogenicities: Some(BTreeSet::from([Stereogenicity::Symmetric])), ..Default::default() }, vec![(None, None)])]
#[case::intersection(StereoSelection { sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(0))])), kinds: Some(BTreeSet::from([StereoKind::CisTrans])), ..Default::default() }, vec![(None, None)])]
#[case::empty_sites(StereoSelection { sites: Some(BTreeSet::new()), ..Default::default() }, vec![(None, None)])]
#[case::empty_kinds(StereoSelection { kinds: Some(BTreeSet::new()), ..Default::default() }, vec![(None, None)])]
#[case::empty_classes(StereoSelection { stereogenicities: Some(BTreeSet::new()), ..Default::default() }, vec![(None, None)])]
#[case::no_entity(StereoSelection { sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(1))])), ..Default::default() }, vec![(None, None)])]
fn test_stereo_enumerator_with_config(
    mixed_stereo: Molecule,
    #[case] selection: StereoSelection,
    #[case] values: Vec<(Option<u32>, Option<u32>)>,
) {
    let transformer = StereoEnumerator::with_config(
        &StereoModel::default(),
        EnumerateStereoConfig {
            selection,
            ..Default::default()
        },
        vec![
            StereoSite::Atom(AtomId(99)),
            StereoSite::Atom(AtomId(0)),
            StereoSite::Atom(AtomId(99)),
            StereoSite::Bond(BondId(4)),
        ],
    );
    let expected: Vec<_> = values
        .into_iter()
        .map(|(atom, bond)| {
            let mut result = mixed_stereo.clone();
            if let Some(value) = atom {
                result
                    .stereo_atom_mut(StereoAtomId(0))
                    .attributes_mut()
                    .configuration = (StereoKind::Tetrahedral, value).into();
            }
            if let Some(value) = bond {
                result
                    .stereo_bond_mut(StereoBondId(0))
                    .attributes_mut()
                    .configuration = (StereoKind::CisTrans, value).into();
            }
            result
        })
        .collect();
    let original = mixed_stereo.clone();
    assert_eq!(
        transformer
            .transform_iter(&mixed_stereo)
            .collect::<Vec<_>>(),
        expected
    );
    assert_eq!(
        transformer.transform(mixed_stereo.clone()),
        Ok(expected[0].clone())
    );
    let mut in_place = mixed_stereo.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected[0]);
    assert_eq!(mixed_stereo, original);
}

#[rstest]
fn test_stereo_enumerator_with_config_model(mut mixed_stereo: Molecule) {
    let mut model = StereoModel {
        kind_models: [None, None, None, None, None, None],
        stereo_bond_minimum_ring_size: 64,
        ..Default::default()
    };
    let transformer = StereoEnumerator::with_config(
        &model,
        EnumerateStereoConfig {
            selection: StereoSelection {
                kinds: Some(BTreeSet::from([StereoKind::Tetrahedral])),
                stereogenicities: Some(BTreeSet::from([Stereogenicity::Stereogenic])),
                ..Default::default()
            },
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(0))],
    );
    model = StereoModel::default();
    model.kind_models[StereoKind::Tetrahedral as usize]
        .as_mut()
        .unwrap()
        .fluxionality = true;
    mixed_stereo.atom_mut(AtomId(1)).attributes_mut().charge = NumForm::Undetermined;
    let expected: Vec<_> = [0, 1]
        .into_iter()
        .map(|value| {
            let mut molecule = mixed_stereo.clone();
            molecule
                .stereo_atom_mut(StereoAtomId(0))
                .attributes_mut()
                .configuration = (StereoKind::Tetrahedral, value).into();
            molecule
        })
        .collect();
    assert_eq!(
        transformer
            .transform_iter(&mixed_stereo)
            .collect::<Vec<_>>(),
        expected
    );
    assert_eq!(
        transformer.transform(mixed_stereo.clone()),
        Ok(expected[0].clone())
    );
    let other = StereoEnumerator::with_config(
        &model,
        EnumerateStereoConfig {
            selection: StereoSelection {
                kinds: Some(BTreeSet::from([StereoKind::Tetrahedral])),
                ..Default::default()
            },
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(0))],
    );
    assert_eq!(
        other.transform(mixed_stereo),
        Err(EnumerateStereoError::Perception(
            StereoPerceptionError::UnsupportedFluxionality {
                kind: StereoKind::Tetrahedral
            }
        ))
    );
}

#[rstest]
#[case::fixed_and_open(vec![StereoCoset::Lit(1), StereoCoset::Undetermined])]
#[case::finite(vec![StereoCoset::lit_set([0, 1]), StereoCoset::lit_set([1])])]
#[case::both_fixed(vec![StereoCoset::Lit(0), StereoCoset::Lit(1)])]
#[case::singletons(vec![StereoCoset::lit_set([1]), StereoCoset::lit_set([0])])]
fn test_stereo_enumerator_transform_domains(
    mut mixed_stereo: Molecule,
    #[case] source: Vec<StereoCoset>,
    #[values(false, true)] reverse: bool,
    #[values(false, true)] classified: bool,
) {
    mixed_stereo
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration =
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, source[0].clone());
    mixed_stereo
        .stereo_bond_mut(StereoBondId(0))
        .attributes_mut()
        .configuration = StereoConfigurationForm::Kinded(StereoKind::CisTrans, source[1].clone());
    let order = if reverse { vec![1, 0] } else { vec![0, 1] };
    let sites = [StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(4))];
    let expected: Vec<_> = absolute_assignments(&source, &BTreeSet::from([0, 1]), &order, |_| true)
        .unwrap()
        .into_iter()
        .map(|values| {
            let mut molecule = mixed_stereo.clone();
            molecule
                .stereo_atom_mut(StereoAtomId(0))
                .attributes_mut()
                .configuration =
                StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, values[0].clone());
            molecule
                .stereo_bond_mut(StereoBondId(0))
                .attributes_mut()
                .configuration =
                StereoConfigurationForm::Kinded(StereoKind::CisTrans, values[1].clone());
            molecule
        })
        .collect();
    let transformer = StereoEnumerator::with_config(
        &StereoModel::default(),
        EnumerateStereoConfig {
            selection: StereoSelection {
                kinds: None,
                stereogenicities: classified.then(|| BTreeSet::from([Stereogenicity::Stereogenic])),
                ..Default::default()
            },
            ..Default::default()
        },
        order.into_iter().map(|index| sites[index]).collect(),
    );
    let original = mixed_stereo.clone();
    assert_eq!(
        transformer
            .transform_iter(&mixed_stereo)
            .collect::<Vec<_>>(),
        expected
    );
    assert_eq!(
        transformer.transform(mixed_stereo.clone()),
        Ok(expected[0].clone())
    );
    let mut in_place = mixed_stereo.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected[0]);
    assert_eq!(mixed_stereo, original);
}

#[rstest]
#[case::swap(StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(0))), 1)]
#[case::mirror(StereoCoset::term(StereoTerm::mirror(StereoTerm::Lit(1))), 0)]
fn test_stereo_enumerator_transform_expression(
    mut mixed_stereo: Molecule,
    #[case] coset: StereoCoset,
    #[case] expected_value: u32,
    #[values(false, true)] classified: bool,
) {
    mixed_stereo
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, coset);
    let mut expected = mixed_stereo.clone();
    expected
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = (StereoKind::Tetrahedral, expected_value).into();
    let transformer = StereoEnumerator::with_config(
        &StereoModel::default(),
        EnumerateStereoConfig {
            selection: StereoSelection {
                sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(0))])),
                stereogenicities: classified.then(|| BTreeSet::from([Stereogenicity::Stereogenic])),
                ..Default::default()
            },
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(0))],
    );
    let original = mixed_stereo.clone();
    assert_eq!(
        transformer
            .transform_iter(&mixed_stereo)
            .collect::<Vec<_>>(),
        vec![expected.clone()]
    );
    let mut in_place = mixed_stereo.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected);
    assert_eq!(mixed_stereo, original);
    assert_eq!(transformer.transform(mixed_stereo), Ok(expected));
}

#[rstest]
#[case::source_assertions(0, 1, vec![(1, 1)])]
#[case::first_assertions(1, 0, vec![(0, 0)])]
fn test_stereo_enumerator_transform_evidence(
    mut mixed_stereo: Molecule,
    #[case] asserted_atom: u32,
    #[case] asserted_bond: u32,
    #[case] values: Vec<(u32, u32)>,
) {
    mixed_stereo
        .atom_mut(AtomId(0))
        .attributes_mut()
        .constraints
        .set(AtomConstraintForm::TetrahedralStereo(
            TetrahedralStereoForm::Stereo(asserted_atom.into()),
        ));
    mixed_stereo
        .bond_mut(BondId(4))
        .attributes_mut()
        .constraints
        .set(BondConstraintForm::CisTransStereo(
            CisTransStereoForm::Stereo(asserted_bond.into()),
        ));
    let original = mixed_stereo.clone();
    let expected: Vec<_> = values
        .into_iter()
        .map(|(atom, bond)| {
            let mut molecule = mixed_stereo.clone();
            molecule
                .stereo_atom_mut(StereoAtomId(0))
                .attributes_mut()
                .configuration = (StereoKind::Tetrahedral, atom).into();
            molecule
                .stereo_bond_mut(StereoBondId(0))
                .attributes_mut()
                .configuration = (StereoKind::CisTrans, bond).into();
            molecule
        })
        .collect();
    let transformer = StereoEnumerator::new(
        &StereoModel::default(),
        vec![StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(4))],
    );
    assert_eq!(
        transformer
            .transform_iter(&mixed_stereo)
            .collect::<Vec<_>>(),
        expected
    );
    let mut in_place = mixed_stereo.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected[0]);
    assert_eq!(mixed_stereo, original);
    assert_eq!(transformer.transform(mixed_stereo), Ok(expected[0].clone()));
}

#[rstest]
fn test_stereo_enumerator_transform_constraint(
    mixed_stereo: Molecule,
    #[values(false, true)] reverse: bool,
) {
    let mut editor = mixed_stereo.edit();
    editor.constraints_mut().push(Constraint::Or(
        [0u32, 1]
            .map(|value| {
                Constraint::And(vec![
                    Constraint::Atom(
                        AtomId(0),
                        AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Stereo(
                            value.into(),
                        )),
                    ),
                    Constraint::Bond(
                        BondId(4),
                        BondConstraintForm::CisTransStereo(CisTransStereoForm::Stereo(
                            value.into(),
                        )),
                    ),
                ])
            })
            .into(),
    ));
    let molecule = editor.finish().unwrap();
    let order = if reverse { vec![1, 0] } else { vec![0, 1] };
    let expected: Vec<_> = absolute_assignments(
        &[StereoCoset::Undetermined, StereoCoset::Undetermined],
        &BTreeSet::from([0, 1]),
        &order,
        |values| values[0] != values[1],
    )
    .unwrap()
    .into_iter()
    .map(|values| {
        let mut result = molecule.clone();
        result
            .stereo_atom_mut(StereoAtomId(0))
            .attributes_mut()
            .configuration =
            StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, values[0].clone());
        result
            .stereo_bond_mut(StereoBondId(0))
            .attributes_mut()
            .configuration =
            StereoConfigurationForm::Kinded(StereoKind::CisTrans, values[1].clone());
        result
    })
    .collect();
    let sites = [StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(4))];
    let transformer = StereoEnumerator::new(
        &StereoModel::default(),
        order.into_iter().map(|index| sites[index]).collect(),
    );
    let original = molecule.clone();
    assert_eq!(
        transformer.transform_iter(&molecule).collect::<Vec<_>>(),
        expected
    );
    let mut in_place = molecule.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected[0]);
    assert_eq!(molecule, original);
    assert_eq!(transformer.transform(molecule), Ok(expected[0].clone()));
}

#[rstest]
#[case::kindless(StereoConfigurationForm::Undetermined)]
#[case::variable(StereoConfigurationForm::Kinded(
    StereoKind::CisTrans,
    StereoCoset::term(StereoTerm::var("retained"))
))]
#[case::other_kind((StereoKind::SquarePlanar, 0).into())]
fn test_stereo_enumerator_transform_scope(
    mut mixed_stereo: Molecule,
    #[case] retained: StereoConfigurationForm,
    #[values(false, true)] classified: bool,
) {
    mixed_stereo
        .stereo_bond_mut(StereoBondId(0))
        .attributes_mut()
        .configuration = retained;
    let transformer = StereoEnumerator::with_config(
        &StereoModel::default(),
        EnumerateStereoConfig {
            selection: StereoSelection {
                kinds: Some(BTreeSet::from([StereoKind::Tetrahedral])),
                stereogenicities: classified.then(|| BTreeSet::from([Stereogenicity::Stereogenic])),
                ..Default::default()
            },
            ..Default::default()
        },
        vec![
            StereoSite::Bond(BondId(4)),
            StereoSite::Atom(AtomId(0)),
            StereoSite::Bond(BondId(4)),
        ],
    );
    let expected: Vec<_> = [0, 1]
        .into_iter()
        .map(|value| {
            let mut molecule = mixed_stereo.clone();
            molecule
                .stereo_atom_mut(StereoAtomId(0))
                .attributes_mut()
                .configuration = (StereoKind::Tetrahedral, value).into();
            molecule
        })
        .collect();
    let original = mixed_stereo.clone();
    assert_eq!(
        transformer
            .transform_iter(&mixed_stereo)
            .collect::<Vec<_>>(),
        expected
    );
    assert_eq!(
        transformer.transform(mixed_stereo.clone()),
        Ok(expected[0].clone())
    );
    let mut in_place = mixed_stereo.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected[0]);
    assert_eq!(mixed_stereo, original);
}

#[rstest]
#[case::stereogenic(mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]}"#), Stereogenicity::Stereogenic)]
#[case::prochiral(mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "H" "H"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]}"#), Stereogenicity::Prochiral)]
#[case::symmetric(mol_dsl_concrete!(r#"{:atoms ["C" "H" "H" "H" "H"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]}"#), Stereogenicity::Symmetric)]
fn test_stereo_enumerator_transform_stereogenicity(
    #[case] mut molecule: Molecule,
    #[case] class: Stereogenicity,
    #[values(
        None,
        Some(Stereogenicity::Stereogenic),
        Some(Stereogenicity::Prochiral),
        Some(Stereogenicity::Symmetric)
    )]
    selected_class: Option<Stereogenicity>,
) {
    molecule
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .constraints
        .set(StereoAtomConstraintForm::Stereogenicity(
            StereogenicityForm::Lit(class),
        ));
    let expected = if selected_class.is_none_or(|selected| selected == class) {
        [0, 1]
            .into_iter()
            .map(|value| {
                let mut result = molecule.clone();
                result
                    .stereo_atom_mut(StereoAtomId(0))
                    .attributes_mut()
                    .configuration = (StereoKind::Tetrahedral, value).into();
                result
            })
            .collect()
    } else {
        vec![molecule.clone()]
    };
    let transformer = StereoEnumerator::with_config(
        &StereoModel::default(),
        EnumerateStereoConfig {
            selection: StereoSelection {
                stereogenicities: selected_class.map(|class| BTreeSet::from([class])),
                ..Default::default()
            },
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(0))],
    );
    let original = molecule.clone();
    assert_eq!(
        transformer.transform_iter(&molecule).collect::<Vec<_>>(),
        expected
    );
    let mut in_place = molecule.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected[0]);
    assert_eq!(molecule, original);
    assert_eq!(transformer.transform(molecule), Ok(expected[0].clone()));
}

#[rstest]
fn test_stereo_enumerator_transform_ring() {
    let molecule = mol_dsl_concrete!(
        r#"{:atoms ["C#h" "C#h" "C#h2" "C#h2" "C#h2" "C#h2"]
        :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
        :stereo-bonds [{:site 0 :ligands [5 [:h 0] 2 [:h 1]] :attrs "Ct*"}]}"#
    );
    let expected: Vec<_> = [0, 1]
        .into_iter()
        .map(|value| {
            let mut result = molecule.clone();
            result
                .stereo_bond_mut(StereoBondId(0))
                .attributes_mut()
                .configuration = (StereoKind::CisTrans, value).into();
            result
        })
        .collect();
    let transformer =
        StereoEnumerator::new(&StereoModel::default(), vec![StereoSite::Bond(BondId(0))]);
    assert_eq!(
        transformer.transform_iter(&molecule).collect::<Vec<_>>(),
        expected
    );
    assert_eq!(transformer.transform(molecule), Ok(expected[0].clone()));
}

#[rstest]
#[case::empty_domain(StereoCoset::lit_set([]), TetrahedralStereoForm::Undetermined)]
#[case::fixed_conflict(
    StereoCoset::Lit(0),
    TetrahedralStereoForm::Stereo(StereoCoset::Lit(0))
)]
#[case::finite_conflict(StereoCoset::lit_set([0]), TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)))]
#[case::negative(StereoCoset::Undetermined, TetrahedralStereoForm::NotStereo)]
#[case::empty_assertion(StereoCoset::Undetermined, TetrahedralStereoForm::Stereo(StereoCoset::lit_set([])))]
fn test_stereo_enumerator_transform_rejection(
    mut mixed_stereo: Molecule,
    #[case] coset: StereoCoset,
    #[case] asserted: TetrahedralStereoForm,
) {
    mixed_stereo
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, coset);
    mixed_stereo
        .atom_mut(AtomId(0))
        .attributes_mut()
        .constraints
        .set(AtomConstraintForm::TetrahedralStereo(asserted));
    let transformer = StereoEnumerator::new(
        &StereoModel::default(),
        vec![StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(4))],
    );
    let original = mixed_stereo.clone();
    let mut in_place = mixed_stereo.clone();
    assert_eq!(
        transformer.transform_into(&mut in_place),
        Err(EnumerateStereoError::NoAssignment)
    );
    assert_eq!(in_place, original);
    let mut iter = transformer.transform_iter(&mixed_stereo);
    assert_eq!(iter.next(), None);
    assert_eq!(iter.next(), None);
    assert_eq!(mixed_stereo, original);
    assert_eq!(
        transformer.transform(mixed_stereo),
        Err(EnumerateStereoError::NoAssignment)
    );
}

#[rstest]
#[case::frame(mol_dsl_concrete!(r#"{:atoms ["C#h" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3] :attrs "Th*"}]}"#), vec![StereoSite::Atom(AtomId(0))],
    StereoPerceptionError::Inconsistency(StereoInconsistency::StereoAtomFailure { stereo_atom: StereoAtomId(0) }))]
#[case::later_bond(mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2"] [5 7 "1"] [6 8 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]
    :stereo-bonds [{:site 4 :ligands [7 [:h 5] 8 [:h 6]] :attrs "Ct?configuration"}]}"#),
    vec![StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(4))],
    StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Bond(BondId(4)) })]
fn test_stereo_enumerator_transform_error(
    #[case] molecule: Molecule,
    #[case] order: Vec<StereoSite>,
    #[case] cause: StereoPerceptionError,
) {
    let transformer = StereoEnumerator::new(&StereoModel::default(), order);
    let error = EnumerateStereoError::Perception(cause);
    let original = molecule.clone();
    let mut in_place = molecule.clone();
    assert_eq!(
        transformer.transform_into(&mut in_place),
        Err(error.clone())
    );
    assert_eq!(in_place, original);
    assert_eq!(
        transformer.transform_iter(&molecule).collect::<Vec<_>>(),
        Vec::<Molecule>::new()
    );
    assert_eq!(molecule, original);
    assert_eq!(transformer.transform(molecule), Err(error));
}

#[rstest]
fn test_stereo_enumerator_transform_constraint_error() {
    let molecule = mol_dsl_concrete!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I" "N" "B"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :dative-bonds [{:donors [5] :acceptor 6 :attrs "1"}]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]}"#
    );
    let mut editor = molecule.edit();
    editor.constraints_mut().push(Constraint::Or(vec![
        Constraint::Atom(
            AtomId(0),
            AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Stereo(0.into())),
        ),
        Constraint::DativeBond(
            DativeBondId(0),
            DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1),
        ),
    ]));
    let molecule = editor.finish().unwrap();
    let original = molecule.clone();
    let transformer =
        StereoEnumerator::new(&StereoModel::default(), vec![StereoSite::Atom(AtomId(0))]);
    let error = EnumerateStereoError::Constraint(
        ConstraintInvariantsError::DativeBondRingMembershipUnsupported {
            bond: DativeBondId(0),
        },
    );
    let mut in_place = molecule.clone();
    assert_eq!(
        transformer.transform_into(&mut in_place),
        Err(error.clone())
    );
    assert_eq!(in_place, original);
    assert_eq!(
        transformer.transform_iter(&molecule).collect::<Vec<_>>(),
        Vec::<Molecule>::new()
    );
    assert_eq!(molecule, original);
    assert_eq!(transformer.transform(molecule), Err(error));
}

#[rstest]
#[case::empty(Molecule::default(), vec![])]
#[case::no_stereo(mol_dsl_concrete!(r#"{:atoms ["C#T1" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#), vec![])]
#[case::fixed(mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [2 1 3 4] :attrs "Th1"}]}"#), vec![StereoSite::Atom(AtomId(0))])]
#[case::excluded_kind(mol_dsl_concrete!(r#"{:atoms ["Pt#c+2" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Sp*"}]}"#), vec![])]
#[case::excluded_kindless(mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "*"}]}"#), vec![])]
fn test_stereo_enumerator_transform_identity(
    #[case] molecule: Molecule,
    #[case] order: Vec<StereoSite>,
) {
    let transformer = StereoEnumerator::new(&StereoModel::default(), order);
    let mut results = transformer.transform_iter(&molecule);
    assert_eq!(results.next(), Some(molecule.clone()));
    assert_eq!(results.next(), None);
    assert_eq!(results.next(), None);
    let mut in_place = molecule.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, molecule);
    assert_eq!(transformer.transform(molecule.clone()), Ok(molecule));
}

#[rstest]
#[case::missing(vec![StereoSite::Atom(AtomId(0))], vec![StereoSite::Bond(BondId(4))], vec![])]
#[case::duplicate(vec![StereoSite::Atom(AtomId(0)), StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(4))], vec![], vec![StereoSite::Atom(AtomId(0))])]
#[case::both(vec![StereoSite::Atom(AtomId(0)), StereoSite::Atom(AtomId(0))], vec![StereoSite::Bond(BondId(4))], vec![StereoSite::Atom(AtomId(0))])]
#[case::duplicates(vec![StereoSite::Bond(BondId(4)), StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(4)), StereoSite::Atom(AtomId(0))], vec![], vec![StereoSite::Bond(BondId(4)), StereoSite::Atom(AtomId(0))])]
fn test_stereo_enumerator_transform_order(
    mixed_stereo: Molecule,
    #[case] order: Vec<StereoSite>,
    #[case] missing: Vec<StereoSite>,
    #[case] duplicates: Vec<StereoSite>,
) {
    let transformer = StereoEnumerator::new(&StereoModel::default(), order);
    let error = EnumerateStereoError::InvalidSiteOrder {
        missing,
        duplicates,
    };
    let original = mixed_stereo.clone();
    let mut in_place = mixed_stereo.clone();
    assert_eq!(
        transformer.transform_into(&mut in_place),
        Err(error.clone())
    );
    assert_eq!(in_place, original);
    let mut iter = transformer.transform_iter(&mixed_stereo);
    assert_eq!(iter.next(), None);
    assert_eq!(iter.next(), None);
    assert_eq!(mixed_stereo, original);
    assert_eq!(transformer.transform(mixed_stereo), Err(error));
}

#[rstest]
#[case::kindless(StereoConfigurationForm::Undetermined, StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Atom(AtomId(0)) })]
#[case::square_planar((StereoKind::SquarePlanar, 0).into(), StereoPerceptionError::UnsupportedKind { kind: StereoKind::SquarePlanar })]
#[case::variable(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::term(StereoTerm::var("configuration"))), StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Atom(AtomId(0)) })]
#[case::finite_variable(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::term(StereoTerm::var_in("configuration", [0, 1]))), StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Atom(AtomId(0)) })]
#[case::empty_variable(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::term(StereoTerm::var_in("configuration", []))), StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Atom(AtomId(0)) })]
#[case::invalid((StereoKind::Tetrahedral, 2).into(), StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Atom(AtomId(0)) })]
#[case::site_kind((StereoKind::CisTrans, 0).into(), StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Atom(AtomId(0)) })]
fn test_stereo_enumerator_transform_configuration(
    mut mixed_stereo: Molecule,
    #[case] configuration: StereoConfigurationForm,
    #[case] cause: StereoPerceptionError,
    #[values(false, true)] classified: bool,
) {
    mixed_stereo
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = configuration;
    let original = mixed_stereo.clone();
    let transformer = StereoEnumerator::with_config(
        &StereoModel::default(),
        EnumerateStereoConfig {
            selection: StereoSelection {
                sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(0))])),
                stereogenicities: classified.then(|| BTreeSet::from([Stereogenicity::Stereogenic])),
                ..Default::default()
            },
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(0))],
    );
    let error = EnumerateStereoError::Perception(cause);
    let mut in_place = mixed_stereo.clone();
    assert_eq!(
        transformer.transform_into(&mut in_place),
        Err(error.clone())
    );
    assert_eq!(in_place, original);
    assert_eq!(
        transformer
            .transform_iter(&mixed_stereo)
            .collect::<Vec<_>>(),
        Vec::<Molecule>::new()
    );
    assert_eq!(mixed_stereo, original);
    assert_eq!(transformer.transform(mixed_stereo), Err(error));
}

#[rstest]
#[case::invalid_atom(StereoSelection { sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(9))])), ..Default::default() }, StereoSite::Atom(AtomId(9)))]
#[case::invalid_bond(StereoSelection { sites: Some(BTreeSet::from([StereoSite::Bond(BondId(7))])), ..Default::default() }, StereoSite::Bond(BondId(7)))]
#[case::invalid_empty(StereoSelection { sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(9))])), kinds: Some(BTreeSet::new()), stereogenicities: Some(BTreeSet::new()) }, StereoSite::Atom(AtomId(9)))]
fn test_stereo_enumerator_transform_selection(
    mixed_stereo: Molecule,
    #[case] selection: StereoSelection,
    #[case] site: StereoSite,
) {
    let transformer = StereoEnumerator::with_config(
        &StereoModel::default(),
        EnumerateStereoConfig {
            selection,
            ..Default::default()
        },
        vec![],
    );
    let error = EnumerateStereoError::Perception(StereoPerceptionError::InvalidSite { site });
    let mut in_place = mixed_stereo.clone();
    assert_eq!(
        transformer.transform_into(&mut in_place),
        Err(error.clone())
    );
    assert_eq!(in_place, mixed_stereo);
    assert_eq!(
        transformer
            .transform_iter(&mixed_stereo)
            .collect::<Vec<_>>(),
        Vec::<Molecule>::new()
    );
    assert_eq!(transformer.transform(mixed_stereo), Err(error));
}

#[rstest]
fn test_stereo_enumerator_transform_reference(
    #[values(0u32, 1, 2, 4, 8)] count: u32,
    #[values(false, true)] reverse: bool,
) {
    let mut entries = MoleculeEntries::default();
    for component in 0..count {
        let base = component * 5;
        entries.atoms.extend(
            [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
                .map(|element| AtomForm::from_element(element).into_concrete()),
        );
        entries.bonds.extend((1..5).map(|ligand| {
            (
                AtomId(base),
                AtomId(base + ligand),
                BondForm::from_order(1).into_concrete(),
            )
        }));
        entries.stereo_atoms.push((
            AtomId(base),
            (1..5)
                .map(|ligand| StereoLigand::new(AtomId(base + ligand), StereoLigandKind::Atom))
                .collect(),
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
        ));
    }
    let molecule = Molecule::from_entries(entries);
    let mut order: Vec<_> = (0..count as usize).collect();
    if reverse {
        order.reverse();
    }
    let expected: Vec<_> = absolute_assignments(
        &vec![StereoCoset::Undetermined; count as usize],
        &(0..count as usize).collect(),
        &order,
        |_| true,
    )
    .unwrap()
    .into_iter()
    .map(|values| {
        let mut result = molecule.clone();
        for (index, coset) in values.into_iter().enumerate() {
            result
                .stereo_atom_mut(StereoAtomId(index as u32))
                .attributes_mut()
                .configuration = StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, coset);
        }
        result
    })
    .collect();
    let transformer = StereoEnumerator::new(
        &StereoModel::default(),
        order
            .into_iter()
            .map(|index| StereoSite::Atom(AtomId(index as u32 * 5)))
            .collect(),
    );
    assert_eq!(
        transformer.transform_iter(&molecule).collect::<Vec<_>>(),
        expected
    );
    assert_eq!(
        transformer.transform(molecule.clone()),
        Ok(expected[0].clone())
    );
    let mut in_place = molecule.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected[0]);
}

#[rstest]
fn test_stereo_enumerator_transform_iter(mixed_stereo: Molecule) {
    let original = mixed_stereo.clone();
    let mut source = mixed_stereo;
    let mut iter: EnumerateStereoIter = {
        let transformer = StereoEnumerator::new(
            &StereoModel::default(),
            vec![StereoSite::Bond(BondId(4)), StereoSite::Atom(AtomId(0))],
        );
        transformer.transform_iter(&source)
    };
    source.atom_mut(AtomId(0)).attributes_mut().charge = NumForm::Lit(1);
    source
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = (StereoKind::Tetrahedral, 1).into();
    drop(source);
    let expected: Vec<_> = [(0, 0), (1, 0), (0, 1), (1, 1)]
        .into_iter()
        .map(|(atom, bond)| {
            let mut result = original.clone();
            result
                .stereo_atom_mut(StereoAtomId(0))
                .attributes_mut()
                .configuration = (StereoKind::Tetrahedral, atom).into();
            result
                .stereo_bond_mut(StereoBondId(0))
                .attributes_mut()
                .configuration = (StereoKind::CisTrans, bond).into();
            result
        })
        .collect();
    let mut first = iter.next().unwrap();
    assert_eq!(first, expected[0]);
    first.atom_mut(AtomId(1)).attributes_mut().charge = NumForm::Lit(-1);
    first
        .stereo_bond_mut(StereoBondId(0))
        .attributes_mut()
        .configuration = (StereoKind::CisTrans, 1).into();
    assert_eq!(iter.by_ref().collect::<Vec<_>>(), expected[1..]);
    assert_eq!(iter.next(), None);
    assert_eq!(iter.next(), None);
}

#[rstest]
fn test_stereo_enumerator_transform_iter_prefix() {
    let mut entries = MoleculeEntries::default();
    for component in 0..128 {
        let base = component * 5;
        entries.atoms.extend(
            [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
                .map(|element| AtomForm::from_element(element).into_concrete()),
        );
        entries.bonds.extend((1..5).map(|ligand| {
            (
                AtomId(base),
                AtomId(base + ligand),
                BondForm::from_order(1).into_concrete(),
            )
        }));
        entries.stereo_atoms.push((
            AtomId(base),
            (1..5)
                .map(|ligand| StereoLigand::new(AtomId(base + ligand), StereoLigandKind::Atom))
                .collect(),
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
        ));
    }
    let source = Molecule::from_entries(entries);
    let transformer = StereoEnumerator::new(
        &StereoModel::default(),
        (0..128)
            .map(|index| StereoSite::Atom(AtomId(index * 5)))
            .collect(),
    );
    let mut first = source.clone();
    for id in first.stereo_atoms().ids().collect::<Vec<_>>() {
        first.stereo_atom_mut(id).attributes_mut().configuration =
            (StereoKind::Tetrahedral, 0).into();
    }
    let mut second = first.clone();
    second
        .stereo_atom_mut(StereoAtomId(127))
        .attributes_mut()
        .configuration = (StereoKind::Tetrahedral, 1).into();
    let mut iter = transformer.transform_iter(&source);
    assert_eq!(iter.next(), Some(first.clone()));
    assert_eq!(iter.next(), Some(second));
    assert_eq!(transformer.transform(source), Ok(first));
}
