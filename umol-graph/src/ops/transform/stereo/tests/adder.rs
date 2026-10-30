use std::collections::BTreeSet;

use rstest::{fixture, rstest};
use umol_graph_ir::ir::{
    AtomConstraintForm, AtomFieldKind, AtomId, BondId, ElementForm, Molecule, NumForm,
    StereoAtomForm, StereoAtomId, StereoBondForm, StereoBondId, StereoConfigurationForm,
    StereoCoset, StereoKind, StereoLigand, StereoLigandKind, StereoTerm, Stereogenicity,
    TetrahedralStereoForm,
};
use umol_graph_ir::mol_dsl_concrete;

use crate::ops::model::{ElementScope, StereoKindModel, StereoModel};
use crate::ops::stereo::{StereoInconsistency, StereoPerceptionError, StereoSelection, StereoSite};
use crate::ops::transform::{AddStereoConfig, AddStereoError, StereoAdder, Transformer};

#[fixture]
fn mixed_sites() -> Molecule {
    mol_dsl_concrete!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2"] [5 7 "1"] [6 8 "1"]]}"#
    )
}

#[rstest]
#[case::atom_disabled(StereoKind::Tetrahedral, None, StereoSite::Bond(BondId(4)))]
#[case::bond_disabled(StereoKind::CisTrans, None, StereoSite::Atom(AtomId(0)))]
#[case::atom_scope(StereoKind::Tetrahedral, Some(StereoKindModel { scope: ElementScope::AllowList(vec![]), fluxionality: false }), StereoSite::Bond(BondId(4)))]
#[case::bond_scope(StereoKind::CisTrans, Some(StereoKindModel { scope: ElementScope::AllowList(vec![]), fluxionality: false }), StereoSite::Atom(AtomId(0)))]
fn test_stereo_adder_new(
    mixed_sites: Molecule,
    #[case] kind: StereoKind,
    #[case] kind_model: Option<StereoKindModel>,
    #[case] site: StereoSite,
) {
    let mut model = StereoModel::default();
    model.kind_models[kind as usize] = kind_model;
    let transformer = StereoAdder::new(&model);
    let mut editor = mixed_sites.clone().edit();
    match site {
        StereoSite::Atom(atom) => {
            editor.add_stereo_atom(
                atom,
                &[1, 2, 3, 4].map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
            );
        }
        StereoSite::Bond(bond) => {
            editor.add_stereo_bond(
                bond,
                &[
                    StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(8), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(6), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined),
            );
        }
    }
    assert_eq!(
        transformer.transform(mixed_sites),
        Ok(editor.finish().unwrap())
    );
}

#[rstest]
#[case::default(AddStereoConfig::default().selection, vec![StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(4))])]
#[case::atom(StereoSelection { sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(0))])), ..AddStereoConfig::default().selection }, vec![StereoSite::Atom(AtomId(0))])]
#[case::bond(StereoSelection { sites: Some(BTreeSet::from([StereoSite::Bond(BondId(4))])), ..AddStereoConfig::default().selection }, vec![StereoSite::Bond(BondId(4))])]
#[case::tetrahedral(StereoSelection { kinds: Some(BTreeSet::from([StereoKind::Tetrahedral])), ..AddStereoConfig::default().selection }, vec![StereoSite::Atom(AtomId(0))])]
#[case::cis_trans(StereoSelection { kinds: Some(BTreeSet::from([StereoKind::CisTrans])), ..AddStereoConfig::default().selection }, vec![StereoSite::Bond(BondId(4))])]
#[case::intersection(StereoSelection { sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(0))])), kinds: Some(BTreeSet::from([StereoKind::CisTrans])), ..AddStereoConfig::default().selection }, vec![])]
#[case::unrestricted(StereoSelection::default(), vec![StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(4))])]
#[case::prochiral(StereoSelection { stereogenicities: Some(BTreeSet::from([Stereogenicity::Prochiral])), ..AddStereoConfig::default().selection }, vec![])]
#[case::empty_sites(StereoSelection { sites: Some(BTreeSet::new()), ..AddStereoConfig::default().selection }, vec![])]
#[case::empty_kinds(StereoSelection { kinds: Some(BTreeSet::new()), ..AddStereoConfig::default().selection }, vec![])]
#[case::empty_classes(StereoSelection { stereogenicities: Some(BTreeSet::new()), ..AddStereoConfig::default().selection }, vec![])]
fn test_stereo_adder_with_config(
    mixed_sites: Molecule,
    #[case] selection: StereoSelection,
    #[case] sites: Vec<StereoSite>,
) {
    let transformer = StereoAdder::with_config(
        &StereoModel::default(),
        AddStereoConfig {
            selection,
            ..Default::default()
        },
    );
    let original = mixed_sites.clone();
    let mut editor = mixed_sites.clone().edit();
    for site in sites {
        match site {
            StereoSite::Atom(atom) => {
                editor.add_stereo_atom(
                    atom,
                    &[1, 2, 3, 4].map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)),
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
                );
            }
            StereoSite::Bond(bond) => {
                editor.add_stereo_bond(
                    bond,
                    &[
                        StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(8), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(6), StereoLigandKind::ImplicitHydrogen),
                    ],
                    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined),
                );
            }
        }
    }
    let expected = editor.finish().unwrap();
    let mut in_place = mixed_sites.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected);
    assert_eq!(
        transformer.transform_iter(&mixed_sites).collect::<Vec<_>>(),
        vec![expected.clone()]
    );
    assert_eq!(mixed_sites, original);
    assert_eq!(transformer.transform(mixed_sites), Ok(expected.clone()));
    assert_eq!(transformer.transform(expected.clone()), Ok(expected));
}

#[rstest]
#[case::tetrahedral(
    mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]}"#),
)]
#[case::implicit_hydrogen(
    mol_dsl_concrete!(r#"{:atoms ["C#h" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C#h" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th*"}]}"#),
)]
#[case::lone_pair(
    mol_dsl_concrete!(r#"{:atoms ["P#n" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#),
    mol_dsl_concrete!(r#"{:atoms ["P#n" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 [:lp 0]] :attrs "Th*"}]}"#),
)]
#[case::cis_trans(
    mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
        :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct*"}]}"#),
)]
#[case::imine(
    mol_dsl_concrete!(r#"{:atoms ["C#h" "N#n" "F" "C#h3"] :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C#h" "N#n" "F" "C#h3"] :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
        :stereo-bonds [{:site 0 :ligands [2 [:h 0] 3 [:lp 1]] :attrs "Ct*"}]}"#),
)]
#[case::source_evidence(
    mol_dsl_concrete!(r#"{:atoms ["C#T1#D4" "C#h" "F" "Cl" "Br" "C#h" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "2#C0"] [5 6 "1"]]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C#T1#D4" "C#h" "F" "Cl" "Br" "C#h" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "2#C0"] [5 6 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]
        :stereo-bonds [{:site 4 :ligands [0 [:h 1] 6 [:h 5]] :attrs "Ct0"}]}"#),
)]
#[case::retained_frame(
    mol_dsl_concrete!(r#"{:atoms ["C#T0" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2"] [5 7 "1"] [6 8 "1"]]
        :stereo-atoms [{:site 0 :ligands [2 1 3 4] :attrs "Th1"}]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C#T0" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2"] [5 7 "1"] [6 8 "1"]]
        :stereo-atoms [{:site 0 :ligands [2 1 3 4] :attrs "Th1"}]
        :stereo-bonds [{:site 4 :ligands [7 [:h 5] 8 [:h 6]] :attrs "Ct*"}]}"#),
)]
#[case::retained_geometry(
    mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I" "Pt#c+2" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "1"] [5 7 "1"] [5 8 "1"] [5 9 "1"]]
        :stereo-atoms [{:site 5 :ligands [6 7 8 9] :attrs "Sp0"}]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I" "Pt#c+2" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "1"] [5 7 "1"] [5 8 "1"] [5 9 "1"]]
        :stereo-atoms [{:site 5 :ligands [6 7 8 9] :attrs "Sp0"} {:site 0 :ligands [1 2 3 4] :attrs "Th*"}]}"#),
)]
#[case::remote_stereo(
    mol_dsl_concrete!(r#"{:atoms ["C" "C#T0" "C#T1" "F" "H" "Cl" "Br" "H" "Cl" "Br" "H"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"] [1 6 "1"] [1 7 "1"] [2 8 "1"] [2 9 "1"] [2 10 "1"]]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C" "C#T0" "C#T1" "F" "H" "Cl" "Br" "H" "Cl" "Br" "H"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"] [1 6 "1"] [1 7 "1"] [2 8 "1"] [2 9 "1"] [2 10 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"} {:site 1 :ligands [0 5 6 7] :attrs "Th0"} {:site 2 :ligands [0 8 9 10] :attrs "Th1"}]}"#),
)]
fn test_stereo_adder_transform(#[case] molecule: Molecule, #[case] expected: Molecule) {
    let transformer = StereoAdder::new(&StereoModel::default());
    let original = molecule.clone();
    let mut in_place = molecule.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected);
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected);
    assert_eq!(
        transformer.transform_iter(&molecule).collect::<Vec<_>>(),
        vec![expected.clone()]
    );
    assert_eq!(molecule, original);
    assert_eq!(transformer.transform(molecule), Ok(expected.clone()));
    assert_eq!(
        transformer.transform(expected.clone()),
        Ok(expected.clone())
    );
    assert_eq!(expected.clone().edit().finish(), Ok(expected));
}

#[rstest]
#[case::empty(Molecule::default())]
#[case::acyclic(mol_dsl_concrete!(r#"{:atoms ["C#h4"]}"#))]
#[case::explicit_hydrogens(mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "H" "H"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#))]
#[case::virtual_hydrogens(mol_dsl_concrete!(r#"{:atoms ["C#h2" "F" "Cl"] :bonds [[0 1 "1"] [0 2 "1"]]}"#))]
#[case::negative_atom(mol_dsl_concrete!(r#"{:atoms ["C#T!" "F" "Cl" "Br" "I"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#))]
#[case::negative_bond(mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "F" "Cl"]
    :bonds [[0 1 "2#C!"] [0 2 "1"] [1 3 "1"]]}"#))]
#[case::symmetric_bond(mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "F" "Br" "I"]
    :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]}"#))]
#[case::existing_atom(mol_dsl_concrete!(r#"{:atoms ["C#T0" "F" "Cl" "Br" "I"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [2 1 3 4] :attrs "Th1"}]}"#))]
#[case::existing_kind(mol_dsl_concrete!(r#"{:atoms ["Pt#c+2" "F" "Cl" "Br" "I"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Sp0"}]}"#))]
#[case::existing_bond(mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "F" "Cl"]
    :bonds [[0 1 "2#C0"] [0 2 "1"] [1 3 "1"]]
    :stereo-bonds [{:site 0 :ligands [[:h 0] 2 3 [:h 1]] :attrs "Ct1"}]}"#))]
#[case::existing_domain(mol_dsl_concrete!(r#"{:atoms ["C#T1" "F" "Cl" "Br" "I"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th{0,1}"}]}"#))]
#[case::existing_symbol(mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th?configuration"}]}"#))]
#[case::unknown_kind(mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "*"}]}"#))]
#[case::small_ring(mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h2" "C#h2" "C#h2" "C#h2"]
    :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#))]
fn test_stereo_adder_transform_identity(#[case] molecule: Molecule) {
    let transformer = StereoAdder::new(&StereoModel::default());
    let mut in_place = molecule.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, molecule);
    assert_eq!(
        transformer.transform_iter(&molecule).collect::<Vec<_>>(),
        vec![molecule.clone()]
    );
    assert_eq!(transformer.transform(molecule.clone()), Ok(molecule));
}

#[rstest]
#[case::prochiral(mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "H" "H"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#), StereoSite::Atom(AtomId(0)), Stereogenicity::Prochiral)]
#[case::symmetric(mol_dsl_concrete!(r#"{:atoms ["C" "H" "H" "H" "H"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#), StereoSite::Atom(AtomId(0)), Stereogenicity::Symmetric)]
#[case::symmetric_bond(mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "F" "Br" "I"]
    :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]}"#), StereoSite::Bond(BondId(0)), Stereogenicity::Symmetric)]
fn test_stereo_adder_transform_selection(
    #[case] molecule: Molecule,
    #[case] site: StereoSite,
    #[case] class: Stereogenicity,
    #[values(false, true)] unrestricted: bool,
) {
    let transformer = StereoAdder::with_config(
        &StereoModel::default(),
        AddStereoConfig {
            selection: StereoSelection {
                stereogenicities: (!unrestricted).then_some(BTreeSet::from([class])),
                ..AddStereoConfig::default().selection
            },
            ..Default::default()
        },
    );
    let mut editor = molecule.clone().edit();
    match site {
        StereoSite::Atom(atom) => {
            editor.add_stereo_atom(
                atom,
                &[1, 2, 3, 4].map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
            );
        }
        StereoSite::Bond(bond) => {
            editor.add_stereo_bond(
                bond,
                &[2, 3, 4, 5].map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)),
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined),
            );
        }
    }
    let expected = editor.finish().unwrap();
    let mut in_place = molecule.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected);
    assert_eq!(
        transformer.transform_iter(&molecule).collect::<Vec<_>>(),
        vec![expected.clone()]
    );
    assert_eq!(transformer.transform(molecule), Ok(expected.clone()));
    assert_eq!(transformer.transform(expected.clone()), Ok(expected));
}

#[rstest]
#[case::open(StereoCoset::Undetermined)]
#[case::first(StereoCoset::Lit(0))]
#[case::second(StereoCoset::Lit(1))]
#[case::finite(StereoCoset::lit_set([0, 1]))]
#[case::symbolic(StereoCoset::term(StereoTerm::var("configuration")))]
fn test_stereo_adder_transform_configuration(#[case] coset: StereoCoset) {
    let mut molecule = mol_dsl_concrete!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#
    );
    molecule.atom_mut(AtomId(0)).attributes_mut().constraints =
        AtomConstraintForm::tetrahedral_stereo(TetrahedralStereoForm::Stereo(coset.clone())).into();
    let original = molecule.clone();
    let mut editor = molecule.clone().edit();
    editor.add_stereo_atom(
        AtomId(0),
        &[1, 2, 3, 4].map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)),
        StereoAtomForm::new(StereoKind::Tetrahedral, coset),
    );
    let expected = editor.finish().unwrap();
    let transformer = StereoAdder::with_config(
        &StereoModel::default(),
        AddStereoConfig {
            selection: StereoSelection {
                stereogenicities: None,
                ..AddStereoConfig::default().selection
            },
            ..Default::default()
        },
    );
    let mut in_place = molecule.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected);
    assert_eq!(
        transformer.transform_iter(&molecule).collect::<Vec<_>>(),
        vec![expected.clone()]
    );
    assert_eq!(molecule, original);
    assert_eq!(transformer.transform(molecule), Ok(expected.clone()));
    assert_eq!(transformer.transform(expected.clone()), Ok(expected));
}

#[rstest]
fn test_stereo_adder_transform_partial() {
    let mut molecule = mol_dsl_concrete!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I" "C"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#
    );
    let mut expected = mol_dsl_concrete!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I" "C"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]}"#
    );
    for input in [&mut molecule, &mut expected] {
        input.atom_mut(AtomId(5)).attributes_mut().element = ElementForm::Undetermined;
        input
            .atom_mut(AtomId(5))
            .attributes_mut()
            .implicit_hydrogens = NumForm::Undetermined;
        input.atom_mut(AtomId(5)).attributes_mut().lone_pairs = NumForm::Undetermined;
    }
    let original = molecule.clone();
    let transformer = StereoAdder::new(&StereoModel::default());
    let mut in_place = molecule.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected);
    assert_eq!(
        transformer.transform_iter(&molecule).collect::<Vec<_>>(),
        vec![expected.clone()]
    );
    assert_eq!(molecule, original);
    assert_eq!(transformer.transform(molecule), Ok(expected));
}

#[rstest]
#[case::invalid_atom(mol_dsl_concrete!(r#"{:atoms ["C#h4"]}"#),
    StereoSelection { sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(1))])), ..AddStereoConfig::default().selection },
    AddStereoError::Perception(StereoPerceptionError::InvalidSite { site: StereoSite::Atom(AtomId(1)) }))]
#[case::invalid_bond(mol_dsl_concrete!(r#"{:atoms ["C#h4"]}"#),
    StereoSelection { sites: Some(BTreeSet::from([StereoSite::Bond(BondId(0))])), ..AddStereoConfig::default().selection },
    AddStereoError::Perception(StereoPerceptionError::InvalidSite { site: StereoSite::Bond(BondId(0)) }))]
#[case::unsupported(mol_dsl_concrete!(r#"{:atoms ["C#h4"]}"#),
    StereoSelection { kinds: Some(BTreeSet::from([StereoKind::Axial])), ..AddStereoConfig::default().selection },
    AddStereoError::Perception(StereoPerceptionError::UnsupportedKind { kind: StereoKind::Axial }))]
#[case::atom_evidence(mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I" "C#h2#T+" "F" "Cl"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "1"] [5 7 "1"]]}"#),
    AddStereoConfig::default().selection,
    AddStereoError::Inconsistency(StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(5) }))]
#[case::bond_evidence(mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I" "C#h2" "C#h" "Cl"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2#C1"] [6 7 "1"]]}"#),
    AddStereoConfig::default().selection,
    AddStereoError::Inconsistency(StereoInconsistency::CisTransStereoFailure { bond: BondId(4) }))]
#[case::atom_mismatch(mol_dsl_concrete!(r#"{:atoms ["C#T0" "F" "Cl" "Br" "I"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#),
    AddStereoConfig::default().selection,
    AddStereoError::Inconsistency(StereoInconsistency::TetrahedralStereoMismatch { atom: AtomId(0), stereo_atom: StereoAtomId(0) }))]
#[case::competing_kind(mol_dsl_concrete!(r#"{:atoms ["C#T0" "F" "Cl" "Br" "I"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Sp0"}]}"#),
    AddStereoConfig::default().selection,
    AddStereoError::Inconsistency(StereoInconsistency::TetrahedralStereoMismatch { atom: AtomId(0), stereo_atom: StereoAtomId(0) }))]
#[case::bond_mismatch(mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "F" "Cl"]
    :bonds [[0 1 "2#C0"] [0 2 "1"] [1 3 "1"]]
    :stereo-bonds [{:site 0 :ligands [2 [:h 0] 3 [:h 1]] :attrs "Ct1"}]}"#),
    AddStereoConfig::default().selection,
    AddStereoError::Inconsistency(StereoInconsistency::CisTransStereoMismatch { bond: BondId(0), stereo_bond: StereoBondId(0) }))]
#[case::symbolic_configuration(mol_dsl_concrete!(r#"{:atoms ["C#T?configuration" "F" "Cl" "Br" "I"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#),
    AddStereoConfig::default().selection,
    AddStereoError::Perception(StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Atom(AtomId(0)) }))]
fn test_stereo_adder_transform_error(
    #[case] molecule: Molecule,
    #[case] selection: StereoSelection,
    #[case] error: AddStereoError,
) {
    let transformer = StereoAdder::with_config(
        &StereoModel::default(),
        AddStereoConfig {
            selection,
            ..Default::default()
        },
    );
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
#[case::element(AtomFieldKind::Element)]
#[case::charge(AtomFieldKind::Charge)]
#[case::hydrogens(AtomFieldKind::ImplicitHydrogens)]
#[case::lone_pairs(AtomFieldKind::LonePairs)]
fn test_stereo_adder_transform_scope(#[case] field: AtomFieldKind) {
    let mut molecule = mol_dsl_concrete!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I" "C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "1"] [5 7 "1"] [5 8 "1"] [5 9 "1"]]}"#
    );
    let mut atom = molecule.atom_mut(AtomId(5));
    let attributes = atom.attributes_mut();
    match field {
        AtomFieldKind::Element => attributes.element = ElementForm::Undetermined,
        AtomFieldKind::Charge => attributes.charge = NumForm::Undetermined,
        AtomFieldKind::ImplicitHydrogens => attributes.implicit_hydrogens = NumForm::Undetermined,
        AtomFieldKind::LonePairs => attributes.lone_pairs = NumForm::Undetermined,
        _ => unreachable!(),
    }
    let original = molecule.clone();
    let error = AddStereoError::Perception(StereoPerceptionError::UnderdeterminedAtomField {
        atom: AtomId(5),
        field,
    });
    let unrestricted = StereoAdder::new(&StereoModel::default());
    let mut in_place = molecule.clone();
    assert_eq!(
        unrestricted.transform_into(&mut in_place),
        Err(error.clone())
    );
    assert_eq!(in_place, original);
    assert_eq!(
        unrestricted.transform_iter(&molecule).collect::<Vec<_>>(),
        Vec::<Molecule>::new()
    );
    assert_eq!(unrestricted.transform(molecule.clone()), Err(error));
    let restricted = StereoAdder::with_config(
        &StereoModel::default(),
        AddStereoConfig {
            selection: StereoSelection {
                sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(0))])),
                ..AddStereoConfig::default().selection
            },
            ..Default::default()
        },
    );
    let mut editor = molecule.clone().edit();
    editor.add_stereo_atom(
        AtomId(0),
        &[1, 2, 3, 4].map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)),
        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
    );
    let expected = editor.finish().unwrap();
    assert_eq!(restricted.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected);
    assert_eq!(
        restricted.transform_iter(&molecule).collect::<Vec<_>>(),
        vec![expected.clone()]
    );
    assert_eq!(molecule, original);
    assert_eq!(restricted.transform(molecule), Ok(expected));
}

#[rstest]
fn test_stereo_adder_transform_order(mut mixed_sites: Molecule) {
    mixed_sites.bond_mut(BondId(4)).attributes_mut().order = NumForm::Undetermined;
    let original = mixed_sites.clone();
    let error = AddStereoError::Perception(StereoPerceptionError::UnderdeterminedBondOrder {
        bond: BondId(4),
    });
    let transformer = StereoAdder::new(&StereoModel::default());
    let mut in_place = mixed_sites.clone();
    assert_eq!(
        transformer.transform_into(&mut in_place),
        Err(error.clone())
    );
    assert_eq!(in_place, original);
    assert_eq!(
        transformer.transform_iter(&mixed_sites).collect::<Vec<_>>(),
        Vec::<Molecule>::new()
    );
    assert_eq!(transformer.transform(mixed_sites.clone()), Err(error));
    assert_eq!(mixed_sites, original);
}

#[rstest]
#[case::atom(StereoKind::Tetrahedral)]
#[case::bond(StereoKind::CisTrans)]
fn test_stereo_adder_transform_model(mixed_sites: Molecule, #[case] kind: StereoKind) {
    let mut model = StereoModel::default();
    model.kind_models[kind as usize]
        .as_mut()
        .unwrap()
        .fluxionality = true;
    let error = AddStereoError::Perception(StereoPerceptionError::UnsupportedFluxionality { kind });
    let transformer = StereoAdder::new(&model);
    let original = mixed_sites.clone();
    let mut in_place = mixed_sites.clone();
    assert_eq!(
        transformer.transform_into(&mut in_place),
        Err(error.clone())
    );
    assert_eq!(in_place, original);
    assert_eq!(
        transformer.transform_iter(&mixed_sites).collect::<Vec<_>>(),
        Vec::<Molecule>::new()
    );
    assert_eq!(mixed_sites, original);
    assert_eq!(transformer.transform(mixed_sites), Err(error));
}

#[rstest]
fn test_stereo_adder_transform_iter(mixed_sites: Molecule) {
    let transformer = StereoAdder::new(&StereoModel::default());
    let original = mixed_sites.clone();
    let expected = mol_dsl_concrete!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2"] [5 7 "1"] [6 8 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]
        :stereo-bonds [{:site 4 :ligands [7 [:h 5] 8 [:h 6]] :attrs "Ct*"}]}"#
    );
    let mut first = transformer.transform_iter(&mixed_sites).next().unwrap();
    let second = transformer.transform_iter(&mixed_sites).collect::<Vec<_>>();
    first.atom_mut(AtomId(0)).attributes_mut().charge = NumForm::Lit(1);
    first
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration =
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1));
    let changed = mol_dsl_concrete!(
        r#"{:atoms ["C#c+" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2"] [5 7 "1"] [6 8 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]
        :stereo-bonds [{:site 4 :ligands [7 [:h 5] 8 [:h 6]] :attrs "Ct*"}]}"#
    );
    assert_eq!(first, changed);
    assert_eq!(second, vec![expected]);
    assert_eq!(mixed_sites, original);
}
