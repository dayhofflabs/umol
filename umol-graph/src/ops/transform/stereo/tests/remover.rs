use std::collections::BTreeSet;

use rstest::{fixture, rstest};
use umol_graph_ir::ir::{
    AtomConstraintForm, AtomId, BondId, Constraint, ElementForm, Molecule, NumForm,
    StereoAtomConstraintForm, StereoAtomId, StereoBondConstraintForm, StereoBondId,
    StereoConfigurationForm, StereoCoset, StereoKind, StereoTerm, Stereogenicity,
    StereogenicityForm,
};
use umol_graph_ir::mol_dsl_concrete;

use crate::ops::model::StereoModel;
use crate::ops::stereo::{StereoInconsistency, StereoPerceptionError, StereoSelection, StereoSite};
use crate::ops::transform::{
    RemoveStereoConfig, RemoveStereoError, RemoveStereoMode, StereoRemover, Transformer,
};

#[fixture]
fn mixed_sites() -> Molecule {
    mol_dsl_concrete!(
        r#"{:atoms ["C#T1#D4" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2#C0#R!"] [5 7 "1"] [6 8 "1"]]
        :stereo-atoms [{:site 0 :ligands [2 1 3 4] :attrs "Th0"}]
        :stereo-bonds [{:site 4 :ligands [7 [:h 5] 8 [:h 6]] :attrs "Ct*"}]}"#
    )
}

#[rstest]
#[case::all(StereoSelection::default(), vec![StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(4))], vec![StereoSite::Bond(BondId(4))])]
#[case::atom(StereoSelection { sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(0))])), ..Default::default() }, vec![StereoSite::Atom(AtomId(0))], vec![])]
#[case::bond(StereoSelection { sites: Some(BTreeSet::from([StereoSite::Bond(BondId(4))])), ..Default::default() }, vec![StereoSite::Bond(BondId(4))], vec![StereoSite::Bond(BondId(4))])]
#[case::tetrahedral(StereoSelection { kinds: Some(BTreeSet::from([StereoKind::Tetrahedral])), ..Default::default() }, vec![StereoSite::Atom(AtomId(0))], vec![])]
#[case::cis_trans(StereoSelection { kinds: Some(BTreeSet::from([StereoKind::CisTrans])), ..Default::default() }, vec![StereoSite::Bond(BondId(4))], vec![StereoSite::Bond(BondId(4))])]
#[case::other_kind(StereoSelection { kinds: Some(BTreeSet::from([StereoKind::SquarePlanar])), ..Default::default() }, vec![], vec![])]
#[case::stereogenic(StereoSelection { stereogenicities: Some(BTreeSet::from([Stereogenicity::Stereogenic])), ..Default::default() }, vec![StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(4))], vec![StereoSite::Bond(BondId(4))])]
#[case::symmetric(StereoSelection { stereogenicities: Some(BTreeSet::from([Stereogenicity::Symmetric])), ..Default::default() }, vec![], vec![])]
#[case::intersection(StereoSelection { sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(0))])), kinds: Some(BTreeSet::from([StereoKind::CisTrans])), ..Default::default() }, vec![], vec![])]
#[case::empty_sites(StereoSelection { sites: Some(BTreeSet::new()), ..Default::default() }, vec![], vec![])]
#[case::empty_kinds(StereoSelection { kinds: Some(BTreeSet::new()), ..Default::default() }, vec![], vec![])]
#[case::empty_classes(StereoSelection { stereogenicities: Some(BTreeSet::new()), ..Default::default() }, vec![], vec![])]
fn test_stereo_remover_with_config(
    mixed_sites: Molecule,
    #[case] selection: StereoSelection,
    #[case] all: Vec<StereoSite>,
    #[case] undetermined: Vec<StereoSite>,
    #[values(RemoveStereoMode::All, RemoveStereoMode::Undetermined)] mode: RemoveStereoMode,
) {
    let sites = match mode {
        RemoveStereoMode::All => all,
        RemoveStereoMode::Undetermined => undetermined,
    };
    let expected = match (
        sites.contains(&StereoSite::Atom(AtomId(0))),
        sites.contains(&StereoSite::Bond(BondId(4))),
    ) {
        (true, true) => mol_dsl_concrete!(
            r#"{:atoms ["C#D4" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2#R!"] [5 7 "1"] [6 8 "1"]]}"#
        ),
        (true, false) => mol_dsl_concrete!(
            r#"{:atoms ["C#D4" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2#C0#R!"] [5 7 "1"] [6 8 "1"]]
            :stereo-bonds [{:site 4 :ligands [7 [:h 5] 8 [:h 6]] :attrs "Ct*"}]}"#
        ),
        (false, true) => mol_dsl_concrete!(
            r#"{:atoms ["C#T1#D4" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2#R!"] [5 7 "1"] [6 8 "1"]]
            :stereo-atoms [{:site 0 :ligands [2 1 3 4] :attrs "Th0"}]}"#
        ),
        (false, false) => mixed_sites.clone(),
    };
    let transformer = StereoRemover::with_config(
        &StereoModel::default(),
        RemoveStereoConfig {
            selection,
            mode,
            ..Default::default()
        },
    );
    let original = mixed_sites.clone();
    let mut in_place = mixed_sites.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected);
    assert_eq!(
        transformer.transform_iter(&mixed_sites).collect::<Vec<_>>(),
        vec![expected.clone()]
    );
    assert_eq!(mixed_sites, original);
    assert_eq!(transformer.transform(mixed_sites), Ok(expected.clone()));
    assert_eq!(expected.clone().edit().finish(), Ok(expected));
}

#[rstest]
#[case::tetrahedral(
    mol_dsl_concrete!(r#"{:atoms ["C#T1#D4" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [2 1 3 4] :attrs "Th0"}]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C#D4" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#),
)]
#[case::cis_trans(
    mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "2#C0#R!"] [0 2 "1"] [1 3 "1"]]
        :stereo-bonds [{:site 0 :ligands [[:h 0] 2 3 [:h 1]] :attrs "Ct1"}]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "2#R!"] [0 2 "1"] [1 3 "1"]]}"#),
)]
#[case::kindless(
    mol_dsl_concrete!(r#"{:atoms ["C#T+" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "*"}]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#),
)]
#[case::symbolic(
    mol_dsl_concrete!(r#"{:atoms ["C#T?configuration" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th?configuration"}]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#),
)]
#[case::other_kind(
    mol_dsl_concrete!(r#"{:atoms ["Pt#c+2#T!" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Sp0"}]}"#),
    mol_dsl_concrete!(r#"{:atoms ["Pt#c+2" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#),
)]
#[case::conflicting_evidence(
    mol_dsl_concrete!(r#"{:atoms ["C#T0" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#),
)]
fn test_stereo_remover_transform(#[case] molecule: Molecule, #[case] expected: Molecule) {
    let transformer = StereoRemover::new(&StereoModel::default());
    let original = molecule.clone();
    let mut in_place = molecule.clone();
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
#[case::no_entities(mol_dsl_concrete!(r#"{:atoms ["C#T1" "F" "Cl" "Br" "I"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#))]
fn test_stereo_remover_transform_identity(#[case] molecule: Molecule) {
    let transformer = StereoRemover::new(&StereoModel::default());
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
#[case::open(Some(StereoCoset::Undetermined))]
#[case::first(Some(StereoCoset::Lit(0)))]
#[case::second(Some(StereoCoset::Lit(1)))]
#[case::finite(Some(StereoCoset::lit_set([0, 1])))]
#[case::symbolic(Some(StereoCoset::term(StereoTerm::var("configuration"))))]
#[case::invalid_domain(Some(StereoCoset::Lit(2)))]
#[case::empty_domain(Some(StereoCoset::lit_set([])))]
#[case::kindless(None)]
fn test_stereo_remover_transform_configuration(
    mut mixed_sites: Molecule,
    #[case] coset: Option<StereoCoset>,
    #[values(StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(4)))] site: StereoSite,
) {
    let expected = match site {
        StereoSite::Atom(_) => {
            mixed_sites
                .stereo_atom_mut(StereoAtomId(0))
                .attributes_mut()
                .configuration = coset.map_or(StereoConfigurationForm::Undetermined, |coset| {
                StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, coset)
            });
            mol_dsl_concrete!(
                r#"{:atoms ["C#D4" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2#C0#R!"] [5 7 "1"] [6 8 "1"]]
                :stereo-bonds [{:site 4 :ligands [7 [:h 5] 8 [:h 6]] :attrs "Ct*"}]}"#
            )
        }
        StereoSite::Bond(_) => {
            mixed_sites
                .stereo_bond_mut(StereoBondId(0))
                .attributes_mut()
                .configuration = coset.map_or(StereoConfigurationForm::Undetermined, |coset| {
                StereoConfigurationForm::Kinded(StereoKind::CisTrans, coset)
            });
            mol_dsl_concrete!(
                r#"{:atoms ["C#T1#D4" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2#R!"] [5 7 "1"] [6 8 "1"]]
                :stereo-atoms [{:site 0 :ligands [2 1 3 4] :attrs "Th0"}]}"#
            )
        }
    };
    let transformer = StereoRemover::with_config(
        &StereoModel::default(),
        RemoveStereoConfig {
            selection: StereoSelection {
                sites: Some(BTreeSet::from([site])),
                ..Default::default()
            },
            ..Default::default()
        },
    );
    let original = mixed_sites.clone();
    let mut in_place = mixed_sites.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected);
    assert_eq!(
        transformer.transform_iter(&mixed_sites).collect::<Vec<_>>(),
        vec![expected.clone()]
    );
    assert_eq!(mixed_sites, original);
    assert_eq!(transformer.transform(mixed_sites), Ok(expected));
}

#[rstest]
#[case::literal(Some(StereoCoset::Lit(0)))]
#[case::finite(Some(StereoCoset::lit_set([0, 1])))]
#[case::symbolic(Some(StereoCoset::term(StereoTerm::var("configuration"))))]
#[case::invalid_domain(Some(StereoCoset::Lit(2)))]
#[case::empty_domain(Some(StereoCoset::lit_set([])))]
#[case::kindless(None)]
fn test_stereo_remover_transform_configuration_identity(
    mut mixed_sites: Molecule,
    #[case] coset: Option<StereoCoset>,
    #[values(StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(4)))] site: StereoSite,
) {
    match site {
        StereoSite::Atom(_) => {
            mixed_sites
                .stereo_atom_mut(StereoAtomId(0))
                .attributes_mut()
                .configuration = coset.map_or(StereoConfigurationForm::Undetermined, |coset| {
                StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, coset)
            })
        }
        StereoSite::Bond(_) => {
            mixed_sites
                .stereo_bond_mut(StereoBondId(0))
                .attributes_mut()
                .configuration = coset.map_or(StereoConfigurationForm::Undetermined, |coset| {
                StereoConfigurationForm::Kinded(StereoKind::CisTrans, coset)
            })
        }
    }
    let transformer = StereoRemover::with_config(
        &StereoModel::default(),
        RemoveStereoConfig {
            selection: StereoSelection {
                sites: Some(BTreeSet::from([site])),
                ..Default::default()
            },
            mode: RemoveStereoMode::Undetermined,
            ..Default::default()
        },
    );
    let mut in_place = mixed_sites.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, mixed_sites);
    assert_eq!(
        transformer.transform_iter(&mixed_sites).collect::<Vec<_>>(),
        vec![mixed_sites.clone()]
    );
    assert_eq!(transformer.transform(mixed_sites.clone()), Ok(mixed_sites));
}

#[rstest]
#[case::open_atom(
    mol_dsl_concrete!(r#"{:atoms ["C#T1" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#),
)]
#[case::other_kind(
    mol_dsl_concrete!(r#"{:atoms ["Pt#c+2" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Sp*"}]}"#),
    mol_dsl_concrete!(r#"{:atoms ["Pt#c+2" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#),
)]
fn test_stereo_remover_transform_mode(#[case] molecule: Molecule, #[case] expected: Molecule) {
    let transformer = StereoRemover::with_config(
        &StereoModel::default(),
        RemoveStereoConfig {
            mode: RemoveStereoMode::Undetermined,
            ..Default::default()
        },
    );
    let original = molecule.clone();
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
#[case::stereogenic_atom(
    mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#),
    StereoSite::Atom(AtomId(0)), Stereogenicity::Stereogenic,
)]
#[case::prochiral_atom(
    mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "H" "H"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "H" "H"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#),
    StereoSite::Atom(AtomId(0)), Stereogenicity::Prochiral,
)]
#[case::symmetric_atom(
    mol_dsl_concrete!(r#"{:atoms ["C" "H" "H" "H" "H"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C" "H" "H" "H" "H"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#),
    StereoSite::Atom(AtomId(0)), Stereogenicity::Symmetric,
)]
#[case::stereogenic_bond(
    mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"] :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
        :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct0"}]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"] :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]}"#),
    StereoSite::Bond(BondId(0)), Stereogenicity::Stereogenic,
)]
#[case::symmetric_bond(
    mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "F" "Br" "I"] :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
        :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct0"}]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "F" "Br" "I"] :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]}"#),
    StereoSite::Bond(BondId(0)), Stereogenicity::Symmetric,
)]
fn test_stereo_remover_transform_stereogenicity(
    #[case] mut molecule: Molecule,
    #[case] expected: Molecule,
    #[case] site: StereoSite,
    #[case] class: Stereogenicity,
    #[values(StereoCoset::Undetermined, StereoCoset::Lit(0), StereoCoset::Lit(1), StereoCoset::lit_set([0, 1]))]
    coset: StereoCoset,
) {
    match site {
        StereoSite::Atom(_) => {
            molecule
                .stereo_atom_mut(StereoAtomId(0))
                .attributes_mut()
                .configuration = StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, coset)
        }
        StereoSite::Bond(_) => {
            molecule
                .stereo_bond_mut(StereoBondId(0))
                .attributes_mut()
                .configuration = StereoConfigurationForm::Kinded(StereoKind::CisTrans, coset)
        }
    }
    let transformer = StereoRemover::with_config(
        &StereoModel::default(),
        RemoveStereoConfig {
            selection: StereoSelection {
                stereogenicities: Some(BTreeSet::from([class])),
                ..Default::default()
            },
            ..Default::default()
        },
    );
    let original = molecule.clone();
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
#[case::symmetric_only(mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "H" "H"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]}"#), Stereogenicity::Symmetric)]
#[case::stereogenic_only(mol_dsl_concrete!(r#"{:atoms ["C" "H" "H" "H" "H"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#), Stereogenicity::Stereogenic)]
fn test_stereo_remover_transform_stereogenicity_identity(
    #[case] molecule: Molecule,
    #[case] class: Stereogenicity,
) {
    let transformer = StereoRemover::with_config(
        &StereoModel::default(),
        RemoveStereoConfig {
            selection: StereoSelection {
                stereogenicities: Some(BTreeSet::from([class])),
                ..Default::default()
            },
            ..Default::default()
        },
    );
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
fn test_stereo_remover_transform_context() {
    let molecule = mol_dsl_concrete!(
        r#"{:atoms ["C#T+#D4" "C" "C" "F" "H" "Cl" "Br" "H" "Cl" "Br" "H"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"] [1 6 "1"] [1 7 "1"] [2 8 "1"] [2 9 "1"] [2 10 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"} {:site 1 :ligands [0 5 6 7] :attrs "Th0"} {:site 2 :ligands [0 8 9 10] :attrs "Th1"}]}"#
    );
    let expected = mol_dsl_concrete!(
        r#"{:atoms ["C#D4" "C" "C" "F" "H" "Cl" "Br" "H" "Cl" "Br" "H"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"] [1 6 "1"] [1 7 "1"] [2 8 "1"] [2 9 "1"] [2 10 "1"]]
        :stereo-atoms [{:site 1 :ligands [0 5 6 7] :attrs "Th0"} {:site 2 :ligands [0 8 9 10] :attrs "Th1"}]}"#
    );
    let transformer = StereoRemover::with_config(
        &StereoModel::default(),
        RemoveStereoConfig {
            selection: StereoSelection {
                sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(0))])),
                stereogenicities: Some(BTreeSet::from([Stereogenicity::Stereogenic])),
                ..Default::default()
            },
            mode: RemoveStereoMode::Undetermined,
            ..Default::default()
        },
    );
    let original = molecule.clone();
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
fn test_stereo_remover_transform_constraints() {
    let source = mol_dsl_concrete!(
        r#"{:atoms ["C#T+#D4" "F" "Cl" "Br" "I" "C#T0" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "1"] [5 7 "1"] [5 8 "1"] [5 9 "1"]
            [10 11 "2#C+"] [10 12 "1"] [11 13 "1"] [14 15 "2#C0"] [14 16 "1"] [15 17 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"} {:site 5 :ligands [7 6 8 9] :attrs "Th1"}]
        :stereo-bonds [{:site 8 :ligands [12 [:h 10] 13 [:h 11]] :attrs "Ct*"} {:site 11 :ligands [[:h 14] 16 17 [:h 15]] :attrs "Ct1"}]}"#
    );
    let expected = mol_dsl_concrete!(
        r#"{:atoms ["C#D4" "F" "Cl" "Br" "I" "C#T0" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "1"] [5 7 "1"] [5 8 "1"] [5 9 "1"]
            [10 11 "2"] [10 12 "1"] [11 13 "1"] [14 15 "2#C0"] [14 16 "1"] [15 17 "1"]]
        :stereo-atoms [{:site 5 :ligands [7 6 8 9] :attrs "Th1"}]
        :stereo-bonds [{:site 11 :ligands [[:h 14] 16 17 [:h 15]] :attrs "Ct1"}]}"#
    );
    let mut source_editor = source.edit();
    source_editor.constraints_mut().extend([
        Constraint::StereoAtom(
            StereoAtomId(0),
            StereoKind::Tetrahedral,
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
        ),
        Constraint::StereoAtom(
            StereoAtomId(1),
            StereoKind::Tetrahedral,
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
                Stereogenicity::Stereogenic,
            )),
        ),
        Constraint::StereoBond(
            StereoBondId(0),
            StereoKind::CisTrans,
            StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
        ),
        Constraint::StereoBond(
            StereoBondId(1),
            StereoKind::CisTrans,
            StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
                Stereogenicity::Stereogenic,
            )),
        ),
        Constraint::Atom(AtomId(0), AtomConstraintForm::Degree(NumForm::Lit(4))),
        Constraint::Not(Box::new(Constraint::StereoBond(
            StereoBondId(1),
            StereoKind::CisTrans,
            StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
                Stereogenicity::Symmetric,
            )),
        ))),
        Constraint::And(vec![
            Constraint::StereoAtom(
                StereoAtomId(0),
                StereoKind::Tetrahedral,
                StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
            ),
            Constraint::StereoAtom(
                StereoAtomId(1),
                StereoKind::Tetrahedral,
                StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
            ),
        ]),
    ]);
    let mut molecule = source_editor.finish().unwrap();
    molecule
        .stereo_atom_mut(StereoAtomId(1))
        .attributes_mut()
        .constraints
        .set(StereoAtomConstraintForm::Stereogenicity(
            StereogenicityForm::Lit(Stereogenicity::Stereogenic),
        ));
    molecule
        .stereo_bond_mut(StereoBondId(1))
        .attributes_mut()
        .constraints
        .set(StereoBondConstraintForm::Stereogenicity(
            StereogenicityForm::Lit(Stereogenicity::Stereogenic),
        ));
    let mut expected_editor = expected.edit();
    expected_editor.constraints_mut().extend([
        Constraint::StereoAtom(
            StereoAtomId(0),
            StereoKind::Tetrahedral,
            StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
                Stereogenicity::Stereogenic,
            )),
        ),
        Constraint::StereoBond(
            StereoBondId(0),
            StereoKind::CisTrans,
            StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
                Stereogenicity::Stereogenic,
            )),
        ),
        Constraint::Atom(AtomId(0), AtomConstraintForm::Degree(NumForm::Lit(4))),
        Constraint::Not(Box::new(Constraint::StereoBond(
            StereoBondId(0),
            StereoKind::CisTrans,
            StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
                Stereogenicity::Symmetric,
            )),
        ))),
    ]);
    let mut expected = expected_editor.finish().unwrap();
    expected
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .constraints
        .set(StereoAtomConstraintForm::Stereogenicity(
            StereogenicityForm::Lit(Stereogenicity::Stereogenic),
        ));
    expected
        .stereo_bond_mut(StereoBondId(0))
        .attributes_mut()
        .constraints
        .set(StereoBondConstraintForm::Stereogenicity(
            StereogenicityForm::Lit(Stereogenicity::Stereogenic),
        ));
    let transformer = StereoRemover::with_config(
        &StereoModel::default(),
        RemoveStereoConfig {
            mode: RemoveStereoMode::Undetermined,
            ..Default::default()
        },
    );
    let original = molecule.clone();
    let mut in_place = molecule.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected);
    assert_eq!(
        transformer.transform_iter(&molecule).collect::<Vec<_>>(),
        vec![expected.clone()]
    );
    assert_eq!(molecule, original);
    assert_eq!(transformer.transform(molecule), Ok(expected.clone()));
    assert_eq!(expected.clone().edit().finish(), Ok(expected));
}

#[rstest]
fn test_stereo_remover_transform_partial(mut mixed_sites: Molecule) {
    let mut expected = mol_dsl_concrete!(
        r#"{:atoms ["C#D4" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2#R!"] [5 7 "1"] [6 8 "1"]]}"#
    );
    for input in [&mut mixed_sites, &mut expected] {
        input.atom_mut(AtomId(0)).attributes_mut().element = ElementForm::Undetermined;
        input
            .atom_mut(AtomId(0))
            .attributes_mut()
            .implicit_hydrogens = NumForm::Undetermined;
        input.atom_mut(AtomId(0)).attributes_mut().lone_pairs = NumForm::Undetermined;
        input.bond_mut(BondId(4)).attributes_mut().order = NumForm::Undetermined;
    }
    let transformer = StereoRemover::new(&StereoModel::default());
    let original = mixed_sites.clone();
    let mut in_place = mixed_sites.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected);
    assert_eq!(
        transformer.transform_iter(&mixed_sites).collect::<Vec<_>>(),
        vec![expected.clone()]
    );
    assert_eq!(mixed_sites, original);
    assert_eq!(transformer.transform(mixed_sites), Ok(expected));
}

#[rstest]
#[case::unrestricted(StereoSelection::default())]
#[case::classified(StereoSelection { stereogenicities: Some(BTreeSet::from([Stereogenicity::Stereogenic])), ..Default::default() })]
fn test_stereo_remover_transform_model(#[case] selection: StereoSelection) {
    let molecule = mol_dsl_concrete!(
        r#"{:atoms ["N#n" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 [:lp 0]] :attrs "Th*"}]}"#
    );
    let expected = mol_dsl_concrete!(
        r#"{:atoms ["N#n" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#
    );
    let model = StereoModel {
        kind_models: [None, None, None, None, None, None],
        ..Default::default()
    };
    let transformer = StereoRemover::with_config(
        &model,
        RemoveStereoConfig {
            selection,
            ..Default::default()
        },
    );
    let original = molecule.clone();
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
#[case::invalid_atom(mol_dsl_concrete!(r#"{:atoms ["C#h4"]}"#), StereoSelection { sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(1))])), ..Default::default() },
    StereoPerceptionError::InvalidSite { site: StereoSite::Atom(AtomId(1)) })]
#[case::invalid_bond(mol_dsl_concrete!(r#"{:atoms ["C#h4"]}"#), StereoSelection { sites: Some(BTreeSet::from([StereoSite::Bond(BondId(0))])), ..Default::default() },
    StereoPerceptionError::InvalidSite { site: StereoSite::Bond(BondId(0)) })]
#[case::invalid_empty(mol_dsl_concrete!(r#"{:atoms ["C#h4"]}"#), StereoSelection { sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(1))])), kinds: Some(BTreeSet::new()), ..Default::default() },
    StereoPerceptionError::InvalidSite { site: StereoSite::Atom(AtomId(1)) })]
#[case::kindless(mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "*"}]}"#),
    StereoSelection { stereogenicities: Some(BTreeSet::from([Stereogenicity::Symmetric])), ..Default::default() },
    StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Atom(AtomId(0)) })]
#[case::symbolic(mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th?configuration"}]}"#),
    StereoSelection { stereogenicities: Some(BTreeSet::from([Stereogenicity::Symmetric])), ..Default::default() },
    StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Atom(AtomId(0)) })]
#[case::geometry(mol_dsl_concrete!(r#"{:atoms ["Pt#c+2" "F" "Cl" "Br" "I"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Sp0"}]}"#),
    StereoSelection { stereogenicities: Some(BTreeSet::from([Stereogenicity::Symmetric])), ..Default::default() },
    StereoPerceptionError::UnsupportedKind { kind: StereoKind::SquarePlanar })]
#[case::frame(mol_dsl_concrete!(r#"{:atoms ["C#h" "F" "Cl" "Br"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3] :attrs "Th0"}]}"#),
    StereoSelection { stereogenicities: Some(BTreeSet::from([Stereogenicity::Symmetric])), ..Default::default() },
    StereoPerceptionError::Inconsistency(StereoInconsistency::StereoAtomFailure { stereo_atom: StereoAtomId(0) }))]
#[case::later_bond(mol_dsl_concrete!(r#"{:atoms ["C#T1#D4" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2#C0#R!"] [5 7 "1"] [6 8 "1"]]
    :stereo-atoms [{:site 0 :ligands [2 1 3 4] :attrs "Th0"}]
    :stereo-bonds [{:site 4 :ligands [7 [:h 5] 8 [:h 6]] :attrs "Ct?configuration"}]}"#),
    StereoSelection { stereogenicities: Some(BTreeSet::from([Stereogenicity::Stereogenic])), ..Default::default() },
    StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Bond(BondId(4)) })]
fn test_stereo_remover_transform_error(
    #[case] molecule: Molecule,
    #[case] selection: StereoSelection,
    #[case] cause: StereoPerceptionError,
) {
    let transformer = StereoRemover::with_config(
        &StereoModel::default(),
        RemoveStereoConfig {
            selection,
            ..Default::default()
        },
    );
    let error = RemoveStereoError::Perception(cause);
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
#[case::site(StereoSelection { sites: Some(BTreeSet::from([StereoSite::Bond(BondId(4))])), stereogenicities: Some(BTreeSet::from([Stereogenicity::Stereogenic])), ..Default::default() })]
#[case::kind(StereoSelection { kinds: Some(BTreeSet::from([StereoKind::CisTrans])), stereogenicities: Some(BTreeSet::from([Stereogenicity::Stereogenic])), ..Default::default() })]
fn test_stereo_remover_transform_scope(
    mut mixed_sites: Molecule,
    #[case] selection: StereoSelection,
) {
    mixed_sites
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = StereoConfigurationForm::Undetermined;
    let expected = mol_dsl_concrete!(
        r#"{:atoms ["C#T1#D4" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2#R!"] [5 7 "1"] [6 8 "1"]]
        :stereo-atoms [{:site 0 :ligands [2 1 3 4] :attrs "*"}]}"#
    );
    let transformer = StereoRemover::with_config(
        &StereoModel::default(),
        RemoveStereoConfig {
            selection,
            ..Default::default()
        },
    );
    let original = mixed_sites.clone();
    let mut in_place = mixed_sites.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected);
    assert_eq!(
        transformer.transform_iter(&mixed_sites).collect::<Vec<_>>(),
        vec![expected.clone()]
    );
    assert_eq!(mixed_sites, original);
    assert_eq!(transformer.transform(mixed_sites), Ok(expected));
}

#[rstest]
#[case::literal(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)))]
#[case::symbolic(StereoConfigurationForm::Kinded(
    StereoKind::Tetrahedral,
    StereoCoset::term(StereoTerm::var("configuration"))
))]
#[case::kindless(StereoConfigurationForm::Undetermined)]
#[case::geometry(StereoConfigurationForm::Kinded(StereoKind::SquarePlanar, StereoCoset::Lit(0)))]
fn test_stereo_remover_transform_mode_scope(
    mut mixed_sites: Molecule,
    #[case] configuration: StereoConfigurationForm,
) {
    mixed_sites
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = configuration.clone();
    let mut expected = mol_dsl_concrete!(
        r#"{:atoms ["C#T1#D4" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2#R!"] [5 7 "1"] [6 8 "1"]]
        :stereo-atoms [{:site 0 :ligands [2 1 3 4] :attrs "Th0"}]}"#
    );
    expected
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = configuration;
    let transformer = StereoRemover::with_config(
        &StereoModel::default(),
        RemoveStereoConfig {
            selection: StereoSelection {
                stereogenicities: Some(BTreeSet::from([Stereogenicity::Stereogenic])),
                ..Default::default()
            },
            mode: RemoveStereoMode::Undetermined,
            ..Default::default()
        },
    );
    let original = mixed_sites.clone();
    let mut in_place = mixed_sites.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected);
    assert_eq!(
        transformer.transform_iter(&mixed_sites).collect::<Vec<_>>(),
        vec![expected.clone()]
    );
    assert_eq!(mixed_sites, original);
    assert_eq!(transformer.transform(mixed_sites), Ok(expected));
}

#[rstest]
fn test_stereo_remover_transform_iter(mixed_sites: Molecule) {
    let transformer = StereoRemover::with_config(
        &StereoModel::default(),
        RemoveStereoConfig {
            selection: StereoSelection {
                sites: Some(BTreeSet::from([StereoSite::Bond(BondId(4))])),
                ..Default::default()
            },
            ..Default::default()
        },
    );
    let original = mixed_sites.clone();
    let expected = mol_dsl_concrete!(
        r#"{:atoms ["C#T1#D4" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2#R!"] [5 7 "1"] [6 8 "1"]]
        :stereo-atoms [{:site 0 :ligands [2 1 3 4] :attrs "Th0"}]}"#
    );
    let mut first = transformer.transform_iter(&mixed_sites).next().unwrap();
    let second = transformer.transform_iter(&mixed_sites).collect::<Vec<_>>();
    first.atom_mut(AtomId(0)).attributes_mut().charge = NumForm::Lit(1);
    first
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration =
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined);
    let changed = mol_dsl_concrete!(
        r#"{:atoms ["C#c+#T1#D4" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2#R!"] [5 7 "1"] [6 8 "1"]]
        :stereo-atoms [{:site 0 :ligands [2 1 3 4] :attrs "Th*"}]}"#
    );
    assert_eq!(first, changed);
    assert_eq!(second, vec![expected]);
    assert_eq!(mixed_sites, original);
}
