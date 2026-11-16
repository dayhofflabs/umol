use std::ops::Range;

use rstest::*;
use umol_chem::element::Element;
use umol_perm::{Orientation, Permutation};

use super::{molecule, virtual_ligand_molecule};
use crate::ir::{
    AtomForm, AtomId, BondForm, GraphSymmetryConfig, LigandPermutation, Molecule, MoleculeEntries,
    OrientedLigandPermutation, StereoAtomConstraintForm, StereoAtomForm, StereoAtomId,
    StereoBondId, StereoConfigurationForm, StereoCoset, StereoKind, StereoLigand, StereoLigandKind,
    StereoLigandPair, StereoLigandPosition, StereoTerm, Stereogenicity, Topicity, TopicityForm,
    TopicityRelationForm,
};
use crate::mol_dsl_concrete;

#[fixture]
fn dichlorofluoromethane() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C#h" "F" "Cl" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]
    }"#
    )
}

#[fixture]
fn difluorodichloroethene() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["C" "C" "F" "Cl" "F" "Cl"]
        :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
        :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct0"}]
    }"#
    )
}

#[rstest]
fn test_stereo_atom_view_with_symmetry_context(dichlorofluoromethane: Molecule) {
    let symmetry = Molecule::default()
        .graph_symmetry(&GraphSymmetryConfig::default())
        .into_graph_symmetry();
    let stereo = dichlorofluoromethane
        .stereo_atom(StereoAtomId(0))
        .with_symmetry(&symmetry);
    let _ = stereo.stereogenicity();
    let _ = stereo.topicity(StereoLigandPosition(1), StereoLigandPosition(2));
    let _ = stereo.ligand_symmetry(OrientedLigandPermutation {
        permutation: LigandPermutation(Permutation::identity(4)),
        orientation: Orientation::Proper,
    });
}

#[rstest]
#[case::dichlorodifluoroethane_chiral(mol_dsl_concrete!(r#"{
        :atoms ["C" "C" "H" "F" "Cl" "H" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"] [1 6 "1"] [1 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                       {:site 1 :ligands [0 5 6 7] :attrs "Th0"}]
    }"#), vec![Stereogenicity::Stereogenic, Stereogenicity::Stereogenic])]
#[case::dichlorodifluoroethane_meso(mol_dsl_concrete!(r#"{
        :atoms ["C" "C" "H" "F" "Cl" "H" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"] [1 6 "1"] [1 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                       {:site 1 :ligands [0 5 6 7] :attrs "Th1"}]
    }"#), vec![Stereogenicity::Stereogenic, Stereogenicity::Stereogenic])]
#[case::dichlorodifluoroethane_chiral_renumbered(mol_dsl_concrete!(r#"{
        :atoms ["C" "C" "H" "F" "Cl" "F" "H" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"] [1 6 "1"] [1 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                       {:site 1 :ligands [0 5 6 7] :attrs "Th1"}]
    }"#), vec![Stereogenicity::Stereogenic, Stereogenicity::Stereogenic])]
#[case::dichlorodifluoroethane_meso_renumbered(mol_dsl_concrete!(r#"{
        :atoms ["C" "C" "H" "F" "Cl" "F" "H" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"] [1 6 "1"] [1 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                       {:site 1 :ligands [0 5 6 7] :attrs "Th0"}]
    }"#), vec![Stereogenicity::Stereogenic, Stereogenicity::Stereogenic])]
#[case::trichloropentane_meso(mol_dsl_concrete!(r#"{
        :atoms ["C#h3" "C" "C" "C" "C#h3" "Cl" "Cl" "Cl" "H" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [1 5 "1"] [2 6 "1"] [3 7 "1"]
                [1 8 "1"] [2 9 "1"] [3 10 "1"]]
        :stereo-atoms [{:site 1 :ligands [0 2 5 8] :attrs "Th0"}
                       {:site 2 :ligands [1 3 6 9] :attrs "Th0"}
                       {:site 3 :ligands [2 4 7 10] :attrs "Th0"}]
    }"#), vec![Stereogenicity::Stereogenic, Stereogenicity::Stereogenic, Stereogenicity::Stereogenic])]
#[case::trichloropentane_meso_inverted_c3(mol_dsl_concrete!(r#"{
        :atoms ["C#h3" "C" "C" "C" "C#h3" "Cl" "Cl" "Cl" "H" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [1 5 "1"] [2 6 "1"] [3 7 "1"]
                [1 8 "1"] [2 9 "1"] [3 10 "1"]]
        :stereo-atoms [{:site 1 :ligands [0 2 5 8] :attrs "Th0"}
                       {:site 2 :ligands [1 3 6 9] :attrs "Th1"}
                       {:site 3 :ligands [2 4 7 10] :attrs "Th0"}]
    }"#), vec![Stereogenicity::Stereogenic, Stereogenicity::Stereogenic, Stereogenicity::Stereogenic])]
#[case::trichloropentane_chiral(mol_dsl_concrete!(r#"{
        :atoms ["C#h3" "C" "C" "C" "C#h3" "Cl" "Cl" "Cl" "H" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [1 5 "1"] [2 6 "1"] [3 7 "1"]
                [1 8 "1"] [2 9 "1"] [3 10 "1"]]
        :stereo-atoms [{:site 1 :ligands [0 2 5 8] :attrs "Th0"}
                       {:site 2 :ligands [1 3 6 9] :attrs "Th0"}
                       {:site 3 :ligands [2 4 7 10] :attrs "Th1"}]
    }"#), vec![Stereogenicity::Stereogenic, Stereogenicity::Prostereogenic, Stereogenicity::Stereogenic])]
#[case::dichlorocyclobutane_cis(mol_dsl_concrete!(r#"{
        :atoms ["C" "C#h2" "C" "C#h2" "Cl" "Cl" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 0 "1"] [0 4 "1"] [2 5 "1"] [0 6 "1"] [2 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 3 4 6] :attrs "Th0"}
                       {:site 2 :ligands [1 3 5 7] :attrs "Th1"}]
    }"#), vec![Stereogenicity::Stereogenic, Stereogenicity::Stereogenic])]
#[case::dichlorocyclobutane_trans(mol_dsl_concrete!(r#"{
        :atoms ["C" "C#h2" "C" "C#h2" "Cl" "Cl" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 0 "1"] [0 4 "1"] [2 5 "1"] [0 6 "1"] [2 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 3 4 6] :attrs "Th0"}
                       {:site 2 :ligands [1 3 5 7] :attrs "Th0"}]
    }"#), vec![Stereogenicity::Stereogenic, Stereogenicity::Stereogenic])]
#[case::chlorofluoropropane(mol_dsl_concrete!(r#"{
        :atoms ["C" "C" "C#h3" "F" "Cl" "H" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [0 3 "1"] [0 4 "1"] [0 5 "1"] [1 6 "1"] [1 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 3 4 5] :attrs "Th0"}
                       {:site 1 :ligands [0 2 6 7] :attrs "Th0"}]
    }"#), vec![Stereogenicity::Stereogenic, Stereogenicity::Prostereogenic])]
#[case::chlorofluoropropane_without_c1(mol_dsl_concrete!(r#"{
        :atoms ["C" "C" "C#h3" "F" "Cl" "H" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [0 3 "1"] [0 4 "1"] [0 5 "1"] [1 6 "1"] [1 7 "1"]]
        :stereo-atoms [{:site 1 :ligands [0 2 6 7] :attrs "Th0"}]
    }"#), vec![Stereogenicity::Prostereogenic])]
#[case::trimethylcyclohexane_all_cis(mol_dsl_concrete!(r#"{
        :atoms ["C" "C#h2" "C" "C#h2" "C" "C#h2" "C#h3" "C#h3" "C#h3" "H" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                [0 6 "1"] [2 7 "1"] [4 8 "1"] [0 9 "1"] [2 10 "1"] [4 11 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 5 6 9] :attrs "Th0"}
                       {:site 2 :ligands [1 3 7 10] :attrs "Th1"}
                       {:site 4 :ligands [3 5 8 11] :attrs "Th1"}]
    }"#), vec![Stereogenicity::Stereogenic, Stereogenicity::Stereogenic, Stereogenicity::Stereogenic])]
#[case::trimethylcyclohexane_cis_trans(mol_dsl_concrete!(r#"{
        :atoms ["C" "C#h2" "C" "C#h2" "C" "C#h2" "C#h3" "C#h3" "C#h3" "H" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                [0 6 "1"] [2 7 "1"] [4 8 "1"] [0 9 "1"] [2 10 "1"] [4 11 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 5 6 9] :attrs "Th0"}
                       {:site 2 :ligands [1 3 7 10] :attrs "Th0"}
                       {:site 4 :ligands [3 5 8 11] :attrs "Th0"}]
    }"#), vec![Stereogenicity::Stereogenic, Stereogenicity::Prostereogenic, Stereogenicity::Prostereogenic])]
fn test_stereo_atom_view_stereogenicity(
    #[case] molecule: Molecule,
    #[case] expected: Vec<Stereogenicity>,
    #[values(false, true)] include_constraints: bool,
) {
    let original = molecule.clone();
    let symmetry = molecule
        .graph_symmetry(&GraphSymmetryConfig {
            include_constraints,
            ..Default::default()
        })
        .into_graph_symmetry();
    let actual = molecule
        .stereo_atoms()
        .iter()
        .map(|stereo| stereo.with_symmetry(&symmetry).stereogenicity())
        .collect::<Vec<_>>();
    assert_eq!(actual, expected);
    assert_eq!(molecule, original);
}

#[rstest]
#[case::square_planar(mol_dsl_concrete!(r#"{
    :atoms ["Pt#c-2" "N#c+#h3" "N#c+#h3" "Cl" "Cl"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Sp0"}]
}"#), 0..3)]
#[case::trigonal_bipyramidal(mol_dsl_concrete!(r#"{
    :atoms ["P" "Cl" "F" "F" "F" "F"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [0 5 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4 5] :attrs "Tb0"}]
}"#), 0..20)]
#[case::octahedral(mol_dsl_concrete!(r#"{
    :atoms ["Pt" "F" "Cl" "Br" "I" "N" "O"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [0 5 "1"] [0 6 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4 5 6] :attrs "Oh0"}]
}"#), 0..30)]
fn test_stereo_atom_view_stereogenicity_cosets(
    #[case] mut molecule: Molecule,
    #[case] cosets: Range<u32>,
) {
    let kind = molecule.stereo_atom(StereoAtomId(0)).kind();
    for coset in cosets {
        molecule
            .stereo_atom_mut(StereoAtomId(0))
            .attributes_mut()
            .configuration = StereoConfigurationForm::kinded(kind, coset);
        let symmetry = molecule
            .graph_symmetry(&GraphSymmetryConfig::default())
            .into_graph_symmetry();
        assert_eq!(
            molecule
                .stereo_atom(StereoAtomId(0))
                .with_symmetry(&symmetry)
                .stereogenicity(),
            Stereogenicity::Stereogenic
        );
    }
}

#[rstest]
#[case::moving_focus(mol_dsl_concrete!(r#"{
    :atoms ["C" "C" "Cl" "Cl" "H" "H" "H" "H"]
    :bonds [[0 1 "1"] [0 2 "1"] [1 3 "1"] [0 4 "1"] [0 5 "1"] [1 6 "1"] [1 7 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 4 5] :attrs "Th0"} {:site 1 :ligands [0 3 6 7] :attrs "Th0"}]
}"#),vec![Stereogenicity::Prostereogenic,Stereogenicity::Prostereogenic])]
#[case::methane(mol_dsl_concrete!(r#"{
    :atoms ["C" "H" "H" "H" "H"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]
}"#),vec![Stereogenicity::Nonstereogenic])]
fn test_stereo_atom_view_stereogenicity_replacement(
    #[case] molecule: Molecule,
    #[case] expected: Vec<Stereogenicity>,
) {
    let symmetry = molecule
        .graph_symmetry(&GraphSymmetryConfig::default())
        .into_graph_symmetry();
    assert_eq!(
        molecule
            .stereo_atoms()
            .iter()
            .map(|stereo| stereo.with_symmetry(&symmetry).stereogenicity())
            .collect::<Vec<_>>(),
        expected
    );
}

#[rstest]
#[case::undetermined(StereoCoset::Undetermined)]
#[case::singleton(StereoCoset::lit_set([0]))]
#[case::term(StereoCoset::term(StereoTerm::var_in("configuration",[])))]
fn test_stereo_atom_view_stereogenicity_raw(
    mut dichlorofluoromethane: Molecule,
    #[case] coset: StereoCoset,
) {
    dichlorofluoromethane
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = StereoConfigurationForm::kinded(StereoKind::Tetrahedral, coset);
    let original = dichlorofluoromethane.clone();
    let symmetry = dichlorofluoromethane
        .graph_symmetry(&GraphSymmetryConfig::default())
        .into_graph_symmetry();
    assert_eq!(
        dichlorofluoromethane
            .stereo_atom(StereoAtomId(0))
            .with_symmetry(&symmetry)
            .stereogenicity(),
        Stereogenicity::Prostereogenic
    );
    assert_eq!(dichlorofluoromethane, original);
}

#[rstest]
#[case::excluded(false, Stereogenicity::Prostereogenic, Topicity::Enantiotopic, true)]
#[case::included(true, Stereogenicity::Stereogenic, Topicity::Distinct, false)]
fn test_stereo_atom_view_stereogenicity_constraints(
    mut dichlorofluoromethane: Molecule,
    #[case] include_constraints: bool,
    #[case] expected: Stereogenicity,
    #[case] topicity: Topicity,
    #[case] action: bool,
) {
    dichlorofluoromethane
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .constraints
        .set(StereoAtomConstraintForm::Topicity(TopicityForm {
            pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(7)),
            relation: TopicityRelationForm::Lit(Topicity::Homotopic),
        }));
    let original = dichlorofluoromethane.clone();
    let symmetry = dichlorofluoromethane
        .graph_symmetry(&GraphSymmetryConfig {
            include_constraints,
            ..Default::default()
        })
        .into_graph_symmetry();
    let stereo = dichlorofluoromethane
        .stereo_atom(StereoAtomId(0))
        .with_symmetry(&symmetry);
    assert_eq!(stereo.stereogenicity(), expected);
    assert_eq!(
        stereo.topicity(StereoLigandPosition(1), StereoLigandPosition(2)),
        topicity
    );
    assert_eq!(
        stereo.ligand_symmetry(OrientedLigandPermutation {
            permutation: LigandPermutation(Permutation::from_image(&[0, 2, 1, 3])),
            orientation: Orientation::Improper,
        }),
        action
    );
    assert_eq!(dichlorofluoromethane, original);
}

#[rstest]
#[case::degree(StereoConfigurationForm::kinded(StereoKind::TrigonalBipyramidal, 0))]
#[case::coset(StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 2))]
fn test_stereo_atom_view_stereogenicity_encoding(
    mut dichlorofluoromethane: Molecule,
    #[case] configuration: StereoConfigurationForm,
) {
    dichlorofluoromethane
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = configuration;
    let symmetry = dichlorofluoromethane
        .graph_symmetry(&GraphSymmetryConfig::default())
        .into_graph_symmetry();
    let stereo = dichlorofluoromethane
        .stereo_atom(StereoAtomId(0))
        .with_symmetry(&symmetry);
    assert_eq!(stereo.stereogenicity(), Stereogenicity::Stereogenic);
    assert_eq!(
        stereo.topicity(StereoLigandPosition(1), StereoLigandPosition(2)),
        Topicity::Distinct
    );
    assert!(!stereo.ligand_symmetry(OrientedLigandPermutation {
        permutation: LigandPermutation(Permutation::identity(4)),
        orientation: Orientation::Proper,
    }));
}

#[rstest]
#[case::missing_data(false)]
#[case::kindless(true)]
#[should_panic(expected = "stereo view has")]
fn test_stereo_atom_view_stereogenicity_error(mut molecule: Molecule, #[case] kindless: bool) {
    if kindless {
        molecule
            .stereo_atom_mut(StereoAtomId(0))
            .attributes_mut()
            .configuration = StereoConfigurationForm::Undetermined;
    }
    let symmetry = molecule
        .graph_symmetry(&GraphSymmetryConfig::default())
        .into_graph_symmetry();
    let stereo = molecule.stereo_atom(StereoAtomId(0));
    let stereo = if kindless {
        stereo.with_symmetry(&symmetry)
    } else {
        stereo
    };
    let _ = stereo.stereogenicity();
}

#[rstest]
#[case::chlorofluoropropane(mol_dsl_concrete!(r#"{
        :atoms ["C" "C" "C#h3" "F" "Cl" "H" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [0 3 "1"] [0 4 "1"] [0 5 "1"] [1 6 "1"] [1 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 3 4 5] :attrs "Th0"}
                       {:site 1 :ligands [0 2 6 7] :attrs "Th0"}]
    }"#), 1, 2, 3, Topicity::Diastereotopic)]
#[case::chlorofluoropropane_without_c1(mol_dsl_concrete!(r#"{
        :atoms ["C" "C" "C#h3" "F" "Cl" "H" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [0 3 "1"] [0 4 "1"] [0 5 "1"] [1 6 "1"] [1 7 "1"]]
        :stereo-atoms [{:site 1 :ligands [0 2 6 7] :attrs "Th0"}]
    }"#), 0, 2, 3, Topicity::Enantiotopic)]
#[case::trichloropentane_meso(mol_dsl_concrete!(r#"{
        :atoms ["C#h3" "C" "C" "C" "C#h3" "Cl" "Cl" "Cl" "H" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [1 5 "1"] [2 6 "1"] [3 7 "1"]
                [1 8 "1"] [2 9 "1"] [3 10 "1"]]
        :stereo-atoms [{:site 1 :ligands [0 2 5 8] :attrs "Th0"}
                       {:site 2 :ligands [1 3 6 9] :attrs "Th0"}
                       {:site 3 :ligands [2 4 7 10] :attrs "Th0"}]
    }"#), 1, 0, 1, Topicity::Enantiotopic)]
#[case::trichloropentane_chiral(mol_dsl_concrete!(r#"{
        :atoms ["C#h3" "C" "C" "C" "C#h3" "Cl" "Cl" "Cl" "H" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [1 5 "1"] [2 6 "1"] [3 7 "1"]
                [1 8 "1"] [2 9 "1"] [3 10 "1"]]
        :stereo-atoms [{:site 1 :ligands [0 2 5 8] :attrs "Th0"}
                       {:site 2 :ligands [1 3 6 9] :attrs "Th0"}
                       {:site 3 :ligands [2 4 7 10] :attrs "Th1"}]
    }"#), 1, 0, 1, Topicity::Diastereotopic)]
#[case::distinct(mol_dsl_concrete!(r#"{
        :atoms ["C" "C" "C#h3" "F" "Cl" "H" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [0 3 "1"] [0 4 "1"] [0 5 "1"] [1 6 "1"] [1 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 3 4 5] :attrs "Th0"}
                       {:site 1 :ligands [0 2 6 7] :attrs "Th0"}]
    }"#), 1, 0, 1, Topicity::Distinct)]
#[case::methane(Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::H),
                    AtomForm::from_element(Element::H), AtomForm::from_element(Element::H), AtomForm::from_element(Element::H)],
        bonds: (1..=4).map(|atom| (AtomId(0), AtomId(atom), BondForm::from_order(1))).collect(),
        stereo_atoms: vec![(AtomId(0), (1..=4).map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom)).collect(),
            StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32))],
        ..Default::default()
    }), 0, 0, 1, Topicity::Homotopic)]
fn test_stereo_atom_view_topicity(
    #[case] molecule: Molecule,
    #[case] record: u32,
    #[case] first: u32,
    #[case] second: u32,
    #[case] expected: Topicity,
    #[values(false, true)] include_constraints: bool,
) {
    let original = molecule.clone();
    let symmetry = molecule
        .graph_symmetry(&GraphSymmetryConfig {
            include_constraints,
            ..Default::default()
        })
        .into_graph_symmetry();
    let stereo = molecule
        .stereo_atom(StereoAtomId(record))
        .with_symmetry(&symmetry);
    assert_eq!(
        stereo.topicity(StereoLigandPosition(first), StereoLigandPosition(second)),
        expected
    );
    assert_eq!(
        stereo.topicity(StereoLigandPosition(second), StereoLigandPosition(first)),
        expected
    );
    assert_eq!(molecule, original);
}

#[rstest]
#[case::same(1, 1, Topicity::Homotopic)]
#[case::different_virtual(1, 2, Topicity::Distinct)]
#[case::atom_virtual(0, 1, Topicity::Distinct)]
fn test_stereo_atom_view_topicity_virtual(
    virtual_ligand_molecule: Molecule,
    #[case] first: u32,
    #[case] second: u32,
    #[case] expected: Topicity,
) {
    let symmetry = virtual_ligand_molecule
        .graph_symmetry(&GraphSymmetryConfig::default())
        .into_graph_symmetry();
    assert_eq!(
        virtual_ligand_molecule
            .stereo_atom(StereoAtomId(0))
            .with_symmetry(&symmetry)
            .topicity(StereoLigandPosition(first), StereoLigandPosition(second)),
        expected
    );
}

#[rstest]
#[case::missing_data(false)]
#[case::kindless(true)]
#[should_panic(expected = "stereo view has")]
fn test_stereo_atom_view_topicity_error(mut molecule: Molecule, #[case] kindless: bool) {
    if kindless {
        molecule
            .stereo_atom_mut(StereoAtomId(0))
            .attributes_mut()
            .configuration = StereoConfigurationForm::Undetermined;
    }
    let symmetry = molecule
        .graph_symmetry(&GraphSymmetryConfig::default())
        .into_graph_symmetry();
    let stereo = molecule.stereo_atom(StereoAtomId(0));
    let stereo = if kindless {
        stereo.with_symmetry(&symmetry)
    } else {
        stereo
    };
    let _ = stereo.topicity(StereoLigandPosition(0), StereoLigandPosition(1));
}

#[rstest]
#[case::first(4, 0)]
#[case::second(0, 4)]
#[should_panic(expected = "ligand id must refer to a ligand of this stereo atom")]
fn test_stereo_atom_view_topicity_position_error(
    dichlorofluoromethane: Molecule,
    #[case] first: u32,
    #[case] second: u32,
) {
    let symmetry = dichlorofluoromethane
        .graph_symmetry(&GraphSymmetryConfig::default())
        .into_graph_symmetry();
    let _ = dichlorofluoromethane
        .stereo_atom(StereoAtomId(0))
        .with_symmetry(&symmetry)
        .topicity(StereoLigandPosition(first), StereoLigandPosition(second));
}

#[rstest]
#[case::identity(vec![0,1,2,3],Orientation::Proper,true)]
#[case::improper_exchange(vec![0,2,1,3],Orientation::Improper,true)]
#[case::proper_exchange(vec![0,2,1,3],Orientation::Proper,false)]
#[case::improper_identity(vec![0,1,2,3],Orientation::Improper,false)]
#[case::virtual_exchange(vec![3,1,2,0],Orientation::Improper,false)]
#[case::degree(vec![0,1,2],Orientation::Proper,false)]
fn test_stereo_atom_view_ligand_symmetry(
    dichlorofluoromethane: Molecule,
    #[case] image: Vec<usize>,
    #[case] orientation: Orientation,
    #[case] expected: bool,
) {
    let original = dichlorofluoromethane.clone();
    let symmetry = dichlorofluoromethane
        .graph_symmetry(&GraphSymmetryConfig::default())
        .into_graph_symmetry();
    let stereo = dichlorofluoromethane
        .stereo_atom(StereoAtomId(0))
        .with_symmetry(&symmetry);
    assert_eq!(
        stereo.ligand_symmetry(OrientedLigandPermutation {
            permutation: LigandPermutation(Permutation::from_image(&image)),
            orientation,
        }),
        expected
    );
    assert_eq!(dichlorofluoromethane, original);
}

#[rstest]
#[case::missing_data(false)]
#[case::kindless(true)]
#[should_panic(expected = "stereo view has")]
fn test_stereo_atom_view_ligand_symmetry_error(mut molecule: Molecule, #[case] kindless: bool) {
    if kindless {
        molecule
            .stereo_atom_mut(StereoAtomId(0))
            .attributes_mut()
            .configuration = StereoConfigurationForm::Undetermined;
    }
    let symmetry = molecule
        .graph_symmetry(&GraphSymmetryConfig::default())
        .into_graph_symmetry();
    let stereo = molecule.stereo_atom(StereoAtomId(0));
    let stereo = if kindless {
        stereo.with_symmetry(&symmetry)
    } else {
        stereo
    };
    let _ = stereo.ligand_symmetry(OrientedLigandPermutation {
        permutation: LigandPermutation(Permutation::identity(4)),
        orientation: Orientation::Proper,
    });
}

#[rstest]
fn test_stereo_bond_view_with_symmetry_context(difluorodichloroethene: Molecule) {
    let symmetry = Molecule::default()
        .graph_symmetry(&GraphSymmetryConfig::default())
        .into_graph_symmetry();
    let stereo = difluorodichloroethene
        .stereo_bond(StereoBondId(0))
        .with_symmetry(&symmetry);
    let _ = stereo.stereogenicity();
    let _ = stereo.topicity(StereoLigandPosition(0), StereoLigandPosition(2));
    let _ = stereo.ligand_symmetry(OrientedLigandPermutation {
        permutation: LigandPermutation(Permutation::identity(4)),
        orientation: Orientation::Proper,
    });
}

#[rstest]
#[case::stereogenic(mol_dsl_concrete!(r#"{
    :atoms ["C" "C" "F" "Cl" "F" "Cl"]
    :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
    :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct0"}]
}"#),Stereogenicity::Stereogenic)]
#[case::prostereogenic(mol_dsl_concrete!(r#"{
    :atoms ["C" "C" "H" "H" "F" "Cl"]
    :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
    :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct0"}]
}"#),Stereogenicity::Prostereogenic)]
#[case::nonstereogenic(mol_dsl_concrete!(r#"{
    :atoms ["C" "C" "H" "H" "H" "H"]
    :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
    :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct0"}]
}"#),Stereogenicity::Nonstereogenic)]
fn test_stereo_bond_view_stereogenicity(
    #[case] molecule: Molecule,
    #[case] expected: Stereogenicity,
) {
    let symmetry = molecule
        .graph_symmetry(&GraphSymmetryConfig::default())
        .into_graph_symmetry();
    assert_eq!(
        molecule
            .stereo_bond(StereoBondId(0))
            .with_symmetry(&symmetry)
            .stereogenicity(),
        expected
    );
}

#[rstest]
#[case::missing_data(false)]
#[case::kindless(true)]
#[should_panic(expected = "stereo view has")]
fn test_stereo_bond_view_stereogenicity_error(mut molecule: Molecule, #[case] kindless: bool) {
    if kindless {
        molecule
            .stereo_bond_mut(StereoBondId(0))
            .attributes_mut()
            .configuration = StereoConfigurationForm::Undetermined;
    }
    let symmetry = molecule
        .graph_symmetry(&GraphSymmetryConfig::default())
        .into_graph_symmetry();
    let stereo = molecule.stereo_bond(StereoBondId(0));
    let stereo = if kindless {
        stereo.with_symmetry(&symmetry)
    } else {
        stereo
    };
    let _ = stereo.stereogenicity();
}

#[rstest]
#[case::same(0, 0, Topicity::Homotopic)]
#[case::opposite_ends(0, 2, Topicity::Homotopic)]
#[case::different_elements(0, 1, Topicity::Distinct)]
fn test_stereo_bond_view_topicity(
    difluorodichloroethene: Molecule,
    #[case] first: u32,
    #[case] second: u32,
    #[case] expected: Topicity,
) {
    let symmetry = difluorodichloroethene
        .graph_symmetry(&GraphSymmetryConfig::default())
        .into_graph_symmetry();
    let stereo = difluorodichloroethene
        .stereo_bond(StereoBondId(0))
        .with_symmetry(&symmetry);
    assert_eq!(
        stereo.topicity(StereoLigandPosition(first), StereoLigandPosition(second)),
        expected
    );
    assert_eq!(
        stereo.topicity(StereoLigandPosition(second), StereoLigandPosition(first)),
        expected
    );
}

#[rstest]
#[case::missing_data(false)]
#[case::kindless(true)]
#[should_panic(expected = "stereo view has")]
fn test_stereo_bond_view_topicity_error(mut molecule: Molecule, #[case] kindless: bool) {
    if kindless {
        molecule
            .stereo_bond_mut(StereoBondId(0))
            .attributes_mut()
            .configuration = StereoConfigurationForm::Undetermined;
    }
    let symmetry = molecule
        .graph_symmetry(&GraphSymmetryConfig::default())
        .into_graph_symmetry();
    let stereo = molecule.stereo_bond(StereoBondId(0));
    let stereo = if kindless {
        stereo.with_symmetry(&symmetry)
    } else {
        stereo
    };
    let _ = stereo.topicity(StereoLigandPosition(0), StereoLigandPosition(1));
}

#[rstest]
#[case::first(4, 0)]
#[case::second(0, 4)]
#[should_panic(expected = "ligand id must refer to a ligand of this stereo bond")]
fn test_stereo_bond_view_topicity_position_error(
    difluorodichloroethene: Molecule,
    #[case] first: u32,
    #[case] second: u32,
) {
    let symmetry = difluorodichloroethene
        .graph_symmetry(&GraphSymmetryConfig::default())
        .into_graph_symmetry();
    let _ = difluorodichloroethene
        .stereo_bond(StereoBondId(0))
        .with_symmetry(&symmetry)
        .topicity(StereoLigandPosition(first), StereoLigandPosition(second));
}

#[rstest]
#[case::identity(vec![0,1,2,3],Orientation::Proper,true)]
#[case::endpoint_exchange(vec![2,3,0,1],Orientation::Proper,true)]
#[case::improper_endpoint_exchange(vec![2,3,0,1],Orientation::Improper,true)]
#[case::incompatible(vec![0,2,1,3],Orientation::Proper,false)]
#[case::degree(vec![0,1,2],Orientation::Proper,false)]
fn test_stereo_bond_view_ligand_symmetry(
    difluorodichloroethene: Molecule,
    #[case] image: Vec<usize>,
    #[case] orientation: Orientation,
    #[case] expected: bool,
) {
    let original = difluorodichloroethene.clone();
    let symmetry = difluorodichloroethene
        .graph_symmetry(&GraphSymmetryConfig::default())
        .into_graph_symmetry();
    assert_eq!(
        difluorodichloroethene
            .stereo_bond(StereoBondId(0))
            .with_symmetry(&symmetry)
            .ligand_symmetry(OrientedLigandPermutation {
                permutation: LigandPermutation(Permutation::from_image(&image)),
                orientation
            }),
        expected
    );
    assert_eq!(difluorodichloroethene, original);
}

#[rstest]
fn test_stereo_bond_view_ligand_symmetry_virtual() {
    let molecule = mol_dsl_concrete!(
        r#"{
        :atoms ["C#h" "C#h" "F" "F"]
        :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
        :stereo-bonds [{:site 0 :ligands [2 [:h 0] 3 [:h 1]] :attrs "Ct0"}]
    }"#
    );
    let symmetry = molecule
        .graph_symmetry(&GraphSymmetryConfig::default())
        .into_graph_symmetry();
    let stereo = molecule
        .stereo_bond(StereoBondId(0))
        .with_symmetry(&symmetry);
    assert!(stereo.ligand_symmetry(OrientedLigandPermutation {
        permutation: LigandPermutation(Permutation::identity(4)),
        orientation: Orientation::Proper,
    }));
    assert!(!stereo.ligand_symmetry(OrientedLigandPermutation {
        permutation: LigandPermutation(Permutation::from_image(&[2, 3, 0, 1])),
        orientation: Orientation::Proper,
    }));
    assert_eq!(
        stereo.topicity(StereoLigandPosition(0), StereoLigandPosition(2)),
        Topicity::Distinct
    );
    assert_eq!(
        stereo.topicity(StereoLigandPosition(1), StereoLigandPosition(1)),
        Topicity::Homotopic
    );
}

#[rstest]
#[case::missing_data(false)]
#[case::kindless(true)]
#[should_panic(expected = "stereo view has")]
fn test_stereo_bond_view_ligand_symmetry_error(mut molecule: Molecule, #[case] kindless: bool) {
    if kindless {
        molecule
            .stereo_bond_mut(StereoBondId(0))
            .attributes_mut()
            .configuration = StereoConfigurationForm::Undetermined;
    }
    let symmetry = molecule
        .graph_symmetry(&GraphSymmetryConfig::default())
        .into_graph_symmetry();
    let stereo = molecule.stereo_bond(StereoBondId(0));
    let stereo = if kindless {
        stereo.with_symmetry(&symmetry)
    } else {
        stereo
    };
    let _ = stereo.ligand_symmetry(OrientedLigandPermutation {
        permutation: LigandPermutation(Permutation::identity(4)),
        orientation: Orientation::Proper,
    });
}
