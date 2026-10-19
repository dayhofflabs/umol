use pretty_assertions::assert_eq;
use rstest::{fixture, rstest};
use umol_chem::element::Element;
use umol_graph_ir::ir::{
    AtomConstraintForm, AtomForm, AtomId, BondConstraintForm, BondForm, BondId, BooleanForm,
    Constraint, ElementForm, Entity, FluxionalityForm, LigandPermutation, Molecule,
    MoleculeConstraint, MoleculeEntries, NumForm, StereoAtomConstraintForm, StereoAtomId,
    StereoBondConstraintForm, StereoBondId, StereoConfigurationForm, StereoCoset, StereoKind,
    StereoTerm, TetrahedralStereoForm,
};
use umol_graph_ir::mol_dsl_concrete;
use umol_perm::Permutation;

use crate::ops::transform::{FoldHydrogensError, HydrogenFolder, Transformer};

#[fixture]
fn methyl() -> Molecule {
    mol_dsl_concrete!(r#"{:atoms ["H" "C#h2"] :bonds [[0 1 "1"]]}"#)
}

#[rstest]
#[case::single(
    mol_dsl_concrete!(r#"{:atoms ["H" "C#h2"] :bonds [[0 1 "1"]]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C#h3"]}"#),
)]
#[case::several_bearers(
    mol_dsl_concrete!(r#"{
        :atoms ["H" "C#h" "H" "O" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"] [3 4 "1"] [1 5 "1"]]
    }"#),
    mol_dsl_concrete!(r#"{:atoms ["C#h4" "O#h"] :bonds [[0 1 "1"]]}"#),
)]
#[case::constraints(
    mol_dsl_concrete!(r#"{
        :atoms ["H" "C#h#v3#D3#X4#H3#R!" "H" "F"]
        :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
    }"#),
    mol_dsl_concrete!(r#"{:atoms ["C#h3#v1#D1#X4#H3#R!" "F"] :bonds [[0 1 "1"]]}"#),
)]
#[case::isotope_and_charge_preserved(
    mol_dsl_concrete!(r#"{
        :atoms ["H" "N#i15#c+#h#n" "H#i2" "H#c-"]
        :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
    }"#),
    mol_dsl_concrete!(r#"{
        :atoms ["N#i15#c+#h2#n" "H#i2" "H#c-"]
        :bonds [[0 1 "1"] [0 2 "1"]]
    }"#),
)]
#[case::tetrahedral_frame(
    mol_dsl_concrete!(r#"{
        :atoms ["H" "C#T0" "F" "Cl" "Br"]
        :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"] [1 4 "1"]]
        :stereo-atoms [{:site 1 :ligands [3 0 4 2] :attrs "Th1"}]
    }"#),
    mol_dsl_concrete!(r#"{
        :atoms ["C#h#T1" "F" "Cl" "Br"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
        :stereo-atoms [{:site 0 :ligands [2 [:h 0] 3 1] :attrs "Th1"}]
    }"#),
)]
#[case::tetrahedral_assertion_without_overlay(
    mol_dsl_concrete!(r#"{
        :atoms ["H" "C#T0" "F" "Cl" "Br"]
        :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"] [1 4 "1"]]
    }"#),
    mol_dsl_concrete!(r#"{
        :atoms ["C#h#T1" "F" "Cl" "Br"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
    }"#),
)]
#[case::last_actual_ligand(
    mol_dsl_concrete!(r#"{
        :atoms ["C#T1" "F" "Cl" "Br" "H"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    }"#),
    mol_dsl_concrete!(r#"{
        :atoms ["C#h#T1" "F" "Cl" "Br"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
    }"#),
)]
#[case::lone_pair_frame(
    mol_dsl_concrete!(r#"{
        :atoms ["H" "N#n#T0" "F" "Cl"]
        :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
        :stereo-atoms [{:site 1 :ligands [0 2 3 [:lp 1]] :attrs "Th0"}]
    }"#),
    mol_dsl_concrete!(r#"{
        :atoms ["N#h#n#T0" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"]]
        :stereo-atoms [{:site 0 :ligands [[:h 0] 1 2 [:lp 0]] :attrs "Th0"}]
    }"#),
)]
#[case::cis_trans_reversed_endpoint_blocks(
    mol_dsl_concrete!(r#"{
        :atoms ["H" "C" "F" "C" "Cl" "H"]
        :bonds [[1 3 "2#C0"] [0 1 "1"] [1 2 "1"] [3 4 "1"] [3 5 "1"]]
        :stereo-bonds [{:site 0 :ligands [5 4 0 2] :attrs "Ct1"}]
    }"#),
    mol_dsl_concrete!(r#"{
        :atoms ["C#h" "F" "C#h" "Cl"]
        :bonds [[0 2 "2#C1"] [0 1 "1"] [2 3 "1"]]
        :stereo-bonds [{:site 0 :ligands [[:h 2] 3 [:h 0] 1] :attrs "Ct1"}]
    }"#),
)]
#[case::cis_trans_both_endpoint_frames(
    mol_dsl_concrete!(r#"{
        :atoms ["H" "H" "C" "C" "F" "Cl"]
        :bonds [[2 3 "2#C0"] [0 2 "1"] [1 3 "1"] [2 4 "1"] [3 5 "1"]]
        :stereo-bonds [{:site 0 :ligands [0 4 1 5] :attrs "Ct0"}]
    }"#),
    mol_dsl_concrete!(r#"{
        :atoms ["C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "2#C0"] [0 2 "1"] [1 3 "1"]]
        :stereo-bonds [{:site 0 :ligands [[:h 0] 2 [:h 1] 3] :attrs "Ct0"}]
    }"#),
)]
#[case::shared_stereo_bearer(
    mol_dsl_concrete!(r#"{
        :atoms ["H" "N#n" "F" "C#h" "Cl"]
        :bonds [[1 3 "2"] [0 1 "1"] [1 2 "1"] [3 4 "1"]]
        :stereo-atoms [{:site 1 :ligands [3 0 2 [:lp 1]] :attrs "Th0"}]
        :stereo-bonds [{:site 0 :ligands [0 2 4 [:h 3]] :attrs "Ct0"}]
    }"#),
    mol_dsl_concrete!(r#"{
        :atoms ["N#h#n" "F" "C#h" "Cl"]
        :bonds [[0 2 "2"] [0 1 "1"] [2 3 "1"]]
        :stereo-atoms [{:site 0 :ligands [2 [:h 0] 1 [:lp 0]] :attrs "Th0"}]
        :stereo-bonds [{:site 0 :ligands [[:h 0] 1 3 [:h 2]] :attrs "Ct0"}]
    }"#),
)]
#[case::retained_relations(
    mol_dsl_concrete!(r#"{
        :atoms ["H" "O" "C" "C" "C"]
        :bonds [[0 1 "1"] [2 3 "1"] [3 4 "1"] [2 4 "1"]]
        :dative-bonds [{:donors [1] :acceptor 2 :attrs "1"}]
        :aromatic-systems [{:atoms [2 3 4] :attrs "[1,1,1]"}]
        :multicenter-bonds [{:atoms [1 3 4] :attrs "[2,0,0]"}]
        :noncovalent-bonds [{:atoms [1 2] :attrs "*"}]
    }"#),
    mol_dsl_concrete!(r#"{
        :atoms ["O#h" "C" "C" "C"]
        :bonds [[1 2 "1"] [2 3 "1"] [1 3 "1"]]
        :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]
        :aromatic-systems [{:atoms [1 2 3] :attrs "[1,1,1]"}]
        :multicenter-bonds [{:atoms [0 2 3] :attrs "[2,0,0]"}]
        :noncovalent-bonds [{:atoms [0 1] :attrs "*"}]
    }"#),
)]
fn test_hydrogen_folder_transform(#[case] input: Molecule, #[case] expected: Molecule) {
    assert_eq!(
        HydrogenFolder.transform(input.clone()),
        Ok(expected.clone())
    );
    let mut in_place = input;
    assert_eq!(HydrogenFolder.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected);
    assert_eq!(HydrogenFolder.transform(in_place), Ok(expected.clone()));
    assert_eq!(expected.clone().edit().finish(), Ok(expected));
}

#[rstest]
#[case::empty(Molecule::new())]
#[case::already_folded(mol_dsl_concrete!(r#"{:atoms ["C#h4"]}"#))]
#[case::hydrogen_molecule(mol_dsl_concrete!(r#"{:atoms ["H" "H"] :bonds [[0 1 "1"]]}"#))]
#[case::isolated(mol_dsl_concrete!(r#"{:atoms ["H" "C"]}"#))]
#[case::two_neighbors(mol_dsl_concrete!(r#"{
    :atoms ["C" "H" "C"] :bonds [[0 1 "1"] [1 2 "1"]]
}"#))]
#[case::negative_without_eligible_hydrogen(mol_dsl_concrete!(r#"{
    :atoms ["C#h-1" "H#i2"] :bonds [[0 1 "1"]]
}"#))]
#[case::dative_donor(mol_dsl_concrete!(r#"{
    :atoms ["H" "C" "B"] :bonds [[0 1 "1"]]
    :dative-bonds [{:donors [0] :acceptor 2 :attrs "1"}]
}"#))]
#[case::dative_acceptor(mol_dsl_concrete!(r#"{
    :atoms ["H" "C" "N"] :bonds [[0 1 "1"]]
    :dative-bonds [{:donors [2] :acceptor 0 :attrs "1"}]
}"#))]
#[case::aromatic_participant(mol_dsl_concrete!(r#"{
    :atoms ["H" "C" "C"] :bonds [[0 1 "1"]]
    :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,1]"}]
}"#))]
#[case::multicenter_participant(mol_dsl_concrete!(r#"{
    :atoms ["H" "B" "B"] :bonds [[0 1 "1"]]
    :multicenter-bonds [{:atoms [0 1 2] :attrs "[0,1,1]"}]
}"#))]
#[case::noncovalent_participant(mol_dsl_concrete!(r#"{
    :atoms ["H" "C" "O"] :bonds [[0 1 "1"]]
    :noncovalent-bonds [{:atoms [0 2] :attrs "*"}]
}"#))]
#[case::stereo_atom_site(mol_dsl_concrete!(r#"{
    :atoms ["H" "C"] :bonds [[0 1 "1"]]
    :stereo-atoms [{:site 0 :ligands [1] :attrs "*"}]
}"#))]
#[case::stereo_bond_site(mol_dsl_concrete!(r#"{
    :atoms ["H" "C#h" "F"] :bonds [[0 1 "1"] [1 2 "1"]]
    :stereo-bonds [{:site 0 :ligands [[:h 0] [:lp 0] 2 [:h 1]] :attrs "Ct0"}]
}"#))]
fn test_hydrogen_folder_transform_identity(#[case] input: Molecule) {
    assert_eq!(HydrogenFolder.transform(input.clone()), Ok(input.clone()));
    let mut in_place = input.clone();
    assert_eq!(HydrogenFolder.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, input);
}

#[rstest]
#[case::deuterium("H#i2")]
#[case::tritium("H#i3")]
#[case::protium_label("H#i1")]
#[case::proton("H#c+")]
#[case::nonliteral_charge("H#c{0,1}")]
#[case::hydride("H#c-")]
#[case::radical("H#u#s2")]
#[case::spin("H#s3")]
#[case::lone_pair("H#n")]
#[case::implicit_hydrogen("H#h")]
#[case::assertion("H#v1")]
fn test_hydrogen_folder_transform_atom_identity(mut methyl: Molecule, #[case] atom: &str) {
    *methyl.atom_mut(AtomId(0)).attributes_mut() = AtomForm::from(atom).into_concrete();
    assert_eq!(HydrogenFolder.transform(methyl.clone()), Ok(methyl));
}

#[rstest]
#[case::double("2")]
#[case::charged("1#c+")]
#[case::radical("1#u#s2")]
#[case::spin("1#s3")]
#[case::assertion("1#a!")]
fn test_hydrogen_folder_transform_bond_identity(mut methyl: Molecule, #[case] bond: &str) {
    *methyl.bond_mut(BondId(0)).attributes_mut() = BondForm::from(bond).into_concrete();
    assert_eq!(HydrogenFolder.transform(methyl.clone()), Ok(methyl));
}

#[rstest]
#[case::unknown(NumForm::Undetermined)]
#[case::set(NumForm::from(vec![0, 1]))]
#[case::variable(NumForm::var("h"))]
fn test_hydrogen_folder_transform_count_identity(mut methyl: Molecule, #[case] count: NumForm) {
    methyl
        .atom_mut(AtomId(1))
        .attributes_mut()
        .implicit_hydrogens = count;
    assert_eq!(HydrogenFolder.transform(methyl.clone()), Ok(methyl));
}

#[rstest]
#[case::unknown(ElementForm::Undetermined)]
#[case::hydrogen_or_carbon(ElementForm::lit_set([Element::H, Element::C]))]
#[case::variable(ElementForm::var("element"))]
fn test_hydrogen_folder_transform_element_identity(
    mut methyl: Molecule,
    #[case] element: ElementForm,
) {
    methyl.atom_mut(AtomId(1)).attributes_mut().element = element;
    assert_eq!(HydrogenFolder.transform(methyl.clone()), Ok(methyl));
}

#[rstest]
#[case::non_hydrogen(ElementForm::not(Element::H))]
#[case::non_hydrogen_set(ElementForm::lit_set([Element::C, Element::N]))]
fn test_hydrogen_folder_transform_element(mut methyl: Molecule, #[case] element: ElementForm) {
    methyl.atom_mut(AtomId(1)).attributes_mut().element = element.clone();
    let expected = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm {
            element,
            ..AtomForm::from("C#h3").into_concrete()
        }],
        ..Default::default()
    });
    assert_eq!(HydrogenFolder.transform(methyl), Ok(expected));
}

#[rstest]
#[case::set(NumForm::from(vec![1, 2]))]
#[case::variable(NumForm::var("v"))]
fn test_hydrogen_folder_transform_numeric_assertion(mut methyl: Molecule, #[case] value: NumForm) {
    methyl
        .atom_mut(AtomId(1))
        .attributes_mut()
        .constraints
        .set(AtomConstraintForm::Valence(value));
    let mut expected = mol_dsl_concrete!(r#"{:atoms ["C#h3"]}"#);
    expected
        .atom_mut(AtomId(0))
        .attributes_mut()
        .constraints
        .set(AtomConstraintForm::Valence(NumForm::Undetermined));
    assert_eq!(HydrogenFolder.transform(methyl), Ok(expected));
}

#[rstest]
#[case::atom(Constraint::Atom(AtomId(0), AtomConstraintForm::Degree(1.into())))]
#[case::bond(Constraint::Bond(BondId(0), BondConstraintForm::Aromatic(BooleanForm::Lit(false))))]
#[case::nested(Constraint::Not(Box::new(Constraint::Or(vec![
    Constraint::Atom(AtomId(0), AtomConstraintForm::Degree(1.into())),
]))))]
#[case::atom_subset(Constraint::Molecule(MoleculeConstraint::ChargeSum {
    atoms: Some(vec![AtomId(0)]), sum: 0.into(),
}))]
#[case::bond_subset(Constraint::Molecule(MoleculeConstraint::BondOrderSum {
    bonds: Some(vec![BondId(0)]), sum: 1.into(),
}))]
fn test_hydrogen_folder_transform_constraint_identity(
    methyl: Molecule,
    #[case] constraint: Constraint,
) {
    let mut editor = methyl.edit();
    editor.constraints_mut().push(constraint);
    let molecule = editor.finish().unwrap();
    assert_eq!(HydrogenFolder.transform(molecule.clone()), Ok(molecule));
}

#[rstest]
fn test_hydrogen_folder_transform_constraints(methyl: Molecule) {
    let constraints = vec![
        Constraint::Molecule(MoleculeConstraint::BondOrderSum {
            bonds: None,
            sum: 1.into(),
        }),
        Constraint::And(vec![Constraint::Atom(
            AtomId(1),
            AtomConstraintForm::Degree(1.into()),
        )]),
    ];
    let mut editor = methyl.edit();
    editor.constraints_mut().extend(constraints);
    let molecule = editor.finish().unwrap();
    let expected = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from("C#h3").into_concrete()],
        constraints: vec![
            Constraint::Molecule(MoleculeConstraint::BondOrderSum {
                bonds: None,
                sum: 1.into(),
            }),
            Constraint::And(vec![Constraint::Atom(
                AtomId(0),
                AtomConstraintForm::Degree(1.into()),
            )]),
        ]
        .into(),
        ..Default::default()
    });
    assert_eq!(HydrogenFolder.transform(molecule), Ok(expected));
}

#[rstest]
#[case::literal(StereoCoset::Lit(0), StereoCoset::Lit(1))]
#[case::set(StereoCoset::lit_set([0]), StereoCoset::lit_set([1]))]
#[case::term(StereoCoset::term(StereoTerm::Lit(0)), StereoCoset::term(StereoTerm::apply(
    StereoTerm::Lit(0), Permutation::from_image(&[1, 2, 3, 0]),
)))]
#[case::open(StereoCoset::Undetermined, StereoCoset::Undetermined)]
fn test_hydrogen_folder_transform_coset(#[case] input: StereoCoset, #[case] expected: StereoCoset) {
    let mut molecule = mol_dsl_concrete!(
        r#"{
        :atoms ["H" "C" "F" "Cl" "Br"]
        :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"] [1 4 "1"]]
    }"#
    );
    molecule
        .atom_mut(AtomId(1))
        .attributes_mut()
        .constraints
        .set(AtomConstraintForm::tetrahedral_stereo(
            TetrahedralStereoForm::Stereo(input),
        ));
    let mut product = mol_dsl_concrete!(
        r#"{
        :atoms ["C#h" "F" "Cl" "Br"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
    }"#
    );
    product
        .atom_mut(AtomId(0))
        .attributes_mut()
        .constraints
        .set(AtomConstraintForm::tetrahedral_stereo(
            TetrahedralStereoForm::Stereo(expected),
        ));
    assert_eq!(HydrogenFolder.transform(molecule), Ok(product));
}

#[rstest]
#[case::determined(
    StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
    StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0))
)]
#[case::open_coset(
    StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined),
    StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Undetermined)
)]
#[case::kindless(
    StereoConfigurationForm::Undetermined,
    StereoConfigurationForm::Undetermined
)]
fn test_hydrogen_folder_transform_stereo_attributes(
    #[case] atom_configuration: StereoConfigurationForm,
    #[case] bond_configuration: StereoConfigurationForm,
) {
    let mut input = mol_dsl_concrete!(
        r#"{
        :atoms ["H" "N#n" "F" "C#h" "Cl"]
        :bonds [[1 3 "2"] [0 1 "1"] [1 2 "1"] [3 4 "1"]]
        :stereo-atoms [{:site 1 :ligands [3 0 2 [:lp 1]] :attrs "Th0"}]
        :stereo-bonds [{:site 0 :ligands [0 2 4 [:h 3]] :attrs "Ct0"}]
    }"#
    );
    let mut expected = mol_dsl_concrete!(
        r#"{
        :atoms ["N#h#n" "F" "C#h" "Cl"]
        :bonds [[0 2 "2"] [0 1 "1"] [2 3 "1"]]
        :stereo-atoms [{:site 0 :ligands [2 [:h 0] 1 [:lp 0]] :attrs "Th0"}]
        :stereo-bonds [{:site 0 :ligands [[:h 0] 1 3 [:h 2]] :attrs "Ct0"}]
    }"#
    );
    let fluxionality = FluxionalityForm {
        permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])),
        active: BooleanForm::Lit(true),
    };
    for molecule in [&mut input, &mut expected] {
        let mut stereo = molecule.stereo_atom_mut(StereoAtomId(0));
        stereo.attributes_mut().configuration = atom_configuration.clone();
        stereo
            .attributes_mut()
            .constraints
            .set(StereoAtomConstraintForm::Fluxionality(fluxionality));
        let mut stereo = molecule.stereo_bond_mut(StereoBondId(0));
        stereo.attributes_mut().configuration = bond_configuration.clone();
        stereo
            .attributes_mut()
            .constraints
            .set(StereoBondConstraintForm::Fluxionality(fluxionality));
    }
    assert_eq!(HydrogenFolder.transform(input), Ok(expected));
}

#[rstest]
#[case::negative_after_eligible(
    mol_dsl_concrete!(r#"{
        :atoms ["C#h2" "H" "N#h-1" "H"]
        :bonds [[0 1 "1"] [2 3 "1"]]
    }"#),
    FoldHydrogensError::NegativeImplicitHydrogenCount { atom: AtomId(2), count: -1 },
)]
#[case::existing_implicit_hydrogen(
    mol_dsl_concrete!(r#"{
        :atoms ["C#h" "H" "F" "Cl"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]
    }"#),
    FoldHydrogensError::StereoHydrogenConflict { entity: Entity::StereoAtom(StereoAtomId(0)), hydrogen: AtomId(1) },
)]
#[case::existing_virtual_with_zero_count(
    mol_dsl_concrete!(r#"{
        :atoms ["C" "H" "F" "Cl"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]
    }"#),
    FoldHydrogensError::StereoHydrogenConflict { entity: Entity::StereoAtom(StereoAtomId(0)), hydrogen: AtomId(1) },
)]
#[case::two_explicit_hydrogens(
    mol_dsl_concrete!(r#"{
        :atoms ["C" "H" "H" "F" "Cl"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]
    }"#),
    FoldHydrogensError::StereoHydrogenConflict { entity: Entity::StereoAtom(StereoAtomId(0)), hydrogen: AtomId(1) },
)]
#[case::stereo_bond_endpoint(
    mol_dsl_concrete!(r#"{
        :atoms ["H" "C#h" "C#h" "F"] :bonds [[1 2 "2"] [0 1 "1"] [2 3 "1"]]
        :stereo-bonds [{:site 0 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct0"}]
    }"#),
    FoldHydrogensError::StereoHydrogenConflict { entity: Entity::StereoBond(StereoBondId(0)), hydrogen: AtomId(0) },
)]
fn test_hydrogen_folder_transform_error(
    #[case] input: Molecule,
    #[case] expected: FoldHydrogensError,
) {
    assert_eq!(
        HydrogenFolder.transform(input.clone()),
        Err(expected.clone())
    );
    let mut in_place = input.clone();
    assert_eq!(HydrogenFolder.transform_into(&mut in_place), Err(expected));
    assert_eq!(in_place, input);
    assert_eq!(HydrogenFolder.transform_iter(&input).next(), None);
}

#[rstest]
fn test_hydrogen_folder_transform_iter(methyl: Molecule) {
    let before = methyl.clone();
    let expected = mol_dsl_concrete!(r#"{:atoms ["C#h3"]}"#);
    let mut results = HydrogenFolder.transform_iter(&methyl);
    let mut product = results.next().unwrap();
    assert_eq!(product, expected);
    assert_eq!(results.next(), None);
    product.atom_mut(AtomId(0)).attributes_mut().charge = 1.into();
    assert_eq!(methyl, before);
}
