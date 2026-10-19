use pretty_assertions::assert_eq;
use rstest::{fixture, rstest};
use umol_graph_ir::ir::{
    AtomConstraintForm, AtomId, BondConstraintForm, BondId, BooleanForm, CisTransStereoForm,
    Constraint, Entity, FluxionalityForm, LigandPermutation, Molecule, MoleculeConstraint, NumForm,
    StereoAtomConstraintForm, StereoAtomId, StereoBondConstraintForm, StereoBondId,
    StereoConfigurationForm, StereoCoset, StereoKind, StereoTerm, TetrahedralStereoForm,
};
use umol_graph_ir::{mol_dsl, mol_dsl_concrete};
use umol_perm::Permutation;

use crate::ops::transform::{HydrogenUnfolder, Transformer, UnfoldHydrogensError};

#[fixture]
fn methyl() -> Molecule {
    mol_dsl_concrete!(r#"{:atoms ["C#h3"]}"#)
}

#[fixture]
fn stereo() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["N#h#n#T1" "F" "C#h" "Cl"]
        :bonds [[0 2 "2#C0"] [0 1 "1"] [2 3 "1"]]
        :stereo-atoms [{:site 0 :ligands [2 [:h 0] 1 [:lp 0]] :attrs "Th1"}]
        :stereo-bonds [{:site 0 :ligands [3 [:h 2] [:h 0] 1] :attrs "Ct0"}]
    }"#
    )
}

#[fixture]
fn explicit_stereo() -> Molecule {
    mol_dsl_concrete!(
        r#"{
        :atoms ["N#n#T1" "F" "C" "Cl" "H" "H"]
        :bonds [[0 2 "2#C0"] [0 1 "1"] [2 3 "1"] [0 4 "1"] [2 5 "1"]]
        :stereo-atoms [{:site 0 :ligands [2 4 1 [:lp 0]] :attrs "Th1"}]
        :stereo-bonds [{:site 0 :ligands [3 5 4 1] :attrs "Ct0"}]
    }"#
    )
}

#[rstest]
#[case::one(
    mol_dsl_concrete!(r#"{:atoms ["O#h"]}"#),
    mol_dsl_concrete!(r#"{:atoms ["O" "H"] :bonds [[0 1 "1"]]}"#),
)]
#[case::several(
    mol_dsl_concrete!(r#"{:atoms ["C#h4"]}"#),
    mol_dsl_concrete!(r#"{
        :atoms ["C" "H" "H" "H" "H"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    }"#),
)]
#[case::several_bearers(
    mol_dsl_concrete!(r#"{:atoms ["C#h3" "O#h"] :bonds [[0 1 "1"]]}"#),
    mol_dsl_concrete!(r#"{
        :atoms ["C" "O" "H" "H" "H" "H"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"]]
    }"#),
)]
#[case::mixed(
    mol_dsl_concrete!(r#"{
        :atoms ["C#h2" "H#i2" "H"] :bonds [[0 1 "1"] [0 2 "1"]]
    }"#),
    mol_dsl_concrete!(r#"{
        :atoms ["C" "H#i2" "H" "H" "H"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    }"#),
)]
#[case::constraints(
    mol_dsl_concrete!(r#"{
        :atoms ["C#h2#v1#D1#V3#X3#H2#R!" "F"] :bonds [[0 1 "1"]]
    }"#),
    mol_dsl_concrete!(r#"{
        :atoms ["C#v3#D3#V3#X3#H2#R!" "F" "H" "H"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
    }"#),
)]
#[case::other_fields(
    mol_dsl_concrete!(r#"{:atoms ["N#i15#c+#h2#n#u#s2"]}"#),
    mol_dsl_concrete!(r#"{
        :atoms ["N#i15#c+#n#u#s2" "H" "H"] :bonds [[0 1 "1"] [0 2 "1"]]
    }"#),
)]
#[case::partial(
    mol_dsl!(r#"{:atoms ["*#h2"]}"#),
    mol_dsl!(r#"{
        :atoms ["*#h0" "H#i=#c0#h0#n0#u0#s" "H#i=#c0#h0#n0#u0#s"]
        :bonds [[0 1 "1#c0#u0#s"] [0 2 "1#c0#u0#s"]]
    }"#),
)]
#[case::tetrahedral(
    mol_dsl_concrete!(r#"{
        :atoms ["C#h#T1" "F" "Cl" "Br"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
        :stereo-atoms [{:site 0 :ligands [2 [:h 0] 3 1] :attrs "Th1"}]
    }"#),
    mol_dsl_concrete!(r#"{
        :atoms ["C#T1" "F" "Cl" "Br" "H"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [2 4 3 1] :attrs "Th1"}]
    }"#),
)]
#[case::assertions_without_overlays(
    mol_dsl_concrete!(r#"{
        :atoms ["N#h#n#T1" "F" "C#h" "Cl"]
        :bonds [[0 2 "2#C0"] [0 1 "1"] [2 3 "1"]]
    }"#),
    mol_dsl_concrete!(r#"{
        :atoms ["N#n#T1" "F" "C" "Cl" "H" "H"]
        :bonds [[0 2 "2#C0"] [0 1 "1"] [2 3 "1"] [0 4 "1"] [2 5 "1"]]
    }"#),
)]
#[case::shared_bearer(stereo(), explicit_stereo())]
#[case::retained_relations(
    mol_dsl_concrete!(r#"{
        :atoms ["O#h" "C" "C" "C"] :bonds [[1 2 "1"] [2 3 "1"] [1 3 "1"]]
        :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]
        :aromatic-systems [{:atoms [1 2 3] :attrs "[1,1,1]"}]
        :multicenter-bonds [{:atoms [0 2 3] :attrs "[2,0,0]"}]
        :noncovalent-bonds [{:atoms [0 1] :attrs "*"}]
    }"#),
    mol_dsl_concrete!(r#"{
        :atoms ["O" "C" "C" "C" "H"]
        :bonds [[1 2 "1"] [2 3 "1"] [1 3 "1"] [0 4 "1"]]
        :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]
        :aromatic-systems [{:atoms [1 2 3] :attrs "[1,1,1]"}]
        :multicenter-bonds [{:atoms [0 2 3] :attrs "[2,0,0]"}]
        :noncovalent-bonds [{:atoms [0 1] :attrs "*"}]
    }"#),
)]
fn test_hydrogen_unfolder_transform(#[case] input: Molecule, #[case] expected: Molecule) {
    assert_eq!(
        HydrogenUnfolder.transform(input.clone()),
        Ok(expected.clone())
    );
    let mut in_place = input;
    assert_eq!(HydrogenUnfolder.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected);
    assert_eq!(HydrogenUnfolder.transform(in_place), Ok(expected.clone()));
    assert_eq!(expected.clone().edit().finish(), Ok(expected));
}

#[rstest]
#[case::empty(Molecule::new())]
#[case::zero(mol_dsl_concrete!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#))]
#[case::explicit(mol_dsl_concrete!(r#"{
    :atoms ["C" "H" "H#i2"] :bonds [[0 1 "1"] [0 2 "1"]]
}"#))]
#[case::stereo(explicit_stereo())]
#[case::nonliteral_stereo_without_virtual(mol_dsl!(r#"{
    :atoms ["C" "F" "Cl" "Br" "I"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]
}"#))]
fn test_hydrogen_unfolder_transform_identity(#[case] input: Molecule) {
    assert_eq!(HydrogenUnfolder.transform(input.clone()), Ok(input.clone()));
    let mut in_place = input.clone();
    assert_eq!(HydrogenUnfolder.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, input);
}

#[rstest]
#[case::unknown(NumForm::Undetermined)]
#[case::set(NumForm::from(vec![0, 1]))]
#[case::variable(NumForm::var("h"))]
fn test_hydrogen_unfolder_transform_partial(#[case] count: NumForm) {
    let mut input = mol_dsl_concrete!(r#"{:atoms ["O#h" "C"] :bonds [[0 1 "1"]]}"#);
    input
        .atom_mut(AtomId(1))
        .attributes_mut()
        .implicit_hydrogens = count.clone();
    let mut expected = mol_dsl_concrete!(
        r#"{
        :atoms ["O" "C" "H"] :bonds [[0 1 "1"] [0 2 "1"]]
    }"#
    );
    expected
        .atom_mut(AtomId(1))
        .attributes_mut()
        .implicit_hydrogens = count;
    assert_eq!(HydrogenUnfolder.transform(input), Ok(expected));
}

#[rstest]
#[case::unknown(NumForm::Undetermined)]
#[case::set(NumForm::from(vec![0, 1]))]
#[case::variable(NumForm::var("v"))]
fn test_hydrogen_unfolder_transform_numeric_assertions(
    mut methyl: Molecule,
    #[case] value: NumForm,
) {
    let mut expected = mol_dsl_concrete!(
        r#"{
        :atoms ["C" "H" "H" "H"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
    }"#
    );
    for assertion in [
        AtomConstraintForm::Valence(value.clone()),
        AtomConstraintForm::Degree(value),
    ] {
        methyl
            .atom_mut(AtomId(0))
            .attributes_mut()
            .constraints
            .set(assertion);
    }
    for assertion in [
        AtomConstraintForm::Valence(NumForm::Undetermined),
        AtomConstraintForm::Degree(NumForm::Undetermined),
    ] {
        expected
            .atom_mut(AtomId(0))
            .attributes_mut()
            .constraints
            .set(assertion);
    }
    assert_eq!(HydrogenUnfolder.transform(methyl), Ok(expected));
}

#[rstest]
fn test_hydrogen_unfolder_transform_constraints(methyl: Molecule) {
    let mut editor = methyl.edit();
    editor.constraints_mut().extend([
        Constraint::Molecule(MoleculeConstraint::BondOrderSum {
            bonds: None,
            sum: 0.into(),
        }),
        Constraint::Not(Box::new(Constraint::Atom(
            AtomId(0),
            AtomConstraintForm::Degree(1.into()),
        ))),
    ]);
    let input = editor.finish().unwrap();
    let mut editor = mol_dsl_concrete!(
        r#"{
        :atoms ["C" "H" "H" "H"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
    }"#
    )
    .edit();
    editor.constraints_mut().extend([
        Constraint::Molecule(MoleculeConstraint::BondOrderSum {
            bonds: None,
            sum: 0.into(),
        }),
        Constraint::Not(Box::new(Constraint::Atom(
            AtomId(0),
            AtomConstraintForm::Degree(1.into()),
        ))),
    ]);
    assert_eq!(
        HydrogenUnfolder.transform(input),
        Ok(editor.finish().unwrap())
    );
}

#[rstest]
#[case::literal(StereoCoset::Lit(1))]
#[case::set(StereoCoset::lit_set([0]))]
#[case::term(StereoCoset::term(StereoTerm::Lit(0)))]
#[case::open(StereoCoset::Undetermined)]
fn test_hydrogen_unfolder_transform_coset(
    mut stereo: Molecule,
    mut explicit_stereo: Molecule,
    #[case] coset: StereoCoset,
) {
    for molecule in [&mut stereo, &mut explicit_stereo] {
        molecule
            .atom_mut(AtomId(0))
            .attributes_mut()
            .constraints
            .set(AtomConstraintForm::tetrahedral_stereo(
                TetrahedralStereoForm::Stereo(coset.clone()),
            ));
        molecule
            .bond_mut(BondId(0))
            .attributes_mut()
            .constraints
            .set(BondConstraintForm::cis_trans_stereo(
                CisTransStereoForm::Stereo(coset.clone()),
            ));
    }
    assert_eq!(HydrogenUnfolder.transform(stereo), Ok(explicit_stereo));
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
fn test_hydrogen_unfolder_transform_stereo_attributes(
    mut stereo: Molecule,
    mut explicit_stereo: Molecule,
    #[case] atom_configuration: StereoConfigurationForm,
    #[case] bond_configuration: StereoConfigurationForm,
) {
    let fluxionality = FluxionalityForm {
        permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])),
        active: BooleanForm::Lit(true),
    };
    for molecule in [&mut stereo, &mut explicit_stereo] {
        let mut atom = molecule.stereo_atom_mut(StereoAtomId(0));
        atom.attributes_mut().configuration = atom_configuration.clone();
        atom.attributes_mut()
            .constraints
            .set(StereoAtomConstraintForm::Fluxionality(fluxionality));
        let mut bond = molecule.stereo_bond_mut(StereoBondId(0));
        bond.attributes_mut().configuration = bond_configuration.clone();
        bond.attributes_mut()
            .constraints
            .set(StereoBondConstraintForm::Fluxionality(fluxionality));
    }
    assert_eq!(HydrogenUnfolder.transform(stereo), Ok(explicit_stereo));
}

#[rstest]
#[case::negative(mol_dsl_concrete!(r#"{:atoms ["C#h-1"]}"#),
    UnfoldHydrogensError::NegativeImplicitHydrogenCount { atom: AtomId(0), count: -1 })]
#[case::negative_after_applicable(mol_dsl_concrete!(r#"{:atoms ["C#h3" "O#h-1"]}"#),
    UnfoldHydrogensError::NegativeImplicitHydrogenCount { atom: AtomId(1), count: -1 })]
#[case::zero_with_virtual(mol_dsl_concrete!(r#"{
    :atoms ["C" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]
}"#), UnfoldHydrogensError::StereoHydrogenCountMismatch {
    entity: Entity::StereoAtom(StereoAtomId(0)), atom: AtomId(0),
    implicit_hydrogens: 0, virtual_hydrogens: 1,
})]
#[case::unknown_with_virtual(mol_dsl!(r#"{
    :atoms ["C" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]
}"#), UnfoldHydrogensError::NonConcreteStereoHydrogenCount {
    entity: Entity::StereoAtom(StereoAtomId(0)), atom: AtomId(0), value: NumForm::Undetermined,
})]
#[case::atom_missing_virtual(mol_dsl_concrete!(r#"{
    :atoms ["C#h" "F" "Cl" "Br" "I"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]
}"#), UnfoldHydrogensError::StereoHydrogenCountMismatch {
    entity: Entity::StereoAtom(StereoAtomId(0)), atom: AtomId(0),
    implicit_hydrogens: 1, virtual_hydrogens: 0,
})]
#[case::atom_several_without_virtual(mol_dsl_concrete!(r#"{
    :atoms ["C#h2" "F" "Cl" "Br" "I"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]
}"#), UnfoldHydrogensError::StereoHydrogenCountMismatch {
    entity: Entity::StereoAtom(StereoAtomId(0)), atom: AtomId(0),
    implicit_hydrogens: 2, virtual_hydrogens: 0,
})]
#[case::bond_missing_virtual(mol_dsl_concrete!(r#"{
    :atoms ["C#h" "C" "F" "Cl" "Br" "I"]
    :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
    :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct0"}]
}"#), UnfoldHydrogensError::StereoHydrogenCountMismatch {
    entity: Entity::StereoBond(StereoBondId(0)), atom: AtomId(0),
    implicit_hydrogens: 1, virtual_hydrogens: 0,
})]
fn test_hydrogen_unfolder_transform_error(
    #[case] input: Molecule,
    #[case] expected: UnfoldHydrogensError,
) {
    assert_eq!(
        HydrogenUnfolder.transform(input.clone()),
        Err(expected.clone())
    );
    let mut in_place = input.clone();
    assert_eq!(
        HydrogenUnfolder.transform_into(&mut in_place),
        Err(expected)
    );
    assert_eq!(in_place, input);
    assert_eq!(HydrogenUnfolder.transform_iter(&input).next(), None);
}

#[rstest]
#[case::zero(NumForm::Lit(0))]
#[case::several(NumForm::Lit(2))]
#[case::unknown(NumForm::Undetermined)]
#[case::set(NumForm::from(vec![0, 1]))]
#[case::variable(NumForm::var("h"))]
fn test_hydrogen_unfolder_transform_stereo_error(
    stereo: Molecule,
    #[case] value: NumForm,
    #[values(AtomId(0), AtomId(2))] atom: AtomId,
) {
    let entity = if atom == AtomId(0) {
        Entity::StereoAtom(StereoAtomId(0))
    } else {
        Entity::StereoBond(StereoBondId(0))
    };
    let expected = match &value {
        NumForm::Lit(count) => UnfoldHydrogensError::StereoHydrogenCountMismatch {
            entity,
            atom,
            implicit_hydrogens: *count,
            virtual_hydrogens: 1,
        },
        _ => UnfoldHydrogensError::NonConcreteStereoHydrogenCount {
            entity,
            atom,
            value: value.clone(),
        },
    };
    let mut input = stereo;
    input.atom_mut(atom).attributes_mut().implicit_hydrogens = value;
    assert_eq!(
        HydrogenUnfolder.transform(input.clone()),
        Err(expected.clone())
    );
    let mut in_place = input.clone();
    assert_eq!(
        HydrogenUnfolder.transform_into(&mut in_place),
        Err(expected)
    );
    assert_eq!(in_place, input);
    assert_eq!(HydrogenUnfolder.transform_iter(&input).next(), None);
}

#[rstest]
fn test_hydrogen_unfolder_transform_iter(methyl: Molecule) {
    let before = methyl.clone();
    let expected = mol_dsl_concrete!(
        r#"{
        :atoms ["C" "H" "H" "H"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
    }"#
    );
    let mut results = HydrogenUnfolder.transform_iter(&methyl);
    let mut product = results.next().unwrap();
    assert_eq!(product, expected);
    assert_eq!(results.next(), None);
    product.atom_mut(AtomId(1)).attributes_mut().charge = 1.into();
    assert_eq!(methyl, before);
}
