use pretty_assertions::assert_eq;
use rstest::rstest;
use umol_graph_ir::ir::{AromaticSystemId, AtomId, Molecule, NumForm};
use umol_graph_ir::{mol_dsl, mol_dsl_concrete};
use umol_utils::solution::Solution;

use crate::ops::transform::{ChargeDelocalizer, Transformer};
use crate::ops::validate::{SpinInvariantsValidator, ValenceInvariantsValidator};

#[rstest]
#[case::cyclopentadienyl(
    mol_dsl_concrete!(r#"{
        :atoms ["C#c-#h#a2" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
        :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                [3 4 :aromatic] [4 0 :aromatic]]
        :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]
    }"#),
    mol_dsl_concrete!(r#"{
        :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
        :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                [3 4 :aromatic] [4 0 :aromatic]]
        :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[1,1,1,1,1]#c-"}]
    }"#)
)]
#[case::tropylium(
    mol_dsl_concrete!(r#"{
        :atoms ["C#c+#h#a0" "C#h#a" "C#h#a" "C#h#a"
                "C#h#a" "C#h#a" "C#h#a"]
        :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                [3 4 :aromatic] [4 5 :aromatic] [5 6 :aromatic]
                [6 0 :aromatic]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5 6]
                            :attrs "[0,1,1,1,1,1,1]"}]
    }"#),
    mol_dsl_concrete!(r#"{
        :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a"
                "C#h#a" "C#h#a" "C#h#a"]
        :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                [3 4 :aromatic] [4 5 :aromatic] [5 6 :aromatic]
                [6 0 :aromatic]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5 6]
                            :attrs "[1,1,1,1,1,1,1]#c+"}]
    }"#)
)]
#[case::multiple_aromatic_systems(
    mol_dsl_concrete!(r#"{
        :atoms ["C#c+#h#a0" "C#h#a" "C#h#a"
                "C#c+#h#a0" "C#h#a" "C#h#a"]
        :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 0 :aromatic]
                [3 4 :aromatic] [4 5 :aromatic] [5 3 :aromatic]]
        :aromatic-systems [
            {:atoms [0 1 2] :attrs "[0,1,1]"}
            {:atoms [3 4 5] :attrs "[0,1,1]"}]
    }"#),
    mol_dsl_concrete!(r#"{
        :atoms ["C#h#a" "C#h#a" "C#h#a"
                "C#h#a" "C#h#a" "C#h#a"]
        :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 0 :aromatic]
                [3 4 :aromatic] [4 5 :aromatic] [5 3 :aromatic]]
        :aromatic-systems [
            {:atoms [0 1 2] :attrs "[1,1,1]#c+"}
            {:atoms [3 4 5] :attrs "[1,1,1]#c+"}]
    }"#)
)]
#[case::radical_anion(
    mol_dsl_concrete!(r#"{
        :atoms ["C#c-#h0#u1#s2#a2" "C#h1#a1" "C#h1#a1"
                "C#h1#a1" "C#h1#a1" "C#h1#a1"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]
                [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[2,1,1,1,1,1]"}]
    }"#),
    mol_dsl_concrete!(r#"{
        :atoms ["C#h0#u1#s2#a1" "C#h1#a1" "C#h1#a1"
                "C#h1#a1" "C#h1#a1" "C#h1#a1"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]
                [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]#c-"}]
    }"#)
)]
#[case::radical_cation(
    mol_dsl_concrete!(r#"{
        :atoms ["C#c+#h0#u1#s2#a0" "C#h1#a1" "C#h1#a1"
                "C#h1#a1" "C#h1#a1" "C#h1#a1"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]
                [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[0,1,1,1,1,1]"}]
    }"#),
    mol_dsl_concrete!(r#"{
        :atoms ["C#h0#u1#s2#a1" "C#h1#a1" "C#h1#a1"
                "C#h1#a1" "C#h1#a1" "C#h1#a1"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]
                [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]#c+"}]
    }"#)
)]
#[case::mixed_charges(
    mol_dsl_concrete!(r#"{
        :atoms ["C#c+#h1#a0" "C#c-#h1#a2" "C#h1#a1"
                "C#h1#a1" "C#h1#a1" "C#h1#a1"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]
                [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[0,2,1,1,1,1]"}]
    }"#),
    mol_dsl_concrete!(r#"{
        :atoms ["C#h1#a1" "C#h1#a1" "C#h1#a1"
                "C#h1#a1" "C#h1#a1" "C#h1#a1"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]
                [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
    }"#)
)]
#[case::existing_system_charge(
    mol_dsl_concrete!(r#"{
        :atoms ["C#c+#h1#a0" "C#h1#a1" "C#h1#a1"
                "C#h1#a1" "C#h1#a1" "C#h1#a1"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]
                [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[0,1,1,1,1,1]#c-"}]
    }"#),
    mol_dsl_concrete!(r#"{
        :atoms ["C#h1#a1" "C#h1#a1" "C#h1#a1"
                "C#h1#a1" "C#h1#a1" "C#h1#a1"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]
                [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
    }"#)
)]
fn test_charge_delocalizer_transform(#[case] input: Molecule, #[case] expected: Molecule) {
    let mut borrowed = input.clone();
    assert_eq!(ChargeDelocalizer.transform_into(&mut borrowed), Ok(()));
    assert_eq!(borrowed, expected);
    for molecule in [&input, &borrowed] {
        assert_eq!(
            ValenceInvariantsValidator.validate(molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(
            SpinInvariantsValidator.validate(molecule),
            Ok(Solution::Determined(()))
        );
    }
    for system in input.aromatic_systems().iter() {
        assert_eq!(
            borrowed.aromatic_system(system.id()).electron_count()
                - borrowed.aromatic_system(system.id()).charge(),
            system.electron_count() - system.charge()
        );
    }
    assert_eq!(ChargeDelocalizer.transform(input), Ok(expected));
}

#[rstest]
#[case::open_occupations(
    mol_dsl!(r#"{
        :atoms ["C#c-#a2" "C#c0#a1" "C#c0#a1"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2] :attrs "[2,1,1]#c0"}]
    }"#),
    mol_dsl!(r#"{
        :atoms ["C#c0#a1" "C#c0#a1" "C#c0#a1"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,1]#c-"}]
    }"#)
)]
fn test_charge_delocalizer_transform_partial(#[case] input: Molecule, #[case] expected: Molecule) {
    let mut borrowed = input.clone();
    assert_eq!(ChargeDelocalizer.transform_into(&mut borrowed), Ok(()));
    assert_eq!(borrowed, expected);
    assert_eq!(ChargeDelocalizer.transform(input), Ok(expected));
}

#[rstest]
#[case::short_counts(mol_dsl!(r#"{:atoms ["C" "C" "C"] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,2]"}]}"#))]
#[case::long_counts(mol_dsl!(r#"{:atoms ["C" "C" "C"] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,2,0,5]"}]}"#))]
#[case::already_delocalized(mol_dsl_concrete!(r#"{
    :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
    :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
            [3 4 :aromatic] [4 0 :aromatic]]
    :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[1,1,1,1,1]#c-"}]
}"#))]
#[case::localized_radical(mol_dsl_concrete!(r#"{
    :atoms ["C#h0#u1#s2#a1" "C#h1#a1" "C#h1#a1"
            "C#h1#a1" "C#h1#a1" "C#h1#a1"]
    :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]
            [3 4 "1"] [4 5 "1"] [5 0 "1"]]
    :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
}"#))]
#[case::undetermined_atom_charge({
    let mut molecule = mol_dsl_concrete!(r#"{
        :atoms ["C#c-#h1#a2" "C#h1#a1" "C#h1#a1"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2] :attrs "[2,1,1]"}]
    }"#);
    molecule.atom_mut(AtomId(1)).attributes_mut().charge = NumForm::Undetermined;
    molecule
})]
#[case::undetermined_system_charge({
    let mut molecule = mol_dsl_concrete!(r#"{
        :atoms ["C#c-#h1#a2" "C#h1#a1" "C#h1#a1"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2] :attrs "[2,1,1]"}]
    }"#);
    molecule.aromatic_system_mut(AromaticSystemId(0)).attributes_mut().charge = NumForm::Undetermined;
    molecule
})]
#[case::heterogeneous(mol_dsl_concrete!(r#"{
    :atoms ["B#c-#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
    :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
            [3 4 :aromatic] [4 5 :aromatic] [5 0 :aromatic]]
    :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
}"#))]
#[case::non_literal(mol_dsl!(r#"{
    :atoms ["C" "C" "C"]
    :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 0 :aromatic]]
    :aromatic-systems [{:atoms [0 1 2] :attrs "*"}]
}"#))]
fn test_charge_delocalizer_transform_identity(#[case] input: Molecule) {
    let mut borrowed = input.clone();
    assert_eq!(ChargeDelocalizer.transform_into(&mut borrowed), Ok(()));
    assert_eq!(borrowed, input);
    assert_eq!(ChargeDelocalizer.transform(input.clone()), Ok(input));
}

#[rstest]
#[case(
    mol_dsl_concrete!(r#"{
        :atoms ["C#c+#h#a0" "C#h#a" "C#h#a"]
        :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 0 :aromatic]]
        :aromatic-systems [{:atoms [0 1 2] :attrs "[0,1,1]"}]
    }"#),
    mol_dsl_concrete!(r#"{
        :atoms ["C#h#a" "C#h#a" "C#h#a"]
        :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 0 :aromatic]]
        :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,1]#c+"}]
    }"#)
)]
fn test_charge_delocalizer_transform_iter(#[case] input: Molecule, #[case] expected: Molecule) {
    let original = input.clone();
    let mut results = ChargeDelocalizer.transform_iter(&input);
    let mut first = results.next().unwrap();
    assert_eq!(first, expected);
    assert_eq!(results.next(), None);
    assert_eq!(results.next(), None);

    first.atom_mut(AtomId(0)).attributes_mut().charge = NumForm::Lit(-1);
    assert_eq!(input, original);
    assert_eq!(
        ChargeDelocalizer.transform_iter(&input).next(),
        Some(expected)
    );
}
