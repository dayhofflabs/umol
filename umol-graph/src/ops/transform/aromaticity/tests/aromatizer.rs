use rstest::rstest;
use umol_graph_core::{
    ConnectedComponentsAlgorithm, MaximumIndependentSetAlgorithm,
    RelevantCycleEnumerationAlgorithm, SimpleCycleEnumerationAlgorithm,
};
use umol_graph_ir::ir::{
    AromaticSystemId, AromaticValenceForm, AtomConstraintForm, AtomId, BondId, Molecule, NumForm,
    RingConfig,
};
use umol_graph_ir::mol_dsl_concrete;

use crate::ops::aromaticity::{
    AromaticityConfig, AromaticityContradiction, AromaticityInconsistency,
};
use crate::ops::invariant::ValenceMismatch;
use crate::ops::model::{AromaticityModel, AromaticityRule};
use crate::ops::transform::{AromatizeError, Aromatizer, Transformer};

#[rstest]
fn test_aromatizer_with_config() {
    let molecule = mol_dsl_concrete!(
        r#"{
        :atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"]
                [3 4 "1"] [4 5 "2"] [5 0 "1"]]
    }"#
    );
    let expected = Aromatizer::new(&AromaticityModel::daylight())
        .transform(molecule.clone())
        .unwrap();
    let configured = Aromatizer::with_config(
        &AromaticityModel::daylight(),
        AromaticityConfig {
            ring_config: RingConfig {
                simple_cycle_algorithm: SimpleCycleEnumerationAlgorithm::ReadTarjan,
                relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
            },
            connected_components_algorithm: ConnectedComponentsAlgorithm::Bfs,
            maximum_independent_set_algorithm: MaximumIndependentSetAlgorithm::BranchAndBound,
        },
    )
    .transform(molecule);
    assert_eq!(configured, Ok(expected));
}

#[rstest]
#[case::ring_double_bonds(
    mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#),
)]
#[case::nitrogen_donor(
    mol_dsl_concrete!(r#"{:atoms ["N#h#n" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"] [3 4 "2"] [4 0 "1"]]}"#),
    mol_dsl_concrete!(r#"{:atoms ["N#h" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#),
)]
#[case::oxygen_donor(
    mol_dsl_concrete!(r#"{:atoms ["O#n2" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"] [3 4 "2"] [4 0 "1"]]}"#),
    mol_dsl_concrete!(r#"{:atoms ["O#n" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#),
)]
#[case::sulfur_donor(
    mol_dsl_concrete!(r#"{:atoms ["S#n2" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"] [3 4 "2"] [4 0 "1"]]}"#),
    mol_dsl_concrete!(r#"{:atoms ["S#n" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#),
)]
#[case::anionic_carbon_donor(
    mol_dsl_concrete!(r#"{:atoms ["C#c-#h#n" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"] [3 4 "2"] [4 0 "1"]]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C#c-#h" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#),
)]
#[case::external_carbonyl(
    mol_dsl_concrete!(r#"{:atoms ["C" "C#h" "C#h" "C#h" "C#h" "C#h" "C#h" "O#n2"]
        :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"] [3 4 "2"] [4 5 "1"] [5 6 "2"] [6 0 "1"] [0 7 "2"]]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C" "C#h" "C#h" "C#h" "C#h" "C#h" "C#h" "O#n2"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 6 "1"] [6 0 "1"] [0 7 "2"]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5 6] :attrs "[0,1,1,1,1,1,1]"}]}"#),
)]
#[case::empty_carbon_orbital(
    mol_dsl_concrete!(r#"{:atoms ["C#c+#h" "C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"] [3 4 "2"] [4 5 "1"] [5 6 "2"] [6 0 "1"]]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C#c+#h" "C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 6 "1"] [6 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5 6] :attrs "[0,1,1,1,1,1,1]"}]}"#),
)]
#[case::aromatic_assertion(
    mol_dsl_concrete!(r#"{:atoms ["C#h#a#D2" "C#h" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C#h#D2" "C#h" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#),
)]
#[case::excluded_peripheral_bond(
    mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C" "C" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "2#a!"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]
                [5 6 "1"] [6 7 "2"] [7 8 "1"] [8 9 "2"] [9 4 "1"]]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C" "C" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "2#a!"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                [5 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 4 "1"]]
        :aromatic-systems [{:atoms [4 5 6 7 8 9] :attrs "[1,1,1,1,1,1]"}]}"#),
)]
#[case::nonliteral_aromatic_assertion(
    {
        let mut molecule = mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]}"#);
        molecule.atom_mut(AtomId(0)).attributes_mut().constraints.set(
            AtomConstraintForm::aromatic_valence(AromaticValenceForm::aromatic(NumForm::lit_set([0, 1])))
        );
        molecule
    },
    mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#),
)]
#[case::unrelated_undetermined_fields(
    {
        let mut molecule = mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h" "*"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]}"#);
        molecule.atom_mut(AtomId(0)).attributes_mut().implicit_hydrogens = NumForm::Undetermined;
        molecule.atom_mut(AtomId(0)).attributes_mut().lone_pairs = NumForm::Undetermined;
        molecule
    },
    {
        let mut molecule = mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h" "*"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#);
        molecule.atom_mut(AtomId(0)).attributes_mut().implicit_hydrogens = NumForm::Undetermined;
        molecule.atom_mut(AtomId(0)).attributes_mut().lone_pairs = NumForm::Undetermined;
        molecule
    },
)]
fn test_aromatizer_transform(#[case] molecule: Molecule, #[case] expected: Molecule) {
    let transformer = Aromatizer::new(&AromaticityModel::daylight());
    let mut in_place = molecule.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, expected);
    assert_eq!(
        transformer.transform_iter(&molecule).collect::<Vec<_>>(),
        vec![expected.clone()]
    );
    assert_eq!(transformer.transform(molecule), Ok(expected));
}

#[rstest]
#[case::acyclic(mol_dsl_concrete!(r#"{:atoms ["C#h4"]}"#), AromaticityModel::daylight())]
#[case::acyclic_donor(
    mol_dsl_concrete!(r#"{:atoms ["N#h3#n"]}"#),
    AromaticityModel { rule: AromaticityRule::Clar, ..AromaticityModel::daylight() },
)]
#[case::excluded_atom(mol_dsl_concrete!(r#"{
    :atoms ["C#h#a!" "C#h" "C#h" "C#h" "C#h" "C#h"]
    :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]
}"#), AromaticityModel::daylight())]
#[case::excluded_undetermined_atom(mol_dsl_concrete!(r#"{
    :atoms ["*#a!" "C#h" "C#h" "C#h" "C#h" "C#h"]
    :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]
}"#), AromaticityModel::daylight())]
#[case::saturated_ring(mol_dsl_concrete!(r#"{
    :atoms ["C#h2" "C#h2" "C#h2" "C#h2" "C#h2" "C#h2"]
    :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
}"#), AromaticityModel::daylight())]
#[case::unavailable_donor_pair(mol_dsl_concrete!(r#"{
    :atoms ["N#c+#h2" "C#h" "C#h" "C#h" "C#h"]
    :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"] [3 4 "2"] [4 0 "1"]]
}"#), AromaticityModel::daylight())]
#[case::cumulated_bonds(mol_dsl_concrete!(r#"{
    :atoms ["C" "C#h" "C#h2" "C#h2" "C#h2" "C#h"]
    :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "2"]]
}"#), AromaticityModel::daylight())]
#[case::already_aromatic(mol_dsl_concrete!(r#"{
    :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
    :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
            [3 4 :aromatic] [4 5 :aromatic] [5 0 :aromatic]]
    :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
}"#), AromaticityModel::daylight())]
fn test_aromatizer_transform_identity(#[case] molecule: Molecule, #[case] model: AromaticityModel) {
    let transformer = Aromatizer::new(&model);
    let mut in_place = molecule.clone();
    assert_eq!(transformer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, molecule);
    assert_eq!(transformer.transform(molecule.clone()), Ok(molecule));
}

#[rstest]
#[case::aromatic_contribution(
    mol_dsl_concrete!(r#"{:atoms ["C#h#a2" "C#h" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]}"#),
    AromaticityModel::daylight(),
    AromatizeError::Contradiction(AromaticityInconsistency::AromaticValenceMismatch {
        atom: AtomId(0), aromatic_system: AromaticSystemId(0),
    }.into()),
)]
#[case::later_aromatic_contribution(
    mol_dsl_concrete!(r#"{
        :atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"
                "C#h" "C#h" "C#h#a2" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]
                [6 7 "2"] [7 8 "1"] [8 9 "2"] [9 10 "1"] [10 11 "2"] [11 6 "1"]]
    }"#),
    AromaticityModel::daylight(),
    AromatizeError::Contradiction(AromaticityInconsistency::AromaticValenceMismatch {
        atom: AtomId(8), aromatic_system: AromaticSystemId(1),
    }.into()),
)]
#[case::nonliteral_aromatic_contribution(
    {
        let mut molecule = mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]}"#);
        molecule.atom_mut(AtomId(0)).attributes_mut().constraints.set(
            AtomConstraintForm::aromatic_valence(AromaticValenceForm::aromatic(NumForm::lit_set([0, 2])))
        );
        molecule
    },
    AromaticityModel::daylight(),
    AromatizeError::Contradiction(AromaticityInconsistency::AromaticValenceMismatch {
        atom: AtomId(0), aromatic_system: AromaticSystemId(0),
    }.into()),
)]
#[case::unclaimed_ring_atom(
    mol_dsl_concrete!(r#"{:atoms ["C#h2#a" "C#h2" "C#h2" "C#h2" "C#h2" "C#h2"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
    AromaticityModel::daylight(),
    AromatizeError::Contradiction(AromaticityInconsistency::AromaticValenceFailure {
        atom: AtomId(0),
    }.into()),
)]
#[case::unclaimed_ring_bond(
    mol_dsl_concrete!(r#"{:atoms ["C#h2" "C#h2" "C#h2" "C#h2" "C#h2" "C#h2"]
        :bonds [[0 1 "1#a"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
    AromaticityModel::daylight(),
    AromatizeError::Contradiction(AromaticityInconsistency::AromaticBondConstraintFailure {
        bond: BondId(0),
    }.into()),
)]
#[case::unclaimed_acyclic_atom(
    mol_dsl_concrete!(r#"{:atoms ["C#h4#a"]}"#),
    AromaticityModel::daylight(),
    AromatizeError::Contradiction(AromaticityInconsistency::AromaticValenceFailure {
        atom: AtomId(0),
    }.into()),
)]
#[case::unclaimed_acyclic_bond(
    mol_dsl_concrete!(r#"{:atoms ["C#h3" "C#h3"] :bonds [[0 1 "1#a"]]}"#),
    AromaticityModel::daylight(),
    AromatizeError::Contradiction(AromaticityInconsistency::AromaticBondConstraintFailure {
        bond: BondId(0),
    }.into()),
)]
#[case::conflicting_bond_evidence(
    mol_dsl_concrete!(r#"{:atoms ["C#h#a!" "C#h" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "2#a"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]}"#),
    AromaticityModel::daylight(),
    AromatizeError::Contradiction(AromaticityInconsistency::AromaticValenceBondMismatch {
        atom: AtomId(0), bond: BondId(0),
    }.into()),
)]
#[case::clar_heterocycle(
    mol_dsl_concrete!(r#"{
        :atoms ["N#n" "C#h" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]
    }"#),
    AromaticityModel { rule: AromaticityRule::Clar, ..AromaticityModel::daylight() },
    AromatizeError::Contradiction(AromaticityContradiction::ClarNonBenzenoid(
        "Clar model requires benzenoid input but non-carbon aromatic atoms are present".into()
    ))
)]
#[case::undetermined_donor_pair(
    {
        let mut molecule = mol_dsl_concrete!(r#"{:atoms ["N#h#n" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"] [3 4 "2"] [4 0 "1"]]}"#);
        molecule.atom_mut(AtomId(0)).attributes_mut().lone_pairs = NumForm::Undetermined;
        molecule
    },
    AromaticityModel::daylight(),
    AromatizeError::Underdetermined,
)]
#[case::undetermined_bond_order(
    mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "*"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]}"#),
    AromaticityModel::daylight(),
    AromatizeError::Underdetermined,
)]
#[case::contradicted_valence(
    mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h#v3" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]]}"#),
    AromaticityModel::daylight(),
    AromatizeError::ContradictedAssertion {
        aromatic_system: AromaticSystemId(0),
        atom: AtomId(2),
        constraint: AtomConstraintForm::valence(3),
    },
)]
#[case::later_aromatic_system_mismatch(
    mol_dsl_concrete!(r#"{
        :atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"
                "C#h" "C#h" "C#h#n" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [5 0 "1"]
                [6 7 "2"] [7 8 "1"] [8 9 "2"] [9 10 "1"] [10 11 "2"] [11 6 "1"]]
    }"#),
    AromaticityModel::daylight(),
    AromatizeError::ValenceInvariant {
        aromatic_system: AromaticSystemId(1),
        mismatch: ValenceMismatch::OrbitalCount {
            atom_id: AtomId(8), orbital_count: 10, electron_count: 8,
        },
    },
)]
fn test_aromatizer_transform_error(
    #[case] molecule: Molecule,
    #[case] model: AromaticityModel,
    #[case] expected: AromatizeError,
) {
    let transformer = Aromatizer::new(&model);
    let mut in_place = molecule.clone();
    assert_eq!(
        transformer.transform_into(&mut in_place),
        Err(expected.clone())
    );
    assert_eq!(in_place, molecule);
    assert_eq!(transformer.transform_iter(&molecule).next(), None);
    assert_eq!(transformer.transform(molecule), Err(expected));
}
