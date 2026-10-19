use rstest::rstest;
use umol_graph_ir::ir::{
    AromaticSystemId, AtomConstraintForm, AtomConstraintKey, AtomId, BondId, Molecule, NumForm,
};
use umol_graph_ir::{mol_dsl, mol_dsl_concrete};

use crate::ops::invariant::ValenceMismatch;
use crate::ops::transform::aromaticity::{
    AromaticSystemPlan, ExposedAtom, MatchingInput, MatchingInputMode,
};
use crate::ops::transform::{
    KekulizeConfig, KekulizeError, Kekulizer, MaximumMatchingAlgorithm, Transformer,
};

#[rstest]
#[case::neutral_all_covered(
    mol_dsl_concrete!(r#"{:atoms ["C" "C" "C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1 2 3] :attrs "[1,1,1,1]"}]}"#),
    MatchingInput {
        required_covered: vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)],
        required_exposed: vec![],
        exposed_count: 0,
        mode: MatchingInputMode::Prescribed,
    }
)]
#[case::prescribed_donor(
    mol_dsl_concrete!(r#"{:atoms ["N" "C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1 2] :attrs "[2,1,1]"}]}"#),
    MatchingInput {
        required_covered: vec![AtomId(1), AtomId(2)],
        required_exposed: vec![ExposedAtom { atom: AtomId(0), electrons: 2 }],
        exposed_count: 1,
        mode: MatchingInputMode::Prescribed,
    }
)]
#[case::prescribed_acceptor_at_positional_atom(
    mol_dsl_concrete!(r#"{:atoms ["C" "B" "C"] :bonds [] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,0,1]"}]}"#),
    MatchingInput {
        required_covered: vec![AtomId(0), AtomId(2)],
        required_exposed: vec![ExposedAtom { atom: AtomId(1), electrons: 0 }],
        exposed_count: 1,
        mode: MatchingInputMode::Prescribed,
    }
)]
#[case::multiple_prescribed_exposures(
    mol_dsl_concrete!(r#"{:atoms ["B" "C" "N"] :bonds [] :aromatic-systems [{:atoms [0 1 2] :attrs "[0,1,2]"}]}"#),
    MatchingInput {
        required_covered: vec![AtomId(1)],
        required_exposed: vec![ExposedAtom { atom: AtomId(0), electrons: 0 }, ExposedAtom { atom: AtomId(2), electrons: 2 }],
        exposed_count: 2,
        mode: MatchingInputMode::Prescribed,
    }
)]
#[case::mobile_positive_charge(
    mol_dsl_concrete!(r#"{:atoms ["C" "C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,1]#c+"}]}"#),
    MatchingInput {
        required_covered: vec![],
        required_exposed: vec![],
        exposed_count: 1,
        mode: MatchingInputMode::OneMobileExposure,
    }
)]
#[case::mobile_negative_charge(
    mol_dsl_concrete!(r#"{:atoms ["C" "C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,1]#c-"}]}"#),
    MatchingInput {
        required_covered: vec![],
        required_exposed: vec![],
        exposed_count: 1,
        mode: MatchingInputMode::OneMobileExposure,
    }
)]
fn test_matching_input_from_aromatic_system(
    #[case] input: Molecule,
    #[case] expected: MatchingInput,
) {
    let aromatic_system = input.aromatic_systems().ids().next().unwrap();
    assert_eq!(
        MatchingInput::from_aromatic_system(&input, aromatic_system),
        Ok(expected)
    );
}

#[rstest]
#[case::short_counts(
    mol_dsl!(r#"{:atoms ["C" "C" "C"] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,2]"}]}"#),
    KekulizeError::ElectronCountMismatch { aromatic_system: AromaticSystemId(0), member_count: 3, electron_count: 2 }
)]
#[case::long_counts(
    mol_dsl!(r#"{:atoms ["C" "C" "C"] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,2,0,5]"}]}"#),
    KekulizeError::ElectronCountMismatch { aromatic_system: AromaticSystemId(0), member_count: 3, electron_count: 4 }
)]
#[case::undetermined_electrons(
    mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :attrs "*#c0#u0#s1"}]}"#),
    KekulizeError::UndeterminedElectrons(AromaticSystemId(0))
)]
#[case::undetermined_charge(
    mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :attrs "[1,1]#u0#s1"}]}"#),
    KekulizeError::UndeterminedCharge(AromaticSystemId(0))
)]
#[case::undetermined_unpaired_electrons(
    mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :attrs "[1,1]#c0#u0"}]}"#),
    KekulizeError::UndeterminedUnpairedElectrons(AromaticSystemId(0))
)]
#[case::unsupported_contribution_at_positional_atom(
    mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :attrs "[1,3]#c0#u0#s1"}]}"#),
    KekulizeError::UnsupportedElectronContribution {
        aromatic_system: AromaticSystemId(0),
        atom: AtomId(1),
        contribution: 3,
    }
)]
#[case::open_shell(
    mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :attrs "[1,1]#c0#u1#s2"}]}"#),
    KekulizeError::OpenShell(AromaticSystemId(0))
)]
#[case::unsupported_charge(
    mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :attrs "[1,1]#c+2#u0#s1"}]}"#),
    KekulizeError::UnsupportedCharge {
        aromatic_system: AromaticSystemId(0),
        charge: 2,
    }
)]
#[case::mixed_prescribed_and_mobile_demand(
    mol_dsl!(r#"{:atoms ["B" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :attrs "[0,1]#c+#u0#s1"}]}"#),
    KekulizeError::MixedExposureDemand(AromaticSystemId(0))
)]
fn test_matching_input_from_aromatic_system_error(
    #[case] input: Molecule,
    #[case] expected: KekulizeError,
) {
    let aromatic_system = input.aromatic_systems().ids().next().unwrap();
    assert_eq!(
        MatchingInput::from_aromatic_system(&input, aromatic_system),
        Err(expected)
    );
}

#[rstest]
#[case::edmonds(
    MaximumMatchingAlgorithm::Edmonds,
    KekulizeConfig {
        algorithm: MaximumMatchingAlgorithm::Edmonds,
    }
)]
#[case::hopcroft_karp(
    MaximumMatchingAlgorithm::HopcroftKarp,
    KekulizeConfig {
        algorithm: MaximumMatchingAlgorithm::HopcroftKarp,
    }
)]
fn test_kekulize_config_new(
    #[case] algorithm: MaximumMatchingAlgorithm,
    #[case] expected: KekulizeConfig,
) {
    assert_eq!(KekulizeConfig::new(algorithm), expected);
}

#[rstest]
fn test_kekulize_config_default() {
    assert_eq!(
        KekulizeConfig::default(),
        KekulizeConfig {
            algorithm: MaximumMatchingAlgorithm::Edmonds,
        }
    );
}

#[rstest]
#[case::benzene(
    mol_dsl_concrete!(r#"{:atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [0 5 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#),
    (0..6).map(AtomId).collect(),
    mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"] :bonds [[0 1 :double] [1 2 :single] [2 3 :double] [3 4 :single] [4 5 :double] [0 5 :single]]}"#)
)]
#[case::embedded_benzene_nonidentity_correspondence(
    mol_dsl_concrete!(r#"{:atoms ["O#h#n2" "H" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"] :bonds [[0 1 :single] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [5 6 :aromatic] [6 7 :aromatic] [2 7 :aromatic]] :aromatic-systems [{:atoms [2 3 4 5 6 7] :attrs "[1,1,1,1,1,1]"}]}"#),
    vec![AtomId(0), AtomId(2), AtomId(3), AtomId(4), AtomId(5), AtomId(6), AtomId(7), AtomId(1)],
    mol_dsl_concrete!(r#"{:atoms ["O#h#n2" "H" "C#h" "C#h" "C#h" "C#h" "C#h" "C#h"] :bonds [[0 1 :single] [2 3 :double] [3 4 :single] [4 5 :double] [5 6 :single] [6 7 :double] [2 7 :single]]}"#)
)]
#[case::multiple_disjoint_aromatic_systems(
    mol_dsl_concrete!(r#"{:atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 0 :aromatic] [4 5 :aromatic] [5 6 :aromatic] [6 7 :aromatic] [7 4 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3] :attrs "[1,1,1,1]"} {:atoms [4 5 6 7] :attrs "[1,1,1,1]"}]}"#),
    vec![AtomId(0), AtomId(4), AtomId(1), AtomId(5), AtomId(2), AtomId(6), AtomId(3), AtomId(7)],
    mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h" "C#h" "C#h"] :bonds [[0 1 :double] [1 2 :single] [2 3 :double] [3 0 :single] [4 5 :double] [5 6 :single] [6 7 :double] [7 4 :single]]}"#)
)]
#[case::unchanged_spin(
    mol_dsl_concrete!(r#"{:atoms ["N#a#u2#s2" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
        :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [0 5 :aromatic]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#),
    (0..6).map(AtomId).collect(),
    mol_dsl_concrete!(r#"{:atoms ["N#u2#s2" "C#h" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [0 5 "1"]]}"#)
)]
#[case::unknown_external_order(
    mol_dsl_concrete!(r#"{:atoms ["C#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "*"]
        :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [0 5 :aromatic] [0 6 "*"]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#),
    (0..7).map(AtomId).collect(),
    mol_dsl_concrete!(r#"{:atoms ["C" "C#h" "C#h" "C#h" "C#h" "C#h" "*"]
        :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [0 5 "1"] [0 6 "*"]]}"#)
)]
#[case::retained_incidence(
    mol_dsl_concrete!(r#"{:atoms ["N#a#v3#d1" "C#h#a" "C#h#a" "N#a#m2" "C#h#a" "C#h#a" "B#h3" "B#h3"]
        :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [0 5 :aromatic]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        :dative-bonds [{:donors [0] :acceptor 6 :attrs "1"}]
        :multicenter-bonds [{:atoms [3 7] :attrs "[2,0]"}]}"#),
    (0..8).map(AtomId).collect(),
    mol_dsl_concrete!(r#"{:atoms ["N#v3#d1" "C#h" "C#h" "N#m2" "C#h" "C#h" "B#h3" "B#h3"]
        :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [0 5 "1"]]
        :dative-bonds [{:donors [0] :acceptor 6 :attrs "1"}]
        :multicenter-bonds [{:atoms [3 7] :attrs "[2,0]"}]}"#)
)]
#[case::normalized_external_order(
    {
        let mut input = mol_dsl_concrete!(r#"{:atoms ["C#a#v4" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "F#n3"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [0 5 :aromatic] [0 6 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#);
        input.bond_mut(BondId(6)).attributes_mut().order = NumForm::lit_set([1]);
        input
    },
    (0..7).map(AtomId).collect(),
    {
        let mut expected = mol_dsl_concrete!(r#"{:atoms ["C#v4" "C#h" "C#h" "C#h" "C#h" "C#h" "F#n3"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [0 5 "1"] [0 6 "1"]]}"#);
        expected.bond_mut(BondId(6)).attributes_mut().order = NumForm::lit_set([1]);
        expected
    }
)]
fn test_kekulizer_transform_into(
    #[case] input: Molecule,
    #[case] node_order: Vec<AtomId>,
    #[case] expected: Molecule,
) {
    let transformer = Kekulizer::new(KekulizeConfig::default(), node_order);
    assert_eq!(transformer.transform(input.clone()), Ok(expected.clone()));
    let mut molecule = input;
    transformer.transform_into(&mut molecule).unwrap();
    assert_eq!(molecule, expected);
}

#[rstest]
#[case::kekule_benzene( mol_dsl_concrete!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 :double] [1 2 :single] [2 3 :double] [3 4 :single] [4 5 :double] [0 5 :single]]}"#))]
fn test_kekulizer_transform_into_identity(#[case] input: Molecule) {
    let transformer = Kekulizer::new(KekulizeConfig::default(), (0..6).map(AtomId).collect());
    assert_eq!(transformer.transform(input.clone()), Ok(input.clone()));
    let mut molecule = input.clone();
    transformer.transform_into(&mut molecule).unwrap();
    assert_eq!(molecule, input);
}

#[rstest]
#[case::no_matching(
    mol_dsl_concrete!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [0 4 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[1,1,1,1,1]"}]}"#),
    (0..5).map(AtomId).collect(),
    KekulizeError::NoMatching(AromaticSystemId(0))
)]
#[case::missing_aromatic_system_atom(
    mol_dsl_concrete!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 0 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3] :attrs "[1,1,1,1]"}]}"#),
    vec![AtomId(0), AtomId(1), AtomId(2)],
    KekulizeError::InvalidNodeOrder {
        aromatic_system: AromaticSystemId(0),
        missing: vec![AtomId(3)],
        duplicates: vec![],
    }
)]
#[case::duplicate_aromatic_system_atom(
    mol_dsl_concrete!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 0 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3] :attrs "[1,1,1,1]"}]}"#),
    vec![AtomId(0), AtomId(1), AtomId(1), AtomId(2), AtomId(3)],
    KekulizeError::InvalidNodeOrder {
        aromatic_system: AromaticSystemId(0),
        missing: vec![],
        duplicates: vec![AtomId(1)],
    }
)]
#[case::contradicted_valence(
    mol_dsl_concrete!(r#"{:atoms ["C#h#a" "C#h#a" "C#h#a#v2" "C#h#a" "C#h#a" "C#h#a"]
        :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [0 5 :aromatic]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#),
    (0..6).map(AtomId).collect(),
    KekulizeError::ContradictedAssertion {
        aromatic_system: AromaticSystemId(0),
        atom: AtomId(2),
        constraint: AtomConstraintForm::valence(2),
    }
)]
#[case::later_aromatic_system_mismatch(
    mol_dsl_concrete!(r#"{:atoms ["*#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"
                                "C#h#a" "C#h#a" "C#h#a#n" "C#h#a" "C#h#a" "C#h#a"]
        :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [0 5 :aromatic]
                [6 7 :aromatic] [7 8 :aromatic] [8 9 :aromatic] [9 10 :aromatic] [10 11 :aromatic] [6 11 :aromatic]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}
                           {:atoms [6 7 8 9 10 11] :attrs "[1,1,1,1,1,1]"}]}"#),
    (0..12).map(AtomId).collect(),
    KekulizeError::ValenceInvariant {
        aromatic_system: AromaticSystemId(1),
        mismatch: ValenceMismatch::OrbitalCount {
            atom_id: AtomId(8), orbital_count: 10, electron_count: 8,
        },
    }
)]
#[case::contradicted_donation(
    mol_dsl_concrete!(r#"{:atoms ["C#h#a#d1" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
        :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [0 5 :aromatic]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#),
    (0..6).map(AtomId).collect(),
    KekulizeError::ContradictedAssertion {
        aromatic_system: AromaticSystemId(0),
        atom: AtomId(0),
        constraint: AtomConstraintForm::DonatedPairs(NumForm::Lit(1)),
    }
)]
#[case::contradictory_external_order(
    {
        let mut input = mol_dsl_concrete!(r#"{:atoms ["C#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "F#n3"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [0 5 :aromatic] [0 6 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#);
        input.bond_mut(BondId(6)).attributes_mut().order = NumForm::lit_set([]);
        input
    },
    (0..7).map(AtomId).collect(),
    KekulizeError::ValenceInvariant {
        aromatic_system: AromaticSystemId(0),
        mismatch: ValenceMismatch::ContradictoryConstraint {
            atom_id: AtomId(0), key: AtomConstraintKey::Valence,
        },
    }
)]
fn test_kekulizer_transform_into_error(
    #[case] input: Molecule,
    #[case] node_order: Vec<AtomId>,
    #[case] expected: KekulizeError,
) {
    let transformer = Kekulizer::new(KekulizeConfig::default(), node_order);
    assert_eq!(transformer.transform(input.clone()), Err(expected.clone()));
    let mut molecule = input.clone();
    let result = transformer.transform_into(&mut molecule);
    assert_eq!(result, Err(expected));
    assert_eq!(molecule, input);
}

#[rstest]
#[case::non_bipartite(
    mol_dsl_concrete!(r#"{:atoms ["C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [0 2 :aromatic]] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,1]"}]}"#),
    (0..3).map(AtomId).collect(),
    KekulizeError::NonBipartiteMatching(AromaticSystemId(0)),
)]
#[case::mobile_exposure(
    mol_dsl_concrete!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 0 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[1,1,1,1,1]#c-"}]}"#),
    (0..5).map(AtomId).collect(),
    KekulizeError::NonBipartiteMatching(AromaticSystemId(0)),
)]
fn test_kekulizer_transform_into_matching_error(
    #[case] input: Molecule,
    #[case] node_order: Vec<AtomId>,
    #[case] expected: KekulizeError,
) {
    let transformer = Kekulizer::new(
        KekulizeConfig::new(MaximumMatchingAlgorithm::HopcroftKarp),
        node_order,
    );
    assert_eq!(transformer.transform(input.clone()), Err(expected.clone()));
    let mut actual = input.clone();
    let result = transformer.transform_into(&mut actual);

    assert_eq!(result, Err(expected));
    assert_eq!(actual, input);
}

#[rstest]
#[case::benzene(
    mol_dsl_concrete!(r#"{:atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [0 5 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#),
    mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"] :bonds [[0 1 :double] [1 2 :single] [2 3 :double] [3 4 :single] [4 5 :double] [0 5 :single]]}"#)
)]
fn test_kekulizer_transform_iter(#[case] input: Molecule, #[case] expected: Molecule) {
    let original = input.clone();
    let transformer = Kekulizer::new(KekulizeConfig::default(), input.atoms().ids().collect());
    let mut results = transformer.transform_iter(&input);
    assert_eq!(results.next(), Some(expected));
    assert_eq!(results.next(), None);
    assert_eq!(results.next(), None);
    assert_eq!(input, original);
}

#[rstest]
#[case::no_matching(
    mol_dsl_concrete!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [0 4 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[1,1,1,1,1]"}]}"#),
    KekulizeError::NoMatching(AromaticSystemId(0))
)]
fn test_kekulizer_transform_iter_error(#[case] input: Molecule, #[case] expected: KekulizeError) {
    let original = input.clone();
    let transformer = Kekulizer::new(KekulizeConfig::default(), input.atoms().ids().collect());
    {
        let mut results = transformer.transform_iter(&input);
        assert_eq!(results.next(), None);
        assert_eq!(results.next(), None);
    }
    assert_eq!(input, original);
    assert_eq!(transformer.transform(input), Err(expected));
}

#[rstest]
#[case::five_membered_ring(
    mol_dsl_concrete!(r#"{:atoms ["N#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 0 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#),
    (0..5).map(AtomId).collect(),
    vec![AromaticSystemPlan {
        aromatic_system_idx: AromaticSystemId(0),
        atoms: (0..5).map(AtomId).collect(),
        exposed_atoms: vec![ExposedAtom { atom: AtomId(0), electrons: 2 }],
        mobile_charge: None,
        matched_bonds: vec![BondId(1), BondId(3)],
        unmatched_bonds: vec![BondId(0), BondId(4), BondId(2)],
    }]
)]
#[case::seven_membered_ring(
    mol_dsl_concrete!(r#"{:atoms ["B#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [5 6 :aromatic] [6 0 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4 5 6] :attrs "[0,1,1,1,1,1,1]"}]}"#),
    (0..7).map(AtomId).collect(),
    vec![AromaticSystemPlan {
        aromatic_system_idx: AromaticSystemId(0),
        atoms: (0..7).map(AtomId).collect(),
        exposed_atoms: vec![ExposedAtom { atom: AtomId(0), electrons: 0 }],
        mobile_charge: None,
        matched_bonds: vec![BondId(1), BondId(3), BondId(5)],
        unmatched_bonds: vec![BondId(0), BondId(6), BondId(2), BondId(4)],
    }]
)]
#[case::fused_heterocycle(
    mol_dsl_concrete!(r#"{:atoms ["N#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 0 :aromatic] [4 5 :aromatic] [5 6 :aromatic] [6 7 :aromatic] [7 8 :aromatic] [8 3 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4 5 6 7 8] :attrs "[2,1,1,1,1,1,1,1,1]"}]}"#),
    (0..9).map(AtomId).collect(),
    vec![AromaticSystemPlan {
        aromatic_system_idx: AromaticSystemId(0),
        atoms: (0..9).map(AtomId).collect(),
        exposed_atoms: vec![ExposedAtom { atom: AtomId(0), electrons: 2 }],
        mobile_charge: None,
        matched_bonds: vec![BondId(1), BondId(3), BondId(6), BondId(8)],
        unmatched_bonds: vec![BondId(0), BondId(4), BondId(2), BondId(9), BondId(5), BondId(7)],
    }]
)]
#[case::nontrivial_host_ids(
    mol_dsl_concrete!(r#"{:atoms ["O" "H" "N#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :single] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [5 6 :aromatic] [6 2 :aromatic]] :aromatic-systems [{:atoms [2 3 4 5 6] :attrs "[2,1,1,1,1]"}]}"#),
    vec![AtomId(0), AtomId(2), AtomId(3), AtomId(4), AtomId(5), AtomId(6), AtomId(1)],
    vec![AromaticSystemPlan {
        aromatic_system_idx: AromaticSystemId(0),
        atoms: (2..7).map(AtomId).collect(),
        exposed_atoms: vec![ExposedAtom { atom: AtomId(2), electrons: 2 }],
        mobile_charge: None,
        matched_bonds: vec![BondId(2), BondId(4)],
        unmatched_bonds: vec![BondId(1), BondId(5), BondId(3)],
    }]
)]
#[case::five_membered_mobile_hole_natural_order(
    mol_dsl_concrete!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 0 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[1,1,1,1,1]#c-"}]}"#),
    (0..5).map(AtomId).collect(),
    vec![AromaticSystemPlan {
        aromatic_system_idx: AromaticSystemId(0),
        atoms: (0..5).map(AtomId).collect(),
        exposed_atoms: vec![ExposedAtom { atom: AtomId(4), electrons: 2 }],
        mobile_charge: Some(-1),
        matched_bonds: vec![BondId(0), BondId(2)],
        unmatched_bonds: vec![BondId(4), BondId(1), BondId(3)],
    }]
)]
#[case::five_membered_mobile_hole_rotated_order(
    mol_dsl_concrete!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 0 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[1,1,1,1,1]#c-"}]}"#),
    vec![AtomId(2), AtomId(3), AtomId(4), AtomId(0), AtomId(1)],
    vec![AromaticSystemPlan {
        aromatic_system_idx: AromaticSystemId(0),
        atoms: (0..5).map(AtomId).collect(),
        exposed_atoms: vec![ExposedAtom { atom: AtomId(0), electrons: 2 }],
        mobile_charge: Some(-1),
        matched_bonds: vec![BondId(1), BondId(3)],
        unmatched_bonds: vec![BondId(0), BondId(4), BondId(2)],
    }]
)]
#[case::seven_membered_mobile_hole_natural_order(
    mol_dsl_concrete!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [5 6 :aromatic] [6 0 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4 5 6] :attrs "[1,1,1,1,1,1,1]#c+"}]}"#),
    (0..7).map(AtomId).collect(),
    vec![AromaticSystemPlan {
        aromatic_system_idx: AromaticSystemId(0),
        atoms: (0..7).map(AtomId).collect(),
        exposed_atoms: vec![ExposedAtom { atom: AtomId(6), electrons: 0 }],
        mobile_charge: Some(1),
        matched_bonds: vec![BondId(0), BondId(2), BondId(4)],
        unmatched_bonds: vec![BondId(6), BondId(1), BondId(3), BondId(5)],
    }]
)]
#[case::seven_membered_mobile_hole_rotated_order(
    mol_dsl_concrete!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [5 6 :aromatic] [6 0 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4 5 6] :attrs "[1,1,1,1,1,1,1]#c+"}]}"#),
    vec![AtomId(3), AtomId(4), AtomId(5), AtomId(6), AtomId(0), AtomId(1), AtomId(2)],
    vec![AromaticSystemPlan {
        aromatic_system_idx: AromaticSystemId(0),
        atoms: (0..7).map(AtomId).collect(),
        exposed_atoms: vec![ExposedAtom { atom: AtomId(1), electrons: 0 }],
        mobile_charge: Some(1),
        matched_bonds: vec![BondId(6), BondId(2), BondId(4)],
        unmatched_bonds: vec![BondId(0), BondId(1), BondId(3), BondId(5)],
    }]
)]
fn test_kekulizer_plan_aromatic_systems(
    #[case] input: Molecule,
    #[case] node_order: Vec<AtomId>,
    #[case] mut expected: Vec<AromaticSystemPlan>,
) {
    let kekulizer = Kekulizer::new(KekulizeConfig::default(), node_order);
    let mut actual = kekulizer.plan_aromatic_systems(&input).unwrap();
    assert_eq!(kekulizer.plan_aromatic_systems(&input), Ok(actual.clone()));
    for plan in actual.iter_mut().chain(&mut expected) {
        plan.matched_bonds.sort_unstable();
        plan.unmatched_bonds.sort_unstable();
    }
    assert_eq!(actual, expected);
}

#[rstest]
#[case::impossible_prescribed_hole(
    mol_dsl_concrete!(r#"{:atoms ["N#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 0 :aromatic] [0 4 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#),
    (0..5).map(AtomId).collect(),
    KekulizeError::NoMatching(AromaticSystemId(0))
)]
#[case::mobile_matching_deficiency(
    mol_dsl_concrete!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[1,1,1,1,1]#c-"}]}"#),
    (0..5).map(AtomId).collect(),
    KekulizeError::MatchingDeficiency {
        aromatic_system: AromaticSystemId(0),
        expected: 1,
        actual: 3,
    }
)]
fn test_kekulizer_plan_aromatic_systems_error(
    #[case] input: Molecule,
    #[case] node_order: Vec<AtomId>,
    #[case] expected: KekulizeError,
) {
    let kekulizer = Kekulizer::new(KekulizeConfig::default(), node_order);
    assert_eq!(kekulizer.plan_aromatic_systems(&input), Err(expected));
}
