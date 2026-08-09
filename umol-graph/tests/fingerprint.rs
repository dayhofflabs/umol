use rstest::{fixture, rstest};
use umol_graph::fingerprint::{
    featurize_reaction, EcfpFeaturizer, Featurizer, FingerprintError, MorganFeaturizer,
    PatternFingerprinter, ReactionCombinator, ReactionFingerprint, ReactionSide,
    SubstructureFeaturizer, WlFeaturizer,
};
use umol_graph::ingest::ingest_smiles;
use umol_graph_core::{
    RefinementRounds, RelevantCycleEnumerationAlgorithm, SimpleCycleEnumerationAlgorithm,
};
use umol_graph_ir::ir::{
    AtomDelta, AtomFieldChange, AtomId, BondDelta, BondId, Delta, Deltas, Molecule, NumForm,
    ReactionAst, RingConfig,
};
use umol_graph_ir::{mol_dsl, mol_dsl_ground};

#[fixture]
fn ethanol() -> Molecule {
    ingest_smiles("CCO").unwrap()
}

#[fixture]
fn benzene() -> Molecule {
    ingest_smiles("c1ccccc1").unwrap()
}

#[fixture]
fn ethanol_deoxygenation(ethanol: Molecule) -> ReactionAst {
    let oxygen = ethanol.atom(AtomId(2)).ast.clone();
    let bond = ethanol.bond(BondId(1)).ast.clone();
    ReactionAst::new(
        ethanol,
        Deltas::from_iter([
            Delta::Atom(AtomDelta::Remove {
                id: AtomId(2),
                ast: oxygen,
            }),
            Delta::Bond(BondDelta::Remove {
                id: BondId(1),
                atoms: [AtomId(1), AtomId(2)],
                ast: bond,
            }),
        ]),
    )
}

#[rstest]
fn test_wl_featurizer_featurize(ethanol: Molecule) {
    let fingerprint = WlFeaturizer::new(RefinementRounds::Fixed(3))
        .featurize(&ethanol)
        .unwrap();
    assert_eq!(
        fingerprint.ids(),
        &[
            2520347590860685079,
            3352603313223549703,
            4152249898001161146,
            5715207763479934940,
            5807737097854608645,
            7542810387455301591,
            11457795998246593156,
            11986000156817227245,
            12895020514073294021,
            13932567567828606490,
            17305796300852423160,
            17417400371411086222,
        ]
    );
}

#[rstest]
fn test_wl_featurizer_featurize_error() {
    assert_eq!(
        WlFeaturizer::new(RefinementRounds::Fixed(3))
            .featurize(&mol_dsl!(r#"{:atoms ["C"] :bonds []}"#))
            .unwrap_err(),
        FingerprintError::NotGround
    );
}

#[rstest]
fn test_wl_featurizer_featurize_counted_error() {
    assert_eq!(
        WlFeaturizer::new(RefinementRounds::Fixed(3))
            .featurize_counted(&mol_dsl!(r#"{:atoms ["C"] :bonds []}"#))
            .unwrap_err(),
        FingerprintError::NotGround
    );
}

#[rstest]
fn test_ecfp_featurizer_featurize(benzene: Molecule) {
    assert_eq!(
        EcfpFeaturizer {
            radius: 2,
            hashing_scheme: Default::default(),
            ring_config: RingConfig {
                simple_cycle_algorithm: SimpleCycleEnumerationAlgorithm::ReadTarjan,
                relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
            },
        }
        .featurize(&benzene)
        .unwrap()
        .ids(),
        &[
            3716727142329830942,
            7364724293986779056,
            16614949630732484927,
        ]
    );
}

#[rstest]
fn test_ecfp_featurizer_featurize_error() {
    assert_eq!(
        EcfpFeaturizer::new(2)
            .featurize(&mol_dsl!(r#"{:atoms ["C"] :bonds []}"#))
            .unwrap_err(),
        FingerprintError::NotGround
    );
}

#[rstest]
fn test_ecfp_featurizer_featurize_counted_error() {
    assert_eq!(
        EcfpFeaturizer::new(2)
            .featurize_counted(&mol_dsl!(r#"{:atoms ["C"] :bonds []}"#))
            .unwrap_err(),
        FingerprintError::NotGround
    );
}

#[rstest]
fn test_morgan_featurizer_featurize(benzene: Molecule) {
    assert_eq!(
        MorganFeaturizer {
            radius: 2,
            ring_config: RingConfig {
                simple_cycle_algorithm: SimpleCycleEnumerationAlgorithm::ReadTarjan,
                relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
            },
        }
        .featurize(&benzene)
        .unwrap()
        .ids(),
        &[98513984, 2763854213, 3218693969]
    );
}

#[rstest]
fn test_morgan_featurizer_featurize_error() {
    assert_eq!(
        MorganFeaturizer::new(2)
            .featurize(&mol_dsl!(r#"{:atoms ["C"] :bonds []}"#))
            .unwrap_err(),
        FingerprintError::NotGround
    );
}

#[rstest]
fn test_morgan_featurizer_featurize_counted_error() {
    assert_eq!(
        MorganFeaturizer::new(2)
            .featurize_counted(&mol_dsl!(r#"{:atoms ["C"] :bonds []}"#))
            .unwrap_err(),
        FingerprintError::NotGround
    );
}

#[rstest]
#[case::wl(Featurizer::Wl(WlFeaturizer::new(RefinementRounds::Fixed(3))))]
#[case::ecfp(Featurizer::Ecfp(EcfpFeaturizer::new(2)))]
#[case::morgan(Featurizer::Morgan(MorganFeaturizer::new(2)))]
fn test_featurizer_featurize_error(#[case] featurizer: Featurizer) {
    assert_eq!(
        featurizer
            .featurize(&mol_dsl!(r#"{:atoms ["C"] :bonds []}"#))
            .unwrap_err(),
        FingerprintError::NotGround
    );
}

#[rstest]
#[case::wl(Featurizer::Wl(WlFeaturizer::new(RefinementRounds::Fixed(3))))]
#[case::ecfp(Featurizer::Ecfp(EcfpFeaturizer::new(2)))]
#[case::morgan(Featurizer::Morgan(MorganFeaturizer::new(2)))]
fn test_featurizer_featurize_counted(ethanol: Molecule, #[case] featurizer: Featurizer) {
    let binary = featurizer.featurize(&ethanol).unwrap();
    let counted = featurizer.featurize_counted(&ethanol).unwrap();
    assert_eq!(
        counted
            .entries()
            .iter()
            .map(|(identifier, _)| *identifier)
            .collect::<Vec<_>>(),
        binary.ids()
    );
}

#[rstest]
#[case::wl(Featurizer::Wl(WlFeaturizer::new(RefinementRounds::Fixed(3))))]
#[case::ecfp(Featurizer::Ecfp(EcfpFeaturizer::new(2)))]
#[case::morgan(Featurizer::Morgan(MorganFeaturizer::new(2)))]
fn test_featurizer_featurize_counted_error(#[case] featurizer: Featurizer) {
    assert_eq!(
        featurizer
            .featurize_counted(&mol_dsl!(r#"{:atoms ["C"] :bonds []}"#))
            .unwrap_err(),
        FingerprintError::NotGround
    );
}

#[rstest]
fn test_pattern_fingerprinter_fingerprint(ethanol: Molecule) {
    let fingerprint = PatternFingerprinter::new().fingerprint(&ethanol).unwrap();
    assert_eq!(
        (0..fingerprint.width())
            .filter(|&bit| fingerprint.get(bit) == Some(true))
            .collect::<Vec<_>>(),
        vec![54, 173, 217, 429, 622, 759, 778, 874, 946, 967, 1022, 1033, 1061, 1236, 1289, 1295]
    );
}

#[rstest]
fn test_pattern_fingerprinter_fingerprint_error() {
    assert_eq!(
        PatternFingerprinter::new()
            .fingerprint(&mol_dsl!(r#"{:atoms ["C"] :bonds []}"#))
            .unwrap_err(),
        FingerprintError::NotGround
    );
}

#[rstest]
fn test_substructure_featurizer_featurize(ethanol: Molecule) {
    assert_eq!(
        SubstructureFeaturizer::new(2)
            .featurize(&ethanol)
            .unwrap()
            .ids(),
        &[
            vec![1, 0, 0, 0, 5, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0],
            vec![1, 0, 0, 0, 5, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0],
            vec![
                3, 0, 0, 0, 5, 0, 0, 0, 0, 6, 0, 0, 0, 5, 0, 0, 0, 0, 6, 0, 0, 0, 3, 0, 0, 0, 1, 1,
                0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0,
            ],
            vec![
                3, 0, 0, 0, 5, 0, 0, 0, 0, 6, 0, 0, 0, 5, 0, 0, 0, 0, 8, 0, 0, 0, 3, 0, 0, 0, 1, 1,
                0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0,
            ],
            vec![
                5, 0, 0, 0, 5, 0, 0, 0, 0, 6, 0, 0, 0, 5, 0, 0, 0, 0, 6, 0, 0, 0, 5, 0, 0, 0, 0, 8,
                0, 0, 0, 3, 0, 0, 0, 1, 1, 0, 3, 0, 0, 0, 1, 1, 0, 4, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0,
                0, 1, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 4, 0, 0, 0, 2, 0, 0, 0, 4, 0, 0, 0,
            ],
        ]
    );
}

#[rstest]
fn test_substructure_featurizer_featurize_error() {
    assert_eq!(
        SubstructureFeaturizer::new(2)
            .featurize(&mol_dsl!(r#"{:atoms ["C"] :bonds []}"#))
            .unwrap_err(),
        FingerprintError::NotGround
    );
}

#[rstest]
fn test_featurize_reaction_difference(ethanol_deoxygenation: ReactionAst) {
    let fingerprint = featurize_reaction(
        &ethanol_deoxygenation,
        &Featurizer::Morgan(MorganFeaturizer {
            radius: 1,
            ring_config: RingConfig {
                simple_cycle_algorithm: SimpleCycleEnumerationAlgorithm::ReadTarjan,
                relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
            },
        }),
        ReactionCombinator::Difference,
    )
    .unwrap();
    match fingerprint {
        ReactionFingerprint::Difference(features) => assert_eq!(
            features.entries(),
            &[
                (864662311, -1),
                (1535166686, -1),
                (2245384272, -1),
                (2246997334, 1),
                (3542456614, -1),
                (3548082732, 1),
                (4018048386, -1),
            ]
        ),
        other => panic!("expected Difference, got {other:?}"),
    }
}

#[rstest]
#[case::non_ground(
    ReactionAst::new(mol_dsl!(r#"{:atoms ["C"] :bonds []}"#), Deltas::new()),
    FingerprintError::NotGround
)]
#[case::inconsistent(
    ReactionAst::new(
        mol_dsl_ground!(r#"{:atoms ["C #h4"] :bonds []}"#),
        Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(0),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(1),
                new: NumForm::Lit(0),
            },
        })]),
    ),
    FingerprintError::Inconsistent
)]
fn test_featurize_reaction_error(
    #[case] reaction: ReactionAst,
    #[case] expected: FingerprintError,
) {
    assert_eq!(
        featurize_reaction(
            &reaction,
            &Featurizer::Morgan(MorganFeaturizer::new(1)),
            ReactionCombinator::Difference,
        )
        .unwrap_err(),
        expected
    );
}

#[rstest]
fn test_featurize_reaction_disjoint_union(ethanol_deoxygenation: ReactionAst) {
    let fingerprint = featurize_reaction(
        &ethanol_deoxygenation,
        &Featurizer::Morgan(MorganFeaturizer::new(1)),
        ReactionCombinator::DisjointUnion,
    )
    .unwrap();
    match fingerprint {
        ReactionFingerprint::DisjointUnion(features) => assert_eq!(
            features.ids(),
            &[
                (ReactionSide::Reactant, 864662311),
                (ReactionSide::Reactant, 1535166686),
                (ReactionSide::Reactant, 2245384272),
                (ReactionSide::Reactant, 2246728737),
                (ReactionSide::Reactant, 3542456614),
                (ReactionSide::Reactant, 4018048386),
                (ReactionSide::Product, 2246728737),
                (ReactionSide::Product, 2246997334),
                (ReactionSide::Product, 3548082732),
            ]
        ),
        other => panic!("expected DisjointUnion, got {other:?}"),
    }
}
