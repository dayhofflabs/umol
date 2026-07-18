use rstest::{fixture, rstest};
use umol_ast::ast::{
    AtomDelta, AtomId, BondDelta, BondId, Delta, Deltas, MoleculeAst, ReactionAst,
};
use umol_graph::fingerprint::{
    featurize_reaction, EcfpFeaturizer, Featurizer, MorganFeaturizer, PatternFingerprinter,
    ReactionCombinator, ReactionFingerprint, Side, SubstructureFeaturizer, WlFeaturizer,
};
use umol_graph::hash::RefinementXxh3Scheme;
use umol_graph::ingest::smiles;
use umol_graph_core::RefinementRounds;

#[fixture]
fn ethanol() -> MoleculeAst {
    smiles("CCO").unwrap()
}

#[fixture]
fn ethanol_deoxygenation(ethanol: MoleculeAst) -> ReactionAst {
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
fn test_wl_featurizer_featurize(ethanol: MoleculeAst) {
    let fingerprint = WlFeaturizer {
        rounds: RefinementRounds::Fixed(3),
        scheme: RefinementXxh3Scheme::albatross(),
    }
    .featurize(&ethanol);
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
fn test_ecfp_featurizer_featurize(ethanol: MoleculeAst) {
    assert_eq!(
        EcfpFeaturizer::new(2).featurize(&ethanol).ids(),
        &[
            63839236075656913,
            1189585227353469813,
            3822471596818936039,
            13652293261850732425,
            15001976065402722634,
            16149328945726899460,
        ]
    );
}

#[rstest]
fn test_morgan_featurizer_featurize(ethanol: MoleculeAst) {
    assert_eq!(
        MorganFeaturizer::new(2).featurize(&ethanol).ids(),
        &[864662311, 1535166686, 2245384272, 2246728737, 3542456614, 4018048386]
    );
}

#[rstest]
fn test_pattern_fingerprinter_fingerprint(ethanol: MoleculeAst) {
    let fingerprint = PatternFingerprinter::new().fingerprint(&ethanol).unwrap();
    assert_eq!(
        (0..fingerprint.width())
            .filter(|&bit| fingerprint.get(bit))
            .collect::<Vec<_>>(),
        vec![54, 173, 217, 429, 622, 759, 778, 874, 946, 967, 1022, 1033, 1061, 1236, 1289, 1295]
    );
}

#[rstest]
fn test_substructure_featurizer_featurize(ethanol: MoleculeAst) {
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
fn test_featurize_reaction_difference(ethanol_deoxygenation: ReactionAst) {
    let fingerprint = featurize_reaction(
        &ethanol_deoxygenation,
        &Featurizer::Morgan(MorganFeaturizer::new(1)),
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
                (Side::Reactant, 864662311),
                (Side::Reactant, 1535166686),
                (Side::Reactant, 2245384272),
                (Side::Reactant, 2246728737),
                (Side::Reactant, 3542456614),
                (Side::Reactant, 4018048386),
                (Side::Product, 2246728737),
                (Side::Product, 2246997334),
                (Side::Product, 3548082732),
            ]
        ),
        other => panic!("expected DisjointUnion, got {other:?}"),
    }
}
