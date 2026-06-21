use proptest::prelude::*;
use umol_ast::ast::SubstructureMatchAlgorithm::{GraphAndOverlays, Incidence};
use umol_ast::ast::{AtomId, MoleculeAst, SubstructureMatchAlgorithm};
use umol_graph_core::SubgraphIsomorphismAlgorithm::{
    ArcMatch, RayKirsch, Ri, Ullmann, Vf2, Vf2Rdkit,
};
use umol_graph_core::{SubgraphIsomorphismAlgorithm, ARCMATCH_DEFAULT_PATH_LENGTH};

use crate::strategies::molecule_ast_strategy;

const SUBISO: [SubgraphIsomorphismAlgorithm; 6] = [
    Vf2,
    Ullmann,
    Ri,
    ArcMatch { path_length: ARCMATCH_DEFAULT_PATH_LENGTH },
    Vf2Rdkit,
    RayKirsch,
];
const STRATEGIES: [SubstructureMatchAlgorithm; 2] = [GraphAndOverlays, Incidence];

fn sorted_matches(
    pattern: &MoleculeAst,
    host: &MoleculeAst,
    strategy: SubstructureMatchAlgorithm,
    subiso: SubgraphIsomorphismAlgorithm,
) -> Vec<Vec<AtomId>> {
    let mut occurrences: Vec<Vec<AtomId>> = pattern
        .substructure_matches(host, strategy, subiso)
        .iter()
        .map(|e| e.host_atoms().to_vec())
        .collect();
    occurrences.sort();
    occurrences
}

proptest! {
    #[test]
    fn test_substructure_cross_validation(
        host in molecule_ast_strategy(),
        pattern in molecule_ast_strategy(),
    ) {
        let reference = sorted_matches(&pattern, &host, GraphAndOverlays, Vf2);
        for strategy in STRATEGIES {
            for subiso in SUBISO {
                prop_assert_eq!(
                    sorted_matches(&pattern, &host, strategy, subiso),
                    reference.clone(),
                    "{:?}/{:?}", strategy, subiso
                );
            }
        }
    }

    #[test]
    fn test_substructure_cross_validation_planted(
        (host, subset) in molecule_ast_strategy()
            .prop_filter("non-empty", |m| m.atoms().count() > 0)
            .prop_flat_map(|m| {
                let n = m.atoms().count();
                let ids: Vec<AtomId> = (0..n as u32).map(AtomId).collect();
                (Just(m), prop::sample::subsequence(ids, 1..=n))
            }),
    ) {
        // A pattern induced from `host` shares its structure, so it exercises the
        // match path far more than independent random pairs. (No `!is_empty` check:
        // the generator may emit constraints inconsistent with the topology — e.g. a
        // stored valence differing from the derived one — so a self-match is not
        // guaranteed; only cross-strategy / cross-algorithm agreement is.)
        let pattern = host.induced_subgraph(&subset).extract();
        let reference = sorted_matches(&pattern, &host, GraphAndOverlays, Vf2);
        for strategy in STRATEGIES {
            for subiso in SUBISO {
                prop_assert_eq!(
                    sorted_matches(&pattern, &host, strategy, subiso),
                    reference.clone(),
                    "{:?}/{:?}", strategy, subiso
                );
            }
        }
    }
}
