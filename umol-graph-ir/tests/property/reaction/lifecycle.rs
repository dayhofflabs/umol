//! Reaction construction and application properties.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_core::{
    Correspondence, RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm,
};
use umol_graph_ir::ir::{SubstructureMatchAlgorithm, SubstructureMatchConfig};

use crate::strategies::*;

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_reaction_try_new_roundtrip(reaction in comprehensive_reaction_strategy()) {
        let republished = Reaction::try_new(reaction.lhs().clone(), reaction.deltas().clone());

        prop_assert_eq!(republished, Ok(reaction));
    }

    #[test]
    fn test_reaction_apply_at(reaction in reaction_strategy()) {
        let atom_count = reaction.lhs().atoms().count();
        let atom_images = (0..atom_count).map(AtomId::from).collect::<Vec<_>>();
        let correspondence = MoleculeCorrespondence::induce(
            reaction.lhs(),
            reaction.lhs(),
            Correspondence::from_images(&atom_images, atom_count),
        )
        .expect("identity atom correspondence induces a molecule correspondence");
        let direct = reaction.apply_at(reaction.lhs(), &correspondence).map_err(|error| {
            TestCaseError::fail(format!("identity application failed: {error}"))
        })?.expect("identity application is applicable");
        let mut applications = reaction
            .apply(
                reaction.lhs(),
                SubstructureMatchConfig {
                    match_algorithm: SubstructureMatchAlgorithm::GraphAndOverlays,
                    subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::Vf2,
                    relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                },
            )
            .map_err(|error| TestCaseError::fail(format!("matching failed: {error}")))?;
        let mut found = false;
        for application in &mut applications {
            let application = application.map_err(|error| {
                TestCaseError::fail(format!("matched application failed: {error}"))
            })?;
            if application == direct {
                found = true;
                break;
            }
        }
        prop_assert!(found);
    }

    #[test]
    fn test_reaction_tracked_apply_roundtrip(reaction in reaction_strategy()) {
        let atom_count = reaction.lhs().atoms().count();
        let atom_images = (0..atom_count).map(AtomId::from).collect::<Vec<_>>();
        let correspondence = MoleculeCorrespondence::induce(
            reaction.lhs(),
            reaction.lhs(),
            Correspondence::from_images(&atom_images, atom_count),
        )
        .expect("identity atom correspondence induces a molecule correspondence");
        let (product, witness) = reaction.tracked_apply(reaction.lhs(), SubstructureMatchConfig {
            match_algorithm: SubstructureMatchAlgorithm::GraphAndOverlays,
            subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::Vf2,
            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
        }).unwrap().next().unwrap().unwrap();

        prop_assert_eq!(witness.reverse().reverse(), witness.clone());
        let recovered = Reaction::new(
            reaction.lhs().clone(),
            reaction.lhs().difference_to(&product, &witness).unwrap(),
        );
        let recovered_product = recovered.apply_at(reaction.lhs(), &correspondence)
            .unwrap().expect("recovered reaction is applicable");
        prop_assert_eq!(&recovered_product, &product);

        let roundtrip = witness.compose(&witness.reverse()).unwrap();
        prop_assert!(roundtrip.is_compatible(reaction.lhs(), reaction.lhs()));
        for &(left, right) in roundtrip.atoms().matched_pairs() {
            prop_assert_eq!(left, right);
        }
    }
}
