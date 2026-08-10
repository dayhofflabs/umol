//! Reaction construction, normalization, and derivation properties.

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
    fn test_reaction_normalize(reaction in reaction_strategy()) {
        if let Ok(normalized) = reaction.normalize() {
            prop_assert_eq!(normalized.clone().normalize(), Ok(normalized));
        }
    }

    #[test]
    fn test_reaction_apply_at(reaction in reaction_strategy()) {
        let atom_count = reaction.lhs.atoms().count();
        let atom_images = (0..atom_count).map(AtomId::from).collect::<Vec<_>>();
        let correspondence = MoleculeCorrespondence::induce(
            &reaction.lhs,
            &reaction.lhs,
            Correspondence::from_images(&atom_images, atom_count),
        )
        .expect("identity atom correspondence induces a molecule correspondence");
        let direct = reaction.apply_at(&reaction.lhs, &correspondence).map_err(|error| {
            TestCaseError::fail(format!("identity application failed: {error}"))
        })?;
        let mut applications = reaction
            .apply(
                &reaction.lhs,
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
            if application.rhs() == direct.rhs() {
                found = true;
                break;
            }
        }
        prop_assert!(found);
    }

    #[test]
    fn test_reaction_derivation_roundtrip(reaction in reaction_strategy()) {
        let atom_count = reaction.lhs.atoms().count();
        let atom_images = (0..atom_count).map(AtomId::from).collect::<Vec<_>>();
        let correspondence = MoleculeCorrespondence::induce(
            &reaction.lhs,
            &reaction.lhs,
            Correspondence::from_images(&atom_images, atom_count),
        )
        .expect("identity atom correspondence induces a molecule correspondence");
        let derivation = reaction.apply_at(&reaction.lhs, &correspondence).map_err(|error| {
            TestCaseError::fail(format!("identity application failed: {error}"))
        })?;

        prop_assert_eq!(derivation.reverse().reverse(), derivation.clone());

        let recovered = derivation.to_reaction();
        let recovered_correspondence = MoleculeCorrespondence::induce(
            derivation.lhs(),
            derivation.lhs(),
            correspondence.atoms().clone(),
        )
        .expect("identity atom correspondence induces a molecule correspondence");
        let recovered_derivation = recovered
            .apply_at(derivation.lhs(), &recovered_correspondence)
            .map_err(|error| {
                TestCaseError::fail(format!("recovered reaction did not apply: {error}"))
            })?;
        prop_assert_eq!(recovered_derivation.rhs(), derivation.rhs());

        let identity = derivation.chain(&derivation.reverse());
        prop_assert_eq!(identity.lhs(), derivation.lhs());
        prop_assert_eq!(identity.rhs(), derivation.lhs());
        prop_assert_eq!(identity.comap(), &derivation.comap().compose(&derivation.comap().reverse()));
    }
}
