//! Reaction construction, canonicalization, and derivation properties.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_ast::ast::{SubstructureMatchAlgorithm, SubstructureMatchConfig};
use umol_graph_core::{
    Correspondence, NodeId, RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm,
};

use crate::strategies::*;

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_reaction_ast_canonicalize(reaction in reaction_strategy()) {
        if let Ok(canonical) = reaction.canonicalize() {
            prop_assert_eq!(canonical.clone().canonicalize(), Ok(canonical));
        }
    }

    #[test]
    fn test_reaction_ast_from_sides(reaction in comprehensive_reaction_strategy()) {
        if let Ok(span) = reaction.to_reaction_span() {
            let rebuilt = ReactionAst::from_sides(
                span.lhs(),
                span.rhs(),
                span.correspondence().atoms().clone(),
            );
            let rebuilt_span = rebuilt.to_reaction_span().map_err(|error| {
                TestCaseError::fail(format!("rebuilt reaction did not materialize: {error}"))
            })?;
            prop_assert_eq!(rebuilt_span, span);
        }
    }

    #[test]
    fn test_reaction_ast_apply_at(reaction in reaction_strategy()) {
        let atom_count = reaction.lhs.atoms().count();
        let atom_images = (0..atom_count).map(NodeId::from).collect::<Vec<_>>();
        let correspondence = MoleculeCorrespondence::induce(
            &reaction.lhs,
            &reaction.lhs,
            Correspondence::from_images(&atom_images, atom_count),
        );
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
        let atom_images = (0..atom_count).map(NodeId::from).collect::<Vec<_>>();
        let correspondence = MoleculeCorrespondence::induce(
            &reaction.lhs,
            &reaction.lhs,
            Correspondence::from_images(&atom_images, atom_count),
        );
        let derivation = reaction.apply_at(&reaction.lhs, &correspondence).map_err(|error| {
            TestCaseError::fail(format!("identity application failed: {error}"))
        })?;

        prop_assert_eq!(derivation.reverse().reverse(), derivation.clone());

        let recovered = derivation.to_reaction();
        let recovered_correspondence = MoleculeCorrespondence::induce(
            derivation.lhs(),
            derivation.lhs(),
            correspondence.atoms().clone(),
        );
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
