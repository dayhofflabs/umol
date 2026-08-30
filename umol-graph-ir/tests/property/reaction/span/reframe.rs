//! Reaction-span frame-transport and reframe properties.
//!
//! Integrity-valid direct spans and materialized reactions exercise the complete action algebra,
//! fused/witness agreement, normalization-prefix law, projection agreement, and the
//! operational-reaction roundtrip. These overlap deliberately: transport isolates an independently
//! supplied action, while reframe validates representative selection across all six span
//! aggregates. Pipeline fixpoint and absorption laws are checked in the canonicalization module.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_ir::ir::{Contradiction, FrameTransport, Normalize, Reframe};

use super::reaction_span_strategy;
use crate::strategies::intrinsic_contradiction_scenario_strategy;

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_reaction_span_reframe_by(span in reaction_span_strategy()) {
        let action = span.representative_action();
        let identity = action.identity();
        let inverse = action.inverse();
        let composite = action
            .compose(&inverse)
            .expect("an action and its inverse have the same domain");

        prop_assert_eq!(span.clone().reframe_by(&identity), Some(span.clone()));
        prop_assert_eq!(
            span.clone()
                .reframe_by(&action)
                .and_then(|transported| transported.reframe_by(&inverse)),
            span.clone().reframe_by(&composite),
        );
        prop_assert_eq!(span.clone().reframe_by(&composite), Some(span));
    }

    #[test]
    fn test_reaction_span_representative_action_contradiction(
        scenario in intrinsic_contradiction_scenario_strategy(),
    ) {
        for span in scenario.spans {
            let action = span.representative_action();

            prop_assert_eq!(
                action.compose(&action.identity()),
                Some(action.clone()),
            );
            prop_assert_eq!(span.clone().normalize(), Err(Contradiction));
            prop_assert_eq!(span.reframe(), Err(Contradiction));
        }
    }

    #[test]
    fn test_reaction_span_reframe_with_action(span in reaction_span_strategy()) {
        let fused = span.clone().reframe().map_err(|_| {
            TestCaseError::fail("generated reaction span is intrinsically contradictory")
        })?;
        let (witnessed, action) = span.clone().reframe_with_action().map_err(|_| {
            TestCaseError::fail("generated reaction span is intrinsically contradictory")
        })?;
        let transported = span
            .normalize()
            .map_err(|_| {
                TestCaseError::fail("generated reaction span is intrinsically contradictory")
            })?
            .reframe_by(&action)
            .ok_or_else(|| TestCaseError::fail("representative action did not cover its source"))?
            .normalize()
            .map_err(|_| {
                TestCaseError::fail("transported reaction span is intrinsically contradictory")
            })?;
        let selected_action = witnessed.representative_action();

        prop_assert_eq!(fused, witnessed.clone());
        prop_assert_eq!(transported, witnessed);
        prop_assert_eq!(selected_action.clone(), selected_action.identity());
    }

    #[test]
    fn test_reaction_span_reframe(span in reaction_span_strategy()) {
        let normalized = span.clone().normalize().map_err(|_| {
            TestCaseError::fail("generated reaction span is intrinsically contradictory")
        })?;
        let once = span.clone().reframe().map_err(|_| {
            TestCaseError::fail("generated reaction span is intrinsically contradictory")
        })?;

        prop_assert!(span.normalized_eq(&normalized));
        prop_assert!(span.framed_eq(&normalized));
        prop_assert!(span.framed_eq(&once));
    }

    #[test]
    fn test_reaction_span_reframe_projection(span in reaction_span_strategy()) {
        let lhs = span.lhs().reframe().map_err(|_| {
            TestCaseError::fail("generated lhs is intrinsically contradictory")
        })?;
        let rhs = span.rhs().reframe().map_err(|_| {
            TestCaseError::fail("generated rhs is intrinsically contradictory")
        })?;
        let reframed = span.reframe().map_err(|_| {
            TestCaseError::fail("generated reaction span is intrinsically contradictory")
        })?;

        prop_assert_eq!(reframed.lhs(), lhs);
        prop_assert_eq!(reframed.rhs(), rhs);
    }

    #[test]
    fn test_reaction_span_reframe_roundtrip(span in reaction_span_strategy()) {
        let reframed = span.reframe().map_err(|_| {
            TestCaseError::fail("generated reaction span is intrinsically contradictory")
        })?;
        let roundtrip = reframed
            .to_reaction()
            .to_reaction_span()
            .map_err(|_| TestCaseError::fail("reframed span did not materialize again"))?;

        prop_assert_eq!(roundtrip, reframed);
    }
}
