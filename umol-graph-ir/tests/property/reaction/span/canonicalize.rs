//! Aggregate canonicalization and reaction-normal-form properties for reaction spans.
//!
//! Integrity-valid generated spans are compared with independently renumbered union frames. The
//! properties cover the complete normalization/reframe/canonicalize fixpoint and absorption matrix,
//! complete equality and canonical-hash relations, LHS anchoring, reaction-normal-form convergence,
//! integrity, and the documented weakened reversal law.

use std::hash::{DefaultHasher, Hash, Hasher};

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{
    Canonicalize, CanonicalizeContext, Contradiction, Normalize, ReactionSpan,
    ReactionSpanCanonicalizeError, Reframe,
};

use super::reaction_span_scenario_strategy;
use crate::strategies::{
    intrinsic_contradiction_scenario_strategy, standardization_scenario_strategy,
};

fn context() -> CanonicalizeContext {
    CanonicalizeContext {
        para_stereo: false,
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    }
}

fn structural_hash<T: Hash>(value: &T) -> u64 {
    let mut hasher = DefaultHasher::new();
    value.hash(&mut hasher);
    hasher.finish()
}

fn reaction_normal_form(span: &ReactionSpan) -> ReactionSpan {
    span.to_reaction()
        .to_reaction_span()
        .expect("a reaction derived from an integrity-valid span materializes")
}

fn reversed(span: &ReactionSpan) -> Result<ReactionSpan, TestCaseError> {
    span.to_reaction()
        .reverse()
        .map_err(|error| TestCaseError::fail(format!("span reversal failed: {error}")))?
        .to_reaction_span()
        .map_err(|error| {
            TestCaseError::fail(format!("reversed reaction did not materialize: {error}"))
        })
}

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_reaction_span_canonicalize(
        scenario in reaction_span_scenario_strategy(),
    ) {
        let context = context();
        let renumbered = scenario.span.remap(&scenario.first);
        let canonical = scenario.span.clone().canonicalize(&context);
        let renumbered_canonical = renumbered.canonicalize(&context);

        prop_assert_eq!(&renumbered_canonical, &canonical);
        if let Ok(canonical) = canonical {
            let (with_correspondence, correspondence) = scenario
                .span
                .clone()
                .canonicalize_with_correspondence(&context)
                .expect("successful canonicalization returns its correspondence");

            prop_assert!(canonical.to_reaction().to_reaction_span().is_ok());
            prop_assert_eq!(&with_correspondence, &canonical);
            prop_assert_eq!(
                scenario
                    .span
                    .remap(&correspondence)
                    .reframe(),
                Ok(canonical.clone()),
            );
            prop_assert_eq!(reaction_normal_form(&canonical), canonical.clone());
        }
    }

    #[test]
    fn test_reaction_span_canonicalize_standardization(
        scenario in standardization_scenario_strategy(),
    ) {
        let context = context();
        let source = scenario.span;
        let identical = source.clone();
        let normalized = source.clone().normalize().map_err(|_| {
            TestCaseError::fail("generated reaction span is intrinsically contradictory")
        })?;
        let normalized_again = normalized.clone().normalize().map_err(|_| {
            TestCaseError::fail("normalized reaction span became contradictory")
        })?;
        let reframed = source.clone().reframe().map_err(|_| {
            TestCaseError::fail("generated reaction span is intrinsically contradictory")
        })?;
        let canonical = source.clone().canonicalize(&context).map_err(|error| {
            TestCaseError::fail(format!("generated reaction span did not canonicalize: {error}"))
        })?;
        let renumbered = source.remap(&scenario.correspondence);

        prop_assert_eq!(normalized.clone().normalize(), Ok(normalized.clone()));
        prop_assert_eq!(reframed.clone().reframe(), Ok(reframed.clone()));
        prop_assert_eq!(canonical.clone().canonicalize(&context), Ok(canonical.clone()));
        prop_assert_eq!(normalized.clone().reframe(), Ok(reframed.clone()));
        prop_assert_eq!(reframed.clone().normalize(), Ok(reframed.clone()));
        prop_assert_eq!(normalized.clone().canonicalize(&context), Ok(canonical.clone()));
        prop_assert_eq!(reframed.clone().canonicalize(&context), Ok(canonical.clone()));
        prop_assert_eq!(canonical.clone().normalize(), Ok(canonical.clone()));
        prop_assert_eq!(canonical.clone().reframe(), Ok(canonical.clone()));

        prop_assert_eq!(&source, &identical);
        prop_assert!(source.normalized_eq(&identical));
        prop_assert!(source.framed_eq(&identical));
        prop_assert!(source.canonical_eq(&identical, &context));
        prop_assert!(source.normalized_eq(&normalized));
        prop_assert_eq!(
            source.normalized_eq(&normalized),
            normalized.normalized_eq(&source),
        );
        prop_assert!(normalized.normalized_eq(&normalized_again));
        prop_assert!(source.normalized_eq(&normalized_again));
        prop_assert!(source.framed_eq(&normalized));
        prop_assert!(normalized.framed_eq(&reframed));
        prop_assert!(source.framed_eq(&reframed));
        prop_assert_eq!(source.framed_eq(&reframed), reframed.framed_eq(&source));
        prop_assert!(source.canonical_eq(&reframed, &context));
        prop_assert_eq!(
            source.canonical_eq(&reframed, &context),
            reframed.canonical_eq(&source, &context),
        );
        prop_assert!(reframed.canonical_eq(&renumbered, &context));
        prop_assert!(source.canonical_eq(&renumbered, &context));
        prop_assert!(renumbered.canonical_eq(&canonical, &context));
    }

    #[test]
    fn test_reaction_span_canonical_hash(
        scenario in reaction_span_scenario_strategy(),
    ) {
        let context = context();
        let renumbered = scenario.span.remap(&scenario.first);

        prop_assert_eq!(
            scenario.span.clone().canonical_hash(&context),
            renumbered.clone().canonical_hash(&context),
        );
        if let Ok(canonical) = scenario.span.clone().canonicalize(&context) {
            prop_assert_eq!(
                scenario.span.canonical_hash(&context),
                Ok(structural_hash(&canonical)),
            );
        }
    }

    #[test]
    fn test_reaction_span_canonical_eq(
        scenario in reaction_span_scenario_strategy(),
    ) {
        let context = context();
        let renumbered = scenario.span.remap(&scenario.first);
        let canonical = scenario.span.clone().canonicalize(&context);

        prop_assert!(scenario.span.canonical_eq(&scenario.span, &context));
        prop_assert!(scenario.span.canonical_eq(&renumbered, &context));
        prop_assert_eq!(
            scenario.span.canonical_eq(&renumbered, &context),
            renumbered.canonical_eq(&scenario.span, &context),
        );
        if let Ok(canonical) = canonical {
            prop_assert!(renumbered.canonical_eq(&canonical, &context));
            prop_assert!(scenario.span.canonical_eq(&canonical, &context));
        }
    }

    #[test]
    fn test_reaction_span_canonicalize_roundtrip(
        scenario in reaction_span_scenario_strategy(),
    ) {
        let context = context();
        let reordered = scenario.span.remap(&scenario.first);
        let normalized = reaction_normal_form(&reordered);
        let normalized_twice = reaction_normal_form(&normalized);

        prop_assert_eq!(normalized_twice, normalized.clone());
        prop_assert_eq!(
            normalized.canonicalize(&context),
            reordered.canonicalize(&context),
        );
    }

    #[test]
    fn test_reaction_span_canonicalize_reversal(
        scenario in reaction_span_scenario_strategy(),
    ) {
        let context = context();
        let canonical = scenario.span.clone().canonicalize(&context).map_err(|error| {
            TestCaseError::fail(format!("generated span did not canonicalize: {error}"))
        })?;
        let canonical_reversed = reversed(&canonical)?.canonicalize(&context);
        let reversed_canonical = reversed(&scenario.span)?.canonicalize(&context);

        prop_assert_eq!(canonical_reversed, reversed_canonical);
    }

    #[test]
    fn test_reaction_span_canonical_eq_contradiction(
        scenario in intrinsic_contradiction_scenario_strategy(),
    ) {
        let context = context();
        let [first, second, third] = scenario.spans;

        prop_assert_ne!(&first, &second);
        prop_assert_ne!(&second, &third);
        prop_assert!(first.normalized_eq(&second));
        prop_assert!(second.normalized_eq(&third));
        prop_assert!(first.normalized_eq(&third));
        prop_assert!(first.framed_eq(&second));
        prop_assert!(second.framed_eq(&third));
        prop_assert!(first.framed_eq(&third));
        prop_assert!(first.canonical_eq(&second, &context));
        prop_assert!(second.canonical_eq(&third, &context));
        prop_assert!(first.canonical_eq(&third, &context));
        prop_assert_eq!(
            first.canonicalize_with_correspondence(&context),
            Err(ReactionSpanCanonicalizeError::Contradiction(Contradiction)),
        );
    }
}
