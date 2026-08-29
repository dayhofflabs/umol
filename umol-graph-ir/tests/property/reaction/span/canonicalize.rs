//! Aggregate canonicalization and reaction-normal-form properties for reaction spans.
//!
//! Integrity-valid generated spans are compared with independently renumbered union frames. The
//! properties cover exact complete-form idempotence, all level-specific equality and canonical-hash
//! relations, LHS anchoring, reaction-normal-form convergence, integrity, and the documented
//! weakened reversal law.

use std::hash::{DefaultHasher, Hash, Hasher};

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{Canonicalize, CanonicalizeContext, DescriptionLevel, ReactionSpan};

use super::reaction_span_scenario_strategy;

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
        prop_assert_eq!(
            scenario
                .span
                .clone()
                .canonicalize_by(DescriptionLevel::Full, &context),
            canonical.clone(),
        );
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
                    .canonicalize(&context),
                Ok(canonical.clone()),
            );
            prop_assert_eq!(reaction_normal_form(&canonical), canonical.clone());
            prop_assert_eq!(canonical.clone().canonicalize(&context), Ok(canonical));
        }
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
        for level in [
            DescriptionLevel::Topology,
            DescriptionLevel::Constitution,
            DescriptionLevel::Structure,
            DescriptionLevel::Full,
        ] {
            prop_assert_eq!(
                scenario.span.clone().canonical_hash_by(level, &context),
                renumbered.clone().canonical_hash_by(level, &context),
            );
        }
        prop_assert_eq!(
            scenario
                .span
                .clone()
                .canonical_hash_by(DescriptionLevel::Full, &context),
            scenario.span.clone().canonical_hash(&context),
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
    fn test_reaction_span_canonical_eq_by(
        scenario in reaction_span_scenario_strategy(),
    ) {
        let context = context();
        let renumbered = scenario.span.remap(&scenario.first);

        for level in [
            DescriptionLevel::Topology,
            DescriptionLevel::Constitution,
            DescriptionLevel::Structure,
            DescriptionLevel::Full,
        ] {
            let canonical = scenario.span.clone().canonicalize_by(level, &context);

            prop_assert!(scenario.span.canonical_eq_by(&scenario.span, level, &context));
            prop_assert!(scenario.span.canonical_eq_by(&renumbered, level, &context));
            prop_assert_eq!(
                scenario.span.canonical_eq_by(&renumbered, level, &context),
                renumbered.canonical_eq_by(&scenario.span, level, &context),
            );
            if let Ok(canonical) = canonical {
                prop_assert!(renumbered.canonical_eq_by(&canonical, level, &context));
                prop_assert!(scenario.span.canonical_eq_by(&canonical, level, &context));
            }
        }
        prop_assert_eq!(
            scenario.span.canonical_eq_by(
                &renumbered,
                DescriptionLevel::Full,
                &context,
            ),
            scenario.span.canonical_eq(&renumbered, &context),
        );
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
}
