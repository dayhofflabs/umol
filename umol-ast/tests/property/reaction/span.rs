//! Property tests for reaction spans and reaction reversal.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};

use crate::strategies::*;

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    /// Cross-validate the two span constructions: the direct `superimpose` (Strategy A) reproduces
    /// the span the delta path (`to_reaction_span`) builds. Recover `(L, R, C)` from the delta-path
    /// span and reassemble; a mismatch flags a diff-completeness or frame gap between the paths.
    #[test]
    fn test_reaction_span_ast_superimpose_matches_delta_path(reaction in reaction_strategy()) {
        if let Ok(span) = reaction.to_reaction_span() {
            let rebuilt =
                ReactionSpanAst::superimpose(&span.lhs(), &span.rhs(), &span.correspondence());
            prop_assert_eq!(rebuilt, span);
        }
    }

    /// `reverse` swaps the span's sides and reverses its correspondence. Constructing that span
    /// directly must reproduce the span obtained by reversing the reaction, including the union
    /// frame chosen for entities unmatched on only one side.
    #[test]
    fn test_reaction_ast_reverse_swaps_sides(reaction in reaction_strategy()) {
        if let (Ok(span), Ok(reverse)) = (reaction.to_reaction_span(), reaction.reverse()) {
            if let Ok(reverse_span) = reverse.to_reaction_span() {
                let expected = ReactionSpanAst::superimpose(
                    &span.rhs(),
                    &span.lhs(),
                    &span.correspondence().reverse(),
                );
                prop_assert_eq!(reverse_span, expected);
            }
        }
    }

    /// Cross-validate the two span constructions with overlays present: the direct `superimpose`
    /// reassembles the delta-path span across all overlay families, not just atoms/bonds.
    #[test]
    fn test_reaction_span_ast_superimpose_matches_delta_path_overlay(
        reaction in overlay_reaction_strategy(),
    ) {
        if let Ok(span) = reaction.to_reaction_span() {
            let rebuilt =
                ReactionSpanAst::superimpose(&span.lhs(), &span.rhs(), &span.correspondence());
            prop_assert_eq!(rebuilt, span);
        }
    }

    /// Reaction ↔ span roundtrip fidelity: recovering the reaction from a span and re-materializing
    /// reproduces the span (`to_reaction` then `to_reaction_span` is the identity on spans),
    /// exercising the overlay `EntitySpan` columns and the span→delta recovery in both directions.
    #[test]
    fn test_reaction_ast_span_roundtrip(reaction in overlay_reaction_strategy()) {
        if let Ok(span) = reaction.to_reaction_span() {
            if let Ok(rebuilt) = span.to_reaction().to_reaction_span() {
                prop_assert_eq!(rebuilt, span);
            }
        }
    }
}
