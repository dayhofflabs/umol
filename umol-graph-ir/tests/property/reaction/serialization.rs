//! Property tests for reaction serialization.

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

    /// The reaction round-trips through the EDN surface: render → parse reaches a
    /// fixpoint, exercising the atom/bond add / remove / modify-field delta ops
    /// (`Reaction::to_edn` then `from_edn`, twice, must agree).
    #[test]
    fn test_reaction_edn_roundtrip_stable(reaction in reaction_strategy()) {
        let once = Reaction::from_edn(&reaction.to_edn())
            .map_err(|e| TestCaseError::fail(format!("first reparse failed: {e}")))?;
        let twice = Reaction::from_edn(&once.to_edn())
            .map_err(|e| TestCaseError::fail(format!("second reparse failed: {e}")))?;
        prop_assert_eq!(once, twice);
    }

    #[test]
    fn test_reaction_dsl_to_edn_from_edn_roundtrip(dsl in reaction_dsl_strategy()) {
        let via_tree = ReactionDsl::from_edn(&dsl.to_edn())
            .map_err(|error| TestCaseError::fail(format!("tree parse failed: {error}")))?;
        let rendered = dsl.to_edn().to_string();
        let via_stream = ReactionDsl::from_edn_str(&rendered)
            .map_err(|error| TestCaseError::fail(format!("streaming parse failed: {error}")))?;
        prop_assert_eq!(&via_tree, &dsl);
        prop_assert_eq!(via_stream, dsl);
    }

    #[test]
    fn test_reaction_dsl_parser_parity(input in any::<String>()) {
        let via_stream = ReactionDsl::from_edn_str(&input).ok();
        let via_tree = read_string(&input)
            .ok()
            .and_then(|edn| ReactionDsl::from_edn(&edn).ok());
        prop_assert_eq!(via_stream, via_tree);
    }

    #[test]
    fn test_reaction_defaults_roundtrip(reaction in comprehensive_reaction_strategy()) {
        let defaults = ReactionDefaults::new();
        let rebuilt = ReactionDsl::from_ir(&reaction, &defaults).into_ir(&defaults);
        prop_assert_eq!(rebuilt, reaction);
    }

    #[test]
    fn test_reaction_defaults_roundtrip_ground(reaction in comprehensive_reaction_strategy()) {
        let required = ReactionDefaults::new();
        let ground = ReactionDefaults::concrete();
        let grounded = ReactionDsl::from_ir(&reaction, &required).into_ir(&ground);
        let rebuilt = ReactionDsl::from_ir(&grounded, &ground).into_ir(&ground);
        prop_assert_eq!(rebuilt, grounded);
    }

    #[test]
    fn test_reaction_span_dsl_to_edn_from_edn_roundtrip(dsl in reaction_span_dsl_strategy()) {
        let via_tree = ReactionSpanDsl::from_edn(&dsl.to_edn())
            .map_err(|error| TestCaseError::fail(format!("tree parse failed: {error}")))?;
        let rendered = dsl.to_edn().to_string();
        let via_stream = ReactionSpanDsl::from_edn_str(&rendered)
            .map_err(|error| TestCaseError::fail(format!("streaming parse failed: {error}")))?;
        prop_assert_eq!(&via_tree, &dsl);
        prop_assert_eq!(via_stream, dsl);
    }

    #[test]
    fn test_reaction_span_dsl_parser_parity(input in any::<String>()) {
        let via_stream = ReactionSpanDsl::from_edn_str(&input).ok();
        let via_tree = read_string(&input)
            .ok()
            .and_then(|edn| ReactionSpanDsl::from_edn(&edn).ok());
        prop_assert_eq!(via_stream, via_tree);
    }

    #[test]
    fn test_reaction_span_defaults_roundtrip_ground(
        reaction in materializable_reaction_strategy(),
    ) {
        let span = reaction
            .to_reaction_span()
            .expect("generated reaction materializes a span");
        let required = MoleculeDefaults::new();
        let ground = MoleculeDefaults::concrete();
        let grounded = ReactionSpanDsl::from_ir(&span, &required).into_ir(&ground);
        let rebuilt = ReactionSpanDsl::from_ir(&grounded, &ground).into_ir(&ground);
        prop_assert_eq!(rebuilt, grounded);
    }
}
