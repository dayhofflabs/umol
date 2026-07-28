//! Reaction and reaction-span metadata properties.
//!
//! These properties intentionally overlap concrete metadata and checked-
//! constructor unit tests: they state the same contracts over generated
//! reactions, scopes, aliases, and span frames.

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

    #[test]
    fn test_reaction_metadata_lookup_roundtrip(dsl in reaction_dsl_strategy()) {
        let metadata = dsl.metadata();

        for (entity, keyword) in metadata.iter_keywords() {
            prop_assert_eq!(metadata.keyword(entity), Some(keyword));
            prop_assert_eq!(metadata.entity(keyword), Some(entity));
        }
        for (entity, keyword) in metadata.iter_delta_keywords() {
            prop_assert_eq!(metadata.delta_keyword(entity), Some(keyword));
            prop_assert_eq!(metadata.delta_entity(keyword), Some(entity));
        }
        for (name, atom) in metadata.lhs().iter_atom_aliases() {
            prop_assert_eq!(metadata.atom_alias(name), Some(atom));
            prop_assert_eq!(metadata.atom_alias_name(atom), Some(name));
        }
        for (name, atom) in metadata.iter_reaction_atom_aliases() {
            prop_assert_eq!(metadata.atom_alias(name), Some(atom));
            prop_assert_eq!(metadata.atom_alias_name(atom), Some(name));
        }
    }

    #[test]
    fn test_reaction_metadata_insertion_atomicity(dsl in reaction_dsl_strategy()) {
        let metadata = dsl.metadata();

        let mut keyword_collision = metadata.clone();
        let expected = keyword_collision.clone();
        prop_assert_eq!(
            keyword_collision
                .set_delta_keyword(Entity::Atom(AtomId(u32::MAX)), "reaction_alias"),
            Err(MetadataError::DuplicateKeyword("reaction_alias".to_string())),
        );
        prop_assert_eq!(keyword_collision, expected);

        let (alias, atom) = metadata
            .lhs()
            .iter_atom_aliases()
            .next()
            .expect("generated lhs metadata has aliases");
        let mut alias_collision = metadata.clone();
        let expected = alias_collision.clone();
        prop_assert_eq!(
            alias_collision.add_atom_alias("duplicate_alias", atom.clone()),
            Err(MetadataError::DuplicateAtomAlias(alias.to_string())),
        );
        prop_assert_eq!(alias_collision, expected);
    }

    #[test]
    fn test_reaction_metadata_namespace(dsl in reaction_dsl_strategy()) {
        let metadata = dsl.metadata();
        let keywords = metadata
            .iter_keywords()
            .map(|(_, keyword)| keyword)
            .collect::<HashSet<_>>();
        let lhs_aliases = metadata
            .lhs()
            .iter_atom_aliases()
            .map(|(name, _)| name)
            .collect::<HashSet<_>>();
        let reaction_aliases = metadata
            .iter_reaction_atom_aliases()
            .map(|(name, _)| name)
            .collect::<HashSet<_>>();

        prop_assert_eq!(keywords.len(), metadata.iter_keywords().len());
        prop_assert_eq!(lhs_aliases.len(), metadata.lhs().iter_atom_aliases().len());
        prop_assert_eq!(
            reaction_aliases.len(),
            metadata.iter_reaction_atom_aliases().len(),
        );
        prop_assert!(keywords.is_disjoint(&lhs_aliases));
        prop_assert!(keywords.is_disjoint(&reaction_aliases));
        prop_assert!(lhs_aliases.is_disjoint(&reaction_aliases));
    }

    #[test]
    fn test_reaction_dsl_context_metadata(dsl in reaction_dsl_strategy()) {
        let rendered = dsl.to_edn().to_string();
        let parsed = ReactionDsl::from_edn_str(&rendered)
            .map_err(|error| TestCaseError::fail(format!("parse failed: {error}")))?;

        prop_assert_eq!(parsed.metadata(), dsl.metadata());
    }

    #[test]
    fn test_reaction_dsl_new_parsed(dsl in reaction_dsl_strategy()) {
        let rendered = dsl.to_edn().to_string();
        let parsed = ReactionDsl::from_edn_str(&rendered)
            .map_err(|error| TestCaseError::fail(format!("parse failed: {error}")))?;
        let expected = parsed.clone();
        let (ast, metadata) = parsed.into_parts();

        prop_assert_eq!(ReactionDsl::new(ast, metadata), Ok(expected));
    }

    #[test]
    fn test_reaction_dsl_new_error(
        (ast, metadata, expected) in invalid_reaction_dsl_parts_strategy(),
    ) {
        prop_assert_eq!(ReactionDsl::new(ast, metadata), Err(expected));
    }

    #[test]
    fn test_reaction_span_dsl_context_metadata(dsl in reaction_span_dsl_strategy()) {
        let rendered = dsl.to_edn().to_string();
        let parsed = ReactionSpanDsl::from_edn_str(&rendered)
            .map_err(|error| TestCaseError::fail(format!("parse failed: {error}")))?;

        prop_assert_eq!(parsed.metadata(), dsl.metadata());
    }

    #[test]
    fn test_reaction_span_dsl_new_parsed(dsl in reaction_span_dsl_strategy()) {
        let rendered = dsl.to_edn().to_string();
        let parsed = ReactionSpanDsl::from_edn_str(&rendered)
            .map_err(|error| TestCaseError::fail(format!("parse failed: {error}")))?;
        let expected = parsed.clone();
        let (ast, metadata) = parsed.into_parts();

        prop_assert_eq!(ReactionSpanDsl::new(ast, metadata), Ok(expected));
    }

    #[test]
    fn test_reaction_span_dsl_new_error(
        (ast, metadata, entity) in invalid_reaction_span_dsl_parts_strategy(),
    ) {
        prop_assert_eq!(
            ReactionSpanDsl::new(ast, metadata),
            Err(MetadataError::EntityOutOfRange(entity)),
        );
    }
}
