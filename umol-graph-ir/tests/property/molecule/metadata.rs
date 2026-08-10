//! Molecule metadata properties.
//!
//! These properties intentionally overlap concrete metadata unit tests: the
//! unit tests pin examples and error variants, while this module states the
//! lookup, construction, namespace, and remapping laws over generated values.

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
    fn test_molecule_metadata_lookup_roundtrip(dsl in molecule_dsl_strategy()) {
        let metadata = dsl.metadata();

        for (entity, keyword) in metadata.iter_keywords() {
            prop_assert_eq!(metadata.keyword(entity), Some(keyword));
            prop_assert_eq!(metadata.entity(keyword), Some(entity));
        }
        for (name, atom) in metadata.iter_atom_aliases() {
            prop_assert_eq!(metadata.atom_alias(name), Some(atom));
            prop_assert_eq!(metadata.atom_alias_name(atom), Some(name));
        }
    }

    #[test]
    fn test_molecule_metadata_insertion_atomicity(dsl in molecule_dsl_strategy()) {
        let metadata = dsl.metadata();
        let (alias, atom) = metadata
            .iter_atom_aliases()
            .next()
            .expect("generated metadata has aliases");

        let mut keyword_collision = metadata.clone();
        let expected = keyword_collision.clone();
        prop_assert_eq!(
            keyword_collision.set_keyword(Entity::Atom(AtomId(0)), alias),
            Err(MetadataError::DuplicateKeyword(alias.to_string())),
        );
        prop_assert_eq!(keyword_collision, expected);

        let mut alias_collision = metadata.clone();
        let expected = alias_collision.clone();
        prop_assert_eq!(
            alias_collision.add_atom_alias(format!("{alias}_duplicate"), atom.clone()),
            Err(MetadataError::DuplicateAtomAlias(alias.to_string())),
        );
        prop_assert_eq!(alias_collision, expected);
    }

    #[test]
    fn test_molecule_metadata_namespace(dsl in molecule_dsl_strategy()) {
        let metadata = dsl.metadata();
        let keywords = metadata
            .iter_keywords()
            .map(|(_, keyword)| keyword)
            .collect::<HashSet<_>>();
        let aliases = metadata
            .iter_atom_aliases()
            .map(|(name, _)| name)
            .collect::<HashSet<_>>();

        prop_assert_eq!(keywords.len(), metadata.iter_keywords().len());
        prop_assert_eq!(aliases.len(), metadata.iter_atom_aliases().len());
        prop_assert!(keywords.is_disjoint(&aliases));
    }

    #[test]
    fn test_molecule_dsl_context_metadata(dsl in molecule_dsl_strategy()) {
        let rendered = dsl.to_edn().to_string();
        let parsed = MoleculeDsl::from_edn_str(&rendered)
            .map_err(|error| TestCaseError::fail(format!("parse failed: {error}")))?;

        prop_assert_eq!(parsed.metadata(), dsl.metadata());
    }

    #[test]
    fn test_molecule_dsl_new_parsed(dsl in molecule_dsl_strategy()) {
        let rendered = dsl.to_edn().to_string();
        let parsed = MoleculeDsl::from_edn_str(&rendered)
            .map_err(|error| TestCaseError::fail(format!("parse failed: {error}")))?;
        let expected = parsed.clone();
        let (molecule, metadata) = parsed.into_parts();

        prop_assert_eq!(MoleculeDsl::new(molecule, metadata), Ok(expected));
    }

    #[test]
    fn test_molecule_dsl_new_error(
        (molecule, metadata, entity) in invalid_molecule_dsl_parts_strategy(),
    ) {
        prop_assert_eq!(
            MoleculeDsl::new(molecule, metadata),
            Err(MetadataError::EntityOutOfRange(entity)),
        );
    }

    #[test]
    fn test_molecule_metadata_remap_identity(dsl in molecule_dsl_strategy()) {
        let atoms = dsl.molecule().atoms().ids().collect::<Vec<_>>();
        let identity = dsl.molecule().induced_subgraph(&atoms).reverse();

        prop_assert_eq!(
            dsl.metadata().clone().remap(&identity),
            dsl.metadata().clone(),
        );
    }

    #[test]
    fn test_molecule_metadata_remap_composition(
        (molecule, metadata, atoms) in molecule_metadata_with_atom_subset_strategy(),
    ) {
        let sub_to_host = molecule.induced_subgraph(&atoms);
        let host_to_sub = sub_to_host.reverse();
        let direct = metadata
            .clone()
            .remap(&host_to_sub.compose(&sub_to_host));
        let sequential = metadata
            .clone()
            .remap(&host_to_sub)
            .remap(&sub_to_host);

        prop_assert_eq!(direct, sequential);
    }

    #[test]
    fn test_molecule_metadata_remap_roundtrip(dsl in molecule_dsl_strategy()) {
        let mut atoms = dsl.molecule().atoms().ids().collect::<Vec<_>>();
        atoms.reverse();
        let copy_to_host = dsl.molecule().induced_subgraph(&atoms);
        let host_to_copy = copy_to_host.reverse();
        let roundtrip = dsl
            .metadata()
            .clone()
            .remap(&host_to_copy)
            .remap(&copy_to_host);

        prop_assert_eq!(roundtrip, dsl.metadata().clone());
    }

    #[test]
    fn test_molecule_metadata_remap_partial(
        (molecule, metadata, atoms) in molecule_metadata_with_atom_subset_strategy(),
    ) {
        let host_to_sub = molecule.induced_subgraph(&atoms).reverse();
        let restricted = metadata.clone().remap(&host_to_sub);

        for (entity, keyword) in metadata.iter_keywords() {
            match host_to_sub.right_of(entity) {
                Some(right) => prop_assert_eq!(restricted.keyword(right), Some(keyword)),
                None => prop_assert_eq!(restricted.entity(keyword), None),
            }
        }
    }

    #[test]
    fn test_molecule_metadata_remap_aliases(
        (molecule, metadata, atoms) in molecule_metadata_with_atom_subset_strategy(),
    ) {
        let host_to_sub = molecule.induced_subgraph(&atoms).reverse();
        let expected = metadata
            .iter_atom_aliases()
            .map(|(name, atom)| (name.to_string(), atom.clone()))
            .collect::<Vec<_>>();
        let remapped = metadata.remap(&host_to_sub);

        prop_assert_eq!(
            remapped
                .iter_atom_aliases()
                .map(|(name, atom)| (name.to_string(), atom.clone()))
                .collect::<Vec<_>>(),
            expected,
        );
    }
}
