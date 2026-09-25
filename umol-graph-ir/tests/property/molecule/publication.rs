//! Representation-integrity preservation by trusted molecule publishers.
//!
//! The generated inputs contain every entity kind and frame-sensitive payloads. Republishing an
//! output through the checked editor boundary must reproduce it exactly; a trusted transform may
//! not emit a value outside the same aggregate contract.
//! Constraint-view writes are also compared with checked editor publication, using independent
//! raw constraint trees. Rejection must leave the entire source molecule unchanged.

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
    fn test_molecule_remap_integrity_preservation(
        (source, correspondence) in molecule_dense_renumbering_strategy(),
    ) {
        let published = source.remap(&correspondence);

        prop_assert_eq!(published.edit().try_build(), Ok(published));
    }

    #[test]
    fn test_molecule_extract_integrity_preservation(
        (source, atoms) in molecule_with_atom_subset_strategy(),
    ) {
        let correspondence = source.induced_subgraph(&atoms);
        let published = source.extract(&correspondence);

        prop_assert_eq!(published.edit().try_build(), Ok(published));
    }

    #[test]
    fn test_molecule_combine_all_integrity_preservation(
        sources in prop::collection::vec(molecule_with_constraints_strategy(), 0..5),
    ) {
        let published = Molecule::combine_all(&sources);

        prop_assert_eq!(published.edit().try_build(), Ok(published));
    }

    #[test]
    fn test_molecule_split_integrity_preservation(
        source in molecule_with_constraints_strategy(),
    ) {
        for published in source.split() {
            prop_assert_eq!(published.edit().try_build(), Ok(published));
        }
    }

    #[test]
    fn test_molecule_apply_integrity_preservation(
        (source, edits) in transaction_edits_strategy(),
    ) {
        let published = source.apply(edits).map_err(|error| {
            TestCaseError::fail(format!("generated edit application failed: {error}"))
        })?;

        prop_assert_eq!(published.edit().try_build(), Ok(published));
    }

    #[test]
    fn test_constraints_view_mut_push_publication(
        source in molecule_with_constraints_strategy(),
        incoming in raw_constraint_tree_strategy(),
    ) {
        let mut editor = source.edit();
        editor.constraints_mut().push(incoming.clone());
        let expected = editor.try_build();
        let mut actual = source.clone();
        let result = actual.constraints_mut().push(incoming);
        match expected {
            Ok(expected) => {
                prop_assert_eq!(result, Ok(()));
                prop_assert_eq!(actual, expected);
            }
            Err(error) => {
                prop_assert_eq!(result, Err(error));
                prop_assert_eq!(actual, source);
            }
        }
    }

    #[test]
    fn test_constraints_view_mut_extend_publication(
        source in molecule_with_constraints_strategy(),
        incoming in prop::collection::vec(raw_constraint_tree_strategy(), 0..8),
    ) {
        let mut editor = source.edit();
        for constraint in &incoming {
            editor.constraints_mut().push(constraint.clone());
        }
        let expected = editor.try_build();
        let mut actual = source.clone();
        let result = actual.constraints_mut().extend(incoming.into());
        match expected {
            Ok(expected) => {
                prop_assert_eq!(result, Ok(()));
                prop_assert_eq!(actual, expected);
            }
            Err(error) => {
                prop_assert_eq!(result, Err(error));
                prop_assert_eq!(actual, source);
            }
        }
    }

    #[test]
    fn test_constraints_view_mut_replace_publication(
        source in molecule_with_constraints_strategy(),
        incoming in prop::collection::vec(raw_constraint_tree_strategy(), 0..8),
    ) {
        let mut editor = source.edit();
        *editor.constraints_mut() = incoming.clone().into();
        let expected = editor.try_build();
        let mut actual = source.clone();
        let result = actual.constraints_mut().replace(incoming.into());
        match expected {
            Ok(expected) => {
                prop_assert_eq!(result, Ok(()));
                prop_assert_eq!(actual, expected);
            }
            Err(error) => {
                prop_assert_eq!(result, Err(error));
                prop_assert_eq!(actual, source);
            }
        }
    }
}
