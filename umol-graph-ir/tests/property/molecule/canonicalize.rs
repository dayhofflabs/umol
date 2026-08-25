//! Aggregate-canonicalization properties.
//!
//! The generated domain contains integrity-valid molecules with every entity family and optional
//! constraints. Independently shuffled complete permutations in every entity namespace supply the
//! dense-remapping action. Exact fixtures and bounded exhaustive minima remain in the unit suite;
//! this module asserts the public idempotence, remapping-invariance, equality, and canonical-hash
//! laws without selecting a particular symmetry-equivalent correspondence.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{Canonicalize, CanonicalizeContext, CanonicalizeLevel};

use crate::strategies::*;

fn context() -> CanonicalizeContext {
    CanonicalizeContext {
        para_stereo: false,
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    }
}

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_molecule_canonicalize(
        (molecule, renumbering) in molecule_dense_renumbering_strategy(),
    ) {
        let context = context();
        let renumbered = molecule.remap(&renumbering);
        let canonical = molecule.clone().canonicalize(&context);
        let renumbered_canonical = renumbered.canonicalize(&context);

        prop_assert_eq!(&renumbered_canonical, &canonical);
        if let Ok(canonical) = canonical {
            let (_, correspondence) = molecule
                .clone()
                .canonicalize_with_correspondence(&context)
                .expect("successful canonicalization returns its correspondence");

            prop_assert!(molecule.equiv_under(&canonical, &correspondence));
            prop_assert_eq!(canonical.clone().canonicalize(&context), Ok(canonical));
        }
    }

    #[test]
    fn test_molecule_canonical_hash(
        (molecule, renumbering) in molecule_dense_renumbering_strategy(),
    ) {
        let context = context();
        let renumbered = molecule.remap(&renumbering);

        prop_assert_eq!(
            molecule.clone().canonical_hash(&context),
            renumbered.clone().canonical_hash(&context),
        );
        for level in [
            CanonicalizeLevel::Topology,
            CanonicalizeLevel::Constitution,
            CanonicalizeLevel::Structure,
            CanonicalizeLevel::Full,
        ] {
            prop_assert_eq!(
                molecule.clone().canonical_hash_by(level, &context),
                renumbered.clone().canonical_hash_by(level, &context),
            );
        }
        prop_assert_eq!(
            molecule.clone().canonical_hash_by(CanonicalizeLevel::Full, &context),
            molecule.canonical_hash(&context),
        );
    }

    #[test]
    fn test_molecule_canonical_eq(
        (molecule, renumbering) in molecule_dense_renumbering_strategy(),
    ) {
        let context = context();
        let renumbered = molecule.remap(&renumbering);
        let canonical = molecule.clone().canonicalize(&context);

        prop_assert!(molecule.canonical_eq(&molecule, &context));
        prop_assert!(molecule.canonical_eq(&renumbered, &context));
        prop_assert_eq!(
            molecule.canonical_eq(&renumbered, &context),
            renumbered.canonical_eq(&molecule, &context),
        );
        if let Ok(canonical) = canonical {
            prop_assert!(renumbered.canonical_eq(&canonical, &context));
            prop_assert!(molecule.canonical_eq(&canonical, &context));
        }
    }

    #[test]
    fn test_molecule_canonical_eq_by(
        (molecule, renumbering) in molecule_dense_renumbering_strategy(),
    ) {
        let context = context();
        let renumbered = molecule.remap(&renumbering);

        for level in [
            CanonicalizeLevel::Topology,
            CanonicalizeLevel::Constitution,
            CanonicalizeLevel::Structure,
            CanonicalizeLevel::Full,
        ] {
            let canonical = molecule.clone().canonicalize_by(level, &context);

            prop_assert!(molecule.canonical_eq_by(&molecule, level, &context));
            prop_assert!(molecule.canonical_eq_by(&renumbered, level, &context));
            prop_assert_eq!(
                molecule.canonical_eq_by(&renumbered, level, &context),
                renumbered.canonical_eq_by(&molecule, level, &context),
            );
            if let Ok(canonical) = canonical {
                prop_assert!(renumbered.canonical_eq_by(&canonical, level, &context));
                prop_assert!(molecule.canonical_eq_by(&canonical, level, &context));
            }
        }
        prop_assert_eq!(
            molecule.canonical_eq_by(&renumbered, CanonicalizeLevel::Full, &context),
            molecule.canonical_eq(&renumbered, &context),
        );
    }
}
