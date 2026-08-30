//! Aggregate-canonicalization properties.
//!
//! The generated domain contains integrity-valid molecules with every entity kind and optional
//! constraints. Independently shuffled complete permutations in every entity namespace supply the
//! dense-remapping action. Exact fixtures and bounded exhaustive minima remain in the unit suite;
//! this module asserts the full normalization/reframe/canonicalize fixpoint and absorption matrix,
//! remapping invariance, equality, and canonical-hash laws without selecting a particular
//! symmetry-equivalent correspondence.

use std::hash::{DefaultHasher, Hash, Hasher};

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{Canonicalize, CanonicalizeContext, Normalize, Reframe};

use crate::strategies::*;

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
            let (with_correspondence, correspondence) = molecule
                .clone()
                .canonicalize_with_correspondence(&context)
                .expect("successful canonicalization returns its correspondence");
            let reframed = molecule
                .remap(&correspondence)
                .reframe()
                .expect("a canonical correspondence preserves molecule integrity");

            prop_assert_eq!(&with_correspondence, &canonical);
            prop_assert_eq!(&reframed, &canonical);
            prop_assert!(molecule.framed_eq_under(&canonical, &correspondence));
        }
    }

    #[test]
    fn test_molecule_canonicalize_standardization(
        scenario in standardization_scenario_strategy(),
    ) {
        let context = context();
        let source = scenario.molecule;
        let normalized = source.clone().normalize().map_err(|_| {
            TestCaseError::fail("generated molecule is intrinsically contradictory")
        })?;
        let reframed = source.clone().reframe().map_err(|_| {
            TestCaseError::fail("generated molecule is intrinsically contradictory")
        })?;
        let canonical = source.clone().canonicalize(&context).map_err(|error| {
            TestCaseError::fail(format!("generated molecule did not canonicalize: {error}"))
        })?;

        prop_assert_eq!(normalized.clone().normalize(), Ok(normalized.clone()));
        prop_assert_eq!(reframed.clone().reframe(), Ok(reframed.clone()));
        prop_assert_eq!(canonical.clone().canonicalize(&context), Ok(canonical.clone()));
        prop_assert_eq!(normalized.clone().reframe(), Ok(reframed.clone()));
        prop_assert_eq!(reframed.clone().normalize(), Ok(reframed.clone()));
        prop_assert_eq!(normalized.canonicalize(&context), Ok(canonical.clone()));
        prop_assert_eq!(reframed.canonicalize(&context), Ok(canonical.clone()));
        prop_assert_eq!(canonical.clone().normalize(), Ok(canonical.clone()));
        prop_assert_eq!(canonical.clone().reframe(), Ok(canonical));
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
        if let Ok(canonical) = molecule.clone().canonicalize(&context) {
            prop_assert_eq!(
                molecule.canonical_hash(&context),
                Ok(structural_hash(&canonical)),
            );
        }
    }

    #[test]
    fn test_molecule_canonicalize_reframed(
        (left, right) in stereo_reframed_molecule_pair_strategy(),
    ) {
        let context = context();

        prop_assert_eq!(
            right.clone().canonicalize(&context),
            left.clone().canonicalize(&context),
        );
        prop_assert!(left.canonical_eq(&right, &context));
        prop_assert_eq!(right.canonical_hash(&context), left.canonical_hash(&context));
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

}
