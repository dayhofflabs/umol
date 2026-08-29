//! Molecule description-level properties.
//!
//! Integrity-valid molecules may populate any entity kind and any inline or molecule-level
//! constraint store. Independently shuffled dense bijections in every entity namespace preserve
//! the molecule's projected description level.

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
    fn test_molecule_description_level(
        (molecule, renumbering) in molecule_dense_renumbering_strategy(),
    ) {
        prop_assert_eq!(
            molecule.description_level(),
            molecule.remap(&renumbering).description_level(),
        );
    }
}
