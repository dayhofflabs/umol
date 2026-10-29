//! Stereo serialization/semantic properties and bounded coset-action references.
//!
//! The semantic and serialization properties share one regression file through
//! `super::REGRESSION_FILE`, so reorganizing them does not orphan minimized failures.
//! The reference checks both binary kinds over every four-position permutation.

const REGRESSION_FILE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/property/stereo.proptest-regressions"
);

#[path = "stereo/reference.rs"]
mod reference;
#[path = "stereo/semantics.rs"]
mod semantics;
#[path = "stereo/serialization.rs"]
mod serialization;
