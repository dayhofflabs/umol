//! Stereo serialization/semantic properties and bounded exhaustive checks.
//!
//! The semantic and serialization properties share one regression file through
//! `super::REGRESSION_FILE`, so reorganizing them does not orphan minimized failures.
//! The checks exhaust binary coset actions, single-carrier ligand-equality domains and
//! assignments of the selected literature molecules.

const REGRESSION_FILE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/property/stereo.proptest-regressions"
);

#[path = "stereo/exhaustive.rs"]
mod exhaustive;
#[path = "stereo/semantics.rs"]
mod semantics;
#[path = "stereo/serialization.rs"]
mod serialization;
