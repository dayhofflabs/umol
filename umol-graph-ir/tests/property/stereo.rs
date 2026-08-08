//! Stereo properties separated into DSL serialization and AST semantics.
//!
//! The parent owns one regression file for the whole subject. Both children use
//! that file through `super::REGRESSION_FILE`, so reorganizing properties by
//! operation does not orphan minimized failures.

const REGRESSION_FILE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/property/stereo.proptest-regressions"
);

#[path = "stereo/semantics.rs"]
mod semantics;
#[path = "stereo/serialization.rs"]
mod serialization;
