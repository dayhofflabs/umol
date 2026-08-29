//! Reaction properties grouped by operation.
//!
//! The parent owns one regression file for the whole subject. Every child uses
//! that file through `super::REGRESSION_FILE`, so reorganizing properties by
//! operation does not orphan minimized failures.

const REGRESSION_FILE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/property/reaction.proptest-regressions"
);

#[path = "reaction/application.rs"]
mod application;
#[path = "reaction/canonicalize.rs"]
mod canonicalize;
#[path = "reaction/composition.rs"]
mod composition;
#[path = "reaction/lifecycle.rs"]
mod lifecycle;
#[path = "reaction/malformed.rs"]
mod malformed;
#[path = "reaction/metadata.rs"]
mod metadata;
#[path = "reaction/reframe.rs"]
mod reframe;
#[path = "reaction/serialization.rs"]
mod serialization;
#[path = "reaction/span.rs"]
mod span;
