//! Property-based tests for `umol-graph`.
//!
//! The suite is organized by subject and operation; shared generators live in
//! `strategies`. This test target and `cargo test --test property -- --list`
//! are the authoritative inventory.

const REGRESSION_FILE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/property.proptest-regressions"
);

mod strategies;

mod isotope;
mod publication;
mod resolve;

mod project;

mod transform;

mod atom_typing;

mod stereo;
