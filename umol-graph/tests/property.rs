//! Property-based tests for `umol-graph`.
//!
//! The suite is organized by subject and operation; shared generators live in
//! `strategies`. This test target and `cargo test --test property -- --list`
//! are the authoritative inventory.

const REGRESSION_FILE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/property.proptest-regressions"
);

#[path = "property/strategies.rs"]
mod strategies;

#[path = "property/isotope.rs"]
mod isotope;
#[path = "property/publication.rs"]
mod publication;
#[path = "property/resolve.rs"]
mod resolve;

#[path = "property/project.rs"]
mod project;

#[path = "property/atom_typing.rs"]
mod atom_typing;
