//! Property-based tests for `umol-graph`.
//!
//! The suite is organized by subject and operation; shared generators live in
//! `strategies`. This test target and `cargo test --test property -- --list`
//! are the authoritative inventory.

#[path = "property/strategies.rs"]
mod strategies;

#[path = "property/resolve.rs"]
mod resolve;
