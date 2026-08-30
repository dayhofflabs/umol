//! Property-based tests for `umol-graph-ir`.
//!
//! The suite is organized first by subject and, for the larger molecule,
//! reaction, and stereo families, then by operation. Uniform law and surface
//! families remain flat: splitting lattice laws, entity DSL serialization,
//! deltas, or edits solely by file size would obscure the invariant shared by
//! the tests.
//!
//! Shared generators live in `strategies`. This test target and `cargo test
//! --test property -- --list` are the authoritative inventory; a separate
//! README would duplicate coverage information and inevitably drift.

#[path = "property/strategies.rs"]
mod strategies;

#[path = "property/constraint.rs"]
mod constraint;
#[path = "property/delta.rs"]
mod delta;
#[path = "property/edit.rs"]
mod edit;
#[path = "property/entity.rs"]
mod entity;
#[path = "property/frame.rs"]
mod frame;
#[path = "property/lattice.rs"]
mod lattice;
#[path = "property/molecule.rs"]
mod molecule;
#[path = "property/num.rs"]
mod num;
#[path = "property/reaction.rs"]
mod reaction;
#[path = "property/stereo.rs"]
mod stereo;
#[path = "property/substructure.rs"]
mod substructure;
