//! Property-based tests for umol-ast, split by area. Shared generators live in
//! the `strategies` module; each sibling module holds one area's proptests.

#[path = "property/strategies.rs"]
mod strategies;

#[path = "property/delta.rs"]
mod delta;
#[path = "property/edit.rs"]
mod edit;
#[path = "property/entity.rs"]
mod entity;
#[path = "property/lattice.rs"]
mod lattice;
#[path = "property/molecule.rs"]
mod molecule;
#[path = "property/stereo.rs"]
mod stereo;
#[path = "property/substructure.rs"]
mod substructure;
#[path = "property/value.rs"]
mod value;
