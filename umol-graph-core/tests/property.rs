//! Property tests organized by graph operation.
//!
//! Checked-in graph collections live in `corpus`, generated inputs live in
//! `strategy`, and operation-specific validation facilities remain under their
//! operation module.

#[path = "property/corpus.rs"]
mod corpus;
#[path = "property/cycles.rs"]
mod cycles;
#[path = "property/strategy.rs"]
mod strategy;
#[path = "property/subgraph_isomorphism.rs"]
mod subgraph_isomorphism;
