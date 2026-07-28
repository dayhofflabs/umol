//! Property tests organized by graph operation.
//!
//! Checked-in graph collections live in `corpus`, shared generated inputs live
//! in `strategy`, and operation-specific validation facilities remain under
//! their operation module.

#[path = "property/common_subgraph.rs"]
mod common_subgraph;
#[path = "property/corpus.rs"]
mod corpus;
#[path = "property/correspondence.rs"]
mod correspondence;
#[path = "property/cycles.rs"]
mod cycles;
#[path = "property/graph.rs"]
mod graph;
#[path = "property/relation.rs"]
mod relation;
#[path = "property/strategy.rs"]
mod strategy;
#[path = "property/subgraph_isomorphism.rs"]
mod subgraph_isomorphism;
