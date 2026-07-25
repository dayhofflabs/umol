//! Graph algorithms grouped by graph problem.
//!
//! Each problem module owns its operation entry points, algorithm selectors,
//! and result types. A `visit_*` operation emits results incrementally,
//! `enumerate_*` collects them eagerly, and an eventual `iter_*` operation
//! returns resumable iteration state; these forms stay together in the problem
//! module.

pub mod automorphism;
pub mod bipartite;
pub mod common_subgraph;
pub mod connected_subgraphs;
pub mod connectivity;
pub mod cycles;
pub mod independent_set;
pub mod matching;
pub mod paths;
pub mod refinement;
pub mod subgraph_isomorphism;
pub mod topological_sort;
pub mod traversal;
