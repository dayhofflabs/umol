//! Core graph data structures and algorithms for umol.
//!
//! Provides a CSR-based undirected `Graph` (topology only, Arc-shared with
//! copy-on-write), `FixedRelationSet` and `VarRelationSet` for N-ary relations
//! over graph nodes, and graph algorithms (connected components, biconnected
//! components, cycle enumeration, maximum independent set).

pub mod algorithms;
pub mod graph;
pub mod relation;

pub use graph::{EdgeId, Graph, Neighbor, NodeId, Remapping};
pub use relation::{FixedRelationSet, RelationId, VarRelationSet};
