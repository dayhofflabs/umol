//! Core graph data structures and algorithms for umol.
//!
//! Provides a CSR-based undirected `Graph` (topology only, Arc-shared with
//! copy-on-write), `FixedRelationSet` and `VarRelationSet` for N-ary relations
//! over graph nodes, and graph algorithms (connected components, biconnected
//! components, cycle enumeration, maximum independent set).

pub mod algorithms;
pub mod graph;
pub mod relation;
pub mod union_find;

pub use algorithms::vf2::subgraph_isomorphisms;
pub use graph::{EdgeId, Graph, Neighbor, NodeId, Remapping, Subgraph};
pub use relation::{FixedRelationSet, RelationId, VarRelationSet};
pub use union_find::UnionFind;
