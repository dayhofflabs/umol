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

pub use algorithms::auto::{AutoGroupOrder, Automorphism, AutomorphismAlgorithm};
pub use algorithms::bcc::BiconnectedComponentsAlgorithm;
pub use algorithms::connected::ConnectedComponentsAlgorithm;
pub use algorithms::cycles::{CycleEnumerationAlgorithm, ShortestCycleAlgorithm};
pub use algorithms::matching::{MatchingEnumerationAlgorithm, MaxMatchingAlgorithm, Matching};
pub use algorithms::mis::MaxIndependentSetAlgorithm;
pub use algorithms::subiso::SubgraphIsomorphismAlgorithm;
pub use graph::{EdgeId, Graph, Neighbor, NodeId, Remapping, Subgraph};
pub use relation::{FixedRelationSet, RelationId, VarRelationSet};
pub use union_find::UnionFind;
