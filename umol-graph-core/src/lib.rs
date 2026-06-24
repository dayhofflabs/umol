//! Core graph data structures and algorithms for umol.
//!
//! Provides a CSR-based undirected `Graph` (topology only, Arc-shared with
//! copy-on-write), `FixedRelationSet` and `VarRelationSet` for N-ary relations
//! over graph nodes, and graph algorithms (connected components, biconnected
//! components, cycle enumeration, maximum independent set).

pub mod algorithms;
pub(crate) mod graph;
pub(crate) mod relation;
pub(crate) mod union_find;

pub use algorithms::auto::{AutoGroupOrder, Automorphism, AutomorphismAlgorithm};
pub use algorithms::bcc::BiconnectedComponentsAlgorithm;
pub use algorithms::coloring::BipartitionAlgorithm;
pub use algorithms::connected::ConnectedComponentsAlgorithm;
pub use algorithms::cycles::{CycleEnumerationAlgorithm, ShortestCycleAlgorithm};
pub use algorithms::enumeration::{PathEnumerationAlgorithm, SubgraphEnumerationAlgorithm};
pub use algorithms::matching::{
    Matching, MatchingEnumerationAlgorithm, MaxMatchingAlgorithm, PerfectMatchingAlgorithm,
};
pub use algorithms::mcs::{CommonSubgraph, McesAlgorithm, McisAlgorithm, McsConnectivity};
pub use algorithms::mis::MaxIndependentSetAlgorithm;
pub use algorithms::refine::{
    CircularRefinementAlgorithm, CircularRefinementHash, Refinement, RefinementAlgorithm,
    RefinementHash, RefinementRounds,
};
pub use algorithms::subiso::{SubgraphIsomorphismAlgorithm, ARCMATCH_DEFAULT_PATH_LENGTH};
pub use algorithms::traversal::TraversalAlgorithm;
pub use graph::{EdgeId, Embedding, Graph, Neighbor, NodeId, Remapping};
pub use relation::{
    FactorOrdering, FixedFixedBirelationSet, FixedRelationSet, FixedVarBirelationSet, Ordered,
    ParticipantRefs, RelationId, RelationParticipant, Unordered, VarRelationSet,
    VarVarBirelationSet,
};
pub use union_find::UnionFind;
