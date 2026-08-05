//! Core graph data structures and algorithms for umol.
//!
//! Provides a CSR-based undirected `Graph` (topology only, Arc-shared with
//! copy-on-write), `FixedRelationSet` and `VarRelationSet` for N-ary relations
//! over graph nodes, and graph algorithms (connected components, biconnected
//! components, cycle enumeration, maximum independent set).

pub mod algorithms;
pub(crate) mod correspondence;
pub(crate) mod digraph;
pub(crate) mod graph;
pub(crate) mod relation;
pub(crate) mod rewriting;
pub(crate) mod union_find;

pub use algorithms::automorphism::{
    AutomorphismAlgorithm, AutomorphismGroupOrder, AutomorphismOutput,
};
pub use algorithms::bipartite::{BipartitionAlgorithm, NonBipartiteGraphError};
pub use algorithms::common_subgraph::{
    CommonSubgraphEnumerationAlgorithm, EmbeddingKind, MaximalCommonSubgraphAlgorithm,
    McesAlgorithm, McisAlgorithm, McsConnectivity,
};
pub use algorithms::connected_subgraphs::SubgraphEnumerationAlgorithm;
pub use algorithms::connectivity::{BiconnectedComponentsAlgorithm, ConnectedComponentsAlgorithm};
pub use algorithms::cycles::{
    Cycle, MinimumCycleBasis, MinimumCycleBasisAlgorithm, NonSimpleGraphError, RelevantCycleCount,
    RelevantCycleEnumerationAlgorithm, ShortestCycleAlgorithm, SimpleCycleEnumerationAlgorithm,
    UniqueRingFamilies, UniqueRingFamily, UniqueRingFamilyAlgorithm, UniqueRingFamilyId,
};
pub use algorithms::independent_set::MaximumIndependentSetAlgorithm;
pub use algorithms::matching::{
    BipartiteMaximumMatchingAlgorithm, FaceBoundary, GeneralMaximumMatchingAlgorithm, Matching,
    MatchingEnumerationAlgorithm, PerfectMatchingAlgorithm, PlanarEmbedding, PlanarEmbeddingError,
    PlanarMatchingCountError,
};
pub use algorithms::paths::PathEnumerationAlgorithm;
pub use algorithms::refinement::{
    CircularRefinementAlgorithm, CircularRefinementHash, Refinement, RefinementAlgorithm,
    RefinementHash, RefinementRounds,
};
pub use algorithms::subgraph_isomorphism::{
    SubgraphIsomorphismAlgorithm, ARCMATCH_DEFAULT_PATH_LENGTH,
};
pub use algorithms::topological_sort::TopologicalSortAlgorithm;
pub use algorithms::traversal::TraversalAlgorithm;
pub use correspondence::{Correspondence, CorrespondenceError, GraphCorrespondence};
pub use digraph::DiGraph;
pub use graph::{
    compact_edge_vec, compact_node_vec, Compaction, EdgeId, Graph, Neighbor, NodeId, Remapping,
    SubdividedGraph, SubdivisionNodeSource,
};
pub use relation::{
    BiRelationData, FactorOrdering, FixedFixedBirelationSet, FixedRelationSet,
    FixedVarBirelationSet, Ordered, ParticipantAnchor, ParticipantPosition, ParticipantRefs,
    RelationData, RelationId, RelationParticipant, RelationPullback, RelationPushout, Unordered,
    VarRelationSet, VarVarBirelationSet,
};
pub use rewriting::{Pullback, Pushout, PushoutComplement};
pub use union_find::UnionFind;
