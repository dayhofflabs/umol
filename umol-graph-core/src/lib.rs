//! Core graph data structures and algorithms for umol.
//!
//! Provides a CSR-based undirected `Graph` (topology only, Arc-shared with
//! copy-on-write), `FixedRelationSet` and `VarRelationSet` for N-ary relations
//! over graph nodes, and graph algorithms (connected components, biconnected
//! components, cycle enumeration, maximum independent set).

pub mod algorithms;
pub(crate) mod compact;
pub(crate) mod correspondence;
pub(crate) mod digraph;
pub(crate) mod graph;
pub(crate) mod relation;
pub(crate) mod remap;
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
pub use compact::{Compaction, CompactionError, GraphCompaction};
pub use correspondence::{
    Correspondence, CorrespondenceComposeError, CorrespondenceError, GraphCorrespondence,
    GraphCorrespondenceComposeError,
};
pub use digraph::DiGraph;
pub use graph::{EdgeId, Graph, Neighbor, NodeId, SubdividedGraph, SubdivisionNodeSource};
pub use relation::{
    FixedFixedBirelationSet, FixedRelationSet, FixedVarBirelationSet, ParticipantPosition,
    ParticipantRefs, RelationId, RelationParticipant, RelationPullbackCorrespondence,
    RelationPushoutCorrespondence, VarRelationSet, VarVarBirelationSet,
};
pub use remap::{GraphRemapping, Remapping, RemappingError};
pub use rewriting::{
    GraphPushoutCorrespondence, PullbackCorrespondence, PushoutComplementCorrespondence,
};
pub use union_find::UnionFind;
