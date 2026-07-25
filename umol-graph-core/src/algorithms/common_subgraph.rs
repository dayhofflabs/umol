//! Common subgraphs of two graphs.

mod enumeration;
mod maximum;

pub use enumeration::{
    CommonSubgraphEnumerationAlgorithm, EmbeddingKind, MaximalCommonSubgraphAlgorithm,
};
pub use maximum::{McesAlgorithm, McisAlgorithm, McsConnectivity};
