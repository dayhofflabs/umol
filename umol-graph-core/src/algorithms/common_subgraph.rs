//! Common subgraphs of two graphs.
//!
//! The current operations enumerate all or all maximal common subgraphs through
//! modular-product clique search and find maximum common induced or edge
//! subgraphs through McGregor backtracking. See
//! [Bron and Kerbosch (1973)](https://doi.org/10.1145/362342.362367) and
//! [McGregor (1982)](https://doi.org/10.1002/spe.4380120103).

mod enumeration;
mod maximum;

pub use enumeration::{
    CommonSubgraphEnumerationAlgorithm, EmbeddingKind, MaximalCommonSubgraphAlgorithm,
};
pub use maximum::{McesAlgorithm, McisAlgorithm, McsConnectivity};
