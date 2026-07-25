//! Graph matchings.
//!
//! Current operations find maximum and perfect matchings, visit or collect
//! perfect and maximum matchings, and count perfect matchings in a supplied
//! planar embedding. Maximum matching uses Edmonds for general graphs or
//! Hopcroft--Karp for bipartite graphs; planar counting uses Kasteleyn signing
//! and exact Pfaffian evaluation. See [Edmonds
//! (1965)](https://doi.org/10.4153/CJM-1965-045-4),
//! [Hopcroft and Karp (1973)](https://doi.org/10.1137/0202019), and
//! [Kasteleyn (1961)](https://doi.org/10.1016/0031-8914(61)90063-5).

mod count;
mod enumeration;
mod maximum;

pub use count::{FaceBoundary, PlanarEmbedding, PlanarEmbeddingError, PlanarMatchingCountError};
pub use enumeration::MatchingEnumerationAlgorithm;
pub use maximum::{
    Matching, MaximumMatchingAlgorithm, MaximumMatchingError, PerfectMatchingAlgorithm,
};
