//! Graph matchings.

mod count;
mod enumeration;
mod maximum;

pub use count::{FaceBoundary, PlanarEmbedding, PlanarEmbeddingError, PlanarMatchingCountError};
pub use enumeration::MatchingEnumerationAlgorithm;
pub use maximum::{
    Matching, MaximumMatchingAlgorithm, MaximumMatchingError, PerfectMatchingAlgorithm,
};
