//! Molecular fingerprints: featurizers over molecules plus similarity ops.

mod bit_fp;
mod ecfp;
mod feature_set;
mod featurizer;
mod morgan;
mod pattern;
mod reaction;
mod substructure;
mod wl;

pub use bit_fp::BitFp;
pub use ecfp::EcfpFeaturizer;
pub use feature_set::{CountedFeatureSet, FeatureSet, SignedFeatureSet};
pub use featurizer::{Featurizer, FingerprintError};
pub use morgan::MorganFeaturizer;
pub use pattern::{PatternFingerprinter, PATTERN_FP_WIDTH};
pub use reaction::{featurize_reaction, ReactionCombinator, ReactionFingerprint, Side};
pub use substructure::SubstructureFeaturizer;
pub use wl::WlFeaturizer;
