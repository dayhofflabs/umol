//! Molecular fingerprints: featurizers over molecules plus similarity ops.

mod bit_fp;
mod ecfp;
mod feature_set;
mod featurizer;
mod morgan;
mod wl;

pub use bit_fp::BitFp;
pub use ecfp::{EcfpFeaturizer, ECFP_SEED};
pub use feature_set::FeatureSet;
pub use featurizer::{Featurizer, FingerprintError};
pub use morgan::MorganFeaturizer;
pub use wl::WlFeaturizer;
