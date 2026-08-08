//! The featurizer enum: dispatches to a concrete featurizer struct.

use umol_graph_ir::ir::MoleculeAst;

use super::ecfp::EcfpFeaturizer;
use super::feature_set::{CountedFeatureSet, FeatureSet};
use super::morgan::MorganFeaturizer;
use super::wl::WlFeaturizer;

/// A named fingerprint algorithm, dispatching to the concrete featurizer it wraps.
/// The output is an unfolded [`FeatureSet`].
#[derive(Clone, Copy, Debug)]
pub enum Featurizer {
    Wl(WlFeaturizer),
    Ecfp(EcfpFeaturizer),
    Morgan(MorganFeaturizer),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FingerprintError {
    /// The molecule contains values that are not ground terms.
    NotGround,
    /// The reaction deltas cannot be applied consistently.
    Inconsistent,
    /// A fixed-width fingerprint cannot have zero width.
    ZeroWidth,
    /// An operation requires equal-width fingerprints.
    WidthMismatch { left: usize, right: usize },
}

impl Featurizer {
    /// Featurize `mol`, yielding an unfolded [`FeatureSet`].
    pub fn featurize(&self, mol: &MoleculeAst) -> Result<FeatureSet<u64>, FingerprintError> {
        match self {
            Featurizer::Wl(featurizer) => featurizer.featurize(mol),
            Featurizer::Ecfp(featurizer) => featurizer.featurize(mol),
            Featurizer::Morgan(featurizer) => featurizer.featurize(mol),
        }
    }

    /// Featurize `mol`, keeping per-identifier counts.
    pub fn featurize_counted(
        &self,
        mol: &MoleculeAst,
    ) -> Result<CountedFeatureSet<u64>, FingerprintError> {
        match self {
            Featurizer::Wl(featurizer) => featurizer.featurize_counted(mol),
            Featurizer::Ecfp(featurizer) => featurizer.featurize_counted(mol),
            Featurizer::Morgan(featurizer) => featurizer.featurize_counted(mol),
        }
    }
}
