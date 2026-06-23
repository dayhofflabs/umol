//! The featurizer enum: dispatches to a concrete featurizer struct.

use umol_ast::ast::MoleculeAst;

use super::ecfp::EcfpFeaturizer;
use super::feature_set::FeatureSet;
use super::wl::WlFeaturizer;

/// A named fingerprint algorithm, dispatching to the concrete featurizer it wraps.
/// The output is an unfolded [`FeatureSet`].
#[derive(Clone, Copy, Debug)]
pub enum Featurizer {
    Wl(WlFeaturizer),
    Ecfp(EcfpFeaturizer),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FingerprintError {
    NotGround,
}

impl Featurizer {
    /// Featurize `mol`, yielding an unfolded [`FeatureSet`]. The molecule must be
    /// ground — the precondition shared by every featurizer, checked once here.
    pub fn featurize(&self, mol: &MoleculeAst) -> Result<FeatureSet<u64>, FingerprintError> {
        if !mol.is_ground() {
            return Err(FingerprintError::NotGround);
        }
        Ok(match self {
            Featurizer::Wl(featurizer) => featurizer.featurize(mol),
            Featurizer::Ecfp(featurizer) => featurizer.featurize(mol),
        })
    }
}
