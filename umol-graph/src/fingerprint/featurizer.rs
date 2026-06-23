//! The featurizer enum: dispatches to a concrete featurizer struct.

use umol_ast::ast::MoleculeAst;

use super::feature_set::FeatureSet;
use super::wl::WlFeaturizer;

/// A named fingerprint algorithm, dispatching to the concrete featurizer it wraps.
/// The output is an unfolded [`FeatureSet`].
#[derive(Clone, Copy, Debug)]
pub enum Featurizer {
    Wl(WlFeaturizer),
}

/// A molecule could not be featurized. Fingerprints are defined only on ground
/// molecules; the featurizer never coerces a non-ground field to a default.
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
        })
    }
}

#[cfg(test)]
mod tests {
    use rstest::{fixture, rstest};
    use umol_ast::{mol, mol_ground};
    use umol_graph_core::{RefinementRounds, RefinementXxh3Scheme};

    use super::*;

    #[fixture]
    fn wl() -> Featurizer {
        Featurizer::Wl(WlFeaturizer {
            rounds: RefinementRounds::Fixed(3),
            scheme: RefinementXxh3Scheme::albatross(),
        })
    }

    #[rstest]
    #[case::ethane(r#"{:atoms ["C #h3" "C #h3"] :bonds [[0 1 "1"]]}"#)]
    fn test_featurizer_featurize(wl: Featurizer, #[case] edn: &str) {
        let mol = mol_ground!(edn);
        let expected = match wl {
            Featurizer::Wl(inner) => inner.featurize(&mol),
        };
        assert_eq!(wl.featurize(&mol), Ok(expected));
    }

    #[rstest]
    #[case::non_ground_atom(r#"{:atoms ["C"] :bonds []}"#)]
    fn test_featurizer_featurize_error(wl: Featurizer, #[case] edn: &str) {
        assert_eq!(wl.featurize(&mol!(edn)), Err(FingerprintError::NotGround));
    }
}
