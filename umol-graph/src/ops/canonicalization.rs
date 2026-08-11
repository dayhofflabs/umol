//! Aggregate canonicalization configuration.

use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::CanonicalizationContext;

use super::model::StereoModel;

/// Operational configuration for aggregate canonicalization.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CanonicalizationConfig {
    /// Graph automorphism algorithm used during canonical-frame search.
    pub automorphism_algorithm: AutomorphismAlgorithm,
}

impl Default for CanonicalizationConfig {
    fn default() -> Self {
        Self {
            automorphism_algorithm: AutomorphismAlgorithm::Nauty,
        }
    }
}

impl CanonicalizationConfig {
    /// Combine the operation config with the canonicalization-relevant part of the stereo model.
    pub fn context(&self, model: &StereoModel) -> CanonicalizationContext {
        CanonicalizationContext {
            para_stereo: model.para_stereo,
            automorphism_algorithm: self.automorphism_algorithm,
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    fn test_canonicalization_config_default() {
        assert_eq!(
            CanonicalizationConfig::default(),
            CanonicalizationConfig {
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
            }
        );
    }

    #[rstest]
    #[case::without_para_stereo(false)]
    #[case::with_para_stereo(true)]
    fn test_canonicalization_config_context(#[case] para_stereo: bool) {
        let model = StereoModel {
            para_stereo,
            ..StereoModel::default()
        };

        assert_eq!(
            CanonicalizationConfig::default().context(&model),
            CanonicalizationContext {
                para_stereo,
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
            }
        );
    }
}
