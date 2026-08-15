//! Aggregate canonicalization configuration.

use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::CanonicalizeContext;

use super::model::StereoModel;

/// Operational configuration for aggregate canonicalization.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CanonicalizeConfig {
    /// Graph automorphism algorithm used during canonical-frame search.
    pub automorphism_algorithm: AutomorphismAlgorithm,
}

impl Default for CanonicalizeConfig {
    fn default() -> Self {
        Self {
            automorphism_algorithm: AutomorphismAlgorithm::Nauty,
        }
    }
}

impl CanonicalizeConfig {
    /// Combine the operation config with the canonicalization-relevant part of the stereo model.
    pub fn context(&self, model: &StereoModel) -> CanonicalizeContext {
        CanonicalizeContext {
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
    fn test_canonicalize_config_default() {
        assert_eq!(
            CanonicalizeConfig::default(),
            CanonicalizeConfig {
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
            }
        );
    }

    #[rstest]
    #[case::without_para_stereo(false)]
    #[case::with_para_stereo(true)]
    fn test_canonicalize_config_context(#[case] para_stereo: bool) {
        let model = StereoModel {
            para_stereo,
            ..StereoModel::default()
        };

        assert_eq!(
            CanonicalizeConfig::default().context(&model),
            CanonicalizeContext {
                para_stereo,
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
            }
        );
    }
}
