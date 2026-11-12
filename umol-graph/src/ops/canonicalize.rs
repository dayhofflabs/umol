//! Aggregate canonicalization configuration.

use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::CanonicalizeContext;

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
    /// The graph-IR canonicalization context of this config.
    pub fn context(&self) -> CanonicalizeContext {
        CanonicalizeContext {
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
    fn test_canonicalize_config_context() {
        assert_eq!(
            CanonicalizeConfig::default().context(),
            CanonicalizeContext {
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
            }
        );
    }
}
