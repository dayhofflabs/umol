//! Python configuration values for molecular transformations.

use pyo3::prelude::*;
use umol_graph::ops::transform::KekulizeConfig as GraphKekulizeConfig;

use crate::algorithm::MaximumMatchingAlgorithm;

/// Operational configuration for kekulization.
///
/// The algorithm is immutable and defaults to Edmonds, matching the Rust default.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct KekulizeConfig(GraphKekulizeConfig);

#[pymethods]
impl KekulizeConfig {
    #[new]
    #[pyo3(signature = (*, algorithm=MaximumMatchingAlgorithm::Edmonds()))]
    fn new(algorithm: MaximumMatchingAlgorithm) -> Self {
        Self::from_rust(GraphKekulizeConfig::new(algorithm.to_rust()))
    }

    #[getter]
    fn algorithm(&self) -> MaximumMatchingAlgorithm {
        MaximumMatchingAlgorithm::from_rust(self.0.algorithm)
    }

    fn __repr__(&self) -> String {
        format!("KekulizeConfig(algorithm={})", self.algorithm().repr())
    }
}

impl KekulizeConfig {
    pub(crate) fn from_rust(config: GraphKekulizeConfig) -> Self {
        Self(config)
    }

    pub(crate) fn to_rust(self) -> GraphKekulizeConfig {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph::ops::transform::MaximumMatchingAlgorithm as GraphMaximumMatchingAlgorithm;

    use super::*;

    #[rstest]
    #[case::edmonds(GraphKekulizeConfig::default(), MaximumMatchingAlgorithm::Edmonds())]
    #[case::hopcroft_karp(
        GraphKekulizeConfig::new(GraphMaximumMatchingAlgorithm::HopcroftKarp),
        MaximumMatchingAlgorithm::HopcroftKarp()
    )]
    fn test_kekulize_config_from_rust(
        #[case] config: GraphKekulizeConfig,
        #[case] expected_algorithm: MaximumMatchingAlgorithm,
    ) {
        let bound = KekulizeConfig::from_rust(config);
        assert_eq!(bound.algorithm(), expected_algorithm);
        assert_eq!(bound.to_rust(), config);
    }
}
