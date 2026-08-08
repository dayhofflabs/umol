//! Python bindings for ring-enumeration configuration.

use pyo3::prelude::*;
use umol_graph_ir::ir::RingConfig as GraphIrRingConfig;

use crate::algorithm::{RelevantCycleEnumerationAlgorithm, SimpleCycleEnumerationAlgorithm};

/// Algorithms used to compute each supported ring-set kind.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RingConfig(GraphIrRingConfig);

#[pymethods]
impl RingConfig {
    #[new]
    #[pyo3(signature = (*, simple_cycle_algorithm=None, relevant_cycle_algorithm=None))]
    fn new(
        simple_cycle_algorithm: Option<SimpleCycleEnumerationAlgorithm>,
        relevant_cycle_algorithm: Option<RelevantCycleEnumerationAlgorithm>,
    ) -> Self {
        let defaults = Self::default().0;
        Self(GraphIrRingConfig {
            simple_cycle_algorithm: simple_cycle_algorithm
                .map_or(defaults.simple_cycle_algorithm, |algorithm| {
                    algorithm.to_rust()
                }),
            relevant_cycle_algorithm: relevant_cycle_algorithm
                .map_or(defaults.relevant_cycle_algorithm, |algorithm| {
                    algorithm.to_rust()
                }),
        })
    }

    #[getter]
    fn simple_cycle_algorithm(&self) -> SimpleCycleEnumerationAlgorithm {
        SimpleCycleEnumerationAlgorithm::from_rust(self.0.simple_cycle_algorithm)
    }

    #[getter]
    fn relevant_cycle_algorithm(&self) -> RelevantCycleEnumerationAlgorithm {
        RelevantCycleEnumerationAlgorithm::from_rust(self.0.relevant_cycle_algorithm)
    }

    pub(crate) fn __repr__(&self) -> String {
        format!(
            "RingConfig(simple_cycle_algorithm={}, relevant_cycle_algorithm={})",
            self.simple_cycle_algorithm().repr(),
            self.relevant_cycle_algorithm().repr(),
        )
    }
}

impl RingConfig {
    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion API for ring-consuming configurations"
    )]
    pub(crate) fn from_rust(config: GraphIrRingConfig) -> Self {
        Self(config)
    }

    #[allow(
        dead_code,
        reason = "Python-to-Rust conversion API for ring-consuming configurations"
    )]
    pub(crate) fn to_rust(self) -> GraphIrRingConfig {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_core::{
        RelevantCycleEnumerationAlgorithm as GraphCoreRelevantCycleEnumerationAlgorithm,
        SimpleCycleEnumerationAlgorithm as GraphCoreSimpleCycleEnumerationAlgorithm,
    };

    use super::*;

    #[rstest]
    #[case::default(
        None,
        None,
        GraphIrRingConfig::default(),
        "RingConfig(simple_cycle_algorithm=SimpleCycleEnumerationAlgorithm.ReadTarjan(), relevant_cycle_algorithm=RelevantCycleEnumerationAlgorithm.Vismara())",
    )]
    #[case::explicit(
        Some(SimpleCycleEnumerationAlgorithm::ReadTarjan()),
        Some(RelevantCycleEnumerationAlgorithm::Vismara()),
        GraphIrRingConfig {
            simple_cycle_algorithm: GraphCoreSimpleCycleEnumerationAlgorithm::ReadTarjan,
            relevant_cycle_algorithm: GraphCoreRelevantCycleEnumerationAlgorithm::Vismara,
        },
        "RingConfig(simple_cycle_algorithm=SimpleCycleEnumerationAlgorithm.ReadTarjan(), relevant_cycle_algorithm=RelevantCycleEnumerationAlgorithm.Vismara())",
    )]
    fn test_ring_config_new(
        #[case] simple_cycle_algorithm: Option<SimpleCycleEnumerationAlgorithm>,
        #[case] relevant_cycle_algorithm: Option<RelevantCycleEnumerationAlgorithm>,
        #[case] expected: GraphIrRingConfig,
        #[case] expected_repr: &str,
    ) {
        let config = RingConfig::new(simple_cycle_algorithm, relevant_cycle_algorithm);

        assert_eq!(config.to_rust(), expected);
        assert_eq!(config.__repr__(), expected_repr);
    }

    #[rstest]
    #[case::default(GraphIrRingConfig::default())]
    fn test_ring_config_from_rust(#[case] config: GraphIrRingConfig) {
        assert_eq!(RingConfig::from_rust(config).to_rust(), config);
    }
}
