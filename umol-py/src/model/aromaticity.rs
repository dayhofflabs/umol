//! Python bindings for aromaticity model and operation-configuration values.

use pyo3::prelude::*;
use umol_graph::ops::aromaticity::AromaticityConfig as GraphAromaticityConfig;
use umol_graph::ops::model::{
    AromaticityModel as GraphAromaticityModel, RingLimits as GraphRingLimits,
};

use super::ElementScope;
use crate::algorithm::{ConnectedComponentsAlgorithm, MaximumIndependentSetAlgorithm};
use crate::ring::RingConfig;

/// Algorithms used to perform aromaticity perception.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct AromaticityConfig(GraphAromaticityConfig);

#[pymethods]
impl AromaticityConfig {
    #[new]
    #[pyo3(signature = (*, ring_config=RingConfig::default(), connected_components_algorithm=ConnectedComponentsAlgorithm::Bfs(), maximum_independent_set_algorithm=MaximumIndependentSetAlgorithm::BranchAndBound()))]
    fn new(
        ring_config: RingConfig,
        connected_components_algorithm: ConnectedComponentsAlgorithm,
        maximum_independent_set_algorithm: MaximumIndependentSetAlgorithm,
    ) -> Self {
        Self(GraphAromaticityConfig {
            ring_config: ring_config.to_rust(),
            connected_components_algorithm: connected_components_algorithm.to_rust(),
            maximum_independent_set_algorithm: maximum_independent_set_algorithm.to_rust(),
        })
    }

    #[getter]
    fn ring_config(&self) -> RingConfig {
        RingConfig::from_rust(self.0.ring_config)
    }

    #[getter]
    fn connected_components_algorithm(&self) -> ConnectedComponentsAlgorithm {
        ConnectedComponentsAlgorithm::from_rust(self.0.connected_components_algorithm)
    }

    #[getter]
    fn maximum_independent_set_algorithm(&self) -> MaximumIndependentSetAlgorithm {
        MaximumIndependentSetAlgorithm::from_rust(self.0.maximum_independent_set_algorithm)
    }

    pub(crate) fn __repr__(&self) -> String {
        format!(
            "AromaticityConfig(ring_config={}, connected_components_algorithm={}, maximum_independent_set_algorithm={})",
            self.ring_config().__repr__(),
            self.connected_components_algorithm().repr(),
            self.maximum_independent_set_algorithm().repr(),
        )
    }
}

impl AromaticityConfig {
    pub(crate) fn from_rust(config: GraphAromaticityConfig) -> Self {
        Self(config)
    }

    pub(crate) fn to_rust(self) -> GraphAromaticityConfig {
        self.0
    }
}

/// Ring-size and fused-ring search bounds for aromaticity perception.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RingLimits(GraphRingLimits);

#[pymethods]
impl RingLimits {
    #[new]
    #[pyo3(signature = (*, min_ring_size=3, max_ring_size=22, include_fused=true, max_fused_combination=6, max_fused_search=10_000))]
    fn new(
        min_ring_size: usize,
        max_ring_size: usize,
        include_fused: bool,
        max_fused_combination: usize,
        max_fused_search: usize,
    ) -> Self {
        Self(GraphRingLimits {
            min_ring_size,
            max_ring_size,
            include_fused,
            max_fused_combination,
            max_fused_search,
        })
    }

    #[getter]
    fn min_ring_size(&self) -> usize {
        self.0.min_ring_size
    }

    #[getter]
    fn max_ring_size(&self) -> usize {
        self.0.max_ring_size
    }

    #[getter]
    fn include_fused(&self) -> bool {
        self.0.include_fused
    }

    #[getter]
    fn max_fused_combination(&self) -> usize {
        self.0.max_fused_combination
    }

    #[getter]
    fn max_fused_search(&self) -> usize {
        self.0.max_fused_search
    }

    fn __repr__(&self) -> String {
        format!(
            "RingLimits(min_ring_size={}, max_ring_size={}, include_fused={}, max_fused_combination={}, max_fused_search={})",
            self.0.min_ring_size,
            self.0.max_ring_size,
            if self.0.include_fused { "True" } else { "False" },
            self.0.max_fused_combination,
            self.0.max_fused_search,
        )
    }
}

impl RingLimits {
    pub(crate) fn from_rust(limits: &GraphRingLimits) -> Self {
        Self(limits.clone())
    }

    #[allow(
        dead_code,
        reason = "Python-to-Rust conversion API for aggregate model configuration"
    )]
    pub(crate) fn to_rust(&self) -> &GraphRingLimits {
        &self.0
    }
}

/// Aromaticity perception model and its model parameters.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq)]
pub enum AromaticityModel {
    /// Hückel-rule aromaticity over rings within the configured limits.
    #[pyo3(constructor = (*, scope, ring_limits))]
    HueckelRule {
        scope: ElementScope,
        ring_limits: RingLimits,
    },
    /// Hückel molecular-orbital aromaticity with a stabilization threshold.
    #[pyo3(constructor = (*, scope, stabilization_threshold))]
    Hmo {
        scope: ElementScope,
        stabilization_threshold: f64,
    },
    /// Clar aromaticity over rings within the configured limits.
    #[pyo3(constructor = (*, scope, ring_limits))]
    Clar {
        scope: ElementScope,
        ring_limits: RingLimits,
    },
}

#[pymethods]
impl AromaticityModel {
    #[staticmethod]
    fn daylight() -> Self {
        Self::from_rust(&GraphAromaticityModel::daylight())
    }

    #[staticmethod]
    fn mdl() -> Self {
        Self::from_rust(&GraphAromaticityModel::mdl())
    }

    #[staticmethod]
    fn permissive() -> Self {
        Self::from_rust(&GraphAromaticityModel::permissive())
    }

    pub(crate) fn __repr__(&self) -> String {
        match self {
            Self::HueckelRule { scope, ring_limits } => format!(
                "AromaticityModel.HueckelRule(scope={}, ring_limits={})",
                scope.__repr__(),
                ring_limits.__repr__(),
            ),
            Self::Hmo {
                scope,
                stabilization_threshold,
            } => format!(
                "AromaticityModel.Hmo(scope={}, stabilization_threshold={stabilization_threshold})",
                scope.__repr__(),
            ),
            Self::Clar { scope, ring_limits } => format!(
                "AromaticityModel.Clar(scope={}, ring_limits={})",
                scope.__repr__(),
                ring_limits.__repr__(),
            ),
        }
    }
}

impl AromaticityModel {
    pub(crate) fn from_rust(model: &GraphAromaticityModel) -> Self {
        match model {
            GraphAromaticityModel::HueckelRule { scope, ring_limits } => Self::HueckelRule {
                scope: ElementScope::from_rust(scope),
                ring_limits: RingLimits::from_rust(ring_limits),
            },
            GraphAromaticityModel::Hmo {
                scope,
                stabilization_threshold,
            } => Self::Hmo {
                scope: ElementScope::from_rust(scope),
                stabilization_threshold: *stabilization_threshold,
            },
            GraphAromaticityModel::Clar { scope, ring_limits } => Self::Clar {
                scope: ElementScope::from_rust(scope),
                ring_limits: RingLimits::from_rust(ring_limits),
            },
        }
    }

    #[allow(
        dead_code,
        reason = "Python-to-Rust conversion API for ChemistryModel configuration"
    )]
    pub(crate) fn to_rust(&self) -> GraphAromaticityModel {
        match self {
            Self::HueckelRule { scope, ring_limits } => GraphAromaticityModel::HueckelRule {
                scope: scope.to_rust(),
                ring_limits: ring_limits.to_rust().clone(),
            },
            Self::Hmo {
                scope,
                stabilization_threshold,
            } => GraphAromaticityModel::Hmo {
                scope: scope.to_rust(),
                stabilization_threshold: *stabilization_threshold,
            },
            Self::Clar { scope, ring_limits } => GraphAromaticityModel::Clar {
                scope: scope.to_rust(),
                ring_limits: ring_limits.to_rust().clone(),
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_chem::element::Element as ChemElement;
    use umol_graph::ops::aromaticity::AromaticityConfig as GraphAromaticityConfig;
    use umol_graph::ops::model::ElementScope as GraphElementScope;

    use super::*;

    #[rstest]
    #[case::default(
        RingConfig::default(),
        ConnectedComponentsAlgorithm::Bfs(),
        MaximumIndependentSetAlgorithm::BranchAndBound(),
        GraphAromaticityConfig::default()
    )]
    fn test_aromaticity_config_new(
        #[case] ring_config: RingConfig,
        #[case] connected_components_algorithm: ConnectedComponentsAlgorithm,
        #[case] maximum_independent_set_algorithm: MaximumIndependentSetAlgorithm,
        #[case] expected: GraphAromaticityConfig,
    ) {
        assert_eq!(
            AromaticityConfig::new(
                ring_config,
                connected_components_algorithm,
                maximum_independent_set_algorithm,
            )
            .to_rust(),
            expected
        );
    }

    #[rstest]
    #[case::default(
        AromaticityConfig::default(),
        "AromaticityConfig(ring_config=RingConfig(simple_cycle_algorithm=SimpleCycleEnumerationAlgorithm.ReadTarjan(), relevant_cycle_algorithm=RelevantCycleEnumerationAlgorithm.Vismara()), connected_components_algorithm=ConnectedComponentsAlgorithm.Bfs(), maximum_independent_set_algorithm=MaximumIndependentSetAlgorithm.BranchAndBound())"
    )]
    fn test_aromaticity_config_repr(#[case] config: AromaticityConfig, #[case] expected: &str) {
        assert_eq!(config.__repr__(), expected);
    }

    #[rstest]
    #[case::default(GraphAromaticityConfig::default())]
    fn test_aromaticity_config_from_rust(#[case] config: GraphAromaticityConfig) {
        assert_eq!(AromaticityConfig::from_rust(config).to_rust(), config);
    }

    #[rstest]
    #[case::default(3, 22, true, 6, 10_000, GraphRingLimits::default())]
    #[case::zero(0, 0, false, 0, 0, GraphRingLimits {
        min_ring_size: 0,
        max_ring_size: 0,
        include_fused: false,
        max_fused_combination: 0,
        max_fused_search: 0,
    })]
    #[case::nondefault(5, 18, false, 4, 2_500, GraphRingLimits {
        min_ring_size: 5,
        max_ring_size: 18,
        include_fused: false,
        max_fused_combination: 4,
        max_fused_search: 2_500,
    })]
    fn test_ring_limits_new(
        #[case] min_ring_size: usize,
        #[case] max_ring_size: usize,
        #[case] include_fused: bool,
        #[case] max_fused_combination: usize,
        #[case] max_fused_search: usize,
        #[case] expected: GraphRingLimits,
    ) {
        assert_eq!(
            RingLimits::new(
                min_ring_size,
                max_ring_size,
                include_fused,
                max_fused_combination,
                max_fused_search,
            )
            .0,
            expected
        );
    }

    #[rstest]
    #[case::default(
        RingLimits::new(3, 22, true, 6, 10_000),
        "RingLimits(min_ring_size=3, max_ring_size=22, include_fused=True, max_fused_combination=6, max_fused_search=10000)"
    )]
    #[case::nondefault(
        RingLimits::new(5, 18, false, 4, 2_500),
        "RingLimits(min_ring_size=5, max_ring_size=18, include_fused=False, max_fused_combination=4, max_fused_search=2500)"
    )]
    fn test_ring_limits_repr(#[case] limits: RingLimits, #[case] expected: &str) {
        assert_eq!(limits.__repr__(), expected);
    }

    #[rstest]
    #[case::default(GraphRingLimits::default())]
    #[case::nondefault(GraphRingLimits {
        min_ring_size: 5,
        max_ring_size: 18,
        include_fused: false,
        max_fused_combination: 4,
        max_fused_search: 2_500,
    })]
    fn test_ring_limits_from_rust(#[case] limits: GraphRingLimits) {
        assert_eq!(RingLimits::from_rust(&limits).0, limits);
    }

    #[rstest]
    #[case::default(RingLimits::new(3, 22, true, 6, 10_000))]
    #[case::nondefault(RingLimits::new(5, 18, false, 4, 2_500))]
    fn test_ring_limits_to_rust(#[case] limits: RingLimits) {
        assert_eq!(limits.to_rust(), &limits.0);
    }

    #[rstest]
    fn test_aromaticity_model_daylight() {
        assert_eq!(
            AromaticityModel::daylight(),
            AromaticityModel::HueckelRule {
                scope: ElementScope::AllowList {
                    elements: vec![
                        ChemElement::C.into(),
                        ChemElement::N.into(),
                        ChemElement::O.into(),
                        ChemElement::S.into(),
                        ChemElement::Se.into(),
                        ChemElement::As.into(),
                    ],
                },
                ring_limits: RingLimits::new(3, 22, true, 6, 10_000),
            }
        );
    }

    #[rstest]
    fn test_aromaticity_model_mdl() {
        assert_eq!(
            AromaticityModel::mdl(),
            AromaticityModel::HueckelRule {
                scope: ElementScope::AllowList {
                    elements: vec![ChemElement::C.into(), ChemElement::N.into()],
                },
                ring_limits: RingLimits::new(6, 22, true, 6, 10_000),
            }
        );
    }

    #[rstest]
    fn test_aromaticity_model_permissive() {
        assert_eq!(
            AromaticityModel::permissive(),
            AromaticityModel::HueckelRule {
                scope: ElementScope::Any {},
                ring_limits: RingLimits::new(3, 22, true, 6, 10_000),
            }
        );
    }

    #[rstest]
    #[case::hueckel_rule(
        AromaticityModel::HueckelRule {
            scope: ElementScope::Any {},
            ring_limits: RingLimits::new(4, 18, false, 3, 2_000),
        },
        "AromaticityModel.HueckelRule(scope=ElementScope.Any(), ring_limits=RingLimits(min_ring_size=4, max_ring_size=18, include_fused=False, max_fused_combination=3, max_fused_search=2000))"
    )]
    #[case::hmo(
        AromaticityModel::Hmo {
            scope: ElementScope::AllowList {
                elements: vec![ChemElement::C.into(), ChemElement::N.into()],
            },
            stabilization_threshold: 0.375,
        },
        "AromaticityModel.Hmo(scope=ElementScope.AllowList([Element('C'), Element('N')]), stabilization_threshold=0.375)"
    )]
    #[case::clar(
        AromaticityModel::Clar {
            scope: ElementScope::AllowList {
                elements: vec![ChemElement::C.into()],
            },
            ring_limits: RingLimits::new(6, 14, true, 4, 1_500),
        },
        "AromaticityModel.Clar(scope=ElementScope.AllowList([Element('C')]), ring_limits=RingLimits(min_ring_size=6, max_ring_size=14, include_fused=True, max_fused_combination=4, max_fused_search=1500))"
    )]
    fn test_aromaticity_model_repr(#[case] model: AromaticityModel, #[case] expected: &str) {
        assert_eq!(model.__repr__(), expected);
    }

    #[rstest]
    #[case::hueckel_rule(
        GraphAromaticityModel::HueckelRule {
            scope: GraphElementScope::Any,
            ring_limits: GraphRingLimits {
                min_ring_size: 4,
                max_ring_size: 18,
                include_fused: false,
                max_fused_combination: 3,
                max_fused_search: 2_000,
            },
        },
        AromaticityModel::HueckelRule {
            scope: ElementScope::Any {},
            ring_limits: RingLimits::new(4, 18, false, 3, 2_000),
        }
    )]
    #[case::hmo(
        GraphAromaticityModel::Hmo {
            scope: GraphElementScope::AllowList(vec![ChemElement::C, ChemElement::N]),
            stabilization_threshold: 0.375,
        },
        AromaticityModel::Hmo {
            scope: ElementScope::AllowList {
                elements: vec![ChemElement::C.into(), ChemElement::N.into()],
            },
            stabilization_threshold: 0.375,
        }
    )]
    #[case::clar(
        GraphAromaticityModel::Clar {
            scope: GraphElementScope::AllowList(vec![ChemElement::C]),
            ring_limits: GraphRingLimits {
                min_ring_size: 6,
                max_ring_size: 14,
                include_fused: true,
                max_fused_combination: 4,
                max_fused_search: 1_500,
            },
        },
        AromaticityModel::Clar {
            scope: ElementScope::AllowList {
                elements: vec![ChemElement::C.into()],
            },
            ring_limits: RingLimits::new(6, 14, true, 4, 1_500),
        }
    )]
    fn test_aromaticity_model_from_rust(
        #[case] model: GraphAromaticityModel,
        #[case] expected: AromaticityModel,
    ) {
        assert_eq!(AromaticityModel::from_rust(&model), expected);
    }

    #[rstest]
    #[case::hueckel_rule(
        AromaticityModel::HueckelRule {
            scope: ElementScope::Any {},
            ring_limits: RingLimits::new(4, 18, false, 3, 2_000),
        },
        GraphAromaticityModel::HueckelRule {
            scope: GraphElementScope::Any,
            ring_limits: GraphRingLimits {
                min_ring_size: 4,
                max_ring_size: 18,
                include_fused: false,
                max_fused_combination: 3,
                max_fused_search: 2_000,
            },
        }
    )]
    #[case::hmo(
        AromaticityModel::Hmo {
            scope: ElementScope::AllowList {
                elements: vec![ChemElement::C.into(), ChemElement::N.into()],
            },
            stabilization_threshold: 0.375,
        },
        GraphAromaticityModel::Hmo {
            scope: GraphElementScope::AllowList(vec![ChemElement::C, ChemElement::N]),
            stabilization_threshold: 0.375,
        }
    )]
    #[case::clar(
        AromaticityModel::Clar {
            scope: ElementScope::AllowList {
                elements: vec![ChemElement::C.into()],
            },
            ring_limits: RingLimits::new(6, 14, true, 4, 1_500),
        },
        GraphAromaticityModel::Clar {
            scope: GraphElementScope::AllowList(vec![ChemElement::C]),
            ring_limits: GraphRingLimits {
                min_ring_size: 6,
                max_ring_size: 14,
                include_fused: true,
                max_fused_combination: 4,
                max_fused_search: 1_500,
            },
        }
    )]
    fn test_aromaticity_model_to_rust(
        #[case] model: AromaticityModel,
        #[case] expected: GraphAromaticityModel,
    ) {
        assert_eq!(model.to_rust(), expected);
    }
}
