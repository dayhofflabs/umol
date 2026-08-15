//! Python bindings for aromaticity model and operation-configuration values.

use pyo3::prelude::*;
use umol_graph::ops::aromaticity::AromaticityConfig as GraphAromaticityConfig;
use umol_graph::ops::model::{
    AromaticityModel as GraphAromaticityModel, AromaticityRule as GraphAromaticityRule,
    AromaticityTieBreak as GraphAromaticityTieBreak, RingLimits as GraphRingLimits,
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

/// The aromaticity perception rule and its parameters.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq)]
pub enum AromaticityRule {
    /// Hückel-rule aromaticity over rings within the configured limits.
    #[pyo3(constructor = (*, ring_limits))]
    Hueckel { ring_limits: RingLimits },
    /// Hückel molecular-orbital aromaticity with a stabilization threshold.
    #[pyo3(constructor = (*, stabilization_threshold))]
    Hmo { stabilization_threshold: f64 },
    /// Clar aromaticity over rings within the configured limits.
    #[pyo3(constructor = (*, ring_limits))]
    Clar { ring_limits: RingLimits },
}

#[pymethods]
impl AromaticityRule {
    pub(crate) fn __repr__(&self) -> String {
        match self {
            Self::Hueckel { ring_limits } => format!(
                "AromaticityRule.Hueckel(ring_limits={})",
                ring_limits.__repr__(),
            ),
            Self::Hmo {
                stabilization_threshold,
            } => format!("AromaticityRule.Hmo(stabilization_threshold={stabilization_threshold})"),
            Self::Clar { ring_limits } => format!(
                "AromaticityRule.Clar(ring_limits={})",
                ring_limits.__repr__(),
            ),
        }
    }
}

impl AromaticityRule {
    pub(crate) fn from_rust(rule: &GraphAromaticityRule) -> Self {
        match rule {
            GraphAromaticityRule::Hueckel { ring_limits } => Self::Hueckel {
                ring_limits: RingLimits::from_rust(ring_limits),
            },
            GraphAromaticityRule::Hmo {
                stabilization_threshold,
            } => Self::Hmo {
                stabilization_threshold: *stabilization_threshold,
            },
            GraphAromaticityRule::Clar { ring_limits } => Self::Clar {
                ring_limits: RingLimits::from_rust(ring_limits),
            },
        }
    }

    pub(crate) fn to_rust(&self) -> GraphAromaticityRule {
        match self {
            Self::Hueckel { ring_limits } => GraphAromaticityRule::Hueckel {
                ring_limits: ring_limits.to_rust().clone(),
            },
            Self::Hmo {
                stabilization_threshold,
            } => GraphAromaticityRule::Hmo {
                stabilization_threshold: *stabilization_threshold,
            },
            Self::Clar { ring_limits } => GraphAromaticityRule::Clar {
                ring_limits: ring_limits.to_rust().clone(),
            },
        }
    }
}

/// Disposal policy for structurally distinct valid aromatic assignments.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AromaticityTieBreak {
    /// No structural preference: structurally distinct survivors stay plural.
    Strict,
    /// Maximal claimed-atom count; the perception decides what a system is,
    /// the policy only orders how much of the evidence is realized.
    MaxAtomCount,
}

#[pymethods]
impl AromaticityTieBreak {
    pub(crate) fn __repr__(&self) -> String {
        match self {
            Self::Strict => "AromaticityTieBreak.Strict".to_owned(),
            Self::MaxAtomCount => "AromaticityTieBreak.MaxAtomCount".to_owned(),
        }
    }
}

impl AromaticityTieBreak {
    pub(crate) fn from_rust(tie_break: GraphAromaticityTieBreak) -> Self {
        match tie_break {
            GraphAromaticityTieBreak::Strict => Self::Strict,
            GraphAromaticityTieBreak::MaxAtomCount => Self::MaxAtomCount,
        }
    }

    pub(crate) fn to_rust(self) -> GraphAromaticityTieBreak {
        match self {
            Self::Strict => GraphAromaticityTieBreak::Strict,
            Self::MaxAtomCount => GraphAromaticityTieBreak::MaxAtomCount,
        }
    }
}

/// Aromaticity model: the participating elements, the perception rule, and
/// how structurally distinct valid assignments are disposed of.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq)]
pub struct AromaticityModel {
    #[pyo3(get)]
    pub(crate) scope: ElementScope,
    #[pyo3(get)]
    pub(crate) rule: AromaticityRule,
    #[pyo3(get)]
    pub(crate) tie_break: AromaticityTieBreak,
}

#[pymethods]
impl AromaticityModel {
    #[new]
    #[pyo3(signature = (*, scope, rule, tie_break=AromaticityTieBreak::Strict))]
    fn new(scope: ElementScope, rule: AromaticityRule, tie_break: AromaticityTieBreak) -> Self {
        Self {
            scope,
            rule,
            tie_break,
        }
    }

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
        format!(
            "AromaticityModel(scope={}, rule={}, tie_break={})",
            self.scope.__repr__(),
            self.rule.__repr__(),
            self.tie_break.__repr__(),
        )
    }
}

impl AromaticityModel {
    pub(crate) fn from_rust(model: &GraphAromaticityModel) -> Self {
        Self {
            scope: ElementScope::from_rust(&model.scope),
            rule: AromaticityRule::from_rust(&model.rule),
            tie_break: AromaticityTieBreak::from_rust(model.tie_break),
        }
    }

    #[allow(
        dead_code,
        reason = "Python-to-Rust conversion API for ChemistryModel configuration"
    )]
    pub(crate) fn to_rust(&self) -> GraphAromaticityModel {
        GraphAromaticityModel {
            scope: self.scope.to_rust(),
            rule: self.rule.to_rust(),
            tie_break: self.tie_break.to_rust(),
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_chem::element::Element as ChemElement;
    use umol_graph::ops::aromaticity::AromaticityConfig as GraphAromaticityConfig;
    use umol_graph::ops::model::{
        AromaticityRule as GraphAromaticityRule, ElementScope as GraphElementScope,
    };

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
            AromaticityModel {
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
                rule: AromaticityRule::Hueckel {
                    ring_limits: RingLimits::new(3, 22, true, 6, 10_000)
                },
                tie_break: AromaticityTieBreak::Strict,
            }
        );
    }

    #[rstest]
    fn test_aromaticity_model_mdl() {
        assert_eq!(
            AromaticityModel::mdl(),
            AromaticityModel {
                scope: ElementScope::AllowList {
                    elements: vec![ChemElement::C.into(), ChemElement::N.into()],
                },
                rule: AromaticityRule::Hueckel {
                    ring_limits: RingLimits::new(6, 22, true, 6, 10_000)
                },
                tie_break: AromaticityTieBreak::Strict,
            }
        );
    }

    #[rstest]
    fn test_aromaticity_model_permissive() {
        assert_eq!(
            AromaticityModel::permissive(),
            AromaticityModel {
                scope: ElementScope::Any {},
                rule: AromaticityRule::Hueckel {
                    ring_limits: RingLimits::new(3, 22, true, 6, 10_000)
                },
                tie_break: AromaticityTieBreak::Strict,
            }
        );
    }

    #[rstest]
    #[case::hueckel_rule(
        AromaticityModel { scope: ElementScope::Any {}, rule: AromaticityRule::Hueckel { ring_limits: RingLimits::new(4, 18, false, 3, 2_000) }, tie_break: AromaticityTieBreak::Strict },
        "AromaticityModel(scope=ElementScope.Any(), rule=AromaticityRule.Hueckel(ring_limits=RingLimits(min_ring_size=4, max_ring_size=18, include_fused=False, max_fused_combination=3, max_fused_search=2000)), tie_break=AromaticityTieBreak.Strict)"
    )]
    #[case::hmo(
        AromaticityModel { scope: ElementScope::AllowList {
                elements: vec![ChemElement::C.into(), ChemElement::N.into()],
            }, rule: AromaticityRule::Hmo { stabilization_threshold: 0.375 }, tie_break: AromaticityTieBreak::MaxAtomCount },
        "AromaticityModel(scope=ElementScope.AllowList([Element('C'), Element('N')]), rule=AromaticityRule.Hmo(stabilization_threshold=0.375), tie_break=AromaticityTieBreak.MaxAtomCount)"
    )]
    #[case::clar(
        AromaticityModel { scope: ElementScope::AllowList {
                elements: vec![ChemElement::C.into()],
            }, rule: AromaticityRule::Clar { ring_limits: RingLimits::new(6, 14, true, 4, 1_500) }, tie_break: AromaticityTieBreak::Strict },
        "AromaticityModel(scope=ElementScope.AllowList([Element('C')]), rule=AromaticityRule.Clar(ring_limits=RingLimits(min_ring_size=6, max_ring_size=14, include_fused=True, max_fused_combination=4, max_fused_search=1500)), tie_break=AromaticityTieBreak.Strict)"
    )]
    fn test_aromaticity_model_repr(#[case] model: AromaticityModel, #[case] expected: &str) {
        assert_eq!(model.__repr__(), expected);
    }

    #[rstest]
    #[case::hueckel_rule(
        GraphAromaticityModel { scope: GraphElementScope::Any, rule: GraphAromaticityRule::Hueckel { ring_limits: GraphRingLimits {
                min_ring_size: 4,
                max_ring_size: 18,
                include_fused: false,
                max_fused_combination: 3,
                max_fused_search: 2_000,
            } }, tie_break: GraphAromaticityTieBreak::Strict },
        AromaticityModel { scope: ElementScope::Any {}, rule: AromaticityRule::Hueckel { ring_limits: RingLimits::new(4, 18, false, 3, 2_000) }, tie_break: AromaticityTieBreak::Strict }
    )]
    #[case::hmo(
        GraphAromaticityModel { scope: GraphElementScope::AllowList(vec![ChemElement::C, ChemElement::N]), rule: GraphAromaticityRule::Hmo { stabilization_threshold: 0.375 }, tie_break: GraphAromaticityTieBreak::MaxAtomCount },
        AromaticityModel { scope: ElementScope::AllowList {
                elements: vec![ChemElement::C.into(), ChemElement::N.into()],
            }, rule: AromaticityRule::Hmo { stabilization_threshold: 0.375 }, tie_break: AromaticityTieBreak::MaxAtomCount }
    )]
    #[case::clar(
        GraphAromaticityModel { scope: GraphElementScope::AllowList(vec![ChemElement::C]), rule: GraphAromaticityRule::Clar { ring_limits: GraphRingLimits {
                min_ring_size: 6,
                max_ring_size: 14,
                include_fused: true,
                max_fused_combination: 4,
                max_fused_search: 1_500,
            } }, tie_break: GraphAromaticityTieBreak::Strict },
        AromaticityModel { scope: ElementScope::AllowList {
                elements: vec![ChemElement::C.into()],
            }, rule: AromaticityRule::Clar { ring_limits: RingLimits::new(6, 14, true, 4, 1_500) }, tie_break: AromaticityTieBreak::Strict }
    )]
    fn test_aromaticity_model_from_rust(
        #[case] model: GraphAromaticityModel,
        #[case] expected: AromaticityModel,
    ) {
        assert_eq!(AromaticityModel::from_rust(&model), expected);
    }

    #[rstest]
    #[case::hueckel_rule(
        AromaticityModel { scope: ElementScope::Any {}, rule: AromaticityRule::Hueckel { ring_limits: RingLimits::new(4, 18, false, 3, 2_000) }, tie_break: AromaticityTieBreak::Strict },
        GraphAromaticityModel { scope: GraphElementScope::Any, rule: GraphAromaticityRule::Hueckel { ring_limits: GraphRingLimits {
                min_ring_size: 4,
                max_ring_size: 18,
                include_fused: false,
                max_fused_combination: 3,
                max_fused_search: 2_000,
            } }, tie_break: GraphAromaticityTieBreak::Strict }
    )]
    #[case::hmo(
        AromaticityModel { scope: ElementScope::AllowList {
                elements: vec![ChemElement::C.into(), ChemElement::N.into()],
            }, rule: AromaticityRule::Hmo { stabilization_threshold: 0.375 }, tie_break: AromaticityTieBreak::MaxAtomCount },
        GraphAromaticityModel { scope: GraphElementScope::AllowList(vec![ChemElement::C, ChemElement::N]), rule: GraphAromaticityRule::Hmo { stabilization_threshold: 0.375 }, tie_break: GraphAromaticityTieBreak::MaxAtomCount }
    )]
    #[case::clar(
        AromaticityModel { scope: ElementScope::AllowList {
                elements: vec![ChemElement::C.into()],
            }, rule: AromaticityRule::Clar { ring_limits: RingLimits::new(6, 14, true, 4, 1_500) }, tie_break: AromaticityTieBreak::Strict },
        GraphAromaticityModel { scope: GraphElementScope::AllowList(vec![ChemElement::C]), rule: GraphAromaticityRule::Clar { ring_limits: GraphRingLimits {
                min_ring_size: 6,
                max_ring_size: 14,
                include_fused: true,
                max_fused_combination: 4,
                max_fused_search: 1_500,
            } }, tie_break: GraphAromaticityTieBreak::Strict }
    )]
    fn test_aromaticity_model_to_rust(
        #[case] model: AromaticityModel,
        #[case] expected: GraphAromaticityModel,
    ) {
        assert_eq!(model.to_rust(), expected);
    }
}
