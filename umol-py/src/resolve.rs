//! Python bindings for molecule-resolution configuration.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use pyo3::prelude::*;
use umol_graph::ops::resolve::{
    AromaticityInconsistencyPolicy as GraphAromaticityInconsistencyPolicy,
    AromaticityResolveConfig as GraphAromaticityResolveConfig, ResolveConfig as GraphResolveConfig,
    StereoInconsistencyPolicy as GraphStereoInconsistencyPolicy,
    StereoResolveConfig as GraphStereoResolveConfig,
};

use crate::model::aromaticity::AromaticityConfig;

/// Policy for aromaticity assertions that disagree with perception.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum AromaticityInconsistencyPolicy {
    Keep,
    Error,
}

impl AromaticityInconsistencyPolicy {
    pub(crate) fn from_rust(policy: GraphAromaticityInconsistencyPolicy) -> Self {
        match policy {
            GraphAromaticityInconsistencyPolicy::Keep => Self::Keep,
            GraphAromaticityInconsistencyPolicy::Error => Self::Error,
        }
    }

    pub(crate) fn to_rust(self) -> GraphAromaticityInconsistencyPolicy {
        match self {
            Self::Keep => GraphAromaticityInconsistencyPolicy::Keep,
            Self::Error => GraphAromaticityInconsistencyPolicy::Error,
        }
    }
}

/// Policy for stereo assertions that cannot be fully realized.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum StereoInconsistencyPolicy {
    Keep,
    Strip,
    Error,
}

impl StereoInconsistencyPolicy {
    pub(crate) fn from_rust(policy: GraphStereoInconsistencyPolicy) -> Self {
        match policy {
            GraphStereoInconsistencyPolicy::Keep => Self::Keep,
            GraphStereoInconsistencyPolicy::Strip => Self::Strip,
            GraphStereoInconsistencyPolicy::Error => Self::Error,
        }
    }

    pub(crate) fn to_rust(self) -> GraphStereoInconsistencyPolicy {
        match self {
            Self::Keep => GraphStereoInconsistencyPolicy::Keep,
            Self::Strip => GraphStereoInconsistencyPolicy::Strip,
            Self::Error => GraphStereoInconsistencyPolicy::Error,
        }
    }
}

/// Operational policy for aromaticity resolution.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AromaticityResolveConfig(GraphAromaticityResolveConfig);

#[pymethods]
impl AromaticityResolveConfig {
    #[new]
    #[pyo3(signature = (*, perception=AromaticityConfig::default(), inconsistency=AromaticityInconsistencyPolicy::Error, reset_aromatic_valence=false))]
    fn new(
        perception: AromaticityConfig,
        inconsistency: AromaticityInconsistencyPolicy,
        reset_aromatic_valence: bool,
    ) -> Self {
        Self(GraphAromaticityResolveConfig {
            perception: perception.to_rust(),
            inconsistency: inconsistency.to_rust(),
            reset_aromatic_valence,
        })
    }

    #[getter]
    fn perception(&self) -> AromaticityConfig {
        AromaticityConfig::from_rust(self.0.perception)
    }

    #[getter]
    fn inconsistency(&self) -> AromaticityInconsistencyPolicy {
        AromaticityInconsistencyPolicy::from_rust(self.0.inconsistency)
    }

    #[getter]
    fn reset_aromatic_valence(&self) -> bool {
        self.0.reset_aromatic_valence
    }

    fn __repr__(&self) -> String {
        format!(
            "AromaticityResolveConfig(perception={}, inconsistency=AromaticityInconsistencyPolicy.{:?}, reset_aromatic_valence={})",
            self.perception().__repr__(),
            self.inconsistency(),
            if self.0.reset_aromatic_valence {
                "True"
            } else {
                "False"
            },
        )
    }
}

impl AromaticityResolveConfig {
    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion API for ResolveConfig composition"
    )]
    pub(crate) fn from_rust(config: GraphAromaticityResolveConfig) -> Self {
        Self(config)
    }

    #[allow(
        dead_code,
        reason = "Python-to-Rust conversion API for ResolveConfig composition"
    )]
    pub(crate) fn to_rust(self) -> GraphAromaticityResolveConfig {
        self.0
    }
}

/// Operational policy for stereo resolution.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct StereoResolveConfig(GraphStereoResolveConfig);

#[pymethods]
impl StereoResolveConfig {
    #[new]
    #[pyo3(signature = (*, reset_stereo_constraints=false, inconsistency=StereoInconsistencyPolicy::Error))]
    fn new(reset_stereo_constraints: bool, inconsistency: StereoInconsistencyPolicy) -> Self {
        Self(GraphStereoResolveConfig {
            reset_stereo_constraints,
            inconsistency: inconsistency.to_rust(),
        })
    }

    #[getter]
    fn reset_stereo_constraints(&self) -> bool {
        self.0.reset_stereo_constraints
    }

    #[getter]
    fn inconsistency(&self) -> StereoInconsistencyPolicy {
        StereoInconsistencyPolicy::from_rust(self.0.inconsistency)
    }

    fn __repr__(&self) -> String {
        format!(
            "StereoResolveConfig(reset_stereo_constraints={}, inconsistency=StereoInconsistencyPolicy.{:?})",
            if self.0.reset_stereo_constraints {
                "True"
            } else {
                "False"
            },
            self.inconsistency(),
        )
    }
}

impl StereoResolveConfig {
    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion API for ResolveConfig composition"
    )]
    pub(crate) fn from_rust(config: GraphStereoResolveConfig) -> Self {
        Self(config)
    }

    #[allow(
        dead_code,
        reason = "Python-to-Rust conversion API for ResolveConfig composition"
    )]
    pub(crate) fn to_rust(self) -> GraphStereoResolveConfig {
        self.0
    }
}

/// Operational policy for molecule resolution.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ResolveConfig(GraphResolveConfig);

#[pymethods]
impl ResolveConfig {
    #[new]
    #[pyo3(signature = (*, aromaticity, stereo))]
    fn new(aromaticity: AromaticityResolveConfig, stereo: StereoResolveConfig) -> Self {
        Self(GraphResolveConfig {
            aromaticity: aromaticity.to_rust(),
            stereo: stereo.to_rust(),
        })
    }

    #[staticmethod]
    fn default() -> Self {
        Self::from_rust(GraphResolveConfig::default())
    }

    #[getter]
    fn aromaticity(&self) -> AromaticityResolveConfig {
        AromaticityResolveConfig::from_rust(self.0.aromaticity)
    }

    #[getter]
    fn stereo(&self) -> StereoResolveConfig {
        StereoResolveConfig::from_rust(self.0.stereo)
    }

    fn __repr__(&self) -> String {
        if self.0 == GraphResolveConfig::default() {
            return "ResolveConfig.default()".to_owned();
        }
        format!(
            "ResolveConfig(aromaticity={}, stereo={})",
            self.aromaticity().__repr__(),
            self.stereo().__repr__(),
        )
    }
}

impl ResolveConfig {
    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion API for configured molecule ingestion"
    )]
    pub(crate) fn from_rust(config: GraphResolveConfig) -> Self {
        Self(config)
    }

    pub(crate) fn to_rust(self) -> GraphResolveConfig {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::keep(
        GraphAromaticityInconsistencyPolicy::Keep,
        AromaticityInconsistencyPolicy::Keep
    )]
    #[case::error(
        GraphAromaticityInconsistencyPolicy::Error,
        AromaticityInconsistencyPolicy::Error
    )]
    fn test_aromaticity_inconsistency_policy_from_rust(
        #[case] policy: GraphAromaticityInconsistencyPolicy,
        #[case] expected: AromaticityInconsistencyPolicy,
    ) {
        assert_eq!(AromaticityInconsistencyPolicy::from_rust(policy), expected);
    }

    #[rstest]
    #[case::keep(
        AromaticityInconsistencyPolicy::Keep,
        GraphAromaticityInconsistencyPolicy::Keep
    )]
    #[case::error(
        AromaticityInconsistencyPolicy::Error,
        GraphAromaticityInconsistencyPolicy::Error
    )]
    fn test_aromaticity_inconsistency_policy_to_rust(
        #[case] policy: AromaticityInconsistencyPolicy,
        #[case] expected: GraphAromaticityInconsistencyPolicy,
    ) {
        assert_eq!(policy.to_rust(), expected);
    }

    #[rstest]
    #[case::keep(GraphStereoInconsistencyPolicy::Keep, StereoInconsistencyPolicy::Keep)]
    #[case::strip(
        GraphStereoInconsistencyPolicy::Strip,
        StereoInconsistencyPolicy::Strip
    )]
    #[case::error(
        GraphStereoInconsistencyPolicy::Error,
        StereoInconsistencyPolicy::Error
    )]
    fn test_stereo_inconsistency_policy_from_rust(
        #[case] policy: GraphStereoInconsistencyPolicy,
        #[case] expected: StereoInconsistencyPolicy,
    ) {
        assert_eq!(StereoInconsistencyPolicy::from_rust(policy), expected);
    }

    #[rstest]
    #[case::keep(StereoInconsistencyPolicy::Keep, GraphStereoInconsistencyPolicy::Keep)]
    #[case::strip(
        StereoInconsistencyPolicy::Strip,
        GraphStereoInconsistencyPolicy::Strip
    )]
    #[case::error(
        StereoInconsistencyPolicy::Error,
        GraphStereoInconsistencyPolicy::Error
    )]
    fn test_stereo_inconsistency_policy_to_rust(
        #[case] policy: StereoInconsistencyPolicy,
        #[case] expected: GraphStereoInconsistencyPolicy,
    ) {
        assert_eq!(policy.to_rust(), expected);
    }

    #[rstest]
    #[case::default(
        AromaticityConfig::default(),
        AromaticityInconsistencyPolicy::Error,
        false,
        GraphAromaticityResolveConfig::default()
    )]
    #[case::keep(AromaticityConfig::default(), AromaticityInconsistencyPolicy::Keep, false, GraphAromaticityResolveConfig {
        perception: Default::default(),
        inconsistency: GraphAromaticityInconsistencyPolicy::Keep,
        reset_aromatic_valence: false,
    })]
    #[case::reset_valence(AromaticityConfig::default(), AromaticityInconsistencyPolicy::Error, true, GraphAromaticityResolveConfig {
        perception: Default::default(),
        inconsistency: GraphAromaticityInconsistencyPolicy::Error,
        reset_aromatic_valence: true,
    })]
    #[case::both(AromaticityConfig::default(), AromaticityInconsistencyPolicy::Keep, true, GraphAromaticityResolveConfig {
        perception: Default::default(),
        inconsistency: GraphAromaticityInconsistencyPolicy::Keep,
        reset_aromatic_valence: true,
    })]
    fn test_aromaticity_resolve_config_new(
        #[case] perception: AromaticityConfig,
        #[case] inconsistency: AromaticityInconsistencyPolicy,
        #[case] reset_aromatic_valence: bool,
        #[case] expected: GraphAromaticityResolveConfig,
    ) {
        assert_eq!(
            AromaticityResolveConfig::new(perception, inconsistency, reset_aromatic_valence).0,
            expected
        );
    }

    #[rstest]
    #[case::default(
        AromaticityResolveConfig::new(AromaticityConfig::default(), AromaticityInconsistencyPolicy::Error, false),
        "AromaticityResolveConfig(perception=AromaticityConfig(ring_config=RingConfig(simple_cycle_algorithm=SimpleCycleEnumerationAlgorithm.ReadTarjan(), relevant_cycle_algorithm=RelevantCycleEnumerationAlgorithm.Vismara()), connected_components_algorithm=ConnectedComponentsAlgorithm.Bfs(), maximum_independent_set_algorithm=MaximumIndependentSetAlgorithm.BranchAndBound()), inconsistency=AromaticityInconsistencyPolicy.Error, reset_aromatic_valence=False)"
    )]
    #[case::nondefault(
        AromaticityResolveConfig::new(AromaticityConfig::default(), AromaticityInconsistencyPolicy::Keep, true),
        "AromaticityResolveConfig(perception=AromaticityConfig(ring_config=RingConfig(simple_cycle_algorithm=SimpleCycleEnumerationAlgorithm.ReadTarjan(), relevant_cycle_algorithm=RelevantCycleEnumerationAlgorithm.Vismara()), connected_components_algorithm=ConnectedComponentsAlgorithm.Bfs(), maximum_independent_set_algorithm=MaximumIndependentSetAlgorithm.BranchAndBound()), inconsistency=AromaticityInconsistencyPolicy.Keep, reset_aromatic_valence=True)"
    )]
    fn test_aromaticity_resolve_config_repr(
        #[case] config: AromaticityResolveConfig,
        #[case] expected: &str,
    ) {
        assert_eq!(config.__repr__(), expected);
    }

    #[rstest]
    #[case::default(GraphAromaticityResolveConfig::default())]
    #[case::nondefault(GraphAromaticityResolveConfig {
        perception: Default::default(),
        inconsistency: GraphAromaticityInconsistencyPolicy::Keep,
        reset_aromatic_valence: true,
    })]
    fn test_aromaticity_resolve_config_from_rust(#[case] config: GraphAromaticityResolveConfig) {
        assert_eq!(AromaticityResolveConfig::from_rust(config).0, config);
    }

    #[rstest]
    #[case::default(AromaticityResolveConfig::new(
        AromaticityConfig::default(),
        AromaticityInconsistencyPolicy::Error,
        false
    ))]
    #[case::nondefault(AromaticityResolveConfig::new(
        AromaticityConfig::default(),
        AromaticityInconsistencyPolicy::Keep,
        true
    ))]
    fn test_aromaticity_resolve_config_to_rust(#[case] config: AromaticityResolveConfig) {
        assert_eq!(config.to_rust(), config.0);
    }

    #[rstest]
    #[case::default(
        false,
        StereoInconsistencyPolicy::Error,
        GraphStereoResolveConfig::default()
    )]
    #[case::keep(false, StereoInconsistencyPolicy::Keep, GraphStereoResolveConfig {
        reset_stereo_constraints: false,
        inconsistency: GraphStereoInconsistencyPolicy::Keep,
    })]
    #[case::strip(true, StereoInconsistencyPolicy::Strip, GraphStereoResolveConfig {
        reset_stereo_constraints: true,
        inconsistency: GraphStereoInconsistencyPolicy::Strip,
    })]
    fn test_stereo_resolve_config_new(
        #[case] reset_stereo_constraints: bool,
        #[case] inconsistency: StereoInconsistencyPolicy,
        #[case] expected: GraphStereoResolveConfig,
    ) {
        assert_eq!(
            StereoResolveConfig::new(reset_stereo_constraints, inconsistency).0,
            expected
        );
    }

    #[rstest]
    #[case::default(
        StereoResolveConfig::new(false, StereoInconsistencyPolicy::Error),
        "StereoResolveConfig(reset_stereo_constraints=False, inconsistency=StereoInconsistencyPolicy.Error)"
    )]
    #[case::configured(
        StereoResolveConfig::new(true, StereoInconsistencyPolicy::Strip),
        "StereoResolveConfig(reset_stereo_constraints=True, inconsistency=StereoInconsistencyPolicy.Strip)"
    )]
    fn test_stereo_resolve_config_repr(
        #[case] config: StereoResolveConfig,
        #[case] expected: &str,
    ) {
        assert_eq!(config.__repr__(), expected);
    }

    #[rstest]
    #[case::default(GraphStereoResolveConfig::default())]
    #[case::configured(GraphStereoResolveConfig {
        reset_stereo_constraints: true,
        inconsistency: GraphStereoInconsistencyPolicy::Strip,
    })]
    fn test_stereo_resolve_config_from_rust(#[case] config: GraphStereoResolveConfig) {
        assert_eq!(StereoResolveConfig::from_rust(config).0, config);
    }

    #[rstest]
    #[case::default(StereoResolveConfig::new(false, StereoInconsistencyPolicy::Error))]
    #[case::configured(StereoResolveConfig::new(true, StereoInconsistencyPolicy::Strip))]
    fn test_stereo_resolve_config_to_rust(#[case] config: StereoResolveConfig) {
        assert_eq!(config.to_rust(), config.0);
    }

    #[rstest]
    #[case::default(
        AromaticityResolveConfig::new(
            AromaticityConfig::default(),
            AromaticityInconsistencyPolicy::Error,
            false
        ),
        StereoResolveConfig::new(false, StereoInconsistencyPolicy::Error),
        GraphResolveConfig::default()
    )]
    #[case::aromaticity(
        AromaticityResolveConfig::new(AromaticityConfig::default(), AromaticityInconsistencyPolicy::Keep, true),
        StereoResolveConfig::new(false, StereoInconsistencyPolicy::Error),
        GraphResolveConfig {
            aromaticity: GraphAromaticityResolveConfig {
                perception: Default::default(),
                inconsistency: GraphAromaticityInconsistencyPolicy::Keep,
                reset_aromatic_valence: true,
            },
            stereo: GraphStereoResolveConfig::default(),
        },
    )]
    #[case::stereo(
        AromaticityResolveConfig::new(AromaticityConfig::default(), AromaticityInconsistencyPolicy::Error, false),
        StereoResolveConfig::new(true, StereoInconsistencyPolicy::Strip),
        GraphResolveConfig {
            aromaticity: GraphAromaticityResolveConfig::default(),
            stereo: GraphStereoResolveConfig {
                reset_stereo_constraints: true,
                inconsistency: GraphStereoInconsistencyPolicy::Strip,
            },
        },
    )]
    fn test_resolve_config_new(
        #[case] aromaticity: AromaticityResolveConfig,
        #[case] stereo: StereoResolveConfig,
        #[case] expected: GraphResolveConfig,
    ) {
        assert_eq!(ResolveConfig::new(aromaticity, stereo).0, expected);
    }

    #[rstest]
    fn test_resolve_config_default() {
        assert_eq!(
            ResolveConfig::default(),
            ResolveConfig(GraphResolveConfig {
                aromaticity: GraphAromaticityResolveConfig::default(),
                stereo: GraphStereoResolveConfig::default(),
            })
        );
    }

    #[rstest]
    #[case::default(ResolveConfig::default(), "ResolveConfig.default()")]
    #[case::configured(
        ResolveConfig::new(
            AromaticityResolveConfig::new(AromaticityConfig::default(), AromaticityInconsistencyPolicy::Keep, true),
            StereoResolveConfig::new(true, StereoInconsistencyPolicy::Strip),
        ),
        "ResolveConfig(aromaticity=AromaticityResolveConfig(perception=AromaticityConfig(ring_config=RingConfig(simple_cycle_algorithm=SimpleCycleEnumerationAlgorithm.ReadTarjan(), relevant_cycle_algorithm=RelevantCycleEnumerationAlgorithm.Vismara()), connected_components_algorithm=ConnectedComponentsAlgorithm.Bfs(), maximum_independent_set_algorithm=MaximumIndependentSetAlgorithm.BranchAndBound()), inconsistency=AromaticityInconsistencyPolicy.Keep, reset_aromatic_valence=True), stereo=StereoResolveConfig(reset_stereo_constraints=True, inconsistency=StereoInconsistencyPolicy.Strip))",
    )]
    fn test_resolve_config_repr(#[case] config: ResolveConfig, #[case] expected: &str) {
        assert_eq!(config.__repr__(), expected);
    }

    #[rstest]
    #[case::default(GraphResolveConfig::default())]
    #[case::configured(GraphResolveConfig {
        aromaticity: GraphAromaticityResolveConfig {
            perception: Default::default(),
            inconsistency: GraphAromaticityInconsistencyPolicy::Keep,
            reset_aromatic_valence: true,
        },
        stereo: GraphStereoResolveConfig {
            reset_stereo_constraints: true,
            inconsistency: GraphStereoInconsistencyPolicy::Strip,
        },
    })]
    fn test_resolve_config_from_rust(#[case] config: GraphResolveConfig) {
        assert_eq!(ResolveConfig::from_rust(config).0, config);
    }

    #[rstest]
    #[case::default(ResolveConfig::default())]
    #[case::configured(ResolveConfig::new(
        AromaticityResolveConfig::new(
            AromaticityConfig::default(),
            AromaticityInconsistencyPolicy::Keep,
            true
        ),
        StereoResolveConfig::new(true, StereoInconsistencyPolicy::Strip),
    ))]
    fn test_resolve_config_to_rust(#[case] config: ResolveConfig) {
        assert_eq!(config.to_rust(), config.0);
    }
}
