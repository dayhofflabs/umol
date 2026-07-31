//! Python bindings for molecule-resolution configuration.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use pyo3::prelude::*;
use umol_graph::ops::resolve::{
    AromaticityInconsistencyPolicy as GraphAromaticityInconsistencyPolicy,
    AromaticityResolveConfig as GraphAromaticityResolveConfig,
    InconsistencyPolicy as GraphInconsistencyPolicy, ResolveConfig as GraphResolveConfig,
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
pub enum InconsistencyPolicy {
    Keep,
    Strip,
    Error,
}

impl InconsistencyPolicy {
    pub(crate) fn from_rust(policy: GraphInconsistencyPolicy) -> Self {
        match policy {
            GraphInconsistencyPolicy::Keep => Self::Keep,
            GraphInconsistencyPolicy::Strip => Self::Strip,
            GraphInconsistencyPolicy::Error => Self::Error,
        }
    }

    pub(crate) fn to_rust(self) -> GraphInconsistencyPolicy {
        match self {
            Self::Keep => GraphInconsistencyPolicy::Keep,
            Self::Strip => GraphInconsistencyPolicy::Strip,
            Self::Error => GraphInconsistencyPolicy::Error,
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
    #[pyo3(signature = (*, reset_stereo_constraints=false, inconsistency=InconsistencyPolicy::Error))]
    fn new(reset_stereo_constraints: bool, inconsistency: InconsistencyPolicy) -> Self {
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
    fn inconsistency(&self) -> InconsistencyPolicy {
        InconsistencyPolicy::from_rust(self.0.inconsistency)
    }

    fn __repr__(&self) -> String {
        format!(
            "StereoResolveConfig(reset_stereo_constraints={}, inconsistency=InconsistencyPolicy.{:?})",
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
    #[case::keep(GraphInconsistencyPolicy::Keep, InconsistencyPolicy::Keep)]
    #[case::strip(GraphInconsistencyPolicy::Strip, InconsistencyPolicy::Strip)]
    #[case::error(GraphInconsistencyPolicy::Error, InconsistencyPolicy::Error)]
    fn test_inconsistency_policy_from_rust(
        #[case] policy: GraphInconsistencyPolicy,
        #[case] expected: InconsistencyPolicy,
    ) {
        assert_eq!(InconsistencyPolicy::from_rust(policy), expected);
    }

    #[rstest]
    #[case::keep(InconsistencyPolicy::Keep, GraphInconsistencyPolicy::Keep)]
    #[case::strip(InconsistencyPolicy::Strip, GraphInconsistencyPolicy::Strip)]
    #[case::error(InconsistencyPolicy::Error, GraphInconsistencyPolicy::Error)]
    fn test_inconsistency_policy_to_rust(
        #[case] policy: InconsistencyPolicy,
        #[case] expected: GraphInconsistencyPolicy,
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
    #[case::default(false, InconsistencyPolicy::Error, GraphStereoResolveConfig::default())]
    #[case::keep(false, InconsistencyPolicy::Keep, GraphStereoResolveConfig {
        reset_stereo_constraints: false,
        inconsistency: GraphInconsistencyPolicy::Keep,
    })]
    #[case::strip(true, InconsistencyPolicy::Strip, GraphStereoResolveConfig {
        reset_stereo_constraints: true,
        inconsistency: GraphInconsistencyPolicy::Strip,
    })]
    fn test_stereo_resolve_config_new(
        #[case] reset_stereo_constraints: bool,
        #[case] inconsistency: InconsistencyPolicy,
        #[case] expected: GraphStereoResolveConfig,
    ) {
        assert_eq!(
            StereoResolveConfig::new(reset_stereo_constraints, inconsistency).0,
            expected
        );
    }

    #[rstest]
    #[case::default(
        StereoResolveConfig::new(false, InconsistencyPolicy::Error),
        "StereoResolveConfig(reset_stereo_constraints=False, inconsistency=InconsistencyPolicy.Error)"
    )]
    #[case::configured(
        StereoResolveConfig::new(true, InconsistencyPolicy::Strip),
        "StereoResolveConfig(reset_stereo_constraints=True, inconsistency=InconsistencyPolicy.Strip)"
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
        inconsistency: GraphInconsistencyPolicy::Strip,
    })]
    fn test_stereo_resolve_config_from_rust(#[case] config: GraphStereoResolveConfig) {
        assert_eq!(StereoResolveConfig::from_rust(config).0, config);
    }

    #[rstest]
    #[case::default(StereoResolveConfig::new(false, InconsistencyPolicy::Error))]
    #[case::configured(StereoResolveConfig::new(true, InconsistencyPolicy::Strip))]
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
        StereoResolveConfig::new(false, InconsistencyPolicy::Error),
        GraphResolveConfig::default()
    )]
    #[case::aromaticity(
        AromaticityResolveConfig::new(AromaticityConfig::default(), AromaticityInconsistencyPolicy::Keep, true),
        StereoResolveConfig::new(false, InconsistencyPolicy::Error),
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
        StereoResolveConfig::new(true, InconsistencyPolicy::Strip),
        GraphResolveConfig {
            aromaticity: GraphAromaticityResolveConfig::default(),
            stereo: GraphStereoResolveConfig {
                reset_stereo_constraints: true,
                inconsistency: GraphInconsistencyPolicy::Strip,
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
            StereoResolveConfig::new(true, InconsistencyPolicy::Strip),
        ),
        "ResolveConfig(aromaticity=AromaticityResolveConfig(perception=AromaticityConfig(ring_config=RingConfig(simple_cycle_algorithm=SimpleCycleEnumerationAlgorithm.ReadTarjan(), relevant_cycle_algorithm=RelevantCycleEnumerationAlgorithm.Vismara()), connected_components_algorithm=ConnectedComponentsAlgorithm.Bfs(), maximum_independent_set_algorithm=MaximumIndependentSetAlgorithm.BranchAndBound()), inconsistency=AromaticityInconsistencyPolicy.Keep, reset_aromatic_valence=True), stereo=StereoResolveConfig(reset_stereo_constraints=True, inconsistency=InconsistencyPolicy.Strip))",
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
            inconsistency: GraphInconsistencyPolicy::Strip,
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
        StereoResolveConfig::new(true, InconsistencyPolicy::Strip),
    ))]
    fn test_resolve_config_to_rust(#[case] config: ResolveConfig) {
        assert_eq!(config.to_rust(), config.0);
    }
}
