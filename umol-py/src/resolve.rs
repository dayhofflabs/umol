//! Python bindings for molecule-resolution configuration.

use pyo3::prelude::*;
use umol_graph::ops::resolve::{
    AromaticityResolveConfig as GraphAromaticityResolveConfig, ResolveConfig as GraphResolveConfig,
    StereoResolveConfig as GraphStereoResolveConfig,
};

/// Operational policy for aromaticity resolution.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AromaticityResolveConfig(GraphAromaticityResolveConfig);

#[pymethods]
impl AromaticityResolveConfig {
    #[new]
    #[pyo3(signature = (*, delocalize_charge=true, reset_aromatic_valence=false))]
    fn new(delocalize_charge: bool, reset_aromatic_valence: bool) -> Self {
        Self(GraphAromaticityResolveConfig {
            perception: Default::default(),
            delocalize_charge,
            reset_aromatic_valence,
        })
    }

    #[getter]
    fn delocalize_charge(&self) -> bool {
        self.0.delocalize_charge
    }

    #[getter]
    fn reset_aromatic_valence(&self) -> bool {
        self.0.reset_aromatic_valence
    }

    fn __repr__(&self) -> String {
        format!(
            "AromaticityResolveConfig(delocalize_charge={}, reset_aromatic_valence={})",
            if self.0.delocalize_charge {
                "True"
            } else {
                "False"
            },
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
    #[pyo3(signature = (*, reset_stereo_constraints=false))]
    fn new(reset_stereo_constraints: bool) -> Self {
        Self(GraphStereoResolveConfig {
            reset_stereo_constraints,
        })
    }

    #[getter]
    fn reset_stereo_constraints(&self) -> bool {
        self.0.reset_stereo_constraints
    }

    fn __repr__(&self) -> String {
        format!(
            "StereoResolveConfig(reset_stereo_constraints={})",
            if self.0.reset_stereo_constraints {
                "True"
            } else {
                "False"
            },
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
    #[case::default(true, false, GraphAromaticityResolveConfig::default())]
    #[case::retain_charge(false, false, GraphAromaticityResolveConfig {
        perception: Default::default(),
        delocalize_charge: false,
        reset_aromatic_valence: false,
    })]
    #[case::reset_valence(true, true, GraphAromaticityResolveConfig {
        perception: Default::default(),
        delocalize_charge: true,
        reset_aromatic_valence: true,
    })]
    #[case::both(false, true, GraphAromaticityResolveConfig {
        perception: Default::default(),
        delocalize_charge: false,
        reset_aromatic_valence: true,
    })]
    fn test_aromaticity_resolve_config_new(
        #[case] delocalize_charge: bool,
        #[case] reset_aromatic_valence: bool,
        #[case] expected: GraphAromaticityResolveConfig,
    ) {
        assert_eq!(
            AromaticityResolveConfig::new(delocalize_charge, reset_aromatic_valence).0,
            expected
        );
    }

    #[rstest]
    #[case::default(
        AromaticityResolveConfig::new(true, false),
        "AromaticityResolveConfig(delocalize_charge=True, reset_aromatic_valence=False)"
    )]
    #[case::nondefault(
        AromaticityResolveConfig::new(false, true),
        "AromaticityResolveConfig(delocalize_charge=False, reset_aromatic_valence=True)"
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
        delocalize_charge: false,
        reset_aromatic_valence: true,
    })]
    fn test_aromaticity_resolve_config_from_rust(#[case] config: GraphAromaticityResolveConfig) {
        assert_eq!(AromaticityResolveConfig::from_rust(config).0, config);
    }

    #[rstest]
    #[case::default(AromaticityResolveConfig::new(true, false))]
    #[case::nondefault(AromaticityResolveConfig::new(false, true))]
    fn test_aromaticity_resolve_config_to_rust(#[case] config: AromaticityResolveConfig) {
        assert_eq!(config.to_rust(), config.0);
    }

    #[rstest]
    #[case::default(false, GraphStereoResolveConfig::default())]
    #[case::reset(true, GraphStereoResolveConfig {
        reset_stereo_constraints: true,
    })]
    fn test_stereo_resolve_config_new(
        #[case] reset_stereo_constraints: bool,
        #[case] expected: GraphStereoResolveConfig,
    ) {
        assert_eq!(
            StereoResolveConfig::new(reset_stereo_constraints).0,
            expected
        );
    }

    #[rstest]
    #[case::default(
        StereoResolveConfig::new(false),
        "StereoResolveConfig(reset_stereo_constraints=False)"
    )]
    #[case::reset(
        StereoResolveConfig::new(true),
        "StereoResolveConfig(reset_stereo_constraints=True)"
    )]
    fn test_stereo_resolve_config_repr(
        #[case] config: StereoResolveConfig,
        #[case] expected: &str,
    ) {
        assert_eq!(config.__repr__(), expected);
    }

    #[rstest]
    #[case::default(GraphStereoResolveConfig::default())]
    #[case::reset(GraphStereoResolveConfig {
        reset_stereo_constraints: true,
    })]
    fn test_stereo_resolve_config_from_rust(#[case] config: GraphStereoResolveConfig) {
        assert_eq!(StereoResolveConfig::from_rust(config).0, config);
    }

    #[rstest]
    #[case::default(StereoResolveConfig::new(false))]
    #[case::reset(StereoResolveConfig::new(true))]
    fn test_stereo_resolve_config_to_rust(#[case] config: StereoResolveConfig) {
        assert_eq!(config.to_rust(), config.0);
    }

    #[rstest]
    #[case::default(
        AromaticityResolveConfig::new(true, false),
        StereoResolveConfig::new(false),
        GraphResolveConfig::default()
    )]
    #[case::aromaticity(
        AromaticityResolveConfig::new(false, true),
        StereoResolveConfig::new(false),
        GraphResolveConfig {
            aromaticity: GraphAromaticityResolveConfig {
                perception: Default::default(),
                delocalize_charge: false,
                reset_aromatic_valence: true,
            },
            stereo: GraphStereoResolveConfig::default(),
        },
    )]
    #[case::stereo(
        AromaticityResolveConfig::new(true, false),
        StereoResolveConfig::new(true),
        GraphResolveConfig {
            aromaticity: GraphAromaticityResolveConfig::default(),
            stereo: GraphStereoResolveConfig {
                reset_stereo_constraints: true,
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
            AromaticityResolveConfig::new(false, true),
            StereoResolveConfig::new(true),
        ),
        "ResolveConfig(aromaticity=AromaticityResolveConfig(delocalize_charge=False, reset_aromatic_valence=True), stereo=StereoResolveConfig(reset_stereo_constraints=True))",
    )]
    fn test_resolve_config_repr(#[case] config: ResolveConfig, #[case] expected: &str) {
        assert_eq!(config.__repr__(), expected);
    }

    #[rstest]
    #[case::default(GraphResolveConfig::default())]
    #[case::configured(GraphResolveConfig {
        aromaticity: GraphAromaticityResolveConfig {
            perception: Default::default(),
            delocalize_charge: false,
            reset_aromatic_valence: true,
        },
        stereo: GraphStereoResolveConfig {
            reset_stereo_constraints: true,
        },
    })]
    fn test_resolve_config_from_rust(#[case] config: GraphResolveConfig) {
        assert_eq!(ResolveConfig::from_rust(config).0, config);
    }

    #[rstest]
    #[case::default(ResolveConfig::default())]
    #[case::configured(ResolveConfig::new(
        AromaticityResolveConfig::new(false, true),
        StereoResolveConfig::new(true),
    ))]
    fn test_resolve_config_to_rust(#[case] config: ResolveConfig) {
        assert_eq!(config.to_rust(), config.0);
    }
}
