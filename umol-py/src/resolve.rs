//! Python bindings for molecule-resolution configuration.

use pyo3::prelude::*;
use umol_graph::ops::resolve::AromaticityResolveConfig as GraphAromaticityResolveConfig;

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

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::default(true, false, GraphAromaticityResolveConfig::default())]
    #[case::retain_charge(false, false, GraphAromaticityResolveConfig {
        delocalize_charge: false,
        reset_aromatic_valence: false,
    })]
    #[case::reset_valence(true, true, GraphAromaticityResolveConfig {
        delocalize_charge: true,
        reset_aromatic_valence: true,
    })]
    #[case::both(false, true, GraphAromaticityResolveConfig {
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
}
