//! Transformation configuration values and the owned stereo iterator.

use pyo3::prelude::*;
use umol_graph::ops::transform::{
    AddStereoConfig as GraphAddStereoConfig, EnumerateStereoConfig as GraphEnumerateStereoConfig,
    EnumerateStereoIter as GraphEnumerateStereoIter, KekulizeConfig as GraphKekulizeConfig,
    RemoveStereoConfig as GraphRemoveStereoConfig, RemoveStereoMode as GraphRemoveStereoMode,
};

use crate::algorithm::{AutomorphismAlgorithm, MaximumMatchingAlgorithm};
use crate::molecule::Molecule;
use crate::ring::RingConfig;
use crate::stereo::{StereoConfigurationForm, StereoSelection};

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

/// Which stored configuration states stereo removal accepts.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RemoveStereoMode {
    All,
    Undetermined,
}

#[pymethods]
impl RemoveStereoMode {
    /// Test the stored configuration without classification or resolution.
    fn matches(&self, py: Python<'_>, configuration: &StereoConfigurationForm) -> bool {
        self.to_rust().matches(&configuration.to_rust(py))
    }

    fn __repr__(&self) -> String {
        format!("RemoveStereoMode.{self:?}")
    }
}

impl RemoveStereoMode {
    pub(crate) fn from_rust(mode: GraphRemoveStereoMode) -> Self {
        match mode {
            GraphRemoveStereoMode::All => Self::All,
            GraphRemoveStereoMode::Undetermined => Self::Undetermined,
        }
    }
    pub(crate) fn to_rust(self) -> GraphRemoveStereoMode {
        match self {
            Self::All => GraphRemoveStereoMode::All,
            Self::Undetermined => GraphRemoveStereoMode::Undetermined,
        }
    }
}

/// Immutable operational configuration corresponding to Rust AddStereoConfig.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AddStereoConfig(GraphAddStereoConfig);

#[pymethods]
impl AddStereoConfig {
    #[new]
    #[pyo3(signature = (*, selection=None, automorphism_algorithm=AutomorphismAlgorithm::Nauty(), ring_config=None))]
    fn new(
        selection: Option<&StereoSelection>,
        automorphism_algorithm: AutomorphismAlgorithm,
        ring_config: Option<RingConfig>,
    ) -> Self {
        Self::from_rust(GraphAddStereoConfig {
            selection: selection.map_or_else(
                || GraphAddStereoConfig::default().selection,
                |value| value.to_rust().clone(),
            ),
            automorphism_algorithm: automorphism_algorithm.to_rust(),
            ring_config: ring_config.unwrap_or_default().to_rust(),
        })
    }

    #[getter]
    fn selection(&self) -> StereoSelection {
        StereoSelection::from_rust(self.to_rust().selection.clone())
    }
    #[getter]
    fn automorphism_algorithm(&self) -> AutomorphismAlgorithm {
        AutomorphismAlgorithm::from_rust(self.to_rust().automorphism_algorithm)
    }

    #[getter]
    fn ring_config(&self) -> RingConfig {
        RingConfig::from_rust(self.to_rust().ring_config)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "AddStereoConfig(selection={}, automorphism_algorithm={}, ring_config={})",
            self.selection().__repr__(py)?,
            self.automorphism_algorithm().repr(),
            self.ring_config().__repr__()
        ))
    }
}

impl AddStereoConfig {
    pub(crate) fn from_rust(config: GraphAddStereoConfig) -> Self {
        Self(config)
    }
    pub(crate) fn to_rust(&self) -> &GraphAddStereoConfig {
        &self.0
    }
}

/// Immutable operational configuration corresponding to Rust RemoveStereoConfig.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RemoveStereoConfig(GraphRemoveStereoConfig);

#[pymethods]
impl RemoveStereoConfig {
    #[new]
    #[pyo3(signature = (*, selection=None, mode=RemoveStereoMode::All, automorphism_algorithm=AutomorphismAlgorithm::Nauty()))]
    fn new(
        selection: Option<&StereoSelection>,
        mode: RemoveStereoMode,
        automorphism_algorithm: AutomorphismAlgorithm,
    ) -> Self {
        Self::from_rust(GraphRemoveStereoConfig {
            selection: selection.map_or_else(
                || GraphRemoveStereoConfig::default().selection,
                |value| value.to_rust().clone(),
            ),
            mode: mode.to_rust(),
            automorphism_algorithm: automorphism_algorithm.to_rust(),
        })
    }

    #[getter]
    fn selection(&self) -> StereoSelection {
        StereoSelection::from_rust(self.to_rust().selection.clone())
    }
    #[getter]
    fn automorphism_algorithm(&self) -> AutomorphismAlgorithm {
        AutomorphismAlgorithm::from_rust(self.to_rust().automorphism_algorithm)
    }

    #[getter]
    fn mode(&self) -> RemoveStereoMode {
        RemoveStereoMode::from_rust(self.to_rust().mode)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!("RemoveStereoConfig(selection={}, mode=RemoveStereoMode.{:?}, automorphism_algorithm={})", self.selection().__repr__(py)?, self.mode(), self.automorphism_algorithm().repr()))
    }
}

impl RemoveStereoConfig {
    pub(crate) fn from_rust(config: GraphRemoveStereoConfig) -> Self {
        Self(config)
    }
    pub(crate) fn to_rust(&self) -> &GraphRemoveStereoConfig {
        &self.0
    }
}

/// Immutable operational configuration corresponding to Rust EnumerateStereoConfig.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct EnumerateStereoConfig(GraphEnumerateStereoConfig);

#[pymethods]
impl EnumerateStereoConfig {
    #[new]
    #[pyo3(signature = (*, selection=None, automorphism_algorithm=AutomorphismAlgorithm::Nauty()))]
    fn new(
        selection: Option<&StereoSelection>,
        automorphism_algorithm: AutomorphismAlgorithm,
    ) -> Self {
        Self::from_rust(GraphEnumerateStereoConfig {
            selection: selection.map_or_else(
                || GraphEnumerateStereoConfig::default().selection,
                |value| value.to_rust().clone(),
            ),
            automorphism_algorithm: automorphism_algorithm.to_rust(),
        })
    }

    #[getter]
    fn selection(&self) -> StereoSelection {
        StereoSelection::from_rust(self.to_rust().selection.clone())
    }
    #[getter]
    fn automorphism_algorithm(&self) -> AutomorphismAlgorithm {
        AutomorphismAlgorithm::from_rust(self.to_rust().automorphism_algorithm)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "EnumerateStereoConfig(selection={}, automorphism_algorithm={})",
            self.selection().__repr__(py)?,
            self.automorphism_algorithm().repr()
        ))
    }
}

impl EnumerateStereoConfig {
    pub(crate) fn from_rust(config: GraphEnumerateStereoConfig) -> Self {
        Self(config)
    }
    pub(crate) fn to_rust(&self) -> &GraphEnumerateStereoConfig {
        &self.0
    }
}

#[pyclass(skip_from_py_object)]
pub(crate) struct EnumerateStereoIter {
    inner: GraphEnumerateStereoIter,
}

impl EnumerateStereoIter {
    pub(crate) fn from_rust(inner: GraphEnumerateStereoIter) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl EnumerateStereoIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    fn __next__(&mut self) -> Option<Molecule> {
        self.inner.next().map(Molecule::from_rust)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use rstest::rstest;
    use umol_graph::ops::stereo::{
        StereoSelection as GraphStereoSelection, StereoSite as GraphStereoSite,
    };
    use umol_graph::ops::transform::MaximumMatchingAlgorithm as GraphMaximumMatchingAlgorithm;
    use umol_graph_ir::ir::{
        BondId as GraphIrBondId, StereoKind as GraphIrStereoKind,
        Stereogenicity as GraphIrStereogenicity,
    };

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

    #[rstest]
    #[case::unrestricted(GraphStereoSelection::default())]
    #[case::empty(GraphStereoSelection { sites: Some(BTreeSet::new()), kinds: Some(BTreeSet::new()), stereogenicities: Some(BTreeSet::new()) })]
    #[case::selected(GraphStereoSelection {
        sites: Some(BTreeSet::from([GraphStereoSite::Bond(GraphIrBondId(3))])),
        kinds: Some(BTreeSet::from([GraphIrStereoKind::SquarePlanar])),
        stereogenicities: Some(BTreeSet::from([GraphIrStereogenicity::Prochiral])),
    })]
    fn test_add_stereo_config_from_rust(#[case] selection: GraphStereoSelection) {
        let config = GraphAddStereoConfig {
            selection: selection.clone(),
            ..Default::default()
        };
        let bound = AddStereoConfig::from_rust(config.clone());
        assert_eq!(bound.to_rust(), &config);
        assert_eq!(bound.selection().to_rust(), &selection);
        assert_eq!(
            bound.automorphism_algorithm().to_rust(),
            config.automorphism_algorithm
        );
        assert_eq!(bound.ring_config().to_rust(), config.ring_config);
    }

    #[rstest]
    #[case::unrestricted(GraphStereoSelection::default())]
    #[case::empty(GraphStereoSelection { sites: Some(BTreeSet::new()), kinds: Some(BTreeSet::new()), stereogenicities: Some(BTreeSet::new()) })]
    #[case::selected(GraphStereoSelection {
        sites: Some(BTreeSet::from([GraphStereoSite::Bond(GraphIrBondId(3))])),
        kinds: Some(BTreeSet::from([GraphIrStereoKind::SquarePlanar])),
        stereogenicities: Some(BTreeSet::from([GraphIrStereogenicity::Prochiral])),
    })]
    fn test_remove_stereo_config_from_rust(#[case] selection: GraphStereoSelection) {
        let config = GraphRemoveStereoConfig {
            selection: selection.clone(),
            mode: GraphRemoveStereoMode::Undetermined,
            ..Default::default()
        };
        let bound = RemoveStereoConfig::from_rust(config.clone());
        assert_eq!(bound.to_rust(), &config);
        assert_eq!(bound.selection().to_rust(), &selection);
        assert_eq!(
            bound.automorphism_algorithm().to_rust(),
            config.automorphism_algorithm
        );
        assert_eq!(bound.mode().to_rust(), config.mode);
    }

    #[rstest]
    #[case::unrestricted(GraphStereoSelection::default())]
    #[case::empty(GraphStereoSelection { sites: Some(BTreeSet::new()), kinds: Some(BTreeSet::new()), stereogenicities: Some(BTreeSet::new()) })]
    #[case::selected(GraphStereoSelection {
        sites: Some(BTreeSet::from([GraphStereoSite::Bond(GraphIrBondId(3))])),
        kinds: Some(BTreeSet::from([GraphIrStereoKind::SquarePlanar])),
        stereogenicities: Some(BTreeSet::from([GraphIrStereogenicity::Prochiral])),
    })]
    fn test_enumerate_stereo_config_from_rust(#[case] selection: GraphStereoSelection) {
        let config = GraphEnumerateStereoConfig {
            selection: selection.clone(),
            ..Default::default()
        };
        let bound = EnumerateStereoConfig::from_rust(config.clone());
        assert_eq!(bound.to_rust(), &config);
        assert_eq!(bound.selection().to_rust(), &selection);
        assert_eq!(
            bound.automorphism_algorithm().to_rust(),
            config.automorphism_algorithm
        );
    }
}
