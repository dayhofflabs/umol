//! Python configuration and aggregate operations for canonicalization.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use pyo3::prelude::*;
use umol_graph::ops::canonicalize::CanonicalizeConfig as GraphCanonicalizeConfig;
use umol_graph::ops::model::StereoModel as GraphStereoModel;
use umol_graph_ir::ir::{
    Canonicalize, CanonicalizeContext as GraphIrCanonicalizeContext,
    MoleculeCanonicalizeError as GraphIrMoleculeCanonicalizeError,
    ReactionCanonicalizeError as GraphIrReactionCanonicalizeError,
    ReactionSpanCanonicalizeError as GraphIrReactionSpanCanonicalizeError,
};

use crate::algorithm::AutomorphismAlgorithm;
use crate::error::ContradictionError;
use crate::model::stereo::StereoModel;
use crate::molecule::Molecule;
use crate::reaction::Reaction;
use crate::reaction_span::ReactionSpan;
use crate::remap::MoleculeRemapping;

/// Operational configuration for aggregate canonicalization.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CanonicalizeConfig {
    automorphism_algorithm: AutomorphismAlgorithm,
}

impl Default for CanonicalizeConfig {
    fn default() -> Self {
        Self::from_rust(GraphCanonicalizeConfig::default())
    }
}

#[pymethods]
impl CanonicalizeConfig {
    #[new]
    #[pyo3(signature = (*, automorphism_algorithm=AutomorphismAlgorithm::Nauty()))]
    fn new(automorphism_algorithm: AutomorphismAlgorithm) -> Self {
        Self {
            automorphism_algorithm,
        }
    }

    #[staticmethod]
    fn default() -> Self {
        Default::default()
    }

    #[getter]
    fn automorphism_algorithm(&self) -> AutomorphismAlgorithm {
        self.automorphism_algorithm
    }

    fn __repr__(&self) -> String {
        if self == &Self::default() {
            return "CanonicalizeConfig.default()".to_owned();
        }
        format!(
            "CanonicalizeConfig(automorphism_algorithm={})",
            self.automorphism_algorithm.repr(),
        )
    }
}

impl CanonicalizeConfig {
    pub(crate) fn from_rust(config: GraphCanonicalizeConfig) -> Self {
        Self {
            automorphism_algorithm: AutomorphismAlgorithm::from_rust(config.automorphism_algorithm),
        }
    }

    pub(crate) fn to_rust(self) -> GraphCanonicalizeConfig {
        GraphCanonicalizeConfig {
            automorphism_algorithm: self.automorphism_algorithm.to_rust(),
        }
    }
}

fn canonicalize_context(
    stereo_model: Option<StereoModel>,
    config: Option<CanonicalizeConfig>,
) -> GraphIrCanonicalizeContext {
    let model = stereo_model.map_or_else(GraphStereoModel::default, |model| model.to_rust());
    config.unwrap_or_default().to_rust().context(&model)
}

fn molecule_canonicalization_error(error: GraphIrMoleculeCanonicalizeError) -> PyErr {
    match error {
        GraphIrMoleculeCanonicalizeError::Contradiction(error) => {
            ContradictionError::new_err(error.to_string())
        }
    }
}

fn reaction_span_canonicalization_error(error: GraphIrReactionSpanCanonicalizeError) -> PyErr {
    match error {
        GraphIrReactionSpanCanonicalizeError::Contradiction(error) => {
            ContradictionError::new_err(error.to_string())
        }
    }
}

fn reaction_canonicalization_error(error: GraphIrReactionCanonicalizeError) -> PyErr {
    match error {
        GraphIrReactionCanonicalizeError::Contradiction(error) => {
            ContradictionError::new_err(error.to_string())
        }
    }
}

#[pymethods]
impl Molecule {
    /// Return the complete canonical form without changing this molecule.
    ///
    /// Canonical representatives may change between umol 0.x releases and are not persistent ids.
    #[pyo3(signature = (*, stereo_model=None, config=None))]
    fn canonicalize(
        &self,
        stereo_model: Option<StereoModel>,
        config: Option<CanonicalizeConfig>,
    ) -> PyResult<Self> {
        self.to_rust()
            .clone()
            .canonicalize(&canonicalize_context(stereo_model, config))
            .map(Self::from_rust)
            .map_err(molecule_canonicalization_error)
    }

    /// Return the complete canonical form and its source-to-canonical remapping.
    ///
    /// The remapping is total across every entity kind and maps entity ids; participant
    /// frames are selected internally and are not encoded in it. Canonical representatives may
    /// change between umol 0.x releases and are not persistent ids.
    #[pyo3(signature = (*, stereo_model=None, config=None))]
    fn canonicalize_with_remapping(
        &self,
        stereo_model: Option<StereoModel>,
        config: Option<CanonicalizeConfig>,
    ) -> PyResult<(Self, MoleculeRemapping)> {
        self.to_rust()
            .clone()
            .canonicalize_with_remapping(&canonicalize_context(stereo_model, config))
            .map(|(canonical, remapping)| {
                (
                    Self::from_rust(canonical),
                    MoleculeRemapping::from_rust(remapping),
                )
            })
            .map_err(molecule_canonicalization_error)
    }

    /// Compare complete canonical forms under the same model and config.
    #[pyo3(signature = (other, *, stereo_model=None, config=None))]
    fn canonical_eq(
        &self,
        other: &Self,
        stereo_model: Option<StereoModel>,
        config: Option<CanonicalizeConfig>,
    ) -> bool {
        self.to_rust()
            .canonical_eq(other.to_rust(), &canonicalize_context(stereo_model, config))
    }
}

#[pymethods]
impl ReactionSpan {
    /// Return the complete canonical form without changing this reaction span.
    ///
    /// Canonical representatives may change between umol 0.x releases and are not persistent ids.
    #[pyo3(signature = (*, stereo_model=None, config=None))]
    fn canonicalize(
        &self,
        stereo_model: Option<StereoModel>,
        config: Option<CanonicalizeConfig>,
    ) -> PyResult<Self> {
        self.to_rust()
            .clone()
            .canonicalize(&canonicalize_context(stereo_model, config))
            .map(Self::from_rust)
            .map_err(reaction_span_canonicalization_error)
    }

    /// Return the complete canonical form and its source-to-canonical remapping.
    ///
    /// The remapping is total across every union-frame entity kind and maps entity ids;
    /// participant frames are selected internally and are not encoded in it. Canonical
    /// representatives may change between umol 0.x releases and are not persistent ids.
    #[pyo3(signature = (*, stereo_model=None, config=None))]
    fn canonicalize_with_remapping(
        &self,
        stereo_model: Option<StereoModel>,
        config: Option<CanonicalizeConfig>,
    ) -> PyResult<(Self, MoleculeRemapping)> {
        self.to_rust()
            .clone()
            .canonicalize_with_remapping(&canonicalize_context(stereo_model, config))
            .map(|(canonical, remapping)| {
                (
                    Self::from_rust(canonical),
                    MoleculeRemapping::from_rust(remapping),
                )
            })
            .map_err(reaction_span_canonicalization_error)
    }

    /// Compare complete canonical forms under the same model and config.
    #[pyo3(signature = (other, *, stereo_model=None, config=None))]
    fn canonical_eq(
        &self,
        other: &Self,
        stereo_model: Option<StereoModel>,
        config: Option<CanonicalizeConfig>,
    ) -> bool {
        self.to_rust()
            .canonical_eq(other.to_rust(), &canonicalize_context(stereo_model, config))
    }
}

#[pymethods]
impl Reaction {
    /// Return the complete canonical form without changing this reaction.
    ///
    /// Canonical representatives may change between umol 0.x releases and are not persistent ids.
    #[pyo3(signature = (*, stereo_model=None, config=None))]
    fn canonicalize(
        &self,
        py: Python<'_>,
        stereo_model: Option<StereoModel>,
        config: Option<CanonicalizeConfig>,
    ) -> PyResult<Self> {
        let canonical = self
            .to_rust(py)?
            .canonicalize(&canonicalize_context(stereo_model, config))
            .map_err(reaction_canonicalization_error)?;
        Self::from_rust(py, canonical)
    }

    /// Return the complete canonical form and its source-to-canonical remapping.
    ///
    /// The remapping is total across every materialized union-frame entity kind and maps
    /// entity ids; participant frames are selected internally and are not encoded in it. Canonical
    /// representatives may change between umol 0.x releases and are not persistent ids.
    #[pyo3(signature = (*, stereo_model=None, config=None))]
    fn canonicalize_with_remapping(
        &self,
        py: Python<'_>,
        stereo_model: Option<StereoModel>,
        config: Option<CanonicalizeConfig>,
    ) -> PyResult<(Self, MoleculeRemapping)> {
        let (canonical, remapping) = self
            .to_rust(py)?
            .canonicalize_with_remapping(&canonicalize_context(stereo_model, config))
            .map_err(reaction_canonicalization_error)?;
        Ok((
            Self::from_rust(py, canonical)?,
            MoleculeRemapping::from_rust(remapping),
        ))
    }

    /// Compare complete canonical forms under the same model and config.
    #[pyo3(signature = (other, *, stereo_model=None, config=None))]
    fn canonical_eq(
        &self,
        other: &Self,
        py: Python<'_>,
        stereo_model: Option<StereoModel>,
        config: Option<CanonicalizeConfig>,
    ) -> PyResult<bool> {
        Ok(self.to_rust(py)?.canonical_eq(
            &other.to_rust(py)?,
            &canonicalize_context(stereo_model, config),
        ))
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_core::AutomorphismAlgorithm as GraphCoreAutomorphismAlgorithm;

    use super::*;

    #[rstest]
    #[case::default(CanonicalizeConfig::default())]
    fn test_canonicalize_config_new(#[case] expected: CanonicalizeConfig) {
        assert_eq!(
            CanonicalizeConfig::new(AutomorphismAlgorithm::Nauty()),
            expected,
        );
    }

    #[rstest]
    #[case::default(CanonicalizeConfig::default(), "CanonicalizeConfig.default()")]
    fn test_canonicalize_config_repr(#[case] config: CanonicalizeConfig, #[case] expected: &str) {
        assert_eq!(config.__repr__(), expected);
    }

    #[rstest]
    #[case::default(GraphCanonicalizeConfig::default())]
    fn test_canonicalize_config_conversion(#[case] config: GraphCanonicalizeConfig) {
        assert_eq!(CanonicalizeConfig::from_rust(config).to_rust(), config);
    }

    #[rstest]
    #[case::without_para_stereo(false)]
    #[case::with_para_stereo(true)]
    fn test_canonicalize_context(#[case] para_stereo: bool) {
        let model = StereoModel::from_rust(&GraphStereoModel {
            para_stereo,
            ..GraphStereoModel::default()
        });

        assert_eq!(
            canonicalize_context(Some(model), Some(CanonicalizeConfig::default())),
            GraphIrCanonicalizeContext {
                para_stereo,
                automorphism_algorithm: GraphCoreAutomorphismAlgorithm::Nauty,
            },
        );
    }
}
