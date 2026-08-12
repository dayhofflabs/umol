//! Python configuration and aggregate operations for canonicalization.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use pyo3::prelude::*;
use umol_graph::ops::canonicalize::CanonicalizationConfig as GraphCanonicalizationConfig;
use umol_graph::ops::model::StereoModel as GraphStereoModel;
use umol_graph_ir::ir::{
    CanonicalizationContext as GraphIrCanonicalizationContext,
    CanonicalizationLevel as GraphIrCanonicalizationLevel, Canonicalize,
    MoleculeCanonicalizationError as GraphIrMoleculeCanonicalizationError,
    ReactionCanonicalizationError as GraphIrReactionCanonicalizationError,
    ReactionSpanCanonicalizationError as GraphIrReactionSpanCanonicalizationError,
};

use crate::algorithm::AutomorphismAlgorithm;
use crate::error::{ContradictionError, InvalidStructureError};
use crate::model::stereo::StereoModel;
use crate::molecule::Molecule;
use crate::reaction::Reaction;
use crate::reaction_span::ReactionSpan;

/// Nested structural layer used to select or compare a canonical entity frame.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum CanonicalizationLevel {
    Topology,
    Constitution,
    Structure,
    Full,
}

impl CanonicalizationLevel {
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "Rust-to-Python conversion is part of the canonicalization binding contract"
        )
    )]
    pub(crate) fn from_rust(level: GraphIrCanonicalizationLevel) -> Self {
        match level {
            GraphIrCanonicalizationLevel::Topology => Self::Topology,
            GraphIrCanonicalizationLevel::Constitution => Self::Constitution,
            GraphIrCanonicalizationLevel::Structure => Self::Structure,
            GraphIrCanonicalizationLevel::Full => Self::Full,
        }
    }

    pub(crate) fn to_rust(self) -> GraphIrCanonicalizationLevel {
        match self {
            Self::Topology => GraphIrCanonicalizationLevel::Topology,
            Self::Constitution => GraphIrCanonicalizationLevel::Constitution,
            Self::Structure => GraphIrCanonicalizationLevel::Structure,
            Self::Full => GraphIrCanonicalizationLevel::Full,
        }
    }
}

/// Operational configuration for aggregate canonicalization.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CanonicalizationConfig {
    automorphism_algorithm: AutomorphismAlgorithm,
}

impl Default for CanonicalizationConfig {
    fn default() -> Self {
        Self::from_rust(GraphCanonicalizationConfig::default())
    }
}

#[pymethods]
impl CanonicalizationConfig {
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
            return "CanonicalizationConfig.default()".to_owned();
        }
        format!(
            "CanonicalizationConfig(automorphism_algorithm={})",
            self.automorphism_algorithm.repr(),
        )
    }
}

impl CanonicalizationConfig {
    pub(crate) fn from_rust(config: GraphCanonicalizationConfig) -> Self {
        Self {
            automorphism_algorithm: AutomorphismAlgorithm::from_rust(config.automorphism_algorithm),
        }
    }

    pub(crate) fn to_rust(self) -> GraphCanonicalizationConfig {
        GraphCanonicalizationConfig {
            automorphism_algorithm: self.automorphism_algorithm.to_rust(),
        }
    }
}

fn canonicalization_context(
    stereo_model: Option<StereoModel>,
    config: Option<CanonicalizationConfig>,
) -> GraphIrCanonicalizationContext {
    let model = stereo_model.map_or_else(GraphStereoModel::default, |model| model.to_rust());
    config.unwrap_or_default().to_rust().context(&model)
}

fn molecule_canonicalization_error(error: GraphIrMoleculeCanonicalizationError) -> PyErr {
    match error {
        GraphIrMoleculeCanonicalizationError::Integrity(error) => {
            InvalidStructureError::new_err(error.to_string())
        }
        GraphIrMoleculeCanonicalizationError::Contradiction(error) => {
            ContradictionError::new_err(error.to_string())
        }
    }
}

fn reaction_span_canonicalization_error(error: GraphIrReactionSpanCanonicalizationError) -> PyErr {
    match error {
        GraphIrReactionSpanCanonicalizationError::Integrity(error) => {
            InvalidStructureError::new_err(error.to_string())
        }
        GraphIrReactionSpanCanonicalizationError::Contradiction(error) => {
            ContradictionError::new_err(error.to_string())
        }
    }
}

fn reaction_canonicalization_error(error: GraphIrReactionCanonicalizationError) -> PyErr {
    match error {
        GraphIrReactionCanonicalizationError::Integrity(error) => {
            InvalidStructureError::new_err(error.to_string())
        }
        GraphIrReactionCanonicalizationError::Contradiction(error) => {
            ContradictionError::new_err(error.to_string())
        }
    }
}

#[pymethods]
impl Molecule {
    /// Return the complete canonical form without changing this molecule.
    #[pyo3(signature = (*, stereo_model=None, config=None))]
    fn canonicalize(
        &self,
        stereo_model: Option<StereoModel>,
        config: Option<CanonicalizationConfig>,
    ) -> PyResult<Self> {
        self.to_rust()
            .clone()
            .canonicalize(&canonicalization_context(stereo_model, config))
            .map(Self::from_rust)
            .map_err(molecule_canonicalization_error)
    }

    /// Return this molecule in the canonical frame selected at `level`.
    #[pyo3(signature = (level, *, stereo_model=None, config=None))]
    fn canonicalize_by(
        &self,
        level: CanonicalizationLevel,
        stereo_model: Option<StereoModel>,
        config: Option<CanonicalizationConfig>,
    ) -> PyResult<Self> {
        self.to_rust()
            .clone()
            .canonicalize_by(
                level.to_rust(),
                &canonicalization_context(stereo_model, config),
            )
            .map(Self::from_rust)
            .map_err(molecule_canonicalization_error)
    }

    /// Compare complete canonical forms under the same model and config.
    #[pyo3(signature = (other, *, stereo_model=None, config=None))]
    fn canonical_eq(
        &self,
        other: &Self,
        stereo_model: Option<StereoModel>,
        config: Option<CanonicalizationConfig>,
    ) -> bool {
        self.to_rust().canonical_eq(
            other.to_rust(),
            &canonicalization_context(stereo_model, config),
        )
    }

    /// Compare canonical forms at `level` under the same model and config.
    #[pyo3(signature = (other, level, *, stereo_model=None, config=None))]
    fn canonical_eq_by(
        &self,
        other: &Self,
        level: CanonicalizationLevel,
        stereo_model: Option<StereoModel>,
        config: Option<CanonicalizationConfig>,
    ) -> bool {
        self.to_rust().canonical_eq_by(
            other.to_rust(),
            level.to_rust(),
            &canonicalization_context(stereo_model, config),
        )
    }
}

#[pymethods]
impl ReactionSpan {
    /// Return the complete canonical form without changing this reaction span.
    #[pyo3(signature = (*, stereo_model=None, config=None))]
    fn canonicalize(
        &self,
        stereo_model: Option<StereoModel>,
        config: Option<CanonicalizationConfig>,
    ) -> PyResult<Self> {
        self.to_rust()
            .clone()
            .canonicalize(&canonicalization_context(stereo_model, config))
            .map(Self::from_rust)
            .map_err(reaction_span_canonicalization_error)
    }

    /// Return this reaction span in the canonical frame selected at `level`.
    #[pyo3(signature = (level, *, stereo_model=None, config=None))]
    fn canonicalize_by(
        &self,
        level: CanonicalizationLevel,
        stereo_model: Option<StereoModel>,
        config: Option<CanonicalizationConfig>,
    ) -> PyResult<Self> {
        self.to_rust()
            .clone()
            .canonicalize_by(
                level.to_rust(),
                &canonicalization_context(stereo_model, config),
            )
            .map(Self::from_rust)
            .map_err(reaction_span_canonicalization_error)
    }

    /// Compare complete canonical forms under the same model and config.
    #[pyo3(signature = (other, *, stereo_model=None, config=None))]
    fn canonical_eq(
        &self,
        other: &Self,
        stereo_model: Option<StereoModel>,
        config: Option<CanonicalizationConfig>,
    ) -> bool {
        self.to_rust().canonical_eq(
            other.to_rust(),
            &canonicalization_context(stereo_model, config),
        )
    }

    /// Compare canonical forms at `level` under the same model and config.
    #[pyo3(signature = (other, level, *, stereo_model=None, config=None))]
    fn canonical_eq_by(
        &self,
        other: &Self,
        level: CanonicalizationLevel,
        stereo_model: Option<StereoModel>,
        config: Option<CanonicalizationConfig>,
    ) -> bool {
        self.to_rust().canonical_eq_by(
            other.to_rust(),
            level.to_rust(),
            &canonicalization_context(stereo_model, config),
        )
    }
}

#[pymethods]
impl Reaction {
    /// Return the complete canonical form without changing this reaction.
    #[pyo3(signature = (*, stereo_model=None, config=None))]
    fn canonicalize(
        &self,
        py: Python<'_>,
        stereo_model: Option<StereoModel>,
        config: Option<CanonicalizationConfig>,
    ) -> PyResult<Self> {
        let canonical = self
            .to_rust(py)
            .canonicalize(&canonicalization_context(stereo_model, config))
            .map_err(reaction_canonicalization_error)?;
        Self::from_rust(py, canonical)
    }

    /// Return this reaction in the canonical frame selected at `level`.
    #[pyo3(signature = (level, *, stereo_model=None, config=None))]
    fn canonicalize_by(
        &self,
        py: Python<'_>,
        level: CanonicalizationLevel,
        stereo_model: Option<StereoModel>,
        config: Option<CanonicalizationConfig>,
    ) -> PyResult<Self> {
        let canonical = self
            .to_rust(py)
            .canonicalize_by(
                level.to_rust(),
                &canonicalization_context(stereo_model, config),
            )
            .map_err(reaction_canonicalization_error)?;
        Self::from_rust(py, canonical)
    }

    /// Compare complete canonical forms under the same model and config.
    #[pyo3(signature = (other, *, stereo_model=None, config=None))]
    fn canonical_eq(
        &self,
        other: &Self,
        py: Python<'_>,
        stereo_model: Option<StereoModel>,
        config: Option<CanonicalizationConfig>,
    ) -> bool {
        self.to_rust(py).canonical_eq(
            &other.to_rust(py),
            &canonicalization_context(stereo_model, config),
        )
    }

    /// Compare canonical forms at `level` under the same model and config.
    #[pyo3(signature = (other, level, *, stereo_model=None, config=None))]
    fn canonical_eq_by(
        &self,
        other: &Self,
        py: Python<'_>,
        level: CanonicalizationLevel,
        stereo_model: Option<StereoModel>,
        config: Option<CanonicalizationConfig>,
    ) -> bool {
        self.to_rust(py).canonical_eq_by(
            &other.to_rust(py),
            level.to_rust(),
            &canonicalization_context(stereo_model, config),
        )
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_core::AutomorphismAlgorithm as GraphCoreAutomorphismAlgorithm;

    use super::*;

    #[rstest]
    #[case::topology(
        GraphIrCanonicalizationLevel::Topology,
        CanonicalizationLevel::Topology
    )]
    #[case::constitution(
        GraphIrCanonicalizationLevel::Constitution,
        CanonicalizationLevel::Constitution
    )]
    #[case::structure(
        GraphIrCanonicalizationLevel::Structure,
        CanonicalizationLevel::Structure
    )]
    #[case::full(GraphIrCanonicalizationLevel::Full, CanonicalizationLevel::Full)]
    fn test_canonicalization_level_conversion(
        #[case] rust: GraphIrCanonicalizationLevel,
        #[case] python: CanonicalizationLevel,
    ) {
        assert_eq!(CanonicalizationLevel::from_rust(rust), python);
        assert_eq!(python.to_rust(), rust);
    }

    #[rstest]
    #[case::default(CanonicalizationConfig::default())]
    fn test_canonicalization_config_new(#[case] expected: CanonicalizationConfig) {
        assert_eq!(
            CanonicalizationConfig::new(AutomorphismAlgorithm::Nauty()),
            expected,
        );
    }

    #[rstest]
    #[case::default(CanonicalizationConfig::default(), "CanonicalizationConfig.default()")]
    fn test_canonicalization_config_repr(
        #[case] config: CanonicalizationConfig,
        #[case] expected: &str,
    ) {
        assert_eq!(config.__repr__(), expected);
    }

    #[rstest]
    #[case::default(GraphCanonicalizationConfig::default())]
    fn test_canonicalization_config_conversion(#[case] config: GraphCanonicalizationConfig) {
        assert_eq!(CanonicalizationConfig::from_rust(config).to_rust(), config);
    }

    #[rstest]
    #[case::without_para_stereo(false)]
    #[case::with_para_stereo(true)]
    fn test_canonicalization_context(#[case] para_stereo: bool) {
        let model = StereoModel::from_rust(&GraphStereoModel {
            para_stereo,
            ..GraphStereoModel::default()
        });

        assert_eq!(
            canonicalization_context(Some(model), Some(CanonicalizationConfig::default())),
            GraphIrCanonicalizationContext {
                para_stereo,
                automorphism_algorithm: GraphCoreAutomorphismAlgorithm::Nauty,
            },
        );
    }
}
