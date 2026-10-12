//! Molecule-resolution configuration, outcomes, and reports.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use pyo3::prelude::*;
use umol_graph::ops::resolve::{
    AromaticityFailurePolicy as GraphAromaticityFailurePolicy,
    AromaticityResolveConfig as GraphAromaticityResolveConfig, IsotopePolicy as GraphIsotopePolicy,
    ResolveConfig as GraphResolveConfig, ResolveContradiction as GraphResolveContradiction,
    StereoFailurePolicy as GraphStereoFailurePolicy,
    StereoResolveConfig as GraphStereoResolveConfig,
};
use umol_graph::ops::valence::{
    AtomCompletions as GraphAtomCompletions, ResolveReport as GraphResolveReport,
};
use umol_graph_ir::ir::AtomId as GraphIrAtomId;

use crate::atom::AtomForm;
use crate::model::aromaticity::AromaticityConfig;
use crate::molecule::Molecule;

/// Completion of unspecified isotope composition, independently of valence policy.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum IsotopePolicy {
    Strict,
    Natural,
}

impl IsotopePolicy {
    pub(crate) fn from_rust(policy: GraphIsotopePolicy) -> Self {
        match policy {
            GraphIsotopePolicy::Strict => Self::Strict,
            GraphIsotopePolicy::Natural => Self::Natural,
        }
    }

    pub(crate) fn to_rust(self) -> GraphIsotopePolicy {
        match self {
            Self::Strict => GraphIsotopePolicy::Strict,
            Self::Natural => GraphIsotopePolicy::Natural,
        }
    }
}

/// Policy for an independently invalid aromatic constraint or entity.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum AromaticityFailurePolicy {
    Error,
    Keep,
}

impl AromaticityFailurePolicy {
    pub(crate) fn from_rust(policy: GraphAromaticityFailurePolicy) -> Self {
        match policy {
            GraphAromaticityFailurePolicy::Error => Self::Error,
            GraphAromaticityFailurePolicy::Keep => Self::Keep,
        }
    }

    pub(crate) fn to_rust(self) -> GraphAromaticityFailurePolicy {
        match self {
            Self::Error => GraphAromaticityFailurePolicy::Error,
            Self::Keep => GraphAromaticityFailurePolicy::Keep,
        }
    }
}

/// Policy for an independently invalid stereo constraint or entity.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum StereoFailurePolicy {
    Error,
    Keep,
    Remove,
}

impl StereoFailurePolicy {
    pub(crate) fn from_rust(policy: GraphStereoFailurePolicy) -> Self {
        match policy {
            GraphStereoFailurePolicy::Error => Self::Error,
            GraphStereoFailurePolicy::Keep => Self::Keep,
            GraphStereoFailurePolicy::Remove => Self::Remove,
        }
    }

    pub(crate) fn to_rust(self) -> GraphStereoFailurePolicy {
        match self {
            Self::Error => GraphStereoFailurePolicy::Error,
            Self::Keep => GraphStereoFailurePolicy::Keep,
            Self::Remove => GraphStereoFailurePolicy::Remove,
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
    #[pyo3(signature = (*, perception=AromaticityConfig::default(), aromatic_valence_failure=AromaticityFailurePolicy::Error, aromatic_system_failure=AromaticityFailurePolicy::Error, reset_aromatic_valence=false))]
    fn new(
        perception: AromaticityConfig,
        aromatic_valence_failure: AromaticityFailurePolicy,
        aromatic_system_failure: AromaticityFailurePolicy,
        reset_aromatic_valence: bool,
    ) -> Self {
        Self(GraphAromaticityResolveConfig {
            perception: perception.to_rust(),
            aromatic_valence_failure: aromatic_valence_failure.to_rust(),
            aromatic_system_failure: aromatic_system_failure.to_rust(),
            reset_aromatic_valence,
        })
    }

    #[getter]
    fn perception(&self) -> AromaticityConfig {
        AromaticityConfig::from_rust(self.0.perception)
    }

    #[getter]
    fn aromatic_valence_failure(&self) -> AromaticityFailurePolicy {
        AromaticityFailurePolicy::from_rust(self.0.aromatic_valence_failure)
    }

    #[getter]
    fn aromatic_system_failure(&self) -> AromaticityFailurePolicy {
        AromaticityFailurePolicy::from_rust(self.0.aromatic_system_failure)
    }

    #[getter]
    fn reset_aromatic_valence(&self) -> bool {
        self.0.reset_aromatic_valence
    }

    fn __repr__(&self) -> String {
        format!(
            "AromaticityResolveConfig(perception={}, aromatic_valence_failure=AromaticityFailurePolicy.{:?}, aromatic_system_failure=AromaticityFailurePolicy.{:?}, reset_aromatic_valence={})",
            self.perception().__repr__(),
            self.aromatic_valence_failure(),
            self.aromatic_system_failure(),
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
    #[pyo3(signature = (*, tetrahedral_stereo_failure=StereoFailurePolicy::Error, stereo_atom_failure=StereoFailurePolicy::Error, cis_trans_stereo_failure=StereoFailurePolicy::Error, stereo_bond_failure=StereoFailurePolicy::Error, reset_stereo_constraints=false))]
    fn new(
        tetrahedral_stereo_failure: StereoFailurePolicy,
        stereo_atom_failure: StereoFailurePolicy,
        cis_trans_stereo_failure: StereoFailurePolicy,
        stereo_bond_failure: StereoFailurePolicy,
        reset_stereo_constraints: bool,
    ) -> Self {
        Self(GraphStereoResolveConfig {
            tetrahedral_stereo_failure: tetrahedral_stereo_failure.to_rust(),
            stereo_atom_failure: stereo_atom_failure.to_rust(),
            cis_trans_stereo_failure: cis_trans_stereo_failure.to_rust(),
            stereo_bond_failure: stereo_bond_failure.to_rust(),
            reset_stereo_constraints,
        })
    }

    #[getter]
    fn tetrahedral_stereo_failure(&self) -> StereoFailurePolicy {
        StereoFailurePolicy::from_rust(self.0.tetrahedral_stereo_failure)
    }

    #[getter]
    fn stereo_atom_failure(&self) -> StereoFailurePolicy {
        StereoFailurePolicy::from_rust(self.0.stereo_atom_failure)
    }

    #[getter]
    fn cis_trans_stereo_failure(&self) -> StereoFailurePolicy {
        StereoFailurePolicy::from_rust(self.0.cis_trans_stereo_failure)
    }

    #[getter]
    fn stereo_bond_failure(&self) -> StereoFailurePolicy {
        StereoFailurePolicy::from_rust(self.0.stereo_bond_failure)
    }

    #[getter]
    fn reset_stereo_constraints(&self) -> bool {
        self.0.reset_stereo_constraints
    }

    fn __repr__(&self) -> String {
        format!(
            "StereoResolveConfig(tetrahedral_stereo_failure=StereoFailurePolicy.{:?}, stereo_atom_failure=StereoFailurePolicy.{:?}, cis_trans_stereo_failure=StereoFailurePolicy.{:?}, stereo_bond_failure=StereoFailurePolicy.{:?}, reset_stereo_constraints={})",
            self.tetrahedral_stereo_failure(),
            self.stereo_atom_failure(),
            self.cis_trans_stereo_failure(),
            self.stereo_bond_failure(),
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

/// Operational policy for molecule resolution. Isotope completion defaults to Strict.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ResolveConfig(GraphResolveConfig);

#[pymethods]
impl ResolveConfig {
    #[new]
    #[pyo3(signature = (*, aromaticity, stereo, isotope=IsotopePolicy::Strict))]
    fn new(
        aromaticity: AromaticityResolveConfig,
        stereo: StereoResolveConfig,
        isotope: IsotopePolicy,
    ) -> Self {
        Self(GraphResolveConfig {
            isotope: isotope.to_rust(),
            aromaticity: aromaticity.to_rust(),
            stereo: stereo.to_rust(),
        })
    }

    #[staticmethod]
    fn default() -> Self {
        Self::from_rust(GraphResolveConfig::default())
    }

    #[getter]
    fn isotope(&self) -> IsotopePolicy {
        IsotopePolicy::from_rust(self.0.isotope)
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
            "ResolveConfig(isotope=IsotopePolicy.{:?}, aromaticity={}, stereo={})",
            self.isotope(),
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
    use umol_graph::ops::aromaticity::AromaticityContradiction as GraphAromaticityContradiction;

    use super::*;

    #[rstest]
    #[case::strict(GraphIsotopePolicy::Strict, IsotopePolicy::Strict)]
    #[case::natural(GraphIsotopePolicy::Natural, IsotopePolicy::Natural)]
    fn test_isotope_policy_from_rust(
        #[case] policy: GraphIsotopePolicy,
        #[case] expected: IsotopePolicy,
    ) {
        assert_eq!(IsotopePolicy::from_rust(policy), expected);
    }

    #[rstest]
    #[case::strict(IsotopePolicy::Strict, GraphIsotopePolicy::Strict)]
    #[case::natural(IsotopePolicy::Natural, GraphIsotopePolicy::Natural)]
    fn test_isotope_policy_to_rust(
        #[case] policy: IsotopePolicy,
        #[case] expected: GraphIsotopePolicy,
    ) {
        assert_eq!(policy.to_rust(), expected);
    }

    #[rstest]
    #[case::error(GraphAromaticityFailurePolicy::Error, AromaticityFailurePolicy::Error)]
    #[case::keep(GraphAromaticityFailurePolicy::Keep, AromaticityFailurePolicy::Keep)]
    fn test_aromaticity_failure_policy_from_rust(
        #[case] policy: GraphAromaticityFailurePolicy,
        #[case] expected: AromaticityFailurePolicy,
    ) {
        assert_eq!(AromaticityFailurePolicy::from_rust(policy), expected);
    }

    #[rstest]
    #[case::error(AromaticityFailurePolicy::Error, GraphAromaticityFailurePolicy::Error)]
    #[case::keep(AromaticityFailurePolicy::Keep, GraphAromaticityFailurePolicy::Keep)]
    fn test_aromaticity_failure_policy_to_rust(
        #[case] policy: AromaticityFailurePolicy,
        #[case] expected: GraphAromaticityFailurePolicy,
    ) {
        assert_eq!(policy.to_rust(), expected);
    }

    #[rstest]
    #[case::error(GraphStereoFailurePolicy::Error, StereoFailurePolicy::Error)]
    #[case::keep(GraphStereoFailurePolicy::Keep, StereoFailurePolicy::Keep)]
    #[case::remove(GraphStereoFailurePolicy::Remove, StereoFailurePolicy::Remove)]
    fn test_stereo_failure_policy_from_rust(
        #[case] policy: GraphStereoFailurePolicy,
        #[case] expected: StereoFailurePolicy,
    ) {
        assert_eq!(StereoFailurePolicy::from_rust(policy), expected);
    }

    #[rstest]
    #[case::error(StereoFailurePolicy::Error, GraphStereoFailurePolicy::Error)]
    #[case::keep(StereoFailurePolicy::Keep, GraphStereoFailurePolicy::Keep)]
    #[case::remove(StereoFailurePolicy::Remove, GraphStereoFailurePolicy::Remove)]
    fn test_stereo_failure_policy_to_rust(
        #[case] policy: StereoFailurePolicy,
        #[case] expected: GraphStereoFailurePolicy,
    ) {
        assert_eq!(policy.to_rust(), expected);
    }

    #[rstest]
    #[case::default(
        AromaticityConfig::default(),
        AromaticityFailurePolicy::Error,
        AromaticityFailurePolicy::Error,
        false,
        GraphAromaticityResolveConfig::default()
    )]
    #[case::configured(AromaticityConfig::default(), AromaticityFailurePolicy::Keep, AromaticityFailurePolicy::Keep, true, GraphAromaticityResolveConfig {
        perception: Default::default(),
        aromatic_valence_failure: GraphAromaticityFailurePolicy::Keep,
        aromatic_system_failure: GraphAromaticityFailurePolicy::Keep,
        reset_aromatic_valence: true,
    })]
    fn test_aromaticity_resolve_config_new(
        #[case] perception: AromaticityConfig,
        #[case] aromatic_valence_failure: AromaticityFailurePolicy,
        #[case] aromatic_system_failure: AromaticityFailurePolicy,
        #[case] reset_aromatic_valence: bool,
        #[case] expected: GraphAromaticityResolveConfig,
    ) {
        assert_eq!(
            AromaticityResolveConfig::new(
                perception,
                aromatic_valence_failure,
                aromatic_system_failure,
                reset_aromatic_valence,
            )
            .0,
            expected
        );
    }

    #[rstest]
    #[case::default(
        AromaticityResolveConfig::new(AromaticityConfig::default(), AromaticityFailurePolicy::Error, AromaticityFailurePolicy::Error, false),
        "AromaticityResolveConfig(perception=AromaticityConfig(ring_config=RingConfig(simple_cycle_algorithm=SimpleCycleEnumerationAlgorithm.ReadTarjan(), relevant_cycle_algorithm=RelevantCycleEnumerationAlgorithm.Vismara()), connected_components_algorithm=ConnectedComponentsAlgorithm.Bfs(), maximum_independent_set_algorithm=MaximumIndependentSetAlgorithm.BranchAndBound()), aromatic_valence_failure=AromaticityFailurePolicy.Error, aromatic_system_failure=AromaticityFailurePolicy.Error, reset_aromatic_valence=False)"
    )]
    #[case::nondefault(
        AromaticityResolveConfig::new(AromaticityConfig::default(), AromaticityFailurePolicy::Keep, AromaticityFailurePolicy::Keep, true),
        "AromaticityResolveConfig(perception=AromaticityConfig(ring_config=RingConfig(simple_cycle_algorithm=SimpleCycleEnumerationAlgorithm.ReadTarjan(), relevant_cycle_algorithm=RelevantCycleEnumerationAlgorithm.Vismara()), connected_components_algorithm=ConnectedComponentsAlgorithm.Bfs(), maximum_independent_set_algorithm=MaximumIndependentSetAlgorithm.BranchAndBound()), aromatic_valence_failure=AromaticityFailurePolicy.Keep, aromatic_system_failure=AromaticityFailurePolicy.Keep, reset_aromatic_valence=True)"
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
        aromatic_valence_failure: GraphAromaticityFailurePolicy::Keep,
        aromatic_system_failure: GraphAromaticityFailurePolicy::Keep,
        reset_aromatic_valence: true,
    })]
    fn test_aromaticity_resolve_config_from_rust(#[case] config: GraphAromaticityResolveConfig) {
        assert_eq!(AromaticityResolveConfig::from_rust(config).0, config);
    }

    #[rstest]
    #[case::default(AromaticityResolveConfig::new(
        AromaticityConfig::default(),
        AromaticityFailurePolicy::Error,
        AromaticityFailurePolicy::Error,
        false
    ))]
    #[case::nondefault(AromaticityResolveConfig::new(
        AromaticityConfig::default(),
        AromaticityFailurePolicy::Keep,
        AromaticityFailurePolicy::Keep,
        true
    ))]
    fn test_aromaticity_resolve_config_to_rust(#[case] config: AromaticityResolveConfig) {
        assert_eq!(config.to_rust(), config.0);
    }

    #[rstest]
    #[case::default(
        StereoResolveConfig::new(
            StereoFailurePolicy::Error,
            StereoFailurePolicy::Error,
            StereoFailurePolicy::Error,
            StereoFailurePolicy::Error,
            false
        ),
        GraphStereoResolveConfig::default()
    )]
    #[case::configured(
        StereoResolveConfig::new(StereoFailurePolicy::Keep, StereoFailurePolicy::Remove, StereoFailurePolicy::Remove, StereoFailurePolicy::Keep, true),
        GraphStereoResolveConfig {
            tetrahedral_stereo_failure: GraphStereoFailurePolicy::Keep,
            stereo_atom_failure: GraphStereoFailurePolicy::Remove,
            cis_trans_stereo_failure: GraphStereoFailurePolicy::Remove,
            stereo_bond_failure: GraphStereoFailurePolicy::Keep,
            reset_stereo_constraints: true,
        }
    )]
    fn test_stereo_resolve_config_new(
        #[case] config: StereoResolveConfig,
        #[case] expected: GraphStereoResolveConfig,
    ) {
        assert_eq!(config.0, expected);
    }

    #[rstest]
    #[case::default(
        StereoResolveConfig::new(StereoFailurePolicy::Error, StereoFailurePolicy::Error, StereoFailurePolicy::Error, StereoFailurePolicy::Error, false),
        "StereoResolveConfig(tetrahedral_stereo_failure=StereoFailurePolicy.Error, stereo_atom_failure=StereoFailurePolicy.Error, cis_trans_stereo_failure=StereoFailurePolicy.Error, stereo_bond_failure=StereoFailurePolicy.Error, reset_stereo_constraints=False)"
    )]
    #[case::configured(
        StereoResolveConfig::new(StereoFailurePolicy::Keep, StereoFailurePolicy::Remove, StereoFailurePolicy::Remove, StereoFailurePolicy::Keep, true),
        "StereoResolveConfig(tetrahedral_stereo_failure=StereoFailurePolicy.Keep, stereo_atom_failure=StereoFailurePolicy.Remove, cis_trans_stereo_failure=StereoFailurePolicy.Remove, stereo_bond_failure=StereoFailurePolicy.Keep, reset_stereo_constraints=True)"
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
        tetrahedral_stereo_failure: GraphStereoFailurePolicy::Keep,
        stereo_atom_failure: GraphStereoFailurePolicy::Remove,
        cis_trans_stereo_failure: GraphStereoFailurePolicy::Remove,
        stereo_bond_failure: GraphStereoFailurePolicy::Keep,
        reset_stereo_constraints: true,
    })]
    fn test_stereo_resolve_config_from_rust(#[case] config: GraphStereoResolveConfig) {
        assert_eq!(StereoResolveConfig::from_rust(config).0, config);
    }

    #[rstest]
    #[case::default(StereoResolveConfig::new(
        StereoFailurePolicy::Error,
        StereoFailurePolicy::Error,
        StereoFailurePolicy::Error,
        StereoFailurePolicy::Error,
        false
    ))]
    #[case::configured(StereoResolveConfig::new(
        StereoFailurePolicy::Keep,
        StereoFailurePolicy::Remove,
        StereoFailurePolicy::Remove,
        StereoFailurePolicy::Keep,
        true
    ))]
    fn test_stereo_resolve_config_to_rust(#[case] config: StereoResolveConfig) {
        assert_eq!(config.to_rust(), config.0);
    }

    #[rstest]
    #[case::default(
        AromaticityResolveConfig::new(
            AromaticityConfig::default(),
            AromaticityFailurePolicy::Error,
            AromaticityFailurePolicy::Error,
            false
        ),
        StereoResolveConfig::new(
            StereoFailurePolicy::Error,
            StereoFailurePolicy::Error,
            StereoFailurePolicy::Error,
            StereoFailurePolicy::Error,
            false
        ),
        GraphResolveConfig::default()
    )]
    #[case::aromaticity(
        AromaticityResolveConfig::new(AromaticityConfig::default(), AromaticityFailurePolicy::Keep, AromaticityFailurePolicy::Keep, true),
        StereoResolveConfig::new(StereoFailurePolicy::Error, StereoFailurePolicy::Error, StereoFailurePolicy::Error, StereoFailurePolicy::Error, false),
        GraphResolveConfig {
            isotope: GraphIsotopePolicy::Strict,
            aromaticity: GraphAromaticityResolveConfig {
                perception: Default::default(),
                aromatic_valence_failure: GraphAromaticityFailurePolicy::Keep,
                aromatic_system_failure: GraphAromaticityFailurePolicy::Keep,
                reset_aromatic_valence: true,
            },
            stereo: GraphStereoResolveConfig::default(),
        },
    )]
    #[case::stereo(
        AromaticityResolveConfig::new(AromaticityConfig::default(), AromaticityFailurePolicy::Error, AromaticityFailurePolicy::Error, false),
        StereoResolveConfig::new(StereoFailurePolicy::Keep, StereoFailurePolicy::Remove, StereoFailurePolicy::Remove, StereoFailurePolicy::Keep, true),
        GraphResolveConfig {
            isotope: GraphIsotopePolicy::Strict,
            aromaticity: GraphAromaticityResolveConfig::default(),
            stereo: GraphStereoResolveConfig {
                tetrahedral_stereo_failure: GraphStereoFailurePolicy::Keep,
                stereo_atom_failure: GraphStereoFailurePolicy::Remove,
                cis_trans_stereo_failure: GraphStereoFailurePolicy::Remove,
                stereo_bond_failure: GraphStereoFailurePolicy::Keep,
                reset_stereo_constraints: true,
            },
        },
    )]
    fn test_resolve_config_new(
        #[case] aromaticity: AromaticityResolveConfig,
        #[case] stereo: StereoResolveConfig,
        #[case] expected: GraphResolveConfig,
    ) {
        assert_eq!(
            ResolveConfig::new(aromaticity, stereo, IsotopePolicy::Strict).0,
            expected
        );
    }

    #[rstest]
    fn test_resolve_config_default() {
        assert_eq!(
            ResolveConfig::default(),
            ResolveConfig(GraphResolveConfig {
                isotope: GraphIsotopePolicy::Strict,
                aromaticity: GraphAromaticityResolveConfig::default(),
                stereo: GraphStereoResolveConfig::default(),
            })
        );
    }

    #[rstest]
    #[case::default(ResolveConfig::default(), "ResolveConfig.default()")]
    #[case::configured(
        ResolveConfig::new(
            AromaticityResolveConfig::new(AromaticityConfig::default(), AromaticityFailurePolicy::Keep, AromaticityFailurePolicy::Keep, true),
            StereoResolveConfig::new(StereoFailurePolicy::Keep, StereoFailurePolicy::Remove, StereoFailurePolicy::Remove, StereoFailurePolicy::Keep, true),
            IsotopePolicy::Strict,
        ),
        "ResolveConfig(isotope=IsotopePolicy.Strict, aromaticity=AromaticityResolveConfig(perception=AromaticityConfig(ring_config=RingConfig(simple_cycle_algorithm=SimpleCycleEnumerationAlgorithm.ReadTarjan(), relevant_cycle_algorithm=RelevantCycleEnumerationAlgorithm.Vismara()), connected_components_algorithm=ConnectedComponentsAlgorithm.Bfs(), maximum_independent_set_algorithm=MaximumIndependentSetAlgorithm.BranchAndBound()), aromatic_valence_failure=AromaticityFailurePolicy.Keep, aromatic_system_failure=AromaticityFailurePolicy.Keep, reset_aromatic_valence=True), stereo=StereoResolveConfig(tetrahedral_stereo_failure=StereoFailurePolicy.Keep, stereo_atom_failure=StereoFailurePolicy.Remove, cis_trans_stereo_failure=StereoFailurePolicy.Remove, stereo_bond_failure=StereoFailurePolicy.Keep, reset_stereo_constraints=True))",
    )]
    fn test_resolve_config_repr(#[case] config: ResolveConfig, #[case] expected: &str) {
        assert_eq!(config.__repr__(), expected);
    }

    #[rstest]
    #[case::default(GraphResolveConfig::default())]
    #[case::configured(GraphResolveConfig {
        isotope: GraphIsotopePolicy::Strict,
        aromaticity: GraphAromaticityResolveConfig {
            perception: Default::default(),
            aromatic_valence_failure: GraphAromaticityFailurePolicy::Keep,
            aromatic_system_failure: GraphAromaticityFailurePolicy::Keep,
            reset_aromatic_valence: true,
        },
        stereo: GraphStereoResolveConfig {
            tetrahedral_stereo_failure: GraphStereoFailurePolicy::Keep,
            stereo_atom_failure: GraphStereoFailurePolicy::Remove,
            cis_trans_stereo_failure: GraphStereoFailurePolicy::Remove,
            stereo_bond_failure: GraphStereoFailurePolicy::Keep,
            reset_stereo_constraints: true,
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
            AromaticityFailurePolicy::Keep,
            AromaticityFailurePolicy::Keep,
            true
        ),
        StereoResolveConfig::new(
            StereoFailurePolicy::Keep,
            StereoFailurePolicy::Remove,
            StereoFailurePolicy::Remove,
            StereoFailurePolicy::Keep,
            true
        ),
        IsotopePolicy::Strict
    ))]
    fn test_resolve_config_to_rust(#[case] config: ResolveConfig) {
        assert_eq!(config.to_rust(), config.0);
    }

    #[rstest]
    #[case::hmo(
        GraphResolveContradiction::Aromaticity(GraphAromaticityContradiction::HmoInvalidInput(
            String::from("odd component"),
        )),
        "hmo: invalid input: odd component"
    )]
    fn test_resolve_contradiction_str(
        #[case] contradiction: GraphResolveContradiction,
        #[case] expected: &str,
    ) {
        assert_eq!(
            ResolveContradiction::from_rust(contradiction).__str__(),
            expected
        );
    }

    #[rstest]
    #[case::hmo(
        GraphResolveContradiction::Aromaticity(GraphAromaticityContradiction::HmoInvalidInput(
            String::from("odd component"),
        )),
        "ResolveContradiction(\"hmo: invalid input: odd component\")"
    )]
    fn test_resolve_contradiction_repr(
        #[case] contradiction: GraphResolveContradiction,
        #[case] expected: &str,
    ) {
        assert_eq!(
            ResolveContradiction::from_rust(contradiction).__repr__(),
            expected
        );
    }

    #[rstest]
    fn test_resolve_contradiction_eq() {
        let contradiction = || {
            ResolveContradiction::from_rust(GraphResolveContradiction::Aromaticity(
                GraphAromaticityContradiction::HmoInvalidInput(String::from("odd component")),
            ))
        };
        let different = ResolveContradiction::from_rust(GraphResolveContradiction::Aromaticity(
            GraphAromaticityContradiction::HmoInvalidInput(String::from("other")),
        ));

        assert_eq!(contradiction(), contradiction());
        assert_ne!(contradiction(), different);
    }

    #[rstest]
    #[case::determined(
        Solution::Determined { molecule: None },
        Solution::Determined { molecule: None },
        true
    )]
    #[case::underdetermined(
        Solution::Underdetermined {}, Solution::Underdetermined {}, true
    )]
    #[case::different(
        Solution::Determined { molecule: None }, Solution::Underdetermined {}, false
    )]
    #[case::contradictory(
        Solution::Contradictory {
            contradiction: ResolveContradiction::from_rust(GraphResolveContradiction::Aromaticity(
                GraphAromaticityContradiction::HmoInvalidInput(String::from("odd component")),
            )),
        },
        Solution::Contradictory {
            contradiction: ResolveContradiction::from_rust(GraphResolveContradiction::Aromaticity(
                GraphAromaticityContradiction::HmoInvalidInput(String::from("odd component")),
            )),
        },
        true
    )]
    fn test_solution_eq(#[case] left: Solution, #[case] right: Solution, #[case] expected: bool) {
        Python::attach(|py| {
            assert_eq!(left.__eq__(&right, py).unwrap(), expected);
        });
    }

    #[rstest]
    #[case::determined(
        Solution::Determined { molecule: None }, "Solution.Determined(molecule=None)"
    )]
    #[case::underdetermined(Solution::Underdetermined {}, "Solution.Underdetermined()")]
    #[case::contradictory(
        Solution::Contradictory {
            contradiction: ResolveContradiction::from_rust(GraphResolveContradiction::Aromaticity(
                GraphAromaticityContradiction::HmoInvalidInput(String::from("odd component")),
            )),
        },
        "Solution.Contradictory(contradiction=ResolveContradiction(\"hmo: invalid input: odd component\"))"
    )]
    fn test_solution_repr(#[case] solution: Solution, #[case] expected: &str) {
        Python::attach(|py| {
            assert_eq!(solution.__repr__(py).unwrap(), expected);
        });
    }
}

/// Read-only per-atom candidate sets from a resolution run: each entry maps
/// an atom index to its surviving completions.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AtomCompletions(GraphAtomCompletions);

#[pymethods]
impl AtomCompletions {
    fn __len__(&self) -> usize {
        self.0.len()
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// The surviving completions of an atom, or `None`.
    pub(crate) fn get(&self, atom: u32) -> Option<Vec<AtomForm>> {
        self.0
            .get(GraphIrAtomId(atom))
            .map(|forms| forms.iter().cloned().map(AtomForm::from_rust).collect())
    }

    /// Every entry, in ascending atom order.
    fn items(&self) -> Vec<(u32, Vec<AtomForm>)> {
        self.0
            .iter()
            .map(|(atom, forms)| {
                (
                    atom.0,
                    forms.iter().cloned().map(AtomForm::from_rust).collect(),
                )
            })
            .collect()
    }

    pub(crate) fn __repr__(&self) -> String {
        let entries = self
            .0
            .iter()
            .map(|(atom, forms)| {
                let rendered = forms
                    .iter()
                    .map(|form| format!("{:?}", form.to_string()))
                    .collect::<Vec<_>>()
                    .join(", ");
                format!("{}: [{rendered}]", atom.0)
            })
            .collect::<Vec<_>>()
            .join(", ");
        format!("AtomCompletions({{{entries}}})")
    }
}

impl AtomCompletions {
    pub(crate) fn from_rust(completions: &GraphAtomCompletions) -> Self {
        Self(completions.clone())
    }
}

/// Read-only resolution report: the plural survivors and the recorded
/// tie-break uses.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ResolveReport(GraphResolveReport);

#[pymethods]
impl ResolveReport {
    #[getter]
    pub(crate) fn unresolved(&self) -> AtomCompletions {
        AtomCompletions::from_rust(&self.0.unresolved)
    }

    #[getter]
    pub(crate) fn tie_breaks(&self) -> Vec<u32> {
        self.0.tie_breaks.iter().map(|atom| atom.0).collect()
    }

    pub(crate) fn __repr__(&self) -> String {
        format!(
            "ResolveReport(unresolved={}, tie_breaks={:?})",
            self.unresolved().__repr__(),
            self.tie_breaks(),
        )
    }
}

impl ResolveReport {
    pub(crate) fn from_rust(report: GraphResolveReport) -> Self {
        Self(report)
    }
}

/// Resolution contradiction: the chemistry model rejected the input.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ResolveContradiction(GraphResolveContradiction);

#[pymethods]
impl ResolveContradiction {
    pub(crate) fn __str__(&self) -> String {
        self.0.to_string()
    }

    pub(crate) fn __repr__(&self) -> String {
        format!("ResolveContradiction({:?})", self.0.to_string())
    }
}

impl ResolveContradiction {
    pub(crate) fn from_rust(contradiction: GraphResolveContradiction) -> Self {
        Self(contradiction)
    }
}

/// The resolution outcome. Reports are returned separately by with_report methods.
#[pyclass(frozen)]
#[derive(Debug)]
pub enum Solution {
    /// Resolution completed. Consuming resolution returns a molecule; in-place resolution
    /// returns None. The molecule getter retains the same Python object on every access.
    #[pyo3(constructor = (*, molecule=None))]
    Determined { molecule: Option<Py<Molecule>> },
    /// Resolution could not determine the molecule completely.
    Underdetermined {},
    /// The chemistry model rejected the input.
    #[pyo3(constructor = (*, contradiction))]
    Contradictory { contradiction: ResolveContradiction },
}

#[pymethods]
impl Solution {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> PyResult<bool> {
        Ok(match (self, other) {
            (Self::Determined { molecule: left }, Self::Determined { molecule: right }) => {
                let left = left
                    .as_ref()
                    .map(|molecule| molecule.try_borrow(py))
                    .transpose()?;
                let right = right
                    .as_ref()
                    .map(|molecule| molecule.try_borrow(py))
                    .transpose()?;
                left.as_ref()
                    .map(|molecule| molecule.to_rust())
                    .transpose()?
                    == right
                        .as_ref()
                        .map(|molecule| molecule.to_rust())
                        .transpose()?
            }
            (Self::Underdetermined {}, Self::Underdetermined {}) => true,
            (
                Self::Contradictory {
                    contradiction: left,
                },
                Self::Contradictory {
                    contradiction: right,
                },
            ) => left == right,
            _ => false,
        })
    }

    pub(crate) fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(match self {
            Self::Determined { molecule } => format!(
                "Solution.Determined(molecule={})",
                match molecule {
                    Some(molecule) => molecule.try_borrow(py)?.__repr__()?,
                    None => String::from("None"),
                },
            ),
            Self::Underdetermined {} => String::from("Solution.Underdetermined()"),
            Self::Contradictory { contradiction } => format!(
                "Solution.Contradictory(contradiction={})",
                contradiction.__repr__(),
            ),
        })
    }
}
