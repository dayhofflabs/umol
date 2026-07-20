//! Python bindings for fingerprint configuration values.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use umol_graph::fingerprint::{
    EcfpFeaturizer as GraphEcfpFeaturizer, Featurizer as GraphFeaturizer,
    MorganFeaturizer as GraphMorganFeaturizer, PatternFingerprinter as GraphPatternFingerprinter,
    ReactionCombinator as GraphReactionCombinator,
    SubstructureFeaturizer as GraphSubstructureFeaturizer, WlFeaturizer as GraphWlFeaturizer,
};
use umol_graph::hash::{EcfpHashScheme as GraphEcfpHashScheme, WlHashScheme as GraphWlHashScheme};
use umol_graph_core::RefinementRounds as GraphCoreRefinementRounds;

/// Number of graph-refinement rounds: fixed or until stabilization.
#[pyclass(from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RefinementRounds {
    Fixed { rounds: u32 },
    ToFixpoint(),
}

#[pymethods]
impl RefinementRounds {
    fn __eq__(&self, other: &Self) -> bool {
        self == other
    }

    fn __repr__(&self) -> String {
        match self {
            Self::Fixed { rounds } => format!("RefinementRounds.Fixed(rounds={rounds})"),
            Self::ToFixpoint() => "RefinementRounds.ToFixpoint()".to_owned(),
        }
    }
}

impl RefinementRounds {
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "boundary conversion is part of the binding contract without a production caller"
        )
    )]
    pub(crate) fn from_rust(rounds: GraphCoreRefinementRounds) -> Self {
        match rounds {
            GraphCoreRefinementRounds::Fixed(rounds) => Self::Fixed { rounds },
            GraphCoreRefinementRounds::ToFixpoint => Self::ToFixpoint(),
        }
    }

    pub(crate) fn to_rust(self) -> GraphCoreRefinementRounds {
        match self {
            Self::Fixed { rounds } => GraphCoreRefinementRounds::Fixed(rounds),
            Self::ToFixpoint() => GraphCoreRefinementRounds::ToFixpoint,
        }
    }
}

/// Frozen hashing recipe for Weisfeiler-Lehman refinement fingerprints.
#[pyclass(from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WlHashScheme {
    Xxh3SortedWidth64V1(),
}

#[pymethods]
impl WlHashScheme {
    fn __eq__(&self, other: &Self) -> bool {
        self == other
    }

    fn __repr__(&self) -> String {
        match self {
            Self::Xxh3SortedWidth64V1() => "WlHashScheme.Xxh3SortedWidth64V1()".to_owned(),
        }
    }
}

impl WlHashScheme {
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "boundary conversion is part of the binding contract without a production caller"
        )
    )]
    pub(crate) fn from_rust(scheme: GraphWlHashScheme) -> Self {
        match scheme {
            GraphWlHashScheme::Xxh3SortedWidth64V1 => Self::Xxh3SortedWidth64V1(),
        }
    }

    pub(crate) fn to_rust(self) -> GraphWlHashScheme {
        match self {
            Self::Xxh3SortedWidth64V1() => GraphWlHashScheme::Xxh3SortedWidth64V1,
        }
    }
}

/// Frozen hashing recipe for extended-connectivity fingerprints.
#[pyclass(from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EcfpHashScheme {
    Xxh3Width64V1(),
}

#[pymethods]
impl EcfpHashScheme {
    fn __eq__(&self, other: &Self) -> bool {
        self == other
    }

    fn __repr__(&self) -> String {
        match self {
            Self::Xxh3Width64V1() => "EcfpHashScheme.Xxh3Width64V1()".to_owned(),
        }
    }
}

impl EcfpHashScheme {
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "boundary conversion is part of the binding contract without a production caller"
        )
    )]
    pub(crate) fn from_rust(scheme: GraphEcfpHashScheme) -> Self {
        match scheme {
            GraphEcfpHashScheme::Xxh3Width64V1 => Self::Xxh3Width64V1(),
        }
    }

    pub(crate) fn to_rust(self) -> GraphEcfpHashScheme {
        match self {
            Self::Xxh3Width64V1() => GraphEcfpHashScheme::Xxh3Width64V1,
        }
    }
}

/// Configuration for hashed molecular fingerprint generation.
#[pyclass(from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HashedFingerprintConfig {
    #[pyo3(constructor = (*, radius=2))]
    Morgan { radius: u32 },
    #[pyo3(constructor = (*, radius=2, scheme=EcfpHashScheme::Xxh3Width64V1()))]
    Ecfp { radius: u32, scheme: EcfpHashScheme },
    #[pyo3(constructor = (*, rounds, scheme=WlHashScheme::Xxh3SortedWidth64V1()))]
    Wl {
        rounds: RefinementRounds,
        scheme: WlHashScheme,
    },
}

#[pymethods]
impl HashedFingerprintConfig {
    fn __eq__(&self, other: &Self) -> bool {
        self == other
    }

    fn __repr__(&self) -> String {
        match self {
            Self::Morgan { radius } => {
                format!("HashedFingerprintConfig.Morgan(radius={radius})")
            }
            Self::Ecfp { radius, scheme } => format!(
                "HashedFingerprintConfig.Ecfp(radius={radius}, scheme={})",
                scheme.__repr__()
            ),
            Self::Wl { rounds, scheme } => format!(
                "HashedFingerprintConfig.Wl(rounds={}, scheme={})",
                rounds.__repr__(),
                scheme.__repr__()
            ),
        }
    }
}

impl HashedFingerprintConfig {
    pub(crate) fn to_rust(self) -> GraphFeaturizer {
        match self {
            Self::Morgan { radius } => GraphFeaturizer::Morgan(GraphMorganFeaturizer::new(radius)),
            Self::Ecfp { radius, scheme } => GraphFeaturizer::Ecfp(GraphEcfpFeaturizer {
                radius,
                scheme: scheme.to_rust(),
            }),
            Self::Wl { rounds, scheme } => GraphFeaturizer::Wl(GraphWlFeaturizer {
                rounds: rounds.to_rust(),
                scheme: scheme.to_rust(),
            }),
        }
    }
}

/// Configuration for the fixed-width pattern fingerprint.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PatternFingerprintConfig {
    width: usize,
}

#[pymethods]
impl PatternFingerprintConfig {
    #[new]
    #[pyo3(signature = (*, width=2048))]
    fn new(width: isize) -> PyResult<Self> {
        let width = usize::try_from(width)
            .ok()
            .filter(|width| *width > 0)
            .ok_or_else(|| PyValueError::new_err("width must be positive"))?;
        Ok(Self { width })
    }

    #[getter]
    fn width(&self) -> usize {
        self.width
    }

    fn __repr__(&self) -> String {
        format!("PatternFingerprintConfig(width={})", self.width)
    }
}

impl PatternFingerprintConfig {
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "boundary conversion is part of the binding contract without a production caller"
        )
    )]
    pub(crate) fn from_rust(fingerprinter: GraphPatternFingerprinter) -> Self {
        Self {
            width: fingerprinter.width,
        }
    }

    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "boundary conversion is part of the binding contract without a production caller"
        )
    )]
    pub(crate) fn to_rust(self) -> GraphPatternFingerprinter {
        GraphPatternFingerprinter { width: self.width }
    }
}

/// Configuration for exact structural-feature generation.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct StructuralFingerprintConfig {
    max_bonds: u32,
}

#[pymethods]
impl StructuralFingerprintConfig {
    #[new]
    #[pyo3(signature = (*, max_bonds))]
    fn new(max_bonds: u32) -> Self {
        Self { max_bonds }
    }

    #[getter]
    fn max_bonds(&self) -> u32 {
        self.max_bonds
    }

    fn __repr__(&self) -> String {
        format!("StructuralFingerprintConfig(max_bonds={})", self.max_bonds)
    }
}

impl StructuralFingerprintConfig {
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "boundary conversion is part of the binding contract without a production caller"
        )
    )]
    pub(crate) fn from_rust(featurizer: GraphSubstructureFeaturizer) -> Self {
        Self {
            max_bonds: featurizer.max_bonds,
        }
    }

    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "boundary conversion is part of the binding contract without a production caller"
        )
    )]
    pub(crate) fn to_rust(self) -> GraphSubstructureFeaturizer {
        GraphSubstructureFeaturizer::new(self.max_bonds)
    }
}

/// Configuration for combining molecular features across a reaction.
#[pyclass(from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReactionCombinedFingerprintConfig {
    #[pyo3(constructor = (*, molecule))]
    Difference { molecule: HashedFingerprintConfig },
    #[pyo3(constructor = (*, molecule))]
    DisjointUnion { molecule: HashedFingerprintConfig },
}

#[pymethods]
impl ReactionCombinedFingerprintConfig {
    fn __eq__(&self, other: &Self) -> bool {
        self == other
    }

    fn __repr__(&self) -> String {
        match self {
            Self::Difference { molecule } => format!(
                "ReactionCombinedFingerprintConfig.Difference(molecule={})",
                molecule.__repr__()
            ),
            Self::DisjointUnion { molecule } => format!(
                "ReactionCombinedFingerprintConfig.DisjointUnion(molecule={})",
                molecule.__repr__()
            ),
        }
    }
}

impl ReactionCombinedFingerprintConfig {
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "boundary conversion is part of the binding contract without a production caller"
        )
    )]
    pub(crate) fn to_rust(self) -> (GraphFeaturizer, GraphReactionCombinator) {
        match self {
            Self::Difference { molecule } => {
                (molecule.to_rust(), GraphReactionCombinator::Difference)
            }
            Self::DisjointUnion { molecule } => {
                (molecule.to_rust(), GraphReactionCombinator::DisjointUnion)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;
    use crate::convert::into_py_variant;

    #[rstest]
    #[case::zero(
        GraphCoreRefinementRounds::Fixed(0),
        RefinementRounds::Fixed { rounds: 0 }
    )]
    #[case::fixed(
        GraphCoreRefinementRounds::Fixed(3),
        RefinementRounds::Fixed { rounds: 3 }
    )]
    #[case::fixpoint(GraphCoreRefinementRounds::ToFixpoint, RefinementRounds::ToFixpoint())]
    fn test_refinement_rounds_from_rust(
        #[case] rounds: GraphCoreRefinementRounds,
        #[case] expected: RefinementRounds,
    ) {
        assert_eq!(RefinementRounds::from_rust(rounds), expected);
    }

    #[rstest]
    #[case::zero(
        RefinementRounds::Fixed { rounds: 0 },
        GraphCoreRefinementRounds::Fixed(0)
    )]
    #[case::fixed(
        RefinementRounds::Fixed { rounds: 3 },
        GraphCoreRefinementRounds::Fixed(3)
    )]
    #[case::fixpoint(RefinementRounds::ToFixpoint(), GraphCoreRefinementRounds::ToFixpoint)]
    fn test_refinement_rounds_to_rust(
        #[case] rounds: RefinementRounds,
        #[case] expected: GraphCoreRefinementRounds,
    ) {
        assert_eq!(rounds.to_rust(), expected);
    }

    #[rstest]
    #[case::zero(
        RefinementRounds::Fixed { rounds: 0 },
        "RefinementRounds.Fixed(rounds=0)",
        Some(0)
    )]
    #[case::fixed(
        RefinementRounds::Fixed { rounds: 3 },
        "RefinementRounds.Fixed(rounds=3)",
        Some(3)
    )]
    #[case::fixpoint(RefinementRounds::ToFixpoint(), "RefinementRounds.ToFixpoint()", None)]
    fn test_refinement_rounds_value(
        #[case] rounds: RefinementRounds,
        #[case] expected_repr: &str,
        #[case] expected_rounds: Option<u32>,
    ) {
        Python::attach(|py| {
            let expected =
                into_py_variant(py, RefinementRounds::from_rust(rounds.to_rust())).unwrap();
            let rounds = into_py_variant(py, rounds).unwrap();
            let expected = expected.bind(py).as_any();
            let rounds = rounds.bind(py).as_any();

            assert!(rounds.eq(expected).unwrap());
            assert_eq!(
                rounds.repr().unwrap().extract::<String>().unwrap(),
                expected_repr
            );
            match expected_rounds {
                Some(expected_rounds) => assert_eq!(
                    rounds.getattr("rounds").unwrap().extract::<u32>().unwrap(),
                    expected_rounds
                ),
                None => assert!(rounds.getattr("rounds").is_err()),
            }
        });
    }

    #[rstest]
    #[case::xxh3_sorted_width64_v1(
        GraphWlHashScheme::Xxh3SortedWidth64V1,
        WlHashScheme::Xxh3SortedWidth64V1()
    )]
    fn test_wl_hash_scheme_from_rust(
        #[case] scheme: GraphWlHashScheme,
        #[case] expected: WlHashScheme,
    ) {
        assert_eq!(WlHashScheme::from_rust(scheme), expected);
    }

    #[rstest]
    #[case::xxh3_sorted_width64_v1(
        WlHashScheme::Xxh3SortedWidth64V1(),
        GraphWlHashScheme::Xxh3SortedWidth64V1,
        1,
        64
    )]
    fn test_wl_hash_scheme_to_rust(
        #[case] scheme: WlHashScheme,
        #[case] expected: GraphWlHashScheme,
        #[case] expected_version: u16,
        #[case] expected_width: u16,
    ) {
        let scheme = scheme.to_rust();

        assert_eq!(scheme, expected);
        assert_eq!(scheme.version(), expected_version);
        assert_eq!(scheme.identifier_width(), expected_width);
    }

    #[rstest]
    #[case::xxh3_sorted_width64_v1(
        WlHashScheme::Xxh3SortedWidth64V1(),
        "WlHashScheme.Xxh3SortedWidth64V1()"
    )]
    fn test_wl_hash_scheme_value(#[case] scheme: WlHashScheme, #[case] expected_repr: &str) {
        Python::attach(|py| {
            let expected = into_py_variant(
                py,
                WlHashScheme::from_rust(GraphWlHashScheme::Xxh3SortedWidth64V1),
            )
            .unwrap();
            let scheme = into_py_variant(py, scheme).unwrap();
            let expected = expected.bind(py).as_any();
            let scheme = scheme.bind(py).as_any();

            assert!(scheme.eq(expected).unwrap());
            assert_eq!(
                scheme.repr().unwrap().extract::<String>().unwrap(),
                expected_repr
            );
        });
    }

    #[rstest]
    #[case::xxh3_width64_v1(GraphEcfpHashScheme::Xxh3Width64V1, EcfpHashScheme::Xxh3Width64V1())]
    fn test_ecfp_hash_scheme_from_rust(
        #[case] scheme: GraphEcfpHashScheme,
        #[case] expected: EcfpHashScheme,
    ) {
        assert_eq!(EcfpHashScheme::from_rust(scheme), expected);
    }

    #[rstest]
    #[case::xxh3_width64_v1(
        EcfpHashScheme::Xxh3Width64V1(),
        GraphEcfpHashScheme::Xxh3Width64V1,
        1,
        64
    )]
    fn test_ecfp_hash_scheme_to_rust(
        #[case] scheme: EcfpHashScheme,
        #[case] expected: GraphEcfpHashScheme,
        #[case] expected_version: u16,
        #[case] expected_width: u16,
    ) {
        let scheme = scheme.to_rust();

        assert_eq!(scheme, expected);
        assert_eq!(scheme.version(), expected_version);
        assert_eq!(scheme.identifier_width(), expected_width);
    }

    #[rstest]
    #[case::xxh3_width64_v1(EcfpHashScheme::Xxh3Width64V1(), "EcfpHashScheme.Xxh3Width64V1()")]
    fn test_ecfp_hash_scheme_value(#[case] scheme: EcfpHashScheme, #[case] expected_repr: &str) {
        Python::attach(|py| {
            let expected = into_py_variant(
                py,
                EcfpHashScheme::from_rust(GraphEcfpHashScheme::Xxh3Width64V1),
            )
            .unwrap();
            let scheme = into_py_variant(py, scheme).unwrap();
            let expected = expected.bind(py).as_any();
            let scheme = scheme.bind(py).as_any();

            assert!(scheme.eq(expected).unwrap());
            assert_eq!(
                scheme.repr().unwrap().extract::<String>().unwrap(),
                expected_repr
            );
        });
    }

    #[rstest]
    #[case::morgan_default(HashedFingerprintConfig::Morgan { radius: 2 })]
    #[case::morgan_explicit(HashedFingerprintConfig::Morgan { radius: 3 })]
    #[case::ecfp_default(HashedFingerprintConfig::Ecfp {
        radius: 2,
        scheme: EcfpHashScheme::Xxh3Width64V1(),
    })]
    #[case::ecfp_explicit(HashedFingerprintConfig::Ecfp {
        radius: 3,
        scheme: EcfpHashScheme::Xxh3Width64V1(),
    })]
    #[case::wl_fixed(HashedFingerprintConfig::Wl {
        rounds: RefinementRounds::Fixed { rounds: 3 },
        scheme: WlHashScheme::Xxh3SortedWidth64V1(),
    })]
    #[case::wl_fixpoint(HashedFingerprintConfig::Wl {
        rounds: RefinementRounds::ToFixpoint(),
        scheme: WlHashScheme::Xxh3SortedWidth64V1(),
    })]
    fn test_hashed_fingerprint_config_to_rust(#[case] config: HashedFingerprintConfig) {
        let featurizer = config.to_rust();

        match (config, featurizer) {
            (HashedFingerprintConfig::Morgan { radius }, GraphFeaturizer::Morgan(featurizer)) => {
                assert_eq!(featurizer.radius, radius)
            }
            (
                HashedFingerprintConfig::Ecfp { radius, scheme },
                GraphFeaturizer::Ecfp(featurizer),
            ) => {
                assert_eq!(featurizer.radius, radius);
                assert_eq!(featurizer.scheme, scheme.to_rust());
            }
            (HashedFingerprintConfig::Wl { rounds, scheme }, GraphFeaturizer::Wl(featurizer)) => {
                assert_eq!(featurizer.rounds, rounds.to_rust());
                assert_eq!(featurizer.scheme, scheme.to_rust());
            }
            (config, featurizer) => {
                panic!("config {config:?} lowered to mismatched featurizer {featurizer:?}")
            }
        }
    }

    #[rstest]
    #[case::morgan(
        HashedFingerprintConfig::Morgan { radius: 3 },
        "HashedFingerprintConfig.Morgan(radius=3)"
    )]
    #[case::ecfp(
        HashedFingerprintConfig::Ecfp {
            radius: 3,
            scheme: EcfpHashScheme::Xxh3Width64V1(),
        },
        "HashedFingerprintConfig.Ecfp(radius=3, scheme=EcfpHashScheme.Xxh3Width64V1())"
    )]
    #[case::wl(
        HashedFingerprintConfig::Wl {
            rounds: RefinementRounds::Fixed { rounds: 3 },
            scheme: WlHashScheme::Xxh3SortedWidth64V1(),
        },
        "HashedFingerprintConfig.Wl(rounds=RefinementRounds.Fixed(rounds=3), scheme=WlHashScheme.Xxh3SortedWidth64V1())"
    )]
    fn test_hashed_fingerprint_config_value(
        #[case] config: HashedFingerprintConfig,
        #[case] expected_repr: &str,
    ) {
        Python::attach(|py| {
            let expected = into_py_variant(py, config).unwrap();
            let config = into_py_variant(py, config).unwrap();
            let expected = expected.bind(py).as_any();
            let config = config.bind(py).as_any();

            assert!(config.eq(expected).unwrap());
            assert_eq!(
                config.repr().unwrap().extract::<String>().unwrap(),
                expected_repr
            );
        });
    }

    #[rstest]
    #[case::default(2048, PatternFingerprintConfig { width: 2048 })]
    #[case::custom(512, PatternFingerprintConfig { width: 512 })]
    fn test_pattern_fingerprint_config_new(
        #[case] width: isize,
        #[case] expected: PatternFingerprintConfig,
    ) {
        assert_eq!(PatternFingerprintConfig::new(width).unwrap(), expected);
    }

    #[rstest]
    #[case::zero(0)]
    #[case::negative(-1)]
    fn test_pattern_fingerprint_config_new_error(#[case] width: isize) {
        Python::attach(|py| {
            let error = PatternFingerprintConfig::new(width).unwrap_err();

            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "width must be positive"
            );
        });
    }

    #[rstest]
    #[case::default(
        GraphPatternFingerprinter::new(),
        PatternFingerprintConfig { width: 2048 }
    )]
    #[case::custom(
        GraphPatternFingerprinter { width: 512 },
        PatternFingerprintConfig { width: 512 }
    )]
    fn test_pattern_fingerprint_config_from_rust(
        #[case] fingerprinter: GraphPatternFingerprinter,
        #[case] expected: PatternFingerprintConfig,
    ) {
        assert_eq!(PatternFingerprintConfig::from_rust(fingerprinter), expected);
    }

    #[rstest]
    #[case::default(PatternFingerprintConfig { width: 2048 }, 2048)]
    #[case::custom(PatternFingerprintConfig { width: 512 }, 512)]
    fn test_pattern_fingerprint_config_to_rust(
        #[case] config: PatternFingerprintConfig,
        #[case] expected_width: usize,
    ) {
        assert_eq!(config.to_rust().width, expected_width);
    }

    #[rstest]
    #[case::default(
        PatternFingerprintConfig { width: 2048 },
        2048,
        "PatternFingerprintConfig(width=2048)"
    )]
    #[case::custom(
        PatternFingerprintConfig { width: 512 },
        512,
        "PatternFingerprintConfig(width=512)"
    )]
    fn test_pattern_fingerprint_config_value(
        #[case] config: PatternFingerprintConfig,
        #[case] expected_width: usize,
        #[case] expected_repr: &str,
    ) {
        Python::attach(|py| {
            let expected = into_py_variant(py, config).unwrap();
            let config = into_py_variant(py, config).unwrap();
            let expected = expected.bind(py).as_any();
            let config = config.bind(py).as_any();

            assert!(config.eq(expected).unwrap());
            assert_eq!(
                config.getattr("width").unwrap().extract::<usize>().unwrap(),
                expected_width
            );
            assert_eq!(
                config.repr().unwrap().extract::<String>().unwrap(),
                expected_repr
            );
        });
    }

    #[rstest]
    #[case::zero(
        GraphSubstructureFeaturizer::new(0),
        StructuralFingerprintConfig { max_bonds: 0 }
    )]
    #[case::positive(
        GraphSubstructureFeaturizer::new(3),
        StructuralFingerprintConfig { max_bonds: 3 }
    )]
    fn test_structural_fingerprint_config_from_rust(
        #[case] featurizer: GraphSubstructureFeaturizer,
        #[case] expected: StructuralFingerprintConfig,
    ) {
        assert_eq!(StructuralFingerprintConfig::from_rust(featurizer), expected);
    }

    #[rstest]
    #[case::zero(StructuralFingerprintConfig { max_bonds: 0 }, 0)]
    #[case::positive(StructuralFingerprintConfig { max_bonds: 3 }, 3)]
    fn test_structural_fingerprint_config_to_rust(
        #[case] config: StructuralFingerprintConfig,
        #[case] expected_max_bonds: u32,
    ) {
        assert_eq!(config.to_rust().max_bonds, expected_max_bonds);
    }

    #[rstest]
    #[case::zero(
        StructuralFingerprintConfig { max_bonds: 0 },
        0,
        "StructuralFingerprintConfig(max_bonds=0)"
    )]
    #[case::positive(
        StructuralFingerprintConfig { max_bonds: 3 },
        3,
        "StructuralFingerprintConfig(max_bonds=3)"
    )]
    fn test_structural_fingerprint_config_value(
        #[case] config: StructuralFingerprintConfig,
        #[case] expected_max_bonds: u32,
        #[case] expected_repr: &str,
    ) {
        Python::attach(|py| {
            let expected = into_py_variant(py, config).unwrap();
            let config = into_py_variant(py, config).unwrap();
            let expected = expected.bind(py).as_any();
            let config = config.bind(py).as_any();

            assert!(config.eq(expected).unwrap());
            assert_eq!(
                config
                    .getattr("max_bonds")
                    .unwrap()
                    .extract::<u32>()
                    .unwrap(),
                expected_max_bonds
            );
            assert_eq!(
                config.repr().unwrap().extract::<String>().unwrap(),
                expected_repr
            );
        });
    }

    #[rstest]
    #[case::difference_morgan(
        ReactionCombinedFingerprintConfig::Difference {
            molecule: HashedFingerprintConfig::Morgan { radius: 2 },
        },
        GraphReactionCombinator::Difference
    )]
    #[case::difference_ecfp(
        ReactionCombinedFingerprintConfig::Difference {
            molecule: HashedFingerprintConfig::Ecfp {
                radius: 2,
                scheme: EcfpHashScheme::Xxh3Width64V1(),
            },
        },
        GraphReactionCombinator::Difference
    )]
    #[case::difference_wl(
        ReactionCombinedFingerprintConfig::Difference {
            molecule: HashedFingerprintConfig::Wl {
                rounds: RefinementRounds::Fixed { rounds: 3 },
                scheme: WlHashScheme::Xxh3SortedWidth64V1(),
            },
        },
        GraphReactionCombinator::Difference
    )]
    #[case::disjoint_union_morgan(
        ReactionCombinedFingerprintConfig::DisjointUnion {
            molecule: HashedFingerprintConfig::Morgan { radius: 2 },
        },
        GraphReactionCombinator::DisjointUnion
    )]
    #[case::disjoint_union_ecfp(
        ReactionCombinedFingerprintConfig::DisjointUnion {
            molecule: HashedFingerprintConfig::Ecfp {
                radius: 2,
                scheme: EcfpHashScheme::Xxh3Width64V1(),
            },
        },
        GraphReactionCombinator::DisjointUnion
    )]
    #[case::disjoint_union_wl(
        ReactionCombinedFingerprintConfig::DisjointUnion {
            molecule: HashedFingerprintConfig::Wl {
                rounds: RefinementRounds::Fixed { rounds: 3 },
                scheme: WlHashScheme::Xxh3SortedWidth64V1(),
            },
        },
        GraphReactionCombinator::DisjointUnion
    )]
    fn test_reaction_combined_fingerprint_config_to_rust(
        #[case] config: ReactionCombinedFingerprintConfig,
        #[case] expected_combinator: GraphReactionCombinator,
    ) {
        let molecule = match config {
            ReactionCombinedFingerprintConfig::Difference { molecule }
            | ReactionCombinedFingerprintConfig::DisjointUnion { molecule } => molecule,
        };
        let (featurizer, combinator) = config.to_rust();

        assert_eq!(combinator, expected_combinator);
        match (molecule, featurizer) {
            (HashedFingerprintConfig::Morgan { radius }, GraphFeaturizer::Morgan(featurizer)) => {
                assert_eq!(featurizer.radius, radius)
            }
            (
                HashedFingerprintConfig::Ecfp { radius, scheme },
                GraphFeaturizer::Ecfp(featurizer),
            ) => {
                assert_eq!(featurizer.radius, radius);
                assert_eq!(featurizer.scheme, scheme.to_rust());
            }
            (HashedFingerprintConfig::Wl { rounds, scheme }, GraphFeaturizer::Wl(featurizer)) => {
                assert_eq!(featurizer.rounds, rounds.to_rust());
                assert_eq!(featurizer.scheme, scheme.to_rust());
            }
            (molecule, featurizer) => {
                panic!("config {molecule:?} lowered to mismatched featurizer {featurizer:?}")
            }
        }
    }
}
