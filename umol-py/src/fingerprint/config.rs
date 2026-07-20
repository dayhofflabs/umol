//! Python bindings for fingerprint configuration values.

use pyo3::prelude::*;
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

    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "boundary conversion is part of the binding contract without a production caller"
        )
    )]
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

    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "boundary conversion is part of the binding contract without a production caller"
        )
    )]
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

    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "boundary conversion is part of the binding contract without a production caller"
        )
    )]
    pub(crate) fn to_rust(self) -> GraphEcfpHashScheme {
        match self {
            Self::Xxh3Width64V1() => GraphEcfpHashScheme::Xxh3Width64V1,
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
}
