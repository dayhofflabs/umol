//! Python bindings for stereo-model values.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use pyo3::prelude::*;
use umol_graph::ops::model::InconsistencyPolicy as GraphInconsistencyPolicy;

/// Policy for stereo assertions that cannot be fully realized.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum InconsistencyPolicy {
    Keep,
    Strip,
    Error,
}

impl InconsistencyPolicy {
    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion API for StereoModel configuration"
    )]
    pub(crate) fn from_rust(policy: GraphInconsistencyPolicy) -> Self {
        match policy {
            GraphInconsistencyPolicy::Keep => Self::Keep,
            GraphInconsistencyPolicy::Strip => Self::Strip,
            GraphInconsistencyPolicy::Error => Self::Error,
        }
    }

    #[allow(
        dead_code,
        reason = "Python-to-Rust conversion API for StereoModel configuration"
    )]
    pub(crate) fn to_rust(self) -> GraphInconsistencyPolicy {
        match self {
            Self::Keep => GraphInconsistencyPolicy::Keep,
            Self::Strip => GraphInconsistencyPolicy::Strip,
            Self::Error => GraphInconsistencyPolicy::Error,
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

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
}
