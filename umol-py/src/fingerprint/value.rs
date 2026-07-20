//! Python values for molecular fingerprint results.

use std::vec::IntoIter;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use umol_graph::fingerprint::FeatureSet as GraphFeatureSet;

#[derive(Clone, Debug, PartialEq, Eq)]
enum HashedFeatureSetData {
    U32(GraphFeatureSet<u32>),
    U64(GraphFeatureSet<u64>),
    U128(GraphFeatureSet<u128>),
}

/// Immutable sparse set of hashed feature identifiers.
#[pyclass(eq, frozen, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HashedFeatureSet {
    data: HashedFeatureSetData,
}

#[pymethods]
impl HashedFeatureSet {
    /// Hash-identifier width in bits.
    #[getter]
    fn id_width(&self) -> u16 {
        match self.data {
            HashedFeatureSetData::U32(_) => 32,
            HashedFeatureSetData::U64(_) => 64,
            HashedFeatureSetData::U128(_) => 128,
        }
    }

    /// Sorted identifier snapshot as ordinary Python integers.
    #[getter]
    fn ids(&self) -> Vec<u128> {
        match &self.data {
            HashedFeatureSetData::U32(features) => {
                features.ids().iter().copied().map(u128::from).collect()
            }
            HashedFeatureSetData::U64(features) => {
                features.ids().iter().copied().map(u128::from).collect()
            }
            HashedFeatureSetData::U128(features) => features.ids().to_vec(),
        }
    }

    fn tanimoto(&self, other: &Self) -> PyResult<f64> {
        match (&self.data, &other.data) {
            (HashedFeatureSetData::U32(left), HashedFeatureSetData::U32(right)) => {
                Ok(left.tanimoto(right))
            }
            (HashedFeatureSetData::U64(left), HashedFeatureSetData::U64(right)) => {
                Ok(left.tanimoto(right))
            }
            (HashedFeatureSetData::U128(left), HashedFeatureSetData::U128(right)) => {
                Ok(left.tanimoto(right))
            }
            _ => Err(self.width_mismatch(other)),
        }
    }

    fn dice(&self, other: &Self) -> PyResult<f64> {
        match (&self.data, &other.data) {
            (HashedFeatureSetData::U32(left), HashedFeatureSetData::U32(right)) => {
                Ok(left.dice(right))
            }
            (HashedFeatureSetData::U64(left), HashedFeatureSetData::U64(right)) => {
                Ok(left.dice(right))
            }
            (HashedFeatureSetData::U128(left), HashedFeatureSetData::U128(right)) => {
                Ok(left.dice(right))
            }
            _ => Err(self.width_mismatch(other)),
        }
    }

    fn is_subset(&self, other: &Self) -> PyResult<bool> {
        match (&self.data, &other.data) {
            (HashedFeatureSetData::U32(left), HashedFeatureSetData::U32(right)) => {
                Ok(left.is_subset(right))
            }
            (HashedFeatureSetData::U64(left), HashedFeatureSetData::U64(right)) => {
                Ok(left.is_subset(right))
            }
            (HashedFeatureSetData::U128(left), HashedFeatureSetData::U128(right)) => {
                Ok(left.is_subset(right))
            }
            _ => Err(self.width_mismatch(other)),
        }
    }

    fn __len__(&self) -> usize {
        match &self.data {
            HashedFeatureSetData::U32(features) => features.len(),
            HashedFeatureSetData::U64(features) => features.len(),
            HashedFeatureSetData::U128(features) => features.len(),
        }
    }

    fn __iter__(&self) -> HashedFeatureSetIter {
        HashedFeatureSetIter {
            ids: self.ids().into_iter(),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "HashedFeatureSet(ids={:?}, id_width={})",
            self.ids(),
            self.id_width()
        )
    }
}

impl HashedFeatureSet {
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "Rust-to-Python conversion is used by molecular fingerprint operations"
        )
    )]
    pub(crate) fn from_rust<Id>(features: GraphFeatureSet<Id>) -> Self
    where
        Self: From<GraphFeatureSet<Id>>,
    {
        Self::from(features)
    }

    fn width_mismatch(&self, other: &Self) -> PyErr {
        PyValueError::new_err(format!(
            "identifier width mismatch: {} != {}",
            self.id_width(),
            other.id_width()
        ))
    }
}

impl From<GraphFeatureSet<u32>> for HashedFeatureSet {
    fn from(features: GraphFeatureSet<u32>) -> Self {
        Self {
            data: HashedFeatureSetData::U32(features),
        }
    }
}

impl From<GraphFeatureSet<u64>> for HashedFeatureSet {
    fn from(features: GraphFeatureSet<u64>) -> Self {
        Self {
            data: HashedFeatureSetData::U64(features),
        }
    }
}

impl From<GraphFeatureSet<u128>> for HashedFeatureSet {
    fn from(features: GraphFeatureSet<u128>) -> Self {
        Self {
            data: HashedFeatureSetData::U128(features),
        }
    }
}

#[pyclass]
struct HashedFeatureSetIter {
    ids: IntoIter<u128>,
}

#[pymethods]
impl HashedFeatureSetIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<u128> {
        self.ids.next()
    }
}

#[cfg(test)]
mod tests {
    use pyo3::types::PyList;
    use rstest::rstest;

    use super::*;
    use crate::convert::into_py_variant;

    #[rstest]
    #[case::width32(
        HashedFeatureSet::from_rust(GraphFeatureSet::from_features([3u32, 1, 2, 1])),
        vec![1, 2, 3],
        32,
        "HashedFeatureSet(ids=[1, 2, 3], id_width=32)"
    )]
    #[case::width64(
        HashedFeatureSet::from_rust(GraphFeatureSet::from_features([3u64, 1, 2, 1])),
        vec![1, 2, 3],
        64,
        "HashedFeatureSet(ids=[1, 2, 3], id_width=64)"
    )]
    #[case::width128(
        HashedFeatureSet::from_rust(GraphFeatureSet::from_features([
            (u64::MAX as u128) + 1,
            1,
        ])),
        vec![1, (u64::MAX as u128) + 1],
        128,
        "HashedFeatureSet(ids=[1, 18446744073709551616], id_width=128)"
    )]
    fn test_hashed_feature_set_value(
        #[case] features: HashedFeatureSet,
        #[case] expected_ids: Vec<u128>,
        #[case] expected_width: u16,
        #[case] expected_repr: &str,
    ) {
        Python::attach(|py| {
            let expected = into_py_variant(py, features.clone()).unwrap();
            let features = into_py_variant(py, features).unwrap();
            let expected = expected.bind(py).as_any();
            let features = features.bind(py).as_any();

            assert!(features.eq(expected).unwrap());
            assert_eq!(
                features
                    .getattr("ids")
                    .unwrap()
                    .extract::<Vec<u128>>()
                    .unwrap(),
                expected_ids
            );
            assert_eq!(
                features
                    .getattr("id_width")
                    .unwrap()
                    .extract::<u16>()
                    .unwrap(),
                expected_width
            );
            assert_eq!(features.len().unwrap(), expected_ids.len());
            assert_eq!(
                features
                    .call_method0("__iter__")
                    .unwrap()
                    .try_iter()
                    .unwrap()
                    .map(|item| item.unwrap().extract::<u128>().unwrap())
                    .collect::<Vec<_>>(),
                expected_ids
            );
            assert_eq!(
                features.repr().unwrap().extract::<String>().unwrap(),
                expected_repr
            );
        });
    }

    #[rstest]
    #[case::width32(
        HashedFeatureSet::from_rust(GraphFeatureSet::from_features([1u32, 2, 3])),
        HashedFeatureSet::from_rust(GraphFeatureSet::from_features([2u32, 3, 4])),
        HashedFeatureSet::from_rust(GraphFeatureSet::from_features([2u32, 3]))
    )]
    #[case::width64(
        HashedFeatureSet::from_rust(GraphFeatureSet::from_features([1u64, 2, 3])),
        HashedFeatureSet::from_rust(GraphFeatureSet::from_features([2u64, 3, 4])),
        HashedFeatureSet::from_rust(GraphFeatureSet::from_features([2u64, 3]))
    )]
    #[case::width128(
        HashedFeatureSet::from_rust(GraphFeatureSet::from_features([1u128, 2, 3])),
        HashedFeatureSet::from_rust(GraphFeatureSet::from_features([2u128, 3, 4])),
        HashedFeatureSet::from_rust(GraphFeatureSet::from_features([2u128, 3]))
    )]
    fn test_hashed_feature_set_operations(
        #[case] left: HashedFeatureSet,
        #[case] right: HashedFeatureSet,
        #[case] subset: HashedFeatureSet,
    ) {
        assert_eq!(left.tanimoto(&right).unwrap(), 0.5);
        assert_eq!(left.dice(&right).unwrap(), 4.0 / 6.0);
        assert!(!left.is_subset(&right).unwrap());
        assert!(subset.is_subset(&right).unwrap());
    }

    #[rstest]
    #[case::width32_width64(
        HashedFeatureSet::from_rust(GraphFeatureSet::from_features([1u32])),
        HashedFeatureSet::from_rust(GraphFeatureSet::from_features([1u64]))
    )]
    fn test_hashed_feature_set_operations_error(
        #[case] left: HashedFeatureSet,
        #[case] right: HashedFeatureSet,
    ) {
        Python::attach(|py| {
            for error in [
                left.tanimoto(&right).unwrap_err(),
                left.dice(&right).unwrap_err(),
                left.is_subset(&right).unwrap_err(),
            ] {
                assert!(error.is_instance_of::<PyValueError>(py));
                assert_eq!(
                    error.value(py).str().unwrap().extract::<String>().unwrap(),
                    "identifier width mismatch: 32 != 64"
                );
            }
        });
    }

    #[rstest]
    #[case::width64(HashedFeatureSet::from_rust(GraphFeatureSet::from_features([
        1u64, 2, 3,
    ])))]
    fn test_hashed_feature_set_ids_detached(#[case] features: HashedFeatureSet) {
        Python::attach(|py| {
            let features = into_py_variant(py, features).unwrap();
            let features = features.bind(py).as_any();
            let ids = features.getattr("ids").unwrap();

            ids.cast::<PyList>().unwrap().append(9u128).unwrap();

            assert_eq!(
                features
                    .getattr("ids")
                    .unwrap()
                    .extract::<Vec<u128>>()
                    .unwrap(),
                vec![1, 2, 3]
            );
        });
    }

    #[rstest]
    #[case::same_width(
        HashedFeatureSet::from_rust(GraphFeatureSet::from_features([1u64, 2])),
        HashedFeatureSet::from_rust(GraphFeatureSet::from_features([1u64, 2])),
        true
    )]
    #[case::different_ids(
        HashedFeatureSet::from_rust(GraphFeatureSet::from_features([1u64, 2])),
        HashedFeatureSet::from_rust(GraphFeatureSet::from_features([1u64, 3])),
        false
    )]
    #[case::different_widths(
        HashedFeatureSet::from_rust(GraphFeatureSet::from_features([1u32, 2])),
        HashedFeatureSet::from_rust(GraphFeatureSet::from_features([1u64, 2])),
        false
    )]
    fn test_hashed_feature_set_eq(
        #[case] left: HashedFeatureSet,
        #[case] right: HashedFeatureSet,
        #[case] expected: bool,
    ) {
        assert_eq!(left == right, expected);
    }
}
