//! Python values for molecular fingerprint results.

use std::vec::IntoIter;

use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyList};
use umol_graph::fingerprint::{
    BitFp as GraphBitFp, CountedFeatureSet as GraphCountedFeatureSet, FeatureSet as GraphFeatureSet,
};

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

    fn fold(&self, width: usize) -> PyResult<BitFp> {
        let folded = match &self.data {
            HashedFeatureSetData::U32(features) => features.fold(width),
            HashedFeatureSetData::U64(features) => features.fold(width),
            HashedFeatureSetData::U128(features) => features.fold(width),
        };
        folded
            .map(BitFp::from_rust)
            .map_err(|_| PyValueError::new_err("width must be positive"))
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

#[derive(Clone, Debug, PartialEq, Eq)]
enum CountedHashedFeatureSetData {
    U32(GraphCountedFeatureSet<u32>),
    U64(GraphCountedFeatureSet<u64>),
    U128(GraphCountedFeatureSet<u128>),
}

/// Immutable sparse set of hashed feature identifiers and occurrence counts.
#[pyclass(eq, frozen, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CountedHashedFeatureSet {
    data: CountedHashedFeatureSetData,
}

#[pymethods]
impl CountedHashedFeatureSet {
    /// Hash-identifier width in bits.
    #[getter]
    fn id_width(&self) -> u16 {
        match self.data {
            CountedHashedFeatureSetData::U32(_) => 32,
            CountedHashedFeatureSetData::U64(_) => 64,
            CountedHashedFeatureSetData::U128(_) => 128,
        }
    }

    /// Sorted `(identifier, count)` snapshot using ordinary Python integers.
    #[getter]
    fn entries(&self) -> Vec<(u128, u32)> {
        match &self.data {
            CountedHashedFeatureSetData::U32(features) => features
                .entries()
                .iter()
                .map(|(id, count)| (u128::from(*id), *count))
                .collect(),
            CountedHashedFeatureSetData::U64(features) => features
                .entries()
                .iter()
                .map(|(id, count)| (u128::from(*id), *count))
                .collect(),
            CountedHashedFeatureSetData::U128(features) => features.entries().to_vec(),
        }
    }

    /// Occurrence count of `id`, or zero when the identifier is absent.
    fn count(&self, id: u128) -> u32 {
        match &self.data {
            CountedHashedFeatureSetData::U32(features) => {
                u32::try_from(id).map(|id| features.count(&id)).unwrap_or(0)
            }
            CountedHashedFeatureSetData::U64(features) => {
                u64::try_from(id).map(|id| features.count(&id)).unwrap_or(0)
            }
            CountedHashedFeatureSetData::U128(features) => features.count(&id),
        }
    }

    fn __len__(&self) -> usize {
        match &self.data {
            CountedHashedFeatureSetData::U32(features) => features.len(),
            CountedHashedFeatureSetData::U64(features) => features.len(),
            CountedHashedFeatureSetData::U128(features) => features.len(),
        }
    }

    fn __iter__(&self) -> CountedHashedFeatureSetIter {
        CountedHashedFeatureSetIter {
            entries: self.entries().into_iter(),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "CountedHashedFeatureSet(entries={:?}, id_width={})",
            self.entries(),
            self.id_width()
        )
    }
}

impl CountedHashedFeatureSet {
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "Rust-to-Python conversion is used by molecular fingerprint operations"
        )
    )]
    pub(crate) fn from_rust<Id>(features: GraphCountedFeatureSet<Id>) -> Self
    where
        Self: From<GraphCountedFeatureSet<Id>>,
    {
        Self::from(features)
    }
}

impl From<GraphCountedFeatureSet<u32>> for CountedHashedFeatureSet {
    fn from(features: GraphCountedFeatureSet<u32>) -> Self {
        Self {
            data: CountedHashedFeatureSetData::U32(features),
        }
    }
}

impl From<GraphCountedFeatureSet<u64>> for CountedHashedFeatureSet {
    fn from(features: GraphCountedFeatureSet<u64>) -> Self {
        Self {
            data: CountedHashedFeatureSetData::U64(features),
        }
    }
}

impl From<GraphCountedFeatureSet<u128>> for CountedHashedFeatureSet {
    fn from(features: GraphCountedFeatureSet<u128>) -> Self {
        Self {
            data: CountedHashedFeatureSetData::U128(features),
        }
    }
}

#[pyclass]
struct CountedHashedFeatureSetIter {
    entries: IntoIter<(u128, u32)>,
}

#[pymethods]
impl CountedHashedFeatureSetIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<(u128, u32)> {
        self.entries.next()
    }
}

/// Immutable fixed-width bit fingerprint.
#[pyclass(eq, frozen, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BitFp {
    inner: GraphBitFp,
}

#[pymethods]
impl BitFp {
    /// Number of bits in the fingerprint.
    #[getter]
    fn width(&self) -> usize {
        self.inner.width()
    }

    /// Number of set bits.
    fn count_ones(&self) -> usize {
        self.inner.count_ones()
    }

    fn tanimoto(&self, other: &Self) -> PyResult<f64> {
        self.inner
            .tanimoto(&other.inner)
            .map_err(|_| self.width_mismatch(other))
    }

    fn dice(&self, other: &Self) -> PyResult<f64> {
        self.inner
            .dice(&other.inner)
            .map_err(|_| self.width_mismatch(other))
    }

    fn is_subset(&self, other: &Self) -> PyResult<bool> {
        self.inner
            .is_subset(&other.inner)
            .map_err(|_| self.width_mismatch(other))
    }

    fn __getitem__(&self, index: isize) -> PyResult<bool> {
        let index = if index < 0 {
            self.inner.width().checked_add_signed(index)
        } else {
            usize::try_from(index).ok()
        };
        index
            .and_then(|index| self.inner.get(index))
            .ok_or_else(|| PyIndexError::new_err("bit index out of range"))
    }

    fn __repr__(&self) -> String {
        format!(
            "BitFp(width={}, count_ones={})",
            self.width(),
            self.count_ones()
        )
    }
}

impl BitFp {
    pub(crate) fn from_rust(inner: GraphBitFp) -> Self {
        Self { inner }
    }

    fn width_mismatch(&self, other: &Self) -> PyErr {
        PyValueError::new_err(format!(
            "fingerprint width mismatch: {} != {}",
            self.width(),
            other.width()
        ))
    }
}

/// Immutable set of exact canonical structural feature keys.
#[pyclass(eq, frozen, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StructuralFeatureSet {
    inner: GraphFeatureSet<Vec<u8>>,
}

#[pymethods]
impl StructuralFeatureSet {
    /// Sorted key snapshot as Python bytes.
    #[getter]
    fn keys(&self, py: Python<'_>) -> Vec<Py<PyBytes>> {
        self.inner
            .ids()
            .iter()
            .map(|key| PyBytes::new(py, key).unbind())
            .collect()
    }

    fn is_subset(&self, other: &Self) -> bool {
        self.inner.is_subset(&other.inner)
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __iter__(&self) -> StructuralFeatureSetIter {
        StructuralFeatureSetIter {
            keys: self.inner.ids().to_vec().into_iter(),
        }
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let keys = PyList::new(py, self.keys(py))?
            .repr()?
            .extract::<String>()?;
        Ok(format!("StructuralFeatureSet(keys={keys})"))
    }
}

impl StructuralFeatureSet {
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "Rust-to-Python conversion is used by structural fingerprint operations"
        )
    )]
    pub(crate) fn from_rust(inner: GraphFeatureSet<Vec<u8>>) -> Self {
        Self { inner }
    }
}

#[pyclass]
struct StructuralFeatureSetIter {
    keys: IntoIter<Vec<u8>>,
}

#[pymethods]
impl StructuralFeatureSetIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> Option<Py<PyBytes>> {
        self.keys.next().map(|key| PyBytes::new(py, &key).unbind())
    }
}

#[cfg(test)]
mod tests {
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
    #[case::width32(
        HashedFeatureSet::from_rust(GraphFeatureSet::from_features([1u32, 2, 5, 9])),
        vec![false, true, true, false, false, true, false, false]
    )]
    #[case::width64(
        HashedFeatureSet::from_rust(GraphFeatureSet::from_features([1u64, 2, 5, 9])),
        vec![false, true, true, false, false, true, false, false]
    )]
    #[case::width128(
        HashedFeatureSet::from_rust(GraphFeatureSet::from_features([
            (u64::MAX as u128) + 2,
        ])),
        vec![false, true, false, false, false, false, false, false]
    )]
    fn test_hashed_feature_set_fold(
        #[case] features: HashedFeatureSet,
        #[case] expected: Vec<bool>,
    ) {
        Python::attach(|py| {
            let features = into_py_variant(py, features).unwrap();
            let folded = features.bind(py).call_method1("fold", (8,)).unwrap();

            assert_eq!(
                (0..8)
                    .map(|index| folded.get_item(index).unwrap().extract::<bool>().unwrap())
                    .collect::<Vec<_>>(),
                expected
            );
        });
    }

    #[rstest]
    #[case::zero_width(HashedFeatureSet::from_rust(GraphFeatureSet::from_features([
        1u64, 2,
    ])))]
    fn test_hashed_feature_set_fold_error(#[case] features: HashedFeatureSet) {
        Python::attach(|py| {
            let features = into_py_variant(py, features).unwrap();
            let error = features.bind(py).call_method1("fold", (0,)).unwrap_err();

            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "width must be positive"
            );
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

    #[rstest]
    #[case::width32(
        CountedHashedFeatureSet::from_rust(GraphCountedFeatureSet::from_counts([
            (3u32, 2),
            (1, 4),
        ])),
        vec![(1, 4), (3, 2)],
        32,
        "CountedHashedFeatureSet(entries=[(1, 4), (3, 2)], id_width=32)"
    )]
    #[case::width64(
        CountedHashedFeatureSet::from_rust(GraphCountedFeatureSet::from_counts([
            (3u64, 2),
            (1, 4),
        ])),
        vec![(1, 4), (3, 2)],
        64,
        "CountedHashedFeatureSet(entries=[(1, 4), (3, 2)], id_width=64)"
    )]
    #[case::width128(
        CountedHashedFeatureSet::from_rust(GraphCountedFeatureSet::from_counts([
            ((u64::MAX as u128) + 1, 2),
            (1, 4),
        ])),
        vec![(1, 4), ((u64::MAX as u128) + 1, 2)],
        128,
        "CountedHashedFeatureSet(entries=[(1, 4), (18446744073709551616, 2)], id_width=128)"
    )]
    fn test_counted_hashed_feature_set_value(
        #[case] features: CountedHashedFeatureSet,
        #[case] expected_entries: Vec<(u128, u32)>,
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
                    .getattr("entries")
                    .unwrap()
                    .extract::<Vec<(u128, u32)>>()
                    .unwrap(),
                expected_entries
            );
            assert_eq!(
                features
                    .getattr("id_width")
                    .unwrap()
                    .extract::<u16>()
                    .unwrap(),
                expected_width
            );
            assert_eq!(features.len().unwrap(), expected_entries.len());
            assert_eq!(
                features
                    .call_method0("__iter__")
                    .unwrap()
                    .try_iter()
                    .unwrap()
                    .map(|item| item.unwrap().extract::<(u128, u32)>().unwrap())
                    .collect::<Vec<_>>(),
                expected_entries
            );
            assert_eq!(
                features.repr().unwrap().extract::<String>().unwrap(),
                expected_repr
            );
        });
    }

    #[rstest]
    #[case::width32(
        CountedHashedFeatureSet::from_rust(GraphCountedFeatureSet::from_counts([(3u32, 2)])),
        3,
        2,
        u128::from(u32::MAX) + 1
    )]
    #[case::width64(
        CountedHashedFeatureSet::from_rust(GraphCountedFeatureSet::from_counts([(3u64, 2)])),
        3,
        2,
        u128::from(u64::MAX) + 1
    )]
    #[case::width128(
        CountedHashedFeatureSet::from_rust(GraphCountedFeatureSet::from_counts([
            ((u64::MAX as u128) + 1, 2),
        ])),
        (u64::MAX as u128) + 1,
        2,
        3
    )]
    fn test_counted_hashed_feature_set_count(
        #[case] features: CountedHashedFeatureSet,
        #[case] present: u128,
        #[case] expected: u32,
        #[case] absent: u128,
    ) {
        Python::attach(|py| {
            let features = into_py_variant(py, features).unwrap();
            let features = features.bind(py).as_any();

            assert_eq!(
                features
                    .call_method1("count", (present,))
                    .unwrap()
                    .extract::<u32>()
                    .unwrap(),
                expected
            );
            assert_eq!(
                features
                    .call_method1("count", (absent,))
                    .unwrap()
                    .extract::<u32>()
                    .unwrap(),
                0
            );
        });
    }

    #[rstest]
    #[case::width64(CountedHashedFeatureSet::from_rust(
        GraphCountedFeatureSet::from_counts([(1u64, 4), (3, 2)])
    ))]
    fn test_counted_hashed_feature_set_entries_detached(#[case] features: CountedHashedFeatureSet) {
        Python::attach(|py| {
            let features = into_py_variant(py, features).unwrap();
            let features = features.bind(py).as_any();
            let entries = features.getattr("entries").unwrap();

            entries
                .cast::<PyList>()
                .unwrap()
                .append((9u128, 6u32))
                .unwrap();

            assert_eq!(
                features
                    .getattr("entries")
                    .unwrap()
                    .extract::<Vec<(u128, u32)>>()
                    .unwrap(),
                vec![(1, 4), (3, 2)]
            );
        });
    }

    #[rstest]
    #[case::same_width(
        CountedHashedFeatureSet::from_rust(GraphCountedFeatureSet::from_counts([(1u64, 4)])),
        CountedHashedFeatureSet::from_rust(GraphCountedFeatureSet::from_counts([(1u64, 4)])),
        true
    )]
    #[case::different_counts(
        CountedHashedFeatureSet::from_rust(GraphCountedFeatureSet::from_counts([(1u64, 4)])),
        CountedHashedFeatureSet::from_rust(GraphCountedFeatureSet::from_counts([(1u64, 3)])),
        false
    )]
    #[case::different_widths(
        CountedHashedFeatureSet::from_rust(GraphCountedFeatureSet::from_counts([(1u32, 4)])),
        CountedHashedFeatureSet::from_rust(GraphCountedFeatureSet::from_counts([(1u64, 4)])),
        false
    )]
    fn test_counted_hashed_feature_set_eq(
        #[case] left: CountedHashedFeatureSet,
        #[case] right: CountedHashedFeatureSet,
        #[case] expected: bool,
    ) {
        assert_eq!(left == right, expected);
    }

    #[rstest]
    #[case::nonempty(
        BitFp::from_rust(GraphFeatureSet::from_features([1u64, 2, 5, 9]).fold(8).unwrap()),
        8,
        3,
        "BitFp(width=8, count_ones=3)"
    )]
    #[case::empty(
        BitFp::from_rust(GraphFeatureSet::from_features(Vec::<u64>::new()).fold(4).unwrap()),
        4,
        0,
        "BitFp(width=4, count_ones=0)"
    )]
    fn test_bit_fp_value(
        #[case] fingerprint: BitFp,
        #[case] expected_width: usize,
        #[case] expected_count: usize,
        #[case] expected_repr: &str,
    ) {
        Python::attach(|py| {
            let expected = into_py_variant(py, fingerprint.clone()).unwrap();
            let fingerprint = into_py_variant(py, fingerprint).unwrap();
            let expected = expected.bind(py).as_any();
            let fingerprint = fingerprint.bind(py).as_any();

            assert!(fingerprint.eq(expected).unwrap());
            assert_eq!(
                fingerprint
                    .getattr("width")
                    .unwrap()
                    .extract::<usize>()
                    .unwrap(),
                expected_width
            );
            assert_eq!(
                fingerprint
                    .call_method0("count_ones")
                    .unwrap()
                    .extract::<usize>()
                    .unwrap(),
                expected_count
            );
            assert_eq!(
                fingerprint.repr().unwrap().extract::<String>().unwrap(),
                expected_repr
            );
        });
    }

    #[rstest]
    #[case::first(BitFp::from_rust(
        GraphFeatureSet::from_features([1u64, 7]).fold(8).unwrap()
    ), 0, false)]
    #[case::last(BitFp::from_rust(
        GraphFeatureSet::from_features([1u64, 7]).fold(8).unwrap()
    ), 7, true)]
    #[case::negative_last(BitFp::from_rust(
        GraphFeatureSet::from_features([1u64, 7]).fold(8).unwrap()
    ), -1, true)]
    fn test_bit_fp_getitem(
        #[case] fingerprint: BitFp,
        #[case] index: isize,
        #[case] expected: bool,
    ) {
        assert_eq!(fingerprint.__getitem__(index).unwrap(), expected);
    }

    #[rstest]
    #[case::first_invalid(BitFp::from_rust(
        GraphFeatureSet::from_features([1u64, 7]).fold(8).unwrap()
    ), 8)]
    #[case::negative_invalid(BitFp::from_rust(
        GraphFeatureSet::from_features([1u64, 7]).fold(8).unwrap()
    ), -9)]
    fn test_bit_fp_getitem_error(#[case] fingerprint: BitFp, #[case] index: isize) {
        Python::attach(|py| {
            let error = fingerprint.__getitem__(index).unwrap_err();

            assert!(error.is_instance_of::<PyIndexError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "bit index out of range"
            );
        });
    }

    #[rstest]
    #[case::width8(
        BitFp::from_rust(GraphFeatureSet::from_features([1u64, 2, 3]).fold(8).unwrap()),
        BitFp::from_rust(GraphFeatureSet::from_features([2u64, 3, 4]).fold(8).unwrap()),
        BitFp::from_rust(GraphFeatureSet::from_features([2u64, 3]).fold(8).unwrap())
    )]
    fn test_bit_fp_operations(#[case] left: BitFp, #[case] right: BitFp, #[case] subset: BitFp) {
        assert_eq!(left.tanimoto(&right).unwrap(), 0.5);
        assert_eq!(left.dice(&right).unwrap(), 4.0 / 6.0);
        assert!(!left.is_subset(&right).unwrap());
        assert!(subset.is_subset(&right).unwrap());
    }

    #[rstest]
    #[case::unequal_widths(
        BitFp::from_rust(GraphFeatureSet::from_features([1u64]).fold(8).unwrap()),
        BitFp::from_rust(GraphFeatureSet::from_features([1u64]).fold(4).unwrap())
    )]
    fn test_bit_fp_operations_error(#[case] left: BitFp, #[case] right: BitFp) {
        Python::attach(|py| {
            for error in [
                left.tanimoto(&right).unwrap_err(),
                left.dice(&right).unwrap_err(),
                left.is_subset(&right).unwrap_err(),
            ] {
                assert!(error.is_instance_of::<PyValueError>(py));
                assert_eq!(
                    error.value(py).str().unwrap().extract::<String>().unwrap(),
                    "fingerprint width mismatch: 8 != 4"
                );
            }
        });
    }

    #[rstest]
    #[case::same_bits(
        BitFp::from_rust(GraphFeatureSet::from_features([1u64, 2]).fold(8).unwrap()),
        BitFp::from_rust(GraphFeatureSet::from_features([1u64, 2]).fold(8).unwrap()),
        true
    )]
    #[case::different_bits(
        BitFp::from_rust(GraphFeatureSet::from_features([1u64, 2]).fold(8).unwrap()),
        BitFp::from_rust(GraphFeatureSet::from_features([1u64, 3]).fold(8).unwrap()),
        false
    )]
    #[case::different_widths(
        BitFp::from_rust(GraphFeatureSet::from_features([1u64, 2]).fold(8).unwrap()),
        BitFp::from_rust(GraphFeatureSet::from_features([1u64, 2]).fold(4).unwrap()),
        false
    )]
    fn test_bit_fp_eq(#[case] left: BitFp, #[case] right: BitFp, #[case] expected: bool) {
        assert_eq!(left == right, expected);
    }

    #[rstest]
    #[case::keys(
        StructuralFeatureSet::from_rust(GraphFeatureSet::from_features([
            b"b".to_vec(),
            b"a\0b".to_vec(),
            b"a".to_vec(),
            b"a".to_vec(),
        ])),
        vec![b"a".to_vec(), b"a\0b".to_vec(), b"b".to_vec()],
        "StructuralFeatureSet(keys=[b'a', b'a\\x00b', b'b'])"
    )]
    fn test_structural_feature_set_value(
        #[case] features: StructuralFeatureSet,
        #[case] expected_keys: Vec<Vec<u8>>,
        #[case] expected_repr: &str,
    ) {
        Python::attach(|py| {
            let expected = into_py_variant(py, features.clone()).unwrap();
            let features = into_py_variant(py, features).unwrap();
            let expected = expected.bind(py).as_any();
            let features = features.bind(py).as_any();
            let keys = features.getattr("keys").unwrap();

            assert!(features.eq(expected).unwrap());
            assert_eq!(keys.extract::<Vec<Vec<u8>>>().unwrap(), expected_keys);
            assert!(keys
                .cast::<PyList>()
                .unwrap()
                .iter()
                .all(|key| key.is_instance_of::<PyBytes>()));
            assert_eq!(features.len().unwrap(), expected_keys.len());
            assert_eq!(
                features
                    .call_method0("__iter__")
                    .unwrap()
                    .try_iter()
                    .unwrap()
                    .map(|item| {
                        let item = item.unwrap();
                        assert!(item.is_instance_of::<PyBytes>());
                        item.extract::<Vec<u8>>().unwrap()
                    })
                    .collect::<Vec<_>>(),
                expected_keys
            );
            assert_eq!(
                features.repr().unwrap().extract::<String>().unwrap(),
                expected_repr
            );
        });
    }

    #[rstest]
    #[case::proper_subset(
        StructuralFeatureSet::from_rust(GraphFeatureSet::from_features([
            b"a".to_vec(),
            b"b".to_vec(),
        ])),
        StructuralFeatureSet::from_rust(GraphFeatureSet::from_features([
            b"a".to_vec(),
            b"b".to_vec(),
            b"c".to_vec(),
        ])),
        true
    )]
    #[case::missing_key(
        StructuralFeatureSet::from_rust(GraphFeatureSet::from_features([
            b"a".to_vec(),
            b"d".to_vec(),
        ])),
        StructuralFeatureSet::from_rust(GraphFeatureSet::from_features([
            b"a".to_vec(),
            b"b".to_vec(),
            b"c".to_vec(),
        ])),
        false
    )]
    #[case::empty(
        StructuralFeatureSet::from_rust(GraphFeatureSet::from_features(Vec::<Vec<u8>>::new())),
        StructuralFeatureSet::from_rust(GraphFeatureSet::from_features([b"a".to_vec()])),
        true
    )]
    fn test_structural_feature_set_is_subset(
        #[case] query: StructuralFeatureSet,
        #[case] target: StructuralFeatureSet,
        #[case] expected: bool,
    ) {
        assert_eq!(query.is_subset(&target), expected);
    }

    #[rstest]
    #[case::snapshot(StructuralFeatureSet::from_rust(
        GraphFeatureSet::from_features([b"a".to_vec(), b"b".to_vec()])
    ))]
    fn test_structural_feature_set_keys_snapshot(#[case] features: StructuralFeatureSet) {
        Python::attach(|py| {
            let features = into_py_variant(py, features).unwrap();
            let features = features.bind(py).as_any();
            let keys = features.getattr("keys").unwrap();

            keys.cast::<PyList>()
                .unwrap()
                .append(PyBytes::new(py, b"c"))
                .unwrap();

            assert_eq!(
                features
                    .getattr("keys")
                    .unwrap()
                    .extract::<Vec<Vec<u8>>>()
                    .unwrap(),
                vec![b"a".to_vec(), b"b".to_vec()]
            );
        });
    }

    #[rstest]
    #[case::same_keys(
        StructuralFeatureSet::from_rust(GraphFeatureSet::from_features([
            b"a".to_vec(),
            b"b".to_vec(),
        ])),
        StructuralFeatureSet::from_rust(GraphFeatureSet::from_features([
            b"a".to_vec(),
            b"b".to_vec(),
        ])),
        true
    )]
    #[case::different_keys(
        StructuralFeatureSet::from_rust(GraphFeatureSet::from_features([
            b"a".to_vec(),
            b"b".to_vec(),
        ])),
        StructuralFeatureSet::from_rust(GraphFeatureSet::from_features([
            b"a".to_vec(),
            b"c".to_vec(),
        ])),
        false
    )]
    fn test_structural_feature_set_eq(
        #[case] left: StructuralFeatureSet,
        #[case] right: StructuralFeatureSet,
        #[case] expected: bool,
    ) {
        assert_eq!(left == right, expected);
    }
}
