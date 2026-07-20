//! Python values for reaction fingerprint results.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use std::vec::IntoIter;

use pyo3::prelude::*;
use pyo3::types::PyList;
use umol_graph::fingerprint::{
    FeatureSet as GraphFeatureSet, ReactionSide as GraphReactionSide,
    SignedFeatureSet as GraphSignedFeatureSet,
};

/// Side of a reaction from which a fingerprint feature originates.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ReactionSide {
    Reactant,
    Product,
}

impl ReactionSide {
    pub(crate) fn from_rust(side: GraphReactionSide) -> Self {
        match side {
            GraphReactionSide::Reactant => Self::Reactant,
            GraphReactionSide::Product => Self::Product,
        }
    }

    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "Python-to-Rust conversion is used by reaction fingerprint operations"
        )
    )]
    pub(crate) fn to_rust(self) -> GraphReactionSide {
        match self {
            Self::Reactant => GraphReactionSide::Reactant,
            Self::Product => GraphReactionSide::Product,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum SignedHashedFeatureSetData {
    U32(GraphSignedFeatureSet<u32>),
    U64(GraphSignedFeatureSet<u64>),
    U128(GraphSignedFeatureSet<u128>),
}

/// Immutable sparse set of hashed identifiers and signed reaction differences.
#[pyclass(eq, frozen, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SignedHashedFeatureSet {
    data: SignedHashedFeatureSetData,
}

#[pymethods]
impl SignedHashedFeatureSet {
    /// Hash-identifier width in bits.
    #[getter]
    fn id_width(&self) -> u16 {
        match self.data {
            SignedHashedFeatureSetData::U32(_) => 32,
            SignedHashedFeatureSetData::U64(_) => 64,
            SignedHashedFeatureSetData::U128(_) => 128,
        }
    }

    /// Sorted `(identifier, signed_count)` snapshot using ordinary Python integers.
    #[getter]
    fn entries(&self) -> Vec<(u128, i32)> {
        match &self.data {
            SignedHashedFeatureSetData::U32(features) => features
                .entries()
                .iter()
                .map(|(id, count)| (u128::from(*id), *count))
                .collect(),
            SignedHashedFeatureSetData::U64(features) => features
                .entries()
                .iter()
                .map(|(id, count)| (u128::from(*id), *count))
                .collect(),
            SignedHashedFeatureSetData::U128(features) => features.entries().to_vec(),
        }
    }

    /// Signed count of `id`, or zero when the identifier is absent.
    fn count(&self, id: u128) -> i32 {
        match &self.data {
            SignedHashedFeatureSetData::U32(features) => {
                u32::try_from(id).map(|id| features.count(&id)).unwrap_or(0)
            }
            SignedHashedFeatureSetData::U64(features) => {
                u64::try_from(id).map(|id| features.count(&id)).unwrap_or(0)
            }
            SignedHashedFeatureSetData::U128(features) => features.count(&id),
        }
    }

    fn __len__(&self) -> usize {
        match &self.data {
            SignedHashedFeatureSetData::U32(features) => features.len(),
            SignedHashedFeatureSetData::U64(features) => features.len(),
            SignedHashedFeatureSetData::U128(features) => features.len(),
        }
    }

    fn __iter__(&self) -> SignedHashedFeatureSetIter {
        SignedHashedFeatureSetIter {
            entries: self.entries().into_iter(),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "SignedHashedFeatureSet(entries={:?}, id_width={})",
            self.entries(),
            self.id_width()
        )
    }
}

impl SignedHashedFeatureSet {
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "Rust-to-Python conversion is used by reaction fingerprint operations"
        )
    )]
    pub(crate) fn from_rust<Id>(features: GraphSignedFeatureSet<Id>) -> Self
    where
        Self: From<GraphSignedFeatureSet<Id>>,
    {
        Self::from(features)
    }
}

impl From<GraphSignedFeatureSet<u32>> for SignedHashedFeatureSet {
    fn from(features: GraphSignedFeatureSet<u32>) -> Self {
        Self {
            data: SignedHashedFeatureSetData::U32(features),
        }
    }
}

impl From<GraphSignedFeatureSet<u64>> for SignedHashedFeatureSet {
    fn from(features: GraphSignedFeatureSet<u64>) -> Self {
        Self {
            data: SignedHashedFeatureSetData::U64(features),
        }
    }
}

impl From<GraphSignedFeatureSet<u128>> for SignedHashedFeatureSet {
    fn from(features: GraphSignedFeatureSet<u128>) -> Self {
        Self {
            data: SignedHashedFeatureSetData::U128(features),
        }
    }
}

#[pyclass]
struct SignedHashedFeatureSetIter {
    entries: IntoIter<(u128, i32)>,
}

#[pymethods]
impl SignedHashedFeatureSetIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<(u128, i32)> {
        self.entries.next()
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum RoleTaggedHashedFeatureSetData {
    U32(GraphFeatureSet<(GraphReactionSide, u32)>),
    U64(GraphFeatureSet<(GraphReactionSide, u64)>),
    U128(GraphFeatureSet<(GraphReactionSide, u128)>),
}

/// Immutable sparse set of reaction-side-tagged hashed identifiers.
#[pyclass(eq, frozen, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RoleTaggedHashedFeatureSet {
    data: RoleTaggedHashedFeatureSetData,
}

#[pymethods]
impl RoleTaggedHashedFeatureSet {
    /// Hash-identifier width in bits.
    #[getter]
    fn id_width(&self) -> u16 {
        match self.data {
            RoleTaggedHashedFeatureSetData::U32(_) => 32,
            RoleTaggedHashedFeatureSetData::U64(_) => 64,
            RoleTaggedHashedFeatureSetData::U128(_) => 128,
        }
    }

    /// Sorted `(reaction_side, identifier)` snapshot using ordinary Python integers.
    #[getter]
    fn ids(&self) -> Vec<(ReactionSide, u128)> {
        match &self.data {
            RoleTaggedHashedFeatureSetData::U32(features) => features
                .ids()
                .iter()
                .map(|(side, id)| (ReactionSide::from_rust(*side), u128::from(*id)))
                .collect(),
            RoleTaggedHashedFeatureSetData::U64(features) => features
                .ids()
                .iter()
                .map(|(side, id)| (ReactionSide::from_rust(*side), u128::from(*id)))
                .collect(),
            RoleTaggedHashedFeatureSetData::U128(features) => features
                .ids()
                .iter()
                .map(|(side, id)| (ReactionSide::from_rust(*side), *id))
                .collect(),
        }
    }

    fn __len__(&self) -> usize {
        match &self.data {
            RoleTaggedHashedFeatureSetData::U32(features) => features.len(),
            RoleTaggedHashedFeatureSetData::U64(features) => features.len(),
            RoleTaggedHashedFeatureSetData::U128(features) => features.len(),
        }
    }

    fn __iter__(&self) -> RoleTaggedHashedFeatureSetIter {
        RoleTaggedHashedFeatureSetIter {
            ids: self.ids().into_iter(),
        }
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let ids = PyList::new(py, self.ids())?.repr()?.extract::<String>()?;
        Ok(format!(
            "RoleTaggedHashedFeatureSet(ids={ids}, id_width={})",
            self.id_width()
        ))
    }
}

impl RoleTaggedHashedFeatureSet {
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "Rust-to-Python conversion is used by reaction fingerprint operations"
        )
    )]
    pub(crate) fn from_rust<Id>(features: GraphFeatureSet<(GraphReactionSide, Id)>) -> Self
    where
        Self: From<GraphFeatureSet<(GraphReactionSide, Id)>>,
    {
        Self::from(features)
    }
}

impl From<GraphFeatureSet<(GraphReactionSide, u32)>> for RoleTaggedHashedFeatureSet {
    fn from(features: GraphFeatureSet<(GraphReactionSide, u32)>) -> Self {
        Self {
            data: RoleTaggedHashedFeatureSetData::U32(features),
        }
    }
}

impl From<GraphFeatureSet<(GraphReactionSide, u64)>> for RoleTaggedHashedFeatureSet {
    fn from(features: GraphFeatureSet<(GraphReactionSide, u64)>) -> Self {
        Self {
            data: RoleTaggedHashedFeatureSetData::U64(features),
        }
    }
}

impl From<GraphFeatureSet<(GraphReactionSide, u128)>> for RoleTaggedHashedFeatureSet {
    fn from(features: GraphFeatureSet<(GraphReactionSide, u128)>) -> Self {
        Self {
            data: RoleTaggedHashedFeatureSetData::U128(features),
        }
    }
}

#[pyclass]
struct RoleTaggedHashedFeatureSetIter {
    ids: IntoIter<(ReactionSide, u128)>,
}

#[pymethods]
impl RoleTaggedHashedFeatureSetIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<(ReactionSide, u128)> {
        self.ids.next()
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph::fingerprint::CountedFeatureSet as GraphCountedFeatureSet;

    use super::*;

    #[rstest]
    #[case::reactant(GraphReactionSide::Reactant, ReactionSide::Reactant)]
    #[case::product(GraphReactionSide::Product, ReactionSide::Product)]
    fn test_reaction_side_from_rust(
        #[case] side: GraphReactionSide,
        #[case] expected: ReactionSide,
    ) {
        assert_eq!(ReactionSide::from_rust(side), expected);
    }

    #[rstest]
    #[case::reactant(ReactionSide::Reactant, GraphReactionSide::Reactant)]
    #[case::product(ReactionSide::Product, GraphReactionSide::Product)]
    fn test_reaction_side_to_rust(#[case] side: ReactionSide, #[case] expected: GraphReactionSide) {
        assert_eq!(side.to_rust(), expected);
    }

    #[rstest]
    #[case::reactant(
        ReactionSide::Reactant,
        ReactionSide::Reactant,
        ReactionSide::Product,
        "ReactionSide.Reactant"
    )]
    #[case::product(
        ReactionSide::Product,
        ReactionSide::Product,
        ReactionSide::Reactant,
        "ReactionSide.Product"
    )]
    fn test_reaction_side_python_value(
        #[case] side: ReactionSide,
        #[case] equal: ReactionSide,
        #[case] unequal: ReactionSide,
        #[case] expected_repr: &str,
    ) {
        Python::attach(|py| {
            let side = Py::new(py, side).unwrap();
            let equal = Py::new(py, equal).unwrap();
            let unequal = Py::new(py, unequal).unwrap();

            assert!(side.bind(py).as_any().eq(equal.bind(py).as_any()).unwrap());
            assert!(!side
                .bind(py)
                .as_any()
                .eq(unequal.bind(py).as_any())
                .unwrap());
            assert_eq!(
                side.bind(py).as_any().hash().unwrap(),
                equal.bind(py).as_any().hash().unwrap()
            );
            assert_eq!(
                side.bind(py)
                    .as_any()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                expected_repr
            );
        });
    }

    #[rstest]
    #[case::width32(
        SignedHashedFeatureSet::from_rust(GraphSignedFeatureSet::difference(
            &GraphCountedFeatureSet::from_counts([(3u32, 1), (1, 3), (4, 2)]),
            &GraphCountedFeatureSet::from_counts([(2u32, 4), (1, 1), (4, 2)]),
        )),
        vec![(1, 2), (2, -4), (3, 1)],
        32,
        "SignedHashedFeatureSet(entries=[(1, 2), (2, -4), (3, 1)], id_width=32)"
    )]
    #[case::width64(
        SignedHashedFeatureSet::from_rust(GraphSignedFeatureSet::difference(
            &GraphCountedFeatureSet::from_counts([(3u64, 1), (1, 3), (4, 2)]),
            &GraphCountedFeatureSet::from_counts([(2u64, 4), (1, 1), (4, 2)]),
        )),
        vec![(1, 2), (2, -4), (3, 1)],
        64,
        "SignedHashedFeatureSet(entries=[(1, 2), (2, -4), (3, 1)], id_width=64)"
    )]
    #[case::width128(
        SignedHashedFeatureSet::from_rust(GraphSignedFeatureSet::difference(
            &GraphCountedFeatureSet::from_counts([
                ((u64::MAX as u128) + 1, 2),
                (1, 3),
            ]),
            &GraphCountedFeatureSet::from_counts([(2u128, 4), (1, 1)]),
        )),
        vec![(1, 2), (2, -4), ((u64::MAX as u128) + 1, 2)],
        128,
        "SignedHashedFeatureSet(entries=[(1, 2), (2, -4), (18446744073709551616, 2)], id_width=128)"
    )]
    fn test_signed_hashed_feature_set_value(
        #[case] features: SignedHashedFeatureSet,
        #[case] expected_entries: Vec<(u128, i32)>,
        #[case] expected_width: u16,
        #[case] expected_repr: &str,
    ) {
        Python::attach(|py| {
            let expected = Py::new(py, features.clone()).unwrap();
            let features = Py::new(py, features).unwrap();
            let expected = expected.bind(py).as_any();
            let features = features.bind(py).as_any();

            assert!(features.eq(expected).unwrap());
            assert_eq!(
                features
                    .getattr("entries")
                    .unwrap()
                    .extract::<Vec<(u128, i32)>>()
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
                    .map(|item| item.unwrap().extract::<(u128, i32)>().unwrap())
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
        SignedHashedFeatureSet::from_rust(GraphSignedFeatureSet::difference(
            &GraphCountedFeatureSet::from_counts([(1u32, 3), (3, 1)]),
            &GraphCountedFeatureSet::from_counts([(1u32, 1), (2, 4)]),
        )),
        1,
        2,
        2,
        -4,
        u128::from(u32::MAX) + 1
    )]
    #[case::width64(
        SignedHashedFeatureSet::from_rust(GraphSignedFeatureSet::difference(
            &GraphCountedFeatureSet::from_counts([(1u64, 3), (3, 1)]),
            &GraphCountedFeatureSet::from_counts([(1u64, 1), (2, 4)]),
        )),
        1,
        2,
        2,
        -4,
        u128::from(u64::MAX) + 1
    )]
    #[case::width128(
        SignedHashedFeatureSet::from_rust(GraphSignedFeatureSet::difference(
            &GraphCountedFeatureSet::from_counts([
                ((u64::MAX as u128) + 1, 3),
                (3, 1),
            ]),
            &GraphCountedFeatureSet::from_counts([
                ((u64::MAX as u128) + 1, 1),
                (2, 4),
            ]),
        )),
        (u64::MAX as u128) + 1,
        2,
        2,
        -4,
        9
    )]
    fn test_signed_hashed_feature_set_count(
        #[case] features: SignedHashedFeatureSet,
        #[case] positive_id: u128,
        #[case] positive_count: i32,
        #[case] negative_id: u128,
        #[case] negative_count: i32,
        #[case] absent_id: u128,
    ) {
        Python::attach(|py| {
            let features = Py::new(py, features).unwrap();
            let features = features.bind(py).as_any();

            assert_eq!(
                features
                    .call_method1("count", (positive_id,))
                    .unwrap()
                    .extract::<i32>()
                    .unwrap(),
                positive_count
            );
            assert_eq!(
                features
                    .call_method1("count", (negative_id,))
                    .unwrap()
                    .extract::<i32>()
                    .unwrap(),
                negative_count
            );
            assert_eq!(
                features
                    .call_method1("count", (absent_id,))
                    .unwrap()
                    .extract::<i32>()
                    .unwrap(),
                0
            );
        });
    }

    #[rstest]
    #[case::snapshot(SignedHashedFeatureSet::from_rust(
        GraphSignedFeatureSet::difference(
            &GraphCountedFeatureSet::from_counts([(1u64, 3), (3, 1)]),
            &GraphCountedFeatureSet::from_counts([(1u64, 1), (2, 4)]),
        )
    ))]
    fn test_signed_hashed_feature_set_entries_snapshot(#[case] features: SignedHashedFeatureSet) {
        Python::attach(|py| {
            let features = Py::new(py, features).unwrap();
            let features = features.bind(py).as_any();
            let entries = features.getattr("entries").unwrap();

            entries
                .cast::<PyList>()
                .unwrap()
                .append((9u128, -2i32))
                .unwrap();

            assert_eq!(
                features
                    .getattr("entries")
                    .unwrap()
                    .extract::<Vec<(u128, i32)>>()
                    .unwrap(),
                vec![(1, 2), (2, -4), (3, 1)]
            );
        });
    }

    #[rstest]
    #[case::same_width(
        SignedHashedFeatureSet::from_rust(GraphSignedFeatureSet::difference(
            &GraphCountedFeatureSet::from_counts([(1u64, 3)]),
            &GraphCountedFeatureSet::from_counts([(1u64, 1)]),
        )),
        SignedHashedFeatureSet::from_rust(GraphSignedFeatureSet::difference(
            &GraphCountedFeatureSet::from_counts([(1u64, 3)]),
            &GraphCountedFeatureSet::from_counts([(1u64, 1)]),
        )),
        true
    )]
    #[case::different_counts(
        SignedHashedFeatureSet::from_rust(GraphSignedFeatureSet::difference(
            &GraphCountedFeatureSet::from_counts([(1u64, 3)]),
            &GraphCountedFeatureSet::from_counts([(1u64, 1)]),
        )),
        SignedHashedFeatureSet::from_rust(GraphSignedFeatureSet::difference(
            &GraphCountedFeatureSet::from_counts([(1u64, 4)]),
            &GraphCountedFeatureSet::from_counts([(1u64, 1)]),
        )),
        false
    )]
    #[case::different_widths(
        SignedHashedFeatureSet::from_rust(GraphSignedFeatureSet::difference(
            &GraphCountedFeatureSet::from_counts([(1u32, 3)]),
            &GraphCountedFeatureSet::from_counts([(1u32, 1)]),
        )),
        SignedHashedFeatureSet::from_rust(GraphSignedFeatureSet::difference(
            &GraphCountedFeatureSet::from_counts([(1u64, 3)]),
            &GraphCountedFeatureSet::from_counts([(1u64, 1)]),
        )),
        false
    )]
    fn test_signed_hashed_feature_set_eq(
        #[case] left: SignedHashedFeatureSet,
        #[case] right: SignedHashedFeatureSet,
        #[case] expected: bool,
    ) {
        assert_eq!(left == right, expected);
    }

    #[rstest]
    #[case::width32(
        RoleTaggedHashedFeatureSet::from_rust(GraphFeatureSet::from_features([
            (GraphReactionSide::Product, 2u32),
            (GraphReactionSide::Reactant, 3),
            (GraphReactionSide::Reactant, 1),
            (GraphReactionSide::Product, 1),
        ])),
        vec![
            (ReactionSide::Reactant, 1),
            (ReactionSide::Reactant, 3),
            (ReactionSide::Product, 1),
            (ReactionSide::Product, 2),
        ],
        32,
        "RoleTaggedHashedFeatureSet(ids=[(ReactionSide.Reactant, 1), (ReactionSide.Reactant, 3), (ReactionSide.Product, 1), (ReactionSide.Product, 2)], id_width=32)"
    )]
    #[case::width64(
        RoleTaggedHashedFeatureSet::from_rust(GraphFeatureSet::from_features([
            (GraphReactionSide::Product, 2u64),
            (GraphReactionSide::Reactant, 3),
            (GraphReactionSide::Reactant, 1),
            (GraphReactionSide::Product, 1),
        ])),
        vec![
            (ReactionSide::Reactant, 1),
            (ReactionSide::Reactant, 3),
            (ReactionSide::Product, 1),
            (ReactionSide::Product, 2),
        ],
        64,
        "RoleTaggedHashedFeatureSet(ids=[(ReactionSide.Reactant, 1), (ReactionSide.Reactant, 3), (ReactionSide.Product, 1), (ReactionSide.Product, 2)], id_width=64)"
    )]
    #[case::width128(
        RoleTaggedHashedFeatureSet::from_rust(GraphFeatureSet::from_features([
            (GraphReactionSide::Product, (u64::MAX as u128) + 1),
            (GraphReactionSide::Reactant, 1),
        ])),
        vec![
            (ReactionSide::Reactant, 1),
            (ReactionSide::Product, (u64::MAX as u128) + 1),
        ],
        128,
        "RoleTaggedHashedFeatureSet(ids=[(ReactionSide.Reactant, 1), (ReactionSide.Product, 18446744073709551616)], id_width=128)"
    )]
    fn test_role_tagged_hashed_feature_set_value(
        #[case] features: RoleTaggedHashedFeatureSet,
        #[case] expected_ids: Vec<(ReactionSide, u128)>,
        #[case] expected_width: u16,
        #[case] expected_repr: &str,
    ) {
        Python::attach(|py| {
            let expected = Py::new(py, features.clone()).unwrap();
            let features = Py::new(py, features).unwrap();
            let expected = expected.bind(py).as_any();
            let features = features.bind(py).as_any();

            assert!(features.eq(expected).unwrap());
            assert_eq!(
                features
                    .getattr("ids")
                    .unwrap()
                    .extract::<Vec<(ReactionSide, u128)>>()
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
                    .map(|item| item.unwrap().extract::<(ReactionSide, u128)>().unwrap())
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
    #[case::snapshot(RoleTaggedHashedFeatureSet::from_rust(
        GraphFeatureSet::from_features([
            (GraphReactionSide::Reactant, 1u64),
            (GraphReactionSide::Product, 2),
        ])
    ))]
    fn test_role_tagged_hashed_feature_set_ids_snapshot(
        #[case] features: RoleTaggedHashedFeatureSet,
    ) {
        Python::attach(|py| {
            let features = Py::new(py, features).unwrap();
            let features = features.bind(py).as_any();
            let ids = features.getattr("ids").unwrap();

            ids.cast::<PyList>()
                .unwrap()
                .append((ReactionSide::Product, 9u128))
                .unwrap();

            assert_eq!(
                features
                    .getattr("ids")
                    .unwrap()
                    .extract::<Vec<(ReactionSide, u128)>>()
                    .unwrap(),
                vec![(ReactionSide::Reactant, 1), (ReactionSide::Product, 2),]
            );
        });
    }

    #[rstest]
    #[case::same_width(
        RoleTaggedHashedFeatureSet::from_rust(GraphFeatureSet::from_features([
            (GraphReactionSide::Reactant, 1u64),
            (GraphReactionSide::Product, 2),
        ])),
        RoleTaggedHashedFeatureSet::from_rust(GraphFeatureSet::from_features([
            (GraphReactionSide::Reactant, 1u64),
            (GraphReactionSide::Product, 2),
        ])),
        true
    )]
    #[case::different_roles(
        RoleTaggedHashedFeatureSet::from_rust(GraphFeatureSet::from_features([
            (GraphReactionSide::Reactant, 1u64),
        ])),
        RoleTaggedHashedFeatureSet::from_rust(GraphFeatureSet::from_features([
            (GraphReactionSide::Product, 1u64),
        ])),
        false
    )]
    #[case::different_widths(
        RoleTaggedHashedFeatureSet::from_rust(GraphFeatureSet::from_features([
            (GraphReactionSide::Reactant, 1u32),
        ])),
        RoleTaggedHashedFeatureSet::from_rust(GraphFeatureSet::from_features([
            (GraphReactionSide::Reactant, 1u64),
        ])),
        false
    )]
    fn test_role_tagged_hashed_feature_set_eq(
        #[case] left: RoleTaggedHashedFeatureSet,
        #[case] right: RoleTaggedHashedFeatureSet,
        #[case] expected: bool,
    ) {
        assert_eq!(left == right, expected);
    }
}
