//! Read-only Python values for reaction correspondences.

use std::fmt::Debug;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use umol_graph_core::{Correspondence as GraphCoreCorrespondence, NodeId};
use umol_graph_ir::ir::{
    AromaticSystemId, AtomId, BondId, DativeBondId,
    MoleculeCorrespondence as GraphIrMoleculeCorrespondence, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};

/// An id that can cross the Python boundary as an integer index.
pub(crate) trait CorrespondenceId: Copy + Debug + Ord + From<usize> {
    fn index(self) -> usize;
}

macro_rules! correspondence_ids {
    ($($id:ty),+ $(,)?) => {
        $(
            impl CorrespondenceId for $id {
                fn index(self) -> usize {
                    self.index()
                }
            }
        )+
    };
}

correspondence_ids!(
    NodeId,
    AtomId,
    BondId,
    DativeBondId,
    AromaticSystemId,
    MulticenterBondId,
    NoncovalentBondId,
    StereoAtomId,
    StereoBondId,
);

/// A read-only partial bijection between two integer id spaces.
#[pyclass(eq, frozen, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Correspondence(GraphCoreCorrespondence<usize>);

#[pymethods]
impl Correspondence {
    /// Construct an in-range partial bijection between two integer id spaces.
    ///
    /// Raises `ValueError` when either id is outside its declared space or an id occurs in more than
    /// one pair on the same side. The indices are not bound to any particular molecule.
    #[new]
    fn new(
        matched_pairs: Vec<(usize, usize)>,
        left_count: usize,
        right_count: usize,
    ) -> PyResult<Self> {
        GraphCoreCorrespondence::new(matched_pairs, left_count, right_count)
            .map(Self)
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// Matched `(left, right)` id pairs, ordered by left id.
    #[getter]
    fn matched_pairs(&self) -> Vec<(usize, usize)> {
        self.0.matched_pairs().to_vec()
    }

    /// Size of the left id space.
    #[getter]
    fn left_count(&self) -> usize {
        self.0.left_count()
    }

    /// Size of the right id space.
    #[getter]
    fn right_count(&self) -> usize {
        self.0.right_count()
    }

    /// Left ids without a match.
    #[getter]
    fn left_unmatched(&self) -> Vec<usize> {
        self.0.left_unmatched()
    }

    /// Right ids without a match.
    #[getter]
    fn right_unmatched(&self) -> Vec<usize> {
        self.0.right_unmatched()
    }

    /// Number of matched pairs.
    fn __len__(&self) -> usize {
        self.0.matched_pair_count()
    }

    /// Right id matched to `left`, if any.
    fn right_of(&self, left: usize) -> Option<usize> {
        self.0.right_of(left)
    }

    /// Left id matched to `right`, if any.
    fn left_of(&self, right: usize) -> Option<usize> {
        self.0.left_of(right)
    }

    /// Whether every id on both sides is matched.
    fn is_total(&self) -> bool {
        self.0.is_total()
    }

    /// Invert the correspondence.
    fn reverse(&self) -> Self {
        Self(self.0.reverse())
    }

    /// Relational composition with a following correspondence.
    fn compose(&self, other: &Self) -> Self {
        Self(self.0.compose(&other.0))
    }

    /// Compose an iterable of correspondences in iteration order.
    #[staticmethod]
    fn compose_all(correspondences: &Bound<'_, PyAny>) -> PyResult<Option<Self>> {
        let correspondences = correspondences
            .try_iter()?
            .map(|item| {
                let item = item?.cast_into::<Correspondence>()?;
                Ok(item.borrow().0.clone())
            })
            .collect::<PyResult<Vec<_>>>()?;
        Ok(GraphCoreCorrespondence::compose_all(correspondences).map(Self))
    }

    fn __repr__(&self) -> String {
        self.repr()
    }
}

impl Correspondence {
    pub(crate) fn from_rust<Id: CorrespondenceId>(
        correspondence: &GraphCoreCorrespondence<Id>,
    ) -> Self {
        Self(
            GraphCoreCorrespondence::new(
                correspondence
                    .matched_pairs()
                    .iter()
                    .map(|&(left, right)| (left.index(), right.index()))
                    .collect(),
                correspondence.left_count(),
                correspondence.right_count(),
            )
            .expect("correspondence producer preserves partial-bijection invariants"),
        )
    }

    pub(crate) fn to_rust<Id: CorrespondenceId>(&self) -> GraphCoreCorrespondence<Id> {
        GraphCoreCorrespondence::new(
            self.0
                .matched_pairs()
                .iter()
                .map(|&(left, right)| (Id::from(left), Id::from(right)))
                .collect(),
            self.0.left_count(),
            self.0.right_count(),
        )
        .expect("Python correspondence preserves partial-bijection invariants")
    }

    fn repr(&self) -> String {
        format!(
            "Correspondence(matched_pairs={:?}, left_count={}, right_count={})",
            self.0.matched_pairs(),
            self.0.left_count(),
            self.0.right_count()
        )
    }
}

/// A read-only correspondence across every molecule entity family.
#[pyclass(eq, frozen, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MoleculeCorrespondence(GraphIrMoleculeCorrespondence);

#[pymethods]
impl MoleculeCorrespondence {
    /// Atom correspondence.
    #[getter]
    fn atoms(&self) -> Correspondence {
        Correspondence::from_rust(self.0.atoms())
    }

    /// Bond correspondence.
    #[getter]
    fn bonds(&self) -> Correspondence {
        Correspondence::from_rust(self.0.bonds())
    }

    /// Dative-bond correspondence.
    #[getter]
    fn dative_bonds(&self) -> Correspondence {
        Correspondence::from_rust(self.0.dative_bonds())
    }

    /// Aromatic-system correspondence.
    #[getter]
    fn aromatic_systems(&self) -> Correspondence {
        Correspondence::from_rust(self.0.aromatic_systems())
    }

    /// Multicenter-bond correspondence.
    #[getter]
    fn multicenter_bonds(&self) -> Correspondence {
        Correspondence::from_rust(self.0.multicenter_bonds())
    }

    /// Noncovalent-bond correspondence.
    #[getter]
    fn noncovalent_bonds(&self) -> Correspondence {
        Correspondence::from_rust(self.0.noncovalent_bonds())
    }

    /// Stereo-atom correspondence.
    #[getter]
    fn stereo_atoms(&self) -> Correspondence {
        Correspondence::from_rust(self.0.stereo_atoms())
    }

    /// Stereo-bond correspondence.
    #[getter]
    fn stereo_bonds(&self) -> Correspondence {
        Correspondence::from_rust(self.0.stereo_bonds())
    }

    /// Whether every id in every entity family is matched.
    fn is_total(&self) -> bool {
        self.0.is_total()
    }

    /// Invert every per-family correspondence.
    fn reverse(&self) -> Self {
        Self(self.0.reverse())
    }

    /// Relational composition with a following molecule correspondence.
    fn compose(&self, other: &Self) -> Self {
        Self(self.0.compose(&other.0))
    }

    /// Compose an iterable of molecule correspondences in iteration order.
    #[staticmethod]
    fn compose_all(correspondences: &Bound<'_, PyAny>) -> PyResult<Option<Self>> {
        let correspondences = correspondences
            .try_iter()?
            .map(|item| {
                let item = item?.cast_into::<MoleculeCorrespondence>()?;
                Ok(item.borrow().0.clone())
            })
            .collect::<PyResult<Vec<_>>>()?;
        Ok(GraphIrMoleculeCorrespondence::compose_all(correspondences).map(Self))
    }

    fn __repr__(&self) -> String {
        format!(
            concat!(
                "MoleculeCorrespondence(",
                "atoms={}, bonds={}, dative_bonds={}, aromatic_systems={}, ",
                "multicenter_bonds={}, noncovalent_bonds={}, stereo_atoms={}, ",
                "stereo_bonds={})"
            ),
            self.atoms().repr(),
            self.bonds().repr(),
            self.dative_bonds().repr(),
            self.aromatic_systems().repr(),
            self.multicenter_bonds().repr(),
            self.noncovalent_bonds().repr(),
            self.stereo_atoms().repr(),
            self.stereo_bonds().repr(),
        )
    }
}

impl MoleculeCorrespondence {
    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion API for molecule correspondence values"
    )]
    pub(crate) fn from_rust(correspondence: GraphIrMoleculeCorrespondence) -> Self {
        Self(correspondence)
    }

    pub(crate) fn to_rust(&self) -> &GraphIrMoleculeCorrespondence {
        &self.0
    }
}

#[cfg(test)]
mod tests {
    use pyo3::types::PyList;
    use rstest::{fixture, rstest};

    use super::*;

    #[fixture]
    fn molecule_correspondence() -> GraphIrMoleculeCorrespondence {
        GraphIrMoleculeCorrespondence::new(
            GraphCoreCorrespondence::new(vec![(AtomId(0), AtomId(1))], 2, 3)
                .expect("correspondence producer preserves partial-bijection invariants"),
            GraphCoreCorrespondence::new(vec![(BondId(0), BondId(2))], 1, 3)
                .expect("correspondence producer preserves partial-bijection invariants"),
            GraphCoreCorrespondence::new(vec![(DativeBondId(1), DativeBondId(0))], 2, 1)
                .expect("correspondence producer preserves partial-bijection invariants"),
            GraphCoreCorrespondence::new(vec![], 1, 2)
                .expect("correspondence producer preserves partial-bijection invariants"),
            GraphCoreCorrespondence::new(vec![(MulticenterBondId(0), MulticenterBondId(0))], 1, 1)
                .expect("correspondence producer preserves partial-bijection invariants"),
            GraphCoreCorrespondence::new(vec![(NoncovalentBondId(0), NoncovalentBondId(1))], 2, 2)
                .expect("correspondence producer preserves partial-bijection invariants"),
            GraphCoreCorrespondence::new(
                vec![
                    (StereoAtomId(0), StereoAtomId(0)),
                    (StereoAtomId(1), StereoAtomId(1)),
                ],
                2,
                2,
            )
            .expect("correspondence producer preserves partial-bijection invariants"),
            GraphCoreCorrespondence::new(vec![(StereoBondId(0), StereoBondId(1))], 1, 2)
                .expect("correspondence producer preserves partial-bijection invariants"),
        )
    }

    #[rstest]
    #[case::empty(
        vec![],
        2,
        3,
        Correspondence(GraphCoreCorrespondence::new(vec![], 2, 3).expect("correspondence producer preserves partial-bijection invariants")),
    )]
    #[case::partial(
        vec![(0, 2)],
        2,
        3,
        Correspondence(GraphCoreCorrespondence::new(vec![(0, 2)], 2, 3).expect("correspondence producer preserves partial-bijection invariants")),
    )]
    #[case::unsorted(
        vec![(2, 0), (0, 2)],
        3,
        3,
        Correspondence(GraphCoreCorrespondence::new(vec![(0, 2), (2, 0)], 3, 3).expect("correspondence producer preserves partial-bijection invariants")),
    )]
    fn test_correspondence_new(
        #[case] matched_pairs: Vec<(usize, usize)>,
        #[case] left_count: usize,
        #[case] right_count: usize,
        #[case] expected: Correspondence,
    ) {
        assert_eq!(
            Correspondence::new(matched_pairs, left_count, right_count).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::left_out_of_range(
        vec![(2, 0)],
        2,
        1,
        "left id 2 is out of range for 2 entries",
    )]
    #[case::right_out_of_range(
        vec![(0, 1)],
        1,
        1,
        "right id 1 is out of range for 1 entries",
    )]
    #[case::duplicate_left(
        vec![(0, 0), (0, 1)],
        1,
        2,
        "left id 0 occurs more than once",
    )]
    #[case::duplicate_right(
        vec![(0, 0), (1, 0)],
        2,
        1,
        "right id 0 occurs more than once",
    )]
    fn test_correspondence_new_error(
        #[case] matched_pairs: Vec<(usize, usize)>,
        #[case] left_count: usize,
        #[case] right_count: usize,
        #[case] message: &str,
    ) {
        Python::attach(|py| {
            let error = Correspondence::new(matched_pairs, left_count, right_count).unwrap_err();

            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(error.value(py).to_string(), message);
        });
    }

    #[rstest]
    #[case::empty(Correspondence(GraphCoreCorrespondence::new(vec![], 2, 3).expect("correspondence producer preserves partial-bijection invariants")))]
    #[case::partial(Correspondence(GraphCoreCorrespondence::new(vec![(0, 2)], 2, 3).expect("correspondence producer preserves partial-bijection invariants")))]
    fn test_correspondence_to_rust(#[case] correspondence: Correspondence) {
        assert_eq!(
            correspondence.to_rust::<AtomId>(),
            GraphCoreCorrespondence::new(
                correspondence
                    .matched_pairs()
                    .into_iter()
                    .map(|(left, right)| (AtomId::from(left), AtomId::from(right)))
                    .collect(),
                correspondence.left_count(),
                correspondence.right_count(),
            )
            .expect("correspondence producer preserves partial-bijection invariants")
        );
    }

    #[rstest]
    #[case::empty(
        GraphCoreCorrespondence::new(vec![], 2, 3).expect("correspondence producer preserves partial-bijection invariants"),
        Correspondence(GraphCoreCorrespondence::new(vec![], 2, 3).expect("correspondence producer preserves partial-bijection invariants")),
    )]
    #[case::partial(
        GraphCoreCorrespondence::new(vec![(NodeId(0), NodeId(2))], 2, 3).expect("correspondence producer preserves partial-bijection invariants"),
        Correspondence(GraphCoreCorrespondence::new(vec![(0, 2)], 2, 3).expect("correspondence producer preserves partial-bijection invariants")),
    )]
    #[case::total(
        GraphCoreCorrespondence::new(
            vec![(NodeId(0), NodeId(1)), (NodeId(1), NodeId(0))],
            2,
            2,
        ).expect("correspondence producer preserves partial-bijection invariants"),
        Correspondence(GraphCoreCorrespondence::new(
            vec![(0, 1), (1, 0)],
            2,
            2,
        ).expect("correspondence producer preserves partial-bijection invariants")),
    )]
    #[case::unsorted(
        GraphCoreCorrespondence::new(
            vec![(NodeId(2), NodeId(0)), (NodeId(0), NodeId(2))],
            3,
            3,
        ).expect("correspondence producer preserves partial-bijection invariants"),
        Correspondence(GraphCoreCorrespondence::new(
            vec![(0, 2), (2, 0)],
            3,
            3,
        ).expect("correspondence producer preserves partial-bijection invariants")),
    )]
    fn test_correspondence_from_rust(
        #[case] correspondence: GraphCoreCorrespondence<NodeId>,
        #[case] expected: Correspondence,
    ) {
        assert_eq!(Correspondence::from_rust(&correspondence), expected);
    }

    #[rstest]
    #[case::empty(
        GraphCoreCorrespondence::new(vec![], 2, 3).expect("correspondence producer preserves partial-bijection invariants"),
        vec![],
        2,
        3,
        vec![0, 1],
        vec![0, 1, 2],
        0,
    )]
    #[case::partial(
        GraphCoreCorrespondence::new(vec![(NodeId(0), NodeId(2))], 2, 3).expect("correspondence producer preserves partial-bijection invariants"),
        vec![(0, 2)],
        2,
        3,
        vec![1],
        vec![0, 1],
        1,
    )]
    #[case::total(
        GraphCoreCorrespondence::new(
            vec![(NodeId(0), NodeId(1)), (NodeId(1), NodeId(0))],
            2,
            2,
        ).expect("correspondence producer preserves partial-bijection invariants"),
        vec![(0, 1), (1, 0)],
        2,
        2,
        vec![],
        vec![],
        2,
    )]
    fn test_correspondence_accessors(
        #[case] correspondence: GraphCoreCorrespondence<NodeId>,
        #[case] matched_pairs: Vec<(usize, usize)>,
        #[case] left_count: usize,
        #[case] right_count: usize,
        #[case] left_unmatched: Vec<usize>,
        #[case] right_unmatched: Vec<usize>,
        #[case] matched_pair_count: usize,
    ) {
        let correspondence = Correspondence::from_rust(&correspondence);
        assert_eq!(correspondence.matched_pairs(), matched_pairs);
        assert_eq!(correspondence.left_count(), left_count);
        assert_eq!(correspondence.right_count(), right_count);
        assert_eq!(correspondence.left_unmatched(), left_unmatched);
        assert_eq!(correspondence.right_unmatched(), right_unmatched);
        assert_eq!(correspondence.__len__(), matched_pair_count);
    }

    #[rstest]
    #[case::empty(
        GraphCoreCorrespondence::new(vec![], 2, 3).expect("correspondence producer preserves partial-bijection invariants"),
        0,
        None,
        0,
        None,
    )]
    #[case::partial(
        GraphCoreCorrespondence::new(vec![(NodeId(0), NodeId(2))], 2, 3).expect("correspondence producer preserves partial-bijection invariants"),
        0,
        Some(2),
        2,
        Some(0),
    )]
    #[case::total(
        GraphCoreCorrespondence::new(
            vec![(NodeId(0), NodeId(1)), (NodeId(1), NodeId(0))],
            2,
            2,
        ).expect("correspondence producer preserves partial-bijection invariants"),
        0,
        Some(1),
        1,
        Some(0),
    )]
    fn test_correspondence_lookup(
        #[case] correspondence: GraphCoreCorrespondence<NodeId>,
        #[case] left: usize,
        #[case] right_of_left: Option<usize>,
        #[case] right: usize,
        #[case] left_of_right: Option<usize>,
    ) {
        let correspondence = Correspondence::from_rust(&correspondence);
        assert_eq!(correspondence.right_of(left), right_of_left);
        assert_eq!(correspondence.left_of(right), left_of_right);
    }

    #[rstest]
    #[case::empty(GraphCoreCorrespondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"), true)]
    #[case::partial(
        GraphCoreCorrespondence::new(vec![(NodeId(0), NodeId(2))], 2, 3).expect("correspondence producer preserves partial-bijection invariants"),
        false,
    )]
    #[case::total(
        GraphCoreCorrespondence::new(
            vec![(NodeId(0), NodeId(1)), (NodeId(1), NodeId(0))],
            2,
            2,
        ).expect("correspondence producer preserves partial-bijection invariants"),
        true,
    )]
    fn test_correspondence_is_total(
        #[case] correspondence: GraphCoreCorrespondence<NodeId>,
        #[case] expected: bool,
    ) {
        assert_eq!(
            Correspondence::from_rust(&correspondence).is_total(),
            expected
        );
    }

    #[rstest]
    #[case::partial(
        Correspondence(GraphCoreCorrespondence::new(vec![(0, 2), (2, 0)], 3, 4).expect("correspondence producer preserves partial-bijection invariants")),
        Correspondence(GraphCoreCorrespondence::new(vec![(0, 2), (2, 0)], 4, 3).expect("correspondence producer preserves partial-bijection invariants")),
    )]
    #[case::total(
        Correspondence(GraphCoreCorrespondence::new(vec![(0, 1), (1, 0)], 2, 2).expect("correspondence producer preserves partial-bijection invariants")),
        Correspondence(GraphCoreCorrespondence::new(vec![(0, 1), (1, 0)], 2, 2).expect("correspondence producer preserves partial-bijection invariants")),
    )]
    fn test_correspondence_reverse(
        #[case] correspondence: Correspondence,
        #[case] expected: Correspondence,
    ) {
        assert_eq!(correspondence.reverse(), expected);
    }

    #[rstest]
    #[case::ordinary(
        Correspondence(GraphCoreCorrespondence::new(
            vec![(0, 1), (1, 0)],
            2,
            2,
        ).expect("correspondence producer preserves partial-bijection invariants")),
        Correspondence(GraphCoreCorrespondence::new(
            vec![(0, 2), (1, 1)],
            2,
            3,
        ).expect("correspondence producer preserves partial-bijection invariants")),
        Correspondence(GraphCoreCorrespondence::new(
            vec![(0, 1), (1, 2)],
            2,
            3,
        ).expect("correspondence producer preserves partial-bijection invariants")),
    )]
    #[case::mismatched_intermediate(
        Correspondence(GraphCoreCorrespondence::new(
            vec![(0, 2), (1, 0)],
            2,
            3,
        ).expect("correspondence producer preserves partial-bijection invariants")),
        Correspondence(GraphCoreCorrespondence::new(vec![(0, 4)], 1, 5).expect("correspondence producer preserves partial-bijection invariants")),
        Correspondence(GraphCoreCorrespondence::new(vec![(1, 4)], 2, 5).expect("correspondence producer preserves partial-bijection invariants")),
    )]
    fn test_correspondence_compose(
        #[case] left: Correspondence,
        #[case] right: Correspondence,
        #[case] expected: Correspondence,
    ) {
        assert_eq!(left.compose(&right), expected);
    }

    #[rstest]
    #[case::empty(vec![], None)]
    #[case::singleton(
        vec![Correspondence(GraphCoreCorrespondence::new(vec![(0, 1)], 1, 2).expect("correspondence producer preserves partial-bijection invariants"))],
        Some(Correspondence(GraphCoreCorrespondence::new(vec![(0, 1)], 1, 2).expect("correspondence producer preserves partial-bijection invariants"))),
    )]
    #[case::multiple(
        vec![
            Correspondence(GraphCoreCorrespondence::new(
                vec![(0, 1), (1, 0)],
                2,
                2,
            ).expect("correspondence producer preserves partial-bijection invariants")),
            Correspondence(GraphCoreCorrespondence::new(
                vec![(0, 2), (1, 1)],
                2,
                3,
            ).expect("correspondence producer preserves partial-bijection invariants")),
            Correspondence(GraphCoreCorrespondence::new(
                vec![(1, 0), (2, 1)],
                3,
                2,
            ).expect("correspondence producer preserves partial-bijection invariants")),
        ],
        Some(Correspondence(GraphCoreCorrespondence::new(
            vec![(0, 0), (1, 1)],
            2,
            2,
        ).expect("correspondence producer preserves partial-bijection invariants"))),
    )]
    fn test_correspondence_compose_all(
        #[case] correspondences: Vec<Correspondence>,
        #[case] expected: Option<Correspondence>,
    ) {
        Python::attach(|py| {
            let correspondences = correspondences
                .into_iter()
                .map(|correspondence| Py::new(py, correspondence).unwrap())
                .collect::<Vec<_>>();
            let correspondences = PyList::new(py, correspondences).unwrap();

            assert_eq!(
                Correspondence::compose_all(correspondences.as_any()).unwrap(),
                expected
            );
        });
    }

    #[rstest]
    fn test_correspondence_value() {
        let correspondence = Correspondence::from_rust(
            &GraphCoreCorrespondence::new(vec![(NodeId(0), NodeId(2))], 2, 3)
                .expect("correspondence producer preserves partial-bijection invariants"),
        );
        let mut matched_pairs = correspondence.matched_pairs();
        let mut left_unmatched = correspondence.left_unmatched();
        let mut right_unmatched = correspondence.right_unmatched();

        matched_pairs.push((1, 1));
        left_unmatched.clear();
        right_unmatched.clear();

        assert_eq!(correspondence.matched_pairs(), vec![(0, 2)]);
        assert_eq!(correspondence.left_unmatched(), vec![1]);
        assert_eq!(correspondence.right_unmatched(), vec![0, 1]);
        assert_eq!(
            correspondence.__repr__(),
            "Correspondence(matched_pairs=[(0, 2)], left_count=2, right_count=3)"
        );
        assert_eq!(
            correspondence,
            Correspondence::from_rust(
                &GraphCoreCorrespondence::new(vec![(NodeId(0), NodeId(2))], 2, 3,)
                    .expect("correspondence producer preserves partial-bijection invariants")
            )
        );
        assert_ne!(
            correspondence,
            Correspondence::from_rust(
                &GraphCoreCorrespondence::new(vec![(NodeId(0), NodeId(2))], 2, 4,)
                    .expect("correspondence producer preserves partial-bijection invariants")
            )
        );
    }

    #[rstest]
    fn test_molecule_correspondence_accessors(
        molecule_correspondence: GraphIrMoleculeCorrespondence,
    ) {
        let correspondence = MoleculeCorrespondence::from_rust(molecule_correspondence);

        assert_eq!(
            correspondence.atoms(),
            Correspondence(
                GraphCoreCorrespondence::new(vec![(0, 1)], 2, 3)
                    .expect("correspondence producer preserves partial-bijection invariants")
            )
        );
        assert_eq!(
            correspondence.bonds(),
            Correspondence(
                GraphCoreCorrespondence::new(vec![(0, 2)], 1, 3)
                    .expect("correspondence producer preserves partial-bijection invariants")
            )
        );
        assert_eq!(
            correspondence.dative_bonds(),
            Correspondence(
                GraphCoreCorrespondence::new(vec![(1, 0)], 2, 1)
                    .expect("correspondence producer preserves partial-bijection invariants")
            )
        );
        assert_eq!(
            correspondence.aromatic_systems(),
            Correspondence(
                GraphCoreCorrespondence::new(vec![], 1, 2)
                    .expect("correspondence producer preserves partial-bijection invariants")
            )
        );
        assert_eq!(
            correspondence.multicenter_bonds(),
            Correspondence(
                GraphCoreCorrespondence::new(vec![(0, 0)], 1, 1)
                    .expect("correspondence producer preserves partial-bijection invariants")
            )
        );
        assert_eq!(
            correspondence.noncovalent_bonds(),
            Correspondence(
                GraphCoreCorrespondence::new(vec![(0, 1)], 2, 2)
                    .expect("correspondence producer preserves partial-bijection invariants")
            )
        );
        assert_eq!(
            correspondence.stereo_atoms(),
            Correspondence(
                GraphCoreCorrespondence::new(vec![(0, 0), (1, 1)], 2, 2,)
                    .expect("correspondence producer preserves partial-bijection invariants")
            )
        );
        assert_eq!(
            correspondence.stereo_bonds(),
            Correspondence(
                GraphCoreCorrespondence::new(vec![(0, 1)], 1, 2)
                    .expect("correspondence producer preserves partial-bijection invariants")
            )
        );
    }

    #[rstest]
    #[case::empty_spaces(GraphIrMoleculeCorrespondence::new(
        GraphCoreCorrespondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
        GraphCoreCorrespondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
        GraphCoreCorrespondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
        GraphCoreCorrespondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
        GraphCoreCorrespondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
        GraphCoreCorrespondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
        GraphCoreCorrespondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
        GraphCoreCorrespondence::new(vec![], 0, 0).expect("correspondence producer preserves partial-bijection invariants"),
    ))]
    fn test_molecule_correspondence_is_total(
        #[case] correspondence: GraphIrMoleculeCorrespondence,
    ) {
        assert!(MoleculeCorrespondence::from_rust(correspondence).is_total());
    }

    #[rstest]
    fn test_molecule_correspondence_is_total_partial(
        molecule_correspondence: GraphIrMoleculeCorrespondence,
    ) {
        assert!(!MoleculeCorrespondence::from_rust(molecule_correspondence).is_total());
    }

    #[rstest]
    fn test_molecule_correspondence_reverse(
        molecule_correspondence: GraphIrMoleculeCorrespondence,
    ) {
        assert_eq!(
            MoleculeCorrespondence::from_rust(molecule_correspondence).reverse(),
            MoleculeCorrespondence::from_rust(GraphIrMoleculeCorrespondence::new(
                GraphCoreCorrespondence::new(vec![(AtomId(1), AtomId(0))], 3, 2)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                GraphCoreCorrespondence::new(vec![(BondId(2), BondId(0))], 3, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                GraphCoreCorrespondence::new(vec![(DativeBondId(0), DativeBondId(1))], 1, 2,)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                GraphCoreCorrespondence::new(vec![], 2, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                GraphCoreCorrespondence::new(
                    vec![(MulticenterBondId(0), MulticenterBondId(0))],
                    1,
                    1,
                )
                .expect("correspondence producer preserves partial-bijection invariants"),
                GraphCoreCorrespondence::new(
                    vec![(NoncovalentBondId(1), NoncovalentBondId(0))],
                    2,
                    2,
                )
                .expect("correspondence producer preserves partial-bijection invariants"),
                GraphCoreCorrespondence::new(
                    vec![
                        (StereoAtomId(0), StereoAtomId(0)),
                        (StereoAtomId(1), StereoAtomId(1)),
                    ],
                    2,
                    2,
                )
                .expect("correspondence producer preserves partial-bijection invariants"),
                GraphCoreCorrespondence::new(vec![(StereoBondId(1), StereoBondId(0))], 2, 1,)
                    .expect("correspondence producer preserves partial-bijection invariants"),
            ))
        );
    }

    #[rstest]
    fn test_molecule_correspondence_compose(
        molecule_correspondence: GraphIrMoleculeCorrespondence,
    ) {
        let left = MoleculeCorrespondence::from_rust(molecule_correspondence.clone());
        let right = left.reverse();

        assert_eq!(
            left.compose(&right),
            MoleculeCorrespondence::from_rust(GraphIrMoleculeCorrespondence::new(
                GraphCoreCorrespondence::new(vec![(AtomId(0), AtomId(0))], 2, 2)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                GraphCoreCorrespondence::new(vec![(BondId(0), BondId(0))], 1, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                GraphCoreCorrespondence::new(vec![(DativeBondId(1), DativeBondId(1))], 2, 2,)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                GraphCoreCorrespondence::new(vec![], 1, 1)
                    .expect("correspondence producer preserves partial-bijection invariants"),
                GraphCoreCorrespondence::new(
                    vec![(MulticenterBondId(0), MulticenterBondId(0))],
                    1,
                    1,
                )
                .expect("correspondence producer preserves partial-bijection invariants"),
                GraphCoreCorrespondence::new(
                    vec![(NoncovalentBondId(0), NoncovalentBondId(0))],
                    2,
                    2,
                )
                .expect("correspondence producer preserves partial-bijection invariants"),
                GraphCoreCorrespondence::new(
                    vec![
                        (StereoAtomId(0), StereoAtomId(0)),
                        (StereoAtomId(1), StereoAtomId(1)),
                    ],
                    2,
                    2,
                )
                .expect("correspondence producer preserves partial-bijection invariants"),
                GraphCoreCorrespondence::new(vec![(StereoBondId(0), StereoBondId(0))], 1, 1,)
                    .expect("correspondence producer preserves partial-bijection invariants"),
            ))
        );
    }

    #[rstest]
    fn test_molecule_correspondence_compose_all(
        molecule_correspondence: GraphIrMoleculeCorrespondence,
    ) {
        Python::attach(|py| {
            let left = MoleculeCorrespondence::from_rust(molecule_correspondence);
            let right = left.reverse();
            let expected = left.compose(&right);
            let correspondences = PyList::new(
                py,
                [Py::new(py, left).unwrap(), Py::new(py, right).unwrap()],
            )
            .unwrap();

            assert_eq!(
                MoleculeCorrespondence::compose_all(correspondences.as_any()).unwrap(),
                Some(expected)
            );
        });
    }

    #[rstest]
    fn test_molecule_correspondence_compose_all_empty() {
        Python::attach(|py| {
            let correspondences = PyList::empty(py);

            assert_eq!(
                MoleculeCorrespondence::compose_all(correspondences.as_any()).unwrap(),
                None
            );
        });
    }

    #[rstest]
    fn test_molecule_correspondence_value(molecule_correspondence: GraphIrMoleculeCorrespondence) {
        Python::attach(|py| {
            let correspondence = Py::new(
                py,
                MoleculeCorrespondence::from_rust(molecule_correspondence.clone()),
            )
            .unwrap();
            let first = correspondence.bind(py).getattr("atoms").unwrap();
            let second = correspondence.bind(py).getattr("atoms").unwrap();
            let first_matched_pairs = first
                .getattr("matched_pairs")
                .unwrap()
                .cast_into::<PyList>()
                .unwrap();
            first_matched_pairs.append((1, 2)).unwrap();

            assert!(!first.is(&second));
            assert_eq!(
                second
                    .getattr("matched_pairs")
                    .unwrap()
                    .extract::<Vec<(usize, usize)>>()
                    .unwrap(),
                vec![(0, 1)]
            );
        });

        let correspondence = MoleculeCorrespondence::from_rust(molecule_correspondence.clone());
        assert_eq!(
            correspondence,
            MoleculeCorrespondence::from_rust(molecule_correspondence)
        );
        assert_eq!(
            correspondence.__repr__(),
            concat!(
                "MoleculeCorrespondence(",
                "atoms=Correspondence(matched_pairs=[(0, 1)], left_count=2, right_count=3), ",
                "bonds=Correspondence(matched_pairs=[(0, 2)], left_count=1, right_count=3), ",
                "dative_bonds=Correspondence(matched_pairs=[(1, 0)], left_count=2, right_count=1), ",
                "aromatic_systems=Correspondence(matched_pairs=[], left_count=1, right_count=2), ",
                "multicenter_bonds=Correspondence(matched_pairs=[(0, 0)], left_count=1, right_count=1), ",
                "noncovalent_bonds=Correspondence(matched_pairs=[(0, 1)], left_count=2, right_count=2), ",
                "stereo_atoms=Correspondence(matched_pairs=[(0, 0), (1, 1)], left_count=2, right_count=2), ",
                "stereo_bonds=Correspondence(matched_pairs=[(0, 1)], left_count=1, right_count=2))"
            )
        );
    }
}
