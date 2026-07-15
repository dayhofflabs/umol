//! Read-only Python values for reaction correspondences and matching algorithms.

use pyo3::prelude::*;
use umol_ast::ast::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use umol_graph_core::{Correspondence as RustCorrespondence, NodeId};

/// An id that can cross the Python boundary as an integer index.
#[allow(dead_code, reason = "consumed by the S7a.2 molecule correspondence")]
pub(crate) trait CorrespondenceId: Copy + Ord + From<usize> {
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
pub struct Correspondence {
    mates: Vec<(usize, usize)>,
    left_count: usize,
    right_count: usize,
}

#[pymethods]
impl Correspondence {
    /// Mated `(left, right)` id pairs, ordered by left id.
    #[getter]
    fn mates(&self) -> Vec<(usize, usize)> {
        self.mates.clone()
    }

    /// Size of the left id space.
    #[getter]
    fn left_count(&self) -> usize {
        self.left_count
    }

    /// Size of the right id space.
    #[getter]
    fn right_count(&self) -> usize {
        self.right_count
    }

    /// Left ids without a mate.
    #[getter]
    fn left_exposed(&self) -> Vec<usize> {
        exposed_ids(self.left_count, self.mates.iter().map(|&(left, _)| left))
    }

    /// Right ids without a mate.
    #[getter]
    fn right_exposed(&self) -> Vec<usize> {
        exposed_ids(self.right_count, self.mates.iter().map(|&(_, right)| right))
    }

    fn __repr__(&self) -> String {
        self.repr()
    }
}

impl Correspondence {
    #[allow(dead_code, reason = "consumed by the S7a.2 molecule correspondence")]
    pub(crate) fn from_rust<Id: CorrespondenceId>(correspondence: &RustCorrespondence<Id>) -> Self {
        Self {
            mates: correspondence
                .mates()
                .iter()
                .map(|&(left, right)| (left.index(), right.index()))
                .collect(),
            left_count: correspondence.left_count(),
            right_count: correspondence.right_count(),
        }
    }

    fn repr(&self) -> String {
        format!(
            "Correspondence(mates={:?}, left_count={}, right_count={})",
            self.mates, self.left_count, self.right_count
        )
    }
}

fn exposed_ids(count: usize, mated: impl Iterator<Item = usize>) -> Vec<usize> {
    let mut present = vec![false; count];
    for id in mated {
        present[id] = true;
    }
    present
        .into_iter()
        .enumerate()
        .filter_map(|(id, is_mated)| (!is_mated).then_some(id))
        .collect()
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::empty(
        RustCorrespondence::new(vec![], 2, 3),
        Correspondence { mates: vec![], left_count: 2, right_count: 3 },
    )]
    #[case::partial(
        RustCorrespondence::new(vec![(NodeId(0), NodeId(2))], 2, 3),
        Correspondence { mates: vec![(0, 2)], left_count: 2, right_count: 3 },
    )]
    #[case::total(
        RustCorrespondence::new(
            vec![(NodeId(0), NodeId(1)), (NodeId(1), NodeId(0))],
            2,
            2,
        ),
        Correspondence {
            mates: vec![(0, 1), (1, 0)],
            left_count: 2,
            right_count: 2,
        },
    )]
    #[case::unsorted(
        RustCorrespondence::new(
            vec![(NodeId(2), NodeId(0)), (NodeId(0), NodeId(2))],
            3,
            3,
        ),
        Correspondence {
            mates: vec![(0, 2), (2, 0)],
            left_count: 3,
            right_count: 3,
        },
    )]
    fn test_correspondence_from_rust(
        #[case] correspondence: RustCorrespondence<NodeId>,
        #[case] expected: Correspondence,
    ) {
        assert_eq!(Correspondence::from_rust(&correspondence), expected);
    }

    #[rstest]
    #[case::empty(
        RustCorrespondence::new(vec![], 2, 3),
        vec![],
        2,
        3,
        vec![0, 1],
        vec![0, 1, 2],
    )]
    #[case::partial(
        RustCorrespondence::new(vec![(NodeId(0), NodeId(2))], 2, 3),
        vec![(0, 2)],
        2,
        3,
        vec![1],
        vec![0, 1],
    )]
    #[case::total(
        RustCorrespondence::new(
            vec![(NodeId(0), NodeId(1)), (NodeId(1), NodeId(0))],
            2,
            2,
        ),
        vec![(0, 1), (1, 0)],
        2,
        2,
        vec![],
        vec![],
    )]
    fn test_correspondence_accessors(
        #[case] correspondence: RustCorrespondence<NodeId>,
        #[case] mates: Vec<(usize, usize)>,
        #[case] left_count: usize,
        #[case] right_count: usize,
        #[case] left_exposed: Vec<usize>,
        #[case] right_exposed: Vec<usize>,
    ) {
        let correspondence = Correspondence::from_rust(&correspondence);
        assert_eq!(correspondence.mates(), mates);
        assert_eq!(correspondence.left_count(), left_count);
        assert_eq!(correspondence.right_count(), right_count);
        assert_eq!(correspondence.left_exposed(), left_exposed);
        assert_eq!(correspondence.right_exposed(), right_exposed);
    }

    #[rstest]
    fn test_correspondence_value() {
        let correspondence =
            Correspondence::from_rust(&RustCorrespondence::new(vec![(NodeId(0), NodeId(2))], 2, 3));
        let mut mates = correspondence.mates();
        let mut left_exposed = correspondence.left_exposed();
        let mut right_exposed = correspondence.right_exposed();

        mates.push((1, 1));
        left_exposed.clear();
        right_exposed.clear();

        assert_eq!(correspondence.mates(), vec![(0, 2)]);
        assert_eq!(correspondence.left_exposed(), vec![1]);
        assert_eq!(correspondence.right_exposed(), vec![0, 1]);
        assert_eq!(
            correspondence.__repr__(),
            "Correspondence(mates=[(0, 2)], left_count=2, right_count=3)"
        );
        assert_eq!(
            correspondence,
            Correspondence::from_rust(
                &RustCorrespondence::new(vec![(NodeId(0), NodeId(2))], 2, 3,)
            )
        );
        assert_ne!(
            correspondence,
            Correspondence::from_rust(
                &RustCorrespondence::new(vec![(NodeId(0), NodeId(2))], 2, 4,)
            )
        );
    }
}
