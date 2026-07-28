//! Read-only Python values for reaction correspondences.

use std::collections::HashSet;

use pyo3::prelude::*;
use umol_ast::ast::{
    AromaticSystemId, BondId, DativeBondId, MoleculeCorrespondence as AstMoleculeCorrespondence,
    MulticenterBondId, NoncovalentBondId, StereoAtomId, StereoBondId,
};
use umol_graph_core::{Correspondence as GraphCoreCorrespondence, NodeId};

/// An id that can cross the Python boundary as an integer index.
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
    matched_pairs: Vec<(usize, usize)>,
    left_count: usize,
    right_count: usize,
}

#[pymethods]
impl Correspondence {
    /// Matched `(left, right)` id pairs, ordered by left id.
    #[getter]
    fn matched_pairs(&self) -> Vec<(usize, usize)> {
        self.matched_pairs.clone()
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

    /// Left ids without a match.
    #[getter]
    fn left_unmatched(&self) -> Vec<usize> {
        unmatched_ids(
            self.left_count,
            self.matched_pairs.iter().map(|&(left, _)| left),
        )
    }

    /// Right ids without a match.
    #[getter]
    fn right_unmatched(&self) -> Vec<usize> {
        unmatched_ids(
            self.right_count,
            self.matched_pairs.iter().map(|&(_, right)| right),
        )
    }

    fn __repr__(&self) -> String {
        self.repr()
    }
}

impl Correspondence {
    pub(crate) fn from_rust<Id: CorrespondenceId>(
        correspondence: &GraphCoreCorrespondence<Id>,
    ) -> Self {
        Self {
            matched_pairs: correspondence
                .matched_pairs()
                .iter()
                .map(|&(left, right)| (left.index(), right.index()))
                .collect(),
            left_count: correspondence.left_count(),
            right_count: correspondence.right_count(),
        }
    }

    fn repr(&self) -> String {
        format!(
            "Correspondence(matched_pairs={:?}, left_count={}, right_count={})",
            self.matched_pairs, self.left_count, self.right_count
        )
    }
}

/// A read-only correspondence across every molecule entity family.
#[pyclass(eq, frozen, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MoleculeCorrespondence(AstMoleculeCorrespondence);

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
    pub(crate) fn from_rust(correspondence: AstMoleculeCorrespondence) -> Self {
        Self(correspondence)
    }
}

fn unmatched_ids(count: usize, matched: impl Iterator<Item = usize>) -> Vec<usize> {
    let present: HashSet<_> = matched.collect();
    (0..count).filter(|id| !present.contains(id)).collect()
}

#[cfg(test)]
mod tests {
    use pyo3::types::PyList;
    use rstest::{fixture, rstest};

    use super::*;

    #[fixture]
    fn molecule_correspondence() -> AstMoleculeCorrespondence {
        AstMoleculeCorrespondence::new(
            GraphCoreCorrespondence::new(vec![(NodeId(0), NodeId(1))], 2, 3),
            GraphCoreCorrespondence::new(vec![(BondId(0), BondId(2))], 1, 3),
            GraphCoreCorrespondence::new(vec![(DativeBondId(1), DativeBondId(0))], 2, 1),
            GraphCoreCorrespondence::new(vec![], 1, 2),
            GraphCoreCorrespondence::new(vec![(MulticenterBondId(0), MulticenterBondId(0))], 1, 1),
            GraphCoreCorrespondence::new(vec![(NoncovalentBondId(0), NoncovalentBondId(1))], 2, 2),
            GraphCoreCorrespondence::new(
                vec![
                    (StereoAtomId(0), StereoAtomId(0)),
                    (StereoAtomId(1), StereoAtomId(1)),
                ],
                2,
                2,
            ),
            GraphCoreCorrespondence::new(vec![(StereoBondId(0), StereoBondId(1))], 1, 2),
        )
    }

    #[rstest]
    #[case::empty(
        GraphCoreCorrespondence::new(vec![], 2, 3),
        Correspondence { matched_pairs: vec![], left_count: 2, right_count: 3 },
    )]
    #[case::partial(
        GraphCoreCorrespondence::new(vec![(NodeId(0), NodeId(2))], 2, 3),
        Correspondence { matched_pairs: vec![(0, 2)], left_count: 2, right_count: 3 },
    )]
    #[case::total(
        GraphCoreCorrespondence::new(
            vec![(NodeId(0), NodeId(1)), (NodeId(1), NodeId(0))],
            2,
            2,
        ),
        Correspondence {
            matched_pairs: vec![(0, 1), (1, 0)],
            left_count: 2,
            right_count: 2,
        },
    )]
    #[case::unsorted(
        GraphCoreCorrespondence::new(
            vec![(NodeId(2), NodeId(0)), (NodeId(0), NodeId(2))],
            3,
            3,
        ),
        Correspondence {
            matched_pairs: vec![(0, 2), (2, 0)],
            left_count: 3,
            right_count: 3,
        },
    )]
    fn test_correspondence_from_rust(
        #[case] correspondence: GraphCoreCorrespondence<NodeId>,
        #[case] expected: Correspondence,
    ) {
        assert_eq!(Correspondence::from_rust(&correspondence), expected);
    }

    #[rstest]
    #[case::empty(
        GraphCoreCorrespondence::new(vec![], 2, 3),
        vec![],
        2,
        3,
        vec![0, 1],
        vec![0, 1, 2],
    )]
    #[case::partial(
        GraphCoreCorrespondence::new(vec![(NodeId(0), NodeId(2))], 2, 3),
        vec![(0, 2)],
        2,
        3,
        vec![1],
        vec![0, 1],
    )]
    #[case::total(
        GraphCoreCorrespondence::new(
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
        #[case] correspondence: GraphCoreCorrespondence<NodeId>,
        #[case] matched_pairs: Vec<(usize, usize)>,
        #[case] left_count: usize,
        #[case] right_count: usize,
        #[case] left_unmatched: Vec<usize>,
        #[case] right_unmatched: Vec<usize>,
    ) {
        let correspondence = Correspondence::from_rust(&correspondence);
        assert_eq!(correspondence.matched_pairs(), matched_pairs);
        assert_eq!(correspondence.left_count(), left_count);
        assert_eq!(correspondence.right_count(), right_count);
        assert_eq!(correspondence.left_unmatched(), left_unmatched);
        assert_eq!(correspondence.right_unmatched(), right_unmatched);
    }

    #[rstest]
    fn test_correspondence_value() {
        let correspondence = Correspondence::from_rust(&GraphCoreCorrespondence::new(
            vec![(NodeId(0), NodeId(2))],
            2,
            3,
        ));
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
            Correspondence::from_rust(&GraphCoreCorrespondence::new(
                vec![(NodeId(0), NodeId(2))],
                2,
                3,
            ))
        );
        assert_ne!(
            correspondence,
            Correspondence::from_rust(&GraphCoreCorrespondence::new(
                vec![(NodeId(0), NodeId(2))],
                2,
                4,
            ))
        );
    }

    #[rstest]
    fn test_molecule_correspondence_accessors(molecule_correspondence: AstMoleculeCorrespondence) {
        let correspondence = MoleculeCorrespondence::from_rust(molecule_correspondence);

        assert_eq!(
            correspondence.atoms(),
            Correspondence {
                matched_pairs: vec![(0, 1)],
                left_count: 2,
                right_count: 3,
            }
        );
        assert_eq!(
            correspondence.bonds(),
            Correspondence {
                matched_pairs: vec![(0, 2)],
                left_count: 1,
                right_count: 3,
            }
        );
        assert_eq!(
            correspondence.dative_bonds(),
            Correspondence {
                matched_pairs: vec![(1, 0)],
                left_count: 2,
                right_count: 1,
            }
        );
        assert_eq!(
            correspondence.aromatic_systems(),
            Correspondence {
                matched_pairs: vec![],
                left_count: 1,
                right_count: 2,
            }
        );
        assert_eq!(
            correspondence.multicenter_bonds(),
            Correspondence {
                matched_pairs: vec![(0, 0)],
                left_count: 1,
                right_count: 1,
            }
        );
        assert_eq!(
            correspondence.noncovalent_bonds(),
            Correspondence {
                matched_pairs: vec![(0, 1)],
                left_count: 2,
                right_count: 2,
            }
        );
        assert_eq!(
            correspondence.stereo_atoms(),
            Correspondence {
                matched_pairs: vec![(0, 0), (1, 1)],
                left_count: 2,
                right_count: 2,
            }
        );
        assert_eq!(
            correspondence.stereo_bonds(),
            Correspondence {
                matched_pairs: vec![(0, 1)],
                left_count: 1,
                right_count: 2,
            }
        );
    }

    #[rstest]
    fn test_molecule_correspondence_value(molecule_correspondence: AstMoleculeCorrespondence) {
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
