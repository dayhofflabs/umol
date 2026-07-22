//! Read-only Python values for reaction correspondences and matching algorithms.

use std::collections::HashSet;

use pyo3::prelude::*;
use umol_ast::ast::{
    AromaticSystemId, BondId, DativeBondId, MoleculeCorrespondence as AstMoleculeCorrespondence,
    MulticenterBondId, NoncovalentBondId, StereoAtomId, StereoBondId,
    SubstructureMatchAlgorithm as AstSubstructureMatchAlgorithm,
};
use umol_graph_core::{
    Correspondence as RustCorrespondence, NodeId,
    SubgraphIsomorphismAlgorithm as RustSubgraphIsomorphismAlgorithm,
};

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

/// Strategy used to match molecule structure and overlays.
#[pyclass(from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SubstructureMatchAlgorithm {
    GraphAndOverlays(),
    Incidence(),
}

#[pymethods]
impl SubstructureMatchAlgorithm {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __repr__(&self) -> &'static str {
        match self {
            Self::GraphAndOverlays() => "SubstructureMatchAlgorithm.GraphAndOverlays()",
            Self::Incidence() => "SubstructureMatchAlgorithm.Incidence()",
        }
    }
}

impl SubstructureMatchAlgorithm {
    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion API for substructure match algorithms"
    )]
    pub(crate) fn from_rust(algorithm: AstSubstructureMatchAlgorithm) -> Self {
        match algorithm {
            AstSubstructureMatchAlgorithm::GraphAndOverlays => Self::GraphAndOverlays(),
            AstSubstructureMatchAlgorithm::Incidence => Self::Incidence(),
        }
    }

    pub(crate) fn to_rust(self) -> AstSubstructureMatchAlgorithm {
        match self {
            Self::GraphAndOverlays() => AstSubstructureMatchAlgorithm::GraphAndOverlays,
            Self::Incidence() => AstSubstructureMatchAlgorithm::Incidence,
        }
    }
}

/// Algorithm used to enumerate subgraph-isomorphism matches for reaction application.
#[pyclass(from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SubgraphIsomorphismAlgorithm {
    Vf2(),
    Ullmann(),
    Ri(),
    ArcMatch { path_length: usize },
    Vf2Rdkit(),
    RayKirsch(),
}

#[pymethods]
impl SubgraphIsomorphismAlgorithm {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __repr__(&self) -> String {
        match self {
            Self::Vf2() => "SubgraphIsomorphismAlgorithm.Vf2()".to_owned(),
            Self::Ullmann() => "SubgraphIsomorphismAlgorithm.Ullmann()".to_owned(),
            Self::Ri() => "SubgraphIsomorphismAlgorithm.Ri()".to_owned(),
            Self::ArcMatch { path_length } => {
                format!("SubgraphIsomorphismAlgorithm.ArcMatch(path_length={path_length})")
            }
            Self::Vf2Rdkit() => "SubgraphIsomorphismAlgorithm.Vf2Rdkit()".to_owned(),
            Self::RayKirsch() => "SubgraphIsomorphismAlgorithm.RayKirsch()".to_owned(),
        }
    }
}

impl SubgraphIsomorphismAlgorithm {
    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion API for subgraph algorithms"
    )]
    pub(crate) fn from_rust(algorithm: RustSubgraphIsomorphismAlgorithm) -> Self {
        match algorithm {
            RustSubgraphIsomorphismAlgorithm::Vf2 => Self::Vf2(),
            RustSubgraphIsomorphismAlgorithm::Ullmann => Self::Ullmann(),
            RustSubgraphIsomorphismAlgorithm::Ri => Self::Ri(),
            RustSubgraphIsomorphismAlgorithm::ArcMatch { path_length } => {
                Self::ArcMatch { path_length }
            }
            RustSubgraphIsomorphismAlgorithm::Vf2Rdkit => Self::Vf2Rdkit(),
            RustSubgraphIsomorphismAlgorithm::RayKirsch => Self::RayKirsch(),
        }
    }

    pub(crate) fn to_rust(self) -> RustSubgraphIsomorphismAlgorithm {
        match self {
            Self::Vf2() => RustSubgraphIsomorphismAlgorithm::Vf2,
            Self::Ullmann() => RustSubgraphIsomorphismAlgorithm::Ullmann,
            Self::Ri() => RustSubgraphIsomorphismAlgorithm::Ri,
            Self::ArcMatch { path_length } => {
                RustSubgraphIsomorphismAlgorithm::ArcMatch { path_length }
            }
            Self::Vf2Rdkit() => RustSubgraphIsomorphismAlgorithm::Vf2Rdkit,
            Self::RayKirsch() => RustSubgraphIsomorphismAlgorithm::RayKirsch,
        }
    }
}

fn exposed_ids(count: usize, mated: impl Iterator<Item = usize>) -> Vec<usize> {
    let present: HashSet<_> = mated.collect();
    (0..count).filter(|id| !present.contains(id)).collect()
}

#[cfg(test)]
mod tests {
    use pyo3::types::PyList;
    use rstest::{fixture, rstest};

    use super::*;
    use crate::convert::into_py_variant;

    #[fixture]
    fn molecule_correspondence() -> AstMoleculeCorrespondence {
        AstMoleculeCorrespondence::new(
            RustCorrespondence::new(vec![(NodeId(0), NodeId(1))], 2, 3),
            RustCorrespondence::new(vec![(BondId(0), BondId(2))], 1, 3),
            RustCorrespondence::new(vec![(DativeBondId(1), DativeBondId(0))], 2, 1),
            RustCorrespondence::new(vec![], 1, 2),
            RustCorrespondence::new(vec![(MulticenterBondId(0), MulticenterBondId(0))], 1, 1),
            RustCorrespondence::new(vec![(NoncovalentBondId(0), NoncovalentBondId(1))], 2, 2),
            RustCorrespondence::new(
                vec![
                    (StereoAtomId(0), StereoAtomId(0)),
                    (StereoAtomId(1), StereoAtomId(1)),
                ],
                2,
                2,
            ),
            RustCorrespondence::new(vec![(StereoBondId(0), StereoBondId(1))], 1, 2),
        )
    }

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

    #[rstest]
    fn test_molecule_correspondence_accessors(molecule_correspondence: AstMoleculeCorrespondence) {
        let correspondence = MoleculeCorrespondence::from_rust(molecule_correspondence);

        assert_eq!(
            correspondence.atoms(),
            Correspondence {
                mates: vec![(0, 1)],
                left_count: 2,
                right_count: 3,
            }
        );
        assert_eq!(
            correspondence.bonds(),
            Correspondence {
                mates: vec![(0, 2)],
                left_count: 1,
                right_count: 3,
            }
        );
        assert_eq!(
            correspondence.dative_bonds(),
            Correspondence {
                mates: vec![(1, 0)],
                left_count: 2,
                right_count: 1,
            }
        );
        assert_eq!(
            correspondence.aromatic_systems(),
            Correspondence {
                mates: vec![],
                left_count: 1,
                right_count: 2,
            }
        );
        assert_eq!(
            correspondence.multicenter_bonds(),
            Correspondence {
                mates: vec![(0, 0)],
                left_count: 1,
                right_count: 1,
            }
        );
        assert_eq!(
            correspondence.noncovalent_bonds(),
            Correspondence {
                mates: vec![(0, 1)],
                left_count: 2,
                right_count: 2,
            }
        );
        assert_eq!(
            correspondence.stereo_atoms(),
            Correspondence {
                mates: vec![(0, 0), (1, 1)],
                left_count: 2,
                right_count: 2,
            }
        );
        assert_eq!(
            correspondence.stereo_bonds(),
            Correspondence {
                mates: vec![(0, 1)],
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
            let first_mates = first
                .getattr("mates")
                .unwrap()
                .cast_into::<PyList>()
                .unwrap();
            first_mates.append((1, 2)).unwrap();

            assert!(!first.is(&second));
            assert_eq!(
                second
                    .getattr("mates")
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
                "atoms=Correspondence(mates=[(0, 1)], left_count=2, right_count=3), ",
                "bonds=Correspondence(mates=[(0, 2)], left_count=1, right_count=3), ",
                "dative_bonds=Correspondence(mates=[(1, 0)], left_count=2, right_count=1), ",
                "aromatic_systems=Correspondence(mates=[], left_count=1, right_count=2), ",
                "multicenter_bonds=Correspondence(mates=[(0, 0)], left_count=1, right_count=1), ",
                "noncovalent_bonds=Correspondence(mates=[(0, 1)], left_count=2, right_count=2), ",
                "stereo_atoms=Correspondence(mates=[(0, 0), (1, 1)], left_count=2, right_count=2), ",
                "stereo_bonds=Correspondence(mates=[(0, 1)], left_count=1, right_count=2))"
            )
        );
    }

    #[rstest]
    #[case::graph_and_overlays(
        AstSubstructureMatchAlgorithm::GraphAndOverlays,
        SubstructureMatchAlgorithm::GraphAndOverlays()
    )]
    #[case::incidence(
        AstSubstructureMatchAlgorithm::Incidence,
        SubstructureMatchAlgorithm::Incidence()
    )]
    fn test_substructure_match_algorithm_from_rust(
        #[case] algorithm: AstSubstructureMatchAlgorithm,
        #[case] expected: SubstructureMatchAlgorithm,
    ) {
        assert_eq!(SubstructureMatchAlgorithm::from_rust(algorithm), expected);
    }

    #[rstest]
    #[case::graph_and_overlays(
        SubstructureMatchAlgorithm::GraphAndOverlays(),
        AstSubstructureMatchAlgorithm::GraphAndOverlays
    )]
    #[case::incidence(
        SubstructureMatchAlgorithm::Incidence(),
        AstSubstructureMatchAlgorithm::Incidence
    )]
    fn test_substructure_match_algorithm_to_rust(
        #[case] algorithm: SubstructureMatchAlgorithm,
        #[case] expected: AstSubstructureMatchAlgorithm,
    ) {
        assert_eq!(algorithm.to_rust(), expected);
    }

    #[rstest]
    #[case::graph_and_overlays(
        SubstructureMatchAlgorithm::GraphAndOverlays(),
        "SubstructureMatchAlgorithm.GraphAndOverlays()"
    )]
    #[case::incidence(
        SubstructureMatchAlgorithm::Incidence(),
        "SubstructureMatchAlgorithm.Incidence()"
    )]
    fn test_substructure_match_algorithm_value(
        #[case] algorithm: SubstructureMatchAlgorithm,
        #[case] expected_repr: &str,
    ) {
        Python::attach(|py| {
            let expected = into_py_variant(
                py,
                SubstructureMatchAlgorithm::from_rust(algorithm.to_rust()),
            )
            .unwrap();
            let algorithm = into_py_variant(py, algorithm).unwrap();
            let expected = expected.bind(py).as_any();
            let algorithm = algorithm.bind(py).as_any();

            assert_eq!(
                algorithm.repr().unwrap().extract::<String>().unwrap(),
                expected_repr
            );
            assert!(algorithm.eq(expected).unwrap());
        });
    }

    #[rstest]
    #[case::vf2(
        RustSubgraphIsomorphismAlgorithm::Vf2,
        SubgraphIsomorphismAlgorithm::Vf2()
    )]
    #[case::ullmann(
        RustSubgraphIsomorphismAlgorithm::Ullmann,
        SubgraphIsomorphismAlgorithm::Ullmann()
    )]
    #[case::ri(
        RustSubgraphIsomorphismAlgorithm::Ri,
        SubgraphIsomorphismAlgorithm::Ri()
    )]
    #[case::arc_match(
        RustSubgraphIsomorphismAlgorithm::ArcMatch { path_length: 6 },
        SubgraphIsomorphismAlgorithm::ArcMatch { path_length: 6 },
    )]
    #[case::vf2_rdkit(
        RustSubgraphIsomorphismAlgorithm::Vf2Rdkit,
        SubgraphIsomorphismAlgorithm::Vf2Rdkit()
    )]
    #[case::ray_kirsch(
        RustSubgraphIsomorphismAlgorithm::RayKirsch,
        SubgraphIsomorphismAlgorithm::RayKirsch()
    )]
    fn test_subgraph_isomorphism_algorithm_from_rust(
        #[case] algorithm: RustSubgraphIsomorphismAlgorithm,
        #[case] expected: SubgraphIsomorphismAlgorithm,
    ) {
        assert_eq!(SubgraphIsomorphismAlgorithm::from_rust(algorithm), expected);
    }

    #[rstest]
    #[case::vf2(
        SubgraphIsomorphismAlgorithm::Vf2(),
        RustSubgraphIsomorphismAlgorithm::Vf2
    )]
    #[case::ullmann(
        SubgraphIsomorphismAlgorithm::Ullmann(),
        RustSubgraphIsomorphismAlgorithm::Ullmann
    )]
    #[case::ri(
        SubgraphIsomorphismAlgorithm::Ri(),
        RustSubgraphIsomorphismAlgorithm::Ri
    )]
    #[case::arc_match(
        SubgraphIsomorphismAlgorithm::ArcMatch { path_length: 6 },
        RustSubgraphIsomorphismAlgorithm::ArcMatch { path_length: 6 },
    )]
    #[case::vf2_rdkit(
        SubgraphIsomorphismAlgorithm::Vf2Rdkit(),
        RustSubgraphIsomorphismAlgorithm::Vf2Rdkit
    )]
    #[case::ray_kirsch(
        SubgraphIsomorphismAlgorithm::RayKirsch(),
        RustSubgraphIsomorphismAlgorithm::RayKirsch
    )]
    fn test_subgraph_isomorphism_algorithm_to_rust(
        #[case] algorithm: SubgraphIsomorphismAlgorithm,
        #[case] expected: RustSubgraphIsomorphismAlgorithm,
    ) {
        assert_eq!(algorithm.to_rust(), expected);
    }

    #[rstest]
    #[case::vf2(
        SubgraphIsomorphismAlgorithm::Vf2(),
        "SubgraphIsomorphismAlgorithm.Vf2()",
        None
    )]
    #[case::ullmann(
        SubgraphIsomorphismAlgorithm::Ullmann(),
        "SubgraphIsomorphismAlgorithm.Ullmann()",
        None
    )]
    #[case::ri(
        SubgraphIsomorphismAlgorithm::Ri(),
        "SubgraphIsomorphismAlgorithm.Ri()",
        None
    )]
    #[case::arc_match(
        SubgraphIsomorphismAlgorithm::ArcMatch { path_length: 6 },
        "SubgraphIsomorphismAlgorithm.ArcMatch(path_length=6)",
        Some(6),
    )]
    #[case::vf2_rdkit(
        SubgraphIsomorphismAlgorithm::Vf2Rdkit(),
        "SubgraphIsomorphismAlgorithm.Vf2Rdkit()",
        None
    )]
    #[case::ray_kirsch(
        SubgraphIsomorphismAlgorithm::RayKirsch(),
        "SubgraphIsomorphismAlgorithm.RayKirsch()",
        None
    )]
    fn test_subgraph_isomorphism_algorithm_value(
        #[case] algorithm: SubgraphIsomorphismAlgorithm,
        #[case] expected_repr: &str,
        #[case] expected_path_length: Option<usize>,
    ) {
        Python::attach(|py| {
            let expected = into_py_variant(
                py,
                SubgraphIsomorphismAlgorithm::from_rust(algorithm.to_rust()),
            )
            .unwrap();
            let algorithm = into_py_variant(py, algorithm).unwrap();
            let expected = expected.bind(py).as_any();
            let algorithm = algorithm.bind(py).as_any();

            assert_eq!(
                algorithm.repr().unwrap().extract::<String>().unwrap(),
                expected_repr
            );
            assert!(algorithm.eq(expected).unwrap());
            match expected_path_length {
                Some(path_length) => assert_eq!(
                    algorithm
                        .getattr("path_length")
                        .unwrap()
                        .extract::<usize>()
                        .unwrap(),
                    path_length
                ),
                None => assert!(algorithm.getattr("path_length").is_err()),
            }
        });
    }
}
