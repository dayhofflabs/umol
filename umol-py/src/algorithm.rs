//! Python values selecting graph algorithms.

use pyo3::prelude::*;
use umol_ast::ast::SubstructureMatchAlgorithm as AstSubstructureMatchAlgorithm;
use umol_graph_core::{
    AutomorphismAlgorithm as GraphCoreAutomorphismAlgorithm,
    CommonSubgraphEnumerationAlgorithm as GraphCoreCommonSubgraphEnumerationAlgorithm,
    ConnectedComponentsAlgorithm as GraphCoreConnectedComponentsAlgorithm,
    CycleEnumerationAlgorithm as GraphCoreCycleEnumerationAlgorithm,
    MaximumIndependentSetAlgorithm as GraphCoreMaximumIndependentSetAlgorithm,
    SubgraphEnumerationAlgorithm as GraphCoreSubgraphEnumerationAlgorithm,
    SubgraphIsomorphismAlgorithm as GraphCoreSubgraphIsomorphismAlgorithm,
};

/// Algorithm used to find graph automorphisms.
#[pyclass(from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AutomorphismAlgorithm {
    Nauty(),
}

#[pymethods]
impl AutomorphismAlgorithm {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __repr__(&self) -> &'static str {
        self.repr()
    }
}

impl AutomorphismAlgorithm {
    pub(crate) fn repr(self) -> &'static str {
        match self {
            Self::Nauty() => "AutomorphismAlgorithm.Nauty()",
        }
    }

    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion API for automorphism algorithms"
    )]
    pub(crate) fn from_rust(algorithm: GraphCoreAutomorphismAlgorithm) -> Self {
        match algorithm {
            GraphCoreAutomorphismAlgorithm::Nauty => Self::Nauty(),
        }
    }

    pub(crate) fn to_rust(self) -> GraphCoreAutomorphismAlgorithm {
        match self {
            Self::Nauty() => GraphCoreAutomorphismAlgorithm::Nauty,
        }
    }
}

/// Algorithm used to enumerate every common subgraph.
#[pyclass(from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CommonSubgraphEnumerationAlgorithm {
    Backtracking(),
}

#[pymethods]
impl CommonSubgraphEnumerationAlgorithm {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __repr__(&self) -> &'static str {
        self.repr()
    }
}

impl CommonSubgraphEnumerationAlgorithm {
    pub(crate) fn repr(self) -> &'static str {
        match self {
            Self::Backtracking() => "CommonSubgraphEnumerationAlgorithm.Backtracking()",
        }
    }

    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion API for common-subgraph enumeration algorithms"
    )]
    pub(crate) fn from_rust(algorithm: GraphCoreCommonSubgraphEnumerationAlgorithm) -> Self {
        match algorithm {
            GraphCoreCommonSubgraphEnumerationAlgorithm::Backtracking => Self::Backtracking(),
        }
    }

    pub(crate) fn to_rust(self) -> GraphCoreCommonSubgraphEnumerationAlgorithm {
        match self {
            Self::Backtracking() => GraphCoreCommonSubgraphEnumerationAlgorithm::Backtracking,
        }
    }
}

/// Algorithm used to label connected components.
#[pyclass(from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConnectedComponentsAlgorithm {
    Bfs(),
}

#[pymethods]
impl ConnectedComponentsAlgorithm {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __repr__(&self) -> &'static str {
        self.repr()
    }
}

impl ConnectedComponentsAlgorithm {
    pub(crate) fn repr(self) -> &'static str {
        match self {
            Self::Bfs() => "ConnectedComponentsAlgorithm.Bfs()",
        }
    }

    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion API for connected-components algorithms"
    )]
    pub(crate) fn from_rust(algorithm: GraphCoreConnectedComponentsAlgorithm) -> Self {
        match algorithm {
            GraphCoreConnectedComponentsAlgorithm::Bfs => Self::Bfs(),
        }
    }

    pub(crate) fn to_rust(self) -> GraphCoreConnectedComponentsAlgorithm {
        match self {
            Self::Bfs() => GraphCoreConnectedComponentsAlgorithm::Bfs,
        }
    }
}

/// Algorithm used to enumerate relevant cycles.
#[pyclass(from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CycleEnumerationAlgorithm {
    Vismara(),
}

#[pymethods]
impl CycleEnumerationAlgorithm {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __repr__(&self) -> &'static str {
        self.repr()
    }
}

impl CycleEnumerationAlgorithm {
    pub(crate) fn repr(self) -> &'static str {
        match self {
            Self::Vismara() => "CycleEnumerationAlgorithm.Vismara()",
        }
    }

    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion API for cycle-enumeration algorithms"
    )]
    pub(crate) fn from_rust(algorithm: GraphCoreCycleEnumerationAlgorithm) -> Self {
        match algorithm {
            GraphCoreCycleEnumerationAlgorithm::Vismara => Self::Vismara(),
        }
    }

    pub(crate) fn to_rust(self) -> GraphCoreCycleEnumerationAlgorithm {
        match self {
            Self::Vismara() => GraphCoreCycleEnumerationAlgorithm::Vismara,
        }
    }
}

/// Algorithm used to find a maximum independent set.
#[pyclass(from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MaximumIndependentSetAlgorithm {
    BranchAndBound(),
}

#[pymethods]
impl MaximumIndependentSetAlgorithm {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __repr__(&self) -> &'static str {
        self.repr()
    }
}

impl MaximumIndependentSetAlgorithm {
    pub(crate) fn repr(self) -> &'static str {
        match self {
            Self::BranchAndBound() => "MaximumIndependentSetAlgorithm.BranchAndBound()",
        }
    }

    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion API for maximum-independent-set algorithms"
    )]
    pub(crate) fn from_rust(algorithm: GraphCoreMaximumIndependentSetAlgorithm) -> Self {
        match algorithm {
            GraphCoreMaximumIndependentSetAlgorithm::BranchAndBound => Self::BranchAndBound(),
        }
    }

    pub(crate) fn to_rust(self) -> GraphCoreMaximumIndependentSetAlgorithm {
        match self {
            Self::BranchAndBound() => GraphCoreMaximumIndependentSetAlgorithm::BranchAndBound,
        }
    }
}

/// Algorithm used to enumerate connected subgraphs.
#[pyclass(from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SubgraphEnumerationAlgorithm {
    Esu(),
}

#[pymethods]
impl SubgraphEnumerationAlgorithm {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __repr__(&self) -> &'static str {
        self.repr()
    }
}

impl SubgraphEnumerationAlgorithm {
    pub(crate) fn repr(self) -> &'static str {
        match self {
            Self::Esu() => "SubgraphEnumerationAlgorithm.Esu()",
        }
    }

    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion API for subgraph-enumeration algorithms"
    )]
    pub(crate) fn from_rust(algorithm: GraphCoreSubgraphEnumerationAlgorithm) -> Self {
        match algorithm {
            GraphCoreSubgraphEnumerationAlgorithm::Esu => Self::Esu(),
        }
    }

    pub(crate) fn to_rust(self) -> GraphCoreSubgraphEnumerationAlgorithm {
        match self {
            Self::Esu() => GraphCoreSubgraphEnumerationAlgorithm::Esu,
        }
    }
}

/// Algorithm used to enumerate subgraph-isomorphism matches.
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
        self.repr()
    }
}

impl SubgraphIsomorphismAlgorithm {
    pub(crate) fn repr(self) -> String {
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

    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion API for subgraph algorithms"
    )]
    pub(crate) fn from_rust(algorithm: GraphCoreSubgraphIsomorphismAlgorithm) -> Self {
        match algorithm {
            GraphCoreSubgraphIsomorphismAlgorithm::Vf2 => Self::Vf2(),
            GraphCoreSubgraphIsomorphismAlgorithm::Ullmann => Self::Ullmann(),
            GraphCoreSubgraphIsomorphismAlgorithm::Ri => Self::Ri(),
            GraphCoreSubgraphIsomorphismAlgorithm::ArcMatch { path_length } => {
                Self::ArcMatch { path_length }
            }
            GraphCoreSubgraphIsomorphismAlgorithm::Vf2Rdkit => Self::Vf2Rdkit(),
            GraphCoreSubgraphIsomorphismAlgorithm::RayKirsch => Self::RayKirsch(),
        }
    }

    pub(crate) fn to_rust(self) -> GraphCoreSubgraphIsomorphismAlgorithm {
        match self {
            Self::Vf2() => GraphCoreSubgraphIsomorphismAlgorithm::Vf2,
            Self::Ullmann() => GraphCoreSubgraphIsomorphismAlgorithm::Ullmann,
            Self::Ri() => GraphCoreSubgraphIsomorphismAlgorithm::Ri,
            Self::ArcMatch { path_length } => {
                GraphCoreSubgraphIsomorphismAlgorithm::ArcMatch { path_length }
            }
            Self::Vf2Rdkit() => GraphCoreSubgraphIsomorphismAlgorithm::Vf2Rdkit,
            Self::RayKirsch() => GraphCoreSubgraphIsomorphismAlgorithm::RayKirsch,
        }
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
        self.repr()
    }
}

impl SubstructureMatchAlgorithm {
    pub(crate) fn repr(self) -> &'static str {
        match self {
            Self::GraphAndOverlays() => "SubstructureMatchAlgorithm.GraphAndOverlays()",
            Self::Incidence() => "SubstructureMatchAlgorithm.Incidence()",
        }
    }

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

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;
    use crate::convert::into_py_variant;

    #[rstest]
    #[case::nauty(
        GraphCoreAutomorphismAlgorithm::Nauty,
        AutomorphismAlgorithm::Nauty(),
        "AutomorphismAlgorithm.Nauty()"
    )]
    fn test_automorphism_algorithm_value(
        #[case] rust: GraphCoreAutomorphismAlgorithm,
        #[case] python: AutomorphismAlgorithm,
        #[case] expected_repr: &str,
    ) {
        assert_eq!(AutomorphismAlgorithm::from_rust(rust), python);
        assert_eq!(python.to_rust(), rust);
        Python::attach(|py| {
            let expected = into_py_variant(py, AutomorphismAlgorithm::from_rust(rust)).unwrap();
            let value = into_py_variant(py, python).unwrap();

            assert_eq!(
                value
                    .bind(py)
                    .as_any()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                expected_repr
            );
            assert!(value.bind(py).as_any().eq(expected.bind(py)).unwrap());
        });
    }

    #[rstest]
    #[case::backtracking(
        GraphCoreCommonSubgraphEnumerationAlgorithm::Backtracking,
        CommonSubgraphEnumerationAlgorithm::Backtracking(),
        "CommonSubgraphEnumerationAlgorithm.Backtracking()"
    )]
    fn test_common_subgraph_enumeration_algorithm_value(
        #[case] rust: GraphCoreCommonSubgraphEnumerationAlgorithm,
        #[case] python: CommonSubgraphEnumerationAlgorithm,
        #[case] expected_repr: &str,
    ) {
        assert_eq!(CommonSubgraphEnumerationAlgorithm::from_rust(rust), python);
        assert_eq!(python.to_rust(), rust);
        Python::attach(|py| {
            let expected =
                into_py_variant(py, CommonSubgraphEnumerationAlgorithm::from_rust(rust)).unwrap();
            let value = into_py_variant(py, python).unwrap();

            assert_eq!(
                value
                    .bind(py)
                    .as_any()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                expected_repr
            );
            assert!(value.bind(py).as_any().eq(expected.bind(py)).unwrap());
        });
    }

    #[rstest]
    #[case::bfs(
        GraphCoreConnectedComponentsAlgorithm::Bfs,
        ConnectedComponentsAlgorithm::Bfs(),
        "ConnectedComponentsAlgorithm.Bfs()"
    )]
    fn test_connected_components_algorithm_value(
        #[case] rust: GraphCoreConnectedComponentsAlgorithm,
        #[case] python: ConnectedComponentsAlgorithm,
        #[case] expected_repr: &str,
    ) {
        assert_eq!(ConnectedComponentsAlgorithm::from_rust(rust), python);
        assert_eq!(python.to_rust(), rust);
        Python::attach(|py| {
            let expected =
                into_py_variant(py, ConnectedComponentsAlgorithm::from_rust(rust)).unwrap();
            let value = into_py_variant(py, python).unwrap();

            assert_eq!(
                value
                    .bind(py)
                    .as_any()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                expected_repr
            );
            assert!(value.bind(py).as_any().eq(expected.bind(py)).unwrap());
        });
    }

    #[rstest]
    #[case::vismara(
        GraphCoreCycleEnumerationAlgorithm::Vismara,
        CycleEnumerationAlgorithm::Vismara(),
        "CycleEnumerationAlgorithm.Vismara()"
    )]
    fn test_cycle_enumeration_algorithm_value(
        #[case] rust: GraphCoreCycleEnumerationAlgorithm,
        #[case] python: CycleEnumerationAlgorithm,
        #[case] expected_repr: &str,
    ) {
        assert_eq!(CycleEnumerationAlgorithm::from_rust(rust), python);
        assert_eq!(python.to_rust(), rust);
        Python::attach(|py| {
            let expected = into_py_variant(py, CycleEnumerationAlgorithm::from_rust(rust)).unwrap();
            let value = into_py_variant(py, python).unwrap();

            assert_eq!(
                value
                    .bind(py)
                    .as_any()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                expected_repr
            );
            assert!(value.bind(py).as_any().eq(expected.bind(py)).unwrap());
        });
    }

    #[rstest]
    #[case::branch_and_bound(
        GraphCoreMaximumIndependentSetAlgorithm::BranchAndBound,
        MaximumIndependentSetAlgorithm::BranchAndBound(),
        "MaximumIndependentSetAlgorithm.BranchAndBound()"
    )]
    fn test_maximum_independent_set_algorithm_value(
        #[case] rust: GraphCoreMaximumIndependentSetAlgorithm,
        #[case] python: MaximumIndependentSetAlgorithm,
        #[case] expected_repr: &str,
    ) {
        assert_eq!(MaximumIndependentSetAlgorithm::from_rust(rust), python);
        assert_eq!(python.to_rust(), rust);
        Python::attach(|py| {
            let expected =
                into_py_variant(py, MaximumIndependentSetAlgorithm::from_rust(rust)).unwrap();
            let value = into_py_variant(py, python).unwrap();

            assert_eq!(
                value
                    .bind(py)
                    .as_any()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                expected_repr
            );
            assert!(value.bind(py).as_any().eq(expected.bind(py)).unwrap());
        });
    }

    #[rstest]
    #[case::esu(
        GraphCoreSubgraphEnumerationAlgorithm::Esu,
        SubgraphEnumerationAlgorithm::Esu(),
        "SubgraphEnumerationAlgorithm.Esu()"
    )]
    fn test_subgraph_enumeration_algorithm_value(
        #[case] rust: GraphCoreSubgraphEnumerationAlgorithm,
        #[case] python: SubgraphEnumerationAlgorithm,
        #[case] expected_repr: &str,
    ) {
        assert_eq!(SubgraphEnumerationAlgorithm::from_rust(rust), python);
        assert_eq!(python.to_rust(), rust);
        Python::attach(|py| {
            let expected =
                into_py_variant(py, SubgraphEnumerationAlgorithm::from_rust(rust)).unwrap();
            let value = into_py_variant(py, python).unwrap();

            assert_eq!(
                value
                    .bind(py)
                    .as_any()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                expected_repr
            );
            assert!(value.bind(py).as_any().eq(expected.bind(py)).unwrap());
        });
    }

    #[rstest]
    #[case::vf2(
        GraphCoreSubgraphIsomorphismAlgorithm::Vf2,
        SubgraphIsomorphismAlgorithm::Vf2()
    )]
    #[case::ullmann(
        GraphCoreSubgraphIsomorphismAlgorithm::Ullmann,
        SubgraphIsomorphismAlgorithm::Ullmann()
    )]
    #[case::ri(
        GraphCoreSubgraphIsomorphismAlgorithm::Ri,
        SubgraphIsomorphismAlgorithm::Ri()
    )]
    #[case::arc_match(
        GraphCoreSubgraphIsomorphismAlgorithm::ArcMatch { path_length: 6 },
        SubgraphIsomorphismAlgorithm::ArcMatch { path_length: 6 },
    )]
    #[case::vf2_rdkit(
        GraphCoreSubgraphIsomorphismAlgorithm::Vf2Rdkit,
        SubgraphIsomorphismAlgorithm::Vf2Rdkit()
    )]
    #[case::ray_kirsch(
        GraphCoreSubgraphIsomorphismAlgorithm::RayKirsch,
        SubgraphIsomorphismAlgorithm::RayKirsch()
    )]
    fn test_subgraph_isomorphism_algorithm_from_rust(
        #[case] algorithm: GraphCoreSubgraphIsomorphismAlgorithm,
        #[case] expected: SubgraphIsomorphismAlgorithm,
    ) {
        assert_eq!(SubgraphIsomorphismAlgorithm::from_rust(algorithm), expected);
    }

    #[rstest]
    #[case::vf2(
        SubgraphIsomorphismAlgorithm::Vf2(),
        GraphCoreSubgraphIsomorphismAlgorithm::Vf2
    )]
    #[case::ullmann(
        SubgraphIsomorphismAlgorithm::Ullmann(),
        GraphCoreSubgraphIsomorphismAlgorithm::Ullmann
    )]
    #[case::ri(
        SubgraphIsomorphismAlgorithm::Ri(),
        GraphCoreSubgraphIsomorphismAlgorithm::Ri
    )]
    #[case::arc_match(
        SubgraphIsomorphismAlgorithm::ArcMatch { path_length: 6 },
        GraphCoreSubgraphIsomorphismAlgorithm::ArcMatch { path_length: 6 },
    )]
    #[case::vf2_rdkit(
        SubgraphIsomorphismAlgorithm::Vf2Rdkit(),
        GraphCoreSubgraphIsomorphismAlgorithm::Vf2Rdkit
    )]
    #[case::ray_kirsch(
        SubgraphIsomorphismAlgorithm::RayKirsch(),
        GraphCoreSubgraphIsomorphismAlgorithm::RayKirsch
    )]
    fn test_subgraph_isomorphism_algorithm_to_rust(
        #[case] algorithm: SubgraphIsomorphismAlgorithm,
        #[case] expected: GraphCoreSubgraphIsomorphismAlgorithm,
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
}
