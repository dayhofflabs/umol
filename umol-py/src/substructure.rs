//! Python configuration for molecule substructure search.

use pyo3::prelude::*;
use umol_ast::ast::SubstructureMatchAlgorithm as AstSubstructureMatchAlgorithm;
use umol_graph_core::SubgraphIsomorphismAlgorithm as GraphCoreSubgraphIsomorphismAlgorithm;

use crate::algorithm::{SubgraphIsomorphismAlgorithm, SubstructureMatchAlgorithm};

/// Algorithms used to enumerate substructure matches.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SubstructureSearchConfig {
    match_algorithm: SubstructureMatchAlgorithm,
    subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm,
}

impl Default for SubstructureSearchConfig {
    fn default() -> Self {
        Self {
            match_algorithm: SubstructureMatchAlgorithm::GraphAndOverlays(),
            subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::Vf2Rdkit(),
        }
    }
}

#[pymethods]
impl SubstructureSearchConfig {
    #[new]
    #[pyo3(signature = (
        *,
        match_algorithm=SubstructureMatchAlgorithm::GraphAndOverlays(),
        subgraph_isomorphism_algorithm=SubgraphIsomorphismAlgorithm::Vf2Rdkit(),
    ))]
    fn new(
        match_algorithm: SubstructureMatchAlgorithm,
        subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm,
    ) -> Self {
        Self {
            match_algorithm,
            subgraph_isomorphism_algorithm,
        }
    }

    #[staticmethod]
    fn default() -> Self {
        Default::default()
    }

    #[getter]
    fn match_algorithm(&self) -> SubstructureMatchAlgorithm {
        self.match_algorithm
    }

    #[getter]
    fn subgraph_isomorphism_algorithm(&self) -> SubgraphIsomorphismAlgorithm {
        self.subgraph_isomorphism_algorithm
    }

    fn __repr__(&self) -> String {
        format!(
            "SubstructureSearchConfig(match_algorithm={}, subgraph_isomorphism_algorithm={})",
            self.match_algorithm.repr(),
            self.subgraph_isomorphism_algorithm.repr(),
        )
    }
}

impl SubstructureSearchConfig {
    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion API for configured substructure search"
    )]
    pub(crate) fn from_rust(
        match_algorithm: AstSubstructureMatchAlgorithm,
        subgraph_isomorphism_algorithm: GraphCoreSubgraphIsomorphismAlgorithm,
    ) -> Self {
        Self {
            match_algorithm: SubstructureMatchAlgorithm::from_rust(match_algorithm),
            subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::from_rust(
                subgraph_isomorphism_algorithm,
            ),
        }
    }

    pub(crate) fn to_rust(
        self,
    ) -> (
        AstSubstructureMatchAlgorithm,
        GraphCoreSubgraphIsomorphismAlgorithm,
    ) {
        (
            self.match_algorithm.to_rust(),
            self.subgraph_isomorphism_algorithm.to_rust(),
        )
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    fn test_substructure_search_config_new() {
        let config = SubstructureSearchConfig::new(
            SubstructureMatchAlgorithm::Incidence(),
            SubgraphIsomorphismAlgorithm::ArcMatch { path_length: 6 },
        );

        assert_eq!(
            config.match_algorithm(),
            SubstructureMatchAlgorithm::Incidence()
        );
        assert_eq!(
            config.subgraph_isomorphism_algorithm(),
            SubgraphIsomorphismAlgorithm::ArcMatch { path_length: 6 }
        );
        assert_eq!(
            config.__repr__(),
            concat!(
                "SubstructureSearchConfig(",
                "match_algorithm=SubstructureMatchAlgorithm.Incidence(), ",
                "subgraph_isomorphism_algorithm=",
                "SubgraphIsomorphismAlgorithm.ArcMatch(path_length=6))"
            )
        );
        assert_ne!(config, SubstructureSearchConfig::default());
    }

    #[rstest]
    fn test_substructure_search_config_default() {
        assert_eq!(
            SubstructureSearchConfig::default(),
            SubstructureSearchConfig::new(
                SubstructureMatchAlgorithm::GraphAndOverlays(),
                SubgraphIsomorphismAlgorithm::Vf2Rdkit(),
            )
        );
    }

    #[rstest]
    #[case::default(
        AstSubstructureMatchAlgorithm::GraphAndOverlays,
        GraphCoreSubgraphIsomorphismAlgorithm::Vf2Rdkit,
        SubstructureSearchConfig::default()
    )]
    #[case::incidence_arc_match(
        AstSubstructureMatchAlgorithm::Incidence,
        GraphCoreSubgraphIsomorphismAlgorithm::ArcMatch { path_length: 6 },
        SubstructureSearchConfig::new(
            SubstructureMatchAlgorithm::Incidence(),
            SubgraphIsomorphismAlgorithm::ArcMatch { path_length: 6 },
        ),
    )]
    fn test_substructure_search_config_from_rust(
        #[case] match_algorithm: AstSubstructureMatchAlgorithm,
        #[case] subgraph_isomorphism_algorithm: GraphCoreSubgraphIsomorphismAlgorithm,
        #[case] expected: SubstructureSearchConfig,
    ) {
        assert_eq!(
            SubstructureSearchConfig::from_rust(match_algorithm, subgraph_isomorphism_algorithm,),
            expected
        );
    }

    #[rstest]
    #[case::graph_vf2_rdkit(
        SubstructureSearchConfig::default(),
        (
            AstSubstructureMatchAlgorithm::GraphAndOverlays,
            GraphCoreSubgraphIsomorphismAlgorithm::Vf2Rdkit,
        )
    )]
    #[case::incidence_ullmann(
        SubstructureSearchConfig::new(
            SubstructureMatchAlgorithm::Incidence(),
            SubgraphIsomorphismAlgorithm::Ullmann(),
        ),
        (
            AstSubstructureMatchAlgorithm::Incidence,
            GraphCoreSubgraphIsomorphismAlgorithm::Ullmann,
        )
    )]
    fn test_substructure_search_config_to_rust(
        #[case] config: SubstructureSearchConfig,
        #[case] expected: (
            AstSubstructureMatchAlgorithm,
            GraphCoreSubgraphIsomorphismAlgorithm,
        ),
    ) {
        assert_eq!(config.to_rust(), expected);
    }
}
