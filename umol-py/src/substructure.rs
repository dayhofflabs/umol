//! Python configuration for molecule substructure search.

use pyo3::prelude::*;
#[cfg(test)]
use umol_ast::ast::SubstructureMatchAlgorithm as AstSubstructureMatchAlgorithm;
use umol_ast::ast::SubstructureMatchConfig as AstSubstructureMatchConfig;
#[cfg(test)]
use umol_graph_core::{
    RelevantCycleEnumerationAlgorithm as GraphCoreRelevantCycleEnumerationAlgorithm,
    SubgraphIsomorphismAlgorithm as GraphCoreSubgraphIsomorphismAlgorithm,
};

use crate::algorithm::{
    RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm, SubstructureMatchAlgorithm,
};

/// Algorithms used to enumerate substructure matches.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SubstructureSearchConfig {
    match_algorithm: SubstructureMatchAlgorithm,
    subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm,
    relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
}

impl Default for SubstructureSearchConfig {
    fn default() -> Self {
        Self {
            match_algorithm: SubstructureMatchAlgorithm::GraphAndOverlays(),
            subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::Vf2Rdkit(),
            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara(),
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
        relevant_cycle_algorithm=RelevantCycleEnumerationAlgorithm::Vismara(),
    ))]
    fn new(
        match_algorithm: SubstructureMatchAlgorithm,
        subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm,
        relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
    ) -> Self {
        Self {
            match_algorithm,
            subgraph_isomorphism_algorithm,
            relevant_cycle_algorithm,
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

    #[getter]
    fn relevant_cycle_algorithm(&self) -> RelevantCycleEnumerationAlgorithm {
        self.relevant_cycle_algorithm
    }

    fn __repr__(&self) -> String {
        format!(
            "SubstructureSearchConfig(match_algorithm={}, subgraph_isomorphism_algorithm={}, relevant_cycle_algorithm={})",
            self.match_algorithm.repr(),
            self.subgraph_isomorphism_algorithm.repr(),
            self.relevant_cycle_algorithm.repr(),
        )
    }
}

impl SubstructureSearchConfig {
    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion API for configured substructure search"
    )]
    pub(crate) fn from_rust(config: AstSubstructureMatchConfig) -> Self {
        Self {
            match_algorithm: SubstructureMatchAlgorithm::from_rust(config.match_algorithm),
            subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::from_rust(
                config.subgraph_isomorphism_algorithm,
            ),
            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::from_rust(
                config.relevant_cycle_algorithm,
            ),
        }
    }

    pub(crate) fn to_rust(self) -> AstSubstructureMatchConfig {
        AstSubstructureMatchConfig {
            match_algorithm: self.match_algorithm.to_rust(),
            subgraph_isomorphism_algorithm: self.subgraph_isomorphism_algorithm.to_rust(),
            relevant_cycle_algorithm: self.relevant_cycle_algorithm.to_rust(),
        }
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
            RelevantCycleEnumerationAlgorithm::Vismara(),
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
            config.relevant_cycle_algorithm(),
            RelevantCycleEnumerationAlgorithm::Vismara()
        );
        assert_eq!(
            config.__repr__(),
            concat!(
                "SubstructureSearchConfig(",
                "match_algorithm=SubstructureMatchAlgorithm.Incidence(), ",
                "subgraph_isomorphism_algorithm=",
                "SubgraphIsomorphismAlgorithm.ArcMatch(path_length=6), ",
                "relevant_cycle_algorithm=RelevantCycleEnumerationAlgorithm.Vismara())"
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
                RelevantCycleEnumerationAlgorithm::Vismara(),
            )
        );
    }

    #[rstest]
    #[case::default(
        AstSubstructureMatchAlgorithm::GraphAndOverlays,
        GraphCoreSubgraphIsomorphismAlgorithm::Vf2Rdkit,
        GraphCoreRelevantCycleEnumerationAlgorithm::Vismara,
        SubstructureSearchConfig::default()
    )]
    #[case::incidence_arc_match(
        AstSubstructureMatchAlgorithm::Incidence,
        GraphCoreSubgraphIsomorphismAlgorithm::ArcMatch { path_length: 6 },
        GraphCoreRelevantCycleEnumerationAlgorithm::Vismara,
        SubstructureSearchConfig::new(
            SubstructureMatchAlgorithm::Incidence(),
            SubgraphIsomorphismAlgorithm::ArcMatch { path_length: 6 },
            RelevantCycleEnumerationAlgorithm::Vismara(),
        ),
    )]
    fn test_substructure_search_config_from_rust(
        #[case] match_algorithm: AstSubstructureMatchAlgorithm,
        #[case] subgraph_isomorphism_algorithm: GraphCoreSubgraphIsomorphismAlgorithm,
        #[case] relevant_cycle_algorithm: GraphCoreRelevantCycleEnumerationAlgorithm,
        #[case] expected: SubstructureSearchConfig,
    ) {
        assert_eq!(
            SubstructureSearchConfig::from_rust(AstSubstructureMatchConfig {
                match_algorithm,
                subgraph_isomorphism_algorithm,
                relevant_cycle_algorithm,
            }),
            expected
        );
    }

    #[rstest]
    #[case::graph_vf2_rdkit(
        SubstructureSearchConfig::default(),
        AstSubstructureMatchConfig {
            match_algorithm: AstSubstructureMatchAlgorithm::GraphAndOverlays,
            subgraph_isomorphism_algorithm: GraphCoreSubgraphIsomorphismAlgorithm::Vf2Rdkit,
            relevant_cycle_algorithm: GraphCoreRelevantCycleEnumerationAlgorithm::Vismara,
        }
    )]
    #[case::incidence_ullmann(
        SubstructureSearchConfig::new(
            SubstructureMatchAlgorithm::Incidence(),
            SubgraphIsomorphismAlgorithm::Ullmann(),
            RelevantCycleEnumerationAlgorithm::Vismara(),
        ),
        AstSubstructureMatchConfig {
            match_algorithm: AstSubstructureMatchAlgorithm::Incidence,
            subgraph_isomorphism_algorithm: GraphCoreSubgraphIsomorphismAlgorithm::Ullmann,
            relevant_cycle_algorithm: GraphCoreRelevantCycleEnumerationAlgorithm::Vismara,
        }
    )]
    fn test_substructure_search_config_to_rust(
        #[case] config: SubstructureSearchConfig,
        #[case] expected: AstSubstructureMatchConfig,
    ) {
        assert_eq!(config.to_rust(), expected);
    }
}
