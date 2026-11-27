//! Python configuration for molecule substructure search.

use pyo3::prelude::*;
use umol_graph_core::AutomorphismAlgorithm as GraphCoreAutomorphismAlgorithm;
#[cfg(test)]
use umol_graph_core::{
    RelevantCycleEnumerationAlgorithm as GraphCoreRelevantCycleEnumerationAlgorithm,
    SubgraphIsomorphismAlgorithm as GraphCoreSubgraphIsomorphismAlgorithm,
};
#[cfg(test)]
use umol_graph_ir::ir::SubstructureMatchAlgorithm as GraphIrSubstructureMatchAlgorithm;
use umol_graph_ir::ir::SubstructureMatchConfig as GraphIrSubstructureMatchConfig;

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
    pub(crate) fn from_rust(config: GraphIrSubstructureMatchConfig) -> Self {
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

    pub(crate) fn to_rust(self) -> GraphIrSubstructureMatchConfig {
        GraphIrSubstructureMatchConfig {
            match_algorithm: self.match_algorithm.to_rust(),
            subgraph_isomorphism_algorithm: self.subgraph_isomorphism_algorithm.to_rust(),
            relevant_cycle_algorithm: self.relevant_cycle_algorithm.to_rust(),
            automorphism_algorithm: GraphCoreAutomorphismAlgorithm::TrailIr,
            reduce_pattern_symmetry: false,
            reduce_host_symmetry: false,
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
            SubgraphIsomorphismAlgorithm::BitsetRi(),
            RelevantCycleEnumerationAlgorithm::Vismara(),
        );

        assert_eq!(
            config.match_algorithm(),
            SubstructureMatchAlgorithm::Incidence()
        );
        assert_eq!(
            config.subgraph_isomorphism_algorithm(),
            SubgraphIsomorphismAlgorithm::BitsetRi()
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
                "SubgraphIsomorphismAlgorithm.BitsetRi(), ",
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
        GraphIrSubstructureMatchAlgorithm::GraphAndOverlays,
        GraphCoreSubgraphIsomorphismAlgorithm::Vf2Rdkit,
        GraphCoreRelevantCycleEnumerationAlgorithm::Vismara,
        SubstructureSearchConfig::default()
    )]
    #[case::incidence_bitset_ri(
        GraphIrSubstructureMatchAlgorithm::Incidence,
        GraphCoreSubgraphIsomorphismAlgorithm::BitsetRi,
        GraphCoreRelevantCycleEnumerationAlgorithm::Vismara,
        SubstructureSearchConfig::new(
            SubstructureMatchAlgorithm::Incidence(),
            SubgraphIsomorphismAlgorithm::BitsetRi(),
            RelevantCycleEnumerationAlgorithm::Vismara(),
        )
    )]
    fn test_substructure_search_config_from_rust(
        #[case] match_algorithm: GraphIrSubstructureMatchAlgorithm,
        #[case] subgraph_isomorphism_algorithm: GraphCoreSubgraphIsomorphismAlgorithm,
        #[case] relevant_cycle_algorithm: GraphCoreRelevantCycleEnumerationAlgorithm,
        #[case] expected: SubstructureSearchConfig,
    ) {
        assert_eq!(
            SubstructureSearchConfig::from_rust(GraphIrSubstructureMatchConfig {
                match_algorithm,
                subgraph_isomorphism_algorithm,
                relevant_cycle_algorithm,
                automorphism_algorithm: GraphCoreAutomorphismAlgorithm::TrailIr,
                reduce_pattern_symmetry: false,
                reduce_host_symmetry: false,
            }),
            expected
        );
    }

    #[rstest]
    #[case::graph_vf2_rdkit(
        SubstructureSearchConfig::default(),
        GraphIrSubstructureMatchConfig {
            match_algorithm: GraphIrSubstructureMatchAlgorithm::GraphAndOverlays,
            subgraph_isomorphism_algorithm: GraphCoreSubgraphIsomorphismAlgorithm::Vf2Rdkit,
            relevant_cycle_algorithm: GraphCoreRelevantCycleEnumerationAlgorithm::Vismara,
            automorphism_algorithm: GraphCoreAutomorphismAlgorithm::TrailIr,
            reduce_pattern_symmetry: false,
            reduce_host_symmetry: false,
        }
    )]
    #[case::incidence_glasgow(
        SubstructureSearchConfig::new(
            SubstructureMatchAlgorithm::Incidence(),
            SubgraphIsomorphismAlgorithm::Glasgow(),
            RelevantCycleEnumerationAlgorithm::Vismara(),
        ),
        GraphIrSubstructureMatchConfig {
            match_algorithm: GraphIrSubstructureMatchAlgorithm::Incidence,
            subgraph_isomorphism_algorithm: GraphCoreSubgraphIsomorphismAlgorithm::Glasgow,
            relevant_cycle_algorithm: GraphCoreRelevantCycleEnumerationAlgorithm::Vismara,
            automorphism_algorithm: GraphCoreAutomorphismAlgorithm::TrailIr,
            reduce_pattern_symmetry: false,
            reduce_host_symmetry: false,
        }
    )]
    fn test_substructure_search_config_to_rust(
        #[case] config: SubstructureSearchConfig,
        #[case] expected: GraphIrSubstructureMatchConfig,
    ) {
        assert_eq!(config.to_rust(), expected);
    }
}
