//! Python configuration and prepared carriers for molecule substructure search.

use pyo3::prelude::*;
#[cfg(test)]
use umol_graph_core::{
    AutomorphismAlgorithm as GraphCoreAutomorphismAlgorithm,
    RelevantCycleEnumerationAlgorithm as GraphCoreRelevantCycleEnumerationAlgorithm,
    SubgraphIsomorphismAlgorithm as GraphCoreSubgraphIsomorphismAlgorithm,
};
#[cfg(test)]
use umol_graph_ir::ir::SubstructureMatchAlgorithm as GraphIrSubstructureMatchAlgorithm;
use umol_graph_ir::ir::{
    PreparedSubstructureHost as GraphIrPreparedSubstructureHost,
    PreparedSubstructurePattern as GraphIrPreparedSubstructurePattern, React as GraphIrReact,
    SubstructureMatchConfig as GraphIrSubstructureMatchConfig,
};

use crate::algorithm::{
    AutomorphismAlgorithm, RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm,
    SubstructureMatchAlgorithm,
};
use crate::correspondence::MoleculeCorrespondence;
use crate::error::InvalidStructureError;
use crate::molecule::Molecule;
use crate::reaction::{
    PreparedReaction, Reaction, ReactionApplicationConfig, ReactionProductsIter,
};

/// Algorithms and symmetry-reduction options used to enumerate substructure matches.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SubstructureSearchConfig {
    match_algorithm: SubstructureMatchAlgorithm,
    subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm,
    relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
    automorphism_algorithm: AutomorphismAlgorithm,
    reduce_pattern_symmetry: bool,
    reduce_host_symmetry: bool,
}

impl Default for SubstructureSearchConfig {
    fn default() -> Self {
        Self {
            match_algorithm: SubstructureMatchAlgorithm::GraphAndOverlays(),
            subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::BitsetRi(),
            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara(),
            automorphism_algorithm: AutomorphismAlgorithm::TrailIr(),
            reduce_pattern_symmetry: false,
            reduce_host_symmetry: false,
        }
    }
}

#[pymethods]
impl SubstructureSearchConfig {
    #[new]
    #[pyo3(signature = (
        *,
        match_algorithm=SubstructureMatchAlgorithm::GraphAndOverlays(),
        subgraph_isomorphism_algorithm=SubgraphIsomorphismAlgorithm::BitsetRi(),
        relevant_cycle_algorithm=RelevantCycleEnumerationAlgorithm::Vismara(),
        automorphism_algorithm=AutomorphismAlgorithm::TrailIr(),
        reduce_pattern_symmetry=false,
        reduce_host_symmetry=false,
    ))]
    fn new(
        match_algorithm: SubstructureMatchAlgorithm,
        subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm,
        relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
        automorphism_algorithm: AutomorphismAlgorithm,
        reduce_pattern_symmetry: bool,
        reduce_host_symmetry: bool,
    ) -> Self {
        Self {
            match_algorithm,
            subgraph_isomorphism_algorithm,
            relevant_cycle_algorithm,
            automorphism_algorithm,
            reduce_pattern_symmetry,
            reduce_host_symmetry,
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

    #[getter]
    fn automorphism_algorithm(&self) -> AutomorphismAlgorithm {
        self.automorphism_algorithm
    }

    #[getter]
    fn reduce_pattern_symmetry(&self) -> bool {
        self.reduce_pattern_symmetry
    }

    #[getter]
    fn reduce_host_symmetry(&self) -> bool {
        self.reduce_host_symmetry
    }

    fn __repr__(&self) -> String {
        format!(
            "SubstructureSearchConfig(match_algorithm={}, subgraph_isomorphism_algorithm={}, relevant_cycle_algorithm={}, automorphism_algorithm={}, reduce_pattern_symmetry={}, reduce_host_symmetry={})",
            self.match_algorithm.repr(),
            self.subgraph_isomorphism_algorithm.repr(),
            self.relevant_cycle_algorithm.repr(),
            self.automorphism_algorithm.repr(),
            if self.reduce_pattern_symmetry {
                "True"
            } else {
                "False"
            },
            if self.reduce_host_symmetry {
                "True"
            } else {
                "False"
            },
        )
    }
}

impl SubstructureSearchConfig {
    pub(crate) fn from_rust(config: GraphIrSubstructureMatchConfig) -> Self {
        Self {
            match_algorithm: SubstructureMatchAlgorithm::from_rust(config.match_algorithm),
            subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::from_rust(
                config.subgraph_isomorphism_algorithm,
            ),
            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::from_rust(
                config.relevant_cycle_algorithm,
            ),
            automorphism_algorithm: AutomorphismAlgorithm::from_rust(config.automorphism_algorithm),
            reduce_pattern_symmetry: config.reduce_pattern_symmetry,
            reduce_host_symmetry: config.reduce_host_symmetry,
        }
    }

    pub(crate) fn to_rust(self) -> GraphIrSubstructureMatchConfig {
        GraphIrSubstructureMatchConfig {
            match_algorithm: self.match_algorithm.to_rust(),
            subgraph_isomorphism_algorithm: self.subgraph_isomorphism_algorithm.to_rust(),
            relevant_cycle_algorithm: self.relevant_cycle_algorithm.to_rust(),
            automorphism_algorithm: self.automorphism_algorithm.to_rust(),
            reduce_pattern_symmetry: self.reduce_pattern_symmetry,
            reduce_host_symmetry: self.reduce_host_symmetry,
        }
    }
}

/// Argument adapter accepting either a `Molecule` or a `PreparedSubstructureHost`.
#[derive(FromPyObject)]
pub(crate) enum SubstructureHostLike<'py> {
    Prepared(PyRef<'py, PreparedSubstructureHost>),
    Molecule(PyRef<'py, Molecule>),
}

/// Immutable prepared host molecule for repeated substructure matching and reaction application.
///
/// Captures a detached snapshot of the supplied `Molecule` together with its prebuilt bitset
/// target and lazy ring, incidence, and host-symmetry caches. Mutating the source molecule after
/// construction does not affect this prepared host.
#[pyclass(frozen)]
#[derive(Debug)]
pub struct PreparedSubstructureHost(GraphIrPreparedSubstructureHost);

#[pymethods]
impl PreparedSubstructureHost {
    #[new]
    fn new(molecule: &Molecule) -> Self {
        Self::from_rust(GraphIrPreparedSubstructureHost::new(molecule.to_rust()))
    }

    /// Return a detached snapshot of the prepared host molecule.
    fn molecule(&self) -> Molecule {
        Molecule::from_rust(self.to_rust().molecule().clone())
    }

    /// Apply `reaction` to this prepared host and lazily emit connected product components.
    ///
    /// Reuses this host's prepared bitsets and lazy caches. Reaction-wide precondition failures
    /// raise `InvalidStructureError` here; per-match failures are raised during iteration.
    #[pyo3(signature = (reaction, *, config=None))]
    fn react(
        &self,
        py: Python<'_>,
        reaction: &Reaction,
        config: Option<ReactionApplicationConfig>,
    ) -> PyResult<Py<ReactionProductsIter>> {
        let reaction = reaction.to_rust(py)?;
        let products = GraphIrReact::react(
            self.to_rust(),
            &reaction,
            config.unwrap_or_default().to_rust(),
        )
        .map_err(|error| InvalidStructureError::new_err(error.to_string()))?;

        Py::new(py, ReactionProductsIter::from_rust(products))
    }

    /// Apply a `PreparedReaction` to this prepared host and lazily emit connected product
    /// components.
    ///
    /// Reuses both the prepared reaction and this host's prepared bitsets and lazy caches.
    /// Per-match execution failures are raised during iteration.
    fn react_prepared(
        &self,
        py: Python<'_>,
        reaction: &PreparedReaction,
    ) -> PyResult<Py<ReactionProductsIter>> {
        let products = GraphIrReact::react_prepared(self.to_rust(), reaction.to_rust());
        Py::new(py, ReactionProductsIter::from_rust(products))
    }

    fn __repr__(&self) -> String {
        let molecule = self.to_rust().molecule();
        format!(
            "PreparedSubstructureHost(atoms={}, bonds={})",
            molecule.atoms().count(),
            molecule.bonds().count(),
        )
    }
}

impl PreparedSubstructureHost {
    pub(crate) fn from_rust(host: GraphIrPreparedSubstructureHost) -> Self {
        Self(host)
    }

    pub(crate) fn to_rust(&self) -> &GraphIrPreparedSubstructureHost {
        &self.0
    }
}

/// Immutable prepared pattern molecule for repeated substructure matching.
///
/// Captures a detached snapshot of the supplied pattern `Molecule`, binds its search
/// configuration, and precompiles its query plan, predicate classes, screening masks, and
/// pattern-symmetry constraints. Mutating the source molecule after construction does not affect
/// this prepared pattern.
#[pyclass(frozen)]
#[derive(Debug)]
pub struct PreparedSubstructurePattern(GraphIrPreparedSubstructurePattern);

#[pymethods]
impl PreparedSubstructurePattern {
    #[new]
    #[pyo3(signature = (pattern, *, config=None))]
    fn new(pattern: &Molecule, config: Option<SubstructureSearchConfig>) -> Self {
        let config = config.unwrap_or_default().to_rust();
        Self::from_rust(GraphIrPreparedSubstructurePattern::new(
            pattern.to_rust(),
            config,
        ))
    }

    /// Return a detached snapshot of the prepared pattern molecule.
    fn pattern(&self) -> Molecule {
        Molecule::from_rust(self.to_rust().pattern().clone())
    }

    /// The immutable search configuration bound to this prepared pattern.
    #[getter]
    fn config(&self) -> SubstructureSearchConfig {
        SubstructureSearchConfig::from_rust(self.to_rust().config())
    }

    /// Collect all occurrences of this prepared pattern in a `Molecule` or
    /// `PreparedSubstructureHost`.
    fn matches(&self, host: SubstructureHostLike<'_>) -> Vec<MoleculeCorrespondence> {
        let matches = match host {
            SubstructureHostLike::Prepared(host) => self.to_rust().matches_in_host(host.to_rust()),
            SubstructureHostLike::Molecule(host) => self.to_rust().matches(host.to_rust()),
        };
        matches
            .into_iter()
            .map(MoleculeCorrespondence::from_rust)
            .collect()
    }

    /// Return whether this prepared pattern has at least one occurrence in a `Molecule` or
    /// `PreparedSubstructureHost`.
    fn has_match(&self, host: SubstructureHostLike<'_>) -> bool {
        match host {
            SubstructureHostLike::Prepared(host) => {
                self.to_rust().has_match_in_host(host.to_rust())
            }
            SubstructureHostLike::Molecule(host) => self.to_rust().has_match(host.to_rust()),
        }
    }

    /// Return the first occurrence of this prepared pattern in a `Molecule` or
    /// `PreparedSubstructureHost`, or `None` if no match exists.
    fn first_match(&self, host: SubstructureHostLike<'_>) -> Option<MoleculeCorrespondence> {
        let correspondence = match host {
            SubstructureHostLike::Prepared(host) => {
                self.to_rust().first_match_in_host(host.to_rust())
            }
            SubstructureHostLike::Molecule(host) => self.to_rust().first_match(host.to_rust()),
        };
        correspondence.map(MoleculeCorrespondence::from_rust)
    }

    /// Count all occurrences of this prepared pattern in a `Molecule` or
    /// `PreparedSubstructureHost` without materializing correspondences.
    fn count_matches(&self, host: SubstructureHostLike<'_>) -> u64 {
        match host {
            SubstructureHostLike::Prepared(host) => {
                self.to_rust().count_matches_in_host(host.to_rust())
            }
            SubstructureHostLike::Molecule(host) => self.to_rust().count_matches(host.to_rust()),
        }
    }

    fn __repr__(&self) -> String {
        let pattern = self.to_rust().pattern();
        format!(
            "PreparedSubstructurePattern(atoms={}, bonds={}, config={})",
            pattern.atoms().count(),
            pattern.bonds().count(),
            self.config().__repr__(),
        )
    }
}

impl PreparedSubstructurePattern {
    pub(crate) fn from_rust(pattern: GraphIrPreparedSubstructurePattern) -> Self {
        Self(pattern)
    }

    pub(crate) fn to_rust(&self) -> &GraphIrPreparedSubstructurePattern {
        &self.0
    }
}

#[cfg(test)]
mod tests {
    use std::iter;

    use rstest::rstest;
    use umol_chem::element::Element as ChemElement;
    use umol_graph_ir::ir::{
        AtomDelta as GraphIrAtomDelta, AtomFieldChange as GraphIrAtomFieldChange,
        AtomForm as GraphIrAtomForm, AtomId as GraphIrAtomId, BondForm as GraphIrBondForm,
        Delta as GraphIrDelta, Deltas as GraphIrDeltas, Molecule as GraphIrMolecule,
        MoleculeEntries as GraphIrMoleculeEntries, NumForm as GraphIrNumForm,
        Reaction as GraphIrReaction,
    };

    use super::*;

    #[rstest]
    fn test_substructure_search_config_new() {
        let config = SubstructureSearchConfig::new(
            SubstructureMatchAlgorithm::Incidence(),
            SubgraphIsomorphismAlgorithm::Glasgow(),
            RelevantCycleEnumerationAlgorithm::Vismara(),
            AutomorphismAlgorithm::Nauty(),
            true,
            true,
        );

        assert_eq!(
            config.match_algorithm(),
            SubstructureMatchAlgorithm::Incidence()
        );
        assert_eq!(
            config.subgraph_isomorphism_algorithm(),
            SubgraphIsomorphismAlgorithm::Glasgow()
        );
        assert_eq!(
            config.relevant_cycle_algorithm(),
            RelevantCycleEnumerationAlgorithm::Vismara()
        );
        assert_eq!(
            config.automorphism_algorithm(),
            AutomorphismAlgorithm::Nauty()
        );
        assert!(config.reduce_pattern_symmetry());
        assert!(config.reduce_host_symmetry());
        assert_eq!(
            config.__repr__(),
            concat!(
                "SubstructureSearchConfig(",
                "match_algorithm=SubstructureMatchAlgorithm.Incidence(), ",
                "subgraph_isomorphism_algorithm=SubgraphIsomorphismAlgorithm.Glasgow(), ",
                "relevant_cycle_algorithm=RelevantCycleEnumerationAlgorithm.Vismara(), ",
                "automorphism_algorithm=AutomorphismAlgorithm.Nauty(), ",
                "reduce_pattern_symmetry=True, ",
                "reduce_host_symmetry=True)"
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
                SubgraphIsomorphismAlgorithm::BitsetRi(),
                RelevantCycleEnumerationAlgorithm::Vismara(),
                AutomorphismAlgorithm::TrailIr(),
                false,
                false,
            )
        );
    }

    #[rstest]
    #[case::default(
        GraphIrSubstructureMatchAlgorithm::GraphAndOverlays,
        GraphCoreSubgraphIsomorphismAlgorithm::BitsetRi,
        GraphCoreRelevantCycleEnumerationAlgorithm::Vismara,
        GraphCoreAutomorphismAlgorithm::TrailIr,
        false,
        false,
        SubstructureSearchConfig::default()
    )]
    #[case::incidence_glasgow_sym(
        GraphIrSubstructureMatchAlgorithm::Incidence,
        GraphCoreSubgraphIsomorphismAlgorithm::Glasgow,
        GraphCoreRelevantCycleEnumerationAlgorithm::Vismara,
        GraphCoreAutomorphismAlgorithm::Nauty,
        true,
        true,
        SubstructureSearchConfig::new(
            SubstructureMatchAlgorithm::Incidence(),
            SubgraphIsomorphismAlgorithm::Glasgow(),
            RelevantCycleEnumerationAlgorithm::Vismara(),
            AutomorphismAlgorithm::Nauty(),
            true,
            true,
        )
    )]
    fn test_substructure_search_config_from_rust(
        #[case] match_algorithm: GraphIrSubstructureMatchAlgorithm,
        #[case] subgraph_isomorphism_algorithm: GraphCoreSubgraphIsomorphismAlgorithm,
        #[case] relevant_cycle_algorithm: GraphCoreRelevantCycleEnumerationAlgorithm,
        #[case] automorphism_algorithm: GraphCoreAutomorphismAlgorithm,
        #[case] reduce_pattern_symmetry: bool,
        #[case] reduce_host_symmetry: bool,
        #[case] expected: SubstructureSearchConfig,
    ) {
        assert_eq!(
            SubstructureSearchConfig::from_rust(GraphIrSubstructureMatchConfig {
                match_algorithm,
                subgraph_isomorphism_algorithm,
                relevant_cycle_algorithm,
                automorphism_algorithm,
                reduce_pattern_symmetry,
                reduce_host_symmetry,
            }),
            expected
        );
    }

    #[rstest]
    #[case::default(
        SubstructureSearchConfig::default(),
        GraphIrSubstructureMatchConfig {
            match_algorithm: GraphIrSubstructureMatchAlgorithm::GraphAndOverlays,
            subgraph_isomorphism_algorithm: GraphCoreSubgraphIsomorphismAlgorithm::BitsetRi,
            relevant_cycle_algorithm: GraphCoreRelevantCycleEnumerationAlgorithm::Vismara,
            automorphism_algorithm: GraphCoreAutomorphismAlgorithm::TrailIr,
            reduce_pattern_symmetry: false,
            reduce_host_symmetry: false,
        }
    )]
    #[case::incidence_glasgow_sym(
        SubstructureSearchConfig::new(
            SubstructureMatchAlgorithm::Incidence(),
            SubgraphIsomorphismAlgorithm::Glasgow(),
            RelevantCycleEnumerationAlgorithm::Vismara(),
            AutomorphismAlgorithm::Nauty(),
            true,
            true,
        ),
        GraphIrSubstructureMatchConfig {
            match_algorithm: GraphIrSubstructureMatchAlgorithm::Incidence,
            subgraph_isomorphism_algorithm: GraphCoreSubgraphIsomorphismAlgorithm::Glasgow,
            relevant_cycle_algorithm: GraphCoreRelevantCycleEnumerationAlgorithm::Vismara,
            automorphism_algorithm: GraphCoreAutomorphismAlgorithm::Nauty,
            reduce_pattern_symmetry: true,
            reduce_host_symmetry: true,
        }
    )]
    fn test_substructure_search_config_to_rust(
        #[case] config: SubstructureSearchConfig,
        #[case] expected: GraphIrSubstructureMatchConfig,
    ) {
        assert_eq!(config.to_rust(), expected);
    }

    #[rstest]
    fn test_prepared_substructure_host_and_pattern() {
        let host_rust = GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
            atoms: vec![
                GraphIrAtomForm::from_element(ChemElement::C),
                GraphIrAtomForm::from_element(ChemElement::C),
                GraphIrAtomForm::from_element(ChemElement::O),
            ],
            bonds: vec![
                (
                    GraphIrAtomId(0),
                    GraphIrAtomId(1),
                    GraphIrBondForm::from_order(1),
                ),
                (
                    GraphIrAtomId(1),
                    GraphIrAtomId(2),
                    GraphIrBondForm::from_order(1),
                ),
            ],
            ..Default::default()
        });
        let pattern_rust = GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
            atoms: vec![
                GraphIrAtomForm::from_element(ChemElement::C),
                GraphIrAtomForm::from_element(ChemElement::C),
            ],
            bonds: vec![(
                GraphIrAtomId(0),
                GraphIrAtomId(1),
                GraphIrBondForm::from_order(1),
            )],
            ..Default::default()
        });

        Python::attach(|py| {
            let host = Py::new(py, Molecule::from_rust(host_rust.clone())).unwrap();
            let pattern = Py::new(py, Molecule::from_rust(pattern_rust.clone())).unwrap();
            let prepared_host =
                Py::new(py, PreparedSubstructureHost::new(&host.bind(py).borrow())).unwrap();
            let config = SubstructureSearchConfig::new(
                SubstructureMatchAlgorithm::GraphAndOverlays(),
                SubgraphIsomorphismAlgorithm::BitsetRi(),
                RelevantCycleEnumerationAlgorithm::Vismara(),
                AutomorphismAlgorithm::TrailIr(),
                true,
                false,
            );
            let prepared_pattern =
                PreparedSubstructurePattern::new(&pattern.bind(py).borrow(), Some(config));

            // Mutating source molecules after preparation must not affect prepared snapshots.
            *host.bind(py).borrow_mut().to_rust_mut() = GraphIrMolecule::new();
            *pattern.bind(py).borrow_mut().to_rust_mut() = GraphIrMolecule::new();

            assert_eq!(
                prepared_host.bind(py).borrow().molecule().to_rust(),
                &host_rust
            );
            assert_eq!(prepared_pattern.pattern().to_rust(), &pattern_rust);
            assert_eq!(prepared_pattern.config(), config);
            assert_eq!(
                prepared_host.bind(py).borrow().__repr__(),
                "PreparedSubstructureHost(atoms=3, bonds=2)"
            );
            assert_eq!(
                prepared_pattern.__repr__(),
                format!(
                    "PreparedSubstructurePattern(atoms=2, bonds=1, config={})",
                    config.__repr__()
                )
            );

            let fresh_host = Py::new(py, Molecule::from_rust(host_rust.clone())).unwrap();
            let matches_prepared = prepared_pattern.matches(SubstructureHostLike::Prepared(
                prepared_host.bind(py).borrow(),
            ));
            let matches_mol = prepared_pattern
                .matches(SubstructureHostLike::Molecule(fresh_host.bind(py).borrow()));
            assert_eq!(matches_prepared.len(), 1);
            assert_eq!(matches_prepared, matches_mol);
            assert!(prepared_pattern.has_match(SubstructureHostLike::Prepared(
                prepared_host.bind(py).borrow()
            )));
            assert_eq!(
                prepared_pattern.first_match(SubstructureHostLike::Prepared(
                    prepared_host.bind(py).borrow()
                )),
                matches_prepared.first().cloned()
            );
            assert_eq!(
                prepared_pattern.count_matches(SubstructureHostLike::Prepared(
                    prepared_host.bind(py).borrow()
                )),
                1
            );

            let reaction_rust = GraphIrReaction::new(
                pattern_rust,
                GraphIrDeltas::from_iter([GraphIrDelta::Atom(GraphIrAtomDelta::ModifyField {
                    id: GraphIrAtomId(0),
                    change: GraphIrAtomFieldChange::Charge {
                        old: GraphIrNumForm::Undetermined,
                        new: GraphIrNumForm::Lit(1),
                    },
                })]),
            );
            let reaction = Reaction::from_rust(py, reaction_rust).unwrap();
            let prepared_reaction = PreparedReaction::new(py, &reaction, None).unwrap();
            let products_direct = prepared_host
                .bind(py)
                .borrow()
                .react(py, &reaction, None)
                .unwrap();
            let products_prepared = prepared_host
                .bind(py)
                .borrow()
                .react_prepared(py, &prepared_reaction)
                .unwrap();
            let direct_sets: Vec<Vec<GraphIrMolecule>> = iter::from_fn(|| {
                products_direct
                    .borrow_mut(py)
                    .__next__()
                    .unwrap()
                    .map(|set| set.into_iter().map(|m| m.to_rust().clone()).collect())
            })
            .collect();
            let prepared_sets: Vec<Vec<GraphIrMolecule>> = iter::from_fn(|| {
                products_prepared
                    .borrow_mut(py)
                    .__next__()
                    .unwrap()
                    .map(|set| set.into_iter().map(|m| m.to_rust().clone()).collect())
            })
            .collect();
            assert_eq!(direct_sets.len(), 2);
            assert_eq!(prepared_sets, direct_sets);
        });
    }
}
