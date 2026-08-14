//! `Reaction` — an owned Python facade over `umol_graph_ir::ir::Reaction`.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use std::str::FromStr;

use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use umol_graph::fingerprint::featurize_reaction;
use umol_graph::ingest::ingest_reaction_smiles_with;
use umol_graph::ops::model::{
    ChemistryModel as GraphChemistryModel, ValenceModel as GraphValenceModel,
};
use umol_graph::ops::resolve::ResolveConfig as GraphResolveConfig;
use umol_graph_core::CommonSubgraphEnumerationAlgorithm as GraphCoreCommonSubgraphEnumerationAlgorithm;
#[cfg(test)]
use umol_graph_core::{
    Correspondence,
    RelevantCycleEnumerationAlgorithm as GraphCoreRelevantCycleEnumerationAlgorithm,
    SubgraphIsomorphismAlgorithm as GraphCoreSubgraphIsomorphismAlgorithm,
};
use umol_graph_ir::dsl::ReactionDsl as GraphIrReactionDsl;
#[cfg(test)]
use umol_graph_ir::ir::SubstructureMatchAlgorithm as GraphIrSubstructureMatchAlgorithm;
use umol_graph_ir::ir::{
    ApplyError as GraphIrApplyError, AtomId, FromIr, IntoIr, Reaction as GraphIrReaction,
    ReactionApplicationIter as GraphIrReactionApplicationIter,
    ReactionDerivation as GraphIrReactionDerivation,
    ReactionProductsIter as GraphIrReactionProductsIter,
    SubstructureMatchConfig as GraphIrSubstructureMatchConfig,
};
use umol_io::smiles::SmilesIoConfig as IoSmilesIoConfig;

use crate::algorithm::{
    CommonSubgraphEnumerationAlgorithm, RelevantCycleEnumerationAlgorithm,
    SubgraphIsomorphismAlgorithm, SubstructureMatchAlgorithm,
};
use crate::correspondence::{
    Correspondence as PyCorrespondence, MoleculeCorrespondence as PyMoleculeCorrespondence,
};
use crate::defaults::ReactionDefaults;
use crate::delta::Deltas;
use crate::error::{
    contradiction_error, fingerprint_error, metadata_error, parse_error,
    reaction_smiles_input_error, transaction_error, InvalidStructureError,
};
use crate::fingerprint::config::ReactionCombinedFingerprintConfig;
use crate::fingerprint::reaction::ReactionCombinedFingerprint;
use crate::metadata::ReactionMetadata;
use crate::model::ChemistryModel;
use crate::molecule::Molecule;
use crate::reaction_span::ReactionSpan as PyReactionSpan;
use crate::resolve::ResolveConfig;
use crate::smiles::SmilesIoConfig;

/// Algorithm used to enumerate common subgraphs during reaction composition.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReactionCompositionConfig {
    common_subgraph_enumeration_algorithm: CommonSubgraphEnumerationAlgorithm,
}

impl Default for ReactionCompositionConfig {
    fn default() -> Self {
        Self::from_rust(GraphCoreCommonSubgraphEnumerationAlgorithm::DirectBacktracking)
    }
}

#[pymethods]
impl ReactionCompositionConfig {
    #[new]
    #[pyo3(signature = (
        *,
        common_subgraph_enumeration_algorithm=
            CommonSubgraphEnumerationAlgorithm::DirectBacktracking(),
    ))]
    fn new(common_subgraph_enumeration_algorithm: CommonSubgraphEnumerationAlgorithm) -> Self {
        Self {
            common_subgraph_enumeration_algorithm,
        }
    }

    #[staticmethod]
    fn default() -> Self {
        Default::default()
    }

    #[getter]
    fn common_subgraph_enumeration_algorithm(&self) -> CommonSubgraphEnumerationAlgorithm {
        self.common_subgraph_enumeration_algorithm
    }

    fn __repr__(&self) -> String {
        format!(
            "ReactionCompositionConfig(common_subgraph_enumeration_algorithm={})",
            self.common_subgraph_enumeration_algorithm.repr(),
        )
    }
}

impl ReactionCompositionConfig {
    pub(crate) fn from_rust(
        common_subgraph_enumeration_algorithm: GraphCoreCommonSubgraphEnumerationAlgorithm,
    ) -> Self {
        Self {
            common_subgraph_enumeration_algorithm: CommonSubgraphEnumerationAlgorithm::from_rust(
                common_subgraph_enumeration_algorithm,
            ),
        }
    }

    pub(crate) fn to_rust(self) -> GraphCoreCommonSubgraphEnumerationAlgorithm {
        self.common_subgraph_enumeration_algorithm.to_rust()
    }
}

/// Algorithms used to enumerate matches for reaction application.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReactionApplicationConfig {
    match_algorithm: SubstructureMatchAlgorithm,
    subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm,
    relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
}

impl Default for ReactionApplicationConfig {
    fn default() -> Self {
        Self {
            match_algorithm: SubstructureMatchAlgorithm::GraphAndOverlays(),
            subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::Vf2Rdkit(),
            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara(),
        }
    }
}

#[pymethods]
impl ReactionApplicationConfig {
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
            "ReactionApplicationConfig(match_algorithm={}, subgraph_isomorphism_algorithm={}, relevant_cycle_algorithm={})",
            self.match_algorithm.repr(),
            self.subgraph_isomorphism_algorithm.repr(),
            self.relevant_cycle_algorithm.repr(),
        )
    }
}

impl ReactionApplicationConfig {
    #[allow(
        dead_code,
        reason = "Rust-to-Python conversion API for configured reaction application"
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
        }
    }
}

/// A reaction whose molecule and delta components remain live Python values.
#[pyclass]
pub struct Reaction {
    lhs: Py<Molecule>,
    deltas: Py<Deltas>,
}

#[pymethods]
impl Reaction {
    /// Build a reaction from detached component snapshots.
    #[new]
    #[pyo3(signature = (lhs=None, deltas=None))]
    fn new(
        py: Python<'_>,
        lhs: Option<Py<Molecule>>,
        deltas: Option<Py<Deltas>>,
    ) -> PyResult<Self> {
        Self::from_rust(
            py,
            GraphIrReaction::new(
                lhs.map(|value| value.bind(py).borrow().to_rust().clone())
                    .unwrap_or_default(),
                deltas
                    .map(|value| value.bind(py).borrow().to_rust().clone())
                    .unwrap_or_default(),
            ),
        )
    }

    /// Parse a reaction from its EDN representation.
    #[staticmethod]
    #[pyo3(signature = (text, *, defaults=None))]
    fn parse(py: Python<'_>, text: &str, defaults: Option<ReactionDefaults>) -> PyResult<Self> {
        let defaults = defaults.unwrap_or_else(ReactionDefaults::new);
        let reaction = GraphIrReactionDsl::from_str(text)
            .map_err(parse_error)?
            .into_ir(defaults.to_rust());
        Self::from_rust(py, reaction)
    }

    /// Parse a reaction and return `(reaction, metadata)`, retaining lhs and
    /// delta entity keywords and atom aliases for metadata-preserving
    /// rendering.
    #[staticmethod]
    #[pyo3(signature = (text, *, defaults=None))]
    fn parse_with_metadata(
        py: Python<'_>,
        text: &str,
        defaults: Option<ReactionDefaults>,
    ) -> PyResult<(Self, ReactionMetadata)> {
        let defaults = defaults.unwrap_or_else(ReactionDefaults::new);
        let dsl = GraphIrReactionDsl::from_str(text).map_err(parse_error)?;
        let metadata = ReactionMetadata::from_rust(dsl.metadata().clone());
        Ok((
            Self::from_rust(py, dsl.into_ir(defaults.to_rust()))?,
            metadata,
        ))
    }

    /// Render a canonical positional DSL representation without entity
    /// keywords or atom aliases.
    #[pyo3(signature = (*, defaults=None))]
    fn render(&self, py: Python<'_>, defaults: Option<ReactionDefaults>) -> String {
        let defaults = defaults.unwrap_or_else(ReactionDefaults::new);
        GraphIrReactionDsl::from_ir(&self.to_rust(py), defaults.to_rust()).to_string()
    }

    /// Render a canonical DSL representation with persistent metadata.
    ///
    /// Raises `MetadataError` if the detached lhs or delta metadata is not
    /// coherent with this reaction.
    #[pyo3(signature = (metadata, *, defaults=None))]
    fn render_with_metadata(
        &self,
        py: Python<'_>,
        metadata: &ReactionMetadata,
        defaults: Option<ReactionDefaults>,
    ) -> PyResult<String> {
        let defaults = defaults.unwrap_or_else(ReactionDefaults::new);
        let lowered = GraphIrReactionDsl::from_ir(&self.to_rust(py), defaults.to_rust())
            .into_parts()
            .0;
        GraphIrReactionDsl::new(lowered, metadata.to_rust().clone())
            .map(|dsl| dsl.to_string())
            .map_err(metadata_error)
    }

    /// Construct a reaction by comparing two molecule snapshots under an atom correspondence.
    #[staticmethod]
    fn from_sides(
        py: Python<'_>,
        lhs: Py<Molecule>,
        rhs: Py<Molecule>,
        atom_correspondence: &PyCorrespondence,
    ) -> PyResult<Self> {
        let lhs = lhs.bind(py).borrow().to_rust().clone();
        let rhs = rhs.bind(py).borrow().to_rust().clone();
        let atom_correspondence = atom_correspondence.to_rust::<AtomId>();

        let reaction =
            GraphIrReaction::from_sides(lhs, rhs, atom_correspondence).ok_or_else(|| {
                PyValueError::new_err("atom correspondence is incompatible with the reaction sides")
            })?;
        Self::from_rust(py, reaction)
    }

    /// Ingest a determined reaction from reaction SMILES under explicit IO,
    /// chemistry, and resolution policies.
    #[staticmethod]
    #[pyo3(signature = (
        source,
        *,
        io_config=None,
        chemistry_model=None,
        resolve_config=None,
    ))]
    fn from_reaction_smiles(
        py: Python<'_>,
        source: &str,
        io_config: Option<SmilesIoConfig>,
        chemistry_model: Option<ChemistryModel>,
        resolve_config: Option<ResolveConfig>,
    ) -> PyResult<Self> {
        let io_config =
            io_config.map_or_else(IoSmilesIoConfig::opensmiles, SmilesIoConfig::to_rust);
        let chemistry_model = chemistry_model.map_or_else(
            || GraphChemistryModel {
                valence: GraphValenceModel::smiles(),
                ..GraphChemistryModel::default()
            },
            |model| model.to_rust(),
        );
        let resolve_config =
            resolve_config.map_or_else(GraphResolveConfig::default, ResolveConfig::to_rust);
        let reaction =
            ingest_reaction_smiles_with(source, &io_config, &chemistry_model, &resolve_config)
                .map_err(reaction_smiles_input_error)?;

        Self::from_rust(py, reaction)
    }

    /// The live left-hand molecule component.
    #[getter]
    fn lhs(&self, py: Python<'_>) -> Py<Molecule> {
        self.lhs.clone_ref(py)
    }

    /// Replace the left-hand molecule with a detached snapshot.
    #[setter]
    fn set_lhs(slf: Py<Self>, py: Python<'_>, value: Py<Molecule>) -> PyResult<()> {
        let resolved = Py::new(
            py,
            Molecule::from_rust(value.bind(py).borrow().to_rust().clone()),
        )?;
        slf.borrow_mut(py).lhs = resolved;
        Ok(())
    }

    /// The live delta component.
    #[getter]
    fn deltas(&self, py: Python<'_>) -> Py<Deltas> {
        self.deltas.clone_ref(py)
    }

    /// Replace the deltas with a detached snapshot.
    #[setter]
    fn set_deltas(slf: Py<Self>, py: Python<'_>, value: Py<Deltas>) -> PyResult<()> {
        let resolved = Py::new(
            py,
            Deltas::from_rust(value.bind(py).borrow().to_rust().clone()),
        )?;
        slf.borrow_mut(py).deltas = resolved;
        Ok(())
    }

    /// Materialize the superimposed reaction span.
    ///
    /// Raises `ContradictionError` when the deltas are internally inconsistent or cannot form a
    /// structurally intact right-hand molecule.
    fn to_reaction_span(&self, py: Python<'_>) -> PyResult<PyReactionSpan> {
        self.to_rust(py)
            .to_reaction_span()
            .map(PyReactionSpan::from_rust)
            .map_err(contradiction_error)
    }

    /// Return the reverse reaction in the product's compacted id space.
    fn reverse(&self, py: Python<'_>) -> PyResult<Self> {
        let reaction = self.to_rust(py).reverse().map_err(contradiction_error)?;
        Self::from_rust(py, reaction)
    }

    /// Return the sequential composites with another reaction.
    #[pyo3(signature = (other, *, config=None))]
    fn compose(
        &self,
        py: Python<'_>,
        other: &Self,
        config: Option<ReactionCompositionConfig>,
    ) -> PyResult<Vec<Self>> {
        let first = self.to_rust(py);
        let second = other.to_rust(py);
        let algorithm = config.unwrap_or_default().to_rust();

        first
            .compose(&second, algorithm)
            .into_iter()
            .map(|reaction| Self::from_rust(py, reaction))
            .collect()
    }

    /// Return one derivation per successful match through a one-shot iterator.
    ///
    /// Matching is eager; derivation construction is lazy. The iterator owns snapshots of the
    /// reaction and host, so later mutations do not affect it. Reaction-wide precondition failures
    /// raise `InvalidStructureError` here; failures while realizing a match are raised by
    /// iteration.
    #[pyo3(signature = (host, *, config=None))]
    fn apply(
        &self,
        py: Python<'_>,
        host: Py<Molecule>,
        config: Option<ReactionApplicationConfig>,
    ) -> PyResult<Py<ReactionApplicationIter>> {
        let reaction = self.to_rust(py);
        let host = host.bind(py).borrow().to_rust().clone();
        let config = config.unwrap_or_default().to_rust();
        let application = reaction
            .apply(&host, config)
            .map_err(|error| InvalidStructureError::new_err(error.to_string()))?;

        Py::new(py, ReactionApplicationIter::from_rust(application))
    }

    /// Generate a combined fingerprint over the reactant and product sides.
    #[pyo3(signature = (*, config))]
    fn combined_fingerprint(
        &self,
        py: Python<'_>,
        config: ReactionCombinedFingerprintConfig,
    ) -> PyResult<ReactionCombinedFingerprint> {
        let (featurizer, combinator) = config.to_rust();
        featurize_reaction(&self.to_rust(py), &featurizer, combinator)
            .map(ReactionCombinedFingerprint::from_rust)
            .map_err(fingerprint_error)
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __str__(&self, py: Python<'_>) -> String {
        self.render(py, None)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let lhs = self.lhs.bind(py).repr()?.extract::<String>()?;
        let deltas = self.deltas.bind(py).repr()?.extract::<String>()?;
        Ok(format!("Reaction(lhs={lhs}, deltas={deltas})"))
    }
}

impl Reaction {
    /// Wrap a Rust reaction in fresh Python-owned components.
    pub(crate) fn from_rust(py: Python<'_>, reaction: GraphIrReaction) -> PyResult<Self> {
        Ok(Self {
            lhs: Py::new(py, Molecule::from_rust(reaction.lhs))?,
            deltas: Py::new(py, Deltas::from_rust(reaction.deltas))?,
        })
    }

    /// Snapshot the current Python-owned components as a Rust reaction.
    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrReaction {
        GraphIrReaction::new(
            self.lhs.bind(py).borrow().to_rust().clone(),
            self.deltas.bind(py).borrow().to_rust().clone(),
        )
    }
}

/// One owned firing of a reaction, exposed as an immutable result value.
#[pyclass(eq, frozen, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReactionDerivation(GraphIrReactionDerivation);

#[pymethods]
impl ReactionDerivation {
    /// The molecule matched by the reaction, as a fresh snapshot.
    #[getter]
    fn lhs(&self) -> Molecule {
        Molecule::from_rust(self.0.lhs().clone())
    }

    /// The molecule produced by the reaction, as a fresh snapshot.
    #[getter]
    fn rhs(&self) -> Molecule {
        Molecule::from_rust(self.0.rhs().clone())
    }

    /// The correspondence between the two molecule sides, as a fresh snapshot.
    #[getter]
    fn comap(&self) -> PyMoleculeCorrespondence {
        PyMoleculeCorrespondence::from_rust(self.0.comap().clone())
    }

    /// The atom-level correspondence, as a fresh snapshot.
    #[getter]
    fn atom_correspondence(&self) -> PyCorrespondence {
        PyCorrespondence::from_rust(self.0.atom_correspondence())
    }

    /// Return the reverse derivation with swapped sides and inverted correspondence.
    fn reverse(&self) -> Self {
        Self::from_rust(self.to_rust().reverse())
    }

    /// Chain this derivation onto a compatible following derivation.
    fn chain(&self, next: &Self) -> Self {
        let first = self.to_rust();
        let next = next.to_rust();
        Self::from_rust(first.chain(next))
    }

    /// Recover the reaction rule represented by this concrete firing.
    fn to_reaction(&self, py: Python<'_>) -> PyResult<Reaction> {
        Reaction::from_rust(py, self.to_rust().to_reaction())
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let lhs = Py::new(py, self.lhs())?;
        let rhs = Py::new(py, self.rhs())?;
        let comap = Py::new(py, self.comap())?;
        Ok(format!(
            "ReactionDerivation(lhs={}, rhs={}, comap={})",
            lhs.bind(py).repr()?.extract::<String>()?,
            rhs.bind(py).repr()?.extract::<String>()?,
            comap.bind(py).repr()?.extract::<String>()?,
        ))
    }
}

impl ReactionDerivation {
    pub(crate) fn from_rust(derivation: GraphIrReactionDerivation) -> Self {
        Self(derivation)
    }

    pub(crate) fn to_rust(&self) -> &GraphIrReactionDerivation {
        &self.0
    }
}

/// One-shot application results with eager matching and lazy derivation construction.
#[pyclass(skip_from_py_object)]
pub(crate) struct ReactionApplicationIter {
    inner: GraphIrReactionApplicationIter,
}

impl ReactionApplicationIter {
    pub(crate) fn from_rust(inner: GraphIrReactionApplicationIter) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl ReactionApplicationIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> PyResult<Option<ReactionDerivation>> {
        match self.inner.next() {
            Some(Ok(derivation)) => Ok(Some(ReactionDerivation::from_rust(derivation))),
            Some(Err(GraphIrApplyError::Transaction(error))) => Err(transaction_error(error)),
            Some(Err(error)) => Err(PyRuntimeError::new_err(error.to_string())),
            None => Ok(None),
        }
    }
}

/// One-shot product component collections derived lazily from reaction applications.
#[pyclass(skip_from_py_object)]
pub(crate) struct ReactionProductsIter {
    inner: GraphIrReactionProductsIter,
}

impl ReactionProductsIter {
    pub(crate) fn from_rust(inner: GraphIrReactionProductsIter) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl ReactionProductsIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> PyResult<Option<Vec<Molecule>>> {
        match self.inner.next() {
            Some(Ok(products)) => Ok(Some(
                products.into_iter().map(Molecule::from_rust).collect(),
            )),
            Some(Err(GraphIrApplyError::Transaction(error))) => Err(transaction_error(error)),
            Some(Err(error)) => Err(PyRuntimeError::new_err(error.to_string())),
            None => Ok(None),
        }
    }
}

#[cfg(test)]
mod tests {
    use pyo3::exceptions::{PyStopIteration, PyTypeError};
    use pyo3::types::{PyDict, PyList};
    use rstest::{fixture, rstest};
    use umol_chem::element::Element as ChemElement;
    use umol_graph::ingest::ingest_smiles;
    use umol_graph_ir::dsl::{
        AtomDsl as GraphIrAtomDsl, MoleculeMetadata as GraphIrMoleculeMetadata,
        ReactionMetadata as GraphIrReactionMetadata,
    };
    use umol_graph_ir::ir::{
        AromaticSystemDelta as GraphIrAromaticSystemDelta,
        AromaticSystemForm as GraphIrAromaticSystemForm,
        AromaticSystemId as GraphIrAromaticSystemId, AtomDelta as GraphIrAtomDelta,
        AtomFieldChange as GraphIrAtomFieldChange, AtomForm as GraphIrAtomForm,
        AtomId as GraphIrAtomId, BondDelta as GraphIrBondDelta,
        BondFieldChange as GraphIrBondFieldChange, BondForm as GraphIrBondForm,
        BondId as GraphIrBondId, Constraint as GraphIrConstraint,
        ConstraintDelta as GraphIrConstraintDelta, DativeBondDelta as GraphIrDativeBondDelta,
        DativeBondForm as GraphIrDativeBondForm, DativeBondId as GraphIrDativeBondId,
        Delta as GraphIrDelta, Deltas as GraphIrDeltas, Entity as GraphIrEntity,
        EntitySpan as GraphIrEntitySpan, Molecule as GraphIrMolecule,
        MoleculeConstraint as GraphIrMoleculeConstraint,
        MoleculeCorrespondence as GraphIrMoleculeCorrespondence,
        MoleculeEntries as GraphIrMoleculeEntries,
        MulticenterBondDelta as GraphIrMulticenterBondDelta,
        MulticenterBondForm as GraphIrMulticenterBondForm,
        MulticenterBondId as GraphIrMulticenterBondId,
        NoncovalentBondDelta as GraphIrNoncovalentBondDelta,
        NoncovalentBondForm as GraphIrNoncovalentBondForm,
        NoncovalentBondId as GraphIrNoncovalentBondId,
        NoncovalentBondKind as GraphIrNoncovalentBondKind, Normalize, NumForm as GraphIrNumForm,
        React as GraphIrReact, ReactionSpan as GraphIrReactionSpan,
        ReactionSpanEntries as GraphIrReactionSpanEntries,
        StereoAtomDelta as GraphIrStereoAtomDelta,
        StereoAtomFieldChange as GraphIrStereoAtomFieldChange,
        StereoAtomForm as GraphIrStereoAtomForm, StereoAtomId as GraphIrStereoAtomId,
        StereoBondDelta as GraphIrStereoBondDelta, StereoBondForm as GraphIrStereoBondForm,
        StereoBondId as GraphIrStereoBondId,
        StereoConfigurationForm as GraphIrStereoConfigurationForm,
        StereoCoset as GraphIrStereoCoset, StereoKind as GraphIrStereoKind,
        StereoLigand as GraphIrStereoLigand, StereoLigandKind as GraphIrStereoLigandKind,
    };
    use umol_graph_ir::{mol_dsl, mol_dsl_ground};

    use super::*;
    use crate::convert::into_py_variant;
    use crate::delta::Delta;
    use crate::error::{ContradictionError, MetadataError, ParseError, TransactionError};
    use crate::fingerprint::config::{
        EcfpHashScheme, HashedFingerprintConfig, RefinementRounds, WlHashScheme,
    };
    use crate::fingerprint::reaction::{
        ReactionSide, RoleTaggedHashedFeatureSet, SignedHashedFeatureSet,
    };
    use crate::ring::RingConfig;

    #[rstest]
    #[case::direct(
        CommonSubgraphEnumerationAlgorithm::DirectBacktracking(),
        "CommonSubgraphEnumerationAlgorithm.DirectBacktracking()"
    )]
    #[case::modular_product(
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking(),
        "CommonSubgraphEnumerationAlgorithm.ModularProductBacktracking()"
    )]
    fn test_reaction_composition_config_new(
        #[case] algorithm: CommonSubgraphEnumerationAlgorithm,
        #[case] expected_algorithm_repr: &str,
    ) {
        let config = ReactionCompositionConfig::new(algorithm);

        assert_eq!(config.common_subgraph_enumeration_algorithm(), algorithm);
        assert_eq!(
            config.__repr__(),
            format!(
                "ReactionCompositionConfig(\
                 common_subgraph_enumeration_algorithm={expected_algorithm_repr})"
            )
        );
        assert_eq!(config, ReactionCompositionConfig::new(algorithm));
        assert_ne!(
            config,
            ReactionCompositionConfig::new(match algorithm {
                CommonSubgraphEnumerationAlgorithm::DirectBacktracking() => {
                    CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking()
                }
                CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking() => {
                    CommonSubgraphEnumerationAlgorithm::DirectBacktracking()
                }
            })
        );
    }

    #[rstest]
    fn test_reaction_composition_config_default() {
        assert_eq!(
            ReactionCompositionConfig::default(),
            ReactionCompositionConfig::new(CommonSubgraphEnumerationAlgorithm::DirectBacktracking())
        );
    }

    #[rstest]
    #[case::direct(
        GraphCoreCommonSubgraphEnumerationAlgorithm::DirectBacktracking,
        ReactionCompositionConfig::new(CommonSubgraphEnumerationAlgorithm::DirectBacktracking())
    )]
    #[case::modular_product(
        GraphCoreCommonSubgraphEnumerationAlgorithm::ModularProductBacktracking,
        ReactionCompositionConfig::new(
            CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking()
        )
    )]
    fn test_reaction_composition_config_from_rust(
        #[case] input: GraphCoreCommonSubgraphEnumerationAlgorithm,
        #[case] expected: ReactionCompositionConfig,
    ) {
        assert_eq!(ReactionCompositionConfig::from_rust(input), expected);
    }

    #[rstest]
    #[case::direct(
        ReactionCompositionConfig::new(CommonSubgraphEnumerationAlgorithm::DirectBacktracking()),
        GraphCoreCommonSubgraphEnumerationAlgorithm::DirectBacktracking
    )]
    #[case::modular_product(
        ReactionCompositionConfig::new(
            CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking()
        ),
        GraphCoreCommonSubgraphEnumerationAlgorithm::ModularProductBacktracking
    )]
    fn test_reaction_composition_config_to_rust(
        #[case] input: ReactionCompositionConfig,
        #[case] expected: GraphCoreCommonSubgraphEnumerationAlgorithm,
    ) {
        assert_eq!(input.to_rust(), expected);
    }

    #[rstest]
    fn test_reaction_application_config_default() {
        let expected = ReactionApplicationConfig::new(
            SubstructureMatchAlgorithm::GraphAndOverlays(),
            SubgraphIsomorphismAlgorithm::Vf2Rdkit(),
            RelevantCycleEnumerationAlgorithm::Vismara(),
        );

        assert_eq!(ReactionApplicationConfig::default(), expected);
    }

    #[rstest]
    fn test_reaction_application_config_value() {
        let config = ReactionApplicationConfig::new(
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
                "ReactionApplicationConfig(",
                "match_algorithm=SubstructureMatchAlgorithm.Incidence(), ",
                "subgraph_isomorphism_algorithm=",
                "SubgraphIsomorphismAlgorithm.ArcMatch(path_length=6), ",
                "relevant_cycle_algorithm=RelevantCycleEnumerationAlgorithm.Vismara())"
            )
        );
        assert_ne!(config, ReactionApplicationConfig::default());
    }

    #[rstest]
    #[case::default(
        GraphIrSubstructureMatchAlgorithm::GraphAndOverlays,
        GraphCoreSubgraphIsomorphismAlgorithm::Vf2Rdkit,
        GraphCoreRelevantCycleEnumerationAlgorithm::Vismara,
        ReactionApplicationConfig::default()
    )]
    #[case::incidence_arc_match(
        GraphIrSubstructureMatchAlgorithm::Incidence,
        GraphCoreSubgraphIsomorphismAlgorithm::ArcMatch { path_length: 6 },
        GraphCoreRelevantCycleEnumerationAlgorithm::Vismara,
        ReactionApplicationConfig::new(
            SubstructureMatchAlgorithm::Incidence(),
            SubgraphIsomorphismAlgorithm::ArcMatch { path_length: 6 },
            RelevantCycleEnumerationAlgorithm::Vismara(),
        ),
    )]
    fn test_reaction_application_config_from_rust(
        #[case] match_algorithm: GraphIrSubstructureMatchAlgorithm,
        #[case] subgraph_isomorphism_algorithm: GraphCoreSubgraphIsomorphismAlgorithm,
        #[case] relevant_cycle_algorithm: GraphCoreRelevantCycleEnumerationAlgorithm,
        #[case] expected: ReactionApplicationConfig,
    ) {
        assert_eq!(
            ReactionApplicationConfig::from_rust(GraphIrSubstructureMatchConfig {
                match_algorithm,
                subgraph_isomorphism_algorithm,
                relevant_cycle_algorithm,
            }),
            expected,
        );
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
    fn test_reaction_application_config_to_rust(
        #[case] subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm,
        #[case] expected_subgraph_isomorphism_algorithm: GraphCoreSubgraphIsomorphismAlgorithm,
    ) {
        let config = ReactionApplicationConfig::new(
            SubstructureMatchAlgorithm::Incidence(),
            subgraph_isomorphism_algorithm,
            RelevantCycleEnumerationAlgorithm::Vismara(),
        );

        assert_eq!(
            config.to_rust(),
            GraphIrSubstructureMatchConfig {
                match_algorithm: GraphIrSubstructureMatchAlgorithm::Incidence,
                subgraph_isomorphism_algorithm: expected_subgraph_isomorphism_algorithm,
                relevant_cycle_algorithm: GraphCoreRelevantCycleEnumerationAlgorithm::Vismara,
            }
        );
    }

    #[rstest]
    #[case::empty(None, None, GraphIrReaction::default())]
    #[case::populated(
        Some(GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
            atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
            ..Default::default()
        })),
        Some(vec![GraphIrDelta::Atom(GraphIrAtomDelta::Add {
            id: GraphIrAtomId(1),
            attributes: GraphIrAtomForm::from_element(ChemElement::O),
        })].into_iter().collect()),
        GraphIrReaction::new(
            GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                ..Default::default()
            }),
            vec![GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                id: GraphIrAtomId(1),
                attributes: GraphIrAtomForm::from_element(ChemElement::O),
            })].into_iter().collect(),
        ),
    )]
    fn test_reaction_new(
        #[case] lhs: Option<GraphIrMolecule>,
        #[case] deltas: Option<GraphIrDeltas>,
        #[case] expected: GraphIrReaction,
    ) {
        Python::attach(|py| {
            let lhs = lhs.map(|value| Py::new(py, Molecule::from_rust(value)).unwrap());
            let deltas = deltas.map(|value| Py::new(py, Deltas::from_rust(value)).unwrap());

            let reaction = Reaction::new(py, lhs, deltas).unwrap();

            assert_eq!(reaction.to_rust(py), expected);
        });
    }

    #[rstest]
    fn test_reaction_new_snapshot() {
        Python::attach(|py| {
            let lhs = Py::new(
                py,
                Molecule::from_rust(GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                    ..Default::default()
                })),
            )
            .unwrap();
            let deltas = Py::new(
                py,
                Deltas::from_rust(
                    vec![GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                        id: GraphIrAtomId(1),
                        attributes: GraphIrAtomForm::from_element(ChemElement::O),
                    })]
                    .into_iter()
                    .collect(),
                ),
            )
            .unwrap();
            let expected = GraphIrReaction::new(
                lhs.bind(py).borrow().to_rust().clone(),
                deltas.bind(py).borrow().to_rust().clone(),
            );

            let reaction =
                Reaction::new(py, Some(lhs.clone_ref(py)), Some(deltas.clone_ref(py))).unwrap();
            *lhs.bind(py).borrow_mut().to_rust_mut() = GraphIrMolecule::new();
            let delta = into_py_variant(
                py,
                Delta::from_rust(
                    py,
                    &GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                        id: GraphIrAtomId(2),
                        attributes: GraphIrAtomForm::from_element(ChemElement::N),
                    }),
                )
                .unwrap(),
            )
            .unwrap();
            deltas.bind(py).call_method1("append", (delta,)).unwrap();

            assert_eq!(reaction.to_rust(py), expected);
            assert_ne!(reaction.lhs.as_ptr(), lhs.as_ptr());
            assert_ne!(reaction.deltas.as_ptr(), deltas.as_ptr());
        });
    }

    #[rstest]
    #[case::atom_add_remove(
        r##"{:lhs {:atoms ["C" "O"]} :deltas [{:atom {:add "N"}} {:atom {:remove 1}}]}"##,
        2,
        vec![
            GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                id: GraphIrAtomId(2),
                attributes: GraphIrAtomForm::from_element(ChemElement::N),
            }),
            GraphIrDelta::Atom(GraphIrAtomDelta::Remove {
                id: GraphIrAtomId(1),
                attributes: GraphIrAtomForm::from_element(ChemElement::O),
            }),
        ],
    )]
    #[case::atom_modify(
        r##"{:lhs {:atoms ["Br#c0"]} :deltas [{:atom {:modify [0 "#c-1"]}}]}"##,
        1,
        vec![GraphIrDelta::Atom(GraphIrAtomDelta::ModifyField {
            id: GraphIrAtomId(0),
            change: GraphIrAtomFieldChange::Charge {
                old: GraphIrNumForm::Lit(0),
                new: GraphIrNumForm::Lit(-1),
            },
        })],
    )]
    #[case::stereo_modify(
        r##"{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]} :deltas [{:stereo-atom {:modify [0 "Th0"]}}]}"##,
        5,
        vec![GraphIrDelta::StereoAtom(GraphIrStereoAtomDelta::ModifyField {
            id: GraphIrStereoAtomId(0),
            change: GraphIrStereoAtomFieldChange::Configuration {
                old: GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(1),
                ),
                new: GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(0),
                ),
            },
        })],
    )]
    #[case::molecule_constraint(
        r##"{:lhs {:atoms ["C"]} :deltas [{:constraint {:add {:connected {}}}}]}"##,
        1,
        vec![GraphIrDelta::Constraint(GraphIrConstraintDelta::Add(
            GraphIrConstraint::Molecule(GraphIrMoleculeConstraint::Connected { atoms: None }),
        ))],
    )]
    fn test_reaction_parse(
        #[case] text: &str,
        #[case] atom_count: usize,
        #[case] expected_deltas: Vec<GraphIrDelta>,
    ) {
        Python::attach(|py| {
            let reaction = Reaction::parse(py, text, None).unwrap().to_rust(py);

            assert_eq!(reaction.lhs.atoms().count(), atom_count);
            assert_eq!(reaction.deltas.as_slice(), expected_deltas.as_slice());
        });
    }

    #[rstest]
    fn test_reaction_parse_error() {
        Python::attach(|py| {
            let error = Reaction::parse(py, "not edn", None).err().unwrap();

            assert!(error.is_instance_of::<ParseError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "EDN parse: unexpected token 'n' at byte 0"
            );
        });
    }

    #[rstest]
    #[case::required(
        r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add "O"}}]}"##,
        None,
        r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add "O"}}]}"##
    )]
    #[case::ground(
        r##"{:lhs {:atoms ["C#h4#v0#d0#t0#a!#m!"]} :deltas [{:atom {:add "O#n2#v0#d0#t0#a!#m!"}}]}"##,
        Some(ReactionDefaults::ground()),
        r##"{:lhs {:atoms ["C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"]} :deltas [{:atom {:add "O#i=#c0#h0#n2#u0#s#v0#d0#t0#a!#m!"}}]}"##
    )]
    fn test_reaction_parse_defaults(
        #[case] text: &str,
        #[case] defaults: Option<ReactionDefaults>,
        #[case] expected: &str,
    ) {
        Python::attach(|py| {
            assert_eq!(
                Reaction::parse(py, text, defaults).unwrap().to_rust(py),
                expected.parse::<GraphIrReaction>().unwrap()
            );
        });
    }

    #[rstest]
    fn test_reaction_parse_with_metadata() {
        Python::attach(|py| {
            let (reaction, metadata) = Reaction::parse_with_metadata(
                py,
                concat!(
                    r#"{:lhs {:atoms [[:lhs :lhs-c]] :atom-aliases [:lhs-c "C"]} "#,
                    r#":atom-aliases [:delta-o "O"] "#,
                    r#":deltas [{:atom {:add [:added :delta-o]}}]}"#,
                ),
                None,
            )
            .unwrap();
            let metadata = metadata.to_rust();

            assert_eq!(
                reaction.to_rust(py),
                r#"{:lhs {:atoms ["C"]} :deltas [{:atom {:add "O"}}]}"#
                    .parse()
                    .unwrap()
            );
            assert_eq!(
                metadata
                    .lhs()
                    .keyword(GraphIrEntity::Atom(GraphIrAtomId(0))),
                Some("lhs")
            );
            assert_eq!(
                metadata.delta_keyword(GraphIrEntity::Atom(GraphIrAtomId(1))),
                Some("added")
            );
            assert_eq!(
                metadata.lhs().atom_alias("lhs-c"),
                Some(&GraphIrAtomDsl(GraphIrAtomForm::from_element(
                    ChemElement::C
                )))
            );
            assert_eq!(
                metadata.atom_alias("delta-o"),
                Some(&GraphIrAtomDsl(GraphIrAtomForm::from_element(
                    ChemElement::O
                )))
            );
        });
    }

    #[rstest]
    fn test_reaction_parse_with_metadata_defaults() {
        Python::attach(|py| {
            let (reaction, metadata) = Reaction::parse_with_metadata(
                py,
                concat!(
                    r#"{:lhs {:atoms ["C#h4#v0#d0#t0#a!#m!"]} "#,
                    r#":deltas [{:atom {:add "O#n2#v0#d0#t0#a!#m!"}}]}"#,
                ),
                Some(ReactionDefaults::ground()),
            )
            .unwrap();

            assert_eq!(
                reaction.to_rust(py),
                concat!(
                    r#"{:lhs {:atoms ["C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"]} "#,
                    r#":deltas [{:atom {:add "O#i=#c0#h0#n2#u0#s#v0#d0#t0#a!#m!"}}]}"#,
                )
                .parse()
                .unwrap()
            );
            assert_eq!(
                metadata,
                ReactionMetadata::from_rust(GraphIrReactionMetadata::default())
            );
        });
    }

    #[rstest]
    #[case::required(
        r#"{:lhs {:atoms ["C"]} :deltas [{:atom {:add "O"}}]}"#.parse().unwrap(),
        None,
        r#"{:deltas [{:atom {:add "O"}}] :lhs {:atoms ["C"] :bonds []}}"#
    )]
    #[case::ground(
        concat!(
            r#"{:lhs {:atoms ["C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"]} "#,
            r#":deltas [{:atom {:add "O#i=#c0#h0#n2#u0#s#v0#d0#t0#a!#m!"}}]}"#,
        ).parse().unwrap(),
        Some(ReactionDefaults::ground()),
        concat!(
            r#"{:deltas [{:atom {:add "O#n2#v0#d0#t0#a!#m!"}}] "#,
            r#":lhs {:atoms ["C#h4#v0#d0#t0#a!#m!"] :bonds []}}"#,
        )
    )]
    fn test_reaction_render(
        #[case] reaction: GraphIrReaction,
        #[case] defaults: Option<ReactionDefaults>,
        #[case] expected: &str,
    ) {
        Python::attach(|py| {
            assert_eq!(
                Reaction::from_rust(py, reaction)
                    .unwrap()
                    .render(py, defaults),
                expected
            );
        });
    }

    #[rstest]
    fn test_reaction_render_with_metadata() {
        Python::attach(|py| {
            let reaction = Reaction::from_rust(
                py,
                r#"{:lhs {:atoms ["C"]} :deltas [{:atom {:add "O"}}]}"#
                    .parse()
                    .unwrap(),
            )
            .unwrap();
            let mut lhs = GraphIrMoleculeMetadata::new();
            lhs.set_keyword(GraphIrEntity::Atom(GraphIrAtomId(0)), "lhs")
                .unwrap();
            let mut metadata = GraphIrReactionMetadata::from(lhs);
            metadata
                .set_delta_keyword(GraphIrEntity::Atom(GraphIrAtomId(1)), "added")
                .unwrap();

            assert_eq!(
                reaction
                    .render_with_metadata(py, &ReactionMetadata::from_rust(metadata), None,)
                    .unwrap(),
                concat!(
                    r#"{:deltas [{:atom {:add [:added "O"]}}] "#,
                    r#":lhs {:atoms [[:lhs "C"]] :bonds []}}"#,
                )
            );
        });
    }

    #[rstest]
    fn test_reaction_render_with_metadata_error() {
        Python::attach(|py| {
            let reaction =
                Reaction::from_rust(py, r#"{:lhs {:atoms ["C"]} :deltas []}"#.parse().unwrap())
                    .unwrap();
            let mut metadata = GraphIrReactionMetadata::default();
            metadata
                .set_delta_keyword(GraphIrEntity::Atom(GraphIrAtomId(1)), "absent")
                .unwrap();

            let error = reaction
                .render_with_metadata(py, &ReactionMetadata::from_rust(metadata), None)
                .unwrap_err();

            assert!(error.is_instance_of::<MetadataError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "metadata entity is not introduced by an add delta: atom 1"
            );
        });
    }

    #[rstest]
    #[case::identity(
        GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
            atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
            ..Default::default()
        }),
        GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
            atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
            ..Default::default()
        }),
        vec![(0, 0)],
        GraphIrReaction::new(
            GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                ..Default::default()
            }),
            GraphIrDeltas::default(),
        ),
    )]
    #[case::partial_correspondence(
        GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
            atoms: vec![
                GraphIrAtomForm::from_element(ChemElement::C),
                GraphIrAtomForm::from_element(ChemElement::O),
            ],
            ..Default::default()
        }),
        GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
            atoms: vec![
                GraphIrAtomForm::from_element(ChemElement::C),
                GraphIrAtomForm::from_element(ChemElement::N),
            ],
            ..Default::default()
        }),
        vec![(0, 0)],
        GraphIrReaction::new(
            GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                atoms: vec![
                    GraphIrAtomForm::from_element(ChemElement::C),
                    GraphIrAtomForm::from_element(ChemElement::O),
                ],
                ..Default::default()
            }),
            vec![
                GraphIrDelta::Atom(GraphIrAtomDelta::Remove {
                    id: GraphIrAtomId(1),
                    attributes: GraphIrAtomForm::from_element(ChemElement::O),
                }),
                GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                    id: GraphIrAtomId(2),
                    attributes: GraphIrAtomForm::from_element(ChemElement::N),
                }),
            ]
            .into_iter()
            .collect(),
        ),
    )]
    #[case::bond_order(
        GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
            atoms: vec![
                GraphIrAtomForm::from_element(ChemElement::C),
                GraphIrAtomForm::from_element(ChemElement::C),
            ],
            bonds: vec![(GraphIrAtomId(0), GraphIrAtomId(1), GraphIrBondForm::from_order(1))],
            ..Default::default()
        }),
        GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
            atoms: vec![
                GraphIrAtomForm::from_element(ChemElement::C),
                GraphIrAtomForm::from_element(ChemElement::C),
            ],
            bonds: vec![(GraphIrAtomId(0), GraphIrAtomId(1), GraphIrBondForm::from_order(2))],
            ..Default::default()
        }),
        vec![(0, 0), (1, 1)],
        GraphIrReaction::new(
            GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                atoms: vec![
                    GraphIrAtomForm::from_element(ChemElement::C),
                    GraphIrAtomForm::from_element(ChemElement::C),
                ],
                bonds: vec![(GraphIrAtomId(0), GraphIrAtomId(1), GraphIrBondForm::from_order(1))],
                ..Default::default()
            }),
            vec![GraphIrDelta::Bond(GraphIrBondDelta::ModifyField {
                id: GraphIrBondId(0),
                change: GraphIrBondFieldChange::Order {
                    old: GraphIrNumForm::Lit(1),
                    new: GraphIrNumForm::Lit(2),
                },
            })]
            .into_iter()
            .collect(),
        ),
    )]
    fn test_reaction_from_sides(
        #[case] lhs: GraphIrMolecule,
        #[case] rhs: GraphIrMolecule,
        #[case] atom_pairs: Vec<(usize, usize)>,
        #[case] expected: GraphIrReaction,
    ) {
        Python::attach(|py| {
            let lhs_before = lhs.clone();
            let rhs_before = rhs.clone();
            let atom_correspondence = PyCorrespondence::from_rust(
                &Correspondence::new(
                    atom_pairs
                        .into_iter()
                        .map(|(left, right)| {
                            (GraphIrAtomId::from(left), GraphIrAtomId::from(right))
                        })
                        .collect(),
                    lhs.atoms().count(),
                    rhs.atoms().count(),
                )
                .expect("correspondence producer preserves partial-bijection invariants"),
            );
            let lhs = Py::new(py, Molecule::from_rust(lhs)).unwrap();
            let rhs = Py::new(py, Molecule::from_rust(rhs)).unwrap();

            let reaction = Reaction::from_sides(
                py,
                lhs.clone_ref(py),
                rhs.clone_ref(py),
                &atom_correspondence,
            )
            .unwrap();

            assert_eq!(reaction.to_rust(py), expected);
            assert_eq!(*lhs.bind(py).borrow().to_rust(), lhs_before);
            assert_eq!(*rhs.bind(py).borrow().to_rust(), rhs_before);
            assert_ne!(reaction.lhs.as_ptr(), lhs.as_ptr());
        });
    }

    #[rstest]
    #[case::dative_bond(
        r#"{:atoms ["N" "B"] :bonds []}"#,
        r#"{:atoms ["N" "B"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#,
        vec![(0, 0), (1, 1)],
        vec![GraphIrDelta::DativeBond(GraphIrDativeBondDelta::Add {
            id: GraphIrDativeBondId(0),
            donors: vec![GraphIrAtomId(0)],
            acceptor: GraphIrAtomId(1),
            attributes: GraphIrDativeBondForm::from_order(1),
        })],
    )]
    #[case::aromatic_system(
        r#"{:atoms ["C" "C"] :bonds []}"#,
        r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :attrs "[1,1]"}]}"#,
        vec![(0, 0), (1, 1)],
        vec![GraphIrDelta::AromaticSystem(GraphIrAromaticSystemDelta::Add {
            id: GraphIrAromaticSystemId(0),
            atoms: vec![GraphIrAtomId(0), GraphIrAtomId(1)],
            attributes: GraphIrAromaticSystemForm::from_electrons(vec![1, 1]),
        })],
    )]
    #[case::multicenter_bond(
        r#"{:atoms ["B" "H" "B"] :bonds []}"#,
        r#"{:atoms ["B" "H" "B"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :attrs "[3,5,7]"}]}"#,
        vec![(0, 0), (1, 1), (2, 2)],
        vec![GraphIrDelta::MulticenterBond(GraphIrMulticenterBondDelta::Add {
            id: GraphIrMulticenterBondId(0),
            atoms: vec![GraphIrAtomId(0), GraphIrAtomId(1), GraphIrAtomId(2)],
            attributes: GraphIrMulticenterBondForm::from_electrons(vec![3, 5, 7]),
        })],
    )]
    #[case::noncovalent_bond(
        r#"{:atoms ["O" "O"] :bonds []}"#,
        r#"{:atoms ["O" "O"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :attrs "Hbd"}]}"#,
        vec![(0, 0), (1, 1)],
        vec![GraphIrDelta::NoncovalentBond(GraphIrNoncovalentBondDelta::Add {
            id: GraphIrNoncovalentBondId(0),
            atoms: [GraphIrAtomId(0), GraphIrAtomId(1)],
            attributes: GraphIrNoncovalentBondForm::from_kind(GraphIrNoncovalentBondKind::HydrogenBond),
        })],
    )]
    #[case::stereo_atom(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds []}"#,
        r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#,
        vec![(0, 0), (1, 1), (2, 2), (3, 3), (4, 4)],
        vec![GraphIrDelta::StereoAtom(GraphIrStereoAtomDelta::Add {
            id: GraphIrStereoAtomId(0),
            site: GraphIrAtomId(0),
            ligands: vec![
                GraphIrStereoLigand::new(GraphIrAtomId(1), GraphIrStereoLigandKind::Atom),
                GraphIrStereoLigand::new(GraphIrAtomId(2), GraphIrStereoLigandKind::Atom),
                GraphIrStereoLigand::new(GraphIrAtomId(3), GraphIrStereoLigandKind::Atom),
                GraphIrStereoLigand::new(GraphIrAtomId(4), GraphIrStereoLigandKind::Atom),
            ],
            attributes: GraphIrStereoAtomForm::new(GraphIrStereoKind::Tetrahedral, GraphIrStereoCoset::Lit(1)),
        })],
    )]
    #[case::stereo_bond(
        r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]}"#,
        r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#,
        vec![(0, 0), (1, 1), (2, 2), (3, 3)],
        vec![GraphIrDelta::StereoBond(GraphIrStereoBondDelta::Add {
            id: GraphIrStereoBondId(0),
            site: GraphIrBondId(1),
            ligands: vec![
                GraphIrStereoLigand::new(GraphIrAtomId(0), GraphIrStereoLigandKind::Atom),
                GraphIrStereoLigand::new(
                    GraphIrAtomId(1),
                    GraphIrStereoLigandKind::ImplicitHydrogen,
                ),
                GraphIrStereoLigand::new(GraphIrAtomId(3), GraphIrStereoLigandKind::Atom),
                GraphIrStereoLigand::new(
                    GraphIrAtomId(2),
                    GraphIrStereoLigandKind::ImplicitHydrogen,
                ),
            ],
            attributes: GraphIrStereoBondForm::new(GraphIrStereoKind::CisTrans, GraphIrStereoCoset::Lit(1)),
        })],
    )]
    #[case::molecule_constraint(
        r#"{:atoms ["C"] :bonds []}"#,
        r#"{:atoms ["C"] :bonds [] :constraints [{:connected {}}]}"#,
        vec![(0, 0)],
        vec![GraphIrDelta::Constraint(GraphIrConstraintDelta::Add(
            GraphIrConstraint::Molecule(GraphIrMoleculeConstraint::Connected { atoms: None }),
        ))],
    )]
    fn test_reaction_from_sides_entities(
        #[case] lhs: &str,
        #[case] rhs: &str,
        #[case] atom_pairs: Vec<(usize, usize)>,
        #[case] expected_deltas: Vec<GraphIrDelta>,
    ) {
        Python::attach(|py| {
            let lhs = lhs.parse::<GraphIrMolecule>().unwrap();
            let rhs = rhs.parse::<GraphIrMolecule>().unwrap();
            let atom_correspondence = PyCorrespondence::from_rust(
                &Correspondence::new(
                    atom_pairs
                        .into_iter()
                        .map(|(left, right)| {
                            (GraphIrAtomId::from(left), GraphIrAtomId::from(right))
                        })
                        .collect(),
                    lhs.atoms().count(),
                    rhs.atoms().count(),
                )
                .expect("correspondence producer preserves partial-bijection invariants"),
            );
            let reaction = Reaction::from_sides(
                py,
                Py::new(py, Molecule::from_rust(lhs.clone())).unwrap(),
                Py::new(py, Molecule::from_rust(rhs)).unwrap(),
                &atom_correspondence,
            )
            .unwrap();

            assert_eq!(
                reaction.to_rust(py),
                GraphIrReaction::new(lhs, expected_deltas.into_iter().collect())
            );
        });
    }

    #[rstest]
    fn test_reaction_from_sides_snapshot() {
        Python::attach(|py| {
            let lhs_before = GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                atoms: vec![
                    GraphIrAtomForm::from_element(ChemElement::C),
                    GraphIrAtomForm::from_element(ChemElement::O),
                ],
                ..Default::default()
            });
            let rhs_before = GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                atoms: vec![
                    GraphIrAtomForm::from_element(ChemElement::C),
                    GraphIrAtomForm::from_element(ChemElement::N),
                ],
                ..Default::default()
            });
            let lhs = Py::new(py, Molecule::from_rust(lhs_before.clone())).unwrap();
            let rhs = Py::new(py, Molecule::from_rust(rhs_before.clone())).unwrap();
            let atom_correspondence = PyCorrespondence::from_rust(
                &Correspondence::new(vec![(GraphIrAtomId(0), GraphIrAtomId(0))], 2, 2)
                    .expect("correspondence producer preserves partial-bijection invariants"),
            );
            let reaction = Reaction::from_sides(
                py,
                lhs.clone_ref(py),
                rhs.clone_ref(py),
                &atom_correspondence,
            )
            .unwrap();
            let expected = reaction.to_rust(py);

            *lhs.bind(py).borrow_mut().to_rust_mut() = GraphIrMolecule::new();
            *rhs.bind(py).borrow_mut().to_rust_mut() = GraphIrMolecule::new();

            assert_eq!(reaction.to_rust(py), expected);
            assert_ne!(reaction.lhs.as_ptr(), lhs.as_ptr());

            *reaction.lhs.bind(py).borrow_mut().to_rust_mut() =
                GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![GraphIrAtomForm::from_element(ChemElement::F)],
                    ..Default::default()
                });
            let delta = into_py_variant(
                py,
                Delta::from_rust(
                    py,
                    &GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                        id: GraphIrAtomId(3),
                        attributes: GraphIrAtomForm::from_element(ChemElement::Cl),
                    }),
                )
                .unwrap(),
            )
            .unwrap();
            reaction
                .deltas
                .bind(py)
                .call_method1("append", (delta,))
                .unwrap();
            let changed = reaction.to_rust(py);

            assert_eq!(
                changed.lhs,
                GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![GraphIrAtomForm::from_element(ChemElement::F)],
                    ..Default::default()
                })
            );
            assert_eq!(
                changed.deltas.as_slice().last(),
                Some(&GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                    id: GraphIrAtomId(3),
                    attributes: GraphIrAtomForm::from_element(ChemElement::Cl),
                }))
            );
        });
    }

    #[rstest]
    fn test_reaction_components() {
        Python::attach(|py| {
            let reaction = Py::new(py, Reaction::new(py, None, None).unwrap()).unwrap();
            let first_lhs = reaction.bind(py).borrow().lhs(py);
            let second_lhs = reaction.bind(py).borrow().lhs(py);
            let first_deltas = reaction.bind(py).borrow().deltas(py);
            let second_deltas = reaction.bind(py).borrow().deltas(py);

            *first_lhs.bind(py).borrow_mut().to_rust_mut() =
                GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                    ..Default::default()
                });
            let delta = into_py_variant(
                py,
                Delta::from_rust(
                    py,
                    &GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                        id: GraphIrAtomId(1),
                        attributes: GraphIrAtomForm::from_element(ChemElement::O),
                    }),
                )
                .unwrap(),
            )
            .unwrap();
            first_deltas
                .bind(py)
                .call_method1("append", (delta,))
                .unwrap();

            assert_eq!(first_lhs.as_ptr(), second_lhs.as_ptr());
            assert_eq!(first_deltas.as_ptr(), second_deltas.as_ptr());
            assert_eq!(
                reaction.bind(py).borrow().to_rust(py),
                GraphIrReaction::new(
                    GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                        atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                        ..Default::default()
                    }),
                    vec![GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                        id: GraphIrAtomId(1),
                        attributes: GraphIrAtomForm::from_element(ChemElement::O),
                    })]
                    .into_iter()
                    .collect(),
                )
            );
        });
    }

    #[rstest]
    fn test_reaction_set_components() {
        Python::attach(|py| {
            let reaction = Py::new(py, Reaction::new(py, None, None).unwrap()).unwrap();
            let lhs = Py::new(
                py,
                Molecule::from_rust(GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                    ..Default::default()
                })),
            )
            .unwrap();
            let deltas = Py::new(
                py,
                Deltas::from_rust(
                    vec![GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                        id: GraphIrAtomId(1),
                        attributes: GraphIrAtomForm::from_element(ChemElement::O),
                    })]
                    .into_iter()
                    .collect(),
                ),
            )
            .unwrap();
            let expected = GraphIrReaction::new(
                lhs.bind(py).borrow().to_rust().clone(),
                deltas.bind(py).borrow().to_rust().clone(),
            );

            Reaction::set_lhs(reaction.clone_ref(py), py, lhs.clone_ref(py)).unwrap();
            Reaction::set_deltas(reaction.clone_ref(py), py, deltas.clone_ref(py)).unwrap();
            *lhs.bind(py).borrow_mut().to_rust_mut() = GraphIrMolecule::new();
            let delta = into_py_variant(
                py,
                Delta::from_rust(
                    py,
                    &GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                        id: GraphIrAtomId(2),
                        attributes: GraphIrAtomForm::from_element(ChemElement::N),
                    }),
                )
                .unwrap(),
            )
            .unwrap();
            deltas.bind(py).call_method1("append", (delta,)).unwrap();

            assert_eq!(reaction.bind(py).borrow().to_rust(py), expected);
        });
    }

    #[rstest]
    fn test_reaction_set_components_self() {
        Python::attach(|py| {
            let expected = GraphIrReaction::new(
                GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                    ..Default::default()
                }),
                vec![GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                    id: GraphIrAtomId(1),
                    attributes: GraphIrAtomForm::from_element(ChemElement::O),
                })]
                .into_iter()
                .collect(),
            );
            let reaction = Py::new(py, Reaction::from_rust(py, expected.clone()).unwrap()).unwrap();
            let own_lhs = reaction.bind(py).borrow().lhs(py);
            let own_deltas = reaction.bind(py).borrow().deltas(py);

            Reaction::set_lhs(reaction.clone_ref(py), py, own_lhs).unwrap();
            Reaction::set_deltas(reaction.clone_ref(py), py, own_deltas).unwrap();

            assert_eq!(reaction.bind(py).borrow().to_rust(py), expected);
        });
    }

    #[rstest]
    fn test_reaction_to_reaction_span() {
        Python::attach(|py| {
            let reaction = Reaction::from_rust(
                py,
                GraphIrReaction::new(
                    GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                        atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                        ..Default::default()
                    }),
                    vec![GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                        id: GraphIrAtomId(1),
                        attributes: GraphIrAtomForm::from_element(ChemElement::O),
                    })]
                    .into_iter()
                    .collect(),
                ),
            )
            .unwrap();

            assert_eq!(
                reaction.to_reaction_span(py).unwrap().to_rust(),
                &GraphIrReactionSpan::from_entries(GraphIrReactionSpanEntries {
                    atoms: vec![
                        GraphIrEntitySpan::Unchanged(GraphIrAtomForm::from_element(ChemElement::C)),
                        GraphIrEntitySpan::Added(GraphIrAtomForm::from_element(ChemElement::O)),
                    ],
                    ..Default::default()
                })
            );
        });
    }

    #[rstest]
    fn test_reaction_to_reaction_span_error() {
        Python::attach(|py| {
            let reaction = Reaction::from_rust(
                py,
                GraphIrReaction::new(
                    GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
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
                    }),
                    vec![GraphIrDelta::Bond(GraphIrBondDelta::ModifyField {
                        id: GraphIrBondId(0),
                        change: GraphIrBondFieldChange::Order {
                            old: GraphIrNumForm::Lit(2),
                            new: GraphIrNumForm::Lit(3),
                        },
                    })]
                    .into_iter()
                    .collect(),
                ),
            )
            .unwrap();

            let error = reaction.to_reaction_span(py).unwrap_err();

            assert!(error.is_instance_of::<ContradictionError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "reached a contradiction"
            );
        });
    }

    #[rstest]
    fn test_reaction_reverse() {
        Python::attach(|py| {
            let source = Reaction::parse(
                py,
                r##"{:lhs {:atoms ["C" "O"]} :deltas [{:atom {:add "N"}} {:atom {:remove 1}}]}"##,
                None,
            )
            .unwrap();
            let before = source.to_rust(py);
            let expected_deltas = before.deltas.clone().normalize().unwrap();

            let reversed = source.reverse(py).unwrap();
            let roundtrip = reversed.reverse(py).unwrap();
            let roundtrip = roundtrip.to_rust(py);

            assert_eq!(
                reversed.to_rust(py).lhs,
                GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![
                        GraphIrAtomForm::from_element(ChemElement::C),
                        GraphIrAtomForm::from_element(ChemElement::N),
                    ],
                    ..Default::default()
                })
            );
            assert_eq!(roundtrip.lhs, before.lhs);
            assert_eq!(roundtrip.deltas.normalize().unwrap(), expected_deltas);
            assert_eq!(source.to_rust(py), before);
            assert_ne!(reversed.lhs.as_ptr(), source.lhs.as_ptr());
            assert_ne!(reversed.deltas.as_ptr(), source.deltas.as_ptr());
        });
    }

    #[rstest]
    #[case::no_match(
        r##"{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+"]}}]}"##,
        r##"{:lhs {:atoms ["N#c0"]} :deltas [{:atom {:modify [0 "#c+"]}}]}"##,
        vec![
            r##"{:lhs {:atoms ["C#c0" "N#c0"]} :deltas [{:atom {:modify [0 "#c+"]}} {:atom {:modify [1 "#c+"]}}]}"##
        ]
    )]
    #[case::admissible(
        r##"{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+"]}}]}"##,
        r##"{:lhs {:atoms ["C#c+"]} :deltas [{:atom {:modify [0 "#c+2"]}}]}"##,
        vec![
            r##"{:lhs {:atoms ["C#c0" "C#c+"]} :deltas [{:atom {:modify [0 "#c+"]}} {:atom {:modify [1 "#c+2"]}}]}"##,
            r##"{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+2"]}}]}"##
        ],
    )]
    fn test_reaction_compose(
        #[case] first: &str,
        #[case] second: &str,
        #[case] expected: Vec<&str>,
    ) {
        Python::attach(|py| {
            let first = Reaction::parse(py, first, None).unwrap();
            let second = Reaction::parse(py, second, None).unwrap();
            let expected: Vec<GraphIrReaction> = expected
                .into_iter()
                .map(|reaction| GraphIrReaction::from_str(reaction).unwrap())
                .collect();

            let actual: Vec<GraphIrReaction> = first
                .compose(py, &second, None)
                .unwrap()
                .iter()
                .map(|reaction| reaction.to_rust(py))
                .collect();

            assert_eq!(actual.len(), expected.len());
            for expected in expected {
                assert!(actual.contains(&expected));
            }
        });
    }

    #[rstest]
    #[case::direct(ReactionCompositionConfig::new(
        CommonSubgraphEnumerationAlgorithm::DirectBacktracking()
    ))]
    #[case::modular_product(ReactionCompositionConfig::new(
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking()
    ))]
    fn test_reaction_compose_config(#[case] config: ReactionCompositionConfig) {
        Python::attach(|py| {
            let first = Reaction::parse(
                py,
                r##"{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+"]}}]}"##,
                None,
            )
            .unwrap();
            let second = Reaction::parse(
                py,
                r##"{:lhs {:atoms ["C#c+"]} :deltas [{:atom {:modify [0 "#c+2"]}}]}"##,
                None,
            )
            .unwrap();

            let actual: Vec<_> = first
                .compose(py, &second, Some(config))
                .unwrap()
                .into_iter()
                .map(|reaction| reaction.to_rust(py))
                .collect();
            let expected = vec![
                GraphIrReaction::from_str(
                    r##"{:lhs {:atoms ["C#c0" "C#c+"]} :deltas [{:atom {:modify [0 "#c+"]}} {:atom {:modify [1 "#c+2"]}}]}"##,
                )
                .unwrap(),
                GraphIrReaction::from_str(
                    r##"{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+2"]}}]}"##,
                )
                .unwrap(),
            ];

            assert_eq!(actual.len(), expected.len());
            for expected in expected {
                assert!(actual.contains(&expected));
            }
        });
    }

    #[rstest]
    fn test_reaction_compose_default() {
        Python::attach(|py| {
            let first = Py::new(
                py,
                Reaction::parse(
                    py,
                    r##"{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+"]}}]}"##,
                    None,
                )
                .unwrap(),
            )
            .unwrap();
            let second = Py::new(
                py,
                Reaction::parse(
                    py,
                    r##"{:lhs {:atoms ["C#c+"]} :deltas [{:atom {:modify [0 "#c+2"]}}]}"##,
                    None,
                )
                .unwrap(),
            )
            .unwrap();
            let config = Py::new(
                py,
                ReactionCompositionConfig::new(
                    CommonSubgraphEnumerationAlgorithm::DirectBacktracking(),
                ),
            )
            .unwrap();
            let kwargs = PyDict::new(py);
            kwargs.set_item("config", config).unwrap();

            let omitted: Vec<Py<Reaction>> = first
                .bind(py)
                .call_method1("compose", (second.clone_ref(py),))
                .unwrap()
                .extract()
                .unwrap();
            let explicit: Vec<Py<Reaction>> = first
                .bind(py)
                .call_method("compose", (second,), Some(&kwargs))
                .unwrap()
                .extract()
                .unwrap();
            let omitted: Vec<GraphIrReaction> = omitted
                .iter()
                .map(|reaction| reaction.bind(py).borrow().to_rust(py))
                .collect();
            let explicit: Vec<GraphIrReaction> = explicit
                .iter()
                .map(|reaction| reaction.bind(py).borrow().to_rust(py))
                .collect();

            assert_eq!(omitted, explicit);
            let expected = vec![
                GraphIrReaction::from_str(
                    r##"{:lhs {:atoms ["C#c0" "C#c+"]} :deltas [{:atom {:modify [0 "#c+"]}} {:atom {:modify [1 "#c+2"]}}]}"##,
                )
                .unwrap(),
                GraphIrReaction::from_str(
                    r##"{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+2"]}}]}"##,
                )
                .unwrap(),
            ];

            assert_eq!(omitted.len(), expected.len());
            for expected in expected {
                assert!(omitted.contains(&expected));
            }
        });
    }

    #[rstest]
    fn test_reaction_compose_snapshot() {
        Python::attach(|py| {
            let first = Reaction::parse(
                py,
                r##"{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+"]}}]}"##,
                None,
            )
            .unwrap();
            let second = Reaction::parse(
                py,
                r##"{:lhs {:atoms ["C#c+"]} :deltas [{:atom {:modify [0 "#c+2"]}}]}"##,
                None,
            )
            .unwrap();
            let first_before = first.to_rust(py);
            let second_before = second.to_rust(py);

            let _self_composites = first.compose(py, &first, None).unwrap();
            let composites = first.compose(py, &second, None).unwrap();

            assert_eq!(first.to_rust(py), first_before);
            assert_eq!(second.to_rust(py), second_before);
            assert_eq!(composites.len(), 2);
            assert_ne!(composites[0].lhs.as_ptr(), first.lhs.as_ptr());
            assert_ne!(composites[0].lhs.as_ptr(), second.lhs.as_ptr());
            assert_ne!(composites[0].deltas.as_ptr(), first.deltas.as_ptr());
            assert_ne!(composites[0].deltas.as_ptr(), second.deltas.as_ptr());
            assert_ne!(composites[0].lhs.as_ptr(), composites[1].lhs.as_ptr());
            assert_ne!(composites[0].deltas.as_ptr(), composites[1].deltas.as_ptr());

            for composite in &composites {
                *composite.lhs.bind(py).borrow_mut().to_rust_mut() =
                    GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                        atoms: vec![GraphIrAtomForm::from_element(ChemElement::F)],
                        ..Default::default()
                    });
                let delta = into_py_variant(
                    py,
                    Delta::from_rust(
                        py,
                        &GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                            id: GraphIrAtomId(8),
                            attributes: GraphIrAtomForm::from_element(ChemElement::Cl),
                        }),
                    )
                    .unwrap(),
                )
                .unwrap();
                composite
                    .deltas
                    .bind(py)
                    .call_method1("append", (delta,))
                    .unwrap();

                assert_eq!(
                    composite.to_rust(py).lhs,
                    GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                        atoms: vec![GraphIrAtomForm::from_element(ChemElement::F)],
                        ..Default::default()
                    })
                );
                assert_eq!(
                    composite.to_rust(py).deltas.as_slice().last(),
                    Some(&GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                        id: GraphIrAtomId(8),
                        attributes: GraphIrAtomForm::from_element(ChemElement::Cl),
                    }))
                );
            }

            assert_eq!(first.to_rust(py), first_before);
            assert_eq!(second.to_rust(py), second_before);
        });
    }

    #[rstest]
    fn test_reaction_apply(reaction_application: (GraphIrReaction, GraphIrMolecule)) {
        let (expected_reaction, expected_host) = reaction_application;
        Python::attach(|py| {
            let reaction = Reaction::from_rust(py, expected_reaction.clone()).unwrap();
            let host = Py::new(py, Molecule::from_rust(expected_host.clone())).unwrap();
            let application = reaction.apply(py, host.clone_ref(py), None).unwrap();

            assert_eq!(reaction.to_rust(py), expected_reaction);
            assert_eq!(host.bind(py).borrow().to_rust(), &expected_host);

            let first = application.borrow_mut(py).__next__().unwrap().unwrap();
            let second = application.borrow_mut(py).__next__().unwrap().unwrap();
            assert_eq!(application.borrow_mut(py).__next__().unwrap(), None);
            assert_eq!(
                [
                    first.rhs().to_rust().clone(),
                    second.rhs().to_rust().clone()
                ],
                [
                    GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                        atoms: vec![
                            GraphIrAtomForm::from_element(ChemElement::C).with_charge(1),
                            GraphIrAtomForm::from_element(ChemElement::C),
                        ],
                        ..Default::default()
                    }),
                    GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                        atoms: vec![
                            GraphIrAtomForm::from_element(ChemElement::C),
                            GraphIrAtomForm::from_element(ChemElement::C).with_charge(1),
                        ],
                        ..Default::default()
                    }),
                ]
            );
        });
    }

    #[rstest]
    fn test_reaction_apply_snapshot(reaction_application: (GraphIrReaction, GraphIrMolecule)) {
        let (expected_reaction, expected_host) = reaction_application;
        Python::attach(|py| {
            let mut reaction = Reaction::from_rust(py, expected_reaction).unwrap();
            let host = Py::new(py, Molecule::from_rust(expected_host)).unwrap();
            let application = reaction.apply(py, host.clone_ref(py), None).unwrap();

            *reaction.lhs.bind(py).borrow_mut().to_rust_mut() =
                GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![GraphIrAtomForm::from_element(ChemElement::N)],
                    ..Default::default()
                });
            reaction.deltas = Py::new(py, Deltas::from_rust(GraphIrDeltas::default())).unwrap();
            *host.bind(py).borrow_mut().to_rust_mut() =
                GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![GraphIrAtomForm::from_element(ChemElement::F)],
                    ..Default::default()
                });

            let products: Vec<GraphIrMolecule> = std::iter::from_fn(|| {
                application
                    .borrow_mut(py)
                    .__next__()
                    .unwrap()
                    .map(|derivation| derivation.rhs().to_rust().clone())
            })
            .collect();
            assert_eq!(
                products,
                vec![
                    GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                        atoms: vec![
                            GraphIrAtomForm::from_element(ChemElement::C).with_charge(1),
                            GraphIrAtomForm::from_element(ChemElement::C),
                        ],
                        ..Default::default()
                    }),
                    GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                        atoms: vec![
                            GraphIrAtomForm::from_element(ChemElement::C),
                            GraphIrAtomForm::from_element(ChemElement::C).with_charge(1),
                        ],
                        ..Default::default()
                    }),
                ]
            );
        });
    }

    #[rstest]
    #[case::vf2(ReactionApplicationConfig::new(
        SubstructureMatchAlgorithm::GraphAndOverlays(),
        SubgraphIsomorphismAlgorithm::Vf2(),
        RelevantCycleEnumerationAlgorithm::Vismara(),
    ))]
    #[case::ullmann(ReactionApplicationConfig::new(
        SubstructureMatchAlgorithm::GraphAndOverlays(),
        SubgraphIsomorphismAlgorithm::Ullmann(),
        RelevantCycleEnumerationAlgorithm::Vismara(),
    ))]
    #[case::ri(ReactionApplicationConfig::new(
        SubstructureMatchAlgorithm::GraphAndOverlays(),
        SubgraphIsomorphismAlgorithm::Ri(),
        RelevantCycleEnumerationAlgorithm::Vismara(),
    ))]
    #[case::arc_match(ReactionApplicationConfig::new(
        SubstructureMatchAlgorithm::GraphAndOverlays(),
        SubgraphIsomorphismAlgorithm::ArcMatch { path_length: 6 },
        RelevantCycleEnumerationAlgorithm::Vismara(),
    ))]
    #[case::vf2_rdkit(ReactionApplicationConfig::new(
        SubstructureMatchAlgorithm::GraphAndOverlays(),
        SubgraphIsomorphismAlgorithm::Vf2Rdkit(),
        RelevantCycleEnumerationAlgorithm::Vismara(),
    ))]
    #[case::ray_kirsch(ReactionApplicationConfig::new(
        SubstructureMatchAlgorithm::GraphAndOverlays(),
        SubgraphIsomorphismAlgorithm::RayKirsch(),
        RelevantCycleEnumerationAlgorithm::Vismara(),
    ))]
    #[case::incidence(ReactionApplicationConfig::new(
        SubstructureMatchAlgorithm::Incidence(),
        SubgraphIsomorphismAlgorithm::Vf2Rdkit(),
        RelevantCycleEnumerationAlgorithm::Vismara(),
    ))]
    fn test_reaction_apply_config(
        reaction_application: (GraphIrReaction, GraphIrMolecule),
        #[case] config: ReactionApplicationConfig,
    ) {
        let (reaction, host) = reaction_application;
        Python::attach(|py| {
            let reaction = Reaction::from_rust(py, reaction).unwrap();
            let host = Py::new(py, Molecule::from_rust(host)).unwrap();
            let application = reaction.apply(py, host, Some(config)).unwrap();

            let products: Vec<GraphIrMolecule> = std::iter::from_fn(|| {
                application
                    .borrow_mut(py)
                    .__next__()
                    .unwrap()
                    .map(|derivation| derivation.rhs().to_rust().clone())
            })
            .collect();
            assert_eq!(
                products,
                vec![
                    GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                        atoms: vec![
                            GraphIrAtomForm::from_element(ChemElement::C).with_charge(1),
                            GraphIrAtomForm::from_element(ChemElement::C),
                        ],
                        ..Default::default()
                    }),
                    GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                        atoms: vec![
                            GraphIrAtomForm::from_element(ChemElement::C),
                            GraphIrAtomForm::from_element(ChemElement::C).with_charge(1),
                        ],
                        ..Default::default()
                    }),
                ]
            );
        });
    }

    #[rstest]
    fn test_reaction_apply_error() {
        Python::attach(|py| {
            let reaction = Reaction::from_rust(
                py,
                GraphIrReaction::new(
                    GraphIrMolecule::default(),
                    [GraphIrDelta::Atom(GraphIrAtomDelta::Remove {
                        id: GraphIrAtomId(0),
                        attributes: GraphIrAtomForm::from_element(ChemElement::C),
                    })]
                    .into_iter()
                    .collect(),
                ),
            )
            .unwrap();
            let host = Py::new(py, Molecule::from_rust(GraphIrMolecule::default())).unwrap();

            let error = reaction.apply(py, host, None).err().unwrap();

            assert!(error.is_instance_of::<InvalidStructureError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "reaction references unavailable entity Atom(AtomId(0))"
            );
        });
    }

    #[fixture]
    fn ethanol_deoxygenation() -> GraphIrReaction {
        let ethanol = ingest_smiles("CCO").unwrap();
        let oxygen = ethanol.atom(GraphIrAtomId(2)).attributes.clone();
        let bond = ethanol.bond(GraphIrBondId(1)).attributes.clone();
        GraphIrReaction::new(
            ethanol,
            GraphIrDeltas::from_iter([
                GraphIrDelta::Atom(GraphIrAtomDelta::Remove {
                    id: GraphIrAtomId(2),
                    attributes: oxygen,
                }),
                GraphIrDelta::Bond(GraphIrBondDelta::Remove {
                    id: GraphIrBondId(1),
                    atoms: [GraphIrAtomId(1), GraphIrAtomId(2)],
                    attributes: bond,
                }),
            ]),
        )
    }

    #[fixture]
    fn ethanol_identity() -> GraphIrReaction {
        GraphIrReaction::new(ingest_smiles("CCO").unwrap(), GraphIrDeltas::new())
    }

    #[rstest]
    #[case::morgan(
        ReactionCombinedFingerprintConfig::Difference {
            molecule: HashedFingerprintConfig::Morgan {
                radius: 2,
                ring_config: RingConfig::default(),
            },
        },
        vec![
            (864662311, -1),
            (1535166686, -1),
            (2245384272, -1),
            (2246997334, 1),
            (3542456614, -1),
            (3548082732, 1),
            (4018048386, -1),
        ]
    )]
    #[case::ecfp(
        ReactionCombinedFingerprintConfig::Difference {
            molecule: HashedFingerprintConfig::Ecfp {
                radius: 2,
                hashing_scheme: EcfpHashScheme::Xxh3Width64V1(),
                ring_config: RingConfig::default(),
            },
        },
        vec![
            (63839236075656913, -1),
            (896060437578512973, 1),
            (1189585227353469813, -1),
            (3822471596818936039, -1),
            (13327007941213506523, 1),
            (13652293261850732425, -1),
            (15001976065402722634, -1),
        ]
    )]
    #[case::wl(
        ReactionCombinedFingerprintConfig::Difference {
            molecule: HashedFingerprintConfig::Wl {
                rounds: RefinementRounds::Fixed { rounds: 3 },
                hashing_scheme: WlHashScheme::Xxh3SortedWidth64V1(),
            },
        },
        vec![
            (2520347590860685079, -1),
            (3352603313223549703, -1),
            (4152249898001161146, -1),
            (5807737097854608645, -1),
            (6786829771653353480, 1),
            (7404535559284410087, 1),
            (8754482138526219790, 1),
            (11986000156817227245, -1),
            (12849090138728295812, 1),
            (12895020514073294021, -1),
            (13932567567828606490, -1),
            (16456488943967932267, 1),
            (17305796300852423160, -1),
            (17417400371411086222, -1),
        ]
    )]
    fn test_reaction_combined_fingerprint_difference(
        ethanol_deoxygenation: GraphIrReaction,
        #[case] config: ReactionCombinedFingerprintConfig,
        #[case] expected_entries: Vec<(u128, i32)>,
    ) {
        Python::attach(|py| {
            let reaction = Reaction::from_rust(py, ethanol_deoxygenation).unwrap();
            let fingerprint = reaction.combined_fingerprint(py, config).unwrap();
            let fingerprint = Py::new(py, fingerprint).unwrap();
            let fingerprint = fingerprint.bind(py).as_any();
            let features = fingerprint.getattr("features").unwrap();
            let entries = features.getattr("entries").unwrap();

            assert!(features.is_instance_of::<SignedHashedFeatureSet>());
            assert!(!features.is_instance_of::<RoleTaggedHashedFeatureSet>());
            assert_eq!(
                features
                    .getattr("id_width")
                    .unwrap()
                    .extract::<u16>()
                    .unwrap(),
                64
            );
            assert_eq!(
                entries.extract::<Vec<(u128, i32)>>().unwrap(),
                expected_entries
            );
            entries
                .cast::<PyList>()
                .unwrap()
                .append((9u128, 3i32))
                .unwrap();
            assert_eq!(
                fingerprint
                    .getattr("features")
                    .unwrap()
                    .getattr("entries")
                    .unwrap()
                    .extract::<Vec<(u128, i32)>>()
                    .unwrap(),
                expected_entries
            );
        });
    }

    #[rstest]
    #[case::morgan(
        ReactionCombinedFingerprintConfig::DisjointUnion {
            molecule: HashedFingerprintConfig::Morgan {
                radius: 2,
                ring_config: RingConfig::default(),
            },
        },
        vec![
            (ReactionSide::Reactant, 864662311),
            (ReactionSide::Reactant, 1535166686),
            (ReactionSide::Reactant, 2245384272),
            (ReactionSide::Reactant, 2246728737),
            (ReactionSide::Reactant, 3542456614),
            (ReactionSide::Reactant, 4018048386),
            (ReactionSide::Product, 2246728737),
            (ReactionSide::Product, 2246997334),
            (ReactionSide::Product, 3548082732),
        ]
    )]
    #[case::ecfp(
        ReactionCombinedFingerprintConfig::DisjointUnion {
            molecule: HashedFingerprintConfig::Ecfp {
                radius: 2,
                hashing_scheme: EcfpHashScheme::Xxh3Width64V1(),
                ring_config: RingConfig::default(),
            },
        },
        vec![
            (ReactionSide::Reactant, 63839236075656913),
            (ReactionSide::Reactant, 1189585227353469813),
            (ReactionSide::Reactant, 3822471596818936039),
            (ReactionSide::Reactant, 13652293261850732425),
            (ReactionSide::Reactant, 15001976065402722634),
            (ReactionSide::Reactant, 16149328945726899460),
            (ReactionSide::Product, 896060437578512973),
            (ReactionSide::Product, 13327007941213506523),
            (ReactionSide::Product, 16149328945726899460),
        ]
    )]
    #[case::wl(
        ReactionCombinedFingerprintConfig::DisjointUnion {
            molecule: HashedFingerprintConfig::Wl {
                rounds: RefinementRounds::Fixed { rounds: 3 },
                hashing_scheme: WlHashScheme::Xxh3SortedWidth64V1(),
            },
        },
        vec![
            (ReactionSide::Reactant, 2520347590860685079),
            (ReactionSide::Reactant, 3352603313223549703),
            (ReactionSide::Reactant, 4152249898001161146),
            (ReactionSide::Reactant, 5715207763479934940),
            (ReactionSide::Reactant, 5807737097854608645),
            (ReactionSide::Reactant, 7542810387455301591),
            (ReactionSide::Reactant, 11457795998246593156),
            (ReactionSide::Reactant, 11986000156817227245),
            (ReactionSide::Reactant, 12895020514073294021),
            (ReactionSide::Reactant, 13932567567828606490),
            (ReactionSide::Reactant, 17305796300852423160),
            (ReactionSide::Reactant, 17417400371411086222),
            (ReactionSide::Product, 5715207763479934940),
            (ReactionSide::Product, 6786829771653353480),
            (ReactionSide::Product, 7404535559284410087),
            (ReactionSide::Product, 7542810387455301591),
            (ReactionSide::Product, 8754482138526219790),
            (ReactionSide::Product, 11457795998246593156),
            (ReactionSide::Product, 12849090138728295812),
            (ReactionSide::Product, 16456488943967932267),
        ]
    )]
    fn test_reaction_combined_fingerprint_disjoint_union(
        ethanol_deoxygenation: GraphIrReaction,
        #[case] config: ReactionCombinedFingerprintConfig,
        #[case] expected_ids: Vec<(ReactionSide, u128)>,
    ) {
        Python::attach(|py| {
            let reaction = Reaction::from_rust(py, ethanol_deoxygenation).unwrap();
            let fingerprint = reaction.combined_fingerprint(py, config).unwrap();
            let fingerprint = Py::new(py, fingerprint).unwrap();
            let fingerprint = fingerprint.bind(py).as_any();
            let features = fingerprint.getattr("features").unwrap();
            let ids = features.getattr("ids").unwrap();

            assert!(features.is_instance_of::<RoleTaggedHashedFeatureSet>());
            assert!(!features.is_instance_of::<SignedHashedFeatureSet>());
            assert_eq!(
                features
                    .getattr("id_width")
                    .unwrap()
                    .extract::<u16>()
                    .unwrap(),
                64
            );
            assert_eq!(
                ids.extract::<Vec<(ReactionSide, u128)>>().unwrap(),
                expected_ids
            );
            ids.cast::<PyList>()
                .unwrap()
                .append((ReactionSide::Product, 9u128))
                .unwrap();
            assert_eq!(
                fingerprint
                    .getattr("features")
                    .unwrap()
                    .getattr("ids")
                    .unwrap()
                    .extract::<Vec<(ReactionSide, u128)>>()
                    .unwrap(),
                expected_ids
            );
        });
    }

    #[rstest]
    #[case::difference(
        ReactionCombinedFingerprintConfig::Difference {
            molecule: HashedFingerprintConfig::Morgan {
                radius: 2,
                ring_config: RingConfig::default(),
            },
        }
    )]
    fn test_reaction_combined_fingerprint_difference_identity(
        ethanol_identity: GraphIrReaction,
        #[case] config: ReactionCombinedFingerprintConfig,
    ) {
        Python::attach(|py| {
            let reaction = Reaction::from_rust(py, ethanol_identity).unwrap();
            let fingerprint = reaction.combined_fingerprint(py, config).unwrap();
            let fingerprint = Py::new(py, fingerprint).unwrap();
            let features = fingerprint.bind(py).getattr("features").unwrap();

            assert!(features.is_instance_of::<SignedHashedFeatureSet>());
            assert_eq!(
                features
                    .getattr("entries")
                    .unwrap()
                    .extract::<Vec<(u128, i32)>>()
                    .unwrap(),
                Vec::new()
            );
            assert_eq!(
                features
                    .getattr("id_width")
                    .unwrap()
                    .extract::<u16>()
                    .unwrap(),
                64
            );
        });
    }

    #[rstest]
    #[case::disjoint_union(
        ReactionCombinedFingerprintConfig::DisjointUnion {
            molecule: HashedFingerprintConfig::Morgan {
                radius: 2,
                ring_config: RingConfig::default(),
            },
        },
        vec![
            (ReactionSide::Reactant, 864662311),
            (ReactionSide::Reactant, 1535166686),
            (ReactionSide::Reactant, 2245384272),
            (ReactionSide::Reactant, 2246728737),
            (ReactionSide::Reactant, 3542456614),
            (ReactionSide::Reactant, 4018048386),
            (ReactionSide::Product, 864662311),
            (ReactionSide::Product, 1535166686),
            (ReactionSide::Product, 2245384272),
            (ReactionSide::Product, 2246728737),
            (ReactionSide::Product, 3542456614),
            (ReactionSide::Product, 4018048386),
        ]
    )]
    fn test_reaction_combined_fingerprint_disjoint_union_identity(
        ethanol_identity: GraphIrReaction,
        #[case] config: ReactionCombinedFingerprintConfig,
        #[case] expected_ids: Vec<(ReactionSide, u128)>,
    ) {
        Python::attach(|py| {
            let reaction = Reaction::from_rust(py, ethanol_identity).unwrap();
            let fingerprint = reaction.combined_fingerprint(py, config).unwrap();
            let fingerprint = Py::new(py, fingerprint).unwrap();
            let features = fingerprint.bind(py).getattr("features").unwrap();

            assert!(features.is_instance_of::<RoleTaggedHashedFeatureSet>());
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
                64
            );
        });
    }

    #[rstest]
    #[case::reactant_not_ground(
        GraphIrReaction::new(
            mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
            GraphIrDeltas::new(),
        ),
        ReactionCombinedFingerprintConfig::Difference {
            molecule: HashedFingerprintConfig::Morgan {
                radius: 2,
                ring_config: RingConfig::default(),
            },
        },
        "UnderdeterminedError",
        "fingerprint requires a determined molecule",
    )]
    #[case::product_not_ground(
        GraphIrReaction::new(
            mol_dsl_ground!(r#"{:atoms ["C #h4"] :bonds []}"#),
            GraphIrDeltas::from_iter([GraphIrDelta::Atom(GraphIrAtomDelta::ModifyField {
                id: GraphIrAtomId(0),
                change: GraphIrAtomFieldChange::Charge {
                    old: GraphIrNumForm::Lit(0),
                    new: GraphIrNumForm::Undetermined,
                },
            })]),
        ),
        ReactionCombinedFingerprintConfig::DisjointUnion {
            molecule: HashedFingerprintConfig::Morgan {
                radius: 2,
                ring_config: RingConfig::default(),
            },
        },
        "UnderdeterminedError",
        "fingerprint requires a determined molecule",
    )]
    #[case::inconsistent(
        GraphIrReaction::new(
            mol_dsl_ground!(r#"{:atoms ["C #h4"] :bonds []}"#),
            GraphIrDeltas::from_iter([GraphIrDelta::Atom(GraphIrAtomDelta::ModifyField {
                id: GraphIrAtomId(0),
                change: GraphIrAtomFieldChange::Charge {
                    old: GraphIrNumForm::Lit(1),
                    new: GraphIrNumForm::Lit(0),
                },
            })]),
        ),
        ReactionCombinedFingerprintConfig::Difference {
            molecule: HashedFingerprintConfig::Morgan {
                radius: 2,
                ring_config: RingConfig::default(),
            },
        },
        "ContradictionError",
        "reaction fingerprint input is inconsistent",
    )]
    fn test_reaction_combined_fingerprint_error(
        #[case] input: GraphIrReaction,
        #[case] config: ReactionCombinedFingerprintConfig,
        #[case] expected_type: &str,
        #[case] expected_message: &str,
    ) {
        Python::attach(|py| {
            let reaction = Reaction::from_rust(py, input).unwrap();
            let error = reaction.combined_fingerprint(py, config).unwrap_err();

            assert_eq!(error.get_type(py).name().unwrap(), expected_type);
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                expected_message
            );
        });
    }

    #[rstest]
    fn test_reaction_eq() {
        Python::attach(|py| {
            let empty = Reaction::new(py, None, None).unwrap();
            let other_empty = Reaction::new(py, None, None).unwrap();
            let populated = Reaction::from_rust(
                py,
                GraphIrReaction::new(
                    GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                        atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                        ..Default::default()
                    }),
                    GraphIrDeltas::new(),
                ),
            )
            .unwrap();

            assert!(empty.__eq__(&other_empty, py));
            assert!(!empty.__eq__(&populated, py));
            let empty = Py::new(py, empty).unwrap();
            assert!(empty
                .bind(py)
                .hash()
                .unwrap_err()
                .is_instance_of::<PyTypeError>(py));
        });
    }

    #[rstest]
    #[case::empty(
        GraphIrReaction::default(),
        r##"{:deltas [] :lhs {:atoms [] :bonds []}}"##
    )]
    #[case::populated(
        GraphIrReaction::new(
            GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                ..Default::default()
            }),
            vec![GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                id: GraphIrAtomId(1),
                attributes: GraphIrAtomForm::from_element(ChemElement::O),
            })].into_iter().collect(),
        ),
        r##"{:deltas [{:atom {:add "O"}}] :lhs {:atoms ["C"] :bonds []}}"##,
    )]
    fn test_reaction_str(#[case] input: GraphIrReaction, #[case] expected: &str) {
        Python::attach(|py| {
            let reaction = Reaction::from_rust(py, input).unwrap();

            assert_eq!(reaction.__str__(py), expected);
            assert_eq!(reaction.__str__(py), reaction.render(py, None));
        });
    }

    #[rstest]
    fn test_reaction_str_components() {
        Python::attach(|py| {
            let reaction = Reaction::from_rust(
                py,
                GraphIrReaction::new(
                    GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                        atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                        ..Default::default()
                    }),
                    GraphIrDeltas::new(),
                ),
            )
            .unwrap();

            *reaction.lhs.bind(py).borrow_mut().to_rust_mut() =
                GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![GraphIrAtomForm::from_element(ChemElement::C).with_charge(1)],
                    ..Default::default()
                });
            let delta = into_py_variant(
                py,
                Delta::from_rust(
                    py,
                    &GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                        id: GraphIrAtomId(1),
                        attributes: GraphIrAtomForm::from_element(ChemElement::O),
                    }),
                )
                .unwrap(),
            )
            .unwrap();
            reaction
                .deltas
                .bind(py)
                .call_method1("append", (delta,))
                .unwrap();

            assert_eq!(
                reaction.__str__(py),
                r##"{:deltas [{:atom {:add "O"}}] :lhs {:atoms ["C#c+"] :bonds []}}"##
            );
        });
    }

    #[rstest]
    #[case::atom_add_remove(
        r##"{:lhs {:atoms ["C" "O"]} :deltas [{:atom {:add "N"}} {:atom {:remove 1}}]}"##
    )]
    #[case::atom_modify(r##"{:lhs {:atoms ["Br#c0"]} :deltas [{:atom {:modify [0 "#c-1"]}}]}"##)]
    #[case::stereo_modify(
        r##"{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]} :deltas [{:stereo-atom {:modify [0 "Th0"]}}]}"##
    )]
    #[case::molecule_constraint(
        r##"{:lhs {:atoms ["C"]} :deltas [{:constraint {:add {:connected {}}}}]}"##
    )]
    fn test_reaction_str_roundtrip(#[case] text: &str) {
        Python::attach(|py| {
            let first = Reaction::parse(py, text, None).unwrap();

            let canonical = first.__str__(py);
            let second = Reaction::parse(py, &canonical, None).unwrap();

            assert!(first.__eq__(&second, py));
            assert_eq!(second.__str__(py), canonical);
        });
    }

    #[rstest]
    fn test_reaction_repr() {
        Python::attach(|py| {
            let reaction = Reaction::from_rust(
                py,
                GraphIrReaction::new(
                    GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                        atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                        ..Default::default()
                    }),
                    vec![GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                        id: GraphIrAtomId(1),
                        attributes: GraphIrAtomForm::from_element(ChemElement::O),
                    })]
                    .into_iter()
                    .collect(),
                ),
            )
            .unwrap();

            assert_eq!(
                reaction.__repr__(py).unwrap(),
                "Reaction(lhs=Molecule(atoms=1, bonds=0), deltas=Deltas([Delta.Atom(AtomDelta.Add(id=1, attributes=AtomForm.parse('O')))]))"
            );
        });
    }

    #[rstest]
    #[case::empty(GraphIrReaction::default())]
    #[case::populated(GraphIrReaction::new(
        GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
            atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
            ..Default::default()
        }),
        vec![GraphIrDelta::Atom(GraphIrAtomDelta::Add {
            id: GraphIrAtomId(1),
            attributes: GraphIrAtomForm::from_element(ChemElement::O),
        })]
        .into_iter()
        .collect(),
    ))]
    fn test_reaction_from_rust(#[case] expected: GraphIrReaction) {
        Python::attach(|py| {
            let reaction = Reaction::from_rust(py, expected.clone()).unwrap();

            assert_eq!(reaction.to_rust(py), expected);
        });
    }

    #[rstest]
    fn test_reaction_to_rust() {
        Python::attach(|py| {
            let expected = GraphIrReaction::new(
                GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                    ..Default::default()
                }),
                vec![GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                    id: GraphIrAtomId(1),
                    attributes: GraphIrAtomForm::from_element(ChemElement::O),
                })]
                .into_iter()
                .collect(),
            );
            let reaction = Reaction::from_rust(py, expected.clone()).unwrap();

            let mut snapshot = reaction.to_rust(py);
            snapshot.lhs = GraphIrMolecule::new();
            snapshot.deltas = GraphIrDeltas::new();

            assert_eq!(reaction.to_rust(py), expected);
        });
    }

    #[rstest]
    fn test_reaction_to_rust_roundtrip() {
        Python::attach(|py| {
            let expected = GraphIrReaction::new(
                GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                    ..Default::default()
                }),
                vec![GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                    id: GraphIrAtomId(1),
                    attributes: GraphIrAtomForm::from_element(ChemElement::O),
                })]
                .into_iter()
                .collect(),
            );
            let python = Py::new(py, Reaction::from_rust(py, expected.clone()).unwrap()).unwrap();

            let rust = python.bind(py).borrow().to_rust(py);
            let roundtrip = Py::new(py, Reaction::from_rust(py, rust).unwrap()).unwrap();

            assert_eq!(roundtrip.bind(py).borrow().to_rust(py), expected);
            assert_ne!(python.as_ptr(), roundtrip.as_ptr());
            assert_ne!(
                python.bind(py).borrow().lhs.as_ptr(),
                roundtrip.bind(py).borrow().lhs.as_ptr()
            );
            assert_ne!(
                python.bind(py).borrow().deltas.as_ptr(),
                roundtrip.bind(py).borrow().deltas.as_ptr()
            );
        });
    }

    #[fixture]
    fn derivation_and_host() -> (GraphIrReactionDerivation, GraphIrMolecule) {
        let pattern = GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
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
        let host = pattern.clone();
        let reaction = GraphIrReaction::new(
            pattern.clone(),
            GraphIrDeltas::from_iter([GraphIrDelta::Bond(GraphIrBondDelta::ModifyField {
                id: GraphIrBondId(0),
                change: GraphIrBondFieldChange::Order {
                    old: GraphIrNumForm::Lit(1),
                    new: GraphIrNumForm::Lit(2),
                },
            })]),
        );
        let correspondence = GraphIrMoleculeCorrespondence::induce(
            &pattern,
            &host,
            Correspondence::new(
                vec![
                    (GraphIrAtomId(0), GraphIrAtomId(0)),
                    (GraphIrAtomId(1), GraphIrAtomId(1)),
                ],
                2,
                2,
            )
            .expect("correspondence producer preserves partial-bijection invariants"),
        )
        .expect("the atom correspondence describes the molecule pair");
        let derivation = reaction.apply_at(&host, &correspondence).unwrap();
        (derivation, host)
    }

    #[rstest]
    fn test_reaction_derivation_observations(
        derivation_and_host: (GraphIrReactionDerivation, GraphIrMolecule),
    ) {
        let (expected, mut host) = derivation_and_host;
        let derivation = ReactionDerivation::from_rust(expected.clone());

        assert_eq!(derivation.lhs().to_rust(), expected.lhs());
        assert_eq!(derivation.rhs().to_rust(), expected.rhs());
        assert_eq!(
            derivation.comap(),
            PyMoleculeCorrespondence::from_rust(expected.comap().clone())
        );
        assert_eq!(
            derivation.atom_correspondence(),
            PyCorrespondence::from_rust(expected.atom_correspondence())
        );

        *host.atom_mut(GraphIrAtomId(0)).attributes = GraphIrAtomForm::from_element(ChemElement::F);
        let mut lhs = derivation.lhs();
        *lhs.to_rust_mut().atom_mut(GraphIrAtomId(0)).attributes =
            GraphIrAtomForm::from_element(ChemElement::N);

        assert_eq!(derivation.to_rust(), &expected);
        assert_ne!(derivation.lhs().to_rust(), &host);
        assert_ne!(derivation.lhs().to_rust(), lhs.to_rust());
    }

    #[rstest]
    fn test_reaction_derivation_reverse(
        derivation_and_host: (GraphIrReactionDerivation, GraphIrMolecule),
    ) {
        let (expected, _) = derivation_and_host;
        let derivation = ReactionDerivation::from_rust(expected.clone());
        let reversed = derivation.reverse();
        let mut reversed_lhs = reversed.lhs();
        *reversed_lhs
            .to_rust_mut()
            .atom_mut(GraphIrAtomId(0))
            .attributes = GraphIrAtomForm::from_element(ChemElement::N);

        assert_eq!(reversed.to_rust(), &expected.reverse());
        assert_eq!(derivation.to_rust(), &expected);
        assert_ne!(reversed.lhs().to_rust(), reversed_lhs.to_rust());
    }

    #[rstest]
    fn test_reaction_derivation_chain(
        derivation_and_host: (GraphIrReactionDerivation, GraphIrMolecule),
    ) {
        let (first, _) = derivation_and_host;
        let middle = first.rhs().clone();
        let reaction = GraphIrReaction::new(
            middle.clone(),
            GraphIrDeltas::from_iter([GraphIrDelta::Bond(GraphIrBondDelta::ModifyField {
                id: GraphIrBondId(0),
                change: GraphIrBondFieldChange::Order {
                    old: GraphIrNumForm::Lit(2),
                    new: GraphIrNumForm::Lit(3),
                },
            })]),
        );
        let correspondence = GraphIrMoleculeCorrespondence::induce(
            &middle,
            &middle,
            Correspondence::new(
                vec![
                    (GraphIrAtomId(0), GraphIrAtomId(0)),
                    (GraphIrAtomId(1), GraphIrAtomId(1)),
                ],
                2,
                2,
            )
            .expect("correspondence producer preserves partial-bijection invariants"),
        )
        .expect("the atom correspondence describes the molecule pair");
        let second = reaction.apply_at(&middle, &correspondence).unwrap();
        let first_value = ReactionDerivation::from_rust(first.clone());
        let second_value = ReactionDerivation::from_rust(second.clone());
        let chained = first_value.chain(&second_value);
        let mut chained_rhs = chained.rhs();
        *chained_rhs
            .to_rust_mut()
            .atom_mut(GraphIrAtomId(0))
            .attributes = GraphIrAtomForm::from_element(ChemElement::N);

        assert_eq!(chained.to_rust(), &first.chain(&second));
        assert_eq!(first_value.to_rust(), &first);
        assert_eq!(second_value.to_rust(), &second);
        assert_ne!(chained.rhs().to_rust(), chained_rhs.to_rust());
    }

    #[rstest]
    fn test_reaction_derivation_to_reaction(
        derivation_and_host: (GraphIrReactionDerivation, GraphIrMolecule),
    ) {
        let (expected_derivation, _) = derivation_and_host;
        let expected_reaction = GraphIrReaction::new(
            GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
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
            }),
            GraphIrDeltas::from_iter([GraphIrDelta::Bond(GraphIrBondDelta::ModifyField {
                id: GraphIrBondId(0),
                change: GraphIrBondFieldChange::Order {
                    old: GraphIrNumForm::Lit(1),
                    new: GraphIrNumForm::Lit(2),
                },
            })]),
        );
        let derivation = ReactionDerivation::from_rust(expected_derivation.clone());

        Python::attach(|py| {
            let first = derivation.to_reaction(py).unwrap();
            let second = derivation.to_reaction(py).unwrap();

            assert_eq!(first.to_rust(py), expected_reaction);
            assert_eq!(second.to_rust(py), expected_reaction);
            assert_ne!(first.lhs.as_ptr(), second.lhs.as_ptr());
            assert_ne!(first.deltas.as_ptr(), second.deltas.as_ptr());

            *first.lhs.bind(py).borrow_mut().to_rust_mut() = GraphIrMolecule::new();
            let delta = into_py_variant(
                py,
                Delta::from_rust(
                    py,
                    &GraphIrDelta::Atom(GraphIrAtomDelta::Add {
                        id: GraphIrAtomId(2),
                        attributes: GraphIrAtomForm::from_element(ChemElement::O),
                    }),
                )
                .unwrap(),
            )
            .unwrap();
            first
                .deltas
                .bind(py)
                .call_method1("append", (delta,))
                .unwrap();

            assert_eq!(second.to_rust(py), expected_reaction);
            assert_eq!(derivation.to_rust(), &expected_derivation);
        });
    }

    #[rstest]
    fn test_reaction_derivation_value(
        derivation_and_host: (GraphIrReactionDerivation, GraphIrMolecule),
    ) {
        let (expected, _) = derivation_and_host;
        Python::attach(|py| {
            let derivation = Py::new(py, ReactionDerivation::from_rust(expected.clone())).unwrap();
            let equal = Py::new(py, ReactionDerivation::from_rust(expected.clone())).unwrap();
            let unequal = Py::new(py, ReactionDerivation::from_rust(expected.reverse())).unwrap();
            let first_lhs = derivation.bind(py).getattr("lhs").unwrap();
            let second_lhs = derivation.bind(py).getattr("lhs").unwrap();
            let first_comap = derivation.bind(py).getattr("comap").unwrap();
            let second_comap = derivation.bind(py).getattr("comap").unwrap();

            assert!(derivation
                .bind(py)
                .as_any()
                .eq(equal.bind(py).as_any())
                .unwrap());
            assert!(!derivation
                .bind(py)
                .as_any()
                .eq(unequal.bind(py).as_any())
                .unwrap());
            assert!(!first_lhs.is(&second_lhs));
            assert!(!first_comap.is(&second_comap));
            assert_eq!(
                derivation
                    .bind(py)
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                concat!(
                    "ReactionDerivation(lhs=Molecule(atoms=2, bonds=1), ",
                    "rhs=Molecule(atoms=2, bonds=1), ",
                    "comap=MoleculeCorrespondence(",
                    "atoms=Correspondence(matched_pairs=[(0, 0), (1, 1)], left_count=2, right_count=2), ",
                    "bonds=Correspondence(matched_pairs=[(0, 0)], left_count=1, right_count=1), ",
                    "dative_bonds=Correspondence(matched_pairs=[], left_count=0, right_count=0), ",
                    "aromatic_systems=Correspondence(matched_pairs=[], left_count=0, right_count=0), ",
                    "multicenter_bonds=Correspondence(matched_pairs=[], left_count=0, right_count=0), ",
                    "noncovalent_bonds=Correspondence(matched_pairs=[], left_count=0, right_count=0), ",
                    "stereo_atoms=Correspondence(matched_pairs=[], left_count=0, right_count=0), ",
                    "stereo_bonds=Correspondence(matched_pairs=[], left_count=0, right_count=0)))"
                )
            );
        });
    }

    #[rstest]
    fn test_reaction_derivation_roundtrip(
        derivation_and_host: (GraphIrReactionDerivation, GraphIrMolecule),
    ) {
        let (expected, _) = derivation_and_host;
        assert_eq!(
            ReactionDerivation::from_rust(expected.clone()).to_rust(),
            &expected
        );
    }

    #[fixture]
    fn reaction_application() -> (GraphIrReaction, GraphIrMolecule) {
        let reaction = GraphIrReaction::new(
            GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                ..Default::default()
            }),
            GraphIrDeltas::from_iter([GraphIrDelta::Atom(GraphIrAtomDelta::ModifyField {
                id: GraphIrAtomId(0),
                change: GraphIrAtomFieldChange::Charge {
                    old: GraphIrNumForm::Undetermined,
                    new: GraphIrNumForm::Lit(1),
                },
            })]),
        );
        let host = GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
            atoms: vec![
                GraphIrAtomForm::from_element(ChemElement::C),
                GraphIrAtomForm::from_element(ChemElement::C),
            ],
            ..Default::default()
        });
        (reaction, host)
    }

    #[rstest]
    fn test_reaction_application_iter_identity(
        reaction_application: (GraphIrReaction, GraphIrMolecule),
    ) {
        let (reaction, host) = reaction_application;
        Python::attach(|py| {
            let application = Py::new(
                py,
                ReactionApplicationIter::from_rust(
                    reaction
                        .apply(&host, ReactionApplicationConfig::default().to_rust())
                        .unwrap(),
                ),
            )
            .unwrap();

            let iter = application.bind(py).call_method0("__iter__").unwrap();
            assert!(iter.is(application.bind(py)));
        });
    }

    #[rstest]
    fn test_reaction_application_iter(reaction_application: (GraphIrReaction, GraphIrMolecule)) {
        let (reaction, host) = reaction_application;
        let mut application = ReactionApplicationIter::from_rust(
            reaction
                .apply(&host, ReactionApplicationConfig::default().to_rust())
                .unwrap(),
        );

        let first = application.__next__().unwrap().unwrap();
        let second = application.__next__().unwrap().unwrap();
        assert_eq!(
            [
                first.rhs().to_rust().clone(),
                second.rhs().to_rust().clone()
            ],
            [
                GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![
                        GraphIrAtomForm::from_element(ChemElement::C).with_charge(1),
                        GraphIrAtomForm::from_element(ChemElement::C),
                    ],
                    ..Default::default()
                }),
                GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![
                        GraphIrAtomForm::from_element(ChemElement::C),
                        GraphIrAtomForm::from_element(ChemElement::C).with_charge(1),
                    ],
                    ..Default::default()
                }),
            ]
        );
        assert_eq!(application.__next__().unwrap(), None);
        assert_eq!(application.__next__().unwrap(), None);

        let expected_first = first.to_rust();
        let expected_second = second.to_rust();
        let mut detached = first.rhs();
        *detached.to_rust_mut().atom_mut(GraphIrAtomId(0)).attributes =
            GraphIrAtomForm::from_element(ChemElement::F);
        assert_eq!(first.to_rust(), expected_first);
        assert_eq!(second.to_rust(), expected_second);
    }

    #[rstest]
    fn test_reaction_application_iter_empty() {
        let reaction = GraphIrReaction::new(
            GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                atoms: vec![GraphIrAtomForm::from_element(ChemElement::N)],
                ..Default::default()
            }),
            GraphIrDeltas::new(),
        );
        let host = GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
            atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
            ..Default::default()
        });
        let mut application = ReactionApplicationIter::from_rust(
            reaction
                .apply(&host, ReactionApplicationConfig::default().to_rust())
                .unwrap(),
        );

        assert_eq!(application.__next__().unwrap(), None);
        assert_eq!(application.__next__().unwrap(), None);
    }

    #[rstest]
    #[ignore = "re-enable when matching evaluates molecule-scope pattern constraints"]
    fn test_reaction_application_iter_error() {
        let constraint = GraphIrConstraint::Molecule(GraphIrMoleculeConstraint::ChargeSum {
            atoms: Some(vec![GraphIrAtomId(0)]),
            sum: GraphIrNumForm::Lit(0),
        });
        let reaction = GraphIrReaction::new(
            GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                constraints: constraint.clone().into(),
                ..Default::default()
            }),
            GraphIrDeltas::from_iter([GraphIrDelta::Constraint(GraphIrConstraintDelta::Remove(
                constraint,
            ))]),
        );
        let host = GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
            atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
            ..Default::default()
        });
        let mut application = ReactionApplicationIter::from_rust(
            reaction
                .apply(&host, ReactionApplicationConfig::default().to_rust())
                .unwrap(),
        );

        let error = application.__next__().unwrap_err();

        Python::attach(|py| {
            assert!(error.is_instance_of::<TransactionError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "missing constraint entry on remove"
            );
        });
        assert_eq!(application.__next__().unwrap(), None);
        assert_eq!(application.__next__().unwrap(), None);
    }

    #[rstest]
    fn test_reaction_products_iter_identity() {
        let reaction = GraphIrReaction::default();
        let host = GraphIrMolecule::default();
        Python::attach(|py| {
            let products = Py::new(
                py,
                ReactionProductsIter::from_rust(
                    host.react(&reaction, ReactionApplicationConfig::default().to_rust())
                        .unwrap(),
                ),
            )
            .unwrap();

            let iter = products.bind(py).call_method0("__iter__").unwrap();
            assert!(iter.is(products.bind(py)));
        });
    }

    #[rstest]
    fn test_reaction_products_iter_python() {
        let reaction = GraphIrReaction::default();
        let host = GraphIrMolecule::default();
        Python::attach(|py| {
            let products = Py::new(
                py,
                ReactionProductsIter::from_rust(
                    host.react(&reaction, ReactionApplicationConfig::default().to_rust())
                        .unwrap(),
                ),
            )
            .unwrap();

            let first = products.bind(py).call_method0("__next__").unwrap();
            assert!(first.cast::<PyList>().unwrap().is_empty());
            assert!(products
                .bind(py)
                .call_method0("__next__")
                .unwrap_err()
                .is_instance_of::<PyStopIteration>(py));
            assert!(products
                .bind(py)
                .call_method0("__next__")
                .unwrap_err()
                .is_instance_of::<PyStopIteration>(py));
        });
    }

    #[rstest]
    fn test_reaction_products_iter() {
        let reaction = GraphIrReaction::new(
            GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                ..Default::default()
            }),
            GraphIrDeltas::from_iter([GraphIrDelta::Atom(GraphIrAtomDelta::ModifyField {
                id: GraphIrAtomId(0),
                change: GraphIrAtomFieldChange::Charge {
                    old: GraphIrNumForm::Undetermined,
                    new: GraphIrNumForm::Lit(1),
                },
            })]),
        );
        let host = GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
            atoms: vec![
                GraphIrAtomForm::from_element(ChemElement::C),
                GraphIrAtomForm::from_element(ChemElement::C),
                GraphIrAtomForm::from_element(ChemElement::N),
            ],
            ..Default::default()
        });
        let expected_host = host.clone();
        let mut products = ReactionProductsIter::from_rust(
            host.react(&reaction, ReactionApplicationConfig::default().to_rust())
                .unwrap(),
        );

        let mut first = products.__next__().unwrap().unwrap();
        let second = products.__next__().unwrap().unwrap();

        assert_eq!(
            first
                .iter()
                .map(|molecule| molecule.to_rust().clone())
                .collect::<Vec<_>>(),
            vec![
                GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![GraphIrAtomForm::from_element(ChemElement::C).with_charge(1)],
                    ..Default::default()
                }),
                GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                    ..Default::default()
                }),
                GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![GraphIrAtomForm::from_element(ChemElement::N)],
                    ..Default::default()
                }),
            ]
        );
        assert_eq!(
            second
                .iter()
                .map(|molecule| molecule.to_rust().clone())
                .collect::<Vec<_>>(),
            vec![
                GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                    ..Default::default()
                }),
                GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![GraphIrAtomForm::from_element(ChemElement::C).with_charge(1)],
                    ..Default::default()
                }),
                GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![GraphIrAtomForm::from_element(ChemElement::N)],
                    ..Default::default()
                }),
            ]
        );
        assert_eq!(products.__next__().unwrap(), None);
        assert_eq!(products.__next__().unwrap(), None);

        *first[0].to_rust_mut().atom_mut(GraphIrAtomId(0)).attributes =
            GraphIrAtomForm::from_element(ChemElement::F);
        assert_eq!(host, expected_host);
        assert_eq!(
            second[0].to_rust(),
            &GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                ..Default::default()
            })
        );
    }

    #[rstest]
    fn test_reaction_products_iter_empty() {
        let reaction = GraphIrReaction::new(
            GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                atoms: vec![GraphIrAtomForm::from_element(ChemElement::N)],
                ..Default::default()
            }),
            GraphIrDeltas::new(),
        );
        let host = GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
            atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
            ..Default::default()
        });
        let mut products = ReactionProductsIter::from_rust(
            host.react(&reaction, ReactionApplicationConfig::default().to_rust())
                .unwrap(),
        );

        assert_eq!(products.__next__().unwrap(), None);
        assert_eq!(products.__next__().unwrap(), None);
    }

    #[rstest]
    #[ignore = "re-enable when matching evaluates molecule-scope pattern constraints"]
    fn test_reaction_products_iter_error() {
        let constraint = GraphIrConstraint::Molecule(GraphIrMoleculeConstraint::ChargeSum {
            atoms: Some(vec![GraphIrAtomId(0)]),
            sum: GraphIrNumForm::Lit(0),
        });
        let reaction = GraphIrReaction::new(
            GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                constraints: constraint.clone().into(),
                ..Default::default()
            }),
            GraphIrDeltas::from_iter([GraphIrDelta::Constraint(GraphIrConstraintDelta::Remove(
                constraint,
            ))]),
        );
        let host = GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
            atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
            ..Default::default()
        });
        let mut products = ReactionProductsIter::from_rust(
            host.react(&reaction, ReactionApplicationConfig::default().to_rust())
                .unwrap(),
        );

        let error = products.__next__().unwrap_err();

        Python::attach(|py| {
            assert!(error.is_instance_of::<TransactionError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "missing constraint entry on remove"
            );
        });
        assert_eq!(products.__next__().unwrap(), None);
        assert_eq!(products.__next__().unwrap(), None);
    }
}
