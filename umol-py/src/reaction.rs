//! `ReactionAst` — an owned Python component facade over the Rust reaction AST.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use std::str::FromStr;
use std::vec::IntoIter;

use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
#[cfg(test)]
use umol_ast::ast::SubstructureMatchAlgorithm as AstSubstructureMatchAlgorithm;
use umol_ast::ast::{
    ApplyError as AstApplyError, AtomId, FromAst, IntoAst, MoleculeAst as AstMoleculeAst,
    MoleculeCorrespondence as AstMoleculeCorrespondence, ReactionAst as AstReactionAst,
    ReactionDerivation as AstReactionDerivation,
    SubstructureMatchConfig as AstSubstructureMatchConfig,
};
use umol_ast::dsl::ReactionDsl as AstReactionDsl;
use umol_graph::fingerprint::featurize_reaction;
use umol_graph::ingest::ingest_reaction_smiles_with;
use umol_graph::ops::model::ChemistryModel as GraphChemistryModel;
use umol_graph::ops::resolve::ResolveConfig as GraphResolveConfig;
use umol_graph_core::{
    CommonSubgraphEnumerationAlgorithm as GraphCoreCommonSubgraphEnumerationAlgorithm,
    Correspondence, CorrespondenceError,
};
#[cfg(test)]
use umol_graph_core::{
    RelevantCycleEnumerationAlgorithm as GraphCoreRelevantCycleEnumerationAlgorithm,
    SubgraphIsomorphismAlgorithm as GraphCoreSubgraphIsomorphismAlgorithm,
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
use crate::lattice::impl_py_canonicalize;
use crate::metadata::ReactionMetadata;
use crate::model::ChemistryModel;
use crate::molecule::MoleculeAst;
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

/// Validate atom pairs and construct their partial bijection over the two side sizes.
fn atom_correspondence(
    pairs: Vec<(usize, usize)>,
    lhs_count: usize,
    rhs_count: usize,
) -> PyResult<Correspondence<AtomId>> {
    let matched_pairs = pairs
        .into_iter()
        .map(|(left, right)| (AtomId::from(left), AtomId::from(right)))
        .collect();
    Correspondence::new(matched_pairs, lhs_count, rhs_count).map_err(|error| {
        PyValueError::new_err(match error {
            CorrespondenceError::LeftIdOutOfRange { id, count } => {
                format!("left atom id {} out of range for {count} atoms", id.index())
            }
            CorrespondenceError::RightIdOutOfRange { id, count } => {
                format!(
                    "right atom id {} out of range for {count} atoms",
                    id.index()
                )
            }
            CorrespondenceError::DuplicateLeftId { id } => {
                format!("duplicate left atom id {}", id.index())
            }
            CorrespondenceError::DuplicateRightId { id } => {
                format!("duplicate right atom id {}", id.index())
            }
            CorrespondenceError::LeftCountMismatch { declared, actual } => {
                format!("declared left atom count {declared} does not match actual count {actual}")
            }
            CorrespondenceError::RightCountMismatch { declared, actual } => {
                format!("declared right atom count {declared} does not match actual count {actual}")
            }
        })
    })
}

/// A reaction whose molecule and delta components remain live Python values.
#[pyclass]
pub struct ReactionAst {
    lhs: Py<MoleculeAst>,
    deltas: Py<Deltas>,
}

#[pymethods]
impl ReactionAst {
    /// Build a reaction from detached component snapshots.
    #[new]
    #[pyo3(signature = (lhs=None, deltas=None))]
    fn new(
        py: Python<'_>,
        lhs: Option<Py<MoleculeAst>>,
        deltas: Option<Py<Deltas>>,
    ) -> PyResult<Self> {
        Self::from_rust(
            py,
            AstReactionAst::new(
                lhs.map(|value| value.bind(py).borrow().inner().clone())
                    .unwrap_or_default(),
                deltas
                    .map(|value| value.bind(py).borrow().to_rust())
                    .unwrap_or_default(),
            ),
        )
    }

    /// Parse a reaction from its EDN representation.
    #[staticmethod]
    #[pyo3(signature = (text, *, defaults=None))]
    fn parse(py: Python<'_>, text: &str, defaults: Option<ReactionDefaults>) -> PyResult<Self> {
        let defaults = defaults.unwrap_or_else(ReactionDefaults::new).to_rust();
        let reaction = AstReactionDsl::from_str(text)
            .map_err(parse_error)?
            .into_ast(&defaults);
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
        let defaults = defaults.unwrap_or_else(ReactionDefaults::new).to_rust();
        let dsl = AstReactionDsl::from_str(text).map_err(parse_error)?;
        let metadata = ReactionMetadata::from_rust(dsl.metadata().clone());
        Ok((Self::from_rust(py, dsl.into_ast(&defaults))?, metadata))
    }

    /// Render a canonical positional DSL representation without entity
    /// keywords or atom aliases.
    #[pyo3(signature = (*, defaults=None))]
    fn render(&self, py: Python<'_>, defaults: Option<ReactionDefaults>) -> String {
        let defaults = defaults.unwrap_or_else(ReactionDefaults::new).to_rust();
        AstReactionDsl::from_ast(&self.to_rust(py), &defaults).to_string()
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
        let defaults = defaults.unwrap_or_else(ReactionDefaults::new).to_rust();
        let lowered = AstReactionDsl::from_ast(&self.to_rust(py), &defaults)
            .into_parts()
            .0;
        AstReactionDsl::new(lowered, metadata.to_rust())
            .map(|dsl| dsl.to_string())
            .map_err(metadata_error)
    }

    /// Construct a reaction by comparing two molecule snapshots under an atom correspondence.
    #[staticmethod]
    fn from_sides(
        py: Python<'_>,
        lhs: Py<MoleculeAst>,
        rhs: Py<MoleculeAst>,
        atom_pairs: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        let lhs = lhs.bind(py).borrow().inner().clone();
        let rhs = rhs.bind(py).borrow().inner().clone();
        let atom_pairs = atom_pairs
            .try_iter()?
            .map(|item| item?.extract::<(usize, usize)>())
            .collect::<PyResult<Vec<_>>>()?;
        let atom = atom_correspondence(atom_pairs, lhs.atoms().count(), rhs.atoms().count())?;

        Self::from_rust(py, AstReactionAst::from_sides(lhs, rhs, atom))
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
        let chemistry_model =
            chemistry_model.map_or_else(GraphChemistryModel::default, |model| model.to_rust());
        let resolve_config =
            resolve_config.map_or_else(GraphResolveConfig::default, ResolveConfig::to_rust);
        let reaction =
            ingest_reaction_smiles_with(source, &io_config, &chemistry_model, &resolve_config)
                .map_err(reaction_smiles_input_error)?;

        Self::from_rust(py, reaction)
    }

    /// The live left-hand molecule component.
    #[getter]
    fn lhs(&self, py: Python<'_>) -> Py<MoleculeAst> {
        self.lhs.clone_ref(py)
    }

    /// Replace the left-hand molecule with a detached snapshot.
    #[setter]
    fn set_lhs(slf: Py<Self>, py: Python<'_>, value: Py<MoleculeAst>) -> PyResult<()> {
        let resolved = Py::new(
            py,
            MoleculeAst::from_inner(value.bind(py).borrow().inner().clone()),
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
        let resolved = Py::new(py, Deltas::from_rust(value.bind(py).borrow().to_rust()))?;
        slf.borrow_mut(py).deltas = resolved;
        Ok(())
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

    /// Return one-shot application results with eager matching and lazy derivation construction.
    #[pyo3(signature = (host, *, config=None))]
    fn apply(
        &self,
        py: Python<'_>,
        host: Py<MoleculeAst>,
        config: Option<ReactionApplicationConfig>,
    ) -> PyResult<Py<ReactionApplicationIter>> {
        let reaction = self.to_rust(py);
        let host = host.bind(py).borrow().inner().clone();
        reaction
            .validate_application(&host)
            .map_err(|error| InvalidStructureError::new_err(error.to_string()))?;
        let config = config.unwrap_or_default().to_rust();
        let correspondences = reaction.lhs.substructure_matches(&host, config);

        Py::new(
            py,
            ReactionApplicationIter::new(reaction, host, correspondences),
        )
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
        Ok(format!("ReactionAst(lhs={lhs}, deltas={deltas})"))
    }
}

impl_py_canonicalize!(
    ReactionAst,
    AstReactionAst,
    |value: &ReactionAst, py: Python<'_>| -> PyResult<AstReactionAst> { Ok(value.to_rust(py)) },
    |py: Python<'_>, value: AstReactionAst| -> PyResult<ReactionAst> {
        ReactionAst::from_rust(py, value)
    }
);

impl ReactionAst {
    /// Wrap a Rust reaction in fresh Python-owned components.
    pub(crate) fn from_rust(py: Python<'_>, reaction: AstReactionAst) -> PyResult<Self> {
        Ok(Self {
            lhs: Py::new(py, MoleculeAst::from_inner(reaction.lhs))?,
            deltas: Py::new(py, Deltas::from_rust(reaction.deltas))?,
        })
    }

    /// Snapshot the current Python-owned components as a Rust reaction.
    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstReactionAst {
        AstReactionAst::new(
            self.lhs.bind(py).borrow().inner().clone(),
            self.deltas.bind(py).borrow().to_rust(),
        )
    }
}

/// One owned firing of a reaction, exposed as an immutable result value.
#[pyclass(eq, frozen, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReactionDerivation(AstReactionDerivation);

#[pymethods]
impl ReactionDerivation {
    /// The molecule matched by the reaction, as a fresh snapshot.
    #[getter]
    fn lhs(&self) -> MoleculeAst {
        MoleculeAst::from_inner(self.0.lhs().clone())
    }

    /// The molecule produced by the reaction, as a fresh snapshot.
    #[getter]
    fn rhs(&self) -> MoleculeAst {
        MoleculeAst::from_inner(self.0.rhs().clone())
    }

    /// The correspondence between the two molecule sides, as a fresh snapshot.
    #[getter]
    fn comap(&self) -> PyMoleculeCorrespondence {
        PyMoleculeCorrespondence::from_rust(self.0.comap().clone())
    }

    /// The atom-level correspondence, as a fresh snapshot.
    #[getter]
    fn atom_map(&self) -> PyCorrespondence {
        PyCorrespondence::from_rust(self.0.atom_map())
    }

    /// Return the reverse derivation with swapped sides and inverted correspondence.
    fn reverse(&self) -> Self {
        Self::from_rust(self.to_rust().reverse())
    }

    /// Chain this derivation onto a compatible following derivation.
    fn chain(&self, next: &Self) -> Self {
        let first = self.to_rust();
        let next = next.to_rust();
        Self::from_rust(first.chain(&next))
    }

    /// Recover the reaction rule represented by this concrete firing.
    fn to_reaction(&self, py: Python<'_>) -> PyResult<ReactionAst> {
        ReactionAst::from_rust(py, self.to_rust().to_reaction())
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
    pub(crate) fn from_rust(derivation: AstReactionDerivation) -> Self {
        Self(derivation)
    }

    pub(crate) fn to_rust(&self) -> AstReactionDerivation {
        self.0.clone()
    }
}

/// One-shot application results over an eagerly enumerated correspondence set.
#[pyclass(skip_from_py_object)]
pub(crate) struct ReactionApplicationIter {
    reaction: AstReactionAst,
    host: AstMoleculeAst,
    correspondences: IntoIter<AstMoleculeCorrespondence>,
}

impl ReactionApplicationIter {
    pub(crate) fn new(
        reaction: AstReactionAst,
        host: AstMoleculeAst,
        correspondences: Vec<AstMoleculeCorrespondence>,
    ) -> Self {
        Self {
            reaction,
            host,
            correspondences: correspondences.into_iter(),
        }
    }
}

#[pymethods]
impl ReactionApplicationIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> PyResult<Option<ReactionDerivation>> {
        loop {
            let Some(correspondence) = self.correspondences.next() else {
                return Ok(None);
            };
            match self.reaction.apply_at(&self.host, &correspondence) {
                Ok(derivation) => return Ok(Some(ReactionDerivation::from_rust(derivation))),
                Err(error) if error.is_match_rejection() => {}
                Err(AstApplyError::Transaction(error)) => {
                    self.correspondences = Vec::new().into_iter();
                    return Err(transaction_error(error));
                }
                Err(error) => {
                    self.correspondences = Vec::new().into_iter();
                    return Err(PyRuntimeError::new_err(error.to_string()));
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use pyo3::exceptions::{PyTypeError, PyValueError};
    use pyo3::types::{PyDict, PyList};
    use rstest::{fixture, rstest};
    use umol_ast::ast::{
        AromaticSystemAst as AstAromaticSystemAst, AromaticSystemDelta as AstAromaticSystemDelta,
        AromaticSystemId as AstAromaticSystemId, AtomAst as AstAtomAst, AtomDelta as AstAtomDelta,
        AtomFieldChange as AstAtomFieldChange, AtomId as AstAtomId, BondAst as AstBondAst,
        BondDelta as AstBondDelta, BondFieldChange as AstBondFieldChange, BondId as AstBondId,
        Canonicalize, Constraint as AstConstraint, ConstraintDelta as AstConstraintDelta,
        DativeBondAst as AstDativeBondAst, DativeBondDelta as AstDativeBondDelta,
        DativeBondId as AstDativeBondId, Delta as AstDelta, Deltas as AstDeltas,
        Entity as AstEntity, MoleculeAst as AstMoleculeAst,
        MoleculeConstraint as AstMoleculeConstraint,
        MoleculeCorrespondence as AstMoleculeCorrespondence, MoleculeEntries as AstMoleculeEntries,
        MulticenterBondAst as AstMulticenterBondAst,
        MulticenterBondDelta as AstMulticenterBondDelta, MulticenterBondId as AstMulticenterBondId,
        NoncovalentBondAst as AstNoncovalentBondAst,
        NoncovalentBondDelta as AstNoncovalentBondDelta, NoncovalentBondId as AstNoncovalentBondId,
        NoncovalentBondKind as AstNoncovalentBondKind, StereoAtomAst as AstStereoAtomAst,
        StereoAtomDelta as AstStereoAtomDelta, StereoAtomId as AstStereoAtomId,
        StereoBondAst as AstStereoBondAst, StereoBondDelta as AstStereoBondDelta,
        StereoBondId as AstStereoBondId, StereoCoset as AstStereoCoset,
        StereoKind as AstStereoKind, StereoLigand as AstStereoLigand,
        StereoLigandKind as AstStereoLigandKind, ValueAst as AstValueAst,
    };
    use umol_ast::dsl::{
        AtomDsl as AstAtomDsl, MoleculeMetadata as AstMoleculeMetadata,
        ReactionMetadata as AstReactionMetadata,
    };
    use umol_ast::{mol_dsl, mol_dsl_ground};
    use umol_chem::element::Element as ChemElement;
    use umol_graph::ingest::ingest_smiles;

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
        AstSubstructureMatchAlgorithm::GraphAndOverlays,
        GraphCoreSubgraphIsomorphismAlgorithm::Vf2Rdkit,
        GraphCoreRelevantCycleEnumerationAlgorithm::Vismara,
        ReactionApplicationConfig::default()
    )]
    #[case::incidence_arc_match(
        AstSubstructureMatchAlgorithm::Incidence,
        GraphCoreSubgraphIsomorphismAlgorithm::ArcMatch { path_length: 6 },
        GraphCoreRelevantCycleEnumerationAlgorithm::Vismara,
        ReactionApplicationConfig::new(
            SubstructureMatchAlgorithm::Incidence(),
            SubgraphIsomorphismAlgorithm::ArcMatch { path_length: 6 },
            RelevantCycleEnumerationAlgorithm::Vismara(),
        ),
    )]
    fn test_reaction_application_config_from_rust(
        #[case] match_algorithm: AstSubstructureMatchAlgorithm,
        #[case] subgraph_isomorphism_algorithm: GraphCoreSubgraphIsomorphismAlgorithm,
        #[case] relevant_cycle_algorithm: GraphCoreRelevantCycleEnumerationAlgorithm,
        #[case] expected: ReactionApplicationConfig,
    ) {
        assert_eq!(
            ReactionApplicationConfig::from_rust(AstSubstructureMatchConfig {
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
            AstSubstructureMatchConfig {
                match_algorithm: AstSubstructureMatchAlgorithm::Incidence,
                subgraph_isomorphism_algorithm: expected_subgraph_isomorphism_algorithm,
                relevant_cycle_algorithm: GraphCoreRelevantCycleEnumerationAlgorithm::Vismara,
            }
        );
    }

    #[rstest]
    #[case::empty(Vec::new(), 0, 0, Vec::new())]
    #[case::partial(vec![(1, 2)], 3, 4, vec![(AstAtomId(1), AstAtomId(2))])]
    #[case::total(
        vec![(0, 1), (1, 0)],
        2,
        2,
        vec![(AstAtomId(0), AstAtomId(1)), (AstAtomId(1), AstAtomId(0))],
    )]
    #[case::unsorted(
        vec![(2, 0), (0, 2)],
        3,
        3,
        vec![(AstAtomId(0), AstAtomId(2)), (AstAtomId(2), AstAtomId(0))],
    )]
    fn test_atom_correspondence(
        #[case] pairs: Vec<(usize, usize)>,
        #[case] lhs_count: usize,
        #[case] rhs_count: usize,
        #[case] expected_matched_pairs: Vec<(AstAtomId, AstAtomId)>,
    ) {
        let correspondence = atom_correspondence(pairs, lhs_count, rhs_count).unwrap();

        assert_eq!(
            correspondence.matched_pairs(),
            expected_matched_pairs.as_slice()
        );
        assert_eq!(correspondence.left_count(), lhs_count);
        assert_eq!(correspondence.right_count(), rhs_count);
    }

    #[rstest]
    #[case::duplicate_left(
        vec![(0, 0), (0, 1)],
        2,
        2,
        "duplicate left atom id 0",
    )]
    #[case::duplicate_right(
        vec![(0, 1), (1, 1)],
        2,
        2,
        "duplicate right atom id 1",
    )]
    #[case::left_out_of_range(
        vec![(2, 0)],
        2,
        1,
        "left atom id 2 out of range for 2 atoms",
    )]
    #[case::right_out_of_range(
        vec![(0, 1)],
        1,
        1,
        "right atom id 1 out of range for 1 atoms",
    )]
    fn test_atom_correspondence_error(
        #[case] pairs: Vec<(usize, usize)>,
        #[case] lhs_count: usize,
        #[case] rhs_count: usize,
        #[case] expected: &str,
    ) {
        Python::attach(|py| {
            let error = atom_correspondence(pairs, lhs_count, rhs_count)
                .err()
                .unwrap();

            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                expected
            );
        });
    }

    #[rstest]
    #[case::empty(None, None, AstReactionAst::default())]
    #[case::populated(
        Some(AstMoleculeAst::from_entries(AstMoleculeEntries {
            atoms: vec![AstAtomAst::from_element(ChemElement::C)],
            ..Default::default()
        })),
        Some(vec![AstDelta::Atom(AstAtomDelta::Add {
            id: AstAtomId(1),
            ast: AstAtomAst::from_element(ChemElement::O),
        })].into_iter().collect()),
        AstReactionAst::new(
            AstMoleculeAst::from_entries(AstMoleculeEntries {
                atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                ..Default::default()
            }),
            vec![AstDelta::Atom(AstAtomDelta::Add {
                id: AstAtomId(1),
                ast: AstAtomAst::from_element(ChemElement::O),
            })].into_iter().collect(),
        ),
    )]
    fn test_reaction_ast_new(
        #[case] lhs: Option<AstMoleculeAst>,
        #[case] deltas: Option<AstDeltas>,
        #[case] expected: AstReactionAst,
    ) {
        Python::attach(|py| {
            let lhs = lhs.map(|value| Py::new(py, MoleculeAst::from_inner(value)).unwrap());
            let deltas = deltas.map(|value| Py::new(py, Deltas::from_rust(value)).unwrap());

            let reaction = ReactionAst::new(py, lhs, deltas).unwrap();

            assert_eq!(reaction.to_rust(py), expected);
        });
    }

    #[rstest]
    fn test_reaction_ast_new_snapshot() {
        Python::attach(|py| {
            let lhs = Py::new(
                py,
                MoleculeAst::from_inner(AstMoleculeAst::from_entries(AstMoleculeEntries {
                    atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                    ..Default::default()
                })),
            )
            .unwrap();
            let deltas = Py::new(
                py,
                Deltas::from_rust(
                    vec![AstDelta::Atom(AstAtomDelta::Add {
                        id: AstAtomId(1),
                        ast: AstAtomAst::from_element(ChemElement::O),
                    })]
                    .into_iter()
                    .collect(),
                ),
            )
            .unwrap();
            let expected = AstReactionAst::new(
                lhs.bind(py).borrow().inner().clone(),
                deltas.bind(py).borrow().to_rust(),
            );

            let reaction =
                ReactionAst::new(py, Some(lhs.clone_ref(py)), Some(deltas.clone_ref(py))).unwrap();
            *lhs.bind(py).borrow_mut().inner_mut() = AstMoleculeAst::new();
            let delta = into_py_variant(
                py,
                Delta::from_rust(
                    py,
                    &AstDelta::Atom(AstAtomDelta::Add {
                        id: AstAtomId(2),
                        ast: AstAtomAst::from_element(ChemElement::N),
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
            AstDelta::Atom(AstAtomDelta::Add {
                id: AstAtomId(2),
                ast: AstAtomAst::from_element(ChemElement::N),
            }),
            AstDelta::Atom(AstAtomDelta::Remove {
                id: AstAtomId(1),
                ast: AstAtomAst::from_element(ChemElement::O),
            }),
        ],
    )]
    #[case::atom_modify(
        r##"{:lhs {:atoms ["Br#c0"]} :deltas [{:atom {:modify [0 "#c-1"]}}]}"##,
        1,
        vec![AstDelta::Atom(AstAtomDelta::ModifyField {
            id: AstAtomId(0),
            change: AstAtomFieldChange::Charge {
                old: AstValueAst::Lit(0),
                new: AstValueAst::Lit(-1),
            },
        })],
    )]
    #[case::stereo_mirror(
        r##"{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]} :deltas [{:stereo-atom {:mirror [0 :tetrahedral]}}]}"##,
        5,
        vec![AstDelta::StereoAtom(AstStereoAtomDelta::Mirror {
            id: AstStereoAtomId(0),
            kind: AstStereoKind::Tetrahedral,
        })],
    )]
    #[case::molecule_constraint(
        r##"{:lhs {:atoms ["C"]} :deltas [{:constraint {:add {:connected {}}}}]}"##,
        1,
        vec![AstDelta::Constraint(AstConstraintDelta::Add(
            AstConstraint::Molecule(AstMoleculeConstraint::Connected { atoms: None }),
        ))],
    )]
    fn test_reaction_ast_parse(
        #[case] text: &str,
        #[case] atom_count: usize,
        #[case] expected_deltas: Vec<AstDelta>,
    ) {
        Python::attach(|py| {
            let reaction = ReactionAst::parse(py, text, None).unwrap().to_rust(py);

            assert_eq!(reaction.lhs.atoms().count(), atom_count);
            assert_eq!(reaction.deltas.as_slice(), expected_deltas.as_slice());
        });
    }

    #[rstest]
    fn test_reaction_ast_parse_error() {
        Python::attach(|py| {
            let error = ReactionAst::parse(py, "not edn", None).err().unwrap();

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
    fn test_reaction_ast_parse_defaults(
        #[case] text: &str,
        #[case] defaults: Option<ReactionDefaults>,
        #[case] expected: &str,
    ) {
        Python::attach(|py| {
            assert_eq!(
                ReactionAst::parse(py, text, defaults).unwrap().to_rust(py),
                expected.parse::<AstReactionAst>().unwrap()
            );
        });
    }

    #[rstest]
    fn test_reaction_ast_parse_with_metadata() {
        Python::attach(|py| {
            let (reaction, metadata) = ReactionAst::parse_with_metadata(
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
                metadata.lhs().keyword(AstEntity::Atom(AstAtomId(0))),
                Some("lhs")
            );
            assert_eq!(
                metadata.delta_keyword(AstEntity::Atom(AstAtomId(1))),
                Some("added")
            );
            assert_eq!(
                metadata.lhs().atom_alias("lhs-c"),
                Some(&AstAtomDsl(AstAtomAst::from_element(ChemElement::C)))
            );
            assert_eq!(
                metadata.atom_alias("delta-o"),
                Some(&AstAtomDsl(AstAtomAst::from_element(ChemElement::O)))
            );
        });
    }

    #[rstest]
    fn test_reaction_ast_parse_with_metadata_defaults() {
        Python::attach(|py| {
            let (reaction, metadata) = ReactionAst::parse_with_metadata(
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
                ReactionMetadata::from_rust(AstReactionMetadata::default())
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
    fn test_reaction_ast_render(
        #[case] reaction: AstReactionAst,
        #[case] defaults: Option<ReactionDefaults>,
        #[case] expected: &str,
    ) {
        Python::attach(|py| {
            assert_eq!(
                ReactionAst::from_rust(py, reaction)
                    .unwrap()
                    .render(py, defaults),
                expected
            );
        });
    }

    #[rstest]
    fn test_reaction_ast_render_with_metadata() {
        Python::attach(|py| {
            let reaction = ReactionAst::from_rust(
                py,
                r#"{:lhs {:atoms ["C"]} :deltas [{:atom {:add "O"}}]}"#
                    .parse()
                    .unwrap(),
            )
            .unwrap();
            let mut lhs = AstMoleculeMetadata::new();
            lhs.set_keyword(AstEntity::Atom(AstAtomId(0)), "lhs")
                .unwrap();
            let mut metadata = AstReactionMetadata::from(lhs);
            metadata
                .set_delta_keyword(AstEntity::Atom(AstAtomId(1)), "added")
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
    fn test_reaction_ast_render_with_metadata_error() {
        Python::attach(|py| {
            let reaction =
                ReactionAst::from_rust(py, r#"{:lhs {:atoms ["C"]} :deltas []}"#.parse().unwrap())
                    .unwrap();
            let mut metadata = AstReactionMetadata::default();
            metadata
                .set_delta_keyword(AstEntity::Atom(AstAtomId(1)), "absent")
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
        AstMoleculeAst::from_entries(AstMoleculeEntries {
            atoms: vec![AstAtomAst::from_element(ChemElement::C)],
            ..Default::default()
        }),
        AstMoleculeAst::from_entries(AstMoleculeEntries {
            atoms: vec![AstAtomAst::from_element(ChemElement::C)],
            ..Default::default()
        }),
        vec![(0, 0)],
        AstReactionAst::new(
            AstMoleculeAst::from_entries(AstMoleculeEntries {
                atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                ..Default::default()
            }),
            AstDeltas::default(),
        ),
    )]
    #[case::partial_correspondence(
        AstMoleculeAst::from_entries(AstMoleculeEntries {
            atoms: vec![
                AstAtomAst::from_element(ChemElement::C),
                AstAtomAst::from_element(ChemElement::O),
            ],
            ..Default::default()
        }),
        AstMoleculeAst::from_entries(AstMoleculeEntries {
            atoms: vec![
                AstAtomAst::from_element(ChemElement::C),
                AstAtomAst::from_element(ChemElement::N),
            ],
            ..Default::default()
        }),
        vec![(0, 0)],
        AstReactionAst::new(
            AstMoleculeAst::from_entries(AstMoleculeEntries {
                atoms: vec![
                    AstAtomAst::from_element(ChemElement::C),
                    AstAtomAst::from_element(ChemElement::O),
                ],
                ..Default::default()
            }),
            vec![
                AstDelta::Atom(AstAtomDelta::Remove {
                    id: AstAtomId(1),
                    ast: AstAtomAst::from_element(ChemElement::O),
                }),
                AstDelta::Atom(AstAtomDelta::Add {
                    id: AstAtomId(2),
                    ast: AstAtomAst::from_element(ChemElement::N),
                }),
            ]
            .into_iter()
            .collect(),
        ),
    )]
    #[case::bond_order(
        AstMoleculeAst::from_entries(AstMoleculeEntries {
            atoms: vec![
                AstAtomAst::from_element(ChemElement::C),
                AstAtomAst::from_element(ChemElement::C),
            ],
            bonds: vec![(AstAtomId(0), AstAtomId(1), AstBondAst::from_order(1))],
            ..Default::default()
        }),
        AstMoleculeAst::from_entries(AstMoleculeEntries {
            atoms: vec![
                AstAtomAst::from_element(ChemElement::C),
                AstAtomAst::from_element(ChemElement::C),
            ],
            bonds: vec![(AstAtomId(0), AstAtomId(1), AstBondAst::from_order(2))],
            ..Default::default()
        }),
        vec![(0, 0), (1, 1)],
        AstReactionAst::new(
            AstMoleculeAst::from_entries(AstMoleculeEntries {
                atoms: vec![
                    AstAtomAst::from_element(ChemElement::C),
                    AstAtomAst::from_element(ChemElement::C),
                ],
                bonds: vec![(AstAtomId(0), AstAtomId(1), AstBondAst::from_order(1))],
                ..Default::default()
            }),
            vec![AstDelta::Bond(AstBondDelta::ModifyField {
                id: AstBondId(0),
                change: AstBondFieldChange::Order {
                    old: AstValueAst::Lit(1),
                    new: AstValueAst::Lit(2),
                },
            })]
            .into_iter()
            .collect(),
        ),
    )]
    fn test_reaction_ast_from_sides(
        #[case] lhs: AstMoleculeAst,
        #[case] rhs: AstMoleculeAst,
        #[case] atom_pairs: Vec<(usize, usize)>,
        #[case] expected: AstReactionAst,
    ) {
        Python::attach(|py| {
            let lhs_before = lhs.clone();
            let rhs_before = rhs.clone();
            let lhs = Py::new(py, MoleculeAst::from_inner(lhs)).unwrap();
            let rhs = Py::new(py, MoleculeAst::from_inner(rhs)).unwrap();

            let atom_pairs = PyList::new(py, atom_pairs).unwrap();
            let reaction = ReactionAst::from_sides(
                py,
                lhs.clone_ref(py),
                rhs.clone_ref(py),
                atom_pairs.as_any(),
            )
            .unwrap();

            assert_eq!(reaction.to_rust(py), expected);
            assert_eq!(*lhs.bind(py).borrow().inner(), lhs_before);
            assert_eq!(*rhs.bind(py).borrow().inner(), rhs_before);
            assert_ne!(reaction.lhs.as_ptr(), lhs.as_ptr());
        });
    }

    #[rstest]
    #[case::dative_bond(
        r#"{:atoms ["N" "B"] :bonds []}"#,
        r#"{:atoms ["N" "B"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type "1"}]}"#,
        vec![(0, 0), (1, 1)],
        vec![AstDelta::DativeBond(AstDativeBondDelta::Add {
            id: AstDativeBondId(0),
            donors: vec![AstAtomId(0)],
            acceptor: AstAtomId(1),
            ast: AstDativeBondAst::from_order(1),
        })],
    )]
    #[case::aromatic_system(
        r#"{:atoms ["C" "C"] :bonds []}"#,
        r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "[1,1]"}]}"#,
        vec![(0, 0), (1, 1)],
        vec![AstDelta::AromaticSystem(AstAromaticSystemDelta::Add {
            id: AstAromaticSystemId(0),
            atoms: vec![AstAtomId(0), AstAtomId(1)],
            ast: AstAromaticSystemAst::from_electrons(vec![1, 1]),
        })],
    )]
    #[case::multicenter_bond(
        r#"{:atoms ["B" "H" "B"] :bonds []}"#,
        r#"{:atoms ["B" "H" "B"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :type "[3,5,7]"}]}"#,
        vec![(0, 0), (1, 1), (2, 2)],
        vec![AstDelta::MulticenterBond(AstMulticenterBondDelta::Add {
            id: AstMulticenterBondId(0),
            atoms: vec![AstAtomId(0), AstAtomId(1), AstAtomId(2)],
            ast: AstMulticenterBondAst::from_electrons(vec![3, 5, 7]),
        })],
    )]
    #[case::noncovalent_bond(
        r#"{:atoms ["O" "O"] :bonds []}"#,
        r#"{:atoms ["O" "O"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}]}"#,
        vec![(0, 0), (1, 1)],
        vec![AstDelta::NoncovalentBond(AstNoncovalentBondDelta::Add {
            id: AstNoncovalentBondId(0),
            atoms: [AstAtomId(0), AstAtomId(1)],
            ast: AstNoncovalentBondAst::from_kind(AstNoncovalentBondKind::HydrogenBond),
        })],
    )]
    #[case::stereo_atom(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds []}"#,
        r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]}"#,
        vec![(0, 0), (1, 1), (2, 2), (3, 3), (4, 4)],
        vec![AstDelta::StereoAtom(AstStereoAtomDelta::Add {
            id: AstStereoAtomId(0),
            site: AstAtomId(0),
            ligands: vec![
                AstStereoLigand::new(AstAtomId(1), AstStereoLigandKind::Atom),
                AstStereoLigand::new(AstAtomId(2), AstStereoLigandKind::Atom),
                AstStereoLigand::new(AstAtomId(3), AstStereoLigandKind::Atom),
                AstStereoLigand::new(AstAtomId(4), AstStereoLigandKind::Atom),
            ],
            ast: AstStereoAtomAst::new(AstStereoKind::Tetrahedral, AstStereoCoset::Lit(1)),
        })],
    )]
    #[case::stereo_bond(
        r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]}"#,
        r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3] :type "Ct1"}]}"#,
        vec![(0, 0), (1, 1), (2, 2), (3, 3)],
        vec![AstDelta::StereoBond(AstStereoBondDelta::Add {
            id: AstStereoBondId(0),
            site: AstBondId(1),
            ligands: vec![
                AstStereoLigand::new(AstAtomId(0), AstStereoLigandKind::Atom),
                AstStereoLigand::new(AstAtomId(3), AstStereoLigandKind::Atom),
            ],
            ast: AstStereoBondAst::new(AstStereoKind::CisTrans, AstStereoCoset::Lit(1)),
        })],
    )]
    #[case::molecule_constraint(
        r#"{:atoms ["C"] :bonds []}"#,
        r#"{:atoms ["C"] :bonds [] :constraints [{:connected {}}]}"#,
        vec![(0, 0)],
        vec![AstDelta::Constraint(AstConstraintDelta::Add(
            AstConstraint::Molecule(AstMoleculeConstraint::Connected { atoms: None }),
        ))],
    )]
    fn test_reaction_ast_from_sides_entities(
        #[case] lhs: &str,
        #[case] rhs: &str,
        #[case] atom_pairs: Vec<(usize, usize)>,
        #[case] expected_deltas: Vec<AstDelta>,
    ) {
        Python::attach(|py| {
            let lhs = lhs.parse::<AstMoleculeAst>().unwrap();
            let rhs = rhs.parse::<AstMoleculeAst>().unwrap();
            let atom_pairs = PyList::new(py, atom_pairs).unwrap();
            let reaction = ReactionAst::from_sides(
                py,
                Py::new(py, MoleculeAst::from_inner(lhs.clone())).unwrap(),
                Py::new(py, MoleculeAst::from_inner(rhs)).unwrap(),
                atom_pairs.as_any(),
            )
            .unwrap();

            assert_eq!(
                reaction.to_rust(py),
                AstReactionAst::new(lhs, expected_deltas.into_iter().collect())
            );
        });
    }

    #[rstest]
    fn test_reaction_ast_from_sides_snapshot() {
        Python::attach(|py| {
            let lhs_before = AstMoleculeAst::from_entries(AstMoleculeEntries {
                atoms: vec![
                    AstAtomAst::from_element(ChemElement::C),
                    AstAtomAst::from_element(ChemElement::O),
                ],
                ..Default::default()
            });
            let rhs_before = AstMoleculeAst::from_entries(AstMoleculeEntries {
                atoms: vec![
                    AstAtomAst::from_element(ChemElement::C),
                    AstAtomAst::from_element(ChemElement::N),
                ],
                ..Default::default()
            });
            let lhs = Py::new(py, MoleculeAst::from_inner(lhs_before.clone())).unwrap();
            let rhs = Py::new(py, MoleculeAst::from_inner(rhs_before.clone())).unwrap();
            let atom_pairs = PyList::new(py, [(0, 0)]).unwrap();
            let reaction = ReactionAst::from_sides(
                py,
                lhs.clone_ref(py),
                rhs.clone_ref(py),
                atom_pairs.as_any(),
            )
            .unwrap();
            let expected = reaction.to_rust(py);

            *lhs.bind(py).borrow_mut().inner_mut() = AstMoleculeAst::new();
            *rhs.bind(py).borrow_mut().inner_mut() = AstMoleculeAst::new();

            assert_eq!(reaction.to_rust(py), expected);
            assert_ne!(reaction.lhs.as_ptr(), lhs.as_ptr());

            *reaction.lhs.bind(py).borrow_mut().inner_mut() =
                AstMoleculeAst::from_entries(AstMoleculeEntries {
                    atoms: vec![AstAtomAst::from_element(ChemElement::F)],
                    ..Default::default()
                });
            let delta = into_py_variant(
                py,
                Delta::from_rust(
                    py,
                    &AstDelta::Atom(AstAtomDelta::Add {
                        id: AstAtomId(3),
                        ast: AstAtomAst::from_element(ChemElement::Cl),
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
                AstMoleculeAst::from_entries(AstMoleculeEntries {
                    atoms: vec![AstAtomAst::from_element(ChemElement::F)],
                    ..Default::default()
                })
            );
            assert_eq!(
                changed.deltas.as_slice().last(),
                Some(&AstDelta::Atom(AstAtomDelta::Add {
                    id: AstAtomId(3),
                    ast: AstAtomAst::from_element(ChemElement::Cl),
                }))
            );
        });
    }

    #[rstest]
    fn test_reaction_ast_components() {
        Python::attach(|py| {
            let reaction = Py::new(py, ReactionAst::new(py, None, None).unwrap()).unwrap();
            let first_lhs = reaction.bind(py).borrow().lhs(py);
            let second_lhs = reaction.bind(py).borrow().lhs(py);
            let first_deltas = reaction.bind(py).borrow().deltas(py);
            let second_deltas = reaction.bind(py).borrow().deltas(py);

            *first_lhs.bind(py).borrow_mut().inner_mut() =
                AstMoleculeAst::from_entries(AstMoleculeEntries {
                    atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                    ..Default::default()
                });
            let delta = into_py_variant(
                py,
                Delta::from_rust(
                    py,
                    &AstDelta::Atom(AstAtomDelta::Add {
                        id: AstAtomId(1),
                        ast: AstAtomAst::from_element(ChemElement::O),
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
                AstReactionAst::new(
                    AstMoleculeAst::from_entries(AstMoleculeEntries {
                        atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                        ..Default::default()
                    }),
                    vec![AstDelta::Atom(AstAtomDelta::Add {
                        id: AstAtomId(1),
                        ast: AstAtomAst::from_element(ChemElement::O),
                    })]
                    .into_iter()
                    .collect(),
                )
            );
        });
    }

    #[rstest]
    fn test_reaction_ast_set_components() {
        Python::attach(|py| {
            let reaction = Py::new(py, ReactionAst::new(py, None, None).unwrap()).unwrap();
            let lhs = Py::new(
                py,
                MoleculeAst::from_inner(AstMoleculeAst::from_entries(AstMoleculeEntries {
                    atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                    ..Default::default()
                })),
            )
            .unwrap();
            let deltas = Py::new(
                py,
                Deltas::from_rust(
                    vec![AstDelta::Atom(AstAtomDelta::Add {
                        id: AstAtomId(1),
                        ast: AstAtomAst::from_element(ChemElement::O),
                    })]
                    .into_iter()
                    .collect(),
                ),
            )
            .unwrap();
            let expected = AstReactionAst::new(
                lhs.bind(py).borrow().inner().clone(),
                deltas.bind(py).borrow().to_rust(),
            );

            ReactionAst::set_lhs(reaction.clone_ref(py), py, lhs.clone_ref(py)).unwrap();
            ReactionAst::set_deltas(reaction.clone_ref(py), py, deltas.clone_ref(py)).unwrap();
            *lhs.bind(py).borrow_mut().inner_mut() = AstMoleculeAst::new();
            let delta = into_py_variant(
                py,
                Delta::from_rust(
                    py,
                    &AstDelta::Atom(AstAtomDelta::Add {
                        id: AstAtomId(2),
                        ast: AstAtomAst::from_element(ChemElement::N),
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
    fn test_reaction_ast_set_components_self() {
        Python::attach(|py| {
            let expected = AstReactionAst::new(
                AstMoleculeAst::from_entries(AstMoleculeEntries {
                    atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                    ..Default::default()
                }),
                vec![AstDelta::Atom(AstAtomDelta::Add {
                    id: AstAtomId(1),
                    ast: AstAtomAst::from_element(ChemElement::O),
                })]
                .into_iter()
                .collect(),
            );
            let reaction =
                Py::new(py, ReactionAst::from_rust(py, expected.clone()).unwrap()).unwrap();
            let own_lhs = reaction.bind(py).borrow().lhs(py);
            let own_deltas = reaction.bind(py).borrow().deltas(py);

            ReactionAst::set_lhs(reaction.clone_ref(py), py, own_lhs).unwrap();
            ReactionAst::set_deltas(reaction.clone_ref(py), py, own_deltas).unwrap();

            assert_eq!(reaction.bind(py).borrow().to_rust(py), expected);
        });
    }

    #[rstest]
    fn test_reaction_ast_canonicalize() {
        Python::attach(|py| {
            let source = ReactionAst::from_rust(
                py,
                AstReactionAst::new(
                    AstMoleculeAst::from_entries(AstMoleculeEntries {
                        atoms: vec![AstAtomAst::from_element(ChemElement::C).with_charge(0)],
                        ..Default::default()
                    }),
                    vec![
                        AstDelta::Atom(AstAtomDelta::ModifyField {
                            id: AstAtomId(0),
                            change: AstAtomFieldChange::Charge {
                                old: AstValueAst::Lit(0),
                                new: AstValueAst::Lit(1),
                            },
                        }),
                        AstDelta::Atom(AstAtomDelta::ModifyField {
                            id: AstAtomId(0),
                            change: AstAtomFieldChange::Charge {
                                old: AstValueAst::Lit(1),
                                new: AstValueAst::Lit(2),
                            },
                        }),
                    ]
                    .into_iter()
                    .collect(),
                ),
            )
            .unwrap();
            let before = source.to_rust(py);
            let expected = AstReactionAst::new(
                before.lhs.clone(),
                vec![AstDelta::Atom(AstAtomDelta::ModifyField {
                    id: AstAtomId(0),
                    change: AstAtomFieldChange::Charge {
                        old: AstValueAst::Lit(0),
                        new: AstValueAst::Lit(2),
                    },
                })]
                .into_iter()
                .collect(),
            );

            let canonical = source.canonicalize(py).unwrap();
            let twice = canonical.canonicalize(py).unwrap();

            assert_eq!(canonical.to_rust(py), expected);
            assert_eq!(twice.to_rust(py), expected);
            assert_eq!(source.to_rust(py), before);
            assert_ne!(canonical.lhs.as_ptr(), source.lhs.as_ptr());
            assert_ne!(canonical.deltas.as_ptr(), source.deltas.as_ptr());
        });
    }

    #[rstest]
    fn test_reaction_ast_canonicalize_error() {
        Python::attach(|py| {
            let source = ReactionAst::from_rust(
                py,
                AstReactionAst::new(
                    AstMoleculeAst::from_entries(AstMoleculeEntries {
                        atoms: vec![AstAtomAst::from_element(ChemElement::C).with_charge(0)],
                        ..Default::default()
                    }),
                    vec![
                        AstDelta::Atom(AstAtomDelta::ModifyField {
                            id: AstAtomId(0),
                            change: AstAtomFieldChange::Charge {
                                old: AstValueAst::Lit(0),
                                new: AstValueAst::Lit(1),
                            },
                        }),
                        AstDelta::Atom(AstAtomDelta::ModifyField {
                            id: AstAtomId(0),
                            change: AstAtomFieldChange::Charge {
                                old: AstValueAst::Lit(2),
                                new: AstValueAst::Lit(3),
                            },
                        }),
                    ]
                    .into_iter()
                    .collect(),
                ),
            )
            .unwrap();
            let before = source.to_rust(py);

            let error = source.canonicalize(py).err().unwrap();

            assert!(error.is_instance_of::<ContradictionError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "reached a contradiction"
            );
            assert_eq!(source.to_rust(py), before);
        });
    }

    #[rstest]
    fn test_reaction_ast_reverse() {
        Python::attach(|py| {
            let source = ReactionAst::parse(
                py,
                r##"{:lhs {:atoms ["C" "O"]} :deltas [{:atom {:add "N"}} {:atom {:remove 1}}]}"##,
                None,
            )
            .unwrap();
            let before = source.to_rust(py);
            let expected_roundtrip = before.clone().canonicalize().unwrap();

            let reversed = source.reverse(py).unwrap();
            let roundtrip = reversed.reverse(py).unwrap();

            assert_eq!(
                reversed.to_rust(py).lhs,
                AstMoleculeAst::from_entries(AstMoleculeEntries {
                    atoms: vec![
                        AstAtomAst::from_element(ChemElement::C),
                        AstAtomAst::from_element(ChemElement::N),
                    ],
                    ..Default::default()
                })
            );
            assert_eq!(
                roundtrip.to_rust(py).canonicalize().unwrap(),
                expected_roundtrip
            );
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
    fn test_reaction_ast_compose(
        #[case] first: &str,
        #[case] second: &str,
        #[case] expected: Vec<&str>,
    ) {
        Python::attach(|py| {
            let first = ReactionAst::parse(py, first, None).unwrap();
            let second = ReactionAst::parse(py, second, None).unwrap();
            let expected: Vec<AstReactionAst> = expected
                .into_iter()
                .map(|reaction| AstReactionAst::from_str(reaction).unwrap())
                .collect();

            let actual: Vec<AstReactionAst> = first
                .compose(py, &second, None)
                .unwrap()
                .iter()
                .map(|reaction| reaction.to_rust(py))
                .collect();

            assert_eq!(actual, expected);
        });
    }

    #[rstest]
    #[case::direct(ReactionCompositionConfig::new(
        CommonSubgraphEnumerationAlgorithm::DirectBacktracking()
    ))]
    #[case::modular_product(ReactionCompositionConfig::new(
        CommonSubgraphEnumerationAlgorithm::ModularProductBacktracking()
    ))]
    fn test_reaction_ast_compose_config(#[case] config: ReactionCompositionConfig) {
        Python::attach(|py| {
            let first = ReactionAst::parse(
                py,
                r##"{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+"]}}]}"##,
                None,
            )
            .unwrap();
            let second = ReactionAst::parse(
                py,
                r##"{:lhs {:atoms ["C#c+"]} :deltas [{:atom {:modify [0 "#c+2"]}}]}"##,
                None,
            )
            .unwrap();

            assert_eq!(
                first
                    .compose(py, &second, Some(config))
                    .unwrap()
                    .into_iter()
                    .map(|reaction| reaction.to_rust(py))
                    .collect::<Vec<_>>(),
                vec![
                    AstReactionAst::from_str(
                        r##"{:lhs {:atoms ["C#c0" "C#c+"]} :deltas [{:atom {:modify [0 "#c+"]}} {:atom {:modify [1 "#c+2"]}}]}"##,
                    )
                    .unwrap(),
                    AstReactionAst::from_str(
                        r##"{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+2"]}}]}"##,
                    )
                    .unwrap(),
                ]
            );
        });
    }

    #[rstest]
    fn test_reaction_ast_compose_default() {
        Python::attach(|py| {
            let first = Py::new(
                py,
                ReactionAst::parse(
                    py,
                    r##"{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+"]}}]}"##,
                    None,
                )
                .unwrap(),
            )
            .unwrap();
            let second = Py::new(
                py,
                ReactionAst::parse(
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

            let omitted: Vec<Py<ReactionAst>> = first
                .bind(py)
                .call_method1("compose", (second.clone_ref(py),))
                .unwrap()
                .extract()
                .unwrap();
            let explicit: Vec<Py<ReactionAst>> = first
                .bind(py)
                .call_method("compose", (second,), Some(&kwargs))
                .unwrap()
                .extract()
                .unwrap();
            let omitted: Vec<AstReactionAst> = omitted
                .iter()
                .map(|reaction| reaction.bind(py).borrow().to_rust(py))
                .collect();
            let explicit: Vec<AstReactionAst> = explicit
                .iter()
                .map(|reaction| reaction.bind(py).borrow().to_rust(py))
                .collect();

            assert_eq!(omitted, explicit);
            assert_eq!(
                omitted,
                vec![
                    AstReactionAst::from_str(
                        r##"{:lhs {:atoms ["C#c0" "C#c+"]} :deltas [{:atom {:modify [0 "#c+"]}} {:atom {:modify [1 "#c+2"]}}]}"##,
                    )
                    .unwrap(),
                    AstReactionAst::from_str(
                        r##"{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+2"]}}]}"##,
                    )
                    .unwrap(),
                ]
            );
        });
    }

    #[rstest]
    fn test_reaction_ast_compose_snapshot() {
        Python::attach(|py| {
            let first = ReactionAst::parse(
                py,
                r##"{:lhs {:atoms ["C#c0"]} :deltas [{:atom {:modify [0 "#c+"]}}]}"##,
                None,
            )
            .unwrap();
            let second = ReactionAst::parse(
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
                *composite.lhs.bind(py).borrow_mut().inner_mut() =
                    AstMoleculeAst::from_entries(AstMoleculeEntries {
                        atoms: vec![AstAtomAst::from_element(ChemElement::F)],
                        ..Default::default()
                    });
                let delta = into_py_variant(
                    py,
                    Delta::from_rust(
                        py,
                        &AstDelta::Atom(AstAtomDelta::Add {
                            id: AstAtomId(8),
                            ast: AstAtomAst::from_element(ChemElement::Cl),
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
                    AstMoleculeAst::from_entries(AstMoleculeEntries {
                        atoms: vec![AstAtomAst::from_element(ChemElement::F)],
                        ..Default::default()
                    })
                );
                assert_eq!(
                    composite.to_rust(py).deltas.as_slice().last(),
                    Some(&AstDelta::Atom(AstAtomDelta::Add {
                        id: AstAtomId(8),
                        ast: AstAtomAst::from_element(ChemElement::Cl),
                    }))
                );
            }

            assert_eq!(first.to_rust(py), first_before);
            assert_eq!(second.to_rust(py), second_before);
        });
    }

    #[rstest]
    fn test_reaction_ast_apply(
        reaction_application: (
            AstReactionAst,
            AstMoleculeAst,
            Vec<AstMoleculeCorrespondence>,
        ),
    ) {
        let (expected_reaction, expected_host, _) = reaction_application;
        Python::attach(|py| {
            let reaction = ReactionAst::from_rust(py, expected_reaction.clone()).unwrap();
            let host = Py::new(py, MoleculeAst::from_inner(expected_host.clone())).unwrap();
            let application = reaction.apply(py, host.clone_ref(py), None).unwrap();

            assert_eq!(application.borrow(py).correspondences.len(), 2);
            assert_eq!(reaction.to_rust(py), expected_reaction);
            assert_eq!(host.bind(py).borrow().inner(), &expected_host);

            let first = application.borrow_mut(py).__next__().unwrap().unwrap();
            assert_eq!(application.borrow(py).correspondences.len(), 1);
            let second = application.borrow_mut(py).__next__().unwrap().unwrap();
            assert_eq!(application.borrow(py).correspondences.len(), 0);
            assert_eq!(application.borrow_mut(py).__next__().unwrap(), None);
            assert_eq!(
                [first.rhs().inner().clone(), second.rhs().inner().clone()],
                [
                    AstMoleculeAst::from_entries(AstMoleculeEntries {
                        atoms: vec![
                            AstAtomAst::from_element(ChemElement::C).with_charge(1),
                            AstAtomAst::from_element(ChemElement::C),
                        ],
                        ..Default::default()
                    }),
                    AstMoleculeAst::from_entries(AstMoleculeEntries {
                        atoms: vec![
                            AstAtomAst::from_element(ChemElement::C),
                            AstAtomAst::from_element(ChemElement::C).with_charge(1),
                        ],
                        ..Default::default()
                    }),
                ]
            );
        });
    }

    #[rstest]
    fn test_reaction_ast_apply_snapshot(
        reaction_application: (
            AstReactionAst,
            AstMoleculeAst,
            Vec<AstMoleculeCorrespondence>,
        ),
    ) {
        let (expected_reaction, expected_host, _) = reaction_application;
        Python::attach(|py| {
            let mut reaction = ReactionAst::from_rust(py, expected_reaction).unwrap();
            let host = Py::new(py, MoleculeAst::from_inner(expected_host)).unwrap();
            let application = reaction.apply(py, host.clone_ref(py), None).unwrap();

            *reaction.lhs.bind(py).borrow_mut().inner_mut() =
                AstMoleculeAst::from_entries(AstMoleculeEntries {
                    atoms: vec![AstAtomAst::from_element(ChemElement::N)],
                    ..Default::default()
                });
            reaction.deltas = Py::new(py, Deltas::from_rust(AstDeltas::default())).unwrap();
            *host.bind(py).borrow_mut().inner_mut() =
                AstMoleculeAst::from_entries(AstMoleculeEntries {
                    atoms: vec![AstAtomAst::from_element(ChemElement::F)],
                    ..Default::default()
                });

            let products: Vec<AstMoleculeAst> = std::iter::from_fn(|| {
                application
                    .borrow_mut(py)
                    .__next__()
                    .unwrap()
                    .map(|derivation| derivation.rhs().inner().clone())
            })
            .collect();
            assert_eq!(
                products,
                vec![
                    AstMoleculeAst::from_entries(AstMoleculeEntries {
                        atoms: vec![
                            AstAtomAst::from_element(ChemElement::C).with_charge(1),
                            AstAtomAst::from_element(ChemElement::C),
                        ],
                        ..Default::default()
                    }),
                    AstMoleculeAst::from_entries(AstMoleculeEntries {
                        atoms: vec![
                            AstAtomAst::from_element(ChemElement::C),
                            AstAtomAst::from_element(ChemElement::C).with_charge(1),
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
    fn test_reaction_ast_apply_config(
        reaction_application: (
            AstReactionAst,
            AstMoleculeAst,
            Vec<AstMoleculeCorrespondence>,
        ),
        #[case] config: ReactionApplicationConfig,
    ) {
        let (reaction, host, _) = reaction_application;
        Python::attach(|py| {
            let reaction = ReactionAst::from_rust(py, reaction).unwrap();
            let host = Py::new(py, MoleculeAst::from_inner(host)).unwrap();
            let application = reaction.apply(py, host, Some(config)).unwrap();

            let products: Vec<AstMoleculeAst> = std::iter::from_fn(|| {
                application
                    .borrow_mut(py)
                    .__next__()
                    .unwrap()
                    .map(|derivation| derivation.rhs().inner().clone())
            })
            .collect();
            assert_eq!(
                products,
                vec![
                    AstMoleculeAst::from_entries(AstMoleculeEntries {
                        atoms: vec![
                            AstAtomAst::from_element(ChemElement::C).with_charge(1),
                            AstAtomAst::from_element(ChemElement::C),
                        ],
                        ..Default::default()
                    }),
                    AstMoleculeAst::from_entries(AstMoleculeEntries {
                        atoms: vec![
                            AstAtomAst::from_element(ChemElement::C),
                            AstAtomAst::from_element(ChemElement::C).with_charge(1),
                        ],
                        ..Default::default()
                    }),
                ]
            );
        });
    }

    #[rstest]
    fn test_reaction_ast_apply_error() {
        Python::attach(|py| {
            let reaction = ReactionAst::new(py, None, None).unwrap();
            let host = Py::new(
                py,
                MoleculeAst::from_inner(AstMoleculeAst::from_entries(AstMoleculeEntries {
                    atoms: vec![
                        AstAtomAst::from_element(ChemElement::C),
                        AstAtomAst::from_element(ChemElement::O),
                    ],
                    bonds: vec![
                        (AstAtomId(0), AstAtomId(1), AstBondAst::from_order(1)),
                        (AstAtomId(0), AstAtomId(1), AstBondAst::from_order(2)),
                    ],
                    ..Default::default()
                })),
            )
            .unwrap();

            let error = reaction.apply(py, host, None).err().unwrap();

            assert!(error.is_instance_of::<InvalidStructureError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "invalid host: bond: parallel bonds on atoms [AtomId(0), AtomId(1)]"
            );
        });
    }

    #[fixture]
    fn ethanol_deoxygenation() -> AstReactionAst {
        let ethanol = ingest_smiles("CCO").unwrap();
        let oxygen = ethanol.atom(AstAtomId(2)).ast.clone();
        let bond = ethanol.bond(AstBondId(1)).ast.clone();
        AstReactionAst::new(
            ethanol,
            AstDeltas::from_iter([
                AstDelta::Atom(AstAtomDelta::Remove {
                    id: AstAtomId(2),
                    ast: oxygen,
                }),
                AstDelta::Bond(AstBondDelta::Remove {
                    id: AstBondId(1),
                    atoms: [AstAtomId(1), AstAtomId(2)],
                    ast: bond,
                }),
            ]),
        )
    }

    #[fixture]
    fn ethanol_identity() -> AstReactionAst {
        AstReactionAst::new(ingest_smiles("CCO").unwrap(), AstDeltas::new())
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
    fn test_reaction_ast_combined_fingerprint_difference(
        ethanol_deoxygenation: AstReactionAst,
        #[case] config: ReactionCombinedFingerprintConfig,
        #[case] expected_entries: Vec<(u128, i32)>,
    ) {
        Python::attach(|py| {
            let reaction = ReactionAst::from_rust(py, ethanol_deoxygenation).unwrap();
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
    fn test_reaction_ast_combined_fingerprint_disjoint_union(
        ethanol_deoxygenation: AstReactionAst,
        #[case] config: ReactionCombinedFingerprintConfig,
        #[case] expected_ids: Vec<(ReactionSide, u128)>,
    ) {
        Python::attach(|py| {
            let reaction = ReactionAst::from_rust(py, ethanol_deoxygenation).unwrap();
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
    fn test_reaction_ast_combined_fingerprint_difference_identity(
        ethanol_identity: AstReactionAst,
        #[case] config: ReactionCombinedFingerprintConfig,
    ) {
        Python::attach(|py| {
            let reaction = ReactionAst::from_rust(py, ethanol_identity).unwrap();
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
    fn test_reaction_ast_combined_fingerprint_disjoint_union_identity(
        ethanol_identity: AstReactionAst,
        #[case] config: ReactionCombinedFingerprintConfig,
        #[case] expected_ids: Vec<(ReactionSide, u128)>,
    ) {
        Python::attach(|py| {
            let reaction = ReactionAst::from_rust(py, ethanol_identity).unwrap();
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
        AstReactionAst::new(
            mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
            AstDeltas::new(),
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
        AstReactionAst::new(
            mol_dsl_ground!(r#"{:atoms ["C #h4"] :bonds []}"#),
            AstDeltas::from_iter([AstDelta::Atom(AstAtomDelta::ModifyField {
                id: AstAtomId(0),
                change: AstAtomFieldChange::Charge {
                    old: AstValueAst::Lit(0),
                    new: AstValueAst::Undetermined,
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
        AstReactionAst::new(
            mol_dsl_ground!(r#"{:atoms ["C #h4"] :bonds []}"#),
            AstDeltas::from_iter([AstDelta::Atom(AstAtomDelta::ModifyField {
                id: AstAtomId(0),
                change: AstAtomFieldChange::Charge {
                    old: AstValueAst::Lit(1),
                    new: AstValueAst::Lit(0),
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
    fn test_reaction_ast_combined_fingerprint_error(
        #[case] input: AstReactionAst,
        #[case] config: ReactionCombinedFingerprintConfig,
        #[case] expected_type: &str,
        #[case] expected_message: &str,
    ) {
        Python::attach(|py| {
            let reaction = ReactionAst::from_rust(py, input).unwrap();
            let error = reaction.combined_fingerprint(py, config).unwrap_err();

            assert_eq!(error.get_type(py).name().unwrap(), expected_type);
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                expected_message
            );
        });
    }

    #[rstest]
    fn test_reaction_ast_eq() {
        Python::attach(|py| {
            let empty = ReactionAst::new(py, None, None).unwrap();
            let other_empty = ReactionAst::new(py, None, None).unwrap();
            let populated = ReactionAst::from_rust(
                py,
                AstReactionAst::new(
                    AstMoleculeAst::from_entries(AstMoleculeEntries {
                        atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                        ..Default::default()
                    }),
                    AstDeltas::new(),
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
        AstReactionAst::default(),
        r##"{:deltas [] :lhs {:atoms [] :bonds []}}"##
    )]
    #[case::populated(
        AstReactionAst::new(
            AstMoleculeAst::from_entries(AstMoleculeEntries {
                atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                ..Default::default()
            }),
            vec![AstDelta::Atom(AstAtomDelta::Add {
                id: AstAtomId(1),
                ast: AstAtomAst::from_element(ChemElement::O),
            })].into_iter().collect(),
        ),
        r##"{:deltas [{:atom {:add "O"}}] :lhs {:atoms ["C"] :bonds []}}"##,
    )]
    fn test_reaction_ast_str(#[case] input: AstReactionAst, #[case] expected: &str) {
        Python::attach(|py| {
            let reaction = ReactionAst::from_rust(py, input).unwrap();

            assert_eq!(reaction.__str__(py), expected);
            assert_eq!(reaction.__str__(py), reaction.render(py, None));
        });
    }

    #[rstest]
    fn test_reaction_ast_str_components() {
        Python::attach(|py| {
            let reaction = ReactionAst::from_rust(
                py,
                AstReactionAst::new(
                    AstMoleculeAst::from_entries(AstMoleculeEntries {
                        atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                        ..Default::default()
                    }),
                    AstDeltas::new(),
                ),
            )
            .unwrap();

            *reaction.lhs.bind(py).borrow_mut().inner_mut() =
                AstMoleculeAst::from_entries(AstMoleculeEntries {
                    atoms: vec![AstAtomAst::from_element(ChemElement::C).with_charge(1)],
                    ..Default::default()
                });
            let delta = into_py_variant(
                py,
                Delta::from_rust(
                    py,
                    &AstDelta::Atom(AstAtomDelta::Add {
                        id: AstAtomId(1),
                        ast: AstAtomAst::from_element(ChemElement::O),
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
    #[case::stereo_mirror(
        r##"{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]} :deltas [{:stereo-atom {:mirror [0 :tetrahedral]}}]}"##
    )]
    #[case::molecule_constraint(
        r##"{:lhs {:atoms ["C"]} :deltas [{:constraint {:add {:connected {}}}}]}"##
    )]
    fn test_reaction_ast_str_roundtrip(#[case] text: &str) {
        Python::attach(|py| {
            let first = ReactionAst::parse(py, text, None).unwrap();

            let canonical = first.__str__(py);
            let second = ReactionAst::parse(py, &canonical, None).unwrap();

            assert!(first.__eq__(&second, py));
            assert_eq!(second.__str__(py), canonical);
        });
    }

    #[rstest]
    fn test_reaction_ast_repr() {
        Python::attach(|py| {
            let reaction = ReactionAst::from_rust(
                py,
                AstReactionAst::new(
                    AstMoleculeAst::from_entries(AstMoleculeEntries {
                        atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                        ..Default::default()
                    }),
                    vec![AstDelta::Atom(AstAtomDelta::Add {
                        id: AstAtomId(1),
                        ast: AstAtomAst::from_element(ChemElement::O),
                    })]
                    .into_iter()
                    .collect(),
                ),
            )
            .unwrap();

            assert_eq!(
                reaction.__repr__(py).unwrap(),
                "ReactionAst(lhs=MoleculeAst(atoms=1, bonds=0), deltas=Deltas([Delta.Atom(AtomDelta.Add(id=1, ast=AtomAst.parse('O')))]))"
            );
        });
    }

    #[rstest]
    #[case::empty(AstReactionAst::default())]
    #[case::populated(AstReactionAst::new(
        AstMoleculeAst::from_entries(AstMoleculeEntries {
            atoms: vec![AstAtomAst::from_element(ChemElement::C)],
            ..Default::default()
        }),
        vec![AstDelta::Atom(AstAtomDelta::Add {
            id: AstAtomId(1),
            ast: AstAtomAst::from_element(ChemElement::O),
        })]
        .into_iter()
        .collect(),
    ))]
    fn test_reaction_ast_from_rust(#[case] expected: AstReactionAst) {
        Python::attach(|py| {
            let reaction = ReactionAst::from_rust(py, expected.clone()).unwrap();

            assert_eq!(reaction.to_rust(py), expected);
        });
    }

    #[rstest]
    fn test_reaction_ast_to_rust() {
        Python::attach(|py| {
            let expected = AstReactionAst::new(
                AstMoleculeAst::from_entries(AstMoleculeEntries {
                    atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                    ..Default::default()
                }),
                vec![AstDelta::Atom(AstAtomDelta::Add {
                    id: AstAtomId(1),
                    ast: AstAtomAst::from_element(ChemElement::O),
                })]
                .into_iter()
                .collect(),
            );
            let reaction = ReactionAst::from_rust(py, expected.clone()).unwrap();

            let mut snapshot = reaction.to_rust(py);
            snapshot.lhs = AstMoleculeAst::new();
            snapshot.deltas = AstDeltas::new();

            assert_eq!(reaction.to_rust(py), expected);
        });
    }

    #[rstest]
    fn test_reaction_ast_to_rust_roundtrip() {
        Python::attach(|py| {
            let expected = AstReactionAst::new(
                AstMoleculeAst::from_entries(AstMoleculeEntries {
                    atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                    ..Default::default()
                }),
                vec![AstDelta::Atom(AstAtomDelta::Add {
                    id: AstAtomId(1),
                    ast: AstAtomAst::from_element(ChemElement::O),
                })]
                .into_iter()
                .collect(),
            );
            let python =
                Py::new(py, ReactionAst::from_rust(py, expected.clone()).unwrap()).unwrap();

            let rust = python.bind(py).borrow().to_rust(py);
            let roundtrip = Py::new(py, ReactionAst::from_rust(py, rust).unwrap()).unwrap();

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
    fn derivation_and_host() -> (AstReactionDerivation, AstMoleculeAst) {
        let pattern = AstMoleculeAst::from_entries(AstMoleculeEntries {
            atoms: vec![
                AstAtomAst::from_element(ChemElement::C),
                AstAtomAst::from_element(ChemElement::C),
            ],
            bonds: vec![(AstAtomId(0), AstAtomId(1), AstBondAst::from_order(1))],
            ..Default::default()
        });
        let host = pattern.clone();
        let reaction = AstReactionAst::new(
            pattern.clone(),
            AstDeltas::from_iter([AstDelta::Bond(AstBondDelta::ModifyField {
                id: AstBondId(0),
                change: AstBondFieldChange::Order {
                    old: AstValueAst::Lit(1),
                    new: AstValueAst::Lit(2),
                },
            })]),
        );
        let correspondence = AstMoleculeCorrespondence::induce(
            &pattern,
            &host,
            Correspondence::new(
                vec![(AstAtomId(0), AstAtomId(0)), (AstAtomId(1), AstAtomId(1))],
                2,
                2,
            )
            .expect("correspondence producer preserves partial-bijection invariants"),
        );
        let derivation = reaction.apply_at(&host, &correspondence).unwrap();
        (derivation, host)
    }

    #[rstest]
    fn test_reaction_derivation_observations(
        derivation_and_host: (AstReactionDerivation, AstMoleculeAst),
    ) {
        let (expected, mut host) = derivation_and_host;
        let derivation = ReactionDerivation::from_rust(expected.clone());

        assert_eq!(derivation.lhs().inner(), expected.lhs());
        assert_eq!(derivation.rhs().inner(), expected.rhs());
        assert_eq!(
            derivation.comap(),
            PyMoleculeCorrespondence::from_rust(expected.comap().clone())
        );
        assert_eq!(
            derivation.atom_map(),
            PyCorrespondence::from_rust(expected.atom_map())
        );

        *host.atom_mut(AstAtomId(0)).ast = AstAtomAst::from_element(ChemElement::F);
        let mut lhs = derivation.lhs();
        *lhs.inner_mut().atom_mut(AstAtomId(0)).ast = AstAtomAst::from_element(ChemElement::N);

        assert_eq!(derivation.to_rust(), expected);
        assert_ne!(derivation.lhs().inner(), &host);
        assert_ne!(derivation.lhs().inner(), lhs.inner());
    }

    #[rstest]
    fn test_reaction_derivation_reverse(
        derivation_and_host: (AstReactionDerivation, AstMoleculeAst),
    ) {
        let (expected, _) = derivation_and_host;
        let derivation = ReactionDerivation::from_rust(expected.clone());
        let reversed = derivation.reverse();
        let mut reversed_lhs = reversed.lhs();
        *reversed_lhs.inner_mut().atom_mut(AstAtomId(0)).ast =
            AstAtomAst::from_element(ChemElement::N);

        assert_eq!(reversed.to_rust(), expected.reverse());
        assert_eq!(derivation.to_rust(), expected);
        assert_ne!(reversed.lhs().inner(), reversed_lhs.inner());
    }

    #[rstest]
    fn test_reaction_derivation_chain(
        derivation_and_host: (AstReactionDerivation, AstMoleculeAst),
    ) {
        let (first, _) = derivation_and_host;
        let middle = first.rhs().clone();
        let reaction = AstReactionAst::new(
            middle.clone(),
            AstDeltas::from_iter([AstDelta::Bond(AstBondDelta::ModifyField {
                id: AstBondId(0),
                change: AstBondFieldChange::Order {
                    old: AstValueAst::Lit(2),
                    new: AstValueAst::Lit(3),
                },
            })]),
        );
        let correspondence = AstMoleculeCorrespondence::induce(
            &middle,
            &middle,
            Correspondence::new(
                vec![(AstAtomId(0), AstAtomId(0)), (AstAtomId(1), AstAtomId(1))],
                2,
                2,
            )
            .expect("correspondence producer preserves partial-bijection invariants"),
        );
        let second = reaction.apply_at(&middle, &correspondence).unwrap();
        let first_value = ReactionDerivation::from_rust(first.clone());
        let second_value = ReactionDerivation::from_rust(second.clone());
        let chained = first_value.chain(&second_value);
        let mut chained_rhs = chained.rhs();
        *chained_rhs.inner_mut().atom_mut(AstAtomId(0)).ast =
            AstAtomAst::from_element(ChemElement::N);

        assert_eq!(chained.to_rust(), first.chain(&second));
        assert_eq!(first_value.to_rust(), first);
        assert_eq!(second_value.to_rust(), second);
        assert_ne!(chained.rhs().inner(), chained_rhs.inner());
    }

    #[rstest]
    fn test_reaction_derivation_to_reaction(
        derivation_and_host: (AstReactionDerivation, AstMoleculeAst),
    ) {
        let (expected_derivation, _) = derivation_and_host;
        let expected_reaction = AstReactionAst::new(
            AstMoleculeAst::from_entries(AstMoleculeEntries {
                atoms: vec![
                    AstAtomAst::from_element(ChemElement::C),
                    AstAtomAst::from_element(ChemElement::C),
                ],
                bonds: vec![(AstAtomId(0), AstAtomId(1), AstBondAst::from_order(1))],
                ..Default::default()
            }),
            AstDeltas::from_iter([AstDelta::Bond(AstBondDelta::ModifyField {
                id: AstBondId(0),
                change: AstBondFieldChange::Order {
                    old: AstValueAst::Lit(1),
                    new: AstValueAst::Lit(2),
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

            *first.lhs.bind(py).borrow_mut().inner_mut() = AstMoleculeAst::new();
            let delta = into_py_variant(
                py,
                Delta::from_rust(
                    py,
                    &AstDelta::Atom(AstAtomDelta::Add {
                        id: AstAtomId(2),
                        ast: AstAtomAst::from_element(ChemElement::O),
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
            assert_eq!(derivation.to_rust(), expected_derivation);
        });
    }

    #[rstest]
    fn test_reaction_derivation_value(
        derivation_and_host: (AstReactionDerivation, AstMoleculeAst),
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
                    "ReactionDerivation(lhs=MoleculeAst(atoms=2, bonds=1), ",
                    "rhs=MoleculeAst(atoms=2, bonds=1), ",
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
        derivation_and_host: (AstReactionDerivation, AstMoleculeAst),
    ) {
        let (expected, _) = derivation_and_host;
        assert_eq!(
            ReactionDerivation::from_rust(expected.clone()).to_rust(),
            expected
        );
    }

    #[fixture]
    fn reaction_application() -> (
        AstReactionAst,
        AstMoleculeAst,
        Vec<AstMoleculeCorrespondence>,
    ) {
        let reaction = AstReactionAst::new(
            AstMoleculeAst::from_entries(AstMoleculeEntries {
                atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                ..Default::default()
            }),
            AstDeltas::from_iter([AstDelta::Atom(AstAtomDelta::ModifyField {
                id: AstAtomId(0),
                change: AstAtomFieldChange::Charge {
                    old: AstValueAst::Undetermined,
                    new: AstValueAst::Lit(1),
                },
            })]),
        );
        let host = AstMoleculeAst::from_entries(AstMoleculeEntries {
            atoms: vec![
                AstAtomAst::from_element(ChemElement::C),
                AstAtomAst::from_element(ChemElement::C),
            ],
            ..Default::default()
        });
        let correspondences = [AstAtomId(0), AstAtomId(1)]
            .into_iter()
            .map(|host_atom| {
                AstMoleculeCorrespondence::induce(
                    &reaction.lhs,
                    &host,
                    Correspondence::from_images(&[host_atom], host.atoms().count()),
                )
            })
            .collect();
        (reaction, host, correspondences)
    }

    #[rstest]
    fn test_reaction_application_iter_identity(
        reaction_application: (
            AstReactionAst,
            AstMoleculeAst,
            Vec<AstMoleculeCorrespondence>,
        ),
    ) {
        let (reaction, host, correspondences) = reaction_application;
        Python::attach(|py| {
            let application = Py::new(
                py,
                ReactionApplicationIter::new(reaction, host, correspondences),
            )
            .unwrap();

            let iter = application.bind(py).call_method0("__iter__").unwrap();
            assert!(iter.is(application.bind(py)));
        });
    }

    #[rstest]
    fn test_reaction_application_iter(
        reaction_application: (
            AstReactionAst,
            AstMoleculeAst,
            Vec<AstMoleculeCorrespondence>,
        ),
    ) {
        let (reaction, host, correspondences) = reaction_application;
        let mut application = ReactionApplicationIter::new(reaction, host, correspondences);

        let first = application.__next__().unwrap().unwrap();
        let second = application.__next__().unwrap().unwrap();
        assert_eq!(
            [first.rhs().inner().clone(), second.rhs().inner().clone()],
            [
                AstMoleculeAst::from_entries(AstMoleculeEntries {
                    atoms: vec![
                        AstAtomAst::from_element(ChemElement::C).with_charge(1),
                        AstAtomAst::from_element(ChemElement::C),
                    ],
                    ..Default::default()
                }),
                AstMoleculeAst::from_entries(AstMoleculeEntries {
                    atoms: vec![
                        AstAtomAst::from_element(ChemElement::C),
                        AstAtomAst::from_element(ChemElement::C).with_charge(1),
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
        *detached.inner_mut().atom_mut(AstAtomId(0)).ast = AstAtomAst::from_element(ChemElement::F);
        assert_eq!(first.to_rust(), expected_first);
        assert_eq!(second.to_rust(), expected_second);
    }

    #[rstest]
    fn test_reaction_application_iter_empty() {
        let mut application = ReactionApplicationIter::new(
            AstReactionAst::default(),
            AstMoleculeAst::default(),
            Vec::new(),
        );

        assert_eq!(application.__next__().unwrap(), None);
        assert_eq!(application.__next__().unwrap(), None);
    }

    #[rstest]
    fn test_reaction_application_iter_rejection() {
        let reaction = AstReactionAst::new(
            AstMoleculeAst::from_entries(AstMoleculeEntries {
                atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                ..Default::default()
            }),
            AstDeltas::from_iter([AstDelta::Atom(AstAtomDelta::Remove {
                id: AstAtomId(0),
                ast: AstAtomAst::from_element(ChemElement::C),
            })]),
        );
        let host = AstMoleculeAst::from_entries(AstMoleculeEntries {
            atoms: vec![
                AstAtomAst::from_element(ChemElement::C),
                AstAtomAst::from_element(ChemElement::C),
                AstAtomAst::from_element(ChemElement::C),
                AstAtomAst::from_element(ChemElement::O),
            ],
            bonds: vec![(AstAtomId(1), AstAtomId(3), AstBondAst::from_order(1))],
            ..Default::default()
        });
        let correspondences = [AstAtomId(0), AstAtomId(1), AstAtomId(2)]
            .into_iter()
            .map(|host_atom| {
                AstMoleculeCorrespondence::induce(
                    &reaction.lhs,
                    &host,
                    Correspondence::from_images(&[host_atom], host.atoms().count()),
                )
            })
            .collect();
        let mut application = ReactionApplicationIter::new(reaction, host, correspondences);

        let first = application.__next__().unwrap().unwrap();
        let second = application.__next__().unwrap().unwrap();

        assert_eq!(
            [first.rhs().inner().clone(), second.rhs().inner().clone()],
            [
                AstMoleculeAst::from_entries(AstMoleculeEntries {
                    atoms: vec![
                        AstAtomAst::from_element(ChemElement::C),
                        AstAtomAst::from_element(ChemElement::C),
                        AstAtomAst::from_element(ChemElement::O),
                    ],
                    bonds: vec![(AstAtomId(0), AstAtomId(2), AstBondAst::from_order(1))],
                    ..Default::default()
                }),
                AstMoleculeAst::from_entries(AstMoleculeEntries {
                    atoms: vec![
                        AstAtomAst::from_element(ChemElement::C),
                        AstAtomAst::from_element(ChemElement::C),
                        AstAtomAst::from_element(ChemElement::O),
                    ],
                    bonds: vec![(AstAtomId(1), AstAtomId(2), AstBondAst::from_order(1))],
                    ..Default::default()
                }),
            ]
        );
        assert_eq!(application.__next__().unwrap(), None);
    }

    #[rstest]
    fn test_reaction_application_iter_error() {
        let constraint = AstConstraint::Molecule(AstMoleculeConstraint::ChargeSum {
            atoms: Some(vec![AstAtomId(0)]),
            sum: AstValueAst::Lit(0),
        });
        let reaction = AstReactionAst::new(
            AstMoleculeAst::from_entries(AstMoleculeEntries {
                atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                constraints: constraint.clone().into(),
                ..Default::default()
            }),
            AstDeltas::from_iter([AstDelta::Constraint(AstConstraintDelta::Remove(constraint))]),
        );
        let host = AstMoleculeAst::from_entries(AstMoleculeEntries {
            atoms: vec![AstAtomAst::from_element(ChemElement::C)],
            ..Default::default()
        });
        let correspondence = AstMoleculeCorrespondence::induce(
            &reaction.lhs,
            &host,
            Correspondence::from_images(&[AstAtomId(0)], 1),
        );
        let mut application = ReactionApplicationIter::new(
            reaction,
            host,
            vec![correspondence.clone(), correspondence],
        );

        let error = application.__next__().unwrap_err();

        Python::attach(|py| {
            assert!(error.is_instance_of::<TransactionError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "missing constraint entry on remove"
            );
        });
        assert_eq!(application.correspondences.len(), 0);
        assert_eq!(application.__next__().unwrap(), None);
        assert_eq!(application.__next__().unwrap(), None);
    }
}
