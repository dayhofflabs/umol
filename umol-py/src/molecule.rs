//! `Molecule` — an owned graph-IR root, wrapping
//! `umol_graph_ir::ir::Molecule`.

use std::str::FromStr;

use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use umol_graph::fingerprint::PatternFingerprinter as GraphPatternFingerprinter;
use umol_graph::ingest::ingest_smiles_with;
use umol_graph::ops::model::{
    ChemistryModel as GraphChemistryModel, ValenceModel as GraphValenceModel,
};
use umol_graph::ops::resolve::{
    IsotopePolicy as GraphIsotopePolicy, ResolveConfig as GraphResolveConfig,
    Resolver as GraphResolver,
};
use umol_graph_ir::dsl::MoleculeDsl as GraphIrMoleculeDsl;
use umol_graph_ir::ir::{
    AtomId as GraphIrAtomId, BondId as GraphIrBondId, FromIr, IntoIr, Molecule as GraphIrMolecule,
    MoleculeEntries as GraphIrMoleculeEntries, React as GraphIrReact,
};
use umol_io::smiles::SmilesIoConfig as IoSmilesIoConfig;
use umol_utils::solution::Solution as GraphSolution;

use crate::aromatic::{AromaticSystemForm, AromaticSystemViews};
use crate::atom::{AtomForm, AtomViews};
use crate::bond::{BondForm, BondViews};
use crate::compact::MoleculeCompaction;
use crate::constraint::molecule::{Constraint, ConstraintsView};
use crate::correspondence::MoleculeCorrespondence;
use crate::dative::{DativeBondForm, DativeBondViews};
use crate::defaults::MoleculeDefaults;
use crate::edit::Edits;
use crate::error::{
    fingerprint_error, metadata_error, molecule_apply_error, parse_error, smiles_input_error,
    ConsumedError, InvalidStructureError, InvalidatedViewError,
};
use crate::fingerprint::config::{
    HashedFingerprintConfig, PatternFingerprintConfig, StructuralFingerprintConfig,
};
use crate::fingerprint::value::{
    BitFp, CountedHashedFeatureSet, HashedFeatureSet, StructuralFeatureSet,
};
use crate::metadata::MoleculeMetadata;
use crate::model::ChemistryModel;
use crate::multicenter::{MulticenterBondForm, MulticenterBondViews};
use crate::noncovalent::{NoncovalentBondForm, NoncovalentBondViews};
use crate::reaction::{Reaction, ReactionApplicationConfig, ReactionProductsIter};
use crate::resolve::{ResolveConfig, ResolveContradiction, ResolveReport, Solution};
use crate::smiles::SmilesIoConfig;
use crate::stereo::{
    StereoAtomForm, StereoAtomViews, StereoBondForm, StereoBondViews, StereoLigand,
};
use crate::substructure::SubstructureSearchConfig;
use crate::transaction::MoleculeEditor;

/// A molecule: the owned graph-IR root.
///
/// Consuming operations move its contents. Subsequent access raises ConsumedError;
/// existing views raise InvalidatedViewError. Use copy to retain an independent molecule.
#[pyclass]
#[derive(Debug)]
pub struct Molecule {
    value: Option<GraphIrMolecule>,
    counter: u64,
}

#[pymethods]
impl Molecule {
    /// An empty molecule: zero atoms, zero bonds.
    #[new]
    fn new() -> Self {
        Self::from_rust(GraphIrMolecule::new())
    }

    /// Parse a molecule from its EDN representation under explicit construction defaults.
    #[staticmethod]
    #[pyo3(signature = (text, *, defaults=None))]
    fn parse(text: &str, defaults: Option<MoleculeDefaults>) -> PyResult<Self> {
        let defaults = defaults.unwrap_or_else(MoleculeDefaults::new);
        let molecule = GraphIrMoleculeDsl::from_str(text)
            .map_err(parse_error)?
            .into_ir(defaults.to_rust());
        Ok(Self::from_rust(molecule))
    }

    /// Parse a molecule and return `(molecule, metadata)`, retaining entity
    /// keywords and atom aliases for metadata-preserving rendering.
    #[staticmethod]
    #[pyo3(signature = (text, *, defaults=None))]
    fn parse_with_metadata(
        text: &str,
        defaults: Option<MoleculeDefaults>,
    ) -> PyResult<(Self, MoleculeMetadata)> {
        let defaults = defaults.unwrap_or_else(MoleculeDefaults::new);
        let dsl = GraphIrMoleculeDsl::from_str(text).map_err(parse_error)?;
        let metadata = MoleculeMetadata::from_rust(dsl.metadata().clone());
        Ok((Self::from_rust(dsl.into_ir(defaults.to_rust())), metadata))
    }

    /// Render a canonical positional DSL representation without entity
    /// keywords or atom aliases.
    #[pyo3(signature = (*, defaults=None))]
    fn render(&self, defaults: Option<MoleculeDefaults>) -> PyResult<String> {
        let defaults = defaults.unwrap_or_else(MoleculeDefaults::new);
        Ok(GraphIrMoleculeDsl::from_ir(self.to_rust()?, defaults.to_rust()).to_string())
    }

    /// Render a canonical DSL representation with persistent metadata.
    ///
    /// Raises `MetadataError` if the detached metadata is not coherent with
    /// this molecule.
    #[pyo3(signature = (metadata, *, defaults=None))]
    fn render_with_metadata(
        &self,
        metadata: &MoleculeMetadata,
        defaults: Option<MoleculeDefaults>,
    ) -> PyResult<String> {
        let defaults = defaults.unwrap_or_else(MoleculeDefaults::new);
        let lowered = GraphIrMoleculeDsl::from_ir(self.to_rust()?, defaults.to_rust())
            .into_parts()
            .0;
        GraphIrMoleculeDsl::new(lowered, metadata.to_rust().clone())
            .map(|dsl| dsl.to_string())
            .map_err(metadata_error)
    }

    fn __str__(&self) -> PyResult<String> {
        self.render(None)
    }

    /// A molecule from its entries. Each bond is a `(first, second, bond)` triple:
    /// two atom indices into `atoms` and a `BondForm`. Each dative bond is a
    /// `(donors, acceptor, bond)` triple: a list of donor atom indices, one
    /// acceptor atom index, and a `DativeBondForm`. Each aromatic system is an
    /// `(atoms, system)` pair: a list of member atom indices and an `AromaticSystemForm`.
    /// Each multicenter bond is an `(atoms, bond)` pair: a list of member atom indices
    /// and a `MulticenterBondForm`. Each noncovalent bond is a `([first, second], bond)`
    /// pair: the two (unordered) endpoint atom indices and a `NoncovalentBondForm`. Each
    /// stereo atom / stereo bond is a `(site, ligands, value)` triple: the site atom / bond
    /// index, a list of `StereoLigand`s in frame order, and a `StereoAtomForm` / `StereoBondForm`.
    #[staticmethod]
    #[pyo3(signature = (atoms, *, bonds=Vec::new(), dative_bonds=Vec::new(), aromatic_systems=Vec::new(), multicenter_bonds=Vec::new(), noncovalent_bonds=Vec::new(), stereo_atoms=Vec::new(), stereo_bonds=Vec::new(), constraints=Vec::new()))]
    #[allow(clippy::too_many_arguments)] // one argument per entity kind — the full molecule surface
    fn from_entries(
        py: Python<'_>,
        atoms: Vec<Py<AtomForm>>,
        bonds: Vec<(u32, u32, Py<BondForm>)>,
        dative_bonds: Vec<(Vec<u32>, u32, Py<DativeBondForm>)>,
        aromatic_systems: Vec<(Vec<u32>, Py<AromaticSystemForm>)>,
        multicenter_bonds: Vec<(Vec<u32>, Py<MulticenterBondForm>)>,
        noncovalent_bonds: Vec<([u32; 2], Py<NoncovalentBondForm>)>,
        stereo_atoms: Vec<(u32, Vec<StereoLigand>, Py<StereoAtomForm>)>,
        stereo_bonds: Vec<(u32, Vec<StereoLigand>, Py<StereoBondForm>)>,
        constraints: Vec<Py<Constraint>>,
    ) -> PyResult<Self> {
        let ir_atoms = atoms
            .iter()
            .map(|atom| atom.bind(py).borrow().to_rust().clone())
            .collect();
        let ir_bonds = bonds
            .iter()
            .map(|(first, second, bond)| {
                (
                    GraphIrAtomId(*first),
                    GraphIrAtomId(*second),
                    bond.bind(py).borrow().to_rust().clone(),
                )
            })
            .collect();
        let ir_dative = dative_bonds
            .iter()
            .map(|(donors, acceptor, bond)| {
                (
                    donors.iter().map(|&donor| GraphIrAtomId(donor)).collect(),
                    GraphIrAtomId(*acceptor),
                    bond.bind(py).borrow().to_rust().clone(),
                )
            })
            .collect();
        let ir_aromatic = aromatic_systems
            .iter()
            .map(|(atoms, system)| {
                (
                    atoms.iter().map(|&atom| GraphIrAtomId(atom)).collect(),
                    system.bind(py).borrow().to_rust().clone(),
                )
            })
            .collect();
        let ir_multicenter = multicenter_bonds
            .iter()
            .map(|(atoms, bond)| {
                (
                    atoms.iter().map(|&atom| GraphIrAtomId(atom)).collect(),
                    bond.bind(py).borrow().to_rust().clone(),
                )
            })
            .collect();
        let ir_noncovalent = noncovalent_bonds
            .iter()
            .map(|([first, second], bond)| {
                (
                    [GraphIrAtomId(*first), GraphIrAtomId(*second)],
                    bond.bind(py).borrow().to_rust().clone(),
                )
            })
            .collect();
        let ir_stereo_atoms = stereo_atoms
            .iter()
            .map(|(site, ligands, value)| {
                (
                    GraphIrAtomId(*site),
                    ligands.iter().copied().map(StereoLigand::to_rust).collect(),
                    value.bind(py).borrow().to_rust().clone(),
                )
            })
            .collect();
        let ir_stereo_bonds = stereo_bonds
            .iter()
            .map(|(site, ligands, value)| {
                (
                    GraphIrBondId(*site),
                    ligands.iter().copied().map(StereoLigand::to_rust).collect(),
                    value.bind(py).borrow().to_rust().clone(),
                )
            })
            .collect();
        let ir_constraints = constraints
            .iter()
            .map(|constraint| constraint.bind(py).borrow().to_rust(py))
            .collect::<PyResult<Vec<_>>>()?;
        GraphIrMolecule::try_from_entries(GraphIrMoleculeEntries {
            atoms: ir_atoms,
            bonds: ir_bonds,
            dative: ir_dative,
            aromatic: ir_aromatic,
            multicenter: ir_multicenter,
            noncovalent: ir_noncovalent,
            stereo_atoms: ir_stereo_atoms,
            stereo_bonds: ir_stereo_bonds,
            constraints: ir_constraints.into(),
        })
        .map(Self::from_rust)
        .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// Ingest a determined molecule from SMILES under explicit IO, chemistry,
    /// and resolution policies.
    /// Omitted options select OpenSMILES, SMILES valence, and Natural isotope policy.
    /// An explicitly supplied resolve_config is preserved.
    #[staticmethod]
    #[pyo3(signature = (source, *, io_config=None, chemistry_model=None, resolve_config=None))]
    fn from_smiles(
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
        let resolve_config = resolve_config.map_or_else(
            || GraphResolveConfig {
                isotope: GraphIsotopePolicy::Natural,
                ..Default::default()
            },
            ResolveConfig::to_rust,
        );

        ingest_smiles_with(source, &io_config, &chemistry_model, &resolve_config)
            .map(Self::from_rust)
            .map_err(smiles_input_error)
    }

    /// Copy the molecule into an independent owner.
    fn copy(&self) -> PyResult<Self> {
        Ok(Self::from_rust(self.to_rust()?.clone()))
    }

    /// Consume this molecule and move its contents into an editor.
    fn edit(&mut self) -> PyResult<MoleculeEditor> {
        Ok(MoleculeEditor::from_rust(self.take()?.edit()))
    }

    /// Consume this molecule and the edit batch, returning the changed molecule.
    ///
    /// Raises `TransactionError` when the edits cannot be applied and `InvalidStructureError` when
    /// the modified draft cannot be published as a molecule.
    /// Failure leaves the molecule consumed.
    fn apply(&mut self, py: Python<'_>, edits: Py<Edits>) -> PyResult<Self> {
        self.take()?
            .apply(edits.try_borrow_mut(py)?.take()?)
            .map(Self::from_rust)
            .map_err(molecule_apply_error)
    }

    /// Apply the same checked edit batch and return the source-to-result correspondence.
    fn tracked_apply(
        &mut self,
        py: Python<'_>,
        edits: Py<Edits>,
    ) -> PyResult<(Self, MoleculeCorrespondence)> {
        self.take()?
            .tracked_apply(edits.try_borrow_mut(py)?.take()?)
            .map(|(molecule, correspondence)| {
                (
                    Self::from_rust(molecule),
                    MoleculeCorrespondence::from_rust(correspondence),
                )
            })
            .map_err(molecule_apply_error)
    }

    /// Consume prepared batches in order and apply them in one transaction.
    ///
    /// The input iterable is collected before any batch is consumed.
    /// Each batch has its own handle namespace. Once execution starts, existing views are
    /// invalidated, including on failure. Application and integrity failures restore the molecule.
    /// If a batch cannot be consumed, earlier batches remain consumed and the molecule is unchanged.
    fn transact(slf: Py<Self>, py: Python<'_>, batches: &Bound<'_, PyAny>) -> PyResult<()> {
        let batches = batches
            .try_iter()?
            .map(|item| -> PyResult<Py<Edits>> { Ok(item?.cast_into::<Edits>()?.unbind()) })
            .collect::<PyResult<Vec<_>>>()?;
        let batches = batches
            .into_iter()
            .map(|edits| edits.try_borrow_mut(py)?.take())
            .collect::<PyResult<Vec<_>>>()?;
        let mut molecule = slf.try_borrow_mut(py)?;
        molecule.advance_counter()?;
        molecule
            .to_rust_mut()?
            .transact(batches)
            .map_err(molecule_apply_error)
    }

    /// Apply prepared batches and return the transaction-entry to result correspondence.
    ///
    /// Batches are consumed and existing views are invalidated as in transact.
    fn tracked_transact(
        slf: Py<Self>,
        py: Python<'_>,
        batches: &Bound<'_, PyAny>,
    ) -> PyResult<MoleculeCorrespondence> {
        let batches = batches
            .try_iter()?
            .map(|item| -> PyResult<Py<Edits>> { Ok(item?.cast_into::<Edits>()?.unbind()) })
            .collect::<PyResult<Vec<_>>>()?;
        let batches = batches
            .into_iter()
            .map(|edits| edits.try_borrow_mut(py)?.take())
            .collect::<PyResult<Vec<_>>>()?;
        let mut molecule = slf.try_borrow_mut(py)?;
        molecule.advance_counter()?;
        molecule
            .to_rust_mut()?
            .tracked_transact(batches)
            .map(MoleculeCorrespondence::from_rust)
            .map_err(molecule_apply_error)
    }

    /// Combine by disjoint concatenation. For each entity kind, this molecule's ids remain
    /// the prefix and other follows in its original order.
    fn combine(&self, other: &Self) -> PyResult<Self> {
        Ok(Self::from_rust(self.to_rust()?.combine(other.to_rust()?)))
    }

    /// Append other in place, preserving existing ids and each entity kind's order.
    fn combine_from(slf: Py<Self>, py: Python<'_>, other: Py<Self>) -> PyResult<()> {
        let other = other.try_borrow(py)?.to_rust()?.clone();
        let mut molecule = slf.try_borrow_mut(py)?;
        molecule.advance_counter()?;
        molecule.to_rust_mut()?.combine_from(&other);
        Ok(())
    }

    /// Combine an iterable by disjoint concatenation, in input order for each entity kind.
    #[staticmethod]
    fn combine_all(py: Python<'_>, molecules: &Bound<'_, PyAny>) -> PyResult<Self> {
        let molecules = molecules
            .try_iter()?
            .map(|item| -> PyResult<Py<Molecule>> { Ok(item?.cast_into::<Molecule>()?.unbind()) })
            .collect::<PyResult<Vec<_>>>()?;
        let borrowed = molecules
            .iter()
            .map(|molecule| molecule.try_borrow(py))
            .collect::<Result<Vec<_>, _>>()?;
        let molecules = borrowed
            .iter()
            .map(|molecule| molecule.to_rust())
            .collect::<PyResult<Vec<_>>>()?;
        Ok(Self::from_rust(GraphIrMolecule::combine_all(molecules)))
    }

    /// Resolve under a chemistry model, returning the three-valued solution.
    ///
    /// The receiver is never modified; `Determined` carries the resolved copy
    /// together with the tie-break record, `Underdetermined` the surviving
    /// per-atom candidate lists, and `Contradictory` the model's rejection.
    /// `chemistry_model` defaults to `ChemistryModel.default()` — presets are
    /// reader conventions, and a constructed molecule has no format.
    /// Resolution fills what is open; validating committed structure under a
    /// chemistry model is a separate operation.
    #[pyo3(signature = (*, chemistry_model=None, resolve_config=None))]
    fn resolve(
        &self,
        py: Python<'_>,
        chemistry_model: Option<ChemistryModel>,
        resolve_config: Option<ResolveConfig>,
    ) -> PyResult<Solution> {
        let chemistry_model =
            chemistry_model.map_or_else(GraphChemistryModel::default, |model| model.to_rust());
        let resolve_config =
            resolve_config.map_or_else(GraphResolveConfig::default, ResolveConfig::to_rust);
        let mut molecule = self.to_rust()?.clone();
        let solution = GraphResolver::with_config(&chemistry_model, resolve_config)
            .resolve_into(&mut molecule)
            .map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
        Ok(match solution {
            GraphSolution::Determined(report) => Solution::Determined {
                molecule: Py::new(py, Self::from_rust(molecule))?,
                report: ResolveReport::from_rust(&report),
            },
            GraphSolution::Underdetermined(report) => Solution::Underdetermined {
                report: ResolveReport::from_rust(&report),
            },
            GraphSolution::Contradictory(contradiction) => Solution::Contradictory {
                contradiction: ResolveContradiction::from_rust(contradiction),
            },
        })
    }

    /// Apply `reaction` and lazily emit one connected product-component list per match.
    ///
    /// Matching is eager; product construction and splitting are lazy. The returned one-shot
    /// iterator owns snapshots of this molecule and the reaction. Reaction-wide precondition
    /// failures raise `InvalidStructureError` here; failures while realizing a match are raised by
    /// iteration.
    ///
    /// Example: `product_sets = molecule.react(reaction)`.
    #[pyo3(signature = (reaction, *, config=None))]
    fn react(
        &self,
        py: Python<'_>,
        reaction: &Reaction,
        config: Option<ReactionApplicationConfig>,
    ) -> PyResult<Py<ReactionProductsIter>> {
        let reaction = reaction.to_rust(py)?;
        let products = GraphIrReact::react(
            self.to_rust()?,
            &reaction,
            config.unwrap_or_default().to_rust(),
        )
        .map_err(|error| InvalidStructureError::new_err(error.to_string()))?;

        Py::new(py, ReactionProductsIter::from_rust(products))
    }

    /// Combine `reactants` in iterable order, apply `reaction`, and lazily emit product components.
    ///
    /// Any iterable is accepted, including an empty iterable. Matching is eager; product
    /// construction and splitting are lazy. The returned one-shot iterator owns snapshots of all
    /// inputs. Reaction-wide precondition failures raise `InvalidStructureError` here; failures
    /// while realizing a match are raised by iteration.
    ///
    /// Example: `product_sets = Molecule.react_all([first, second], reaction)`.
    #[staticmethod]
    #[pyo3(signature = (reactants, reaction, *, config=None))]
    fn react_all(
        py: Python<'_>,
        reactants: &Bound<'_, PyAny>,
        reaction: &Reaction,
        config: Option<ReactionApplicationConfig>,
    ) -> PyResult<Py<ReactionProductsIter>> {
        let reactants = reactants
            .try_iter()?
            .map(|item| -> PyResult<Py<Molecule>> { Ok(item?.cast_into::<Molecule>()?.unbind()) })
            .collect::<PyResult<Vec<_>>>()?;
        let reactants = reactants
            .iter()
            .map(|molecule| -> PyResult<_> { Ok(molecule.try_borrow(py)?.to_rust()?.clone()) })
            .collect::<PyResult<Vec<_>>>()?;
        let reaction = reaction.to_rust(py)?;
        let products = GraphIrReact::react(
            reactants.as_slice(),
            &reaction,
            config.unwrap_or_default().to_rust(),
        )
        .map_err(|error| InvalidStructureError::new_err(error.to_string()))?;

        Py::new(py, ReactionProductsIter::from_rust(products))
    }

    /// Decompose into conservatively connected components, ordered by lowest source atom id.
    fn split(&self) -> PyResult<Vec<Self>> {
        Ok(self
            .to_rust()?
            .split()
            .into_iter()
            .map(Self::from_rust)
            .collect())
    }

    /// Return the same components as split, paired with source-to-component correspondences.
    fn tracked_split(&self) -> PyResult<Vec<(Self, MoleculeCorrespondence)>> {
        Ok(self
            .to_rust()?
            .tracked_split()
            .into_iter()
            .map(|(component, correspondence)| {
                (
                    Self::from_rust(component),
                    MoleculeCorrespondence::from_rust(correspondence),
                )
            })
            .collect())
    }

    /// Extract the atoms selected by a sub-to-host correspondence, preserving host order.
    fn extract(&self, selection: &MoleculeCorrespondence) -> PyResult<Self> {
        Ok(Self::from_rust(
            self.to_rust()?.extract(selection.to_rust()),
        ))
    }

    /// Return the same extraction and its host-to-result compaction.
    fn tracked_extract(
        &self,
        selection: &MoleculeCorrespondence,
    ) -> PyResult<(Self, MoleculeCompaction)> {
        let (molecule, compaction) = self.to_rust()?.tracked_extract(selection.to_rust());
        Ok((
            Self::from_rust(molecule),
            MoleculeCompaction::from_rust(compaction),
        ))
    }

    /// Find occurrences of this pattern in `host`.
    #[pyo3(signature = (host, *, config=None))]
    fn substructure_matches(
        &self,
        host: &Self,
        config: Option<SubstructureSearchConfig>,
    ) -> PyResult<Vec<MoleculeCorrespondence>> {
        let config = config.unwrap_or_default().to_rust();
        Ok(self
            .to_rust()?
            .substructure_matches(host.to_rust()?, config)
            .map_err(|error| PyValueError::new_err(error.to_string()))?
            .into_iter()
            .map(MoleculeCorrespondence::from_rust)
            .collect())
    }

    /// Generate an unfolded binary hashed fingerprint.
    #[pyo3(signature = (*, config))]
    fn hashed_fingerprint(&self, config: HashedFingerprintConfig) -> PyResult<HashedFeatureSet> {
        config
            .to_rust()
            .featurize(self.to_rust()?)
            .map(HashedFeatureSet::from_rust)
            .map_err(fingerprint_error)
    }

    /// Generate an unfolded counted hashed fingerprint.
    #[pyo3(signature = (*, config))]
    fn counted_hashed_fingerprint(
        &self,
        config: HashedFingerprintConfig,
    ) -> PyResult<CountedHashedFeatureSet> {
        config
            .to_rust()
            .featurize_counted(self.to_rust()?)
            .map(CountedHashedFeatureSet::from_rust)
            .map_err(fingerprint_error)
    }

    /// Generate a fixed-width pattern fingerprint.
    #[pyo3(signature = (*, config=None))]
    fn pattern_fingerprint(&self, config: Option<PatternFingerprintConfig>) -> PyResult<BitFp> {
        config
            .map_or_else(
                GraphPatternFingerprinter::new,
                PatternFingerprintConfig::to_rust,
            )
            .fingerprint(self.to_rust()?)
            .map(BitFp::from_rust)
            .map_err(fingerprint_error)
    }

    /// Generate exact canonical structural features.
    #[pyo3(signature = (*, config))]
    fn structural_fingerprint(
        &self,
        config: StructuralFingerprintConfig,
    ) -> PyResult<StructuralFeatureSet> {
        config
            .to_rust()
            .featurize(self.to_rust()?)
            .map(StructuralFeatureSet::from_rust)
            .map_err(fingerprint_error)
    }

    /// The atoms, indexed by integer position.
    #[getter]
    fn atoms(slf: Py<Self>, py: Python<'_>) -> PyResult<AtomViews> {
        AtomViews::new(slf, py)
    }

    /// The bonds, indexed by integer position.
    #[getter]
    fn bonds(slf: Py<Self>, py: Python<'_>) -> PyResult<BondViews> {
        BondViews::new(slf, py)
    }

    /// The dative bonds, indexed by integer position.
    #[getter]
    fn dative_bonds(slf: Py<Self>, py: Python<'_>) -> PyResult<DativeBondViews> {
        DativeBondViews::new(slf, py)
    }

    /// The aromatic systems, indexed by integer position.
    #[getter]
    fn aromatic_systems(slf: Py<Self>, py: Python<'_>) -> PyResult<AromaticSystemViews> {
        AromaticSystemViews::new(slf, py)
    }

    /// The multicenter bonds, indexed by integer position.
    #[getter]
    fn multicenter_bonds(slf: Py<Self>, py: Python<'_>) -> PyResult<MulticenterBondViews> {
        MulticenterBondViews::new(slf, py)
    }

    /// The noncovalent bonds, indexed by integer position.
    #[getter]
    fn noncovalent_bonds(slf: Py<Self>, py: Python<'_>) -> PyResult<NoncovalentBondViews> {
        NoncovalentBondViews::new(slf, py)
    }

    /// The stereo atoms, indexed by integer position.
    #[getter]
    fn stereo_atoms(slf: Py<Self>, py: Python<'_>) -> PyResult<StereoAtomViews> {
        StereoAtomViews::new(slf, py)
    }

    /// The stereo bonds, indexed by integer position.
    #[getter]
    fn stereo_bonds(slf: Py<Self>, py: Python<'_>) -> PyResult<StereoBondViews> {
        StereoBondViews::new(slf, py)
    }

    /// The molecule-level constraints in insertion order.
    #[getter]
    fn constraints(slf: Py<Self>, py: Python<'_>) -> PyResult<ConstraintsView> {
        ConstraintsView::new(slf, py)
    }

    fn __eq__(&self, other: &Self) -> PyResult<bool> {
        Ok(self.to_rust()? == other.to_rust()?)
    }

    pub(crate) fn __repr__(&self) -> PyResult<String> {
        let molecule = self.to_rust()?;
        // Atoms and bonds always; the other entity kinds (dative bonds, aromatic systems,
        // multicenter bonds, noncovalent bonds, stereo atoms, stereo bonds) only when present,
        // so a plain covalent molecule stays uncluttered. Names match the `from_entries` kwargs.
        let mut parts = vec![
            format!("atoms={}", molecule.atoms().count()),
            format!("bonds={}", molecule.bonds().count()),
        ];
        for (name, count) in [
            ("dative_bonds", molecule.dative_bonds().count()),
            ("aromatic_systems", molecule.aromatic_systems().count()),
            ("multicenter_bonds", molecule.multicenter_bonds().count()),
            ("noncovalent_bonds", molecule.noncovalent_bonds().count()),
            ("stereo_atoms", molecule.stereo_atoms().count()),
            ("stereo_bonds", molecule.stereo_bonds().count()),
        ] {
            if count > 0 {
                parts.push(format!("{name}={count}"));
            }
        }
        Ok(format!("Molecule({})", parts.join(", ")))
    }
}

impl Molecule {
    pub(crate) fn view_counter(&self) -> PyResult<u64> {
        self.to_rust()?;
        Ok(self.counter)
    }

    pub(crate) fn check_access(&self, expected: u64, accessor: &'static str) -> PyResult<()> {
        if self.value.is_none() || self.counter != expected {
            return Err(InvalidatedViewError::new_err(format!(
                "{accessor} was invalidated by Molecule mutation or consumption"
            )));
        }
        Ok(())
    }

    pub(crate) fn advance_counter(&mut self) -> PyResult<()> {
        self.to_rust()?;
        self.counter += 1;
        Ok(())
    }

    /// Read access to the wrapped IR molecule.
    pub(crate) fn to_rust(&self) -> PyResult<&GraphIrMolecule> {
        self.value
            .as_ref()
            .ok_or_else(|| ConsumedError::new_err("Molecule has been consumed"))
    }

    /// Mutable access to the wrapped IR molecule.
    pub(crate) fn to_rust_mut(&mut self) -> PyResult<&mut GraphIrMolecule> {
        self.value
            .as_mut()
            .ok_or_else(|| ConsumedError::new_err("Molecule has been consumed"))
    }

    pub(crate) fn take(&mut self) -> PyResult<GraphIrMolecule> {
        self.value
            .take()
            .ok_or_else(|| ConsumedError::new_err("Molecule has been consumed"))
    }

    /// Wrap a Rust molecule as a Python molecule value.
    pub(crate) fn from_rust(molecule: GraphIrMolecule) -> Self {
        Self {
            value: Some(molecule),
            counter: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use pyo3::exceptions::PyRuntimeError;
    use pyo3::types::{PyBytes, PyList};
    use rstest::{fixture, rstest};
    use umol_chem::element::Element as ChemElement;
    use umol_graph::fingerprint::{
        CountedFeatureSet as GraphCountedFeatureSet, FeatureSet as GraphFeatureSet,
        SubstructureFeaturizer as GraphSubstructureFeaturizer,
    };
    use umol_graph::ingest::ingest_smiles;
    use umol_graph::ops::valence::ValenceTable as GraphValenceTable;
    use umol_graph_core::{
        Correspondence as GraphCoreCorrespondence,
        RelevantCycleEnumerationAlgorithm as GraphCoreRelevantCycleEnumerationAlgorithm,
        SubgraphIsomorphismAlgorithm as GraphCoreSubgraphIsomorphismAlgorithm,
    };
    use umol_graph_ir::dsl::{
        AtomDsl as GraphIrAtomDsl, MoleculeMetadata as GraphIrMoleculeMetadata,
    };
    use umol_graph_ir::ir::{
        AromaticSystemForm as GraphIrAromaticSystemForm,
        AromaticSystemId as GraphIrAromaticSystemId, AtomFieldChange as GraphIrAtomFieldChange,
        AtomForm as GraphIrAtomForm, AtomHandle as GraphIrAtomHandle,
        AtomUpdate as GraphIrAtomUpdate, BondForm as GraphIrBondForm,
        Constraint as GraphIrConstraint, DativeBondForm as GraphIrDativeBondForm,
        DativeBondId as GraphIrDativeBondId, Edit as GraphIrEdit, Edits as GraphIrEdits,
        Entity as GraphIrEntity, MoleculeConstraint as GraphIrMoleculeConstraint,
        MoleculeCorrespondence as GraphIrMoleculeCorrespondence,
        MulticenterBondForm as GraphIrMulticenterBondForm,
        MulticenterBondId as GraphIrMulticenterBondId,
        NoncovalentBondForm as GraphIrNoncovalentBondForm,
        NoncovalentBondId as GraphIrNoncovalentBondId,
        NoncovalentBondKind as GraphIrNoncovalentBondKind, NumForm as GraphIrNumForm,
        SubstructureMatchAlgorithm as GraphIrSubstructureMatchAlgorithm,
        SubstructureMatchConfig as GraphIrSubstructureMatchConfig,
    };
    use umol_graph_ir::mol_dsl;

    use super::*;
    use crate::atom::AtomForm as PyAtomForm;
    use crate::error::{
        ConsumedError, InvalidStructureError, MetadataError, ParseError, TransactionError,
        UnderdeterminedError,
    };
    use crate::fingerprint::config::{
        EcfpHashScheme, PatternFingerprintConfig, RefinementRounds, StructuralFingerprintConfig,
        WlHashScheme,
    };
    use crate::ring::RingConfig;

    #[fixture]
    fn ethanol() -> Molecule {
        Molecule::from_rust(ingest_smiles("CCO").unwrap())
    }

    #[fixture]
    fn ethane() -> Molecule {
        Molecule::from_rust(ingest_smiles("CC").unwrap())
    }

    #[rstest]
    fn test_molecule_new() {
        assert_eq!(Molecule::new().to_rust().unwrap(), &GraphIrMolecule::new());
    }

    #[rstest]
    #[case::required(
        r#"{:atoms ["C"]}"#,
        None,
        mol_dsl!(r#"{:atoms ["C"]}"#)
    )]
    #[case::ground(
        r#"{:atoms ["C#h4#v0#d0#t0#a!#m!"]}"#,
        Some(MoleculeDefaults::concrete()),
        mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"]}"#)
    )]
    fn test_molecule_parse(
        #[case] text: &str,
        #[case] defaults: Option<MoleculeDefaults>,
        #[case] expected: GraphIrMolecule,
    ) {
        assert_eq!(
            Molecule::parse(text, defaults).unwrap().to_rust().unwrap(),
            &expected
        );
    }

    #[rstest]
    fn test_molecule_parse_error() {
        Python::attach(|py| {
            let error = Molecule::parse("not edn", None).unwrap_err();

            assert!(error.is_instance_of::<ParseError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "EDN parse: unexpected token 'n' at byte 0"
            );
        });
    }

    #[rstest]
    fn test_molecule_parse_with_metadata() {
        let (molecule, metadata) = Molecule::parse_with_metadata(
            r#"{:atoms [[:carbon :x]] :bonds [] :atom-aliases [:x "C"]}"#,
            None,
        )
        .unwrap();
        let metadata = metadata.to_rust();

        assert_eq!(molecule.to_rust().unwrap(), &mol_dsl!(r#"{:atoms ["C"]}"#));
        assert_eq!(
            metadata.keyword(GraphIrEntity::Atom(GraphIrAtomId(0))),
            Some("carbon")
        );
        assert_eq!(
            metadata.atom_alias("x"),
            Some(&GraphIrAtomDsl(GraphIrAtomForm::from_element(
                ChemElement::C
            )))
        );
    }

    #[rstest]
    fn test_molecule_parse_with_metadata_defaults() {
        let (molecule, metadata) = Molecule::parse_with_metadata(
            r#"{:atoms ["C#h4#v0#d0#t0#a!#m!"]}"#,
            Some(MoleculeDefaults::concrete()),
        )
        .unwrap();

        assert_eq!(
            molecule.to_rust().unwrap(),
            &mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"]}"#)
        );
        assert_eq!(
            metadata,
            MoleculeMetadata::from_rust(GraphIrMoleculeMetadata::new())
        );
    }

    #[rstest]
    #[case::required(
        mol_dsl!(r#"{:atoms ["C"]}"#),
        None,
        r#"{:atoms ["C"] :bonds []}"#
    )]
    #[case::ground(
        mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"]}"#),
        Some(MoleculeDefaults::concrete()),
        r#"{:atoms ["C#h4#v0#d0#t0#a!#m!"] :bonds []}"#
    )]
    fn test_molecule_render(
        #[case] molecule: GraphIrMolecule,
        #[case] defaults: Option<MoleculeDefaults>,
        #[case] expected: &str,
    ) {
        assert_eq!(
            Molecule::from_rust(molecule).render(defaults).unwrap(),
            expected
        );
    }

    #[rstest]
    fn test_molecule_render_with_metadata() {
        let molecule = Molecule::from_rust(mol_dsl!(r#"{:atoms ["C"]}"#));
        let mut metadata = GraphIrMoleculeMetadata::new();
        metadata
            .set_keyword(GraphIrEntity::Atom(GraphIrAtomId(0)), "carbon")
            .unwrap();
        metadata
            .add_atom_alias(
                "x",
                GraphIrAtomDsl(GraphIrAtomForm::from_element(ChemElement::C)),
            )
            .unwrap();

        assert_eq!(
            molecule
                .render_with_metadata(&MoleculeMetadata::from_rust(metadata), None)
                .unwrap(),
            r#"{:atom-aliases [:x "C"] :atoms [[:carbon :x]] :bonds []}"#
        );
    }

    #[rstest]
    fn test_molecule_render_with_metadata_error() {
        Python::attach(|py| {
            let molecule = Molecule::from_rust(mol_dsl!(r#"{:atoms ["C"]}"#));
            let mut metadata = GraphIrMoleculeMetadata::new();
            metadata
                .set_keyword(GraphIrEntity::Atom(GraphIrAtomId(1)), "outside")
                .unwrap();

            let error = molecule
                .render_with_metadata(&MoleculeMetadata::from_rust(metadata), None)
                .unwrap_err();

            assert!(error.is_instance_of::<MetadataError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "metadata entity is out of range: atom 1"
            );
        });
    }

    #[rstest]
    fn test_molecule_str() {
        let molecule = Molecule::from_rust(mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#));

        assert_eq!(molecule.__str__().unwrap(), molecule.render(None).unwrap());
    }

    #[rstest]
    fn test_molecule_from_entries() {
        Python::attach(|py| {
            let atoms = vec![
                Py::new(
                    py,
                    PyAtomForm::from_rust(GraphIrAtomForm::from_element(ChemElement::C)),
                )
                .unwrap(),
                Py::new(
                    py,
                    PyAtomForm::from_rust(GraphIrAtomForm::from_element(ChemElement::B)),
                )
                .unwrap(),
                Py::new(
                    py,
                    PyAtomForm::from_rust(GraphIrAtomForm::from_element(ChemElement::N)),
                )
                .unwrap(),
            ];
            let bonds = vec![(
                0,
                1,
                Py::new(py, BondForm::from_rust(GraphIrBondForm::from_order(1))).unwrap(),
            )];
            let dative = vec![(
                vec![2],
                1,
                Py::new(
                    py,
                    DativeBondForm::from_rust(GraphIrDativeBondForm::from_order(1)),
                )
                .unwrap(),
            )];
            let aromatic = vec![(
                vec![0, 1, 2],
                Py::new(
                    py,
                    AromaticSystemForm::from_rust(GraphIrAromaticSystemForm::from_electrons(vec![
                        1, 1, 1,
                    ])),
                )
                .unwrap(),
            )];
            let multicenter = vec![(
                vec![0, 1, 2],
                Py::new(
                    py,
                    MulticenterBondForm::from_rust(GraphIrMulticenterBondForm::from_electrons(
                        vec![1, 1, 1],
                    )),
                )
                .unwrap(),
            )];
            let noncovalent = vec![(
                [0, 2],
                Py::new(
                    py,
                    NoncovalentBondForm::from_rust(GraphIrNoncovalentBondForm::from_kind(
                        GraphIrNoncovalentBondKind::HydrogenBond,
                    )),
                )
                .unwrap(),
            )];
            let constraint = GraphIrConstraint::Molecule(GraphIrMoleculeConstraint::Connected {
                atoms: Some(vec![GraphIrAtomId(0), GraphIrAtomId(2)]),
            });
            let constraints = vec![Constraint::from_rust(py, &constraint).unwrap()];
            let molecule = Molecule::from_entries(
                py,
                atoms,
                bonds,
                dative,
                aromatic,
                multicenter,
                noncovalent,
                Vec::new(),
                Vec::new(),
                constraints,
            )
            .unwrap();
            assert_eq!(molecule.to_rust().unwrap().atoms().count(), 3);
            assert_eq!(molecule.to_rust().unwrap().bonds().count(), 1);
            let dative_bonds = molecule.to_rust().unwrap().dative_bonds();
            assert_eq!(dative_bonds.count(), 1);
            let dative_view = dative_bonds.get(GraphIrDativeBondId(0)).unwrap();
            assert_eq!(dative_view.acceptor_id(), GraphIrAtomId(1));
            assert_eq!(
                dative_view.donor_ids().collect::<Vec<_>>(),
                vec![GraphIrAtomId(2)]
            );
            let aromatic_systems = molecule.to_rust().unwrap().aromatic_systems();
            assert_eq!(aromatic_systems.count(), 1);
            let aromatic_view = aromatic_systems.get(GraphIrAromaticSystemId(0)).unwrap();
            assert_eq!(
                aromatic_view.atom_ids().collect::<Vec<_>>(),
                vec![GraphIrAtomId(0), GraphIrAtomId(1), GraphIrAtomId(2)]
            );
            let multicenter_bonds = molecule.to_rust().unwrap().multicenter_bonds();
            assert_eq!(multicenter_bonds.count(), 1);
            let multicenter_view = multicenter_bonds.get(GraphIrMulticenterBondId(0)).unwrap();
            assert_eq!(
                multicenter_view.atom_ids().collect::<Vec<_>>(),
                vec![GraphIrAtomId(0), GraphIrAtomId(1), GraphIrAtomId(2)]
            );
            let noncovalent_bonds = molecule.to_rust().unwrap().noncovalent_bonds();
            assert_eq!(noncovalent_bonds.count(), 1);
            let noncovalent_view = noncovalent_bonds.get(GraphIrNoncovalentBondId(0)).unwrap();
            assert_eq!(
                noncovalent_view.atom_ids(),
                [GraphIrAtomId(0), GraphIrAtomId(2)]
            );
            assert_eq!(
                molecule.to_rust().unwrap().constraints().as_slice(),
                &[constraint]
            );
        });
    }

    #[rstest]
    #[case::defaults(None, None, None)]
    #[case::explicit(
        Some(SmilesIoConfig::from_rust(&IoSmilesIoConfig::opensmiles())),
        Some(ChemistryModel::from_rust(&GraphChemistryModel {
            valence: GraphValenceModel::smiles(),
            ..GraphChemistryModel::default()
        })),
        Some(ResolveConfig::from_rust(GraphResolveConfig::default())),
    )]
    fn test_molecule_from_smiles(
        #[case] io_config: Option<SmilesIoConfig>,
        #[case] chemistry_model: Option<ChemistryModel>,
        #[case] resolve_config: Option<ResolveConfig>,
    ) {
        assert_eq!(
            Molecule::from_smiles("C", io_config, chemistry_model, resolve_config)
                .unwrap()
                .to_rust()
                .unwrap(),
            &mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s"]}"#)
        );
    }

    #[rstest]
    #[case::syntax(" C", "ParseError", "Leading whitespace")]
    #[case::incomplete_stereo(
        "C[S@]C",
        "ContradictionError",
        "stereo inconsistency: stereo atom StereoAtomId(0) cannot be realized"
    )]
    #[case::underdetermined("*", "UnderdeterminedError", "resolution underdetermined")]
    fn test_molecule_from_smiles_error(
        #[case] source: &str,
        #[case] expected_type: &str,
        #[case] expected_message: &str,
    ) {
        Python::attach(|py| {
            let error = Molecule::from_smiles(source, None, None, None).unwrap_err();
            assert_eq!(error.get_type(py).name().unwrap(), expected_type);
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                expected_message
            );
        });
    }

    #[rstest]
    fn test_molecule_copy() {
        let initial = mol_dsl!(r#"{:atoms ["C"]}"#);
        let mut molecule = Molecule::from_rust(initial.clone());
        molecule.advance_counter().unwrap();
        let mut copied = molecule.copy().unwrap();
        assert_eq!(copied.view_counter().unwrap(), 0);
        assert!(molecule.__eq__(&copied).unwrap());
        copied
            .to_rust_mut()
            .unwrap()
            .atom_mut(GraphIrAtomId(0))
            .attributes_mut()
            .charge = GraphIrNumForm::Lit(1);
        assert_eq!(molecule.take().unwrap(), initial);
        assert_eq!(copied.to_rust().unwrap(), &mol_dsl!(r#"{:atoms ["C#c+"]}"#));
    }

    #[rstest]
    fn test_molecule_edit() {
        Python::attach(|py| {
            let expected = mol_dsl!(r#"{:atoms ["N#h3"]}"#);
            let editor =
                Py::new(py, Molecule::from_rust(expected.clone()).edit().unwrap()).unwrap();
            let published = editor
                .bind(py)
                .call_method0("finish")
                .unwrap()
                .extract::<Py<Molecule>>()
                .unwrap();

            assert_eq!(published.bind(py).borrow().to_rust().unwrap(), &expected);
        });
    }

    #[rstest]
    fn test_molecule_apply() {
        let initial = mol_dsl!(r#"{:atoms ["N#h3"]}"#);
        let mut rust_edits = GraphIrEdits::new();
        rust_edits.update_atom(
            GraphIrAtomHandle::Id(GraphIrAtomId(0)),
            initial.atom(GraphIrAtomId(0)).attributes(),
            &GraphIrAtomUpdate {
                implicit_hydrogens: Some(GraphIrNumForm::Lit(2)),
                ..Default::default()
            },
        );
        let methyl = rust_edits
            .add_atom(GraphIrAtomForm::from_element(ChemElement::C).with_implicit_hydrogens(3_i64));
        rust_edits.add_bond(
            GraphIrAtomHandle::Id(GraphIrAtomId(0)),
            methyl,
            GraphIrBondForm::from_order(1),
        );
        let mut molecule = Molecule::from_rust(initial);

        Python::attach(|py| {
            let edits = Py::new(py, Edits::from_rust(rust_edits)).unwrap();

            let result = molecule.apply(py, edits).unwrap();

            assert_eq!(
                result.to_rust().unwrap(),
                &mol_dsl!(r#"{:atoms ["N#h2" "C#h3"] :bonds [[0 1 "1"]]}"#)
            );
            assert!(molecule
                .to_rust()
                .unwrap_err()
                .is_instance_of::<ConsumedError>(py));
        });
    }

    #[rstest]
    fn test_molecule_apply_error() {
        let initial = mol_dsl!(r#"{:atoms ["C"]}"#);
        let mut molecule = Molecule::from_rust(initial);
        let mut rust_edits = GraphIrEdits::new();
        rust_edits.add_atom(GraphIrAtomForm::from_element(ChemElement::N));
        rust_edits.push(GraphIrEdit::ModifyAtomField {
            id: GraphIrAtomHandle::Id(GraphIrAtomId(7)),
            change: GraphIrAtomFieldChange::Charge {
                old: GraphIrNumForm::Lit(0),
                new: GraphIrNumForm::Lit(1),
            },
        });

        Python::attach(|py| {
            let edits = Py::new(py, Edits::from_rust(rust_edits)).unwrap();

            let error = molecule.apply(py, edits).unwrap_err();

            assert!(error.is_instance_of::<TransactionError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "atom handle 7 is out of range for 1 entries"
            );
            assert!(molecule
                .to_rust()
                .unwrap_err()
                .is_instance_of::<ConsumedError>(py));
        });
    }

    #[rstest]
    fn test_molecule_apply_publication_error() {
        let initial = mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#);
        let mut molecule = Molecule::from_rust(initial);
        let mut rust_edits = GraphIrEdits::new();
        rust_edits.add_bond(
            GraphIrAtomHandle::Id(GraphIrAtomId(0)),
            GraphIrAtomHandle::Id(GraphIrAtomId(1)),
            GraphIrBondForm::from_order(1),
        );

        Python::attach(|py| {
            let edits = Py::new(py, Edits::from_rust(rust_edits)).unwrap();

            let error = molecule.apply(py, edits).unwrap_err();

            assert!(error.is_instance_of::<InvalidStructureError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "bond: parallel bonds on atoms [AtomId(0), AtomId(1)]"
            );
            assert!(molecule
                .to_rust()
                .unwrap_err()
                .is_instance_of::<ConsumedError>(py));
        });
    }

    #[rstest]
    #[case::untracked(false)]
    #[case::tracked(true)]
    fn test_molecule_transact_borrow_error(#[case] tracked: bool) {
        Python::attach(|py| {
            let initial = mol_dsl!(r#"{:atoms ["C"]}"#);
            let molecule = Py::new(py, Molecule::from_rust(initial.clone())).unwrap();
            let edits = Py::new(py, Edits::from_rust(GraphIrEdits::new())).unwrap();
            let borrowed = molecule.borrow(py);
            let batches = PyList::new(py, [edits.clone_ref(py)]).unwrap();

            let result = if tracked {
                Molecule::tracked_transact(molecule.clone_ref(py), py, batches.as_any()).map(drop)
            } else {
                Molecule::transact(molecule.clone_ref(py), py, batches.as_any())
            };

            let error = result.unwrap_err();
            assert!(error.is_instance_of::<PyRuntimeError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "Already borrowed"
            );
            assert_eq!(borrowed.to_rust().unwrap(), &initial);
            assert_eq!(borrowed.view_counter().unwrap(), 0);
            assert!(edits
                .borrow(py)
                .to_rust()
                .unwrap_err()
                .is_instance_of::<ConsumedError>(py));
        });
    }

    #[rstest]
    #[case::untracked(false)]
    #[case::tracked(true)]
    fn test_molecule_transact_batch_borrow_error(#[case] tracked: bool) {
        Python::attach(|py| {
            let initial = mol_dsl!(r#"{:atoms ["C"]}"#);
            let molecule = Py::new(py, Molecule::from_rust(initial.clone())).unwrap();
            let first = Py::new(py, Edits::from_rust(GraphIrEdits::new())).unwrap();
            let second = Py::new(py, Edits::from_rust(GraphIrEdits::new())).unwrap();
            let third = Py::new(py, Edits::from_rust(GraphIrEdits::new())).unwrap();
            let borrowed = second.borrow(py);
            let batches = PyList::new(
                py,
                [
                    first.clone_ref(py),
                    second.clone_ref(py),
                    third.clone_ref(py),
                ],
            )
            .unwrap();

            let result = if tracked {
                Molecule::tracked_transact(molecule.clone_ref(py), py, batches.as_any()).map(drop)
            } else {
                Molecule::transact(molecule.clone_ref(py), py, batches.as_any())
            };

            let error = result.unwrap_err();
            assert!(error.is_instance_of::<PyRuntimeError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "Already borrowed"
            );
            assert_eq!(molecule.borrow(py).to_rust().unwrap(), &initial);
            assert_eq!(molecule.borrow(py).view_counter().unwrap(), 0);
            assert!(first
                .borrow(py)
                .to_rust()
                .unwrap_err()
                .is_instance_of::<ConsumedError>(py));
            assert_eq!(borrowed.to_rust().unwrap(), &GraphIrEdits::new());
            assert_eq!(third.borrow(py).to_rust().unwrap(), &GraphIrEdits::new());
        });
    }

    #[rstest]
    fn test_molecule_substructure_matches() {
        let pattern = Molecule::from_rust(mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#));
        let host = Molecule::from_rust(mol_dsl!(
            r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#
        ));
        let pattern_before = pattern.to_rust().unwrap().clone();
        let host_before = host.to_rust().unwrap().clone();
        let expected = vec![
            MoleculeCorrespondence::from_rust(
                GraphIrMoleculeCorrespondence::induce(
                    pattern.to_rust().unwrap(),
                    host.to_rust().unwrap(),
                    GraphCoreCorrespondence::new(
                        vec![
                            (GraphIrAtomId(0), GraphIrAtomId(0)),
                            (GraphIrAtomId(1), GraphIrAtomId(1)),
                        ],
                        2,
                        3,
                    )
                    .expect("correspondence producer preserves partial-bijection invariants"),
                )
                .expect("the atom correspondence describes the molecule pair"),
            ),
            MoleculeCorrespondence::from_rust(
                GraphIrMoleculeCorrespondence::induce(
                    pattern.to_rust().unwrap(),
                    host.to_rust().unwrap(),
                    GraphCoreCorrespondence::new(
                        vec![
                            (GraphIrAtomId(0), GraphIrAtomId(1)),
                            (GraphIrAtomId(1), GraphIrAtomId(0)),
                        ],
                        2,
                        3,
                    )
                    .expect("correspondence producer preserves partial-bijection invariants"),
                )
                .expect("the atom correspondence describes the molecule pair"),
            ),
        ];

        assert_eq!(pattern.substructure_matches(&host, None).unwrap(), expected);
        assert_eq!(pattern.to_rust().unwrap(), &pattern_before);
        assert_eq!(host.to_rust().unwrap(), &host_before);
    }

    #[rstest]
    fn test_molecule_substructure_matches_overlay() {
        let pattern = Molecule::from_rust(mol_dsl!(
            r#"{
                :atoms ["N" "B"]
                :bonds []
                :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]
            }"#
        ));
        let host = Molecule::from_rust(mol_dsl!(
            r#"{
                :atoms ["N" "B" "C"]
                :bonds []
                :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]
            }"#
        ));
        let expected = vec![MoleculeCorrespondence::from_rust(
            GraphIrMoleculeCorrespondence::induce(
                pattern.to_rust().unwrap(),
                host.to_rust().unwrap(),
                GraphCoreCorrespondence::new(
                    vec![
                        (GraphIrAtomId(0), GraphIrAtomId(0)),
                        (GraphIrAtomId(1), GraphIrAtomId(1)),
                    ],
                    2,
                    3,
                )
                .expect("correspondence producer preserves partial-bijection invariants"),
            )
            .expect("the atom correspondence describes the molecule pair"),
        )];
        let config = SubstructureSearchConfig::from_rust(GraphIrSubstructureMatchConfig {
            match_algorithm: GraphIrSubstructureMatchAlgorithm::Incidence,
            subgraph_isomorphism_algorithm: GraphCoreSubgraphIsomorphismAlgorithm::Ullmann,
            relevant_cycle_algorithm: GraphCoreRelevantCycleEnumerationAlgorithm::Vismara,
        });

        assert_eq!(
            pattern.substructure_matches(&host, Some(config)).unwrap(),
            expected
        );
    }

    #[rstest]
    fn test_molecule_substructure_matches_empty() {
        let pattern = Molecule::from_rust(mol_dsl!(r#"{:atoms ["O"] :bonds []}"#));
        let host = Molecule::from_rust(mol_dsl!(r#"{:atoms ["C"] :bonds []}"#));
        let config = SubstructureSearchConfig::from_rust(GraphIrSubstructureMatchConfig {
            match_algorithm: GraphIrSubstructureMatchAlgorithm::GraphAndOverlays,
            subgraph_isomorphism_algorithm: GraphCoreSubgraphIsomorphismAlgorithm::Vf2,
            relevant_cycle_algorithm: GraphCoreRelevantCycleEnumerationAlgorithm::Vismara,
        });

        assert_eq!(
            pattern.substructure_matches(&host, Some(config)).unwrap(),
            Vec::new()
        );
    }

    #[rstest]
    #[case::morgan_default(
        HashedFingerprintConfig::Morgan {
            radius: 2,
            ring_config: RingConfig::default(),
        },
        &[
            864662311,
            1535166686,
            2245384272,
            2246728737,
            3542456614,
            4018048386,
        ]
    )]
    #[case::morgan_explicit(
        HashedFingerprintConfig::Morgan {
            radius: 0,
            ring_config: RingConfig::default(),
        },
        &[864662311, 2245384272, 2246728737]
    )]
    #[case::ecfp_default(
        HashedFingerprintConfig::Ecfp {
            radius: 2,
            hashing_scheme: EcfpHashScheme::Xxh3Width64V1(),
            ring_config: RingConfig::default(),
        },
        &[
            63839236075656913,
            1189585227353469813,
            3822471596818936039,
            13652293261850732425,
            15001976065402722634,
            16149328945726899460,
        ]
    )]
    #[case::ecfp_explicit(
        HashedFingerprintConfig::Ecfp {
            radius: 0,
            hashing_scheme: EcfpHashScheme::Xxh3Width64V1(),
            ring_config: RingConfig::default(),
        },
        &[
            1189585227353469813,
            3822471596818936039,
            16149328945726899460,
        ]
    )]
    #[case::wl_default_scheme(
        HashedFingerprintConfig::Wl {
            rounds: RefinementRounds::Fixed { rounds: 3 },
            hashing_scheme: WlHashScheme::Xxh3SortedWidth64V1(),
        },
        &[
            2520347590860685079,
            3352603313223549703,
            4152249898001161146,
            5715207763479934940,
            5807737097854608645,
            7542810387455301591,
            11457795998246593156,
            11986000156817227245,
            12895020514073294021,
            13932567567828606490,
            17305796300852423160,
            17417400371411086222,
        ]
    )]
    #[case::wl_explicit_rounds(
        HashedFingerprintConfig::Wl {
            rounds: RefinementRounds::Fixed { rounds: 1 },
            hashing_scheme: WlHashScheme::Xxh3SortedWidth64V1(),
        },
        &[
            5715207763479934940,
            5807737097854608645,
            7542810387455301591,
            11457795998246593156,
            12895020514073294021,
            17417400371411086222,
        ]
    )]
    fn test_molecule_hashed_fingerprint(
        ethanol: Molecule,
        #[case] config: HashedFingerprintConfig,
        #[case] expected_ids: &[u64],
    ) {
        let fingerprint = ethanol.hashed_fingerprint(config).unwrap();
        assert_eq!(
            fingerprint,
            HashedFeatureSet::from_rust(GraphFeatureSet::from_features(
                expected_ids.iter().copied()
            ))
        );

        Python::attach(|py| {
            let fingerprint = Py::new(py, fingerprint).unwrap();
            let fingerprint = fingerprint.bind(py).as_any();
            fingerprint
                .getattr("ids")
                .unwrap()
                .cast::<PyList>()
                .unwrap()
                .append(9u64)
                .unwrap();
            assert_eq!(
                fingerprint
                    .getattr("ids")
                    .unwrap()
                    .extract::<Vec<u64>>()
                    .unwrap(),
                expected_ids
            );
        });
    }

    #[rstest]
    #[case::morgan(HashedFingerprintConfig::Morgan {
        radius: 2,
        ring_config: RingConfig::default(),
    })]
    #[case::ecfp(HashedFingerprintConfig::Ecfp {
        radius: 2,
        hashing_scheme: EcfpHashScheme::Xxh3Width64V1(),
        ring_config: RingConfig::default(),
    })]
    #[case::wl(HashedFingerprintConfig::Wl {
        rounds: RefinementRounds::Fixed { rounds: 3 },
        hashing_scheme: WlHashScheme::Xxh3SortedWidth64V1(),
    })]
    fn test_molecule_hashed_fingerprint_error(#[case] config: HashedFingerprintConfig) {
        Python::attach(|py| {
            let molecule = Molecule::from_rust(mol_dsl!(r#"{:atoms ["C"] :bonds []}"#));
            let error = molecule.hashed_fingerprint(config).unwrap_err();
            assert!(error.is_instance_of::<UnderdeterminedError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "fingerprint requires a determined molecule"
            );
        });
    }

    #[rstest]
    #[case::morgan(
        HashedFingerprintConfig::Morgan {
            radius: 2,
            ring_config: RingConfig::default(),
        },
        &[(2246728737, 2), (3545175291, 1)]
    )]
    #[case::ecfp(
        HashedFingerprintConfig::Ecfp {
            radius: 2,
            hashing_scheme: EcfpHashScheme::Xxh3Width64V1(),
            ring_config: RingConfig::default(),
        },
        &[(5513743581508886362, 1), (16149328945726899460, 2)]
    )]
    #[case::wl(
        HashedFingerprintConfig::Wl {
            rounds: RefinementRounds::Fixed { rounds: 3 },
            hashing_scheme: WlHashScheme::Xxh3SortedWidth64V1(),
        },
        &[
            (2659163409134283895, 2),
            (7542810387455301591, 2),
            (9541344068636876323, 2),
            (12512207080905326651, 2),
        ]
    )]
    fn test_molecule_counted_hashed_fingerprint(
        ethane: Molecule,
        #[case] config: HashedFingerprintConfig,
        #[case] expected_entries: &[(u64, u32)],
    ) {
        let fingerprint = ethane.counted_hashed_fingerprint(config).unwrap();
        assert_eq!(
            fingerprint,
            CountedHashedFeatureSet::from_rust(GraphCountedFeatureSet::from_counts(
                expected_entries.iter().copied()
            ))
        );

        Python::attach(|py| {
            let fingerprint = Py::new(py, fingerprint).unwrap();
            let fingerprint = fingerprint.bind(py).as_any();
            fingerprint
                .getattr("entries")
                .unwrap()
                .cast::<PyList>()
                .unwrap()
                .append((9u64, 3u32))
                .unwrap();
            assert_eq!(
                fingerprint
                    .getattr("entries")
                    .unwrap()
                    .extract::<Vec<(u64, u32)>>()
                    .unwrap(),
                expected_entries
            );
        });
    }

    #[rstest]
    #[case::morgan(HashedFingerprintConfig::Morgan {
        radius: 2,
        ring_config: RingConfig::default(),
    })]
    #[case::ecfp(HashedFingerprintConfig::Ecfp {
        radius: 2,
        hashing_scheme: EcfpHashScheme::Xxh3Width64V1(),
        ring_config: RingConfig::default(),
    })]
    #[case::wl(HashedFingerprintConfig::Wl {
        rounds: RefinementRounds::Fixed { rounds: 3 },
        hashing_scheme: WlHashScheme::Xxh3SortedWidth64V1(),
    })]
    fn test_molecule_counted_hashed_fingerprint_error(#[case] config: HashedFingerprintConfig) {
        Python::attach(|py| {
            let molecule = Molecule::from_rust(mol_dsl!(r#"{:atoms ["C"] :bonds []}"#));
            let error = molecule.counted_hashed_fingerprint(config).unwrap_err();
            assert!(error.is_instance_of::<UnderdeterminedError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "fingerprint requires a determined molecule"
            );
        });
    }

    #[rstest]
    #[case::omitted(
        None,
        2048,
        &[54, 173, 217, 429, 622, 759, 778, 874, 946, 967, 1022, 1033, 1061, 1236, 1289, 1295]
    )]
    #[case::default(
        Some(PatternFingerprintConfig::from_rust(GraphPatternFingerprinter {
            width: 2048,
            ..GraphPatternFingerprinter::new()
        })),
        2048,
        &[54, 173, 217, 429, 622, 759, 778, 874, 946, 967, 1022, 1033, 1061, 1236, 1289, 1295]
    )]
    #[case::custom(
        Some(PatternFingerprintConfig::from_rust(GraphPatternFingerprinter {
            width: 64,
            ..GraphPatternFingerprinter::new()
        })),
        64,
        &[7, 9, 10, 15, 20, 25, 37, 42, 45, 46, 50, 54, 55, 62]
    )]
    fn test_molecule_pattern_fingerprint(
        ethanol: Molecule,
        #[case] config: Option<PatternFingerprintConfig>,
        #[case] width: usize,
        #[case] expected_bits: &[u64],
    ) {
        let fingerprint = ethanol.pattern_fingerprint(config).unwrap();
        let expected = GraphFeatureSet::from_features(expected_bits.iter().copied())
            .fold(width)
            .unwrap();
        assert_eq!(fingerprint, BitFp::from_rust(expected));
    }

    #[rstest]
    #[case::omitted(None)]
    #[case::default(Some(PatternFingerprintConfig::from_rust(GraphPatternFingerprinter {
        width: 2048,
        ..GraphPatternFingerprinter::new()
    })))]
    #[case::custom(Some(PatternFingerprintConfig::from_rust(GraphPatternFingerprinter {
        width: 64,
        ..GraphPatternFingerprinter::new()
    })))]
    fn test_molecule_pattern_fingerprint_error(#[case] config: Option<PatternFingerprintConfig>) {
        Python::attach(|py| {
            let molecule = Molecule::from_rust(mol_dsl!(r#"{:atoms ["C"] :bonds []}"#));
            let error = molecule.pattern_fingerprint(config).unwrap_err();
            assert!(error.is_instance_of::<UnderdeterminedError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "fingerprint requires a determined molecule"
            );
        });
    }

    #[rstest]
    #[case::atoms(
        StructuralFingerprintConfig::from_rust(GraphSubstructureFeaturizer::new(0)),
        vec![
            vec![1, 0, 0, 0, 5, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0],
            vec![1, 0, 0, 0, 5, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0],
        ]
    )]
    #[case::bounded(
        StructuralFingerprintConfig::from_rust(GraphSubstructureFeaturizer::new(2)),
        vec![
            vec![1, 0, 0, 0, 5, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0],
            vec![1, 0, 0, 0, 5, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0],
            vec![
                3, 0, 0, 0, 5, 0, 0, 0, 0, 6, 0, 0, 0, 5, 0, 0, 0, 0, 6, 0, 0, 0, 3, 0, 0, 0, 1, 1,
                0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0,
            ],
            vec![
                3, 0, 0, 0, 5, 0, 0, 0, 0, 6, 0, 0, 0, 5, 0, 0, 0, 0, 8, 0, 0, 0, 3, 0, 0, 0, 1, 1,
                0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0,
            ],
            vec![
                5, 0, 0, 0, 5, 0, 0, 0, 0, 6, 0, 0, 0, 5, 0, 0, 0, 0, 6, 0, 0, 0, 5, 0, 0, 0, 0, 8,
                0, 0, 0, 3, 0, 0, 0, 1, 1, 0, 3, 0, 0, 0, 1, 1, 0, 4, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0,
                0, 1, 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 4, 0, 0, 0, 2, 0, 0, 0, 4, 0, 0, 0,
            ],
        ]
    )]
    fn test_molecule_structural_fingerprint(
        ethanol: Molecule,
        #[case] config: StructuralFingerprintConfig,
        #[case] expected_keys: Vec<Vec<u8>>,
    ) {
        let fingerprint = ethanol.structural_fingerprint(config).unwrap();
        assert_eq!(
            fingerprint,
            StructuralFeatureSet::from_rust(GraphFeatureSet::from_features(
                expected_keys.iter().cloned()
            ))
        );

        Python::attach(|py| {
            let fingerprint = Py::new(py, fingerprint).unwrap();
            let fingerprint = fingerprint.bind(py).as_any();
            fingerprint
                .getattr("keys")
                .unwrap()
                .cast::<PyList>()
                .unwrap()
                .append(PyBytes::new(py, b"detached"))
                .unwrap();
            assert_eq!(
                fingerprint
                    .getattr("keys")
                    .unwrap()
                    .extract::<Vec<Vec<u8>>>()
                    .unwrap(),
                expected_keys
            );
        });
    }

    #[rstest]
    #[case::bounded(StructuralFingerprintConfig::from_rust(GraphSubstructureFeaturizer::new(2)))]
    fn test_molecule_structural_fingerprint_error(#[case] config: StructuralFingerprintConfig) {
        Python::attach(|py| {
            let molecule = Molecule::from_rust(mol_dsl!(r#"{:atoms ["C"] :bonds []}"#));
            let error = molecule.structural_fingerprint(config).unwrap_err();
            assert!(error.is_instance_of::<UnderdeterminedError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "fingerprint requires a determined molecule"
            );
        });
    }

    #[rstest]
    #[case(vec![], 0)]
    #[case(vec![ChemElement::C], 1)]
    #[case(vec![ChemElement::C, ChemElement::O], 2)]
    fn test_molecule_atoms(#[case] elements: Vec<ChemElement>, #[case] expected: usize) {
        let atoms = elements
            .into_iter()
            .map(GraphIrAtomForm::from_element)
            .collect();
        let molecule = Molecule::from_rust(GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
            atoms,
            ..Default::default()
        }));
        assert_eq!(molecule.to_rust().unwrap().atoms().count(), expected);
    }

    #[rstest]
    fn test_molecule_constraints() {
        Python::attach(|py| {
            let molecule = Py::new(py, Molecule::new()).unwrap();
            let view = Molecule::constraints(molecule, py).unwrap();
            assert_eq!(
                view.read(py, |constraints| Ok(constraints.to_vec()))
                    .unwrap(),
                Vec::new()
            );
        });
    }

    #[rstest]
    fn test_molecule_eq() {
        assert!(Molecule::new().__eq__(&Molecule::new()).unwrap());
        let carbon = Molecule::from_rust(GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
            atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
            ..Default::default()
        }));
        assert!(!Molecule::new().__eq__(&carbon).unwrap());
    }

    #[rstest]
    fn test_molecule_eq_error() {
        Python::attach(|py| {
            let mut consumed = Molecule::new();
            consumed.take().unwrap();
            let available = Molecule::new();
            assert!(consumed
                .__eq__(&available)
                .unwrap_err()
                .is_instance_of::<ConsumedError>(py));
            assert!(available
                .__eq__(&consumed)
                .unwrap_err()
                .is_instance_of::<ConsumedError>(py));
            assert!(consumed
                .__eq__(&consumed)
                .unwrap_err()
                .is_instance_of::<ConsumedError>(py));
            assert!(consumed
                .copy()
                .unwrap_err()
                .is_instance_of::<ConsumedError>(py));
        });
    }

    #[rstest]
    #[case::empty(Molecule::new(), "Molecule(atoms=0, bonds=0)")]
    #[case::noncovalent(
        Molecule::from_rust(GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
            atoms: vec![
                GraphIrAtomForm::from_element(ChemElement::O),
                GraphIrAtomForm::from_element(ChemElement::O),
            ],
            noncovalent: vec![([GraphIrAtomId(0), GraphIrAtomId(1)],
                GraphIrNoncovalentBondForm::from_kind(GraphIrNoncovalentBondKind::HydrogenBond),
            )],
            ..Default::default()
        })),
        "Molecule(atoms=2, bonds=0, noncovalent_bonds=1)"
    )]
    fn test_molecule_repr(#[case] molecule: Molecule, #[case] expected: &str) {
        assert_eq!(molecule.__repr__().unwrap(), expected);
    }

    #[rstest]
    fn test_molecule_resolve() {
        Python::attach(|py| {
            let molecule = Molecule::from_rust(mol_dsl!(r#"{:atoms ["C#c0"]}"#));
            let model = ChemistryModel::from_rust(&GraphChemistryModel {
                valence: GraphValenceModel::smiles(),
                ..GraphChemistryModel::default()
            });

            let solution = molecule
                .resolve(
                    py,
                    Some(model),
                    Some(ResolveConfig::from_rust(GraphResolveConfig {
                        isotope: GraphIsotopePolicy::Natural,
                        ..Default::default()
                    })),
                )
                .unwrap();

            let Solution::Determined {
                molecule: resolved,
                report,
            } = solution
            else {
                panic!("expected Determined");
            };
            assert_eq!(
                resolved.borrow(py).to_rust().unwrap(),
                &mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s"]}"#)
            );
            assert_eq!(report.tie_breaks(), vec![0]);
            assert_eq!(
                molecule.to_rust().unwrap(),
                &mol_dsl!(r#"{:atoms ["C#c0"]}"#)
            );
        });
    }

    #[rstest]
    fn test_molecule_resolve_underdetermined() {
        Python::attach(|py| {
            let molecule = Molecule::from_rust(mol_dsl!(r#"{:atoms ["C#c0"]}"#));
            let model = ChemistryModel::from_rust(&GraphChemistryModel {
                valence: GraphValenceModel::counts(Cow::Borrowed(
                    GraphValenceTable::default_table(),
                )),
                ..GraphChemistryModel::default()
            });

            let solution = molecule.resolve(py, Some(model), None).unwrap();

            let Solution::Underdetermined { report } = solution else {
                panic!("expected Underdetermined");
            };
            assert_eq!(report.unresolved().get(0).unwrap().len(), 5);
            assert_eq!(
                molecule.to_rust().unwrap(),
                &mol_dsl!(r#"{:atoms ["C#c0"]}"#)
            );
        });
    }

    #[rstest]
    fn test_molecule_resolve_contradiction() {
        Python::attach(|py| {
            let molecule = Molecule::from_rust(mol_dsl!(r#"{:atoms ["C#c0#h5"]}"#));
            let model = ChemistryModel::from_rust(&GraphChemistryModel {
                valence: GraphValenceModel::smiles(),
                ..GraphChemistryModel::default()
            });

            let solution = molecule.resolve(py, Some(model), None).unwrap();

            let Solution::Contradictory { contradiction } = solution else {
                panic!("expected Contradictory");
            };
            assert_eq!(contradiction.__str__(), "no matching valence state");
            assert_eq!(
                molecule.to_rust().unwrap(),
                &mol_dsl!(r#"{:atoms ["C#c0#h5"]}"#)
            );
        });
    }

    #[rstest]
    fn test_molecule_resolve_default_model() {
        Python::attach(|py| {
            let molecule = Molecule::from_rust(mol_dsl!(r#"{:atoms ["C"]}"#));

            let solution = molecule.resolve(py, None, None).unwrap();

            let Solution::Underdetermined { report } = solution else {
                panic!("expected Underdetermined under the default model");
            };
            // The charge-open atom takes the registry's charge-less lookup: every
            // carbon row is a candidate.
            assert_eq!(report.unresolved().get(0).unwrap().len(), 9);
        });
    }

    #[rstest]
    #[case::empty(GraphIrMolecule::new())]
    #[case::atom(mol_dsl!(r#"{:atoms ["C"]}"#))]
    fn test_molecule_take(#[case] initial: GraphIrMolecule) {
        Python::attach(|py| {
            let mut molecule = Molecule::from_rust(initial.clone());
            assert_eq!(molecule.take().unwrap(), initial);
            assert!(molecule
                .to_rust()
                .unwrap_err()
                .is_instance_of::<ConsumedError>(py));
            assert!(molecule
                .to_rust_mut()
                .unwrap_err()
                .is_instance_of::<ConsumedError>(py));
            assert!(molecule
                .take()
                .unwrap_err()
                .is_instance_of::<ConsumedError>(py));
            assert!(molecule
                .view_counter()
                .unwrap_err()
                .is_instance_of::<ConsumedError>(py));
            assert!(molecule
                .advance_counter()
                .unwrap_err()
                .is_instance_of::<ConsumedError>(py));
            assert!(molecule
                .check_access(0, "AtomView")
                .unwrap_err()
                .is_instance_of::<InvalidatedViewError>(py));
        });
    }
}
