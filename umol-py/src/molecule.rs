//! `MoleculeAst` — a molecule (owned graph-AST root), wrapping
//! `umol_graph_ir::ir::MoleculeAst`.

use std::str::FromStr;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use umol_graph::fingerprint::PatternFingerprinter as GraphPatternFingerprinter;
use umol_graph::ingest::ingest_smiles_with;
use umol_graph::ops::model::ChemistryModel as GraphChemistryModel;
use umol_graph::ops::resolve::ResolveConfig as GraphResolveConfig;
use umol_graph_ir::dsl::MoleculeDsl as GraphIrMoleculeDsl;
use umol_graph_ir::ir::{
    AtomId as GraphIrAtomId, BondId as GraphIrBondId, FromAst, IntoAst,
    MoleculeAst as GraphIrMoleculeAst, MoleculeEntries as GraphIrMoleculeEntries,
};
use umol_io::smiles::SmilesIoConfig as IoSmilesIoConfig;

use crate::aromatic::{AromaticSystemAst, AromaticSystemViews};
use crate::atom::{AtomAst, AtomViews};
use crate::bond::{BondAst, BondViews};
use crate::constraint::molecule::{Constraint, ConstraintsLike, ConstraintsView};
use crate::correspondence::MoleculeCorrespondence;
use crate::dative::{DativeBondAst, DativeBondViews};
use crate::defaults::MoleculeDefaults;
use crate::edit::Edits;
use crate::error::{
    fingerprint_error, metadata_error, parse_error, smiles_input_error, transaction_error,
};
use crate::fingerprint::config::{
    HashedFingerprintConfig, PatternFingerprintConfig, StructuralFingerprintConfig,
};
use crate::fingerprint::value::{
    BitFp, CountedHashedFeatureSet, HashedFeatureSet, StructuralFeatureSet,
};
use crate::metadata::MoleculeMetadata;
use crate::model::ChemistryModel;
use crate::multicenter::{MulticenterBondAst, MulticenterBondViews};
use crate::noncovalent::{NoncovalentBondAst, NoncovalentBondViews};
use crate::resolve::ResolveConfig;
use crate::smiles::SmilesIoConfig;
use crate::stereo::{StereoAtomAst, StereoAtomViews, StereoBondAst, StereoBondViews, StereoLigand};
use crate::substructure::SubstructureSearchConfig;
use crate::transaction::MoleculeEditor;

/// A molecule: the owned graph-AST root.
#[pyclass(eq)]
#[derive(Debug, PartialEq)]
pub struct MoleculeAst(GraphIrMoleculeAst);

#[pymethods]
impl MoleculeAst {
    /// An empty molecule: zero atoms, zero bonds.
    #[new]
    fn new() -> Self {
        Self(GraphIrMoleculeAst::new())
    }

    /// Parse a molecule from its EDN representation under explicit construction defaults.
    #[staticmethod]
    #[pyo3(signature = (text, *, defaults=None))]
    fn parse(text: &str, defaults: Option<MoleculeDefaults>) -> PyResult<Self> {
        let defaults = defaults.unwrap_or_else(MoleculeDefaults::new).to_rust();
        let molecule = GraphIrMoleculeDsl::from_str(text)
            .map_err(parse_error)?
            .into_ast(&defaults);
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
        let defaults = defaults.unwrap_or_else(MoleculeDefaults::new).to_rust();
        let dsl = GraphIrMoleculeDsl::from_str(text).map_err(parse_error)?;
        let metadata = MoleculeMetadata::from_rust(dsl.metadata().clone());
        Ok((Self::from_rust(dsl.into_ast(&defaults)), metadata))
    }

    /// Render a canonical positional DSL representation without entity
    /// keywords or atom aliases.
    #[pyo3(signature = (*, defaults=None))]
    fn render(&self, defaults: Option<MoleculeDefaults>) -> String {
        let defaults = defaults.unwrap_or_else(MoleculeDefaults::new).to_rust();
        GraphIrMoleculeDsl::from_ast(&self.0, &defaults).to_string()
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
        let defaults = defaults.unwrap_or_else(MoleculeDefaults::new).to_rust();
        let lowered = GraphIrMoleculeDsl::from_ast(&self.0, &defaults)
            .into_parts()
            .0;
        GraphIrMoleculeDsl::new(lowered, metadata.to_rust())
            .map(|dsl| dsl.to_string())
            .map_err(metadata_error)
    }

    fn __str__(&self) -> String {
        self.render(None)
    }

    /// A molecule from its entries. Each bond is a `(first, second, bond)` triple:
    /// two atom indices into `atoms` and a `BondAst`. Each dative bond is a
    /// `(donors, acceptor, bond)` triple: a list of donor atom indices, one
    /// acceptor atom index, and a `DativeBondAst`. Each aromatic system is an
    /// `(atoms, system)` pair: a list of member atom indices and an `AromaticSystemAst`.
    /// Each multicenter bond is an `(atoms, bond)` pair: a list of member atom indices
    /// and a `MulticenterBondAst`. Each noncovalent bond is a `([first, second], bond)`
    /// pair: the two (unordered) endpoint atom indices and a `NoncovalentBondAst`. Each
    /// stereo atom / stereo bond is a `(site, ligands, value)` triple: the site atom / bond
    /// index, a list of `StereoLigand`s in frame order, and a `StereoAtomAst` / `StereoBondAst`.
    #[staticmethod]
    #[pyo3(signature = (atoms, *, bonds=Vec::new(), dative_bonds=Vec::new(), aromatic_systems=Vec::new(), multicenter_bonds=Vec::new(), noncovalent_bonds=Vec::new(), stereo_atoms=Vec::new(), stereo_bonds=Vec::new(), constraints=Vec::new()))]
    #[allow(clippy::too_many_arguments)] // one argument per entity family — the full molecule surface
    fn from_entries(
        py: Python<'_>,
        atoms: Vec<Py<AtomAst>>,
        bonds: Vec<(u32, u32, Py<BondAst>)>,
        dative_bonds: Vec<(Vec<u32>, u32, Py<DativeBondAst>)>,
        aromatic_systems: Vec<(Vec<u32>, Py<AromaticSystemAst>)>,
        multicenter_bonds: Vec<(Vec<u32>, Py<MulticenterBondAst>)>,
        noncovalent_bonds: Vec<([u32; 2], Py<NoncovalentBondAst>)>,
        stereo_atoms: Vec<(u32, Vec<StereoLigand>, Py<StereoAtomAst>)>,
        stereo_bonds: Vec<(u32, Vec<StereoLigand>, Py<StereoBondAst>)>,
        constraints: Vec<Py<Constraint>>,
    ) -> PyResult<Self> {
        let ast_atoms = atoms
            .iter()
            .map(|atom| atom.bind(py).borrow().inner().clone())
            .collect();
        let ast_bonds = bonds
            .iter()
            .map(|(first, second, bond)| {
                (
                    GraphIrAtomId(*first),
                    GraphIrAtomId(*second),
                    bond.bind(py).borrow().inner().clone(),
                )
            })
            .collect();
        let ast_dative = dative_bonds
            .iter()
            .map(|(donors, acceptor, bond)| {
                (
                    donors.iter().map(|&donor| GraphIrAtomId(donor)).collect(),
                    GraphIrAtomId(*acceptor),
                    bond.bind(py).borrow().inner().clone(),
                )
            })
            .collect();
        let ast_aromatic = aromatic_systems
            .iter()
            .map(|(atoms, system)| {
                (
                    atoms.iter().map(|&atom| GraphIrAtomId(atom)).collect(),
                    system.bind(py).borrow().inner().clone(),
                )
            })
            .collect();
        let ast_multicenter = multicenter_bonds
            .iter()
            .map(|(atoms, bond)| {
                (
                    atoms.iter().map(|&atom| GraphIrAtomId(atom)).collect(),
                    bond.bind(py).borrow().inner().clone(),
                )
            })
            .collect();
        let ast_noncovalent = noncovalent_bonds
            .iter()
            .map(|([first, second], bond)| {
                (
                    GraphIrAtomId(*first),
                    GraphIrAtomId(*second),
                    bond.bind(py).borrow().inner().clone(),
                )
            })
            .collect();
        let ast_stereo_atoms = stereo_atoms
            .iter()
            .map(|(site, ligands, value)| {
                (
                    GraphIrAtomId(*site),
                    ligands.iter().copied().map(StereoLigand::to_rust).collect(),
                    value.bind(py).borrow().inner().clone(),
                )
            })
            .collect();
        let ast_stereo_bonds = stereo_bonds
            .iter()
            .map(|(site, ligands, value)| {
                (
                    GraphIrBondId(*site),
                    ligands.iter().copied().map(StereoLigand::to_rust).collect(),
                    value.bind(py).borrow().inner().clone(),
                )
            })
            .collect();
        let ast_constraints = constraints
            .iter()
            .map(|constraint| constraint.bind(py).borrow().to_rust(py))
            .collect::<Vec<_>>();
        GraphIrMoleculeAst::try_from_entries(GraphIrMoleculeEntries {
            atoms: ast_atoms,
            bonds: ast_bonds,
            dative: ast_dative,
            aromatic: ast_aromatic,
            multicenter: ast_multicenter,
            noncovalent: ast_noncovalent,
            stereo_atoms: ast_stereo_atoms,
            stereo_bonds: ast_stereo_bonds,
            constraints: ast_constraints.into(),
        })
        .map(MoleculeAst)
        .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// Ingest a determined molecule from SMILES under explicit IO, chemistry,
    /// and resolution policies.
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
        let chemistry_model =
            chemistry_model.map_or_else(GraphChemistryModel::default, |model| model.to_rust());
        let resolve_config =
            resolve_config.map_or_else(GraphResolveConfig::default, ResolveConfig::to_rust);

        ingest_smiles_with(source, &io_config, &chemistry_model, &resolve_config)
            .map(Self::from_rust)
            .map_err(smiles_input_error)
    }

    /// Create a mutable editor initialized from this molecule.
    fn edit(&self) -> MoleculeEditor {
        MoleculeEditor::from_rust(self.0.edit())
    }

    /// Apply a checked edit batch without modifying this molecule.
    fn apply(&self, py: Python<'_>, edits: Py<Edits>) -> PyResult<Self> {
        self.0
            .apply(edits.bind(py).borrow().to_rust())
            .map(Self::from_rust)
            .map_err(transaction_error)
    }

    /// Combine this molecule and `other` by disjoint concatenation without modifying either input.
    /// The correspondence maps `other` into the combined molecule.
    fn combine(&self, other: &Self) -> (Self, MoleculeCorrespondence) {
        let (combined, correspondence) = self.0.combine(&other.0);
        (
            Self::from_rust(combined),
            MoleculeCorrespondence::from_rust(correspondence),
        )
    }

    /// Append `other` by disjoint concatenation. Existing ids in this molecule remain stable; the
    /// returned correspondence maps `other` into this molecule.
    fn combine_from(slf: Py<Self>, py: Python<'_>, other: Py<Self>) -> MoleculeCorrespondence {
        let other = other.bind(py).borrow().inner().clone();
        let correspondence = slf.borrow_mut(py).inner_mut().combine_from(&other);
        MoleculeCorrespondence::from_rust(correspondence)
    }

    /// Combine an iterable of molecules by disjoint concatenation. Returns one correspondence per
    /// input, in input order, mapping that input into the combined molecule.
    #[staticmethod]
    fn combine_all(
        py: Python<'_>,
        molecules: &Bound<'_, PyAny>,
    ) -> PyResult<(Self, Vec<MoleculeCorrespondence>)> {
        let molecules = molecules
            .try_iter()?
            .map(|item| -> PyResult<Py<MoleculeAst>> {
                Ok(item?.cast_into::<MoleculeAst>()?.unbind())
            })
            .collect::<PyResult<Vec<_>>>()?;
        let borrowed = molecules
            .iter()
            .map(|molecule| molecule.bind(py).borrow())
            .collect::<Vec<_>>();
        let (combined, correspondences) =
            GraphIrMoleculeAst::combine_all(borrowed.iter().map(|molecule| molecule.inner()));
        Ok((
            Self::from_rust(combined),
            correspondences
                .into_iter()
                .map(MoleculeCorrespondence::from_rust)
                .collect(),
        ))
    }

    /// Decompose this molecule into components connected by any relation. Each correspondence maps
    /// the returned component into this molecule.
    fn split(&self) -> Vec<(Self, MoleculeCorrespondence)> {
        self.0
            .split()
            .into_iter()
            .map(|(component, correspondence)| {
                (
                    Self::from_rust(component),
                    MoleculeCorrespondence::from_rust(correspondence),
                )
            })
            .collect()
    }

    /// Find occurrences of this pattern in `host`.
    #[pyo3(signature = (host, *, config=None))]
    fn substructure_matches(
        &self,
        host: &Self,
        config: Option<SubstructureSearchConfig>,
    ) -> Vec<MoleculeCorrespondence> {
        let config = config.unwrap_or_default().to_rust();
        self.0
            .substructure_matches(&host.0, config)
            .into_iter()
            .map(MoleculeCorrespondence::from_rust)
            .collect()
    }

    /// Generate an unfolded binary hashed fingerprint.
    #[pyo3(signature = (*, config))]
    fn hashed_fingerprint(&self, config: HashedFingerprintConfig) -> PyResult<HashedFeatureSet> {
        config
            .to_rust()
            .featurize(&self.0)
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
            .featurize_counted(&self.0)
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
            .fingerprint(&self.0)
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
            .featurize(&self.0)
            .map(StructuralFeatureSet::from_rust)
            .map_err(fingerprint_error)
    }

    /// The atoms, indexed by integer position.
    #[getter]
    fn atoms(slf: Py<Self>) -> AtomViews {
        AtomViews::new(slf)
    }

    /// The bonds, indexed by integer position.
    #[getter]
    fn bonds(slf: Py<Self>) -> BondViews {
        BondViews::new(slf)
    }

    /// The dative bonds, indexed by integer position.
    #[getter]
    fn dative_bonds(slf: Py<Self>) -> DativeBondViews {
        DativeBondViews::new(slf)
    }

    /// The aromatic systems, indexed by integer position.
    #[getter]
    fn aromatic_systems(slf: Py<Self>) -> AromaticSystemViews {
        AromaticSystemViews::new(slf)
    }

    /// The multicenter bonds, indexed by integer position.
    #[getter]
    fn multicenter_bonds(slf: Py<Self>) -> MulticenterBondViews {
        MulticenterBondViews::new(slf)
    }

    /// The noncovalent bonds, indexed by integer position.
    #[getter]
    fn noncovalent_bonds(slf: Py<Self>) -> NoncovalentBondViews {
        NoncovalentBondViews::new(slf)
    }

    /// The stereo atoms, indexed by integer position.
    #[getter]
    fn stereo_atoms(slf: Py<Self>) -> StereoAtomViews {
        StereoAtomViews::new(slf)
    }

    /// The stereo bonds, indexed by integer position.
    #[getter]
    fn stereo_bonds(slf: Py<Self>) -> StereoBondViews {
        StereoBondViews::new(slf)
    }

    /// The molecule-level constraints in insertion order.
    #[getter]
    fn constraints(slf: Py<Self>) -> ConstraintsView {
        ConstraintsView::new(slf)
    }

    /// Replace the molecule-level constraints from an owned container or live view.
    #[setter]
    fn set_constraints(slf: Py<Self>, py: Python<'_>, value: ConstraintsLike) -> PyResult<()> {
        let constraints = value.to_rust(py)?;
        *slf.borrow_mut(py).inner_mut().constraints_mut() = constraints;
        Ok(())
    }

    fn __repr__(&self) -> String {
        // Atoms and bonds always; the other entity families (dative bonds, aromatic systems,
        // multicenter bonds, noncovalent bonds, stereo atoms, stereo bonds) only when present,
        // so a plain covalent molecule stays uncluttered. Names match the `from_entries` kwargs.
        let mut parts = vec![
            format!("atoms={}", self.0.atoms().count()),
            format!("bonds={}", self.0.bonds().count()),
        ];
        for (name, count) in [
            ("dative_bonds", self.0.dative_bonds().count()),
            ("aromatic_systems", self.0.aromatic_systems().count()),
            ("multicenter_bonds", self.0.multicenter_bonds().count()),
            ("noncovalent_bonds", self.0.noncovalent_bonds().count()),
            ("stereo_atoms", self.0.stereo_atoms().count()),
            ("stereo_bonds", self.0.stereo_bonds().count()),
        ] {
            if count > 0 {
                parts.push(format!("{name}={count}"));
            }
        }
        format!("MoleculeAst({})", parts.join(", "))
    }
}

impl MoleculeAst {
    /// The wrapped AST molecule — read access for atom views.
    pub(crate) fn inner(&self) -> &GraphIrMoleculeAst {
        &self.0
    }

    /// Mutable access to the wrapped AST molecule — write access for the live
    /// atom and constraint views (copy-on-write through `atom_mut`).
    pub(crate) fn inner_mut(&mut self) -> &mut GraphIrMoleculeAst {
        &mut self.0
    }

    /// Wrap a Rust molecule as a Python molecule value.
    pub(crate) fn from_rust(molecule: GraphIrMoleculeAst) -> Self {
        MoleculeAst(molecule)
    }
}

#[cfg(test)]
mod tests {
    use pyo3::types::{PyBytes, PyList};
    use rstest::{fixture, rstest};
    use umol_chem::element::Element as ChemElement;
    use umol_graph::fingerprint::{
        CountedFeatureSet as GraphCountedFeatureSet, FeatureSet as GraphFeatureSet,
        SubstructureFeaturizer as GraphSubstructureFeaturizer,
    };
    use umol_graph::ingest::ingest_smiles;
    use umol_graph_core::{
        Correspondence as GraphCoreCorrespondence,
        RelevantCycleEnumerationAlgorithm as GraphCoreRelevantCycleEnumerationAlgorithm,
        SubgraphIsomorphismAlgorithm as GraphCoreSubgraphIsomorphismAlgorithm,
    };
    use umol_graph_ir::dsl::{
        AtomDsl as GraphIrAtomDsl, MoleculeMetadata as GraphIrMoleculeMetadata,
    };
    use umol_graph_ir::ir::{
        AromaticSystemAst as GraphIrAromaticSystemAst, AromaticSystemId as GraphIrAromaticSystemId,
        AtomAst as GraphIrAtomAst, AtomFieldChange as GraphIrAtomFieldChange,
        AtomHandle as GraphIrAtomHandle, AtomUpdate as GraphIrAtomUpdate,
        BondAst as GraphIrBondAst, Constraint as GraphIrConstraint,
        Constraints as GraphIrConstraints, DativeBondAst as GraphIrDativeBondAst,
        DativeBondId as GraphIrDativeBondId, Edit as GraphIrEdit, Edits as GraphIrEdits,
        Entity as GraphIrEntity, MoleculeConstraint as GraphIrMoleculeConstraint,
        MoleculeCorrespondence as GraphIrMoleculeCorrespondence,
        MulticenterBondAst as GraphIrMulticenterBondAst,
        MulticenterBondId as GraphIrMulticenterBondId,
        NoncovalentBondAst as GraphIrNoncovalentBondAst,
        NoncovalentBondId as GraphIrNoncovalentBondId,
        NoncovalentBondKind as GraphIrNoncovalentBondKind,
        SubstructureMatchAlgorithm as GraphIrSubstructureMatchAlgorithm,
        SubstructureMatchConfig as GraphIrSubstructureMatchConfig, ValueAst as GraphIrValueAst,
    };
    use umol_graph_ir::mol_dsl;

    use super::*;
    use crate::atom::AtomAst as PyAtomAst;
    use crate::constraint::molecule::Constraints;
    use crate::convert::into_py_variant;
    use crate::error::{MetadataError, ParseError, TransactionError, UnderdeterminedError};
    use crate::fingerprint::config::{
        EcfpHashScheme, PatternFingerprintConfig, RefinementRounds, StructuralFingerprintConfig,
        WlHashScheme,
    };
    use crate::ring::RingConfig;

    #[fixture]
    fn ethanol() -> MoleculeAst {
        MoleculeAst::from_rust(ingest_smiles("CCO").unwrap())
    }

    #[fixture]
    fn ethane() -> MoleculeAst {
        MoleculeAst::from_rust(ingest_smiles("CC").unwrap())
    }

    #[rstest]
    fn test_molecule_ast_new() {
        assert_eq!(MoleculeAst::new().inner().atoms().count(), 0);
    }

    #[rstest]
    #[case::required(
        r#"{:atoms ["C"]}"#,
        None,
        mol_dsl!(r#"{:atoms ["C"]}"#)
    )]
    #[case::ground(
        r#"{:atoms ["C#h4#v0#d0#t0#a!#m!"]}"#,
        Some(MoleculeDefaults::ground()),
        mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"]}"#)
    )]
    fn test_molecule_ast_parse(
        #[case] text: &str,
        #[case] defaults: Option<MoleculeDefaults>,
        #[case] expected: GraphIrMoleculeAst,
    ) {
        assert_eq!(
            MoleculeAst::parse(text, defaults).unwrap().inner(),
            &expected
        );
    }

    #[rstest]
    fn test_molecule_ast_parse_error() {
        Python::attach(|py| {
            let error = MoleculeAst::parse("not edn", None).unwrap_err();

            assert!(error.is_instance_of::<ParseError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "EDN parse: unexpected token 'n' at byte 0"
            );
        });
    }

    #[rstest]
    fn test_molecule_ast_parse_with_metadata() {
        let (molecule, metadata) = MoleculeAst::parse_with_metadata(
            r#"{:atoms [[:carbon :x]] :bonds [] :atom-aliases [:x "C"]}"#,
            None,
        )
        .unwrap();
        let metadata = metadata.to_rust();

        assert_eq!(molecule.inner(), &mol_dsl!(r#"{:atoms ["C"]}"#));
        assert_eq!(
            metadata.keyword(GraphIrEntity::Atom(GraphIrAtomId(0))),
            Some("carbon")
        );
        assert_eq!(
            metadata.atom_alias("x"),
            Some(&GraphIrAtomDsl(GraphIrAtomAst::from_element(
                ChemElement::C
            )))
        );
    }

    #[rstest]
    fn test_molecule_ast_parse_with_metadata_defaults() {
        let (molecule, metadata) = MoleculeAst::parse_with_metadata(
            r#"{:atoms ["C#h4#v0#d0#t0#a!#m!"]}"#,
            Some(MoleculeDefaults::ground()),
        )
        .unwrap();

        assert_eq!(
            molecule.inner(),
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
        Some(MoleculeDefaults::ground()),
        r#"{:atoms ["C#h4#v0#d0#t0#a!#m!"] :bonds []}"#
    )]
    fn test_molecule_ast_render(
        #[case] molecule: GraphIrMoleculeAst,
        #[case] defaults: Option<MoleculeDefaults>,
        #[case] expected: &str,
    ) {
        assert_eq!(MoleculeAst::from_rust(molecule).render(defaults), expected);
    }

    #[rstest]
    fn test_molecule_ast_render_with_metadata() {
        let molecule = MoleculeAst::from_rust(mol_dsl!(r#"{:atoms ["C"]}"#));
        let mut metadata = GraphIrMoleculeMetadata::new();
        metadata
            .set_keyword(GraphIrEntity::Atom(GraphIrAtomId(0)), "carbon")
            .unwrap();
        metadata
            .add_atom_alias(
                "x",
                GraphIrAtomDsl(GraphIrAtomAst::from_element(ChemElement::C)),
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
    fn test_molecule_ast_render_with_metadata_error() {
        Python::attach(|py| {
            let molecule = MoleculeAst::from_rust(mol_dsl!(r#"{:atoms ["C"]}"#));
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
    fn test_molecule_ast_str() {
        let molecule = MoleculeAst::from_rust(mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#));

        assert_eq!(molecule.__str__(), molecule.render(None));
    }

    #[rstest]
    fn test_molecule_ast_from_entries() {
        Python::attach(|py| {
            let atoms = vec![
                Py::new(
                    py,
                    PyAtomAst::from_inner(GraphIrAtomAst::from_element(ChemElement::C)),
                )
                .unwrap(),
                Py::new(
                    py,
                    PyAtomAst::from_inner(GraphIrAtomAst::from_element(ChemElement::B)),
                )
                .unwrap(),
                Py::new(
                    py,
                    PyAtomAst::from_inner(GraphIrAtomAst::from_element(ChemElement::N)),
                )
                .unwrap(),
            ];
            let bonds = vec![(
                0,
                1,
                Py::new(py, BondAst::from_inner(GraphIrBondAst::from_order(1))).unwrap(),
            )];
            let dative = vec![(
                vec![2],
                1,
                Py::new(
                    py,
                    DativeBondAst::from_inner(GraphIrDativeBondAst::from_order(1)),
                )
                .unwrap(),
            )];
            let aromatic = vec![(
                vec![0, 1, 2],
                Py::new(
                    py,
                    AromaticSystemAst::from_inner(GraphIrAromaticSystemAst::from_electrons(vec![
                        1, 1, 1,
                    ])),
                )
                .unwrap(),
            )];
            let multicenter = vec![(
                vec![0, 1, 2],
                Py::new(
                    py,
                    MulticenterBondAst::from_inner(GraphIrMulticenterBondAst::from_electrons(
                        vec![1, 1, 1],
                    )),
                )
                .unwrap(),
            )];
            let noncovalent = vec![(
                [0, 2],
                Py::new(
                    py,
                    NoncovalentBondAst::from_inner(GraphIrNoncovalentBondAst::from_kind(
                        GraphIrNoncovalentBondKind::HydrogenBond,
                    )),
                )
                .unwrap(),
            )];
            let constraint = GraphIrConstraint::Molecule(GraphIrMoleculeConstraint::Connected {
                atoms: Some(vec![GraphIrAtomId(0), GraphIrAtomId(2)]),
            });
            let constraints =
                vec![into_py_variant(py, Constraint::from_rust(py, &constraint).unwrap()).unwrap()];
            let molecule = MoleculeAst::from_entries(
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
            assert_eq!(molecule.inner().atoms().count(), 3);
            assert_eq!(molecule.inner().bonds().count(), 1);
            let dative_bonds = molecule.inner().dative_bonds();
            assert_eq!(dative_bonds.count(), 1);
            let dative_view = dative_bonds.get(GraphIrDativeBondId(0)).unwrap();
            assert_eq!(dative_view.acceptor_id(), GraphIrAtomId(1));
            assert_eq!(
                dative_view.donor_ids().collect::<Vec<_>>(),
                vec![GraphIrAtomId(2)]
            );
            let aromatic_systems = molecule.inner().aromatic_systems();
            assert_eq!(aromatic_systems.count(), 1);
            let aromatic_view = aromatic_systems.get(GraphIrAromaticSystemId(0)).unwrap();
            assert_eq!(
                aromatic_view.atom_ids().collect::<Vec<_>>(),
                vec![GraphIrAtomId(0), GraphIrAtomId(1), GraphIrAtomId(2)]
            );
            let multicenter_bonds = molecule.inner().multicenter_bonds();
            assert_eq!(multicenter_bonds.count(), 1);
            let multicenter_view = multicenter_bonds.get(GraphIrMulticenterBondId(0)).unwrap();
            assert_eq!(
                multicenter_view.atom_ids().collect::<Vec<_>>(),
                vec![GraphIrAtomId(0), GraphIrAtomId(1), GraphIrAtomId(2)]
            );
            let noncovalent_bonds = molecule.inner().noncovalent_bonds();
            assert_eq!(noncovalent_bonds.count(), 1);
            let noncovalent_view = noncovalent_bonds.get(GraphIrNoncovalentBondId(0)).unwrap();
            assert_eq!(
                noncovalent_view.atom_ids(),
                [GraphIrAtomId(0), GraphIrAtomId(2)]
            );
            assert_eq!(molecule.inner().constraints().as_slice(), &[constraint]);
        });
    }

    #[rstest]
    #[case::defaults(None, None, None)]
    #[case::explicit(
        Some(SmilesIoConfig::from_rust(&IoSmilesIoConfig::opensmiles())),
        Some(ChemistryModel::from_rust(&GraphChemistryModel::default())),
        Some(ResolveConfig::from_rust(GraphResolveConfig::default())),
    )]
    fn test_molecule_ast_from_smiles(
        #[case] io_config: Option<SmilesIoConfig>,
        #[case] chemistry_model: Option<ChemistryModel>,
        #[case] resolve_config: Option<ResolveConfig>,
    ) {
        assert_eq!(
            MoleculeAst::from_smiles("C", io_config, chemistry_model, resolve_config).unwrap(),
            MoleculeAst::from_rust(mol_dsl!(
                r#"{:atoms ["C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"]}"#
            ))
        );
    }

    #[rstest]
    #[case::syntax(" C", "ParseError", "Leading whitespace")]
    #[case::model_conversion(
        "C[S@]C",
        "ModelConversionError",
        "tetrahedral stereo at atom 1 with 2 ligands, expected 3 or 4 ligands"
    )]
    #[case::underdetermined("*", "UnderdeterminedError", "resolution underdetermined")]
    fn test_molecule_ast_from_smiles_error(
        #[case] source: &str,
        #[case] expected_type: &str,
        #[case] expected_message: &str,
    ) {
        Python::attach(|py| {
            let error = MoleculeAst::from_smiles(source, None, None, None).unwrap_err();
            assert_eq!(error.get_type(py).name().unwrap(), expected_type);
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                expected_message
            );
        });
    }

    #[rstest]
    fn test_molecule_ast_edit() {
        Python::attach(|py| {
            let expected = mol_dsl!(r#"{:atoms ["N#h3"]}"#);
            let editor = Py::new(py, MoleculeAst::from_rust(expected.clone()).edit()).unwrap();
            let snapshot = editor
                .bind(py)
                .call_method0("snapshot")
                .unwrap()
                .extract::<Py<MoleculeAst>>()
                .unwrap();

            assert_eq!(snapshot.bind(py).borrow().inner(), &expected);
        });
    }

    #[rstest]
    fn test_molecule_ast_apply() {
        let initial = mol_dsl!(r#"{:atoms ["N#h3"]}"#);
        let mut rust_edits = GraphIrEdits::new();
        rust_edits.update_atom(
            GraphIrAtomHandle::Id(GraphIrAtomId(0)),
            initial.atom(GraphIrAtomId(0)).ast,
            &GraphIrAtomUpdate {
                implicit_hydrogens: Some(GraphIrValueAst::Lit(2)),
                ..Default::default()
            },
        );
        let methyl = rust_edits
            .add_atom(GraphIrAtomAst::from_element(ChemElement::C).with_implicit_hydrogens(3_i64));
        rust_edits.add_bond(
            GraphIrAtomHandle::Id(GraphIrAtomId(0)),
            methyl,
            GraphIrBondAst::from_order(1),
        );
        let molecule = MoleculeAst::from_rust(initial.clone());

        Python::attach(|py| {
            let edits = Py::new(py, Edits::from_rust(rust_edits)).unwrap();

            let result = molecule.apply(py, edits).unwrap();

            assert_eq!(
                result.inner(),
                &mol_dsl!(r#"{:atoms ["N#h2" "C#h3"] :bonds [[0 1 "1"]]}"#)
            );
            assert_eq!(molecule.inner(), &initial);
        });
    }

    #[rstest]
    fn test_molecule_ast_apply_error() {
        let initial = mol_dsl!(r#"{:atoms ["C"]}"#);
        let molecule = MoleculeAst::from_rust(initial.clone());
        let mut rust_edits = GraphIrEdits::new();
        rust_edits.add_atom(GraphIrAtomAst::from_element(ChemElement::N));
        rust_edits.push(GraphIrEdit::ModifyAtomField {
            id: GraphIrAtomHandle::Id(GraphIrAtomId(7)),
            change: GraphIrAtomFieldChange::Charge {
                old: GraphIrValueAst::Lit(0),
                new: GraphIrValueAst::Lit(1),
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
            assert_eq!(molecule.inner(), &initial);
        });
    }

    #[rstest]
    fn test_molecule_ast_substructure_matches() {
        let pattern = MoleculeAst::from_rust(mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#));
        let host = MoleculeAst::from_rust(mol_dsl!(
            r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#
        ));
        let pattern_before = pattern.inner().clone();
        let host_before = host.inner().clone();
        let expected = vec![
            MoleculeCorrespondence::from_rust(
                GraphIrMoleculeCorrespondence::induce(
                    pattern.inner(),
                    host.inner(),
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
                    pattern.inner(),
                    host.inner(),
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

        assert_eq!(pattern.substructure_matches(&host, None), expected);
        assert_eq!(pattern.inner(), &pattern_before);
        assert_eq!(host.inner(), &host_before);
    }

    #[rstest]
    fn test_molecule_ast_substructure_matches_overlay() {
        let pattern = MoleculeAst::from_rust(mol_dsl!(
            r#"{
                :atoms ["N" "B"]
                :bonds []
                :dative-bonds [{:donors [0] :acceptor 1 :type "1"}]
            }"#
        ));
        let host = MoleculeAst::from_rust(mol_dsl!(
            r#"{
                :atoms ["N" "B" "C"]
                :bonds []
                :dative-bonds [{:donors [0] :acceptor 1 :type "1"}]
            }"#
        ));
        let expected = vec![MoleculeCorrespondence::from_rust(
            GraphIrMoleculeCorrespondence::induce(
                pattern.inner(),
                host.inner(),
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

        assert_eq!(pattern.substructure_matches(&host, Some(config)), expected);
    }

    #[rstest]
    fn test_molecule_ast_substructure_matches_empty() {
        let pattern = MoleculeAst::from_rust(mol_dsl!(r#"{:atoms ["O"] :bonds []}"#));
        let host = MoleculeAst::from_rust(mol_dsl!(r#"{:atoms ["C"] :bonds []}"#));
        let config = SubstructureSearchConfig::from_rust(GraphIrSubstructureMatchConfig {
            match_algorithm: GraphIrSubstructureMatchAlgorithm::GraphAndOverlays,
            subgraph_isomorphism_algorithm: GraphCoreSubgraphIsomorphismAlgorithm::Vf2,
            relevant_cycle_algorithm: GraphCoreRelevantCycleEnumerationAlgorithm::Vismara,
        });

        assert_eq!(
            pattern.substructure_matches(&host, Some(config)),
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
    fn test_molecule_ast_hashed_fingerprint(
        ethanol: MoleculeAst,
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
    fn test_molecule_ast_hashed_fingerprint_error(#[case] config: HashedFingerprintConfig) {
        Python::attach(|py| {
            let molecule = MoleculeAst::from_rust(mol_dsl!(r#"{:atoms ["C"] :bonds []}"#));
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
    fn test_molecule_ast_counted_hashed_fingerprint(
        ethane: MoleculeAst,
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
    fn test_molecule_ast_counted_hashed_fingerprint_error(#[case] config: HashedFingerprintConfig) {
        Python::attach(|py| {
            let molecule = MoleculeAst::from_rust(mol_dsl!(r#"{:atoms ["C"] :bonds []}"#));
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
    fn test_molecule_ast_pattern_fingerprint(
        ethanol: MoleculeAst,
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
    fn test_molecule_ast_pattern_fingerprint_error(
        #[case] config: Option<PatternFingerprintConfig>,
    ) {
        Python::attach(|py| {
            let molecule = MoleculeAst::from_rust(mol_dsl!(r#"{:atoms ["C"] :bonds []}"#));
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
    fn test_molecule_ast_structural_fingerprint(
        ethanol: MoleculeAst,
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
    fn test_molecule_ast_structural_fingerprint_error(#[case] config: StructuralFingerprintConfig) {
        Python::attach(|py| {
            let molecule = MoleculeAst::from_rust(mol_dsl!(r#"{:atoms ["C"] :bonds []}"#));
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
    fn test_molecule_ast_atoms(#[case] elements: Vec<ChemElement>, #[case] expected: usize) {
        let atoms = elements
            .into_iter()
            .map(GraphIrAtomAst::from_element)
            .collect();
        let molecule = MoleculeAst(GraphIrMoleculeAst::from_entries(GraphIrMoleculeEntries {
            atoms,
            ..Default::default()
        }));
        assert_eq!(molecule.inner().atoms().count(), expected);
    }

    #[rstest]
    fn test_molecule_ast_constraints() {
        Python::attach(|py| {
            let molecule = Py::new(py, MoleculeAst::new()).unwrap();
            let view = MoleculeAst::constraints(molecule.clone_ref(py));
            let constraint =
                GraphIrConstraint::Molecule(GraphIrMoleculeConstraint::Connected { atoms: None });
            view.with_mut(py, |constraints| constraints.push(constraint.clone()));

            assert_eq!(
                molecule.bind(py).borrow().inner().constraints().as_slice(),
                &[constraint]
            );
        });
    }

    #[rstest]
    fn test_molecule_ast_set_constraints() {
        Python::attach(|py| {
            let molecule = Py::new(py, MoleculeAst::new()).unwrap();
            let constraint = GraphIrConstraint::Molecule(GraphIrMoleculeConstraint::Connected {
                atoms: Some(vec![]),
            });
            let constraints = Py::new(
                py,
                Constraints::from_inner(GraphIrConstraints::from(vec![constraint.clone()])),
            )
            .unwrap();

            MoleculeAst::set_constraints(
                molecule.clone_ref(py),
                py,
                ConstraintsLike::Container(constraints),
            )
            .unwrap();

            assert_eq!(
                molecule.bind(py).borrow().inner().constraints().as_slice(),
                &[constraint]
            );
        });
    }

    #[rstest]
    fn test_molecule_ast_set_constraints_self() {
        Python::attach(|py| {
            let constraint =
                GraphIrConstraint::Molecule(GraphIrMoleculeConstraint::Connected { atoms: None });
            let molecule = Py::new(
                py,
                MoleculeAst(GraphIrMoleculeAst::from_entries(GraphIrMoleculeEntries {
                    constraints: GraphIrConstraints::from(vec![constraint.clone()]),
                    ..Default::default()
                })),
            )
            .unwrap();
            let own_view = Py::new(py, MoleculeAst::constraints(molecule.clone_ref(py))).unwrap();

            MoleculeAst::set_constraints(
                molecule.clone_ref(py),
                py,
                ConstraintsLike::View(own_view),
            )
            .unwrap();

            assert_eq!(
                molecule.bind(py).borrow().inner().constraints().as_slice(),
                &[constraint]
            );
        });
    }

    #[rstest]
    fn test_molecule_ast_eq() {
        assert_eq!(MoleculeAst::new(), MoleculeAst::new());
        let carbon = MoleculeAst(GraphIrMoleculeAst::from_entries(GraphIrMoleculeEntries {
            atoms: vec![GraphIrAtomAst::from_element(ChemElement::C)],
            ..Default::default()
        }));
        assert_ne!(MoleculeAst::new(), carbon);
    }

    #[rstest]
    #[case::empty(MoleculeAst::new(), "MoleculeAst(atoms=0, bonds=0)")]
    #[case::noncovalent(
        MoleculeAst(GraphIrMoleculeAst::from_entries(GraphIrMoleculeEntries {
            atoms: vec![
                GraphIrAtomAst::from_element(ChemElement::O),
                GraphIrAtomAst::from_element(ChemElement::O),
            ],
            noncovalent: vec![(
                GraphIrAtomId(0),
                GraphIrAtomId(1),
                GraphIrNoncovalentBondAst::from_kind(GraphIrNoncovalentBondKind::HydrogenBond),
            )],
            ..Default::default()
        })),
        "MoleculeAst(atoms=2, bonds=0, noncovalent_bonds=1)"
    )]
    fn test_molecule_ast_repr(#[case] molecule: MoleculeAst, #[case] expected: &str) {
        assert_eq!(molecule.__repr__(), expected);
    }
}
