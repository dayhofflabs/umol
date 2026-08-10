//! Noncovalent-bond kind values, owned ASTs, and molecule-backed views.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use std::str::FromStr;
use std::vec::IntoIter;

use pyo3::exceptions::{PyIndexError, PyTypeError};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use umol_graph_ir::ir::{
    AsLit, AtomId as GraphIrAtomId, Molecule as GraphIrMolecule,
    NoncovalentBondForm as GraphIrNoncovalentBondForm,
    NoncovalentBondId as GraphIrNoncovalentBondId,
    NoncovalentBondKind as GraphIrNoncovalentBondKind,
    NoncovalentBondKindForm as GraphIrNoncovalentBondKindForm,
    NoncovalentBondUpdate as GraphIrNoncovalentBondUpdate,
    NoncovalentBondView as GraphIrNoncovalentBondView,
};

use crate::constraint::noncovalent::{
    noncovalent_bond_constraints_asdict, NoncovalentBondConstraintsBacking,
    NoncovalentBondConstraintsForm, NoncovalentBondConstraintsLike, NoncovalentBondConstraintsView,
};
#[cfg(test)]
use crate::constraint::noncovalent::{
    NoncovalentBondConstraintForm, NoncovalentBondConstraintKey, NoncovalentBondConstraintsUpdate,
};
use crate::convert::{hash_rust, variant_repr};
use crate::entity::EntityForm;
use crate::error::parse_error;
use crate::lattice::impl_py_lattice;
use crate::molecule::MoleculeAst;

/// A noncovalent interaction kind. A fieldless, hashable value enum whose members
/// value the Rust `NoncovalentBondKind` exactly.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum NoncovalentBondKind {
    HydrogenBond,
    HalogenBond,
    ChalcogenBond,
    Ionic,
    VanDerWaals,
}

impl NoncovalentBondKind {
    pub(crate) fn from_rust(ast: GraphIrNoncovalentBondKind) -> Self {
        match ast {
            GraphIrNoncovalentBondKind::HydrogenBond => Self::HydrogenBond,
            GraphIrNoncovalentBondKind::HalogenBond => Self::HalogenBond,
            GraphIrNoncovalentBondKind::ChalcogenBond => Self::ChalcogenBond,
            GraphIrNoncovalentBondKind::Ionic => Self::Ionic,
            GraphIrNoncovalentBondKind::VanDerWaals => Self::VanDerWaals,
        }
    }

    pub(crate) fn to_rust(self) -> GraphIrNoncovalentBondKind {
        match self {
            Self::HydrogenBond => GraphIrNoncovalentBondKind::HydrogenBond,
            Self::HalogenBond => GraphIrNoncovalentBondKind::HalogenBond,
            Self::ChalcogenBond => GraphIrNoncovalentBondKind::ChalcogenBond,
            Self::Ionic => GraphIrNoncovalentBondKind::Ionic,
            Self::VanDerWaals => GraphIrNoncovalentBondKind::VanDerWaals,
        }
    }
}

/// A noncovalent bond's interaction kind: undetermined, or a concrete
/// `NoncovalentBondKind`. Corresponds to the Rust `NoncovalentBondKindForm`.
#[pyclass]
pub enum NoncovalentBondKindForm {
    Undetermined(),
    Lit(NoncovalentBondKind),
}

#[pymethods]
impl NoncovalentBondKindForm {
    /// The concrete interaction kind, or `None` when undetermined.
    fn as_lit(&self) -> Option<NoncovalentBondKind> {
        self.to_rust().as_lit().map(NoncovalentBondKind::from_rust)
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            NoncovalentBondKindForm::Undetermined() => ("Undetermined", 0),
            NoncovalentBondKindForm::Lit(_) => ("Lit", 1),
        };
        variant_repr(
            slf.bind(py).as_any(),
            "NoncovalentBondKindForm",
            variant,
            arity,
        )
    }
}

impl_py_lattice!(
    NoncovalentBondKindForm,
    GraphIrNoncovalentBondKindForm,
    |value: &NoncovalentBondKindForm,
     _py: Python<'_>|
     -> PyResult<GraphIrNoncovalentBondKindForm> { Ok(value.to_rust()) },
    |_py: Python<'_>, value: GraphIrNoncovalentBondKindForm| -> PyResult<NoncovalentBondKindForm> {
        Ok(NoncovalentBondKindForm::from_rust(&value))
    }
);

impl NoncovalentBondKindForm {
    pub(crate) fn from_rust(ast: &GraphIrNoncovalentBondKindForm) -> Self {
        match ast {
            GraphIrNoncovalentBondKindForm::Undetermined => Self::Undetermined(),
            GraphIrNoncovalentBondKindForm::Lit(k) => Self::Lit(NoncovalentBondKind::from_rust(*k)),
        }
    }

    pub(crate) fn to_rust(&self) -> GraphIrNoncovalentBondKindForm {
        match self {
            Self::Undetermined() => GraphIrNoncovalentBondKindForm::Undetermined,
            Self::Lit(k) => GraphIrNoncovalentBondKindForm::Lit(k.to_rust()),
        }
    }
}

/// Setter coercion for a noncovalent `kind` field: a bare `NoncovalentBondKind` →
/// `Lit`, or a `NoncovalentBondKindForm` passthrough (matching the `Undetermined |
/// Lit` structure).
#[derive(FromPyObject)]
pub(crate) enum NoncovalentBondKindLike {
    Kind(NoncovalentBondKind),
    Ast(Py<NoncovalentBondKindForm>),
}

impl NoncovalentBondKindLike {
    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrNoncovalentBondKindForm {
        match self {
            NoncovalentBondKindLike::Kind(k) => GraphIrNoncovalentBondKindForm::Lit(k.to_rust()),
            NoncovalentBondKindLike::Ast(a) => a.bind(py).borrow().to_rust(),
        }
    }
}

/// Attribute updates for a noncovalent bond.
#[pyclass(frozen, skip_from_py_object)]
#[derive(Clone)]
pub struct NoncovalentBondUpdate(GraphIrNoncovalentBondUpdate);

#[pymethods]
impl NoncovalentBondUpdate {
    #[new]
    #[pyo3(signature = (*, kind=None, constraints=None))]
    fn new(
        py: Python<'_>,
        kind: Option<NoncovalentBondKindLike>,
        constraints: Option<Py<NoncovalentBondConstraintsForm>>,
    ) -> Self {
        Self::from_rust(&GraphIrNoncovalentBondUpdate {
            kind: kind.map(|value| value.to_rust(py)),
            constraints: constraints
                .map(|value| value.bind(py).borrow().to_rust().clone())
                .unwrap_or_default(),
        })
    }

    /// Parse a noncovalent-bond-update DSL string into a `NoncovalentBondUpdate`.
    #[staticmethod]
    fn parse(s: &str) -> PyResult<Self> {
        GraphIrNoncovalentBondUpdate::from_str(s)
            .map(Self)
            .map_err(parse_error)
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("NoncovalentBondUpdate.parse('{}')", self.0)
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    #[getter]
    fn kind(&self) -> Option<NoncovalentBondKindForm> {
        self.0.kind.as_ref().map(NoncovalentBondKindForm::from_rust)
    }

    #[getter]
    fn constraints(&self) -> NoncovalentBondConstraintsForm {
        NoncovalentBondConstraintsForm::from_rust(self.0.constraints.clone())
    }
}

impl NoncovalentBondUpdate {
    pub(crate) fn from_rust(update: &GraphIrNoncovalentBondUpdate) -> Self {
        Self(update.clone())
    }

    pub(crate) fn to_rust(&self) -> &GraphIrNoncovalentBondUpdate {
        &self.0
    }
}

/// A noncovalent bond: an interaction `kind` plus noncovalent-bond-scope constraints.
/// No bond order, charge, or unpaired-electron field — these do not apply to
/// noncovalent interactions.
/// The endpoint atom pair is the owning molecule's relation topology (the view half),
/// not part of the value.
#[pyclass(eq)]
pub struct NoncovalentBondForm {
    value: GraphIrNoncovalentBondForm,
    readonly: bool,
}

impl PartialEq for NoncovalentBondForm {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

#[pymethods]
impl NoncovalentBondForm {
    /// Construct from an interaction kind — a `NoncovalentBondKind` or a
    /// `NoncovalentBondKindForm` — optionally setting constraints.
    #[new]
    #[pyo3(signature = (kind, *, constraints=None))]
    fn new(
        py: Python<'_>,
        kind: NoncovalentBondKindLike,
        constraints: Option<Py<NoncovalentBondConstraintsForm>>,
    ) -> Self {
        let mut bond = GraphIrNoncovalentBondForm::new(kind.to_rust(py));
        if let Some(constraints) = constraints {
            bond.constraints = constraints.bind(py).borrow().to_rust().clone();
        }
        NoncovalentBondForm::from_rust(bond)
    }

    /// Parse a noncovalent-bond-DSL string (e.g. `"Hbd#I"`) into a `NoncovalentBondForm`.
    #[staticmethod]
    fn parse(s: &str) -> PyResult<Self> {
        GraphIrNoncovalentBondForm::from_str(s)
            .map(Self::from_rust)
            .map_err(parse_error)
    }

    fn __str__(&self) -> String {
        self.value.to_string()
    }

    fn __repr__(&self) -> String {
        format!("NoncovalentBondForm.parse('{}')", self.value)
    }

    /// The interaction kind.
    #[getter]
    fn kind(&self) -> NoncovalentBondKindForm {
        NoncovalentBondKindForm::from_rust(&self.value.kind)
    }

    #[setter]
    fn set_kind(&mut self, py: Python<'_>, value: NoncovalentBondKindLike) -> PyResult<()> {
        self.to_rust_mut()?.kind = value.to_rust(py);
        Ok(())
    }

    /// The bond's constraints as a live handle onto this bond: reads borrow the
    /// current state, mutators write through to the bond in place.
    #[getter]
    fn constraints(slf: Py<Self>) -> NoncovalentBondConstraintsView {
        NoncovalentBondConstraintsView {
            backing: NoncovalentBondConstraintsBacking::Noncovalent(slf),
        }
    }

    /// Replace the whole constraint set (wipe-and-set) from a value container or a
    /// live view. Takes `slf` by handle and snapshots `value` *before* the write borrow,
    /// so `bond.constraints = bond.constraints` (a view over the same bond) reads while
    /// the bond is unborrowed instead of self-aliasing into a double-borrow panic.
    #[setter]
    fn set_constraints(
        slf: Py<Self>,
        py: Python<'_>,
        value: NoncovalentBondConstraintsLike,
    ) -> PyResult<()> {
        let snapshot = value.to_rust(py)?;
        slf.borrow_mut(py).to_rust_mut()?.constraints = snapshot;
        Ok(())
    }

    #[getter]
    fn readonly(&self) -> bool {
        self.readonly
    }

    fn copy(&self) -> Self {
        Self::from_rust(self.to_rust().clone())
    }

    /// The fields as a dict keyed by field name; values are Python objects.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("kind", self.kind())?;
        dict.set_item(
            "constraints",
            noncovalent_bond_constraints_asdict(py, &self.value.constraints)?,
        )?;
        Ok(dict)
    }
}

impl NoncovalentBondForm {
    /// The wrapped AST bond — read access for the bond-backed constraints view.
    pub(crate) fn to_rust(&self) -> &GraphIrNoncovalentBondForm {
        &self.value
    }

    /// Mutable access to the wrapped AST bond — write access for the bond-backed
    /// constraints view.
    pub(crate) fn to_rust_mut(&mut self) -> PyResult<&mut GraphIrNoncovalentBondForm> {
        if self.readonly {
            Err(PyTypeError::new_err("read-only entity form"))
        } else {
            Ok(&mut self.value)
        }
    }

    /// Wrap an owned Rust noncovalent-bond AST.
    pub(crate) fn from_rust(bond: GraphIrNoncovalentBondForm) -> Self {
        Self {
            value: bond,
            readonly: false,
        }
    }
}

impl EntityForm for NoncovalentBondForm {
    type RustForm = GraphIrNoncovalentBondForm;
    fn to_rust(&self) -> &Self::RustForm {
        &self.value
    }
    fn new_readonly(py: Python<'_>, value: Self::RustForm) -> PyResult<Py<Self>> {
        Py::new(
            py,
            Self {
                value,
                readonly: true,
            },
        )
    }
}

impl_py_lattice!(
    NoncovalentBondForm,
    GraphIrNoncovalentBondForm,
    |value: &NoncovalentBondForm, _py: Python<'_>| -> PyResult<GraphIrNoncovalentBondForm> {
        Ok(value.to_rust().clone())
    },
    |_py: Python<'_>, value: GraphIrNoncovalentBondForm| -> PyResult<NoncovalentBondForm> {
        Ok(NoncovalentBondForm::from_rust(value))
    }
);

/// A view of one noncovalent bond within a molecule: a handle to the molecule plus the
/// bond's index. Field reads rebuild the transient Rust view; the molecule is never
/// copied. The two endpoint atom indices are read-only topology; the kind and
/// constraints are the mutable bond value.
#[pyclass]
pub struct NoncovalentBondView {
    owner: Py<MoleculeAst>,
    id: GraphIrNoncovalentBondId,
}

impl NoncovalentBondView {
    fn noncovalent_bond<'a>(
        &self,
        molecule: &'a GraphIrMolecule,
    ) -> PyResult<GraphIrNoncovalentBondView<'a>> {
        molecule
            .noncovalent_bonds()
            .get(self.id)
            .ok_or_else(|| PyIndexError::new_err("noncovalent bond id out of range"))
    }
}

#[pymethods]
impl NoncovalentBondView {
    #[getter]
    fn id(&self) -> u32 {
        self.id.0
    }

    /// The two endpoint atom indices (read-only — participants are topology, not part of
    /// the bond value; the pair is unordered).
    #[getter]
    fn atom_ids(&self, py: Python<'_>) -> PyResult<(u32, u32)> {
        let molecule = self.owner.bind(py).borrow();
        let [first, second] = self.noncovalent_bond(molecule.to_rust())?.atom_ids();
        Ok((first.0, second.0))
    }

    fn __repr__(&self) -> String {
        format!("NoncovalentBondView(id={})", self.id.0)
    }

    /// The interaction kind.
    #[getter]
    fn kind(&self, py: Python<'_>) -> PyResult<NoncovalentBondKindForm> {
        let molecule = self.owner.bind(py).borrow();
        Ok(NoncovalentBondKindForm::from_rust(
            &self.noncovalent_bond(molecule.to_rust())?.attributes.kind,
        ))
    }

    #[setter]
    fn set_kind(&self, py: Python<'_>, value: NoncovalentBondKindLike) {
        self.owner
            .borrow_mut(py)
            .to_rust_mut()
            .noncovalent_bond_mut(self.id)
            .attributes
            .kind = value.to_rust(py);
    }

    /// The bond's constraints as a live handle onto the molecule: reads borrow the
    /// current state, mutators write through to the bond in place.
    #[getter]
    fn constraints(&self, py: Python<'_>) -> NoncovalentBondConstraintsView {
        NoncovalentBondConstraintsView {
            backing: NoncovalentBondConstraintsBacking::Molecule {
                owner: self.owner.clone_ref(py),
                id: self.id,
            },
        }
    }

    /// Replace the whole constraint set of the backing bond in place (wipe-and-set)
    /// from a value container or a live view.
    #[setter]
    fn set_constraints(
        &self,
        py: Python<'_>,
        value: NoncovalentBondConstraintsLike,
    ) -> PyResult<()> {
        self.owner
            .borrow_mut(py)
            .to_rust_mut()
            .noncovalent_bond_mut(self.id)
            .attributes
            .constraints = value.to_rust(py)?;
        Ok(())
    }

    /// The value fields as a dict keyed by field name; values are Python objects —
    /// symmetric with `NoncovalentBondForm.asdict`, read through the view.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let molecule = self.owner.bind(py).borrow();
        let bond = self.noncovalent_bond(molecule.to_rust())?.attributes;
        let dict = PyDict::new(py);
        dict.set_item("kind", NoncovalentBondKindForm::from_rust(&bond.kind))?;
        dict.set_item(
            "constraints",
            noncovalent_bond_constraints_asdict(py, &bond.constraints)?,
        )?;
        Ok(dict)
    }
}

/// Resolve a possibly-negative Python index (negative counts from the end) into an
/// existing noncovalent bond id, or `IndexError`. `NoncovalentBondId` is `RelationId`-
/// backed but contiguous for fresh molecules, so integer positions address it directly.
fn resolve_noncovalent_bond_index(
    molecule: &GraphIrMolecule,
    index: isize,
) -> PyResult<GraphIrNoncovalentBondId> {
    let count = molecule.noncovalent_bonds().count();
    let resolved = if index < 0 {
        index + count as isize
    } else {
        index
    };
    if resolved < 0 {
        return Err(PyIndexError::new_err("noncovalent bond id out of range"));
    }
    let id = GraphIrNoncovalentBondId(resolved as u32);
    if molecule.noncovalent_bonds().contains(id) {
        Ok(id)
    } else {
        Err(PyIndexError::new_err("noncovalent bond id out of range"))
    }
}

/// The noncovalent bonds of a molecule, indexed by integer position.
#[pyclass]
pub struct NoncovalentBondViews {
    owner: Py<MoleculeAst>,
}

#[pymethods]
impl NoncovalentBondViews {
    fn __len__(&self, py: Python<'_>) -> usize {
        self.owner
            .bind(py)
            .borrow()
            .to_rust()
            .noncovalent_bonds()
            .count()
    }

    fn __repr__(&self, py: Python<'_>) -> String {
        format!(
            "NoncovalentBondViews(len={})",
            self.owner
                .bind(py)
                .borrow()
                .to_rust()
                .noncovalent_bonds()
                .count()
        )
    }

    fn __getitem__(&self, py: Python<'_>, index: isize) -> PyResult<NoncovalentBondView> {
        let molecule = self.owner.bind(py).borrow();
        let id = resolve_noncovalent_bond_index(molecule.to_rust(), index)?;
        Ok(NoncovalentBondView {
            owner: self.owner.clone_ref(py),
            id,
        })
    }

    /// Replace the whole noncovalent bond value at `index` in place (endpoints unchanged).
    fn __setitem__(
        &self,
        py: Python<'_>,
        index: isize,
        bond: PyRef<'_, NoncovalentBondForm>,
    ) -> PyResult<()> {
        let mut molecule = self.owner.borrow_mut(py);
        let id = resolve_noncovalent_bond_index(molecule.to_rust(), index)?;
        *molecule.to_rust_mut().noncovalent_bond_mut(id).attributes = bond.to_rust().clone();
        Ok(())
    }

    /// The noncovalent bond between atoms `first` and `second`, or `None`.
    fn of(&self, py: Python<'_>, first: u32, second: u32) -> Option<NoncovalentBondView> {
        let molecule = self.owner.bind(py).borrow();
        molecule
            .to_rust()
            .noncovalent_bonds()
            .of_id(GraphIrAtomId(first), GraphIrAtomId(second))
            .map(|id| NoncovalentBondView {
                owner: self.owner.clone_ref(py),
                id,
            })
    }

    /// The noncovalent bonds `atom` is an endpoint of.
    fn incident(&self, py: Python<'_>, atom: u32) -> Vec<NoncovalentBondView> {
        let molecule = self.owner.bind(py).borrow();
        molecule
            .to_rust()
            .noncovalent_bonds()
            .incident_ids(GraphIrAtomId(atom))
            .map(|id| NoncovalentBondView {
                owner: self.owner.clone_ref(py),
                id,
            })
            .collect()
    }

    fn __iter__(&self, py: Python<'_>) -> NoncovalentBondViewIter {
        let ids = self
            .owner
            .bind(py)
            .borrow()
            .to_rust()
            .noncovalent_bonds()
            .ids()
            .collect::<Vec<_>>();
        NoncovalentBondViewIter {
            owner: self.owner.clone_ref(py),
            ids: ids.into_iter(),
        }
    }
}

impl NoncovalentBondViews {
    /// Build the noncovalent-bond-views handle for `owner` (the `.noncovalent_bonds` accessor).
    pub(crate) fn new(owner: Py<MoleculeAst>) -> NoncovalentBondViews {
        NoncovalentBondViews { owner }
    }
}

#[pyclass]
struct NoncovalentBondViewIter {
    owner: Py<MoleculeAst>,
    ids: IntoIter<GraphIrNoncovalentBondId>,
}

#[pymethods]
impl NoncovalentBondViewIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> Option<NoncovalentBondView> {
        self.ids.next().map(|id| NoncovalentBondView {
            owner: self.owner.clone_ref(py),
            id,
        })
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_chem::element::Element as ChemElement;
    use umol_graph_ir::ir::{
        AtomForm as GraphIrAtomForm, BooleanForm as GraphIrBooleanForm, MoleculeEntries,
        NoncovalentBondConstraintForm as GraphIrNoncovalentBondConstraintForm,
        NoncovalentBondConstraintKey as GraphIrNoncovalentBondConstraintKey,
        NoncovalentBondConstraintsForm as GraphIrNoncovalentBondConstraintsForm,
    };

    use super::*;
    use crate::boolean::{BooleanForm, BooleanLike};
    use crate::convert::into_py_variant;

    #[rstest]
    #[case(GraphIrNoncovalentBondKindForm::Undetermined)]
    #[case(GraphIrNoncovalentBondKindForm::Lit(GraphIrNoncovalentBondKind::HydrogenBond))]
    #[case(GraphIrNoncovalentBondKindForm::Lit(GraphIrNoncovalentBondKind::VanDerWaals))]
    fn test_noncovalent_bond_kind_form_roundtrip(#[case] ast: GraphIrNoncovalentBondKindForm) {
        assert_eq!(NoncovalentBondKindForm::from_rust(&ast).to_rust(), ast);
    }

    #[rstest]
    #[case(GraphIrNoncovalentBondKindForm::Undetermined, None)]
    #[case(
        GraphIrNoncovalentBondKindForm::Lit(GraphIrNoncovalentBondKind::Ionic),
        Some(NoncovalentBondKind::Ionic)
    )]
    fn test_noncovalent_bond_kind_form_as_lit(
        #[case] ast: GraphIrNoncovalentBondKindForm,
        #[case] expected: Option<NoncovalentBondKind>,
    ) {
        assert_eq!(NoncovalentBondKindForm::from_rust(&ast).as_lit(), expected);
    }

    #[rstest]
    #[case(GraphIrNoncovalentBondKind::HydrogenBond)]
    #[case(GraphIrNoncovalentBondKind::ChalcogenBond)]
    fn test_noncovalent_bond_kind_roundtrip(#[case] ast: GraphIrNoncovalentBondKind) {
        assert_eq!(NoncovalentBondKind::from_rust(ast).to_rust(), ast);
    }

    #[rstest]
    fn test_noncovalent_bond_kind_like_to_rust() {
        Python::attach(|py| {
            // a bare kind coerces to Lit
            assert_eq!(
                NoncovalentBondKindLike::Kind(NoncovalentBondKind::HydrogenBond).to_rust(py),
                GraphIrNoncovalentBondKindForm::Lit(GraphIrNoncovalentBondKind::HydrogenBond)
            );
            // a NoncovalentBondKindForm passes through
            let ast =
                Py::new(py, NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic)).unwrap();
            assert_eq!(
                NoncovalentBondKindLike::Ast(ast).to_rust(py),
                GraphIrNoncovalentBondKindForm::Lit(GraphIrNoncovalentBondKind::Ionic)
            );
        });
    }

    /// A `Py<NoncovalentBondConstraintForm>` for `Intramolecular(b)`.
    fn intramolecular(py: Python<'_>, b: bool) -> Py<NoncovalentBondConstraintForm> {
        into_py_variant(
            py,
            NoncovalentBondConstraintForm::from_rust(
                py,
                &GraphIrNoncovalentBondConstraintForm::intramolecular(b),
            )
            .unwrap(),
        )
        .unwrap()
    }

    fn intramolecular_key(py: Python<'_>) -> Py<NoncovalentBondConstraintKey> {
        into_py_variant(py, NoncovalentBondConstraintKey::Intramolecular()).unwrap()
    }

    #[rstest]
    fn test_noncovalent_bond_constraint_key_roundtrip() {
        let key = NoncovalentBondConstraintKey::from_rust(
            &GraphIrNoncovalentBondConstraintKey::Intramolecular,
        );
        assert_eq!(
            key.to_rust(),
            GraphIrNoncovalentBondConstraintKey::Intramolecular
        );
    }

    #[rstest]
    fn test_noncovalent_bond_constraint_form_key() {
        Python::attach(|py| {
            let constraint = GraphIrNoncovalentBondConstraintForm::intramolecular(true);
            let key = NoncovalentBondConstraintForm::from_rust(py, &constraint)
                .unwrap()
                .key(py);
            assert_eq!(
                key.to_rust(),
                GraphIrNoncovalentBondConstraintKey::Intramolecular
            );
        });
    }

    #[rstest]
    #[case(GraphIrNoncovalentBondConstraintForm::intramolecular(true))]
    #[case(GraphIrNoncovalentBondConstraintForm::intramolecular(false))]
    #[case(GraphIrNoncovalentBondConstraintForm::Intramolecular(GraphIrBooleanForm::Undetermined))]
    fn test_noncovalent_bond_constraint_form_roundtrip(
        #[case] ast: GraphIrNoncovalentBondConstraintForm,
    ) {
        Python::attach(|py| {
            assert_eq!(
                NoncovalentBondConstraintForm::from_rust(py, &ast)
                    .unwrap()
                    .to_rust(py),
                ast
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_form_new() {
        Python::attach(|py| {
            let constraints =
                NoncovalentBondConstraintsForm::new(py, vec![intramolecular(py, true)]);
            assert_eq!(constraints.__len__(), 1);
            assert_eq!(
                constraints.intramolecular().to_rust(),
                GraphIrBooleanForm::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_form_repr() {
        Python::attach(|py| {
            let constraints =
                NoncovalentBondConstraintsForm::new(py, vec![intramolecular(py, true)]);
            assert_eq!(
                constraints.__repr__(py).unwrap(),
                "NoncovalentBondConstraintsForm([NoncovalentBondConstraintForm.Intramolecular(BooleanForm.Lit(True))])"
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_form_set() {
        Python::attach(|py| {
            let mut constraints = NoncovalentBondConstraintsForm::new(py, vec![]);
            constraints.set(py, intramolecular(py, false));
            assert_eq!(constraints.__len__(), 1);
            assert_eq!(
                constraints.intramolecular().to_rust(),
                GraphIrBooleanForm::Lit(false)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_form_pop() {
        Python::attach(|py| {
            let mut constraints =
                NoncovalentBondConstraintsForm::new(py, vec![intramolecular(py, true)]);
            let removed = constraints.pop(py, intramolecular_key(py)).unwrap();
            match removed {
                Some(NoncovalentBondConstraintForm::Intramolecular(b)) => {
                    assert_eq!(b.bind(py).borrow().to_rust(), GraphIrBooleanForm::Lit(true))
                }
                _ => panic!("expected removed Intramolecular(Lit(true))"),
            }
            assert_eq!(constraints.__len__(), 0);
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_form_update_container() {
        Python::attach(|py| {
            let constraints = Py::new(py, NoncovalentBondConstraintsForm::new(py, vec![])).unwrap();
            let mut other = GraphIrNoncovalentBondConstraintsForm::new();
            other.set(GraphIrNoncovalentBondConstraintForm::intramolecular(true));
            NoncovalentBondConstraintsForm::update(
                constraints.clone_ref(py),
                py,
                NoncovalentBondConstraintsUpdate::Container(
                    Py::new(py, NoncovalentBondConstraintsForm::from_rust(other)).unwrap(),
                ),
            )
            .unwrap();
            assert_eq!(
                constraints.bind(py).borrow().intramolecular().to_rust(),
                GraphIrBooleanForm::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_form_update_entries() {
        Python::attach(|py| {
            let constraints = Py::new(py, NoncovalentBondConstraintsForm::new(py, vec![])).unwrap();
            NoncovalentBondConstraintsForm::update(
                constraints.clone_ref(py),
                py,
                NoncovalentBondConstraintsUpdate::Entries(vec![intramolecular(py, false)]),
            )
            .unwrap();
            assert_eq!(
                constraints.bind(py).borrow().intramolecular().to_rust(),
                GraphIrBooleanForm::Lit(false)
            );
        });
    }

    /// Regression: a container updating itself resolves `other` before the write borrow,
    /// so it is an idempotent no-op, not a RefCell double-borrow panic.
    #[rstest]
    fn test_noncovalent_bond_constraints_form_update_self() {
        Python::attach(|py| {
            let constraints = Py::new(
                py,
                NoncovalentBondConstraintsForm::new(py, vec![intramolecular(py, true)]),
            )
            .unwrap();
            NoncovalentBondConstraintsForm::update(
                constraints.clone_ref(py),
                py,
                NoncovalentBondConstraintsUpdate::Container(constraints.clone_ref(py)),
            )
            .unwrap();
            assert_eq!(
                constraints.bind(py).borrow().intramolecular().to_rust(),
                GraphIrBooleanForm::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_form_len_contains() {
        Python::attach(|py| {
            let constraints =
                NoncovalentBondConstraintsForm::new(py, vec![intramolecular(py, true)]);
            assert_eq!(constraints.__len__(), 1);
            assert!(constraints.__contains__(py, intramolecular_key(py)));
            let empty = NoncovalentBondConstraintsForm::new(py, vec![]);
            assert!(!empty.__contains__(py, intramolecular_key(py)));
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_form_keys_values_items() {
        Python::attach(|py| {
            let constraints =
                NoncovalentBondConstraintsForm::new(py, vec![intramolecular(py, true)]);
            let mut keys = constraints.keys(py).unwrap();
            assert_eq!(
                keys.__next__().unwrap().bind(py).borrow().to_rust(),
                GraphIrNoncovalentBondConstraintKey::Intramolecular
            );
            assert!(keys.__next__().is_none());
            let mut values = constraints.values(py).unwrap();
            assert_eq!(
                values.__next__().unwrap().bind(py).borrow().to_rust(py),
                GraphIrNoncovalentBondConstraintForm::intramolecular(true)
            );
            let mut items = constraints.items(py).unwrap();
            let (k, v) = items.__next__().unwrap();
            assert_eq!(
                k.bind(py).borrow().to_rust(),
                GraphIrNoncovalentBondConstraintKey::Intramolecular
            );
            assert_eq!(
                v.bind(py).borrow().to_rust(py),
                GraphIrNoncovalentBondConstraintForm::intramolecular(true)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_form_get() {
        Python::attach(|py| {
            let constraints =
                NoncovalentBondConstraintsForm::new(py, vec![intramolecular(py, true)]);
            let present = constraints.get(py, intramolecular_key(py), None).unwrap();
            let expected = intramolecular(py, true).into_any();
            assert!(present.bind(py).eq(expected.bind(py)).unwrap());
            // absent → None
            let empty = NoncovalentBondConstraintsForm::new(py, vec![]);
            let absent = empty.get(py, intramolecular_key(py), None).unwrap();
            assert!(absent.bind(py).is_none());
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_form_intramolecular() {
        Python::attach(|py| {
            let present = NoncovalentBondConstraintsForm::new(py, vec![intramolecular(py, true)]);
            assert_eq!(
                present.intramolecular().to_rust(),
                GraphIrBooleanForm::Lit(true)
            );
            // absent → Undetermined (non-optional accessor)
            let empty = NoncovalentBondConstraintsForm::new(py, vec![]);
            assert_eq!(
                empty.intramolecular().to_rust(),
                GraphIrBooleanForm::Undetermined
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_form_set_intramolecular() {
        Python::attach(|py| {
            let mut constraints = NoncovalentBondConstraintsForm::new(py, vec![]);
            constraints.set_intramolecular(py, BooleanLike::Lit(false));
            assert_eq!(
                constraints.intramolecular().to_rust(),
                GraphIrBooleanForm::Lit(false)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_form_getitem_error() {
        Python::attach(|py| {
            let constraints = NoncovalentBondConstraintsForm::new(py, vec![]);
            assert!(constraints.__getitem__(py, intramolecular_key(py)).is_err());
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_form_delitem_error() {
        Python::attach(|py| {
            let mut constraints = NoncovalentBondConstraintsForm::new(py, vec![]);
            assert!(constraints.__delitem__(py, intramolecular_key(py)).is_err());
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_form_asdict() {
        Python::attach(|py| {
            let constraints =
                NoncovalentBondConstraintsForm::new(py, vec![intramolecular(py, true)]);
            let dict = constraints.asdict(py).unwrap();
            let value = dict.get_item("intramolecular").unwrap().unwrap();
            let expected = into_py_variant(py, BooleanForm::Lit(true)).unwrap();
            assert!(value.eq(expected.bind(py)).unwrap());
        });
    }

    /// A standalone `NoncovalentBondForm` value pyclass (hydrogen bond, no constraints).
    fn hbond(py: Python<'_>) -> Py<NoncovalentBondForm> {
        Py::new(
            py,
            NoncovalentBondForm::from_rust(GraphIrNoncovalentBondForm::from_kind(
                GraphIrNoncovalentBondKind::HydrogenBond,
            )),
        )
        .unwrap()
    }

    #[rstest]
    fn test_noncovalent_bond_form_new() {
        Python::attach(|py| {
            let bond = NoncovalentBondForm::new(
                py,
                NoncovalentBondKindLike::Kind(NoncovalentBondKind::HydrogenBond),
                None,
            );
            assert_eq!(
                bond.to_rust().kind,
                GraphIrNoncovalentBondKindForm::Lit(GraphIrNoncovalentBondKind::HydrogenBond)
            );
            assert_eq!(bond.to_rust().constraints.len(), 0);
        });
    }

    #[rstest]
    fn test_noncovalent_bond_form_new_constraints() {
        Python::attach(|py| {
            let constraints = Py::new(
                py,
                NoncovalentBondConstraintsForm::new(py, vec![intramolecular(py, true)]),
            )
            .unwrap();
            let bond = NoncovalentBondForm::new(
                py,
                NoncovalentBondKindLike::Kind(NoncovalentBondKind::HalogenBond),
                Some(constraints),
            );
            assert_eq!(
                bond.to_rust().constraints.intramolecular(),
                GraphIrBooleanForm::Lit(true)
            );
        });
    }

    #[rstest]
    #[case("Hbd")]
    #[case("Hbd#I")]
    #[case("Hbd#I!")]
    #[case("*")]
    fn test_noncovalent_bond_form_parse(#[case] dsl: &str) {
        let bond = NoncovalentBondForm::parse(dsl).unwrap();
        assert_eq!(bond.__str__(), dsl);
        assert_eq!(
            bond.__repr__(),
            format!("NoncovalentBondForm.parse('{dsl}')")
        );
    }

    #[rstest]
    fn test_noncovalent_bond_form_parse_error() {
        assert!(NoncovalentBondForm::parse("z").is_err());
    }

    #[rstest]
    fn test_noncovalent_bond_form_kind() {
        Python::attach(|py| {
            let mut bond = NoncovalentBondForm::from_rust(GraphIrNoncovalentBondForm::from_kind(
                GraphIrNoncovalentBondKind::HydrogenBond,
            ));
            assert_eq!(
                bond.kind().to_rust(),
                GraphIrNoncovalentBondKindForm::Lit(GraphIrNoncovalentBondKind::HydrogenBond)
            );
            bond.set_kind(
                py,
                NoncovalentBondKindLike::Kind(NoncovalentBondKind::Ionic),
            )
            .unwrap();
            assert_eq!(
                bond.kind().to_rust(),
                GraphIrNoncovalentBondKindForm::Lit(GraphIrNoncovalentBondKind::Ionic)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_form_set_constraints() {
        Python::attach(|py| {
            let bond = hbond(py);
            let constraints = Py::new(
                py,
                NoncovalentBondConstraintsForm::new(py, vec![intramolecular(py, true)]),
            )
            .unwrap();
            NoncovalentBondForm::set_constraints(
                bond.clone_ref(py),
                py,
                NoncovalentBondConstraintsLike::Container(constraints),
            )
            .unwrap();
            assert_eq!(
                bond.bind(py)
                    .borrow()
                    .to_rust()
                    .constraints
                    .intramolecular(),
                GraphIrBooleanForm::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_form_set_constraints_from_view() {
        Python::attach(|py| {
            // source bond carrying a constraint, exposed as a live view
            let source = Py::new(
                py,
                NoncovalentBondForm::from_rust(
                    GraphIrNoncovalentBondForm::from_kind(GraphIrNoncovalentBondKind::HydrogenBond)
                        .with_constraint(GraphIrNoncovalentBondConstraintForm::intramolecular(
                            true,
                        )),
                ),
            )
            .unwrap();
            let view = NoncovalentBondForm::constraints(source);
            let dest = hbond(py);
            NoncovalentBondForm::set_constraints(
                dest.clone_ref(py),
                py,
                NoncovalentBondConstraintsLike::View(Py::new(py, view).unwrap()),
            )
            .unwrap();
            assert_eq!(
                dest.bind(py)
                    .borrow()
                    .to_rust()
                    .constraints
                    .intramolecular(),
                GraphIrBooleanForm::Lit(true)
            );
        });
    }

    /// Regression: assigning a bond's own constraints view back to it snapshots before
    /// the write borrow, so it is a no-op, not a RefCell double-borrow panic
    /// (`bond.constraints = bond.constraints`).
    #[rstest]
    fn test_noncovalent_bond_form_set_constraints_self() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                NoncovalentBondForm::from_rust(
                    GraphIrNoncovalentBondForm::from_kind(GraphIrNoncovalentBondKind::HydrogenBond)
                        .with_constraint(GraphIrNoncovalentBondConstraintForm::intramolecular(
                            true,
                        )),
                ),
            )
            .unwrap();
            let own_view = NoncovalentBondForm::constraints(bond.clone_ref(py));
            NoncovalentBondForm::set_constraints(
                bond.clone_ref(py),
                py,
                NoncovalentBondConstraintsLike::View(Py::new(py, own_view).unwrap()),
            )
            .unwrap();
            assert_eq!(
                bond.bind(py)
                    .borrow()
                    .to_rust()
                    .constraints
                    .intramolecular(),
                GraphIrBooleanForm::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_form_constraints_write_through() {
        Python::attach(|py| {
            let bond = hbond(py);
            let view = NoncovalentBondForm::constraints(bond.clone_ref(py));
            view.set(py, intramolecular(py, true)).unwrap();
            // the write hit the standalone bond, not a copy
            assert_eq!(
                bond.bind(py)
                    .borrow()
                    .to_rust()
                    .constraints
                    .intramolecular(),
                GraphIrBooleanForm::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_form_asdict() {
        Python::attach(|py| {
            let bond = NoncovalentBondForm::from_rust(
                GraphIrNoncovalentBondForm::from_kind(GraphIrNoncovalentBondKind::HydrogenBond)
                    .with_constraint(GraphIrNoncovalentBondConstraintForm::intramolecular(true)),
            );
            let dict = bond.asdict(py).unwrap();
            let kind = dict.get_item("kind").unwrap().unwrap();
            let expected_kind = into_py_variant(
                py,
                NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
            )
            .unwrap();
            assert!(kind.eq(expected_kind.bind(py)).unwrap());
            let constraints = dict.get_item("constraints").unwrap().unwrap();
            assert_eq!(constraints.len().unwrap(), 1);
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_view_set() {
        Python::attach(|py| {
            let bond = hbond(py);
            let view = NoncovalentBondConstraintsView {
                backing: NoncovalentBondConstraintsBacking::Noncovalent(bond.clone_ref(py)),
            };
            view.set(py, intramolecular(py, true)).unwrap();
            // a fresh view proves the write hit the standalone bond, not a copy
            let fresh = NoncovalentBondConstraintsView {
                backing: NoncovalentBondConstraintsBacking::Noncovalent(bond),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 1);
            assert_eq!(
                fresh.intramolecular(py).unwrap().to_rust(),
                GraphIrBooleanForm::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_view_pop() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                NoncovalentBondForm::from_rust(
                    GraphIrNoncovalentBondForm::from_kind(GraphIrNoncovalentBondKind::HydrogenBond)
                        .with_constraint(GraphIrNoncovalentBondConstraintForm::intramolecular(
                            true,
                        )),
                ),
            )
            .unwrap();
            let view = NoncovalentBondConstraintsView {
                backing: NoncovalentBondConstraintsBacking::Noncovalent(bond.clone_ref(py)),
            };
            let removed = view.pop(py, intramolecular_key(py)).unwrap();
            match removed {
                Some(NoncovalentBondConstraintForm::Intramolecular(b)) => {
                    assert_eq!(b.bind(py).borrow().to_rust(), GraphIrBooleanForm::Lit(true))
                }
                _ => panic!("expected removed Intramolecular(Lit(true))"),
            }
            let fresh = NoncovalentBondConstraintsView {
                backing: NoncovalentBondConstraintsBacking::Noncovalent(bond),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 0);
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_view_update() {
        Python::attach(|py| {
            let bond = hbond(py);
            let view = NoncovalentBondConstraintsView {
                backing: NoncovalentBondConstraintsBacking::Noncovalent(bond.clone_ref(py)),
            };
            view.update(
                py,
                NoncovalentBondConstraintsUpdate::Entries(vec![intramolecular(py, false)]),
            )
            .unwrap();
            assert_eq!(
                bond.bind(py)
                    .borrow()
                    .to_rust()
                    .constraints
                    .intramolecular(),
                GraphIrBooleanForm::Lit(false)
            );
        });
    }

    /// Regression: a view updating from a view over the same bond resolves `other`
    /// before the write borrow, so it is an idempotent no-op, not a double-borrow panic
    /// (`bond.constraints.update(bond.constraints)`).
    #[rstest]
    fn test_noncovalent_bond_constraints_view_update_self() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                NoncovalentBondForm::from_rust(
                    GraphIrNoncovalentBondForm::from_kind(GraphIrNoncovalentBondKind::HydrogenBond)
                        .with_constraint(GraphIrNoncovalentBondConstraintForm::intramolecular(
                            true,
                        )),
                ),
            )
            .unwrap();
            let view = NoncovalentBondConstraintsView {
                backing: NoncovalentBondConstraintsBacking::Noncovalent(bond.clone_ref(py)),
            };
            let other = NoncovalentBondConstraintsView {
                backing: NoncovalentBondConstraintsBacking::Noncovalent(bond.clone_ref(py)),
            };
            view.update(
                py,
                NoncovalentBondConstraintsUpdate::View(Py::new(py, other).unwrap()),
            )
            .unwrap();
            assert_eq!(
                bond.bind(py)
                    .borrow()
                    .to_rust()
                    .constraints
                    .intramolecular(),
                GraphIrBooleanForm::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_constraints_view_set_intramolecular() {
        Python::attach(|py| {
            let bond = hbond(py);
            let view = NoncovalentBondConstraintsView {
                backing: NoncovalentBondConstraintsBacking::Noncovalent(bond.clone_ref(py)),
            };
            view.set_intramolecular(py, BooleanLike::Lit(true)).unwrap();
            assert_eq!(
                bond.bind(py)
                    .borrow()
                    .to_rust()
                    .constraints
                    .intramolecular(),
                GraphIrBooleanForm::Lit(true)
            );
        });
    }

    /// A molecule of two oxygens with one hydrogen bond over atoms (0, 1), noncovalent id 0.
    fn molecule_with_hbond(py: Python<'_>) -> Py<MoleculeAst> {
        let molecule = GraphIrMolecule::from_entries(MoleculeEntries {
            atoms: vec![
                GraphIrAtomForm::from_element(ChemElement::O),
                GraphIrAtomForm::from_element(ChemElement::O),
            ],
            noncovalent: vec![(
                GraphIrAtomId(0),
                GraphIrAtomId(1),
                GraphIrNoncovalentBondForm::from_kind(GraphIrNoncovalentBondKind::HydrogenBond),
            )],
            ..Default::default()
        });
        Py::new(py, MoleculeAst::from_rust(molecule)).unwrap()
    }

    #[rstest]
    fn test_noncovalent_bond_view_id_atom_ids() {
        Python::attach(|py| {
            let owner = molecule_with_hbond(py);
            let view = NoncovalentBondView {
                owner,
                id: GraphIrNoncovalentBondId(0),
            };
            assert_eq!(view.id(), 0);
            assert_eq!(view.atom_ids(py).unwrap(), (0, 1));
            assert_eq!(view.__repr__(), "NoncovalentBondView(id=0)");
        });
    }

    #[rstest]
    fn test_noncovalent_bond_view_kind() {
        Python::attach(|py| {
            let owner = molecule_with_hbond(py);
            let view = NoncovalentBondView {
                owner: owner.clone_ref(py),
                id: GraphIrNoncovalentBondId(0),
            };
            assert_eq!(
                view.kind(py).unwrap().to_rust(),
                GraphIrNoncovalentBondKindForm::Lit(GraphIrNoncovalentBondKind::HydrogenBond)
            );
            view.set_kind(
                py,
                NoncovalentBondKindLike::Kind(NoncovalentBondKind::Ionic),
            );
            // a fresh read proves the write hit the molecule
            let fresh = NoncovalentBondView {
                owner,
                id: GraphIrNoncovalentBondId(0),
            };
            assert_eq!(
                fresh.kind(py).unwrap().to_rust(),
                GraphIrNoncovalentBondKindForm::Lit(GraphIrNoncovalentBondKind::Ionic)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_view_constraints_write_through() {
        Python::attach(|py| {
            let owner = molecule_with_hbond(py);
            let view = NoncovalentBondView {
                owner: owner.clone_ref(py),
                id: GraphIrNoncovalentBondId(0),
            };
            // the constraints handle is molecule-backed; a write goes through to the bond
            view.constraints(py)
                .set_intramolecular(py, BooleanLike::Lit(true))
                .unwrap();
            assert_eq!(
                owner
                    .bind(py)
                    .borrow()
                    .to_rust()
                    .noncovalent_bond(GraphIrNoncovalentBondId(0))
                    .attributes
                    .constraints
                    .intramolecular(),
                GraphIrBooleanForm::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_view_set_constraints() {
        Python::attach(|py| {
            let owner = molecule_with_hbond(py);
            let view = NoncovalentBondView {
                owner: owner.clone_ref(py),
                id: GraphIrNoncovalentBondId(0),
            };
            let constraints = Py::new(
                py,
                NoncovalentBondConstraintsForm::new(py, vec![intramolecular(py, false)]),
            )
            .unwrap();
            view.set_constraints(py, NoncovalentBondConstraintsLike::Container(constraints))
                .unwrap();
            assert_eq!(
                owner
                    .bind(py)
                    .borrow()
                    .to_rust()
                    .noncovalent_bond(GraphIrNoncovalentBondId(0))
                    .attributes
                    .constraints
                    .intramolecular(),
                GraphIrBooleanForm::Lit(false)
            );
        });
    }

    /// Regression: `mol.noncovalent_bonds[i].constraints.update(same-view)` resolves before
    /// the molecule write borrow, so the molecule-backed self-alias is a no-op, not a panic.
    #[rstest]
    fn test_noncovalent_bond_view_constraints_update_self() {
        Python::attach(|py| {
            let owner = molecule_with_hbond(py);
            owner
                .bind(py)
                .borrow_mut()
                .to_rust_mut()
                .noncovalent_bond_mut(GraphIrNoncovalentBondId(0))
                .attributes
                .constraints
                .set(GraphIrNoncovalentBondConstraintForm::intramolecular(true));
            let view = NoncovalentBondConstraintsView {
                backing: NoncovalentBondConstraintsBacking::Molecule {
                    owner: owner.clone_ref(py),
                    id: GraphIrNoncovalentBondId(0),
                },
            };
            let other = Py::new(
                py,
                NoncovalentBondConstraintsView {
                    backing: NoncovalentBondConstraintsBacking::Molecule {
                        owner: owner.clone_ref(py),
                        id: GraphIrNoncovalentBondId(0),
                    },
                },
            )
            .unwrap();
            view.update(py, NoncovalentBondConstraintsUpdate::View(other))
                .unwrap();
            assert_eq!(
                owner
                    .bind(py)
                    .borrow()
                    .to_rust()
                    .noncovalent_bond(GraphIrNoncovalentBondId(0))
                    .attributes
                    .constraints
                    .intramolecular(),
                GraphIrBooleanForm::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_noncovalent_bond_view_asdict() {
        Python::attach(|py| {
            let owner = molecule_with_hbond(py);
            let view = NoncovalentBondView {
                owner,
                id: GraphIrNoncovalentBondId(0),
            };
            let dict = view.asdict(py).unwrap();
            assert_eq!(
                dict.keys()
                    .iter()
                    .map(|k| k.extract::<String>().unwrap())
                    .collect::<Vec<_>>(),
                vec!["kind".to_string(), "constraints".to_string()]
            );
            let kind = dict.get_item("kind").unwrap().unwrap();
            let expected = into_py_variant(
                py,
                NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
            )
            .unwrap();
            assert!(kind.eq(expected.bind(py)).unwrap());
        });
    }

    /// Three atoms, one hydrogen bond over (0, 1), atom 2 isolated. For the collection
    /// negative cases (`of` / `incident` with no bond).
    fn molecule_with_hbond_and_isolated(py: Python<'_>) -> Py<MoleculeAst> {
        let molecule = GraphIrMolecule::from_entries(MoleculeEntries {
            atoms: vec![
                GraphIrAtomForm::from_element(ChemElement::O),
                GraphIrAtomForm::from_element(ChemElement::O),
                GraphIrAtomForm::from_element(ChemElement::O),
            ],
            noncovalent: vec![(
                GraphIrAtomId(0),
                GraphIrAtomId(1),
                GraphIrNoncovalentBondForm::from_kind(GraphIrNoncovalentBondKind::HydrogenBond),
            )],
            ..Default::default()
        });
        Py::new(py, MoleculeAst::from_rust(molecule)).unwrap()
    }

    #[rstest]
    fn test_noncovalent_bond_views_len_getitem() {
        Python::attach(|py| {
            let views = NoncovalentBondViews {
                owner: molecule_with_hbond(py),
            };
            assert_eq!(views.__len__(py), 1);
            assert_eq!(views.__repr__(py), "NoncovalentBondViews(len=1)");
            assert_eq!(views.__getitem__(py, 0).unwrap().id(), 0);
            // negative index counts from the end
            assert_eq!(views.__getitem__(py, -1).unwrap().id(), 0);
            assert!(views.__getitem__(py, 1).is_err());
            assert!(views.__getitem__(py, -2).is_err());
        });
    }

    #[rstest]
    fn test_noncovalent_bond_views_setitem() {
        Python::attach(|py| {
            let owner = molecule_with_hbond(py);
            let views = NoncovalentBondViews {
                owner: owner.clone_ref(py),
            };
            let replacement = Py::new(
                py,
                NoncovalentBondForm::from_rust(GraphIrNoncovalentBondForm::from_kind(
                    GraphIrNoncovalentBondKind::Ionic,
                )),
            )
            .unwrap();
            views
                .__setitem__(py, 0, replacement.bind(py).borrow())
                .unwrap();
            let view = views.__getitem__(py, 0).unwrap();
            // value replaced, endpoints preserved
            assert_eq!(
                view.kind(py).unwrap().to_rust(),
                GraphIrNoncovalentBondKindForm::Lit(GraphIrNoncovalentBondKind::Ionic)
            );
            assert_eq!(view.atom_ids(py).unwrap(), (0, 1));
        });
    }

    #[rstest]
    fn test_noncovalent_bond_views_setitem_out_of_range() {
        Python::attach(|py| {
            let views = NoncovalentBondViews {
                owner: molecule_with_hbond(py),
            };
            let bond = Py::new(
                py,
                NoncovalentBondForm::from_rust(GraphIrNoncovalentBondForm::from_kind(
                    GraphIrNoncovalentBondKind::Ionic,
                )),
            )
            .unwrap();
            assert!(views.__setitem__(py, 5, bond.bind(py).borrow()).is_err());
        });
    }

    #[rstest]
    fn test_noncovalent_bond_views_iter() {
        Python::attach(|py| {
            let views = NoncovalentBondViews {
                owner: molecule_with_hbond(py),
            };
            let mut iter = views.__iter__(py);
            assert_eq!(iter.__next__(py).unwrap().id(), 0);
            assert!(iter.__next__(py).is_none());
        });
    }

    #[rstest]
    fn test_noncovalent_bond_views_of() {
        Python::attach(|py| {
            let views = NoncovalentBondViews {
                owner: molecule_with_hbond_and_isolated(py),
            };
            // unordered pair — both orders find the same bond
            assert_eq!(views.of(py, 0, 1).unwrap().id(), 0);
            assert_eq!(views.of(py, 1, 0).unwrap().id(), 0);
            // no bond between 0 and the isolated atom 2
            assert!(views.of(py, 0, 2).is_none());
        });
    }

    #[rstest]
    fn test_noncovalent_bond_views_incident() {
        Python::attach(|py| {
            let views = NoncovalentBondViews {
                owner: molecule_with_hbond_and_isolated(py),
            };
            assert_eq!(
                views
                    .incident(py, 0)
                    .iter()
                    .map(|v| v.id())
                    .collect::<Vec<_>>(),
                vec![0]
            );
            assert!(views.incident(py, 2).is_empty());
        });
    }
}
