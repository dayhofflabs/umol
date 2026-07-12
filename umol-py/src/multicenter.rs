//! Multicenter bond value type and multicenter-bond-constraint surface mirroring
//! `umol_ast::ast`: `MulticenterBondAst`, the `MulticenterBondConstraintAst` enum, the
//! `MulticenterBondConstraintsAst` container, and the `MulticenterBondConstraintsView`
//! live handle. A multicenter bond is a single unordered set of member atoms; the
//! value carries a positional per-atom `electrons` vector plus charge, spin, and
//! constraints. The member atoms are the participants of the owning molecule's
//! multicenter relation, so they are topology (the view half) rather than value.

use std::str::FromStr;
use std::vec::IntoIter;

use pyo3::exceptions::{PyIndexError, PyKeyError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict, PyTuple};
use umol_ast::ast::{
    AtomId as AstAtomId, MoleculeAst as AstMoleculeAst,
    MulticenterBondAst as AstMulticenterBondAst,
    MulticenterBondConstraintAst as AstMulticenterBondConstraintAst,
    MulticenterBondConstraintKey as AstMulticenterBondConstraintKey,
    MulticenterBondConstraintsAst as AstMulticenterBondConstraintsAst,
    MulticenterBondId as AstMulticenterBondId, MulticenterBondView as AstMulticenterBondView,
};

use crate::atom::SpinStateAst;
use crate::convert::{hash_ast, into_py_variant, variant_repr};
use crate::electrons::{ElectronCountsArg, ElectronCountsAst};
use crate::error::parse_error;
use crate::molecule::MoleculeAst;
use crate::value::{ValueArg, ValueAst};

/// A multicenter bond: a positional per-member-atom `electrons` vector, charge,
/// spin, and multicenter-bond-scope constraints. The member atoms are the
/// participants of the owning molecule's multicenter relation (the view half); the
/// `electrons` vector is positional, aligned to that atom order.
#[pyclass(eq)]
#[derive(PartialEq)]
pub struct MulticenterBondAst(AstMulticenterBondAst);

#[pymethods]
impl MulticenterBondAst {
    /// Construct from an electron-count vector — a `list[int]` or an
    /// `ElectronCountsAst` — optionally setting fields.
    #[new]
    #[pyo3(signature = (electrons, *, charge=None, spin=None, constraints=None))]
    fn new(
        py: Python<'_>,
        electrons: ElectronCountsArg,
        charge: Option<ValueArg>,
        spin: Option<PyRef<'_, SpinStateAst>>,
        constraints: Option<Py<MulticenterBondConstraintsAst>>,
    ) -> Self {
        let mut bond = AstMulticenterBondAst::new(electrons.to_ast(py));
        if let Some(charge) = charge {
            bond = bond.with_charge(charge.to_ast(py));
        }
        if let Some(spin) = spin {
            bond = bond.with_spin(spin.to_ast(py));
        }
        if let Some(constraints) = constraints {
            bond.constraints = constraints.bind(py).borrow().inner().clone();
        }
        MulticenterBondAst(bond)
    }

    /// Parse a multicenter-bond-DSL string (e.g. `"[1,1,1]#e6"`) into a `MulticenterBondAst`.
    #[staticmethod]
    fn parse(s: &str) -> PyResult<Self> {
        AstMulticenterBondAst::from_str(s)
            .map(Self)
            .map_err(parse_error)
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("MulticenterBondAst.parse('{}')", self.0)
    }

    /// The per-member-atom electron counts (positional, aligned to `atom_ids`).
    #[getter]
    fn electrons(&self) -> ElectronCountsAst {
        ElectronCountsAst::from_ast(&self.0.electrons)
    }

    #[setter]
    fn set_electrons(&mut self, py: Python<'_>, value: ElectronCountsArg) {
        self.0.electrons = value.to_ast(py);
    }

    #[getter]
    fn charge(&self, py: Python<'_>) -> PyResult<ValueAst> {
        ValueAst::from_ast(py, &self.0.charge)
    }

    #[setter]
    fn set_charge(&mut self, py: Python<'_>, value: ValueArg) {
        self.0.charge = value.to_ast(py);
    }

    #[getter]
    fn spin(&self, py: Python<'_>) -> PyResult<SpinStateAst> {
        SpinStateAst::from_ast(py, &self.0.spin)
    }

    #[setter]
    fn set_spin(&mut self, py: Python<'_>, value: PyRef<'_, SpinStateAst>) {
        self.0.spin = value.to_ast(py);
    }

    /// The bond's constraints as a live handle onto this bond: reads borrow the
    /// current state, mutators write through to the bond in place.
    #[getter]
    fn constraints(slf: Py<Self>) -> MulticenterBondConstraintsView {
        MulticenterBondConstraintsView {
            backing: MulticenterBondConstraintsBacking::MulticenterBond(slf),
        }
    }

    /// Replace the whole constraint set (wipe-and-set) from a value container or a
    /// live view.
    #[setter]
    fn set_constraints(
        &mut self,
        py: Python<'_>,
        value: MulticenterBondConstraintsArg,
    ) -> PyResult<()> {
        self.0.constraints = value.to_ast(py)?;
        Ok(())
    }

    /// The fields as a dict keyed by field name; values are the field mirrors.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("electrons", self.electrons())?;
        dict.set_item("charge", self.charge(py)?)?;
        dict.set_item("spin", self.spin(py)?)?;
        dict.set_item(
            "constraints",
            multicenter_bond_constraints_asdict(py, &self.0.constraints)?,
        )?;
        Ok(dict)
    }
}

impl MulticenterBondAst {
    /// The wrapped AST bond — read access for the bond-backed constraints view.
    pub(crate) fn inner(&self) -> &AstMulticenterBondAst {
        &self.0
    }

    /// Mutable access to the wrapped AST bond — write access for the bond-backed
    /// constraints view.
    pub(crate) fn inner_mut(&mut self) -> &mut AstMulticenterBondAst {
        &mut self.0
    }

    /// Wrap an AST bond (the hold-the-value `from_inner` bridge, paired with
    /// `inner`). Test-only — in-crate construction wraps `MulticenterBondAst(..)` directly.
    #[cfg(test)]
    pub(crate) fn from_inner(bond: AstMulticenterBondAst) -> Self {
        MulticenterBondAst(bond)
    }
}

/// A view of one multicenter bond within a molecule: a handle to the molecule plus
/// the bond's index. Field reads rebuild the transient Rust view; the molecule is
/// never copied. The member atom indices are read-only topology; the electrons,
/// charge, spin, and constraints are the mutable bond value.
#[pyclass]
pub struct MulticenterBondView {
    owner: Py<MoleculeAst>,
    id: AstMulticenterBondId,
}

impl MulticenterBondView {
    fn multicenter_bond<'a>(
        &self,
        molecule: &'a AstMoleculeAst,
    ) -> PyResult<AstMulticenterBondView<'a>> {
        molecule
            .multicenter_bonds()
            .get(self.id)
            .ok_or_else(|| PyIndexError::new_err("multicenter bond id out of range"))
    }
}

#[pymethods]
impl MulticenterBondView {
    #[getter]
    fn id(&self) -> u32 {
        self.id.0
    }

    /// The member atom indices (read-only — participants are topology, not part of
    /// the bond value).
    #[getter]
    fn atom_ids<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let molecule = self.owner.bind(py).borrow();
        let atom_ids: Vec<u32> = self
            .multicenter_bond(molecule.inner())?
            .atom_ids()
            .map(|atom| atom.0)
            .collect();
        PyTuple::new(py, atom_ids)
    }

    fn __repr__(&self) -> String {
        format!("MulticenterBondView(id={})", self.id.0)
    }

    /// The per-member-atom electron counts (positional, aligned to `atom_ids`).
    #[getter]
    fn electrons(&self, py: Python<'_>) -> PyResult<ElectronCountsAst> {
        let molecule = self.owner.bind(py).borrow();
        Ok(ElectronCountsAst::from_ast(
            &self.multicenter_bond(molecule.inner())?.ast.electrons,
        ))
    }

    #[setter]
    fn set_electrons(&self, py: Python<'_>, value: ElectronCountsArg) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .multicenter_bond_mut(self.id)
            .ast
            .electrons = value.to_ast(py);
    }

    #[getter]
    fn charge(&self, py: Python<'_>) -> PyResult<ValueAst> {
        let molecule = self.owner.bind(py).borrow();
        ValueAst::from_ast(py, &self.multicenter_bond(molecule.inner())?.ast.charge)
    }

    #[setter]
    fn set_charge(&self, py: Python<'_>, value: ValueArg) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .multicenter_bond_mut(self.id)
            .ast
            .charge = value.to_ast(py);
    }

    #[getter]
    fn spin(&self, py: Python<'_>) -> PyResult<SpinStateAst> {
        let molecule = self.owner.bind(py).borrow();
        SpinStateAst::from_ast(py, &self.multicenter_bond(molecule.inner())?.ast.spin)
    }

    #[setter]
    fn set_spin(&self, py: Python<'_>, value: PyRef<'_, SpinStateAst>) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .multicenter_bond_mut(self.id)
            .ast
            .spin = value.to_ast(py);
    }

    /// The bond's constraints as a live handle onto the molecule: reads borrow the
    /// current state, mutators write through to the bond in place.
    #[getter]
    fn constraints(&self, py: Python<'_>) -> MulticenterBondConstraintsView {
        MulticenterBondConstraintsView {
            backing: MulticenterBondConstraintsBacking::Molecule {
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
        value: MulticenterBondConstraintsArg,
    ) -> PyResult<()> {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .multicenter_bond_mut(self.id)
            .ast
            .constraints = value.to_ast(py)?;
        Ok(())
    }

    /// The value fields as a dict keyed by field name; values are the field mirrors —
    /// symmetric with `MulticenterBondAst.asdict`, read through the view.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let molecule = self.owner.bind(py).borrow();
        let bond = self.multicenter_bond(molecule.inner())?.ast;
        let dict = PyDict::new(py);
        dict.set_item("electrons", ElectronCountsAst::from_ast(&bond.electrons))?;
        dict.set_item("charge", ValueAst::from_ast(py, &bond.charge)?)?;
        dict.set_item("spin", SpinStateAst::from_ast(py, &bond.spin)?)?;
        dict.set_item(
            "constraints",
            multicenter_bond_constraints_asdict(py, &bond.constraints)?,
        )?;
        Ok(dict)
    }
}

/// Resolve a possibly-negative Python index (negative counts from the end) into an
/// existing multicenter bond id, or `IndexError`. `MulticenterBondId` is `RelationId`-
/// backed but contiguous for fresh molecules, so integer positions address it directly.
fn resolve_multicenter_bond_index(
    molecule: &AstMoleculeAst,
    index: isize,
) -> PyResult<AstMulticenterBondId> {
    let count = molecule.multicenter_bonds().count();
    let resolved = if index < 0 {
        index + count as isize
    } else {
        index
    };
    if resolved < 0 {
        return Err(PyIndexError::new_err("multicenter bond id out of range"));
    }
    let id = AstMulticenterBondId(resolved as u32);
    if molecule.multicenter_bonds().contains(id) {
        Ok(id)
    } else {
        Err(PyIndexError::new_err("multicenter bond id out of range"))
    }
}

/// The multicenter bonds of a molecule, indexed by integer position.
#[pyclass]
pub struct MulticenterBondViews {
    owner: Py<MoleculeAst>,
}

#[pymethods]
impl MulticenterBondViews {
    fn __len__(&self, py: Python<'_>) -> usize {
        self.owner
            .bind(py)
            .borrow()
            .inner()
            .multicenter_bonds()
            .count()
    }

    fn __repr__(&self, py: Python<'_>) -> String {
        format!(
            "MulticenterBondViews(len={})",
            self.owner
                .bind(py)
                .borrow()
                .inner()
                .multicenter_bonds()
                .count()
        )
    }

    fn __getitem__(&self, py: Python<'_>, index: isize) -> PyResult<MulticenterBondView> {
        let molecule = self.owner.bind(py).borrow();
        let id = resolve_multicenter_bond_index(molecule.inner(), index)?;
        Ok(MulticenterBondView {
            owner: self.owner.clone_ref(py),
            id,
        })
    }

    /// Replace the whole multicenter bond value at `index` in place (members unchanged).
    fn __setitem__(
        &self,
        py: Python<'_>,
        index: isize,
        bond: PyRef<'_, MulticenterBondAst>,
    ) -> PyResult<()> {
        let mut molecule = self.owner.borrow_mut(py);
        let id = resolve_multicenter_bond_index(molecule.inner(), index)?;
        *molecule.inner_mut().multicenter_bond_mut(id).ast = bond.inner().clone();
        Ok(())
    }

    /// The multicenter bond whose member atom set equals `atoms`, or `None`.
    fn connecting(&self, py: Python<'_>, atoms: Vec<u32>) -> Option<MulticenterBondView> {
        let molecule = self.owner.bind(py).borrow();
        molecule
            .inner()
            .multicenter_bonds()
            .connecting_id(atoms.into_iter().map(AstAtomId))
            .map(|id| MulticenterBondView {
                owner: self.owner.clone_ref(py),
                id,
            })
    }

    /// The multicenter bonds `atom` is a member of.
    fn incident(&self, py: Python<'_>, atom: u32) -> Vec<MulticenterBondView> {
        let molecule = self.owner.bind(py).borrow();
        molecule
            .inner()
            .multicenter_bonds()
            .incident_ids(AstAtomId(atom))
            .map(|id| MulticenterBondView {
                owner: self.owner.clone_ref(py),
                id,
            })
            .collect()
    }

    fn __iter__(&self, py: Python<'_>) -> MulticenterBondViewIter {
        let ids = self
            .owner
            .bind(py)
            .borrow()
            .inner()
            .multicenter_bonds()
            .ids()
            .collect::<Vec<_>>();
        MulticenterBondViewIter {
            owner: self.owner.clone_ref(py),
            ids: ids.into_iter(),
        }
    }
}

impl MulticenterBondViews {
    /// Build the multicenter-bond-views handle for `owner` (the `.multicenter_bonds` accessor).
    pub(crate) fn new(owner: Py<MoleculeAst>) -> MulticenterBondViews {
        MulticenterBondViews { owner }
    }
}

#[pyclass]
struct MulticenterBondViewIter {
    owner: Py<MoleculeAst>,
    ids: IntoIter<AstMulticenterBondId>,
}

#[pymethods]
impl MulticenterBondViewIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> Option<MulticenterBondView> {
        self.ids.next().map(|id| MulticenterBondView {
            owner: self.owner.clone_ref(py),
            id,
        })
    }
}

/// The key (identity) of a multicenter-bond constraint, for keyed lookup. The
/// single key `ElectronCount` is the bare discriminant (no sub-key).
#[pyclass]
pub enum MulticenterBondConstraintKey {
    ElectronCount(),
}

#[pymethods]
impl MulticenterBondConstraintKey {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_ast() == other.to_ast()
    }

    fn __hash__(&self) -> u64 {
        hash_ast(&self.to_ast())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            MulticenterBondConstraintKey::ElectronCount() => ("ElectronCount", 0),
        };
        variant_repr(
            slf.bind(py).as_any(),
            "MulticenterBondConstraintKey",
            variant,
            arity,
        )
    }
}

impl MulticenterBondConstraintKey {
    pub(crate) fn from_ast(ast: &AstMulticenterBondConstraintKey) -> Self {
        match ast {
            AstMulticenterBondConstraintKey::ElectronCount => Self::ElectronCount(),
        }
    }

    pub(crate) fn to_ast(&self) -> AstMulticenterBondConstraintKey {
        match self {
            Self::ElectronCount() => AstMulticenterBondConstraintKey::ElectronCount,
        }
    }
}

/// A multicenter-bond-scope constraint: the asserted total electron count of the
/// bond (cross-checked against `sum(MulticenterBondAst::electrons)`).
#[pyclass]
pub enum MulticenterBondConstraintAst {
    ElectronCount(Py<ValueAst>),
}

#[pymethods]
impl MulticenterBondConstraintAst {
    /// The constraint's key (identity).
    #[getter]
    fn key(&self, py: Python<'_>) -> MulticenterBondConstraintKey {
        MulticenterBondConstraintKey::from_ast(&self.to_ast(py).key())
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let variant = match &*slf.bind(py).borrow() {
            MulticenterBondConstraintAst::ElectronCount(_) => "ElectronCount",
        };
        variant_repr(
            slf.bind(py).as_any(),
            "MulticenterBondConstraintAst",
            variant,
            1,
        )
    }
}

impl MulticenterBondConstraintAst {
    pub(crate) fn from_ast(
        py: Python<'_>,
        ast: &AstMulticenterBondConstraintAst,
    ) -> PyResult<Self> {
        Ok(match ast {
            AstMulticenterBondConstraintAst::ElectronCount(v) => {
                Self::ElectronCount(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstMulticenterBondConstraintAst {
        match self {
            Self::ElectronCount(v) => {
                AstMulticenterBondConstraintAst::ElectronCount(v.bind(py).borrow().to_ast(py))
            }
        }
    }
}

/// The argument to `update`: another constraint container (value or live view) or an
/// iterable of `MulticenterBondConstraintAst` (each `set`, last-wins).
#[derive(FromPyObject)]
enum MulticenterBondConstraintsUpdate {
    Container(Py<MulticenterBondConstraintsAst>),
    View(Py<MulticenterBondConstraintsView>),
    Entries(Vec<Py<MulticenterBondConstraintAst>>),
}

impl MulticenterBondConstraintsUpdate {
    /// Overlay this update onto `target` in place.
    fn apply(&self, py: Python<'_>, target: &mut AstMulticenterBondConstraintsAst) -> PyResult<()> {
        match self {
            MulticenterBondConstraintsUpdate::Container(c) => {
                target.update(c.bind(py).borrow().inner())
            }
            MulticenterBondConstraintsUpdate::View(v) => {
                let snapshot = v.bind(py).borrow().read(py, |cs| Ok(cs.clone()))?;
                target.update(&snapshot);
            }
            MulticenterBondConstraintsUpdate::Entries(entries) => {
                for entry in entries {
                    target.set(entry.bind(py).borrow().to_ast(py));
                }
            }
        }
        Ok(())
    }
}

/// A whole-container argument that snapshots either a value container or a live view
/// — for the multicenter bond `constraints` setter, which accepts either.
#[derive(FromPyObject)]
enum MulticenterBondConstraintsArg {
    Container(Py<MulticenterBondConstraintsAst>),
    View(Py<MulticenterBondConstraintsView>),
}

impl MulticenterBondConstraintsArg {
    fn to_ast(&self, py: Python<'_>) -> PyResult<AstMulticenterBondConstraintsAst> {
        match self {
            MulticenterBondConstraintsArg::Container(c) => Ok(c.bind(py).borrow().inner().clone()),
            MulticenterBondConstraintsArg::View(v) => {
                v.bind(py).borrow().read(py, |cs| Ok(cs.clone()))
            }
        }
    }
}

/// The multicenter-bond-scope constraints on a multicenter bond, in kind-sorted order.
/// Mutable, hence value-equal but unhashable (matching `MulticenterBondAst`).
#[pyclass(eq)]
#[derive(PartialEq)]
pub struct MulticenterBondConstraintsAst(AstMulticenterBondConstraintsAst);

#[pymethods]
impl MulticenterBondConstraintsAst {
    /// Build from a sequence of constraints (a later entry of the same key replaces
    /// an earlier one, last-wins).
    #[new]
    fn new(py: Python<'_>, entries: Vec<Py<MulticenterBondConstraintAst>>) -> Self {
        let mut constraints = AstMulticenterBondConstraintsAst::new();
        constraints.extend(
            entries
                .into_iter()
                .map(|entry| entry.bind(py).borrow().to_ast(py)),
        );
        MulticenterBondConstraintsAst(constraints)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let mut parts = Vec::with_capacity(self.0.len());
        for entry in self.0.iter() {
            let mirror = into_py_variant(py, MulticenterBondConstraintAst::from_ast(py, entry)?)?;
            parts.push(mirror.bind(py).as_any().repr()?.extract::<String>()?);
        }
        Ok(format!(
            "MulticenterBondConstraintsAst([{}])",
            parts.join(", ")
        ))
    }

    /// Insert `c`, replacing any existing entry of the same key (last-wins).
    fn set(&mut self, py: Python<'_>, c: Py<MulticenterBondConstraintAst>) {
        self.0.set(c.bind(py).borrow().to_ast(py));
    }

    /// Remove the entry with the given key, returning it if present (dict `pop`).
    fn pop(
        &mut self,
        py: Python<'_>,
        key: Py<MulticenterBondConstraintKey>,
    ) -> PyResult<Option<MulticenterBondConstraintAst>> {
        self.0
            .remove(key.bind(py).borrow().to_ast())
            .map(|c| MulticenterBondConstraintAst::from_ast(py, &c))
            .transpose()
    }

    /// Overlay `other` onto self in place — another container or an iterable of
    /// `MulticenterBondConstraintAst` (last-wins per key; undetermined entries remove).
    fn update(&mut self, py: Python<'_>, other: MulticenterBondConstraintsUpdate) -> PyResult<()> {
        other.apply(py, &mut self.0)
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    fn __iter__(&self, py: Python<'_>) -> PyResult<MulticenterBondConstraintKeyIter> {
        multicenter_bond_constraint_keys(py, &self.0)
    }

    /// The constraint keys, in canonical order.
    fn keys(&self, py: Python<'_>) -> PyResult<MulticenterBondConstraintKeyIter> {
        multicenter_bond_constraint_keys(py, &self.0)
    }

    /// The constraints, in canonical order.
    fn values(&self, py: Python<'_>) -> PyResult<MulticenterBondConstraintIter> {
        multicenter_bond_constraints_iter(py, &self.0)
    }

    /// The `(key, constraint)` pairs, in canonical order.
    fn items(&self, py: Python<'_>) -> PyResult<MulticenterBondConstraintItemsIter> {
        multicenter_bond_constraint_items(py, &self.0)
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    fn get(
        &self,
        py: Python<'_>,
        key: Py<MulticenterBondConstraintKey>,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        match self.0.get(key.bind(py).borrow().to_ast()) {
            Some(constraint) => Ok(into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(py, constraint)?,
            )?
            .into_any()),
            None => Ok(default.unwrap_or_else(|| py.None())),
        }
    }

    /// The constraint with the given key; raises `KeyError` if absent.
    fn __getitem__(
        &self,
        py: Python<'_>,
        key: Py<MulticenterBondConstraintKey>,
    ) -> PyResult<MulticenterBondConstraintAst> {
        match self.0.get(key.bind(py).borrow().to_ast()) {
            Some(constraint) => MulticenterBondConstraintAst::from_ast(py, constraint),
            None => Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            )),
        }
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    fn __delitem__(
        &mut self,
        py: Python<'_>,
        key: Py<MulticenterBondConstraintKey>,
    ) -> PyResult<()> {
        if self.0.remove(key.bind(py).borrow().to_ast()).is_some() {
            Ok(())
        } else {
            Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            ))
        }
    }

    fn __contains__(&self, py: Python<'_>, key: Py<MulticenterBondConstraintKey>) -> bool {
        self.0.contains(key.bind(py).borrow().to_ast())
    }

    /// The asserted total electron count; `Undetermined` when no `ElectronCount`
    /// constraint is present (mirroring the non-optional Rust accessor).
    #[getter]
    fn electron_count(&self, py: Python<'_>) -> PyResult<ValueAst> {
        ValueAst::from_ast(py, &self.0.electron_count())
    }

    #[setter]
    fn set_electron_count(&mut self, py: Python<'_>, value: ValueArg) {
        self.0.set(AstMulticenterBondConstraintAst::electron_count(
            value.to_ast(py),
        ));
    }

    /// The present constraints as a dict keyed by snake_case name; values are the
    /// inner-value mirrors.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        multicenter_bond_constraints_asdict(py, &self.0)
    }
}

impl MulticenterBondConstraintsAst {
    /// The wrapped AST constraints — read access for multicenter bond construction.
    pub(crate) fn inner(&self) -> &AstMulticenterBondConstraintsAst {
        &self.0
    }

    /// Wrap AST constraints (the hold-the-value `from_inner` bridge). Test-only —
    /// in-crate construction wraps `MulticenterBondConstraintsAst(..)` directly.
    #[cfg(test)]
    pub(crate) fn from_inner(constraints: AstMulticenterBondConstraintsAst) -> Self {
        MulticenterBondConstraintsAst(constraints)
    }
}

/// Build the per-constraint iterator handle from a borrowed container.
fn multicenter_bond_constraints_iter(
    py: Python<'_>,
    constraints: &AstMulticenterBondConstraintsAst,
) -> PyResult<MulticenterBondConstraintIter> {
    let entries = constraints
        .iter()
        .map(|constraint| {
            into_py_variant(py, MulticenterBondConstraintAst::from_ast(py, constraint)?)
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(MulticenterBondConstraintIter {
        entries: entries.into_iter(),
    })
}

/// Build the key iterator handle from a borrowed container (mapping-style keys).
fn multicenter_bond_constraint_keys(
    py: Python<'_>,
    constraints: &AstMulticenterBondConstraintsAst,
) -> PyResult<MulticenterBondConstraintKeyIter> {
    let keys = constraints
        .iter()
        .map(|constraint| {
            into_py_variant(
                py,
                MulticenterBondConstraintKey::from_ast(&constraint.key()),
            )
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(MulticenterBondConstraintKeyIter {
        keys: keys.into_iter(),
    })
}

/// Build the item iterator handle (`(key, constraint)` pairs) from a borrowed container.
fn multicenter_bond_constraint_items(
    py: Python<'_>,
    constraints: &AstMulticenterBondConstraintsAst,
) -> PyResult<MulticenterBondConstraintItemsIter> {
    let items = constraints
        .iter()
        .map(|constraint| {
            Ok((
                into_py_variant(
                    py,
                    MulticenterBondConstraintKey::from_ast(&constraint.key()),
                )?,
                into_py_variant(py, MulticenterBondConstraintAst::from_ast(py, constraint)?)?,
            ))
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(MulticenterBondConstraintItemsIter {
        items: items.into_iter(),
    })
}

/// The present constraints as a dict keyed by snake_case name; values are the
/// inner-value mirrors.
fn multicenter_bond_constraints_asdict<'py>(
    py: Python<'py>,
    constraints: &AstMulticenterBondConstraintsAst,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    for entry in constraints.iter() {
        match entry {
            AstMulticenterBondConstraintAst::ElectronCount(v) => {
                dict.set_item("electron_count", ValueAst::from_ast(py, v)?)?
            }
        }
    }
    Ok(dict)
}

/// What a `MulticenterBondConstraintsView` writes through to: a multicenter bond
/// within a molecule (by index) or a standalone `MulticenterBondAst`.
enum MulticenterBondConstraintsBacking {
    Molecule {
        owner: Py<MoleculeAst>,
        id: AstMulticenterBondId,
    },
    MulticenterBond(Py<MulticenterBondAst>),
}

/// A live handle onto one multicenter bond's constraints, backed by either a
/// molecule-bond or a standalone `MulticenterBondAst`. Reads borrow the constraints
/// and read only the item they need (no whole-container clone); mutators write through
/// to the bond in place, without a clone-and-writeback.
#[pyclass]
pub struct MulticenterBondConstraintsView {
    backing: MulticenterBondConstraintsBacking,
}

impl MulticenterBondConstraintsView {
    /// Borrow the backing bond's constraints and read one item through `f` — no clone.
    fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&AstMulticenterBondConstraintsAst) -> PyResult<R>,
    ) -> PyResult<R> {
        match &self.backing {
            MulticenterBondConstraintsBacking::Molecule { owner, id } => {
                let molecule = owner.bind(py).borrow();
                let view = molecule
                    .inner()
                    .multicenter_bonds()
                    .get(*id)
                    .ok_or_else(|| PyIndexError::new_err("multicenter bond id out of range"))?;
                f(&view.ast.constraints)
            }
            MulticenterBondConstraintsBacking::MulticenterBond(bond) => {
                let bond = bond.bind(py).borrow();
                f(&bond.inner().constraints)
            }
        }
    }

    /// Mutate the backing bond's constraints in place through `f`.
    fn with_mut<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut AstMulticenterBondConstraintsAst) -> R,
    ) -> R {
        match &self.backing {
            MulticenterBondConstraintsBacking::Molecule { owner, id } => f(&mut owner
                .borrow_mut(py)
                .inner_mut()
                .multicenter_bond_mut(*id)
                .ast
                .constraints),
            MulticenterBondConstraintsBacking::MulticenterBond(bond) => {
                f(&mut bond.borrow_mut(py).inner_mut().constraints)
            }
        }
    }

    /// Set one constraint on the backing bond in place (last-wins per key).
    fn set_ast(&self, py: Python<'_>, constraint: AstMulticenterBondConstraintAst) {
        self.with_mut(py, |cs| cs.set(constraint));
    }

    /// Remove one key from the backing bond in place, returning the removed entry.
    fn remove_ast(
        &self,
        py: Python<'_>,
        key: AstMulticenterBondConstraintKey,
    ) -> Option<AstMulticenterBondConstraintAst> {
        self.with_mut(py, |cs| cs.remove(key))
    }
}

#[pymethods]
impl MulticenterBondConstraintsView {
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let count = self.read(py, |cs| Ok(cs.len()))?;
        Ok(format!("MulticenterBondConstraintsView({count} entries)"))
    }

    /// Insert `c` on the bond in place, replacing any existing entry of the same key
    /// (last-wins).
    fn set(&self, py: Python<'_>, c: Py<MulticenterBondConstraintAst>) {
        self.set_ast(py, c.bind(py).borrow().to_ast(py));
    }

    /// Remove the entry with the given key from the bond in place, returning it if
    /// present (dict `pop`).
    fn pop(
        &self,
        py: Python<'_>,
        key: Py<MulticenterBondConstraintKey>,
    ) -> PyResult<Option<MulticenterBondConstraintAst>> {
        self.remove_ast(py, key.bind(py).borrow().to_ast())
            .map(|c| MulticenterBondConstraintAst::from_ast(py, &c))
            .transpose()
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    fn __delitem__(&self, py: Python<'_>, key: Py<MulticenterBondConstraintKey>) -> PyResult<()> {
        if self
            .remove_ast(py, key.bind(py).borrow().to_ast())
            .is_some()
        {
            Ok(())
        } else {
            Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            ))
        }
    }

    /// Overlay `other` onto the bond's constraints in place — another container or an
    /// iterable of `MulticenterBondConstraintAst` (last-wins per key; undetermined entries remove).
    fn update(&self, py: Python<'_>, other: MulticenterBondConstraintsUpdate) -> PyResult<()> {
        self.with_mut(py, |cs| other.apply(py, cs))
    }

    fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        self.read(py, |cs| Ok(cs.len()))
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    fn __iter__(&self, py: Python<'_>) -> PyResult<MulticenterBondConstraintKeyIter> {
        self.read(py, |cs| multicenter_bond_constraint_keys(py, cs))
    }

    /// The constraint keys, in canonical order.
    fn keys(&self, py: Python<'_>) -> PyResult<MulticenterBondConstraintKeyIter> {
        self.read(py, |cs| multicenter_bond_constraint_keys(py, cs))
    }

    /// The constraints, in canonical order.
    fn values(&self, py: Python<'_>) -> PyResult<MulticenterBondConstraintIter> {
        self.read(py, |cs| multicenter_bond_constraints_iter(py, cs))
    }

    /// The `(key, constraint)` pairs, in canonical order.
    fn items(&self, py: Python<'_>) -> PyResult<MulticenterBondConstraintItemsIter> {
        self.read(py, |cs| multicenter_bond_constraint_items(py, cs))
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    fn get(
        &self,
        py: Python<'_>,
        key: Py<MulticenterBondConstraintKey>,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        let key = key.bind(py).borrow().to_ast();
        let found = self.read(py, |cs| {
            cs.get(key)
                .map(|constraint| MulticenterBondConstraintAst::from_ast(py, constraint))
                .transpose()
        })?;
        match found {
            Some(constraint) => Ok(into_py_variant(py, constraint)?.into_any()),
            None => Ok(default.unwrap_or_else(|| py.None())),
        }
    }

    /// The constraint with the given key; raises `KeyError` if absent.
    fn __getitem__(
        &self,
        py: Python<'_>,
        key: Py<MulticenterBondConstraintKey>,
    ) -> PyResult<MulticenterBondConstraintAst> {
        let ast_key = key.bind(py).borrow().to_ast();
        let found = self.read(py, |cs| {
            cs.get(ast_key)
                .map(|constraint| MulticenterBondConstraintAst::from_ast(py, constraint))
                .transpose()
        })?;
        match found {
            Some(constraint) => Ok(constraint),
            None => Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            )),
        }
    }

    fn __contains__(
        &self,
        py: Python<'_>,
        key: Py<MulticenterBondConstraintKey>,
    ) -> PyResult<bool> {
        let key = key.bind(py).borrow().to_ast();
        self.read(py, |cs| Ok(cs.contains(key)))
    }

    /// The asserted total electron count; `Undetermined` when no `ElectronCount`
    /// constraint is present (mirroring the non-optional Rust accessor).
    #[getter]
    fn electron_count(&self, py: Python<'_>) -> PyResult<ValueAst> {
        self.read(py, |cs| ValueAst::from_ast(py, &cs.electron_count()))
    }

    #[setter]
    fn set_electron_count(&self, py: Python<'_>, value: ValueArg) {
        self.set_ast(
            py,
            AstMulticenterBondConstraintAst::electron_count(value.to_ast(py)),
        );
    }

    /// The present constraints as a dict keyed by snake_case name.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        self.read(py, |cs| multicenter_bond_constraints_asdict(py, cs))
    }
}

#[pyclass]
struct MulticenterBondConstraintIter {
    entries: IntoIter<Py<MulticenterBondConstraintAst>>,
}

#[pymethods]
impl MulticenterBondConstraintIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<Py<MulticenterBondConstraintAst>> {
        self.entries.next()
    }
}

#[pyclass]
struct MulticenterBondConstraintKeyIter {
    keys: IntoIter<Py<MulticenterBondConstraintKey>>,
}

#[pymethods]
impl MulticenterBondConstraintKeyIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<Py<MulticenterBondConstraintKey>> {
        self.keys.next()
    }
}

#[pyclass]
struct MulticenterBondConstraintItemsIter {
    items: IntoIter<(
        Py<MulticenterBondConstraintKey>,
        Py<MulticenterBondConstraintAst>,
    )>,
}

#[pymethods]
impl MulticenterBondConstraintItemsIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(
        &mut self,
    ) -> Option<(
        Py<MulticenterBondConstraintKey>,
        Py<MulticenterBondConstraintAst>,
    )> {
        self.items.next()
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::{
        AtomAst as AstAtomAst, AtomId as AstAtomId, ElectronCountsAst as AstElectronCountsAst,
        MoleculeParts, SpinStateAst as AstSpinStateAst, ValueAst as AstValueAst,
    };
    use umol_chem::element::Element as ChemElement;

    use super::*;

    /// Three borons (atom ids 0–2) joined by one 3-center multicenter bond over all
    /// three (electrons `[1,1,1]`), multicenter bond id 0.
    fn three_center_bond(py: Python<'_>) -> Py<MoleculeAst> {
        let molecule = AstMoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AstAtomAst::from_element(ChemElement::B); 3],
            multicenter: vec![(
                (0u32..3).map(AstAtomId).collect(),
                AstMulticenterBondAst::from_electrons(vec![1, 1, 1]),
            )],
            ..Default::default()
        });
        Py::new(py, MoleculeAst::from_inner(molecule)).unwrap()
    }

    #[rstest]
    fn test_multicenter_bond_ast_new() {
        Python::attach(|py| {
            let spin_ast = AstSpinStateAst::from((0_u8, 1_u8));
            let spin = Py::new(py, SpinStateAst::from_ast(py, &spin_ast).unwrap()).unwrap();
            let bond = MulticenterBondAst::new(
                py,
                ElectronCountsArg::Lit(vec![1, 1, 1]),
                Some(ValueArg::Lit(-2)),
                Some(spin.bind(py).borrow()),
                None,
            );
            assert_eq!(
                bond.inner().electrons,
                AstElectronCountsAst::Lit(vec![1, 1, 1])
            );
            assert_eq!(bond.inner().charge, AstValueAst::Lit(-2));
            assert_eq!(bond.inner().spin, spin_ast);
        });
    }

    #[rstest]
    fn test_multicenter_bond_ast_new_constraints() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints =
                Py::new(py, MulticenterBondConstraintsAst::new(py, vec![ec])).unwrap();
            let bond = MulticenterBondAst::new(
                py,
                ElectronCountsArg::Lit(vec![1, 1, 1]),
                None,
                None,
                Some(constraints),
            );
            assert_eq!(
                bond.inner().constraints.electron_count(),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    #[case::undetermined("*")]
    #[case::electron_count("[1,1,1]#e6")]
    #[case::charge("[1,1,1]#c-2")]
    fn test_multicenter_bond_ast_parse(#[case] dsl: &str) {
        let bond = MulticenterBondAst::parse(dsl).unwrap();
        assert_eq!(bond.__str__(), dsl);
        assert_eq!(
            bond.__repr__(),
            format!("MulticenterBondAst.parse('{dsl}')")
        );
    }

    #[rstest]
    fn test_multicenter_bond_ast_parse_error() {
        assert!(MulticenterBondAst::parse("z").is_err());
    }

    #[rstest]
    fn test_multicenter_bond_ast_electrons() {
        Python::attach(|py| {
            let mut bond =
                MulticenterBondAst::from_inner(AstMulticenterBondAst::from_electrons(vec![
                    1, 1, 1,
                ]));
            assert_eq!(
                bond.electrons().to_ast(),
                AstElectronCountsAst::Lit(vec![1, 1, 1])
            );
            bond.set_electrons(py, ElectronCountsArg::Lit(vec![2, 2]));
            assert_eq!(
                bond.electrons().to_ast(),
                AstElectronCountsAst::Lit(vec![2, 2])
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_ast_charge() {
        Python::attach(|py| {
            let mut bond =
                MulticenterBondAst::from_inner(AstMulticenterBondAst::from_electrons(vec![
                    1, 1, 1,
                ]));
            bond.set_charge(py, ValueArg::Lit(-1));
            assert_eq!(bond.charge(py).unwrap().to_ast(py), AstValueAst::Lit(-1));
        });
    }

    #[rstest]
    fn test_multicenter_bond_ast_spin() {
        Python::attach(|py| {
            let spin_ast = AstSpinStateAst::from((0_u8, 1_u8));
            let spin = Py::new(py, SpinStateAst::from_ast(py, &spin_ast).unwrap()).unwrap();
            let mut bond =
                MulticenterBondAst::from_inner(AstMulticenterBondAst::from_electrons(vec![
                    1, 1, 1,
                ]));
            bond.set_spin(py, spin.bind(py).borrow());
            assert_eq!(bond.spin(py).unwrap().to_ast(py), spin_ast);
        });
    }

    #[rstest]
    fn test_multicenter_bond_ast_set_constraints_from_view() {
        Python::attach(|py| {
            let src = Py::new(
                py,
                MulticenterBondAst::from_inner(
                    AstMulticenterBondAst::from_electrons(vec![1, 1, 1])
                        .with_constraint(AstMulticenterBondConstraintAst::electron_count(6)),
                ),
            )
            .unwrap();
            let view = Py::new(
                py,
                MulticenterBondConstraintsView {
                    backing: MulticenterBondConstraintsBacking::MulticenterBond(src),
                },
            )
            .unwrap();
            let mut dst =
                MulticenterBondAst::from_inner(AstMulticenterBondAst::from_electrons(vec![
                    1, 1, 1,
                ]));
            dst.set_constraints(py, MulticenterBondConstraintsArg::View(view))
                .unwrap();
            assert_eq!(
                dst.inner().constraints.electron_count(),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_ast_asdict() {
        Python::attach(|py| {
            let bond = MulticenterBondAst::from_inner(
                AstMulticenterBondAst::from_electrons(vec![1, 1, 1])
                    .with_constraint(AstMulticenterBondConstraintAst::electron_count(6)),
            );
            let dict = bond.asdict(py).unwrap();
            assert_eq!(dict.len(), 4);
            let electrons = dict.get_item("electrons").unwrap().unwrap();
            let expected = into_py_variant(py, ElectronCountsAst::Lit(vec![1, 1, 1])).unwrap();
            assert!(electrons.eq(expected.bind(py)).unwrap());
            assert!(dict.contains("charge").unwrap());
            assert!(dict.contains("spin").unwrap());
            assert!(dict.contains("constraints").unwrap());
        });
    }

    #[rstest]
    fn test_multicenter_bond_view_atom_ids() {
        Python::attach(|py| {
            let view = MulticenterBondView {
                owner: three_center_bond(py),
                id: AstMulticenterBondId(0),
            };
            assert_eq!(view.id(), 0);
            let atom_ids: Vec<u32> = view.atom_ids(py).unwrap().extract().unwrap();
            assert_eq!(atom_ids, vec![0, 1, 2]);
            assert_eq!(view.__repr__(), "MulticenterBondView(id=0)");
        });
    }

    #[rstest]
    fn test_multicenter_bond_view_electrons() {
        Python::attach(|py| {
            let owner = three_center_bond(py);
            let view = MulticenterBondView {
                owner: owner.clone_ref(py),
                id: AstMulticenterBondId(0),
            };
            assert_eq!(
                view.electrons(py).unwrap().to_ast(),
                AstElectronCountsAst::Lit(vec![1, 1, 1])
            );
            view.set_electrons(py, ElectronCountsArg::Lit(vec![2, 2, 2]));
            let fresh = MulticenterBondView {
                owner,
                id: AstMulticenterBondId(0),
            };
            assert_eq!(
                fresh.electrons(py).unwrap().to_ast(),
                AstElectronCountsAst::Lit(vec![2, 2, 2])
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_view_charge() {
        Python::attach(|py| {
            let owner = three_center_bond(py);
            let view = MulticenterBondView {
                owner: owner.clone_ref(py),
                id: AstMulticenterBondId(0),
            };
            view.set_charge(py, ValueArg::Lit(-1));
            let fresh = MulticenterBondView {
                owner,
                id: AstMulticenterBondId(0),
            };
            assert_eq!(fresh.charge(py).unwrap().to_ast(py), AstValueAst::Lit(-1));
        });
    }

    #[rstest]
    fn test_multicenter_bond_view_spin() {
        Python::attach(|py| {
            let spin_ast = AstSpinStateAst::from((0_u8, 1_u8));
            let spin = Py::new(py, SpinStateAst::from_ast(py, &spin_ast).unwrap()).unwrap();
            let owner = three_center_bond(py);
            let view = MulticenterBondView {
                owner: owner.clone_ref(py),
                id: AstMulticenterBondId(0),
            };
            view.set_spin(py, spin.bind(py).borrow());
            let fresh = MulticenterBondView {
                owner,
                id: AstMulticenterBondId(0),
            };
            assert_eq!(fresh.spin(py).unwrap().to_ast(py), spin_ast);
        });
    }

    #[rstest]
    fn test_multicenter_bond_view_constraints() {
        Python::attach(|py| {
            let view = MulticenterBondView {
                owner: three_center_bond(py),
                id: AstMulticenterBondId(0),
            };
            match view.constraints(py).backing {
                MulticenterBondConstraintsBacking::Molecule { id, .. } => {
                    assert_eq!(id, AstMulticenterBondId(0))
                }
                _ => panic!("expected molecule-backed view"),
            }
        });
    }

    #[rstest]
    fn test_multicenter_bond_view_set_constraints() {
        Python::attach(|py| {
            let owner = three_center_bond(py);
            let view = MulticenterBondView {
                owner: owner.clone_ref(py),
                id: AstMulticenterBondId(0),
            };
            let constraints = Py::new(
                py,
                MulticenterBondConstraintsAst::new(
                    py,
                    vec![into_py_variant(
                        py,
                        MulticenterBondConstraintAst::from_ast(
                            py,
                            &AstMulticenterBondConstraintAst::electron_count(6),
                        )
                        .unwrap(),
                    )
                    .unwrap()],
                ),
            )
            .unwrap();
            view.set_constraints(py, MulticenterBondConstraintsArg::Container(constraints))
                .unwrap();
            let fresh = MulticenterBondView {
                owner,
                id: AstMulticenterBondId(0),
            };
            assert_eq!(
                fresh.constraints(py).electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_view_asdict() {
        Python::attach(|py| {
            let view = MulticenterBondView {
                owner: three_center_bond(py),
                id: AstMulticenterBondId(0),
            };
            let dict = view.asdict(py).unwrap();
            assert_eq!(dict.len(), 4);
            let electrons = dict.get_item("electrons").unwrap().unwrap();
            let expected = into_py_variant(py, ElectronCountsAst::Lit(vec![1, 1, 1])).unwrap();
            assert!(electrons.eq(expected.bind(py)).unwrap());
            assert!(dict.contains("charge").unwrap());
            assert!(dict.contains("spin").unwrap());
            assert!(dict.contains("constraints").unwrap());
        });
    }

    #[rstest]
    fn test_multicenter_bond_views_len_and_getitem() {
        Python::attach(|py| {
            let views = MulticenterBondViews {
                owner: three_center_bond(py),
            };
            assert_eq!(views.__len__(py), 1);
            assert_eq!(views.__getitem__(py, 0).unwrap().id(), 0);
            assert_eq!(views.__getitem__(py, -1).unwrap().id(), 0);
            assert!(views.__getitem__(py, 5).is_err());
            assert!(views.__getitem__(py, -2).is_err());
        });
    }

    #[rstest]
    fn test_multicenter_bond_views_repr() {
        Python::attach(|py| {
            let views = MulticenterBondViews {
                owner: three_center_bond(py),
            };
            assert_eq!(views.__repr__(py), "MulticenterBondViews(len=1)");
        });
    }

    #[rstest]
    fn test_multicenter_bond_views_setitem() {
        Python::attach(|py| {
            let owner = three_center_bond(py);
            let views = MulticenterBondViews {
                owner: owner.clone_ref(py),
            };
            let replacement = Py::new(
                py,
                MulticenterBondAst::from_inner(AstMulticenterBondAst::from_electrons(vec![
                    2, 2, 2,
                ])),
            )
            .unwrap();
            views
                .__setitem__(py, 0, replacement.bind(py).borrow())
                .unwrap();
            let view = views.__getitem__(py, 0).unwrap();
            // value replaced, members preserved
            assert_eq!(
                view.electrons(py).unwrap().to_ast(),
                AstElectronCountsAst::Lit(vec![2, 2, 2])
            );
            let atom_ids: Vec<u32> = view.atom_ids(py).unwrap().extract().unwrap();
            assert_eq!(atom_ids, vec![0, 1, 2]);
        });
    }

    #[rstest]
    fn test_multicenter_bond_views_setitem_error() {
        Python::attach(|py| {
            let views = MulticenterBondViews {
                owner: three_center_bond(py),
            };
            let replacement = Py::new(
                py,
                MulticenterBondAst::from_inner(AstMulticenterBondAst::from_electrons(vec![
                    1, 1, 1,
                ])),
            )
            .unwrap();
            assert!(views
                .__setitem__(py, 5, replacement.bind(py).borrow())
                .is_err());
        });
    }

    #[rstest]
    fn test_multicenter_bond_views_connecting() {
        Python::attach(|py| {
            let views = MulticenterBondViews {
                owner: three_center_bond(py),
            };
            assert_eq!(views.connecting(py, vec![0, 1, 2]).unwrap().id(), 0);
            // a subset is not the bond's exact atom set
            assert!(views.connecting(py, vec![0, 1]).is_none());
        });
    }

    #[rstest]
    fn test_multicenter_bond_views_incident() {
        Python::attach(|py| {
            // three borons bonded plus one isolated boron (atom id 3)
            let molecule = AstMoleculeAst::from_parts(MoleculeParts {
                atoms: vec![AstAtomAst::from_element(ChemElement::B); 4],
                multicenter: vec![(
                    (0u32..3).map(AstAtomId).collect(),
                    AstMulticenterBondAst::from_electrons(vec![1, 1, 1]),
                )],
                ..Default::default()
            });
            let views = MulticenterBondViews {
                owner: Py::new(py, MoleculeAst::from_inner(molecule)).unwrap(),
            };
            assert_eq!(
                views
                    .incident(py, 0)
                    .iter()
                    .map(|v| v.id())
                    .collect::<Vec<_>>(),
                vec![0]
            );
            assert!(views.incident(py, 3).is_empty());
        });
    }

    #[rstest]
    fn test_multicenter_bond_views_iter() {
        Python::attach(|py| {
            let views = MulticenterBondViews {
                owner: three_center_bond(py),
            };
            let mut iter = views.__iter__(py);
            assert_eq!(iter.__next__(py).unwrap().id(), 0);
            assert!(iter.__next__(py).is_none());
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraint_key_roundtrip() {
        let key =
            MulticenterBondConstraintKey::from_ast(&AstMulticenterBondConstraintKey::ElectronCount);
        assert_eq!(key.to_ast(), AstMulticenterBondConstraintKey::ElectronCount);
    }

    #[rstest]
    fn test_multicenter_bond_constraint_ast_key() {
        Python::attach(|py| {
            let constraint = AstMulticenterBondConstraintAst::electron_count(6);
            let key = MulticenterBondConstraintAst::from_ast(py, &constraint)
                .unwrap()
                .key(py);
            assert_eq!(key.to_ast(), AstMulticenterBondConstraintKey::ElectronCount);
        });
    }

    #[rstest]
    #[case(AstMulticenterBondConstraintAst::electron_count(6))]
    #[case(AstMulticenterBondConstraintAst::electron_count(AstValueAst::Undetermined))]
    fn test_multicenter_bond_constraint_ast_roundtrip(
        #[case] ast: AstMulticenterBondConstraintAst,
    ) {
        Python::attach(|py| {
            assert_eq!(
                MulticenterBondConstraintAst::from_ast(py, &ast)
                    .unwrap()
                    .to_ast(py),
                ast
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_new() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = MulticenterBondConstraintsAst::new(py, vec![ec]);
            assert_eq!(constraints.__len__(), 1);
            assert_eq!(
                constraints.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_repr() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = MulticenterBondConstraintsAst::new(py, vec![ec]);
            assert_eq!(
                constraints.__repr__(py).unwrap(),
                "MulticenterBondConstraintsAst([MulticenterBondConstraintAst.ElectronCount(ValueAst.Lit(6))])"
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_set() {
        Python::attach(|py| {
            let mut constraints = MulticenterBondConstraintsAst::new(py, vec![]);
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            constraints.set(py, ec);
            assert_eq!(constraints.__len__(), 1);
            assert_eq!(
                constraints.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_pop() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let mut constraints = MulticenterBondConstraintsAst::new(py, vec![ec]);
            let key = into_py_variant(py, MulticenterBondConstraintKey::ElectronCount()).unwrap();
            let removed = constraints.pop(py, key).unwrap();
            match removed {
                Some(MulticenterBondConstraintAst::ElectronCount(v)) => {
                    assert_eq!(v.bind(py).borrow().to_ast(py), AstValueAst::Lit(6))
                }
                _ => panic!("expected removed ElectronCount(Lit(6))"),
            }
            assert_eq!(constraints.__len__(), 0);
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_update() {
        Python::attach(|py| {
            let mut constraints = MulticenterBondConstraintsAst::new(py, vec![]);
            let mut other = AstMulticenterBondConstraintsAst::new();
            other.set(AstMulticenterBondConstraintAst::electron_count(6));
            constraints
                .update(
                    py,
                    MulticenterBondConstraintsUpdate::Container(
                        Py::new(py, MulticenterBondConstraintsAst::from_inner(other)).unwrap(),
                    ),
                )
                .unwrap();
            assert_eq!(constraints.__len__(), 1);
            assert_eq!(
                constraints.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_update_entries() {
        Python::attach(|py| {
            let mut constraints = MulticenterBondConstraintsAst::new(py, vec![]);
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            constraints
                .update(py, MulticenterBondConstraintsUpdate::Entries(vec![ec]))
                .unwrap();
            assert_eq!(constraints.__len__(), 1);
            assert_eq!(
                constraints.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_len_contains() {
        Python::attach(|py| {
            let empty = MulticenterBondConstraintsAst::new(py, vec![]);
            assert_eq!(empty.__len__(), 0);
            assert!(!empty.__contains__(
                py,
                into_py_variant(py, MulticenterBondConstraintKey::ElectronCount()).unwrap()
            ));
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = MulticenterBondConstraintsAst::new(py, vec![ec]);
            assert!(constraints.__contains__(
                py,
                into_py_variant(py, MulticenterBondConstraintKey::ElectronCount()).unwrap()
            ));
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_keys_values_items() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = MulticenterBondConstraintsAst::new(py, vec![ec]);

            let mut keys = constraints.__iter__(py).unwrap();
            assert_eq!(
                keys.__next__().unwrap().bind(py).borrow().to_ast(),
                AstMulticenterBondConstraintKey::ElectronCount
            );
            assert!(keys.__next__().is_none());

            let mut values = constraints.values(py).unwrap();
            assert_eq!(
                values.__next__().unwrap().bind(py).borrow().to_ast(py),
                AstMulticenterBondConstraintAst::electron_count(6)
            );

            let mut items = constraints.items(py).unwrap();
            let (key, value) = items.__next__().unwrap();
            assert_eq!(
                key.bind(py).borrow().to_ast(),
                AstMulticenterBondConstraintKey::ElectronCount
            );
            assert_eq!(
                value.bind(py).borrow().to_ast(py),
                AstMulticenterBondConstraintAst::electron_count(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_get() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = MulticenterBondConstraintsAst::new(py, vec![ec]);
            let present = constraints
                .get(
                    py,
                    into_py_variant(py, MulticenterBondConstraintKey::ElectronCount()).unwrap(),
                    None,
                )
                .unwrap();
            let expected = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap()
            .into_any();
            assert!(present.bind(py).eq(expected.bind(py)).unwrap());

            let empty = MulticenterBondConstraintsAst::new(py, vec![]);
            let absent = empty
                .get(
                    py,
                    into_py_variant(py, MulticenterBondConstraintKey::ElectronCount()).unwrap(),
                    None,
                )
                .unwrap();
            assert!(absent.bind(py).is_none());

            // a caller-supplied default is returned verbatim when the key is absent
            let sentinel = into_py_variant(py, MulticenterBondConstraintKey::ElectronCount())
                .unwrap()
                .into_any();
            let defaulted = empty
                .get(
                    py,
                    into_py_variant(py, MulticenterBondConstraintKey::ElectronCount()).unwrap(),
                    Some(sentinel.clone_ref(py)),
                )
                .unwrap();
            assert_eq!(defaulted.as_ptr(), sentinel.as_ptr());
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_electron_count() {
        Python::attach(|py| {
            let empty = MulticenterBondConstraintsAst::new(py, vec![]);
            assert_eq!(
                empty.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Undetermined
            );
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = MulticenterBondConstraintsAst::new(py, vec![ec]);
            assert_eq!(
                constraints.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_set_electron_count() {
        Python::attach(|py| {
            let mut constraints = MulticenterBondConstraintsAst::new(py, vec![]);
            constraints.set_electron_count(py, ValueArg::Lit(6));
            assert_eq!(
                constraints.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_getitem_error() {
        Python::attach(|py| {
            let constraints = MulticenterBondConstraintsAst::new(py, vec![]);
            let key = into_py_variant(py, MulticenterBondConstraintKey::ElectronCount()).unwrap();
            assert!(constraints.__getitem__(py, key).is_err());
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_delitem_error() {
        Python::attach(|py| {
            let mut constraints = MulticenterBondConstraintsAst::new(py, vec![]);
            let key = into_py_variant(py, MulticenterBondConstraintKey::ElectronCount()).unwrap();
            assert!(constraints.__delitem__(py, key).is_err());
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_asdict() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = MulticenterBondConstraintsAst::new(py, vec![ec]);
            let dict = constraints.asdict(py).unwrap();
            assert_eq!(dict.len(), 1);
            let value = dict.get_item("electron_count").unwrap().unwrap();
            let expected =
                into_py_variant(py, ValueAst::from_ast(py, &AstValueAst::Lit(6)).unwrap()).unwrap();
            assert!(value.eq(expected.bind(py)).unwrap());
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_view_set() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                MulticenterBondAst::from_inner(AstMulticenterBondAst::from_electrons(vec![
                    1, 1, 1,
                ])),
            )
            .unwrap();
            let view = MulticenterBondConstraintsView {
                backing: MulticenterBondConstraintsBacking::MulticenterBond(bond.clone_ref(py)),
            };
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            view.set(py, ec);
            // a fresh view proves the write hit the standalone bond, not a copy
            let fresh = MulticenterBondConstraintsView {
                backing: MulticenterBondConstraintsBacking::MulticenterBond(bond),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 1);
            assert_eq!(
                fresh.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_view_pop() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                MulticenterBondAst::from_inner(
                    AstMulticenterBondAst::from_electrons(vec![1, 1, 1])
                        .with_constraint(AstMulticenterBondConstraintAst::electron_count(6)),
                ),
            )
            .unwrap();
            let view = MulticenterBondConstraintsView {
                backing: MulticenterBondConstraintsBacking::MulticenterBond(bond.clone_ref(py)),
            };
            let removed = view
                .pop(
                    py,
                    into_py_variant(py, MulticenterBondConstraintKey::ElectronCount()).unwrap(),
                )
                .unwrap();
            match removed {
                Some(MulticenterBondConstraintAst::ElectronCount(v)) => {
                    assert_eq!(v.bind(py).borrow().to_ast(py), AstValueAst::Lit(6))
                }
                _ => panic!("expected removed ElectronCount(Lit(6))"),
            }
            let fresh = MulticenterBondConstraintsView {
                backing: MulticenterBondConstraintsBacking::MulticenterBond(bond),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 0);
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_view_update() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                MulticenterBondAst::from_inner(AstMulticenterBondAst::from_electrons(vec![
                    1, 1, 1,
                ])),
            )
            .unwrap();
            let view = MulticenterBondConstraintsView {
                backing: MulticenterBondConstraintsBacking::MulticenterBond(bond.clone_ref(py)),
            };
            let mut other = AstMulticenterBondConstraintsAst::new();
            other.set(AstMulticenterBondConstraintAst::electron_count(6));
            view.update(
                py,
                MulticenterBondConstraintsUpdate::Container(
                    Py::new(py, MulticenterBondConstraintsAst::from_inner(other)).unwrap(),
                ),
            )
            .unwrap();
            let fresh = MulticenterBondConstraintsView {
                backing: MulticenterBondConstraintsBacking::MulticenterBond(bond),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 1);
            assert_eq!(
                fresh.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_view_set_electron_count() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                MulticenterBondAst::from_inner(AstMulticenterBondAst::from_electrons(vec![
                    1, 1, 1,
                ])),
            )
            .unwrap();
            let view = MulticenterBondConstraintsView {
                backing: MulticenterBondConstraintsBacking::MulticenterBond(bond.clone_ref(py)),
            };
            assert_eq!(
                view.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Undetermined
            );
            view.set_electron_count(py, ValueArg::Lit(6));
            let fresh = MulticenterBondConstraintsView {
                backing: MulticenterBondConstraintsBacking::MulticenterBond(bond),
            };
            assert_eq!(
                fresh.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_view_set_molecule_backed() {
        Python::attach(|py| {
            let owner = three_center_bond(py);
            let view = MulticenterBondConstraintsView {
                backing: MulticenterBondConstraintsBacking::Molecule {
                    owner: owner.clone_ref(py),
                    id: AstMulticenterBondId(0),
                },
            };
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_ast(
                    py,
                    &AstMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            view.set(py, ec);
            let fresh = MulticenterBondConstraintsView {
                backing: MulticenterBondConstraintsBacking::Molecule {
                    owner,
                    id: AstMulticenterBondId(0),
                },
            };
            assert_eq!(fresh.__len__(py).unwrap(), 1);
            assert_eq!(
                fresh.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }
}
