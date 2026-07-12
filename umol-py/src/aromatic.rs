//! Aromatic system value type and aromatic-constraint surface mirroring
//! `umol_ast::ast`: `AromaticSystemAst`, the `AromaticSystemConstraintAst` enum, the
//! `AromaticSystemConstraintsAst` container, and the `AromaticSystemConstraintsView`
//! live handle. An aromatic system is a single unordered set of member atoms; the
//! value carries a positional per-atom `electrons` vector plus charge, spin, and
//! constraints. The member atoms are the participants of the owning molecule's
//! aromatic relation, so they are topology (the view half) rather than value.

use std::str::FromStr;
use std::vec::IntoIter;

use pyo3::exceptions::{PyIndexError, PyKeyError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict, PyTuple};
use umol_ast::ast::{
    AromaticSystemAst as AstAromaticSystemAst,
    AromaticSystemConstraintAst as AstAromaticSystemConstraintAst,
    AromaticSystemConstraintKey as AstAromaticSystemConstraintKey,
    AromaticSystemConstraintsAst as AstAromaticSystemConstraintsAst,
    AromaticSystemId as AstAromaticSystemId, AromaticSystemView as AstAromaticSystemView,
    MoleculeAst as AstMoleculeAst,
};

use crate::atom::SpinStateAst;
use crate::convert::{hash_ast, into_py_variant, variant_repr};
use crate::electrons::{ElectronCountsArg, ElectronCountsAst};
use crate::error::parse_error;
use crate::molecule::MoleculeAst;
use crate::value::{ValueArg, ValueAst};

/// An aromatic system: a positional per-member-atom `electrons` vector, charge,
/// spin, and aromatic-system-scope constraints. The member atoms are the
/// participants of the owning molecule's aromatic relation (the view half); the
/// `electrons` vector is positional, aligned to that atom order.
#[pyclass(eq)]
#[derive(PartialEq)]
pub struct AromaticSystemAst(AstAromaticSystemAst);

#[pymethods]
impl AromaticSystemAst {
    /// Construct from an electron-count vector — a `list[int]` or an
    /// `ElectronCountsAst` — optionally setting fields.
    #[new]
    #[pyo3(signature = (electrons, *, charge=None, spin=None, constraints=None))]
    fn new(
        py: Python<'_>,
        electrons: ElectronCountsArg,
        charge: Option<ValueArg>,
        spin: Option<PyRef<'_, SpinStateAst>>,
        constraints: Option<Py<AromaticSystemConstraintsAst>>,
    ) -> Self {
        let mut system = AstAromaticSystemAst::new(electrons.to_ast(py));
        if let Some(charge) = charge {
            system = system.with_charge(charge.to_ast(py));
        }
        if let Some(spin) = spin {
            system = system.with_spin(spin.to_ast(py));
        }
        if let Some(constraints) = constraints {
            system.constraints = constraints.bind(py).borrow().inner().clone();
        }
        AromaticSystemAst(system)
    }

    /// Parse an aromatic-system-DSL string (e.g. `"[1,1,1]#e6"`) into an `AromaticSystemAst`.
    #[staticmethod]
    fn parse(s: &str) -> PyResult<Self> {
        AstAromaticSystemAst::from_str(s)
            .map(Self)
            .map_err(parse_error)
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("AromaticSystemAst.parse('{}')", self.0)
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

    /// The system's constraints as a live handle onto this system: reads borrow the
    /// current state, mutators write through to the system in place.
    #[getter]
    fn constraints(slf: Py<Self>) -> AromaticSystemConstraintsView {
        AromaticSystemConstraintsView {
            backing: AromaticSystemConstraintsBacking::AromaticSystem(slf),
        }
    }

    /// Replace the whole constraint set (wipe-and-set) from a value container or a
    /// live view.
    #[setter]
    fn set_constraints(
        &mut self,
        py: Python<'_>,
        value: AromaticSystemConstraintsArg,
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
            aromatic_system_constraints_asdict(py, &self.0.constraints)?,
        )?;
        Ok(dict)
    }
}

impl AromaticSystemAst {
    /// The wrapped AST system — read access for the system-backed constraints view.
    pub(crate) fn inner(&self) -> &AstAromaticSystemAst {
        &self.0
    }

    /// Mutable access to the wrapped AST system — write access for the system-backed
    /// constraints view.
    pub(crate) fn inner_mut(&mut self) -> &mut AstAromaticSystemAst {
        &mut self.0
    }

    /// Wrap an AST system (the hold-the-value `from_inner` bridge, paired with
    /// `inner`). Test-only — in-crate construction wraps `AromaticSystemAst(..)` directly.
    #[cfg(test)]
    pub(crate) fn from_inner(system: AstAromaticSystemAst) -> Self {
        AromaticSystemAst(system)
    }
}

/// A view of one aromatic system within a molecule: a handle to the molecule plus
/// the system's index. Field reads rebuild the transient Rust view; the molecule is
/// never copied. The member atom indices are read-only topology; the electrons,
/// charge, spin, and constraints are the mutable system value.
#[pyclass]
pub struct AromaticSystemView {
    owner: Py<MoleculeAst>,
    id: AstAromaticSystemId,
}

impl AromaticSystemView {
    fn aromatic_system<'a>(
        &self,
        molecule: &'a AstMoleculeAst,
    ) -> PyResult<AstAromaticSystemView<'a>> {
        molecule
            .aromatic_systems()
            .get(self.id)
            .ok_or_else(|| PyIndexError::new_err("aromatic system id out of range"))
    }
}

#[pymethods]
impl AromaticSystemView {
    #[getter]
    fn id(&self) -> u32 {
        self.id.0
    }

    /// The member atom indices (read-only — participants are topology, not part of
    /// the system value).
    #[getter]
    fn atom_ids<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let molecule = self.owner.bind(py).borrow();
        let atom_ids: Vec<u32> = self
            .aromatic_system(molecule.inner())?
            .atom_ids()
            .map(|atom| atom.0)
            .collect();
        PyTuple::new(py, atom_ids)
    }

    fn __repr__(&self) -> String {
        format!("AromaticSystemView(id={})", self.id.0)
    }

    /// The per-member-atom electron counts (positional, aligned to `atom_ids`).
    #[getter]
    fn electrons(&self, py: Python<'_>) -> PyResult<ElectronCountsAst> {
        let molecule = self.owner.bind(py).borrow();
        Ok(ElectronCountsAst::from_ast(
            &self.aromatic_system(molecule.inner())?.ast.electrons,
        ))
    }

    #[setter]
    fn set_electrons(&self, py: Python<'_>, value: ElectronCountsArg) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .aromatic_system_mut(self.id)
            .ast
            .electrons = value.to_ast(py);
    }

    #[getter]
    fn charge(&self, py: Python<'_>) -> PyResult<ValueAst> {
        let molecule = self.owner.bind(py).borrow();
        ValueAst::from_ast(py, &self.aromatic_system(molecule.inner())?.ast.charge)
    }

    #[setter]
    fn set_charge(&self, py: Python<'_>, value: ValueArg) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .aromatic_system_mut(self.id)
            .ast
            .charge = value.to_ast(py);
    }

    #[getter]
    fn spin(&self, py: Python<'_>) -> PyResult<SpinStateAst> {
        let molecule = self.owner.bind(py).borrow();
        SpinStateAst::from_ast(py, &self.aromatic_system(molecule.inner())?.ast.spin)
    }

    #[setter]
    fn set_spin(&self, py: Python<'_>, value: PyRef<'_, SpinStateAst>) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .aromatic_system_mut(self.id)
            .ast
            .spin = value.to_ast(py);
    }

    /// The system's constraints as a live handle onto the molecule: reads borrow the
    /// current state, mutators write through to the system in place.
    #[getter]
    fn constraints(&self, py: Python<'_>) -> AromaticSystemConstraintsView {
        AromaticSystemConstraintsView {
            backing: AromaticSystemConstraintsBacking::Molecule {
                owner: self.owner.clone_ref(py),
                id: self.id,
            },
        }
    }

    /// Replace the whole constraint set of the backing system in place (wipe-and-set)
    /// from a value container or a live view.
    #[setter]
    fn set_constraints(&self, py: Python<'_>, value: AromaticSystemConstraintsArg) -> PyResult<()> {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .aromatic_system_mut(self.id)
            .ast
            .constraints = value.to_ast(py)?;
        Ok(())
    }

    /// The value fields as a dict keyed by field name; values are the field mirrors —
    /// symmetric with `AromaticSystemAst.asdict`, read through the view.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let molecule = self.owner.bind(py).borrow();
        let system = self.aromatic_system(molecule.inner())?.ast;
        let dict = PyDict::new(py);
        dict.set_item("electrons", ElectronCountsAst::from_ast(&system.electrons))?;
        dict.set_item("charge", ValueAst::from_ast(py, &system.charge)?)?;
        dict.set_item("spin", SpinStateAst::from_ast(py, &system.spin)?)?;
        dict.set_item(
            "constraints",
            aromatic_system_constraints_asdict(py, &system.constraints)?,
        )?;
        Ok(dict)
    }
}

/// The key (identity) of an aromatic-system constraint, for keyed lookup. The
/// single key `ElectronCount` is the bare discriminant (no sub-key).
#[pyclass]
pub enum AromaticSystemConstraintKey {
    ElectronCount(),
}

#[pymethods]
impl AromaticSystemConstraintKey {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_ast() == other.to_ast()
    }

    fn __hash__(&self) -> u64 {
        hash_ast(&self.to_ast())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            AromaticSystemConstraintKey::ElectronCount() => ("ElectronCount", 0),
        };
        variant_repr(
            slf.bind(py).as_any(),
            "AromaticSystemConstraintKey",
            variant,
            arity,
        )
    }
}

impl AromaticSystemConstraintKey {
    pub(crate) fn from_ast(ast: &AstAromaticSystemConstraintKey) -> Self {
        match ast {
            AstAromaticSystemConstraintKey::ElectronCount => Self::ElectronCount(),
        }
    }

    pub(crate) fn to_ast(&self) -> AstAromaticSystemConstraintKey {
        match self {
            Self::ElectronCount() => AstAromaticSystemConstraintKey::ElectronCount,
        }
    }
}

/// An aromatic-system-scope constraint: the asserted total π-electron count of the
/// system (cross-checked against `sum(AromaticSystemAst::electrons)`).
#[pyclass]
pub enum AromaticSystemConstraintAst {
    ElectronCount(Py<ValueAst>),
}

#[pymethods]
impl AromaticSystemConstraintAst {
    /// The constraint's key (identity).
    #[getter]
    fn key(&self, py: Python<'_>) -> AromaticSystemConstraintKey {
        AromaticSystemConstraintKey::from_ast(&self.to_ast(py).key())
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let variant = match &*slf.bind(py).borrow() {
            AromaticSystemConstraintAst::ElectronCount(_) => "ElectronCount",
        };
        variant_repr(
            slf.bind(py).as_any(),
            "AromaticSystemConstraintAst",
            variant,
            1,
        )
    }
}

impl AromaticSystemConstraintAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstAromaticSystemConstraintAst) -> PyResult<Self> {
        Ok(match ast {
            AstAromaticSystemConstraintAst::ElectronCount(v) => {
                Self::ElectronCount(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstAromaticSystemConstraintAst {
        match self {
            Self::ElectronCount(v) => {
                AstAromaticSystemConstraintAst::ElectronCount(v.bind(py).borrow().to_ast(py))
            }
        }
    }
}

/// The argument to `update`: another constraint container (value or live view) or an
/// iterable of `AromaticSystemConstraintAst` (each `set`, last-wins).
#[derive(FromPyObject)]
enum AromaticSystemConstraintsUpdate {
    Container(Py<AromaticSystemConstraintsAst>),
    View(Py<AromaticSystemConstraintsView>),
    Entries(Vec<Py<AromaticSystemConstraintAst>>),
}

impl AromaticSystemConstraintsUpdate {
    /// Overlay this update onto `target` in place.
    fn apply(&self, py: Python<'_>, target: &mut AstAromaticSystemConstraintsAst) -> PyResult<()> {
        match self {
            AromaticSystemConstraintsUpdate::Container(c) => {
                target.update(c.bind(py).borrow().inner())
            }
            AromaticSystemConstraintsUpdate::View(v) => {
                let snapshot = v.bind(py).borrow().read(py, |cs| Ok(cs.clone()))?;
                target.update(&snapshot);
            }
            AromaticSystemConstraintsUpdate::Entries(entries) => {
                for entry in entries {
                    target.set(entry.bind(py).borrow().to_ast(py));
                }
            }
        }
        Ok(())
    }
}

/// A whole-container argument that snapshots either a value container or a live view
/// — for the aromatic system `constraints` setter, which accepts either.
#[derive(FromPyObject)]
enum AromaticSystemConstraintsArg {
    Container(Py<AromaticSystemConstraintsAst>),
    View(Py<AromaticSystemConstraintsView>),
}

impl AromaticSystemConstraintsArg {
    fn to_ast(&self, py: Python<'_>) -> PyResult<AstAromaticSystemConstraintsAst> {
        match self {
            AromaticSystemConstraintsArg::Container(c) => Ok(c.bind(py).borrow().inner().clone()),
            AromaticSystemConstraintsArg::View(v) => {
                v.bind(py).borrow().read(py, |cs| Ok(cs.clone()))
            }
        }
    }
}

/// The aromatic-system-scope constraints on an aromatic system, in kind-sorted order.
/// Mutable, hence value-equal but unhashable (matching `AromaticSystemAst`).
#[pyclass(eq)]
#[derive(PartialEq)]
pub struct AromaticSystemConstraintsAst(AstAromaticSystemConstraintsAst);

#[pymethods]
impl AromaticSystemConstraintsAst {
    /// Build from a sequence of constraints (a later entry of the same key replaces
    /// an earlier one, last-wins).
    #[new]
    fn new(py: Python<'_>, entries: Vec<Py<AromaticSystemConstraintAst>>) -> Self {
        let mut constraints = AstAromaticSystemConstraintsAst::new();
        constraints.extend(
            entries
                .into_iter()
                .map(|entry| entry.bind(py).borrow().to_ast(py)),
        );
        AromaticSystemConstraintsAst(constraints)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let mut parts = Vec::with_capacity(self.0.len());
        for entry in self.0.iter() {
            let mirror = into_py_variant(py, AromaticSystemConstraintAst::from_ast(py, entry)?)?;
            parts.push(mirror.bind(py).as_any().repr()?.extract::<String>()?);
        }
        Ok(format!(
            "AromaticSystemConstraintsAst([{}])",
            parts.join(", ")
        ))
    }

    /// Insert `c`, replacing any existing entry of the same key (last-wins).
    fn set(&mut self, py: Python<'_>, c: Py<AromaticSystemConstraintAst>) {
        self.0.set(c.bind(py).borrow().to_ast(py));
    }

    /// Remove the entry with the given key, returning it if present (dict `pop`).
    fn pop(
        &mut self,
        py: Python<'_>,
        key: Py<AromaticSystemConstraintKey>,
    ) -> PyResult<Option<AromaticSystemConstraintAst>> {
        self.0
            .remove(key.bind(py).borrow().to_ast())
            .map(|c| AromaticSystemConstraintAst::from_ast(py, &c))
            .transpose()
    }

    /// Overlay `other` onto self in place — another container or an iterable of
    /// `AromaticSystemConstraintAst` (last-wins per key; undetermined entries remove).
    fn update(&mut self, py: Python<'_>, other: AromaticSystemConstraintsUpdate) -> PyResult<()> {
        other.apply(py, &mut self.0)
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    fn __iter__(&self, py: Python<'_>) -> PyResult<AromaticSystemConstraintKeyIter> {
        aromatic_system_constraint_keys(py, &self.0)
    }

    /// The constraint keys, in canonical order.
    fn keys(&self, py: Python<'_>) -> PyResult<AromaticSystemConstraintKeyIter> {
        aromatic_system_constraint_keys(py, &self.0)
    }

    /// The constraints, in canonical order.
    fn values(&self, py: Python<'_>) -> PyResult<AromaticSystemConstraintIter> {
        aromatic_system_constraints_iter(py, &self.0)
    }

    /// The `(key, constraint)` pairs, in canonical order.
    fn items(&self, py: Python<'_>) -> PyResult<AromaticSystemConstraintItemsIter> {
        aromatic_system_constraint_items(py, &self.0)
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    fn get(
        &self,
        py: Python<'_>,
        key: Py<AromaticSystemConstraintKey>,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        match self.0.get(key.bind(py).borrow().to_ast()) {
            Some(constraint) => Ok(into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(py, constraint)?,
            )?
            .into_any()),
            None => Ok(default.unwrap_or_else(|| py.None())),
        }
    }

    /// The constraint with the given key; raises `KeyError` if absent.
    fn __getitem__(
        &self,
        py: Python<'_>,
        key: Py<AromaticSystemConstraintKey>,
    ) -> PyResult<AromaticSystemConstraintAst> {
        match self.0.get(key.bind(py).borrow().to_ast()) {
            Some(constraint) => AromaticSystemConstraintAst::from_ast(py, constraint),
            None => Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            )),
        }
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    fn __delitem__(
        &mut self,
        py: Python<'_>,
        key: Py<AromaticSystemConstraintKey>,
    ) -> PyResult<()> {
        if self.0.remove(key.bind(py).borrow().to_ast()).is_some() {
            Ok(())
        } else {
            Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            ))
        }
    }

    fn __contains__(&self, py: Python<'_>, key: Py<AromaticSystemConstraintKey>) -> bool {
        self.0.contains(key.bind(py).borrow().to_ast())
    }

    /// The asserted total π-electron count; `Undetermined` when no `ElectronCount`
    /// constraint is present (mirroring the non-optional Rust accessor).
    #[getter]
    fn electron_count(&self, py: Python<'_>) -> PyResult<ValueAst> {
        ValueAst::from_ast(py, &self.0.electron_count())
    }

    #[setter]
    fn set_electron_count(&mut self, py: Python<'_>, value: ValueArg) {
        self.0.set(AstAromaticSystemConstraintAst::electron_count(
            value.to_ast(py),
        ));
    }

    /// The present constraints as a dict keyed by snake_case name; values are the
    /// inner-value mirrors.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        aromatic_system_constraints_asdict(py, &self.0)
    }
}

impl AromaticSystemConstraintsAst {
    /// The wrapped AST constraints — read access for aromatic system construction.
    pub(crate) fn inner(&self) -> &AstAromaticSystemConstraintsAst {
        &self.0
    }

    /// Wrap AST constraints (the hold-the-value `from_inner` bridge). Test-only —
    /// in-crate construction wraps `AromaticSystemConstraintsAst(..)` directly.
    #[cfg(test)]
    pub(crate) fn from_inner(constraints: AstAromaticSystemConstraintsAst) -> Self {
        AromaticSystemConstraintsAst(constraints)
    }
}

/// Build the per-constraint iterator handle from a borrowed container.
fn aromatic_system_constraints_iter(
    py: Python<'_>,
    constraints: &AstAromaticSystemConstraintsAst,
) -> PyResult<AromaticSystemConstraintIter> {
    let entries = constraints
        .iter()
        .map(|constraint| {
            into_py_variant(py, AromaticSystemConstraintAst::from_ast(py, constraint)?)
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(AromaticSystemConstraintIter {
        entries: entries.into_iter(),
    })
}

/// Build the key iterator handle from a borrowed container (mapping-style keys).
fn aromatic_system_constraint_keys(
    py: Python<'_>,
    constraints: &AstAromaticSystemConstraintsAst,
) -> PyResult<AromaticSystemConstraintKeyIter> {
    let keys = constraints
        .iter()
        .map(|constraint| {
            into_py_variant(py, AromaticSystemConstraintKey::from_ast(&constraint.key()))
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(AromaticSystemConstraintKeyIter {
        keys: keys.into_iter(),
    })
}

/// Build the item iterator handle (`(key, constraint)` pairs) from a borrowed container.
fn aromatic_system_constraint_items(
    py: Python<'_>,
    constraints: &AstAromaticSystemConstraintsAst,
) -> PyResult<AromaticSystemConstraintItemsIter> {
    let items = constraints
        .iter()
        .map(|constraint| {
            Ok((
                into_py_variant(py, AromaticSystemConstraintKey::from_ast(&constraint.key()))?,
                into_py_variant(py, AromaticSystemConstraintAst::from_ast(py, constraint)?)?,
            ))
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(AromaticSystemConstraintItemsIter {
        items: items.into_iter(),
    })
}

/// The present constraints as a dict keyed by snake_case name; values are the
/// inner-value mirrors.
fn aromatic_system_constraints_asdict<'py>(
    py: Python<'py>,
    constraints: &AstAromaticSystemConstraintsAst,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    for entry in constraints.iter() {
        match entry {
            AstAromaticSystemConstraintAst::ElectronCount(v) => {
                dict.set_item("electron_count", ValueAst::from_ast(py, v)?)?
            }
        }
    }
    Ok(dict)
}

/// What an `AromaticSystemConstraintsView` writes through to: an aromatic system
/// within a molecule (by index) or a standalone `AromaticSystemAst`.
enum AromaticSystemConstraintsBacking {
    Molecule {
        owner: Py<MoleculeAst>,
        id: AstAromaticSystemId,
    },
    AromaticSystem(Py<AromaticSystemAst>),
}

/// A live handle onto one aromatic system's constraints, backed by either a
/// molecule-system or a standalone `AromaticSystemAst`. Reads borrow the constraints
/// and read only the item they need (no whole-container clone); mutators write through
/// to the system in place, without a clone-and-writeback.
#[pyclass]
pub struct AromaticSystemConstraintsView {
    backing: AromaticSystemConstraintsBacking,
}

impl AromaticSystemConstraintsView {
    /// Borrow the backing system's constraints and read one item through `f` — no clone.
    fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&AstAromaticSystemConstraintsAst) -> PyResult<R>,
    ) -> PyResult<R> {
        match &self.backing {
            AromaticSystemConstraintsBacking::Molecule { owner, id } => {
                let molecule = owner.bind(py).borrow();
                let view = molecule
                    .inner()
                    .aromatic_systems()
                    .get(*id)
                    .ok_or_else(|| PyIndexError::new_err("aromatic system id out of range"))?;
                f(&view.ast.constraints)
            }
            AromaticSystemConstraintsBacking::AromaticSystem(system) => {
                let system = system.bind(py).borrow();
                f(&system.inner().constraints)
            }
        }
    }

    /// Mutate the backing system's constraints in place through `f`.
    fn with_mut<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut AstAromaticSystemConstraintsAst) -> R,
    ) -> R {
        match &self.backing {
            AromaticSystemConstraintsBacking::Molecule { owner, id } => f(&mut owner
                .borrow_mut(py)
                .inner_mut()
                .aromatic_system_mut(*id)
                .ast
                .constraints),
            AromaticSystemConstraintsBacking::AromaticSystem(system) => {
                f(&mut system.borrow_mut(py).inner_mut().constraints)
            }
        }
    }

    /// Set one constraint on the backing system in place (last-wins per key).
    fn set_ast(&self, py: Python<'_>, constraint: AstAromaticSystemConstraintAst) {
        self.with_mut(py, |cs| cs.set(constraint));
    }

    /// Remove one key from the backing system in place, returning the removed entry.
    fn remove_ast(
        &self,
        py: Python<'_>,
        key: AstAromaticSystemConstraintKey,
    ) -> Option<AstAromaticSystemConstraintAst> {
        self.with_mut(py, |cs| cs.remove(key))
    }
}

#[pymethods]
impl AromaticSystemConstraintsView {
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let count = self.read(py, |cs| Ok(cs.len()))?;
        Ok(format!("AromaticSystemConstraintsView({count} entries)"))
    }

    /// Insert `c` on the system in place, replacing any existing entry of the same key
    /// (last-wins).
    fn set(&self, py: Python<'_>, c: Py<AromaticSystemConstraintAst>) {
        self.set_ast(py, c.bind(py).borrow().to_ast(py));
    }

    /// Remove the entry with the given key from the system in place, returning it if
    /// present (dict `pop`).
    fn pop(
        &self,
        py: Python<'_>,
        key: Py<AromaticSystemConstraintKey>,
    ) -> PyResult<Option<AromaticSystemConstraintAst>> {
        self.remove_ast(py, key.bind(py).borrow().to_ast())
            .map(|c| AromaticSystemConstraintAst::from_ast(py, &c))
            .transpose()
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    fn __delitem__(&self, py: Python<'_>, key: Py<AromaticSystemConstraintKey>) -> PyResult<()> {
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

    /// Overlay `other` onto the system's constraints in place — another container or an
    /// iterable of `AromaticSystemConstraintAst` (last-wins per key; undetermined entries remove).
    fn update(&self, py: Python<'_>, other: AromaticSystemConstraintsUpdate) -> PyResult<()> {
        self.with_mut(py, |cs| other.apply(py, cs))
    }

    fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        self.read(py, |cs| Ok(cs.len()))
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    fn __iter__(&self, py: Python<'_>) -> PyResult<AromaticSystemConstraintKeyIter> {
        self.read(py, |cs| aromatic_system_constraint_keys(py, cs))
    }

    /// The constraint keys, in canonical order.
    fn keys(&self, py: Python<'_>) -> PyResult<AromaticSystemConstraintKeyIter> {
        self.read(py, |cs| aromatic_system_constraint_keys(py, cs))
    }

    /// The constraints, in canonical order.
    fn values(&self, py: Python<'_>) -> PyResult<AromaticSystemConstraintIter> {
        self.read(py, |cs| aromatic_system_constraints_iter(py, cs))
    }

    /// The `(key, constraint)` pairs, in canonical order.
    fn items(&self, py: Python<'_>) -> PyResult<AromaticSystemConstraintItemsIter> {
        self.read(py, |cs| aromatic_system_constraint_items(py, cs))
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    fn get(
        &self,
        py: Python<'_>,
        key: Py<AromaticSystemConstraintKey>,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        let key = key.bind(py).borrow().to_ast();
        let found = self.read(py, |cs| {
            cs.get(key)
                .map(|constraint| AromaticSystemConstraintAst::from_ast(py, constraint))
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
        key: Py<AromaticSystemConstraintKey>,
    ) -> PyResult<AromaticSystemConstraintAst> {
        let ast_key = key.bind(py).borrow().to_ast();
        let found = self.read(py, |cs| {
            cs.get(ast_key)
                .map(|constraint| AromaticSystemConstraintAst::from_ast(py, constraint))
                .transpose()
        })?;
        match found {
            Some(constraint) => Ok(constraint),
            None => Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            )),
        }
    }

    fn __contains__(&self, py: Python<'_>, key: Py<AromaticSystemConstraintKey>) -> PyResult<bool> {
        let key = key.bind(py).borrow().to_ast();
        self.read(py, |cs| Ok(cs.contains(key)))
    }

    /// The asserted total π-electron count; `Undetermined` when no `ElectronCount`
    /// constraint is present (mirroring the non-optional Rust accessor).
    #[getter]
    fn electron_count(&self, py: Python<'_>) -> PyResult<ValueAst> {
        self.read(py, |cs| ValueAst::from_ast(py, &cs.electron_count()))
    }

    #[setter]
    fn set_electron_count(&self, py: Python<'_>, value: ValueArg) {
        self.set_ast(
            py,
            AstAromaticSystemConstraintAst::electron_count(value.to_ast(py)),
        );
    }

    /// The present constraints as a dict keyed by snake_case name.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        self.read(py, |cs| aromatic_system_constraints_asdict(py, cs))
    }
}

#[pyclass]
struct AromaticSystemConstraintIter {
    entries: IntoIter<Py<AromaticSystemConstraintAst>>,
}

#[pymethods]
impl AromaticSystemConstraintIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<Py<AromaticSystemConstraintAst>> {
        self.entries.next()
    }
}

#[pyclass]
struct AromaticSystemConstraintKeyIter {
    keys: IntoIter<Py<AromaticSystemConstraintKey>>,
}

#[pymethods]
impl AromaticSystemConstraintKeyIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<Py<AromaticSystemConstraintKey>> {
        self.keys.next()
    }
}

#[pyclass]
struct AromaticSystemConstraintItemsIter {
    items: IntoIter<(
        Py<AromaticSystemConstraintKey>,
        Py<AromaticSystemConstraintAst>,
    )>,
}

#[pymethods]
impl AromaticSystemConstraintItemsIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(
        &mut self,
    ) -> Option<(
        Py<AromaticSystemConstraintKey>,
        Py<AromaticSystemConstraintAst>,
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

    /// Benzene: six aromatic carbons (atom ids 0–5), one aromatic system over all six
    /// (electrons `[1,1,1,1,1,1]`), aromatic system id 0.
    fn benzene(py: Python<'_>) -> Py<MoleculeAst> {
        let molecule = AstMoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AstAtomAst::from_element(ChemElement::C); 6],
            aromatic: vec![(
                (0u32..6).map(AstAtomId).collect(),
                AstAromaticSystemAst::from_electrons(vec![1, 1, 1, 1, 1, 1]),
            )],
            ..Default::default()
        });
        Py::new(py, MoleculeAst::from_inner(molecule)).unwrap()
    }

    #[rstest]
    fn test_aromatic_system_ast_new() {
        Python::attach(|py| {
            let spin_ast = AstSpinStateAst::from((0_u8, 1_u8));
            let spin = Py::new(py, SpinStateAst::from_ast(py, &spin_ast).unwrap()).unwrap();
            let system = AromaticSystemAst::new(
                py,
                ElectronCountsArg::Lit(vec![1, 1, 1]),
                Some(ValueArg::Lit(-2)),
                Some(spin.bind(py).borrow()),
                None,
            );
            assert_eq!(
                system.inner().electrons,
                AstElectronCountsAst::Lit(vec![1, 1, 1])
            );
            assert_eq!(system.inner().charge, AstValueAst::Lit(-2));
            assert_eq!(system.inner().spin, spin_ast);
        });
    }

    #[rstest]
    fn test_aromatic_system_ast_new_constraints() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(
                    py,
                    &AstAromaticSystemConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = Py::new(py, AromaticSystemConstraintsAst::new(py, vec![ec])).unwrap();
            let system = AromaticSystemAst::new(
                py,
                ElectronCountsArg::Lit(vec![1, 1, 1]),
                None,
                None,
                Some(constraints),
            );
            assert_eq!(
                system.inner().constraints.electron_count(),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    #[case::undetermined("*")]
    #[case::electron_count("[1,1,1]#e6")]
    #[case::charge("[1,1,1]#c-2")]
    fn test_aromatic_system_ast_parse(#[case] dsl: &str) {
        let system = AromaticSystemAst::parse(dsl).unwrap();
        assert_eq!(system.__str__(), dsl);
        assert_eq!(
            system.__repr__(),
            format!("AromaticSystemAst.parse('{dsl}')")
        );
    }

    #[rstest]
    fn test_aromatic_system_ast_parse_error() {
        assert!(AromaticSystemAst::parse("z").is_err());
    }

    #[rstest]
    fn test_aromatic_system_ast_electrons() {
        Python::attach(|py| {
            let mut system =
                AromaticSystemAst::from_inner(AstAromaticSystemAst::from_electrons(vec![1, 1, 1]));
            assert_eq!(
                system.electrons().to_ast(),
                AstElectronCountsAst::Lit(vec![1, 1, 1])
            );
            system.set_electrons(py, ElectronCountsArg::Lit(vec![2, 2]));
            assert_eq!(
                system.electrons().to_ast(),
                AstElectronCountsAst::Lit(vec![2, 2])
            );
        });
    }

    #[rstest]
    fn test_aromatic_system_ast_charge() {
        Python::attach(|py| {
            let mut system =
                AromaticSystemAst::from_inner(AstAromaticSystemAst::from_electrons(vec![1, 1, 1]));
            system.set_charge(py, ValueArg::Lit(-1));
            assert_eq!(system.charge(py).unwrap().to_ast(py), AstValueAst::Lit(-1));
        });
    }

    #[rstest]
    fn test_aromatic_system_ast_spin() {
        Python::attach(|py| {
            let spin_ast = AstSpinStateAst::from((0_u8, 1_u8));
            let spin = Py::new(py, SpinStateAst::from_ast(py, &spin_ast).unwrap()).unwrap();
            let mut system =
                AromaticSystemAst::from_inner(AstAromaticSystemAst::from_electrons(vec![1, 1, 1]));
            system.set_spin(py, spin.bind(py).borrow());
            assert_eq!(system.spin(py).unwrap().to_ast(py), spin_ast);
        });
    }

    #[rstest]
    fn test_aromatic_system_ast_set_constraints_from_view() {
        Python::attach(|py| {
            let src = Py::new(
                py,
                AromaticSystemAst::from_inner(
                    AstAromaticSystemAst::from_electrons(vec![1, 1, 1])
                        .with_constraint(AstAromaticSystemConstraintAst::electron_count(6)),
                ),
            )
            .unwrap();
            let view = Py::new(
                py,
                AromaticSystemConstraintsView {
                    backing: AromaticSystemConstraintsBacking::AromaticSystem(src),
                },
            )
            .unwrap();
            let mut dst =
                AromaticSystemAst::from_inner(AstAromaticSystemAst::from_electrons(vec![1, 1, 1]));
            dst.set_constraints(py, AromaticSystemConstraintsArg::View(view))
                .unwrap();
            assert_eq!(
                dst.inner().constraints.electron_count(),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_aromatic_system_ast_asdict() {
        Python::attach(|py| {
            let system = AromaticSystemAst::from_inner(
                AstAromaticSystemAst::from_electrons(vec![1, 1, 1])
                    .with_constraint(AstAromaticSystemConstraintAst::electron_count(6)),
            );
            let dict = system.asdict(py).unwrap();
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
    fn test_aromatic_system_view_atom_ids() {
        Python::attach(|py| {
            let view = AromaticSystemView {
                owner: benzene(py),
                id: AstAromaticSystemId(0),
            };
            assert_eq!(view.id(), 0);
            let atom_ids: Vec<u32> = view.atom_ids(py).unwrap().extract().unwrap();
            assert_eq!(atom_ids, vec![0, 1, 2, 3, 4, 5]);
            assert_eq!(view.__repr__(), "AromaticSystemView(id=0)");
        });
    }

    #[rstest]
    fn test_aromatic_system_view_electrons() {
        Python::attach(|py| {
            let owner = benzene(py);
            let view = AromaticSystemView {
                owner: owner.clone_ref(py),
                id: AstAromaticSystemId(0),
            };
            assert_eq!(
                view.electrons(py).unwrap().to_ast(),
                AstElectronCountsAst::Lit(vec![1, 1, 1, 1, 1, 1])
            );
            view.set_electrons(py, ElectronCountsArg::Lit(vec![2, 2, 2, 2, 2, 2]));
            let fresh = AromaticSystemView {
                owner,
                id: AstAromaticSystemId(0),
            };
            assert_eq!(
                fresh.electrons(py).unwrap().to_ast(),
                AstElectronCountsAst::Lit(vec![2, 2, 2, 2, 2, 2])
            );
        });
    }

    #[rstest]
    fn test_aromatic_system_view_charge() {
        Python::attach(|py| {
            let owner = benzene(py);
            let view = AromaticSystemView {
                owner: owner.clone_ref(py),
                id: AstAromaticSystemId(0),
            };
            view.set_charge(py, ValueArg::Lit(-1));
            let fresh = AromaticSystemView {
                owner,
                id: AstAromaticSystemId(0),
            };
            assert_eq!(fresh.charge(py).unwrap().to_ast(py), AstValueAst::Lit(-1));
        });
    }

    #[rstest]
    fn test_aromatic_system_view_spin() {
        Python::attach(|py| {
            let spin_ast = AstSpinStateAst::from((0_u8, 1_u8));
            let spin = Py::new(py, SpinStateAst::from_ast(py, &spin_ast).unwrap()).unwrap();
            let owner = benzene(py);
            let view = AromaticSystemView {
                owner: owner.clone_ref(py),
                id: AstAromaticSystemId(0),
            };
            view.set_spin(py, spin.bind(py).borrow());
            let fresh = AromaticSystemView {
                owner,
                id: AstAromaticSystemId(0),
            };
            assert_eq!(fresh.spin(py).unwrap().to_ast(py), spin_ast);
        });
    }

    #[rstest]
    fn test_aromatic_system_view_constraints() {
        Python::attach(|py| {
            let view = AromaticSystemView {
                owner: benzene(py),
                id: AstAromaticSystemId(0),
            };
            match view.constraints(py).backing {
                AromaticSystemConstraintsBacking::Molecule { id, .. } => {
                    assert_eq!(id, AstAromaticSystemId(0))
                }
                _ => panic!("expected molecule-backed view"),
            }
        });
    }

    #[rstest]
    fn test_aromatic_system_view_set_constraints() {
        Python::attach(|py| {
            let owner = benzene(py);
            let view = AromaticSystemView {
                owner: owner.clone_ref(py),
                id: AstAromaticSystemId(0),
            };
            let constraints = Py::new(
                py,
                AromaticSystemConstraintsAst::new(
                    py,
                    vec![into_py_variant(
                        py,
                        AromaticSystemConstraintAst::from_ast(
                            py,
                            &AstAromaticSystemConstraintAst::electron_count(6),
                        )
                        .unwrap(),
                    )
                    .unwrap()],
                ),
            )
            .unwrap();
            view.set_constraints(py, AromaticSystemConstraintsArg::Container(constraints))
                .unwrap();
            let fresh = AromaticSystemView {
                owner,
                id: AstAromaticSystemId(0),
            };
            assert_eq!(
                fresh.constraints(py).electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_aromatic_system_view_asdict() {
        Python::attach(|py| {
            let view = AromaticSystemView {
                owner: benzene(py),
                id: AstAromaticSystemId(0),
            };
            let dict = view.asdict(py).unwrap();
            assert_eq!(dict.len(), 4);
            let electrons = dict.get_item("electrons").unwrap().unwrap();
            let expected =
                into_py_variant(py, ElectronCountsAst::Lit(vec![1, 1, 1, 1, 1, 1])).unwrap();
            assert!(electrons.eq(expected.bind(py)).unwrap());
            assert!(dict.contains("charge").unwrap());
            assert!(dict.contains("spin").unwrap());
            assert!(dict.contains("constraints").unwrap());
        });
    }

    #[rstest]
    fn test_aromatic_system_constraint_key_roundtrip() {
        let key =
            AromaticSystemConstraintKey::from_ast(&AstAromaticSystemConstraintKey::ElectronCount);
        assert_eq!(key.to_ast(), AstAromaticSystemConstraintKey::ElectronCount);
    }

    #[rstest]
    fn test_aromatic_system_constraint_ast_key() {
        Python::attach(|py| {
            let constraint = AstAromaticSystemConstraintAst::electron_count(6);
            let key = AromaticSystemConstraintAst::from_ast(py, &constraint)
                .unwrap()
                .key(py);
            assert_eq!(key.to_ast(), AstAromaticSystemConstraintKey::ElectronCount);
        });
    }

    #[rstest]
    #[case(AstAromaticSystemConstraintAst::electron_count(6))]
    #[case(AstAromaticSystemConstraintAst::electron_count(AstValueAst::Undetermined))]
    fn test_aromatic_system_constraint_ast_roundtrip(#[case] ast: AstAromaticSystemConstraintAst) {
        Python::attach(|py| {
            assert_eq!(
                AromaticSystemConstraintAst::from_ast(py, &ast)
                    .unwrap()
                    .to_ast(py),
                ast
            );
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_ast_new() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(
                    py,
                    &AstAromaticSystemConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = AromaticSystemConstraintsAst::new(py, vec![ec]);
            assert_eq!(constraints.__len__(), 1);
            assert_eq!(
                constraints.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_ast_repr() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(
                    py,
                    &AstAromaticSystemConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = AromaticSystemConstraintsAst::new(py, vec![ec]);
            assert_eq!(
                constraints.__repr__(py).unwrap(),
                "AromaticSystemConstraintsAst([AromaticSystemConstraintAst.ElectronCount(ValueAst.Lit(6))])"
            );
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_ast_set() {
        Python::attach(|py| {
            let mut constraints = AromaticSystemConstraintsAst::new(py, vec![]);
            let ec = into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(
                    py,
                    &AstAromaticSystemConstraintAst::electron_count(6),
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
    fn test_aromatic_system_constraints_ast_pop() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(
                    py,
                    &AstAromaticSystemConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let mut constraints = AromaticSystemConstraintsAst::new(py, vec![ec]);
            let key = into_py_variant(py, AromaticSystemConstraintKey::ElectronCount()).unwrap();
            let removed = constraints.pop(py, key).unwrap();
            match removed {
                Some(AromaticSystemConstraintAst::ElectronCount(v)) => {
                    assert_eq!(v.bind(py).borrow().to_ast(py), AstValueAst::Lit(6))
                }
                _ => panic!("expected removed ElectronCount(Lit(6))"),
            }
            assert_eq!(constraints.__len__(), 0);
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_ast_update() {
        Python::attach(|py| {
            let mut constraints = AromaticSystemConstraintsAst::new(py, vec![]);
            let mut other = AstAromaticSystemConstraintsAst::new();
            other.set(AstAromaticSystemConstraintAst::electron_count(6));
            constraints
                .update(
                    py,
                    AromaticSystemConstraintsUpdate::Container(
                        Py::new(py, AromaticSystemConstraintsAst::from_inner(other)).unwrap(),
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
    fn test_aromatic_system_constraints_ast_update_entries() {
        Python::attach(|py| {
            let mut constraints = AromaticSystemConstraintsAst::new(py, vec![]);
            let ec = into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(
                    py,
                    &AstAromaticSystemConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            constraints
                .update(py, AromaticSystemConstraintsUpdate::Entries(vec![ec]))
                .unwrap();
            assert_eq!(constraints.__len__(), 1);
            assert_eq!(
                constraints.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_ast_len_contains() {
        Python::attach(|py| {
            let empty = AromaticSystemConstraintsAst::new(py, vec![]);
            assert_eq!(empty.__len__(), 0);
            assert!(!empty.__contains__(
                py,
                into_py_variant(py, AromaticSystemConstraintKey::ElectronCount()).unwrap()
            ));
            let ec = into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(
                    py,
                    &AstAromaticSystemConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = AromaticSystemConstraintsAst::new(py, vec![ec]);
            assert!(constraints.__contains__(
                py,
                into_py_variant(py, AromaticSystemConstraintKey::ElectronCount()).unwrap()
            ));
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_ast_keys_values_items() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(
                    py,
                    &AstAromaticSystemConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = AromaticSystemConstraintsAst::new(py, vec![ec]);

            let mut keys = constraints.__iter__(py).unwrap();
            assert_eq!(
                keys.__next__().unwrap().bind(py).borrow().to_ast(),
                AstAromaticSystemConstraintKey::ElectronCount
            );
            assert!(keys.__next__().is_none());

            let mut values = constraints.values(py).unwrap();
            assert_eq!(
                values.__next__().unwrap().bind(py).borrow().to_ast(py),
                AstAromaticSystemConstraintAst::electron_count(6)
            );

            let mut items = constraints.items(py).unwrap();
            let (key, value) = items.__next__().unwrap();
            assert_eq!(
                key.bind(py).borrow().to_ast(),
                AstAromaticSystemConstraintKey::ElectronCount
            );
            assert_eq!(
                value.bind(py).borrow().to_ast(py),
                AstAromaticSystemConstraintAst::electron_count(6)
            );
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_ast_get() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(
                    py,
                    &AstAromaticSystemConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = AromaticSystemConstraintsAst::new(py, vec![ec]);
            let present = constraints
                .get(
                    py,
                    into_py_variant(py, AromaticSystemConstraintKey::ElectronCount()).unwrap(),
                    None,
                )
                .unwrap();
            let expected = into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(
                    py,
                    &AstAromaticSystemConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap()
            .into_any();
            assert!(present.bind(py).eq(expected.bind(py)).unwrap());

            let empty = AromaticSystemConstraintsAst::new(py, vec![]);
            let absent = empty
                .get(
                    py,
                    into_py_variant(py, AromaticSystemConstraintKey::ElectronCount()).unwrap(),
                    None,
                )
                .unwrap();
            assert!(absent.bind(py).is_none());

            // a caller-supplied default is returned verbatim when the key is absent
            let sentinel = into_py_variant(py, AromaticSystemConstraintKey::ElectronCount())
                .unwrap()
                .into_any();
            let defaulted = empty
                .get(
                    py,
                    into_py_variant(py, AromaticSystemConstraintKey::ElectronCount()).unwrap(),
                    Some(sentinel.clone_ref(py)),
                )
                .unwrap();
            assert_eq!(defaulted.as_ptr(), sentinel.as_ptr());
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_ast_electron_count() {
        Python::attach(|py| {
            let empty = AromaticSystemConstraintsAst::new(py, vec![]);
            assert_eq!(
                empty.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Undetermined
            );
            let ec = into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(
                    py,
                    &AstAromaticSystemConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = AromaticSystemConstraintsAst::new(py, vec![ec]);
            assert_eq!(
                constraints.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_ast_set_electron_count() {
        Python::attach(|py| {
            let mut constraints = AromaticSystemConstraintsAst::new(py, vec![]);
            constraints.set_electron_count(py, ValueArg::Lit(6));
            assert_eq!(
                constraints.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_ast_getitem_error() {
        Python::attach(|py| {
            let constraints = AromaticSystemConstraintsAst::new(py, vec![]);
            let key = into_py_variant(py, AromaticSystemConstraintKey::ElectronCount()).unwrap();
            assert!(constraints.__getitem__(py, key).is_err());
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_ast_delitem_error() {
        Python::attach(|py| {
            let mut constraints = AromaticSystemConstraintsAst::new(py, vec![]);
            let key = into_py_variant(py, AromaticSystemConstraintKey::ElectronCount()).unwrap();
            assert!(constraints.__delitem__(py, key).is_err());
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_ast_asdict() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(
                    py,
                    &AstAromaticSystemConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = AromaticSystemConstraintsAst::new(py, vec![ec]);
            let dict = constraints.asdict(py).unwrap();
            assert_eq!(dict.len(), 1);
            let value = dict.get_item("electron_count").unwrap().unwrap();
            let expected =
                into_py_variant(py, ValueAst::from_ast(py, &AstValueAst::Lit(6)).unwrap()).unwrap();
            assert!(value.eq(expected.bind(py)).unwrap());
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_view_set() {
        Python::attach(|py| {
            let system = Py::new(
                py,
                AromaticSystemAst::from_inner(AstAromaticSystemAst::from_electrons(vec![1, 1, 1])),
            )
            .unwrap();
            let view = AromaticSystemConstraintsView {
                backing: AromaticSystemConstraintsBacking::AromaticSystem(system.clone_ref(py)),
            };
            let ec = into_py_variant(
                py,
                AromaticSystemConstraintAst::from_ast(
                    py,
                    &AstAromaticSystemConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            view.set(py, ec);
            // a fresh view proves the write hit the standalone system, not a copy
            let fresh = AromaticSystemConstraintsView {
                backing: AromaticSystemConstraintsBacking::AromaticSystem(system),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 1);
            assert_eq!(
                fresh.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_view_pop() {
        Python::attach(|py| {
            let system = Py::new(
                py,
                AromaticSystemAst::from_inner(
                    AstAromaticSystemAst::from_electrons(vec![1, 1, 1])
                        .with_constraint(AstAromaticSystemConstraintAst::electron_count(6)),
                ),
            )
            .unwrap();
            let view = AromaticSystemConstraintsView {
                backing: AromaticSystemConstraintsBacking::AromaticSystem(system.clone_ref(py)),
            };
            let removed = view
                .pop(
                    py,
                    into_py_variant(py, AromaticSystemConstraintKey::ElectronCount()).unwrap(),
                )
                .unwrap();
            match removed {
                Some(AromaticSystemConstraintAst::ElectronCount(v)) => {
                    assert_eq!(v.bind(py).borrow().to_ast(py), AstValueAst::Lit(6))
                }
                _ => panic!("expected removed ElectronCount(Lit(6))"),
            }
            let fresh = AromaticSystemConstraintsView {
                backing: AromaticSystemConstraintsBacking::AromaticSystem(system),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 0);
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_view_update() {
        Python::attach(|py| {
            let system = Py::new(
                py,
                AromaticSystemAst::from_inner(AstAromaticSystemAst::from_electrons(vec![1, 1, 1])),
            )
            .unwrap();
            let view = AromaticSystemConstraintsView {
                backing: AromaticSystemConstraintsBacking::AromaticSystem(system.clone_ref(py)),
            };
            let mut other = AstAromaticSystemConstraintsAst::new();
            other.set(AstAromaticSystemConstraintAst::electron_count(6));
            view.update(
                py,
                AromaticSystemConstraintsUpdate::Container(
                    Py::new(py, AromaticSystemConstraintsAst::from_inner(other)).unwrap(),
                ),
            )
            .unwrap();
            let fresh = AromaticSystemConstraintsView {
                backing: AromaticSystemConstraintsBacking::AromaticSystem(system),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 1);
            assert_eq!(
                fresh.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_view_set_electron_count() {
        Python::attach(|py| {
            let system = Py::new(
                py,
                AromaticSystemAst::from_inner(AstAromaticSystemAst::from_electrons(vec![1, 1, 1])),
            )
            .unwrap();
            let view = AromaticSystemConstraintsView {
                backing: AromaticSystemConstraintsBacking::AromaticSystem(system.clone_ref(py)),
            };
            assert_eq!(
                view.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Undetermined
            );
            view.set_electron_count(py, ValueArg::Lit(6));
            let fresh = AromaticSystemConstraintsView {
                backing: AromaticSystemConstraintsBacking::AromaticSystem(system),
            };
            assert_eq!(
                fresh.electron_count(py).unwrap().to_ast(py),
                AstValueAst::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_aromatic_system_constraints_view_set_molecule_backed() {
        Python::attach(|py| {
            let owner = benzene(py);
            let view = AromaticSystemConstraintsView {
                backing: AromaticSystemConstraintsBacking::Molecule {
                    owner: owner.clone_ref(py),
                    id: AstAromaticSystemId(0),
                },
            };
            view.set_electron_count(py, ValueArg::Lit(6));
            let fresh = AromaticSystemConstraintsView {
                backing: AromaticSystemConstraintsBacking::Molecule {
                    owner,
                    id: AstAromaticSystemId(0),
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
