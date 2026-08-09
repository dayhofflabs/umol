//! Atom field values, owned atom ASTs, and molecule-backed atom views.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use std::collections::BTreeSet;
use std::str::FromStr;
use std::vec::IntoIter;

use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use umol_chem::element::Element as ChemElement;
use umol_graph_ir::ir::{
    AsLit, AtomForm as GraphIrAtomForm, AtomId as GraphIrAtomId, AtomUpdate as GraphIrAtomUpdate,
    ElementForm as GraphIrElementForm, IsotopeMass as GraphIrIsotopeMass,
    IsotopeMassForm as GraphIrIsotopeMassForm, Molecule as GraphIrMolecule,
};

use crate::constraint::atom::{
    atom_constraints_asdict, AtomConstraintsAst, AtomConstraintsBacking, AtomConstraintsLike,
    AtomConstraintsView,
};
use crate::convert::{hash_rust, variant_repr};
use crate::element::Element;
use crate::error::parse_error;
use crate::lattice::impl_py_lattice;
use crate::molecule::MoleculeAst;
use crate::spin::{UnpairedElectronsForm, UnpairedElectronsUpdate};
use crate::value::{MemOp, NumForm, NumLike};

/// Element expression: undetermined, a single element, a finite element set, a
/// complement set (`!{…}`), or a variable with an optional membership restriction.
#[pyclass]
pub enum ElementForm {
    Undetermined(),
    Lit(Element),
    LitSet(BTreeSet<Element>),
    NotSet(BTreeSet<Element>),
    Var(String, Option<(MemOp, BTreeSet<Element>)>),
}

#[pymethods]
impl ElementForm {
    /// The single element this resolves to, or `None` when it is not a bare
    /// literal (undetermined, a set, a complement, or a variable).
    fn as_lit(&self) -> Option<Element> {
        self.to_rust().as_lit().map(Element::from)
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            ElementForm::Undetermined() => ("Undetermined", 0),
            ElementForm::Lit(_) => ("Lit", 1),
            ElementForm::LitSet(_) => ("LitSet", 1),
            ElementForm::NotSet(_) => ("NotSet", 1),
            ElementForm::Var(_, _) => ("Var", 2),
        };
        variant_repr(slf.bind(py).as_any(), "ElementForm", variant, arity)
    }
}

impl_py_lattice!(
    ElementForm,
    GraphIrElementForm,
    |value: &ElementForm, _py: Python<'_>| -> PyResult<GraphIrElementForm> { Ok(value.to_rust()) },
    |_py: Python<'_>, value: GraphIrElementForm| -> PyResult<ElementForm> {
        Ok(ElementForm::from_rust(&value))
    }
);

impl ElementForm {
    pub(crate) fn from_rust(ast: &GraphIrElementForm) -> ElementForm {
        match ast {
            GraphIrElementForm::Undetermined => ElementForm::Undetermined(),
            GraphIrElementForm::Lit(e) => ElementForm::Lit(Element::from(*e)),
            GraphIrElementForm::LitSet(members) => {
                ElementForm::LitSet(members.iter().copied().map(Element::from).collect())
            }
            GraphIrElementForm::NotSet(members) => {
                ElementForm::NotSet(members.iter().copied().map(Element::from).collect())
            }
            GraphIrElementForm::Var(boxed) => {
                let (name, restriction) = &**boxed;
                ElementForm::Var(
                    name.clone(),
                    restriction.as_ref().map(|(op, members)| {
                        (
                            MemOp::from_rust(*op),
                            members.iter().copied().map(Element::from).collect(),
                        )
                    }),
                )
            }
        }
    }

    pub(crate) fn to_rust(&self) -> GraphIrElementForm {
        match self {
            ElementForm::Undetermined() => GraphIrElementForm::Undetermined,
            ElementForm::Lit(e) => GraphIrElementForm::Lit(ChemElement::from(e)),
            ElementForm::LitSet(members) => GraphIrElementForm::LitSet(Box::new(
                members.iter().map(ChemElement::from).collect(),
            )),
            ElementForm::NotSet(members) => GraphIrElementForm::NotSet(Box::new(
                members.iter().map(ChemElement::from).collect(),
            )),
            ElementForm::Var(name, restriction) => GraphIrElementForm::Var(Box::new((
                name.clone(),
                restriction.as_ref().map(|(op, members)| {
                    (
                        op.to_rust(),
                        members.iter().map(ChemElement::from).collect(),
                    )
                }),
            ))),
        }
    }
}

/// Exact ground isotope mass: the natural isotopic mixture or a specific mass number.
#[pyclass(from_py_object)]
#[derive(Clone, Copy)]
pub enum IsotopeMass {
    Natural(),
    MassNumber(u32),
}

#[pymethods]
impl IsotopeMass {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            Self::Natural() => ("Natural", 0),
            Self::MassNumber(_) => ("MassNumber", 1),
        };
        variant_repr(slf.bind(py).as_any(), "IsotopeMass", variant, arity)
    }
}

impl IsotopeMass {
    pub(crate) fn from_rust(mass: GraphIrIsotopeMass) -> Self {
        match mass {
            GraphIrIsotopeMass::Natural => Self::Natural(),
            GraphIrIsotopeMass::MassNumber(mass) => Self::MassNumber(mass),
        }
    }

    pub(crate) fn to_rust(self) -> GraphIrIsotopeMass {
        match self {
            Self::Natural() => GraphIrIsotopeMass::Natural,
            Self::MassNumber(mass) => GraphIrIsotopeMass::MassNumber(mass),
        }
    }
}

/// Isotope-mass expression: undetermined, the natural isotopic mixture, a single
/// mass number, a finite mass set, or a variable with an optional mass-set restriction.
#[pyclass]
pub enum IsotopeMassForm {
    Undetermined(),
    Natural(),
    Lit(u32),
    LitSet(BTreeSet<u32>),
    Var(String, Option<BTreeSet<u32>>),
}

#[pymethods]
impl IsotopeMassForm {
    /// The exact isotope-mass value, or `None` when this expression is not ground.
    fn as_lit(&self) -> Option<IsotopeMass> {
        self.to_rust().as_lit().map(IsotopeMass::from_rust)
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            IsotopeMassForm::Undetermined() => ("Undetermined", 0),
            IsotopeMassForm::Natural() => ("Natural", 0),
            IsotopeMassForm::Lit(_) => ("Lit", 1),
            IsotopeMassForm::LitSet(_) => ("LitSet", 1),
            IsotopeMassForm::Var(_, _) => ("Var", 2),
        };
        variant_repr(slf.bind(py).as_any(), "IsotopeMassForm", variant, arity)
    }
}

impl IsotopeMassForm {
    pub(crate) fn from_rust(ast: &GraphIrIsotopeMassForm) -> IsotopeMassForm {
        match ast {
            GraphIrIsotopeMassForm::Undetermined => IsotopeMassForm::Undetermined(),
            GraphIrIsotopeMassForm::Natural => IsotopeMassForm::Natural(),
            GraphIrIsotopeMassForm::Lit(mass) => IsotopeMassForm::Lit(*mass),
            GraphIrIsotopeMassForm::LitSet(masses) => IsotopeMassForm::LitSet((**masses).clone()),
            GraphIrIsotopeMassForm::Var(boxed) => {
                let (name, restriction) = &**boxed;
                IsotopeMassForm::Var(name.clone(), restriction.clone())
            }
        }
    }

    pub(crate) fn to_rust(&self) -> GraphIrIsotopeMassForm {
        match self {
            IsotopeMassForm::Undetermined() => GraphIrIsotopeMassForm::Undetermined,
            IsotopeMassForm::Natural() => GraphIrIsotopeMassForm::Natural,
            IsotopeMassForm::Lit(mass) => GraphIrIsotopeMassForm::Lit(*mass),
            IsotopeMassForm::LitSet(masses) => {
                GraphIrIsotopeMassForm::LitSet(Box::new(masses.clone()))
            }
            IsotopeMassForm::Var(name, restriction) => {
                GraphIrIsotopeMassForm::Var(Box::new((name.clone(), restriction.clone())))
            }
        }
    }
}

impl_py_lattice!(
    IsotopeMassForm,
    GraphIrIsotopeMassForm,
    |value: &IsotopeMassForm, _py: Python<'_>| -> PyResult<GraphIrIsotopeMassForm> {
        Ok(value.to_rust())
    },
    |_py: Python<'_>, value: GraphIrIsotopeMassForm| -> PyResult<IsotopeMassForm> {
        Ok(IsotopeMassForm::from_rust(&value))
    }
);

/// Attribute updates for an atom. Omitted scalar fields remain unchanged;
/// constraints form a keyed update where undetermined values remove their key.
#[pyclass(frozen, skip_from_py_object)]
#[derive(Clone)]
pub struct AtomUpdate(GraphIrAtomUpdate);

#[pymethods]
impl AtomUpdate {
    #[new]
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (*, element=None, isotope_mass=None, charge=None, implicit_hydrogens=None, lone_pairs=None, unpaired_electrons=None, constraints=None))]
    fn new(
        py: Python<'_>,
        element: Option<ElementLike>,
        isotope_mass: Option<IsotopeMassLike>,
        charge: Option<NumLike>,
        implicit_hydrogens: Option<NumLike>,
        lone_pairs: Option<NumLike>,
        unpaired_electrons: Option<PyRef<'_, UnpairedElectronsUpdate>>,
        constraints: Option<Py<AtomConstraintsAst>>,
    ) -> Self {
        Self::from_rust(&GraphIrAtomUpdate {
            element: element.map(|value| value.to_rust(py)),
            isotope_mass: isotope_mass.map(|value| value.to_rust(py)),
            charge: charge.map(|value| value.to_rust(py)),
            implicit_hydrogens: implicit_hydrogens.map(|value| value.to_rust(py)),
            lone_pairs: lone_pairs.map(|value| value.to_rust(py)),
            unpaired_electrons: unpaired_electrons
                .map(|value| value.to_rust(py))
                .unwrap_or_default(),
            constraints: constraints
                .map(|value| value.bind(py).borrow().inner().clone())
                .unwrap_or_default(),
        })
    }

    /// Parse an atom-update DSL string into an `AtomUpdate`.
    #[staticmethod]
    fn parse(s: &str) -> PyResult<Self> {
        GraphIrAtomUpdate::from_str(s)
            .map(Self)
            .map_err(parse_error)
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("AtomUpdate.parse('{}')", self.0)
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    #[getter]
    fn element(&self) -> Option<ElementForm> {
        self.0.element.as_ref().map(ElementForm::from_rust)
    }

    #[getter]
    fn isotope_mass(&self) -> Option<IsotopeMassForm> {
        self.0.isotope_mass.as_ref().map(IsotopeMassForm::from_rust)
    }

    #[getter]
    fn charge(&self, py: Python<'_>) -> PyResult<Option<NumForm>> {
        self.0
            .charge
            .as_ref()
            .map(|value| NumForm::from_rust(py, value))
            .transpose()
    }

    #[getter]
    fn implicit_hydrogens(&self, py: Python<'_>) -> PyResult<Option<NumForm>> {
        self.0
            .implicit_hydrogens
            .as_ref()
            .map(|value| NumForm::from_rust(py, value))
            .transpose()
    }

    #[getter]
    fn lone_pairs(&self, py: Python<'_>) -> PyResult<Option<NumForm>> {
        self.0
            .lone_pairs
            .as_ref()
            .map(|value| NumForm::from_rust(py, value))
            .transpose()
    }

    #[getter]
    fn unpaired_electrons(&self, py: Python<'_>) -> PyResult<UnpairedElectronsUpdate> {
        UnpairedElectronsUpdate::from_rust(py, &self.0.unpaired_electrons)
    }

    #[getter]
    fn constraints(&self) -> AtomConstraintsAst {
        AtomConstraintsAst::from_inner(self.0.constraints.clone())
    }
}

impl AtomUpdate {
    pub(crate) fn from_rust(update: &GraphIrAtomUpdate) -> Self {
        Self(update.clone())
    }

    pub(crate) fn to_rust(&self) -> GraphIrAtomUpdate {
        self.0.clone()
    }
}

/// An atom: element, isotope, charge, implicit hydrogens, lone pairs, unpaired
/// electrons, and atom-scope constraints.
#[pyclass(eq)]
#[derive(PartialEq)]
pub struct AtomAst(GraphIrAtomForm);

#[pymethods]
impl AtomAst {
    /// Construct from an element — a single `Element` or an `ElementForm` expression —
    /// optionally setting fields.
    #[new]
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (element, *, isotope_mass=None, charge=None, implicit_hydrogens=None, lone_pairs=None, unpaired_electrons=None, constraints=None))]
    fn new(
        py: Python<'_>,
        element: ElementLike,
        isotope_mass: Option<IsotopeMassLike>,
        charge: Option<NumLike>,
        implicit_hydrogens: Option<NumLike>,
        lone_pairs: Option<NumLike>,
        unpaired_electrons: Option<PyRef<'_, UnpairedElectronsForm>>,
        constraints: Option<Py<AtomConstraintsAst>>,
    ) -> Self {
        let atom = GraphIrAtomForm::new(element.to_rust(py));
        AtomAst(apply_fields(
            atom,
            py,
            isotope_mass,
            charge,
            implicit_hydrogens,
            lone_pairs,
            unpaired_electrons,
            constraints,
        ))
    }

    /// Parse an atom-DSL string (e.g. `"C#c-1#v4"`) into an `AtomAst`.
    #[staticmethod]
    fn parse(s: &str) -> PyResult<Self> {
        GraphIrAtomForm::from_str(s).map(Self).map_err(parse_error)
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("AtomAst.parse('{}')", self.0)
    }

    #[getter]
    fn element(&self) -> ElementForm {
        ElementForm::from_rust(&self.0.element)
    }

    #[setter]
    fn set_element(&mut self, py: Python<'_>, value: ElementLike) {
        self.0.element = value.to_rust(py);
    }

    #[getter]
    fn isotope_mass(&self) -> IsotopeMassForm {
        IsotopeMassForm::from_rust(&self.0.isotope_mass)
    }

    #[setter]
    fn set_isotope_mass(&mut self, py: Python<'_>, value: IsotopeMassLike) {
        self.0.isotope_mass = value.to_rust(py);
    }

    #[getter]
    fn charge(&self, py: Python<'_>) -> PyResult<NumForm> {
        NumForm::from_rust(py, &self.0.charge)
    }

    #[setter]
    fn set_charge(&mut self, py: Python<'_>, value: NumLike) {
        self.0.charge = value.to_rust(py);
    }

    #[getter]
    fn implicit_hydrogens(&self, py: Python<'_>) -> PyResult<NumForm> {
        NumForm::from_rust(py, &self.0.implicit_hydrogens)
    }

    #[setter]
    fn set_implicit_hydrogens(&mut self, py: Python<'_>, value: NumLike) {
        self.0.implicit_hydrogens = value.to_rust(py);
    }

    #[getter]
    fn lone_pairs(&self, py: Python<'_>) -> PyResult<NumForm> {
        NumForm::from_rust(py, &self.0.lone_pairs)
    }

    #[setter]
    fn set_lone_pairs(&mut self, py: Python<'_>, value: NumLike) {
        self.0.lone_pairs = value.to_rust(py);
    }

    #[getter]
    fn unpaired_electrons(&self, py: Python<'_>) -> PyResult<UnpairedElectronsForm> {
        UnpairedElectronsForm::from_rust(py, &self.0.unpaired_electrons)
    }

    #[setter]
    fn set_unpaired_electrons(&mut self, py: Python<'_>, value: PyRef<'_, UnpairedElectronsForm>) {
        self.0.unpaired_electrons = value.to_rust(py);
    }

    /// The atom's constraints as a live handle onto this atom: reads borrow the
    /// current state, mutators write through to the atom in place.
    #[getter]
    fn constraints(slf: Py<Self>) -> AtomConstraintsView {
        AtomConstraintsView {
            backing: AtomConstraintsBacking::Atom(slf),
        }
    }

    /// Replace the whole constraint set (wipe-and-set) from a value container or
    /// a live view. Takes `slf` by handle and snapshots `value` *before* the write
    /// borrow, so `atom.constraints = atom.constraints` (a view over the same atom) reads
    /// while the atom is unborrowed instead of a double-borrow panic.
    #[setter]
    fn set_constraints(slf: Py<Self>, py: Python<'_>, value: AtomConstraintsLike) -> PyResult<()> {
        let snapshot = value.to_rust(py)?;
        slf.borrow_mut(py).0.constraints = snapshot;
        Ok(())
    }

    /// The fields as a dict keyed by field name; values are Python objects.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("element", self.element())?;
        dict.set_item("isotope_mass", self.isotope_mass())?;
        dict.set_item("charge", self.charge(py)?)?;
        dict.set_item("implicit_hydrogens", self.implicit_hydrogens(py)?)?;
        dict.set_item("lone_pairs", self.lone_pairs(py)?)?;
        dict.set_item("unpaired_electrons", self.unpaired_electrons(py)?)?;
        dict.set_item(
            "constraints",
            atom_constraints_asdict(py, &self.0.constraints)?,
        )?;
        Ok(dict)
    }
}

impl_py_lattice!(
    AtomAst,
    GraphIrAtomForm,
    |value: &AtomAst, _py: Python<'_>| -> PyResult<GraphIrAtomForm> { Ok(value.inner().clone()) },
    |_py: Python<'_>, value: GraphIrAtomForm| -> PyResult<AtomAst> {
        Ok(AtomAst::from_inner(value))
    }
);

/// A binding argument that converts a literal or Python value to its Rust value — the `*Like`
/// convention for these inputs (`*Input` is reserved for the DSL side). Extracted as
/// a PyO3 `FromPyObject` union tried in order; variants are `Ast` = the `*Ast`
/// wrapper, `Lit` = the literal, corresponding to `impl Into<..>` on the Rust builders.
///
/// `ElementLike` accepts a concrete `Element` or an `ElementForm`.
#[derive(FromPyObject)]
enum ElementLike {
    Ast(Py<ElementForm>),
    Lit(Element),
}

impl ElementLike {
    fn to_rust(&self, py: Python<'_>) -> GraphIrElementForm {
        match self {
            ElementLike::Ast(expr) => expr.bind(py).borrow().to_rust(),
            ElementLike::Lit(element) => GraphIrElementForm::Lit(ChemElement::from(element)),
        }
    }
}

/// An `IsotopeMassForm` or a Python `int` (→ `IsotopeMassForm::Lit`, a mass number).
#[derive(FromPyObject)]
enum IsotopeMassLike {
    Ast(Py<IsotopeMassForm>),
    Lit(u32),
}

impl IsotopeMassLike {
    fn to_rust(&self, py: Python<'_>) -> GraphIrIsotopeMassForm {
        match self {
            IsotopeMassLike::Ast(mass) => mass.bind(py).borrow().to_rust(),
            IsotopeMassLike::Lit(number) => GraphIrIsotopeMassForm::Lit(*number),
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn apply_fields(
    mut atom: GraphIrAtomForm,
    py: Python<'_>,
    isotope_mass: Option<IsotopeMassLike>,
    charge: Option<NumLike>,
    implicit_hydrogens: Option<NumLike>,
    lone_pairs: Option<NumLike>,
    unpaired_electrons: Option<PyRef<'_, UnpairedElectronsForm>>,
    constraints: Option<Py<AtomConstraintsAst>>,
) -> GraphIrAtomForm {
    if let Some(isotope_mass) = isotope_mass {
        atom = atom.with_isotope_mass(isotope_mass.to_rust(py));
    }
    if let Some(charge) = charge {
        atom = atom.with_charge(charge.to_rust(py));
    }
    if let Some(implicit_hydrogens) = implicit_hydrogens {
        atom = atom.with_implicit_hydrogens(implicit_hydrogens.to_rust(py));
    }
    if let Some(lone_pairs) = lone_pairs {
        atom = atom.with_lone_pairs(lone_pairs.to_rust(py));
    }
    if let Some(unpaired_electrons) = unpaired_electrons {
        atom = atom.with_unpaired_electrons(unpaired_electrons.to_rust(py));
    }
    if let Some(constraints) = constraints {
        atom.constraints = constraints.bind(py).borrow().inner().clone();
    }
    atom
}

impl AtomAst {
    /// The wrapped AST atom — read access for molecule construction.
    pub(crate) fn inner(&self) -> &GraphIrAtomForm {
        &self.0
    }

    /// Mutable access to the wrapped AST atom — write access for the atom-backed
    /// constraints view.
    pub(crate) fn inner_mut(&mut self) -> &mut GraphIrAtomForm {
        &mut self.0
    }

    /// Wrap an AST atom (the hold-the-value `from_inner` bridge, paired with
    /// `inner`).
    pub(crate) fn from_inner(atom: GraphIrAtomForm) -> Self {
        AtomAst(atom)
    }
}

/// A view of one atom within a molecule: a handle to the molecule plus the atom's
/// index. Field reads rebuild the transient Rust view; the molecule is never copied.
#[pyclass]
pub struct AtomView {
    owner: Py<MoleculeAst>,
    id: GraphIrAtomId,
}

impl AtomView {
    fn atom<'a>(&self, molecule: &'a GraphIrMolecule) -> PyResult<&'a GraphIrAtomForm> {
        molecule
            .atoms()
            .get(self.id)
            .map(|view| view.ast)
            .ok_or_else(|| PyIndexError::new_err("atom id out of range"))
    }
}

#[pymethods]
impl AtomView {
    #[getter]
    fn id(&self) -> u32 {
        self.id.0
    }

    fn __repr__(&self) -> String {
        format!("AtomView(id={})", self.id.0)
    }

    #[getter]
    fn element(&self, py: Python<'_>) -> PyResult<ElementForm> {
        let molecule = self.owner.bind(py).borrow();
        Ok(ElementForm::from_rust(
            &self.atom(molecule.inner())?.element,
        ))
    }

    #[setter]
    fn set_element(&self, py: Python<'_>, value: ElementLike) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .atom_mut(self.id)
            .ast
            .element = value.to_rust(py);
    }

    #[getter]
    fn isotope_mass(&self, py: Python<'_>) -> PyResult<IsotopeMassForm> {
        let molecule = self.owner.bind(py).borrow();
        Ok(IsotopeMassForm::from_rust(
            &self.atom(molecule.inner())?.isotope_mass,
        ))
    }

    #[setter]
    fn set_isotope_mass(&self, py: Python<'_>, value: IsotopeMassLike) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .atom_mut(self.id)
            .ast
            .isotope_mass = value.to_rust(py);
    }

    #[getter]
    fn charge(&self, py: Python<'_>) -> PyResult<NumForm> {
        let molecule = self.owner.bind(py).borrow();
        NumForm::from_rust(py, &self.atom(molecule.inner())?.charge)
    }

    #[setter]
    fn set_charge(&self, py: Python<'_>, value: NumLike) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .atom_mut(self.id)
            .ast
            .charge = value.to_rust(py);
    }

    #[getter]
    fn implicit_hydrogens(&self, py: Python<'_>) -> PyResult<NumForm> {
        let molecule = self.owner.bind(py).borrow();
        NumForm::from_rust(py, &self.atom(molecule.inner())?.implicit_hydrogens)
    }

    #[setter]
    fn set_implicit_hydrogens(&self, py: Python<'_>, value: NumLike) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .atom_mut(self.id)
            .ast
            .implicit_hydrogens = value.to_rust(py);
    }

    #[getter]
    fn lone_pairs(&self, py: Python<'_>) -> PyResult<NumForm> {
        let molecule = self.owner.bind(py).borrow();
        NumForm::from_rust(py, &self.atom(molecule.inner())?.lone_pairs)
    }

    #[setter]
    fn set_lone_pairs(&self, py: Python<'_>, value: NumLike) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .atom_mut(self.id)
            .ast
            .lone_pairs = value.to_rust(py);
    }

    #[getter]
    fn unpaired_electrons(&self, py: Python<'_>) -> PyResult<UnpairedElectronsForm> {
        let molecule = self.owner.bind(py).borrow();
        UnpairedElectronsForm::from_rust(py, &self.atom(molecule.inner())?.unpaired_electrons)
    }

    #[setter]
    fn set_unpaired_electrons(&self, py: Python<'_>, value: PyRef<'_, UnpairedElectronsForm>) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .atom_mut(self.id)
            .ast
            .unpaired_electrons = value.to_rust(py);
    }

    /// The atom's constraints as a live handle onto the molecule: reads borrow the
    /// current state, mutators write through to the atom in place.
    #[getter]
    fn constraints(&self, py: Python<'_>) -> AtomConstraintsView {
        AtomConstraintsView {
            backing: AtomConstraintsBacking::Molecule {
                owner: self.owner.clone_ref(py),
                id: self.id,
            },
        }
    }

    /// Replace the whole constraint set of the backing atom in place (wipe-and-set)
    /// from a value container or a live view.
    #[setter]
    fn set_constraints(&self, py: Python<'_>, value: AtomConstraintsLike) -> PyResult<()> {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .atom_mut(self.id)
            .ast
            .constraints = value.to_rust(py)?;
        Ok(())
    }

    /// The fields as a dict keyed by field name; values are Python objects —
    /// symmetric with `AtomAst.asdict`, read through the view.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let molecule = self.owner.bind(py).borrow();
        let atom = self.atom(molecule.inner())?;
        let dict = PyDict::new(py);
        dict.set_item("element", ElementForm::from_rust(&atom.element))?;
        dict.set_item(
            "isotope_mass",
            IsotopeMassForm::from_rust(&atom.isotope_mass),
        )?;
        dict.set_item("charge", NumForm::from_rust(py, &atom.charge)?)?;
        dict.set_item(
            "implicit_hydrogens",
            NumForm::from_rust(py, &atom.implicit_hydrogens)?,
        )?;
        dict.set_item("lone_pairs", NumForm::from_rust(py, &atom.lone_pairs)?)?;
        dict.set_item(
            "unpaired_electrons",
            UnpairedElectronsForm::from_rust(py, &atom.unpaired_electrons)?,
        )?;
        dict.set_item(
            "constraints",
            atom_constraints_asdict(py, &atom.constraints)?,
        )?;
        Ok(dict)
    }
}

/// Resolve a possibly-negative Python index (negative counts from the end) into an
/// existing atom id, or `IndexError`.
fn resolve_atom_index(molecule: &GraphIrMolecule, index: isize) -> PyResult<GraphIrAtomId> {
    let count = molecule.atoms().count();
    let resolved = if index < 0 {
        index + count as isize
    } else {
        index
    };
    if resolved < 0 {
        return Err(PyIndexError::new_err("atom id out of range"));
    }
    let id = GraphIrAtomId(resolved as u32);
    if molecule.atoms().contains(id) {
        Ok(id)
    } else {
        Err(PyIndexError::new_err("atom id out of range"))
    }
}

/// The atoms of a molecule, indexed by integer position.
#[pyclass]
pub struct AtomViews {
    owner: Py<MoleculeAst>,
}

#[pymethods]
impl AtomViews {
    fn __len__(&self, py: Python<'_>) -> usize {
        self.owner.bind(py).borrow().inner().atoms().count()
    }

    fn __repr__(&self, py: Python<'_>) -> String {
        format!(
            "AtomViews(len={})",
            self.owner.bind(py).borrow().inner().atoms().count()
        )
    }

    fn __getitem__(&self, py: Python<'_>, index: isize) -> PyResult<AtomView> {
        let molecule = self.owner.bind(py).borrow();
        let id = resolve_atom_index(molecule.inner(), index)?;
        Ok(AtomView {
            owner: self.owner.clone_ref(py),
            id,
        })
    }

    /// Replace the whole atom at `index` in place.
    fn __setitem__(&self, py: Python<'_>, index: isize, atom: PyRef<'_, AtomAst>) -> PyResult<()> {
        let mut molecule = self.owner.borrow_mut(py);
        let id = resolve_atom_index(molecule.inner(), index)?;
        *molecule.inner_mut().atom_mut(id).ast = atom.inner().clone();
        Ok(())
    }

    fn __iter__(&self, py: Python<'_>) -> AtomViewIter {
        let ids = self
            .owner
            .bind(py)
            .borrow()
            .inner()
            .atoms()
            .ids()
            .collect::<Vec<_>>();
        AtomViewIter {
            owner: self.owner.clone_ref(py),
            ids: ids.into_iter(),
        }
    }
}

impl AtomViews {
    /// Build the atom-views handle for `owner` (the `.atoms` accessor on the molecule).
    pub(crate) fn new(owner: Py<MoleculeAst>) -> AtomViews {
        AtomViews { owner }
    }
}

#[pyclass]
struct AtomViewIter {
    owner: Py<MoleculeAst>,
    ids: IntoIter<GraphIrAtomId>,
}

#[pymethods]
impl AtomViewIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> Option<AtomView> {
        self.ids.next().map(|id| AtomView {
            owner: self.owner.clone_ref(py),
            id,
        })
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_ir::ir::{
        AromaticValenceForm as GraphIrAromaticValenceForm,
        AtomConstraintForm as GraphIrAtomConstraintForm,
        AtomConstraintKey as GraphIrAtomConstraintKey,
        AtomConstraintsForm as GraphIrAtomConstraintsForm, MemOp as GraphIrMemOp,
        MoleculeEntries as GraphIrMoleculeEntries, NumForm as GraphIrNumForm,
        RingScope as GraphIrRingScope, StereoCoset as GraphIrStereoCoset,
        TetrahedralStereoForm as GraphIrTetrahedralStereoForm,
        UnpairedElectronsForm as GraphIrUnpairedElectronsForm,
    };

    use super::*;
    use crate::constraint::atom::{
        AromaticValenceAst, AromaticValenceLike, AtomConstraintAst, AtomConstraintKey,
        AtomConstraintsUpdate, TetrahedralStereoLike,
    };
    use crate::convert::into_py_variant;
    use crate::stereo::{TetrahedralConfiguration, TetrahedralStereoForm};

    #[rstest]
    #[case(GraphIrElementForm::Undetermined)]
    #[case(GraphIrElementForm::Lit(ChemElement::C))]
    #[case(GraphIrElementForm::LitSet(Box::new(BTreeSet::from([ChemElement::C, ChemElement::N]))))]
    #[case(GraphIrElementForm::NotSet(Box::new(BTreeSet::from([ChemElement::O]))))]
    #[case(GraphIrElementForm::Var(Box::new(("x".to_string(), None))))]
    #[case(GraphIrElementForm::Var(Box::new((
        "y".to_string(),
        Some((GraphIrMemOp::In, BTreeSet::from([ChemElement::C, ChemElement::N]))),
    ))))]
    fn test_element_ast_roundtrip(#[case] ast: GraphIrElementForm) {
        assert_eq!(ElementForm::from_rust(&ast).to_rust(), ast);
    }

    #[rstest]
    #[case(GraphIrElementForm::Lit(ChemElement::C), Some(ChemElement::C))]
    #[case(GraphIrElementForm::Undetermined, None)]
    #[case(GraphIrElementForm::LitSet(Box::new(BTreeSet::from([ChemElement::C, ChemElement::N]))), None)]
    fn test_element_ast_as_lit(
        #[case] ast: GraphIrElementForm,
        #[case] expected: Option<ChemElement>,
    ) {
        let got = ElementForm::from_rust(&ast)
            .as_lit()
            .map(|e| ChemElement::from(&e));
        assert_eq!(got, expected);
    }

    #[rstest]
    #[case(GraphIrIsotopeMassForm::Undetermined)]
    #[case(GraphIrIsotopeMassForm::Natural)]
    #[case(GraphIrIsotopeMassForm::Lit(13))]
    #[case(GraphIrIsotopeMassForm::LitSet(Box::new(BTreeSet::from([12, 13, 14]))))]
    #[case(GraphIrIsotopeMassForm::Var(Box::new(("x".to_string(), None))))]
    #[case(GraphIrIsotopeMassForm::Var(Box::new((
        "y".to_string(),
        Some(BTreeSet::from([12, 13])),
    ))))]
    fn test_isotope_mass_ast_roundtrip(#[case] ast: GraphIrIsotopeMassForm) {
        assert_eq!(IsotopeMassForm::from_rust(&ast).to_rust(), ast);
    }

    #[rstest]
    #[case(
        GraphIrIsotopeMassForm::Lit(13),
        Some(GraphIrIsotopeMass::MassNumber(13))
    )]
    #[case(GraphIrIsotopeMassForm::Natural, Some(GraphIrIsotopeMass::Natural))]
    #[case(GraphIrIsotopeMassForm::Undetermined, None)]
    fn test_isotope_mass_ast_as_lit(
        #[case] ast: GraphIrIsotopeMassForm,
        #[case] expected: Option<GraphIrIsotopeMass>,
    ) {
        assert_eq!(
            IsotopeMassForm::from_rust(&ast)
                .as_lit()
                .map(IsotopeMass::to_rust),
            expected
        );
    }

    fn carbon_oxygen(py: Python<'_>) -> Py<MoleculeAst> {
        let molecule = GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
            atoms: vec![
                GraphIrAtomForm::from_element(ChemElement::C),
                GraphIrAtomForm::from_element(ChemElement::O),
            ],
            ..Default::default()
        });
        Py::new(py, MoleculeAst::from_rust(molecule)).unwrap()
    }

    #[rstest]
    fn test_atom_view_element() {
        Python::attach(|py| {
            let view = AtomView {
                owner: carbon_oxygen(py),
                id: GraphIrAtomId(1),
            };
            assert_eq!(view.id(), 1);
            match view.element(py).unwrap() {
                ElementForm::Lit(e) => assert_eq!(ChemElement::from(&e), ChemElement::O),
                _ => panic!("expected Lit"),
            }
        });
    }

    #[rstest]
    fn test_atom_view_set_charge() {
        Python::attach(|py| {
            let owner = carbon_oxygen(py);
            let view = AtomView {
                owner: owner.clone_ref(py),
                id: GraphIrAtomId(0),
            };
            view.set_charge(py, NumLike::Lit(-1));
            let fresh = AtomView {
                owner,
                id: GraphIrAtomId(0),
            };
            match fresh.charge(py).unwrap() {
                NumForm::Lit(n) => assert_eq!(n, -1),
                _ => panic!("expected Lit"),
            }
        });
    }

    #[rstest]
    fn test_atom_view_set_element() {
        Python::attach(|py| {
            let owner = carbon_oxygen(py);
            let view = AtomView {
                owner: owner.clone_ref(py),
                id: GraphIrAtomId(0),
            };
            view.set_element(py, ElementLike::Lit(Element::from(ChemElement::N)));
            let fresh = AtomView {
                owner,
                id: GraphIrAtomId(0),
            };
            match fresh.element(py).unwrap() {
                ElementForm::Lit(e) => assert_eq!(ChemElement::from(&e), ChemElement::N),
                _ => panic!("expected Lit"),
            }
        });
    }

    #[rstest]
    fn test_atom_view_set_isotope_mass() {
        Python::attach(|py| {
            let owner = carbon_oxygen(py);
            let view = AtomView {
                owner: owner.clone_ref(py),
                id: GraphIrAtomId(0),
            };
            view.set_isotope_mass(py, IsotopeMassLike::Lit(13));
            let fresh = AtomView {
                owner,
                id: GraphIrAtomId(0),
            };
            match fresh.isotope_mass(py).unwrap() {
                IsotopeMassForm::Lit(mass) => assert_eq!(mass, 13),
                _ => panic!("expected Lit"),
            }
        });
    }

    #[rstest]
    fn test_atom_view_set_unpaired_electrons() {
        Python::attach(|py| {
            let owner = carbon_oxygen(py);
            let view = AtomView {
                owner: owner.clone_ref(py),
                id: GraphIrAtomId(0),
            };
            let unpaired_electrons = Py::new(
                py,
                UnpairedElectronsForm::from_rust(
                    py,
                    &GraphIrUnpairedElectronsForm {
                        count: GraphIrNumForm::Lit(1),
                        multiplicity: GraphIrNumForm::Lit(2),
                    },
                )
                .unwrap(),
            )
            .unwrap();
            view.set_unpaired_electrons(py, unpaired_electrons.bind(py).borrow());
            let fresh = AtomView {
                owner,
                id: GraphIrAtomId(0),
            };
            assert_eq!(
                fresh.unpaired_electrons(py).unwrap().to_rust(py),
                GraphIrUnpairedElectronsForm {
                    count: GraphIrNumForm::Lit(1),
                    multiplicity: GraphIrNumForm::Lit(2),
                }
            );
        });
    }

    #[rstest]
    fn test_atom_views_len_and_getitem() {
        Python::attach(|py| {
            let views = AtomViews {
                owner: carbon_oxygen(py),
            };
            assert_eq!(views.__len__(py), 2);
            assert_eq!(views.__getitem__(py, 0).unwrap().id(), 0);
            assert_eq!(views.__getitem__(py, -1).unwrap().id(), 1);
            assert_eq!(views.__getitem__(py, -2).unwrap().id(), 0);
            assert!(views.__getitem__(py, 5).is_err());
            assert!(views.__getitem__(py, -3).is_err());
        });
    }

    #[rstest]
    fn test_atom_views_setitem() {
        Python::attach(|py| {
            let views = AtomViews {
                owner: carbon_oxygen(py),
            };
            let nitrogen = Py::new(
                py,
                AtomAst::from_inner(GraphIrAtomForm::from_element(ChemElement::N)),
            )
            .unwrap();
            views
                .__setitem__(py, 0, nitrogen.bind(py).borrow())
                .unwrap();
            match views.__getitem__(py, 0).unwrap().element(py).unwrap() {
                ElementForm::Lit(e) => assert_eq!(ChemElement::from(&e), ChemElement::N),
                _ => panic!("expected Lit"),
            }
        });
    }

    #[rstest]
    fn test_atom_views_setitem_error() {
        Python::attach(|py| {
            let views = AtomViews {
                owner: carbon_oxygen(py),
            };
            let nitrogen = Py::new(
                py,
                AtomAst::from_inner(GraphIrAtomForm::from_element(ChemElement::N)),
            )
            .unwrap();
            assert!(views
                .__setitem__(py, 5, nitrogen.bind(py).borrow())
                .is_err());
        });
    }

    #[rstest]
    fn test_atom_ast_constraints() {
        let atom = AtomAst(
            GraphIrAtomForm::from_element(ChemElement::C)
                .with_constraint(GraphIrAtomConstraintForm::valence(4)),
        );
        assert_eq!(atom.inner().constraints.len(), 1);
    }

    #[rstest]
    #[case::bare("C")]
    #[case::charge("N#c+")]
    #[case::valence("C#v4")]
    #[case::lone_pairs("O#n2")]
    #[case::ring_size("C#R(6)")]
    fn test_atom_ast_parse(#[case] dsl: &str) {
        let atom = AtomAst::parse(dsl).unwrap();
        assert_eq!(atom.__str__(), dsl);
        assert_eq!(atom.__repr__(), format!("AtomAst.parse('{dsl}')"));
    }

    #[rstest]
    fn test_atom_ast_parse_error() {
        assert!(AtomAst::parse("Zz##").is_err());
    }

    #[rstest]
    fn test_atom_ast_set_constraints_from_view() {
        Python::attach(|py| {
            let src = Py::new(
                py,
                AtomAst::from_inner(
                    GraphIrAtomForm::from_element(ChemElement::C)
                        .with_constraint(GraphIrAtomConstraintForm::valence(4)),
                ),
            )
            .unwrap();
            let view = Py::new(
                py,
                AtomConstraintsView {
                    backing: AtomConstraintsBacking::Atom(src),
                },
            )
            .unwrap();
            let dst = Py::new(
                py,
                AtomAst::from_inner(GraphIrAtomForm::from_element(ChemElement::N)),
            )
            .unwrap();
            AtomAst::set_constraints(dst.clone_ref(py), py, AtomConstraintsLike::View(view))
                .unwrap();
            assert_eq!(
                dst.bind(py)
                    .borrow()
                    .inner()
                    .constraints
                    .valence()
                    .unwrap()
                    .clone(),
                GraphIrNumForm::Lit(4)
            );
        });
    }

    /// Regression: assigning an atom's own constraints view back to it snapshots before
    /// the write borrow, so it is a no-op, not a RefCell double-borrow panic
    /// (`atom.constraints = atom.constraints`).
    #[rstest]
    fn test_atom_ast_set_constraints_self() {
        Python::attach(|py| {
            let atom = Py::new(
                py,
                AtomAst::from_inner(
                    GraphIrAtomForm::from_element(ChemElement::C)
                        .with_constraint(GraphIrAtomConstraintForm::valence(4)),
                ),
            )
            .unwrap();
            let own_view = Py::new(
                py,
                AtomConstraintsView {
                    backing: AtomConstraintsBacking::Atom(atom.clone_ref(py)),
                },
            )
            .unwrap();
            AtomAst::set_constraints(atom.clone_ref(py), py, AtomConstraintsLike::View(own_view))
                .unwrap();
            assert_eq!(
                atom.bind(py)
                    .borrow()
                    .inner()
                    .constraints
                    .valence()
                    .unwrap()
                    .clone(),
                GraphIrNumForm::Lit(4)
            );
        });
    }

    #[rstest]
    fn test_atom_view_constraints() {
        Python::attach(|py| {
            let view = AtomView {
                owner: carbon_oxygen(py),
                id: GraphIrAtomId(1),
            };
            match view.constraints(py).backing {
                AtomConstraintsBacking::Molecule { id, .. } => assert_eq!(id, GraphIrAtomId(1)),
                _ => panic!("expected molecule-backed view"),
            }
        });
    }

    #[rstest]
    #[case(GraphIrAtomConstraintForm::valence(4))]
    #[case(GraphIrAtomConstraintForm::aromatic_valence(GraphIrAromaticValenceForm::aromatic(1)))]
    #[case(GraphIrAtomConstraintForm::ring_membership(GraphIrRingScope::All, 2))]
    #[case(GraphIrAtomConstraintForm::tetrahedral_stereo(
        GraphIrTetrahedralStereoForm::not_stereo()
    ))]
    fn test_atom_constraint_roundtrip(#[case] ast: GraphIrAtomConstraintForm) {
        Python::attach(|py| {
            assert_eq!(
                AtomConstraintAst::from_rust(py, &ast).unwrap().to_rust(py),
                ast
            );
        });
    }

    #[rstest]
    fn test_atom_constraints_len_contains() {
        Python::attach(|py| {
            let valence = into_py_variant(
                py,
                AtomConstraintAst::from_rust(py, &GraphIrAtomConstraintForm::valence(4)).unwrap(),
            )
            .unwrap();
            let degree = into_py_variant(
                py,
                AtomConstraintAst::from_rust(py, &GraphIrAtomConstraintForm::degree(3)).unwrap(),
            )
            .unwrap();
            let constraints = AtomConstraintsAst::new(py, vec![valence, degree]);
            assert_eq!(constraints.__len__(), 2);
            assert!(constraints.__contains__(
                py,
                into_py_variant(py, AtomConstraintKey::Valence()).unwrap()
            ));
            assert!(constraints.__contains__(
                py,
                into_py_variant(py, AtomConstraintKey::Degree()).unwrap()
            ));
            assert!(!constraints.__contains__(
                py,
                into_py_variant(py, AtomConstraintKey::TotalHydrogens()).unwrap()
            ));
        });
    }

    #[rstest]
    fn test_atom_constraints_keys_values_items() {
        Python::attach(|py| {
            let valence = into_py_variant(
                py,
                AtomConstraintAst::from_rust(py, &GraphIrAtomConstraintForm::valence(4)).unwrap(),
            )
            .unwrap();
            let degree = into_py_variant(
                py,
                AtomConstraintAst::from_rust(py, &GraphIrAtomConstraintForm::degree(3)).unwrap(),
            )
            .unwrap();
            let constraints = AtomConstraintsAst::new(py, vec![valence, degree]);

            let mut keys = constraints.__iter__(py).unwrap();
            assert_eq!(
                keys.__next__().unwrap().bind(py).borrow().to_rust(py),
                GraphIrAtomConstraintKey::Valence
            );
            assert_eq!(
                keys.__next__().unwrap().bind(py).borrow().to_rust(py),
                GraphIrAtomConstraintKey::Degree
            );
            assert!(keys.__next__().is_none());

            let mut values = constraints.values(py).unwrap();
            assert_eq!(
                values.__next__().unwrap().bind(py).borrow().to_rust(py),
                GraphIrAtomConstraintForm::valence(4)
            );

            let mut items = constraints.items(py).unwrap();
            let (key, value) = items.__next__().unwrap();
            assert_eq!(
                key.bind(py).borrow().to_rust(py),
                GraphIrAtomConstraintKey::Valence
            );
            assert_eq!(
                value.bind(py).borrow().to_rust(py),
                GraphIrAtomConstraintForm::valence(4)
            );
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_get() {
        Python::attach(|py| {
            let valence = into_py_variant(
                py,
                AtomConstraintAst::from_rust(py, &GraphIrAtomConstraintForm::valence(4)).unwrap(),
            )
            .unwrap();
            let constraints = AtomConstraintsAst::new(py, vec![valence]);
            let present = constraints
                .get(
                    py,
                    into_py_variant(py, AtomConstraintKey::Valence()).unwrap(),
                    None,
                )
                .unwrap();
            let expected = into_py_variant(
                py,
                AtomConstraintAst::from_rust(py, &GraphIrAtomConstraintForm::valence(4)).unwrap(),
            )
            .unwrap()
            .into_any();
            assert!(present.bind(py).eq(expected.bind(py)).unwrap());
            let absent = constraints
                .get(
                    py,
                    into_py_variant(py, AtomConstraintKey::Degree()).unwrap(),
                    None,
                )
                .unwrap();
            assert!(absent.bind(py).is_none());
            let sentinel = into_py_variant(py, AtomConstraintKey::Degree())
                .unwrap()
                .into_any();
            let defaulted = constraints
                .get(
                    py,
                    into_py_variant(py, AtomConstraintKey::Degree()).unwrap(),
                    Some(sentinel.clone_ref(py)),
                )
                .unwrap();
            assert_eq!(defaulted.as_ptr(), sentinel.as_ptr());
        });
    }

    #[rstest]
    fn test_atom_constraints_valence() {
        Python::attach(|py| {
            let valence = into_py_variant(
                py,
                AtomConstraintAst::from_rust(py, &GraphIrAtomConstraintForm::valence(4)).unwrap(),
            )
            .unwrap();
            let degree = into_py_variant(
                py,
                AtomConstraintAst::from_rust(py, &GraphIrAtomConstraintForm::degree(3)).unwrap(),
            )
            .unwrap();
            let constraints = AtomConstraintsAst::new(py, vec![valence, degree]);
            assert_eq!(
                constraints.valence(py).unwrap().unwrap().to_rust(py),
                GraphIrNumForm::Lit(4)
            );
            assert_eq!(
                constraints.degree(py).unwrap().unwrap().to_rust(py),
                GraphIrNumForm::Lit(3)
            );
            assert!(constraints.total_valence(py).unwrap().is_none());
            assert!(constraints.aromatic_valence(py).unwrap().is_none());
        });
    }

    #[rstest]
    fn test_atom_constraints_ring_size_count() {
        Python::attach(|py| {
            let membership = into_py_variant(
                py,
                AtomConstraintAst::from_rust(
                    py,
                    &GraphIrAtomConstraintForm::ring_membership(GraphIrRingScope::Size(6), 1),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = Py::new(py, AtomConstraintsAst::new(py, vec![membership])).unwrap();
            let proxy = AtomConstraintsAst::ring_size_count(constraints.clone_ref(py));
            assert_eq!(
                proxy.__getitem__(py, 6).unwrap().unwrap().to_rust(py),
                GraphIrNumForm::Lit(1)
            );
            assert!(proxy.__getitem__(py, 5).unwrap().is_none());
            assert!(constraints
                .bind(py)
                .borrow()
                .ring_count(py)
                .unwrap()
                .is_none());
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_set() {
        Python::attach(|py| {
            let mut constraints = AtomConstraintsAst::new(py, vec![]);
            let valence = into_py_variant(
                py,
                AtomConstraintAst::from_rust(py, &GraphIrAtomConstraintForm::valence(4)).unwrap(),
            )
            .unwrap();
            constraints.set(py, valence);
            assert_eq!(constraints.__len__(), 1);
            assert_eq!(
                constraints.valence(py).unwrap().unwrap().to_rust(py),
                GraphIrNumForm::Lit(4)
            );
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_pop() {
        Python::attach(|py| {
            let valence = into_py_variant(
                py,
                AtomConstraintAst::from_rust(py, &GraphIrAtomConstraintForm::valence(4)).unwrap(),
            )
            .unwrap();
            let mut constraints = AtomConstraintsAst::new(py, vec![valence]);
            let removed = constraints
                .pop(
                    py,
                    into_py_variant(py, AtomConstraintKey::Valence()).unwrap(),
                )
                .unwrap();
            match removed {
                Some(AtomConstraintAst::Valence(v)) => {
                    assert_eq!(v.bind(py).borrow().to_rust(py), GraphIrNumForm::Lit(4))
                }
                _ => panic!("expected removed Valence(Lit(4))"),
            }
            assert_eq!(constraints.__len__(), 0);
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_update() {
        Python::attach(|py| {
            let constraints = Py::new(py, AtomConstraintsAst::new(py, vec![])).unwrap();
            let mut other = GraphIrAtomConstraintsForm::new();
            other.set(GraphIrAtomConstraintForm::valence(4));
            other.set(GraphIrAtomConstraintForm::degree(3));
            AtomConstraintsAst::update(
                constraints.clone_ref(py),
                py,
                AtomConstraintsUpdate::Container(
                    Py::new(py, AtomConstraintsAst::from_inner(other)).unwrap(),
                ),
            )
            .unwrap();
            let c = constraints.bind(py).borrow();
            assert_eq!(c.__len__(), 2);
            assert_eq!(
                c.valence(py).unwrap().unwrap().to_rust(py),
                GraphIrNumForm::Lit(4)
            );
            assert_eq!(
                c.degree(py).unwrap().unwrap().to_rust(py),
                GraphIrNumForm::Lit(3)
            );
        });
    }

    #[rstest]
    fn test_atom_constraints_view_set() {
        Python::attach(|py| {
            let owner = Py::new(
                py,
                MoleculeAst::from_rust(GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                    ..Default::default()
                })),
            )
            .unwrap();
            let view = AtomConstraintsView {
                backing: AtomConstraintsBacking::Molecule {
                    owner: owner.clone_ref(py),
                    id: GraphIrAtomId(0),
                },
            };
            let valence = into_py_variant(
                py,
                AtomConstraintAst::from_rust(py, &GraphIrAtomConstraintForm::valence(4)).unwrap(),
            )
            .unwrap();
            view.set(py, valence);
            let fresh = AtomConstraintsView {
                backing: AtomConstraintsBacking::Molecule {
                    owner,
                    id: GraphIrAtomId(0),
                },
            };
            assert_eq!(fresh.__len__(py).unwrap(), 1);
            match fresh
                .__getitem__(
                    py,
                    into_py_variant(py, AtomConstraintKey::Valence()).unwrap(),
                )
                .unwrap()
            {
                AtomConstraintAst::Valence(v) => {
                    assert_eq!(v.bind(py).borrow().to_rust(py), GraphIrNumForm::Lit(4))
                }
                _ => panic!("expected Valence(Lit(4))"),
            }
        });
    }

    #[rstest]
    fn test_atom_constraints_view_pop() {
        Python::attach(|py| {
            let atom = GraphIrAtomForm::from_element(ChemElement::C)
                .with_constraint(GraphIrAtomConstraintForm::valence(4));
            let owner = Py::new(
                py,
                MoleculeAst::from_rust(GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![atom],
                    ..Default::default()
                })),
            )
            .unwrap();
            let view = AtomConstraintsView {
                backing: AtomConstraintsBacking::Molecule {
                    owner: owner.clone_ref(py),
                    id: GraphIrAtomId(0),
                },
            };
            let removed = view
                .pop(
                    py,
                    into_py_variant(py, AtomConstraintKey::Valence()).unwrap(),
                )
                .unwrap();
            match removed {
                Some(AtomConstraintAst::Valence(v)) => {
                    assert_eq!(v.bind(py).borrow().to_rust(py), GraphIrNumForm::Lit(4))
                }
                _ => panic!("expected removed Valence(Lit(4))"),
            }
            let fresh = AtomConstraintsView {
                backing: AtomConstraintsBacking::Molecule {
                    owner,
                    id: GraphIrAtomId(0),
                },
            };
            assert_eq!(fresh.__len__(py).unwrap(), 0);
        });
    }

    #[rstest]
    fn test_atom_constraints_view_update() {
        Python::attach(|py| {
            let owner = Py::new(
                py,
                MoleculeAst::from_rust(GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                    ..Default::default()
                })),
            )
            .unwrap();
            let view = AtomConstraintsView {
                backing: AtomConstraintsBacking::Molecule {
                    owner: owner.clone_ref(py),
                    id: GraphIrAtomId(0),
                },
            };
            let mut other = GraphIrAtomConstraintsForm::new();
            other.set(GraphIrAtomConstraintForm::valence(4));
            other.set(GraphIrAtomConstraintForm::degree(3));
            view.update(
                py,
                AtomConstraintsUpdate::Container(
                    Py::new(py, AtomConstraintsAst::from_inner(other)).unwrap(),
                ),
            )
            .unwrap();
            let fresh = AtomConstraintsView {
                backing: AtomConstraintsBacking::Molecule {
                    owner,
                    id: GraphIrAtomId(0),
                },
            };
            assert_eq!(fresh.__len__(py).unwrap(), 2);
        });
    }

    #[rstest]
    fn test_atom_constraints_view_set_atom_backed() {
        Python::attach(|py| {
            let atom = Py::new(
                py,
                AtomAst::from_inner(GraphIrAtomForm::from_element(ChemElement::C)),
            )
            .unwrap();
            let view = AtomConstraintsView {
                backing: AtomConstraintsBacking::Atom(atom.clone_ref(py)),
            };
            let valence = into_py_variant(
                py,
                AtomConstraintAst::from_rust(py, &GraphIrAtomConstraintForm::valence(4)).unwrap(),
            )
            .unwrap();
            view.set(py, valence);
            // a fresh view proves the write hit the standalone atom, not a copy
            let fresh = AtomConstraintsView {
                backing: AtomConstraintsBacking::Atom(atom),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 1);
            match fresh
                .__getitem__(
                    py,
                    into_py_variant(py, AtomConstraintKey::Valence()).unwrap(),
                )
                .unwrap()
            {
                AtomConstraintAst::Valence(v) => {
                    assert_eq!(v.bind(py).borrow().to_rust(py), GraphIrNumForm::Lit(4))
                }
                _ => panic!("expected Valence(Lit(4))"),
            }
        });
    }

    #[rstest]
    fn test_atom_constraints_view_pop_atom_backed() {
        Python::attach(|py| {
            let atom = Py::new(
                py,
                AtomAst::from_inner(
                    GraphIrAtomForm::from_element(ChemElement::C)
                        .with_constraint(GraphIrAtomConstraintForm::valence(4)),
                ),
            )
            .unwrap();
            let view = AtomConstraintsView {
                backing: AtomConstraintsBacking::Atom(atom.clone_ref(py)),
            };
            let removed = view
                .pop(
                    py,
                    into_py_variant(py, AtomConstraintKey::Valence()).unwrap(),
                )
                .unwrap();
            match removed {
                Some(AtomConstraintAst::Valence(v)) => {
                    assert_eq!(v.bind(py).borrow().to_rust(py), GraphIrNumForm::Lit(4))
                }
                _ => panic!("expected removed Valence(Lit(4))"),
            }
            let fresh = AtomConstraintsView {
                backing: AtomConstraintsBacking::Atom(atom),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 0);
        });
    }

    #[rstest]
    fn test_atom_constraints_view_update_atom_backed() {
        Python::attach(|py| {
            let atom = Py::new(
                py,
                AtomAst::from_inner(GraphIrAtomForm::from_element(ChemElement::C)),
            )
            .unwrap();
            let view = AtomConstraintsView {
                backing: AtomConstraintsBacking::Atom(atom.clone_ref(py)),
            };
            let mut other = GraphIrAtomConstraintsForm::new();
            other.set(GraphIrAtomConstraintForm::valence(4));
            other.set(GraphIrAtomConstraintForm::degree(3));
            view.update(
                py,
                AtomConstraintsUpdate::Container(
                    Py::new(py, AtomConstraintsAst::from_inner(other)).unwrap(),
                ),
            )
            .unwrap();
            let fresh = AtomConstraintsView {
                backing: AtomConstraintsBacking::Atom(atom),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 2);
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_set_valence() {
        Python::attach(|py| {
            let mut constraints = AtomConstraintsAst::new(py, vec![]);
            constraints.set_valence(py, NumLike::Lit(4));
            assert_eq!(
                constraints.valence(py).unwrap().unwrap().to_rust(py),
                GraphIrNumForm::Lit(4)
            );
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_set_ring_count() {
        Python::attach(|py| {
            let mut constraints = AtomConstraintsAst::new(py, vec![]);
            constraints.set_ring_count(py, NumLike::Lit(2));
            assert_eq!(
                constraints.ring_count(py).unwrap().unwrap().to_rust(py),
                GraphIrNumForm::Lit(2)
            );
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_set_aromatic_valence() {
        Python::attach(|py| {
            let mut constraints = AtomConstraintsAst::new(py, vec![]);
            constraints
                .set_aromatic_valence(py, AromaticValenceLike::Value(NumLike::Lit(1)))
                .unwrap();
            match constraints.aromatic_valence(py).unwrap().unwrap() {
                AromaticValenceAst::Aromatic(v) => {
                    assert_eq!(v.to_rust(py), GraphIrNumForm::Lit(1))
                }
                _ => panic!("expected Aromatic"),
            }
            constraints
                .set_aromatic_valence(py, AromaticValenceLike::Flag(false))
                .unwrap();
            match constraints.aromatic_valence(py).unwrap().unwrap() {
                AromaticValenceAst::NotAromatic() => {}
                _ => panic!("expected NotAromatic"),
            }
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_set_aromatic_valence_error() {
        Python::attach(|py| {
            let mut constraints = AtomConstraintsAst::new(py, vec![]);
            assert!(constraints
                .set_aromatic_valence(py, AromaticValenceLike::Flag(true))
                .is_err());
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_set_tetrahedral_stereo() {
        Python::attach(|py| {
            let mut constraints = AtomConstraintsAst::new(py, vec![]);
            constraints
                .set_tetrahedral_stereo(
                    py,
                    TetrahedralStereoLike::Config(TetrahedralConfiguration::Cw),
                )
                .unwrap();
            match constraints.tetrahedral_stereo(py).unwrap().unwrap() {
                TetrahedralStereoForm::Stereo(coset) => {
                    assert_eq!(
                        coset.bind(py).borrow().to_rust(py),
                        GraphIrStereoCoset::Lit(1)
                    )
                }
                _ => panic!("expected Stereo"),
            }
        });
    }

    #[rstest]
    fn test_atom_constraints_view_set_aromatic_valence() {
        Python::attach(|py| {
            let owner = Py::new(
                py,
                MoleculeAst::from_rust(GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                    ..Default::default()
                })),
            )
            .unwrap();
            let view = AtomConstraintsView {
                backing: AtomConstraintsBacking::Molecule {
                    owner: owner.clone_ref(py),
                    id: GraphIrAtomId(0),
                },
            };
            view.set_aromatic_valence(py, AromaticValenceLike::Value(NumLike::Lit(1)))
                .unwrap();
            let fresh = AtomConstraintsView {
                backing: AtomConstraintsBacking::Molecule {
                    owner,
                    id: GraphIrAtomId(0),
                },
            };
            match fresh.aromatic_valence(py).unwrap().unwrap() {
                AromaticValenceAst::Aromatic(v) => {
                    assert_eq!(v.to_rust(py), GraphIrNumForm::Lit(1))
                }
                _ => panic!("expected Aromatic"),
            }
        });
    }

    #[rstest]
    fn test_ring_size_counts_value_backed() {
        Python::attach(|py| {
            let constraints = Py::new(py, AtomConstraintsAst::new(py, vec![])).unwrap();
            let proxy = AtomConstraintsAst::ring_size_count(constraints.clone_ref(py));
            proxy.__setitem__(py, 6, NumLike::Lit(3));
            assert_eq!(
                proxy.__getitem__(py, 6).unwrap().unwrap().to_rust(py),
                GraphIrNumForm::Lit(3)
            );
            proxy.__delitem__(py, 6);
            assert!(proxy.__getitem__(py, 6).unwrap().is_none());
        });
    }

    #[rstest]
    fn test_ring_size_counts_molecule_backed() {
        Python::attach(|py| {
            let owner = Py::new(
                py,
                MoleculeAst::from_rust(GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
                    atoms: vec![GraphIrAtomForm::from_element(ChemElement::C)],
                    ..Default::default()
                })),
            )
            .unwrap();
            let view = AtomConstraintsView {
                backing: AtomConstraintsBacking::Molecule {
                    owner: owner.clone_ref(py),
                    id: GraphIrAtomId(0),
                },
            };
            view.ring_size_count(py).__setitem__(py, 5, NumLike::Lit(1));
            let fresh = AtomConstraintsView {
                backing: AtomConstraintsBacking::Molecule {
                    owner,
                    id: GraphIrAtomId(0),
                },
            };
            assert_eq!(
                fresh
                    .ring_size_count(py)
                    .__getitem__(py, 5)
                    .unwrap()
                    .unwrap()
                    .to_rust(py),
                GraphIrNumForm::Lit(1)
            );
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_update_entries() {
        Python::attach(|py| {
            let constraints = Py::new(py, AtomConstraintsAst::new(py, vec![])).unwrap();
            let valence = into_py_variant(
                py,
                AtomConstraintAst::from_rust(py, &GraphIrAtomConstraintForm::valence(4)).unwrap(),
            )
            .unwrap();
            let degree = into_py_variant(
                py,
                AtomConstraintAst::from_rust(py, &GraphIrAtomConstraintForm::degree(3)).unwrap(),
            )
            .unwrap();
            AtomConstraintsAst::update(
                constraints.clone_ref(py),
                py,
                AtomConstraintsUpdate::Entries(vec![valence, degree]),
            )
            .unwrap();
            assert_eq!(constraints.bind(py).borrow().__len__(), 2);
        });
    }

    /// Regression: a container updating itself resolves `other` before the write borrow,
    /// so it is an idempotent no-op, not a RefCell double-borrow panic.
    #[rstest]
    fn test_atom_constraints_ast_update_self() {
        Python::attach(|py| {
            let valence = into_py_variant(
                py,
                AtomConstraintAst::from_rust(py, &GraphIrAtomConstraintForm::valence(4)).unwrap(),
            )
            .unwrap();
            let constraints = Py::new(py, AtomConstraintsAst::new(py, vec![valence])).unwrap();
            AtomConstraintsAst::update(
                constraints.clone_ref(py),
                py,
                AtomConstraintsUpdate::Container(constraints.clone_ref(py)),
            )
            .unwrap();
            assert_eq!(
                constraints
                    .bind(py)
                    .borrow()
                    .valence(py)
                    .unwrap()
                    .unwrap()
                    .to_rust(py),
                GraphIrNumForm::Lit(4)
            );
        });
    }

    /// Regression: a view updating from a view over the same atom resolves `other`
    /// before the write borrow, so it is an idempotent no-op, not a double-borrow panic
    /// (`atom.constraints.update(atom.constraints)`).
    #[rstest]
    fn test_atom_constraints_view_update_self() {
        Python::attach(|py| {
            let atom = Py::new(
                py,
                AtomAst::from_inner(
                    GraphIrAtomForm::from_element(ChemElement::C)
                        .with_constraint(GraphIrAtomConstraintForm::valence(4)),
                ),
            )
            .unwrap();
            let view = AtomConstraintsView {
                backing: AtomConstraintsBacking::Atom(atom.clone_ref(py)),
            };
            let other = Py::new(
                py,
                AtomConstraintsView {
                    backing: AtomConstraintsBacking::Atom(atom.clone_ref(py)),
                },
            )
            .unwrap();
            view.update(py, AtomConstraintsUpdate::View(other)).unwrap();
            assert_eq!(
                atom.bind(py)
                    .borrow()
                    .inner()
                    .constraints
                    .valence()
                    .unwrap()
                    .clone(),
                GraphIrNumForm::Lit(4)
            );
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_getitem_error() {
        Python::attach(|py| {
            let constraints = AtomConstraintsAst::new(py, vec![]);
            let key = into_py_variant(py, AtomConstraintKey::Valence()).unwrap();
            assert!(constraints.__getitem__(py, key).is_err());
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_delitem_error() {
        Python::attach(|py| {
            let mut constraints = AtomConstraintsAst::new(py, vec![]);
            let key = into_py_variant(py, AtomConstraintKey::Valence()).unwrap();
            assert!(constraints.__delitem__(py, key).is_err());
        });
    }

    #[rstest]
    fn test_ring_size_counts_len_iter_contains() {
        Python::attach(|py| {
            let constraints = Py::new(py, AtomConstraintsAst::new(py, vec![])).unwrap();
            let proxy = AtomConstraintsAst::ring_size_count(constraints.clone_ref(py));
            proxy.__setitem__(py, 6, NumLike::Lit(3));
            proxy.__setitem__(py, 5, NumLike::Lit(1));
            assert_eq!(proxy.__len__(py).unwrap(), 2);
            assert!(proxy.__contains__(py, 6).unwrap());
            assert!(!proxy.__contains__(py, 4).unwrap());
            let mut iter = proxy.__iter__(py).unwrap();
            let mut sizes = Vec::new();
            while let Some(size) = iter.__next__() {
                sizes.push(size);
            }
            sizes.sort_unstable();
            assert_eq!(sizes, vec![5, 6]);
        });
    }
}
