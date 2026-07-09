//! Atom-field value types and the atom read surface mirroring `umol_ast::ast`:
//! `ElementAst`, `IsotopeMassAst`, `SpinStateAst`, `AtomAst`, and the
//! `AtomView`/`AtomViews` handle views.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use std::collections::BTreeSet;
use std::str::FromStr;
use std::vec::IntoIter;

use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use umol_ast::ast::{
    AsLit, AtomAst as AstAtomAst, AtomId as AstAtomId, ElementAst as AstElementAst,
    IsotopeMassAst as AstIsotopeMassAst, MoleculeAst as AstMoleculeAst,
    SpinStateAst as AstSpinStateAst,
};
use umol_chem::element::Element as ChemElement;

use crate::constraint::{
    atom_constraints_asdict, AtomConstraintsAst, AtomConstraintsView, ConstraintsArg,
    ConstraintsBacking,
};
use crate::convert::{hash_ast, into_py_variant, variant_repr};
use crate::element::Element;
use crate::error::parse_error;
use crate::molecule::MoleculeAst;
use crate::value::{MemOp, ValueArg, ValueAst};

/// Element expression: undetermined, a single element, a finite element set, a
/// complement set (`!{…}`), or a variable with an optional membership restriction.
#[pyclass]
pub enum ElementAst {
    Undetermined(),
    Lit(Element),
    LitSet(BTreeSet<Element>),
    NotSet(BTreeSet<Element>),
    Var(String, Option<(MemOp, BTreeSet<Element>)>),
}

#[pymethods]
impl ElementAst {
    /// The single element this resolves to, or `None` when it is not a bare
    /// literal (undetermined, a set, a complement, or a variable).
    fn as_lit(&self) -> Option<Element> {
        self.to_ast().as_lit().map(Element::from)
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.to_ast() == other.to_ast()
    }

    fn __hash__(&self) -> u64 {
        hash_ast(&self.to_ast())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            ElementAst::Undetermined() => ("Undetermined", 0),
            ElementAst::Lit(_) => ("Lit", 1),
            ElementAst::LitSet(_) => ("LitSet", 1),
            ElementAst::NotSet(_) => ("NotSet", 1),
            ElementAst::Var(_, _) => ("Var", 2),
        };
        variant_repr(slf.bind(py).as_any(), "ElementAst", variant, arity)
    }
}

impl ElementAst {
    pub(crate) fn from_ast(ast: &AstElementAst) -> ElementAst {
        match ast {
            AstElementAst::Undetermined => ElementAst::Undetermined(),
            AstElementAst::Lit(e) => ElementAst::Lit(Element::from(*e)),
            AstElementAst::LitSet(members) => {
                ElementAst::LitSet(members.iter().copied().map(Element::from).collect())
            }
            AstElementAst::NotSet(members) => {
                ElementAst::NotSet(members.iter().copied().map(Element::from).collect())
            }
            AstElementAst::Var(boxed) => {
                let (name, restriction) = &**boxed;
                ElementAst::Var(
                    name.clone(),
                    restriction.as_ref().map(|(op, members)| {
                        (
                            MemOp::from_ast(*op),
                            members.iter().copied().map(Element::from).collect(),
                        )
                    }),
                )
            }
        }
    }

    pub(crate) fn to_ast(&self) -> AstElementAst {
        match self {
            ElementAst::Undetermined() => AstElementAst::Undetermined,
            ElementAst::Lit(e) => AstElementAst::Lit(ChemElement::from(e)),
            ElementAst::LitSet(members) => {
                AstElementAst::LitSet(Box::new(members.iter().map(ChemElement::from).collect()))
            }
            ElementAst::NotSet(members) => {
                AstElementAst::NotSet(Box::new(members.iter().map(ChemElement::from).collect()))
            }
            ElementAst::Var(name, restriction) => AstElementAst::Var(Box::new((
                name.clone(),
                restriction.as_ref().map(|(op, members)| {
                    (op.to_ast(), members.iter().map(ChemElement::from).collect())
                }),
            ))),
        }
    }
}

/// Isotope-mass expression: undetermined, the natural isotopic mixture, a single
/// mass number, a finite mass set, or a variable with an optional mass-set restriction.
#[pyclass]
pub enum IsotopeMassAst {
    Undetermined(),
    Natural(),
    Lit(u32),
    LitSet(BTreeSet<u32>),
    Var(String, Option<BTreeSet<u32>>),
}

#[pymethods]
impl IsotopeMassAst {
    /// The single mass number this resolves to, or `None` when it is not a bare
    /// literal (undetermined, the natural mixture, a set, or a variable).
    fn as_lit(&self) -> Option<u32> {
        self.to_ast().as_lit()
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.to_ast() == other.to_ast()
    }

    fn __hash__(&self) -> u64 {
        hash_ast(&self.to_ast())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            IsotopeMassAst::Undetermined() => ("Undetermined", 0),
            IsotopeMassAst::Natural() => ("Natural", 0),
            IsotopeMassAst::Lit(_) => ("Lit", 1),
            IsotopeMassAst::LitSet(_) => ("LitSet", 1),
            IsotopeMassAst::Var(_, _) => ("Var", 2),
        };
        variant_repr(slf.bind(py).as_any(), "IsotopeMassAst", variant, arity)
    }
}

impl IsotopeMassAst {
    pub(crate) fn from_ast(ast: &AstIsotopeMassAst) -> IsotopeMassAst {
        match ast {
            AstIsotopeMassAst::Undetermined => IsotopeMassAst::Undetermined(),
            AstIsotopeMassAst::Natural => IsotopeMassAst::Natural(),
            AstIsotopeMassAst::Lit(mass) => IsotopeMassAst::Lit(*mass),
            AstIsotopeMassAst::LitSet(masses) => IsotopeMassAst::LitSet((**masses).clone()),
            AstIsotopeMassAst::Var(boxed) => {
                let (name, restriction) = &**boxed;
                IsotopeMassAst::Var(name.clone(), restriction.clone())
            }
        }
    }

    pub(crate) fn to_ast(&self) -> AstIsotopeMassAst {
        match self {
            IsotopeMassAst::Undetermined() => AstIsotopeMassAst::Undetermined,
            IsotopeMassAst::Natural() => AstIsotopeMassAst::Natural,
            IsotopeMassAst::Lit(mass) => AstIsotopeMassAst::Lit(*mass),
            IsotopeMassAst::LitSet(masses) => AstIsotopeMassAst::LitSet(Box::new(masses.clone())),
            IsotopeMassAst::Var(name, restriction) => {
                AstIsotopeMassAst::Var(Box::new((name.clone(), restriction.clone())))
            }
        }
    }
}

/// Spin state: unpaired-electron count and multiplicity as independent value fields.
#[pyclass]
pub struct SpinStateAst {
    #[pyo3(get)]
    unpaired: Py<ValueAst>,
    #[pyo3(get)]
    multiplicity: Py<ValueAst>,
}

#[pymethods]
impl SpinStateAst {
    #[new]
    fn new(py: Python<'_>, unpaired: ValueArg, multiplicity: ValueArg) -> PyResult<Self> {
        Ok(SpinStateAst {
            unpaired: unpaired.to_py(py)?,
            multiplicity: multiplicity.to_py(py)?,
        })
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "SpinStateAst({}, {})",
            self.unpaired.bind(py).as_any().repr()?.extract::<String>()?,
            self.multiplicity
                .bind(py)
                .as_any()
                .repr()?
                .extract::<String>()?,
        ))
    }
}

impl SpinStateAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstSpinStateAst) -> PyResult<SpinStateAst> {
        Ok(SpinStateAst {
            unpaired: into_py_variant(py, ValueAst::from_ast(py, &ast.unpaired)?)?,
            multiplicity: into_py_variant(py, ValueAst::from_ast(py, &ast.multiplicity)?)?,
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstSpinStateAst {
        AstSpinStateAst {
            unpaired: self.unpaired.bind(py).borrow().to_ast(py),
            multiplicity: self.multiplicity.bind(py).borrow().to_ast(py),
        }
    }
}

/// An atom: element, isotope, charge, implicit hydrogens, lone pairs, spin, and
/// atom-scope constraints.
#[pyclass(eq)]
#[derive(PartialEq)]
pub struct AtomAst(AstAtomAst);

#[pymethods]
impl AtomAst {
    /// Construct from an element — a single `Element` or an `ElementAst` expression —
    /// optionally setting fields.
    #[new]
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (element, *, isotope_mass=None, charge=None, implicit_hydrogens=None, lone_pairs=None, spin=None, constraints=None))]
    fn new(
        py: Python<'_>,
        element: ElementArg,
        isotope_mass: Option<IsotopeMassArg>,
        charge: Option<ValueArg>,
        implicit_hydrogens: Option<ValueArg>,
        lone_pairs: Option<ValueArg>,
        spin: Option<PyRef<'_, SpinStateAst>>,
        constraints: Option<Py<AtomConstraintsAst>>,
    ) -> Self {
        let atom = AstAtomAst::new(element.to_ast(py));
        AtomAst(apply_fields(
            atom,
            py,
            isotope_mass,
            charge,
            implicit_hydrogens,
            lone_pairs,
            spin,
            constraints,
        ))
    }

    /// Parse an atom-DSL string (e.g. `"C#c-1#v4"`) into an `AtomAst`.
    #[staticmethod]
    fn parse(s: &str) -> PyResult<Self> {
        AstAtomAst::from_str(s).map(Self).map_err(parse_error)
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("AtomAst.parse('{}')", self.0)
    }

    #[getter]
    fn element(&self) -> ElementAst {
        ElementAst::from_ast(&self.0.element)
    }

    #[setter]
    fn set_element(&mut self, py: Python<'_>, value: ElementArg) {
        self.0.element = value.to_ast(py);
    }

    #[getter]
    fn isotope_mass(&self) -> IsotopeMassAst {
        IsotopeMassAst::from_ast(&self.0.isotope_mass)
    }

    #[setter]
    fn set_isotope_mass(&mut self, py: Python<'_>, value: IsotopeMassArg) {
        self.0.isotope_mass = value.to_ast(py);
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
    fn implicit_hydrogens(&self, py: Python<'_>) -> PyResult<ValueAst> {
        ValueAst::from_ast(py, &self.0.implicit_hydrogens)
    }

    #[setter]
    fn set_implicit_hydrogens(&mut self, py: Python<'_>, value: ValueArg) {
        self.0.implicit_hydrogens = value.to_ast(py);
    }


    #[getter]
    fn lone_pairs(&self, py: Python<'_>) -> PyResult<ValueAst> {
        ValueAst::from_ast(py, &self.0.lone_pairs)
    }

    #[setter]
    fn set_lone_pairs(&mut self, py: Python<'_>, value: ValueArg) {
        self.0.lone_pairs = value.to_ast(py);
    }

    #[getter]
    fn spin(&self, py: Python<'_>) -> PyResult<SpinStateAst> {
        SpinStateAst::from_ast(py, &self.0.spin)
    }

    #[setter]
    fn set_spin(&mut self, py: Python<'_>, value: PyRef<'_, SpinStateAst>) {
        self.0.spin = value.to_ast(py);
    }

    /// The atom's constraints as a live handle onto this atom: reads borrow the
    /// current state, mutators write through to the atom in place.
    #[getter]
    fn constraints(slf: Py<Self>) -> AtomConstraintsView {
        AtomConstraintsView {
            backing: ConstraintsBacking::Atom(slf),
        }
    }

    /// Replace the whole constraint set (wipe-and-set) from a value container or
    /// a live view.
    #[setter]
    fn set_constraints(&mut self, py: Python<'_>, value: ConstraintsArg) -> PyResult<()> {
        self.0.constraints = value.to_ast(py)?;
        Ok(())
    }

    /// The fields as a dict keyed by field name; values are the field mirrors.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("element", self.element())?;
        dict.set_item("isotope_mass", self.isotope_mass())?;
        dict.set_item("charge", self.charge(py)?)?;
        dict.set_item("implicit_hydrogens", self.implicit_hydrogens(py)?)?;
        dict.set_item("lone_pairs", self.lone_pairs(py)?)?;
        dict.set_item("spin", self.spin(py)?)?;
        dict.set_item(
            "constraints",
            atom_constraints_asdict(py, &self.0.constraints)?,
        )?;
        Ok(dict)
    }
}

/// A binding argument that coerces a literal *or* a mirror to the AST — the `*Arg`
/// convention for these inputs (`*Input` is reserved for the DSL side). Extracted as
/// a PyO3 `FromPyObject` union tried in order; variants are `Ast` = the `*Ast`
/// wrapper, `Lit` = the literal, mirroring `impl Into<..>` on the Rust builders.
///
/// `ElementArg` accepts a concrete `Element` or an `ElementAst`.
#[derive(FromPyObject)]
enum ElementArg {
    Ast(Py<ElementAst>),
    Lit(Element),
}

impl ElementArg {
    fn to_ast(&self, py: Python<'_>) -> AstElementAst {
        match self {
            ElementArg::Ast(expr) => expr.bind(py).borrow().to_ast(),
            ElementArg::Lit(element) => AstElementAst::Lit(ChemElement::from(element)),
        }
    }
}

/// An `IsotopeMassAst` or a Python `int` (→ `IsotopeMassAst::Lit`, a mass number).
#[derive(FromPyObject)]
enum IsotopeMassArg {
    Ast(Py<IsotopeMassAst>),
    Lit(u32),
}

impl IsotopeMassArg {
    fn to_ast(&self, py: Python<'_>) -> AstIsotopeMassAst {
        match self {
            IsotopeMassArg::Ast(mass) => mass.bind(py).borrow().to_ast(),
            IsotopeMassArg::Lit(number) => AstIsotopeMassAst::Lit(*number),
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn apply_fields(
    mut atom: AstAtomAst,
    py: Python<'_>,
    isotope_mass: Option<IsotopeMassArg>,
    charge: Option<ValueArg>,
    implicit_hydrogens: Option<ValueArg>,
    lone_pairs: Option<ValueArg>,
    spin: Option<PyRef<'_, SpinStateAst>>,
    constraints: Option<Py<AtomConstraintsAst>>,
) -> AstAtomAst {
    if let Some(isotope_mass) = isotope_mass {
        atom = atom.with_isotope_mass(isotope_mass.to_ast(py));
    }
    if let Some(charge) = charge {
        atom = atom.with_charge(charge.to_ast(py));
    }
    if let Some(implicit_hydrogens) = implicit_hydrogens {
        atom = atom.with_implicit_hydrogens(implicit_hydrogens.to_ast(py));
    }
    if let Some(lone_pairs) = lone_pairs {
        atom = atom.with_lone_pairs(lone_pairs.to_ast(py));
    }
    if let Some(spin) = spin {
        atom = atom.with_spin(spin.to_ast(py));
    }
    if let Some(constraints) = constraints {
        atom.constraints = constraints.bind(py).borrow().inner().clone();
    }
    atom
}

impl AtomAst {
    /// The wrapped AST atom — read access for molecule construction.
    pub(crate) fn inner(&self) -> &AstAtomAst {
        &self.0
    }

    /// Mutable access to the wrapped AST atom — write access for the atom-backed
    /// constraints view.
    pub(crate) fn inner_mut(&mut self) -> &mut AstAtomAst {
        &mut self.0
    }

    /// Wrap an AST atom (the hold-the-value `from_inner` bridge, paired with
    /// `inner`). Test-only — in-crate construction wraps `AtomAst(..)` directly.
    #[cfg(test)]
    pub(crate) fn from_inner(atom: AstAtomAst) -> Self {
        AtomAst(atom)
    }
}

/// A view of one atom within a molecule: a handle to the molecule plus the atom's
/// index. Field reads rebuild the transient Rust view; the molecule is never copied.
#[pyclass]
pub struct AtomView {
    owner: Py<MoleculeAst>,
    id: AstAtomId,
}

impl AtomView {
    fn atom<'a>(&self, molecule: &'a AstMoleculeAst) -> PyResult<&'a AstAtomAst> {
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
    fn element(&self, py: Python<'_>) -> PyResult<ElementAst> {
        let molecule = self.owner.bind(py).borrow();
        Ok(ElementAst::from_ast(&self.atom(molecule.inner())?.element))
    }

    #[setter]
    fn set_element(&self, py: Python<'_>, value: ElementArg) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .atom_mut(self.id)
            .ast
            .element = value.to_ast(py);
    }

    #[getter]
    fn isotope_mass(&self, py: Python<'_>) -> PyResult<IsotopeMassAst> {
        let molecule = self.owner.bind(py).borrow();
        Ok(IsotopeMassAst::from_ast(
            &self.atom(molecule.inner())?.isotope_mass,
        ))
    }

    #[setter]
    fn set_isotope_mass(&self, py: Python<'_>, value: IsotopeMassArg) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .atom_mut(self.id)
            .ast
            .isotope_mass = value.to_ast(py);
    }

    #[getter]
    fn charge(&self, py: Python<'_>) -> PyResult<ValueAst> {
        let molecule = self.owner.bind(py).borrow();
        ValueAst::from_ast(py, &self.atom(molecule.inner())?.charge)
    }

    #[setter]
    fn set_charge(&self, py: Python<'_>, value: ValueArg) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .atom_mut(self.id)
            .ast
            .charge = value.to_ast(py);
    }

    #[getter]
    fn implicit_hydrogens(&self, py: Python<'_>) -> PyResult<ValueAst> {
        let molecule = self.owner.bind(py).borrow();
        ValueAst::from_ast(py, &self.atom(molecule.inner())?.implicit_hydrogens)
    }

    #[setter]
    fn set_implicit_hydrogens(&self, py: Python<'_>, value: ValueArg) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .atom_mut(self.id)
            .ast
            .implicit_hydrogens = value.to_ast(py);
    }

    #[getter]
    fn lone_pairs(&self, py: Python<'_>) -> PyResult<ValueAst> {
        let molecule = self.owner.bind(py).borrow();
        ValueAst::from_ast(py, &self.atom(molecule.inner())?.lone_pairs)
    }

    #[setter]
    fn set_lone_pairs(&self, py: Python<'_>, value: ValueArg) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .atom_mut(self.id)
            .ast
            .lone_pairs = value.to_ast(py);
    }

    #[getter]
    fn spin(&self, py: Python<'_>) -> PyResult<SpinStateAst> {
        let molecule = self.owner.bind(py).borrow();
        SpinStateAst::from_ast(py, &self.atom(molecule.inner())?.spin)
    }

    #[setter]
    fn set_spin(&self, py: Python<'_>, value: PyRef<'_, SpinStateAst>) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .atom_mut(self.id)
            .ast
            .spin = value.to_ast(py);
    }

    /// The atom's constraints as a live handle onto the molecule: reads borrow the
    /// current state, mutators write through to the atom in place.
    #[getter]
    fn constraints(&self, py: Python<'_>) -> AtomConstraintsView {
        AtomConstraintsView {
            backing: ConstraintsBacking::Molecule {
                owner: self.owner.clone_ref(py),
                id: self.id,
            },
        }
    }

    /// Replace the whole constraint set of the backing atom in place (wipe-and-set)
    /// from a value container or a live view.
    #[setter]
    fn set_constraints(&self, py: Python<'_>, value: ConstraintsArg) -> PyResult<()> {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .atom_mut(self.id)
            .ast
            .constraints = value.to_ast(py)?;
        Ok(())
    }

    /// The fields as a dict keyed by field name; values are the field mirrors —
    /// symmetric with `AtomAst.asdict`, read through the view.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let molecule = self.owner.bind(py).borrow();
        let atom = self.atom(molecule.inner())?;
        let dict = PyDict::new(py);
        dict.set_item("element", ElementAst::from_ast(&atom.element))?;
        dict.set_item("isotope_mass", IsotopeMassAst::from_ast(&atom.isotope_mass))?;
        dict.set_item("charge", ValueAst::from_ast(py, &atom.charge)?)?;
        dict.set_item(
            "implicit_hydrogens",
            ValueAst::from_ast(py, &atom.implicit_hydrogens)?,
        )?;
        dict.set_item("lone_pairs", ValueAst::from_ast(py, &atom.lone_pairs)?)?;
        dict.set_item("spin", SpinStateAst::from_ast(py, &atom.spin)?)?;
        dict.set_item(
            "constraints",
            atom_constraints_asdict(py, &atom.constraints)?,
        )?;
        Ok(dict)
    }
}

/// Resolve a possibly-negative Python index (negative counts from the end) into an
/// existing atom id, or `IndexError`.
fn resolve_atom_index(molecule: &AstMoleculeAst, index: isize) -> PyResult<AstAtomId> {
    let count = molecule.atoms().count();
    let resolved = if index < 0 { index + count as isize } else { index };
    if resolved < 0 {
        return Err(PyIndexError::new_err("atom id out of range"));
    }
    let id = AstAtomId(resolved as u32);
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
    ids: IntoIter<AstAtomId>,
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
    use umol_ast::ast::{
        AtomConstraintAst as AstAtomConstraintAst, MemOp as AstMemOp, ValueAst as AstValueAst,
    };

    use super::*;

    #[rstest]
    #[case(AstElementAst::Undetermined)]
    #[case(AstElementAst::Lit(ChemElement::C))]
    #[case(AstElementAst::LitSet(Box::new(BTreeSet::from([ChemElement::C, ChemElement::N]))))]
    #[case(AstElementAst::NotSet(Box::new(BTreeSet::from([ChemElement::O]))))]
    #[case(AstElementAst::Var(Box::new(("x".to_string(), None))))]
    #[case(AstElementAst::Var(Box::new((
        "y".to_string(),
        Some((AstMemOp::In, BTreeSet::from([ChemElement::C, ChemElement::N]))),
    ))))]
    fn test_element_ast_roundtrip(#[case] ast: AstElementAst) {
        assert_eq!(ElementAst::from_ast(&ast).to_ast(), ast);
    }

    #[rstest]
    #[case(AstElementAst::Lit(ChemElement::C), Some(ChemElement::C))]
    #[case(AstElementAst::Undetermined, None)]
    #[case(AstElementAst::LitSet(Box::new(BTreeSet::from([ChemElement::C, ChemElement::N]))), None)]
    fn test_element_ast_as_lit(#[case] ast: AstElementAst, #[case] expected: Option<ChemElement>) {
        let got = ElementAst::from_ast(&ast)
            .as_lit()
            .map(|e| ChemElement::from(&e));
        assert_eq!(got, expected);
    }

    #[rstest]
    #[case(AstIsotopeMassAst::Undetermined)]
    #[case(AstIsotopeMassAst::Natural)]
    #[case(AstIsotopeMassAst::Lit(13))]
    #[case(AstIsotopeMassAst::LitSet(Box::new(BTreeSet::from([12, 13, 14]))))]
    #[case(AstIsotopeMassAst::Var(Box::new(("x".to_string(), None))))]
    #[case(AstIsotopeMassAst::Var(Box::new((
        "y".to_string(),
        Some(BTreeSet::from([12, 13])),
    ))))]
    fn test_isotope_mass_ast_roundtrip(#[case] ast: AstIsotopeMassAst) {
        assert_eq!(IsotopeMassAst::from_ast(&ast).to_ast(), ast);
    }

    #[rstest]
    #[case(AstIsotopeMassAst::Lit(13), Some(13))]
    #[case(AstIsotopeMassAst::Natural, None)]
    #[case(AstIsotopeMassAst::Undetermined, None)]
    fn test_isotope_mass_ast_as_lit(#[case] ast: AstIsotopeMassAst, #[case] expected: Option<u32>) {
        assert_eq!(IsotopeMassAst::from_ast(&ast).as_lit(), expected);
    }

    #[rstest]
    #[case(AstSpinStateAst { unpaired: AstValueAst::Lit(1), multiplicity: AstValueAst::Lit(2) })]
    #[case(AstSpinStateAst {
        unpaired: AstValueAst::Undetermined,
        multiplicity: AstValueAst::Undetermined,
    })]
    fn test_spin_state_ast_roundtrip(#[case] ast: AstSpinStateAst) {
        Python::attach(|py| {
            assert_eq!(SpinStateAst::from_ast(py, &ast).unwrap().to_ast(py), ast);
        });
    }

    fn carbon_oxygen(py: Python<'_>) -> Py<MoleculeAst> {
        let molecule = AstMoleculeAst::from_atoms_and_bonds(
            vec![
                AstAtomAst::from_element(ChemElement::C),
                AstAtomAst::from_element(ChemElement::O),
            ],
            vec![],
        );
        Py::new(py, MoleculeAst::from_inner(molecule)).unwrap()
    }

    #[rstest]
    fn test_atom_view_element() {
        Python::attach(|py| {
            let view = AtomView {
                owner: carbon_oxygen(py),
                id: AstAtomId(1),
            };
            assert_eq!(view.id(), 1);
            match view.element(py).unwrap() {
                ElementAst::Lit(e) => assert_eq!(ChemElement::from(&e), ChemElement::O),
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
                id: AstAtomId(0),
            };
            view.set_charge(py, ValueArg::Lit(-1));
            let fresh = AtomView {
                owner,
                id: AstAtomId(0),
            };
            match fresh.charge(py).unwrap() {
                ValueAst::Lit(n) => assert_eq!(n, -1),
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
                id: AstAtomId(0),
            };
            view.set_element(py, ElementArg::Lit(Element::from(ChemElement::N)));
            let fresh = AtomView {
                owner,
                id: AstAtomId(0),
            };
            match fresh.element(py).unwrap() {
                ElementAst::Lit(e) => assert_eq!(ChemElement::from(&e), ChemElement::N),
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
                id: AstAtomId(0),
            };
            view.set_isotope_mass(py, IsotopeMassArg::Lit(13));
            let fresh = AtomView {
                owner,
                id: AstAtomId(0),
            };
            match fresh.isotope_mass(py).unwrap() {
                IsotopeMassAst::Lit(mass) => assert_eq!(mass, 13),
                _ => panic!("expected Lit"),
            }
        });
    }

    #[rstest]
    fn test_atom_view_set_spin() {
        Python::attach(|py| {
            let owner = carbon_oxygen(py);
            let view = AtomView {
                owner: owner.clone_ref(py),
                id: AstAtomId(0),
            };
            let spin = Py::new(
                py,
                SpinStateAst::from_ast(
                    py,
                    &AstSpinStateAst {
                        unpaired: AstValueAst::Lit(1),
                        multiplicity: AstValueAst::Lit(2),
                    },
                )
                .unwrap(),
            )
            .unwrap();
            view.set_spin(py, spin.bind(py).borrow());
            let fresh = AtomView {
                owner,
                id: AstAtomId(0),
            };
            assert_eq!(
                fresh.spin(py).unwrap().to_ast(py),
                AstSpinStateAst {
                    unpaired: AstValueAst::Lit(1),
                    multiplicity: AstValueAst::Lit(2),
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
            let nitrogen =
                Py::new(py, AtomAst::from_inner(AstAtomAst::from_element(ChemElement::N))).unwrap();
            views.__setitem__(py, 0, nitrogen.bind(py).borrow()).unwrap();
            match views.__getitem__(py, 0).unwrap().element(py).unwrap() {
                ElementAst::Lit(e) => assert_eq!(ChemElement::from(&e), ChemElement::N),
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
            let nitrogen =
                Py::new(py, AtomAst::from_inner(AstAtomAst::from_element(ChemElement::N))).unwrap();
            assert!(views
                .__setitem__(py, 5, nitrogen.bind(py).borrow())
                .is_err());
        });
    }

    #[rstest]
    fn test_atom_ast_constraints() {
        let atom = AtomAst(
            AstAtomAst::from_element(ChemElement::C)
                .with_constraint(AstAtomConstraintAst::valence(4)),
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
                    AstAtomAst::from_element(ChemElement::C)
                        .with_constraint(AstAtomConstraintAst::valence(4)),
                ),
            )
            .unwrap();
            let view = Py::new(
                py,
                AtomConstraintsView {
                    backing: ConstraintsBacking::Atom(src),
                },
            )
            .unwrap();
            let mut dst = AtomAst::from_inner(AstAtomAst::from_element(ChemElement::N));
            dst.set_constraints(py, ConstraintsArg::View(view)).unwrap();
            assert_eq!(
                dst.inner().constraints.valence().unwrap().clone(),
                AstValueAst::Lit(4)
            );
        });
    }

    #[rstest]
    fn test_atom_view_constraints() {
        Python::attach(|py| {
            let view = AtomView {
                owner: carbon_oxygen(py),
                id: AstAtomId(1),
            };
            match view.constraints(py).backing {
                ConstraintsBacking::Molecule { id, .. } => assert_eq!(id, AstAtomId(1)),
                _ => panic!("expected molecule-backed view"),
            }
        });
    }
}
