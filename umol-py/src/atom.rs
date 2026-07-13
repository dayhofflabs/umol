//! Atom-field value types and the atom read surface mirroring `umol_ast::ast`:
//! `ElementAst`, `IsotopeMassAst`, `SpinStateAst`, `AtomAst`, the `AtomView`/`AtomViews`
//! handle views, and the atom-constraint surface (`AtomConstraintKey`,
//! `AtomConstraintAst`, the `AtomConstraintsAst` container, the `AtomConstraintsView`
//! live handle, and `AtomRingSizeCounts`). The shared constraint value/scope leaves
//! (aromatic/multicenter valence, ring scope, ring membership) live in `constraint`.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use std::collections::BTreeSet;
use std::str::FromStr;
use std::vec::IntoIter;

use pyo3::exceptions::{PyIndexError, PyKeyError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};
#[cfg(test)]
use umol_ast::ast::MoleculeParts as AstMoleculeParts;
use umol_ast::ast::{
    AsLit, AtomAst as AstAtomAst, AtomConstraintAst as AstAtomConstraintAst,
    AtomConstraintKey as AstAtomConstraintKey, AtomConstraintsAst as AstAtomConstraintsAst,
    AtomId as AstAtomId, ElementAst as AstElementAst, IsotopeMassAst as AstIsotopeMassAst,
    MoleculeAst as AstMoleculeAst, RingScope as AstRingScope, SpinStateAst as AstSpinStateAst,
};
use umol_chem::element::Element as ChemElement;

use crate::constraint::{
    AromaticValenceArg, AromaticValenceAst, MulticenterValenceArg, MulticenterValenceAst,
    RingMembershipAst, RingScope, TetrahedralStereoArg,
};
use crate::convert::{hash_ast, into_py_variant, variant_repr};
use crate::element::Element;
use crate::error::parse_error;
use crate::molecule::MoleculeAst;
use crate::stereo::TetrahedralStereoAst;
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
            self.unpaired
                .bind(py)
                .as_any()
                .repr()?
                .extract::<String>()?,
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
            backing: AtomConstraintsBacking::Atom(slf),
        }
    }

    /// Replace the whole constraint set (wipe-and-set) from a value container or
    /// a live view. Takes `slf` by handle and snapshots `value` *before* the write
    /// borrow, so `atom.constraints = atom.constraints` (a view over the same atom) reads
    /// while the atom is unborrowed instead of a double-borrow panic.
    #[setter]
    fn set_constraints(slf: Py<Self>, py: Python<'_>, value: AtomConstraintsArg) -> PyResult<()> {
        let snapshot = value.to_ast(py)?;
        slf.borrow_mut(py).0.constraints = snapshot;
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
            backing: AtomConstraintsBacking::Molecule {
                owner: self.owner.clone_ref(py),
                id: self.id,
            },
        }
    }

    /// Replace the whole constraint set of the backing atom in place (wipe-and-set)
    /// from a value container or a live view.
    #[setter]
    fn set_constraints(&self, py: Python<'_>, value: AtomConstraintsArg) -> PyResult<()> {
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
    let resolved = if index < 0 {
        index + count as isize
    } else {
        index
    };
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

/// The key (identity) of an atom constraint, for keyed lookup. The ring-membership
/// key carries its ring scope; all other keys are the bare discriminant.
#[pyclass]
pub enum AtomConstraintKey {
    Valence(),
    DonatedPairs(),
    AcceptedPairs(),
    AromaticValence(),
    MulticenterValence(),
    TetrahedralStereo(),
    Degree(),
    TotalDegree(),
    TotalValence(),
    RingDegree(),
    RingValence(),
    TotalHydrogens(),
    RingMembership(Py<RingScope>),
}

#[pymethods]
impl AtomConstraintKey {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            AtomConstraintKey::Valence() => ("Valence", 0),
            AtomConstraintKey::DonatedPairs() => ("DonatedPairs", 0),
            AtomConstraintKey::AcceptedPairs() => ("AcceptedPairs", 0),
            AtomConstraintKey::AromaticValence() => ("AromaticValence", 0),
            AtomConstraintKey::MulticenterValence() => ("MulticenterValence", 0),
            AtomConstraintKey::TetrahedralStereo() => ("TetrahedralStereo", 0),
            AtomConstraintKey::Degree() => ("Degree", 0),
            AtomConstraintKey::TotalDegree() => ("TotalDegree", 0),
            AtomConstraintKey::TotalValence() => ("TotalValence", 0),
            AtomConstraintKey::RingDegree() => ("RingDegree", 0),
            AtomConstraintKey::RingValence() => ("RingValence", 0),
            AtomConstraintKey::TotalHydrogens() => ("TotalHydrogens", 0),
            AtomConstraintKey::RingMembership(_) => ("RingMembership", 1),
        };
        variant_repr(slf.bind(py).as_any(), "AtomConstraintKey", variant, arity)
    }
}

impl AtomConstraintKey {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstAtomConstraintKey) -> PyResult<Self> {
        Ok(match ast {
            AstAtomConstraintKey::Valence => Self::Valence(),
            AstAtomConstraintKey::DonatedPairs => Self::DonatedPairs(),
            AstAtomConstraintKey::AcceptedPairs => Self::AcceptedPairs(),
            AstAtomConstraintKey::AromaticValence => Self::AromaticValence(),
            AstAtomConstraintKey::MulticenterValence => Self::MulticenterValence(),
            AstAtomConstraintKey::TetrahedralStereo => Self::TetrahedralStereo(),
            AstAtomConstraintKey::Degree => Self::Degree(),
            AstAtomConstraintKey::TotalDegree => Self::TotalDegree(),
            AstAtomConstraintKey::TotalValence => Self::TotalValence(),
            AstAtomConstraintKey::RingDegree => Self::RingDegree(),
            AstAtomConstraintKey::RingValence => Self::RingValence(),
            AstAtomConstraintKey::TotalHydrogens => Self::TotalHydrogens(),
            AstAtomConstraintKey::RingMembership(scope) => {
                Self::RingMembership(into_py_variant(py, RingScope::from_ast(scope))?)
            }
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstAtomConstraintKey {
        match self {
            Self::Valence() => AstAtomConstraintKey::Valence,
            Self::DonatedPairs() => AstAtomConstraintKey::DonatedPairs,
            Self::AcceptedPairs() => AstAtomConstraintKey::AcceptedPairs,
            Self::AromaticValence() => AstAtomConstraintKey::AromaticValence,
            Self::MulticenterValence() => AstAtomConstraintKey::MulticenterValence,
            Self::TetrahedralStereo() => AstAtomConstraintKey::TetrahedralStereo,
            Self::Degree() => AstAtomConstraintKey::Degree,
            Self::TotalDegree() => AstAtomConstraintKey::TotalDegree,
            Self::TotalValence() => AstAtomConstraintKey::TotalValence,
            Self::RingDegree() => AstAtomConstraintKey::RingDegree,
            Self::RingValence() => AstAtomConstraintKey::RingValence,
            Self::TotalHydrogens() => AstAtomConstraintKey::TotalHydrogens,
            Self::RingMembership(scope) => {
                AstAtomConstraintKey::RingMembership(scope.bind(py).borrow().to_ast())
            }
        }
    }
}

/// An atom-scope constraint: a predicate on a valence, degree, ring, or stereo
/// property of a single atom.
#[pyclass]
pub enum AtomConstraintAst {
    Valence(Py<ValueAst>),
    TotalValence(Py<ValueAst>),
    AromaticValence(Py<AromaticValenceAst>),
    MulticenterValence(Py<MulticenterValenceAst>),
    DonatedPairs(Py<ValueAst>),
    AcceptedPairs(Py<ValueAst>),
    Degree(Py<ValueAst>),
    TotalDegree(Py<ValueAst>),
    RingDegree(Py<ValueAst>),
    RingValence(Py<ValueAst>),
    TotalHydrogens(Py<ValueAst>),
    RingMembership(Py<RingMembershipAst>),
    TetrahedralStereo(Py<TetrahedralStereoAst>),
}

#[pymethods]
impl AtomConstraintAst {
    /// The constraint's key (identity).
    #[getter]
    fn key(&self, py: Python<'_>) -> PyResult<AtomConstraintKey> {
        AtomConstraintKey::from_ast(py, &self.to_ast(py).key())
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let variant = match &*slf.bind(py).borrow() {
            AtomConstraintAst::Valence(_) => "Valence",
            AtomConstraintAst::TotalValence(_) => "TotalValence",
            AtomConstraintAst::AromaticValence(_) => "AromaticValence",
            AtomConstraintAst::MulticenterValence(_) => "MulticenterValence",
            AtomConstraintAst::DonatedPairs(_) => "DonatedPairs",
            AtomConstraintAst::AcceptedPairs(_) => "AcceptedPairs",
            AtomConstraintAst::Degree(_) => "Degree",
            AtomConstraintAst::TotalDegree(_) => "TotalDegree",
            AtomConstraintAst::RingDegree(_) => "RingDegree",
            AtomConstraintAst::RingValence(_) => "RingValence",
            AtomConstraintAst::TotalHydrogens(_) => "TotalHydrogens",
            AtomConstraintAst::RingMembership(_) => "RingMembership",
            AtomConstraintAst::TetrahedralStereo(_) => "TetrahedralStereo",
        };
        variant_repr(slf.bind(py).as_any(), "AtomConstraintAst", variant, 1)
    }
}

impl AtomConstraintAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstAtomConstraintAst) -> PyResult<Self> {
        Ok(match ast {
            AstAtomConstraintAst::Valence(v) => {
                Self::Valence(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraintAst::TotalValence(v) => {
                Self::TotalValence(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraintAst::AromaticValence(c) => {
                Self::AromaticValence(into_py_variant(py, AromaticValenceAst::from_ast(py, c)?)?)
            }
            AstAtomConstraintAst::MulticenterValence(c) => Self::MulticenterValence(
                into_py_variant(py, MulticenterValenceAst::from_ast(py, c)?)?,
            ),
            AstAtomConstraintAst::DonatedPairs(v) => {
                Self::DonatedPairs(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraintAst::AcceptedPairs(v) => {
                Self::AcceptedPairs(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraintAst::Degree(v) => {
                Self::Degree(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraintAst::TotalDegree(v) => {
                Self::TotalDegree(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraintAst::RingDegree(v) => {
                Self::RingDegree(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraintAst::RingValence(v) => {
                Self::RingValence(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraintAst::TotalHydrogens(v) => {
                Self::TotalHydrogens(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraintAst::RingMembership(m) => {
                Self::RingMembership(into_py_variant(py, RingMembershipAst::from_ast(py, m)?)?)
            }
            AstAtomConstraintAst::TetrahedralStereo(c) => Self::TetrahedralStereo(into_py_variant(
                py,
                TetrahedralStereoAst::from_ast(py, c)?,
            )?),
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstAtomConstraintAst {
        match self {
            Self::Valence(v) => AstAtomConstraintAst::Valence(v.bind(py).borrow().to_ast(py)),
            Self::TotalValence(v) => {
                AstAtomConstraintAst::TotalValence(v.bind(py).borrow().to_ast(py))
            }
            Self::AromaticValence(c) => {
                AstAtomConstraintAst::AromaticValence(c.bind(py).borrow().to_ast(py))
            }
            Self::MulticenterValence(c) => {
                AstAtomConstraintAst::MulticenterValence(c.bind(py).borrow().to_ast(py))
            }
            Self::DonatedPairs(v) => {
                AstAtomConstraintAst::DonatedPairs(v.bind(py).borrow().to_ast(py))
            }
            Self::AcceptedPairs(v) => {
                AstAtomConstraintAst::AcceptedPairs(v.bind(py).borrow().to_ast(py))
            }
            Self::Degree(v) => AstAtomConstraintAst::Degree(v.bind(py).borrow().to_ast(py)),
            Self::TotalDegree(v) => {
                AstAtomConstraintAst::TotalDegree(v.bind(py).borrow().to_ast(py))
            }
            Self::RingDegree(v) => AstAtomConstraintAst::RingDegree(v.bind(py).borrow().to_ast(py)),
            Self::RingValence(v) => {
                AstAtomConstraintAst::RingValence(v.bind(py).borrow().to_ast(py))
            }
            Self::TotalHydrogens(v) => {
                AstAtomConstraintAst::TotalHydrogens(v.bind(py).borrow().to_ast(py))
            }
            Self::RingMembership(m) => {
                AstAtomConstraintAst::RingMembership(m.bind(py).borrow().to_ast(py))
            }
            Self::TetrahedralStereo(c) => {
                AstAtomConstraintAst::TetrahedralStereo(c.bind(py).borrow().to_ast(py))
            }
        }
    }
}

/// The argument to `update`: another constraint container (value or live view) or
/// an iterable of `AtomConstraintAst` (each `set`, last-wins).
#[derive(FromPyObject)]
enum AtomConstraintsUpdate {
    Container(Py<AtomConstraintsAst>),
    View(Py<AtomConstraintsView>),
    Entries(Vec<Py<AtomConstraintAst>>),
}

impl AtomConstraintsUpdate {
    /// Read every Python object into owned data — no write target is touched. Callers
    /// resolve *before* taking the write borrow so a view (or container) that aliases the
    /// same atom is read while nothing is borrowed (otherwise
    /// `atom.constraints.update(atom.constraints)` self-aliases into a double-borrow panic).
    fn resolve(&self, py: Python<'_>) -> PyResult<ResolvedAtomConstraintsUpdate> {
        Ok(match self {
            AtomConstraintsUpdate::Container(c) => {
                ResolvedAtomConstraintsUpdate::Overlay(c.bind(py).borrow().inner().clone())
            }
            AtomConstraintsUpdate::View(v) => ResolvedAtomConstraintsUpdate::Overlay(
                v.bind(py).borrow().read(py, |cs| Ok(cs.clone()))?,
            ),
            AtomConstraintsUpdate::Entries(entries) => ResolvedAtomConstraintsUpdate::Entries(
                entries
                    .iter()
                    .map(|entry| entry.bind(py).borrow().to_ast(py))
                    .collect(),
            ),
        })
    }
}

/// A `AtomConstraintsUpdate` with all Python-object reads already done, so it can be applied
/// under a write borrow without re-entering Python.
enum ResolvedAtomConstraintsUpdate {
    /// A whole container (from another container or a live view): overlaid via `update`
    /// (last-wins per key; undetermined entries remove).
    Overlay(AstAtomConstraintsAst),
    /// Loose entries: `set` each (last-wins; undetermined entries stored, not removed).
    Entries(Vec<AstAtomConstraintAst>),
}

impl ResolvedAtomConstraintsUpdate {
    /// Overlay onto `target` in place. No Python reads.
    fn apply(self, target: &mut AstAtomConstraintsAst) {
        match self {
            ResolvedAtomConstraintsUpdate::Overlay(overlay) => target.update(&overlay),
            ResolvedAtomConstraintsUpdate::Entries(entries) => {
                for entry in entries {
                    target.set(entry);
                }
            }
        }
    }
}

/// A whole-container argument that snapshots either a value container or a live
/// view — for the atom `constraints` setter, which accepts either.
#[derive(FromPyObject)]
pub(crate) enum AtomConstraintsArg {
    Container(Py<AtomConstraintsAst>),
    View(Py<AtomConstraintsView>),
}

impl AtomConstraintsArg {
    pub(crate) fn to_ast(&self, py: Python<'_>) -> PyResult<AstAtomConstraintsAst> {
        match self {
            AtomConstraintsArg::Container(c) => Ok(c.bind(py).borrow().inner().clone()),
            AtomConstraintsArg::View(v) => v.bind(py).borrow().read(py, |cs| Ok(cs.clone())),
        }
    }
}

/// The atom-scope constraints on an atom, in kind-sorted order. Mutable, hence
/// value-equal but unhashable (matching `AtomAst`).
#[pyclass(eq)]
#[derive(PartialEq)]
pub struct AtomConstraintsAst(AstAtomConstraintsAst);

#[pymethods]
impl AtomConstraintsAst {
    /// Build from a sequence of constraints (kind-sorted; a unique kind replaces
    /// an earlier one, ring memberships accumulate per scope).
    #[new]
    fn new(py: Python<'_>, entries: Vec<Py<AtomConstraintAst>>) -> Self {
        let mut constraints = AstAtomConstraintsAst::new();
        constraints.extend(
            entries
                .into_iter()
                .map(|entry| entry.bind(py).borrow().to_ast(py)),
        );
        AtomConstraintsAst(constraints)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let mut parts = Vec::with_capacity(self.0.len());
        for entry in self.0.iter() {
            let mirror = into_py_variant(py, AtomConstraintAst::from_ast(py, entry)?)?;
            parts.push(mirror.bind(py).as_any().repr()?.extract::<String>()?);
        }
        Ok(format!("AtomConstraintsAst([{}])", parts.join(", ")))
    }

    /// Insert `c`, replacing any existing entry of the same key (last-wins).
    fn set(&mut self, py: Python<'_>, c: Py<AtomConstraintAst>) {
        self.0.set(c.bind(py).borrow().to_ast(py));
    }

    /// Remove the entry with the given key, returning it if present (dict `pop`).
    fn pop(
        &mut self,
        py: Python<'_>,
        key: Py<AtomConstraintKey>,
    ) -> PyResult<Option<AtomConstraintAst>> {
        self.0
            .remove(key.bind(py).borrow().to_ast(py))
            .map(|c| AtomConstraintAst::from_ast(py, &c))
            .transpose()
    }

    /// Overlay `other` onto self in place — another container, a live view, or an
    /// iterable of `AtomConstraintAst` (last-wins per key; undetermined entries remove).
    /// Takes `slf` by handle so `other` is fully read *before* the write borrow —
    /// `cs.update(cs)` on the same container is then a no-op, not a double-borrow panic.
    fn update(slf: Py<Self>, py: Python<'_>, other: AtomConstraintsUpdate) -> PyResult<()> {
        let resolved = other.resolve(py)?;
        resolved.apply(&mut slf.borrow_mut(py).0);
        Ok(())
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    fn __iter__(&self, py: Python<'_>) -> PyResult<AtomConstraintKeyIter> {
        atom_constraint_keys(py, &self.0)
    }

    /// The constraint keys, in canonical order.
    fn keys(&self, py: Python<'_>) -> PyResult<AtomConstraintKeyIter> {
        atom_constraint_keys(py, &self.0)
    }

    /// The constraints, in canonical order.
    fn values(&self, py: Python<'_>) -> PyResult<AtomConstraintIter> {
        atom_constraints_iter(py, &self.0)
    }

    /// The `(key, constraint)` pairs, in canonical order.
    fn items(&self, py: Python<'_>) -> PyResult<AtomConstraintItemsIter> {
        atom_constraint_items(py, &self.0)
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    fn get(
        &self,
        py: Python<'_>,
        key: Py<AtomConstraintKey>,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        match self.0.get(key.bind(py).borrow().to_ast(py)) {
            Some(constraint) => {
                Ok(into_py_variant(py, AtomConstraintAst::from_ast(py, constraint)?)?.into_any())
            }
            None => Ok(default.unwrap_or_else(|| py.None())),
        }
    }

    /// The constraint with the given key; raises `KeyError` if absent.
    fn __getitem__(
        &self,
        py: Python<'_>,
        key: Py<AtomConstraintKey>,
    ) -> PyResult<AtomConstraintAst> {
        match self.0.get(key.bind(py).borrow().to_ast(py)) {
            Some(constraint) => AtomConstraintAst::from_ast(py, constraint),
            None => Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            )),
        }
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    fn __delitem__(&mut self, py: Python<'_>, key: Py<AtomConstraintKey>) -> PyResult<()> {
        if self.0.remove(key.bind(py).borrow().to_ast(py)).is_some() {
            Ok(())
        } else {
            Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            ))
        }
    }

    fn __contains__(&self, py: Python<'_>, key: Py<AtomConstraintKey>) -> bool {
        self.0.contains(key.bind(py).borrow().to_ast(py))
    }

    /// The valence value, or `None`.
    #[getter]
    fn valence(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .valence()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    #[setter]
    fn set_valence(&mut self, py: Python<'_>, value: ValueArg) {
        self.0.set(AstAtomConstraintAst::valence(value.to_ast(py)));
    }

    /// The donated-pairs value, or `None`.
    #[getter]
    fn donated_pairs(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .donated_pairs()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    #[setter]
    fn set_donated_pairs(&mut self, py: Python<'_>, value: ValueArg) {
        self.0
            .set(AstAtomConstraintAst::donated_pairs(value.to_ast(py)));
    }

    /// The accepted-pairs value, or `None`.
    #[getter]
    fn accepted_pairs(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .accepted_pairs()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    #[setter]
    fn set_accepted_pairs(&mut self, py: Python<'_>, value: ValueArg) {
        self.0
            .set(AstAtomConstraintAst::accepted_pairs(value.to_ast(py)));
    }

    /// The aromatic-valence state, or `None`.
    #[getter]
    fn aromatic_valence(&self, py: Python<'_>) -> PyResult<Option<AromaticValenceAst>> {
        self.0
            .aromatic_valence()
            .map(|c| AromaticValenceAst::from_ast(py, c))
            .transpose()
    }

    #[setter]
    fn set_aromatic_valence(&mut self, py: Python<'_>, value: AromaticValenceArg) -> PyResult<()> {
        self.0
            .set(AstAtomConstraintAst::aromatic_valence(value.to_ast(py)?));
        Ok(())
    }

    /// The multicenter-valence state, or `None`.
    #[getter]
    fn multicenter_valence(&self, py: Python<'_>) -> PyResult<Option<MulticenterValenceAst>> {
        self.0
            .multicenter_valence()
            .map(|c| MulticenterValenceAst::from_ast(py, c))
            .transpose()
    }

    #[setter]
    fn set_multicenter_valence(
        &mut self,
        py: Python<'_>,
        value: MulticenterValenceArg,
    ) -> PyResult<()> {
        self.0
            .set(AstAtomConstraintAst::multicenter_valence(value.to_ast(py)?));
        Ok(())
    }

    /// The tetrahedral-stereo state, or `None`.
    #[getter]
    fn tetrahedral_stereo(&self, py: Python<'_>) -> PyResult<Option<TetrahedralStereoAst>> {
        self.0
            .tetrahedral_stereo()
            .map(|c| TetrahedralStereoAst::from_ast(py, c))
            .transpose()
    }

    #[setter]
    fn set_tetrahedral_stereo(
        &mut self,
        py: Python<'_>,
        value: TetrahedralStereoArg,
    ) -> PyResult<()> {
        self.0
            .set(AstAtomConstraintAst::tetrahedral_stereo(value.to_ast(py)?));
        Ok(())
    }

    /// The degree value, or `None`.
    #[getter]
    fn degree(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .degree()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    #[setter]
    fn set_degree(&mut self, py: Python<'_>, value: ValueArg) {
        self.0.set(AstAtomConstraintAst::degree(value.to_ast(py)));
    }

    /// The total-degree value, or `None`.
    #[getter]
    fn total_degree(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .total_degree()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    #[setter]
    fn set_total_degree(&mut self, py: Python<'_>, value: ValueArg) {
        self.0
            .set(AstAtomConstraintAst::total_degree(value.to_ast(py)));
    }

    /// The total-valence value, or `None`.
    #[getter]
    fn total_valence(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .total_valence()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    #[setter]
    fn set_total_valence(&mut self, py: Python<'_>, value: ValueArg) {
        self.0
            .set(AstAtomConstraintAst::total_valence(value.to_ast(py)));
    }

    /// The ring-degree value, or `None`.
    #[getter]
    fn ring_degree(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .ring_degree()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    #[setter]
    fn set_ring_degree(&mut self, py: Python<'_>, value: ValueArg) {
        self.0
            .set(AstAtomConstraintAst::ring_degree(value.to_ast(py)));
    }

    /// The ring-valence value, or `None`.
    #[getter]
    fn ring_valence(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .ring_valence()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    #[setter]
    fn set_ring_valence(&mut self, py: Python<'_>, value: ValueArg) {
        self.0
            .set(AstAtomConstraintAst::ring_valence(value.to_ast(py)));
    }

    /// The total-hydrogens value, or `None`.
    #[getter]
    fn total_hydrogens(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .total_hydrogens()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    #[setter]
    fn set_total_hydrogens(&mut self, py: Python<'_>, value: ValueArg) {
        self.0
            .set(AstAtomConstraintAst::total_hydrogens(value.to_ast(py)));
    }

    /// The all-rings membership count, or `None`.
    #[getter]
    fn ring_count(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .ring_count()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    #[setter]
    fn set_ring_count(&mut self, py: Python<'_>, value: ValueArg) {
        self.0.set(AstAtomConstraintAst::ring_membership(
            AstRingScope::All,
            value.to_ast(py),
        ));
    }

    /// The sized-ring membership counts, as a subscriptable proxy keyed by ring
    /// size: `constraints.ring_size_count[6]`, `[6] = 3`, `del [6]`.
    #[getter]
    fn ring_size_count(slf: Py<Self>) -> AtomRingSizeCounts {
        AtomRingSizeCounts {
            backing: AtomRingSizeBacking::Value(slf),
        }
    }

    /// The present constraints as a dict keyed by snake_case name; values are the
    /// inner-value mirrors. Ring memberships key by scope: `ring_count` for the
    /// all-rings scope, `ring_size_count_<n>` for a specific ring size.
    pub(crate) fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        atom_constraints_asdict(py, &self.0)
    }
}

impl AtomConstraintsAst {
    /// The wrapped AST constraints — read access for atom construction.
    pub(crate) fn inner(&self) -> &AstAtomConstraintsAst {
        &self.0
    }

    /// Mutable access to the wrapped AST constraints — for the value-backed proxy.
    pub(crate) fn inner_mut(&mut self) -> &mut AstAtomConstraintsAst {
        &mut self.0
    }

    /// Wrap AST constraints (the hold-the-value `from_inner` bridge). Test-only —
    /// in-crate construction wraps `AtomConstraintsAst(..)` directly.
    #[cfg(test)]
    pub(crate) fn from_inner(constraints: AstAtomConstraintsAst) -> Self {
        AtomConstraintsAst(constraints)
    }
}

/// Build the per-constraint iterator handle from a borrowed container.
fn atom_constraints_iter(
    py: Python<'_>,
    constraints: &AstAtomConstraintsAst,
) -> PyResult<AtomConstraintIter> {
    let entries = constraints
        .iter()
        .map(|constraint| into_py_variant(py, AtomConstraintAst::from_ast(py, constraint)?))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(AtomConstraintIter {
        entries: entries.into_iter(),
    })
}

/// Build the key iterator handle from a borrowed container (mapping-style keys).
fn atom_constraint_keys(
    py: Python<'_>,
    constraints: &AstAtomConstraintsAst,
) -> PyResult<AtomConstraintKeyIter> {
    let keys = constraints
        .iter()
        .map(|constraint| into_py_variant(py, AtomConstraintKey::from_ast(py, &constraint.key())?))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(AtomConstraintKeyIter {
        keys: keys.into_iter(),
    })
}

/// Build the item iterator handle (`(key, constraint)` pairs) from a borrowed container.
fn atom_constraint_items(
    py: Python<'_>,
    constraints: &AstAtomConstraintsAst,
) -> PyResult<AtomConstraintItemsIter> {
    let items = constraints
        .iter()
        .map(|constraint| {
            Ok((
                into_py_variant(py, AtomConstraintKey::from_ast(py, &constraint.key())?)?,
                into_py_variant(py, AtomConstraintAst::from_ast(py, constraint)?)?,
            ))
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(AtomConstraintItemsIter {
        items: items.into_iter(),
    })
}

/// The present constraints as a dict keyed by snake_case name; values are the
/// inner-value mirrors. Ring memberships key by scope: `ring_count` for the
/// all-rings scope, `ring_size_count_<n>` for a specific ring size.
pub(crate) fn atom_constraints_asdict<'py>(
    py: Python<'py>,
    constraints: &AstAtomConstraintsAst,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    for entry in constraints.iter() {
        match entry {
            AstAtomConstraintAst::Valence(v) => {
                dict.set_item("valence", ValueAst::from_ast(py, v)?)?
            }
            AstAtomConstraintAst::DonatedPairs(v) => {
                dict.set_item("donated_pairs", ValueAst::from_ast(py, v)?)?
            }
            AstAtomConstraintAst::AcceptedPairs(v) => {
                dict.set_item("accepted_pairs", ValueAst::from_ast(py, v)?)?
            }
            AstAtomConstraintAst::AromaticValence(c) => {
                dict.set_item("aromatic_valence", AromaticValenceAst::from_ast(py, c)?)?
            }
            AstAtomConstraintAst::MulticenterValence(c) => dict.set_item(
                "multicenter_valence",
                MulticenterValenceAst::from_ast(py, c)?,
            )?,
            AstAtomConstraintAst::TetrahedralStereo(c) => {
                dict.set_item("tetrahedral_stereo", TetrahedralStereoAst::from_ast(py, c)?)?
            }
            AstAtomConstraintAst::Degree(v) => {
                dict.set_item("degree", ValueAst::from_ast(py, v)?)?
            }
            AstAtomConstraintAst::TotalDegree(v) => {
                dict.set_item("total_degree", ValueAst::from_ast(py, v)?)?
            }
            AstAtomConstraintAst::TotalValence(v) => {
                dict.set_item("total_valence", ValueAst::from_ast(py, v)?)?
            }
            AstAtomConstraintAst::RingDegree(v) => {
                dict.set_item("ring_degree", ValueAst::from_ast(py, v)?)?
            }
            AstAtomConstraintAst::RingValence(v) => {
                dict.set_item("ring_valence", ValueAst::from_ast(py, v)?)?
            }
            AstAtomConstraintAst::TotalHydrogens(v) => {
                dict.set_item("total_hydrogens", ValueAst::from_ast(py, v)?)?
            }
            AstAtomConstraintAst::RingMembership(m) => {
                let key = match m.scope {
                    AstRingScope::All => "ring_count".to_string(),
                    AstRingScope::Size(size) => format!("ring_size_count_{size}"),
                };
                dict.set_item(key, ValueAst::from_ast(py, &m.count)?)?
            }
        }
    }
    Ok(dict)
}

/// What an `AtomConstraintsView` writes through to: an atom within a molecule
/// (by index) or a standalone `AtomAst`.
pub(crate) enum AtomConstraintsBacking {
    Molecule {
        owner: Py<MoleculeAst>,
        id: AstAtomId,
    },
    Atom(Py<AtomAst>),
}

/// A live handle onto one atom's constraints, backed by either a molecule-atom or
/// a standalone `AtomAst`. Reads borrow the atom's constraints and read only the
/// item they need (no whole-container clone); mutators write through to the atom in
/// place, without a clone-and-writeback.
#[pyclass]
pub struct AtomConstraintsView {
    pub(crate) backing: AtomConstraintsBacking,
}

impl AtomConstraintsView {
    /// Borrow the backing atom's constraints and read one item through `f` — no clone.
    fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&AstAtomConstraintsAst) -> PyResult<R>,
    ) -> PyResult<R> {
        match &self.backing {
            AtomConstraintsBacking::Molecule { owner, id } => {
                let molecule = owner.bind(py).borrow();
                let view = molecule
                    .inner()
                    .atoms()
                    .get(*id)
                    .ok_or_else(|| PyIndexError::new_err("atom id out of range"))?;
                f(&view.ast.constraints)
            }
            AtomConstraintsBacking::Atom(atom) => {
                let atom = atom.bind(py).borrow();
                f(&atom.inner().constraints)
            }
        }
    }

    /// Mutate the backing atom's constraints in place through `f`.
    fn with_mut<R>(&self, py: Python<'_>, f: impl FnOnce(&mut AstAtomConstraintsAst) -> R) -> R {
        match &self.backing {
            AtomConstraintsBacking::Molecule { owner, id } => f(&mut owner
                .borrow_mut(py)
                .inner_mut()
                .atom_mut(*id)
                .ast
                .constraints),
            AtomConstraintsBacking::Atom(atom) => {
                f(&mut atom.borrow_mut(py).inner_mut().constraints)
            }
        }
    }

    /// Set one constraint on the backing atom in place (last-wins per key).
    fn set_ast(&self, py: Python<'_>, constraint: AstAtomConstraintAst) {
        self.with_mut(py, |cs| cs.set(constraint));
    }

    /// Remove one key from the backing atom in place, returning the removed entry.
    fn remove_ast(
        &self,
        py: Python<'_>,
        key: AstAtomConstraintKey,
    ) -> Option<AstAtomConstraintAst> {
        self.with_mut(py, |cs| cs.remove(key))
    }
}

#[pymethods]
impl AtomConstraintsView {
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let count = self.read(py, |cs| Ok(cs.len()))?;
        Ok(format!("AtomConstraintsView({count} entries)"))
    }

    /// Insert `c` on the atom in place, replacing any existing entry of the same
    /// key (last-wins).
    fn set(&self, py: Python<'_>, c: Py<AtomConstraintAst>) {
        self.set_ast(py, c.bind(py).borrow().to_ast(py));
    }

    /// Remove the entry with the given key from the atom in place, returning it if
    /// present (dict `pop`).
    fn pop(
        &self,
        py: Python<'_>,
        key: Py<AtomConstraintKey>,
    ) -> PyResult<Option<AtomConstraintAst>> {
        self.remove_ast(py, key.bind(py).borrow().to_ast(py))
            .map(|c| AtomConstraintAst::from_ast(py, &c))
            .transpose()
    }

    /// Remove the entry with the given key; raises `KeyError` if absent.
    fn __delitem__(&self, py: Python<'_>, key: Py<AtomConstraintKey>) -> PyResult<()> {
        if self
            .remove_ast(py, key.bind(py).borrow().to_ast(py))
            .is_some()
        {
            Ok(())
        } else {
            Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            ))
        }
    }

    /// Overlay `other` onto the atom's constraints in place — another container, a live
    /// view, or an iterable of `AtomConstraintAst` (last-wins per key; undetermined
    /// entries remove). Resolves `other` to owned data *before* the write borrow, so a
    /// view aliasing the same atom is not a double-borrow panic.
    fn update(&self, py: Python<'_>, other: AtomConstraintsUpdate) -> PyResult<()> {
        let resolved = other.resolve(py)?;
        self.with_mut(py, |cs| resolved.apply(cs));
        Ok(())
    }

    fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        self.read(py, |cs| Ok(cs.len()))
    }

    /// Iterate the constraint keys (mapping-style, canonical order).
    fn __iter__(&self, py: Python<'_>) -> PyResult<AtomConstraintKeyIter> {
        self.read(py, |cs| atom_constraint_keys(py, cs))
    }

    /// The constraint keys, in canonical order.
    fn keys(&self, py: Python<'_>) -> PyResult<AtomConstraintKeyIter> {
        self.read(py, |cs| atom_constraint_keys(py, cs))
    }

    /// The constraints, in canonical order.
    fn values(&self, py: Python<'_>) -> PyResult<AtomConstraintIter> {
        self.read(py, |cs| atom_constraints_iter(py, cs))
    }

    /// The `(key, constraint)` pairs, in canonical order.
    fn items(&self, py: Python<'_>) -> PyResult<AtomConstraintItemsIter> {
        self.read(py, |cs| atom_constraint_items(py, cs))
    }

    /// The constraint with the given key, or `default` (`None`) if absent.
    #[pyo3(signature = (key, default=None))]
    fn get(
        &self,
        py: Python<'_>,
        key: Py<AtomConstraintKey>,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        let key = key.bind(py).borrow().to_ast(py);
        let found = self.read(py, |cs| {
            cs.get(key)
                .map(|constraint| AtomConstraintAst::from_ast(py, constraint))
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
        key: Py<AtomConstraintKey>,
    ) -> PyResult<AtomConstraintAst> {
        let ast_key = key.bind(py).borrow().to_ast(py);
        let found = self.read(py, |cs| {
            cs.get(ast_key)
                .map(|constraint| AtomConstraintAst::from_ast(py, constraint))
                .transpose()
        })?;
        match found {
            Some(constraint) => Ok(constraint),
            None => Err(PyKeyError::new_err(
                key.bind(py).as_any().repr()?.extract::<String>()?,
            )),
        }
    }

    fn __contains__(&self, py: Python<'_>, key: Py<AtomConstraintKey>) -> PyResult<bool> {
        let key = key.bind(py).borrow().to_ast(py);
        self.read(py, |cs| Ok(cs.contains(key)))
    }

    /// The valence value, or `None`.
    #[getter]
    fn valence(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.valence().map(|v| ValueAst::from_ast(py, v)).transpose()
        })
    }

    #[setter]
    fn set_valence(&self, py: Python<'_>, value: ValueArg) {
        self.set_ast(py, AstAtomConstraintAst::valence(value.to_ast(py)));
    }

    /// The donated-pairs value, or `None`.
    #[getter]
    fn donated_pairs(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.donated_pairs()
                .map(|v| ValueAst::from_ast(py, v))
                .transpose()
        })
    }

    #[setter]
    fn set_donated_pairs(&self, py: Python<'_>, value: ValueArg) {
        self.set_ast(py, AstAtomConstraintAst::donated_pairs(value.to_ast(py)));
    }

    /// The accepted-pairs value, or `None`.
    #[getter]
    fn accepted_pairs(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.accepted_pairs()
                .map(|v| ValueAst::from_ast(py, v))
                .transpose()
        })
    }

    #[setter]
    fn set_accepted_pairs(&self, py: Python<'_>, value: ValueArg) {
        self.set_ast(py, AstAtomConstraintAst::accepted_pairs(value.to_ast(py)));
    }

    /// The aromatic-valence state, or `None`.
    #[getter]
    fn aromatic_valence(&self, py: Python<'_>) -> PyResult<Option<AromaticValenceAst>> {
        self.read(py, |cs| {
            cs.aromatic_valence()
                .map(|c| AromaticValenceAst::from_ast(py, c))
                .transpose()
        })
    }

    #[setter]
    fn set_aromatic_valence(&self, py: Python<'_>, value: AromaticValenceArg) -> PyResult<()> {
        self.set_ast(
            py,
            AstAtomConstraintAst::aromatic_valence(value.to_ast(py)?),
        );
        Ok(())
    }

    /// The multicenter-valence state, or `None`.
    #[getter]
    fn multicenter_valence(&self, py: Python<'_>) -> PyResult<Option<MulticenterValenceAst>> {
        self.read(py, |cs| {
            cs.multicenter_valence()
                .map(|c| MulticenterValenceAst::from_ast(py, c))
                .transpose()
        })
    }

    #[setter]
    fn set_multicenter_valence(
        &self,
        py: Python<'_>,
        value: MulticenterValenceArg,
    ) -> PyResult<()> {
        self.set_ast(
            py,
            AstAtomConstraintAst::multicenter_valence(value.to_ast(py)?),
        );
        Ok(())
    }

    /// The tetrahedral-stereo state, or `None`.
    #[getter]
    fn tetrahedral_stereo(&self, py: Python<'_>) -> PyResult<Option<TetrahedralStereoAst>> {
        self.read(py, |cs| {
            cs.tetrahedral_stereo()
                .map(|c| TetrahedralStereoAst::from_ast(py, c))
                .transpose()
        })
    }

    #[setter]
    fn set_tetrahedral_stereo(&self, py: Python<'_>, value: TetrahedralStereoArg) -> PyResult<()> {
        self.set_ast(
            py,
            AstAtomConstraintAst::tetrahedral_stereo(value.to_ast(py)?),
        );
        Ok(())
    }

    /// The degree value, or `None`.
    #[getter]
    fn degree(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.degree().map(|v| ValueAst::from_ast(py, v)).transpose()
        })
    }

    #[setter]
    fn set_degree(&self, py: Python<'_>, value: ValueArg) {
        self.set_ast(py, AstAtomConstraintAst::degree(value.to_ast(py)));
    }

    /// The total-degree value, or `None`.
    #[getter]
    fn total_degree(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.total_degree()
                .map(|v| ValueAst::from_ast(py, v))
                .transpose()
        })
    }

    #[setter]
    fn set_total_degree(&self, py: Python<'_>, value: ValueArg) {
        self.set_ast(py, AstAtomConstraintAst::total_degree(value.to_ast(py)));
    }

    /// The total-valence value, or `None`.
    #[getter]
    fn total_valence(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.total_valence()
                .map(|v| ValueAst::from_ast(py, v))
                .transpose()
        })
    }

    #[setter]
    fn set_total_valence(&self, py: Python<'_>, value: ValueArg) {
        self.set_ast(py, AstAtomConstraintAst::total_valence(value.to_ast(py)));
    }

    /// The ring-degree value, or `None`.
    #[getter]
    fn ring_degree(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.ring_degree()
                .map(|v| ValueAst::from_ast(py, v))
                .transpose()
        })
    }

    #[setter]
    fn set_ring_degree(&self, py: Python<'_>, value: ValueArg) {
        self.set_ast(py, AstAtomConstraintAst::ring_degree(value.to_ast(py)));
    }

    /// The ring-valence value, or `None`.
    #[getter]
    fn ring_valence(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.ring_valence()
                .map(|v| ValueAst::from_ast(py, v))
                .transpose()
        })
    }

    #[setter]
    fn set_ring_valence(&self, py: Python<'_>, value: ValueArg) {
        self.set_ast(py, AstAtomConstraintAst::ring_valence(value.to_ast(py)));
    }

    /// The total-hydrogens value, or `None`.
    #[getter]
    fn total_hydrogens(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.total_hydrogens()
                .map(|v| ValueAst::from_ast(py, v))
                .transpose()
        })
    }

    #[setter]
    fn set_total_hydrogens(&self, py: Python<'_>, value: ValueArg) {
        self.set_ast(py, AstAtomConstraintAst::total_hydrogens(value.to_ast(py)));
    }

    /// The all-rings membership count, or `None`.
    #[getter]
    fn ring_count(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.ring_count()
                .map(|v| ValueAst::from_ast(py, v))
                .transpose()
        })
    }

    #[setter]
    fn set_ring_count(&self, py: Python<'_>, value: ValueArg) {
        self.set_ast(
            py,
            AstAtomConstraintAst::ring_membership(AstRingScope::All, value.to_ast(py)),
        );
    }

    /// The sized-ring membership counts, as a subscriptable proxy keyed by ring
    /// size: `constraints.ring_size_count[6]`, `[6] = 3`, `del [6]`.
    #[getter]
    fn ring_size_count(&self, py: Python<'_>) -> AtomRingSizeCounts {
        let backing = match &self.backing {
            AtomConstraintsBacking::Molecule { owner, id } => AtomRingSizeBacking::Molecule {
                owner: owner.clone_ref(py),
                id: *id,
            },
            AtomConstraintsBacking::Atom(atom) => AtomRingSizeBacking::Atom(atom.clone_ref(py)),
        };
        AtomRingSizeCounts { backing }
    }

    /// The present constraints as a dict keyed by snake_case name.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        self.read(py, |cs| atom_constraints_asdict(py, cs))
    }
}

/// What a `AtomRingSizeCounts` proxy reads/writes through to: an atom within a molecule,
/// a standalone `AtomAst`, or a standalone `AtomConstraintsAst` value.
pub(crate) enum AtomRingSizeBacking {
    Molecule {
        owner: Py<MoleculeAst>,
        id: AstAtomId,
    },
    Atom(Py<AtomAst>),
    Value(Py<AtomConstraintsAst>),
}

/// A subscriptable proxy over the sized-ring membership counts of an atom, keyed by
/// ring size: `proxy[size]` reads, `proxy[size] = count` sets, `del proxy[size]`
/// removes. Backs onto whichever container produced it (dual-backing, like
/// `AtomConstraintsView`).
#[pyclass]
pub struct AtomRingSizeCounts {
    backing: AtomRingSizeBacking,
}

impl AtomRingSizeCounts {
    /// Borrow the backing constraints and read through `f` — no clone.
    fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&AstAtomConstraintsAst) -> PyResult<R>,
    ) -> PyResult<R> {
        match &self.backing {
            AtomRingSizeBacking::Molecule { owner, id } => {
                let molecule = owner.bind(py).borrow();
                let view = molecule
                    .inner()
                    .atoms()
                    .get(*id)
                    .ok_or_else(|| PyIndexError::new_err("atom id out of range"))?;
                f(&view.ast.constraints)
            }
            AtomRingSizeBacking::Atom(atom) => f(&atom.bind(py).borrow().inner().constraints),
            AtomRingSizeBacking::Value(value) => f(value.bind(py).borrow().inner()),
        }
    }

    /// Mutate the backing constraints in place through `f`.
    fn write(&self, py: Python<'_>, f: impl FnOnce(&mut AstAtomConstraintsAst)) {
        match &self.backing {
            AtomRingSizeBacking::Molecule { owner, id } => f(&mut owner
                .borrow_mut(py)
                .inner_mut()
                .atom_mut(*id)
                .ast
                .constraints),
            AtomRingSizeBacking::Atom(atom) => f(&mut atom.borrow_mut(py).inner_mut().constraints),
            AtomRingSizeBacking::Value(value) => f(value.borrow_mut(py).inner_mut()),
        }
    }
}

#[pymethods]
impl AtomRingSizeCounts {
    /// The membership count for rings of `size`, or `None`.
    fn __getitem__(&self, py: Python<'_>, size: u8) -> PyResult<Option<ValueAst>> {
        self.read(py, |cs| {
            cs.ring_size_count(size)
                .map(|v| ValueAst::from_ast(py, v))
                .transpose()
        })
    }

    /// The number of distinct ring sizes with a membership constraint.
    fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        self.read(py, |cs| Ok(ring_sizes(cs).count()))
    }

    fn __contains__(&self, py: Python<'_>, size: u8) -> PyResult<bool> {
        self.read(py, |cs| Ok(cs.ring_size_count(size).is_some()))
    }

    /// Iterate the present ring sizes (as ints).
    fn __iter__(&self, py: Python<'_>) -> PyResult<AtomRingSizeIter> {
        let sizes = self.read(py, |cs| Ok(ring_sizes(cs).collect::<Vec<u8>>()))?;
        Ok(AtomRingSizeIter {
            sizes: sizes.into_iter(),
        })
    }

    /// Set the membership count for rings of `size` in place.
    fn __setitem__(&self, py: Python<'_>, size: u8, count: ValueArg) {
        let constraint =
            AstAtomConstraintAst::ring_membership(AstRingScope::Size(size), count.to_ast(py));
        self.write(py, |cs| cs.set(constraint));
    }

    /// Remove the sized-ring membership for `size` in place.
    fn __delitem__(&self, py: Python<'_>, size: u8) {
        self.write(py, |cs| {
            cs.remove(AstAtomConstraintKey::RingMembership(AstRingScope::Size(
                size,
            )));
        });
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        self.read(py, |cs| {
            let mut parts = Vec::new();
            for entry in cs.iter() {
                if let AstAtomConstraintAst::RingMembership(m) = entry {
                    if let AstRingScope::Size(size) = m.scope {
                        let count = into_py_variant(py, ValueAst::from_ast(py, &m.count)?)?;
                        parts.push(format!(
                            "{size}: {}",
                            count.bind(py).as_any().repr()?.extract::<String>()?
                        ));
                    }
                }
            }
            Ok(format!("AtomRingSizeCounts({{{}}})", parts.join(", ")))
        })
    }
}

/// The ring sizes with a membership constraint, in kind-sorted order.
fn ring_sizes(constraints: &AstAtomConstraintsAst) -> impl Iterator<Item = u8> + '_ {
    constraints.iter().filter_map(|entry| match entry {
        AstAtomConstraintAst::RingMembership(m) => match m.scope {
            AstRingScope::Size(size) => Some(size),
            AstRingScope::All => None,
        },
        _ => None,
    })
}

#[pyclass]
struct AtomRingSizeIter {
    sizes: IntoIter<u8>,
}

#[pymethods]
impl AtomRingSizeIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<u8> {
        self.sizes.next()
    }
}

#[pyclass]
struct AtomConstraintIter {
    entries: IntoIter<Py<AtomConstraintAst>>,
}

#[pymethods]
impl AtomConstraintIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<Py<AtomConstraintAst>> {
        self.entries.next()
    }
}

#[pyclass]
struct AtomConstraintKeyIter {
    keys: IntoIter<Py<AtomConstraintKey>>,
}

#[pymethods]
impl AtomConstraintKeyIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<Py<AtomConstraintKey>> {
        self.keys.next()
    }
}

#[pyclass]
struct AtomConstraintItemsIter {
    items: IntoIter<(Py<AtomConstraintKey>, Py<AtomConstraintAst>)>,
}

#[pymethods]
impl AtomConstraintItemsIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<(Py<AtomConstraintKey>, Py<AtomConstraintAst>)> {
        self.items.next()
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::{
        AromaticValenceAst as AstAromaticValenceAst, AtomConstraintAst as AstAtomConstraintAst,
        MemOp as AstMemOp, StereoCosetAst as AstStereoCosetAst,
        TetrahedralStereoAst as AstTetrahedralStereoAst, ValueAst as AstValueAst,
    };

    use super::*;
    use crate::stereo::TetrahedralStereo;

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
        let molecule = AstMoleculeAst::from_parts(AstMoleculeParts {
            atoms: vec![
                AstAtomAst::from_element(ChemElement::C),
                AstAtomAst::from_element(ChemElement::O),
            ],
            ..Default::default()
        });
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
            let nitrogen = Py::new(
                py,
                AtomAst::from_inner(AstAtomAst::from_element(ChemElement::N)),
            )
            .unwrap();
            views
                .__setitem__(py, 0, nitrogen.bind(py).borrow())
                .unwrap();
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
            let nitrogen = Py::new(
                py,
                AtomAst::from_inner(AstAtomAst::from_element(ChemElement::N)),
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
                    backing: AtomConstraintsBacking::Atom(src),
                },
            )
            .unwrap();
            let dst = Py::new(
                py,
                AtomAst::from_inner(AstAtomAst::from_element(ChemElement::N)),
            )
            .unwrap();
            AtomAst::set_constraints(dst.clone_ref(py), py, AtomConstraintsArg::View(view))
                .unwrap();
            assert_eq!(
                dst.bind(py)
                    .borrow()
                    .inner()
                    .constraints
                    .valence()
                    .unwrap()
                    .clone(),
                AstValueAst::Lit(4)
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
                    AstAtomAst::from_element(ChemElement::C)
                        .with_constraint(AstAtomConstraintAst::valence(4)),
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
            AtomAst::set_constraints(atom.clone_ref(py), py, AtomConstraintsArg::View(own_view))
                .unwrap();
            assert_eq!(
                atom.bind(py)
                    .borrow()
                    .inner()
                    .constraints
                    .valence()
                    .unwrap()
                    .clone(),
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
                AtomConstraintsBacking::Molecule { id, .. } => assert_eq!(id, AstAtomId(1)),
                _ => panic!("expected molecule-backed view"),
            }
        });
    }

    #[rstest]
    #[case(AstAtomConstraintAst::valence(4))]
    #[case(AstAtomConstraintAst::aromatic_valence(AstAromaticValenceAst::aromatic(1)))]
    #[case(AstAtomConstraintAst::ring_membership(AstRingScope::All, 2))]
    #[case(AstAtomConstraintAst::tetrahedral_stereo(AstTetrahedralStereoAst::not_stereo()))]
    fn test_atom_constraint_roundtrip(#[case] ast: AstAtomConstraintAst) {
        Python::attach(|py| {
            assert_eq!(
                AtomConstraintAst::from_ast(py, &ast).unwrap().to_ast(py),
                ast
            );
        });
    }

    #[rstest]
    fn test_atom_constraints_len_contains() {
        Python::attach(|py| {
            let valence = into_py_variant(
                py,
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::valence(4)).unwrap(),
            )
            .unwrap();
            let degree = into_py_variant(
                py,
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::degree(3)).unwrap(),
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
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::valence(4)).unwrap(),
            )
            .unwrap();
            let degree = into_py_variant(
                py,
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::degree(3)).unwrap(),
            )
            .unwrap();
            let constraints = AtomConstraintsAst::new(py, vec![valence, degree]);

            let mut keys = constraints.__iter__(py).unwrap();
            assert_eq!(
                keys.__next__().unwrap().bind(py).borrow().to_ast(py),
                AstAtomConstraintKey::Valence
            );
            assert_eq!(
                keys.__next__().unwrap().bind(py).borrow().to_ast(py),
                AstAtomConstraintKey::Degree
            );
            assert!(keys.__next__().is_none());

            let mut values = constraints.values(py).unwrap();
            assert_eq!(
                values.__next__().unwrap().bind(py).borrow().to_ast(py),
                AstAtomConstraintAst::valence(4)
            );

            let mut items = constraints.items(py).unwrap();
            let (key, value) = items.__next__().unwrap();
            assert_eq!(
                key.bind(py).borrow().to_ast(py),
                AstAtomConstraintKey::Valence
            );
            assert_eq!(
                value.bind(py).borrow().to_ast(py),
                AstAtomConstraintAst::valence(4)
            );
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_get() {
        Python::attach(|py| {
            let valence = into_py_variant(
                py,
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::valence(4)).unwrap(),
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
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::valence(4)).unwrap(),
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
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::valence(4)).unwrap(),
            )
            .unwrap();
            let degree = into_py_variant(
                py,
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::degree(3)).unwrap(),
            )
            .unwrap();
            let constraints = AtomConstraintsAst::new(py, vec![valence, degree]);
            assert_eq!(
                constraints.valence(py).unwrap().unwrap().to_ast(py),
                AstValueAst::Lit(4)
            );
            assert_eq!(
                constraints.degree(py).unwrap().unwrap().to_ast(py),
                AstValueAst::Lit(3)
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
                AtomConstraintAst::from_ast(
                    py,
                    &AstAtomConstraintAst::ring_membership(AstRingScope::Size(6), 1),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = Py::new(py, AtomConstraintsAst::new(py, vec![membership])).unwrap();
            let proxy = AtomConstraintsAst::ring_size_count(constraints.clone_ref(py));
            assert_eq!(
                proxy.__getitem__(py, 6).unwrap().unwrap().to_ast(py),
                AstValueAst::Lit(1)
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
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::valence(4)).unwrap(),
            )
            .unwrap();
            constraints.set(py, valence);
            assert_eq!(constraints.__len__(), 1);
            assert_eq!(
                constraints.valence(py).unwrap().unwrap().to_ast(py),
                AstValueAst::Lit(4)
            );
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_pop() {
        Python::attach(|py| {
            let valence = into_py_variant(
                py,
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::valence(4)).unwrap(),
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
                    assert_eq!(v.bind(py).borrow().to_ast(py), AstValueAst::Lit(4))
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
            let mut other = AstAtomConstraintsAst::new();
            other.set(AstAtomConstraintAst::valence(4));
            other.set(AstAtomConstraintAst::degree(3));
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
                c.valence(py).unwrap().unwrap().to_ast(py),
                AstValueAst::Lit(4)
            );
            assert_eq!(
                c.degree(py).unwrap().unwrap().to_ast(py),
                AstValueAst::Lit(3)
            );
        });
    }

    #[rstest]
    fn test_atom_constraints_view_set() {
        Python::attach(|py| {
            let owner = Py::new(
                py,
                MoleculeAst::from_inner(AstMoleculeAst::from_parts(AstMoleculeParts {
                    atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                    ..Default::default()
                })),
            )
            .unwrap();
            let view = AtomConstraintsView {
                backing: AtomConstraintsBacking::Molecule {
                    owner: owner.clone_ref(py),
                    id: AstAtomId(0),
                },
            };
            let valence = into_py_variant(
                py,
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::valence(4)).unwrap(),
            )
            .unwrap();
            view.set(py, valence);
            let fresh = AtomConstraintsView {
                backing: AtomConstraintsBacking::Molecule {
                    owner,
                    id: AstAtomId(0),
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
                    assert_eq!(v.bind(py).borrow().to_ast(py), AstValueAst::Lit(4))
                }
                _ => panic!("expected Valence(Lit(4))"),
            }
        });
    }

    #[rstest]
    fn test_atom_constraints_view_pop() {
        Python::attach(|py| {
            let atom = AstAtomAst::from_element(ChemElement::C)
                .with_constraint(AstAtomConstraintAst::valence(4));
            let owner = Py::new(
                py,
                MoleculeAst::from_inner(AstMoleculeAst::from_parts(AstMoleculeParts {
                    atoms: vec![atom],
                    ..Default::default()
                })),
            )
            .unwrap();
            let view = AtomConstraintsView {
                backing: AtomConstraintsBacking::Molecule {
                    owner: owner.clone_ref(py),
                    id: AstAtomId(0),
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
                    assert_eq!(v.bind(py).borrow().to_ast(py), AstValueAst::Lit(4))
                }
                _ => panic!("expected removed Valence(Lit(4))"),
            }
            let fresh = AtomConstraintsView {
                backing: AtomConstraintsBacking::Molecule {
                    owner,
                    id: AstAtomId(0),
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
                MoleculeAst::from_inner(AstMoleculeAst::from_parts(AstMoleculeParts {
                    atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                    ..Default::default()
                })),
            )
            .unwrap();
            let view = AtomConstraintsView {
                backing: AtomConstraintsBacking::Molecule {
                    owner: owner.clone_ref(py),
                    id: AstAtomId(0),
                },
            };
            let mut other = AstAtomConstraintsAst::new();
            other.set(AstAtomConstraintAst::valence(4));
            other.set(AstAtomConstraintAst::degree(3));
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
                    id: AstAtomId(0),
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
                AtomAst::from_inner(AstAtomAst::from_element(ChemElement::C)),
            )
            .unwrap();
            let view = AtomConstraintsView {
                backing: AtomConstraintsBacking::Atom(atom.clone_ref(py)),
            };
            let valence = into_py_variant(
                py,
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::valence(4)).unwrap(),
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
                    assert_eq!(v.bind(py).borrow().to_ast(py), AstValueAst::Lit(4))
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
                    AstAtomAst::from_element(ChemElement::C)
                        .with_constraint(AstAtomConstraintAst::valence(4)),
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
                    assert_eq!(v.bind(py).borrow().to_ast(py), AstValueAst::Lit(4))
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
                AtomAst::from_inner(AstAtomAst::from_element(ChemElement::C)),
            )
            .unwrap();
            let view = AtomConstraintsView {
                backing: AtomConstraintsBacking::Atom(atom.clone_ref(py)),
            };
            let mut other = AstAtomConstraintsAst::new();
            other.set(AstAtomConstraintAst::valence(4));
            other.set(AstAtomConstraintAst::degree(3));
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
            constraints.set_valence(py, ValueArg::Lit(4));
            assert_eq!(
                constraints.valence(py).unwrap().unwrap().to_ast(py),
                AstValueAst::Lit(4)
            );
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_set_ring_count() {
        Python::attach(|py| {
            let mut constraints = AtomConstraintsAst::new(py, vec![]);
            constraints.set_ring_count(py, ValueArg::Lit(2));
            assert_eq!(
                constraints.ring_count(py).unwrap().unwrap().to_ast(py),
                AstValueAst::Lit(2)
            );
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_set_aromatic_valence() {
        Python::attach(|py| {
            let mut constraints = AtomConstraintsAst::new(py, vec![]);
            constraints
                .set_aromatic_valence(py, AromaticValenceArg::Value(ValueArg::Lit(1)))
                .unwrap();
            match constraints.aromatic_valence(py).unwrap().unwrap() {
                AromaticValenceAst::Aromatic(v) => assert_eq!(v.to_ast(py), AstValueAst::Lit(1)),
                _ => panic!("expected Aromatic"),
            }
            constraints
                .set_aromatic_valence(py, AromaticValenceArg::Flag(false))
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
                .set_aromatic_valence(py, AromaticValenceArg::Flag(true))
                .is_err());
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_set_tetrahedral_stereo() {
        Python::attach(|py| {
            let mut constraints = AtomConstraintsAst::new(py, vec![]);
            constraints
                .set_tetrahedral_stereo(py, TetrahedralStereoArg::Config(TetrahedralStereo::Cw))
                .unwrap();
            match constraints.tetrahedral_stereo(py).unwrap().unwrap() {
                TetrahedralStereoAst::Stereo(coset) => {
                    assert_eq!(
                        coset.bind(py).borrow().to_ast(py),
                        AstStereoCosetAst::Lit(1)
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
                MoleculeAst::from_inner(AstMoleculeAst::from_parts(AstMoleculeParts {
                    atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                    ..Default::default()
                })),
            )
            .unwrap();
            let view = AtomConstraintsView {
                backing: AtomConstraintsBacking::Molecule {
                    owner: owner.clone_ref(py),
                    id: AstAtomId(0),
                },
            };
            view.set_aromatic_valence(py, AromaticValenceArg::Value(ValueArg::Lit(1)))
                .unwrap();
            let fresh = AtomConstraintsView {
                backing: AtomConstraintsBacking::Molecule {
                    owner,
                    id: AstAtomId(0),
                },
            };
            match fresh.aromatic_valence(py).unwrap().unwrap() {
                AromaticValenceAst::Aromatic(v) => assert_eq!(v.to_ast(py), AstValueAst::Lit(1)),
                _ => panic!("expected Aromatic"),
            }
        });
    }

    #[rstest]
    fn test_ring_size_counts_value_backed() {
        Python::attach(|py| {
            let constraints = Py::new(py, AtomConstraintsAst::new(py, vec![])).unwrap();
            let proxy = AtomConstraintsAst::ring_size_count(constraints.clone_ref(py));
            proxy.__setitem__(py, 6, ValueArg::Lit(3));
            assert_eq!(
                proxy.__getitem__(py, 6).unwrap().unwrap().to_ast(py),
                AstValueAst::Lit(3)
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
                MoleculeAst::from_inner(AstMoleculeAst::from_parts(AstMoleculeParts {
                    atoms: vec![AstAtomAst::from_element(ChemElement::C)],
                    ..Default::default()
                })),
            )
            .unwrap();
            let view = AtomConstraintsView {
                backing: AtomConstraintsBacking::Molecule {
                    owner: owner.clone_ref(py),
                    id: AstAtomId(0),
                },
            };
            view.ring_size_count(py)
                .__setitem__(py, 5, ValueArg::Lit(1));
            let fresh = AtomConstraintsView {
                backing: AtomConstraintsBacking::Molecule {
                    owner,
                    id: AstAtomId(0),
                },
            };
            assert_eq!(
                fresh
                    .ring_size_count(py)
                    .__getitem__(py, 5)
                    .unwrap()
                    .unwrap()
                    .to_ast(py),
                AstValueAst::Lit(1)
            );
        });
    }

    #[rstest]
    fn test_atom_constraints_ast_update_entries() {
        Python::attach(|py| {
            let constraints = Py::new(py, AtomConstraintsAst::new(py, vec![])).unwrap();
            let valence = into_py_variant(
                py,
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::valence(4)).unwrap(),
            )
            .unwrap();
            let degree = into_py_variant(
                py,
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::degree(3)).unwrap(),
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
                AtomConstraintAst::from_ast(py, &AstAtomConstraintAst::valence(4)).unwrap(),
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
                    .to_ast(py),
                AstValueAst::Lit(4)
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
                    AstAtomAst::from_element(ChemElement::C)
                        .with_constraint(AstAtomConstraintAst::valence(4)),
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
                AstValueAst::Lit(4)
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
            proxy.__setitem__(py, 6, ValueArg::Lit(3));
            proxy.__setitem__(py, 5, ValueArg::Lit(1));
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
