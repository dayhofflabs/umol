//! Atom-field value types mirroring `umol_ast::ast` (`ElementAst`; `IsotopeMassAst`
//! and `SpinStateAst` follow at S2b/c, `AtomAst` at S3).

use std::collections::BTreeSet;
use std::vec::IntoIter;

use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use umol_ast::ast::{
    AtomAst as AstAtomAst, AtomId as AstAtomId, ElementAst as AstElementAst,
    IsotopeMassAst as AstIsotopeMassAst, MoleculeAst as AstMoleculeAst,
    SpinStateAst as AstSpinStateAst,
};
use umol_chem::element::Element as ChemElement;

use crate::convert::into_py_variant;
use crate::element::Element;
use crate::molecule::MoleculeAst;
use crate::value::{MemOp, ValueAst};

fn atom_id_error() -> PyErr {
    PyIndexError::new_err("atom id out of range")
}

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
    fn new(unpaired: Py<ValueAst>, multiplicity: Py<ValueAst>) -> Self {
        SpinStateAst {
            unpaired,
            multiplicity,
        }
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

/// An atom: element, isotope, charge, implicit hydrogens, lone pairs, and spin.
/// (Constraints are not yet exposed — S5.)
#[pyclass(eq)]
#[derive(PartialEq)]
pub struct AtomAst(AstAtomAst);

#[pymethods]
impl AtomAst {
    /// Construct from an element — a single `Element` or an `ElementAst` expression —
    /// optionally setting fields.
    #[new]
    #[pyo3(signature = (element, *, isotope_mass=None, charge=None, implicit_hydrogens=None, lone_pairs=None, spin=None))]
    fn new(
        py: Python<'_>,
        element: ElementInput,
        isotope_mass: Option<PyRef<'_, IsotopeMassAst>>,
        charge: Option<PyRef<'_, ValueAst>>,
        implicit_hydrogens: Option<PyRef<'_, ValueAst>>,
        lone_pairs: Option<PyRef<'_, ValueAst>>,
        spin: Option<PyRef<'_, SpinStateAst>>,
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
        ))
    }

    #[getter]
    fn element(&self) -> ElementAst {
        ElementAst::from_ast(&self.0.element)
    }

    #[getter]
    fn isotope_mass(&self) -> IsotopeMassAst {
        IsotopeMassAst::from_ast(&self.0.isotope_mass)
    }

    #[getter]
    fn charge(&self, py: Python<'_>) -> PyResult<ValueAst> {
        ValueAst::from_ast(py, &self.0.charge)
    }

    #[getter]
    fn implicit_hydrogens(&self, py: Python<'_>) -> PyResult<ValueAst> {
        ValueAst::from_ast(py, &self.0.implicit_hydrogens)
    }

    #[getter]
    fn lone_pairs(&self, py: Python<'_>) -> PyResult<ValueAst> {
        ValueAst::from_ast(py, &self.0.lone_pairs)
    }

    #[getter]
    fn spin(&self, py: Python<'_>) -> PyResult<SpinStateAst> {
        SpinStateAst::from_ast(py, &self.0.spin)
    }

    /// A copy with the given fields replaced — the Python-idiomatic single-copy
    /// builder (one clone per call, vs a `with_*` chain that copies at every step).
    #[pyo3(signature = (*, element=None, isotope_mass=None, charge=None, implicit_hydrogens=None, lone_pairs=None, spin=None))]
    fn replace(
        &self,
        py: Python<'_>,
        element: Option<ElementInput>,
        isotope_mass: Option<PyRef<'_, IsotopeMassAst>>,
        charge: Option<PyRef<'_, ValueAst>>,
        implicit_hydrogens: Option<PyRef<'_, ValueAst>>,
        lone_pairs: Option<PyRef<'_, ValueAst>>,
        spin: Option<PyRef<'_, SpinStateAst>>,
    ) -> Self {
        let mut atom = self.0.clone();
        if let Some(element) = element {
            atom = atom.with_element(element.to_ast(py));
        }
        AtomAst(apply_fields(
            atom,
            py,
            isotope_mass,
            charge,
            implicit_hydrogens,
            lone_pairs,
            spin,
        ))
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
        Ok(dict)
    }
}

/// The `element` argument to `AtomAst` — a concrete `Element` or an `ElementAst`
/// expression, extracted as a PyO3 `FromPyObject` union (tried in order). The
/// `Ast`/`Lit` variant names are the convention for these wrapper-or-literal inputs:
/// `Ast` = the `*Ast` wrapper, `Lit` = the concrete value (mirrors `ElementAst::Lit`).
#[derive(FromPyObject)]
enum ElementInput {
    Ast(Py<ElementAst>),
    Lit(Element),
}

impl ElementInput {
    fn to_ast(&self, py: Python<'_>) -> AstElementAst {
        match self {
            ElementInput::Ast(expr) => expr.bind(py).borrow().to_ast(),
            ElementInput::Lit(element) => AstElementAst::Lit(ChemElement::from(element)),
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn apply_fields(
    mut atom: AstAtomAst,
    py: Python<'_>,
    isotope_mass: Option<PyRef<'_, IsotopeMassAst>>,
    charge: Option<PyRef<'_, ValueAst>>,
    implicit_hydrogens: Option<PyRef<'_, ValueAst>>,
    lone_pairs: Option<PyRef<'_, ValueAst>>,
    spin: Option<PyRef<'_, SpinStateAst>>,
) -> AstAtomAst {
    if let Some(isotope_mass) = isotope_mass {
        atom = atom.with_isotope_mass(isotope_mass.to_ast());
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
    atom
}

/// An atom index into a molecule.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct AtomId(AstAtomId);

impl From<AstAtomId> for AtomId {
    fn from(id: AstAtomId) -> Self {
        AtomId(id)
    }
}

#[pymethods]
impl AtomId {
    #[new]
    fn new(index: u32) -> Self {
        AtomId(AstAtomId(index))
    }

    #[getter]
    fn index(&self) -> u32 {
        self.0 .0
    }

    fn __repr__(&self) -> String {
        format!("AtomId({})", self.0 .0)
    }
}

impl AtomAst {
    /// The wrapped AST atom — read access for molecule construction.
    pub(crate) fn inner(&self) -> &AstAtomAst {
        &self.0
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
            .ok_or_else(atom_id_error)
    }
}

#[pymethods]
impl AtomView {
    #[getter]
    fn id(&self) -> AtomId {
        AtomId(self.id)
    }

    #[getter]
    fn element(&self, py: Python<'_>) -> PyResult<ElementAst> {
        let molecule = self.owner.bind(py).borrow();
        Ok(ElementAst::from_ast(&self.atom(molecule.inner())?.element))
    }

    #[getter]
    fn isotope_mass(&self, py: Python<'_>) -> PyResult<IsotopeMassAst> {
        let molecule = self.owner.bind(py).borrow();
        Ok(IsotopeMassAst::from_ast(
            &self.atom(molecule.inner())?.isotope_mass,
        ))
    }

    #[getter]
    fn charge(&self, py: Python<'_>) -> PyResult<ValueAst> {
        let molecule = self.owner.bind(py).borrow();
        ValueAst::from_ast(py, &self.atom(molecule.inner())?.charge)
    }

    #[getter]
    fn implicit_hydrogens(&self, py: Python<'_>) -> PyResult<ValueAst> {
        let molecule = self.owner.bind(py).borrow();
        ValueAst::from_ast(py, &self.atom(molecule.inner())?.implicit_hydrogens)
    }

    #[getter]
    fn lone_pairs(&self, py: Python<'_>) -> PyResult<ValueAst> {
        let molecule = self.owner.bind(py).borrow();
        ValueAst::from_ast(py, &self.atom(molecule.inner())?.lone_pairs)
    }

    #[getter]
    fn spin(&self, py: Python<'_>) -> PyResult<SpinStateAst> {
        let molecule = self.owner.bind(py).borrow();
        SpinStateAst::from_ast(py, &self.atom(molecule.inner())?.spin)
    }
}

/// The atoms of a molecule, indexed by `AtomId`.
#[pyclass]
pub struct AtomViews {
    owner: Py<MoleculeAst>,
}

#[pymethods]
impl AtomViews {
    fn __len__(&self, py: Python<'_>) -> usize {
        self.owner.bind(py).borrow().inner().atoms().count()
    }

    fn __getitem__(&self, py: Python<'_>, id: AtomId) -> PyResult<AtomView> {
        let molecule = self.owner.bind(py).borrow();
        if molecule.inner().atoms().contains(id.0) {
            Ok(AtomView {
                owner: self.owner.clone_ref(py),
                id: id.0,
            })
        } else {
            Err(atom_id_error())
        }
    }

    /// The atom at `id`, or `None` if out of range.
    fn get(&self, py: Python<'_>, id: AtomId) -> Option<AtomView> {
        let molecule = self.owner.bind(py).borrow();
        molecule.inner().atoms().contains(id.0).then(|| AtomView {
            owner: self.owner.clone_ref(py),
            id: id.0,
        })
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
    use umol_ast::ast::{MemOp as AstMemOp, ValueAst as AstValueAst};

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
            assert_eq!(view.id().index(), 1);
            match view.element(py).unwrap() {
                ElementAst::Lit(e) => assert_eq!(ChemElement::from(&e), ChemElement::O),
                _ => panic!("expected Lit"),
            }
        });
    }

    #[rstest]
    fn test_atom_views_len_and_getitem() {
        Python::attach(|py| {
            let views = AtomViews {
                owner: carbon_oxygen(py),
            };
            assert_eq!(views.__len__(py), 2);
            assert_eq!(
                views.__getitem__(py, AtomId::new(0)).unwrap().id().index(),
                0
            );
            assert!(views.__getitem__(py, AtomId::new(5)).is_err());
        });
    }
}
