//! Owned multicenter-bond ASTs and molecule-backed multicenter-bond views.

use std::str::FromStr;
use std::vec::IntoIter;

use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use umol_graph_ir::ir::{
    AtomId as GraphIrAtomId, MoleculeAst as GraphIrMoleculeAst,
    MulticenterBondForm as GraphIrMulticenterBondForm,
    MulticenterBondId as GraphIrMulticenterBondId,
    MulticenterBondUpdate as GraphIrMulticenterBondUpdate,
    MulticenterBondView as GraphIrMulticenterBondView,
};

use crate::constraint::multicenter::{
    multicenter_bond_constraints_asdict, MulticenterBondConstraintsAst,
    MulticenterBondConstraintsBacking, MulticenterBondConstraintsLike,
    MulticenterBondConstraintsView,
};
#[cfg(test)]
use crate::constraint::multicenter::{
    MulticenterBondConstraintAst, MulticenterBondConstraintKey, MulticenterBondConstraintsUpdate,
};
use crate::convert::hash_rust;
use crate::electrons::{ElectronCountsAst, ElectronCountsLike};
use crate::error::parse_error;
use crate::lattice::impl_py_lattice;
use crate::molecule::MoleculeAst;
use crate::spin::{UnpairedElectronsAst, UnpairedElectronsUpdate};
use crate::value::{NumLike, ValueAst};

/// Attribute updates for a multicenter bond.
#[pyclass(frozen, skip_from_py_object)]
#[derive(Clone)]
pub struct MulticenterBondUpdate(GraphIrMulticenterBondUpdate);

#[pymethods]
impl MulticenterBondUpdate {
    #[new]
    #[pyo3(signature = (*, electrons=None, charge=None, unpaired_electrons=None, constraints=None))]
    fn new(
        py: Python<'_>,
        electrons: Option<ElectronCountsLike>,
        charge: Option<NumLike>,
        unpaired_electrons: Option<PyRef<'_, UnpairedElectronsUpdate>>,
        constraints: Option<Py<MulticenterBondConstraintsAst>>,
    ) -> Self {
        Self::from_rust(&GraphIrMulticenterBondUpdate {
            electrons: electrons.map(|value| value.to_rust(py)),
            charge: charge.map(|value| value.to_rust(py)),
            unpaired_electrons: unpaired_electrons
                .map(|value| value.to_rust(py))
                .unwrap_or_default(),
            constraints: constraints
                .map(|value| value.bind(py).borrow().inner().clone())
                .unwrap_or_default(),
        })
    }

    /// Parse a multicenter-bond-update DSL string into a `MulticenterBondUpdate`.
    #[staticmethod]
    fn parse(s: &str) -> PyResult<Self> {
        GraphIrMulticenterBondUpdate::from_str(s)
            .map(Self)
            .map_err(parse_error)
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("MulticenterBondUpdate.parse('{}')", self.0)
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    #[getter]
    fn electrons(&self) -> Option<ElectronCountsAst> {
        self.0.electrons.as_ref().map(ElectronCountsAst::from_rust)
    }

    #[getter]
    fn charge(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .charge
            .as_ref()
            .map(|value| ValueAst::from_rust(py, value))
            .transpose()
    }

    #[getter]
    fn unpaired_electrons(&self, py: Python<'_>) -> PyResult<UnpairedElectronsUpdate> {
        UnpairedElectronsUpdate::from_rust(py, &self.0.unpaired_electrons)
    }

    #[getter]
    fn constraints(&self) -> MulticenterBondConstraintsAst {
        MulticenterBondConstraintsAst::from_inner(self.0.constraints.clone())
    }
}

impl MulticenterBondUpdate {
    pub(crate) fn from_rust(update: &GraphIrMulticenterBondUpdate) -> Self {
        Self(update.clone())
    }

    pub(crate) fn to_rust(&self) -> GraphIrMulticenterBondUpdate {
        self.0.clone()
    }
}

/// A multicenter bond: a positional per-member-atom `electrons` vector, charge,
/// unpaired electrons, and multicenter-bond-scope constraints. The member atoms are
/// the participants of the owning molecule's multicenter relation (the view half); the
/// `electrons` vector is positional, aligned to that atom order.
#[pyclass(eq)]
#[derive(PartialEq)]
pub struct MulticenterBondAst(GraphIrMulticenterBondForm);

#[pymethods]
impl MulticenterBondAst {
    /// Construct from an electron-count vector — a `list[int]` or an
    /// `ElectronCountsAst` — optionally setting fields.
    #[new]
    #[pyo3(signature = (electrons, *, charge=None, unpaired_electrons=None, constraints=None))]
    fn new(
        py: Python<'_>,
        electrons: ElectronCountsLike,
        charge: Option<NumLike>,
        unpaired_electrons: Option<PyRef<'_, UnpairedElectronsAst>>,
        constraints: Option<Py<MulticenterBondConstraintsAst>>,
    ) -> Self {
        let mut bond = GraphIrMulticenterBondForm::new(electrons.to_rust(py));
        if let Some(charge) = charge {
            bond = bond.with_charge(charge.to_rust(py));
        }
        if let Some(unpaired_electrons) = unpaired_electrons {
            bond = bond.with_unpaired_electrons(unpaired_electrons.to_rust(py));
        }
        if let Some(constraints) = constraints {
            bond.constraints = constraints.bind(py).borrow().inner().clone();
        }
        MulticenterBondAst(bond)
    }

    /// Parse a multicenter-bond-DSL string (e.g. `"[1,1,1]#e6"`) into a `MulticenterBondAst`.
    #[staticmethod]
    fn parse(s: &str) -> PyResult<Self> {
        GraphIrMulticenterBondForm::from_str(s)
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
        ElectronCountsAst::from_rust(&self.0.electrons)
    }

    #[setter]
    fn set_electrons(&mut self, py: Python<'_>, value: ElectronCountsLike) {
        self.0.electrons = value.to_rust(py);
    }

    #[getter]
    fn charge(&self, py: Python<'_>) -> PyResult<ValueAst> {
        ValueAst::from_rust(py, &self.0.charge)
    }

    #[setter]
    fn set_charge(&mut self, py: Python<'_>, value: NumLike) {
        self.0.charge = value.to_rust(py);
    }

    #[getter]
    fn unpaired_electrons(&self, py: Python<'_>) -> PyResult<UnpairedElectronsAst> {
        UnpairedElectronsAst::from_rust(py, &self.0.unpaired_electrons)
    }

    #[setter]
    fn set_unpaired_electrons(&mut self, py: Python<'_>, value: PyRef<'_, UnpairedElectronsAst>) {
        self.0.unpaired_electrons = value.to_rust(py);
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
        slf: Py<Self>,
        py: Python<'_>,
        value: MulticenterBondConstraintsLike,
    ) -> PyResult<()> {
        let snapshot = value.to_rust(py)?;
        slf.borrow_mut(py).0.constraints = snapshot;
        Ok(())
    }

    /// The fields as a dict keyed by field name; values are Python objects.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("electrons", self.electrons())?;
        dict.set_item("charge", self.charge(py)?)?;
        dict.set_item("unpaired_electrons", self.unpaired_electrons(py)?)?;
        dict.set_item(
            "constraints",
            multicenter_bond_constraints_asdict(py, &self.0.constraints)?,
        )?;
        Ok(dict)
    }
}

impl MulticenterBondAst {
    /// The wrapped AST bond — read access for the bond-backed constraints view.
    pub(crate) fn inner(&self) -> &GraphIrMulticenterBondForm {
        &self.0
    }

    /// Mutable access to the wrapped AST bond — write access for the bond-backed
    /// constraints view.
    pub(crate) fn inner_mut(&mut self) -> &mut GraphIrMulticenterBondForm {
        &mut self.0
    }

    /// Wrap an owned Rust multicenter-bond AST.
    pub(crate) fn from_inner(bond: GraphIrMulticenterBondForm) -> Self {
        MulticenterBondAst(bond)
    }
}

impl_py_lattice!(
    MulticenterBondAst,
    GraphIrMulticenterBondForm,
    |value: &MulticenterBondAst, _py: Python<'_>| -> PyResult<GraphIrMulticenterBondForm> {
        Ok(value.inner().clone())
    },
    |_py: Python<'_>, value: GraphIrMulticenterBondForm| -> PyResult<MulticenterBondAst> {
        Ok(MulticenterBondAst::from_inner(value))
    }
);

/// A view of one multicenter bond within a molecule: a handle to the molecule plus
/// the bond's index. Field reads rebuild the transient Rust view; the molecule is
/// never copied. The member atom indices are read-only topology; the electrons,
/// charge, unpaired electrons, and constraints are the mutable bond value.
#[pyclass]
pub struct MulticenterBondView {
    owner: Py<MoleculeAst>,
    id: GraphIrMulticenterBondId,
}

impl MulticenterBondView {
    fn multicenter_bond<'a>(
        &self,
        molecule: &'a GraphIrMoleculeAst,
    ) -> PyResult<GraphIrMulticenterBondView<'a>> {
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
        Ok(ElectronCountsAst::from_rust(
            &self.multicenter_bond(molecule.inner())?.ast.electrons,
        ))
    }

    #[setter]
    fn set_electrons(&self, py: Python<'_>, value: ElectronCountsLike) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .multicenter_bond_mut(self.id)
            .ast
            .electrons = value.to_rust(py);
    }

    #[getter]
    fn charge(&self, py: Python<'_>) -> PyResult<ValueAst> {
        let molecule = self.owner.bind(py).borrow();
        ValueAst::from_rust(py, &self.multicenter_bond(molecule.inner())?.ast.charge)
    }

    #[setter]
    fn set_charge(&self, py: Python<'_>, value: NumLike) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .multicenter_bond_mut(self.id)
            .ast
            .charge = value.to_rust(py);
    }

    #[getter]
    fn unpaired_electrons(&self, py: Python<'_>) -> PyResult<UnpairedElectronsAst> {
        let molecule = self.owner.bind(py).borrow();
        UnpairedElectronsAst::from_rust(
            py,
            &self
                .multicenter_bond(molecule.inner())?
                .ast
                .unpaired_electrons,
        )
    }

    #[setter]
    fn set_unpaired_electrons(&self, py: Python<'_>, value: PyRef<'_, UnpairedElectronsAst>) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .multicenter_bond_mut(self.id)
            .ast
            .unpaired_electrons = value.to_rust(py);
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
        value: MulticenterBondConstraintsLike,
    ) -> PyResult<()> {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .multicenter_bond_mut(self.id)
            .ast
            .constraints = value.to_rust(py)?;
        Ok(())
    }

    /// The value fields as a dict keyed by field name; values are Python objects —
    /// symmetric with `MulticenterBondAst.asdict`, read through the view.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let molecule = self.owner.bind(py).borrow();
        let bond = self.multicenter_bond(molecule.inner())?.ast;
        let dict = PyDict::new(py);
        dict.set_item("electrons", ElectronCountsAst::from_rust(&bond.electrons))?;
        dict.set_item("charge", ValueAst::from_rust(py, &bond.charge)?)?;
        dict.set_item(
            "unpaired_electrons",
            UnpairedElectronsAst::from_rust(py, &bond.unpaired_electrons)?,
        )?;
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
    molecule: &GraphIrMoleculeAst,
    index: isize,
) -> PyResult<GraphIrMulticenterBondId> {
    let count = molecule.multicenter_bonds().count();
    let resolved = if index < 0 {
        index + count as isize
    } else {
        index
    };
    if resolved < 0 {
        return Err(PyIndexError::new_err("multicenter bond id out of range"));
    }
    let id = GraphIrMulticenterBondId(resolved as u32);
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
    fn of(&self, py: Python<'_>, atoms: Vec<u32>) -> Option<MulticenterBondView> {
        let molecule = self.owner.bind(py).borrow();
        molecule
            .inner()
            .multicenter_bonds()
            .of_id(atoms.into_iter().map(GraphIrAtomId))
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
            .incident_ids(GraphIrAtomId(atom))
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
    ids: IntoIter<GraphIrMulticenterBondId>,
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

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_chem::element::Element as ChemElement;
    use umol_graph_ir::ir::{
        AtomForm as GraphIrAtomForm, AtomId as GraphIrAtomId,
        ElectronCountsForm as GraphIrElectronCountsForm, MoleculeEntries,
        MulticenterBondConstraintAst as GraphIrMulticenterBondConstraintAst,
        MulticenterBondConstraintKey as GraphIrMulticenterBondConstraintKey,
        MulticenterBondConstraintsAst as GraphIrMulticenterBondConstraintsAst,
        NumForm as GraphIrNumForm, UnpairedElectronsForm as GraphIrUnpairedElectronsForm,
    };

    use super::*;
    use crate::convert::into_py_variant;

    /// Three borons (atom ids 0–2) joined by one 3-center multicenter bond over all
    /// three (electrons `[1,1,1]`), multicenter bond id 0.
    fn three_center_bond(py: Python<'_>) -> Py<MoleculeAst> {
        let molecule = GraphIrMoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![GraphIrAtomForm::from_element(ChemElement::B); 3],
            multicenter: vec![(
                (0u32..3).map(GraphIrAtomId).collect(),
                GraphIrMulticenterBondForm::from_electrons(vec![1, 1, 1]),
            )],
            ..Default::default()
        });
        Py::new(py, MoleculeAst::from_rust(molecule)).unwrap()
    }

    #[rstest]
    fn test_multicenter_bond_ast_new() {
        Python::attach(|py| {
            let unpaired_electrons_ast = GraphIrUnpairedElectronsForm::from((0_u8, 1_u8));
            let unpaired_electrons = Py::new(
                py,
                UnpairedElectronsAst::from_rust(py, &unpaired_electrons_ast).unwrap(),
            )
            .unwrap();
            let bond = MulticenterBondAst::new(
                py,
                ElectronCountsLike::Lit(vec![1, 1, 1]),
                Some(NumLike::Lit(-2)),
                Some(unpaired_electrons.bind(py).borrow()),
                None,
            );
            assert_eq!(
                bond.inner().electrons,
                GraphIrElectronCountsForm::Lit(vec![1, 1, 1])
            );
            assert_eq!(bond.inner().charge, GraphIrNumForm::Lit(-2));
            assert_eq!(bond.inner().unpaired_electrons, unpaired_electrons_ast);
        });
    }

    #[rstest]
    fn test_multicenter_bond_ast_new_constraints() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_rust(
                    py,
                    &GraphIrMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints =
                Py::new(py, MulticenterBondConstraintsAst::new(py, vec![ec])).unwrap();
            let bond = MulticenterBondAst::new(
                py,
                ElectronCountsLike::Lit(vec![1, 1, 1]),
                None,
                None,
                Some(constraints),
            );
            assert_eq!(
                bond.inner().constraints.electron_count(),
                GraphIrNumForm::Lit(6)
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
                MulticenterBondAst::from_inner(GraphIrMulticenterBondForm::from_electrons(vec![
                    1, 1, 1,
                ]));
            assert_eq!(
                bond.electrons().to_rust(),
                GraphIrElectronCountsForm::Lit(vec![1, 1, 1])
            );
            bond.set_electrons(py, ElectronCountsLike::Lit(vec![2, 2]));
            assert_eq!(
                bond.electrons().to_rust(),
                GraphIrElectronCountsForm::Lit(vec![2, 2])
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_ast_charge() {
        Python::attach(|py| {
            let mut bond =
                MulticenterBondAst::from_inner(GraphIrMulticenterBondForm::from_electrons(vec![
                    1, 1, 1,
                ]));
            bond.set_charge(py, NumLike::Lit(-1));
            assert_eq!(
                bond.charge(py).unwrap().to_rust(py),
                GraphIrNumForm::Lit(-1)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_ast_unpaired_electrons() {
        Python::attach(|py| {
            let unpaired_electrons_ast = GraphIrUnpairedElectronsForm::from((0_u8, 1_u8));
            let unpaired_electrons = Py::new(
                py,
                UnpairedElectronsAst::from_rust(py, &unpaired_electrons_ast).unwrap(),
            )
            .unwrap();
            let mut bond =
                MulticenterBondAst::from_inner(GraphIrMulticenterBondForm::from_electrons(vec![
                    1, 1, 1,
                ]));
            bond.set_unpaired_electrons(py, unpaired_electrons.bind(py).borrow());
            assert_eq!(
                bond.unpaired_electrons(py).unwrap().to_rust(py),
                unpaired_electrons_ast
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_ast_set_constraints_from_view() {
        Python::attach(|py| {
            let src = Py::new(
                py,
                MulticenterBondAst::from_inner(
                    GraphIrMulticenterBondForm::from_electrons(vec![1, 1, 1])
                        .with_constraint(GraphIrMulticenterBondConstraintAst::electron_count(6)),
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
            let dst = Py::new(
                py,
                MulticenterBondAst::from_inner(GraphIrMulticenterBondForm::from_electrons(vec![
                    1, 1, 1,
                ])),
            )
            .unwrap();
            MulticenterBondAst::set_constraints(
                dst.clone_ref(py),
                py,
                MulticenterBondConstraintsLike::View(view),
            )
            .unwrap();
            assert_eq!(
                dst.bind(py).borrow().inner().constraints.electron_count(),
                GraphIrNumForm::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_ast_asdict() {
        Python::attach(|py| {
            let bond = MulticenterBondAst::from_inner(
                GraphIrMulticenterBondForm::from_electrons(vec![1, 1, 1])
                    .with_constraint(GraphIrMulticenterBondConstraintAst::electron_count(6)),
            );
            let dict = bond.asdict(py).unwrap();
            assert_eq!(dict.len(), 4);
            let electrons = dict.get_item("electrons").unwrap().unwrap();
            let expected = into_py_variant(py, ElectronCountsAst::Lit(vec![1, 1, 1])).unwrap();
            assert!(electrons.eq(expected.bind(py)).unwrap());
            assert!(dict.contains("charge").unwrap());
            assert!(dict.contains("unpaired_electrons").unwrap());
            assert!(dict.contains("constraints").unwrap());
        });
    }

    #[rstest]
    fn test_multicenter_bond_view_atom_ids() {
        Python::attach(|py| {
            let view = MulticenterBondView {
                owner: three_center_bond(py),
                id: GraphIrMulticenterBondId(0),
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
                id: GraphIrMulticenterBondId(0),
            };
            assert_eq!(
                view.electrons(py).unwrap().to_rust(),
                GraphIrElectronCountsForm::Lit(vec![1, 1, 1])
            );
            view.set_electrons(py, ElectronCountsLike::Lit(vec![2, 2, 2]));
            let fresh = MulticenterBondView {
                owner,
                id: GraphIrMulticenterBondId(0),
            };
            assert_eq!(
                fresh.electrons(py).unwrap().to_rust(),
                GraphIrElectronCountsForm::Lit(vec![2, 2, 2])
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_view_charge() {
        Python::attach(|py| {
            let owner = three_center_bond(py);
            let view = MulticenterBondView {
                owner: owner.clone_ref(py),
                id: GraphIrMulticenterBondId(0),
            };
            view.set_charge(py, NumLike::Lit(-1));
            let fresh = MulticenterBondView {
                owner,
                id: GraphIrMulticenterBondId(0),
            };
            assert_eq!(
                fresh.charge(py).unwrap().to_rust(py),
                GraphIrNumForm::Lit(-1)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_view_unpaired_electrons() {
        Python::attach(|py| {
            let unpaired_electrons_ast = GraphIrUnpairedElectronsForm::from((0_u8, 1_u8));
            let unpaired_electrons = Py::new(
                py,
                UnpairedElectronsAst::from_rust(py, &unpaired_electrons_ast).unwrap(),
            )
            .unwrap();
            let owner = three_center_bond(py);
            let view = MulticenterBondView {
                owner: owner.clone_ref(py),
                id: GraphIrMulticenterBondId(0),
            };
            view.set_unpaired_electrons(py, unpaired_electrons.bind(py).borrow());
            let fresh = MulticenterBondView {
                owner,
                id: GraphIrMulticenterBondId(0),
            };
            assert_eq!(
                fresh.unpaired_electrons(py).unwrap().to_rust(py),
                unpaired_electrons_ast
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_view_constraints() {
        Python::attach(|py| {
            let view = MulticenterBondView {
                owner: three_center_bond(py),
                id: GraphIrMulticenterBondId(0),
            };
            match view.constraints(py).backing {
                MulticenterBondConstraintsBacking::Molecule { id, .. } => {
                    assert_eq!(id, GraphIrMulticenterBondId(0))
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
                id: GraphIrMulticenterBondId(0),
            };
            let constraints = Py::new(
                py,
                MulticenterBondConstraintsAst::new(
                    py,
                    vec![into_py_variant(
                        py,
                        MulticenterBondConstraintAst::from_rust(
                            py,
                            &GraphIrMulticenterBondConstraintAst::electron_count(6),
                        )
                        .unwrap(),
                    )
                    .unwrap()],
                ),
            )
            .unwrap();
            view.set_constraints(py, MulticenterBondConstraintsLike::Container(constraints))
                .unwrap();
            let fresh = MulticenterBondView {
                owner,
                id: GraphIrMulticenterBondId(0),
            };
            assert_eq!(
                fresh
                    .constraints(py)
                    .electron_count(py)
                    .unwrap()
                    .to_rust(py),
                GraphIrNumForm::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_view_asdict() {
        Python::attach(|py| {
            let view = MulticenterBondView {
                owner: three_center_bond(py),
                id: GraphIrMulticenterBondId(0),
            };
            let dict = view.asdict(py).unwrap();
            assert_eq!(dict.len(), 4);
            let electrons = dict.get_item("electrons").unwrap().unwrap();
            let expected = into_py_variant(py, ElectronCountsAst::Lit(vec![1, 1, 1])).unwrap();
            assert!(electrons.eq(expected.bind(py)).unwrap());
            assert!(dict.contains("charge").unwrap());
            assert!(dict.contains("unpaired_electrons").unwrap());
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
                MulticenterBondAst::from_inner(GraphIrMulticenterBondForm::from_electrons(vec![
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
                view.electrons(py).unwrap().to_rust(),
                GraphIrElectronCountsForm::Lit(vec![2, 2, 2])
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
                MulticenterBondAst::from_inner(GraphIrMulticenterBondForm::from_electrons(vec![
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
    fn test_multicenter_bond_views_of() {
        Python::attach(|py| {
            let views = MulticenterBondViews {
                owner: three_center_bond(py),
            };
            assert_eq!(views.of(py, vec![0, 1, 2]).unwrap().id(), 0);
            // a subset is not the bond's exact atom set
            assert!(views.of(py, vec![0, 1]).is_none());
        });
    }

    #[rstest]
    fn test_multicenter_bond_views_incident() {
        Python::attach(|py| {
            // three borons bonded plus one isolated boron (atom id 3)
            let molecule = GraphIrMoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![GraphIrAtomForm::from_element(ChemElement::B); 4],
                multicenter: vec![(
                    (0u32..3).map(GraphIrAtomId).collect(),
                    GraphIrMulticenterBondForm::from_electrons(vec![1, 1, 1]),
                )],
                ..Default::default()
            });
            let views = MulticenterBondViews {
                owner: Py::new(py, MoleculeAst::from_rust(molecule)).unwrap(),
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
        let key = MulticenterBondConstraintKey::from_rust(
            &GraphIrMulticenterBondConstraintKey::ElectronCount,
        );
        assert_eq!(
            key.to_rust(),
            GraphIrMulticenterBondConstraintKey::ElectronCount
        );
    }

    #[rstest]
    fn test_multicenter_bond_constraint_ast_key() {
        Python::attach(|py| {
            let constraint = GraphIrMulticenterBondConstraintAst::electron_count(6);
            let key = MulticenterBondConstraintAst::from_rust(py, &constraint)
                .unwrap()
                .key(py);
            assert_eq!(
                key.to_rust(),
                GraphIrMulticenterBondConstraintKey::ElectronCount
            );
        });
    }

    #[rstest]
    #[case(GraphIrMulticenterBondConstraintAst::electron_count(6))]
    #[case(GraphIrMulticenterBondConstraintAst::electron_count(GraphIrNumForm::Undetermined))]
    fn test_multicenter_bond_constraint_ast_roundtrip(
        #[case] ast: GraphIrMulticenterBondConstraintAst,
    ) {
        Python::attach(|py| {
            assert_eq!(
                MulticenterBondConstraintAst::from_rust(py, &ast)
                    .unwrap()
                    .to_rust(py),
                ast
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_new() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_rust(
                    py,
                    &GraphIrMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = MulticenterBondConstraintsAst::new(py, vec![ec]);
            assert_eq!(constraints.__len__(), 1);
            assert_eq!(
                constraints.electron_count(py).unwrap().to_rust(py),
                GraphIrNumForm::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_repr() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_rust(
                    py,
                    &GraphIrMulticenterBondConstraintAst::electron_count(6),
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
                MulticenterBondConstraintAst::from_rust(
                    py,
                    &GraphIrMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            constraints.set(py, ec);
            assert_eq!(constraints.__len__(), 1);
            assert_eq!(
                constraints.electron_count(py).unwrap().to_rust(py),
                GraphIrNumForm::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_pop() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_rust(
                    py,
                    &GraphIrMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let mut constraints = MulticenterBondConstraintsAst::new(py, vec![ec]);
            let key = into_py_variant(py, MulticenterBondConstraintKey::ElectronCount()).unwrap();
            let removed = constraints.pop(py, key).unwrap();
            match removed {
                Some(MulticenterBondConstraintAst::ElectronCount(v)) => {
                    assert_eq!(v.bind(py).borrow().to_rust(py), GraphIrNumForm::Lit(6))
                }
                _ => panic!("expected removed ElectronCount(Lit(6))"),
            }
            assert_eq!(constraints.__len__(), 0);
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_update() {
        Python::attach(|py| {
            let constraints = Py::new(py, MulticenterBondConstraintsAst::new(py, vec![])).unwrap();
            let mut other = GraphIrMulticenterBondConstraintsAst::new();
            other.set(GraphIrMulticenterBondConstraintAst::electron_count(6));
            MulticenterBondConstraintsAst::update(
                constraints.clone_ref(py),
                py,
                MulticenterBondConstraintsUpdate::Container(
                    Py::new(py, MulticenterBondConstraintsAst::from_inner(other)).unwrap(),
                ),
            )
            .unwrap();
            let c = constraints.bind(py).borrow();
            assert_eq!(c.__len__(), 1);
            assert_eq!(
                c.electron_count(py).unwrap().to_rust(py),
                GraphIrNumForm::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_update_entries() {
        Python::attach(|py| {
            let constraints = Py::new(py, MulticenterBondConstraintsAst::new(py, vec![])).unwrap();
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_rust(
                    py,
                    &GraphIrMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            MulticenterBondConstraintsAst::update(
                constraints.clone_ref(py),
                py,
                MulticenterBondConstraintsUpdate::Entries(vec![ec]),
            )
            .unwrap();
            let c = constraints.bind(py).borrow();
            assert_eq!(c.__len__(), 1);
            assert_eq!(
                c.electron_count(py).unwrap().to_rust(py),
                GraphIrNumForm::Lit(6)
            );
        });
    }

    /// Regression: a container updating itself resolves `other` before the write borrow,
    /// so it is an idempotent no-op, not a RefCell double-borrow panic.
    #[rstest]
    fn test_multicenter_bond_constraints_ast_update_self() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_rust(
                    py,
                    &GraphIrMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints =
                Py::new(py, MulticenterBondConstraintsAst::new(py, vec![ec])).unwrap();
            MulticenterBondConstraintsAst::update(
                constraints.clone_ref(py),
                py,
                MulticenterBondConstraintsUpdate::Container(constraints.clone_ref(py)),
            )
            .unwrap();
            assert_eq!(
                constraints
                    .bind(py)
                    .borrow()
                    .electron_count(py)
                    .unwrap()
                    .to_rust(py),
                GraphIrNumForm::Lit(6)
            );
        });
    }

    /// Regression: assigning a bond's own constraints view back to it snapshots before
    /// the write borrow, so it is a no-op, not a double-borrow panic.
    #[rstest]
    fn test_multicenter_bond_ast_set_constraints_self() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                MulticenterBondAst::from_inner(
                    GraphIrMulticenterBondForm::from_electrons(vec![1, 1, 1])
                        .with_constraint(GraphIrMulticenterBondConstraintAst::electron_count(6)),
                ),
            )
            .unwrap();
            let own_view = Py::new(
                py,
                MulticenterBondConstraintsView {
                    backing: MulticenterBondConstraintsBacking::MulticenterBond(bond.clone_ref(py)),
                },
            )
            .unwrap();
            MulticenterBondAst::set_constraints(
                bond.clone_ref(py),
                py,
                MulticenterBondConstraintsLike::View(own_view),
            )
            .unwrap();
            assert_eq!(
                bond.bind(py).borrow().inner().constraints.electron_count(),
                GraphIrNumForm::Lit(6)
            );
        });
    }

    /// Regression: a view updating from a view over the same bond resolves `other`
    /// before the write borrow, so it is an idempotent no-op, not a double-borrow panic.
    #[rstest]
    fn test_multicenter_bond_constraints_view_update_self() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                MulticenterBondAst::from_inner(
                    GraphIrMulticenterBondForm::from_electrons(vec![1, 1, 1])
                        .with_constraint(GraphIrMulticenterBondConstraintAst::electron_count(6)),
                ),
            )
            .unwrap();
            let view = MulticenterBondConstraintsView {
                backing: MulticenterBondConstraintsBacking::MulticenterBond(bond.clone_ref(py)),
            };
            let other = Py::new(
                py,
                MulticenterBondConstraintsView {
                    backing: MulticenterBondConstraintsBacking::MulticenterBond(bond.clone_ref(py)),
                },
            )
            .unwrap();
            view.update(py, MulticenterBondConstraintsUpdate::View(other))
                .unwrap();
            assert_eq!(
                bond.bind(py).borrow().inner().constraints.electron_count(),
                GraphIrNumForm::Lit(6)
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
                MulticenterBondConstraintAst::from_rust(
                    py,
                    &GraphIrMulticenterBondConstraintAst::electron_count(6),
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
                MulticenterBondConstraintAst::from_rust(
                    py,
                    &GraphIrMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = MulticenterBondConstraintsAst::new(py, vec![ec]);

            let mut keys = constraints.__iter__(py).unwrap();
            assert_eq!(
                keys.__next__().unwrap().bind(py).borrow().to_rust(),
                GraphIrMulticenterBondConstraintKey::ElectronCount
            );
            assert!(keys.__next__().is_none());

            let mut values = constraints.values(py).unwrap();
            assert_eq!(
                values.__next__().unwrap().bind(py).borrow().to_rust(py),
                GraphIrMulticenterBondConstraintAst::electron_count(6)
            );

            let mut items = constraints.items(py).unwrap();
            let (key, value) = items.__next__().unwrap();
            assert_eq!(
                key.bind(py).borrow().to_rust(),
                GraphIrMulticenterBondConstraintKey::ElectronCount
            );
            assert_eq!(
                value.bind(py).borrow().to_rust(py),
                GraphIrMulticenterBondConstraintAst::electron_count(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_get() {
        Python::attach(|py| {
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_rust(
                    py,
                    &GraphIrMulticenterBondConstraintAst::electron_count(6),
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
                MulticenterBondConstraintAst::from_rust(
                    py,
                    &GraphIrMulticenterBondConstraintAst::electron_count(6),
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
                empty.electron_count(py).unwrap().to_rust(py),
                GraphIrNumForm::Undetermined
            );
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_rust(
                    py,
                    &GraphIrMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = MulticenterBondConstraintsAst::new(py, vec![ec]);
            assert_eq!(
                constraints.electron_count(py).unwrap().to_rust(py),
                GraphIrNumForm::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_ast_set_electron_count() {
        Python::attach(|py| {
            let mut constraints = MulticenterBondConstraintsAst::new(py, vec![]);
            constraints.set_electron_count(py, NumLike::Lit(6));
            assert_eq!(
                constraints.electron_count(py).unwrap().to_rust(py),
                GraphIrNumForm::Lit(6)
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
                MulticenterBondConstraintAst::from_rust(
                    py,
                    &GraphIrMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = MulticenterBondConstraintsAst::new(py, vec![ec]);
            let dict = constraints.asdict(py).unwrap();
            assert_eq!(dict.len(), 1);
            let value = dict.get_item("electron_count").unwrap().unwrap();
            let expected = into_py_variant(
                py,
                ValueAst::from_rust(py, &GraphIrNumForm::Lit(6)).unwrap(),
            )
            .unwrap();
            assert!(value.eq(expected.bind(py)).unwrap());
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_view_set() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                MulticenterBondAst::from_inner(GraphIrMulticenterBondForm::from_electrons(vec![
                    1, 1, 1,
                ])),
            )
            .unwrap();
            let view = MulticenterBondConstraintsView {
                backing: MulticenterBondConstraintsBacking::MulticenterBond(bond.clone_ref(py)),
            };
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_rust(
                    py,
                    &GraphIrMulticenterBondConstraintAst::electron_count(6),
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
                fresh.electron_count(py).unwrap().to_rust(py),
                GraphIrNumForm::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_view_pop() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                MulticenterBondAst::from_inner(
                    GraphIrMulticenterBondForm::from_electrons(vec![1, 1, 1])
                        .with_constraint(GraphIrMulticenterBondConstraintAst::electron_count(6)),
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
                    assert_eq!(v.bind(py).borrow().to_rust(py), GraphIrNumForm::Lit(6))
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
                MulticenterBondAst::from_inner(GraphIrMulticenterBondForm::from_electrons(vec![
                    1, 1, 1,
                ])),
            )
            .unwrap();
            let view = MulticenterBondConstraintsView {
                backing: MulticenterBondConstraintsBacking::MulticenterBond(bond.clone_ref(py)),
            };
            let mut other = GraphIrMulticenterBondConstraintsAst::new();
            other.set(GraphIrMulticenterBondConstraintAst::electron_count(6));
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
                fresh.electron_count(py).unwrap().to_rust(py),
                GraphIrNumForm::Lit(6)
            );
        });
    }

    #[rstest]
    fn test_multicenter_bond_constraints_view_set_electron_count() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                MulticenterBondAst::from_inner(GraphIrMulticenterBondForm::from_electrons(vec![
                    1, 1, 1,
                ])),
            )
            .unwrap();
            let view = MulticenterBondConstraintsView {
                backing: MulticenterBondConstraintsBacking::MulticenterBond(bond.clone_ref(py)),
            };
            assert_eq!(
                view.electron_count(py).unwrap().to_rust(py),
                GraphIrNumForm::Undetermined
            );
            view.set_electron_count(py, NumLike::Lit(6));
            let fresh = MulticenterBondConstraintsView {
                backing: MulticenterBondConstraintsBacking::MulticenterBond(bond),
            };
            assert_eq!(
                fresh.electron_count(py).unwrap().to_rust(py),
                GraphIrNumForm::Lit(6)
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
                    id: GraphIrMulticenterBondId(0),
                },
            };
            let ec = into_py_variant(
                py,
                MulticenterBondConstraintAst::from_rust(
                    py,
                    &GraphIrMulticenterBondConstraintAst::electron_count(6),
                )
                .unwrap(),
            )
            .unwrap();
            view.set(py, ec);
            let fresh = MulticenterBondConstraintsView {
                backing: MulticenterBondConstraintsBacking::Molecule {
                    owner,
                    id: GraphIrMulticenterBondId(0),
                },
            };
            assert_eq!(fresh.__len__(py).unwrap(), 1);
            assert_eq!(
                fresh.electron_count(py).unwrap().to_rust(py),
                GraphIrNumForm::Lit(6)
            );
        });
    }
}
