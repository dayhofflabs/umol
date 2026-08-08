//! Owned bond ASTs and molecule-backed bond views.

use std::str::FromStr;
use std::vec::IntoIter;

use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use umol_graph_ir::ir::{
    AtomId as GraphIrAtomId, BondAst as GraphIrBondAst,
    BondConstraintAst as GraphIrBondConstraintAst, BondId as GraphIrBondId,
    BondUpdate as GraphIrBondUpdate, MoleculeAst as GraphIrMoleculeAst,
};

use crate::constraint::bond::{
    bond_constraints_asdict, BondConstraintsAst, BondConstraintsBacking, BondConstraintsLike,
    BondConstraintsView,
};
use crate::convert::hash_rust;
use crate::error::parse_error;
use crate::lattice::impl_py_lattice;
use crate::molecule::MoleculeAst;
use crate::spin::{UnpairedElectronsAst, UnpairedElectronsUpdate};
use crate::value::{ValueAst, ValueLike};

/// Attribute updates for a localized bond.
#[pyclass(frozen, skip_from_py_object)]
#[derive(Clone)]
pub struct BondUpdate(GraphIrBondUpdate);

#[pymethods]
impl BondUpdate {
    #[new]
    #[pyo3(signature = (*, order=None, charge=None, unpaired_electrons=None, constraints=None))]
    fn new(
        py: Python<'_>,
        order: Option<ValueLike>,
        charge: Option<ValueLike>,
        unpaired_electrons: Option<PyRef<'_, UnpairedElectronsUpdate>>,
        constraints: Option<Py<BondConstraintsAst>>,
    ) -> Self {
        Self::from_rust(&GraphIrBondUpdate {
            order: order.map(|value| value.to_rust(py)),
            charge: charge.map(|value| value.to_rust(py)),
            unpaired_electrons: unpaired_electrons
                .map(|value| value.to_rust(py))
                .unwrap_or_default(),
            constraints: constraints
                .map(|value| value.bind(py).borrow().inner().clone())
                .unwrap_or_default(),
        })
    }

    /// Parse a bond-update DSL string into a `BondUpdate`.
    #[staticmethod]
    fn parse(s: &str) -> PyResult<Self> {
        GraphIrBondUpdate::from_str(s)
            .map(Self)
            .map_err(parse_error)
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("BondUpdate.parse('{}')", self.0)
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    #[getter]
    fn order(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .order
            .as_ref()
            .map(|value| ValueAst::from_rust(py, value))
            .transpose()
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
    fn constraints(&self) -> BondConstraintsAst {
        BondConstraintsAst::from_inner(self.0.constraints.clone())
    }
}

impl BondUpdate {
    pub(crate) fn from_rust(update: &GraphIrBondUpdate) -> Self {
        Self(update.clone())
    }

    pub(crate) fn to_rust(&self) -> GraphIrBondUpdate {
        self.0.clone()
    }
}

/// A bond: order, charge, unpaired electrons, and bond-scope constraints.
#[pyclass(eq)]
#[derive(PartialEq)]
pub struct BondAst(GraphIrBondAst);

#[pymethods]
impl BondAst {
    /// Construct from an order — an `int` or a `ValueAst` expression — optionally
    /// setting fields.
    #[new]
    #[pyo3(signature = (order, *, charge=None, unpaired_electrons=None, constraints=None))]
    fn new(
        py: Python<'_>,
        order: ValueLike,
        charge: Option<ValueLike>,
        unpaired_electrons: Option<PyRef<'_, UnpairedElectronsAst>>,
        constraints: Option<Py<BondConstraintsAst>>,
    ) -> Self {
        let mut bond = GraphIrBondAst::new(order.to_rust(py));
        if let Some(charge) = charge {
            bond = bond.with_charge(charge.to_rust(py));
        }
        if let Some(unpaired_electrons) = unpaired_electrons {
            bond = bond.with_unpaired_electrons(unpaired_electrons.to_rust(py));
        }
        if let Some(constraints) = constraints {
            bond.constraints = constraints.bind(py).borrow().inner().clone();
        }
        BondAst(bond)
    }

    /// Construct the canonical `:single` bond shape.
    #[staticmethod]
    fn single() -> Self {
        Self(GraphIrBondAst::from_order(1))
    }

    /// Construct the canonical `:double` bond shape.
    #[staticmethod]
    fn double() -> Self {
        Self(GraphIrBondAst::from_order(2))
    }

    /// Construct the canonical `:triple` bond shape.
    #[staticmethod]
    fn triple() -> Self {
        Self(GraphIrBondAst::from_order(3))
    }

    /// Construct the canonical `:quadruple` bond shape.
    #[staticmethod]
    fn quadruple() -> Self {
        Self(GraphIrBondAst::from_order(4))
    }

    /// Construct the canonical `:aromatic` shape: an order-1 localized bond
    /// carrying the aromatic constraint, not an aromatic bond order.
    #[staticmethod]
    fn aromatic() -> Self {
        Self(
            GraphIrBondAst::from_order(1).with_constraint(GraphIrBondConstraintAst::aromatic(true)),
        )
    }

    /// Parse a bond-DSL string (e.g. `"2#c-1"`) into a `BondAst`.
    #[staticmethod]
    fn parse(s: &str) -> PyResult<Self> {
        GraphIrBondAst::from_str(s).map(Self).map_err(parse_error)
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("BondAst.parse('{}')", self.0)
    }

    #[getter]
    fn order(&self, py: Python<'_>) -> PyResult<ValueAst> {
        ValueAst::from_rust(py, &self.0.order)
    }

    #[setter]
    fn set_order(&mut self, py: Python<'_>, value: ValueLike) {
        self.0.order = value.to_rust(py);
    }

    #[getter]
    fn charge(&self, py: Python<'_>) -> PyResult<ValueAst> {
        ValueAst::from_rust(py, &self.0.charge)
    }

    #[setter]
    fn set_charge(&mut self, py: Python<'_>, value: ValueLike) {
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
    fn constraints(slf: Py<Self>) -> BondConstraintsView {
        BondConstraintsView {
            backing: BondConstraintsBacking::Bond(slf),
        }
    }

    /// Replace the whole constraint set (wipe-and-set) from a value container or
    /// a live view. Takes `slf` by handle and snapshots `value` *before* the write
    /// borrow, so `bond.constraints = bond.constraints` (a view over the same bond)
    /// reads while the bond is unborrowed instead of a double-borrow panic.
    #[setter]
    fn set_constraints(slf: Py<Self>, py: Python<'_>, value: BondConstraintsLike) -> PyResult<()> {
        let snapshot = value.to_rust(py)?;
        slf.borrow_mut(py).0.constraints = snapshot;
        Ok(())
    }

    /// The fields as a dict keyed by field name; values are Python objects.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("order", self.order(py)?)?;
        dict.set_item("charge", self.charge(py)?)?;
        dict.set_item("unpaired_electrons", self.unpaired_electrons(py)?)?;
        dict.set_item(
            "constraints",
            bond_constraints_asdict(py, &self.0.constraints)?,
        )?;
        Ok(dict)
    }
}

impl BondAst {
    /// The wrapped AST bond — read access for the bond-backed constraints view.
    pub(crate) fn inner(&self) -> &GraphIrBondAst {
        &self.0
    }

    /// Mutable access to the wrapped AST bond — write access for the bond-backed
    /// constraints view.
    pub(crate) fn inner_mut(&mut self) -> &mut GraphIrBondAst {
        &mut self.0
    }

    /// Wrap an AST bond (the hold-the-value `from_inner` bridge, paired with
    /// `inner`).
    pub(crate) fn from_inner(bond: GraphIrBondAst) -> Self {
        BondAst(bond)
    }
}

impl_py_lattice!(
    BondAst,
    GraphIrBondAst,
    |value: &BondAst, _py: Python<'_>| -> PyResult<GraphIrBondAst> { Ok(value.inner().clone()) },
    |_py: Python<'_>, value: GraphIrBondAst| -> PyResult<BondAst> {
        Ok(BondAst::from_inner(value))
    }
);

/// A view of one bond within a molecule: a handle to the molecule plus the bond's
/// index. Field reads rebuild the transient Rust view; the molecule is never copied.
#[pyclass]
pub struct BondView {
    owner: Py<MoleculeAst>,
    id: GraphIrBondId,
}

impl BondView {
    fn bond<'a>(&self, molecule: &'a GraphIrMoleculeAst) -> PyResult<&'a GraphIrBondAst> {
        molecule
            .bonds()
            .get(self.id)
            .map(|view| view.ast)
            .ok_or_else(|| PyIndexError::new_err("bond id out of range"))
    }
}

#[pymethods]
impl BondView {
    #[getter]
    fn id(&self) -> u32 {
        self.id.0
    }

    /// The two atom indices incident to this bond (read-only — endpoints are
    /// topology, not part of the bond value).
    #[getter]
    fn atom_ids(&self, py: Python<'_>) -> PyResult<(u32, u32)> {
        let molecule = self.owner.bind(py).borrow();
        let view = molecule
            .inner()
            .bonds()
            .get(self.id)
            .ok_or_else(|| PyIndexError::new_err("bond id out of range"))?;
        let [first, second] = view.atom_ids();
        Ok((first.0, second.0))
    }

    fn __repr__(&self) -> String {
        format!("BondView(id={})", self.id.0)
    }

    #[getter]
    fn order(&self, py: Python<'_>) -> PyResult<ValueAst> {
        let molecule = self.owner.bind(py).borrow();
        ValueAst::from_rust(py, &self.bond(molecule.inner())?.order)
    }

    #[setter]
    fn set_order(&self, py: Python<'_>, value: ValueLike) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .bond_mut(self.id)
            .ast
            .order = value.to_rust(py);
    }

    #[getter]
    fn charge(&self, py: Python<'_>) -> PyResult<ValueAst> {
        let molecule = self.owner.bind(py).borrow();
        ValueAst::from_rust(py, &self.bond(molecule.inner())?.charge)
    }

    #[setter]
    fn set_charge(&self, py: Python<'_>, value: ValueLike) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .bond_mut(self.id)
            .ast
            .charge = value.to_rust(py);
    }

    #[getter]
    fn unpaired_electrons(&self, py: Python<'_>) -> PyResult<UnpairedElectronsAst> {
        let molecule = self.owner.bind(py).borrow();
        UnpairedElectronsAst::from_rust(py, &self.bond(molecule.inner())?.unpaired_electrons)
    }

    #[setter]
    fn set_unpaired_electrons(&self, py: Python<'_>, value: PyRef<'_, UnpairedElectronsAst>) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .bond_mut(self.id)
            .ast
            .unpaired_electrons = value.to_rust(py);
    }

    /// The bond's constraints as a live handle onto the molecule: reads borrow the
    /// current state, mutators write through to the bond in place.
    #[getter]
    fn constraints(&self, py: Python<'_>) -> BondConstraintsView {
        BondConstraintsView {
            backing: BondConstraintsBacking::Molecule {
                owner: self.owner.clone_ref(py),
                id: self.id,
            },
        }
    }

    /// Replace the whole constraint set of the backing bond in place (wipe-and-set)
    /// from a value container or a live view.
    #[setter]
    fn set_constraints(&self, py: Python<'_>, value: BondConstraintsLike) -> PyResult<()> {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .bond_mut(self.id)
            .ast
            .constraints = value.to_rust(py)?;
        Ok(())
    }

    /// The fields as a dict keyed by field name; values are Python objects —
    /// symmetric with `BondAst.asdict`, read through the view.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let molecule = self.owner.bind(py).borrow();
        let bond = self.bond(molecule.inner())?;
        let dict = PyDict::new(py);
        dict.set_item("order", ValueAst::from_rust(py, &bond.order)?)?;
        dict.set_item("charge", ValueAst::from_rust(py, &bond.charge)?)?;
        dict.set_item(
            "unpaired_electrons",
            UnpairedElectronsAst::from_rust(py, &bond.unpaired_electrons)?,
        )?;
        dict.set_item(
            "constraints",
            bond_constraints_asdict(py, &bond.constraints)?,
        )?;
        Ok(dict)
    }
}

/// Resolve a possibly-negative Python index (negative counts from the end) into an
/// existing bond id, or `IndexError`.
fn resolve_bond_index(molecule: &GraphIrMoleculeAst, index: isize) -> PyResult<GraphIrBondId> {
    let count = molecule.bonds().count();
    let resolved = if index < 0 {
        index + count as isize
    } else {
        index
    };
    if resolved < 0 {
        return Err(PyIndexError::new_err("bond id out of range"));
    }
    let id = GraphIrBondId(resolved as u32);
    if molecule.bonds().contains(id) {
        Ok(id)
    } else {
        Err(PyIndexError::new_err("bond id out of range"))
    }
}

/// The bonds of a molecule, indexed by integer position.
#[pyclass]
pub struct BondViews {
    owner: Py<MoleculeAst>,
}

#[pymethods]
impl BondViews {
    fn __len__(&self, py: Python<'_>) -> usize {
        self.owner.bind(py).borrow().inner().bonds().count()
    }

    fn __repr__(&self, py: Python<'_>) -> String {
        format!(
            "BondViews(len={})",
            self.owner.bind(py).borrow().inner().bonds().count()
        )
    }

    fn __getitem__(&self, py: Python<'_>, index: isize) -> PyResult<BondView> {
        let molecule = self.owner.bind(py).borrow();
        let id = resolve_bond_index(molecule.inner(), index)?;
        Ok(BondView {
            owner: self.owner.clone_ref(py),
            id,
        })
    }

    /// Replace the whole bond value at `index` in place (endpoints unchanged).
    fn __setitem__(&self, py: Python<'_>, index: isize, bond: PyRef<'_, BondAst>) -> PyResult<()> {
        let mut molecule = self.owner.borrow_mut(py);
        let id = resolve_bond_index(molecule.inner(), index)?;
        *molecule.inner_mut().bond_mut(id).ast = bond.inner().clone();
        Ok(())
    }

    /// The bond between atoms `first` and `second`, or `None`.
    fn of(&self, py: Python<'_>, first: u32, second: u32) -> Option<BondView> {
        let molecule = self.owner.bind(py).borrow();
        molecule
            .inner()
            .bonds()
            .of_id(GraphIrAtomId(first), GraphIrAtomId(second))
            .map(|id| BondView {
                owner: self.owner.clone_ref(py),
                id,
            })
    }

    fn __iter__(&self, py: Python<'_>) -> BondViewIter {
        let ids = self
            .owner
            .bind(py)
            .borrow()
            .inner()
            .bonds()
            .ids()
            .collect::<Vec<_>>();
        BondViewIter {
            owner: self.owner.clone_ref(py),
            ids: ids.into_iter(),
        }
    }
}

impl BondViews {
    /// Build the bond-views handle for `owner` (the `.bonds` accessor on the molecule).
    pub(crate) fn new(owner: Py<MoleculeAst>) -> BondViews {
        BondViews { owner }
    }
}

#[pyclass]
struct BondViewIter {
    owner: Py<MoleculeAst>,
    ids: IntoIter<GraphIrBondId>,
}

#[pymethods]
impl BondViewIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> Option<BondView> {
        self.ids.next().map(|id| BondView {
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
        AtomAst as GraphIrAtomAst, BondConstraintAst as GraphIrBondConstraintAst,
        BondConstraintKey as GraphIrBondConstraintKey,
        BondConstraintsAst as GraphIrBondConstraintsAst, BooleanAst as GraphIrBooleanAst,
        CisTransStereoAst as GraphIrCisTransStereoAst, MoleculeEntries as GraphIrMoleculeEntries,
        RingScope as GraphIrRingScope, StereoCoset as GraphIrStereoCoset,
        ValueAst as GraphIrValueAst,
    };

    use super::*;
    use crate::boolean::BooleanLike;
    use crate::constraint::bond::{BondConstraintAst, BondConstraintKey, BondConstraintsUpdate};
    use crate::convert::into_py_variant;
    use crate::stereo::{CisTransConfiguration, CisTransStereoAst, CisTransStereoLike};

    /// A two-carbon molecule joined by one double bond (bond id 0, atoms 0–1).
    fn ethene(py: Python<'_>) -> Py<MoleculeAst> {
        let molecule = GraphIrMoleculeAst::from_entries(GraphIrMoleculeEntries {
            atoms: vec![
                GraphIrAtomAst::from_element(ChemElement::C),
                GraphIrAtomAst::from_element(ChemElement::C),
            ],
            bonds: vec![(
                GraphIrAtomId(0),
                GraphIrAtomId(1),
                GraphIrBondAst::from_order(2),
            )],
            ..Default::default()
        });
        Py::new(py, MoleculeAst::from_rust(molecule)).unwrap()
    }

    #[rstest]
    #[case::single("1")]
    #[case::charge("2#c-")]
    #[case::aromatic("1#a")]
    #[case::ring_size("1#R(6)")]
    fn test_bond_ast_parse(#[case] dsl: &str) {
        let bond = BondAst::parse(dsl).unwrap();
        assert_eq!(bond.__str__(), dsl);
        assert_eq!(bond.__repr__(), format!("BondAst.parse('{dsl}')"));
    }

    #[rstest]
    fn test_bond_ast_parse_error() {
        assert!(BondAst::parse("x#").is_err());
    }

    #[rstest]
    fn test_bond_ast_constraints() {
        let bond = BondAst(GraphIrBondAst::from_order(1).with_constraint(
            GraphIrBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
        ));
        assert_eq!(bond.inner().constraints.len(), 1);
    }

    #[rstest]
    fn test_bond_ast_set_constraints_from_view() {
        Python::attach(|py| {
            let src = Py::new(
                py,
                BondAst::from_inner(GraphIrBondAst::from_order(1).with_constraint(
                    GraphIrBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )),
            )
            .unwrap();
            let view = Py::new(
                py,
                BondConstraintsView {
                    backing: BondConstraintsBacking::Bond(src),
                },
            )
            .unwrap();
            let dst = Py::new(py, BondAst::from_inner(GraphIrBondAst::from_order(2))).unwrap();
            BondAst::set_constraints(dst.clone_ref(py), py, BondConstraintsLike::View(view))
                .unwrap();
            assert_eq!(
                dst.bind(py).borrow().inner().constraints.aromatic(),
                GraphIrBooleanAst::Lit(true)
            );
        });
    }

    #[rstest]
    #[case(GraphIrBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)))]
    #[case(GraphIrBondConstraintAst::cis_trans_stereo(GraphIrCisTransStereoAst::NotStereo))]
    #[case(
        GraphIrBondConstraintAst::cis_trans_stereo(GraphIrCisTransStereoAst::Stereo(
            GraphIrStereoCoset::Lit(1)
        ))
    )]
    #[case(GraphIrBondConstraintAst::ring_membership(GraphIrRingScope::All, 2))]
    #[case(GraphIrBondConstraintAst::ring_membership(GraphIrRingScope::Size(6), 1))]
    fn test_bond_constraint_ast_roundtrip(#[case] ast: GraphIrBondConstraintAst) {
        Python::attach(|py| {
            assert_eq!(
                BondConstraintAst::from_rust(py, &ast).unwrap().to_rust(py),
                ast
            );
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_len_contains() {
        Python::attach(|py| {
            let aromatic = into_py_variant(
                py,
                BondConstraintAst::from_rust(
                    py,
                    &GraphIrBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let ring = into_py_variant(
                py,
                BondConstraintAst::from_rust(
                    py,
                    &GraphIrBondConstraintAst::ring_membership(GraphIrRingScope::All, 2),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = BondConstraintsAst::new(py, vec![aromatic, ring]);
            assert_eq!(constraints.__len__(), 2);
            assert!(constraints.__contains__(
                py,
                into_py_variant(py, BondConstraintKey::Aromatic()).unwrap()
            ));
            assert!(!constraints.__contains__(
                py,
                into_py_variant(py, BondConstraintKey::CisTransStereo()).unwrap()
            ));
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_keys_values_items() {
        Python::attach(|py| {
            let aromatic = into_py_variant(
                py,
                BondConstraintAst::from_rust(
                    py,
                    &GraphIrBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let ring = into_py_variant(
                py,
                BondConstraintAst::from_rust(
                    py,
                    &GraphIrBondConstraintAst::ring_membership(GraphIrRingScope::All, 2),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = BondConstraintsAst::new(py, vec![aromatic, ring]);

            let mut keys = constraints.__iter__(py).unwrap();
            assert_eq!(
                keys.__next__().unwrap().bind(py).borrow().to_rust(py),
                GraphIrBondConstraintKey::Aromatic
            );
            assert_eq!(
                keys.__next__().unwrap().bind(py).borrow().to_rust(py),
                GraphIrBondConstraintKey::RingMembership(GraphIrRingScope::All)
            );
            assert!(keys.__next__().is_none());

            let mut values = constraints.values(py).unwrap();
            assert_eq!(
                values.__next__().unwrap().bind(py).borrow().to_rust(py),
                GraphIrBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true))
            );

            let mut items = constraints.items(py).unwrap();
            let (key, value) = items.__next__().unwrap();
            assert_eq!(
                key.bind(py).borrow().to_rust(py),
                GraphIrBondConstraintKey::Aromatic
            );
            assert_eq!(
                value.bind(py).borrow().to_rust(py),
                GraphIrBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true))
            );
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_get() {
        Python::attach(|py| {
            let aromatic = into_py_variant(
                py,
                BondConstraintAst::from_rust(
                    py,
                    &GraphIrBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = BondConstraintsAst::new(py, vec![aromatic]);
            let present = constraints
                .get(
                    py,
                    into_py_variant(py, BondConstraintKey::Aromatic()).unwrap(),
                    None,
                )
                .unwrap();
            let expected = into_py_variant(
                py,
                BondConstraintAst::from_rust(
                    py,
                    &GraphIrBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap()
            .into_any();
            assert!(present.bind(py).eq(expected.bind(py)).unwrap());
            let absent = constraints
                .get(
                    py,
                    into_py_variant(py, BondConstraintKey::CisTransStereo()).unwrap(),
                    None,
                )
                .unwrap();
            assert!(absent.bind(py).is_none());
            let sentinel = into_py_variant(py, BondConstraintKey::CisTransStereo())
                .unwrap()
                .into_any();
            let defaulted = constraints
                .get(
                    py,
                    into_py_variant(py, BondConstraintKey::CisTransStereo()).unwrap(),
                    Some(sentinel.clone_ref(py)),
                )
                .unwrap();
            assert_eq!(defaulted.as_ptr(), sentinel.as_ptr());
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_aromatic() {
        Python::attach(|py| {
            let empty = BondConstraintsAst::new(py, vec![]);
            assert_eq!(empty.aromatic().to_rust(), GraphIrBooleanAst::Undetermined);
            assert!(empty.cis_trans_stereo(py).unwrap().is_none());
            assert!(empty.ring_count(py).unwrap().is_none());
            let aromatic = into_py_variant(
                py,
                BondConstraintAst::from_rust(
                    py,
                    &GraphIrBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = BondConstraintsAst::new(py, vec![aromatic]);
            assert_eq!(
                constraints.aromatic().to_rust(),
                GraphIrBooleanAst::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_ring_size_count() {
        Python::attach(|py| {
            let membership = into_py_variant(
                py,
                BondConstraintAst::from_rust(
                    py,
                    &GraphIrBondConstraintAst::ring_membership(GraphIrRingScope::Size(6), 1),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = Py::new(py, BondConstraintsAst::new(py, vec![membership])).unwrap();
            let proxy = BondConstraintsAst::ring_size_count(constraints.clone_ref(py));
            assert_eq!(
                proxy.__getitem__(py, 6).unwrap().unwrap().to_rust(py),
                GraphIrValueAst::Lit(1)
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
    fn test_bond_constraints_ast_set() {
        Python::attach(|py| {
            let mut constraints = BondConstraintsAst::new(py, vec![]);
            let aromatic = into_py_variant(
                py,
                BondConstraintAst::from_rust(
                    py,
                    &GraphIrBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            constraints.set(py, aromatic);
            assert_eq!(constraints.__len__(), 1);
            assert_eq!(
                constraints.aromatic().to_rust(),
                GraphIrBooleanAst::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_pop() {
        Python::attach(|py| {
            let aromatic = into_py_variant(
                py,
                BondConstraintAst::from_rust(
                    py,
                    &GraphIrBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let mut constraints = BondConstraintsAst::new(py, vec![aromatic]);
            let removed = constraints
                .pop(
                    py,
                    into_py_variant(py, BondConstraintKey::Aromatic()).unwrap(),
                )
                .unwrap();
            match removed {
                Some(BondConstraintAst::Aromatic(b)) => {
                    assert_eq!(b.bind(py).borrow().to_rust(), GraphIrBooleanAst::Lit(true))
                }
                _ => panic!("expected removed Aromatic(Lit(true))"),
            }
            assert_eq!(constraints.__len__(), 0);
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_update() {
        Python::attach(|py| {
            let constraints = Py::new(py, BondConstraintsAst::new(py, vec![])).unwrap();
            let mut other = GraphIrBondConstraintsAst::new();
            other.set(GraphIrBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(
                true,
            )));
            other.set(GraphIrBondConstraintAst::ring_membership(
                GraphIrRingScope::All,
                2,
            ));
            BondConstraintsAst::update(
                constraints.clone_ref(py),
                py,
                BondConstraintsUpdate::Container(
                    Py::new(py, BondConstraintsAst::from_inner(other)).unwrap(),
                ),
            )
            .unwrap();
            let c = constraints.bind(py).borrow();
            assert_eq!(c.__len__(), 2);
            assert_eq!(c.aromatic().to_rust(), GraphIrBooleanAst::Lit(true));
            assert_eq!(
                c.ring_count(py).unwrap().unwrap().to_rust(py),
                GraphIrValueAst::Lit(2)
            );
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_update_entries() {
        Python::attach(|py| {
            let constraints = Py::new(py, BondConstraintsAst::new(py, vec![])).unwrap();
            let aromatic = into_py_variant(
                py,
                BondConstraintAst::from_rust(
                    py,
                    &GraphIrBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let ring = into_py_variant(
                py,
                BondConstraintAst::from_rust(
                    py,
                    &GraphIrBondConstraintAst::ring_membership(GraphIrRingScope::All, 2),
                )
                .unwrap(),
            )
            .unwrap();
            BondConstraintsAst::update(
                constraints.clone_ref(py),
                py,
                BondConstraintsUpdate::Entries(vec![aromatic, ring]),
            )
            .unwrap();
            assert_eq!(constraints.bind(py).borrow().__len__(), 2);
        });
    }

    /// Regression: a container updating itself resolves `other` before the write borrow,
    /// so it is an idempotent no-op, not a RefCell double-borrow panic.
    #[rstest]
    fn test_bond_constraints_ast_update_self() {
        Python::attach(|py| {
            let aromatic = into_py_variant(
                py,
                BondConstraintAst::from_rust(
                    py,
                    &GraphIrBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = Py::new(py, BondConstraintsAst::new(py, vec![aromatic])).unwrap();
            BondConstraintsAst::update(
                constraints.clone_ref(py),
                py,
                BondConstraintsUpdate::Container(constraints.clone_ref(py)),
            )
            .unwrap();
            assert_eq!(
                constraints.bind(py).borrow().aromatic().to_rust(),
                GraphIrBooleanAst::Lit(true)
            );
        });
    }

    /// Regression: assigning a bond's own constraints view back to it snapshots before
    /// the write borrow, so it is a no-op, not a double-borrow panic
    /// (`bond.constraints = bond.constraints`).
    #[rstest]
    fn test_bond_ast_set_constraints_self() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                BondAst::from_inner(GraphIrBondAst::from_order(1).with_constraint(
                    GraphIrBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )),
            )
            .unwrap();
            let own_view = Py::new(
                py,
                BondConstraintsView {
                    backing: BondConstraintsBacking::Bond(bond.clone_ref(py)),
                },
            )
            .unwrap();
            BondAst::set_constraints(bond.clone_ref(py), py, BondConstraintsLike::View(own_view))
                .unwrap();
            assert_eq!(
                bond.bind(py).borrow().inner().constraints.aromatic(),
                GraphIrBooleanAst::Lit(true)
            );
        });
    }

    /// Regression: a view updating from a view over the same bond resolves `other`
    /// before the write borrow, so it is an idempotent no-op, not a double-borrow panic
    /// (`bond.constraints.update(bond.constraints)`).
    #[rstest]
    fn test_bond_constraints_view_update_self() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                BondAst::from_inner(GraphIrBondAst::from_order(1).with_constraint(
                    GraphIrBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )),
            )
            .unwrap();
            let view = BondConstraintsView {
                backing: BondConstraintsBacking::Bond(bond.clone_ref(py)),
            };
            let other = Py::new(
                py,
                BondConstraintsView {
                    backing: BondConstraintsBacking::Bond(bond.clone_ref(py)),
                },
            )
            .unwrap();
            view.update(py, BondConstraintsUpdate::View(other)).unwrap();
            assert_eq!(
                bond.bind(py).borrow().inner().constraints.aromatic(),
                GraphIrBooleanAst::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_set_aromatic() {
        Python::attach(|py| {
            let mut constraints = BondConstraintsAst::new(py, vec![]);
            constraints.set_aromatic(py, BooleanLike::Lit(true));
            assert_eq!(
                constraints.aromatic().to_rust(),
                GraphIrBooleanAst::Lit(true)
            );
            constraints.set_aromatic(py, BooleanLike::Lit(false));
            assert_eq!(
                constraints.aromatic().to_rust(),
                GraphIrBooleanAst::Lit(false)
            );
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_set_cis_trans_stereo() {
        Python::attach(|py| {
            let mut constraints = BondConstraintsAst::new(py, vec![]);
            constraints
                .set_cis_trans_stereo(py, CisTransStereoLike::Config(CisTransConfiguration::E))
                .unwrap();
            match constraints.cis_trans_stereo(py).unwrap().unwrap() {
                CisTransStereoAst::Stereo(coset) => {
                    assert_eq!(
                        coset.bind(py).borrow().to_rust(py),
                        GraphIrStereoCoset::Lit(1)
                    )
                }
                _ => panic!("expected Stereo"),
            }
            constraints
                .set_cis_trans_stereo(py, CisTransStereoLike::Flag(false))
                .unwrap();
            match constraints.cis_trans_stereo(py).unwrap().unwrap() {
                CisTransStereoAst::NotStereo() => {}
                _ => panic!("expected NotStereo"),
            }
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_set_cis_trans_stereo_error() {
        Python::attach(|py| {
            let mut constraints = BondConstraintsAst::new(py, vec![]);
            assert!(constraints
                .set_cis_trans_stereo(py, CisTransStereoLike::Flag(true))
                .is_err());
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_set_ring_count() {
        Python::attach(|py| {
            let mut constraints = BondConstraintsAst::new(py, vec![]);
            constraints.set_ring_count(py, ValueLike::Lit(2));
            assert_eq!(
                constraints.ring_count(py).unwrap().unwrap().to_rust(py),
                GraphIrValueAst::Lit(2)
            );
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_getitem_error() {
        Python::attach(|py| {
            let constraints = BondConstraintsAst::new(py, vec![]);
            let key = into_py_variant(py, BondConstraintKey::Aromatic()).unwrap();
            assert!(constraints.__getitem__(py, key).is_err());
        });
    }

    #[rstest]
    fn test_bond_constraints_ast_delitem_error() {
        Python::attach(|py| {
            let mut constraints = BondConstraintsAst::new(py, vec![]);
            let key = into_py_variant(py, BondConstraintKey::Aromatic()).unwrap();
            assert!(constraints.__delitem__(py, key).is_err());
        });
    }

    #[rstest]
    fn test_bond_constraints_view_set() {
        Python::attach(|py| {
            let bond = Py::new(py, BondAst::from_inner(GraphIrBondAst::from_order(1))).unwrap();
            let view = BondConstraintsView {
                backing: BondConstraintsBacking::Bond(bond.clone_ref(py)),
            };
            let aromatic = into_py_variant(
                py,
                BondConstraintAst::from_rust(
                    py,
                    &GraphIrBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            view.set(py, aromatic);
            // a fresh view proves the write hit the standalone bond, not a copy
            let fresh = BondConstraintsView {
                backing: BondConstraintsBacking::Bond(bond),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 1);
            match fresh
                .__getitem__(
                    py,
                    into_py_variant(py, BondConstraintKey::Aromatic()).unwrap(),
                )
                .unwrap()
            {
                BondConstraintAst::Aromatic(b) => {
                    assert_eq!(b.bind(py).borrow().to_rust(), GraphIrBooleanAst::Lit(true))
                }
                _ => panic!("expected Aromatic(Lit(true))"),
            }
        });
    }

    #[rstest]
    fn test_bond_constraints_view_pop() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                BondAst::from_inner(GraphIrBondAst::from_order(1).with_constraint(
                    GraphIrBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )),
            )
            .unwrap();
            let view = BondConstraintsView {
                backing: BondConstraintsBacking::Bond(bond.clone_ref(py)),
            };
            let removed = view
                .pop(
                    py,
                    into_py_variant(py, BondConstraintKey::Aromatic()).unwrap(),
                )
                .unwrap();
            match removed {
                Some(BondConstraintAst::Aromatic(b)) => {
                    assert_eq!(b.bind(py).borrow().to_rust(), GraphIrBooleanAst::Lit(true))
                }
                _ => panic!("expected removed Aromatic(Lit(true))"),
            }
            let fresh = BondConstraintsView {
                backing: BondConstraintsBacking::Bond(bond),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 0);
        });
    }

    #[rstest]
    fn test_bond_constraints_view_update() {
        Python::attach(|py| {
            let bond = Py::new(py, BondAst::from_inner(GraphIrBondAst::from_order(1))).unwrap();
            let view = BondConstraintsView {
                backing: BondConstraintsBacking::Bond(bond.clone_ref(py)),
            };
            let mut other = GraphIrBondConstraintsAst::new();
            other.set(GraphIrBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(
                true,
            )));
            other.set(GraphIrBondConstraintAst::ring_membership(
                GraphIrRingScope::All,
                2,
            ));
            view.update(
                py,
                BondConstraintsUpdate::Container(
                    Py::new(py, BondConstraintsAst::from_inner(other)).unwrap(),
                ),
            )
            .unwrap();
            let fresh = BondConstraintsView {
                backing: BondConstraintsBacking::Bond(bond),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 2);
        });
    }

    #[rstest]
    fn test_bond_constraints_view_set_aromatic() {
        Python::attach(|py| {
            let bond = Py::new(py, BondAst::from_inner(GraphIrBondAst::from_order(1))).unwrap();
            let view = BondConstraintsView {
                backing: BondConstraintsBacking::Bond(bond.clone_ref(py)),
            };
            assert_eq!(
                view.aromatic(py).unwrap().to_rust(),
                GraphIrBooleanAst::Undetermined
            );
            view.set_aromatic(py, BooleanLike::Lit(true));
            let fresh = BondConstraintsView {
                backing: BondConstraintsBacking::Bond(bond),
            };
            assert_eq!(
                fresh.aromatic(py).unwrap().to_rust(),
                GraphIrBooleanAst::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_bond_ring_size_counts_value_backed() {
        Python::attach(|py| {
            let constraints = Py::new(py, BondConstraintsAst::new(py, vec![])).unwrap();
            let proxy = BondConstraintsAst::ring_size_count(constraints.clone_ref(py));
            proxy.__setitem__(py, 6, ValueLike::Lit(3));
            assert_eq!(
                proxy.__getitem__(py, 6).unwrap().unwrap().to_rust(py),
                GraphIrValueAst::Lit(3)
            );
            proxy.__delitem__(py, 6);
            assert!(proxy.__getitem__(py, 6).unwrap().is_none());
        });
    }

    #[rstest]
    fn test_bond_ring_size_counts_bond_backed() {
        Python::attach(|py| {
            let bond = Py::new(py, BondAst::from_inner(GraphIrBondAst::from_order(1))).unwrap();
            let view = BondConstraintsView {
                backing: BondConstraintsBacking::Bond(bond.clone_ref(py)),
            };
            view.ring_size_count(py)
                .__setitem__(py, 5, ValueLike::Lit(1));
            let fresh = BondConstraintsView {
                backing: BondConstraintsBacking::Bond(bond),
            };
            assert_eq!(
                fresh
                    .ring_size_count(py)
                    .__getitem__(py, 5)
                    .unwrap()
                    .unwrap()
                    .to_rust(py),
                GraphIrValueAst::Lit(1)
            );
        });
    }

    #[rstest]
    fn test_bond_ring_size_counts_len_iter_contains() {
        Python::attach(|py| {
            let constraints = Py::new(py, BondConstraintsAst::new(py, vec![])).unwrap();
            let proxy = BondConstraintsAst::ring_size_count(constraints.clone_ref(py));
            proxy.__setitem__(py, 6, ValueLike::Lit(3));
            proxy.__setitem__(py, 5, ValueLike::Lit(1));
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

    #[rstest]
    fn test_bond_view_order() {
        Python::attach(|py| {
            let view = BondView {
                owner: ethene(py),
                id: GraphIrBondId(0),
            };
            assert_eq!(view.id(), 0);
            assert_eq!(view.order(py).unwrap().to_rust(py), GraphIrValueAst::Lit(2));
        });
    }

    #[rstest]
    fn test_bond_view_atom_ids() {
        Python::attach(|py| {
            let view = BondView {
                owner: ethene(py),
                id: GraphIrBondId(0),
            };
            assert_eq!(view.atom_ids(py).unwrap(), (0, 1));
        });
    }

    #[rstest]
    fn test_bond_view_set_order() {
        Python::attach(|py| {
            let owner = ethene(py);
            let view = BondView {
                owner: owner.clone_ref(py),
                id: GraphIrBondId(0),
            };
            view.set_order(py, ValueLike::Lit(1));
            let fresh = BondView {
                owner,
                id: GraphIrBondId(0),
            };
            assert_eq!(
                fresh.order(py).unwrap().to_rust(py),
                GraphIrValueAst::Lit(1)
            );
        });
    }

    #[rstest]
    fn test_bond_view_set_charge() {
        Python::attach(|py| {
            let owner = ethene(py);
            let view = BondView {
                owner: owner.clone_ref(py),
                id: GraphIrBondId(0),
            };
            view.set_charge(py, ValueLike::Lit(-1));
            let fresh = BondView {
                owner,
                id: GraphIrBondId(0),
            };
            assert_eq!(
                fresh.charge(py).unwrap().to_rust(py),
                GraphIrValueAst::Lit(-1)
            );
        });
    }

    #[rstest]
    fn test_bond_view_constraints() {
        Python::attach(|py| {
            let view = BondView {
                owner: ethene(py),
                id: GraphIrBondId(0),
            };
            match view.constraints(py).backing {
                BondConstraintsBacking::Molecule { id, .. } => assert_eq!(id, GraphIrBondId(0)),
                _ => panic!("expected molecule-backed view"),
            }
        });
    }

    #[rstest]
    fn test_bond_constraints_view_set_molecule_backed() {
        Python::attach(|py| {
            let owner = ethene(py);
            let view = BondConstraintsView {
                backing: BondConstraintsBacking::Molecule {
                    owner: owner.clone_ref(py),
                    id: GraphIrBondId(0),
                },
            };
            let aromatic = into_py_variant(
                py,
                BondConstraintAst::from_rust(
                    py,
                    &GraphIrBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            view.set(py, aromatic);
            let fresh = BondConstraintsView {
                backing: BondConstraintsBacking::Molecule {
                    owner,
                    id: GraphIrBondId(0),
                },
            };
            assert_eq!(fresh.__len__(py).unwrap(), 1);
            assert_eq!(
                fresh.aromatic(py).unwrap().to_rust(),
                GraphIrBooleanAst::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_bond_ring_size_counts_molecule_backed() {
        Python::attach(|py| {
            let owner = ethene(py);
            let view = BondConstraintsView {
                backing: BondConstraintsBacking::Molecule {
                    owner: owner.clone_ref(py),
                    id: GraphIrBondId(0),
                },
            };
            view.ring_size_count(py)
                .__setitem__(py, 6, ValueLike::Lit(1));
            let fresh = BondConstraintsView {
                backing: BondConstraintsBacking::Molecule {
                    owner,
                    id: GraphIrBondId(0),
                },
            };
            assert_eq!(
                fresh
                    .ring_size_count(py)
                    .__getitem__(py, 6)
                    .unwrap()
                    .unwrap()
                    .to_rust(py),
                GraphIrValueAst::Lit(1)
            );
        });
    }

    #[rstest]
    fn test_bond_views_len_and_getitem() {
        Python::attach(|py| {
            let views = BondViews { owner: ethene(py) };
            assert_eq!(views.__len__(py), 1);
            assert_eq!(views.__getitem__(py, 0).unwrap().id(), 0);
            assert_eq!(views.__getitem__(py, -1).unwrap().id(), 0);
            assert!(views.__getitem__(py, 5).is_err());
            assert!(views.__getitem__(py, -2).is_err());
        });
    }

    #[rstest]
    fn test_bond_views_setitem() {
        Python::attach(|py| {
            let owner = ethene(py);
            let views = BondViews {
                owner: owner.clone_ref(py),
            };
            let single = Py::new(py, BondAst::from_inner(GraphIrBondAst::from_order(1))).unwrap();
            views.__setitem__(py, 0, single.bind(py).borrow()).unwrap();
            let view = views.__getitem__(py, 0).unwrap();
            // value replaced, endpoints preserved
            assert_eq!(view.order(py).unwrap().to_rust(py), GraphIrValueAst::Lit(1));
            assert_eq!(view.atom_ids(py).unwrap(), (0, 1));
        });
    }

    #[rstest]
    fn test_bond_views_setitem_error() {
        Python::attach(|py| {
            let views = BondViews { owner: ethene(py) };
            let single = Py::new(py, BondAst::from_inner(GraphIrBondAst::from_order(1))).unwrap();
            assert!(views.__setitem__(py, 5, single.bind(py).borrow()).is_err());
        });
    }

    #[rstest]
    fn test_bond_views_of() {
        Python::attach(|py| {
            // three carbons, one bond 0–1; atom 2 is isolated
            let molecule = GraphIrMoleculeAst::from_entries(GraphIrMoleculeEntries {
                atoms: vec![
                    GraphIrAtomAst::from_element(ChemElement::C),
                    GraphIrAtomAst::from_element(ChemElement::C),
                    GraphIrAtomAst::from_element(ChemElement::C),
                ],
                bonds: vec![(
                    GraphIrAtomId(0),
                    GraphIrAtomId(1),
                    GraphIrBondAst::from_order(1),
                )],
                ..Default::default()
            });
            let owner = Py::new(py, MoleculeAst::from_rust(molecule)).unwrap();
            let views = BondViews { owner };
            assert_eq!(views.of(py, 0, 1).unwrap().id(), 0);
            assert_eq!(views.of(py, 1, 0).unwrap().id(), 0);
            assert!(views.of(py, 1, 2).is_none());
        });
    }
}
