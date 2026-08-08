//! Owned dative-bond ASTs and molecule-backed dative-bond views.

use std::str::FromStr;
use std::vec::IntoIter;

use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use umol_graph_ir::ir::{
    AtomId as GraphIrAtomId, DativeBondAst as GraphIrDativeBondAst,
    DativeBondId as GraphIrDativeBondId, DativeBondUpdate as GraphIrDativeBondUpdate,
    DativeBondView as GraphIrDativeBondView, MoleculeAst as GraphIrMoleculeAst,
};

use crate::constraint::dative::{
    dative_bond_constraints_asdict, DativeBondConstraintsAst, DativeBondConstraintsBacking,
    DativeBondConstraintsLike, DativeBondConstraintsView,
};
#[cfg(test)]
use crate::constraint::dative::{
    DativeBondConstraintAst, DativeBondConstraintKey, DativeBondConstraintsUpdate,
};
use crate::convert::hash_rust;
use crate::error::parse_error;
use crate::lattice::impl_py_lattice;
use crate::molecule::MoleculeAst;
use crate::value::{ValueAst, ValueLike};

/// Attribute updates for a dative bond.
#[pyclass(frozen, skip_from_py_object)]
#[derive(Clone)]
pub struct DativeBondUpdate(GraphIrDativeBondUpdate);

#[pymethods]
impl DativeBondUpdate {
    #[new]
    #[pyo3(signature = (*, order=None, constraints=None))]
    fn new(
        py: Python<'_>,
        order: Option<ValueLike>,
        constraints: Option<Py<DativeBondConstraintsAst>>,
    ) -> Self {
        Self::from_rust(&GraphIrDativeBondUpdate {
            order: order.map(|value| value.to_rust(py)),
            constraints: constraints
                .map(|value| value.bind(py).borrow().inner().clone())
                .unwrap_or_default(),
        })
    }

    /// Parse a dative-bond-update DSL string into a `DativeBondUpdate`.
    #[staticmethod]
    fn parse(s: &str) -> PyResult<Self> {
        GraphIrDativeBondUpdate::from_str(s)
            .map(Self)
            .map_err(parse_error)
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("DativeBondUpdate.parse('{}')", self.0)
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
    fn constraints(&self) -> DativeBondConstraintsAst {
        DativeBondConstraintsAst::from_inner(self.0.constraints.clone())
    }
}

impl DativeBondUpdate {
    pub(crate) fn from_rust(update: &GraphIrDativeBondUpdate) -> Self {
        Self(update.clone())
    }

    pub(crate) fn to_rust(&self) -> GraphIrDativeBondUpdate {
        self.0.clone()
    }
}

/// A dative bond: order and bond-scope constraints.
#[pyclass(eq)]
#[derive(PartialEq)]
pub struct DativeBondAst(GraphIrDativeBondAst);

#[pymethods]
impl DativeBondAst {
    /// Construct from an order — an `int` or a `ValueAst` expression — optionally
    /// setting constraints.
    #[new]
    #[pyo3(signature = (order, *, constraints=None))]
    fn new(
        py: Python<'_>,
        order: ValueLike,
        constraints: Option<Py<DativeBondConstraintsAst>>,
    ) -> Self {
        let mut bond = GraphIrDativeBondAst::new(order.to_rust(py));
        if let Some(constraints) = constraints {
            bond.constraints = constraints.bind(py).borrow().inner().clone();
        }
        DativeBondAst(bond)
    }

    /// Parse a dative-bond-DSL string (e.g. `"1#R(6)"`) into a `DativeBondAst`.
    #[staticmethod]
    fn parse(s: &str) -> PyResult<Self> {
        GraphIrDativeBondAst::from_str(s)
            .map(Self)
            .map_err(parse_error)
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("DativeBondAst.parse('{}')", self.0)
    }

    #[getter]
    fn order(&self, py: Python<'_>) -> PyResult<ValueAst> {
        ValueAst::from_rust(py, &self.0.order)
    }

    #[setter]
    fn set_order(&mut self, py: Python<'_>, value: ValueLike) {
        self.0.order = value.to_rust(py);
    }

    /// The dative bond's constraints as a live handle onto this bond: reads borrow
    /// the current state, mutators write through to the bond in place.
    #[getter]
    fn constraints(slf: Py<Self>) -> DativeBondConstraintsView {
        DativeBondConstraintsView {
            backing: DativeBondConstraintsBacking::DativeBond(slf),
        }
    }

    /// Replace the whole constraint set (wipe-and-set) from a value container or
    /// a live view.
    #[setter]
    fn set_constraints(
        slf: Py<Self>,
        py: Python<'_>,
        value: DativeBondConstraintsLike,
    ) -> PyResult<()> {
        let snapshot = value.to_rust(py)?;
        slf.borrow_mut(py).0.constraints = snapshot;
        Ok(())
    }

    /// The fields as a dict keyed by field name; values are Python objects.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("order", self.order(py)?)?;
        dict.set_item(
            "constraints",
            dative_bond_constraints_asdict(py, &self.0.constraints)?,
        )?;
        Ok(dict)
    }
}

impl DativeBondAst {
    /// The wrapped AST bond — read access for the bond-backed constraints view.
    pub(crate) fn inner(&self) -> &GraphIrDativeBondAst {
        &self.0
    }

    /// Mutable access to the wrapped AST bond — write access for the bond-backed
    /// constraints view.
    pub(crate) fn inner_mut(&mut self) -> &mut GraphIrDativeBondAst {
        &mut self.0
    }

    /// Wrap an owned Rust dative-bond AST.
    pub(crate) fn from_inner(bond: GraphIrDativeBondAst) -> Self {
        DativeBondAst(bond)
    }
}

impl_py_lattice!(
    DativeBondAst,
    GraphIrDativeBondAst,
    |value: &DativeBondAst, _py: Python<'_>| -> PyResult<GraphIrDativeBondAst> {
        Ok(value.inner().clone())
    },
    |_py: Python<'_>, value: GraphIrDativeBondAst| -> PyResult<DativeBondAst> {
        Ok(DativeBondAst::from_inner(value))
    }
);

/// A view of one dative bond within a molecule: a handle to the molecule plus the
/// bond's index. Field reads rebuild the transient Rust view; the molecule is never
/// copied. The acceptor and donor atom indices are read-only topology; the order
/// and constraints are the mutable bond value.
#[pyclass]
pub struct DativeBondView {
    owner: Py<MoleculeAst>,
    id: GraphIrDativeBondId,
}

impl DativeBondView {
    fn dative_bond<'a>(
        &self,
        molecule: &'a GraphIrMoleculeAst,
    ) -> PyResult<GraphIrDativeBondView<'a>> {
        molecule
            .dative_bonds()
            .get(self.id)
            .ok_or_else(|| PyIndexError::new_err("dative bond id out of range"))
    }
}

#[pymethods]
impl DativeBondView {
    #[getter]
    fn id(&self) -> u32 {
        self.id.0
    }

    /// The acceptor atom index (read-only — participants are topology, not part of
    /// the bond value).
    #[getter]
    fn acceptor(&self, py: Python<'_>) -> PyResult<u32> {
        let molecule = self.owner.bind(py).borrow();
        Ok(self.dative_bond(molecule.inner())?.acceptor_id().0)
    }

    /// The donor atom indices (read-only).
    #[getter]
    fn donors(&self, py: Python<'_>) -> PyResult<Vec<u32>> {
        let molecule = self.owner.bind(py).borrow();
        Ok(self
            .dative_bond(molecule.inner())?
            .donor_ids()
            .map(|donor| donor.0)
            .collect())
    }

    /// All atom indices incident to this dative bond — the donors followed by the
    /// acceptor (read-only).
    #[getter]
    fn atom_ids<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let molecule = self.owner.bind(py).borrow();
        let atom_ids: Vec<u32> = self
            .dative_bond(molecule.inner())?
            .atom_ids()
            .map(|atom| atom.0)
            .collect();
        PyTuple::new(py, atom_ids)
    }

    fn __repr__(&self) -> String {
        format!("DativeBondView(id={})", self.id.0)
    }

    #[getter]
    fn order(&self, py: Python<'_>) -> PyResult<ValueAst> {
        let molecule = self.owner.bind(py).borrow();
        ValueAst::from_rust(py, &self.dative_bond(molecule.inner())?.ast.order)
    }

    #[setter]
    fn set_order(&self, py: Python<'_>, value: ValueLike) {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .dative_bond_mut(self.id)
            .ast
            .order = value.to_rust(py);
    }

    /// The dative bond's constraints as a live handle onto the molecule: reads borrow
    /// the current state, mutators write through to the bond in place.
    #[getter]
    fn constraints(&self, py: Python<'_>) -> DativeBondConstraintsView {
        DativeBondConstraintsView {
            backing: DativeBondConstraintsBacking::Molecule {
                owner: self.owner.clone_ref(py),
                id: self.id,
            },
        }
    }

    /// Replace the whole constraint set of the backing bond in place (wipe-and-set)
    /// from a value container or a live view.
    #[setter]
    fn set_constraints(&self, py: Python<'_>, value: DativeBondConstraintsLike) -> PyResult<()> {
        self.owner
            .borrow_mut(py)
            .inner_mut()
            .dative_bond_mut(self.id)
            .ast
            .constraints = value.to_rust(py)?;
        Ok(())
    }

    /// The value fields as a dict keyed by field name; values are Python objects —
    /// symmetric with `DativeBondAst.asdict`, read through the view.
    fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let molecule = self.owner.bind(py).borrow();
        let bond = self.dative_bond(molecule.inner())?.ast;
        let dict = PyDict::new(py);
        dict.set_item("order", ValueAst::from_rust(py, &bond.order)?)?;
        dict.set_item(
            "constraints",
            dative_bond_constraints_asdict(py, &bond.constraints)?,
        )?;
        Ok(dict)
    }
}

/// Resolve a possibly-negative Python index (negative counts from the end) into an
/// existing dative bond id, or `IndexError`. `DativeBondId` is `RelationId`-backed
/// but contiguous for fresh molecules, so integer positions address it directly.
fn resolve_dative_bond_index(
    molecule: &GraphIrMoleculeAst,
    index: isize,
) -> PyResult<GraphIrDativeBondId> {
    let count = molecule.dative_bonds().count();
    let resolved = if index < 0 {
        index + count as isize
    } else {
        index
    };
    if resolved < 0 {
        return Err(PyIndexError::new_err("dative bond id out of range"));
    }
    let id = GraphIrDativeBondId(resolved as u32);
    if molecule.dative_bonds().contains(id) {
        Ok(id)
    } else {
        Err(PyIndexError::new_err("dative bond id out of range"))
    }
}

/// The dative bonds of a molecule, indexed by integer position.
#[pyclass]
pub struct DativeBondViews {
    owner: Py<MoleculeAst>,
}

#[pymethods]
impl DativeBondViews {
    fn __len__(&self, py: Python<'_>) -> usize {
        self.owner.bind(py).borrow().inner().dative_bonds().count()
    }

    fn __repr__(&self, py: Python<'_>) -> String {
        format!(
            "DativeBondViews(len={})",
            self.owner.bind(py).borrow().inner().dative_bonds().count()
        )
    }

    fn __getitem__(&self, py: Python<'_>, index: isize) -> PyResult<DativeBondView> {
        let molecule = self.owner.bind(py).borrow();
        let id = resolve_dative_bond_index(molecule.inner(), index)?;
        Ok(DativeBondView {
            owner: self.owner.clone_ref(py),
            id,
        })
    }

    /// Replace the whole dative bond value at `index` in place (participants unchanged).
    fn __setitem__(
        &self,
        py: Python<'_>,
        index: isize,
        bond: PyRef<'_, DativeBondAst>,
    ) -> PyResult<()> {
        let mut molecule = self.owner.borrow_mut(py);
        let id = resolve_dative_bond_index(molecule.inner(), index)?;
        *molecule.inner_mut().dative_bond_mut(id).ast = bond.inner().clone();
        Ok(())
    }

    /// The dative bond with exactly this acceptor and donor set, or `None`.
    fn of(&self, py: Python<'_>, donors: Vec<u32>, acceptor: u32) -> Option<DativeBondView> {
        let molecule = self.owner.bind(py).borrow();
        let donor_ids: Vec<GraphIrAtomId> = donors.into_iter().map(GraphIrAtomId).collect();
        molecule
            .inner()
            .dative_bonds()
            .of_id(GraphIrAtomId(acceptor), &donor_ids)
            .map(|id| DativeBondView {
                owner: self.owner.clone_ref(py),
                id,
            })
    }

    /// The dative bonds incident on `atom` (as acceptor or donor).
    fn incident(&self, py: Python<'_>, atom: u32) -> Vec<DativeBondView> {
        let molecule = self.owner.bind(py).borrow();
        molecule
            .inner()
            .dative_bonds()
            .incident_ids(GraphIrAtomId(atom))
            .map(|id| DativeBondView {
                owner: self.owner.clone_ref(py),
                id,
            })
            .collect()
    }

    fn __iter__(&self, py: Python<'_>) -> DativeBondViewIter {
        let ids = self
            .owner
            .bind(py)
            .borrow()
            .inner()
            .dative_bonds()
            .ids()
            .collect::<Vec<_>>();
        DativeBondViewIter {
            owner: self.owner.clone_ref(py),
            ids: ids.into_iter(),
        }
    }
}

impl DativeBondViews {
    /// Build the dative-bond-views handle for `owner` (the `.dative_bonds` accessor).
    pub(crate) fn new(owner: Py<MoleculeAst>) -> DativeBondViews {
        DativeBondViews { owner }
    }
}

#[pyclass]
struct DativeBondViewIter {
    owner: Py<MoleculeAst>,
    ids: IntoIter<GraphIrDativeBondId>,
}

#[pymethods]
impl DativeBondViewIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> Option<DativeBondView> {
        self.ids.next().map(|id| DativeBondView {
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
        AtomAst as GraphIrAtomAst, BooleanAst as GraphIrBooleanAst,
        DativeBondConstraintAst as GraphIrDativeBondConstraintAst,
        DativeBondConstraintKey as GraphIrDativeBondConstraintKey,
        DativeBondConstraintsAst as GraphIrDativeBondConstraintsAst, MoleculeEntries,
        RingScope as GraphIrRingScope, ValueAst as GraphIrValueAst,
    };

    use super::*;
    use crate::boolean::BooleanLike;
    use crate::constraint::ring::RingScope;
    use crate::convert::into_py_variant;

    /// An ammonia-borane adduct: borane B (id 0) accepts from ammonia N (id 1),
    /// dative bond id 0 (acceptor B, donor N, order 1).
    fn ammonia_borane(py: Python<'_>) -> Py<MoleculeAst> {
        let molecule = GraphIrMoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![
                GraphIrAtomAst::from_element(ChemElement::B),
                GraphIrAtomAst::from_element(ChemElement::N),
            ],
            dative: vec![(
                vec![GraphIrAtomId(1)],
                GraphIrAtomId(0),
                GraphIrDativeBondAst::from_order(1),
            )],
            ..Default::default()
        });
        Py::new(py, MoleculeAst::from_rust(molecule)).unwrap()
    }

    #[rstest]
    #[case::single("1")]
    #[case::aromatic("1#a")]
    #[case::ring_size("1#R(6)")]
    fn test_dative_bond_ast_parse(#[case] dsl: &str) {
        let bond = DativeBondAst::parse(dsl).unwrap();
        assert_eq!(bond.__str__(), dsl);
        assert_eq!(bond.__repr__(), format!("DativeBondAst.parse('{dsl}')"));
    }

    #[rstest]
    fn test_dative_bond_ast_parse_error() {
        assert!(DativeBondAst::parse("x#").is_err());
    }

    #[rstest]
    fn test_dative_bond_ast_constraints() {
        let bond = DativeBondAst(GraphIrDativeBondAst::from_order(1).with_constraint(
            GraphIrDativeBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
        ));
        assert_eq!(bond.inner().constraints.len(), 1);
    }

    #[rstest]
    fn test_dative_bond_ast_set_constraints_from_view() {
        Python::attach(|py| {
            let src = Py::new(
                py,
                DativeBondAst::from_inner(GraphIrDativeBondAst::from_order(1).with_constraint(
                    GraphIrDativeBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )),
            )
            .unwrap();
            let view = Py::new(
                py,
                DativeBondConstraintsView {
                    backing: DativeBondConstraintsBacking::DativeBond(src),
                },
            )
            .unwrap();
            let dst = Py::new(
                py,
                DativeBondAst::from_inner(GraphIrDativeBondAst::from_order(2)),
            )
            .unwrap();
            DativeBondAst::set_constraints(
                dst.clone_ref(py),
                py,
                DativeBondConstraintsLike::View(view),
            )
            .unwrap();
            assert_eq!(
                dst.bind(py).borrow().inner().constraints.aromatic(),
                GraphIrBooleanAst::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_dative_bond_view_order() {
        Python::attach(|py| {
            let view = DativeBondView {
                owner: ammonia_borane(py),
                id: GraphIrDativeBondId(0),
            };
            assert_eq!(view.id(), 0);
            assert_eq!(view.order(py).unwrap().to_rust(py), GraphIrValueAst::Lit(1));
        });
    }

    #[rstest]
    fn test_dative_bond_view_participants() {
        Python::attach(|py| {
            let view = DativeBondView {
                owner: ammonia_borane(py),
                id: GraphIrDativeBondId(0),
            };
            assert_eq!(view.acceptor(py).unwrap(), 0);
            assert_eq!(view.donors(py).unwrap(), vec![1]);
            // atom_ids is donors-then-acceptor
            let atom_ids: Vec<u32> = view.atom_ids(py).unwrap().extract().unwrap();
            assert_eq!(atom_ids, vec![1, 0]);
        });
    }

    #[rstest]
    fn test_dative_bond_view_set_order() {
        Python::attach(|py| {
            let owner = ammonia_borane(py);
            let view = DativeBondView {
                owner: owner.clone_ref(py),
                id: GraphIrDativeBondId(0),
            };
            view.set_order(py, ValueLike::Lit(2));
            let fresh = DativeBondView {
                owner,
                id: GraphIrDativeBondId(0),
            };
            assert_eq!(
                fresh.order(py).unwrap().to_rust(py),
                GraphIrValueAst::Lit(2)
            );
        });
    }

    #[rstest]
    fn test_dative_bond_view_constraints() {
        Python::attach(|py| {
            let view = DativeBondView {
                owner: ammonia_borane(py),
                id: GraphIrDativeBondId(0),
            };
            match view.constraints(py).backing {
                DativeBondConstraintsBacking::Molecule { id, .. } => {
                    assert_eq!(id, GraphIrDativeBondId(0))
                }
                _ => panic!("expected molecule-backed view"),
            }
        });
    }

    #[rstest]
    fn test_dative_bond_views_len_and_getitem() {
        Python::attach(|py| {
            let views = DativeBondViews {
                owner: ammonia_borane(py),
            };
            assert_eq!(views.__len__(py), 1);
            assert_eq!(views.__getitem__(py, 0).unwrap().id(), 0);
            assert_eq!(views.__getitem__(py, -1).unwrap().id(), 0);
            assert!(views.__getitem__(py, 5).is_err());
            assert!(views.__getitem__(py, -2).is_err());
        });
    }

    #[rstest]
    fn test_dative_bond_views_setitem() {
        Python::attach(|py| {
            let owner = ammonia_borane(py);
            let views = DativeBondViews {
                owner: owner.clone_ref(py),
            };
            let single = Py::new(
                py,
                DativeBondAst::from_inner(GraphIrDativeBondAst::from_order(2)),
            )
            .unwrap();
            views.__setitem__(py, 0, single.bind(py).borrow()).unwrap();
            let view = views.__getitem__(py, 0).unwrap();
            // value replaced, participants preserved
            assert_eq!(view.order(py).unwrap().to_rust(py), GraphIrValueAst::Lit(2));
            assert_eq!(view.acceptor(py).unwrap(), 0);
            assert_eq!(view.donors(py).unwrap(), vec![1]);
        });
    }

    #[rstest]
    fn test_dative_bond_views_setitem_error() {
        Python::attach(|py| {
            let views = DativeBondViews {
                owner: ammonia_borane(py),
            };
            let single = Py::new(
                py,
                DativeBondAst::from_inner(GraphIrDativeBondAst::from_order(2)),
            )
            .unwrap();
            assert!(views.__setitem__(py, 5, single.bind(py).borrow()).is_err());
        });
    }

    #[rstest]
    fn test_dative_bond_views_of() {
        Python::attach(|py| {
            let views = DativeBondViews {
                owner: ammonia_borane(py),
            };
            // acceptor B(0), donor N(1)
            assert_eq!(views.of(py, vec![1], 0).unwrap().id(), 0);
            // roles swapped: no such dative bond
            assert!(views.of(py, vec![0], 1).is_none());
        });
    }

    #[rstest]
    fn test_dative_bond_views_incident() {
        Python::attach(|py| {
            // B(0) accepts from N(1); C(2) isolated
            let molecule = GraphIrMoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![
                    GraphIrAtomAst::from_element(ChemElement::B),
                    GraphIrAtomAst::from_element(ChemElement::N),
                    GraphIrAtomAst::from_element(ChemElement::C),
                ],
                dative: vec![(
                    vec![GraphIrAtomId(1)],
                    GraphIrAtomId(0),
                    GraphIrDativeBondAst::from_order(1),
                )],
                ..Default::default()
            });
            let views = DativeBondViews {
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
            assert_eq!(
                views
                    .incident(py, 1)
                    .iter()
                    .map(|v| v.id())
                    .collect::<Vec<_>>(),
                vec![0]
            );
            assert!(views.incident(py, 2).is_empty());
        });
    }

    #[rstest]
    #[case(GraphIrDativeBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)))]
    #[case(GraphIrDativeBondConstraintAst::ring_membership(GraphIrRingScope::All, 2))]
    #[case(GraphIrDativeBondConstraintAst::ring_membership(GraphIrRingScope::Size(6), 1))]
    fn test_dative_bond_constraint_ast_roundtrip(#[case] ast: GraphIrDativeBondConstraintAst) {
        Python::attach(|py| {
            assert_eq!(
                DativeBondConstraintAst::from_rust(py, &ast)
                    .unwrap()
                    .to_rust(py),
                ast
            );
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_len_contains() {
        Python::attach(|py| {
            let aromatic = into_py_variant(
                py,
                DativeBondConstraintAst::from_rust(
                    py,
                    &GraphIrDativeBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let ring = into_py_variant(
                py,
                DativeBondConstraintAst::from_rust(
                    py,
                    &GraphIrDativeBondConstraintAst::ring_membership(GraphIrRingScope::All, 2),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = DativeBondConstraintsAst::new(py, vec![aromatic, ring]);
            assert_eq!(constraints.__len__(), 2);
            assert!(constraints.__contains__(
                py,
                into_py_variant(py, DativeBondConstraintKey::Aromatic()).unwrap()
            ));
            assert!(!constraints.__contains__(
                py,
                into_py_variant(
                    py,
                    DativeBondConstraintKey::RingMembership(
                        into_py_variant(py, RingScope::Size(5)).unwrap()
                    ),
                )
                .unwrap()
            ));
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_keys_values_items() {
        Python::attach(|py| {
            let aromatic = into_py_variant(
                py,
                DativeBondConstraintAst::from_rust(
                    py,
                    &GraphIrDativeBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let ring = into_py_variant(
                py,
                DativeBondConstraintAst::from_rust(
                    py,
                    &GraphIrDativeBondConstraintAst::ring_membership(GraphIrRingScope::All, 2),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = DativeBondConstraintsAst::new(py, vec![aromatic, ring]);

            let mut keys = constraints.__iter__(py).unwrap();
            assert_eq!(
                keys.__next__().unwrap().bind(py).borrow().to_rust(py),
                GraphIrDativeBondConstraintKey::Aromatic
            );
            assert_eq!(
                keys.__next__().unwrap().bind(py).borrow().to_rust(py),
                GraphIrDativeBondConstraintKey::RingMembership(GraphIrRingScope::All)
            );
            assert!(keys.__next__().is_none());

            let mut values = constraints.values(py).unwrap();
            assert_eq!(
                values.__next__().unwrap().bind(py).borrow().to_rust(py),
                GraphIrDativeBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true))
            );

            let mut items = constraints.items(py).unwrap();
            let (key, value) = items.__next__().unwrap();
            assert_eq!(
                key.bind(py).borrow().to_rust(py),
                GraphIrDativeBondConstraintKey::Aromatic
            );
            assert_eq!(
                value.bind(py).borrow().to_rust(py),
                GraphIrDativeBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true))
            );
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_get() {
        Python::attach(|py| {
            let aromatic = into_py_variant(
                py,
                DativeBondConstraintAst::from_rust(
                    py,
                    &GraphIrDativeBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = DativeBondConstraintsAst::new(py, vec![aromatic]);
            let present = constraints
                .get(
                    py,
                    into_py_variant(py, DativeBondConstraintKey::Aromatic()).unwrap(),
                    None,
                )
                .unwrap();
            let expected = into_py_variant(
                py,
                DativeBondConstraintAst::from_rust(
                    py,
                    &GraphIrDativeBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap()
            .into_any();
            assert!(present.bind(py).eq(expected.bind(py)).unwrap());
            let sentinel_key = into_py_variant(
                py,
                DativeBondConstraintKey::RingMembership(
                    into_py_variant(py, RingScope::All()).unwrap(),
                ),
            )
            .unwrap();
            let absent = constraints
                .get(py, sentinel_key.clone_ref(py), None)
                .unwrap();
            assert!(absent.bind(py).is_none());
            let sentinel = sentinel_key.into_any();
            let defaulted = constraints
                .get(
                    py,
                    into_py_variant(
                        py,
                        DativeBondConstraintKey::RingMembership(
                            into_py_variant(py, RingScope::All()).unwrap(),
                        ),
                    )
                    .unwrap(),
                    Some(sentinel.clone_ref(py)),
                )
                .unwrap();
            assert_eq!(defaulted.as_ptr(), sentinel.as_ptr());
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_aromatic() {
        Python::attach(|py| {
            let empty = DativeBondConstraintsAst::new(py, vec![]);
            assert_eq!(empty.aromatic().to_rust(), GraphIrBooleanAst::Undetermined);
            assert!(empty.ring_count(py).unwrap().is_none());
            let aromatic = into_py_variant(
                py,
                DativeBondConstraintAst::from_rust(
                    py,
                    &GraphIrDativeBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints = DativeBondConstraintsAst::new(py, vec![aromatic]);
            assert_eq!(
                constraints.aromatic().to_rust(),
                GraphIrBooleanAst::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_ring_size_count() {
        Python::attach(|py| {
            let membership = into_py_variant(
                py,
                DativeBondConstraintAst::from_rust(
                    py,
                    &GraphIrDativeBondConstraintAst::ring_membership(GraphIrRingScope::Size(6), 1),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints =
                Py::new(py, DativeBondConstraintsAst::new(py, vec![membership])).unwrap();
            let proxy = DativeBondConstraintsAst::ring_size_count(constraints.clone_ref(py));
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
    fn test_dative_bond_constraints_ast_set() {
        Python::attach(|py| {
            let mut constraints = DativeBondConstraintsAst::new(py, vec![]);
            let aromatic = into_py_variant(
                py,
                DativeBondConstraintAst::from_rust(
                    py,
                    &GraphIrDativeBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
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
    fn test_dative_bond_constraints_ast_pop() {
        Python::attach(|py| {
            let aromatic = into_py_variant(
                py,
                DativeBondConstraintAst::from_rust(
                    py,
                    &GraphIrDativeBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let mut constraints = DativeBondConstraintsAst::new(py, vec![aromatic]);
            let removed = constraints
                .pop(
                    py,
                    into_py_variant(py, DativeBondConstraintKey::Aromatic()).unwrap(),
                )
                .unwrap();
            match removed {
                Some(DativeBondConstraintAst::Aromatic(b)) => {
                    assert_eq!(b.bind(py).borrow().to_rust(), GraphIrBooleanAst::Lit(true))
                }
                _ => panic!("expected removed Aromatic(Lit(true))"),
            }
            assert_eq!(constraints.__len__(), 0);
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_update() {
        Python::attach(|py| {
            let constraints = Py::new(py, DativeBondConstraintsAst::new(py, vec![])).unwrap();
            let mut other = GraphIrDativeBondConstraintsAst::new();
            other.set(GraphIrDativeBondConstraintAst::aromatic(
                GraphIrBooleanAst::Lit(true),
            ));
            other.set(GraphIrDativeBondConstraintAst::ring_membership(
                GraphIrRingScope::All,
                2,
            ));
            DativeBondConstraintsAst::update(
                constraints.clone_ref(py),
                py,
                DativeBondConstraintsUpdate::Container(
                    Py::new(py, DativeBondConstraintsAst::from_inner(other)).unwrap(),
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
    fn test_dative_bond_constraints_ast_update_entries() {
        Python::attach(|py| {
            let constraints = Py::new(py, DativeBondConstraintsAst::new(py, vec![])).unwrap();
            let aromatic = into_py_variant(
                py,
                DativeBondConstraintAst::from_rust(
                    py,
                    &GraphIrDativeBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let ring = into_py_variant(
                py,
                DativeBondConstraintAst::from_rust(
                    py,
                    &GraphIrDativeBondConstraintAst::ring_membership(GraphIrRingScope::All, 2),
                )
                .unwrap(),
            )
            .unwrap();
            DativeBondConstraintsAst::update(
                constraints.clone_ref(py),
                py,
                DativeBondConstraintsUpdate::Entries(vec![aromatic, ring]),
            )
            .unwrap();
            assert_eq!(constraints.bind(py).borrow().__len__(), 2);
        });
    }

    /// Regression: a container updating itself resolves `other` before the write borrow,
    /// so it is an idempotent no-op, not a RefCell double-borrow panic.
    #[rstest]
    fn test_dative_bond_constraints_ast_update_self() {
        Python::attach(|py| {
            let aromatic = into_py_variant(
                py,
                DativeBondConstraintAst::from_rust(
                    py,
                    &GraphIrDativeBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            let constraints =
                Py::new(py, DativeBondConstraintsAst::new(py, vec![aromatic])).unwrap();
            DativeBondConstraintsAst::update(
                constraints.clone_ref(py),
                py,
                DativeBondConstraintsUpdate::Container(constraints.clone_ref(py)),
            )
            .unwrap();
            assert_eq!(
                constraints.bind(py).borrow().aromatic().to_rust(),
                GraphIrBooleanAst::Lit(true)
            );
        });
    }

    /// Regression: assigning a bond's own constraints view back to it snapshots before
    /// the write borrow, so it is a no-op, not a double-borrow panic.
    #[rstest]
    fn test_dative_bond_ast_set_constraints_self() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                DativeBondAst::from_inner(GraphIrDativeBondAst::from_order(1).with_constraint(
                    GraphIrDativeBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )),
            )
            .unwrap();
            let own_view = Py::new(
                py,
                DativeBondConstraintsView {
                    backing: DativeBondConstraintsBacking::DativeBond(bond.clone_ref(py)),
                },
            )
            .unwrap();
            DativeBondAst::set_constraints(
                bond.clone_ref(py),
                py,
                DativeBondConstraintsLike::View(own_view),
            )
            .unwrap();
            assert_eq!(
                bond.bind(py).borrow().inner().constraints.aromatic(),
                GraphIrBooleanAst::Lit(true)
            );
        });
    }

    /// Regression: a view updating from a view over the same bond resolves `other`
    /// before the write borrow, so it is an idempotent no-op, not a double-borrow panic.
    #[rstest]
    fn test_dative_bond_constraints_view_update_self() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                DativeBondAst::from_inner(GraphIrDativeBondAst::from_order(1).with_constraint(
                    GraphIrDativeBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )),
            )
            .unwrap();
            let view = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::DativeBond(bond.clone_ref(py)),
            };
            let other = Py::new(
                py,
                DativeBondConstraintsView {
                    backing: DativeBondConstraintsBacking::DativeBond(bond.clone_ref(py)),
                },
            )
            .unwrap();
            view.update(py, DativeBondConstraintsUpdate::View(other))
                .unwrap();
            assert_eq!(
                bond.bind(py).borrow().inner().constraints.aromatic(),
                GraphIrBooleanAst::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_set_aromatic() {
        Python::attach(|py| {
            let mut constraints = DativeBondConstraintsAst::new(py, vec![]);
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
    fn test_dative_bond_constraints_ast_set_ring_count() {
        Python::attach(|py| {
            let mut constraints = DativeBondConstraintsAst::new(py, vec![]);
            constraints.set_ring_count(py, ValueLike::Lit(2));
            assert_eq!(
                constraints.ring_count(py).unwrap().unwrap().to_rust(py),
                GraphIrValueAst::Lit(2)
            );
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_getitem_error() {
        Python::attach(|py| {
            let constraints = DativeBondConstraintsAst::new(py, vec![]);
            let key = into_py_variant(py, DativeBondConstraintKey::Aromatic()).unwrap();
            assert!(constraints.__getitem__(py, key).is_err());
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_ast_delitem_error() {
        Python::attach(|py| {
            let mut constraints = DativeBondConstraintsAst::new(py, vec![]);
            let key = into_py_variant(py, DativeBondConstraintKey::Aromatic()).unwrap();
            assert!(constraints.__delitem__(py, key).is_err());
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_view_set() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                DativeBondAst::from_inner(GraphIrDativeBondAst::from_order(1)),
            )
            .unwrap();
            let view = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::DativeBond(bond.clone_ref(py)),
            };
            let aromatic = into_py_variant(
                py,
                DativeBondConstraintAst::from_rust(
                    py,
                    &GraphIrDativeBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            view.set(py, aromatic);
            // a fresh view proves the write hit the standalone bond, not a copy
            let fresh = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::DativeBond(bond),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 1);
            match fresh
                .__getitem__(
                    py,
                    into_py_variant(py, DativeBondConstraintKey::Aromatic()).unwrap(),
                )
                .unwrap()
            {
                DativeBondConstraintAst::Aromatic(b) => {
                    assert_eq!(b.bind(py).borrow().to_rust(), GraphIrBooleanAst::Lit(true))
                }
                _ => panic!("expected Aromatic(Lit(true))"),
            }
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_view_pop() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                DativeBondAst::from_inner(GraphIrDativeBondAst::from_order(1).with_constraint(
                    GraphIrDativeBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )),
            )
            .unwrap();
            let view = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::DativeBond(bond.clone_ref(py)),
            };
            let removed = view
                .pop(
                    py,
                    into_py_variant(py, DativeBondConstraintKey::Aromatic()).unwrap(),
                )
                .unwrap();
            match removed {
                Some(DativeBondConstraintAst::Aromatic(b)) => {
                    assert_eq!(b.bind(py).borrow().to_rust(), GraphIrBooleanAst::Lit(true))
                }
                _ => panic!("expected removed Aromatic(Lit(true))"),
            }
            let fresh = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::DativeBond(bond),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 0);
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_view_update() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                DativeBondAst::from_inner(GraphIrDativeBondAst::from_order(1)),
            )
            .unwrap();
            let view = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::DativeBond(bond.clone_ref(py)),
            };
            let mut other = GraphIrDativeBondConstraintsAst::new();
            other.set(GraphIrDativeBondConstraintAst::aromatic(
                GraphIrBooleanAst::Lit(true),
            ));
            other.set(GraphIrDativeBondConstraintAst::ring_membership(
                GraphIrRingScope::All,
                2,
            ));
            view.update(
                py,
                DativeBondConstraintsUpdate::Container(
                    Py::new(py, DativeBondConstraintsAst::from_inner(other)).unwrap(),
                ),
            )
            .unwrap();
            let fresh = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::DativeBond(bond),
            };
            assert_eq!(fresh.__len__(py).unwrap(), 2);
        });
    }

    #[rstest]
    fn test_dative_bond_constraints_view_set_aromatic() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                DativeBondAst::from_inner(GraphIrDativeBondAst::from_order(1)),
            )
            .unwrap();
            let view = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::DativeBond(bond.clone_ref(py)),
            };
            assert_eq!(
                view.aromatic(py).unwrap().to_rust(),
                GraphIrBooleanAst::Undetermined
            );
            view.set_aromatic(py, BooleanLike::Lit(true));
            let fresh = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::DativeBond(bond),
            };
            assert_eq!(
                fresh.aromatic(py).unwrap().to_rust(),
                GraphIrBooleanAst::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_dative_bond_ring_size_counts_value_backed() {
        Python::attach(|py| {
            let constraints = Py::new(py, DativeBondConstraintsAst::new(py, vec![])).unwrap();
            let proxy = DativeBondConstraintsAst::ring_size_count(constraints.clone_ref(py));
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
    fn test_dative_bond_ring_size_counts_bond_backed() {
        Python::attach(|py| {
            let bond = Py::new(
                py,
                DativeBondAst::from_inner(GraphIrDativeBondAst::from_order(1)),
            )
            .unwrap();
            let view = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::DativeBond(bond.clone_ref(py)),
            };
            view.ring_size_count(py)
                .__setitem__(py, 5, ValueLike::Lit(1));
            let fresh = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::DativeBond(bond),
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
    fn test_dative_bond_ring_size_counts_len_iter_contains() {
        Python::attach(|py| {
            let constraints = Py::new(py, DativeBondConstraintsAst::new(py, vec![])).unwrap();
            let proxy = DativeBondConstraintsAst::ring_size_count(constraints.clone_ref(py));
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
    fn test_dative_bond_constraints_view_set_molecule_backed() {
        Python::attach(|py| {
            let owner = ammonia_borane(py);
            let view = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::Molecule {
                    owner: owner.clone_ref(py),
                    id: GraphIrDativeBondId(0),
                },
            };
            let aromatic = into_py_variant(
                py,
                DativeBondConstraintAst::from_rust(
                    py,
                    &GraphIrDativeBondConstraintAst::aromatic(GraphIrBooleanAst::Lit(true)),
                )
                .unwrap(),
            )
            .unwrap();
            view.set(py, aromatic);
            let fresh = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::Molecule {
                    owner,
                    id: GraphIrDativeBondId(0),
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
    fn test_dative_bond_ring_size_counts_molecule_backed() {
        Python::attach(|py| {
            let owner = ammonia_borane(py);
            let view = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::Molecule {
                    owner: owner.clone_ref(py),
                    id: GraphIrDativeBondId(0),
                },
            };
            view.ring_size_count(py)
                .__setitem__(py, 6, ValueLike::Lit(1));
            let fresh = DativeBondConstraintsView {
                backing: DativeBondConstraintsBacking::Molecule {
                    owner,
                    id: GraphIrDativeBondId(0),
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
}
