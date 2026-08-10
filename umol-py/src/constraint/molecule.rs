//! Molecule-level constraint payloads matching `umol_graph_ir::ir::constraint`.

use std::vec::IntoIter;

use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use umol_graph_ir::ir::{
    AromaticSystemId as GraphIrAromaticSystemId, AtomId as GraphIrAtomId, BondId as GraphIrBondId,
    Constraint as GraphIrConstraint, Constraints as GraphIrConstraints,
    DativeBondId as GraphIrDativeBondId, MoleculeConstraint as GraphIrMoleculeConstraint,
    MulticenterBondId as GraphIrMulticenterBondId, NoncovalentBondId as GraphIrNoncovalentBondId,
    RelationalConstraint as GraphIrRelationalConstraint, StereoAtomId as GraphIrStereoAtomId,
    StereoBondId as GraphIrStereoBondId,
};

use super::aromatic::AromaticSystemConstraintForm;
use super::atom::AtomConstraintForm;
use super::bond::BondConstraintForm;
use super::dative::DativeBondConstraintForm;
use super::multicenter::MulticenterBondConstraintForm;
use super::noncovalent::NoncovalentBondConstraintForm;
use super::stereo::{StereoAtomConstraintForm, StereoBondConstraintForm};
use crate::convert::{into_py_variant, variant_repr};
use crate::lattice::impl_py_canonicalize;
use crate::molecule::MoleculeAst;
use crate::spin::UnpairedElectronsForm;
use crate::stereo::StereoKind;
use crate::value::NumForm;

/// A cross-entity molecule constraint covering dative bonds, aromatic systems,
/// multicenter bonds, noncovalent bonds, stereo atoms, and stereo bonds.
#[pyclass(frozen)]
pub enum RelationalConstraint {
    DativeBondDonors(u32, Vec<u32>),
    DativeBondDonor(u32, u32),
    DativeBondContainsAllDonors(u32, Vec<u32>),
    DativeBondAllDonors(u32, Py<AtomConstraintForm>),
    DativeBondAnyDonor(u32, Py<AtomConstraintForm>),
    DativeBondAcceptor(u32, u32),
    DativeBondAcceptorSatisfies(u32, Py<AtomConstraintForm>),
    DativeBondParallels(u32, u32),
    AromaticSystemAtoms(u32, Vec<u32>),
    AromaticSystemContains(u32, u32),
    AromaticSystemContainsAll(u32, Vec<u32>),
    AromaticSystemAllAtoms(u32, Py<AtomConstraintForm>),
    AromaticSystemAnyAtom(u32, Py<AtomConstraintForm>),
    MulticenterBondAtoms(u32, Vec<u32>),
    MulticenterBondContains(u32, u32),
    MulticenterBondContainsAll(u32, Vec<u32>),
    MulticenterBondAllAtoms(u32, Py<AtomConstraintForm>),
    MulticenterBondAnyAtom(u32, Py<AtomConstraintForm>),
    NoncovalentBondEnds(u32, [u32; 2]),
    NoncovalentBondContains(u32, u32),
    NoncovalentBondEndsSatisfy(u32, [Py<AtomConstraintForm>; 2]),
    StereoAtomSite(u32, u32),
    StereoAtomContains(u32, u32),
    StereoAtomLigands(u32, Vec<u32>),
    StereoAtomAllLigands(u32, Py<AtomConstraintForm>),
    StereoAtomAnyLigand(u32, Py<AtomConstraintForm>),
    StereoBondSite(u32, u32),
    StereoBondContains(u32, u32),
    StereoBondLigands(u32, Vec<u32>),
    StereoBondAllLigands(u32, Py<AtomConstraintForm>),
    StereoBondAnyLigand(u32, Py<AtomConstraintForm>),
}

/// A molecule-scope predicate over values or connectivity.
#[pyclass(frozen)]
pub enum MoleculeConstraint {
    ChargeSum(Option<Vec<u32>>, Py<NumForm>),
    #[pyo3(constructor = (atoms, unpaired_electrons))]
    UnpairedElectronCoupling {
        atoms: Option<Vec<u32>>,
        unpaired_electrons: Py<UnpairedElectronsForm>,
    },
    BondOrderSum(Option<Vec<u32>>, Py<NumForm>),
    Connected(Option<Vec<u32>>),
}

/// A recursive molecule-constraint tree containing entity leaves, aggregate
/// leaves, and Boolean combinators.
#[pyclass(frozen)]
pub enum Constraint {
    Atom(u32, Py<AtomConstraintForm>),
    Bond(u32, Py<BondConstraintForm>),
    DativeBond(u32, Py<DativeBondConstraintForm>),
    AromaticSystem(u32, Py<AromaticSystemConstraintForm>),
    MulticenterBond(u32, Py<MulticenterBondConstraintForm>),
    NoncovalentBond(u32, Py<NoncovalentBondConstraintForm>),
    StereoAtom(u32, StereoKind, Py<StereoAtomConstraintForm>),
    StereoBond(u32, StereoKind, Py<StereoBondConstraintForm>),
    Relational(Py<RelationalConstraint>),
    Molecule(Py<MoleculeConstraint>),
    And(Vec<Py<Constraint>>),
    Or(Vec<Py<Constraint>>),
    Not(Py<Constraint>),
}

#[pymethods]
impl Constraint {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            Self::Atom(_, _) => ("Atom", 2),
            Self::Bond(_, _) => ("Bond", 2),
            Self::DativeBond(_, _) => ("DativeBond", 2),
            Self::AromaticSystem(_, _) => ("AromaticSystem", 2),
            Self::MulticenterBond(_, _) => ("MulticenterBond", 2),
            Self::NoncovalentBond(_, _) => ("NoncovalentBond", 2),
            Self::StereoAtom(_, _, _) => ("StereoAtom", 3),
            Self::StereoBond(_, _, _) => ("StereoBond", 3),
            Self::Relational(_) => ("Relational", 1),
            Self::Molecule(_) => ("Molecule", 1),
            Self::And(_) => ("And", 1),
            Self::Or(_) => ("Or", 1),
            Self::Not(_) => ("Not", 1),
        };
        variant_repr(slf.bind(py).as_any(), "Constraint", variant, arity)
    }
}

impl_py_canonicalize!(
    Constraint,
    GraphIrConstraint,
    |value: &Constraint, py: Python<'_>| -> PyResult<GraphIrConstraint> { Ok(value.to_rust(py)) },
    |py: Python<'_>, value: GraphIrConstraint| -> PyResult<Constraint> {
        Constraint::from_rust(py, &value)
    }
);

impl Constraint {
    pub(crate) fn from_rust(py: Python<'_>, constraint: &GraphIrConstraint) -> PyResult<Self> {
        Ok(match constraint {
            GraphIrConstraint::Atom(id, child) => Self::Atom(
                id.0,
                into_py_variant(py, AtomConstraintForm::from_rust(py, child)?)?,
            ),
            GraphIrConstraint::Bond(id, child) => Self::Bond(
                id.0,
                into_py_variant(py, BondConstraintForm::from_rust(py, child)?)?,
            ),
            GraphIrConstraint::DativeBond(id, child) => Self::DativeBond(
                id.0,
                into_py_variant(py, DativeBondConstraintForm::from_rust(py, child)?)?,
            ),
            GraphIrConstraint::AromaticSystem(id, child) => Self::AromaticSystem(
                id.0,
                into_py_variant(py, AromaticSystemConstraintForm::from_rust(py, child)?)?,
            ),
            GraphIrConstraint::MulticenterBond(id, child) => Self::MulticenterBond(
                id.0,
                into_py_variant(py, MulticenterBondConstraintForm::from_rust(py, child)?)?,
            ),
            GraphIrConstraint::NoncovalentBond(id, child) => Self::NoncovalentBond(
                id.0,
                into_py_variant(py, NoncovalentBondConstraintForm::from_rust(py, child)?)?,
            ),
            GraphIrConstraint::StereoAtom(id, kind, child) => Self::StereoAtom(
                id.0,
                StereoKind::from_rust(*kind),
                into_py_variant(py, StereoAtomConstraintForm::from_rust(py, child)?)?,
            ),
            GraphIrConstraint::StereoBond(id, kind, child) => Self::StereoBond(
                id.0,
                StereoKind::from_rust(*kind),
                into_py_variant(py, StereoBondConstraintForm::from_rust(py, child)?)?,
            ),
            GraphIrConstraint::Relational(child) => Self::Relational(into_py_variant(
                py,
                RelationalConstraint::from_rust(py, child)?,
            )?),
            GraphIrConstraint::Molecule(child) => Self::Molecule(into_py_variant(
                py,
                MoleculeConstraint::from_rust(py, child)?,
            )?),
            GraphIrConstraint::And(children) => Self::And(
                children
                    .iter()
                    .map(|child| into_py_variant(py, Self::from_rust(py, child)?))
                    .collect::<PyResult<_>>()?,
            ),
            GraphIrConstraint::Or(children) => Self::Or(
                children
                    .iter()
                    .map(|child| into_py_variant(py, Self::from_rust(py, child)?))
                    .collect::<PyResult<_>>()?,
            ),
            GraphIrConstraint::Not(child) => {
                Self::Not(into_py_variant(py, Self::from_rust(py, child)?)?)
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrConstraint {
        match self {
            Self::Atom(id, child) => {
                GraphIrConstraint::Atom(GraphIrAtomId(*id), child.bind(py).borrow().to_rust(py))
            }
            Self::Bond(id, child) => {
                GraphIrConstraint::Bond(GraphIrBondId(*id), child.bind(py).borrow().to_rust(py))
            }
            Self::DativeBond(id, child) => GraphIrConstraint::DativeBond(
                GraphIrDativeBondId(*id),
                child.bind(py).borrow().to_rust(py),
            ),
            Self::AromaticSystem(id, child) => GraphIrConstraint::AromaticSystem(
                GraphIrAromaticSystemId(*id),
                child.bind(py).borrow().to_rust(py),
            ),
            Self::MulticenterBond(id, child) => GraphIrConstraint::MulticenterBond(
                GraphIrMulticenterBondId(*id),
                child.bind(py).borrow().to_rust(py),
            ),
            Self::NoncovalentBond(id, child) => GraphIrConstraint::NoncovalentBond(
                GraphIrNoncovalentBondId(*id),
                child.bind(py).borrow().to_rust(py),
            ),
            Self::StereoAtom(id, kind, child) => GraphIrConstraint::StereoAtom(
                GraphIrStereoAtomId(*id),
                kind.to_rust(),
                child.bind(py).borrow().to_rust(py),
            ),
            Self::StereoBond(id, kind, child) => GraphIrConstraint::StereoBond(
                GraphIrStereoBondId(*id),
                kind.to_rust(),
                child.bind(py).borrow().to_rust(py),
            ),
            Self::Relational(child) => {
                GraphIrConstraint::Relational(child.bind(py).borrow().to_rust(py))
            }
            Self::Molecule(child) => {
                GraphIrConstraint::Molecule(child.bind(py).borrow().to_rust(py))
            }
            Self::And(children) => GraphIrConstraint::And(
                children
                    .iter()
                    .map(|child| child.bind(py).borrow().to_rust(py))
                    .collect(),
            ),
            Self::Or(children) => GraphIrConstraint::Or(
                children
                    .iter()
                    .map(|child| child.bind(py).borrow().to_rust(py))
                    .collect(),
            ),
            Self::Not(child) => {
                GraphIrConstraint::Not(Box::new(child.bind(py).borrow().to_rust(py)))
            }
        }
    }
}

/// Resolve a possibly-negative Python index into an existing constraint position.
fn resolve_constraint_index(len: usize, index: isize) -> PyResult<usize> {
    let resolved = if index < 0 {
        index + len as isize
    } else {
        index
    };
    if resolved < 0 || resolved as usize >= len {
        Err(PyIndexError::new_err("constraint index out of range"))
    } else {
        Ok(resolved as usize)
    }
}

/// Build a detached iterator of concrete Python constraint variants.
fn constraint_iter(py: Python<'_>, constraints: &GraphIrConstraints) -> PyResult<ConstraintIter> {
    let entries = constraints
        .iter()
        .map(|constraint| into_py_variant(py, Constraint::from_rust(py, constraint)?))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(ConstraintIter {
        entries: entries.into_iter(),
    })
}

/// A snapshot iterator over molecule-level constraints.
#[pyclass]
pub(crate) struct ConstraintIter {
    entries: IntoIter<Py<Constraint>>,
}

#[pymethods]
impl ConstraintIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<Py<Constraint>> {
        self.entries.next()
    }
}

/// The argument to `update`: another value container, a live view, or constraint entries.
#[derive(FromPyObject)]
pub(crate) enum ConstraintsUpdate {
    Container(Py<Constraints>),
    View(Py<ConstraintsView>),
    Entries(Vec<Py<Constraint>>),
}

impl ConstraintsUpdate {
    /// Snapshot every Python input before the target takes a write borrow.
    pub(crate) fn resolve(&self, py: Python<'_>) -> PyResult<ResolvedConstraintsUpdate> {
        Ok(match self {
            Self::Container(container) => {
                ResolvedConstraintsUpdate::Overlay(container.bind(py).borrow().inner().clone())
            }
            Self::View(view) => ResolvedConstraintsUpdate::Overlay(
                view.bind(py)
                    .borrow()
                    .read(py, |constraints| Ok(constraints.clone()))?,
            ),
            Self::Entries(entries) => ResolvedConstraintsUpdate::Entries(
                entries
                    .iter()
                    .map(|entry| entry.bind(py).borrow().to_rust(py))
                    .collect(),
            ),
        })
    }
}

/// A resolved update containing no Python references that need to be read.
pub(crate) enum ResolvedConstraintsUpdate {
    Overlay(GraphIrConstraints),
    Entries(Vec<GraphIrConstraint>),
}

impl ResolvedConstraintsUpdate {
    /// Append the resolved entries in order, preserving duplicates.
    pub(crate) fn apply(self, target: &mut GraphIrConstraints) {
        match self {
            Self::Overlay(overlay) => {
                for entry in overlay {
                    target.push(entry);
                }
            }
            Self::Entries(entries) => {
                for entry in entries {
                    target.push(entry);
                }
            }
        }
    }
}

/// A whole-container input that snapshots either a value container or a live view.
#[derive(FromPyObject)]
pub(crate) enum ConstraintsLike {
    Container(Py<Constraints>),
    View(Py<ConstraintsView>),
}

impl ConstraintsLike {
    pub(crate) fn to_rust(&self, py: Python<'_>) -> PyResult<GraphIrConstraints> {
        match self {
            Self::Container(container) => Ok(container.bind(py).borrow().inner().clone()),
            Self::View(view) => view
                .bind(py)
                .borrow()
                .read(py, |constraints| Ok(constraints.clone())),
        }
    }
}

/// The molecule-level constraints in insertion order. Mutable, value-equal,
/// and unhashable.
#[pyclass(eq)]
#[derive(Debug, PartialEq)]
pub struct Constraints(GraphIrConstraints);

#[pymethods]
impl Constraints {
    /// Build an owned container from constraint entries, preserving order and duplicates.
    #[new]
    fn new(py: Python<'_>, entries: Vec<Py<Constraint>>) -> Self {
        Self(GraphIrConstraints::from(
            entries
                .into_iter()
                .map(|entry| entry.bind(py).borrow().to_rust(py))
                .collect::<Vec<_>>(),
        ))
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let mut parts = Vec::with_capacity(self.0.len());
        for entry in self.0.iter() {
            let value = into_py_variant(py, Constraint::from_rust(py, entry)?)?;
            parts.push(value.bind(py).as_any().repr()?.extract::<String>()?);
        }
        Ok(format!("Constraints([{}])", parts.join(", ")))
    }

    /// Append one constraint, preserving existing entries and duplicates.
    fn append(&mut self, py: Python<'_>, constraint: Py<Constraint>) {
        self.0.push(constraint.bind(py).borrow().to_rust(py));
    }

    fn clear(&mut self) {
        self.0.clear();
    }

    /// Append another container, live view, or iterable after snapshotting the RHS.
    fn update(slf: Py<Self>, py: Python<'_>, other: ConstraintsUpdate) -> PyResult<()> {
        let resolved = other.resolve(py)?;
        resolved.apply(slf.borrow_mut(py).inner_mut());
        Ok(())
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    fn __getitem__(&self, py: Python<'_>, index: isize) -> PyResult<Constraint> {
        let index = resolve_constraint_index(self.0.len(), index)?;
        Constraint::from_rust(py, &self.0.as_slice()[index])
    }

    fn __iter__(&self, py: Python<'_>) -> PyResult<ConstraintIter> {
        constraint_iter(py, &self.0)
    }
}

impl Constraints {
    pub(crate) fn inner(&self) -> &GraphIrConstraints {
        &self.0
    }

    pub(crate) fn inner_mut(&mut self) -> &mut GraphIrConstraints {
        &mut self.0
    }

    pub(crate) fn from_inner(constraints: GraphIrConstraints) -> Self {
        Self(constraints)
    }
}

impl_py_canonicalize!(
    Constraints,
    GraphIrConstraints,
    |value: &Constraints, _py: Python<'_>| -> PyResult<GraphIrConstraints> {
        Ok(value.inner().clone())
    },
    |_py: Python<'_>, value: GraphIrConstraints| -> PyResult<Constraints> {
        Ok(Constraints::from_inner(value))
    }
);

/// A live handle onto the molecule-level constraints of one `MoleculeAst`.
#[pyclass]
pub struct ConstraintsView {
    pub(crate) owner: Py<MoleculeAst>,
}

impl ConstraintsView {
    pub(crate) fn new(owner: Py<MoleculeAst>) -> Self {
        Self { owner }
    }

    /// Borrow the current constraints and read through `f` without cloning the store.
    pub(crate) fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&GraphIrConstraints) -> PyResult<R>,
    ) -> PyResult<R> {
        let molecule = self.owner.bind(py).borrow();
        f(molecule.inner().constraints())
    }

    /// Mutate the molecule's constraint store in place through `f`.
    pub(crate) fn with_mut<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut GraphIrConstraints) -> R,
    ) -> R {
        f(self.owner.borrow_mut(py).inner_mut().constraints_mut())
    }
}

#[pymethods]
impl ConstraintsView {
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let count = self.read(py, |constraints| Ok(constraints.len()))?;
        Ok(format!("ConstraintsView({count} entries)"))
    }

    /// Append one constraint to the molecule, preserving existing entries and duplicates.
    fn append(&self, py: Python<'_>, constraint: Py<Constraint>) {
        let constraint = constraint.bind(py).borrow().to_rust(py);
        self.with_mut(py, |constraints| constraints.push(constraint));
    }

    fn clear(&self, py: Python<'_>) {
        self.with_mut(py, GraphIrConstraints::clear);
    }

    /// Append another container, live view, or iterable after snapshotting the RHS.
    fn update(&self, py: Python<'_>, other: ConstraintsUpdate) -> PyResult<()> {
        let resolved = other.resolve(py)?;
        self.with_mut(py, |constraints| resolved.apply(constraints));
        Ok(())
    }

    fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        self.read(py, |constraints| Ok(constraints.len()))
    }

    fn __getitem__(&self, py: Python<'_>, index: isize) -> PyResult<Constraint> {
        self.read(py, |constraints| {
            let index = resolve_constraint_index(constraints.len(), index)?;
            Constraint::from_rust(py, &constraints.as_slice()[index])
        })
    }

    fn __iter__(&self, py: Python<'_>) -> PyResult<ConstraintIter> {
        self.read(py, |constraints| constraint_iter(py, constraints))
    }
}

#[pymethods]
impl MoleculeConstraint {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        match &*slf.bind(py).borrow() {
            Self::ChargeSum(_, _) => {
                variant_repr(slf.bind(py).as_any(), "MoleculeConstraint", "ChargeSum", 2)
            }
            Self::UnpairedElectronCoupling { .. } => {
                let object = slf.bind(py).as_any();
                let atoms = object.getattr("atoms")?.repr()?.extract::<String>()?;
                let unpaired_electrons = object
                    .getattr("unpaired_electrons")?
                    .repr()?
                    .extract::<String>()?;
                Ok(format!(
                    "MoleculeConstraint.UnpairedElectronCoupling(atoms={atoms}, \
                     unpaired_electrons={unpaired_electrons})"
                ))
            }
            Self::BondOrderSum(_, _) => variant_repr(
                slf.bind(py).as_any(),
                "MoleculeConstraint",
                "BondOrderSum",
                2,
            ),
            Self::Connected(_) => {
                variant_repr(slf.bind(py).as_any(), "MoleculeConstraint", "Connected", 1)
            }
        }
    }
}

impl_py_canonicalize!(
    MoleculeConstraint,
    GraphIrMoleculeConstraint,
    |value: &MoleculeConstraint, py: Python<'_>| -> PyResult<GraphIrMoleculeConstraint> {
        Ok(value.to_rust(py))
    },
    |py: Python<'_>, value: GraphIrMoleculeConstraint| -> PyResult<MoleculeConstraint> {
        MoleculeConstraint::from_rust(py, &value)
    }
);

impl MoleculeConstraint {
    pub(crate) fn from_rust(
        py: Python<'_>,
        constraint: &GraphIrMoleculeConstraint,
    ) -> PyResult<Self> {
        Ok(match constraint {
            GraphIrMoleculeConstraint::ChargeSum { atoms, sum } => Self::ChargeSum(
                atoms
                    .as_ref()
                    .map(|atoms| atoms.iter().map(|atom| atom.0).collect()),
                into_py_variant(py, NumForm::from_rust(py, sum)?)?,
            ),
            GraphIrMoleculeConstraint::UnpairedElectronCoupling {
                atoms,
                unpaired_electrons,
            } => Self::UnpairedElectronCoupling {
                atoms: atoms
                    .as_ref()
                    .map(|atoms| atoms.iter().map(|atom| atom.0).collect()),
                unpaired_electrons: Py::new(
                    py,
                    UnpairedElectronsForm::from_rust(py, unpaired_electrons)?,
                )?,
            },
            GraphIrMoleculeConstraint::BondOrderSum { bonds, sum } => Self::BondOrderSum(
                bonds
                    .as_ref()
                    .map(|bonds| bonds.iter().map(|bond| bond.0).collect()),
                into_py_variant(py, NumForm::from_rust(py, sum)?)?,
            ),
            GraphIrMoleculeConstraint::Connected { atoms } => Self::Connected(
                atoms
                    .as_ref()
                    .map(|atoms| atoms.iter().map(|atom| atom.0).collect()),
            ),
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrMoleculeConstraint {
        match self {
            Self::ChargeSum(atoms, sum) => GraphIrMoleculeConstraint::ChargeSum {
                atoms: atoms
                    .as_ref()
                    .map(|atoms| atoms.iter().copied().map(GraphIrAtomId).collect()),
                sum: sum.bind(py).borrow().to_rust(py),
            },
            Self::UnpairedElectronCoupling {
                atoms,
                unpaired_electrons,
            } => GraphIrMoleculeConstraint::UnpairedElectronCoupling {
                atoms: atoms
                    .as_ref()
                    .map(|atoms| atoms.iter().copied().map(GraphIrAtomId).collect()),
                unpaired_electrons: unpaired_electrons.bind(py).borrow().to_rust(py),
            },
            Self::BondOrderSum(bonds, sum) => GraphIrMoleculeConstraint::BondOrderSum {
                bonds: bonds
                    .as_ref()
                    .map(|bonds| bonds.iter().copied().map(GraphIrBondId).collect()),
                sum: sum.bind(py).borrow().to_rust(py),
            },
            Self::Connected(atoms) => GraphIrMoleculeConstraint::Connected {
                atoms: atoms
                    .as_ref()
                    .map(|atoms| atoms.iter().copied().map(GraphIrAtomId).collect()),
            },
        }
    }
}

#[pymethods]
impl RelationalConstraint {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let variant = match &*slf.bind(py).borrow() {
            Self::DativeBondDonors(_, _) => "DativeBondDonors",
            Self::DativeBondDonor(_, _) => "DativeBondDonor",
            Self::DativeBondContainsAllDonors(_, _) => "DativeBondContainsAllDonors",
            Self::DativeBondAllDonors(_, _) => "DativeBondAllDonors",
            Self::DativeBondAnyDonor(_, _) => "DativeBondAnyDonor",
            Self::DativeBondAcceptor(_, _) => "DativeBondAcceptor",
            Self::DativeBondAcceptorSatisfies(_, _) => "DativeBondAcceptorSatisfies",
            Self::DativeBondParallels(_, _) => "DativeBondParallels",
            Self::AromaticSystemAtoms(_, _) => "AromaticSystemAtoms",
            Self::AromaticSystemContains(_, _) => "AromaticSystemContains",
            Self::AromaticSystemContainsAll(_, _) => "AromaticSystemContainsAll",
            Self::AromaticSystemAllAtoms(_, _) => "AromaticSystemAllAtoms",
            Self::AromaticSystemAnyAtom(_, _) => "AromaticSystemAnyAtom",
            Self::MulticenterBondAtoms(_, _) => "MulticenterBondAtoms",
            Self::MulticenterBondContains(_, _) => "MulticenterBondContains",
            Self::MulticenterBondContainsAll(_, _) => "MulticenterBondContainsAll",
            Self::MulticenterBondAllAtoms(_, _) => "MulticenterBondAllAtoms",
            Self::MulticenterBondAnyAtom(_, _) => "MulticenterBondAnyAtom",
            Self::NoncovalentBondEnds(_, _) => "NoncovalentBondEnds",
            Self::NoncovalentBondContains(_, _) => "NoncovalentBondContains",
            Self::NoncovalentBondEndsSatisfy(_, _) => "NoncovalentBondEndsSatisfy",
            Self::StereoAtomSite(_, _) => "StereoAtomSite",
            Self::StereoAtomContains(_, _) => "StereoAtomContains",
            Self::StereoAtomLigands(_, _) => "StereoAtomLigands",
            Self::StereoAtomAllLigands(_, _) => "StereoAtomAllLigands",
            Self::StereoAtomAnyLigand(_, _) => "StereoAtomAnyLigand",
            Self::StereoBondSite(_, _) => "StereoBondSite",
            Self::StereoBondContains(_, _) => "StereoBondContains",
            Self::StereoBondLigands(_, _) => "StereoBondLigands",
            Self::StereoBondAllLigands(_, _) => "StereoBondAllLigands",
            Self::StereoBondAnyLigand(_, _) => "StereoBondAnyLigand",
        };
        variant_repr(slf.bind(py).as_any(), "RelationalConstraint", variant, 2)
    }
}

impl RelationalConstraint {
    /// Convert any relational constraint into its Python value.
    pub(crate) fn from_rust(
        py: Python<'_>,
        constraint: &GraphIrRelationalConstraint,
    ) -> PyResult<Self> {
        Ok(match constraint {
            GraphIrRelationalConstraint::DativeBondDonors { bond, atoms } => {
                Self::DativeBondDonors(bond.0, atoms.iter().map(|atom| atom.0).collect())
            }
            GraphIrRelationalConstraint::DativeBondDonor { bond, atom } => {
                Self::DativeBondDonor(bond.0, atom.0)
            }
            GraphIrRelationalConstraint::DativeBondContainsAllDonors { bond, atoms } => {
                Self::DativeBondContainsAllDonors(bond.0, atoms.iter().map(|atom| atom.0).collect())
            }
            GraphIrRelationalConstraint::DativeBondAllDonors { bond, predicate } => {
                Self::DativeBondAllDonors(
                    bond.0,
                    into_py_variant(py, AtomConstraintForm::from_rust(py, predicate)?)?,
                )
            }
            GraphIrRelationalConstraint::DativeBondAnyDonor { bond, predicate } => {
                Self::DativeBondAnyDonor(
                    bond.0,
                    into_py_variant(py, AtomConstraintForm::from_rust(py, predicate)?)?,
                )
            }
            GraphIrRelationalConstraint::DativeBondAcceptor { bond, atom } => {
                Self::DativeBondAcceptor(bond.0, atom.0)
            }
            GraphIrRelationalConstraint::DativeBondAcceptorSatisfies { bond, predicate } => {
                Self::DativeBondAcceptorSatisfies(
                    bond.0,
                    into_py_variant(py, AtomConstraintForm::from_rust(py, predicate)?)?,
                )
            }
            GraphIrRelationalConstraint::DativeBondParallels { dative, parallel } => {
                Self::DativeBondParallels(dative.0, parallel.0)
            }
            GraphIrRelationalConstraint::AromaticSystemAtoms { system, atoms } => {
                Self::AromaticSystemAtoms(system.0, atoms.iter().map(|atom| atom.0).collect())
            }
            GraphIrRelationalConstraint::AromaticSystemContains { system, atom } => {
                Self::AromaticSystemContains(system.0, atom.0)
            }
            GraphIrRelationalConstraint::AromaticSystemContainsAll { system, atoms } => {
                Self::AromaticSystemContainsAll(system.0, atoms.iter().map(|atom| atom.0).collect())
            }
            GraphIrRelationalConstraint::AromaticSystemAllAtoms { system, predicate } => {
                Self::AromaticSystemAllAtoms(
                    system.0,
                    into_py_variant(py, AtomConstraintForm::from_rust(py, predicate)?)?,
                )
            }
            GraphIrRelationalConstraint::AromaticSystemAnyAtom { system, predicate } => {
                Self::AromaticSystemAnyAtom(
                    system.0,
                    into_py_variant(py, AtomConstraintForm::from_rust(py, predicate)?)?,
                )
            }
            GraphIrRelationalConstraint::MulticenterBondAtoms { bond, atoms } => {
                Self::MulticenterBondAtoms(bond.0, atoms.iter().map(|atom| atom.0).collect())
            }
            GraphIrRelationalConstraint::MulticenterBondContains { bond, atom } => {
                Self::MulticenterBondContains(bond.0, atom.0)
            }
            GraphIrRelationalConstraint::MulticenterBondContainsAll { bond, atoms } => {
                Self::MulticenterBondContainsAll(bond.0, atoms.iter().map(|atom| atom.0).collect())
            }
            GraphIrRelationalConstraint::MulticenterBondAllAtoms { bond, predicate } => {
                Self::MulticenterBondAllAtoms(
                    bond.0,
                    into_py_variant(py, AtomConstraintForm::from_rust(py, predicate)?)?,
                )
            }
            GraphIrRelationalConstraint::MulticenterBondAnyAtom { bond, predicate } => {
                Self::MulticenterBondAnyAtom(
                    bond.0,
                    into_py_variant(py, AtomConstraintForm::from_rust(py, predicate)?)?,
                )
            }
            GraphIrRelationalConstraint::NoncovalentBondEnds { bond, atoms } => {
                Self::NoncovalentBondEnds(bond.0, [atoms[0].0, atoms[1].0])
            }
            GraphIrRelationalConstraint::NoncovalentBondContains { bond, atom } => {
                Self::NoncovalentBondContains(bond.0, atom.0)
            }
            GraphIrRelationalConstraint::NoncovalentBondEndsSatisfy { bond, predicates } => {
                Self::NoncovalentBondEndsSatisfy(
                    bond.0,
                    [
                        into_py_variant(py, AtomConstraintForm::from_rust(py, &predicates[0])?)?,
                        into_py_variant(py, AtomConstraintForm::from_rust(py, &predicates[1])?)?,
                    ],
                )
            }
            GraphIrRelationalConstraint::StereoAtomSite { stereo_atom, atom } => {
                Self::StereoAtomSite(stereo_atom.0, atom.0)
            }
            GraphIrRelationalConstraint::StereoAtomContains { stereo_atom, atom } => {
                Self::StereoAtomContains(stereo_atom.0, atom.0)
            }
            GraphIrRelationalConstraint::StereoAtomLigands { stereo_atom, atoms } => {
                Self::StereoAtomLigands(stereo_atom.0, atoms.iter().map(|atom| atom.0).collect())
            }
            GraphIrRelationalConstraint::StereoAtomAllLigands {
                stereo_atom,
                predicate,
            } => Self::StereoAtomAllLigands(
                stereo_atom.0,
                into_py_variant(py, AtomConstraintForm::from_rust(py, predicate)?)?,
            ),
            GraphIrRelationalConstraint::StereoAtomAnyLigand {
                stereo_atom,
                predicate,
            } => Self::StereoAtomAnyLigand(
                stereo_atom.0,
                into_py_variant(py, AtomConstraintForm::from_rust(py, predicate)?)?,
            ),
            GraphIrRelationalConstraint::StereoBondSite { stereo_bond, bond } => {
                Self::StereoBondSite(stereo_bond.0, bond.0)
            }
            GraphIrRelationalConstraint::StereoBondContains { stereo_bond, atom } => {
                Self::StereoBondContains(stereo_bond.0, atom.0)
            }
            GraphIrRelationalConstraint::StereoBondLigands { stereo_bond, atoms } => {
                Self::StereoBondLigands(stereo_bond.0, atoms.iter().map(|atom| atom.0).collect())
            }
            GraphIrRelationalConstraint::StereoBondAllLigands {
                stereo_bond,
                predicate,
            } => Self::StereoBondAllLigands(
                stereo_bond.0,
                into_py_variant(py, AtomConstraintForm::from_rust(py, predicate)?)?,
            ),
            GraphIrRelationalConstraint::StereoBondAnyLigand {
                stereo_bond,
                predicate,
            } => Self::StereoBondAnyLigand(
                stereo_bond.0,
                into_py_variant(py, AtomConstraintForm::from_rust(py, predicate)?)?,
            ),
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrRelationalConstraint {
        match self {
            Self::DativeBondDonors(bond, atoms) => GraphIrRelationalConstraint::DativeBondDonors {
                bond: GraphIrDativeBondId(*bond),
                atoms: atoms.iter().copied().map(GraphIrAtomId).collect(),
            },
            Self::DativeBondDonor(bond, atom) => GraphIrRelationalConstraint::DativeBondDonor {
                bond: GraphIrDativeBondId(*bond),
                atom: GraphIrAtomId(*atom),
            },
            Self::DativeBondContainsAllDonors(bond, atoms) => {
                GraphIrRelationalConstraint::DativeBondContainsAllDonors {
                    bond: GraphIrDativeBondId(*bond),
                    atoms: atoms.iter().copied().map(GraphIrAtomId).collect(),
                }
            }
            Self::DativeBondAllDonors(bond, predicate) => {
                GraphIrRelationalConstraint::DativeBondAllDonors {
                    bond: GraphIrDativeBondId(*bond),
                    predicate: Box::new(predicate.bind(py).borrow().to_rust(py)),
                }
            }
            Self::DativeBondAnyDonor(bond, predicate) => {
                GraphIrRelationalConstraint::DativeBondAnyDonor {
                    bond: GraphIrDativeBondId(*bond),
                    predicate: Box::new(predicate.bind(py).borrow().to_rust(py)),
                }
            }
            Self::DativeBondAcceptor(bond, atom) => {
                GraphIrRelationalConstraint::DativeBondAcceptor {
                    bond: GraphIrDativeBondId(*bond),
                    atom: GraphIrAtomId(*atom),
                }
            }
            Self::DativeBondAcceptorSatisfies(bond, predicate) => {
                GraphIrRelationalConstraint::DativeBondAcceptorSatisfies {
                    bond: GraphIrDativeBondId(*bond),
                    predicate: Box::new(predicate.bind(py).borrow().to_rust(py)),
                }
            }
            Self::DativeBondParallels(dative, parallel) => {
                GraphIrRelationalConstraint::DativeBondParallels {
                    dative: GraphIrDativeBondId(*dative),
                    parallel: GraphIrBondId(*parallel),
                }
            }
            Self::AromaticSystemAtoms(system, atoms) => {
                GraphIrRelationalConstraint::AromaticSystemAtoms {
                    system: GraphIrAromaticSystemId(*system),
                    atoms: atoms.iter().copied().map(GraphIrAtomId).collect(),
                }
            }
            Self::AromaticSystemContains(system, atom) => {
                GraphIrRelationalConstraint::AromaticSystemContains {
                    system: GraphIrAromaticSystemId(*system),
                    atom: GraphIrAtomId(*atom),
                }
            }
            Self::AromaticSystemContainsAll(system, atoms) => {
                GraphIrRelationalConstraint::AromaticSystemContainsAll {
                    system: GraphIrAromaticSystemId(*system),
                    atoms: atoms.iter().copied().map(GraphIrAtomId).collect(),
                }
            }
            Self::AromaticSystemAllAtoms(system, predicate) => {
                GraphIrRelationalConstraint::AromaticSystemAllAtoms {
                    system: GraphIrAromaticSystemId(*system),
                    predicate: Box::new(predicate.bind(py).borrow().to_rust(py)),
                }
            }
            Self::AromaticSystemAnyAtom(system, predicate) => {
                GraphIrRelationalConstraint::AromaticSystemAnyAtom {
                    system: GraphIrAromaticSystemId(*system),
                    predicate: Box::new(predicate.bind(py).borrow().to_rust(py)),
                }
            }
            Self::MulticenterBondAtoms(bond, atoms) => {
                GraphIrRelationalConstraint::MulticenterBondAtoms {
                    bond: GraphIrMulticenterBondId(*bond),
                    atoms: atoms.iter().copied().map(GraphIrAtomId).collect(),
                }
            }
            Self::MulticenterBondContains(bond, atom) => {
                GraphIrRelationalConstraint::MulticenterBondContains {
                    bond: GraphIrMulticenterBondId(*bond),
                    atom: GraphIrAtomId(*atom),
                }
            }
            Self::MulticenterBondContainsAll(bond, atoms) => {
                GraphIrRelationalConstraint::MulticenterBondContainsAll {
                    bond: GraphIrMulticenterBondId(*bond),
                    atoms: atoms.iter().copied().map(GraphIrAtomId).collect(),
                }
            }
            Self::MulticenterBondAllAtoms(bond, predicate) => {
                GraphIrRelationalConstraint::MulticenterBondAllAtoms {
                    bond: GraphIrMulticenterBondId(*bond),
                    predicate: Box::new(predicate.bind(py).borrow().to_rust(py)),
                }
            }
            Self::MulticenterBondAnyAtom(bond, predicate) => {
                GraphIrRelationalConstraint::MulticenterBondAnyAtom {
                    bond: GraphIrMulticenterBondId(*bond),
                    predicate: Box::new(predicate.bind(py).borrow().to_rust(py)),
                }
            }
            Self::NoncovalentBondEnds(bond, atoms) => {
                GraphIrRelationalConstraint::NoncovalentBondEnds {
                    bond: GraphIrNoncovalentBondId(*bond),
                    atoms: [GraphIrAtomId(atoms[0]), GraphIrAtomId(atoms[1])],
                }
            }
            Self::NoncovalentBondContains(bond, atom) => {
                GraphIrRelationalConstraint::NoncovalentBondContains {
                    bond: GraphIrNoncovalentBondId(*bond),
                    atom: GraphIrAtomId(*atom),
                }
            }
            Self::NoncovalentBondEndsSatisfy(bond, predicates) => {
                GraphIrRelationalConstraint::NoncovalentBondEndsSatisfy {
                    bond: GraphIrNoncovalentBondId(*bond),
                    predicates: [
                        Box::new(predicates[0].bind(py).borrow().to_rust(py)),
                        Box::new(predicates[1].bind(py).borrow().to_rust(py)),
                    ],
                }
            }
            Self::StereoAtomSite(stereo_atom, atom) => {
                GraphIrRelationalConstraint::StereoAtomSite {
                    stereo_atom: GraphIrStereoAtomId(*stereo_atom),
                    atom: GraphIrAtomId(*atom),
                }
            }
            Self::StereoAtomContains(stereo_atom, atom) => {
                GraphIrRelationalConstraint::StereoAtomContains {
                    stereo_atom: GraphIrStereoAtomId(*stereo_atom),
                    atom: GraphIrAtomId(*atom),
                }
            }
            Self::StereoAtomLigands(stereo_atom, atoms) => {
                GraphIrRelationalConstraint::StereoAtomLigands {
                    stereo_atom: GraphIrStereoAtomId(*stereo_atom),
                    atoms: atoms.iter().copied().map(GraphIrAtomId).collect(),
                }
            }
            Self::StereoAtomAllLigands(stereo_atom, predicate) => {
                GraphIrRelationalConstraint::StereoAtomAllLigands {
                    stereo_atom: GraphIrStereoAtomId(*stereo_atom),
                    predicate: Box::new(predicate.bind(py).borrow().to_rust(py)),
                }
            }
            Self::StereoAtomAnyLigand(stereo_atom, predicate) => {
                GraphIrRelationalConstraint::StereoAtomAnyLigand {
                    stereo_atom: GraphIrStereoAtomId(*stereo_atom),
                    predicate: Box::new(predicate.bind(py).borrow().to_rust(py)),
                }
            }
            Self::StereoBondSite(stereo_bond, bond) => {
                GraphIrRelationalConstraint::StereoBondSite {
                    stereo_bond: GraphIrStereoBondId(*stereo_bond),
                    bond: GraphIrBondId(*bond),
                }
            }
            Self::StereoBondContains(stereo_bond, atom) => {
                GraphIrRelationalConstraint::StereoBondContains {
                    stereo_bond: GraphIrStereoBondId(*stereo_bond),
                    atom: GraphIrAtomId(*atom),
                }
            }
            Self::StereoBondLigands(stereo_bond, atoms) => {
                GraphIrRelationalConstraint::StereoBondLigands {
                    stereo_bond: GraphIrStereoBondId(*stereo_bond),
                    atoms: atoms.iter().copied().map(GraphIrAtomId).collect(),
                }
            }
            Self::StereoBondAllLigands(stereo_bond, predicate) => {
                GraphIrRelationalConstraint::StereoBondAllLigands {
                    stereo_bond: GraphIrStereoBondId(*stereo_bond),
                    predicate: Box::new(predicate.bind(py).borrow().to_rust(py)),
                }
            }
            Self::StereoBondAnyLigand(stereo_bond, predicate) => {
                GraphIrRelationalConstraint::StereoBondAnyLigand {
                    stereo_bond: GraphIrStereoBondId(*stereo_bond),
                    predicate: Box::new(predicate.bind(py).borrow().to_rust(py)),
                }
            }
        }
    }
}

impl_py_canonicalize!(
    RelationalConstraint,
    GraphIrRelationalConstraint,
    |value: &RelationalConstraint, py: Python<'_>| -> PyResult<GraphIrRelationalConstraint> {
        Ok(value.to_rust(py))
    },
    |py: Python<'_>, value: GraphIrRelationalConstraint| -> PyResult<RelationalConstraint> {
        RelationalConstraint::from_rust(py, &value)
    }
);

#[cfg(test)]
mod tests {
    use std::ffi::CString;

    use pyo3::types::PyDict;
    use rstest::rstest;
    use umol_graph_ir::ir::{
        AromaticSystemConstraintForm as GraphIrAromaticSystemConstraintForm,
        AtomConstraintForm as GraphIrAtomConstraintForm,
        BondConstraintForm as GraphIrBondConstraintForm,
        DativeBondConstraintForm as GraphIrDativeBondConstraintForm, Molecule as GraphIrMolecule,
        MulticenterBondConstraintForm as GraphIrMulticenterBondConstraintForm,
        NoncovalentBondConstraintForm as GraphIrNoncovalentBondConstraintForm,
        NumForm as GraphIrNumForm, StereoAtomConstraintForm as GraphIrStereoAtomConstraintForm,
        StereoBondConstraintForm as GraphIrStereoBondConstraintForm,
        StereoKind as GraphIrStereoKind, Stereogenicity as GraphIrStereogenicity,
        StereogenicityForm as GraphIrStereogenicityForm,
        UnpairedElectronsForm as GraphIrUnpairedElectronsForm,
    };

    use super::*;

    #[rstest]
    #[case::donors(GraphIrRelationalConstraint::DativeBondDonors {
        bond: GraphIrDativeBondId(1),
        atoms: vec![GraphIrAtomId(2), GraphIrAtomId(3)],
    })]
    #[case::donor(GraphIrRelationalConstraint::DativeBondDonor {
        bond: GraphIrDativeBondId(4),
        atom: GraphIrAtomId(5),
    })]
    #[case::contains_all_donors(GraphIrRelationalConstraint::DativeBondContainsAllDonors {
        bond: GraphIrDativeBondId(6),
        atoms: vec![GraphIrAtomId(7), GraphIrAtomId(8)],
    })]
    #[case::all_donors(GraphIrRelationalConstraint::DativeBondAllDonors {
        bond: GraphIrDativeBondId(9),
        predicate: Box::new(GraphIrAtomConstraintForm::degree(2)),
    })]
    #[case::any_donor(GraphIrRelationalConstraint::DativeBondAnyDonor {
        bond: GraphIrDativeBondId(10),
        predicate: Box::new(GraphIrAtomConstraintForm::valence(3)),
    })]
    #[case::acceptor(GraphIrRelationalConstraint::DativeBondAcceptor {
        bond: GraphIrDativeBondId(11),
        atom: GraphIrAtomId(12),
    })]
    #[case::acceptor_satisfies(GraphIrRelationalConstraint::DativeBondAcceptorSatisfies {
        bond: GraphIrDativeBondId(13),
        predicate: Box::new(GraphIrAtomConstraintForm::total_degree(4)),
    })]
    #[case::parallels(GraphIrRelationalConstraint::DativeBondParallels {
        dative: GraphIrDativeBondId(14),
        parallel: GraphIrBondId(15),
    })]
    #[case::aromatic_atoms(GraphIrRelationalConstraint::AromaticSystemAtoms {
        system: GraphIrAromaticSystemId(16),
        atoms: vec![GraphIrAtomId(17), GraphIrAtomId(18)],
    })]
    #[case::aromatic_contains(GraphIrRelationalConstraint::AromaticSystemContains {
        system: GraphIrAromaticSystemId(19),
        atom: GraphIrAtomId(20),
    })]
    #[case::aromatic_contains_all(GraphIrRelationalConstraint::AromaticSystemContainsAll {
        system: GraphIrAromaticSystemId(21),
        atoms: vec![GraphIrAtomId(22), GraphIrAtomId(23)],
    })]
    #[case::aromatic_all_atoms(GraphIrRelationalConstraint::AromaticSystemAllAtoms {
        system: GraphIrAromaticSystemId(24),
        predicate: Box::new(GraphIrAtomConstraintForm::degree(5)),
    })]
    #[case::aromatic_any_atom(GraphIrRelationalConstraint::AromaticSystemAnyAtom {
        system: GraphIrAromaticSystemId(25),
        predicate: Box::new(GraphIrAtomConstraintForm::valence(6)),
    })]
    #[case::multicenter_atoms(GraphIrRelationalConstraint::MulticenterBondAtoms {
        bond: GraphIrMulticenterBondId(26),
        atoms: vec![GraphIrAtomId(27), GraphIrAtomId(28)],
    })]
    #[case::multicenter_contains(GraphIrRelationalConstraint::MulticenterBondContains {
        bond: GraphIrMulticenterBondId(29),
        atom: GraphIrAtomId(30),
    })]
    #[case::multicenter_contains_all(GraphIrRelationalConstraint::MulticenterBondContainsAll {
        bond: GraphIrMulticenterBondId(31),
        atoms: vec![GraphIrAtomId(32), GraphIrAtomId(33)],
    })]
    #[case::multicenter_all_atoms(GraphIrRelationalConstraint::MulticenterBondAllAtoms {
        bond: GraphIrMulticenterBondId(34),
        predicate: Box::new(GraphIrAtomConstraintForm::degree(7)),
    })]
    #[case::multicenter_any_atom(GraphIrRelationalConstraint::MulticenterBondAnyAtom {
        bond: GraphIrMulticenterBondId(35),
        predicate: Box::new(GraphIrAtomConstraintForm::valence(8)),
    })]
    #[case::noncovalent_ends(GraphIrRelationalConstraint::NoncovalentBondEnds {
        bond: GraphIrNoncovalentBondId(36),
        atoms: [GraphIrAtomId(37), GraphIrAtomId(38)],
    })]
    #[case::noncovalent_contains(GraphIrRelationalConstraint::NoncovalentBondContains {
        bond: GraphIrNoncovalentBondId(39),
        atom: GraphIrAtomId(40),
    })]
    #[case::noncovalent_ends_satisfy(GraphIrRelationalConstraint::NoncovalentBondEndsSatisfy {
        bond: GraphIrNoncovalentBondId(41),
        predicates: [
            Box::new(GraphIrAtomConstraintForm::degree(9)),
            Box::new(GraphIrAtomConstraintForm::valence(10)),
        ],
    })]
    #[case::stereo_atom_site(GraphIrRelationalConstraint::StereoAtomSite {
        stereo_atom: GraphIrStereoAtomId(42),
        atom: GraphIrAtomId(43),
    })]
    #[case::stereo_atom_contains(GraphIrRelationalConstraint::StereoAtomContains {
        stereo_atom: GraphIrStereoAtomId(44),
        atom: GraphIrAtomId(45),
    })]
    #[case::stereo_atom_ligands(GraphIrRelationalConstraint::StereoAtomLigands {
        stereo_atom: GraphIrStereoAtomId(46),
        atoms: vec![GraphIrAtomId(47), GraphIrAtomId(48)],
    })]
    #[case::stereo_atom_all_ligands(GraphIrRelationalConstraint::StereoAtomAllLigands {
        stereo_atom: GraphIrStereoAtomId(49),
        predicate: Box::new(GraphIrAtomConstraintForm::degree(11)),
    })]
    #[case::stereo_atom_any_ligand(GraphIrRelationalConstraint::StereoAtomAnyLigand {
        stereo_atom: GraphIrStereoAtomId(50),
        predicate: Box::new(GraphIrAtomConstraintForm::valence(12)),
    })]
    #[case::stereo_bond_site(GraphIrRelationalConstraint::StereoBondSite {
        stereo_bond: GraphIrStereoBondId(51),
        bond: GraphIrBondId(52),
    })]
    #[case::stereo_bond_contains(GraphIrRelationalConstraint::StereoBondContains {
        stereo_bond: GraphIrStereoBondId(53),
        atom: GraphIrAtomId(54),
    })]
    #[case::stereo_bond_ligands(GraphIrRelationalConstraint::StereoBondLigands {
        stereo_bond: GraphIrStereoBondId(55),
        atoms: vec![GraphIrAtomId(56), GraphIrAtomId(57)],
    })]
    #[case::stereo_bond_all_ligands(GraphIrRelationalConstraint::StereoBondAllLigands {
        stereo_bond: GraphIrStereoBondId(58),
        predicate: Box::new(GraphIrAtomConstraintForm::degree(13)),
    })]
    #[case::stereo_bond_any_ligand(GraphIrRelationalConstraint::StereoBondAnyLigand {
        stereo_bond: GraphIrStereoBondId(59),
        predicate: Box::new(GraphIrAtomConstraintForm::valence(14)),
    })]
    fn test_relational_constraint_roundtrip(#[case] constraint: GraphIrRelationalConstraint) {
        Python::attach(|py| {
            let value = RelationalConstraint::from_rust(py, &constraint).unwrap();
            assert_eq!(value.to_rust(py), constraint);
        });
    }

    #[rstest]
    #[case::charge_sum_whole(GraphIrMoleculeConstraint::ChargeSum {
        atoms: None,
        sum: GraphIrNumForm::Lit(1),
    })]
    #[case::charge_sum_empty_subset(GraphIrMoleculeConstraint::ChargeSum {
        atoms: Some(Vec::new()),
        sum: GraphIrNumForm::Lit(2),
    })]
    #[case::unpaired_electron_coupling(GraphIrMoleculeConstraint::UnpairedElectronCoupling {
        atoms: Some(vec![GraphIrAtomId(3), GraphIrAtomId(4)]),
        unpaired_electrons: GraphIrUnpairedElectronsForm::from((1, 2)),
    })]
    #[case::bond_order_sum(GraphIrMoleculeConstraint::BondOrderSum {
        bonds: Some(vec![GraphIrBondId(5), GraphIrBondId(6)]),
        sum: GraphIrNumForm::Lit(3),
    })]
    #[case::connected(GraphIrMoleculeConstraint::Connected {
        atoms: None,
    })]
    fn test_molecule_constraint_roundtrip(#[case] constraint: GraphIrMoleculeConstraint) {
        Python::attach(|py| {
            let value = MoleculeConstraint::from_rust(py, &constraint).unwrap();
            assert_eq!(value.to_rust(py), constraint);
        });
    }

    #[rstest]
    #[case::atom(GraphIrConstraint::Atom(GraphIrAtomId(1), GraphIrAtomConstraintForm::degree(2),))]
    #[case::bond(GraphIrConstraint::Bond(
        GraphIrBondId(3),
        GraphIrBondConstraintForm::aromatic(true),
    ))]
    #[case::dative_bond(GraphIrConstraint::DativeBond(
        GraphIrDativeBondId(4),
        GraphIrDativeBondConstraintForm::aromatic(false),
    ))]
    #[case::aromatic_system(GraphIrConstraint::AromaticSystem(
        GraphIrAromaticSystemId(5),
        GraphIrAromaticSystemConstraintForm::electron_count(6),
    ))]
    #[case::multicenter_bond(GraphIrConstraint::MulticenterBond(
        GraphIrMulticenterBondId(7),
        GraphIrMulticenterBondConstraintForm::electron_count(8),
    ))]
    #[case::noncovalent_bond(GraphIrConstraint::NoncovalentBond(
        GraphIrNoncovalentBondId(9),
        GraphIrNoncovalentBondConstraintForm::intramolecular(true),
    ))]
    #[case::stereo_atom(GraphIrConstraint::StereoAtom(
        GraphIrStereoAtomId(10),
        GraphIrStereoKind::Tetrahedral,
        GraphIrStereoAtomConstraintForm::Stereogenicity(GraphIrStereogenicityForm::Lit(
            GraphIrStereogenicity::Stereogenic,
        )),
    ))]
    #[case::stereo_bond(GraphIrConstraint::StereoBond(
        GraphIrStereoBondId(11),
        GraphIrStereoKind::CisTrans,
        GraphIrStereoBondConstraintForm::Stereogenicity(GraphIrStereogenicityForm::Lit(
            GraphIrStereogenicity::Prochiral,
        )),
    ))]
    #[case::relational(GraphIrConstraint::Relational(
        GraphIrRelationalConstraint::DativeBondDonor {
            bond: GraphIrDativeBondId(12),
            atom: GraphIrAtomId(13),
        },
    ))]
    #[case::molecule(GraphIrConstraint::Molecule(
        GraphIrMoleculeConstraint::Connected {
            atoms: Some(vec![GraphIrAtomId(14), GraphIrAtomId(15)]),
        },
    ))]
    #[case::and(GraphIrConstraint::And(Vec::new()))]
    #[case::or(GraphIrConstraint::Or(Vec::new()))]
    #[case::not(GraphIrConstraint::Not(Box::new(GraphIrConstraint::Atom(
        GraphIrAtomId(16),
        GraphIrAtomConstraintForm::degree(3),
    ))))]
    fn test_constraint_roundtrip(#[case] constraint: GraphIrConstraint) {
        Python::attach(|py| {
            let value = Constraint::from_rust(py, &constraint).unwrap();
            assert_eq!(value.to_rust(py), constraint);
        });
    }

    #[rstest]
    fn test_constraint_roundtrip_recursive() {
        let constraint = GraphIrConstraint::And(vec![
            GraphIrConstraint::Atom(GraphIrAtomId(17), GraphIrAtomConstraintForm::valence(4)),
            GraphIrConstraint::Or(vec![
                GraphIrConstraint::Relational(GraphIrRelationalConstraint::DativeBondDonor {
                    bond: GraphIrDativeBondId(18),
                    atom: GraphIrAtomId(19),
                }),
                GraphIrConstraint::Not(Box::new(GraphIrConstraint::Molecule(
                    GraphIrMoleculeConstraint::Connected {
                        atoms: Some(vec![GraphIrAtomId(20), GraphIrAtomId(21)]),
                    },
                ))),
            ]),
        ]);

        Python::attach(|py| {
            let value =
                into_py_variant(py, Constraint::from_rust(py, &constraint).unwrap()).unwrap();
            let equal =
                into_py_variant(py, Constraint::from_rust(py, &constraint).unwrap()).unwrap();

            assert_eq!(value.bind(py).borrow().to_rust(py), constraint);
            assert!(value.bind(py).as_any().eq(equal.bind(py).as_any()).unwrap());
            assert_eq!(
                value
                    .bind(py)
                    .as_any()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                "Constraint.And([Constraint.Atom(17, AtomConstraintForm.Valence(NumForm.Lit(4))), Constraint.Or([Constraint.Relational(RelationalConstraint.DativeBondDonor(18, 19)), Constraint.Not(Constraint.Molecule(MoleculeConstraint.Connected([20, 21])))])])"
            );

            let children = value.bind(py).as_any().getattr("_0").unwrap();
            assert_eq!(children.len().unwrap(), 2);
            assert_eq!(
                children
                    .get_item(0)
                    .unwrap()
                    .getattr("_0")
                    .unwrap()
                    .extract::<u32>()
                    .unwrap(),
                17
            );

            let locals = PyDict::new(py);
            locals.set_item("node", &value).unwrap();
            let source = CString::new(
                r#"
And = type(node)
Atom = type(node._0[0])
Or = type(node._0[1])
Relational = type(node._0[1]._0[0])
Not = type(node._0[1]._0[1])
Molecule = type(node._0[1]._0[1]._0)
match node:
    case And([Atom(atom_id, _), Or([Relational(_), Not(Molecule(_))])]):
        matched_atom_id = atom_id
    case _:
        matched_atom_id = None
"#,
            )
            .unwrap();
            py.run(source.as_c_str(), None, Some(&locals)).unwrap();
            assert_eq!(
                locals
                    .get_item("matched_atom_id")
                    .unwrap()
                    .unwrap()
                    .extract::<u32>()
                    .unwrap(),
                17
            );
        });
    }

    #[rstest]
    #[case::positive_first(3, 0, 0)]
    #[case::positive_last(3, 2, 2)]
    #[case::negative_last(3, -1, 2)]
    #[case::negative_first(3, -3, 0)]
    fn test_resolve_constraint_index(
        #[case] len: usize,
        #[case] index: isize,
        #[case] expected: usize,
    ) {
        assert_eq!(resolve_constraint_index(len, index).unwrap(), expected);
    }

    #[rstest]
    #[case::empty(0, 0)]
    #[case::positive(3, 3)]
    #[case::negative(3, -4)]
    fn test_resolve_constraint_index_error(#[case] len: usize, #[case] index: isize) {
        assert_eq!(
            resolve_constraint_index(len, index)
                .unwrap_err()
                .to_string(),
            "IndexError: constraint index out of range"
        );
    }

    #[rstest]
    fn test_constraint_iter() {
        let first = GraphIrConstraint::Atom(GraphIrAtomId(1), GraphIrAtomConstraintForm::degree(2));
        let second =
            GraphIrConstraint::Molecule(GraphIrMoleculeConstraint::Connected { atoms: None });
        let mut constraints = GraphIrConstraints::from(vec![first.clone(), second.clone()]);

        Python::attach(|py| {
            let mut iter = constraint_iter(py, &constraints).unwrap();
            constraints.push(GraphIrConstraint::Or(Vec::new()));

            let first_mirror = iter.__next__().unwrap();
            assert_eq!(
                first_mirror
                    .bind(py)
                    .as_any()
                    .getattr("_0")
                    .unwrap()
                    .extract::<u32>()
                    .unwrap(),
                1
            );
            assert_eq!(first_mirror.bind(py).borrow().to_rust(py), first);
            assert_eq!(
                iter.__next__().unwrap().bind(py).borrow().to_rust(py),
                second
            );
            assert!(iter.__next__().is_none());
        });
    }

    #[rstest]
    fn test_constraints_like_to_rust_container() {
        let expected = GraphIrConstraints::from(vec![
            GraphIrConstraint::And(Vec::new()),
            GraphIrConstraint::And(Vec::new()),
        ]);

        Python::attach(|py| {
            let container = Py::new(py, Constraints::from_inner(expected.clone())).unwrap();
            let arg = ConstraintsLike::Container(container);

            assert_eq!(arg.to_rust(py).unwrap(), expected);
        });
    }

    #[rstest]
    fn test_constraints_like_to_rust_view() {
        let expected = GraphIrConstraints::from(vec![
            GraphIrConstraint::Or(Vec::new()),
            GraphIrConstraint::Or(Vec::new()),
        ]);
        let mut molecule = GraphIrMolecule::new();
        *molecule.constraints_mut() = expected.clone();

        Python::attach(|py| {
            let owner = Py::new(py, MoleculeAst::from_rust(molecule)).unwrap();
            let view = Py::new(py, ConstraintsView::new(owner)).unwrap();
            let arg = ConstraintsLike::View(view);

            assert_eq!(arg.to_rust(py).unwrap(), expected);
        });
    }

    #[rstest]
    #[case::empty(Vec::new())]
    #[case::populated(vec![
        GraphIrConstraint::Atom(GraphIrAtomId(1), GraphIrAtomConstraintForm::degree(2)),
        GraphIrConstraint::Molecule(GraphIrMoleculeConstraint::Connected { atoms: None }),
    ])]
    fn test_constraints_new(#[case] entries: Vec<GraphIrConstraint>) {
        Python::attach(|py| {
            let values = entries
                .iter()
                .map(|entry| {
                    into_py_variant(py, Constraint::from_rust(py, entry).unwrap()).unwrap()
                })
                .collect();
            let constraints = Constraints::new(py, values);

            assert_eq!(constraints.inner().as_slice(), entries.as_slice());
        });
    }

    #[rstest]
    #[case::equal(
        GraphIrConstraints::from(vec![GraphIrConstraint::And(Vec::new())]),
        GraphIrConstraints::from(vec![GraphIrConstraint::And(Vec::new())]),
        true,
    )]
    #[case::different(
        GraphIrConstraints::from(vec![GraphIrConstraint::And(Vec::new())]),
        GraphIrConstraints::from(vec![GraphIrConstraint::Or(Vec::new())]),
        false,
    )]
    fn test_constraints_eq(
        #[case] left: GraphIrConstraints,
        #[case] right: GraphIrConstraints,
        #[case] expected: bool,
    ) {
        assert_eq!(
            Constraints::from_inner(left) == Constraints::from_inner(right),
            expected
        );
    }

    #[rstest]
    fn test_constraints_repr() {
        Python::attach(|py| {
            let constraints = Constraints::from_inner(GraphIrConstraints::from(vec![
                GraphIrConstraint::Atom(GraphIrAtomId(1), GraphIrAtomConstraintForm::degree(2)),
                GraphIrConstraint::Or(Vec::new()),
            ]));

            assert_eq!(
                constraints.__repr__(py).unwrap(),
                "Constraints([Constraint.Atom(1, AtomConstraintForm.Degree(NumForm.Lit(2))), Constraint.Or([])])"
            );
        });
    }

    #[rstest]
    fn test_constraints_append() {
        let constraint =
            GraphIrConstraint::Molecule(GraphIrMoleculeConstraint::Connected { atoms: None });
        Python::attach(|py| {
            let mut constraints =
                Constraints::from_inner(GraphIrConstraints::from(vec![constraint.clone()]));
            let value =
                into_py_variant(py, Constraint::from_rust(py, &constraint).unwrap()).unwrap();

            constraints.append(py, value);

            assert_eq!(
                constraints.inner().as_slice(),
                &[constraint.clone(), constraint]
            );
        });
    }

    #[rstest]
    fn test_constraints_clear() {
        let mut constraints =
            Constraints::from_inner(GraphIrConstraints::from(vec![GraphIrConstraint::And(
                Vec::new(),
            )]));

        constraints.clear();

        assert_eq!(constraints.inner(), &GraphIrConstraints::new());
    }

    #[rstest]
    fn test_constraints_update() {
        let initial = GraphIrConstraint::And(Vec::new());
        let from_container = GraphIrConstraint::Or(Vec::new());
        let from_view = GraphIrConstraint::Not(Box::new(GraphIrConstraint::And(Vec::new())));
        let from_entries =
            GraphIrConstraint::Molecule(GraphIrMoleculeConstraint::Connected { atoms: None });

        Python::attach(|py| {
            let target = Py::new(
                py,
                Constraints::from_inner(GraphIrConstraints::from(vec![initial.clone()])),
            )
            .unwrap();
            let container = Py::new(
                py,
                Constraints::from_inner(GraphIrConstraints::from(vec![from_container.clone()])),
            )
            .unwrap();
            let mut molecule = GraphIrMolecule::new();
            molecule.constraints_mut().push(from_view.clone());
            let view = Py::new(
                py,
                ConstraintsView::new(Py::new(py, MoleculeAst::from_rust(molecule)).unwrap()),
            )
            .unwrap();
            let entry =
                into_py_variant(py, Constraint::from_rust(py, &from_entries).unwrap()).unwrap();

            Constraints::update(
                target.clone_ref(py),
                py,
                ConstraintsUpdate::Container(container),
            )
            .unwrap();
            Constraints::update(target.clone_ref(py), py, ConstraintsUpdate::View(view)).unwrap();
            Constraints::update(
                target.clone_ref(py),
                py,
                ConstraintsUpdate::Entries(vec![entry]),
            )
            .unwrap();

            assert_eq!(
                target.bind(py).borrow().inner().as_slice(),
                &[initial, from_container, from_view, from_entries]
            );
        });
    }

    #[rstest]
    fn test_constraints_update_self() {
        let entry = GraphIrConstraint::And(Vec::new());

        Python::attach(|py| {
            let target = Py::new(
                py,
                Constraints::from_inner(GraphIrConstraints::from(vec![entry.clone()])),
            )
            .unwrap();

            Constraints::update(
                target.clone_ref(py),
                py,
                ConstraintsUpdate::Container(target.clone_ref(py)),
            )
            .unwrap();

            assert_eq!(
                target.bind(py).borrow().inner().as_slice(),
                &[entry.clone(), entry]
            );
        });
    }

    #[rstest]
    #[case::empty(GraphIrConstraints::new(), 0)]
    #[case::populated(GraphIrConstraints::from(vec![
        GraphIrConstraint::And(Vec::new()),
        GraphIrConstraint::Or(Vec::new()),
    ]), 2)]
    fn test_constraints_len(#[case] inner: GraphIrConstraints, #[case] expected: usize) {
        assert_eq!(Constraints::from_inner(inner).__len__(), expected);
    }

    #[rstest]
    #[case::positive(
        0,
        GraphIrConstraint::Atom(GraphIrAtomId(1), GraphIrAtomConstraintForm::degree(2),)
    )]
    #[case::negative(-1, GraphIrConstraint::Molecule(
        GraphIrMoleculeConstraint::Connected { atoms: None },
    ))]
    fn test_constraints_getitem(#[case] index: isize, #[case] expected: GraphIrConstraint) {
        Python::attach(|py| {
            let constraints = Constraints::from_inner(GraphIrConstraints::from(vec![
                GraphIrConstraint::Atom(GraphIrAtomId(1), GraphIrAtomConstraintForm::degree(2)),
                GraphIrConstraint::Molecule(GraphIrMoleculeConstraint::Connected { atoms: None }),
            ]));
            let actual = constraints.__getitem__(py, index).unwrap();

            assert_eq!(actual.to_rust(py), expected);
        });
    }

    #[rstest]
    #[case::positive(1)]
    #[case::negative(-2)]
    fn test_constraints_getitem_error(#[case] index: isize) {
        Python::attach(|py| {
            let constraints =
                Constraints::from_inner(GraphIrConstraints::from(vec![GraphIrConstraint::And(
                    Vec::new(),
                )]));

            assert_eq!(
                constraints
                    .__getitem__(py, index)
                    .err()
                    .unwrap()
                    .to_string(),
                "IndexError: constraint index out of range"
            );
        });
    }

    #[rstest]
    fn test_constraints_iter() {
        let first = GraphIrConstraint::And(Vec::new());
        let second = GraphIrConstraint::Or(Vec::new());
        let mut constraints = Constraints::from_inner(GraphIrConstraints::from(vec![
            first.clone(),
            second.clone(),
        ]));

        Python::attach(|py| {
            let mut iter = constraints.__iter__(py).unwrap();
            constraints
                .inner_mut()
                .push(GraphIrConstraint::Not(Box::new(GraphIrConstraint::And(
                    Vec::new(),
                ))));

            assert_eq!(
                iter.__next__().unwrap().bind(py).borrow().to_rust(py),
                first
            );
            assert_eq!(
                iter.__next__().unwrap().bind(py).borrow().to_rust(py),
                second
            );
            assert!(iter.__next__().is_none());
        });
    }

    #[rstest]
    fn test_constraints_from_inner() {
        let entries = vec![GraphIrConstraint::Or(Vec::new())];

        assert_eq!(
            Constraints::from_inner(GraphIrConstraints::from(entries.clone()))
                .inner()
                .as_slice(),
            entries.as_slice()
        );
    }

    #[rstest]
    fn test_constraints_view_repr() {
        let mut molecule = GraphIrMolecule::new();
        molecule
            .constraints_mut()
            .push(GraphIrConstraint::And(Vec::new()));
        molecule
            .constraints_mut()
            .push(GraphIrConstraint::Or(Vec::new()));

        Python::attach(|py| {
            let owner = Py::new(py, MoleculeAst::from_rust(molecule)).unwrap();
            let view = ConstraintsView::new(owner);

            assert_eq!(view.__repr__(py).unwrap(), "ConstraintsView(2 entries)");
        });
    }

    #[rstest]
    fn test_constraints_view_append() {
        let constraint =
            GraphIrConstraint::Molecule(GraphIrMoleculeConstraint::Connected { atoms: None });
        let mut molecule = GraphIrMolecule::new();
        molecule.constraints_mut().push(constraint.clone());

        Python::attach(|py| {
            let owner = Py::new(py, MoleculeAst::from_rust(molecule)).unwrap();
            let view = ConstraintsView::new(owner.clone_ref(py));
            let value =
                into_py_variant(py, Constraint::from_rust(py, &constraint).unwrap()).unwrap();

            view.append(py, value);

            assert_eq!(
                owner.bind(py).borrow().inner().constraints().as_slice(),
                &[constraint.clone(), constraint]
            );
        });
    }

    #[rstest]
    fn test_constraints_view_clear() {
        let mut molecule = GraphIrMolecule::new();
        molecule
            .constraints_mut()
            .push(GraphIrConstraint::And(Vec::new()));

        Python::attach(|py| {
            let owner = Py::new(py, MoleculeAst::from_rust(molecule)).unwrap();
            let view = ConstraintsView::new(owner.clone_ref(py));

            view.clear(py);

            assert_eq!(
                owner.bind(py).borrow().inner().constraints(),
                &GraphIrConstraints::new()
            );
        });
    }

    #[rstest]
    fn test_constraints_view_update() {
        let initial = GraphIrConstraint::And(Vec::new());
        let from_container = GraphIrConstraint::Or(Vec::new());
        let from_view = GraphIrConstraint::Not(Box::new(GraphIrConstraint::And(Vec::new())));
        let from_entries =
            GraphIrConstraint::Molecule(GraphIrMoleculeConstraint::Connected { atoms: None });
        let mut target_molecule = GraphIrMolecule::new();
        target_molecule.constraints_mut().push(initial.clone());

        Python::attach(|py| {
            let target_owner = Py::new(py, MoleculeAst::from_rust(target_molecule)).unwrap();
            let target = ConstraintsView::new(target_owner.clone_ref(py));
            let container = Py::new(
                py,
                Constraints::from_inner(GraphIrConstraints::from(vec![from_container.clone()])),
            )
            .unwrap();
            let mut source_molecule = GraphIrMolecule::new();
            source_molecule.constraints_mut().push(from_view.clone());
            let source_view = Py::new(
                py,
                ConstraintsView::new(Py::new(py, MoleculeAst::from_rust(source_molecule)).unwrap()),
            )
            .unwrap();
            let entry =
                into_py_variant(py, Constraint::from_rust(py, &from_entries).unwrap()).unwrap();

            target
                .update(py, ConstraintsUpdate::Container(container))
                .unwrap();
            target
                .update(py, ConstraintsUpdate::View(source_view))
                .unwrap();
            target
                .update(py, ConstraintsUpdate::Entries(vec![entry]))
                .unwrap();

            assert_eq!(
                target_owner
                    .bind(py)
                    .borrow()
                    .inner()
                    .constraints()
                    .as_slice(),
                &[initial, from_container, from_view, from_entries]
            );
        });
    }

    #[rstest]
    fn test_constraints_view_update_self() {
        let entry = GraphIrConstraint::Or(Vec::new());
        let mut molecule = GraphIrMolecule::new();
        molecule.constraints_mut().push(entry.clone());

        Python::attach(|py| {
            let owner = Py::new(py, MoleculeAst::from_rust(molecule)).unwrap();
            let view = Py::new(py, ConstraintsView::new(owner.clone_ref(py))).unwrap();

            view.bind(py)
                .borrow()
                .update(py, ConstraintsUpdate::View(view.clone_ref(py)))
                .unwrap();

            assert_eq!(
                owner.bind(py).borrow().inner().constraints().as_slice(),
                &[entry.clone(), entry]
            );
        });
    }

    #[rstest]
    fn test_constraints_view_len() {
        Python::attach(|py| {
            let owner = Py::new(py, MoleculeAst::from_rust(GraphIrMolecule::new())).unwrap();
            let view = ConstraintsView::new(owner.clone_ref(py));
            assert_eq!(view.__len__(py).unwrap(), 0);

            owner
                .borrow_mut(py)
                .inner_mut()
                .constraints_mut()
                .push(GraphIrConstraint::And(Vec::new()));

            assert_eq!(view.__len__(py).unwrap(), 1);
        });
    }

    #[rstest]
    #[case::positive(
        0,
        GraphIrConstraint::Atom(GraphIrAtomId(1), GraphIrAtomConstraintForm::degree(2),)
    )]
    #[case::negative(-1, GraphIrConstraint::Molecule(
        GraphIrMoleculeConstraint::Connected { atoms: None },
    ))]
    fn test_constraints_view_getitem(#[case] index: isize, #[case] expected: GraphIrConstraint) {
        let mut molecule = GraphIrMolecule::new();
        molecule.constraints_mut().push(GraphIrConstraint::Atom(
            GraphIrAtomId(1),
            GraphIrAtomConstraintForm::degree(2),
        ));
        molecule.constraints_mut().push(GraphIrConstraint::Molecule(
            GraphIrMoleculeConstraint::Connected { atoms: None },
        ));

        Python::attach(|py| {
            let owner = Py::new(py, MoleculeAst::from_rust(molecule)).unwrap();
            let view = ConstraintsView::new(owner);

            assert_eq!(view.__getitem__(py, index).unwrap().to_rust(py), expected);
        });
    }

    #[rstest]
    #[case::positive(1)]
    #[case::negative(-2)]
    fn test_constraints_view_getitem_error(#[case] index: isize) {
        let mut molecule = GraphIrMolecule::new();
        molecule
            .constraints_mut()
            .push(GraphIrConstraint::And(Vec::new()));

        Python::attach(|py| {
            let owner = Py::new(py, MoleculeAst::from_rust(molecule)).unwrap();
            let view = ConstraintsView::new(owner);

            assert_eq!(
                view.__getitem__(py, index).err().unwrap().to_string(),
                "IndexError: constraint index out of range"
            );
        });
    }

    #[rstest]
    fn test_constraints_view_iter() {
        let first = GraphIrConstraint::And(Vec::new());
        let second = GraphIrConstraint::Or(Vec::new());
        let mut molecule = GraphIrMolecule::new();
        molecule.constraints_mut().push(first.clone());
        molecule.constraints_mut().push(second.clone());

        Python::attach(|py| {
            let owner = Py::new(py, MoleculeAst::from_rust(molecule)).unwrap();
            let view = ConstraintsView::new(owner.clone_ref(py));
            let mut iter = view.__iter__(py).unwrap();
            owner
                .borrow_mut(py)
                .inner_mut()
                .constraints_mut()
                .push(GraphIrConstraint::Not(Box::new(GraphIrConstraint::And(
                    Vec::new(),
                ))));

            assert_eq!(
                iter.__next__().unwrap().bind(py).borrow().to_rust(py),
                first
            );
            assert_eq!(
                iter.__next__().unwrap().bind(py).borrow().to_rust(py),
                second
            );
            assert!(iter.__next__().is_none());
        });
    }
}
