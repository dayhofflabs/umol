//! Molecule-level constraint payloads matching `umol_ast::ast::constraint`.

use std::vec::IntoIter;

use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use umol_ast::ast::{
    AromaticSystemId as AstAromaticSystemId, AtomId as AstAtomId, BondId as AstBondId,
    Constraint as AstConstraint, Constraints as AstConstraints, DativeBondId as AstDativeBondId,
    MoleculeConstraint as AstMoleculeConstraint, MulticenterBondId as AstMulticenterBondId,
    NoncovalentBondId as AstNoncovalentBondId, RelationalConstraint as AstRelationalConstraint,
    StereoAtomId as AstStereoAtomId, StereoBondId as AstStereoBondId,
    SubPatternAnchor as AstSubPatternAnchor,
};

use super::aromatic::AromaticSystemConstraintAst;
use super::atom::AtomConstraintAst;
use super::bond::BondConstraintAst;
use super::dative::DativeBondConstraintAst;
use super::multicenter::MulticenterBondConstraintAst;
use super::noncovalent::NoncovalentBondConstraintAst;
use super::stereo::{StereoAtomConstraintAst, StereoBondConstraintAst};
use crate::convert::{into_py_variant, variant_repr};
use crate::molecule::MoleculeAst;
use crate::spin::UnpairedElectronsAst;
use crate::stereo::StereoKind;
use crate::value::ValueAst;

/// Entity correspondences anchoring a subpattern match. Each collection holds
/// `(target, pattern)` integer-id pairs for one molecule entity family.
#[pyclass(eq, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SubPatternAnchor {
    atoms: Vec<(u32, u32)>,
    bonds: Vec<(u32, u32)>,
    dative_bonds: Vec<(u32, u32)>,
    aromatic_systems: Vec<(u32, u32)>,
    multicenter_bonds: Vec<(u32, u32)>,
    noncovalent_bonds: Vec<(u32, u32)>,
    stereo_atoms: Vec<(u32, u32)>,
    stereo_bonds: Vec<(u32, u32)>,
}

/// A cross-entity molecule constraint covering dative bonds, aromatic systems,
/// multicenter bonds, noncovalent bonds, stereo atoms, and stereo bonds.
#[pyclass]
pub enum RelationalConstraint {
    DativeBondDonors(u32, Vec<u32>),
    DativeBondDonor(u32, u32),
    DativeBondContainsAllDonors(u32, Vec<u32>),
    DativeBondAllDonors(u32, Py<AtomConstraintAst>),
    DativeBondAnyDonor(u32, Py<AtomConstraintAst>),
    DativeBondAcceptor(u32, u32),
    DativeBondAcceptorSatisfies(u32, Py<AtomConstraintAst>),
    DativeBondParallels(u32, u32),
    AromaticSystemAtoms(u32, Vec<u32>),
    AromaticSystemContains(u32, u32),
    AromaticSystemContainsAll(u32, Vec<u32>),
    AromaticSystemAllAtoms(u32, Py<AtomConstraintAst>),
    AromaticSystemAnyAtom(u32, Py<AtomConstraintAst>),
    MulticenterBondAtoms(u32, Vec<u32>),
    MulticenterBondContains(u32, u32),
    MulticenterBondContainsAll(u32, Vec<u32>),
    MulticenterBondAllAtoms(u32, Py<AtomConstraintAst>),
    MulticenterBondAnyAtom(u32, Py<AtomConstraintAst>),
    NoncovalentBondEnds(u32, [u32; 2]),
    NoncovalentBondContains(u32, u32),
    NoncovalentBondEndsSatisfy(u32, [Py<AtomConstraintAst>; 2]),
    StereoAtomSite(u32, u32),
    StereoAtomContains(u32, u32),
    StereoAtomLigands(u32, Vec<u32>),
    StereoAtomAllLigands(u32, Py<AtomConstraintAst>),
    StereoAtomAnyLigand(u32, Py<AtomConstraintAst>),
    StereoBondSite(u32, u32),
    StereoBondContains(u32, u32),
    StereoBondLigands(u32, Vec<u32>),
    StereoBondAllLigands(u32, Py<AtomConstraintAst>),
    StereoBondAnyLigand(u32, Py<AtomConstraintAst>),
}

/// A molecule-scope predicate over values, connectivity, or a nested pattern.
#[pyclass]
pub enum MoleculeConstraint {
    ChargeSum(Option<Vec<u32>>, Py<ValueAst>),
    #[pyo3(constructor = (atoms, unpaired_electrons))]
    UnpairedElectronCoupling {
        atoms: Option<Vec<u32>>,
        unpaired_electrons: Py<UnpairedElectronsAst>,
    },
    BondOrderSum(Option<Vec<u32>>, Py<ValueAst>),
    Connected(Option<Vec<u32>>),
    SubPattern(Py<SubPatternAnchor>, Py<MoleculeAst>),
}

/// A recursive molecule-constraint tree containing entity leaves, aggregate
/// leaves, and Boolean combinators.
#[pyclass]
pub enum Constraint {
    Atom(u32, Py<AtomConstraintAst>),
    Bond(u32, Py<BondConstraintAst>),
    DativeBond(u32, Py<DativeBondConstraintAst>),
    AromaticSystem(u32, Py<AromaticSystemConstraintAst>),
    MulticenterBond(u32, Py<MulticenterBondConstraintAst>),
    NoncovalentBond(u32, Py<NoncovalentBondConstraintAst>),
    StereoAtom(u32, StereoKind, Py<StereoAtomConstraintAst>),
    StereoBond(u32, StereoKind, Py<StereoBondConstraintAst>),
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

impl Constraint {
    pub(crate) fn from_rust(py: Python<'_>, constraint: &AstConstraint) -> PyResult<Self> {
        Ok(match constraint {
            AstConstraint::Atom(id, child) => Self::Atom(
                id.0,
                into_py_variant(py, AtomConstraintAst::from_rust(py, child)?)?,
            ),
            AstConstraint::Bond(id, child) => Self::Bond(
                id.0,
                into_py_variant(py, BondConstraintAst::from_rust(py, child)?)?,
            ),
            AstConstraint::DativeBond(id, child) => Self::DativeBond(
                id.0,
                into_py_variant(py, DativeBondConstraintAst::from_rust(py, child)?)?,
            ),
            AstConstraint::AromaticSystem(id, child) => Self::AromaticSystem(
                id.0,
                into_py_variant(py, AromaticSystemConstraintAst::from_rust(py, child)?)?,
            ),
            AstConstraint::MulticenterBond(id, child) => Self::MulticenterBond(
                id.0,
                into_py_variant(py, MulticenterBondConstraintAst::from_rust(py, child)?)?,
            ),
            AstConstraint::NoncovalentBond(id, child) => Self::NoncovalentBond(
                id.0,
                into_py_variant(py, NoncovalentBondConstraintAst::from_rust(py, child)?)?,
            ),
            AstConstraint::StereoAtom(id, kind, child) => Self::StereoAtom(
                id.0,
                StereoKind::from_rust(*kind),
                into_py_variant(py, StereoAtomConstraintAst::from_rust(py, child)?)?,
            ),
            AstConstraint::StereoBond(id, kind, child) => Self::StereoBond(
                id.0,
                StereoKind::from_rust(*kind),
                into_py_variant(py, StereoBondConstraintAst::from_rust(py, child)?)?,
            ),
            AstConstraint::Relational(child) => Self::Relational(into_py_variant(
                py,
                RelationalConstraint::from_rust(py, child)?,
            )?),
            AstConstraint::Molecule(child) => Self::Molecule(into_py_variant(
                py,
                MoleculeConstraint::from_rust(py, child)?,
            )?),
            AstConstraint::And(children) => Self::And(
                children
                    .iter()
                    .map(|child| into_py_variant(py, Self::from_rust(py, child)?))
                    .collect::<PyResult<_>>()?,
            ),
            AstConstraint::Or(children) => Self::Or(
                children
                    .iter()
                    .map(|child| into_py_variant(py, Self::from_rust(py, child)?))
                    .collect::<PyResult<_>>()?,
            ),
            AstConstraint::Not(child) => {
                Self::Not(into_py_variant(py, Self::from_rust(py, child)?)?)
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstConstraint {
        match self {
            Self::Atom(id, child) => {
                AstConstraint::Atom(AstAtomId(*id), child.bind(py).borrow().to_rust(py))
            }
            Self::Bond(id, child) => {
                AstConstraint::Bond(AstBondId(*id), child.bind(py).borrow().to_rust(py))
            }
            Self::DativeBond(id, child) => {
                AstConstraint::DativeBond(AstDativeBondId(*id), child.bind(py).borrow().to_rust(py))
            }
            Self::AromaticSystem(id, child) => AstConstraint::AromaticSystem(
                AstAromaticSystemId(*id),
                child.bind(py).borrow().to_rust(py),
            ),
            Self::MulticenterBond(id, child) => AstConstraint::MulticenterBond(
                AstMulticenterBondId(*id),
                child.bind(py).borrow().to_rust(py),
            ),
            Self::NoncovalentBond(id, child) => AstConstraint::NoncovalentBond(
                AstNoncovalentBondId(*id),
                child.bind(py).borrow().to_rust(py),
            ),
            Self::StereoAtom(id, kind, child) => AstConstraint::StereoAtom(
                AstStereoAtomId(*id),
                kind.to_rust(),
                child.bind(py).borrow().to_rust(py),
            ),
            Self::StereoBond(id, kind, child) => AstConstraint::StereoBond(
                AstStereoBondId(*id),
                kind.to_rust(),
                child.bind(py).borrow().to_rust(py),
            ),
            Self::Relational(child) => {
                AstConstraint::Relational(child.bind(py).borrow().to_rust(py))
            }
            Self::Molecule(child) => AstConstraint::Molecule(child.bind(py).borrow().to_rust(py)),
            Self::And(children) => AstConstraint::And(
                children
                    .iter()
                    .map(|child| child.bind(py).borrow().to_rust(py))
                    .collect(),
            ),
            Self::Or(children) => AstConstraint::Or(
                children
                    .iter()
                    .map(|child| child.bind(py).borrow().to_rust(py))
                    .collect(),
            ),
            Self::Not(child) => AstConstraint::Not(Box::new(child.bind(py).borrow().to_rust(py))),
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
fn constraint_iter(py: Python<'_>, constraints: &AstConstraints) -> PyResult<ConstraintIter> {
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
    Overlay(AstConstraints),
    Entries(Vec<AstConstraint>),
}

impl ResolvedConstraintsUpdate {
    /// Append the resolved entries in order, preserving duplicates.
    pub(crate) fn apply(self, target: &mut AstConstraints) {
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
pub(crate) enum ConstraintsArg {
    Container(Py<Constraints>),
    View(Py<ConstraintsView>),
}

impl ConstraintsArg {
    pub(crate) fn to_rust(&self, py: Python<'_>) -> PyResult<AstConstraints> {
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
pub struct Constraints(AstConstraints);

#[pymethods]
impl Constraints {
    /// Build an owned container from constraint entries, preserving order and duplicates.
    #[new]
    fn new(py: Python<'_>, entries: Vec<Py<Constraint>>) -> Self {
        Self(AstConstraints::from(
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
    pub(crate) fn inner(&self) -> &AstConstraints {
        &self.0
    }

    pub(crate) fn inner_mut(&mut self) -> &mut AstConstraints {
        &mut self.0
    }

    #[allow(
        dead_code,
        reason = "AST conversion for the unregistered Constraints pyclass"
    )]
    pub(crate) fn from_inner(constraints: AstConstraints) -> Self {
        Self(constraints)
    }
}

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
        f: impl FnOnce(&AstConstraints) -> PyResult<R>,
    ) -> PyResult<R> {
        let molecule = self.owner.bind(py).borrow();
        f(molecule.inner().constraints())
    }

    /// Mutate the molecule's constraint store in place through `f`.
    pub(crate) fn with_mut<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut AstConstraints) -> R,
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
        self.with_mut(py, AstConstraints::clear);
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
            Self::SubPattern(_, _) => {
                variant_repr(slf.bind(py).as_any(), "MoleculeConstraint", "SubPattern", 2)
            }
        }
    }
}

impl MoleculeConstraint {
    pub(crate) fn from_rust(py: Python<'_>, constraint: &AstMoleculeConstraint) -> PyResult<Self> {
        Ok(match constraint {
            AstMoleculeConstraint::ChargeSum { atoms, sum } => Self::ChargeSum(
                atoms
                    .as_ref()
                    .map(|atoms| atoms.iter().map(|atom| atom.0).collect()),
                into_py_variant(py, ValueAst::from_rust(py, sum)?)?,
            ),
            AstMoleculeConstraint::UnpairedElectronCoupling {
                atoms,
                unpaired_electrons,
            } => Self::UnpairedElectronCoupling {
                atoms: atoms
                    .as_ref()
                    .map(|atoms| atoms.iter().map(|atom| atom.0).collect()),
                unpaired_electrons: Py::new(
                    py,
                    UnpairedElectronsAst::from_rust(py, unpaired_electrons)?,
                )?,
            },
            AstMoleculeConstraint::BondOrderSum { bonds, sum } => Self::BondOrderSum(
                bonds
                    .as_ref()
                    .map(|bonds| bonds.iter().map(|bond| bond.0).collect()),
                into_py_variant(py, ValueAst::from_rust(py, sum)?)?,
            ),
            AstMoleculeConstraint::Connected { atoms } => Self::Connected(
                atoms
                    .as_ref()
                    .map(|atoms| atoms.iter().map(|atom| atom.0).collect()),
            ),
            AstMoleculeConstraint::SubPattern { anchor, pattern } => Self::SubPattern(
                Py::new(py, SubPatternAnchor::from_rust(anchor))?,
                Py::new(py, MoleculeAst::from_inner((**pattern).clone()))?,
            ),
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstMoleculeConstraint {
        match self {
            Self::ChargeSum(atoms, sum) => AstMoleculeConstraint::ChargeSum {
                atoms: atoms
                    .as_ref()
                    .map(|atoms| atoms.iter().copied().map(AstAtomId).collect()),
                sum: sum.bind(py).borrow().to_rust(py),
            },
            Self::UnpairedElectronCoupling {
                atoms,
                unpaired_electrons,
            } => AstMoleculeConstraint::UnpairedElectronCoupling {
                atoms: atoms
                    .as_ref()
                    .map(|atoms| atoms.iter().copied().map(AstAtomId).collect()),
                unpaired_electrons: unpaired_electrons.bind(py).borrow().to_rust(py),
            },
            Self::BondOrderSum(bonds, sum) => AstMoleculeConstraint::BondOrderSum {
                bonds: bonds
                    .as_ref()
                    .map(|bonds| bonds.iter().copied().map(AstBondId).collect()),
                sum: sum.bind(py).borrow().to_rust(py),
            },
            Self::Connected(atoms) => AstMoleculeConstraint::Connected {
                atoms: atoms
                    .as_ref()
                    .map(|atoms| atoms.iter().copied().map(AstAtomId).collect()),
            },
            Self::SubPattern(anchor, pattern) => AstMoleculeConstraint::SubPattern {
                anchor: anchor.bind(py).borrow().to_rust(),
                pattern: Box::new(pattern.bind(py).borrow().inner().clone()),
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
        constraint: &AstRelationalConstraint,
    ) -> PyResult<Self> {
        Ok(match constraint {
            AstRelationalConstraint::DativeBondDonors { bond, atoms } => {
                Self::DativeBondDonors(bond.0, atoms.iter().map(|atom| atom.0).collect())
            }
            AstRelationalConstraint::DativeBondDonor { bond, atom } => {
                Self::DativeBondDonor(bond.0, atom.0)
            }
            AstRelationalConstraint::DativeBondContainsAllDonors { bond, atoms } => {
                Self::DativeBondContainsAllDonors(bond.0, atoms.iter().map(|atom| atom.0).collect())
            }
            AstRelationalConstraint::DativeBondAllDonors { bond, predicate } => {
                Self::DativeBondAllDonors(
                    bond.0,
                    into_py_variant(py, AtomConstraintAst::from_rust(py, predicate)?)?,
                )
            }
            AstRelationalConstraint::DativeBondAnyDonor { bond, predicate } => {
                Self::DativeBondAnyDonor(
                    bond.0,
                    into_py_variant(py, AtomConstraintAst::from_rust(py, predicate)?)?,
                )
            }
            AstRelationalConstraint::DativeBondAcceptor { bond, atom } => {
                Self::DativeBondAcceptor(bond.0, atom.0)
            }
            AstRelationalConstraint::DativeBondAcceptorSatisfies { bond, predicate } => {
                Self::DativeBondAcceptorSatisfies(
                    bond.0,
                    into_py_variant(py, AtomConstraintAst::from_rust(py, predicate)?)?,
                )
            }
            AstRelationalConstraint::DativeBondParallels { dative, parallel } => {
                Self::DativeBondParallels(dative.0, parallel.0)
            }
            AstRelationalConstraint::AromaticSystemAtoms { system, atoms } => {
                Self::AromaticSystemAtoms(system.0, atoms.iter().map(|atom| atom.0).collect())
            }
            AstRelationalConstraint::AromaticSystemContains { system, atom } => {
                Self::AromaticSystemContains(system.0, atom.0)
            }
            AstRelationalConstraint::AromaticSystemContainsAll { system, atoms } => {
                Self::AromaticSystemContainsAll(system.0, atoms.iter().map(|atom| atom.0).collect())
            }
            AstRelationalConstraint::AromaticSystemAllAtoms { system, predicate } => {
                Self::AromaticSystemAllAtoms(
                    system.0,
                    into_py_variant(py, AtomConstraintAst::from_rust(py, predicate)?)?,
                )
            }
            AstRelationalConstraint::AromaticSystemAnyAtom { system, predicate } => {
                Self::AromaticSystemAnyAtom(
                    system.0,
                    into_py_variant(py, AtomConstraintAst::from_rust(py, predicate)?)?,
                )
            }
            AstRelationalConstraint::MulticenterBondAtoms { bond, atoms } => {
                Self::MulticenterBondAtoms(bond.0, atoms.iter().map(|atom| atom.0).collect())
            }
            AstRelationalConstraint::MulticenterBondContains { bond, atom } => {
                Self::MulticenterBondContains(bond.0, atom.0)
            }
            AstRelationalConstraint::MulticenterBondContainsAll { bond, atoms } => {
                Self::MulticenterBondContainsAll(bond.0, atoms.iter().map(|atom| atom.0).collect())
            }
            AstRelationalConstraint::MulticenterBondAllAtoms { bond, predicate } => {
                Self::MulticenterBondAllAtoms(
                    bond.0,
                    into_py_variant(py, AtomConstraintAst::from_rust(py, predicate)?)?,
                )
            }
            AstRelationalConstraint::MulticenterBondAnyAtom { bond, predicate } => {
                Self::MulticenterBondAnyAtom(
                    bond.0,
                    into_py_variant(py, AtomConstraintAst::from_rust(py, predicate)?)?,
                )
            }
            AstRelationalConstraint::NoncovalentBondEnds { bond, atoms } => {
                Self::NoncovalentBondEnds(bond.0, [atoms[0].0, atoms[1].0])
            }
            AstRelationalConstraint::NoncovalentBondContains { bond, atom } => {
                Self::NoncovalentBondContains(bond.0, atom.0)
            }
            AstRelationalConstraint::NoncovalentBondEndsSatisfy { bond, predicates } => {
                Self::NoncovalentBondEndsSatisfy(
                    bond.0,
                    [
                        into_py_variant(py, AtomConstraintAst::from_rust(py, &predicates[0])?)?,
                        into_py_variant(py, AtomConstraintAst::from_rust(py, &predicates[1])?)?,
                    ],
                )
            }
            AstRelationalConstraint::StereoAtomSite { stereo_atom, atom } => {
                Self::StereoAtomSite(stereo_atom.0, atom.0)
            }
            AstRelationalConstraint::StereoAtomContains { stereo_atom, atom } => {
                Self::StereoAtomContains(stereo_atom.0, atom.0)
            }
            AstRelationalConstraint::StereoAtomLigands { stereo_atom, atoms } => {
                Self::StereoAtomLigands(stereo_atom.0, atoms.iter().map(|atom| atom.0).collect())
            }
            AstRelationalConstraint::StereoAtomAllLigands {
                stereo_atom,
                predicate,
            } => Self::StereoAtomAllLigands(
                stereo_atom.0,
                into_py_variant(py, AtomConstraintAst::from_rust(py, predicate)?)?,
            ),
            AstRelationalConstraint::StereoAtomAnyLigand {
                stereo_atom,
                predicate,
            } => Self::StereoAtomAnyLigand(
                stereo_atom.0,
                into_py_variant(py, AtomConstraintAst::from_rust(py, predicate)?)?,
            ),
            AstRelationalConstraint::StereoBondSite { stereo_bond, bond } => {
                Self::StereoBondSite(stereo_bond.0, bond.0)
            }
            AstRelationalConstraint::StereoBondContains { stereo_bond, atom } => {
                Self::StereoBondContains(stereo_bond.0, atom.0)
            }
            AstRelationalConstraint::StereoBondLigands { stereo_bond, atoms } => {
                Self::StereoBondLigands(stereo_bond.0, atoms.iter().map(|atom| atom.0).collect())
            }
            AstRelationalConstraint::StereoBondAllLigands {
                stereo_bond,
                predicate,
            } => Self::StereoBondAllLigands(
                stereo_bond.0,
                into_py_variant(py, AtomConstraintAst::from_rust(py, predicate)?)?,
            ),
            AstRelationalConstraint::StereoBondAnyLigand {
                stereo_bond,
                predicate,
            } => Self::StereoBondAnyLigand(
                stereo_bond.0,
                into_py_variant(py, AtomConstraintAst::from_rust(py, predicate)?)?,
            ),
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstRelationalConstraint {
        match self {
            Self::DativeBondDonors(bond, atoms) => AstRelationalConstraint::DativeBondDonors {
                bond: AstDativeBondId(*bond),
                atoms: atoms.iter().copied().map(AstAtomId).collect(),
            },
            Self::DativeBondDonor(bond, atom) => AstRelationalConstraint::DativeBondDonor {
                bond: AstDativeBondId(*bond),
                atom: AstAtomId(*atom),
            },
            Self::DativeBondContainsAllDonors(bond, atoms) => {
                AstRelationalConstraint::DativeBondContainsAllDonors {
                    bond: AstDativeBondId(*bond),
                    atoms: atoms.iter().copied().map(AstAtomId).collect(),
                }
            }
            Self::DativeBondAllDonors(bond, predicate) => {
                AstRelationalConstraint::DativeBondAllDonors {
                    bond: AstDativeBondId(*bond),
                    predicate: Box::new(predicate.bind(py).borrow().to_rust(py)),
                }
            }
            Self::DativeBondAnyDonor(bond, predicate) => {
                AstRelationalConstraint::DativeBondAnyDonor {
                    bond: AstDativeBondId(*bond),
                    predicate: Box::new(predicate.bind(py).borrow().to_rust(py)),
                }
            }
            Self::DativeBondAcceptor(bond, atom) => AstRelationalConstraint::DativeBondAcceptor {
                bond: AstDativeBondId(*bond),
                atom: AstAtomId(*atom),
            },
            Self::DativeBondAcceptorSatisfies(bond, predicate) => {
                AstRelationalConstraint::DativeBondAcceptorSatisfies {
                    bond: AstDativeBondId(*bond),
                    predicate: Box::new(predicate.bind(py).borrow().to_rust(py)),
                }
            }
            Self::DativeBondParallels(dative, parallel) => {
                AstRelationalConstraint::DativeBondParallels {
                    dative: AstDativeBondId(*dative),
                    parallel: AstBondId(*parallel),
                }
            }
            Self::AromaticSystemAtoms(system, atoms) => {
                AstRelationalConstraint::AromaticSystemAtoms {
                    system: AstAromaticSystemId(*system),
                    atoms: atoms.iter().copied().map(AstAtomId).collect(),
                }
            }
            Self::AromaticSystemContains(system, atom) => {
                AstRelationalConstraint::AromaticSystemContains {
                    system: AstAromaticSystemId(*system),
                    atom: AstAtomId(*atom),
                }
            }
            Self::AromaticSystemContainsAll(system, atoms) => {
                AstRelationalConstraint::AromaticSystemContainsAll {
                    system: AstAromaticSystemId(*system),
                    atoms: atoms.iter().copied().map(AstAtomId).collect(),
                }
            }
            Self::AromaticSystemAllAtoms(system, predicate) => {
                AstRelationalConstraint::AromaticSystemAllAtoms {
                    system: AstAromaticSystemId(*system),
                    predicate: Box::new(predicate.bind(py).borrow().to_rust(py)),
                }
            }
            Self::AromaticSystemAnyAtom(system, predicate) => {
                AstRelationalConstraint::AromaticSystemAnyAtom {
                    system: AstAromaticSystemId(*system),
                    predicate: Box::new(predicate.bind(py).borrow().to_rust(py)),
                }
            }
            Self::MulticenterBondAtoms(bond, atoms) => {
                AstRelationalConstraint::MulticenterBondAtoms {
                    bond: AstMulticenterBondId(*bond),
                    atoms: atoms.iter().copied().map(AstAtomId).collect(),
                }
            }
            Self::MulticenterBondContains(bond, atom) => {
                AstRelationalConstraint::MulticenterBondContains {
                    bond: AstMulticenterBondId(*bond),
                    atom: AstAtomId(*atom),
                }
            }
            Self::MulticenterBondContainsAll(bond, atoms) => {
                AstRelationalConstraint::MulticenterBondContainsAll {
                    bond: AstMulticenterBondId(*bond),
                    atoms: atoms.iter().copied().map(AstAtomId).collect(),
                }
            }
            Self::MulticenterBondAllAtoms(bond, predicate) => {
                AstRelationalConstraint::MulticenterBondAllAtoms {
                    bond: AstMulticenterBondId(*bond),
                    predicate: Box::new(predicate.bind(py).borrow().to_rust(py)),
                }
            }
            Self::MulticenterBondAnyAtom(bond, predicate) => {
                AstRelationalConstraint::MulticenterBondAnyAtom {
                    bond: AstMulticenterBondId(*bond),
                    predicate: Box::new(predicate.bind(py).borrow().to_rust(py)),
                }
            }
            Self::NoncovalentBondEnds(bond, atoms) => {
                AstRelationalConstraint::NoncovalentBondEnds {
                    bond: AstNoncovalentBondId(*bond),
                    atoms: [AstAtomId(atoms[0]), AstAtomId(atoms[1])],
                }
            }
            Self::NoncovalentBondContains(bond, atom) => {
                AstRelationalConstraint::NoncovalentBondContains {
                    bond: AstNoncovalentBondId(*bond),
                    atom: AstAtomId(*atom),
                }
            }
            Self::NoncovalentBondEndsSatisfy(bond, predicates) => {
                AstRelationalConstraint::NoncovalentBondEndsSatisfy {
                    bond: AstNoncovalentBondId(*bond),
                    predicates: [
                        Box::new(predicates[0].bind(py).borrow().to_rust(py)),
                        Box::new(predicates[1].bind(py).borrow().to_rust(py)),
                    ],
                }
            }
            Self::StereoAtomSite(stereo_atom, atom) => AstRelationalConstraint::StereoAtomSite {
                stereo_atom: AstStereoAtomId(*stereo_atom),
                atom: AstAtomId(*atom),
            },
            Self::StereoAtomContains(stereo_atom, atom) => {
                AstRelationalConstraint::StereoAtomContains {
                    stereo_atom: AstStereoAtomId(*stereo_atom),
                    atom: AstAtomId(*atom),
                }
            }
            Self::StereoAtomLigands(stereo_atom, atoms) => {
                AstRelationalConstraint::StereoAtomLigands {
                    stereo_atom: AstStereoAtomId(*stereo_atom),
                    atoms: atoms.iter().copied().map(AstAtomId).collect(),
                }
            }
            Self::StereoAtomAllLigands(stereo_atom, predicate) => {
                AstRelationalConstraint::StereoAtomAllLigands {
                    stereo_atom: AstStereoAtomId(*stereo_atom),
                    predicate: Box::new(predicate.bind(py).borrow().to_rust(py)),
                }
            }
            Self::StereoAtomAnyLigand(stereo_atom, predicate) => {
                AstRelationalConstraint::StereoAtomAnyLigand {
                    stereo_atom: AstStereoAtomId(*stereo_atom),
                    predicate: Box::new(predicate.bind(py).borrow().to_rust(py)),
                }
            }
            Self::StereoBondSite(stereo_bond, bond) => AstRelationalConstraint::StereoBondSite {
                stereo_bond: AstStereoBondId(*stereo_bond),
                bond: AstBondId(*bond),
            },
            Self::StereoBondContains(stereo_bond, atom) => {
                AstRelationalConstraint::StereoBondContains {
                    stereo_bond: AstStereoBondId(*stereo_bond),
                    atom: AstAtomId(*atom),
                }
            }
            Self::StereoBondLigands(stereo_bond, atoms) => {
                AstRelationalConstraint::StereoBondLigands {
                    stereo_bond: AstStereoBondId(*stereo_bond),
                    atoms: atoms.iter().copied().map(AstAtomId).collect(),
                }
            }
            Self::StereoBondAllLigands(stereo_bond, predicate) => {
                AstRelationalConstraint::StereoBondAllLigands {
                    stereo_bond: AstStereoBondId(*stereo_bond),
                    predicate: Box::new(predicate.bind(py).borrow().to_rust(py)),
                }
            }
            Self::StereoBondAnyLigand(stereo_bond, predicate) => {
                AstRelationalConstraint::StereoBondAnyLigand {
                    stereo_bond: AstStereoBondId(*stereo_bond),
                    predicate: Box::new(predicate.bind(py).borrow().to_rust(py)),
                }
            }
        }
    }
}

#[pymethods]
impl SubPatternAnchor {
    #[new]
    #[pyo3(signature = (*, atoms=Vec::new(), bonds=Vec::new(), dative_bonds=Vec::new(), aromatic_systems=Vec::new(), multicenter_bonds=Vec::new(), noncovalent_bonds=Vec::new(), stereo_atoms=Vec::new(), stereo_bonds=Vec::new()))]
    #[allow(clippy::too_many_arguments)] // one collection per molecule entity family
    fn new(
        atoms: Vec<(u32, u32)>,
        bonds: Vec<(u32, u32)>,
        dative_bonds: Vec<(u32, u32)>,
        aromatic_systems: Vec<(u32, u32)>,
        multicenter_bonds: Vec<(u32, u32)>,
        noncovalent_bonds: Vec<(u32, u32)>,
        stereo_atoms: Vec<(u32, u32)>,
        stereo_bonds: Vec<(u32, u32)>,
    ) -> Self {
        Self {
            atoms,
            bonds,
            dative_bonds,
            aromatic_systems,
            multicenter_bonds,
            noncovalent_bonds,
            stereo_atoms,
            stereo_bonds,
        }
    }

    #[getter]
    fn atoms(&self) -> Vec<(u32, u32)> {
        self.atoms.clone()
    }

    #[getter]
    fn bonds(&self) -> Vec<(u32, u32)> {
        self.bonds.clone()
    }

    #[getter]
    fn dative_bonds(&self) -> Vec<(u32, u32)> {
        self.dative_bonds.clone()
    }

    #[getter]
    fn aromatic_systems(&self) -> Vec<(u32, u32)> {
        self.aromatic_systems.clone()
    }

    #[getter]
    fn multicenter_bonds(&self) -> Vec<(u32, u32)> {
        self.multicenter_bonds.clone()
    }

    #[getter]
    fn noncovalent_bonds(&self) -> Vec<(u32, u32)> {
        self.noncovalent_bonds.clone()
    }

    #[getter]
    fn stereo_atoms(&self) -> Vec<(u32, u32)> {
        self.stereo_atoms.clone()
    }

    #[getter]
    fn stereo_bonds(&self) -> Vec<(u32, u32)> {
        self.stereo_bonds.clone()
    }
}

impl SubPatternAnchor {
    pub(crate) fn from_rust(anchor: &AstSubPatternAnchor) -> Self {
        Self {
            atoms: anchor
                .atoms()
                .iter()
                .map(|(target, pattern)| (target.0, pattern.0))
                .collect(),
            bonds: anchor
                .bonds()
                .iter()
                .map(|(target, pattern)| (target.0, pattern.0))
                .collect(),
            dative_bonds: anchor
                .dative_bonds()
                .iter()
                .map(|(target, pattern)| (target.0, pattern.0))
                .collect(),
            aromatic_systems: anchor
                .aromatic_systems()
                .iter()
                .map(|(target, pattern)| (target.0, pattern.0))
                .collect(),
            multicenter_bonds: anchor
                .multicenter_bonds()
                .iter()
                .map(|(target, pattern)| (target.0, pattern.0))
                .collect(),
            noncovalent_bonds: anchor
                .noncovalent_bonds()
                .iter()
                .map(|(target, pattern)| (target.0, pattern.0))
                .collect(),
            stereo_atoms: anchor
                .stereo_atoms()
                .iter()
                .map(|(target, pattern)| (target.0, pattern.0))
                .collect(),
            stereo_bonds: anchor
                .stereo_bonds()
                .iter()
                .map(|(target, pattern)| (target.0, pattern.0))
                .collect(),
        }
    }

    pub(crate) fn to_rust(&self) -> AstSubPatternAnchor {
        let mut anchor = AstSubPatternAnchor::new();
        for &(target, pattern) in &self.atoms {
            anchor.push_atom(AstAtomId(target), AstAtomId(pattern));
        }
        for &(target, pattern) in &self.bonds {
            anchor.push_bond(AstBondId(target), AstBondId(pattern));
        }
        for &(target, pattern) in &self.dative_bonds {
            anchor.push_dative_bond(AstDativeBondId(target), AstDativeBondId(pattern));
        }
        for &(target, pattern) in &self.aromatic_systems {
            anchor.push_aromatic_system(AstAromaticSystemId(target), AstAromaticSystemId(pattern));
        }
        for &(target, pattern) in &self.multicenter_bonds {
            anchor
                .push_multicenter_bond(AstMulticenterBondId(target), AstMulticenterBondId(pattern));
        }
        for &(target, pattern) in &self.noncovalent_bonds {
            anchor
                .push_noncovalent_bond(AstNoncovalentBondId(target), AstNoncovalentBondId(pattern));
        }
        for &(target, pattern) in &self.stereo_atoms {
            anchor.push_stereo_atom(AstStereoAtomId(target), AstStereoAtomId(pattern));
        }
        for &(target, pattern) in &self.stereo_bonds {
            anchor.push_stereo_bond(AstStereoBondId(target), AstStereoBondId(pattern));
        }
        anchor
    }
}

#[cfg(test)]
mod tests {
    use std::ffi::CString;

    use pyo3::types::PyDict;
    use rstest::rstest;
    use umol_ast::ast::{
        AromaticSystemConstraintAst as AstAromaticSystemConstraintAst,
        AtomConstraintAst as AstAtomConstraintAst, BondConstraintAst as AstBondConstraintAst,
        DativeBondConstraintAst as AstDativeBondConstraintAst, MoleculeAst as AstMoleculeAst,
        MulticenterBondConstraintAst as AstMulticenterBondConstraintAst,
        NoncovalentBondConstraintAst as AstNoncovalentBondConstraintAst,
        StereoAtomConstraintAst as AstStereoAtomConstraintAst,
        StereoBondConstraintAst as AstStereoBondConstraintAst, StereoKind as AstStereoKind,
        Stereogenicity as AstStereogenicity, StereogenicityAst as AstStereogenicityAst,
        UnpairedElectronsAst as AstUnpairedElectronsAst, ValueAst as AstValueAst,
    };

    use super::*;

    #[rstest]
    fn test_sub_pattern_anchor_new() {
        let anchor = SubPatternAnchor::new(
            vec![(1, 2)],
            vec![(3, 4)],
            vec![(5, 6)],
            vec![(7, 8)],
            vec![(9, 10)],
            vec![(11, 12)],
            vec![(13, 14)],
            vec![(15, 16)],
        );

        assert_eq!(anchor.atoms(), vec![(1, 2)]);
        assert_eq!(anchor.bonds(), vec![(3, 4)]);
        assert_eq!(anchor.dative_bonds(), vec![(5, 6)]);
        assert_eq!(anchor.aromatic_systems(), vec![(7, 8)]);
        assert_eq!(anchor.multicenter_bonds(), vec![(9, 10)]);
        assert_eq!(anchor.noncovalent_bonds(), vec![(11, 12)]);
        assert_eq!(anchor.stereo_atoms(), vec![(13, 14)]);
        assert_eq!(anchor.stereo_bonds(), vec![(15, 16)]);
    }

    #[rstest]
    fn test_sub_pattern_anchor_roundtrip() {
        let mut anchor = AstSubPatternAnchor::new();
        anchor.push_atom(AstAtomId(1), AstAtomId(2));
        anchor.push_bond(AstBondId(3), AstBondId(4));
        anchor.push_dative_bond(AstDativeBondId(5), AstDativeBondId(6));
        anchor.push_aromatic_system(AstAromaticSystemId(7), AstAromaticSystemId(8));
        anchor.push_multicenter_bond(AstMulticenterBondId(9), AstMulticenterBondId(10));
        anchor.push_noncovalent_bond(AstNoncovalentBondId(11), AstNoncovalentBondId(12));
        anchor.push_stereo_atom(AstStereoAtomId(13), AstStereoAtomId(14));
        anchor.push_stereo_bond(AstStereoBondId(15), AstStereoBondId(16));

        assert_eq!(SubPatternAnchor::from_rust(&anchor).to_rust(), anchor);
    }

    #[rstest]
    #[case::donors(AstRelationalConstraint::DativeBondDonors {
        bond: AstDativeBondId(1),
        atoms: vec![AstAtomId(2), AstAtomId(3)],
    })]
    #[case::donor(AstRelationalConstraint::DativeBondDonor {
        bond: AstDativeBondId(4),
        atom: AstAtomId(5),
    })]
    #[case::contains_all_donors(AstRelationalConstraint::DativeBondContainsAllDonors {
        bond: AstDativeBondId(6),
        atoms: vec![AstAtomId(7), AstAtomId(8)],
    })]
    #[case::all_donors(AstRelationalConstraint::DativeBondAllDonors {
        bond: AstDativeBondId(9),
        predicate: Box::new(AstAtomConstraintAst::degree(2)),
    })]
    #[case::any_donor(AstRelationalConstraint::DativeBondAnyDonor {
        bond: AstDativeBondId(10),
        predicate: Box::new(AstAtomConstraintAst::valence(3)),
    })]
    #[case::acceptor(AstRelationalConstraint::DativeBondAcceptor {
        bond: AstDativeBondId(11),
        atom: AstAtomId(12),
    })]
    #[case::acceptor_satisfies(AstRelationalConstraint::DativeBondAcceptorSatisfies {
        bond: AstDativeBondId(13),
        predicate: Box::new(AstAtomConstraintAst::total_degree(4)),
    })]
    #[case::parallels(AstRelationalConstraint::DativeBondParallels {
        dative: AstDativeBondId(14),
        parallel: AstBondId(15),
    })]
    #[case::aromatic_atoms(AstRelationalConstraint::AromaticSystemAtoms {
        system: AstAromaticSystemId(16),
        atoms: vec![AstAtomId(17), AstAtomId(18)],
    })]
    #[case::aromatic_contains(AstRelationalConstraint::AromaticSystemContains {
        system: AstAromaticSystemId(19),
        atom: AstAtomId(20),
    })]
    #[case::aromatic_contains_all(AstRelationalConstraint::AromaticSystemContainsAll {
        system: AstAromaticSystemId(21),
        atoms: vec![AstAtomId(22), AstAtomId(23)],
    })]
    #[case::aromatic_all_atoms(AstRelationalConstraint::AromaticSystemAllAtoms {
        system: AstAromaticSystemId(24),
        predicate: Box::new(AstAtomConstraintAst::degree(5)),
    })]
    #[case::aromatic_any_atom(AstRelationalConstraint::AromaticSystemAnyAtom {
        system: AstAromaticSystemId(25),
        predicate: Box::new(AstAtomConstraintAst::valence(6)),
    })]
    #[case::multicenter_atoms(AstRelationalConstraint::MulticenterBondAtoms {
        bond: AstMulticenterBondId(26),
        atoms: vec![AstAtomId(27), AstAtomId(28)],
    })]
    #[case::multicenter_contains(AstRelationalConstraint::MulticenterBondContains {
        bond: AstMulticenterBondId(29),
        atom: AstAtomId(30),
    })]
    #[case::multicenter_contains_all(AstRelationalConstraint::MulticenterBondContainsAll {
        bond: AstMulticenterBondId(31),
        atoms: vec![AstAtomId(32), AstAtomId(33)],
    })]
    #[case::multicenter_all_atoms(AstRelationalConstraint::MulticenterBondAllAtoms {
        bond: AstMulticenterBondId(34),
        predicate: Box::new(AstAtomConstraintAst::degree(7)),
    })]
    #[case::multicenter_any_atom(AstRelationalConstraint::MulticenterBondAnyAtom {
        bond: AstMulticenterBondId(35),
        predicate: Box::new(AstAtomConstraintAst::valence(8)),
    })]
    #[case::noncovalent_ends(AstRelationalConstraint::NoncovalentBondEnds {
        bond: AstNoncovalentBondId(36),
        atoms: [AstAtomId(37), AstAtomId(38)],
    })]
    #[case::noncovalent_contains(AstRelationalConstraint::NoncovalentBondContains {
        bond: AstNoncovalentBondId(39),
        atom: AstAtomId(40),
    })]
    #[case::noncovalent_ends_satisfy(AstRelationalConstraint::NoncovalentBondEndsSatisfy {
        bond: AstNoncovalentBondId(41),
        predicates: [
            Box::new(AstAtomConstraintAst::degree(9)),
            Box::new(AstAtomConstraintAst::valence(10)),
        ],
    })]
    #[case::stereo_atom_site(AstRelationalConstraint::StereoAtomSite {
        stereo_atom: AstStereoAtomId(42),
        atom: AstAtomId(43),
    })]
    #[case::stereo_atom_contains(AstRelationalConstraint::StereoAtomContains {
        stereo_atom: AstStereoAtomId(44),
        atom: AstAtomId(45),
    })]
    #[case::stereo_atom_ligands(AstRelationalConstraint::StereoAtomLigands {
        stereo_atom: AstStereoAtomId(46),
        atoms: vec![AstAtomId(47), AstAtomId(48)],
    })]
    #[case::stereo_atom_all_ligands(AstRelationalConstraint::StereoAtomAllLigands {
        stereo_atom: AstStereoAtomId(49),
        predicate: Box::new(AstAtomConstraintAst::degree(11)),
    })]
    #[case::stereo_atom_any_ligand(AstRelationalConstraint::StereoAtomAnyLigand {
        stereo_atom: AstStereoAtomId(50),
        predicate: Box::new(AstAtomConstraintAst::valence(12)),
    })]
    #[case::stereo_bond_site(AstRelationalConstraint::StereoBondSite {
        stereo_bond: AstStereoBondId(51),
        bond: AstBondId(52),
    })]
    #[case::stereo_bond_contains(AstRelationalConstraint::StereoBondContains {
        stereo_bond: AstStereoBondId(53),
        atom: AstAtomId(54),
    })]
    #[case::stereo_bond_ligands(AstRelationalConstraint::StereoBondLigands {
        stereo_bond: AstStereoBondId(55),
        atoms: vec![AstAtomId(56), AstAtomId(57)],
    })]
    #[case::stereo_bond_all_ligands(AstRelationalConstraint::StereoBondAllLigands {
        stereo_bond: AstStereoBondId(58),
        predicate: Box::new(AstAtomConstraintAst::degree(13)),
    })]
    #[case::stereo_bond_any_ligand(AstRelationalConstraint::StereoBondAnyLigand {
        stereo_bond: AstStereoBondId(59),
        predicate: Box::new(AstAtomConstraintAst::valence(14)),
    })]
    fn test_relational_constraint_roundtrip(#[case] constraint: AstRelationalConstraint) {
        Python::attach(|py| {
            let value = RelationalConstraint::from_rust(py, &constraint).unwrap();
            assert_eq!(value.to_rust(py), constraint);
        });
    }

    #[rstest]
    #[case::charge_sum_whole(AstMoleculeConstraint::ChargeSum {
        atoms: None,
        sum: AstValueAst::Lit(1),
    })]
    #[case::charge_sum_empty_subset(AstMoleculeConstraint::ChargeSum {
        atoms: Some(Vec::new()),
        sum: AstValueAst::Lit(2),
    })]
    #[case::unpaired_electron_coupling(AstMoleculeConstraint::UnpairedElectronCoupling {
        atoms: Some(vec![AstAtomId(3), AstAtomId(4)]),
        unpaired_electrons: AstUnpairedElectronsAst::from((1, 2)),
    })]
    #[case::bond_order_sum(AstMoleculeConstraint::BondOrderSum {
        bonds: Some(vec![AstBondId(5), AstBondId(6)]),
        sum: AstValueAst::Lit(3),
    })]
    #[case::connected(AstMoleculeConstraint::Connected {
        atoms: None,
    })]
    #[case::sub_pattern(AstMoleculeConstraint::SubPattern {
        anchor: {
            let mut anchor = AstSubPatternAnchor::new();
            anchor.push_atom(AstAtomId(7), AstAtomId(0));
            anchor
        },
        pattern: Box::new(AstMoleculeAst::new()),
    })]
    fn test_molecule_constraint_roundtrip(#[case] constraint: AstMoleculeConstraint) {
        Python::attach(|py| {
            let value = MoleculeConstraint::from_rust(py, &constraint).unwrap();
            assert_eq!(value.to_rust(py), constraint);
        });
    }

    #[rstest]
    #[case::atom(AstConstraint::Atom(AstAtomId(1), AstAtomConstraintAst::degree(2),))]
    #[case::bond(AstConstraint::Bond(AstBondId(3), AstBondConstraintAst::aromatic(true),))]
    #[case::dative_bond(AstConstraint::DativeBond(
        AstDativeBondId(4),
        AstDativeBondConstraintAst::aromatic(false),
    ))]
    #[case::aromatic_system(AstConstraint::AromaticSystem(
        AstAromaticSystemId(5),
        AstAromaticSystemConstraintAst::electron_count(6),
    ))]
    #[case::multicenter_bond(AstConstraint::MulticenterBond(
        AstMulticenterBondId(7),
        AstMulticenterBondConstraintAst::electron_count(8),
    ))]
    #[case::noncovalent_bond(AstConstraint::NoncovalentBond(
        AstNoncovalentBondId(9),
        AstNoncovalentBondConstraintAst::intramolecular(true),
    ))]
    #[case::stereo_atom(AstConstraint::StereoAtom(
        AstStereoAtomId(10),
        AstStereoKind::Tetrahedral,
        AstStereoAtomConstraintAst::Stereogenicity(AstStereogenicityAst::Lit(
            AstStereogenicity::Stereogenic,
        )),
    ))]
    #[case::stereo_bond(AstConstraint::StereoBond(
        AstStereoBondId(11),
        AstStereoKind::CisTrans,
        AstStereoBondConstraintAst::Stereogenicity(AstStereogenicityAst::Lit(
            AstStereogenicity::Prochiral,
        )),
    ))]
    #[case::relational(AstConstraint::Relational(
        AstRelationalConstraint::DativeBondDonor {
            bond: AstDativeBondId(12),
            atom: AstAtomId(13),
        },
    ))]
    #[case::molecule(AstConstraint::Molecule(
        AstMoleculeConstraint::Connected {
            atoms: Some(vec![AstAtomId(14), AstAtomId(15)]),
        },
    ))]
    #[case::and(AstConstraint::And(Vec::new()))]
    #[case::or(AstConstraint::Or(Vec::new()))]
    #[case::not(AstConstraint::Not(Box::new(AstConstraint::Atom(
        AstAtomId(16),
        AstAtomConstraintAst::degree(3),
    ))))]
    fn test_constraint_roundtrip(#[case] constraint: AstConstraint) {
        Python::attach(|py| {
            let value = Constraint::from_rust(py, &constraint).unwrap();
            assert_eq!(value.to_rust(py), constraint);
        });
    }

    #[rstest]
    fn test_constraint_roundtrip_recursive() {
        let constraint = AstConstraint::And(vec![
            AstConstraint::Atom(AstAtomId(17), AstAtomConstraintAst::valence(4)),
            AstConstraint::Or(vec![
                AstConstraint::Relational(AstRelationalConstraint::DativeBondDonor {
                    bond: AstDativeBondId(18),
                    atom: AstAtomId(19),
                }),
                AstConstraint::Not(Box::new(AstConstraint::Molecule(
                    AstMoleculeConstraint::Connected {
                        atoms: Some(vec![AstAtomId(20), AstAtomId(21)]),
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
                "Constraint.And([Constraint.Atom(17, AtomConstraintAst.Valence(ValueAst.Lit(4))), Constraint.Or([Constraint.Relational(RelationalConstraint.DativeBondDonor(18, 19)), Constraint.Not(Constraint.Molecule(MoleculeConstraint.Connected([20, 21])))])])"
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
        let first = AstConstraint::Atom(AstAtomId(1), AstAtomConstraintAst::degree(2));
        let second = AstConstraint::Molecule(AstMoleculeConstraint::Connected { atoms: None });
        let mut constraints = AstConstraints::from(vec![first.clone(), second.clone()]);

        Python::attach(|py| {
            let mut iter = constraint_iter(py, &constraints).unwrap();
            constraints.push(AstConstraint::Or(Vec::new()));

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
    fn test_constraints_arg_to_rust_container() {
        let expected = AstConstraints::from(vec![
            AstConstraint::And(Vec::new()),
            AstConstraint::And(Vec::new()),
        ]);

        Python::attach(|py| {
            let container = Py::new(py, Constraints::from_inner(expected.clone())).unwrap();
            let arg = ConstraintsArg::Container(container);

            assert_eq!(arg.to_rust(py).unwrap(), expected);
        });
    }

    #[rstest]
    fn test_constraints_arg_to_rust_view() {
        let expected = AstConstraints::from(vec![
            AstConstraint::Or(Vec::new()),
            AstConstraint::Or(Vec::new()),
        ]);
        let mut molecule = AstMoleculeAst::new();
        *molecule.constraints_mut() = expected.clone();

        Python::attach(|py| {
            let owner = Py::new(py, MoleculeAst::from_inner(molecule)).unwrap();
            let view = Py::new(py, ConstraintsView::new(owner)).unwrap();
            let arg = ConstraintsArg::View(view);

            assert_eq!(arg.to_rust(py).unwrap(), expected);
        });
    }

    #[rstest]
    #[case::empty(Vec::new())]
    #[case::populated(vec![
        AstConstraint::Atom(AstAtomId(1), AstAtomConstraintAst::degree(2)),
        AstConstraint::Molecule(AstMoleculeConstraint::Connected { atoms: None }),
    ])]
    fn test_constraints_new(#[case] entries: Vec<AstConstraint>) {
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
        AstConstraints::from(vec![AstConstraint::And(Vec::new())]),
        AstConstraints::from(vec![AstConstraint::And(Vec::new())]),
        true,
    )]
    #[case::different(
        AstConstraints::from(vec![AstConstraint::And(Vec::new())]),
        AstConstraints::from(vec![AstConstraint::Or(Vec::new())]),
        false,
    )]
    fn test_constraints_eq(
        #[case] left: AstConstraints,
        #[case] right: AstConstraints,
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
            let constraints = Constraints::from_inner(AstConstraints::from(vec![
                AstConstraint::Atom(AstAtomId(1), AstAtomConstraintAst::degree(2)),
                AstConstraint::Or(Vec::new()),
            ]));

            assert_eq!(
                constraints.__repr__(py).unwrap(),
                "Constraints([Constraint.Atom(1, AtomConstraintAst.Degree(ValueAst.Lit(2))), Constraint.Or([])])"
            );
        });
    }

    #[rstest]
    fn test_constraints_append() {
        let constraint = AstConstraint::Molecule(AstMoleculeConstraint::Connected { atoms: None });
        Python::attach(|py| {
            let mut constraints =
                Constraints::from_inner(AstConstraints::from(vec![constraint.clone()]));
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
            Constraints::from_inner(AstConstraints::from(vec![AstConstraint::And(Vec::new())]));

        constraints.clear();

        assert_eq!(constraints.inner(), &AstConstraints::new());
    }

    #[rstest]
    fn test_constraints_update() {
        let initial = AstConstraint::And(Vec::new());
        let from_container = AstConstraint::Or(Vec::new());
        let from_view = AstConstraint::Not(Box::new(AstConstraint::And(Vec::new())));
        let from_entries =
            AstConstraint::Molecule(AstMoleculeConstraint::Connected { atoms: None });

        Python::attach(|py| {
            let target = Py::new(
                py,
                Constraints::from_inner(AstConstraints::from(vec![initial.clone()])),
            )
            .unwrap();
            let container = Py::new(
                py,
                Constraints::from_inner(AstConstraints::from(vec![from_container.clone()])),
            )
            .unwrap();
            let mut molecule = AstMoleculeAst::new();
            molecule.constraints_mut().push(from_view.clone());
            let view = Py::new(
                py,
                ConstraintsView::new(Py::new(py, MoleculeAst::from_inner(molecule)).unwrap()),
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
        let entry = AstConstraint::And(Vec::new());

        Python::attach(|py| {
            let target = Py::new(
                py,
                Constraints::from_inner(AstConstraints::from(vec![entry.clone()])),
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
    #[case::empty(AstConstraints::new(), 0)]
    #[case::populated(AstConstraints::from(vec![
        AstConstraint::And(Vec::new()),
        AstConstraint::Or(Vec::new()),
    ]), 2)]
    fn test_constraints_len(#[case] inner: AstConstraints, #[case] expected: usize) {
        assert_eq!(Constraints::from_inner(inner).__len__(), expected);
    }

    #[rstest]
    #[case::positive(0, AstConstraint::Atom(AstAtomId(1), AstAtomConstraintAst::degree(2),))]
    #[case::negative(-1, AstConstraint::Molecule(
        AstMoleculeConstraint::Connected { atoms: None },
    ))]
    fn test_constraints_getitem(#[case] index: isize, #[case] expected: AstConstraint) {
        Python::attach(|py| {
            let constraints = Constraints::from_inner(AstConstraints::from(vec![
                AstConstraint::Atom(AstAtomId(1), AstAtomConstraintAst::degree(2)),
                AstConstraint::Molecule(AstMoleculeConstraint::Connected { atoms: None }),
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
                Constraints::from_inner(AstConstraints::from(vec![AstConstraint::And(Vec::new())]));

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
        let first = AstConstraint::And(Vec::new());
        let second = AstConstraint::Or(Vec::new());
        let mut constraints =
            Constraints::from_inner(AstConstraints::from(vec![first.clone(), second.clone()]));

        Python::attach(|py| {
            let mut iter = constraints.__iter__(py).unwrap();
            constraints
                .inner_mut()
                .push(AstConstraint::Not(Box::new(AstConstraint::And(Vec::new()))));

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
        let entries = vec![AstConstraint::Or(Vec::new())];

        assert_eq!(
            Constraints::from_inner(AstConstraints::from(entries.clone()))
                .inner()
                .as_slice(),
            entries.as_slice()
        );
    }

    #[rstest]
    fn test_constraints_view_repr() {
        let mut molecule = AstMoleculeAst::new();
        molecule
            .constraints_mut()
            .push(AstConstraint::And(Vec::new()));
        molecule
            .constraints_mut()
            .push(AstConstraint::Or(Vec::new()));

        Python::attach(|py| {
            let owner = Py::new(py, MoleculeAst::from_inner(molecule)).unwrap();
            let view = ConstraintsView::new(owner);

            assert_eq!(view.__repr__(py).unwrap(), "ConstraintsView(2 entries)");
        });
    }

    #[rstest]
    fn test_constraints_view_append() {
        let constraint = AstConstraint::Molecule(AstMoleculeConstraint::Connected { atoms: None });
        let mut molecule = AstMoleculeAst::new();
        molecule.constraints_mut().push(constraint.clone());

        Python::attach(|py| {
            let owner = Py::new(py, MoleculeAst::from_inner(molecule)).unwrap();
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
        let mut molecule = AstMoleculeAst::new();
        molecule
            .constraints_mut()
            .push(AstConstraint::And(Vec::new()));

        Python::attach(|py| {
            let owner = Py::new(py, MoleculeAst::from_inner(molecule)).unwrap();
            let view = ConstraintsView::new(owner.clone_ref(py));

            view.clear(py);

            assert_eq!(
                owner.bind(py).borrow().inner().constraints(),
                &AstConstraints::new()
            );
        });
    }

    #[rstest]
    fn test_constraints_view_update() {
        let initial = AstConstraint::And(Vec::new());
        let from_container = AstConstraint::Or(Vec::new());
        let from_view = AstConstraint::Not(Box::new(AstConstraint::And(Vec::new())));
        let from_entries =
            AstConstraint::Molecule(AstMoleculeConstraint::Connected { atoms: None });
        let mut target_molecule = AstMoleculeAst::new();
        target_molecule.constraints_mut().push(initial.clone());

        Python::attach(|py| {
            let target_owner = Py::new(py, MoleculeAst::from_inner(target_molecule)).unwrap();
            let target = ConstraintsView::new(target_owner.clone_ref(py));
            let container = Py::new(
                py,
                Constraints::from_inner(AstConstraints::from(vec![from_container.clone()])),
            )
            .unwrap();
            let mut source_molecule = AstMoleculeAst::new();
            source_molecule.constraints_mut().push(from_view.clone());
            let source_view = Py::new(
                py,
                ConstraintsView::new(
                    Py::new(py, MoleculeAst::from_inner(source_molecule)).unwrap(),
                ),
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
        let entry = AstConstraint::Or(Vec::new());
        let mut molecule = AstMoleculeAst::new();
        molecule.constraints_mut().push(entry.clone());

        Python::attach(|py| {
            let owner = Py::new(py, MoleculeAst::from_inner(molecule)).unwrap();
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
            let owner = Py::new(py, MoleculeAst::from_inner(AstMoleculeAst::new())).unwrap();
            let view = ConstraintsView::new(owner.clone_ref(py));
            assert_eq!(view.__len__(py).unwrap(), 0);

            owner
                .borrow_mut(py)
                .inner_mut()
                .constraints_mut()
                .push(AstConstraint::And(Vec::new()));

            assert_eq!(view.__len__(py).unwrap(), 1);
        });
    }

    #[rstest]
    #[case::positive(0, AstConstraint::Atom(AstAtomId(1), AstAtomConstraintAst::degree(2),))]
    #[case::negative(-1, AstConstraint::Molecule(
        AstMoleculeConstraint::Connected { atoms: None },
    ))]
    fn test_constraints_view_getitem(#[case] index: isize, #[case] expected: AstConstraint) {
        let mut molecule = AstMoleculeAst::new();
        molecule.constraints_mut().push(AstConstraint::Atom(
            AstAtomId(1),
            AstAtomConstraintAst::degree(2),
        ));
        molecule.constraints_mut().push(AstConstraint::Molecule(
            AstMoleculeConstraint::Connected { atoms: None },
        ));

        Python::attach(|py| {
            let owner = Py::new(py, MoleculeAst::from_inner(molecule)).unwrap();
            let view = ConstraintsView::new(owner);

            assert_eq!(view.__getitem__(py, index).unwrap().to_rust(py), expected);
        });
    }

    #[rstest]
    #[case::positive(1)]
    #[case::negative(-2)]
    fn test_constraints_view_getitem_error(#[case] index: isize) {
        let mut molecule = AstMoleculeAst::new();
        molecule
            .constraints_mut()
            .push(AstConstraint::And(Vec::new()));

        Python::attach(|py| {
            let owner = Py::new(py, MoleculeAst::from_inner(molecule)).unwrap();
            let view = ConstraintsView::new(owner);

            assert_eq!(
                view.__getitem__(py, index).err().unwrap().to_string(),
                "IndexError: constraint index out of range"
            );
        });
    }

    #[rstest]
    fn test_constraints_view_iter() {
        let first = AstConstraint::And(Vec::new());
        let second = AstConstraint::Or(Vec::new());
        let mut molecule = AstMoleculeAst::new();
        molecule.constraints_mut().push(first.clone());
        molecule.constraints_mut().push(second.clone());

        Python::attach(|py| {
            let owner = Py::new(py, MoleculeAst::from_inner(molecule)).unwrap();
            let view = ConstraintsView::new(owner.clone_ref(py));
            let mut iter = view.__iter__(py).unwrap();
            owner
                .borrow_mut(py)
                .inner_mut()
                .constraints_mut()
                .push(AstConstraint::Not(Box::new(AstConstraint::And(Vec::new()))));

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
