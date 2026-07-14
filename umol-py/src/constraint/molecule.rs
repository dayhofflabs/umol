//! Molecule-level constraint payloads mirroring `umol_ast::ast::constraint`.

use pyo3::prelude::*;
use umol_ast::ast::{
    AromaticSystemId as AstAromaticSystemId, AtomId as AstAtomId, BondId as AstBondId,
    Constraint as AstConstraint, DativeBondId as AstDativeBondId,
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
use crate::spin::SpinStateAst;
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
    SpinSum(Option<Vec<u32>>, Py<SpinStateAst>),
    BondOrderSum(Option<Vec<u32>>, Py<ValueAst>),
    Connected(Option<Vec<u32>>),
    SubPattern(Py<SubPatternAnchor>, Py<MoleculeAst>),
}

/// A recursive molecule-constraint tree. This slice contains the ordinary and
/// stereo entity leaves; aggregate and combinator variants follow in S1b.
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
}

#[pymethods]
impl Constraint {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
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
        };
        variant_repr(slf.bind(py).as_any(), "Constraint", variant, arity)
    }
}

#[allow(dead_code, reason = "consumed by the constraint container later in S1")]
impl Constraint {
    /// Convert an entity leaf. Other variants return `None` until the remaining
    /// S1b slices extend this mirror.
    pub(crate) fn from_ast(py: Python<'_>, constraint: &AstConstraint) -> PyResult<Option<Self>> {
        Ok(Some(match constraint {
            AstConstraint::Atom(id, child) => Self::Atom(
                id.0,
                into_py_variant(py, AtomConstraintAst::from_ast(py, child)?)?,
            ),
            AstConstraint::Bond(id, child) => Self::Bond(
                id.0,
                into_py_variant(py, BondConstraintAst::from_ast(py, child)?)?,
            ),
            AstConstraint::DativeBond(id, child) => Self::DativeBond(
                id.0,
                into_py_variant(py, DativeBondConstraintAst::from_ast(py, child)?)?,
            ),
            AstConstraint::AromaticSystem(id, child) => Self::AromaticSystem(
                id.0,
                into_py_variant(py, AromaticSystemConstraintAst::from_ast(py, child)?)?,
            ),
            AstConstraint::MulticenterBond(id, child) => Self::MulticenterBond(
                id.0,
                into_py_variant(py, MulticenterBondConstraintAst::from_ast(py, child)?)?,
            ),
            AstConstraint::NoncovalentBond(id, child) => Self::NoncovalentBond(
                id.0,
                into_py_variant(py, NoncovalentBondConstraintAst::from_ast(py, child)?)?,
            ),
            AstConstraint::StereoAtom(id, kind, child) => Self::StereoAtom(
                id.0,
                StereoKind::from_ast(*kind),
                into_py_variant(py, StereoAtomConstraintAst::from_ast(py, child)?)?,
            ),
            AstConstraint::StereoBond(id, kind, child) => Self::StereoBond(
                id.0,
                StereoKind::from_ast(*kind),
                into_py_variant(py, StereoBondConstraintAst::from_ast(py, child)?)?,
            ),
            _ => return Ok(None),
        }))
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstConstraint {
        match self {
            Self::Atom(id, child) => {
                AstConstraint::Atom(AstAtomId(*id), child.bind(py).borrow().to_ast(py))
            }
            Self::Bond(id, child) => {
                AstConstraint::Bond(AstBondId(*id), child.bind(py).borrow().to_ast(py))
            }
            Self::DativeBond(id, child) => {
                AstConstraint::DativeBond(AstDativeBondId(*id), child.bind(py).borrow().to_ast(py))
            }
            Self::AromaticSystem(id, child) => AstConstraint::AromaticSystem(
                AstAromaticSystemId(*id),
                child.bind(py).borrow().to_ast(py),
            ),
            Self::MulticenterBond(id, child) => AstConstraint::MulticenterBond(
                AstMulticenterBondId(*id),
                child.bind(py).borrow().to_ast(py),
            ),
            Self::NoncovalentBond(id, child) => AstConstraint::NoncovalentBond(
                AstNoncovalentBondId(*id),
                child.bind(py).borrow().to_ast(py),
            ),
            Self::StereoAtom(id, kind, child) => AstConstraint::StereoAtom(
                AstStereoAtomId(*id),
                kind.to_ast(),
                child.bind(py).borrow().to_ast(py),
            ),
            Self::StereoBond(id, kind, child) => AstConstraint::StereoBond(
                AstStereoBondId(*id),
                kind.to_ast(),
                child.bind(py).borrow().to_ast(py),
            ),
        }
    }
}

#[pymethods]
impl MoleculeConstraint {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            Self::ChargeSum(_, _) => ("ChargeSum", 2),
            Self::SpinSum(_, _) => ("SpinSum", 2),
            Self::BondOrderSum(_, _) => ("BondOrderSum", 2),
            Self::Connected(_) => ("Connected", 1),
            Self::SubPattern(_, _) => ("SubPattern", 2),
        };
        variant_repr(slf.bind(py).as_any(), "MoleculeConstraint", variant, arity)
    }
}

#[allow(
    dead_code,
    reason = "consumed by the recursive constraint mirror later in S1a"
)]
impl MoleculeConstraint {
    pub(crate) fn from_ast(py: Python<'_>, constraint: &AstMoleculeConstraint) -> PyResult<Self> {
        Ok(match constraint {
            AstMoleculeConstraint::ChargeSum { atoms, sum } => Self::ChargeSum(
                atoms
                    .as_ref()
                    .map(|atoms| atoms.iter().map(|atom| atom.0).collect()),
                into_py_variant(py, ValueAst::from_ast(py, sum)?)?,
            ),
            AstMoleculeConstraint::SpinSum { atoms, spin } => Self::SpinSum(
                atoms
                    .as_ref()
                    .map(|atoms| atoms.iter().map(|atom| atom.0).collect()),
                Py::new(py, SpinStateAst::from_ast(py, spin)?)?,
            ),
            AstMoleculeConstraint::BondOrderSum { bonds, sum } => Self::BondOrderSum(
                bonds
                    .as_ref()
                    .map(|bonds| bonds.iter().map(|bond| bond.0).collect()),
                into_py_variant(py, ValueAst::from_ast(py, sum)?)?,
            ),
            AstMoleculeConstraint::Connected { atoms } => Self::Connected(
                atoms
                    .as_ref()
                    .map(|atoms| atoms.iter().map(|atom| atom.0).collect()),
            ),
            AstMoleculeConstraint::SubPattern { anchor, pattern } => Self::SubPattern(
                Py::new(py, SubPatternAnchor::from_ast(anchor))?,
                Py::new(py, MoleculeAst::from_inner((**pattern).clone()))?,
            ),
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstMoleculeConstraint {
        match self {
            Self::ChargeSum(atoms, sum) => AstMoleculeConstraint::ChargeSum {
                atoms: atoms
                    .as_ref()
                    .map(|atoms| atoms.iter().copied().map(AstAtomId).collect()),
                sum: sum.bind(py).borrow().to_ast(py),
            },
            Self::SpinSum(atoms, spin) => AstMoleculeConstraint::SpinSum {
                atoms: atoms
                    .as_ref()
                    .map(|atoms| atoms.iter().copied().map(AstAtomId).collect()),
                spin: spin.bind(py).borrow().to_ast(py),
            },
            Self::BondOrderSum(bonds, sum) => AstMoleculeConstraint::BondOrderSum {
                bonds: bonds
                    .as_ref()
                    .map(|bonds| bonds.iter().copied().map(AstBondId).collect()),
                sum: sum.bind(py).borrow().to_ast(py),
            },
            Self::Connected(atoms) => AstMoleculeConstraint::Connected {
                atoms: atoms
                    .as_ref()
                    .map(|atoms| atoms.iter().copied().map(AstAtomId).collect()),
            },
            Self::SubPattern(anchor, pattern) => AstMoleculeConstraint::SubPattern {
                anchor: anchor.bind(py).borrow().to_ast(),
                pattern: Box::new(pattern.bind(py).borrow().inner().clone()),
            },
        }
    }
}

#[pymethods]
impl RelationalConstraint {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
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

#[allow(
    dead_code,
    reason = "consumed by the recursive constraint mirror later in S1a"
)]
impl RelationalConstraint {
    /// Convert any relational constraint into its Python mirror.
    pub(crate) fn from_ast(py: Python<'_>, constraint: &AstRelationalConstraint) -> PyResult<Self> {
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
                    into_py_variant(py, AtomConstraintAst::from_ast(py, predicate)?)?,
                )
            }
            AstRelationalConstraint::DativeBondAnyDonor { bond, predicate } => {
                Self::DativeBondAnyDonor(
                    bond.0,
                    into_py_variant(py, AtomConstraintAst::from_ast(py, predicate)?)?,
                )
            }
            AstRelationalConstraint::DativeBondAcceptor { bond, atom } => {
                Self::DativeBondAcceptor(bond.0, atom.0)
            }
            AstRelationalConstraint::DativeBondAcceptorSatisfies { bond, predicate } => {
                Self::DativeBondAcceptorSatisfies(
                    bond.0,
                    into_py_variant(py, AtomConstraintAst::from_ast(py, predicate)?)?,
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
                    into_py_variant(py, AtomConstraintAst::from_ast(py, predicate)?)?,
                )
            }
            AstRelationalConstraint::AromaticSystemAnyAtom { system, predicate } => {
                Self::AromaticSystemAnyAtom(
                    system.0,
                    into_py_variant(py, AtomConstraintAst::from_ast(py, predicate)?)?,
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
                    into_py_variant(py, AtomConstraintAst::from_ast(py, predicate)?)?,
                )
            }
            AstRelationalConstraint::MulticenterBondAnyAtom { bond, predicate } => {
                Self::MulticenterBondAnyAtom(
                    bond.0,
                    into_py_variant(py, AtomConstraintAst::from_ast(py, predicate)?)?,
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
                        into_py_variant(py, AtomConstraintAst::from_ast(py, &predicates[0])?)?,
                        into_py_variant(py, AtomConstraintAst::from_ast(py, &predicates[1])?)?,
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
                into_py_variant(py, AtomConstraintAst::from_ast(py, predicate)?)?,
            ),
            AstRelationalConstraint::StereoAtomAnyLigand {
                stereo_atom,
                predicate,
            } => Self::StereoAtomAnyLigand(
                stereo_atom.0,
                into_py_variant(py, AtomConstraintAst::from_ast(py, predicate)?)?,
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
                into_py_variant(py, AtomConstraintAst::from_ast(py, predicate)?)?,
            ),
            AstRelationalConstraint::StereoBondAnyLigand {
                stereo_bond,
                predicate,
            } => Self::StereoBondAnyLigand(
                stereo_bond.0,
                into_py_variant(py, AtomConstraintAst::from_ast(py, predicate)?)?,
            ),
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstRelationalConstraint {
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
                    predicate: Box::new(predicate.bind(py).borrow().to_ast(py)),
                }
            }
            Self::DativeBondAnyDonor(bond, predicate) => {
                AstRelationalConstraint::DativeBondAnyDonor {
                    bond: AstDativeBondId(*bond),
                    predicate: Box::new(predicate.bind(py).borrow().to_ast(py)),
                }
            }
            Self::DativeBondAcceptor(bond, atom) => AstRelationalConstraint::DativeBondAcceptor {
                bond: AstDativeBondId(*bond),
                atom: AstAtomId(*atom),
            },
            Self::DativeBondAcceptorSatisfies(bond, predicate) => {
                AstRelationalConstraint::DativeBondAcceptorSatisfies {
                    bond: AstDativeBondId(*bond),
                    predicate: Box::new(predicate.bind(py).borrow().to_ast(py)),
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
                    predicate: Box::new(predicate.bind(py).borrow().to_ast(py)),
                }
            }
            Self::AromaticSystemAnyAtom(system, predicate) => {
                AstRelationalConstraint::AromaticSystemAnyAtom {
                    system: AstAromaticSystemId(*system),
                    predicate: Box::new(predicate.bind(py).borrow().to_ast(py)),
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
                    predicate: Box::new(predicate.bind(py).borrow().to_ast(py)),
                }
            }
            Self::MulticenterBondAnyAtom(bond, predicate) => {
                AstRelationalConstraint::MulticenterBondAnyAtom {
                    bond: AstMulticenterBondId(*bond),
                    predicate: Box::new(predicate.bind(py).borrow().to_ast(py)),
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
                        Box::new(predicates[0].bind(py).borrow().to_ast(py)),
                        Box::new(predicates[1].bind(py).borrow().to_ast(py)),
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
                    predicate: Box::new(predicate.bind(py).borrow().to_ast(py)),
                }
            }
            Self::StereoAtomAnyLigand(stereo_atom, predicate) => {
                AstRelationalConstraint::StereoAtomAnyLigand {
                    stereo_atom: AstStereoAtomId(*stereo_atom),
                    predicate: Box::new(predicate.bind(py).borrow().to_ast(py)),
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
                    predicate: Box::new(predicate.bind(py).borrow().to_ast(py)),
                }
            }
            Self::StereoBondAnyLigand(stereo_bond, predicate) => {
                AstRelationalConstraint::StereoBondAnyLigand {
                    stereo_bond: AstStereoBondId(*stereo_bond),
                    predicate: Box::new(predicate.bind(py).borrow().to_ast(py)),
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

#[allow(
    dead_code,
    reason = "consumed by the molecule constraint mirror later in S1a"
)]
impl SubPatternAnchor {
    pub(crate) fn from_ast(anchor: &AstSubPatternAnchor) -> Self {
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

    pub(crate) fn to_ast(&self) -> AstSubPatternAnchor {
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
    use rstest::rstest;
    use umol_ast::ast::{
        AromaticSystemConstraintAst as AstAromaticSystemConstraintAst,
        AtomConstraintAst as AstAtomConstraintAst, BondConstraintAst as AstBondConstraintAst,
        DativeBondConstraintAst as AstDativeBondConstraintAst, MoleculeAst as AstMoleculeAst,
        MulticenterBondConstraintAst as AstMulticenterBondConstraintAst,
        NoncovalentBondConstraintAst as AstNoncovalentBondConstraintAst,
        SpinStateAst as AstSpinStateAst, StereoAtomConstraintAst as AstStereoAtomConstraintAst,
        StereoBondConstraintAst as AstStereoBondConstraintAst, StereoKind as AstStereoKind,
        Stereogenicity as AstStereogenicity, StereogenicityAst as AstStereogenicityAst,
        ValueAst as AstValueAst,
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

        assert_eq!(SubPatternAnchor::from_ast(&anchor).to_ast(), anchor);
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
            let mirror = RelationalConstraint::from_ast(py, &constraint).unwrap();
            assert_eq!(mirror.to_ast(py), constraint);
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
    #[case::spin_sum(AstMoleculeConstraint::SpinSum {
        atoms: Some(vec![AstAtomId(3), AstAtomId(4)]),
        spin: AstSpinStateAst::from((1, 2)),
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
            let mirror = MoleculeConstraint::from_ast(py, &constraint).unwrap();
            assert_eq!(mirror.to_ast(py), constraint);
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
    fn test_constraint_roundtrip(#[case] constraint: AstConstraint) {
        Python::attach(|py| {
            let mirror = Constraint::from_ast(py, &constraint)
                .unwrap()
                .expect("ordinary entity leaf");
            assert_eq!(mirror.to_ast(py), constraint);
        });
    }
}
