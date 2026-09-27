//! Molecule-level constraint payloads matching `umol_graph_ir::ir::constraint`.

use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use pyo3::types::{PyList, PySequence};
use pyo3::PyClassInitializer;
use umol_graph_ir::ir::{
    AromaticSystemId as GraphIrAromaticSystemId, AtomId as GraphIrAtomId, BondId as GraphIrBondId,
    Constraint as GraphIrConstraint, Constraints as GraphIrConstraints,
    DativeBondId as GraphIrDativeBondId, MoleculeConstraint as GraphIrMoleculeConstraint,
    MulticenterBondId as GraphIrMulticenterBondId, NoncovalentBondId as GraphIrNoncovalentBondId,
    Normalize, RelationalConstraint as GraphIrRelationalConstraint,
    StereoAtomId as GraphIrStereoAtomId, StereoBondId as GraphIrStereoBondId,
};

use super::aromatic::AromaticSystemConstraintForm;
use super::atom::AtomConstraintForm;
use super::bond::BondConstraintForm;
use super::dative::DativeBondConstraintForm;
use super::multicenter::MulticenterBondConstraintForm;
use super::noncovalent::NoncovalentBondConstraintForm;
use super::stereo::{StereoAtomConstraintForm, StereoBondConstraintForm};
use crate::convert::{into_py_variant, variant_repr};
use crate::error::{contradiction_error, InvalidatedViewError};
use crate::lattice::impl_py_normalize;
use crate::molecule::Molecule;
use crate::num::NumForm;
use crate::spin::UnpairedElectronsForm;
use crate::stereo::StereoKind;

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

/// A constraint value or read-only access to a stored constraint.
#[pyclass(frozen, subclass)]
pub struct Constraint {
    storage: ConstraintStorage,
}

enum ConstraintStorage {
    Owned(GraphIrConstraint),
    Molecule {
        owner: Py<Molecule>,
        counter: u64,
        position: usize,
        path: Box<[usize]>,
    },
}

#[pymethods]
impl Constraint {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> PyResult<bool> {
        self.read(py, |lhs| other.read(py, |rhs| Ok(lhs == rhs)))
    }

    fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        self.read(py, |value| {
            Ok(match value {
                GraphIrConstraint::Atom(..)
                | GraphIrConstraint::Bond(..)
                | GraphIrConstraint::DativeBond(..)
                | GraphIrConstraint::AromaticSystem(..)
                | GraphIrConstraint::MulticenterBond(..)
                | GraphIrConstraint::NoncovalentBond(..) => 2,
                GraphIrConstraint::StereoAtom(..) | GraphIrConstraint::StereoBond(..) => 3,
                GraphIrConstraint::Relational(..)
                | GraphIrConstraint::Molecule(..)
                | GraphIrConstraint::And(..)
                | GraphIrConstraint::Or(..)
                | GraphIrConstraint::Not(..) => 1,
            })
        })
    }

    fn __getitem__(slf: Py<Self>, py: Python<'_>, index: usize) -> PyResult<Py<PyAny>> {
        if index >= slf.try_borrow(py)?.__len__(py)? {
            return Err(PyIndexError::new_err("tuple index out of range"));
        }
        Ok(slf.bind(py).getattr(format!("_{index}"))?.unbind())
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        self.read(py, |value| {
            let (name, fields) = constraint_fields(py, value)?;
            Ok(format!("Constraint.{name}({})", fields.join(", ")))
        })
    }

    /// Copy the selected constraint into an independent value.
    fn copy(&self, py: Python<'_>) -> PyResult<Py<Self>> {
        self.read(py, |value| Self::from_rust(py, value))
    }

    fn normalize(&self, py: Python<'_>) -> PyResult<Py<Self>> {
        let value = self.to_rust(py)?.normalize().map_err(contradiction_error)?;
        Self::from_storage(py, ConstraintStorage::Owned(value))
    }

    fn normalized_eq(&self, py: Python<'_>, other: &Self) -> PyResult<bool> {
        self.read(py, |lhs| other.read(py, |rhs| Ok(lhs.normalized_eq(rhs))))
    }
}

impl Constraint {
    fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&GraphIrConstraint) -> PyResult<R>,
    ) -> PyResult<R> {
        match &self.storage {
            ConstraintStorage::Owned(value) => f(value),
            ConstraintStorage::Molecule {
                owner,
                counter,
                position,
                path,
            } => {
                let owner = owner.try_borrow(py)?;
                owner.check_access(*counter, "Constraint")?;
                let mut value = owner
                    .to_rust()
                    .constraints()
                    .as_slice()
                    .get(*position)
                    .ok_or_else(|| {
                        InvalidatedViewError::new_err("Constraint is no longer available")
                    })?;
                for index in path {
                    value = match value {
                        GraphIrConstraint::And(children) | GraphIrConstraint::Or(children) => {
                            children.get(*index)
                        }
                        GraphIrConstraint::Not(child) if *index == 0 => Some(child.as_ref()),
                        _ => None,
                    }
                    .ok_or_else(|| {
                        InvalidatedViewError::new_err("Constraint is no longer available")
                    })?;
                }
                f(value)
            }
        }
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> PyResult<GraphIrConstraint> {
        self.read(py, |value| Ok(value.clone()))
    }

    pub(crate) fn from_rust(py: Python<'_>, value: &GraphIrConstraint) -> PyResult<Py<Self>> {
        Self::from_storage(py, ConstraintStorage::Owned(value.clone()))
    }

    fn child(&self, py: Python<'_>, index: usize) -> PyResult<Py<Self>> {
        self.read(py, |value| {
            let child = match value {
                GraphIrConstraint::And(children) | GraphIrConstraint::Or(children) => {
                    children.get(index)
                }
                GraphIrConstraint::Not(child) if index == 0 => Some(child.as_ref()),
                _ => None,
            }
            .ok_or_else(|| PyIndexError::new_err("constraint index out of range"))?;
            match &self.storage {
                ConstraintStorage::Owned(_) => Self::from_rust(py, child),
                ConstraintStorage::Molecule {
                    owner,
                    counter,
                    position,
                    path,
                } => {
                    let mut path = path.to_vec();
                    path.push(index);
                    Self::from_storage(
                        py,
                        ConstraintStorage::Molecule {
                            owner: owner.clone_ref(py),
                            counter: *counter,
                            position: *position,
                            path: path.into_boxed_slice(),
                        },
                    )
                }
            }
        })
    }

    fn from_storage(py: Python<'_>, storage: ConstraintStorage) -> PyResult<Py<Self>> {
        let base = Self { storage };
        let initialize = base.read(py, |value| {
            let initialize: fn(Python<'_>, Self) -> PyResult<Py<Self>> = match value {
                GraphIrConstraint::Atom(..) => |py, base| {
                    Ok(Py::new(
                        py,
                        PyClassInitializer::from(base).add_subclass(variants::Atom),
                    )?
                    .into_bound(py)
                    .into_super()
                    .unbind())
                },
                GraphIrConstraint::Bond(..) => |py, base| {
                    Ok(Py::new(
                        py,
                        PyClassInitializer::from(base).add_subclass(variants::Bond),
                    )?
                    .into_bound(py)
                    .into_super()
                    .unbind())
                },
                GraphIrConstraint::DativeBond(..) => |py, base| {
                    Ok(Py::new(
                        py,
                        PyClassInitializer::from(base).add_subclass(variants::DativeBond),
                    )?
                    .into_bound(py)
                    .into_super()
                    .unbind())
                },
                GraphIrConstraint::AromaticSystem(..) => |py, base| {
                    Ok(Py::new(
                        py,
                        PyClassInitializer::from(base).add_subclass(variants::AromaticSystem),
                    )?
                    .into_bound(py)
                    .into_super()
                    .unbind())
                },
                GraphIrConstraint::MulticenterBond(..) => |py, base| {
                    Ok(Py::new(
                        py,
                        PyClassInitializer::from(base).add_subclass(variants::MulticenterBond),
                    )?
                    .into_bound(py)
                    .into_super()
                    .unbind())
                },
                GraphIrConstraint::NoncovalentBond(..) => |py, base| {
                    Ok(Py::new(
                        py,
                        PyClassInitializer::from(base).add_subclass(variants::NoncovalentBond),
                    )?
                    .into_bound(py)
                    .into_super()
                    .unbind())
                },
                GraphIrConstraint::StereoAtom(..) => |py, base| {
                    Ok(Py::new(
                        py,
                        PyClassInitializer::from(base).add_subclass(variants::StereoAtom),
                    )?
                    .into_bound(py)
                    .into_super()
                    .unbind())
                },
                GraphIrConstraint::StereoBond(..) => |py, base| {
                    Ok(Py::new(
                        py,
                        PyClassInitializer::from(base).add_subclass(variants::StereoBond),
                    )?
                    .into_bound(py)
                    .into_super()
                    .unbind())
                },
                GraphIrConstraint::Relational(..) => |py, base| {
                    Ok(Py::new(
                        py,
                        PyClassInitializer::from(base).add_subclass(variants::Relational),
                    )?
                    .into_bound(py)
                    .into_super()
                    .unbind())
                },
                GraphIrConstraint::Molecule(..) => |py, base| {
                    Ok(Py::new(
                        py,
                        PyClassInitializer::from(base).add_subclass(variants::Molecule),
                    )?
                    .into_bound(py)
                    .into_super()
                    .unbind())
                },
                GraphIrConstraint::And(..) => |py, base| {
                    Ok(Py::new(
                        py,
                        PyClassInitializer::from(base).add_subclass(variants::And),
                    )?
                    .into_bound(py)
                    .into_super()
                    .unbind())
                },
                GraphIrConstraint::Or(..) => |py, base| {
                    Ok(Py::new(
                        py,
                        PyClassInitializer::from(base).add_subclass(variants::Or),
                    )?
                    .into_bound(py)
                    .into_super()
                    .unbind())
                },
                GraphIrConstraint::Not(..) => |py, base| {
                    Ok(Py::new(
                        py,
                        PyClassInitializer::from(base).add_subclass(variants::Not),
                    )?
                    .into_bound(py)
                    .into_super()
                    .unbind())
                },
            };
            Ok(initialize)
        })?;
        initialize(py, base)
    }
}

#[allow(
    clippy::just_underscores_and_digits,
    reason = "Python tuple-variant constructor keywords are _0, _1, and _2."
)]
mod variants {
    use super::*;

    macro_rules! entity_variant {
        ($name:ident, $form:ident, $id:ident) => {
            #[pyclass(extends = Constraint, frozen, module = "umol")]
            pub struct $name;
            #[pymethods]
            impl $name {
                #[classattr]
                fn __qualname__() -> &'static str {
                    concat!("Constraint.", stringify!($name))
                }

                #[new]
                fn new(py: Python<'_>, _0: u32, _1: Py<$form>) -> PyClassInitializer<Self> {
                    PyClassInitializer::from(Constraint {
                        storage: ConstraintStorage::Owned(GraphIrConstraint::$name(
                            $id(_0),
                            _1.borrow(py).to_rust(py),
                        )),
                    })
                    .add_subclass(Self)
                }
                #[classattr]
                fn __match_args__() -> (&'static str, &'static str) {
                    ("_0", "_1")
                }
                #[getter]
                fn _0(slf: PyRef<'_, Self>, py: Python<'_>) -> PyResult<u32> {
                    slf.as_super().read(py, |value| match value {
                        GraphIrConstraint::$name(id, _) => Ok(id.0),
                        _ => unreachable!(),
                    })
                }
                #[getter]
                fn _1(slf: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<$form>> {
                    slf.as_super().read(py, |value| match value {
                        GraphIrConstraint::$name(_, form) => {
                            into_py_variant(py, $form::from_rust(py, form)?)
                        }
                        _ => unreachable!(),
                    })
                }
            }
        };
    }
    entity_variant!(Atom, AtomConstraintForm, GraphIrAtomId);
    entity_variant!(Bond, BondConstraintForm, GraphIrBondId);
    entity_variant!(DativeBond, DativeBondConstraintForm, GraphIrDativeBondId);
    entity_variant!(
        AromaticSystem,
        AromaticSystemConstraintForm,
        GraphIrAromaticSystemId
    );
    entity_variant!(
        MulticenterBond,
        MulticenterBondConstraintForm,
        GraphIrMulticenterBondId
    );
    entity_variant!(
        NoncovalentBond,
        NoncovalentBondConstraintForm,
        GraphIrNoncovalentBondId
    );

    #[pyclass(extends = Constraint, frozen, module = "umol")]
    pub struct StereoAtom;
    #[pymethods]
    impl StereoAtom {
        #[classattr]
        fn __qualname__() -> &'static str {
            "Constraint.StereoAtom"
        }

        #[new]
        fn new(
            py: Python<'_>,
            _0: u32,
            _1: StereoKind,
            _2: Py<StereoAtomConstraintForm>,
        ) -> PyClassInitializer<Self> {
            PyClassInitializer::from(Constraint {
                storage: ConstraintStorage::Owned(GraphIrConstraint::StereoAtom(
                    GraphIrStereoAtomId(_0),
                    _1.to_rust(),
                    _2.borrow(py).to_rust(py),
                )),
            })
            .add_subclass(Self)
        }
        #[classattr]
        fn __match_args__() -> (&'static str, &'static str, &'static str) {
            ("_0", "_1", "_2")
        }
        #[getter]
        fn _0(slf: PyRef<'_, Self>, py: Python<'_>) -> PyResult<u32> {
            slf.as_super().read(py, |value| match value {
                GraphIrConstraint::StereoAtom(id, _, _) => Ok(id.0),
                _ => unreachable!(),
            })
        }
        #[getter]
        fn _1(slf: PyRef<'_, Self>, py: Python<'_>) -> PyResult<StereoKind> {
            slf.as_super().read(py, |value| match value {
                GraphIrConstraint::StereoAtom(_, kind, _) => Ok(StereoKind::from_rust(*kind)),
                _ => unreachable!(),
            })
        }
        #[getter]
        fn _2(slf: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<StereoAtomConstraintForm>> {
            slf.as_super().read(py, |value| match value {
                GraphIrConstraint::StereoAtom(_, _, form) => {
                    into_py_variant(py, StereoAtomConstraintForm::from_rust(py, form)?)
                }
                _ => unreachable!(),
            })
        }
    }

    #[pyclass(extends = Constraint, frozen, module = "umol")]
    pub struct StereoBond;
    #[pymethods]
    impl StereoBond {
        #[classattr]
        fn __qualname__() -> &'static str {
            "Constraint.StereoBond"
        }

        #[new]
        fn new(
            py: Python<'_>,
            _0: u32,
            _1: StereoKind,
            _2: Py<StereoBondConstraintForm>,
        ) -> PyClassInitializer<Self> {
            PyClassInitializer::from(Constraint {
                storage: ConstraintStorage::Owned(GraphIrConstraint::StereoBond(
                    GraphIrStereoBondId(_0),
                    _1.to_rust(),
                    _2.borrow(py).to_rust(py),
                )),
            })
            .add_subclass(Self)
        }
        #[classattr]
        fn __match_args__() -> (&'static str, &'static str, &'static str) {
            ("_0", "_1", "_2")
        }
        #[getter]
        fn _0(slf: PyRef<'_, Self>, py: Python<'_>) -> PyResult<u32> {
            slf.as_super().read(py, |value| match value {
                GraphIrConstraint::StereoBond(id, _, _) => Ok(id.0),
                _ => unreachable!(),
            })
        }
        #[getter]
        fn _1(slf: PyRef<'_, Self>, py: Python<'_>) -> PyResult<StereoKind> {
            slf.as_super().read(py, |value| match value {
                GraphIrConstraint::StereoBond(_, kind, _) => Ok(StereoKind::from_rust(*kind)),
                _ => unreachable!(),
            })
        }
        #[getter]
        fn _2(slf: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<StereoBondConstraintForm>> {
            slf.as_super().read(py, |value| match value {
                GraphIrConstraint::StereoBond(_, _, form) => {
                    into_py_variant(py, StereoBondConstraintForm::from_rust(py, form)?)
                }
                _ => unreachable!(),
            })
        }
    }

    #[pyclass(extends = Constraint, frozen, module = "umol")]
    pub struct Relational;
    #[pymethods]
    impl Relational {
        #[classattr]
        fn __qualname__() -> &'static str {
            "Constraint.Relational"
        }

        #[new]
        fn new(py: Python<'_>, _0: Py<RelationalConstraint>) -> PyClassInitializer<Self> {
            PyClassInitializer::from(Constraint {
                storage: ConstraintStorage::Owned(GraphIrConstraint::Relational(
                    _0.borrow(py).to_rust(py),
                )),
            })
            .add_subclass(Self)
        }
        #[classattr]
        fn __match_args__() -> (&'static str,) {
            ("_0",)
        }
        #[getter]
        fn _0(slf: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<RelationalConstraint>> {
            slf.as_super().read(py, |value| match value {
                GraphIrConstraint::Relational(form) => {
                    into_py_variant(py, RelationalConstraint::from_rust(py, form)?)
                }
                _ => unreachable!(),
            })
        }
    }

    #[pyclass(extends = Constraint, frozen, module = "umol")]
    pub struct Molecule;
    #[pymethods]
    impl Molecule {
        #[classattr]
        fn __qualname__() -> &'static str {
            "Constraint.Molecule"
        }

        #[new]
        fn new(py: Python<'_>, _0: Py<MoleculeConstraint>) -> PyClassInitializer<Self> {
            PyClassInitializer::from(Constraint {
                storage: ConstraintStorage::Owned(GraphIrConstraint::Molecule(
                    _0.borrow(py).to_rust(py),
                )),
            })
            .add_subclass(Self)
        }
        #[classattr]
        fn __match_args__() -> (&'static str,) {
            ("_0",)
        }
        #[getter]
        fn _0(slf: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<MoleculeConstraint>> {
            slf.as_super().read(py, |value| match value {
                GraphIrConstraint::Molecule(form) => {
                    into_py_variant(py, MoleculeConstraint::from_rust(py, form)?)
                }
                _ => unreachable!(),
            })
        }
    }

    #[pyclass(extends = Constraint, frozen, module = "umol")]
    pub struct And;
    #[pymethods]
    impl And {
        #[classattr]
        fn __qualname__() -> &'static str {
            "Constraint.And"
        }

        #[new]
        fn new(py: Python<'_>, _0: Vec<Py<Constraint>>) -> PyResult<PyClassInitializer<Self>> {
            let children = _0
                .iter()
                .map(|child| child.borrow(py).to_rust(py))
                .collect::<PyResult<Vec<_>>>()?;
            Ok(PyClassInitializer::from(Constraint {
                storage: ConstraintStorage::Owned(GraphIrConstraint::And(children)),
            })
            .add_subclass(Self))
        }
        #[classattr]
        fn __match_args__() -> (&'static str,) {
            ("_0",)
        }
        #[getter]
        fn _0(slf: PyRef<'_, Self>, py: Python<'_>) -> PyResult<ConstraintsView> {
            slf.as_super().read(py, |_| Ok(()))?;
            Ok(ConstraintsView {
                storage: ConstraintsStorage::Children(slf.into_super().into()),
            })
        }
    }

    #[pyclass(extends = Constraint, frozen, module = "umol")]
    pub struct Or;
    #[pymethods]
    impl Or {
        #[classattr]
        fn __qualname__() -> &'static str {
            "Constraint.Or"
        }

        #[new]
        fn new(py: Python<'_>, _0: Vec<Py<Constraint>>) -> PyResult<PyClassInitializer<Self>> {
            let children = _0
                .iter()
                .map(|child| child.borrow(py).to_rust(py))
                .collect::<PyResult<Vec<_>>>()?;
            Ok(PyClassInitializer::from(Constraint {
                storage: ConstraintStorage::Owned(GraphIrConstraint::Or(children)),
            })
            .add_subclass(Self))
        }
        #[classattr]
        fn __match_args__() -> (&'static str,) {
            ("_0",)
        }
        #[getter]
        fn _0(slf: PyRef<'_, Self>, py: Python<'_>) -> PyResult<ConstraintsView> {
            slf.as_super().read(py, |_| Ok(()))?;
            Ok(ConstraintsView {
                storage: ConstraintsStorage::Children(slf.into_super().into()),
            })
        }
    }

    #[pyclass(extends = Constraint, frozen, module = "umol")]
    pub struct Not;
    #[pymethods]
    impl Not {
        #[classattr]
        fn __qualname__() -> &'static str {
            "Constraint.Not"
        }

        #[new]
        fn new(py: Python<'_>, _0: Py<Constraint>) -> PyResult<PyClassInitializer<Self>> {
            Ok(PyClassInitializer::from(Constraint {
                storage: ConstraintStorage::Owned(GraphIrConstraint::Not(Box::new(
                    _0.borrow(py).to_rust(py)?,
                ))),
            })
            .add_subclass(Self))
        }
        #[classattr]
        fn __match_args__() -> (&'static str,) {
            ("_0",)
        }
        #[getter]
        fn _0(slf: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Constraint>> {
            slf.as_super().child(py, 0)
        }
    }
}

pub(crate) fn register_constraint(py: Python<'_>) -> PyResult<()> {
    PySequence::register::<ConstraintsView>(py)?;
    let base = py.get_type::<Constraint>();
    base.setattr("Atom", py.get_type::<variants::Atom>())?;
    base.setattr("Bond", py.get_type::<variants::Bond>())?;
    base.setattr("DativeBond", py.get_type::<variants::DativeBond>())?;
    base.setattr("AromaticSystem", py.get_type::<variants::AromaticSystem>())?;
    base.setattr(
        "MulticenterBond",
        py.get_type::<variants::MulticenterBond>(),
    )?;
    base.setattr(
        "NoncovalentBond",
        py.get_type::<variants::NoncovalentBond>(),
    )?;
    base.setattr("StereoAtom", py.get_type::<variants::StereoAtom>())?;
    base.setattr("StereoBond", py.get_type::<variants::StereoBond>())?;
    base.setattr("Relational", py.get_type::<variants::Relational>())?;
    base.setattr("Molecule", py.get_type::<variants::Molecule>())?;
    base.setattr("And", py.get_type::<variants::And>())?;
    base.setattr("Or", py.get_type::<variants::Or>())?;
    base.setattr("Not", py.get_type::<variants::Not>())?;
    Ok(())
}

fn constraint_fields(
    py: Python<'_>,
    value: &GraphIrConstraint,
) -> PyResult<(&'static str, Vec<String>)> {
    Ok(match value {
        GraphIrConstraint::Atom(id, form) => (
            "Atom",
            vec![
                id.0.to_string(),
                into_py_variant(py, AtomConstraintForm::from_rust(py, form)?)?
                    .bind(py)
                    .as_any()
                    .repr()?
                    .extract()?,
            ],
        ),
        GraphIrConstraint::Bond(id, form) => (
            "Bond",
            vec![
                id.0.to_string(),
                into_py_variant(py, BondConstraintForm::from_rust(py, form)?)?
                    .bind(py)
                    .as_any()
                    .repr()?
                    .extract()?,
            ],
        ),
        GraphIrConstraint::DativeBond(id, form) => (
            "DativeBond",
            vec![
                id.0.to_string(),
                into_py_variant(py, DativeBondConstraintForm::from_rust(py, form)?)?
                    .bind(py)
                    .as_any()
                    .repr()?
                    .extract()?,
            ],
        ),
        GraphIrConstraint::AromaticSystem(id, form) => (
            "AromaticSystem",
            vec![
                id.0.to_string(),
                into_py_variant(py, AromaticSystemConstraintForm::from_rust(py, form)?)?
                    .bind(py)
                    .as_any()
                    .repr()?
                    .extract()?,
            ],
        ),
        GraphIrConstraint::MulticenterBond(id, form) => (
            "MulticenterBond",
            vec![
                id.0.to_string(),
                into_py_variant(py, MulticenterBondConstraintForm::from_rust(py, form)?)?
                    .bind(py)
                    .as_any()
                    .repr()?
                    .extract()?,
            ],
        ),
        GraphIrConstraint::NoncovalentBond(id, form) => (
            "NoncovalentBond",
            vec![
                id.0.to_string(),
                into_py_variant(py, NoncovalentBondConstraintForm::from_rust(py, form)?)?
                    .bind(py)
                    .as_any()
                    .repr()?
                    .extract()?,
            ],
        ),
        GraphIrConstraint::StereoAtom(id, kind, form) => (
            "StereoAtom",
            vec![
                id.0.to_string(),
                StereoKind::from_rust(*kind)
                    .into_pyobject(py)?
                    .as_any()
                    .repr()?
                    .extract()?,
                into_py_variant(py, StereoAtomConstraintForm::from_rust(py, form)?)?
                    .bind(py)
                    .as_any()
                    .repr()?
                    .extract()?,
            ],
        ),
        GraphIrConstraint::StereoBond(id, kind, form) => (
            "StereoBond",
            vec![
                id.0.to_string(),
                StereoKind::from_rust(*kind)
                    .into_pyobject(py)?
                    .as_any()
                    .repr()?
                    .extract()?,
                into_py_variant(py, StereoBondConstraintForm::from_rust(py, form)?)?
                    .bind(py)
                    .as_any()
                    .repr()?
                    .extract()?,
            ],
        ),
        GraphIrConstraint::Relational(form) => (
            "Relational",
            vec![
                into_py_variant(py, RelationalConstraint::from_rust(py, form)?)?
                    .bind(py)
                    .as_any()
                    .repr()?
                    .extract()?,
            ],
        ),
        GraphIrConstraint::Molecule(form) => (
            "Molecule",
            vec![
                into_py_variant(py, MoleculeConstraint::from_rust(py, form)?)?
                    .bind(py)
                    .as_any()
                    .repr()?
                    .extract()?,
            ],
        ),
        GraphIrConstraint::And(children) | GraphIrConstraint::Or(children) => {
            let children = children
                .iter()
                .map(|child| {
                    let (name, fields) = constraint_fields(py, child)?;
                    Ok(format!("Constraint.{name}({})", fields.join(", ")))
                })
                .collect::<PyResult<Vec<_>>>()?;
            (
                if matches!(value, GraphIrConstraint::And(_)) {
                    "And"
                } else {
                    "Or"
                },
                vec![format!("[{}]", children.join(", "))],
            )
        }
        GraphIrConstraint::Not(child) => {
            let (name, fields) = constraint_fields(py, child)?;
            (
                "Not",
                vec![format!("Constraint.{name}({})", fields.join(", "))],
            )
        }
    })
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

fn constraint_iter(py: Python<'_>, constraints: &GraphIrConstraints) -> PyResult<Py<PyAny>> {
    let entries = constraints
        .iter()
        .map(|value| Constraint::from_rust(py, value))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(PyList::new(py, entries)?.try_iter()?.into_any().unbind())
}

#[pyclass]
pub(crate) struct ConstraintIter {
    collection: Py<ConstraintsView>,
    position: usize,
    end: usize,
}

#[pymethods]
impl ConstraintIter {
    fn __iter__<'py>(slf: PyRef<'py, Self>, py: Python<'py>) -> PyResult<PyRef<'py, Self>> {
        slf.collection.try_borrow(py)?.read(py, |_| Ok(()))?;
        Ok(slf)
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<Py<Constraint>>> {
        let collection = self.collection.try_borrow(py)?;
        collection.read(py, |_| Ok(()))?;
        if self.position == self.end {
            return Ok(None);
        }
        let value = collection.__getitem__(py, self.position as isize)?;
        self.position += 1;
        Ok(Some(value))
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
                ResolvedConstraintsUpdate::Overlay(container.bind(py).borrow().to_rust().clone())
            }
            Self::View(view) => ResolvedConstraintsUpdate::Overlay(
                view.bind(py).borrow().read(py, |constraints| {
                    Ok(GraphIrConstraints::from(constraints.to_vec()))
                })?,
            ),
            Self::Entries(entries) => ResolvedConstraintsUpdate::Entries(
                entries
                    .iter()
                    .map(|entry| entry.bind(py).borrow().to_rust(py))
                    .collect::<PyResult<Vec<_>>>()?,
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

/// The molecule-level constraints in insertion order. Mutable, value-equal,
/// and unhashable.
#[pyclass(eq)]
#[derive(Debug, PartialEq)]
pub struct Constraints(GraphIrConstraints);

#[pymethods]
impl Constraints {
    /// Build an owned container from constraint entries, preserving order and duplicates.
    #[new]
    fn new(py: Python<'_>, entries: Vec<Py<Constraint>>) -> PyResult<Self> {
        Ok(Self(GraphIrConstraints::from(
            entries
                .into_iter()
                .map(|entry| entry.bind(py).borrow().to_rust(py))
                .collect::<PyResult<Vec<_>>>()?,
        )))
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let mut parts = Vec::with_capacity(self.0.len());
        for entry in self.0.iter() {
            let value = Constraint::from_rust(py, entry)?;
            parts.push(value.bind(py).as_any().repr()?.extract::<String>()?);
        }
        Ok(format!("Constraints([{}])", parts.join(", ")))
    }

    /// Append one constraint, preserving existing entries and duplicates.
    fn append(&mut self, py: Python<'_>, constraint: Py<Constraint>) -> PyResult<()> {
        self.0.push(constraint.bind(py).borrow().to_rust(py)?);
        Ok(())
    }

    fn clear(&mut self) {
        self.0.clear();
    }

    /// Append another container, live view, or iterable after snapshotting the RHS.
    fn update(slf: Py<Self>, py: Python<'_>, other: ConstraintsUpdate) -> PyResult<()> {
        let resolved = other.resolve(py)?;
        resolved.apply(slf.borrow_mut(py).to_rust_mut());
        Ok(())
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    fn __getitem__(&self, py: Python<'_>, index: isize) -> PyResult<Py<Constraint>> {
        let index = resolve_constraint_index(self.0.len(), index)?;
        Constraint::from_rust(py, &self.0.as_slice()[index])
    }

    fn __iter__(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        constraint_iter(py, &self.0)
    }
}

impl Constraints {
    pub(crate) fn to_rust(&self) -> &GraphIrConstraints {
        &self.0
    }

    pub(crate) fn to_rust_mut(&mut self) -> &mut GraphIrConstraints {
        &mut self.0
    }

    pub(crate) fn from_rust(constraints: GraphIrConstraints) -> Self {
        Self(constraints)
    }
}

impl_py_normalize!(
    Constraints,
    GraphIrConstraints,
    |value: &Constraints, _py: Python<'_>| -> PyResult<GraphIrConstraints> {
        Ok(value.to_rust().clone())
    },
    |_py: Python<'_>, value: GraphIrConstraints| -> PyResult<Constraints> {
        Ok(Constraints::from_rust(value))
    }
);

/// Read-only access to a molecule's constraints or a composition's children.
#[pyclass(sequence)]
pub struct ConstraintsView {
    storage: ConstraintsStorage,
}

enum ConstraintsStorage {
    Molecule { owner: Py<Molecule>, counter: u64 },
    Children(Py<Constraint>),
}

impl ConstraintsView {
    pub(crate) fn new(owner: Py<Molecule>, py: Python<'_>) -> PyResult<Self> {
        let counter = owner.try_borrow(py)?.view_counter()?;
        Ok(Self {
            storage: ConstraintsStorage::Molecule { owner, counter },
        })
    }

    pub(crate) fn read<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&[GraphIrConstraint]) -> PyResult<R>,
    ) -> PyResult<R> {
        match &self.storage {
            ConstraintsStorage::Molecule { owner, counter } => {
                let owner = owner.try_borrow(py)?;
                owner.check_access(*counter, "ConstraintsView")?;
                f(owner.to_rust().constraints().as_slice())
            }
            ConstraintsStorage::Children(parent) => {
                parent.try_borrow(py)?.read(py, |value| match value {
                    GraphIrConstraint::And(children) | GraphIrConstraint::Or(children) => {
                        f(children)
                    }
                    _ => unreachable!(),
                })
            }
        }
    }
}

#[pymethods]
impl ConstraintsView {
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        self.read(py, |values| {
            Ok(format!("ConstraintsView({} entries)", values.len()))
        })
    }

    fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        self.read(py, |values| Ok(values.len()))
    }

    fn __getitem__(&self, py: Python<'_>, index: isize) -> PyResult<Py<Constraint>> {
        let index = self.read(py, |values| resolve_constraint_index(values.len(), index))?;
        match &self.storage {
            ConstraintsStorage::Molecule { owner, counter } => Constraint::from_storage(
                py,
                ConstraintStorage::Molecule {
                    owner: owner.clone_ref(py),
                    counter: *counter,
                    position: index,
                    path: Box::default(),
                },
            ),
            ConstraintsStorage::Children(parent) => parent.try_borrow(py)?.child(py, index),
        }
    }

    fn __iter__(slf: Py<Self>, py: Python<'_>) -> PyResult<ConstraintIter> {
        let end = slf.try_borrow(py)?.__len__(py)?;
        Ok(ConstraintIter {
            collection: slf,
            position: 0,
            end,
        })
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

impl_py_normalize!(
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

impl_py_normalize!(
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
        AtomConstraintForm as GraphIrAtomConstraintForm, AtomForm as GraphIrAtomForm,
        BondConstraintForm as GraphIrBondConstraintForm,
        DativeBondConstraintForm as GraphIrDativeBondConstraintForm, Molecule as GraphIrMolecule,
        MoleculeEntries as GraphIrMoleculeEntries,
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
            assert_eq!(value.borrow(py).to_rust(py).unwrap(), constraint);
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
            register_constraint(py).unwrap();
            let value = Constraint::from_rust(py, &constraint).unwrap();
            let equal = Constraint::from_rust(py, &constraint).unwrap();

            assert_eq!(value.bind(py).borrow().to_rust(py).unwrap(), constraint);
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
            let iterator = constraint_iter(py, &constraints).unwrap();
            let mut iter = iterator.bind(py).try_iter().unwrap();
            constraints.push(GraphIrConstraint::Or(Vec::new()));

            let first_mirror = iter
                .next()
                .unwrap()
                .unwrap()
                .extract::<Py<Constraint>>()
                .unwrap();
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
            assert_eq!(first_mirror.bind(py).borrow().to_rust(py).unwrap(), first);
            assert_eq!(
                iter.next()
                    .unwrap()
                    .unwrap()
                    .extract::<Py<Constraint>>()
                    .unwrap()
                    .borrow(py)
                    .to_rust(py)
                    .unwrap(),
                second
            );
            assert!(iter.next().is_none());
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
                .map(|entry| Constraint::from_rust(py, entry).unwrap())
                .collect();
            let constraints = Constraints::new(py, values).unwrap();

            assert_eq!(constraints.to_rust().as_slice(), entries.as_slice());
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
            Constraints::from_rust(left) == Constraints::from_rust(right),
            expected
        );
    }

    #[rstest]
    fn test_constraints_repr() {
        Python::attach(|py| {
            let constraints = Constraints::from_rust(GraphIrConstraints::from(vec![
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
                Constraints::from_rust(GraphIrConstraints::from(vec![constraint.clone()]));
            let value = Constraint::from_rust(py, &constraint).unwrap();

            constraints.append(py, value).unwrap();

            assert_eq!(
                constraints.to_rust().as_slice(),
                &[constraint.clone(), constraint]
            );
        });
    }

    #[rstest]
    fn test_constraints_clear() {
        let mut constraints =
            Constraints::from_rust(GraphIrConstraints::from(vec![GraphIrConstraint::And(
                Vec::new(),
            )]));

        constraints.clear();

        assert_eq!(constraints.to_rust(), &GraphIrConstraints::new());
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
                Constraints::from_rust(GraphIrConstraints::from(vec![initial.clone()])),
            )
            .unwrap();
            let container = Py::new(
                py,
                Constraints::from_rust(GraphIrConstraints::from(vec![from_container.clone()])),
            )
            .unwrap();
            let mut molecule = GraphIrMolecule::new();
            let mut editor = molecule.edit();
            editor.constraints_mut().push(from_view.clone());
            molecule = editor.build();
            let view = Py::new(
                py,
                ConstraintsView::new(Py::new(py, Molecule::from_rust(molecule)).unwrap(), py)
                    .unwrap(),
            )
            .unwrap();
            let entry = Constraint::from_rust(py, &from_entries).unwrap();

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
                target.bind(py).borrow().to_rust().as_slice(),
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
                Constraints::from_rust(GraphIrConstraints::from(vec![entry.clone()])),
            )
            .unwrap();

            Constraints::update(
                target.clone_ref(py),
                py,
                ConstraintsUpdate::Container(target.clone_ref(py)),
            )
            .unwrap();

            assert_eq!(
                target.bind(py).borrow().to_rust().as_slice(),
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
    fn test_constraints_len(#[case] constraints: GraphIrConstraints, #[case] expected: usize) {
        assert_eq!(Constraints::from_rust(constraints).__len__(), expected);
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
            let constraints = Constraints::from_rust(GraphIrConstraints::from(vec![
                GraphIrConstraint::Atom(GraphIrAtomId(1), GraphIrAtomConstraintForm::degree(2)),
                GraphIrConstraint::Molecule(GraphIrMoleculeConstraint::Connected { atoms: None }),
            ]));
            let actual = constraints.__getitem__(py, index).unwrap();

            assert_eq!(actual.borrow(py).to_rust(py).unwrap(), expected);
        });
    }

    #[rstest]
    #[case::positive(1)]
    #[case::negative(-2)]
    fn test_constraints_getitem_error(#[case] index: isize) {
        Python::attach(|py| {
            let constraints =
                Constraints::from_rust(GraphIrConstraints::from(vec![GraphIrConstraint::And(
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
        let mut constraints = Constraints::from_rust(GraphIrConstraints::from(vec![
            first.clone(),
            second.clone(),
        ]));

        Python::attach(|py| {
            let iterator = constraints.__iter__(py).unwrap();
            let mut iter = iterator.bind(py).try_iter().unwrap();
            constraints
                .to_rust_mut()
                .push(GraphIrConstraint::Not(Box::new(GraphIrConstraint::And(
                    Vec::new(),
                ))));

            assert_eq!(
                iter.next()
                    .unwrap()
                    .unwrap()
                    .extract::<Py<Constraint>>()
                    .unwrap()
                    .borrow(py)
                    .to_rust(py)
                    .unwrap(),
                first
            );
            assert_eq!(
                iter.next()
                    .unwrap()
                    .unwrap()
                    .extract::<Py<Constraint>>()
                    .unwrap()
                    .borrow(py)
                    .to_rust(py)
                    .unwrap(),
                second
            );
            assert!(iter.next().is_none());
        });
    }

    #[rstest]
    fn test_constraints_from_rust() {
        let entries = vec![GraphIrConstraint::Or(Vec::new())];

        assert_eq!(
            Constraints::from_rust(GraphIrConstraints::from(entries.clone()))
                .to_rust()
                .as_slice(),
            entries.as_slice()
        );
    }

    #[rstest]
    fn test_constraints_view_repr() {
        let mut molecule = GraphIrMolecule::new();
        let mut editor = molecule.edit();
        editor
            .constraints_mut()
            .push(GraphIrConstraint::And(Vec::new()));
        editor
            .constraints_mut()
            .push(GraphIrConstraint::Or(Vec::new()));

        molecule = editor.build();

        Python::attach(|py| {
            let owner = Py::new(py, Molecule::from_rust(molecule)).unwrap();
            let view = ConstraintsView::new(owner, py).unwrap();

            assert_eq!(view.__repr__(py).unwrap(), "ConstraintsView(2 entries)");
        });
    }

    #[rstest]
    #[case::empty(Vec::new(), 0)]
    #[case::nonempty(vec![GraphIrConstraint::And(Vec::new())], 1)]
    fn test_constraints_view_len(#[case] entries: Vec<GraphIrConstraint>, #[case] expected: usize) {
        let mut editor = GraphIrMolecule::new().edit();
        *editor.constraints_mut() = entries.into();
        Python::attach(|py| {
            let owner = Py::new(py, Molecule::from_rust(editor.build())).unwrap();
            let view = ConstraintsView::new(owner, py).unwrap();
            assert_eq!(view.__len__(py).unwrap(), expected);
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
        let mut molecule = GraphIrMolecule::from_entries(GraphIrMoleculeEntries {
            atoms: vec![GraphIrAtomForm::default(), GraphIrAtomForm::default()],
            ..Default::default()
        });
        let mut editor = molecule.edit();
        editor.constraints_mut().push(GraphIrConstraint::Atom(
            GraphIrAtomId(1),
            GraphIrAtomConstraintForm::degree(2),
        ));
        editor.constraints_mut().push(GraphIrConstraint::Molecule(
            GraphIrMoleculeConstraint::Connected { atoms: None },
        ));
        molecule = editor.build();

        Python::attach(|py| {
            let owner = Py::new(py, Molecule::from_rust(molecule)).unwrap();
            let view = ConstraintsView::new(owner, py).unwrap();

            assert_eq!(
                view.__getitem__(py, index)
                    .unwrap()
                    .borrow(py)
                    .to_rust(py)
                    .unwrap(),
                expected
            );
        });
    }

    #[rstest]
    #[case::positive(1)]
    #[case::negative(-2)]
    fn test_constraints_view_getitem_error(#[case] index: isize) {
        let mut molecule = GraphIrMolecule::new();
        let mut editor = molecule.edit();
        editor
            .constraints_mut()
            .push(GraphIrConstraint::And(Vec::new()));

        molecule = editor.build();

        Python::attach(|py| {
            let owner = Py::new(py, Molecule::from_rust(molecule)).unwrap();
            let view = ConstraintsView::new(owner, py).unwrap();

            assert_eq!(
                view.__getitem__(py, index).err().unwrap().to_string(),
                "IndexError: constraint index out of range"
            );
        });
    }

    #[rstest]
    fn test_constraints_view_iter() {
        let first = GraphIrConstraint::And(vec![GraphIrConstraint::Or(Vec::new())]);
        let second = GraphIrConstraint::Not(Box::new(GraphIrConstraint::And(Vec::new())));
        let mut editor = GraphIrMolecule::new().edit();
        editor.constraints_mut().push(first.clone());
        editor.constraints_mut().push(second.clone());
        let molecule = editor.build();
        Python::attach(|py| {
            let owner = Py::new(py, Molecule::from_rust(molecule)).unwrap();
            let view = Py::new(py, ConstraintsView::new(owner, py).unwrap()).unwrap();
            let mut iter = ConstraintsView::__iter__(view, py).unwrap();
            assert_eq!(
                iter.__next__(py)
                    .unwrap()
                    .unwrap()
                    .borrow(py)
                    .to_rust(py)
                    .unwrap(),
                first
            );
            assert_eq!(
                iter.__next__(py)
                    .unwrap()
                    .unwrap()
                    .borrow(py)
                    .to_rust(py)
                    .unwrap(),
                second
            );
            assert!(iter.__next__(py).unwrap().is_none());
        });
    }
}
