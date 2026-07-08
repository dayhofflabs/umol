//! Atom-constraint sub-ASTs mirroring `umol_ast::ast::constraint` (S5a): the
//! aromatic/multicenter valence states, ring scope, and ring membership. The
//! `AtomConstraint` enum and `AtomConstraints` container follow at S5b.

use std::vec::IntoIter;

use pyo3::prelude::*;
use umol_ast::ast::{
    AromaticValenceAst as AstAromaticValenceAst, AtomConstraint as AstAtomConstraint,
    AtomConstraintKey as AstAtomConstraintKey, AtomConstraints as AstAtomConstraints,
    MulticenterValenceAst as AstMulticenterValenceAst, RingMembershipAst as AstRingMembershipAst,
    RingScope as AstRingScope,
};

use crate::convert::into_py_variant;
use crate::stereo::TetrahedralStereoAst;
use crate::value::{ValueArg, ValueAst};

/// Aromatic-valence state: undetermined, explicitly not aromatic, or aromatic with
/// an aromatic-valence count.
#[pyclass]
pub enum AromaticValenceAst {
    Undetermined(),
    NotAromatic(),
    Aromatic(Py<ValueAst>),
}

impl AromaticValenceAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstAromaticValenceAst) -> PyResult<Self> {
        Ok(match ast {
            AstAromaticValenceAst::Undetermined => Self::Undetermined(),
            AstAromaticValenceAst::NotAromatic => Self::NotAromatic(),
            AstAromaticValenceAst::Aromatic(v) => {
                Self::Aromatic(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstAromaticValenceAst {
        match self {
            Self::Undetermined() => AstAromaticValenceAst::Undetermined,
            Self::NotAromatic() => AstAromaticValenceAst::NotAromatic,
            Self::Aromatic(v) => AstAromaticValenceAst::Aromatic(v.bind(py).borrow().to_ast(py)),
        }
    }
}

/// Multicenter-valence state: undetermined, explicitly not multicenter, or
/// multicenter with a multicenter-valence count.
#[pyclass]
pub enum MulticenterValenceAst {
    Undetermined(),
    NotMulticenter(),
    Multicenter(Py<ValueAst>),
}

impl MulticenterValenceAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstMulticenterValenceAst) -> PyResult<Self> {
        Ok(match ast {
            AstMulticenterValenceAst::Undetermined => Self::Undetermined(),
            AstMulticenterValenceAst::NotMulticenter => Self::NotMulticenter(),
            AstMulticenterValenceAst::Multicenter(v) => {
                Self::Multicenter(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstMulticenterValenceAst {
        match self {
            Self::Undetermined() => AstMulticenterValenceAst::Undetermined,
            Self::NotMulticenter() => AstMulticenterValenceAst::NotMulticenter,
            Self::Multicenter(v) => {
                AstMulticenterValenceAst::Multicenter(v.bind(py).borrow().to_ast(py))
            }
        }
    }
}

/// Ring scope: all rings, or rings of a given size.
#[pyclass]
pub enum RingScope {
    All(),
    Size(u8),
}

impl RingScope {
    pub(crate) fn from_ast(ast: &AstRingScope) -> Self {
        match ast {
            AstRingScope::All => Self::All(),
            AstRingScope::Size(size) => Self::Size(*size),
        }
    }

    pub(crate) fn to_ast(&self) -> AstRingScope {
        match self {
            Self::All() => AstRingScope::All,
            Self::Size(size) => AstRingScope::Size(*size),
        }
    }
}

/// Ring-membership fact: a ring scope and a membership count.
#[pyclass]
pub struct RingMembershipAst {
    #[pyo3(get)]
    scope: Py<RingScope>,
    #[pyo3(get)]
    count: Py<ValueAst>,
}

#[pymethods]
impl RingMembershipAst {
    #[new]
    fn new(py: Python<'_>, scope: Py<RingScope>, count: ValueArg) -> PyResult<Self> {
        Ok(RingMembershipAst {
            scope,
            count: count.to_py(py)?,
        })
    }
}

impl RingMembershipAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstRingMembershipAst) -> PyResult<Self> {
        Ok(RingMembershipAst {
            scope: into_py_variant(py, RingScope::from_ast(&ast.scope))?,
            count: into_py_variant(py, ValueAst::from_ast(py, &ast.count)?)?,
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstRingMembershipAst {
        AstRingMembershipAst::new(
            self.scope.bind(py).borrow().to_ast(),
            self.count.bind(py).borrow().to_ast(py),
        )
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
pub enum AtomConstraint {
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
impl AtomConstraint {
    /// The constraint's key (identity).
    #[getter]
    fn key(&self, py: Python<'_>) -> PyResult<AtomConstraintKey> {
        AtomConstraintKey::from_ast(py, &self.to_ast(py).key())
    }
}

impl AtomConstraint {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstAtomConstraint) -> PyResult<Self> {
        Ok(match ast {
            AstAtomConstraint::Valence(v) => {
                Self::Valence(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraint::TotalValence(v) => {
                Self::TotalValence(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraint::AromaticValence(c) => {
                Self::AromaticValence(into_py_variant(py, AromaticValenceAst::from_ast(py, c)?)?)
            }
            AstAtomConstraint::MulticenterValence(c) => Self::MulticenterValence(into_py_variant(
                py,
                MulticenterValenceAst::from_ast(py, c)?,
            )?),
            AstAtomConstraint::DonatedPairs(v) => {
                Self::DonatedPairs(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraint::AcceptedPairs(v) => {
                Self::AcceptedPairs(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraint::Degree(v) => {
                Self::Degree(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraint::TotalDegree(v) => {
                Self::TotalDegree(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraint::RingDegree(v) => {
                Self::RingDegree(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraint::RingValence(v) => {
                Self::RingValence(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraint::TotalHydrogens(v) => {
                Self::TotalHydrogens(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
            AstAtomConstraint::RingMembership(m) => {
                Self::RingMembership(into_py_variant(py, RingMembershipAst::from_ast(py, m)?)?)
            }
            AstAtomConstraint::TetrahedralStereo(c) => Self::TetrahedralStereo(into_py_variant(
                py,
                TetrahedralStereoAst::from_ast(py, c)?,
            )?),
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstAtomConstraint {
        match self {
            Self::Valence(v) => AstAtomConstraint::Valence(v.bind(py).borrow().to_ast(py)),
            Self::TotalValence(v) => {
                AstAtomConstraint::TotalValence(v.bind(py).borrow().to_ast(py))
            }
            Self::AromaticValence(c) => {
                AstAtomConstraint::AromaticValence(c.bind(py).borrow().to_ast(py))
            }
            Self::MulticenterValence(c) => {
                AstAtomConstraint::MulticenterValence(c.bind(py).borrow().to_ast(py))
            }
            Self::DonatedPairs(v) => {
                AstAtomConstraint::DonatedPairs(v.bind(py).borrow().to_ast(py))
            }
            Self::AcceptedPairs(v) => {
                AstAtomConstraint::AcceptedPairs(v.bind(py).borrow().to_ast(py))
            }
            Self::Degree(v) => AstAtomConstraint::Degree(v.bind(py).borrow().to_ast(py)),
            Self::TotalDegree(v) => AstAtomConstraint::TotalDegree(v.bind(py).borrow().to_ast(py)),
            Self::RingDegree(v) => AstAtomConstraint::RingDegree(v.bind(py).borrow().to_ast(py)),
            Self::RingValence(v) => AstAtomConstraint::RingValence(v.bind(py).borrow().to_ast(py)),
            Self::TotalHydrogens(v) => {
                AstAtomConstraint::TotalHydrogens(v.bind(py).borrow().to_ast(py))
            }
            Self::RingMembership(m) => {
                AstAtomConstraint::RingMembership(m.bind(py).borrow().to_ast(py))
            }
            Self::TetrahedralStereo(c) => {
                AstAtomConstraint::TetrahedralStereo(c.bind(py).borrow().to_ast(py))
            }
        }
    }
}

/// The atom-scope constraints on an atom, in kind-sorted order.
#[pyclass]
pub struct AtomConstraints(AstAtomConstraints);

#[pymethods]
impl AtomConstraints {
    /// Build from a sequence of constraints (kind-sorted; a unique kind replaces
    /// an earlier one, ring memberships accumulate per scope).
    #[new]
    fn new(py: Python<'_>, entries: Vec<Py<AtomConstraint>>) -> Self {
        let mut constraints = AstAtomConstraints::new();
        constraints.extend(entries.into_iter().map(|entry| entry.bind(py).borrow().to_ast(py)));
        AtomConstraints(constraints)
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn __iter__(&self, py: Python<'_>) -> PyResult<AtomConstraintIter> {
        let entries = self
            .0
            .iter()
            .map(|constraint| into_py_variant(py, AtomConstraint::from_ast(py, constraint)?))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(AtomConstraintIter {
            entries: entries.into_iter(),
        })
    }

    /// The constraint with the given key, or `None`.
    fn get(&self, py: Python<'_>, key: Py<AtomConstraintKey>) -> PyResult<Option<AtomConstraint>> {
        self.0
            .get(key.bind(py).borrow().to_ast(py))
            .map(|constraint| AtomConstraint::from_ast(py, constraint))
            .transpose()
    }

    fn contains(&self, py: Python<'_>, key: Py<AtomConstraintKey>) -> bool {
        self.0.contains(key.bind(py).borrow().to_ast(py))
    }
}

impl AtomConstraints {
    /// The wrapped AST constraints — read access for atom construction.
    pub(crate) fn inner(&self) -> &AstAtomConstraints {
        &self.0
    }

    /// Wrap AST constraints (the hold-the-value `from_inner` bridge).
    pub(crate) fn from_inner(constraints: AstAtomConstraints) -> Self {
        AtomConstraints(constraints)
    }
}

#[pyclass]
struct AtomConstraintIter {
    entries: IntoIter<Py<AtomConstraint>>,
}

#[pymethods]
impl AtomConstraintIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<Py<AtomConstraint>> {
        self.entries.next()
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::TetrahedralStereoAst as AstTetrahedralStereoAst;

    use super::*;

    #[rstest]
    #[case(AstAromaticValenceAst::Undetermined)]
    #[case(AstAromaticValenceAst::NotAromatic)]
    #[case(AstAromaticValenceAst::aromatic(1))]
    fn test_aromatic_valence_ast_roundtrip(#[case] ast: AstAromaticValenceAst) {
        Python::attach(|py| {
            assert_eq!(
                AromaticValenceAst::from_ast(py, &ast).unwrap().to_ast(py),
                ast
            );
        });
    }

    #[rstest]
    #[case(AstMulticenterValenceAst::Undetermined)]
    #[case(AstMulticenterValenceAst::NotMulticenter)]
    #[case(AstMulticenterValenceAst::multicenter(2))]
    fn test_multicenter_valence_ast_roundtrip(#[case] ast: AstMulticenterValenceAst) {
        Python::attach(|py| {
            assert_eq!(
                MulticenterValenceAst::from_ast(py, &ast)
                    .unwrap()
                    .to_ast(py),
                ast
            );
        });
    }

    #[rstest]
    #[case(AstRingScope::All)]
    #[case(AstRingScope::Size(6))]
    fn test_ring_scope_roundtrip(#[case] ast: AstRingScope) {
        assert_eq!(RingScope::from_ast(&ast).to_ast(), ast);
    }

    #[rstest]
    #[case(AstRingMembershipAst::new(AstRingScope::All, 2))]
    #[case(AstRingMembershipAst::new(AstRingScope::Size(6), 1))]
    fn test_ring_membership_ast_roundtrip(#[case] ast: AstRingMembershipAst) {
        Python::attach(|py| {
            assert_eq!(
                RingMembershipAst::from_ast(py, &ast).unwrap().to_ast(py),
                ast
            );
        });
    }

    #[rstest]
    #[case(AstAtomConstraint::valence(4))]
    #[case(AstAtomConstraint::aromatic_valence(AstAromaticValenceAst::aromatic(1)))]
    #[case(AstAtomConstraint::ring_membership(AstRingScope::All, 2))]
    #[case(AstAtomConstraint::tetrahedral_stereo(AstTetrahedralStereoAst::not_stereo()))]
    fn test_atom_constraint_roundtrip(#[case] ast: AstAtomConstraint) {
        Python::attach(|py| {
            assert_eq!(AtomConstraint::from_ast(py, &ast).unwrap().to_ast(py), ast);
        });
    }

    #[rstest]
    fn test_atom_constraints_len_contains() {
        Python::attach(|py| {
            let valence = into_py_variant(
                py,
                AtomConstraint::from_ast(py, &AstAtomConstraint::valence(4)).unwrap(),
            )
            .unwrap();
            let degree = into_py_variant(
                py,
                AtomConstraint::from_ast(py, &AstAtomConstraint::degree(3)).unwrap(),
            )
            .unwrap();
            let constraints = AtomConstraints::new(py, vec![valence, degree]);
            assert_eq!(constraints.__len__(), 2);
            assert!(constraints.contains(
                py,
                into_py_variant(py, AtomConstraintKey::Valence()).unwrap()
            ));
            assert!(constraints.contains(
                py,
                into_py_variant(py, AtomConstraintKey::Degree()).unwrap()
            ));
            assert!(!constraints.contains(
                py,
                into_py_variant(py, AtomConstraintKey::TotalHydrogens()).unwrap()
            ));
        });
    }
}
