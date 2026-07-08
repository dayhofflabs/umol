//! Atom-constraint sub-ASTs mirroring `umol_ast::ast::constraint` (S5a): the
//! aromatic/multicenter valence states, ring scope, and ring membership. The
//! `AtomConstraintAst` enum and `AtomConstraintsAst` container follow at S5b.

use std::vec::IntoIter;

use pyo3::prelude::*;
use pyo3::types::PyDict;
use umol_ast::ast::{
    AromaticValenceAst as AstAromaticValenceAst, AtomConstraintAst as AstAtomConstraintAst,
    AtomConstraintKey as AstAtomConstraintKey, AtomConstraintsAst as AstAtomConstraintsAst,
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

/// The atom-scope constraints on an atom, in kind-sorted order.
#[pyclass]
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
            .map(|constraint| into_py_variant(py, AtomConstraintAst::from_ast(py, constraint)?))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(AtomConstraintIter {
            entries: entries.into_iter(),
        })
    }

    /// The constraint with the given key, or `None`.
    fn get(
        &self,
        py: Python<'_>,
        key: Py<AtomConstraintKey>,
    ) -> PyResult<Option<AtomConstraintAst>> {
        self.0
            .get(key.bind(py).borrow().to_ast(py))
            .map(|constraint| AtomConstraintAst::from_ast(py, constraint))
            .transpose()
    }

    fn contains(&self, py: Python<'_>, key: Py<AtomConstraintKey>) -> bool {
        self.0.contains(key.bind(py).borrow().to_ast(py))
    }

    /// The valence value, or `None`.
    fn valence(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .valence()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    /// The aromatic-valence state, or `None`.
    fn aromatic_valence(&self, py: Python<'_>) -> PyResult<Option<AromaticValenceAst>> {
        self.0
            .aromatic_valence()
            .map(|c| AromaticValenceAst::from_ast(py, c))
            .transpose()
    }

    /// The multicenter-valence state, or `None`.
    fn multicenter_valence(&self, py: Python<'_>) -> PyResult<Option<MulticenterValenceAst>> {
        self.0
            .multicenter_valence()
            .map(|c| MulticenterValenceAst::from_ast(py, c))
            .transpose()
    }

    /// The tetrahedral-stereo state, or `None`.
    fn tetrahedral_stereo(&self, py: Python<'_>) -> PyResult<Option<TetrahedralStereoAst>> {
        self.0
            .tetrahedral_stereo()
            .map(|c| TetrahedralStereoAst::from_ast(py, c))
            .transpose()
    }

    /// The degree value, or `None`.
    fn degree(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .degree()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    /// The total-degree value, or `None`.
    fn total_degree(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .total_degree()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    /// The total-valence value, or `None`.
    fn total_valence(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .total_valence()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    /// The ring-degree value, or `None`.
    fn ring_degree(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .ring_degree()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    /// The ring-valence value, or `None`.
    fn ring_valence(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .ring_valence()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    /// The total-hydrogens value, or `None`.
    fn total_hydrogens(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .total_hydrogens()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    /// The donated-pairs value, or `None`.
    fn donated_pairs(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .donated_pairs()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    /// The accepted-pairs value, or `None`.
    fn accepted_pairs(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .accepted_pairs()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    /// The all-rings membership count, or `None`.
    fn ring_count(&self, py: Python<'_>) -> PyResult<Option<ValueAst>> {
        self.0
            .ring_count()
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    /// The membership count for rings of the given size, or `None`.
    fn ring_size_count(&self, py: Python<'_>, size: u8) -> PyResult<Option<ValueAst>> {
        self.0
            .ring_size_count(size)
            .map(|v| ValueAst::from_ast(py, v))
            .transpose()
    }

    /// The present constraints as a dict keyed by snake_case name; values are the
    /// inner-value mirrors. Ring memberships key by scope: `ring_count` for the
    /// all-rings scope, `ring_size_count_<n>` for a specific ring size.
    pub(crate) fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        for entry in self.0.iter() {
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
}

impl AtomConstraintsAst {
    /// The wrapped AST constraints — read access for atom construction.
    pub(crate) fn inner(&self) -> &AstAtomConstraintsAst {
        &self.0
    }

    /// Wrap AST constraints (the hold-the-value `from_inner` bridge).
    pub(crate) fn from_inner(constraints: AstAtomConstraintsAst) -> Self {
        AtomConstraintsAst(constraints)
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

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::{TetrahedralStereoAst as AstTetrahedralStereoAst, ValueAst as AstValueAst};

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
            let constraints = AtomConstraintsAst::new(py, vec![membership]);
            assert_eq!(
                constraints
                    .ring_size_count(py, 6)
                    .unwrap()
                    .unwrap()
                    .to_ast(py),
                AstValueAst::Lit(1)
            );
            assert!(constraints.ring_size_count(py, 5).unwrap().is_none());
            assert!(constraints.ring_count(py).unwrap().is_none());
        });
    }
}
