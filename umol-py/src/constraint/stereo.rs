//! Stereo-atom and stereo-bond constraint values, containers, and live views.
#![allow(clippy::absolute_paths)]

use std::collections::BTreeSet;
use std::vec::IntoIter;

use pyo3::exceptions::{PyIndexError, PyKeyError};
use pyo3::prelude::*;
use umol_ast::ast::{
    FluxionalityAst as AstFluxionalityAst, LigandSymmetryAst as AstLigandSymmetryAst,
    StereoAtomConstraintAst as AstStereoAtomConstraintAst,
    StereoAtomConstraintKey as AstStereoAtomConstraintKey,
    StereoAtomConstraintsAst as AstStereoAtomConstraintsAst, StereoAtomId as AstStereoAtomId,
    StereoBondConstraintAst as AstStereoBondConstraintAst,
    StereoBondConstraintKey as AstStereoBondConstraintKey,
    StereoBondConstraintsAst as AstStereoBondConstraintsAst, StereoBondId as AstStereoBondId,
    StereogenicityAst as AstStereogenicityAst, TopicityAst as AstTopicityAst,
    TopicityRelationAst as AstTopicityRelationAst,
};

use crate::boolean::{BooleanArg, BooleanAst};
use crate::convert::{hash_rust, into_py_variant, variant_repr};
use crate::lattice::impl_py_lattice;
use crate::molecule::MoleculeAst;
use crate::stereo::{
    LigandPermutation, OrientedLigandPermutation, StereoAtomAst, StereoBondAst, StereoLigandPair,
    Stereogenicity, Topicity,
};

/// A topicity relation constraint value: the undetermined wildcard, a single topicity, a set
/// of admissible topicities, or the complement of a set. A finite-domain subset lattice over
/// `Topicity`. Corresponds to the Rust `TopicityRelationAst`.
#[pyclass]
pub enum TopicityRelationAst {
    Undetermined(),
    Lit(Topicity),
    LitSet(BTreeSet<Topicity>),
    NotSet(BTreeSet<Topicity>),
}

#[pymethods]
impl TopicityRelationAst {
    /// The single topicity this resolves to, or `None` when it is not a bare literal.
    pub(crate) fn as_lit(&self) -> Option<Topicity> {
        match self {
            Self::Lit(topicity) => Some(*topicity),
            _ => None,
        }
    }

    pub(crate) fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    pub(crate) fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    pub(crate) fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            TopicityRelationAst::Undetermined() => ("Undetermined", 0),
            TopicityRelationAst::Lit(_) => ("Lit", 1),
            TopicityRelationAst::LitSet(_) => ("LitSet", 1),
            TopicityRelationAst::NotSet(_) => ("NotSet", 1),
        };
        variant_repr(slf.bind(py).as_any(), "TopicityRelationAst", variant, arity)
    }
}

impl_py_lattice!(
    TopicityRelationAst,
    AstTopicityRelationAst,
    |value: &TopicityRelationAst, _py: Python<'_>| -> PyResult<AstTopicityRelationAst> {
        Ok(value.to_rust())
    },
    |_py: Python<'_>, value: AstTopicityRelationAst| -> PyResult<TopicityRelationAst> {
        Ok(TopicityRelationAst::from_rust(&value))
    }
);

impl TopicityRelationAst {
    pub(crate) fn from_rust(ast: &AstTopicityRelationAst) -> Self {
        match ast {
            AstTopicityRelationAst::Undetermined => Self::Undetermined(),
            AstTopicityRelationAst::Lit(topicity) => Self::Lit(Topicity::from_rust(*topicity)),
            AstTopicityRelationAst::LitSet(topicities) => {
                Self::LitSet(topicities.iter().map(|t| Topicity::from_rust(*t)).collect())
            }
            AstTopicityRelationAst::NotSet(topicities) => {
                Self::NotSet(topicities.iter().map(|t| Topicity::from_rust(*t)).collect())
            }
        }
    }

    pub(crate) fn to_rust(&self) -> AstTopicityRelationAst {
        match self {
            Self::Undetermined() => AstTopicityRelationAst::Undetermined,
            Self::Lit(topicity) => AstTopicityRelationAst::Lit(topicity.to_rust()),
            Self::LitSet(topicities) => {
                AstTopicityRelationAst::LitSet(topicities.iter().map(|t| t.to_rust()).collect())
            }
            Self::NotSet(topicities) => {
                AstTopicityRelationAst::NotSet(topicities.iter().map(|t| t.to_rust()).collect())
            }
        }
    }
}

/// Setter coercion for a topicity relation: a `Topicity` literal (→ `Lit`) or a
/// `TopicityRelationAst` passthrough (matching `impl From<Topicity>`).
#[derive(FromPyObject)]
pub(crate) enum TopicityRelationArg {
    Lit(Topicity),
    Ast(Py<TopicityRelationAst>),
}

impl TopicityRelationArg {
    /// Coerce to a `Py<TopicityRelationAst>` (for the `TopicityAst.relation` field).
    pub(crate) fn to_py(&self, py: Python<'_>) -> PyResult<Py<TopicityRelationAst>> {
        match self {
            TopicityRelationArg::Lit(topicity) => {
                into_py_variant(py, TopicityRelationAst::Lit(*topicity))
            }
            TopicityRelationArg::Ast(relation) => Ok(relation.clone_ref(py)),
        }
    }
}

/// A stereogenicity constraint value: the undetermined wildcard, a single classification, a
/// set of admissible classifications, or the complement of a set. A finite-domain subset
/// lattice over `Stereogenicity`. Corresponds to the Rust `StereogenicityAst`.
#[pyclass]
pub enum StereogenicityAst {
    Undetermined(),
    Lit(Stereogenicity),
    LitSet(BTreeSet<Stereogenicity>),
    NotSet(BTreeSet<Stereogenicity>),
}

#[pymethods]
impl StereogenicityAst {
    /// The single classification this resolves to, or `None` when it is not a bare literal.
    pub(crate) fn as_lit(&self) -> Option<Stereogenicity> {
        match self {
            Self::Lit(stereogenicity) => Some(*stereogenicity),
            _ => None,
        }
    }

    pub(crate) fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    pub(crate) fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    pub(crate) fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            StereogenicityAst::Undetermined() => ("Undetermined", 0),
            StereogenicityAst::Lit(_) => ("Lit", 1),
            StereogenicityAst::LitSet(_) => ("LitSet", 1),
            StereogenicityAst::NotSet(_) => ("NotSet", 1),
        };
        variant_repr(slf.bind(py).as_any(), "StereogenicityAst", variant, arity)
    }
}

impl_py_lattice!(
    StereogenicityAst,
    AstStereogenicityAst,
    |value: &StereogenicityAst, _py: Python<'_>| -> PyResult<AstStereogenicityAst> {
        Ok(value.to_rust())
    },
    |_py: Python<'_>, value: AstStereogenicityAst| -> PyResult<StereogenicityAst> {
        Ok(StereogenicityAst::from_rust(&value))
    }
);

impl StereogenicityAst {
    pub(crate) fn from_rust(ast: &AstStereogenicityAst) -> Self {
        match ast {
            AstStereogenicityAst::Undetermined => Self::Undetermined(),
            AstStereogenicityAst::Lit(stereogenicity) => {
                Self::Lit(Stereogenicity::from_rust(*stereogenicity))
            }
            AstStereogenicityAst::LitSet(stereogenicities) => Self::LitSet(
                stereogenicities
                    .iter()
                    .map(|g| Stereogenicity::from_rust(*g))
                    .collect(),
            ),
            AstStereogenicityAst::NotSet(stereogenicities) => Self::NotSet(
                stereogenicities
                    .iter()
                    .map(|g| Stereogenicity::from_rust(*g))
                    .collect(),
            ),
        }
    }

    pub(crate) fn to_rust(&self) -> AstStereogenicityAst {
        match self {
            Self::Undetermined() => AstStereogenicityAst::Undetermined,
            Self::Lit(stereogenicity) => AstStereogenicityAst::Lit(stereogenicity.to_rust()),
            Self::LitSet(stereogenicities) => {
                AstStereogenicityAst::LitSet(stereogenicities.iter().map(|g| g.to_rust()).collect())
            }
            Self::NotSet(stereogenicities) => {
                AstStereogenicityAst::NotSet(stereogenicities.iter().map(|g| g.to_rust()).collect())
            }
        }
    }
}

/// A ligand-symmetry constraint value: an oriented ligand permutation with a presence
/// assertion (whether the permutation is a ligand symmetry). Corresponds to the Rust
/// `LigandSymmetryAst`.
#[pyclass]
pub struct LigandSymmetryAst {
    pub(crate) permutation: OrientedLigandPermutation,
    pub(crate) invariant: Py<BooleanAst>,
}

#[pymethods]
impl LigandSymmetryAst {
    #[new]
    pub(crate) fn new(
        py: Python<'_>,
        permutation: OrientedLigandPermutation,
        invariant: BooleanArg,
    ) -> PyResult<Self> {
        Ok(LigandSymmetryAst {
            permutation,
            invariant: into_py_variant(py, BooleanAst::from_rust(&invariant.to_rust(py)))?,
        })
    }

    #[getter]
    pub(crate) fn permutation(&self) -> OrientedLigandPermutation {
        self.permutation
    }

    #[getter]
    pub(crate) fn invariant(&self, py: Python<'_>) -> Py<BooleanAst> {
        self.invariant.clone_ref(py)
    }

    pub(crate) fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    pub(crate) fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
    }

    pub(crate) fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "LigandSymmetryAst({}, {})",
            self.permutation.__repr__(),
            self.invariant
                .bind(py)
                .as_any()
                .repr()?
                .extract::<String>()?,
        ))
    }
}

impl_py_lattice!(
    LigandSymmetryAst,
    AstLigandSymmetryAst,
    |value: &LigandSymmetryAst, py: Python<'_>| -> PyResult<AstLigandSymmetryAst> {
        Ok(value.to_rust(py))
    },
    |py: Python<'_>, value: AstLigandSymmetryAst| -> PyResult<LigandSymmetryAst> {
        LigandSymmetryAst::from_rust(py, &value)
    }
);

impl LigandSymmetryAst {
    pub(crate) fn from_rust(py: Python<'_>, ast: &AstLigandSymmetryAst) -> PyResult<Self> {
        Ok(LigandSymmetryAst {
            permutation: OrientedLigandPermutation::from_rust(ast.permutation),
            invariant: into_py_variant(py, BooleanAst::from_rust(&ast.invariant))?,
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstLigandSymmetryAst {
        AstLigandSymmetryAst {
            permutation: self.permutation.to_rust(),
            invariant: self.invariant.bind(py).borrow().to_rust(),
        }
    }
}

/// A fluxionality constraint value: a proper ligand permutation realized by dynamics, with an
/// assertion of whether the move is `active`. Corresponds to the Rust `FluxionalityAst`.
#[pyclass]
pub struct FluxionalityAst {
    pub(crate) permutation: LigandPermutation,
    pub(crate) active: Py<BooleanAst>,
}

#[pymethods]
impl FluxionalityAst {
    #[new]
    pub(crate) fn new(
        py: Python<'_>,
        permutation: LigandPermutation,
        active: BooleanArg,
    ) -> PyResult<Self> {
        Ok(FluxionalityAst {
            permutation,
            active: into_py_variant(py, BooleanAst::from_rust(&active.to_rust(py)))?,
        })
    }

    #[getter]
    pub(crate) fn permutation(&self) -> LigandPermutation {
        self.permutation
    }

    #[getter]
    pub(crate) fn active(&self, py: Python<'_>) -> Py<BooleanAst> {
        self.active.clone_ref(py)
    }

    pub(crate) fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    pub(crate) fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
    }

    pub(crate) fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "FluxionalityAst({}, {})",
            self.permutation.__repr__(),
            self.active.bind(py).as_any().repr()?.extract::<String>()?,
        ))
    }
}

impl_py_lattice!(
    FluxionalityAst,
    AstFluxionalityAst,
    |value: &FluxionalityAst, py: Python<'_>| -> PyResult<AstFluxionalityAst> {
        Ok(value.to_rust(py))
    },
    |py: Python<'_>, value: AstFluxionalityAst| -> PyResult<FluxionalityAst> {
        FluxionalityAst::from_rust(py, &value)
    }
);

impl FluxionalityAst {
    pub(crate) fn from_rust(py: Python<'_>, ast: &AstFluxionalityAst) -> PyResult<Self> {
        Ok(FluxionalityAst {
            permutation: LigandPermutation::from_rust(ast.permutation),
            active: into_py_variant(py, BooleanAst::from_rust(&ast.active))?,
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstFluxionalityAst {
        AstFluxionalityAst {
            permutation: self.permutation.to_rust(),
            active: self.active.bind(py).borrow().to_rust(),
        }
    }
}

/// A per-pair topicity constraint value: a relation between a pair of ligand positions.
/// Corresponds to the Rust `TopicityAst`.
#[pyclass]
pub struct TopicityAst {
    pub(crate) pair: StereoLigandPair,
    pub(crate) relation: Py<TopicityRelationAst>,
}

#[pymethods]
impl TopicityAst {
    #[new]
    pub(crate) fn new(
        py: Python<'_>,
        pair: StereoLigandPair,
        relation: TopicityRelationArg,
    ) -> PyResult<Self> {
        Ok(TopicityAst {
            pair,
            relation: relation.to_py(py)?,
        })
    }

    #[getter]
    pub(crate) fn pair(&self) -> StereoLigandPair {
        self.pair
    }

    #[getter]
    pub(crate) fn relation(&self, py: Python<'_>) -> Py<TopicityRelationAst> {
        self.relation.clone_ref(py)
    }

    pub(crate) fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    pub(crate) fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
    }

    pub(crate) fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "TopicityAst({}, {})",
            self.pair.__repr__(),
            self.relation
                .bind(py)
                .as_any()
                .repr()?
                .extract::<String>()?,
        ))
    }
}

impl_py_lattice!(
    TopicityAst,
    AstTopicityAst,
    |value: &TopicityAst, py: Python<'_>| -> PyResult<AstTopicityAst> { Ok(value.to_rust(py)) },
    |py: Python<'_>, value: AstTopicityAst| -> PyResult<TopicityAst> {
        TopicityAst::from_rust(py, &value)
    }
);

impl TopicityAst {
    pub(crate) fn from_rust(py: Python<'_>, ast: &AstTopicityAst) -> PyResult<Self> {
        Ok(TopicityAst {
            pair: StereoLigandPair::from_rust(ast.pair),
            relation: into_py_variant(py, TopicityRelationAst::from_rust(&ast.relation))?,
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstTopicityAst {
        AstTopicityAst {
            pair: self.pair.to_rust(),
            relation: self.relation.bind(py).borrow().to_rust(),
        }
    }
}

/// Per-entity stereo constraint surface — key + constraint enum + container + args —
/// macro-generated for the two stereo entities (`StereoAtom`, `StereoBond`), which share the
/// value types (`LigandSymmetryAst`/`FluxionalityAst`/`TopicityAst`/`StereogenicityAst`) and
/// key sub-types (`OrientedLigandPermutation`/`LigandPermutation`/`StereoLigandPair`); only
/// the enum/container/key names and their AST peers differ.
macro_rules! stereo_constraints {
    (
        $key:ident, $constraint:ident, $constraints:ident,
        $update:ident, $resolved:ident, $arg:ident,
        $key_iter:ident, $iter:ident, $items_iter:ident,
        $ast_key:ident, $ast_constraint:ident, $ast_constraints:ident,
        $value:ident, $view:ident, $backing:ident,
        $ast_id:ident, $namespace:ident, $entity_mut:ident, $id_error:literal $(,)?
    ) => {
        /// The key (identity) of a stereo constraint: the sub-keyed oriented/ligand
        /// permutation or ligand pair for the per-permutation / per-pair constraints; the
        /// bare discriminant for stereogenicity.
        #[pyclass]
        pub enum $key {
            LigandSymmetry(Py<OrientedLigandPermutation>),
            Fluxionality(Py<LigandPermutation>),
            Topicity(Py<StereoLigandPair>),
            Stereogenicity(),
        }

        #[pymethods]
        impl $key {
            pub(crate) fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
                self.to_rust(py) == other.to_rust(py)
            }

            pub(crate) fn __hash__(&self, py: Python<'_>) -> u64 {
                hash_rust(&self.to_rust(py))
            }

            pub(crate) fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
                let (variant, arity) = match &*slf.bind(py).borrow() {
                    $key::LigandSymmetry(_) => ("LigandSymmetry", 1),
                    $key::Fluxionality(_) => ("Fluxionality", 1),
                    $key::Topicity(_) => ("Topicity", 1),
                    $key::Stereogenicity() => ("Stereogenicity", 0),
                };
                variant_repr(slf.bind(py).as_any(), stringify!($key), variant, arity)
            }
        }

        impl $key {
            pub(crate) fn from_rust(py: Python<'_>, ast: &$ast_key) -> PyResult<Self> {
                Ok(match ast {
                    $ast_key::LigandSymmetry(permutation) => Self::LigandSymmetry(into_py_variant(
                        py,
                        OrientedLigandPermutation::from_rust(*permutation),
                    )?),
                    $ast_key::Fluxionality(permutation) => Self::Fluxionality(into_py_variant(
                        py,
                        LigandPermutation::from_rust(*permutation),
                    )?),
                    $ast_key::Topicity(pair) => {
                        Self::Topicity(into_py_variant(py, StereoLigandPair::from_rust(*pair))?)
                    }
                    $ast_key::Stereogenicity => Self::Stereogenicity(),
                })
            }

            pub(crate) fn to_rust(&self, py: Python<'_>) -> $ast_key {
                match self {
                    Self::LigandSymmetry(permutation) => {
                        $ast_key::LigandSymmetry(permutation.bind(py).borrow().to_rust())
                    }
                    Self::Fluxionality(permutation) => {
                        $ast_key::Fluxionality(permutation.bind(py).borrow().to_rust())
                    }
                    Self::Topicity(pair) => $ast_key::Topicity(pair.bind(py).borrow().to_rust()),
                    Self::Stereogenicity() => $ast_key::Stereogenicity,
                }
            }
        }

        /// A stereo constraint: a ligand-symmetry, fluxionality, topicity, or stereogenicity
        /// predicate on a stereo atom / bond.
        #[pyclass]
        pub enum $constraint {
            LigandSymmetry(Py<LigandSymmetryAst>),
            Fluxionality(Py<FluxionalityAst>),
            Topicity(Py<TopicityAst>),
            Stereogenicity(Py<StereogenicityAst>),
        }

        #[pymethods]
        impl $constraint {
            /// The constraint's key (identity).
            #[getter]
            pub(crate) fn key(&self, py: Python<'_>) -> PyResult<$key> {
                $key::from_rust(py, &self.to_rust(py).key())
            }

            pub(crate) fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
                self.to_rust(py) == other.to_rust(py)
            }

            pub(crate) fn __hash__(&self, py: Python<'_>) -> u64 {
                hash_rust(&self.to_rust(py))
            }

            pub(crate) fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
                let variant = match &*slf.bind(py).borrow() {
                    $constraint::LigandSymmetry(_) => "LigandSymmetry",
                    $constraint::Fluxionality(_) => "Fluxionality",
                    $constraint::Topicity(_) => "Topicity",
                    $constraint::Stereogenicity(_) => "Stereogenicity",
                };
                variant_repr(slf.bind(py).as_any(), stringify!($constraint), variant, 1)
            }
        }

        impl_py_lattice!(
            $constraint,
            $ast_constraint,
            |value: &$constraint, py: Python<'_>| -> PyResult<$ast_constraint> {
                Ok(value.to_rust(py))
            },
            |py: Python<'_>, value: $ast_constraint| -> PyResult<$constraint> {
                $constraint::from_rust(py, &value)
            }
        );

        impl $constraint {
            pub(crate) fn from_rust(py: Python<'_>, ast: &$ast_constraint) -> PyResult<Self> {
                Ok(match ast {
                    $ast_constraint::LigandSymmetry(value) => Self::LigandSymmetry(
                        into_py_variant(py, LigandSymmetryAst::from_rust(py, value)?)?,
                    ),
                    $ast_constraint::Fluxionality(value) => Self::Fluxionality(into_py_variant(
                        py,
                        FluxionalityAst::from_rust(py, value)?,
                    )?),
                    $ast_constraint::Topicity(value) => {
                        Self::Topicity(into_py_variant(py, TopicityAst::from_rust(py, value)?)?)
                    }
                    $ast_constraint::Stereogenicity(value) => Self::Stereogenicity(
                        into_py_variant(py, StereogenicityAst::from_rust(value))?,
                    ),
                })
            }

            pub(crate) fn to_rust(&self, py: Python<'_>) -> $ast_constraint {
                match self {
                    Self::LigandSymmetry(value) => {
                        $ast_constraint::LigandSymmetry(value.bind(py).borrow().to_rust(py))
                    }
                    Self::Fluxionality(value) => {
                        $ast_constraint::Fluxionality(value.bind(py).borrow().to_rust(py))
                    }
                    Self::Topicity(value) => {
                        $ast_constraint::Topicity(value.bind(py).borrow().to_rust(py))
                    }
                    Self::Stereogenicity(value) => {
                        $ast_constraint::Stereogenicity(value.bind(py).borrow().to_rust())
                    }
                }
            }
        }

        /// Argument to the container's `update`: another container, a live view, or a loose
        /// iterable of constraints.
        #[derive(FromPyObject)]
        pub(crate) enum $update {
            Container(Py<$constraints>),
            View(Py<$view>),
            Entries(Vec<Py<$constraint>>),
        }

        impl $update {
            /// Read every Python object into owned data before any write borrow is taken, so a
            /// container or view that aliases the same entity is read while nothing is borrowed
            /// (otherwise `cs.update(cs)` self-aliases into a double-borrow panic).
            pub(crate) fn resolve(&self, py: Python<'_>) -> PyResult<$resolved> {
                Ok(match self {
                    $update::Container(c) => {
                        $resolved::Overlay(c.bind(py).borrow().inner().clone())
                    }
                    $update::View(v) => {
                        $resolved::Overlay(v.bind(py).borrow().read(py, |cs| Ok(cs.clone()))?)
                    }
                    $update::Entries(entries) => $resolved::Entries(
                        entries
                            .iter()
                            .map(|entry| entry.bind(py).borrow().to_rust(py))
                            .collect(),
                    ),
                })
            }
        }

        /// A `$update` with all Python reads done, applicable under a write borrow.
        pub(crate) enum $resolved {
            Overlay($ast_constraints),
            Entries(Vec<$ast_constraint>),
        }

        impl $resolved {
            pub(crate) fn apply(self, target: &mut $ast_constraints) {
                match self {
                    $resolved::Overlay(overlay) => target.update(&overlay),
                    $resolved::Entries(entries) => {
                        for entry in entries {
                            target.set(entry);
                        }
                    }
                }
            }
        }

        /// A whole-container argument for the entity `constraints` setter: a value container
        /// or a live view (which is read while unborrowed, self-alias safe).
        #[derive(FromPyObject)]
        pub(crate) enum $arg {
            Container(Py<$constraints>),
            View(Py<$view>),
        }

        impl $arg {
            pub(crate) fn to_rust(&self, py: Python<'_>) -> PyResult<$ast_constraints> {
                match self {
                    $arg::Container(c) => Ok(c.bind(py).borrow().inner().clone()),
                    $arg::View(v) => v.bind(py).borrow().read(py, |cs| Ok(cs.clone())),
                }
            }
        }

        /// The stereo constraints on a stereo atom / bond, in kind-sorted order. Mutable,
        /// hence value-equal but unhashable.
        #[pyclass(eq)]
        #[derive(PartialEq)]
        pub struct $constraints($ast_constraints);

        #[pymethods]
        impl $constraints {
            /// Build from a sequence of constraints (kind-sorted; a unique key replaces an
            /// earlier one; per-permutation / per-pair entries accumulate).
            #[new]
            pub(crate) fn new(py: Python<'_>, entries: Vec<Py<$constraint>>) -> Self {
                let mut constraints = $ast_constraints::new();
                constraints.extend(
                    entries
                        .into_iter()
                        .map(|entry| entry.bind(py).borrow().to_rust(py)),
                );
                $constraints(constraints)
            }

            pub(crate) fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
                let mut parts = Vec::with_capacity(self.0.len());
                for entry in self.0.iter() {
                    let value = into_py_variant(py, $constraint::from_rust(py, entry)?)?;
                    parts.push(value.bind(py).as_any().repr()?.extract::<String>()?);
                }
                Ok(format!(
                    "{}([{}])",
                    stringify!($constraints),
                    parts.join(", ")
                ))
            }

            /// Insert `c`, replacing any existing entry of the same key (last-wins).
            pub(crate) fn set(&mut self, py: Python<'_>, c: Py<$constraint>) {
                self.0.set(c.bind(py).borrow().to_rust(py));
            }

            /// Remove the entry with the given key, returning it if present (dict `pop`).
            pub(crate) fn pop(
                &mut self,
                py: Python<'_>,
                key: Py<$key>,
            ) -> PyResult<Option<$constraint>> {
                self.0
                    .remove(key.bind(py).borrow().to_rust(py))
                    .map(|c| $constraint::from_rust(py, &c))
                    .transpose()
            }

            /// Overlay `other` onto self in place — another container or an iterable of
            /// constraints (last-wins per key; undetermined entries remove). Takes `slf` by
            /// handle so `other` is fully read before the write borrow (`cs.update(cs)` is a
            /// no-op, not a double-borrow panic).
            pub(crate) fn update(slf: Py<Self>, py: Python<'_>, other: $update) -> PyResult<()> {
                let resolved = other.resolve(py)?;
                resolved.apply(&mut slf.borrow_mut(py).0);
                Ok(())
            }

            pub(crate) fn __len__(&self) -> usize {
                self.0.len()
            }

            /// Iterate the constraint keys (mapping-style, canonical order).
            pub(crate) fn __iter__(&self, py: Python<'_>) -> PyResult<$key_iter> {
                self.keys(py)
            }

            /// The constraint keys, in canonical order.
            pub(crate) fn keys(&self, py: Python<'_>) -> PyResult<$key_iter> {
                let keys = self
                    .0
                    .iter()
                    .map(|c| into_py_variant(py, $key::from_rust(py, &c.key())?))
                    .collect::<PyResult<Vec<_>>>()?;
                Ok($key_iter {
                    keys: keys.into_iter(),
                })
            }

            /// The constraints, in canonical order.
            pub(crate) fn values(&self, py: Python<'_>) -> PyResult<$iter> {
                let entries = self
                    .0
                    .iter()
                    .map(|c| into_py_variant(py, $constraint::from_rust(py, c)?))
                    .collect::<PyResult<Vec<_>>>()?;
                Ok($iter {
                    entries: entries.into_iter(),
                })
            }

            /// The `(key, constraint)` pairs, in canonical order.
            pub(crate) fn items(&self, py: Python<'_>) -> PyResult<$items_iter> {
                let items = self
                    .0
                    .iter()
                    .map(|c| {
                        Ok((
                            into_py_variant(py, $key::from_rust(py, &c.key())?)?,
                            into_py_variant(py, $constraint::from_rust(py, c)?)?,
                        ))
                    })
                    .collect::<PyResult<Vec<_>>>()?;
                Ok($items_iter {
                    items: items.into_iter(),
                })
            }

            /// The constraint with the given key, or `default` (`None`) if absent.
            #[pyo3(signature = (key, default=None))]
            pub(crate) fn get(
                &self,
                py: Python<'_>,
                key: Py<$key>,
                default: Option<Py<PyAny>>,
            ) -> PyResult<Py<PyAny>> {
                match self.0.get(key.bind(py).borrow().to_rust(py)) {
                    Some(constraint) => Ok(into_py_variant(
                        py,
                        $constraint::from_rust(py, constraint)?,
                    )?
                    .into_any()),
                    None => Ok(default.unwrap_or_else(|| py.None())),
                }
            }

            /// The constraint with the given key; raises `KeyError` if absent.
            pub(crate) fn __getitem__(
                &self,
                py: Python<'_>,
                key: Py<$key>,
            ) -> PyResult<$constraint> {
                match self.0.get(key.bind(py).borrow().to_rust(py)) {
                    Some(constraint) => $constraint::from_rust(py, constraint),
                    None => Err(PyKeyError::new_err(
                        key.bind(py).as_any().repr()?.extract::<String>()?,
                    )),
                }
            }

            /// Remove the entry with the given key; raises `KeyError` if absent.
            pub(crate) fn __delitem__(&mut self, py: Python<'_>, key: Py<$key>) -> PyResult<()> {
                if self.0.remove(key.bind(py).borrow().to_rust(py)).is_some() {
                    Ok(())
                } else {
                    Err(PyKeyError::new_err(
                        key.bind(py).as_any().repr()?.extract::<String>()?,
                    ))
                }
            }

            pub(crate) fn __contains__(&self, py: Python<'_>, key: Py<$key>) -> bool {
                self.0.contains(key.bind(py).borrow().to_rust(py))
            }

            /// The ligand-symmetry constraints.
            pub(crate) fn ligand_symmetries(
                &self,
                py: Python<'_>,
            ) -> PyResult<Vec<LigandSymmetryAst>> {
                self.0
                    .ligand_symmetries()
                    .map(|ls| LigandSymmetryAst::from_rust(py, ls))
                    .collect()
            }

            /// The ligand-symmetry constraint at `permutation` (undetermined if absent).
            pub(crate) fn ligand_symmetry(
                &self,
                py: Python<'_>,
                permutation: OrientedLigandPermutation,
            ) -> PyResult<LigandSymmetryAst> {
                LigandSymmetryAst::from_rust(py, &self.0.ligand_symmetry(permutation.to_rust()))
            }

            /// The fluxionality constraints.
            pub(crate) fn fluxionalities(&self, py: Python<'_>) -> PyResult<Vec<FluxionalityAst>> {
                self.0
                    .fluxionalities()
                    .map(|f| FluxionalityAst::from_rust(py, f))
                    .collect()
            }

            /// The fluxionality constraint at `permutation` (undetermined if absent).
            pub(crate) fn fluxionality(
                &self,
                py: Python<'_>,
                permutation: LigandPermutation,
            ) -> PyResult<FluxionalityAst> {
                FluxionalityAst::from_rust(py, &self.0.fluxionality(permutation.to_rust()))
            }

            /// The topicity constraints.
            pub(crate) fn topicities(&self, py: Python<'_>) -> PyResult<Vec<TopicityAst>> {
                self.0
                    .topicities()
                    .map(|t| TopicityAst::from_rust(py, t))
                    .collect()
            }

            /// The topicity relation at ligand `pair` (undetermined if absent).
            pub(crate) fn topicity(&self, pair: StereoLigandPair) -> TopicityRelationAst {
                TopicityRelationAst::from_rust(&self.0.topicity(pair.to_rust()))
            }

            /// The stereogenicity constraint (undetermined if absent).
            pub(crate) fn stereogenicity(&self) -> StereogenicityAst {
                StereogenicityAst::from_rust(&self.0.stereogenicity())
            }
        }

        impl $constraints {
            pub(crate) fn inner(&self) -> &$ast_constraints {
                &self.0
            }

            #[cfg(test)]
            pub(crate) fn from_inner(constraints: $ast_constraints) -> Self {
                $constraints(constraints)
            }
        }

        impl_py_lattice!(
            $constraints,
            $ast_constraints,
            |value: &$constraints, _py: Python<'_>| -> PyResult<$ast_constraints> {
                Ok(value.inner().clone())
            },
            |_py: Python<'_>, value: $ast_constraints| -> PyResult<$constraints> {
                Ok($constraints(value))
            }
        );

        #[pyclass]
        pub struct $key_iter {
            pub(crate) keys: IntoIter<Py<$key>>,
        }

        #[pymethods]
        impl $key_iter {
            pub(crate) fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
                slf
            }

            pub(crate) fn __next__(&mut self) -> Option<Py<$key>> {
                self.keys.next()
            }
        }

        #[pyclass]
        pub struct $iter {
            pub(crate) entries: IntoIter<Py<$constraint>>,
        }

        #[pymethods]
        impl $iter {
            pub(crate) fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
                slf
            }

            pub(crate) fn __next__(&mut self) -> Option<Py<$constraint>> {
                self.entries.next()
            }
        }

        #[pyclass]
        pub struct $items_iter {
            pub(crate) items: IntoIter<(Py<$key>, Py<$constraint>)>,
        }

        #[pymethods]
        impl $items_iter {
            pub(crate) fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
                slf
            }

            pub(crate) fn __next__(&mut self) -> Option<(Py<$key>, Py<$constraint>)> {
                self.items.next()
            }
        }

        /// What a `$view` writes through to: a stereo entity within a molecule (by id) or a
        /// standalone own-value stereo entity (`Py<$value>`).
        pub(crate) enum $backing {
            Molecule { owner: Py<MoleculeAst>, id: $ast_id },
            Value(Py<$value>),
        }

        /// A live handle onto one stereo entity's constraints, backed by either a
        /// molecule-embedded entity or a standalone value. Reads borrow the entity and read
        /// only what they need; mutators write through in place, without a clone-and-writeback.
        #[pyclass]
        pub struct $view {
            pub(crate) backing: $backing,
        }

        impl $view {
            /// Borrow the backing entity's constraints and read through `f` — no clone.
            pub(crate) fn read<R>(
                &self,
                py: Python<'_>,
                f: impl FnOnce(&$ast_constraints) -> PyResult<R>,
            ) -> PyResult<R> {
                match &self.backing {
                    $backing::Molecule { owner, id } => {
                        let molecule = owner.bind(py).borrow();
                        let view = molecule
                            .inner()
                            .$namespace()
                            .get(*id)
                            .ok_or_else(|| PyIndexError::new_err($id_error))?;
                        f(&view.ast.constraints)
                    }
                    $backing::Value(entity) => {
                        let entity = entity.bind(py).borrow();
                        f(&entity.inner().constraints)
                    }
                }
            }

            /// Mutate the backing entity's constraints in place through `f`.
            pub(crate) fn with_mut<R>(
                &self,
                py: Python<'_>,
                f: impl FnOnce(&mut $ast_constraints) -> R,
            ) -> R {
                match &self.backing {
                    $backing::Molecule { owner, id } => f(&mut owner
                        .borrow_mut(py)
                        .inner_mut()
                        .$entity_mut(*id)
                        .ast
                        .constraints),
                    $backing::Value(entity) => {
                        f(&mut entity.borrow_mut(py).inner_mut().constraints)
                    }
                }
            }
        }

        #[pymethods]
        impl $view {
            pub(crate) fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
                let count = self.read(py, |cs| Ok(cs.len()))?;
                Ok(format!("{}({count} entries)", stringify!($view)))
            }

            /// Insert `c` on the entity in place, replacing any existing entry of the same key
            /// (last-wins).
            pub(crate) fn set(&self, py: Python<'_>, c: Py<$constraint>) {
                let constraint = c.bind(py).borrow().to_rust(py);
                self.with_mut(py, |cs| cs.set(constraint));
            }

            /// Remove the entry with the given key, returning it if present (dict `pop`).
            pub(crate) fn pop(
                &self,
                py: Python<'_>,
                key: Py<$key>,
            ) -> PyResult<Option<$constraint>> {
                let ast_key = key.bind(py).borrow().to_rust(py);
                self.with_mut(py, |cs| cs.remove(ast_key))
                    .map(|c| $constraint::from_rust(py, &c))
                    .transpose()
            }

            /// Remove the entry with the given key; raises `KeyError` if absent.
            pub(crate) fn __delitem__(&self, py: Python<'_>, key: Py<$key>) -> PyResult<()> {
                let ast_key = key.bind(py).borrow().to_rust(py);
                if self.with_mut(py, |cs| cs.remove(ast_key)).is_some() {
                    Ok(())
                } else {
                    Err(PyKeyError::new_err(
                        key.bind(py).as_any().repr()?.extract::<String>()?,
                    ))
                }
            }

            /// Overlay `other` onto the entity's constraints in place — another container, a
            /// live view, or an iterable of constraints (last-wins per key; undetermined
            /// entries remove). Resolves `other` before the write borrow (self-alias safe).
            pub(crate) fn update(&self, py: Python<'_>, other: $update) -> PyResult<()> {
                let resolved = other.resolve(py)?;
                self.with_mut(py, |cs| resolved.apply(cs));
                Ok(())
            }

            pub(crate) fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
                self.read(py, |cs| Ok(cs.len()))
            }

            /// Iterate the constraint keys (mapping-style, canonical order).
            pub(crate) fn __iter__(&self, py: Python<'_>) -> PyResult<$key_iter> {
                self.keys(py)
            }

            /// The constraint keys, in canonical order.
            pub(crate) fn keys(&self, py: Python<'_>) -> PyResult<$key_iter> {
                let keys = self.read(py, |cs| {
                    cs.iter()
                        .map(|c| into_py_variant(py, $key::from_rust(py, &c.key())?))
                        .collect::<PyResult<Vec<_>>>()
                })?;
                Ok($key_iter {
                    keys: keys.into_iter(),
                })
            }

            /// The constraints, in canonical order.
            pub(crate) fn values(&self, py: Python<'_>) -> PyResult<$iter> {
                let entries = self.read(py, |cs| {
                    cs.iter()
                        .map(|c| into_py_variant(py, $constraint::from_rust(py, c)?))
                        .collect::<PyResult<Vec<_>>>()
                })?;
                Ok($iter {
                    entries: entries.into_iter(),
                })
            }

            /// The `(key, constraint)` pairs, in canonical order.
            pub(crate) fn items(&self, py: Python<'_>) -> PyResult<$items_iter> {
                let items = self.read(py, |cs| {
                    cs.iter()
                        .map(|c| {
                            Ok((
                                into_py_variant(py, $key::from_rust(py, &c.key())?)?,
                                into_py_variant(py, $constraint::from_rust(py, c)?)?,
                            ))
                        })
                        .collect::<PyResult<Vec<_>>>()
                })?;
                Ok($items_iter {
                    items: items.into_iter(),
                })
            }

            /// The constraint with the given key, or `default` (`None`) if absent.
            #[pyo3(signature = (key, default=None))]
            pub(crate) fn get(
                &self,
                py: Python<'_>,
                key: Py<$key>,
                default: Option<Py<PyAny>>,
            ) -> PyResult<Py<PyAny>> {
                let ast_key = key.bind(py).borrow().to_rust(py);
                let found = self.read(py, |cs| {
                    cs.get(ast_key)
                        .map(|constraint| $constraint::from_rust(py, constraint))
                        .transpose()
                })?;
                match found {
                    Some(constraint) => Ok(into_py_variant(py, constraint)?.into_any()),
                    None => Ok(default.unwrap_or_else(|| py.None())),
                }
            }

            /// The constraint with the given key; raises `KeyError` if absent.
            pub(crate) fn __getitem__(
                &self,
                py: Python<'_>,
                key: Py<$key>,
            ) -> PyResult<$constraint> {
                let ast_key = key.bind(py).borrow().to_rust(py);
                let found = self.read(py, |cs| {
                    cs.get(ast_key)
                        .map(|constraint| $constraint::from_rust(py, constraint))
                        .transpose()
                })?;
                match found {
                    Some(constraint) => Ok(constraint),
                    None => Err(PyKeyError::new_err(
                        key.bind(py).as_any().repr()?.extract::<String>()?,
                    )),
                }
            }

            pub(crate) fn __contains__(&self, py: Python<'_>, key: Py<$key>) -> PyResult<bool> {
                let ast_key = key.bind(py).borrow().to_rust(py);
                self.read(py, |cs| Ok(cs.contains(ast_key)))
            }

            /// The ligand-symmetry constraints.
            pub(crate) fn ligand_symmetries(
                &self,
                py: Python<'_>,
            ) -> PyResult<Vec<LigandSymmetryAst>> {
                self.read(py, |cs| {
                    cs.ligand_symmetries()
                        .map(|ls| LigandSymmetryAst::from_rust(py, ls))
                        .collect()
                })
            }

            /// The ligand-symmetry constraint at `permutation` (undetermined if absent).
            pub(crate) fn ligand_symmetry(
                &self,
                py: Python<'_>,
                permutation: OrientedLigandPermutation,
            ) -> PyResult<LigandSymmetryAst> {
                self.read(py, |cs| {
                    LigandSymmetryAst::from_rust(py, &cs.ligand_symmetry(permutation.to_rust()))
                })
            }

            /// The fluxionality constraints.
            pub(crate) fn fluxionalities(&self, py: Python<'_>) -> PyResult<Vec<FluxionalityAst>> {
                self.read(py, |cs| {
                    cs.fluxionalities()
                        .map(|f| FluxionalityAst::from_rust(py, f))
                        .collect()
                })
            }

            /// The fluxionality constraint at `permutation` (undetermined if absent).
            pub(crate) fn fluxionality(
                &self,
                py: Python<'_>,
                permutation: LigandPermutation,
            ) -> PyResult<FluxionalityAst> {
                self.read(py, |cs| {
                    FluxionalityAst::from_rust(py, &cs.fluxionality(permutation.to_rust()))
                })
            }

            /// The topicity constraints.
            pub(crate) fn topicities(&self, py: Python<'_>) -> PyResult<Vec<TopicityAst>> {
                self.read(py, |cs| {
                    cs.topicities()
                        .map(|t| TopicityAst::from_rust(py, t))
                        .collect()
                })
            }

            /// The topicity relation at ligand `pair` (undetermined if absent).
            pub(crate) fn topicity(
                &self,
                py: Python<'_>,
                pair: StereoLigandPair,
            ) -> PyResult<TopicityRelationAst> {
                self.read(py, |cs| {
                    Ok(TopicityRelationAst::from_rust(&cs.topicity(pair.to_rust())))
                })
            }

            /// The stereogenicity constraint (undetermined if absent).
            pub(crate) fn stereogenicity(&self, py: Python<'_>) -> PyResult<StereogenicityAst> {
                self.read(py, |cs| {
                    Ok(StereogenicityAst::from_rust(&cs.stereogenicity()))
                })
            }
        }
    };
}

stereo_constraints! {
    StereoAtomConstraintKey, StereoAtomConstraintAst, StereoAtomConstraintsAst,
    StereoAtomConstraintsUpdate, ResolvedStereoAtomConstraintsUpdate, StereoAtomConstraintsArg,
    StereoAtomConstraintKeyIter, StereoAtomConstraintIter, StereoAtomConstraintItemsIter,
    AstStereoAtomConstraintKey, AstStereoAtomConstraintAst, AstStereoAtomConstraintsAst,
    StereoAtomAst, StereoAtomConstraintsView, StereoAtomConstraintsBacking,
    AstStereoAtomId, stereo_atoms, stereo_atom_mut, "stereo atom id out of range",
}

stereo_constraints! {
    StereoBondConstraintKey, StereoBondConstraintAst, StereoBondConstraintsAst,
    StereoBondConstraintsUpdate, ResolvedStereoBondConstraintsUpdate, StereoBondConstraintsArg,
    StereoBondConstraintKeyIter, StereoBondConstraintIter, StereoBondConstraintItemsIter,
    AstStereoBondConstraintKey, AstStereoBondConstraintAst, AstStereoBondConstraintsAst,
    StereoBondAst, StereoBondConstraintsView, StereoBondConstraintsBacking,
    AstStereoBondId, stereo_bonds, stereo_bond_mut, "stereo bond id out of range",
}
