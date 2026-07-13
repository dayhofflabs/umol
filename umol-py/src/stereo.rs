//! Stereo sub-ASTs mirroring `umol_ast::ast::stereo` and `umol_perm` (S5a): the
//! `Permutation` value, the recursive `StereoTerm` transformation algebra, the
//! `StereoCosetAst` coset, and the `TetrahedralStereoAst` atom-stereo state.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use std::collections::BTreeSet;
use std::vec::IntoIter;

use pyo3::exceptions::{PyKeyError, PyValueError};
use pyo3::prelude::*;
// A few `from_ast` mirrors stay `#[cfg(test)]` until S2/S4a wire them into the constraint
// enum and the entity view; their `to_ast` peers are already live (eq/hash).
#[cfg(test)]
use umol_ast::ast::{
    AtomId as AstAtomId, BooleanAst as AstBooleanAst, StereoLigand as AstStereoLigand,
    StereoLigandKind as AstStereoLigandKind,
};
use umol_ast::ast::{
    CisTransStereoAst as AstCisTransStereoAst, FluxionalityAst as AstFluxionalityAst, Lattice,
    LigandPermutation as AstLigandPermutation, LigandSymmetryAst as AstLigandSymmetryAst,
    OrientedLigandPermutation as AstOrientedLigandPermutation,
    StereoAtomConstraintAst as AstStereoAtomConstraintAst,
    StereoAtomConstraintKey as AstStereoAtomConstraintKey,
    StereoAtomConstraintsAst as AstStereoAtomConstraintsAst,
    StereoBondConstraintAst as AstStereoBondConstraintAst,
    StereoBondConstraintKey as AstStereoBondConstraintKey,
    StereoBondConstraintsAst as AstStereoBondConstraintsAst,
    StereoConfigurationAst as AstStereoConfigurationAst, StereoCosetAst as AstStereoCosetAst,
    StereoKind as AstStereoKind, StereoLigandPair as AstStereoLigandPair,
    StereoLigandPosition as AstStereoLigandPosition, StereoTerm as AstStereoTerm,
    Stereogenicity as AstStereogenicity, StereogenicityAst as AstStereogenicityAst,
    TetrahedralStereoAst as AstTetrahedralStereoAst, Topicity as AstTopicity,
    TopicityAst as AstTopicityAst, TopicityRelationAst as AstTopicityRelationAst,
};
use umol_perm::{Orientation as PermOrientation, Permutation as PermPermutation};

use crate::boolean::{BooleanArg, BooleanAst};
use crate::convert::{hash_ast, into_py_variant, variant_repr};

/// A permutation of `0..degree` in one-line (image) notation.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Permutation(PermPermutation);

#[pymethods]
impl Permutation {
    /// Construct from the image (one-line notation); the degree is the image length.
    #[new]
    fn new(image: Vec<u32>) -> Self {
        let image: Vec<u8> = image.iter().map(|&index| index as u8).collect();
        Permutation(PermPermutation::from_image(image.len(), &image))
    }

    /// The identity permutation on `0..degree`.
    #[staticmethod]
    fn identity(degree: usize) -> Self {
        Permutation(PermPermutation::identity(degree))
    }

    #[getter]
    fn degree(&self) -> usize {
        self.0.degree()
    }

    /// The image in one-line notation.
    fn image(&self) -> Vec<u32> {
        (0..self.0.degree())
            .map(|i| self.0.apply(i) as u32)
            .collect()
    }

    fn __repr__(&self) -> String {
        format!("Permutation({:?})", self.image())
    }
}

impl Permutation {
    pub(crate) fn inner(&self) -> PermPermutation {
        self.0
    }

    pub(crate) fn from_inner(permutation: PermPermutation) -> Self {
        Permutation(permutation)
    }
}

/// Stereo coset transformation term: a literal/set coset, a variable, or a
/// swap/mirror/permutation applied to a sub-term.
#[pyclass]
pub enum StereoTerm {
    Var(String, Option<BTreeSet<u32>>),
    Lit(u32),
    LitSet(BTreeSet<u32>),
    Swap(Py<StereoTerm>),
    Mirror(Py<StereoTerm>),
    Apply(Py<StereoTerm>, Permutation),
}

#[pymethods]
impl StereoTerm {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            StereoTerm::Var(_, _) => ("Var", 2),
            StereoTerm::Lit(_) => ("Lit", 1),
            StereoTerm::LitSet(_) => ("LitSet", 1),
            StereoTerm::Swap(_) => ("Swap", 1),
            StereoTerm::Mirror(_) => ("Mirror", 1),
            StereoTerm::Apply(_, _) => ("Apply", 2),
        };
        variant_repr(slf.bind(py).as_any(), "StereoTerm", variant, arity)
    }
}

impl StereoTerm {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstStereoTerm) -> PyResult<Self> {
        Ok(match ast {
            AstStereoTerm::Var(boxed) => {
                let (name, restriction) = &**boxed;
                StereoTerm::Var(name.clone(), restriction.clone())
            }
            AstStereoTerm::Lit(index) => StereoTerm::Lit(*index),
            AstStereoTerm::LitSet(members) => StereoTerm::LitSet(members.clone()),
            AstStereoTerm::Swap(inner) => {
                StereoTerm::Swap(into_py_variant(py, StereoTerm::from_ast(py, inner)?)?)
            }
            AstStereoTerm::Mirror(inner) => {
                StereoTerm::Mirror(into_py_variant(py, StereoTerm::from_ast(py, inner)?)?)
            }
            AstStereoTerm::Apply(inner, permutation) => StereoTerm::Apply(
                into_py_variant(py, StereoTerm::from_ast(py, inner)?)?,
                Permutation::from_inner(*permutation),
            ),
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstStereoTerm {
        match self {
            StereoTerm::Var(name, restriction) => {
                AstStereoTerm::Var(Box::new((name.clone(), restriction.clone())))
            }
            StereoTerm::Lit(index) => AstStereoTerm::Lit(*index),
            StereoTerm::LitSet(members) => AstStereoTerm::LitSet(members.clone()),
            StereoTerm::Swap(inner) => {
                AstStereoTerm::Swap(Box::new(inner.bind(py).borrow().to_ast(py)))
            }
            StereoTerm::Mirror(inner) => {
                AstStereoTerm::Mirror(Box::new(inner.bind(py).borrow().to_ast(py)))
            }
            StereoTerm::Apply(inner, permutation) => AstStereoTerm::Apply(
                Box::new(inner.bind(py).borrow().to_ast(py)),
                permutation.inner(),
            ),
        }
    }
}

/// Stereo coset: undetermined, a literal coset index, a set of indices, or a
/// transformation term.
#[pyclass]
pub enum StereoCosetAst {
    Undetermined(),
    Lit(u32),
    LitSet(BTreeSet<u32>),
    Term(Py<StereoTerm>),
}

#[pymethods]
impl StereoCosetAst {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            StereoCosetAst::Undetermined() => ("Undetermined", 0),
            StereoCosetAst::Lit(_) => ("Lit", 1),
            StereoCosetAst::LitSet(_) => ("LitSet", 1),
            StereoCosetAst::Term(_) => ("Term", 1),
        };
        variant_repr(slf.bind(py).as_any(), "StereoCosetAst", variant, arity)
    }
}

impl StereoCosetAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstStereoCosetAst) -> PyResult<Self> {
        Ok(match ast {
            AstStereoCosetAst::Undetermined => Self::Undetermined(),
            AstStereoCosetAst::Lit(index) => Self::Lit(*index),
            AstStereoCosetAst::LitSet(members) => Self::LitSet(members.clone()),
            AstStereoCosetAst::Term(inner) => {
                Self::Term(into_py_variant(py, StereoTerm::from_ast(py, inner)?)?)
            }
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstStereoCosetAst {
        match self {
            Self::Undetermined() => AstStereoCosetAst::Undetermined,
            Self::Lit(index) => AstStereoCosetAst::Lit(*index),
            Self::LitSet(members) => AstStereoCosetAst::LitSet(members.clone()),
            Self::Term(inner) => {
                AstStereoCosetAst::Term(Box::new(inner.bind(py).borrow().to_ast(py)))
            }
        }
    }
}

/// Tetrahedral atom stereo: undetermined, explicitly not stereogenic, or a
/// stereo coset.
#[pyclass]
pub enum TetrahedralStereoAst {
    Undetermined(),
    NotStereo(),
    Stereo(Py<StereoCosetAst>),
}

#[pymethods]
impl TetrahedralStereoAst {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            TetrahedralStereoAst::Undetermined() => ("Undetermined", 0),
            TetrahedralStereoAst::NotStereo() => ("NotStereo", 0),
            TetrahedralStereoAst::Stereo(_) => ("Stereo", 1),
        };
        variant_repr(
            slf.bind(py).as_any(),
            "TetrahedralStereoAst",
            variant,
            arity,
        )
    }
}

impl TetrahedralStereoAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstTetrahedralStereoAst) -> PyResult<Self> {
        Ok(match ast {
            AstTetrahedralStereoAst::Undetermined => Self::Undetermined(),
            AstTetrahedralStereoAst::NotStereo => Self::NotStereo(),
            AstTetrahedralStereoAst::Stereo(coset) => {
                Self::Stereo(into_py_variant(py, StereoCosetAst::from_ast(py, coset)?)?)
            }
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstTetrahedralStereoAst {
        match self {
            Self::Undetermined() => AstTetrahedralStereoAst::Undetermined,
            Self::NotStereo() => AstTetrahedralStereoAst::NotStereo,
            Self::Stereo(coset) => {
                AstTetrahedralStereoAst::Stereo(coset.bind(py).borrow().to_ast(py))
            }
        }
    }
}

/// Tetrahedral stereo configuration shorthand: counterclockwise (`Ccw`, coset
/// `Th0`) or clockwise (`Cw`, coset `Th1`).
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub enum TetrahedralStereo {
    Ccw,
    Cw,
}

impl TetrahedralStereo {
    /// The tetrahedral-stereo AST for this configuration (a literal coset).
    pub(crate) fn to_ast(self) -> AstTetrahedralStereoAst {
        let coset = match self {
            TetrahedralStereo::Ccw => AstStereoCosetAst::Lit(0),
            TetrahedralStereo::Cw => AstStereoCosetAst::Lit(1),
        };
        AstTetrahedralStereoAst::Stereo(coset)
    }
}

/// Cis/trans bond stereo: undetermined, explicitly not stereogenic, or a stereo coset.
#[pyclass]
pub enum CisTransStereoAst {
    Undetermined(),
    NotStereo(),
    Stereo(Py<StereoCosetAst>),
}

#[pymethods]
impl CisTransStereoAst {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            CisTransStereoAst::Undetermined() => ("Undetermined", 0),
            CisTransStereoAst::NotStereo() => ("NotStereo", 0),
            CisTransStereoAst::Stereo(_) => ("Stereo", 1),
        };
        variant_repr(slf.bind(py).as_any(), "CisTransStereoAst", variant, arity)
    }
}

impl CisTransStereoAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstCisTransStereoAst) -> PyResult<Self> {
        Ok(match ast {
            AstCisTransStereoAst::Undetermined => Self::Undetermined(),
            AstCisTransStereoAst::NotStereo => Self::NotStereo(),
            AstCisTransStereoAst::Stereo(coset) => {
                Self::Stereo(into_py_variant(py, StereoCosetAst::from_ast(py, coset)?)?)
            }
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstCisTransStereoAst {
        match self {
            Self::Undetermined() => AstCisTransStereoAst::Undetermined,
            Self::NotStereo() => AstCisTransStereoAst::NotStereo,
            Self::Stereo(coset) => AstCisTransStereoAst::Stereo(coset.bind(py).borrow().to_ast(py)),
        }
    }
}

/// Cis/trans stereo configuration shorthand: `Z` (coset `Ct0`) or `E` (coset `Ct1`),
/// named for the chemistry keywords.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub enum CisTransStereo {
    Z,
    E,
}

impl CisTransStereo {
    /// The cis/trans-stereo AST for this configuration (a literal coset).
    pub(crate) fn to_ast(self) -> AstCisTransStereoAst {
        let coset = match self {
            CisTransStereo::Z => AstStereoCosetAst::Lit(0),
            CisTransStereo::E => AstStereoCosetAst::Lit(1),
        };
        AstCisTransStereoAst::Stereo(coset)
    }
}

/// Setter coercion for `cis_trans_stereo`: `False` → not stereogenic, a
/// `CisTransStereo` (`Z`/`E`) → that coset, or a `CisTransStereoAst` passthrough.
#[derive(FromPyObject)]
pub(crate) enum CisTransStereoArg {
    Flag(bool),
    Config(CisTransStereo),
    Ast(Py<CisTransStereoAst>),
}

impl CisTransStereoArg {
    pub(crate) fn to_ast(&self, py: Python<'_>) -> PyResult<AstCisTransStereoAst> {
        Ok(match self {
            CisTransStereoArg::Flag(false) => AstCisTransStereoAst::NotStereo,
            CisTransStereoArg::Flag(true) => {
                return Err(PyValueError::new_err(
                    "cis_trans_stereo = True is not meaningful; use CisTransStereo.Z/E or False",
                ))
            }
            CisTransStereoArg::Config(cts) => cts.to_ast(),
            CisTransStereoArg::Ast(a) => a.bind(py).borrow().to_ast(py),
        })
    }
}

/// The coordination geometry of a stereo site. A fieldless, hashable value enum whose
/// members mirror the Rust `StereoKind` exactly.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum StereoKind {
    Tetrahedral,
    CisTrans,
    Axial,
    SquarePlanar,
    TrigonalBipyramidal,
    Octahedral,
}

impl StereoKind {
    #[cfg(test)]
    pub(crate) fn from_ast(ast: AstStereoKind) -> Self {
        match ast {
            AstStereoKind::Tetrahedral => Self::Tetrahedral,
            AstStereoKind::CisTrans => Self::CisTrans,
            AstStereoKind::Axial => Self::Axial,
            AstStereoKind::SquarePlanar => Self::SquarePlanar,
            AstStereoKind::TrigonalBipyramidal => Self::TrigonalBipyramidal,
            AstStereoKind::Octahedral => Self::Octahedral,
        }
    }

    pub(crate) fn to_ast(self) -> AstStereoKind {
        match self {
            Self::Tetrahedral => AstStereoKind::Tetrahedral,
            Self::CisTrans => AstStereoKind::CisTrans,
            Self::Axial => AstStereoKind::Axial,
            Self::SquarePlanar => AstStereoKind::SquarePlanar,
            Self::TrigonalBipyramidal => AstStereoKind::TrigonalBipyramidal,
            Self::Octahedral => AstStereoKind::Octahedral,
        }
    }
}

/// The kind of a stereo ligand: a real atom, or a virtual ligand (implicit hydrogen or
/// lone pair). A fieldless, hashable value enum mirroring the Rust `StereoLigandKind`.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum StereoLigandKind {
    Atom,
    ImplicitHydrogen,
    LonePair,
}

impl StereoLigandKind {
    #[cfg(test)]
    pub(crate) fn from_ast(ast: AstStereoLigandKind) -> Self {
        match ast {
            AstStereoLigandKind::Atom => Self::Atom,
            AstStereoLigandKind::ImplicitHydrogen => Self::ImplicitHydrogen,
            AstStereoLigandKind::LonePair => Self::LonePair,
        }
    }

    #[cfg(test)]
    pub(crate) fn to_ast(self) -> AstStereoLigandKind {
        match self {
            Self::Atom => AstStereoLigandKind::Atom,
            Self::ImplicitHydrogen => AstStereoLigandKind::ImplicitHydrogen,
            Self::LonePair => AstStereoLigandKind::LonePair,
        }
    }
}

/// Topicity of two ligand positions of a stereo carrier (a derived ground classification).
/// A fieldless, hashable value enum mirroring the Rust `Topicity`. `Ord` lets it key the
/// `BTreeSet` in the `TopicityRelationAst` set variants.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Topicity {
    Homotopic,
    Enantiotopic,
    Diastereotopic,
}

impl Topicity {
    pub(crate) fn from_ast(ast: AstTopicity) -> Self {
        match ast {
            AstTopicity::Homotopic => Self::Homotopic,
            AstTopicity::Enantiotopic => Self::Enantiotopic,
            AstTopicity::Diastereotopic => Self::Diastereotopic,
        }
    }

    pub(crate) fn to_ast(self) -> AstTopicity {
        match self {
            Self::Homotopic => AstTopicity::Homotopic,
            Self::Enantiotopic => AstTopicity::Enantiotopic,
            Self::Diastereotopic => AstTopicity::Diastereotopic,
        }
    }
}

/// Stereogenicity classification of a stereo carrier (a derived ground classification).
/// A fieldless, hashable value enum mirroring the Rust `Stereogenicity`. `Ord` lets it key
/// the `BTreeSet` in the `StereogenicityAst` set variants.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Stereogenicity {
    Symmetric,
    Prochiral,
    Stereogenic,
}

impl Stereogenicity {
    pub(crate) fn from_ast(ast: AstStereogenicity) -> Self {
        match ast {
            AstStereogenicity::Symmetric => Self::Symmetric,
            AstStereogenicity::Prochiral => Self::Prochiral,
            AstStereogenicity::Stereogenic => Self::Stereogenic,
        }
    }

    pub(crate) fn to_ast(self) -> AstStereogenicity {
        match self {
            Self::Symmetric => AstStereogenicity::Symmetric,
            Self::Prochiral => AstStereogenicity::Prochiral,
            Self::Stereogenic => AstStereogenicity::Stereogenic,
        }
    }
}

/// A stereo ligand occupying a coordination position of a stereo site: the ligand's atom
/// id and its kind. For a virtual ligand (`ImplicitHydrogen`/`LonePair`) the `atom_id` is
/// the bearing atom; the `kind` disambiguates. Immutable value, hashable. Mirrors the Rust
/// `StereoLigand`.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct StereoLigand {
    #[pyo3(get)]
    atom_id: u32,
    #[pyo3(get)]
    kind: StereoLigandKind,
}

#[pymethods]
impl StereoLigand {
    #[new]
    fn new(atom_id: u32, kind: StereoLigandKind) -> Self {
        StereoLigand { atom_id, kind }
    }

    fn __repr__(&self) -> String {
        format!(
            "StereoLigand(atom_id={}, kind=StereoLigandKind.{:?})",
            self.atom_id, self.kind
        )
    }
}

impl StereoLigand {
    #[cfg(test)]
    pub(crate) fn from_ast(ast: AstStereoLigand) -> Self {
        StereoLigand {
            atom_id: ast.atom_id.0,
            kind: StereoLigandKind::from_ast(ast.kind),
        }
    }

    #[cfg(test)]
    pub(crate) fn to_ast(self) -> AstStereoLigand {
        AstStereoLigand::new(AstAtomId(self.atom_id), self.kind.to_ast())
    }
}

/// A stereo configuration: undetermined (geometry not yet known, so no coset), or `Kinded`
/// — a concrete coordination geometry bound to a coset that may still be open. Mirrors the
/// Rust `StereoConfigurationAst`; `Undetermined` and `Kinded(Tetrahedral, Undetermined)` are
/// distinct.
#[pyclass]
pub enum StereoConfigurationAst {
    Undetermined(),
    Kinded(StereoKind, Py<StereoCosetAst>),
}

#[pymethods]
impl StereoConfigurationAst {
    /// The coordination-geometry kind, or `None` when undetermined.
    #[getter]
    fn kind(&self) -> Option<StereoKind> {
        match self {
            Self::Kinded(kind, _) => Some(*kind),
            Self::Undetermined() => None,
        }
    }

    /// The coset, or `None` when undetermined.
    #[getter]
    fn coset(&self, py: Python<'_>) -> Option<Py<StereoCosetAst>> {
        match self {
            Self::Kinded(_, coset) => Some(coset.clone_ref(py)),
            Self::Undetermined() => None,
        }
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            StereoConfigurationAst::Undetermined() => ("Undetermined", 0),
            StereoConfigurationAst::Kinded(_, _) => ("Kinded", 2),
        };
        variant_repr(
            slf.bind(py).as_any(),
            "StereoConfigurationAst",
            variant,
            arity,
        )
    }
}

impl StereoConfigurationAst {
    #[cfg(test)]
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstStereoConfigurationAst) -> PyResult<Self> {
        Ok(match ast {
            AstStereoConfigurationAst::Undetermined => Self::Undetermined(),
            AstStereoConfigurationAst::Kinded(kind, coset) => Self::Kinded(
                StereoKind::from_ast(*kind),
                into_py_variant(py, StereoCosetAst::from_ast(py, coset)?)?,
            ),
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstStereoConfigurationAst {
        match self {
            Self::Undetermined() => AstStereoConfigurationAst::Undetermined,
            Self::Kinded(kind, coset) => {
                AstStereoConfigurationAst::Kinded(kind.to_ast(), coset.bind(py).borrow().to_ast(py))
            }
        }
    }
}

/// Setter coercion for a stereo `configuration` field: the `TetrahedralStereo` (`Ccw`/`Cw`)
/// or `CisTransStereo` (`Z`/`E`) per-kind coset shorthand, or a `StereoConfigurationAst`
/// passthrough. Axial/square-planar/etc. have no shorthand — use the full `Kinded` form.
#[cfg(test)]
#[derive(FromPyObject)]
pub(crate) enum StereoConfigurationArg {
    Tetrahedral(TetrahedralStereo),
    CisTrans(CisTransStereo),
    Ast(Py<StereoConfigurationAst>),
}

#[cfg(test)]
impl StereoConfigurationArg {
    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstStereoConfigurationAst {
        match self {
            StereoConfigurationArg::Tetrahedral(t) => match t.to_ast() {
                AstTetrahedralStereoAst::Stereo(coset) => {
                    AstStereoConfigurationAst::Kinded(AstStereoKind::Tetrahedral, coset)
                }
                _ => unreachable!("TetrahedralStereo shorthand is always a Stereo coset"),
            },
            StereoConfigurationArg::CisTrans(c) => match c.to_ast() {
                AstCisTransStereoAst::Stereo(coset) => {
                    AstStereoConfigurationAst::Kinded(AstStereoKind::CisTrans, coset)
                }
                _ => unreachable!("CisTransStereo shorthand is always a Stereo coset"),
            },
            StereoConfigurationArg::Ast(a) => a.bind(py).borrow().to_ast(py),
        }
    }
}

/// Orientation grade of a ligand permutation: a proper rotation, or an improper (mirror)
/// operation. A fieldless, hashable value enum mirroring `umol_perm::Orientation`.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Orientation {
    Proper,
    Improper,
}

impl Orientation {
    pub(crate) fn from_ast(orientation: PermOrientation) -> Self {
        match orientation {
            PermOrientation::Proper => Self::Proper,
            PermOrientation::Improper => Self::Improper,
        }
    }

    pub(crate) fn to_ast(self) -> PermOrientation {
        match self {
            Self::Proper => PermOrientation::Proper,
            Self::Improper => PermOrientation::Improper,
        }
    }
}

/// A permutation of a stereo site's ligand positions (frame-relative). Immutable value,
/// hashable. Mirrors the Rust `LigandPermutation`.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct LigandPermutation {
    #[pyo3(get)]
    permutation: Permutation,
}

#[pymethods]
impl LigandPermutation {
    #[new]
    fn new(permutation: Permutation) -> Self {
        LigandPermutation { permutation }
    }

    /// A pattern matches a target iff they are the same permutation.
    fn matches(&self, target: &Self) -> bool {
        self.permutation == target.permutation
    }

    fn __repr__(&self) -> String {
        format!("LigandPermutation({:?})", self.permutation.image())
    }
}

impl LigandPermutation {
    pub(crate) fn from_ast(ast: AstLigandPermutation) -> Self {
        LigandPermutation {
            permutation: Permutation::from_inner(ast.0),
        }
    }

    pub(crate) fn to_ast(self) -> AstLigandPermutation {
        AstLigandPermutation(self.permutation.inner())
    }
}

/// A ligand permutation carrying a proper/improper grade. Immutable value, hashable. Mirrors
/// the Rust `OrientedLigandPermutation`.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct OrientedLigandPermutation {
    #[pyo3(get)]
    permutation: LigandPermutation,
    #[pyo3(get)]
    orientation: Orientation,
}

#[pymethods]
impl OrientedLigandPermutation {
    #[new]
    fn new(permutation: LigandPermutation, orientation: Orientation) -> Self {
        OrientedLigandPermutation {
            permutation,
            orientation,
        }
    }

    fn matches(&self, target: &Self) -> bool {
        self.permutation.matches(&target.permutation) && self.orientation == target.orientation
    }

    fn __repr__(&self) -> String {
        format!(
            "OrientedLigandPermutation(permutation={}, orientation=Orientation.{:?})",
            self.permutation.__repr__(),
            self.orientation
        )
    }
}

impl OrientedLigandPermutation {
    pub(crate) fn from_ast(ast: AstOrientedLigandPermutation) -> Self {
        OrientedLigandPermutation {
            permutation: LigandPermutation::from_ast(ast.permutation),
            orientation: Orientation::from_ast(ast.orientation),
        }
    }

    pub(crate) fn to_ast(self) -> AstOrientedLigandPermutation {
        AstOrientedLigandPermutation {
            permutation: self.permutation.to_ast(),
            orientation: self.orientation.to_ast(),
        }
    }
}

/// An unordered pair of ligand positions of a stereo site, normalized so the lower position
/// is `first`. Keys a per-pair topicity constraint. Immutable value, hashable. Mirrors the
/// Rust `StereoLigandPair`.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct StereoLigandPair {
    #[pyo3(get)]
    first: u32,
    #[pyo3(get)]
    second: u32,
}

#[pymethods]
impl StereoLigandPair {
    /// Normalizes the pair so the lower position is `first`.
    #[new]
    fn new(a: u32, b: u32) -> Self {
        StereoLigandPair::from_ast(AstStereoLigandPair::new(
            AstStereoLigandPosition(a),
            AstStereoLigandPosition(b),
        ))
    }

    fn __repr__(&self) -> String {
        format!("StereoLigandPair({}, {})", self.first, self.second)
    }
}

impl StereoLigandPair {
    fn from_ast(ast: AstStereoLigandPair) -> Self {
        StereoLigandPair {
            first: ast.first().0,
            second: ast.second().0,
        }
    }

    pub(crate) fn to_ast(self) -> AstStereoLigandPair {
        AstStereoLigandPair::new(
            AstStereoLigandPosition(self.first),
            AstStereoLigandPosition(self.second),
        )
    }
}

/// A topicity relation constraint value: the undetermined wildcard, a single topicity, a set
/// of admissible topicities, or the complement of a set. A finite-domain subset lattice over
/// `Topicity`. Mirrors the Rust `TopicityRelationAst`.
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
    fn as_lit(&self) -> Option<Topicity> {
        match self {
            Self::Lit(topicity) => Some(*topicity),
            _ => None,
        }
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.to_ast() == other.to_ast()
    }

    fn __hash__(&self) -> u64 {
        hash_ast(&self.to_ast())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            TopicityRelationAst::Undetermined() => ("Undetermined", 0),
            TopicityRelationAst::Lit(_) => ("Lit", 1),
            TopicityRelationAst::LitSet(_) => ("LitSet", 1),
            TopicityRelationAst::NotSet(_) => ("NotSet", 1),
        };
        variant_repr(slf.bind(py).as_any(), "TopicityRelationAst", variant, arity)
    }
}

impl TopicityRelationAst {
    pub(crate) fn from_ast(ast: &AstTopicityRelationAst) -> Self {
        match ast {
            AstTopicityRelationAst::Undetermined => Self::Undetermined(),
            AstTopicityRelationAst::Lit(topicity) => Self::Lit(Topicity::from_ast(*topicity)),
            AstTopicityRelationAst::LitSet(topicities) => {
                Self::LitSet(topicities.iter().map(|t| Topicity::from_ast(*t)).collect())
            }
            AstTopicityRelationAst::NotSet(topicities) => {
                Self::NotSet(topicities.iter().map(|t| Topicity::from_ast(*t)).collect())
            }
        }
    }

    pub(crate) fn to_ast(&self) -> AstTopicityRelationAst {
        match self {
            Self::Undetermined() => AstTopicityRelationAst::Undetermined,
            Self::Lit(topicity) => AstTopicityRelationAst::Lit(topicity.to_ast()),
            Self::LitSet(topicities) => {
                AstTopicityRelationAst::LitSet(topicities.iter().map(|t| t.to_ast()).collect())
            }
            Self::NotSet(topicities) => {
                AstTopicityRelationAst::NotSet(topicities.iter().map(|t| t.to_ast()).collect())
            }
        }
    }
}

/// Setter coercion for a topicity relation: a `Topicity` literal (→ `Lit`) or a
/// `TopicityRelationAst` passthrough (mirroring `impl From<Topicity>`).
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
/// lattice over `Stereogenicity`. Mirrors the Rust `StereogenicityAst`.
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
    fn as_lit(&self) -> Option<Stereogenicity> {
        match self {
            Self::Lit(stereogenicity) => Some(*stereogenicity),
            _ => None,
        }
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.to_ast() == other.to_ast()
    }

    fn __hash__(&self) -> u64 {
        hash_ast(&self.to_ast())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            StereogenicityAst::Undetermined() => ("Undetermined", 0),
            StereogenicityAst::Lit(_) => ("Lit", 1),
            StereogenicityAst::LitSet(_) => ("LitSet", 1),
            StereogenicityAst::NotSet(_) => ("NotSet", 1),
        };
        variant_repr(slf.bind(py).as_any(), "StereogenicityAst", variant, arity)
    }
}

impl StereogenicityAst {
    pub(crate) fn from_ast(ast: &AstStereogenicityAst) -> Self {
        match ast {
            AstStereogenicityAst::Undetermined => Self::Undetermined(),
            AstStereogenicityAst::Lit(stereogenicity) => {
                Self::Lit(Stereogenicity::from_ast(*stereogenicity))
            }
            AstStereogenicityAst::LitSet(stereogenicities) => Self::LitSet(
                stereogenicities
                    .iter()
                    .map(|g| Stereogenicity::from_ast(*g))
                    .collect(),
            ),
            AstStereogenicityAst::NotSet(stereogenicities) => Self::NotSet(
                stereogenicities
                    .iter()
                    .map(|g| Stereogenicity::from_ast(*g))
                    .collect(),
            ),
        }
    }

    pub(crate) fn to_ast(&self) -> AstStereogenicityAst {
        match self {
            Self::Undetermined() => AstStereogenicityAst::Undetermined,
            Self::Lit(stereogenicity) => AstStereogenicityAst::Lit(stereogenicity.to_ast()),
            Self::LitSet(stereogenicities) => {
                AstStereogenicityAst::LitSet(stereogenicities.iter().map(|g| g.to_ast()).collect())
            }
            Self::NotSet(stereogenicities) => {
                AstStereogenicityAst::NotSet(stereogenicities.iter().map(|g| g.to_ast()).collect())
            }
        }
    }
}

/// A ligand-symmetry constraint value: an oriented ligand permutation with a presence
/// assertion (whether the permutation is a ligand symmetry). Mirrors the Rust
/// `LigandSymmetryAst`.
#[pyclass]
pub struct LigandSymmetryAst {
    permutation: OrientedLigandPermutation,
    invariant: Py<BooleanAst>,
}

#[pymethods]
impl LigandSymmetryAst {
    #[new]
    fn new(
        py: Python<'_>,
        permutation: OrientedLigandPermutation,
        invariant: BooleanArg,
    ) -> PyResult<Self> {
        Ok(LigandSymmetryAst {
            permutation,
            invariant: into_py_variant(py, BooleanAst::from_ast(&invariant.to_ast(py)))?,
        })
    }

    #[getter]
    fn permutation(&self) -> OrientedLigandPermutation {
        self.permutation
    }

    #[getter]
    fn invariant(&self, py: Python<'_>) -> Py<BooleanAst> {
        self.invariant.clone_ref(py)
    }

    /// Matches iff the permutations are equal and the presence assertions match.
    fn matches(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py).matches(&other.to_ast(py))
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
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

impl LigandSymmetryAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstLigandSymmetryAst) -> PyResult<Self> {
        Ok(LigandSymmetryAst {
            permutation: OrientedLigandPermutation::from_ast(ast.permutation),
            invariant: into_py_variant(py, BooleanAst::from_ast(&ast.invariant))?,
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstLigandSymmetryAst {
        AstLigandSymmetryAst {
            permutation: self.permutation.to_ast(),
            invariant: self.invariant.bind(py).borrow().to_ast(),
        }
    }
}

/// A fluxionality constraint value: a proper ligand permutation realized by dynamics, with an
/// assertion of whether the move is `active`. Mirrors the Rust `FluxionalityAst`.
#[pyclass]
pub struct FluxionalityAst {
    permutation: LigandPermutation,
    active: Py<BooleanAst>,
}

#[pymethods]
impl FluxionalityAst {
    #[new]
    fn new(py: Python<'_>, permutation: LigandPermutation, active: BooleanArg) -> PyResult<Self> {
        Ok(FluxionalityAst {
            permutation,
            active: into_py_variant(py, BooleanAst::from_ast(&active.to_ast(py)))?,
        })
    }

    #[getter]
    fn permutation(&self) -> LigandPermutation {
        self.permutation
    }

    #[getter]
    fn active(&self, py: Python<'_>) -> Py<BooleanAst> {
        self.active.clone_ref(py)
    }

    /// Matches iff the permutations are equal and the presence assertions match.
    fn matches(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py).matches(&other.to_ast(py))
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "FluxionalityAst({}, {})",
            self.permutation.__repr__(),
            self.active.bind(py).as_any().repr()?.extract::<String>()?,
        ))
    }
}

impl FluxionalityAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstFluxionalityAst) -> PyResult<Self> {
        Ok(FluxionalityAst {
            permutation: LigandPermutation::from_ast(ast.permutation),
            active: into_py_variant(py, BooleanAst::from_ast(&ast.active))?,
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstFluxionalityAst {
        AstFluxionalityAst {
            permutation: self.permutation.to_ast(),
            active: self.active.bind(py).borrow().to_ast(),
        }
    }
}

/// A per-pair topicity constraint value: a relation between a pair of ligand positions.
/// Mirrors the Rust `TopicityAst`.
#[pyclass]
pub struct TopicityAst {
    pair: StereoLigandPair,
    relation: Py<TopicityRelationAst>,
}

#[pymethods]
impl TopicityAst {
    #[new]
    fn new(
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
    fn pair(&self) -> StereoLigandPair {
        self.pair
    }

    #[getter]
    fn relation(&self, py: Python<'_>) -> Py<TopicityRelationAst> {
        self.relation.clone_ref(py)
    }

    /// Matches iff the pairs are equal and the per-pair relations match.
    fn matches(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py).matches(&other.to_ast(py))
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
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

impl TopicityAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstTopicityAst) -> PyResult<Self> {
        Ok(TopicityAst {
            pair: StereoLigandPair::from_ast(ast.pair),
            relation: into_py_variant(py, TopicityRelationAst::from_ast(&ast.relation))?,
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstTopicityAst {
        AstTopicityAst {
            pair: self.pair.to_ast(),
            relation: self.relation.bind(py).borrow().to_ast(),
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
        $ast_key:ident, $ast_constraint:ident, $ast_constraints:ident $(,)?
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
            fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
                self.to_ast(py) == other.to_ast(py)
            }

            fn __hash__(&self, py: Python<'_>) -> u64 {
                hash_ast(&self.to_ast(py))
            }

            fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
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
            pub(crate) fn from_ast(py: Python<'_>, ast: &$ast_key) -> PyResult<Self> {
                Ok(match ast {
                    $ast_key::LigandSymmetry(permutation) => Self::LigandSymmetry(into_py_variant(
                        py,
                        OrientedLigandPermutation::from_ast(*permutation),
                    )?),
                    $ast_key::Fluxionality(permutation) => Self::Fluxionality(into_py_variant(
                        py,
                        LigandPermutation::from_ast(*permutation),
                    )?),
                    $ast_key::Topicity(pair) => {
                        Self::Topicity(into_py_variant(py, StereoLigandPair::from_ast(*pair))?)
                    }
                    $ast_key::Stereogenicity => Self::Stereogenicity(),
                })
            }

            pub(crate) fn to_ast(&self, py: Python<'_>) -> $ast_key {
                match self {
                    Self::LigandSymmetry(permutation) => {
                        $ast_key::LigandSymmetry(permutation.bind(py).borrow().to_ast())
                    }
                    Self::Fluxionality(permutation) => {
                        $ast_key::Fluxionality(permutation.bind(py).borrow().to_ast())
                    }
                    Self::Topicity(pair) => $ast_key::Topicity(pair.bind(py).borrow().to_ast()),
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
            fn key(&self, py: Python<'_>) -> PyResult<$key> {
                $key::from_ast(py, &self.to_ast(py).key())
            }

            fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
                self.to_ast(py) == other.to_ast(py)
            }

            fn __hash__(&self, py: Python<'_>) -> u64 {
                hash_ast(&self.to_ast(py))
            }

            fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
                let variant = match &*slf.bind(py).borrow() {
                    $constraint::LigandSymmetry(_) => "LigandSymmetry",
                    $constraint::Fluxionality(_) => "Fluxionality",
                    $constraint::Topicity(_) => "Topicity",
                    $constraint::Stereogenicity(_) => "Stereogenicity",
                };
                variant_repr(slf.bind(py).as_any(), stringify!($constraint), variant, 1)
            }
        }

        impl $constraint {
            pub(crate) fn from_ast(py: Python<'_>, ast: &$ast_constraint) -> PyResult<Self> {
                Ok(match ast {
                    $ast_constraint::LigandSymmetry(value) => Self::LigandSymmetry(
                        into_py_variant(py, LigandSymmetryAst::from_ast(py, value)?)?,
                    ),
                    $ast_constraint::Fluxionality(value) => Self::Fluxionality(into_py_variant(
                        py,
                        FluxionalityAst::from_ast(py, value)?,
                    )?),
                    $ast_constraint::Topicity(value) => {
                        Self::Topicity(into_py_variant(py, TopicityAst::from_ast(py, value)?)?)
                    }
                    $ast_constraint::Stereogenicity(value) => Self::Stereogenicity(
                        into_py_variant(py, StereogenicityAst::from_ast(value))?,
                    ),
                })
            }

            pub(crate) fn to_ast(&self, py: Python<'_>) -> $ast_constraint {
                match self {
                    Self::LigandSymmetry(value) => {
                        $ast_constraint::LigandSymmetry(value.bind(py).borrow().to_ast(py))
                    }
                    Self::Fluxionality(value) => {
                        $ast_constraint::Fluxionality(value.bind(py).borrow().to_ast(py))
                    }
                    Self::Topicity(value) => {
                        $ast_constraint::Topicity(value.bind(py).borrow().to_ast(py))
                    }
                    Self::Stereogenicity(value) => {
                        $ast_constraint::Stereogenicity(value.bind(py).borrow().to_ast())
                    }
                }
            }
        }

        /// Argument to the container's `update`: another container or a loose iterable of
        /// constraints. (A live-view arm ships with the entity view in S2b.)
        #[derive(FromPyObject)]
        enum $update {
            Container(Py<$constraints>),
            Entries(Vec<Py<$constraint>>),
        }

        impl $update {
            /// Read every Python object into owned data before any write borrow is taken, so a
            /// container that aliases the same entity is read while nothing is borrowed
            /// (otherwise `cs.update(cs)` self-aliases into a double-borrow panic).
            fn resolve(&self, py: Python<'_>) -> PyResult<$resolved> {
                Ok(match self {
                    $update::Container(c) => {
                        $resolved::Overlay(c.bind(py).borrow().inner().clone())
                    }
                    $update::Entries(entries) => $resolved::Entries(
                        entries
                            .iter()
                            .map(|entry| entry.bind(py).borrow().to_ast(py))
                            .collect(),
                    ),
                })
            }
        }

        /// A `$update` with all Python reads done, applicable under a write borrow.
        enum $resolved {
            Overlay($ast_constraints),
            Entries(Vec<$ast_constraint>),
        }

        impl $resolved {
            fn apply(self, target: &mut $ast_constraints) {
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

        /// A whole-container argument for the entity `constraints` setter. (A live-view arm
        /// ships in S2b; the value pyclass that consumes this lands in S3 — gated until then.)
        #[cfg(test)]
        #[derive(FromPyObject)]
        pub(crate) enum $arg {
            Container(Py<$constraints>),
        }

        #[cfg(test)]
        impl $arg {
            pub(crate) fn to_ast(&self, py: Python<'_>) -> $ast_constraints {
                match self {
                    $arg::Container(c) => c.bind(py).borrow().inner().clone(),
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
            fn new(py: Python<'_>, entries: Vec<Py<$constraint>>) -> Self {
                let mut constraints = $ast_constraints::new();
                constraints.extend(
                    entries
                        .into_iter()
                        .map(|entry| entry.bind(py).borrow().to_ast(py)),
                );
                $constraints(constraints)
            }

            fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
                let mut parts = Vec::with_capacity(self.0.len());
                for entry in self.0.iter() {
                    let mirror = into_py_variant(py, $constraint::from_ast(py, entry)?)?;
                    parts.push(mirror.bind(py).as_any().repr()?.extract::<String>()?);
                }
                Ok(format!(
                    "{}([{}])",
                    stringify!($constraints),
                    parts.join(", ")
                ))
            }

            /// Insert `c`, replacing any existing entry of the same key (last-wins).
            fn set(&mut self, py: Python<'_>, c: Py<$constraint>) {
                self.0.set(c.bind(py).borrow().to_ast(py));
            }

            /// Remove the entry with the given key, returning it if present (dict `pop`).
            fn pop(&mut self, py: Python<'_>, key: Py<$key>) -> PyResult<Option<$constraint>> {
                self.0
                    .remove(key.bind(py).borrow().to_ast(py))
                    .map(|c| $constraint::from_ast(py, &c))
                    .transpose()
            }

            /// Overlay `other` onto self in place — another container or an iterable of
            /// constraints (last-wins per key; undetermined entries remove). Takes `slf` by
            /// handle so `other` is fully read before the write borrow (`cs.update(cs)` is a
            /// no-op, not a double-borrow panic).
            fn update(slf: Py<Self>, py: Python<'_>, other: $update) -> PyResult<()> {
                let resolved = other.resolve(py)?;
                resolved.apply(&mut slf.borrow_mut(py).0);
                Ok(())
            }

            fn __len__(&self) -> usize {
                self.0.len()
            }

            /// Iterate the constraint keys (mapping-style, canonical order).
            fn __iter__(&self, py: Python<'_>) -> PyResult<$key_iter> {
                self.keys(py)
            }

            /// The constraint keys, in canonical order.
            fn keys(&self, py: Python<'_>) -> PyResult<$key_iter> {
                let keys = self
                    .0
                    .iter()
                    .map(|c| into_py_variant(py, $key::from_ast(py, &c.key())?))
                    .collect::<PyResult<Vec<_>>>()?;
                Ok($key_iter {
                    keys: keys.into_iter(),
                })
            }

            /// The constraints, in canonical order.
            fn values(&self, py: Python<'_>) -> PyResult<$iter> {
                let entries = self
                    .0
                    .iter()
                    .map(|c| into_py_variant(py, $constraint::from_ast(py, c)?))
                    .collect::<PyResult<Vec<_>>>()?;
                Ok($iter {
                    entries: entries.into_iter(),
                })
            }

            /// The `(key, constraint)` pairs, in canonical order.
            fn items(&self, py: Python<'_>) -> PyResult<$items_iter> {
                let items = self
                    .0
                    .iter()
                    .map(|c| {
                        Ok((
                            into_py_variant(py, $key::from_ast(py, &c.key())?)?,
                            into_py_variant(py, $constraint::from_ast(py, c)?)?,
                        ))
                    })
                    .collect::<PyResult<Vec<_>>>()?;
                Ok($items_iter {
                    items: items.into_iter(),
                })
            }

            /// The constraint with the given key, or `default` (`None`) if absent.
            #[pyo3(signature = (key, default=None))]
            fn get(
                &self,
                py: Python<'_>,
                key: Py<$key>,
                default: Option<Py<PyAny>>,
            ) -> PyResult<Py<PyAny>> {
                match self.0.get(key.bind(py).borrow().to_ast(py)) {
                    Some(constraint) => {
                        Ok(into_py_variant(py, $constraint::from_ast(py, constraint)?)?.into_any())
                    }
                    None => Ok(default.unwrap_or_else(|| py.None())),
                }
            }

            /// The constraint with the given key; raises `KeyError` if absent.
            fn __getitem__(&self, py: Python<'_>, key: Py<$key>) -> PyResult<$constraint> {
                match self.0.get(key.bind(py).borrow().to_ast(py)) {
                    Some(constraint) => $constraint::from_ast(py, constraint),
                    None => Err(PyKeyError::new_err(
                        key.bind(py).as_any().repr()?.extract::<String>()?,
                    )),
                }
            }

            /// Remove the entry with the given key; raises `KeyError` if absent.
            fn __delitem__(&mut self, py: Python<'_>, key: Py<$key>) -> PyResult<()> {
                if self.0.remove(key.bind(py).borrow().to_ast(py)).is_some() {
                    Ok(())
                } else {
                    Err(PyKeyError::new_err(
                        key.bind(py).as_any().repr()?.extract::<String>()?,
                    ))
                }
            }

            fn __contains__(&self, py: Python<'_>, key: Py<$key>) -> bool {
                self.0.contains(key.bind(py).borrow().to_ast(py))
            }

            /// The ligand-symmetry constraints.
            fn ligand_symmetries(&self, py: Python<'_>) -> PyResult<Vec<LigandSymmetryAst>> {
                self.0
                    .ligand_symmetries()
                    .map(|ls| LigandSymmetryAst::from_ast(py, ls))
                    .collect()
            }

            /// The ligand-symmetry constraint at `permutation` (undetermined if absent).
            fn ligand_symmetry(
                &self,
                py: Python<'_>,
                permutation: OrientedLigandPermutation,
            ) -> PyResult<LigandSymmetryAst> {
                LigandSymmetryAst::from_ast(py, &self.0.ligand_symmetry(permutation.to_ast()))
            }

            /// The fluxionality constraints.
            fn fluxionalities(&self, py: Python<'_>) -> PyResult<Vec<FluxionalityAst>> {
                self.0
                    .fluxionalities()
                    .map(|f| FluxionalityAst::from_ast(py, f))
                    .collect()
            }

            /// The fluxionality constraint at `permutation` (undetermined if absent).
            fn fluxionality(
                &self,
                py: Python<'_>,
                permutation: LigandPermutation,
            ) -> PyResult<FluxionalityAst> {
                FluxionalityAst::from_ast(py, &self.0.fluxionality(permutation.to_ast()))
            }

            /// The topicity constraints.
            fn topicities(&self, py: Python<'_>) -> PyResult<Vec<TopicityAst>> {
                self.0
                    .topicities()
                    .map(|t| TopicityAst::from_ast(py, t))
                    .collect()
            }

            /// The topicity relation at ligand `pair` (undetermined if absent).
            fn topicity(&self, pair: StereoLigandPair) -> TopicityRelationAst {
                TopicityRelationAst::from_ast(&self.0.topicity(pair.to_ast()))
            }

            /// The stereogenicity constraint (undetermined if absent).
            fn stereogenicity(&self) -> StereogenicityAst {
                StereogenicityAst::from_ast(&self.0.stereogenicity())
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

        #[pyclass]
        pub struct $key_iter {
            keys: IntoIter<Py<$key>>,
        }

        #[pymethods]
        impl $key_iter {
            fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
                slf
            }

            fn __next__(&mut self) -> Option<Py<$key>> {
                self.keys.next()
            }
        }

        #[pyclass]
        pub struct $iter {
            entries: IntoIter<Py<$constraint>>,
        }

        #[pymethods]
        impl $iter {
            fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
                slf
            }

            fn __next__(&mut self) -> Option<Py<$constraint>> {
                self.entries.next()
            }
        }

        #[pyclass]
        pub struct $items_iter {
            items: IntoIter<(Py<$key>, Py<$constraint>)>,
        }

        #[pymethods]
        impl $items_iter {
            fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
                slf
            }

            fn __next__(&mut self) -> Option<(Py<$key>, Py<$constraint>)> {
                self.items.next()
            }
        }
    };
}

stereo_constraints! {
    StereoAtomConstraintKey, StereoAtomConstraintAst, StereoAtomConstraintsAst,
    StereoAtomConstraintsUpdate, ResolvedStereoAtomConstraintsUpdate, StereoAtomConstraintsArg,
    StereoAtomConstraintKeyIter, StereoAtomConstraintIter, StereoAtomConstraintItemsIter,
    AstStereoAtomConstraintKey, AstStereoAtomConstraintAst, AstStereoAtomConstraintsAst,
}

stereo_constraints! {
    StereoBondConstraintKey, StereoBondConstraintAst, StereoBondConstraintsAst,
    StereoBondConstraintsUpdate, ResolvedStereoBondConstraintsUpdate, StereoBondConstraintsArg,
    StereoBondConstraintKeyIter, StereoBondConstraintIter, StereoBondConstraintItemsIter,
    AstStereoBondConstraintKey, AstStereoBondConstraintAst, AstStereoBondConstraintsAst,
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case(vec![0, 1, 2, 3])]
    #[case(vec![1, 0, 2, 3])]
    #[case(vec![2, 0, 1])]
    fn test_permutation_image(#[case] image: Vec<u32>) {
        let permutation = Permutation::new(image.clone());
        assert_eq!(permutation.image(), image);
        assert_eq!(permutation.degree(), image.len());
    }

    #[rstest]
    fn test_permutation_identity() {
        assert_eq!(Permutation::identity(4).image(), vec![0, 1, 2, 3]);
    }

    #[rstest]
    #[case(AstStereoTerm::Lit(1))]
    #[case(AstStereoTerm::LitSet(BTreeSet::from([0, 2])))]
    #[case(AstStereoTerm::Var(Box::new(("x".to_string(), None))))]
    #[case(AstStereoTerm::Var(Box::new(("y".to_string(), Some(BTreeSet::from([0, 1]))))))]
    #[case(AstStereoTerm::Swap(Box::new(AstStereoTerm::Lit(0))))]
    #[case(AstStereoTerm::Mirror(Box::new(AstStereoTerm::Lit(0))))]
    #[case(AstStereoTerm::Apply(Box::new(AstStereoTerm::Lit(0)), PermPermutation::from_image(4, &[1, 0, 2, 3])))]
    fn test_stereo_term_roundtrip(#[case] ast: AstStereoTerm) {
        Python::attach(|py| {
            assert_eq!(StereoTerm::from_ast(py, &ast).unwrap().to_ast(py), ast);
        });
    }

    #[rstest]
    #[case(AstStereoCosetAst::Undetermined)]
    #[case(AstStereoCosetAst::Lit(1))]
    #[case(AstStereoCosetAst::LitSet(BTreeSet::from([0, 1])))]
    #[case(AstStereoCosetAst::Term(Box::new(AstStereoTerm::Lit(1))))]
    fn test_stereo_coset_ast_roundtrip(#[case] ast: AstStereoCosetAst) {
        Python::attach(|py| {
            assert_eq!(StereoCosetAst::from_ast(py, &ast).unwrap().to_ast(py), ast);
        });
    }

    #[rstest]
    #[case(AstTetrahedralStereoAst::Undetermined)]
    #[case(AstTetrahedralStereoAst::NotStereo)]
    #[case(AstTetrahedralStereoAst::Stereo(AstStereoCosetAst::Lit(1)))]
    #[case(AstTetrahedralStereoAst::Stereo(AstStereoCosetAst::Term(Box::new(
        AstStereoTerm::Lit(0)
    ))))]
    fn test_tetrahedral_stereo_ast_roundtrip(#[case] ast: AstTetrahedralStereoAst) {
        Python::attach(|py| {
            assert_eq!(
                TetrahedralStereoAst::from_ast(py, &ast).unwrap().to_ast(py),
                ast
            );
        });
    }

    #[rstest]
    #[case(TetrahedralStereo::Ccw, AstStereoCosetAst::Lit(0))]
    #[case(TetrahedralStereo::Cw, AstStereoCosetAst::Lit(1))]
    fn test_tetrahedral_stereo_to_ast(
        #[case] config: TetrahedralStereo,
        #[case] coset: AstStereoCosetAst,
    ) {
        assert_eq!(config.to_ast(), AstTetrahedralStereoAst::Stereo(coset));
    }

    #[rstest]
    #[case(AstCisTransStereoAst::Undetermined)]
    #[case(AstCisTransStereoAst::NotStereo)]
    #[case(AstCisTransStereoAst::Stereo(AstStereoCosetAst::Lit(1)))]
    #[case(AstCisTransStereoAst::Stereo(AstStereoCosetAst::Term(Box::new(AstStereoTerm::Lit(0)))))]
    fn test_cis_trans_stereo_ast_roundtrip(#[case] ast: AstCisTransStereoAst) {
        Python::attach(|py| {
            assert_eq!(
                CisTransStereoAst::from_ast(py, &ast).unwrap().to_ast(py),
                ast
            );
        });
    }

    #[rstest]
    #[case(CisTransStereo::Z, AstStereoCosetAst::Lit(0))]
    #[case(CisTransStereo::E, AstStereoCosetAst::Lit(1))]
    fn test_cis_trans_stereo_to_ast(
        #[case] config: CisTransStereo,
        #[case] coset: AstStereoCosetAst,
    ) {
        assert_eq!(config.to_ast(), AstCisTransStereoAst::Stereo(coset));
    }

    #[rstest]
    #[case(AstStereoKind::Tetrahedral)]
    #[case(AstStereoKind::CisTrans)]
    #[case(AstStereoKind::Axial)]
    #[case(AstStereoKind::SquarePlanar)]
    #[case(AstStereoKind::TrigonalBipyramidal)]
    #[case(AstStereoKind::Octahedral)]
    fn test_stereo_kind_roundtrip(#[case] ast: AstStereoKind) {
        assert_eq!(StereoKind::from_ast(ast).to_ast(), ast);
    }

    #[rstest]
    #[case(AstStereoLigandKind::Atom)]
    #[case(AstStereoLigandKind::ImplicitHydrogen)]
    #[case(AstStereoLigandKind::LonePair)]
    fn test_stereo_ligand_kind_roundtrip(#[case] ast: AstStereoLigandKind) {
        assert_eq!(StereoLigandKind::from_ast(ast).to_ast(), ast);
    }

    #[rstest]
    #[case(AstTopicity::Homotopic)]
    #[case(AstTopicity::Enantiotopic)]
    #[case(AstTopicity::Diastereotopic)]
    fn test_topicity_roundtrip(#[case] ast: AstTopicity) {
        assert_eq!(Topicity::from_ast(ast).to_ast(), ast);
    }

    #[rstest]
    #[case(AstStereogenicity::Symmetric)]
    #[case(AstStereogenicity::Prochiral)]
    #[case(AstStereogenicity::Stereogenic)]
    fn test_stereogenicity_roundtrip(#[case] ast: AstStereogenicity) {
        assert_eq!(Stereogenicity::from_ast(ast).to_ast(), ast);
    }

    #[rstest]
    fn test_stereo_ligand_new() {
        let ligand = StereoLigand::new(3, StereoLigandKind::ImplicitHydrogen);
        assert_eq!(ligand.atom_id, 3);
        assert_eq!(ligand.kind, StereoLigandKind::ImplicitHydrogen);
        assert_eq!(
            ligand.__repr__(),
            "StereoLigand(atom_id=3, kind=StereoLigandKind.ImplicitHydrogen)"
        );
    }

    #[rstest]
    #[case(AstStereoLigand::new(AstAtomId(0), AstStereoLigandKind::Atom))]
    #[case(AstStereoLigand::new(AstAtomId(5), AstStereoLigandKind::LonePair))]
    fn test_stereo_ligand_roundtrip(#[case] ast: AstStereoLigand) {
        assert_eq!(StereoLigand::from_ast(ast).to_ast(), ast);
    }

    #[rstest]
    fn test_stereo_configuration_ast_roundtrip() {
        Python::attach(|py| {
            for ast in [
                AstStereoConfigurationAst::Undetermined,
                AstStereoConfigurationAst::kinded(
                    AstStereoKind::Tetrahedral,
                    AstStereoCosetAst::Lit(1),
                ),
                AstStereoConfigurationAst::kinded(
                    AstStereoKind::Octahedral,
                    AstStereoCosetAst::Undetermined,
                ),
            ] {
                assert_eq!(
                    StereoConfigurationAst::from_ast(py, &ast)
                        .unwrap()
                        .to_ast(py),
                    ast
                );
            }
        });
    }

    #[rstest]
    fn test_stereo_configuration_ast_kind_coset() {
        Python::attach(|py| {
            let coset = into_py_variant(py, StereoCosetAst::Lit(1)).unwrap();
            let config = StereoConfigurationAst::Kinded(StereoKind::Tetrahedral, coset);
            assert_eq!(config.kind(), Some(StereoKind::Tetrahedral));
            assert_eq!(
                config.coset(py).unwrap().bind(py).borrow().to_ast(py),
                AstStereoCosetAst::Lit(1)
            );
            let undetermined = StereoConfigurationAst::Undetermined();
            assert_eq!(undetermined.kind(), None);
            assert!(undetermined.coset(py).is_none());
        });
    }

    #[rstest]
    fn test_stereo_configuration_arg_to_ast() {
        Python::attach(|py| {
            // the Th shorthand → Kinded(Tetrahedral, coset)
            assert_eq!(
                StereoConfigurationArg::Tetrahedral(TetrahedralStereo::Cw).to_ast(py),
                AstStereoConfigurationAst::kinded(
                    AstStereoKind::Tetrahedral,
                    AstStereoCosetAst::Lit(1)
                )
            );
            // the Ct shorthand → Kinded(CisTrans, coset)
            assert_eq!(
                StereoConfigurationArg::CisTrans(CisTransStereo::E).to_ast(py),
                AstStereoConfigurationAst::kinded(
                    AstStereoKind::CisTrans,
                    AstStereoCosetAst::Lit(1)
                )
            );
            // a StereoConfigurationAst passes through
            let config = Py::new(py, StereoConfigurationAst::Undetermined()).unwrap();
            assert_eq!(
                StereoConfigurationArg::Ast(config).to_ast(py),
                AstStereoConfigurationAst::Undetermined
            );
        });
    }

    #[rstest]
    #[case(PermOrientation::Proper)]
    #[case(PermOrientation::Improper)]
    fn test_orientation_roundtrip(#[case] ast: PermOrientation) {
        assert_eq!(Orientation::from_ast(ast).to_ast(), ast);
    }

    #[rstest]
    fn test_ligand_permutation_new() {
        let ligand_permutation = LigandPermutation::new(Permutation::new(vec![1, 0, 2, 3]));
        assert_eq!(ligand_permutation.permutation.image(), vec![1, 0, 2, 3]);
        assert_eq!(
            ligand_permutation.__repr__(),
            "LigandPermutation([1, 0, 2, 3])"
        );
    }

    #[rstest]
    #[case(AstLigandPermutation(PermPermutation::identity(4)))]
    #[case(AstLigandPermutation(PermPermutation::from_image(4, &[1, 0, 2, 3])))]
    fn test_ligand_permutation_roundtrip(#[case] ast: AstLigandPermutation) {
        assert_eq!(LigandPermutation::from_ast(ast).to_ast(), ast);
    }

    #[rstest]
    #[case::equal(vec![1, 0, 2, 3], vec![1, 0, 2, 3], true)]
    #[case::different(vec![1, 0, 2, 3], vec![0, 1, 2, 3], false)]
    fn test_ligand_permutation_matches(
        #[case] a: Vec<u32>,
        #[case] b: Vec<u32>,
        #[case] expected: bool,
    ) {
        let a = LigandPermutation::new(Permutation::new(a));
        let b = LigandPermutation::new(Permutation::new(b));
        assert_eq!(a.matches(&b), expected);
    }

    #[rstest]
    #[case::equal(vec![1, 0, 2, 3], Orientation::Proper, vec![1, 0, 2, 3], Orientation::Proper, true)]
    #[case::different_orientation(vec![1, 0, 2, 3], Orientation::Proper, vec![1, 0, 2, 3], Orientation::Improper, false)]
    #[case::different_permutation(vec![1, 0, 2, 3], Orientation::Proper, vec![0, 1, 2, 3], Orientation::Proper, false)]
    fn test_oriented_ligand_permutation_matches(
        #[case] a_permutation: Vec<u32>,
        #[case] a_orientation: Orientation,
        #[case] b_permutation: Vec<u32>,
        #[case] b_orientation: Orientation,
        #[case] expected: bool,
    ) {
        let a = OrientedLigandPermutation::new(
            LigandPermutation::new(Permutation::new(a_permutation)),
            a_orientation,
        );
        let b = OrientedLigandPermutation::new(
            LigandPermutation::new(Permutation::new(b_permutation)),
            b_orientation,
        );
        assert_eq!(a.matches(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case(AstOrientedLigandPermutation { permutation: AstLigandPermutation(PermPermutation::from_image(4, &[1, 0, 2, 3])), orientation: PermOrientation::Proper })]
    #[case(AstOrientedLigandPermutation { permutation: AstLigandPermutation(PermPermutation::identity(4)), orientation: PermOrientation::Improper })]
    fn test_oriented_ligand_permutation_roundtrip(#[case] ast: AstOrientedLigandPermutation) {
        assert_eq!(OrientedLigandPermutation::from_ast(ast).to_ast(), ast);
    }

    #[rstest]
    #[case::ordered(1, 2, 1, 2)]
    #[case::reversed(2, 1, 1, 2)]
    #[case::equal(3, 3, 3, 3)]
    fn test_stereo_ligand_pair_new(
        #[case] a: u32,
        #[case] b: u32,
        #[case] first: u32,
        #[case] second: u32,
    ) {
        let pair = StereoLigandPair::new(a, b);
        assert_eq!(pair.first, first);
        assert_eq!(pair.second, second);
    }

    #[rstest]
    #[case(AstStereoLigandPair::new(AstStereoLigandPosition(0), AstStereoLigandPosition(3)))]
    #[case(AstStereoLigandPair::new(AstStereoLigandPosition(2), AstStereoLigandPosition(1)))]
    fn test_stereo_ligand_pair_roundtrip(#[case] ast: AstStereoLigandPair) {
        assert_eq!(StereoLigandPair::from_ast(ast).to_ast(), ast);
    }

    #[rstest]
    #[case::lit(
        TopicityRelationAst::Lit(Topicity::Homotopic),
        Some(Topicity::Homotopic)
    )]
    #[case::undetermined(TopicityRelationAst::Undetermined(), None)]
    #[case::set(
        TopicityRelationAst::LitSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic])),
        None
    )]
    fn test_topicity_relation_ast_as_lit(
        #[case] relation: TopicityRelationAst,
        #[case] expected: Option<Topicity>,
    ) {
        assert_eq!(relation.as_lit(), expected);
    }

    #[rstest]
    #[case(AstTopicityRelationAst::Undetermined)]
    #[case(AstTopicityRelationAst::Lit(AstTopicity::Homotopic))]
    #[case(AstTopicityRelationAst::LitSet(BTreeSet::from([
        AstTopicity::Homotopic,
        AstTopicity::Enantiotopic,
    ])))]
    #[case(AstTopicityRelationAst::NotSet(BTreeSet::from([AstTopicity::Diastereotopic])))]
    fn test_topicity_relation_ast_roundtrip(#[case] ast: AstTopicityRelationAst) {
        assert_eq!(TopicityRelationAst::from_ast(&ast).to_ast(), ast);
    }

    #[rstest]
    #[case::lit(
        StereogenicityAst::Lit(Stereogenicity::Prochiral),
        Some(Stereogenicity::Prochiral)
    )]
    #[case::undetermined(StereogenicityAst::Undetermined(), None)]
    #[case::set(
        StereogenicityAst::LitSet(BTreeSet::from([Stereogenicity::Symmetric])),
        None
    )]
    fn test_stereogenicity_ast_as_lit(
        #[case] relation: StereogenicityAst,
        #[case] expected: Option<Stereogenicity>,
    ) {
        assert_eq!(relation.as_lit(), expected);
    }

    #[rstest]
    #[case(AstStereogenicityAst::Undetermined)]
    #[case(AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic))]
    #[case(AstStereogenicityAst::LitSet(BTreeSet::from([
        AstStereogenicity::Symmetric,
        AstStereogenicity::Prochiral,
    ])))]
    #[case(AstStereogenicityAst::NotSet(BTreeSet::from([AstStereogenicity::Stereogenic])))]
    fn test_stereogenicity_ast_roundtrip(#[case] ast: AstStereogenicityAst) {
        assert_eq!(StereogenicityAst::from_ast(&ast).to_ast(), ast);
    }

    #[rstest]
    fn test_ligand_symmetry_ast_new() {
        Python::attach(|py| {
            let permutation = OrientedLigandPermutation::new(
                LigandPermutation::new(Permutation::new(vec![1, 0, 2, 3])),
                Orientation::Proper,
            );
            let value = LigandSymmetryAst::new(py, permutation, BooleanArg::Lit(true)).unwrap();
            assert!(value.permutation() == permutation);
            assert_eq!(
                value.invariant.bind(py).borrow().to_ast(),
                AstBooleanAst::Lit(true)
            );
            assert_eq!(
                value.__repr__(py).unwrap(),
                "LigandSymmetryAst(OrientedLigandPermutation(permutation=LigandPermutation([1, 0, 2, 3]), orientation=Orientation.Proper), BooleanAst.Lit(True))"
            );
        });
    }

    #[rstest]
    fn test_ligand_symmetry_ast_matches() {
        Python::attach(|py| {
            let permutation = OrientedLigandPermutation::new(
                LigandPermutation::new(Permutation::new(vec![1, 0, 2, 3])),
                Orientation::Proper,
            );
            let other_permutation = OrientedLigandPermutation::new(
                LigandPermutation::new(Permutation::new(vec![0, 1, 2, 3])),
                Orientation::Proper,
            );
            let wildcard = LigandSymmetryAst {
                permutation,
                invariant: into_py_variant(py, BooleanAst::Undetermined()).unwrap(),
            };
            let invariant_true = LigandSymmetryAst {
                permutation,
                invariant: into_py_variant(py, BooleanAst::Lit(true)).unwrap(),
            };
            let invariant_false = LigandSymmetryAst {
                permutation,
                invariant: into_py_variant(py, BooleanAst::Lit(false)).unwrap(),
            };
            let other = LigandSymmetryAst {
                permutation: other_permutation,
                invariant: into_py_variant(py, BooleanAst::Lit(true)).unwrap(),
            };
            assert!(wildcard.matches(&invariant_true, py));
            assert!(!invariant_true.matches(&invariant_false, py));
            assert!(!invariant_true.matches(&other, py));
        });
    }

    #[rstest]
    fn test_ligand_symmetry_ast_roundtrip() {
        Python::attach(|py| {
            for ast in [
                AstLigandSymmetryAst {
                    permutation: AstOrientedLigandPermutation {
                        permutation: AstLigandPermutation(PermPermutation::from_image(
                            4,
                            &[1, 0, 2, 3],
                        )),
                        orientation: PermOrientation::Proper,
                    },
                    invariant: AstBooleanAst::Lit(true),
                },
                AstLigandSymmetryAst {
                    permutation: AstOrientedLigandPermutation {
                        permutation: AstLigandPermutation(PermPermutation::identity(4)),
                        orientation: PermOrientation::Improper,
                    },
                    invariant: AstBooleanAst::Undetermined,
                },
            ] {
                assert_eq!(
                    LigandSymmetryAst::from_ast(py, &ast).unwrap().to_ast(py),
                    ast
                );
            }
        });
    }

    #[rstest]
    fn test_fluxionality_ast_new() {
        Python::attach(|py| {
            let permutation = LigandPermutation::new(Permutation::new(vec![1, 0, 2, 3]));
            let value = FluxionalityAst::new(py, permutation, BooleanArg::Lit(false)).unwrap();
            assert!(value.permutation() == permutation);
            assert_eq!(
                value.active.bind(py).borrow().to_ast(),
                AstBooleanAst::Lit(false)
            );
            assert_eq!(
                value.__repr__(py).unwrap(),
                "FluxionalityAst(LigandPermutation([1, 0, 2, 3]), BooleanAst.Lit(False))"
            );
        });
    }

    #[rstest]
    fn test_fluxionality_ast_matches() {
        Python::attach(|py| {
            let permutation = LigandPermutation::new(Permutation::new(vec![1, 0, 2, 3]));
            let other_permutation = LigandPermutation::new(Permutation::new(vec![0, 1, 2, 3]));
            let wildcard = FluxionalityAst {
                permutation,
                active: into_py_variant(py, BooleanAst::Undetermined()).unwrap(),
            };
            let active_true = FluxionalityAst {
                permutation,
                active: into_py_variant(py, BooleanAst::Lit(true)).unwrap(),
            };
            let active_false = FluxionalityAst {
                permutation,
                active: into_py_variant(py, BooleanAst::Lit(false)).unwrap(),
            };
            let other = FluxionalityAst {
                permutation: other_permutation,
                active: into_py_variant(py, BooleanAst::Lit(true)).unwrap(),
            };
            assert!(wildcard.matches(&active_true, py));
            assert!(!active_true.matches(&active_false, py));
            assert!(!active_true.matches(&other, py));
        });
    }

    #[rstest]
    fn test_fluxionality_ast_roundtrip() {
        Python::attach(|py| {
            for ast in [
                AstFluxionalityAst {
                    permutation: AstLigandPermutation(PermPermutation::from_image(
                        4,
                        &[1, 0, 2, 3],
                    )),
                    active: AstBooleanAst::Lit(false),
                },
                AstFluxionalityAst {
                    permutation: AstLigandPermutation(PermPermutation::identity(4)),
                    active: AstBooleanAst::Undetermined,
                },
            ] {
                assert_eq!(FluxionalityAst::from_ast(py, &ast).unwrap().to_ast(py), ast);
            }
        });
    }

    #[rstest]
    fn test_topicity_ast_new() {
        Python::attach(|py| {
            let pair = StereoLigandPair::new(0, 2);
            let value =
                TopicityAst::new(py, pair, TopicityRelationArg::Lit(Topicity::Homotopic)).unwrap();
            assert!(value.pair() == pair);
            assert_eq!(
                value.relation.bind(py).borrow().to_ast(),
                AstTopicityRelationAst::Lit(AstTopicity::Homotopic)
            );
            assert_eq!(
                value.__repr__(py).unwrap(),
                "TopicityAst(StereoLigandPair(0, 2), TopicityRelationAst.Lit(Topicity.Homotopic))"
            );
        });
    }

    #[rstest]
    fn test_topicity_ast_matches() {
        Python::attach(|py| {
            let pair = StereoLigandPair::new(0, 2);
            let other_pair = StereoLigandPair::new(1, 3);
            let wildcard = TopicityAst {
                pair,
                relation: into_py_variant(py, TopicityRelationAst::Undetermined()).unwrap(),
            };
            let homotopic = TopicityAst {
                pair,
                relation: into_py_variant(py, TopicityRelationAst::Lit(Topicity::Homotopic))
                    .unwrap(),
            };
            let enantiotopic = TopicityAst {
                pair,
                relation: into_py_variant(py, TopicityRelationAst::Lit(Topicity::Enantiotopic))
                    .unwrap(),
            };
            let other = TopicityAst {
                pair: other_pair,
                relation: into_py_variant(py, TopicityRelationAst::Lit(Topicity::Homotopic))
                    .unwrap(),
            };
            assert!(wildcard.matches(&homotopic, py));
            assert!(!homotopic.matches(&enantiotopic, py));
            assert!(!homotopic.matches(&other, py));
        });
    }

    #[rstest]
    fn test_topicity_ast_roundtrip() {
        Python::attach(|py| {
            for ast in [
                AstTopicityAst {
                    pair: AstStereoLigandPair::new(
                        AstStereoLigandPosition(0),
                        AstStereoLigandPosition(2),
                    ),
                    relation: AstTopicityRelationAst::Lit(AstTopicity::Homotopic),
                },
                AstTopicityAst {
                    pair: AstStereoLigandPair::new(
                        AstStereoLigandPosition(1),
                        AstStereoLigandPosition(3),
                    ),
                    relation: AstTopicityRelationAst::Undetermined,
                },
            ] {
                assert_eq!(TopicityAst::from_ast(py, &ast).unwrap().to_ast(py), ast);
            }
        });
    }

    #[rustfmt::skip]
    #[rstest]
    #[case(AstStereoAtomConstraintAst::LigandSymmetry(AstLigandSymmetryAst { permutation: AstOrientedLigandPermutation { permutation: AstLigandPermutation(PermPermutation::from_image(4, &[1, 0, 2, 3])), orientation: PermOrientation::Proper }, invariant: AstBooleanAst::Lit(true) }))]
    #[case(AstStereoAtomConstraintAst::Fluxionality(AstFluxionalityAst { permutation: AstLigandPermutation(PermPermutation::identity(4)), active: AstBooleanAst::Lit(false) }))]
    #[case(AstStereoAtomConstraintAst::Topicity(AstTopicityAst { pair: AstStereoLigandPair::new(AstStereoLigandPosition(0), AstStereoLigandPosition(1)), relation: AstTopicityRelationAst::Lit(AstTopicity::Homotopic) }))]
    #[case(AstStereoAtomConstraintAst::Stereogenicity(AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic)))]
    fn test_stereo_atom_constraint_ast_roundtrip(#[case] ast: AstStereoAtomConstraintAst) {
        Python::attach(|py| {
            assert_eq!(
                StereoAtomConstraintAst::from_ast(py, &ast).unwrap().to_ast(py),
                ast
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_constraint_ast_key() {
        Python::attach(|py| {
            let ast = AstStereoAtomConstraintAst::Topicity(AstTopicityAst {
                pair: AstStereoLigandPair::new(
                    AstStereoLigandPosition(0),
                    AstStereoLigandPosition(1),
                ),
                relation: AstTopicityRelationAst::Lit(AstTopicity::Homotopic),
            });
            let key = StereoAtomConstraintAst::from_ast(py, &ast)
                .unwrap()
                .key(py)
                .unwrap();
            assert_eq!(
                key.to_ast(py),
                AstStereoAtomConstraintKey::Topicity(AstStereoLigandPair::new(
                    AstStereoLigandPosition(0),
                    AstStereoLigandPosition(1),
                ))
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_constraints_ast_get() {
        Python::attach(|py| {
            let stereogenicity = AstStereoAtomConstraintAst::Stereogenicity(
                AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic),
            );
            let mut ast_cs = AstStereoAtomConstraintsAst::new();
            ast_cs.extend([stereogenicity.clone()]);
            let constraints = StereoAtomConstraintsAst::from_inner(ast_cs);
            assert_eq!(constraints.__len__(), 1);

            let present = into_py_variant(py, StereoAtomConstraintKey::Stereogenicity()).unwrap();
            assert!(constraints.__contains__(py, present.clone_ref(py)));
            assert_eq!(
                constraints
                    .__getitem__(py, present.clone_ref(py))
                    .unwrap()
                    .to_ast(py),
                stereogenicity
            );

            let absent = into_py_variant(
                py,
                StereoAtomConstraintKey::Topicity(
                    into_py_variant(py, StereoLigandPair::new(0, 1)).unwrap(),
                ),
            )
            .unwrap();
            assert!(!constraints.__contains__(py, absent.clone_ref(py)));
            assert!(constraints.__getitem__(py, absent).is_err());
        });
    }

    #[rstest]
    fn test_stereo_atom_constraints_ast_set_pop() {
        Python::attach(|py| {
            let stereogenicity = into_py_variant(
                py,
                StereoAtomConstraintAst::from_ast(
                    py,
                    &AstStereoAtomConstraintAst::Stereogenicity(AstStereogenicityAst::Lit(
                        AstStereogenicity::Stereogenic,
                    )),
                )
                .unwrap(),
            )
            .unwrap();
            let mut constraints = StereoAtomConstraintsAst::new(py, Vec::new());
            constraints.set(py, stereogenicity);
            assert_eq!(constraints.__len__(), 1);

            let key = into_py_variant(py, StereoAtomConstraintKey::Stereogenicity()).unwrap();
            let popped = constraints.pop(py, key.clone_ref(py)).unwrap();
            assert_eq!(
                popped.unwrap().to_ast(py),
                AstStereoAtomConstraintAst::Stereogenicity(AstStereogenicityAst::Lit(
                    AstStereogenicity::Stereogenic
                ))
            );
            assert_eq!(constraints.__len__(), 0);
            assert!(constraints.pop(py, key).unwrap().is_none());
        });
    }

    #[rstest]
    fn test_stereo_atom_constraints_ast_accessors() {
        Python::attach(|py| {
            let mut ast_cs = AstStereoAtomConstraintsAst::new();
            ast_cs.extend([
                AstStereoAtomConstraintAst::LigandSymmetry(AstLigandSymmetryAst {
                    permutation: AstOrientedLigandPermutation {
                        permutation: AstLigandPermutation(PermPermutation::from_image(
                            4,
                            &[1, 0, 2, 3],
                        )),
                        orientation: PermOrientation::Proper,
                    },
                    invariant: AstBooleanAst::Lit(true),
                }),
                AstStereoAtomConstraintAst::Topicity(AstTopicityAst {
                    pair: AstStereoLigandPair::new(
                        AstStereoLigandPosition(0),
                        AstStereoLigandPosition(1),
                    ),
                    relation: AstTopicityRelationAst::Lit(AstTopicity::Homotopic),
                }),
                AstStereoAtomConstraintAst::Stereogenicity(AstStereogenicityAst::Lit(
                    AstStereogenicity::Stereogenic,
                )),
            ]);
            let constraints = StereoAtomConstraintsAst::from_inner(ast_cs);

            assert_eq!(
                constraints.stereogenicity().to_ast(),
                AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic)
            );
            assert_eq!(
                constraints.topicity(StereoLigandPair::new(0, 1)).to_ast(),
                AstTopicityRelationAst::Lit(AstTopicity::Homotopic)
            );
            let ligand_symmetries = constraints.ligand_symmetries(py).unwrap();
            assert_eq!(ligand_symmetries.len(), 1);
            assert_eq!(
                ligand_symmetries[0].to_ast(py).invariant,
                AstBooleanAst::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_constraints_ast_iter() {
        Python::attach(|py| {
            let mut ast_cs = AstStereoAtomConstraintsAst::new();
            ast_cs.extend([
                AstStereoAtomConstraintAst::Topicity(AstTopicityAst {
                    pair: AstStereoLigandPair::new(
                        AstStereoLigandPosition(0),
                        AstStereoLigandPosition(1),
                    ),
                    relation: AstTopicityRelationAst::Lit(AstTopicity::Homotopic),
                }),
                AstStereoAtomConstraintAst::Stereogenicity(AstStereogenicityAst::Lit(
                    AstStereogenicity::Stereogenic,
                )),
            ]);
            let constraints = StereoAtomConstraintsAst::from_inner(ast_cs);

            let keys: Vec<AstStereoAtomConstraintKey> = constraints
                .keys(py)
                .unwrap()
                .keys
                .map(|k| k.bind(py).borrow().to_ast(py))
                .collect();
            assert_eq!(
                keys,
                vec![
                    AstStereoAtomConstraintKey::Topicity(AstStereoLigandPair::new(
                        AstStereoLigandPosition(0),
                        AstStereoLigandPosition(1),
                    )),
                    AstStereoAtomConstraintKey::Stereogenicity,
                ]
            );
            let values: Vec<AstStereoAtomConstraintAst> = constraints
                .values(py)
                .unwrap()
                .entries
                .map(|c| c.bind(py).borrow().to_ast(py))
                .collect();
            assert_eq!(values.len(), 2);
            assert_eq!(constraints.items(py).unwrap().items.count(), 2);
        });
    }

    #[rstest]
    fn test_stereo_atom_constraints_ast_update() {
        Python::attach(|py| {
            let base = Py::new(py, StereoAtomConstraintsAst::new(py, Vec::new())).unwrap();
            let entry = into_py_variant(
                py,
                StereoAtomConstraintAst::from_ast(
                    py,
                    &AstStereoAtomConstraintAst::Stereogenicity(AstStereogenicityAst::Lit(
                        AstStereogenicity::Stereogenic,
                    )),
                )
                .unwrap(),
            )
            .unwrap();
            StereoAtomConstraintsAst::update(
                base.clone_ref(py),
                py,
                StereoAtomConstraintsUpdate::Entries(vec![entry]),
            )
            .unwrap();
            assert_eq!(base.borrow(py).__len__(), 1);

            let overlay = Py::new(py, StereoAtomConstraintsAst::new(py, Vec::new())).unwrap();
            StereoAtomConstraintsAst::update(
                overlay.clone_ref(py),
                py,
                StereoAtomConstraintsUpdate::Container(base),
            )
            .unwrap();
            assert_eq!(overlay.borrow(py).__len__(), 1);
        });
    }

    #[rstest]
    fn test_stereo_atom_constraints_arg_to_ast() {
        Python::attach(|py| {
            let entry = into_py_variant(
                py,
                StereoAtomConstraintAst::from_ast(
                    py,
                    &AstStereoAtomConstraintAst::Stereogenicity(AstStereogenicityAst::Lit(
                        AstStereogenicity::Stereogenic,
                    )),
                )
                .unwrap(),
            )
            .unwrap();
            let container = Py::new(py, StereoAtomConstraintsAst::new(py, vec![entry])).unwrap();
            let arg = StereoAtomConstraintsArg::Container(container);
            let mut expected = AstStereoAtomConstraintsAst::new();
            expected.extend([AstStereoAtomConstraintAst::Stereogenicity(
                AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic),
            )]);
            assert_eq!(arg.to_ast(py), expected);
        });
    }

    // `StereoBondConstraintsAst` is the second `stereo_constraints!` instantiation; the shared
    // macro is covered by the `StereoAtom` tests above. This confirms the bond instantiation
    // and exercises its `from_inner` / `Arg::to_ast`.
    #[rstest]
    fn test_stereo_bond_constraints_ast() {
        Python::attach(|py| {
            let stereogenicity = AstStereoBondConstraintAst::Stereogenicity(
                AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic),
            );
            let mut ast_cs = AstStereoBondConstraintsAst::new();
            ast_cs.extend([stereogenicity.clone()]);
            let constraints = StereoBondConstraintsAst::from_inner(ast_cs);
            assert_eq!(constraints.__len__(), 1);
            assert_eq!(
                constraints.stereogenicity().to_ast(),
                AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic)
            );

            let mut container_ast = AstStereoBondConstraintsAst::new();
            container_ast.extend([stereogenicity.clone()]);
            let container =
                Py::new(py, StereoBondConstraintsAst::from_inner(container_ast)).unwrap();
            let arg = StereoBondConstraintsArg::Container(container);
            let mut expected = AstStereoBondConstraintsAst::new();
            expected.extend([stereogenicity]);
            assert_eq!(arg.to_ast(py), expected);
        });
    }
}
