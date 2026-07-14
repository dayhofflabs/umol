//! Stereo values, configurations, owned entity ASTs, and molecule-backed views.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use std::collections::BTreeSet;
use std::str::FromStr;
use std::vec::IntoIter;

use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;
// The `BooleanAst` Rust value is still `#[cfg(test)]` (only tests build it directly); its `to_rust`
// peer is already live.
#[cfg(test)]
use umol_ast::ast::BooleanAst as AstBooleanAst;
use umol_ast::ast::{
    AtomId as AstAtomId, BondId as AstBondId, CisTransStereoAst as AstCisTransStereoAst,
    LigandPermutation as AstLigandPermutation, MoleculeAst as AstMoleculeAst,
    OrientedLigandPermutation as AstOrientedLigandPermutation, StereoAtomAst as AstStereoAtomAst,
    StereoAtomId as AstStereoAtomId, StereoAtomView as AstStereoAtomView,
    StereoBondAst as AstStereoBondAst, StereoBondId as AstStereoBondId,
    StereoBondView as AstStereoBondView, StereoConfigurationAst as AstStereoConfigurationAst,
    StereoCosetAst as AstStereoCosetAst, StereoKind as AstStereoKind,
    StereoLigand as AstStereoLigand, StereoLigandKind as AstStereoLigandKind,
    StereoLigandPair as AstStereoLigandPair, StereoLigandPosition as AstStereoLigandPosition,
    StereoTerm as AstStereoTerm, Stereogenicity as AstStereogenicity,
    TetrahedralStereoAst as AstTetrahedralStereoAst, Topicity as AstTopicity,
};
use umol_perm::{Orientation as PermOrientation, Permutation as PermPermutation};

use crate::convert::{hash_rust, into_py_variant, variant_repr};
use crate::error::parse_error;
use crate::molecule::MoleculeAst;

/// A permutation of `0..degree` in one-line (image) notation.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Permutation(PermPermutation);

#[pymethods]
impl Permutation {
    /// Construct from the image (one-line notation); the degree is the image length.
    #[new]
    fn new(image: Vec<u32>) -> Self {
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
        self.to_rust(py) == other.to_rust(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
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
    pub(crate) fn from_rust(py: Python<'_>, ast: &AstStereoTerm) -> PyResult<Self> {
        Ok(match ast {
            AstStereoTerm::Var(boxed) => {
                let (name, restriction) = &**boxed;
                StereoTerm::Var(name.clone(), restriction.clone())
            }
            AstStereoTerm::Lit(index) => StereoTerm::Lit(*index),
            AstStereoTerm::LitSet(members) => StereoTerm::LitSet(members.clone()),
            AstStereoTerm::Swap(inner) => {
                StereoTerm::Swap(into_py_variant(py, StereoTerm::from_rust(py, inner)?)?)
            }
            AstStereoTerm::Mirror(inner) => {
                StereoTerm::Mirror(into_py_variant(py, StereoTerm::from_rust(py, inner)?)?)
            }
            AstStereoTerm::Apply(inner, permutation) => StereoTerm::Apply(
                into_py_variant(py, StereoTerm::from_rust(py, inner)?)?,
                Permutation::from_inner(*permutation),
            ),
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstStereoTerm {
        match self {
            StereoTerm::Var(name, restriction) => {
                AstStereoTerm::Var(Box::new((name.clone(), restriction.clone())))
            }
            StereoTerm::Lit(index) => AstStereoTerm::Lit(*index),
            StereoTerm::LitSet(members) => AstStereoTerm::LitSet(members.clone()),
            StereoTerm::Swap(inner) => {
                AstStereoTerm::Swap(Box::new(inner.bind(py).borrow().to_rust(py)))
            }
            StereoTerm::Mirror(inner) => {
                AstStereoTerm::Mirror(Box::new(inner.bind(py).borrow().to_rust(py)))
            }
            StereoTerm::Apply(inner, permutation) => AstStereoTerm::Apply(
                Box::new(inner.bind(py).borrow().to_rust(py)),
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
        self.to_rust(py) == other.to_rust(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
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
    pub(crate) fn from_rust(py: Python<'_>, ast: &AstStereoCosetAst) -> PyResult<Self> {
        Ok(match ast {
            AstStereoCosetAst::Undetermined => Self::Undetermined(),
            AstStereoCosetAst::Lit(index) => Self::Lit(*index),
            AstStereoCosetAst::LitSet(members) => Self::LitSet(members.clone()),
            AstStereoCosetAst::Term(inner) => {
                Self::Term(into_py_variant(py, StereoTerm::from_rust(py, inner)?)?)
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstStereoCosetAst {
        match self {
            Self::Undetermined() => AstStereoCosetAst::Undetermined,
            Self::Lit(index) => AstStereoCosetAst::Lit(*index),
            Self::LitSet(members) => AstStereoCosetAst::LitSet(members.clone()),
            Self::Term(inner) => {
                AstStereoCosetAst::Term(Box::new(inner.bind(py).borrow().to_rust(py)))
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
        self.to_rust(py) == other.to_rust(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
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
    pub(crate) fn from_rust(py: Python<'_>, ast: &AstTetrahedralStereoAst) -> PyResult<Self> {
        Ok(match ast {
            AstTetrahedralStereoAst::Undetermined => Self::Undetermined(),
            AstTetrahedralStereoAst::NotStereo => Self::NotStereo(),
            AstTetrahedralStereoAst::Stereo(coset) => {
                Self::Stereo(into_py_variant(py, StereoCosetAst::from_rust(py, coset)?)?)
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstTetrahedralStereoAst {
        match self {
            Self::Undetermined() => AstTetrahedralStereoAst::Undetermined,
            Self::NotStereo() => AstTetrahedralStereoAst::NotStereo,
            Self::Stereo(coset) => {
                AstTetrahedralStereoAst::Stereo(coset.bind(py).borrow().to_rust(py))
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
    pub(crate) fn to_rust(self) -> AstTetrahedralStereoAst {
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
        self.to_rust(py) == other.to_rust(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
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
    pub(crate) fn from_rust(py: Python<'_>, ast: &AstCisTransStereoAst) -> PyResult<Self> {
        Ok(match ast {
            AstCisTransStereoAst::Undetermined => Self::Undetermined(),
            AstCisTransStereoAst::NotStereo => Self::NotStereo(),
            AstCisTransStereoAst::Stereo(coset) => {
                Self::Stereo(into_py_variant(py, StereoCosetAst::from_rust(py, coset)?)?)
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstCisTransStereoAst {
        match self {
            Self::Undetermined() => AstCisTransStereoAst::Undetermined,
            Self::NotStereo() => AstCisTransStereoAst::NotStereo,
            Self::Stereo(coset) => {
                AstCisTransStereoAst::Stereo(coset.bind(py).borrow().to_rust(py))
            }
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
    pub(crate) fn to_rust(self) -> AstCisTransStereoAst {
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
    pub(crate) fn to_rust(&self, py: Python<'_>) -> PyResult<AstCisTransStereoAst> {
        Ok(match self {
            CisTransStereoArg::Flag(false) => AstCisTransStereoAst::NotStereo,
            CisTransStereoArg::Flag(true) => {
                return Err(PyValueError::new_err(
                    "cis_trans_stereo = True is not meaningful; use CisTransStereo.Z/E or False",
                ))
            }
            CisTransStereoArg::Config(cts) => cts.to_rust(),
            CisTransStereoArg::Ast(a) => a.bind(py).borrow().to_rust(py),
        })
    }
}

/// The coordination geometry of a stereo site. A fieldless, hashable value enum whose
/// members correspond exactly to the Rust `StereoKind`.
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
    pub(crate) fn from_rust(ast: AstStereoKind) -> Self {
        match ast {
            AstStereoKind::Tetrahedral => Self::Tetrahedral,
            AstStereoKind::CisTrans => Self::CisTrans,
            AstStereoKind::Axial => Self::Axial,
            AstStereoKind::SquarePlanar => Self::SquarePlanar,
            AstStereoKind::TrigonalBipyramidal => Self::TrigonalBipyramidal,
            AstStereoKind::Octahedral => Self::Octahedral,
        }
    }

    pub(crate) fn to_rust(self) -> AstStereoKind {
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
/// lone pair). A fieldless, hashable value enum corresponding to the Rust `StereoLigandKind`.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum StereoLigandKind {
    Atom,
    ImplicitHydrogen,
    LonePair,
}

impl StereoLigandKind {
    pub(crate) fn from_rust(ast: AstStereoLigandKind) -> Self {
        match ast {
            AstStereoLigandKind::Atom => Self::Atom,
            AstStereoLigandKind::ImplicitHydrogen => Self::ImplicitHydrogen,
            AstStereoLigandKind::LonePair => Self::LonePair,
        }
    }

    pub(crate) fn to_rust(self) -> AstStereoLigandKind {
        match self {
            Self::Atom => AstStereoLigandKind::Atom,
            Self::ImplicitHydrogen => AstStereoLigandKind::ImplicitHydrogen,
            Self::LonePair => AstStereoLigandKind::LonePair,
        }
    }
}

/// Topicity of two ligand positions of a stereo carrier (a derived ground classification).
/// A fieldless, hashable value enum corresponding to the Rust `Topicity`. `Ord` lets it key the
/// `BTreeSet` in the `TopicityRelationAst` set variants.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Topicity {
    Homotopic,
    Enantiotopic,
    Diastereotopic,
}

impl Topicity {
    pub(crate) fn from_rust(ast: AstTopicity) -> Self {
        match ast {
            AstTopicity::Homotopic => Self::Homotopic,
            AstTopicity::Enantiotopic => Self::Enantiotopic,
            AstTopicity::Diastereotopic => Self::Diastereotopic,
        }
    }

    pub(crate) fn to_rust(self) -> AstTopicity {
        match self {
            Self::Homotopic => AstTopicity::Homotopic,
            Self::Enantiotopic => AstTopicity::Enantiotopic,
            Self::Diastereotopic => AstTopicity::Diastereotopic,
        }
    }
}

/// Stereogenicity classification of a stereo carrier (a derived ground classification).
/// A fieldless, hashable value enum corresponding to the Rust `Stereogenicity`. `Ord` lets it key
/// the `BTreeSet` in the `StereogenicityAst` set variants.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Stereogenicity {
    Symmetric,
    Prochiral,
    Stereogenic,
}

impl Stereogenicity {
    pub(crate) fn from_rust(ast: AstStereogenicity) -> Self {
        match ast {
            AstStereogenicity::Symmetric => Self::Symmetric,
            AstStereogenicity::Prochiral => Self::Prochiral,
            AstStereogenicity::Stereogenic => Self::Stereogenic,
        }
    }

    pub(crate) fn to_rust(self) -> AstStereogenicity {
        match self {
            Self::Symmetric => AstStereogenicity::Symmetric,
            Self::Prochiral => AstStereogenicity::Prochiral,
            Self::Stereogenic => AstStereogenicity::Stereogenic,
        }
    }
}

/// A stereo ligand occupying a coordination position of a stereo site: the ligand's atom
/// id and its kind. For a virtual ligand (`ImplicitHydrogen`/`LonePair`) the `atom_id` is
/// the bearing atom; the `kind` disambiguates. Immutable value, hashable. Corresponds to the Rust
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
    pub(crate) fn from_rust(ast: AstStereoLigand) -> Self {
        StereoLigand {
            atom_id: ast.atom_id.0,
            kind: StereoLigandKind::from_rust(ast.kind),
        }
    }

    pub(crate) fn to_rust(self) -> AstStereoLigand {
        AstStereoLigand::new(AstAtomId(self.atom_id), self.kind.to_rust())
    }
}

/// A stereo configuration: undetermined (geometry not yet known, so no coset), or `Kinded`
/// — a concrete coordination geometry bound to a coset that may still be open. Corresponds to the
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
        self.to_rust(py) == other.to_rust(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
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
    pub(crate) fn from_rust(py: Python<'_>, ast: &AstStereoConfigurationAst) -> PyResult<Self> {
        Ok(match ast {
            AstStereoConfigurationAst::Undetermined => Self::Undetermined(),
            AstStereoConfigurationAst::Kinded(kind, coset) => Self::Kinded(
                StereoKind::from_rust(*kind),
                into_py_variant(py, StereoCosetAst::from_rust(py, coset)?)?,
            ),
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstStereoConfigurationAst {
        match self {
            Self::Undetermined() => AstStereoConfigurationAst::Undetermined,
            Self::Kinded(kind, coset) => AstStereoConfigurationAst::Kinded(
                kind.to_rust(),
                coset.bind(py).borrow().to_rust(py),
            ),
        }
    }
}

/// Setter coercion for a stereo `configuration` field: the `TetrahedralStereo` (`Ccw`/`Cw`)
/// or `CisTransStereo` (`Z`/`E`) per-kind coset shorthand, or a `StereoConfigurationAst`
/// passthrough. Axial/square-planar/etc. have no shorthand — use the full `Kinded` form.
#[derive(FromPyObject)]
pub(crate) enum StereoConfigurationArg {
    Tetrahedral(TetrahedralStereo),
    CisTrans(CisTransStereo),
    Ast(Py<StereoConfigurationAst>),
}

impl StereoConfigurationArg {
    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstStereoConfigurationAst {
        match self {
            StereoConfigurationArg::Tetrahedral(t) => match t.to_rust() {
                AstTetrahedralStereoAst::Stereo(coset) => {
                    AstStereoConfigurationAst::Kinded(AstStereoKind::Tetrahedral, coset)
                }
                _ => unreachable!("TetrahedralStereo shorthand is always a Stereo coset"),
            },
            StereoConfigurationArg::CisTrans(c) => match c.to_rust() {
                AstCisTransStereoAst::Stereo(coset) => {
                    AstStereoConfigurationAst::Kinded(AstStereoKind::CisTrans, coset)
                }
                _ => unreachable!("CisTransStereo shorthand is always a Stereo coset"),
            },
            StereoConfigurationArg::Ast(a) => a.bind(py).borrow().to_rust(py),
        }
    }
}

/// Orientation grade of a ligand permutation: a proper rotation, or an improper (mirror)
/// operation. A fieldless, hashable value enum corresponding to `umol_perm::Orientation`.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Orientation {
    Proper,
    Improper,
}

impl Orientation {
    pub(crate) fn from_rust(orientation: PermOrientation) -> Self {
        match orientation {
            PermOrientation::Proper => Self::Proper,
            PermOrientation::Improper => Self::Improper,
        }
    }

    pub(crate) fn to_rust(self) -> PermOrientation {
        match self {
            Self::Proper => PermOrientation::Proper,
            Self::Improper => PermOrientation::Improper,
        }
    }
}

/// A permutation of a stereo site's ligand positions (frame-relative). Immutable value,
/// hashable. Corresponds to the Rust `LigandPermutation`.
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

    pub(crate) fn __repr__(&self) -> String {
        format!("LigandPermutation({:?})", self.permutation.image())
    }
}

impl LigandPermutation {
    pub(crate) fn from_rust(ast: AstLigandPermutation) -> Self {
        LigandPermutation {
            permutation: Permutation::from_inner(ast.0),
        }
    }

    pub(crate) fn to_rust(self) -> AstLigandPermutation {
        AstLigandPermutation(self.permutation.inner())
    }
}

/// A ligand permutation carrying a proper/improper grade. Immutable value, hashable. Corresponds
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

    pub(crate) fn __repr__(&self) -> String {
        format!(
            "OrientedLigandPermutation(permutation={}, orientation=Orientation.{:?})",
            self.permutation.__repr__(),
            self.orientation
        )
    }
}

impl OrientedLigandPermutation {
    pub(crate) fn from_rust(ast: AstOrientedLigandPermutation) -> Self {
        OrientedLigandPermutation {
            permutation: LigandPermutation::from_rust(ast.permutation),
            orientation: Orientation::from_rust(ast.orientation),
        }
    }

    pub(crate) fn to_rust(self) -> AstOrientedLigandPermutation {
        AstOrientedLigandPermutation {
            permutation: self.permutation.to_rust(),
            orientation: self.orientation.to_rust(),
        }
    }
}

/// An unordered pair of ligand positions of a stereo site, normalized so the lower position
/// is `first`. Keys a per-pair topicity constraint. Immutable value, hashable. Corresponds to the
/// Rust `StereoLigandPair`.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct StereoLigandPair {
    #[pyo3(get)]
    pub(crate) first: u32,
    #[pyo3(get)]
    pub(crate) second: u32,
}

#[pymethods]
impl StereoLigandPair {
    /// Normalizes the pair so the lower position is `first`.
    #[new]
    pub(crate) fn new(a: u32, b: u32) -> Self {
        StereoLigandPair::from_rust(AstStereoLigandPair::new(
            AstStereoLigandPosition(a),
            AstStereoLigandPosition(b),
        ))
    }

    pub(crate) fn __repr__(&self) -> String {
        format!("StereoLigandPair({}, {})", self.first, self.second)
    }
}

impl StereoLigandPair {
    pub(crate) fn from_rust(ast: AstStereoLigandPair) -> Self {
        StereoLigandPair {
            first: ast.first().0,
            second: ast.second().0,
        }
    }

    pub(crate) fn to_rust(self) -> AstStereoLigandPair {
        AstStereoLigandPair::new(
            AstStereoLigandPosition(self.first),
            AstStereoLigandPosition(self.second),
        )
    }
}

#[cfg(test)]
use crate::constraint::stereo::{
    FluxionalityAst, LigandSymmetryAst, StereoAtomConstraintKey, StereoAtomConstraintsUpdate,
    StereogenicityAst, TopicityAst, TopicityRelationArg, TopicityRelationAst,
};
use crate::constraint::stereo::{
    StereoAtomConstraintAst, StereoAtomConstraintsArg, StereoAtomConstraintsAst,
    StereoAtomConstraintsBacking, StereoAtomConstraintsView, StereoBondConstraintAst,
    StereoBondConstraintsArg, StereoBondConstraintsAst, StereoBondConstraintsBacking,
    StereoBondConstraintsView,
};

/// Per-entity stereo element value pyclass — `StereoAtomAst` / `StereoBondAst`
/// `{configuration, constraints}` — macro-generated for the two stereo entities.
macro_rules! stereo_value {
    (@from_inner production, $value:ident, $ast_value:ident) => {
        /// Wrap an owned Rust stereo-entity AST.
        pub(crate) fn from_inner(value: $ast_value) -> Self {
            $value(value)
        }
    };
    (@from_inner test, $value:ident, $ast_value:ident) => {
        #[cfg(test)]
        pub(crate) fn from_inner(value: $ast_value) -> Self {
            $value(value)
        }
    };
    (
        $value:ident, $ast_value:ident, $constraint:ident, $constraints:ident, $arg:ident,
        $view:ident, $backing:ident, $from_inner:ident $(,)?
    ) => {
        #[pyclass]
        pub struct $value($ast_value);

        #[pymethods]
        impl $value {
            /// Construct from a stereo configuration — a `TetrahedralStereo` / `CisTransStereo`
            /// per-kind shorthand or a `StereoConfigurationAst` — optionally setting constraints.
            #[new]
            #[pyo3(signature = (configuration, *, constraints=None))]
            fn new(
                py: Python<'_>,
                configuration: StereoConfigurationArg,
                constraints: Option<Py<$constraints>>,
            ) -> Self {
                let constraints = constraints
                    .map(|c| c.bind(py).borrow().inner().clone())
                    .unwrap_or_default();
                $value($ast_value {
                    configuration: configuration.to_rust(py),
                    constraints,
                })
            }

            /// Parse a stereo-DSL string (e.g. `"Th0"`) into the value.
            #[staticmethod]
            fn parse(s: &str) -> PyResult<Self> {
                $ast_value::from_str(s).map(Self).map_err(parse_error)
            }

            fn __str__(&self) -> String {
                self.0.to_string()
            }

            fn __repr__(&self) -> String {
                format!("{}.parse('{}')", stringify!($value), self.0)
            }

            /// The stereo configuration (geometry + coset).
            #[getter]
            fn configuration(&self, py: Python<'_>) -> PyResult<StereoConfigurationAst> {
                StereoConfigurationAst::from_rust(py, &self.0.configuration)
            }

            #[setter]
            fn set_configuration(&mut self, py: Python<'_>, value: StereoConfigurationArg) {
                self.0.configuration = value.to_rust(py);
            }

            /// The entity's constraints as a live handle onto this entity: reads borrow the
            /// current state, mutators write through to the entity in place.
            #[getter]
            fn constraints(slf: Py<Self>) -> $view {
                $view {
                    backing: $backing::Value(slf),
                }
            }

            /// Replace the whole constraint set (wipe-and-set) from a value container or a live
            /// view. Snapshots `value` *before* the write borrow, so `x.constraints =
            /// x.constraints` (a view over the same entity) reads while the entity is unborrowed
            /// instead of self-aliasing into a double-borrow panic.
            #[setter]
            fn set_constraints(slf: Py<Self>, py: Python<'_>, value: $arg) -> PyResult<()> {
                let snapshot = value.to_rust(py)?;
                slf.borrow_mut(py).0.constraints = snapshot;
                Ok(())
            }

            /// The fields as a dict: `configuration` plus a `constraints` list of the entries.
            fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
                let dict = PyDict::new(py);
                dict.set_item("configuration", self.configuration(py)?)?;
                let constraints = self
                    .0
                    .constraints
                    .iter()
                    .map(|c| into_py_variant(py, $constraint::from_rust(py, c)?))
                    .collect::<PyResult<Vec<_>>>()?;
                dict.set_item("constraints", constraints)?;
                Ok(dict)
            }
        }

        impl $value {
            /// The wrapped AST entity — read access for the entity-backed constraints view.
            pub(crate) fn inner(&self) -> &$ast_value {
                &self.0
            }

            /// Mutable access to the wrapped AST entity — write access for the view.
            pub(crate) fn inner_mut(&mut self) -> &mut $ast_value {
                &mut self.0
            }

            stereo_value!(@from_inner $from_inner, $value, $ast_value);
        }
    };
}

stereo_value! {
    StereoAtomAst, AstStereoAtomAst, StereoAtomConstraintAst, StereoAtomConstraintsAst,
    StereoAtomConstraintsArg, StereoAtomConstraintsView, StereoAtomConstraintsBacking, production,
}

stereo_value! {
    StereoBondAst, AstStereoBondAst, StereoBondConstraintAst, StereoBondConstraintsAst,
    StereoBondConstraintsArg, StereoBondConstraintsView, StereoBondConstraintsBacking, test,
}

/// Per-entity molecule-embedded stereo view — `StereoAtomView` / `StereoBondView` — a handle
/// to the molecule plus the entity's id. Field reads rebuild the transient Rust view; the
/// molecule is never copied. The site atom/bond and ligands are read-only topology; the
/// configuration and constraints are the mutable value.
macro_rules! stereo_view {
    (
        $view:ident, $ast_view:ident, $ast_id:ident, $namespace:ident, $entity_mut:ident,
        $id_error:literal, $constraint:ident, $constraints_view:ident, $constraints_backing:ident,
        $arg:ident $(,)?
    ) => {
        #[pyclass]
        pub struct $view {
            owner: Py<MoleculeAst>,
            id: $ast_id,
        }

        impl $view {
            /// Rebuild the transient AST view for this entity, or `IndexError` if the id is
            /// no longer present.
            fn view<'a>(&self, molecule: &'a AstMoleculeAst) -> PyResult<$ast_view<'a>> {
                molecule
                    .$namespace()
                    .get(self.id)
                    .ok_or_else(|| PyIndexError::new_err($id_error))
            }
        }

        #[pymethods]
        impl $view {
            #[getter]
            fn id(&self) -> u32 {
                self.id.0
            }

            fn __repr__(&self) -> String {
                format!("{}(id={})", stringify!($view), self.id.0)
            }

            /// The site atom/bond index this stereo entity sits on (read-only topology).
            #[getter]
            fn site_id(&self, py: Python<'_>) -> PyResult<u32> {
                let molecule = self.owner.bind(py).borrow();
                Ok(self.view(molecule.inner())?.site_id().0)
            }

            /// The ligands in frame order (read-only topology).
            #[getter]
            fn ligands(&self, py: Python<'_>) -> PyResult<Vec<StereoLigand>> {
                let molecule = self.owner.bind(py).borrow();
                Ok(self
                    .view(molecule.inner())?
                    .ligand_frame()
                    .into_iter()
                    .map(StereoLigand::from_rust)
                    .collect())
            }

            /// The coordination-geometry kind (from the configuration).
            #[getter]
            fn kind(&self, py: Python<'_>) -> PyResult<StereoKind> {
                let molecule = self.owner.bind(py).borrow();
                Ok(StereoKind::from_rust(self.view(molecule.inner())?.kind()))
            }

            /// The coset (from the configuration).
            #[getter]
            fn coset(&self, py: Python<'_>) -> PyResult<StereoCosetAst> {
                let molecule = self.owner.bind(py).borrow();
                StereoCosetAst::from_rust(py, self.view(molecule.inner())?.coset())
            }

            /// The stereo configuration (geometry + coset).
            #[getter]
            fn configuration(&self, py: Python<'_>) -> PyResult<StereoConfigurationAst> {
                let molecule = self.owner.bind(py).borrow();
                StereoConfigurationAst::from_rust(
                    py,
                    &self.view(molecule.inner())?.ast.configuration,
                )
            }

            #[setter]
            fn set_configuration(&self, py: Python<'_>, value: StereoConfigurationArg) {
                self.owner
                    .borrow_mut(py)
                    .inner_mut()
                    .$entity_mut(self.id)
                    .ast
                    .configuration = value.to_rust(py);
            }

            /// The entity's constraints as a live handle onto the molecule: reads borrow the
            /// current state, mutators write through to the entity in place.
            #[getter]
            fn constraints(&self, py: Python<'_>) -> $constraints_view {
                $constraints_view {
                    backing: $constraints_backing::Molecule {
                        owner: self.owner.clone_ref(py),
                        id: self.id,
                    },
                }
            }

            /// Replace the whole constraint set of the backing entity in place (wipe-and-set)
            /// from a value container or a live view.
            #[setter]
            fn set_constraints(&self, py: Python<'_>, value: $arg) -> PyResult<()> {
                self.owner
                    .borrow_mut(py)
                    .inner_mut()
                    .$entity_mut(self.id)
                    .ast
                    .constraints = value.to_rust(py)?;
                Ok(())
            }

            /// The value fields as a dict: `configuration` plus a `constraints` list of the
            /// entries — symmetric with the value pyclass's `asdict`, read through the view.
            fn asdict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
                let molecule = self.owner.bind(py).borrow();
                let ast = self.view(molecule.inner())?.ast;
                let dict = PyDict::new(py);
                dict.set_item(
                    "configuration",
                    StereoConfigurationAst::from_rust(py, &ast.configuration)?,
                )?;
                let constraints = ast
                    .constraints
                    .iter()
                    .map(|c| into_py_variant(py, $constraint::from_rust(py, c)?))
                    .collect::<PyResult<Vec<_>>>()?;
                dict.set_item("constraints", constraints)?;
                Ok(dict)
            }
        }
    };
}

stereo_view! {
    StereoAtomView, AstStereoAtomView, AstStereoAtomId, stereo_atoms, stereo_atom_mut,
    "stereo atom id out of range", StereoAtomConstraintAst, StereoAtomConstraintsView,
    StereoAtomConstraintsBacking, StereoAtomConstraintsArg,
}

stereo_view! {
    StereoBondView, AstStereoBondView, AstStereoBondId, stereo_bonds, stereo_bond_mut,
    "stereo bond id out of range", StereoBondConstraintAst, StereoBondConstraintsView,
    StereoBondConstraintsBacking, StereoBondConstraintsArg,
}

/// Per-entity molecule-level stereo collection — `StereoAtomViews` / `StereoBondViews` — the
/// stereo atoms/bonds of a molecule, indexed by integer position, plus content lookups by
/// site (`at`) and by site + ligand multiset (`of`).
macro_rules! stereo_views {
    (
        $views:ident, $view:ident, $iter:ident, $ast_id:ident, $site_id:ident, $namespace:ident,
        $entity_mut:ident, $value:ident, $resolve_index:ident, $id_error:literal $(,)?
    ) => {
        /// Resolve a possibly-negative Python index (negative counts from the end) into an
        /// existing stereo entity id, or `IndexError`.
        fn $resolve_index(molecule: &AstMoleculeAst, index: isize) -> PyResult<$ast_id> {
            let count = molecule.$namespace().count();
            let resolved = if index < 0 {
                index + count as isize
            } else {
                index
            };
            if resolved < 0 {
                return Err(PyIndexError::new_err($id_error));
            }
            let id = $ast_id(resolved as u32);
            if molecule.$namespace().contains(id) {
                Ok(id)
            } else {
                Err(PyIndexError::new_err($id_error))
            }
        }

        #[pyclass]
        pub struct $views {
            owner: Py<MoleculeAst>,
        }

        #[pymethods]
        impl $views {
            fn __len__(&self, py: Python<'_>) -> usize {
                self.owner.bind(py).borrow().inner().$namespace().count()
            }

            fn __repr__(&self, py: Python<'_>) -> String {
                format!(
                    "{}(len={})",
                    stringify!($views),
                    self.owner.bind(py).borrow().inner().$namespace().count()
                )
            }

            fn __getitem__(&self, py: Python<'_>, index: isize) -> PyResult<$view> {
                let molecule = self.owner.bind(py).borrow();
                let id = $resolve_index(molecule.inner(), index)?;
                Ok($view {
                    owner: self.owner.clone_ref(py),
                    id,
                })
            }

            /// Replace the whole stereo entity value at `index` in place (site and ligands
            /// unchanged).
            fn __setitem__(
                &self,
                py: Python<'_>,
                index: isize,
                value: PyRef<'_, $value>,
            ) -> PyResult<()> {
                let mut molecule = self.owner.borrow_mut(py);
                let id = $resolve_index(molecule.inner(), index)?;
                *molecule.inner_mut().$entity_mut(id).ast = value.inner().clone();
                Ok(())
            }

            /// The stereo entity sitting on the atom/bond with id `site`, or `None`. Keyed by
            /// site id, *not* by position — use `views[i]` to index by position.
            fn at(&self, py: Python<'_>, site: u32) -> Option<$view> {
                let molecule = self.owner.bind(py).borrow();
                molecule
                    .inner()
                    .$namespace()
                    .at_id($site_id(site))
                    .map(|id| $view {
                        owner: self.owner.clone_ref(py),
                        id,
                    })
            }

            /// The stereo entity on `site` with exactly `ligands` (order-independent), or `None`.
            fn of(&self, py: Python<'_>, site: u32, ligands: Vec<StereoLigand>) -> Option<$view> {
                let ligands: Vec<AstStereoLigand> =
                    ligands.into_iter().map(StereoLigand::to_rust).collect();
                let molecule = self.owner.bind(py).borrow();
                molecule
                    .inner()
                    .$namespace()
                    .of_id($site_id(site), &ligands)
                    .map(|id| $view {
                        owner: self.owner.clone_ref(py),
                        id,
                    })
            }

            fn __iter__(&self, py: Python<'_>) -> $iter {
                let ids = self
                    .owner
                    .bind(py)
                    .borrow()
                    .inner()
                    .$namespace()
                    .ids()
                    .collect::<Vec<_>>();
                $iter {
                    owner: self.owner.clone_ref(py),
                    ids: ids.into_iter(),
                }
            }
        }

        impl $views {
            /// Build the stereo-views handle for `owner` (the `mol.stereo_{atoms,bonds}` accessor).
            pub(crate) fn new(owner: Py<MoleculeAst>) -> $views {
                $views { owner }
            }
        }

        #[pyclass]
        struct $iter {
            owner: Py<MoleculeAst>,
            ids: IntoIter<$ast_id>,
        }

        #[pymethods]
        impl $iter {
            fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
                slf
            }

            fn __next__(&mut self, py: Python<'_>) -> Option<$view> {
                self.ids.next().map(|id| $view {
                    owner: self.owner.clone_ref(py),
                    id,
                })
            }
        }
    };
}

stereo_views! {
    StereoAtomViews, StereoAtomView, StereoAtomViewIter, AstStereoAtomId, AstAtomId, stereo_atoms,
    stereo_atom_mut, StereoAtomAst, resolve_stereo_atom_index, "stereo atom id out of range",
}

stereo_views! {
    StereoBondViews, StereoBondView, StereoBondViewIter, AstStereoBondId, AstBondId, stereo_bonds,
    stereo_bond_mut, StereoBondAst, resolve_stereo_bond_index, "stereo bond id out of range",
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::{
        AtomAst as AstAtomAst, BondAst as AstBondAst, FluxionalityAst as AstFluxionalityAst,
        LigandSymmetryAst as AstLigandSymmetryAst, MoleculeParts,
        StereoAtomConstraintAst as AstStereoAtomConstraintAst,
        StereoAtomConstraintKey as AstStereoAtomConstraintKey,
        StereoAtomConstraintsAst as AstStereoAtomConstraintsAst,
        StereoBondConstraintAst as AstStereoBondConstraintAst,
        StereoBondConstraintsAst as AstStereoBondConstraintsAst,
        StereoLigandPair as AstStereoLigandPair, StereoLigandPosition as AstStereoLigandPosition,
        StereogenicityAst as AstStereogenicityAst, TopicityAst as AstTopicityAst,
        TopicityRelationAst as AstTopicityRelationAst,
    };
    use umol_chem::element::Element as ChemElement;

    use super::*;
    use crate::boolean::{BooleanArg, BooleanAst};

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
            assert_eq!(StereoTerm::from_rust(py, &ast).unwrap().to_rust(py), ast);
        });
    }

    #[rstest]
    #[case(AstStereoCosetAst::Undetermined)]
    #[case(AstStereoCosetAst::Lit(1))]
    #[case(AstStereoCosetAst::LitSet(BTreeSet::from([0, 1])))]
    #[case(AstStereoCosetAst::Term(Box::new(AstStereoTerm::Lit(1))))]
    fn test_stereo_coset_ast_roundtrip(#[case] ast: AstStereoCosetAst) {
        Python::attach(|py| {
            assert_eq!(
                StereoCosetAst::from_rust(py, &ast).unwrap().to_rust(py),
                ast
            );
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
                TetrahedralStereoAst::from_rust(py, &ast)
                    .unwrap()
                    .to_rust(py),
                ast
            );
        });
    }

    #[rstest]
    #[case(TetrahedralStereo::Ccw, AstStereoCosetAst::Lit(0))]
    #[case(TetrahedralStereo::Cw, AstStereoCosetAst::Lit(1))]
    fn test_tetrahedral_stereo_to_rust(
        #[case] config: TetrahedralStereo,
        #[case] coset: AstStereoCosetAst,
    ) {
        assert_eq!(config.to_rust(), AstTetrahedralStereoAst::Stereo(coset));
    }

    #[rstest]
    #[case(AstCisTransStereoAst::Undetermined)]
    #[case(AstCisTransStereoAst::NotStereo)]
    #[case(AstCisTransStereoAst::Stereo(AstStereoCosetAst::Lit(1)))]
    #[case(AstCisTransStereoAst::Stereo(AstStereoCosetAst::Term(Box::new(AstStereoTerm::Lit(0)))))]
    fn test_cis_trans_stereo_ast_roundtrip(#[case] ast: AstCisTransStereoAst) {
        Python::attach(|py| {
            assert_eq!(
                CisTransStereoAst::from_rust(py, &ast).unwrap().to_rust(py),
                ast
            );
        });
    }

    #[rstest]
    #[case(CisTransStereo::Z, AstStereoCosetAst::Lit(0))]
    #[case(CisTransStereo::E, AstStereoCosetAst::Lit(1))]
    fn test_cis_trans_stereo_to_rust(
        #[case] config: CisTransStereo,
        #[case] coset: AstStereoCosetAst,
    ) {
        assert_eq!(config.to_rust(), AstCisTransStereoAst::Stereo(coset));
    }

    #[rstest]
    #[case(AstStereoKind::Tetrahedral)]
    #[case(AstStereoKind::CisTrans)]
    #[case(AstStereoKind::Axial)]
    #[case(AstStereoKind::SquarePlanar)]
    #[case(AstStereoKind::TrigonalBipyramidal)]
    #[case(AstStereoKind::Octahedral)]
    fn test_stereo_kind_roundtrip(#[case] ast: AstStereoKind) {
        assert_eq!(StereoKind::from_rust(ast).to_rust(), ast);
    }

    #[rstest]
    #[case(AstStereoLigandKind::Atom)]
    #[case(AstStereoLigandKind::ImplicitHydrogen)]
    #[case(AstStereoLigandKind::LonePair)]
    fn test_stereo_ligand_kind_roundtrip(#[case] ast: AstStereoLigandKind) {
        assert_eq!(StereoLigandKind::from_rust(ast).to_rust(), ast);
    }

    #[rstest]
    #[case(AstTopicity::Homotopic)]
    #[case(AstTopicity::Enantiotopic)]
    #[case(AstTopicity::Diastereotopic)]
    fn test_topicity_roundtrip(#[case] ast: AstTopicity) {
        assert_eq!(Topicity::from_rust(ast).to_rust(), ast);
    }

    #[rstest]
    #[case(AstStereogenicity::Symmetric)]
    #[case(AstStereogenicity::Prochiral)]
    #[case(AstStereogenicity::Stereogenic)]
    fn test_stereogenicity_roundtrip(#[case] ast: AstStereogenicity) {
        assert_eq!(Stereogenicity::from_rust(ast).to_rust(), ast);
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
        assert_eq!(StereoLigand::from_rust(ast).to_rust(), ast);
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
                    StereoConfigurationAst::from_rust(py, &ast)
                        .unwrap()
                        .to_rust(py),
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
                config.coset(py).unwrap().bind(py).borrow().to_rust(py),
                AstStereoCosetAst::Lit(1)
            );
            let undetermined = StereoConfigurationAst::Undetermined();
            assert_eq!(undetermined.kind(), None);
            assert!(undetermined.coset(py).is_none());
        });
    }

    #[rstest]
    fn test_stereo_configuration_arg_to_rust() {
        Python::attach(|py| {
            // the Th shorthand → Kinded(Tetrahedral, coset)
            assert_eq!(
                StereoConfigurationArg::Tetrahedral(TetrahedralStereo::Cw).to_rust(py),
                AstStereoConfigurationAst::kinded(
                    AstStereoKind::Tetrahedral,
                    AstStereoCosetAst::Lit(1)
                )
            );
            // the Ct shorthand → Kinded(CisTrans, coset)
            assert_eq!(
                StereoConfigurationArg::CisTrans(CisTransStereo::E).to_rust(py),
                AstStereoConfigurationAst::kinded(
                    AstStereoKind::CisTrans,
                    AstStereoCosetAst::Lit(1)
                )
            );
            // a StereoConfigurationAst passes through
            let config = Py::new(py, StereoConfigurationAst::Undetermined()).unwrap();
            assert_eq!(
                StereoConfigurationArg::Ast(config).to_rust(py),
                AstStereoConfigurationAst::Undetermined
            );
        });
    }

    #[rstest]
    #[case(PermOrientation::Proper)]
    #[case(PermOrientation::Improper)]
    fn test_orientation_roundtrip(#[case] ast: PermOrientation) {
        assert_eq!(Orientation::from_rust(ast).to_rust(), ast);
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
        assert_eq!(LigandPermutation::from_rust(ast).to_rust(), ast);
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
        assert_eq!(OrientedLigandPermutation::from_rust(ast).to_rust(), ast);
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
        assert_eq!(StereoLigandPair::from_rust(ast).to_rust(), ast);
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
        assert_eq!(TopicityRelationAst::from_rust(&ast).to_rust(), ast);
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
        assert_eq!(StereogenicityAst::from_rust(&ast).to_rust(), ast);
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
                value.invariant.bind(py).borrow().to_rust(),
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
                    LigandSymmetryAst::from_rust(py, &ast).unwrap().to_rust(py),
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
                value.active.bind(py).borrow().to_rust(),
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
                assert_eq!(
                    FluxionalityAst::from_rust(py, &ast).unwrap().to_rust(py),
                    ast
                );
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
                value.relation.bind(py).borrow().to_rust(),
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
                assert_eq!(TopicityAst::from_rust(py, &ast).unwrap().to_rust(py), ast);
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
                StereoAtomConstraintAst::from_rust(py, &ast).unwrap().to_rust(py),
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
            let key = StereoAtomConstraintAst::from_rust(py, &ast)
                .unwrap()
                .key(py)
                .unwrap();
            assert_eq!(
                key.to_rust(py),
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
                    .to_rust(py),
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
                StereoAtomConstraintAst::from_rust(
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
                popped.unwrap().to_rust(py),
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
                constraints.stereogenicity().to_rust(),
                AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic)
            );
            assert_eq!(
                constraints.topicity(StereoLigandPair::new(0, 1)).to_rust(),
                AstTopicityRelationAst::Lit(AstTopicity::Homotopic)
            );
            let ligand_symmetries = constraints.ligand_symmetries(py).unwrap();
            assert_eq!(ligand_symmetries.len(), 1);
            assert_eq!(
                ligand_symmetries[0].to_rust(py).invariant,
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
                .map(|k| k.bind(py).borrow().to_rust(py))
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
                .map(|c| c.bind(py).borrow().to_rust(py))
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
                StereoAtomConstraintAst::from_rust(
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
    fn test_stereo_atom_constraints_arg_to_rust() {
        Python::attach(|py| {
            let entry = into_py_variant(
                py,
                StereoAtomConstraintAst::from_rust(
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
            assert_eq!(arg.to_rust(py).unwrap(), expected);
        });
    }

    // `StereoBondConstraintsAst` is the second `stereo_constraints!` instantiation; the shared
    // macro is covered by the `StereoAtom` tests above. This confirms the bond instantiation
    // and exercises its `from_inner` / `Arg::to_rust`.
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
                constraints.stereogenicity().to_rust(),
                AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic)
            );

            let mut container_ast = AstStereoBondConstraintsAst::new();
            container_ast.extend([stereogenicity.clone()]);
            let container =
                Py::new(py, StereoBondConstraintsAst::from_inner(container_ast)).unwrap();
            let arg = StereoBondConstraintsArg::Container(container);
            let mut expected = AstStereoBondConstraintsAst::new();
            expected.extend([stereogenicity]);
            assert_eq!(arg.to_rust(py).unwrap(), expected);
        });
    }

    #[rstest]
    fn test_stereo_atom_constraints_view_set() {
        Python::attach(|py| {
            let value = Py::new(
                py,
                StereoAtomAst::from_inner(AstStereoAtomAst::new(
                    AstStereoKind::Tetrahedral,
                    AstStereoCosetAst::Lit(0),
                )),
            )
            .unwrap();
            let view = StereoAtomConstraintsView {
                backing: StereoAtomConstraintsBacking::Value(value.clone_ref(py)),
            };
            let stereogenicity = into_py_variant(
                py,
                StereoAtomConstraintAst::from_rust(
                    py,
                    &AstStereoAtomConstraintAst::Stereogenicity(AstStereogenicityAst::Lit(
                        AstStereogenicity::Stereogenic,
                    )),
                )
                .unwrap(),
            )
            .unwrap();
            view.set(py, stereogenicity);
            assert_eq!(
                value.borrow(py).inner().constraints.stereogenicity(),
                AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic)
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_constraints_view_pop() {
        Python::attach(|py| {
            let mut ast_cs = AstStereoAtomConstraintsAst::new();
            ast_cs.extend([AstStereoAtomConstraintAst::Stereogenicity(
                AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic),
            )]);
            let value = Py::new(
                py,
                StereoAtomAst::from_inner(AstStereoAtomAst {
                    configuration: AstStereoConfigurationAst::Kinded(
                        AstStereoKind::Tetrahedral,
                        AstStereoCosetAst::Lit(0),
                    ),
                    constraints: ast_cs,
                }),
            )
            .unwrap();
            let view = StereoAtomConstraintsView {
                backing: StereoAtomConstraintsBacking::Value(value.clone_ref(py)),
            };
            let key = into_py_variant(py, StereoAtomConstraintKey::Stereogenicity()).unwrap();
            let popped = view.pop(py, key).unwrap();
            assert_eq!(
                popped.unwrap().to_rust(py),
                AstStereoAtomConstraintAst::Stereogenicity(AstStereogenicityAst::Lit(
                    AstStereogenicity::Stereogenic
                ))
            );
            assert_eq!(value.borrow(py).inner().constraints.len(), 0);
        });
    }

    #[rstest]
    fn test_stereo_atom_constraints_view_getitem() {
        Python::attach(|py| {
            let mut ast_cs = AstStereoAtomConstraintsAst::new();
            ast_cs.extend([AstStereoAtomConstraintAst::Stereogenicity(
                AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic),
            )]);
            let value = Py::new(
                py,
                StereoAtomAst::from_inner(AstStereoAtomAst {
                    configuration: AstStereoConfigurationAst::Kinded(
                        AstStereoKind::Tetrahedral,
                        AstStereoCosetAst::Lit(0),
                    ),
                    constraints: ast_cs,
                }),
            )
            .unwrap();
            let view = StereoAtomConstraintsView {
                backing: StereoAtomConstraintsBacking::Value(value.clone_ref(py)),
            };
            assert_eq!(view.__len__(py).unwrap(), 1);
            let present = into_py_variant(py, StereoAtomConstraintKey::Stereogenicity()).unwrap();
            assert!(view.__contains__(py, present.clone_ref(py)).unwrap());
            assert_eq!(
                view.__getitem__(py, present).unwrap().to_rust(py),
                AstStereoAtomConstraintAst::Stereogenicity(AstStereogenicityAst::Lit(
                    AstStereogenicity::Stereogenic
                ))
            );
            let absent = into_py_variant(
                py,
                StereoAtomConstraintKey::Topicity(
                    into_py_variant(py, StereoLigandPair::new(0, 1)).unwrap(),
                ),
            )
            .unwrap();
            assert!(!view.__contains__(py, absent.clone_ref(py)).unwrap());
            assert!(view.__getitem__(py, absent).is_err());
        });
    }

    #[rstest]
    fn test_stereo_atom_constraints_view_items() {
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
            let value = Py::new(
                py,
                StereoAtomAst::from_inner(AstStereoAtomAst {
                    configuration: AstStereoConfigurationAst::Kinded(
                        AstStereoKind::Tetrahedral,
                        AstStereoCosetAst::Lit(0),
                    ),
                    constraints: ast_cs,
                }),
            )
            .unwrap();
            let view = StereoAtomConstraintsView {
                backing: StereoAtomConstraintsBacking::Value(value.clone_ref(py)),
            };
            let keys: Vec<AstStereoAtomConstraintKey> = view
                .keys(py)
                .unwrap()
                .keys
                .map(|k| k.bind(py).borrow().to_rust(py))
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
            assert_eq!(view.values(py).unwrap().entries.count(), 2);
            assert_eq!(view.items(py).unwrap().items.count(), 2);
        });
    }

    #[rstest]
    fn test_stereo_atom_constraints_view_update() {
        Python::attach(|py| {
            let value = Py::new(
                py,
                StereoAtomAst::from_inner(AstStereoAtomAst::new(
                    AstStereoKind::Tetrahedral,
                    AstStereoCosetAst::Lit(0),
                )),
            )
            .unwrap();
            let view = StereoAtomConstraintsView {
                backing: StereoAtomConstraintsBacking::Value(value.clone_ref(py)),
            };
            let entry = into_py_variant(
                py,
                StereoAtomConstraintAst::from_rust(
                    py,
                    &AstStereoAtomConstraintAst::Stereogenicity(AstStereogenicityAst::Lit(
                        AstStereogenicity::Stereogenic,
                    )),
                )
                .unwrap(),
            )
            .unwrap();
            view.update(py, StereoAtomConstraintsUpdate::Entries(vec![entry]))
                .unwrap();
            assert_eq!(
                value.borrow(py).inner().constraints.stereogenicity(),
                AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic)
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_constraints_view_accessors() {
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
            let value = Py::new(
                py,
                StereoAtomAst::from_inner(AstStereoAtomAst {
                    configuration: AstStereoConfigurationAst::Kinded(
                        AstStereoKind::Tetrahedral,
                        AstStereoCosetAst::Lit(0),
                    ),
                    constraints: ast_cs,
                }),
            )
            .unwrap();
            let view = StereoAtomConstraintsView {
                backing: StereoAtomConstraintsBacking::Value(value.clone_ref(py)),
            };
            assert_eq!(
                view.stereogenicity(py).unwrap().to_rust(),
                AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic)
            );
            assert_eq!(
                view.topicity(py, StereoLigandPair::new(0, 1))
                    .unwrap()
                    .to_rust(),
                AstTopicityRelationAst::Lit(AstTopicity::Homotopic)
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_ast_constraints() {
        Python::attach(|py| {
            let value = Py::new(
                py,
                StereoAtomAst::from_inner(AstStereoAtomAst::new(
                    AstStereoKind::Tetrahedral,
                    AstStereoCosetAst::Lit(0),
                )),
            )
            .unwrap();
            let view = StereoAtomAst::constraints(value.clone_ref(py));
            let stereogenicity = into_py_variant(
                py,
                StereoAtomConstraintAst::from_rust(
                    py,
                    &AstStereoAtomConstraintAst::Stereogenicity(AstStereogenicityAst::Lit(
                        AstStereogenicity::Stereogenic,
                    )),
                )
                .unwrap(),
            )
            .unwrap();
            view.set(py, stereogenicity);
            assert_eq!(
                value.borrow(py).inner().constraints.stereogenicity(),
                AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic)
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_ast_set_constraints_self() {
        Python::attach(|py| {
            let mut ast_cs = AstStereoAtomConstraintsAst::new();
            ast_cs.extend([AstStereoAtomConstraintAst::Stereogenicity(
                AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic),
            )]);
            let value = Py::new(
                py,
                StereoAtomAst::from_inner(AstStereoAtomAst {
                    configuration: AstStereoConfigurationAst::Kinded(
                        AstStereoKind::Tetrahedral,
                        AstStereoCosetAst::Lit(0),
                    ),
                    constraints: ast_cs,
                }),
            )
            .unwrap();
            let own_view = StereoAtomAst::constraints(value.clone_ref(py));
            StereoAtomAst::set_constraints(
                value.clone_ref(py),
                py,
                StereoAtomConstraintsArg::View(Py::new(py, own_view).unwrap()),
            )
            .unwrap();
            assert_eq!(
                value.borrow(py).inner().constraints.stereogenicity(),
                AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic)
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_constraints_view_update_self() {
        Python::attach(|py| {
            let mut ast_cs = AstStereoAtomConstraintsAst::new();
            ast_cs.extend([AstStereoAtomConstraintAst::Stereogenicity(
                AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic),
            )]);
            let value = Py::new(
                py,
                StereoAtomAst::from_inner(AstStereoAtomAst {
                    configuration: AstStereoConfigurationAst::Kinded(
                        AstStereoKind::Tetrahedral,
                        AstStereoCosetAst::Lit(0),
                    ),
                    constraints: ast_cs,
                }),
            )
            .unwrap();
            let view = StereoAtomConstraintsView {
                backing: StereoAtomConstraintsBacking::Value(value.clone_ref(py)),
            };
            let own = StereoAtomAst::constraints(value.clone_ref(py));
            view.update(
                py,
                StereoAtomConstraintsUpdate::View(Py::new(py, own).unwrap()),
            )
            .unwrap();
            assert_eq!(value.borrow(py).inner().constraints.len(), 1);
        });
    }

    #[rstest]
    fn test_stereo_bond_constraints_view_set() {
        Python::attach(|py| {
            let value = Py::new(
                py,
                StereoBondAst::from_inner(AstStereoBondAst::new(
                    AstStereoKind::CisTrans,
                    AstStereoCosetAst::Lit(0),
                )),
            )
            .unwrap();
            let view = StereoBondConstraintsView {
                backing: StereoBondConstraintsBacking::Value(value.clone_ref(py)),
            };
            let stereogenicity = into_py_variant(
                py,
                StereoBondConstraintAst::from_rust(
                    py,
                    &AstStereoBondConstraintAst::Stereogenicity(AstStereogenicityAst::Lit(
                        AstStereogenicity::Stereogenic,
                    )),
                )
                .unwrap(),
            )
            .unwrap();
            view.set(py, stereogenicity);
            assert_eq!(
                value.borrow(py).inner().constraints.stereogenicity(),
                AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic)
            );
        });
    }

    #[rstest]
    #[case::ccw(
        StereoConfigurationArg::Tetrahedral(TetrahedralStereo::Ccw),
        AstStereoAtomAst::new(AstStereoKind::Tetrahedral, AstStereoCosetAst::Lit(0))
    )]
    #[case::cw(
        StereoConfigurationArg::Tetrahedral(TetrahedralStereo::Cw),
        AstStereoAtomAst::new(AstStereoKind::Tetrahedral, AstStereoCosetAst::Lit(1))
    )]
    fn test_stereo_atom_ast_new(
        #[case] configuration: StereoConfigurationArg,
        #[case] expected: AstStereoAtomAst,
    ) {
        Python::attach(|py| {
            let value = StereoAtomAst::new(py, configuration, None);
            assert_eq!(*value.inner(), expected);
        });
    }

    #[rstest]
    fn test_stereo_atom_ast_new_constraints() {
        Python::attach(|py| {
            let stereogenicity = AstStereoAtomConstraintAst::Stereogenicity(
                AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic),
            );
            let mut ast_cs = AstStereoAtomConstraintsAst::new();
            ast_cs.extend([stereogenicity.clone()]);
            let container = Py::new(py, StereoAtomConstraintsAst::from_inner(ast_cs)).unwrap();
            let value = StereoAtomAst::new(
                py,
                StereoConfigurationArg::Tetrahedral(TetrahedralStereo::Ccw),
                Some(container),
            );
            let mut expected_cs = AstStereoAtomConstraintsAst::new();
            expected_cs.extend([stereogenicity]);
            assert_eq!(
                *value.inner(),
                AstStereoAtomAst {
                    configuration: AstStereoConfigurationAst::Kinded(
                        AstStereoKind::Tetrahedral,
                        AstStereoCosetAst::Lit(0)
                    ),
                    constraints: expected_cs,
                }
            );
        });
    }

    #[rstest]
    #[case::ccw(
        "Th0",
        AstStereoConfigurationAst::Kinded(AstStereoKind::Tetrahedral, AstStereoCosetAst::Lit(0))
    )]
    #[case::undetermined_coset(
        "Th*",
        AstStereoConfigurationAst::Kinded(
            AstStereoKind::Tetrahedral,
            AstStereoCosetAst::Undetermined
        )
    )]
    fn test_stereo_atom_ast_parse(
        #[case] input: &str,
        #[case] expected: AstStereoConfigurationAst,
    ) {
        let value = StereoAtomAst::parse(input).unwrap();
        assert_eq!(value.inner().configuration, expected);
    }

    #[rstest]
    fn test_stereo_atom_ast_parse_error() {
        assert!(StereoAtomAst::parse("not-a-stereo-atom").is_err());
    }

    #[rstest]
    #[case::ccw(
        AstStereoAtomAst::new(AstStereoKind::Tetrahedral, AstStereoCosetAst::Lit(0)),
        "Th0"
    )]
    #[case::square_planar(
        AstStereoAtomAst::new(AstStereoKind::SquarePlanar, AstStereoCosetAst::Lit(2)),
        "Sp2"
    )]
    fn test_stereo_atom_ast_str(#[case] ast: AstStereoAtomAst, #[case] expected: &str) {
        let value = StereoAtomAst::from_inner(ast);
        assert_eq!(value.__str__(), expected);
    }

    #[rstest]
    fn test_stereo_atom_ast_repr() {
        let value = StereoAtomAst::from_inner(AstStereoAtomAst::new(
            AstStereoKind::Tetrahedral,
            AstStereoCosetAst::Lit(0),
        ));
        assert_eq!(value.__repr__(), "StereoAtomAst.parse('Th0')");
    }

    #[rstest]
    fn test_stereo_atom_ast_configuration() {
        Python::attach(|py| {
            let value = StereoAtomAst::from_inner(AstStereoAtomAst::new(
                AstStereoKind::Tetrahedral,
                AstStereoCosetAst::Lit(0),
            ));
            assert_eq!(
                value.configuration(py).unwrap().to_rust(py),
                AstStereoConfigurationAst::Kinded(
                    AstStereoKind::Tetrahedral,
                    AstStereoCosetAst::Lit(0)
                )
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_ast_set_configuration() {
        Python::attach(|py| {
            let mut value = StereoAtomAst::from_inner(AstStereoAtomAst::new(
                AstStereoKind::Tetrahedral,
                AstStereoCosetAst::Lit(0),
            ));
            value.set_configuration(
                py,
                StereoConfigurationArg::Tetrahedral(TetrahedralStereo::Cw),
            );
            assert_eq!(
                value.inner().configuration,
                AstStereoConfigurationAst::Kinded(
                    AstStereoKind::Tetrahedral,
                    AstStereoCosetAst::Lit(1)
                )
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_ast_set_constraints() {
        Python::attach(|py| {
            let value = Py::new(
                py,
                StereoAtomAst::from_inner(AstStereoAtomAst::new(
                    AstStereoKind::Tetrahedral,
                    AstStereoCosetAst::Lit(0),
                )),
            )
            .unwrap();
            let stereogenicity = AstStereoAtomConstraintAst::Stereogenicity(
                AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic),
            );
            let mut ast_cs = AstStereoAtomConstraintsAst::new();
            ast_cs.extend([stereogenicity.clone()]);
            let container = Py::new(py, StereoAtomConstraintsAst::from_inner(ast_cs)).unwrap();
            StereoAtomAst::set_constraints(
                value.clone_ref(py),
                py,
                StereoAtomConstraintsArg::Container(container),
            )
            .unwrap();
            let mut expected_cs = AstStereoAtomConstraintsAst::new();
            expected_cs.extend([stereogenicity]);
            assert_eq!(value.borrow(py).inner().constraints, expected_cs);
        });
    }

    #[rstest]
    fn test_stereo_atom_ast_asdict() {
        Python::attach(|py| {
            let stereogenicity = AstStereoAtomConstraintAst::Stereogenicity(
                AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic),
            );
            let mut ast_cs = AstStereoAtomConstraintsAst::new();
            ast_cs.extend([stereogenicity]);
            let value = StereoAtomAst::from_inner(AstStereoAtomAst {
                configuration: AstStereoConfigurationAst::Kinded(
                    AstStereoKind::Tetrahedral,
                    AstStereoCosetAst::Lit(0),
                ),
                constraints: ast_cs,
            });
            let dict = value.asdict(py).unwrap();
            let configuration = dict.get_item("configuration").unwrap().unwrap();
            let expected = into_py_variant(
                py,
                StereoConfigurationAst::from_rust(
                    py,
                    &AstStereoConfigurationAst::Kinded(
                        AstStereoKind::Tetrahedral,
                        AstStereoCosetAst::Lit(0),
                    ),
                )
                .unwrap(),
            )
            .unwrap();
            assert!(configuration.eq(expected.bind(py)).unwrap());
            let constraints = dict.get_item("constraints").unwrap().unwrap();
            assert_eq!(constraints.len().unwrap(), 1);
        });
    }

    #[rstest]
    fn test_stereo_bond_ast_new() {
        Python::attach(|py| {
            let value = StereoBondAst::new(
                py,
                StereoConfigurationArg::CisTrans(CisTransStereo::Z),
                None,
            );
            assert_eq!(
                *value.inner(),
                AstStereoBondAst::new(AstStereoKind::CisTrans, AstStereoCosetAst::Lit(0))
            );
            assert_eq!(value.__str__(), "Ct0");
        });
    }

    #[rstest]
    #[case::z(
        "Ct0",
        AstStereoBondAst::new(AstStereoKind::CisTrans, AstStereoCosetAst::Lit(0))
    )]
    #[case::e(
        "Ct1",
        AstStereoBondAst::new(AstStereoKind::CisTrans, AstStereoCosetAst::Lit(1))
    )]
    fn test_stereo_bond_ast_parse(#[case] input: &str, #[case] expected: AstStereoBondAst) {
        let value = StereoBondAst::parse(input).unwrap();
        assert_eq!(*value.inner(), expected);
    }

    #[rstest]
    fn test_stereo_bond_ast_str() {
        let value = StereoBondAst::from_inner(AstStereoBondAst::new(
            AstStereoKind::CisTrans,
            AstStereoCosetAst::Lit(1),
        ));
        assert_eq!(value.__str__(), "Ct1");
    }

    fn stereo_atom_molecule(py: Python<'_>) -> Py<MoleculeAst> {
        let molecule = AstMoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AstAtomAst::from_element(ChemElement::C); 5],
            stereo_atoms: vec![(
                AstAtomId(0),
                vec![
                    AstStereoLigand::new(AstAtomId(1), AstStereoLigandKind::Atom),
                    AstStereoLigand::new(AstAtomId(2), AstStereoLigandKind::Atom),
                    AstStereoLigand::new(AstAtomId(3), AstStereoLigandKind::Atom),
                    AstStereoLigand::new(AstAtomId(4), AstStereoLigandKind::Atom),
                ],
                AstStereoAtomAst::new(AstStereoKind::Tetrahedral, AstStereoCosetAst::Lit(0)),
            )],
            ..Default::default()
        });
        Py::new(py, MoleculeAst::from_inner(molecule)).unwrap()
    }

    fn stereo_bond_molecule(py: Python<'_>) -> Py<MoleculeAst> {
        let molecule = AstMoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AstAtomAst::from_element(ChemElement::C); 4],
            bonds: vec![
                (AstAtomId(0), AstAtomId(1), AstBondAst::from_order(2)),
                (AstAtomId(0), AstAtomId(2), AstBondAst::from_order(1)),
                (AstAtomId(1), AstAtomId(3), AstBondAst::from_order(1)),
            ],
            stereo_bonds: vec![(
                AstBondId(0),
                vec![
                    AstStereoLigand::new(AstAtomId(2), AstStereoLigandKind::Atom),
                    AstStereoLigand::new(AstAtomId(3), AstStereoLigandKind::Atom),
                ],
                AstStereoBondAst::new(AstStereoKind::CisTrans, AstStereoCosetAst::Lit(0)),
            )],
            ..Default::default()
        });
        Py::new(py, MoleculeAst::from_inner(molecule)).unwrap()
    }

    #[rstest]
    fn test_stereo_atom_view_id() {
        Python::attach(|py| {
            let view = StereoAtomView {
                owner: stereo_atom_molecule(py),
                id: AstStereoAtomId(0),
            };
            assert_eq!(view.id(), 0);
            assert_eq!(view.__repr__(), "StereoAtomView(id=0)");
        });
    }

    #[rstest]
    fn test_stereo_atom_view_site_id() {
        Python::attach(|py| {
            let view = StereoAtomView {
                owner: stereo_atom_molecule(py),
                id: AstStereoAtomId(0),
            };
            assert_eq!(view.site_id(py).unwrap(), 0);
        });
    }

    #[rstest]
    fn test_stereo_atom_view_ligands() {
        Python::attach(|py| {
            let view = StereoAtomView {
                owner: stereo_atom_molecule(py),
                id: AstStereoAtomId(0),
            };
            assert_eq!(
                view.ligands(py)
                    .unwrap()
                    .iter()
                    .map(|l| (l.atom_id, l.kind))
                    .collect::<Vec<_>>(),
                vec![
                    (1, StereoLigandKind::Atom),
                    (2, StereoLigandKind::Atom),
                    (3, StereoLigandKind::Atom),
                    (4, StereoLigandKind::Atom),
                ]
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_view_kind() {
        Python::attach(|py| {
            let view = StereoAtomView {
                owner: stereo_atom_molecule(py),
                id: AstStereoAtomId(0),
            };
            assert_eq!(view.kind(py).unwrap(), StereoKind::Tetrahedral);
        });
    }

    #[rstest]
    fn test_stereo_atom_view_coset() {
        Python::attach(|py| {
            let view = StereoAtomView {
                owner: stereo_atom_molecule(py),
                id: AstStereoAtomId(0),
            };
            assert_eq!(
                view.coset(py).unwrap().to_rust(py),
                AstStereoCosetAst::Lit(0)
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_view_configuration() {
        Python::attach(|py| {
            let view = StereoAtomView {
                owner: stereo_atom_molecule(py),
                id: AstStereoAtomId(0),
            };
            assert_eq!(
                view.configuration(py).unwrap().to_rust(py),
                AstStereoConfigurationAst::Kinded(
                    AstStereoKind::Tetrahedral,
                    AstStereoCosetAst::Lit(0)
                )
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_view_set_configuration() {
        Python::attach(|py| {
            let view = StereoAtomView {
                owner: stereo_atom_molecule(py),
                id: AstStereoAtomId(0),
            };
            view.set_configuration(
                py,
                StereoConfigurationArg::Tetrahedral(TetrahedralStereo::Cw),
            );
            assert_eq!(
                view.configuration(py).unwrap().to_rust(py),
                AstStereoConfigurationAst::Kinded(
                    AstStereoKind::Tetrahedral,
                    AstStereoCosetAst::Lit(1)
                )
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_view_constraints() {
        Python::attach(|py| {
            let view = StereoAtomView {
                owner: stereo_atom_molecule(py),
                id: AstStereoAtomId(0),
            };
            let stereogenicity = into_py_variant(
                py,
                StereoAtomConstraintAst::from_rust(
                    py,
                    &AstStereoAtomConstraintAst::Stereogenicity(AstStereogenicityAst::Lit(
                        AstStereogenicity::Stereogenic,
                    )),
                )
                .unwrap(),
            )
            .unwrap();
            view.constraints(py).set(py, stereogenicity);
            // a fresh molecule-backed handle proves the write hit the molecule
            assert_eq!(
                view.constraints(py).stereogenicity(py).unwrap().to_rust(),
                AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic)
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_view_set_constraints() {
        Python::attach(|py| {
            let view = StereoAtomView {
                owner: stereo_atom_molecule(py),
                id: AstStereoAtomId(0),
            };
            let mut ast_cs = AstStereoAtomConstraintsAst::new();
            ast_cs.extend([AstStereoAtomConstraintAst::Stereogenicity(
                AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic),
            )]);
            let container = Py::new(py, StereoAtomConstraintsAst::from_inner(ast_cs)).unwrap();
            view.set_constraints(py, StereoAtomConstraintsArg::Container(container))
                .unwrap();
            assert_eq!(
                view.constraints(py).stereogenicity(py).unwrap().to_rust(),
                AstStereogenicityAst::Lit(AstStereogenicity::Stereogenic)
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_view_asdict() {
        Python::attach(|py| {
            let view = StereoAtomView {
                owner: stereo_atom_molecule(py),
                id: AstStereoAtomId(0),
            };
            let dict = view.asdict(py).unwrap();
            let configuration = dict.get_item("configuration").unwrap().unwrap();
            let expected = into_py_variant(
                py,
                StereoConfigurationAst::from_rust(
                    py,
                    &AstStereoConfigurationAst::Kinded(
                        AstStereoKind::Tetrahedral,
                        AstStereoCosetAst::Lit(0),
                    ),
                )
                .unwrap(),
            )
            .unwrap();
            assert!(configuration.eq(expected.bind(py)).unwrap());
            assert_eq!(
                dict.get_item("constraints")
                    .unwrap()
                    .unwrap()
                    .len()
                    .unwrap(),
                0
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_view_id_out_of_range() {
        Python::attach(|py| {
            let view = StereoAtomView {
                owner: stereo_atom_molecule(py),
                id: AstStereoAtomId(5),
            };
            assert!(view.site_id(py).is_err());
        });
    }

    #[rstest]
    fn test_stereo_bond_view() {
        Python::attach(|py| {
            let view = StereoBondView {
                owner: stereo_bond_molecule(py),
                id: AstStereoBondId(0),
            };
            assert_eq!(view.id(), 0);
            assert_eq!(view.site_id(py).unwrap(), 0);
            assert_eq!(
                view.configuration(py).unwrap().to_rust(py),
                AstStereoConfigurationAst::Kinded(
                    AstStereoKind::CisTrans,
                    AstStereoCosetAst::Lit(0)
                )
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_views_len() {
        Python::attach(|py| {
            let views = StereoAtomViews {
                owner: stereo_atom_molecule(py),
            };
            assert_eq!(views.__len__(py), 1);
        });
    }

    #[rstest]
    #[case::first(0, 0)]
    #[case::negative(-1, 0)]
    fn test_stereo_atom_views_getitem(#[case] index: isize, #[case] expected_id: u32) {
        Python::attach(|py| {
            let views = StereoAtomViews {
                owner: stereo_atom_molecule(py),
            };
            assert_eq!(views.__getitem__(py, index).unwrap().id(), expected_id);
        });
    }

    #[rstest]
    #[case::past_end(1)]
    #[case::far_past_end(5)]
    fn test_stereo_atom_views_getitem_error(#[case] index: isize) {
        Python::attach(|py| {
            let views = StereoAtomViews {
                owner: stereo_atom_molecule(py),
            };
            assert!(views.__getitem__(py, index).is_err());
        });
    }

    #[rstest]
    fn test_stereo_atom_views_setitem() {
        Python::attach(|py| {
            let views = StereoAtomViews {
                owner: stereo_atom_molecule(py),
            };
            let replacement = Py::new(
                py,
                StereoAtomAst::from_inner(AstStereoAtomAst::new(
                    AstStereoKind::Tetrahedral,
                    AstStereoCosetAst::Lit(1),
                )),
            )
            .unwrap();
            views.__setitem__(py, 0, replacement.borrow(py)).unwrap();
            let view = views.__getitem__(py, 0).unwrap();
            // value replaced
            assert_eq!(
                view.configuration(py).unwrap().to_rust(py),
                AstStereoConfigurationAst::Kinded(
                    AstStereoKind::Tetrahedral,
                    AstStereoCosetAst::Lit(1)
                )
            );
            // site topology unchanged
            assert_eq!(view.site_id(py).unwrap(), 0);
        });
    }

    #[rstest]
    #[case::has_stereo(0, Some(0))]
    #[case::no_stereo(1, None)]
    fn test_stereo_atom_views_at(#[case] site: u32, #[case] expected_id: Option<u32>) {
        Python::attach(|py| {
            let views = StereoAtomViews {
                owner: stereo_atom_molecule(py),
            };
            assert_eq!(views.at(py, site).map(|v| v.id()), expected_id);
        });
    }

    #[rstest]
    fn test_stereo_atom_views_of() {
        Python::attach(|py| {
            let views = StereoAtomViews {
                owner: stereo_atom_molecule(py),
            };
            // order-independent full-ligand-set match
            let matched = views.of(
                py,
                0,
                vec![
                    StereoLigand::new(4, StereoLigandKind::Atom),
                    StereoLigand::new(3, StereoLigandKind::Atom),
                    StereoLigand::new(2, StereoLigandKind::Atom),
                    StereoLigand::new(1, StereoLigandKind::Atom),
                ],
            );
            assert_eq!(matched.map(|v| v.id()), Some(0));
            // a partial ligand set does not match
            let missed = views.of(
                py,
                0,
                vec![
                    StereoLigand::new(1, StereoLigandKind::Atom),
                    StereoLigand::new(2, StereoLigandKind::Atom),
                ],
            );
            assert!(missed.is_none());
        });
    }

    #[rstest]
    fn test_stereo_atom_views_iter() {
        Python::attach(|py| {
            let views = StereoAtomViews {
                owner: stereo_atom_molecule(py),
            };
            let mut iter = views.__iter__(py);
            let mut ids = Vec::new();
            while let Some(view) = iter.__next__(py) {
                ids.push(view.id());
            }
            assert_eq!(ids, vec![0]);
        });
    }

    #[rstest]
    fn test_stereo_bond_views() {
        Python::attach(|py| {
            let views = StereoBondViews {
                owner: stereo_bond_molecule(py),
            };
            assert_eq!(views.__len__(py), 1);
            assert_eq!(views.__getitem__(py, 0).unwrap().id(), 0);
            assert_eq!(views.at(py, 0).map(|v| v.id()), Some(0));
            assert!(views.at(py, 2).is_none());
        });
    }
}
