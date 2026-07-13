//! Stereo sub-ASTs mirroring `umol_ast::ast::stereo` and `umol_perm` (S5a): the
//! `Permutation` value, the recursive `StereoTerm` transformation algebra, the
//! `StereoCosetAst` coset, and the `TetrahedralStereoAst` atom-stereo state.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use std::collections::BTreeSet;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
// The stereo value-leaf AST mirrors are consumed only by the `#[cfg(test)]` `from_ast`/`to_ast`
// until S1b/S2/S4a wire them into the constraint values, container, and entity view.
#[cfg(test)]
use umol_ast::ast::{
    AtomId as AstAtomId, LigandPermutation as AstLigandPermutation,
    OrientedLigandPermutation as AstOrientedLigandPermutation, StereoLigand as AstStereoLigand,
    StereoLigandKind as AstStereoLigandKind, Stereogenicity as AstStereogenicity,
};
use umol_ast::ast::{
    CisTransStereoAst as AstCisTransStereoAst, StereoConfigurationAst as AstStereoConfigurationAst,
    StereoCosetAst as AstStereoCosetAst, StereoKind as AstStereoKind,
    StereoLigandPair as AstStereoLigandPair, StereoLigandPosition as AstStereoLigandPosition,
    StereoTerm as AstStereoTerm, TetrahedralStereoAst as AstTetrahedralStereoAst,
    Topicity as AstTopicity, TopicityRelationAst as AstTopicityRelationAst,
};
#[cfg(test)]
use umol_perm::Orientation as PermOrientation;
use umol_perm::Permutation as PermPermutation;

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
    #[cfg(test)]
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
/// A fieldless, hashable value enum mirroring the Rust `Stereogenicity`.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Stereogenicity {
    Symmetric,
    Prochiral,
    Stereogenic,
}

impl Stereogenicity {
    #[cfg(test)]
    pub(crate) fn from_ast(ast: AstStereogenicity) -> Self {
        match ast {
            AstStereogenicity::Symmetric => Self::Symmetric,
            AstStereogenicity::Prochiral => Self::Prochiral,
            AstStereogenicity::Stereogenic => Self::Stereogenic,
        }
    }

    #[cfg(test)]
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
    #[cfg(test)]
    pub(crate) fn from_ast(orientation: PermOrientation) -> Self {
        match orientation {
            PermOrientation::Proper => Self::Proper,
            PermOrientation::Improper => Self::Improper,
        }
    }

    #[cfg(test)]
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
    #[cfg(test)]
    pub(crate) fn from_ast(ast: AstLigandPermutation) -> Self {
        LigandPermutation {
            permutation: Permutation::from_inner(ast.0),
        }
    }

    #[cfg(test)]
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
    #[cfg(test)]
    pub(crate) fn from_ast(ast: AstOrientedLigandPermutation) -> Self {
        OrientedLigandPermutation {
            permutation: LigandPermutation::from_ast(ast.permutation),
            orientation: Orientation::from_ast(ast.orientation),
        }
    }

    #[cfg(test)]
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

    #[cfg(test)]
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
    #[cfg(test)]
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
}
