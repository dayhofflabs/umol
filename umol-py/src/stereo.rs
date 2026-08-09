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
use umol_graph_ir::ir::BooleanForm as GraphIrBooleanForm;
use umol_graph_ir::ir::{
    AsLit, AtomId as GraphIrAtomId, BondId as GraphIrBondId,
    CisTransConfiguration as GraphIrCisTransConfiguration, CisTransStereo as GraphIrCisTransStereo,
    CisTransStereoForm as GraphIrCisTransStereoForm, LigandPermutation as GraphIrLigandPermutation,
    Molecule as GraphIrMolecule, OrientedLigandPermutation as GraphIrOrientedLigandPermutation,
    StereoAtomForm as GraphIrStereoAtomForm, StereoAtomId as GraphIrStereoAtomId,
    StereoAtomUpdate as GraphIrStereoAtomUpdate, StereoAtomView as GraphIrStereoAtomView,
    StereoBondForm as GraphIrStereoBondForm, StereoBondId as GraphIrStereoBondId,
    StereoBondUpdate as GraphIrStereoBondUpdate, StereoBondView as GraphIrStereoBondView,
    StereoConfigurationForm as GraphIrStereoConfigurationForm,
    StereoConfigurationUpdate as GraphIrStereoConfigurationUpdate,
    StereoCoset as GraphIrStereoCoset, StereoKind as GraphIrStereoKind,
    StereoLigand as GraphIrStereoLigand, StereoLigandKind as GraphIrStereoLigandKind,
    StereoLigandPair as GraphIrStereoLigandPair,
    StereoLigandPosition as GraphIrStereoLigandPosition, StereoTerm as GraphIrStereoTerm,
    Stereogenicity as GraphIrStereogenicity,
    TetrahedralConfiguration as GraphIrTetrahedralConfiguration,
    TetrahedralStereo as GraphIrTetrahedralStereo,
    TetrahedralStereoForm as GraphIrTetrahedralStereoForm, Topicity as GraphIrTopicity,
};
use umol_perm::{Orientation as PermOrientation, Permutation as PermPermutation};

use crate::convert::{hash_rust, into_py_variant, variant_repr};
use crate::error::parse_error;
use crate::lattice::impl_py_lattice;
use crate::molecule::MoleculeAst;

/// A permutation of `0..degree` in one-line (image) notation.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Permutation(PermPermutation);

#[pymethods]
impl Permutation {
    /// Construct from the image (one-line notation); the degree is the image length.
    #[new]
    fn new(image: Vec<usize>) -> PyResult<Self> {
        PermPermutation::try_from(image.as_slice())
            .map(Permutation)
            .map_err(|error| PyValueError::new_err(error.to_string()))
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
    fn image(&self) -> Vec<usize> {
        (0..self.0.degree()).map(|i| self.0.apply(i)).collect()
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
    pub(crate) fn from_rust(py: Python<'_>, ast: &GraphIrStereoTerm) -> PyResult<Self> {
        Ok(match ast {
            GraphIrStereoTerm::Var(boxed) => {
                let (name, restriction) = &**boxed;
                StereoTerm::Var(name.clone(), restriction.clone())
            }
            GraphIrStereoTerm::Lit(index) => StereoTerm::Lit(*index),
            GraphIrStereoTerm::LitSet(members) => StereoTerm::LitSet(members.clone()),
            GraphIrStereoTerm::Swap(inner) => {
                StereoTerm::Swap(into_py_variant(py, StereoTerm::from_rust(py, inner)?)?)
            }
            GraphIrStereoTerm::Mirror(inner) => {
                StereoTerm::Mirror(into_py_variant(py, StereoTerm::from_rust(py, inner)?)?)
            }
            GraphIrStereoTerm::Apply(inner, permutation) => StereoTerm::Apply(
                into_py_variant(py, StereoTerm::from_rust(py, inner)?)?,
                Permutation::from_inner(*permutation),
            ),
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrStereoTerm {
        match self {
            StereoTerm::Var(name, restriction) => {
                GraphIrStereoTerm::Var(Box::new((name.clone(), restriction.clone())))
            }
            StereoTerm::Lit(index) => GraphIrStereoTerm::Lit(*index),
            StereoTerm::LitSet(members) => GraphIrStereoTerm::LitSet(members.clone()),
            StereoTerm::Swap(inner) => {
                GraphIrStereoTerm::Swap(Box::new(inner.bind(py).borrow().to_rust(py)))
            }
            StereoTerm::Mirror(inner) => {
                GraphIrStereoTerm::Mirror(Box::new(inner.bind(py).borrow().to_rust(py)))
            }
            StereoTerm::Apply(inner, permutation) => GraphIrStereoTerm::Apply(
                Box::new(inner.bind(py).borrow().to_rust(py)),
                permutation.inner(),
            ),
        }
    }
}

/// Stereo coset: undetermined, a literal coset index, a set of indices, or a
/// transformation term.
#[pyclass]
pub enum StereoCoset {
    Undetermined(),
    Lit(u32),
    LitSet(BTreeSet<u32>),
    Term(Py<StereoTerm>),
}

#[pymethods]
impl StereoCoset {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            StereoCoset::Undetermined() => ("Undetermined", 0),
            StereoCoset::Lit(_) => ("Lit", 1),
            StereoCoset::LitSet(_) => ("LitSet", 1),
            StereoCoset::Term(_) => ("Term", 1),
        };
        variant_repr(slf.bind(py).as_any(), "StereoCoset", variant, arity)
    }
}

impl StereoCoset {
    pub(crate) fn from_rust(py: Python<'_>, ast: &GraphIrStereoCoset) -> PyResult<Self> {
        Ok(match ast {
            GraphIrStereoCoset::Undetermined => Self::Undetermined(),
            GraphIrStereoCoset::Lit(index) => Self::Lit(*index),
            GraphIrStereoCoset::LitSet(members) => Self::LitSet(members.clone()),
            GraphIrStereoCoset::Term(inner) => {
                Self::Term(into_py_variant(py, StereoTerm::from_rust(py, inner)?)?)
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrStereoCoset {
        match self {
            Self::Undetermined() => GraphIrStereoCoset::Undetermined,
            Self::Lit(index) => GraphIrStereoCoset::Lit(*index),
            Self::LitSet(members) => GraphIrStereoCoset::LitSet(members.clone()),
            Self::Term(inner) => {
                GraphIrStereoCoset::Term(Box::new(inner.bind(py).borrow().to_rust(py)))
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
    Stereo(Py<StereoCoset>),
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

    /// The exact absence or stereo-coset value, or `None` when this expression is not ground.
    fn as_lit(&self, py: Python<'_>) -> Option<TetrahedralStereo> {
        self.to_rust(py).as_lit().map(TetrahedralStereo::from_rust)
    }
}

impl_py_lattice!(
    TetrahedralStereoAst,
    GraphIrTetrahedralStereoForm,
    |value: &TetrahedralStereoAst, py: Python<'_>| -> PyResult<GraphIrTetrahedralStereoForm> {
        Ok(value.to_rust(py))
    },
    |py: Python<'_>, value: GraphIrTetrahedralStereoForm| -> PyResult<TetrahedralStereoAst> {
        TetrahedralStereoAst::from_rust(py, &value)
    }
);

impl TetrahedralStereoAst {
    pub(crate) fn from_rust(py: Python<'_>, ast: &GraphIrTetrahedralStereoForm) -> PyResult<Self> {
        Ok(match ast {
            GraphIrTetrahedralStereoForm::Undetermined => Self::Undetermined(),
            GraphIrTetrahedralStereoForm::NotStereo => Self::NotStereo(),
            GraphIrTetrahedralStereoForm::Stereo(coset) => {
                Self::Stereo(into_py_variant(py, StereoCoset::from_rust(py, coset)?)?)
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrTetrahedralStereoForm {
        match self {
            Self::Undetermined() => GraphIrTetrahedralStereoForm::Undetermined,
            Self::NotStereo() => GraphIrTetrahedralStereoForm::NotStereo,
            Self::Stereo(coset) => {
                GraphIrTetrahedralStereoForm::Stereo(coset.bind(py).borrow().to_rust(py))
            }
        }
    }
}

/// Exact ground tetrahedral stereo: explicitly not stereogenic, or a literal coset.
#[pyclass(from_py_object)]
#[derive(Clone, Copy)]
pub enum TetrahedralStereo {
    NotStereo(),
    Stereo(u32),
}

#[pymethods]
impl TetrahedralStereo {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            Self::NotStereo() => ("NotStereo", 0),
            Self::Stereo(_) => ("Stereo", 1),
        };
        variant_repr(slf.bind(py).as_any(), "TetrahedralStereo", variant, arity)
    }
}

impl TetrahedralStereo {
    pub(crate) fn from_rust(stereo: GraphIrTetrahedralStereo) -> Self {
        match stereo {
            GraphIrTetrahedralStereo::NotStereo => Self::NotStereo(),
            GraphIrTetrahedralStereo::Stereo(coset) => Self::Stereo(coset),
        }
    }

    pub(crate) fn to_rust(self) -> GraphIrTetrahedralStereo {
        match self {
            Self::NotStereo() => GraphIrTetrahedralStereo::NotStereo,
            Self::Stereo(coset) => GraphIrTetrahedralStereo::Stereo(coset),
        }
    }
}

/// Named tetrahedral configuration shorthand: counterclockwise (`Ccw`, coset
/// `Th0`) or clockwise (`Cw`, coset `Th1`).
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub enum TetrahedralConfiguration {
    Ccw,
    Cw,
}

impl TetrahedralConfiguration {
    pub(crate) fn to_rust(self) -> GraphIrTetrahedralConfiguration {
        match self {
            Self::Ccw => GraphIrTetrahedralConfiguration::Ccw,
            Self::Cw => GraphIrTetrahedralConfiguration::Cw,
        }
    }
}

/// Cis/trans bond stereo: undetermined, explicitly not stereogenic, or a stereo coset.
#[pyclass]
pub enum CisTransStereoAst {
    Undetermined(),
    NotStereo(),
    Stereo(Py<StereoCoset>),
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

    /// The exact absence or stereo-coset value, or `None` when this expression is not ground.
    fn as_lit(&self, py: Python<'_>) -> Option<CisTransStereo> {
        self.to_rust(py).as_lit().map(CisTransStereo::from_rust)
    }
}

impl_py_lattice!(
    CisTransStereoAst,
    GraphIrCisTransStereoForm,
    |value: &CisTransStereoAst, py: Python<'_>| -> PyResult<GraphIrCisTransStereoForm> {
        Ok(value.to_rust(py))
    },
    |py: Python<'_>, value: GraphIrCisTransStereoForm| -> PyResult<CisTransStereoAst> {
        CisTransStereoAst::from_rust(py, &value)
    }
);

impl CisTransStereoAst {
    pub(crate) fn from_rust(py: Python<'_>, ast: &GraphIrCisTransStereoForm) -> PyResult<Self> {
        Ok(match ast {
            GraphIrCisTransStereoForm::Undetermined => Self::Undetermined(),
            GraphIrCisTransStereoForm::NotStereo => Self::NotStereo(),
            GraphIrCisTransStereoForm::Stereo(coset) => {
                Self::Stereo(into_py_variant(py, StereoCoset::from_rust(py, coset)?)?)
            }
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrCisTransStereoForm {
        match self {
            Self::Undetermined() => GraphIrCisTransStereoForm::Undetermined,
            Self::NotStereo() => GraphIrCisTransStereoForm::NotStereo,
            Self::Stereo(coset) => {
                GraphIrCisTransStereoForm::Stereo(coset.bind(py).borrow().to_rust(py))
            }
        }
    }
}

/// Exact ground cis/trans stereo: explicitly not stereogenic, or a literal coset.
#[pyclass(from_py_object)]
#[derive(Clone, Copy)]
pub enum CisTransStereo {
    NotStereo(),
    Stereo(u32),
}

#[pymethods]
impl CisTransStereo {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            Self::NotStereo() => ("NotStereo", 0),
            Self::Stereo(_) => ("Stereo", 1),
        };
        variant_repr(slf.bind(py).as_any(), "CisTransStereo", variant, arity)
    }
}

impl CisTransStereo {
    pub(crate) fn from_rust(stereo: GraphIrCisTransStereo) -> Self {
        match stereo {
            GraphIrCisTransStereo::NotStereo => Self::NotStereo(),
            GraphIrCisTransStereo::Stereo(coset) => Self::Stereo(coset),
        }
    }

    pub(crate) fn to_rust(self) -> GraphIrCisTransStereo {
        match self {
            Self::NotStereo() => GraphIrCisTransStereo::NotStereo,
            Self::Stereo(coset) => GraphIrCisTransStereo::Stereo(coset),
        }
    }
}

/// Named cis/trans configuration shorthand: `Z` (coset `Ct0`) or `E` (coset `Ct1`),
/// named for the chemistry keywords.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub enum CisTransConfiguration {
    Z,
    E,
}

impl CisTransConfiguration {
    pub(crate) fn to_rust(self) -> GraphIrCisTransConfiguration {
        match self {
            Self::Z => GraphIrCisTransConfiguration::Z,
            Self::E => GraphIrCisTransConfiguration::E,
        }
    }
}

/// Setter coercion for `cis_trans_stereo`: `False` → not stereogenic, a
/// `CisTransConfiguration` (`Z`/`E`) → that coset, or a `CisTransStereoAst` passthrough.
#[derive(FromPyObject)]
pub(crate) enum CisTransStereoLike {
    Flag(bool),
    Config(CisTransConfiguration),
    Ast(Py<CisTransStereoAst>),
}

impl CisTransStereoLike {
    pub(crate) fn to_rust(&self, py: Python<'_>) -> PyResult<GraphIrCisTransStereoForm> {
        Ok(match self {
            CisTransStereoLike::Flag(false) => GraphIrCisTransStereoForm::NotStereo,
            CisTransStereoLike::Flag(true) => return Err(PyValueError::new_err(
                "cis_trans_stereo = True is not meaningful; use CisTransConfiguration.Z/E or False",
            )),
            CisTransStereoLike::Config(cts) => cts.to_rust().into(),
            CisTransStereoLike::Ast(a) => a.bind(py).borrow().to_rust(py),
        })
    }
}

/// The coordination geometry of a stereo site. A fieldless, hashable value enum whose
/// members correspond exactly to the Rust `StereoKind`.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum StereoKind {
    Tetrahedral,
    CisTrans,
    Axial,
    SquarePlanar,
    TrigonalBipyramidal,
    Octahedral,
}

impl StereoKind {
    pub(crate) fn from_rust(ast: GraphIrStereoKind) -> Self {
        match ast {
            GraphIrStereoKind::Tetrahedral => Self::Tetrahedral,
            GraphIrStereoKind::CisTrans => Self::CisTrans,
            GraphIrStereoKind::Axial => Self::Axial,
            GraphIrStereoKind::SquarePlanar => Self::SquarePlanar,
            GraphIrStereoKind::TrigonalBipyramidal => Self::TrigonalBipyramidal,
            GraphIrStereoKind::Octahedral => Self::Octahedral,
        }
    }

    pub(crate) fn to_rust(self) -> GraphIrStereoKind {
        match self {
            Self::Tetrahedral => GraphIrStereoKind::Tetrahedral,
            Self::CisTrans => GraphIrStereoKind::CisTrans,
            Self::Axial => GraphIrStereoKind::Axial,
            Self::SquarePlanar => GraphIrStereoKind::SquarePlanar,
            Self::TrigonalBipyramidal => GraphIrStereoKind::TrigonalBipyramidal,
            Self::Octahedral => GraphIrStereoKind::Octahedral,
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
    pub(crate) fn from_rust(ast: GraphIrStereoLigandKind) -> Self {
        match ast {
            GraphIrStereoLigandKind::Atom => Self::Atom,
            GraphIrStereoLigandKind::ImplicitHydrogen => Self::ImplicitHydrogen,
            GraphIrStereoLigandKind::LonePair => Self::LonePair,
        }
    }

    pub(crate) fn to_rust(self) -> GraphIrStereoLigandKind {
        match self {
            Self::Atom => GraphIrStereoLigandKind::Atom,
            Self::ImplicitHydrogen => GraphIrStereoLigandKind::ImplicitHydrogen,
            Self::LonePair => GraphIrStereoLigandKind::LonePair,
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
    pub(crate) fn from_rust(ast: GraphIrTopicity) -> Self {
        match ast {
            GraphIrTopicity::Homotopic => Self::Homotopic,
            GraphIrTopicity::Enantiotopic => Self::Enantiotopic,
            GraphIrTopicity::Diastereotopic => Self::Diastereotopic,
        }
    }

    pub(crate) fn to_rust(self) -> GraphIrTopicity {
        match self {
            Self::Homotopic => GraphIrTopicity::Homotopic,
            Self::Enantiotopic => GraphIrTopicity::Enantiotopic,
            Self::Diastereotopic => GraphIrTopicity::Diastereotopic,
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
    pub(crate) fn from_rust(ast: GraphIrStereogenicity) -> Self {
        match ast {
            GraphIrStereogenicity::Symmetric => Self::Symmetric,
            GraphIrStereogenicity::Prochiral => Self::Prochiral,
            GraphIrStereogenicity::Stereogenic => Self::Stereogenic,
        }
    }

    pub(crate) fn to_rust(self) -> GraphIrStereogenicity {
        match self {
            Self::Symmetric => GraphIrStereogenicity::Symmetric,
            Self::Prochiral => GraphIrStereogenicity::Prochiral,
            Self::Stereogenic => GraphIrStereogenicity::Stereogenic,
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
    pub(crate) fn from_rust(ast: GraphIrStereoLigand) -> Self {
        StereoLigand {
            atom_id: ast.atom_id.0,
            kind: StereoLigandKind::from_rust(ast.kind),
        }
    }

    pub(crate) fn to_rust(self) -> GraphIrStereoLigand {
        GraphIrStereoLigand::new(GraphIrAtomId(self.atom_id), self.kind.to_rust())
    }
}

/// A stereo configuration: undetermined (geometry not yet known, so no coset), or `Kinded`
/// — a concrete coordination geometry bound to a coset that may still be open. Corresponds to the
/// Rust `StereoConfigurationForm`; `Undetermined` and `Kinded(Tetrahedral, Undetermined)` are
/// distinct.
#[pyclass]
pub enum StereoConfigurationAst {
    Undetermined(),
    Kinded(StereoKind, Py<StereoCoset>),
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
    fn coset(&self, py: Python<'_>) -> Option<Py<StereoCoset>> {
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

/// An update to a stereo configuration.
///
/// `Unchanged` omits the field, `Undetermined` clears it, and `Kinded` sets the
/// geometry while optionally leaving its coset unchanged.
#[pyclass]
pub enum StereoConfigurationUpdate {
    Unchanged(),
    Undetermined(),
    Kinded(StereoKind, Option<Py<StereoCoset>>),
}

#[pymethods]
impl StereoConfigurationUpdate {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            Self::Unchanged() => ("Unchanged", 0),
            Self::Undetermined() => ("Undetermined", 0),
            Self::Kinded(_, _) => ("Kinded", 2),
        };
        variant_repr(
            slf.bind(py).as_any(),
            "StereoConfigurationUpdate",
            variant,
            arity,
        )
    }
}

/// Attribute updates for a stereo atom.
#[pyclass(frozen, skip_from_py_object)]
#[derive(Clone)]
pub struct StereoAtomUpdate(GraphIrStereoAtomUpdate);

#[pymethods]
impl StereoAtomUpdate {
    #[new]
    #[pyo3(signature = (*, configuration=None, constraints=None))]
    fn new(
        py: Python<'_>,
        configuration: Option<PyRef<'_, StereoConfigurationUpdate>>,
        constraints: Option<Py<StereoAtomConstraintsAst>>,
    ) -> Self {
        Self::from_rust(&GraphIrStereoAtomUpdate {
            configuration: configuration
                .map(|value| value.to_rust(py))
                .unwrap_or_default(),
            constraints: constraints
                .map(|value| value.bind(py).borrow().inner().clone())
                .unwrap_or_default(),
        })
    }

    /// Parse a stereo-atom-update DSL string into a `StereoAtomUpdate`.
    #[staticmethod]
    fn parse(s: &str) -> PyResult<Self> {
        GraphIrStereoAtomUpdate::from_str(s)
            .map(Self)
            .map_err(parse_error)
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("StereoAtomUpdate.parse('{}')", self.0)
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    #[getter]
    fn configuration(&self, py: Python<'_>) -> PyResult<StereoConfigurationUpdate> {
        StereoConfigurationUpdate::from_rust(py, &self.0.configuration)
    }

    #[getter]
    fn constraints(&self) -> StereoAtomConstraintsAst {
        StereoAtomConstraintsAst::from_inner(self.0.constraints.clone())
    }
}

impl StereoAtomUpdate {
    pub(crate) fn from_rust(update: &GraphIrStereoAtomUpdate) -> Self {
        Self(update.clone())
    }

    pub(crate) fn to_rust(&self) -> GraphIrStereoAtomUpdate {
        self.0.clone()
    }
}

/// Attribute updates for a stereo bond.
#[pyclass(frozen, skip_from_py_object)]
#[derive(Clone)]
pub struct StereoBondUpdate(GraphIrStereoBondUpdate);

#[pymethods]
impl StereoBondUpdate {
    #[new]
    #[pyo3(signature = (*, configuration=None, constraints=None))]
    fn new(
        py: Python<'_>,
        configuration: Option<PyRef<'_, StereoConfigurationUpdate>>,
        constraints: Option<Py<StereoBondConstraintsAst>>,
    ) -> Self {
        Self::from_rust(&GraphIrStereoBondUpdate {
            configuration: configuration
                .map(|value| value.to_rust(py))
                .unwrap_or_default(),
            constraints: constraints
                .map(|value| value.bind(py).borrow().inner().clone())
                .unwrap_or_default(),
        })
    }

    /// Parse a stereo-bond-update DSL string into a `StereoBondUpdate`.
    #[staticmethod]
    fn parse(s: &str) -> PyResult<Self> {
        GraphIrStereoBondUpdate::from_str(s)
            .map(Self)
            .map_err(parse_error)
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("StereoBondUpdate.parse('{}')", self.0)
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    #[getter]
    fn configuration(&self, py: Python<'_>) -> PyResult<StereoConfigurationUpdate> {
        StereoConfigurationUpdate::from_rust(py, &self.0.configuration)
    }

    #[getter]
    fn constraints(&self) -> StereoBondConstraintsAst {
        StereoBondConstraintsAst::from_inner(self.0.constraints.clone())
    }
}

impl StereoBondUpdate {
    pub(crate) fn from_rust(update: &GraphIrStereoBondUpdate) -> Self {
        Self(update.clone())
    }

    pub(crate) fn to_rust(&self) -> GraphIrStereoBondUpdate {
        self.0.clone()
    }
}

impl StereoConfigurationUpdate {
    pub(crate) fn from_rust(
        py: Python<'_>,
        update: &GraphIrStereoConfigurationUpdate,
    ) -> PyResult<Self> {
        Ok(match update {
            GraphIrStereoConfigurationUpdate::Unchanged => Self::Unchanged(),
            GraphIrStereoConfigurationUpdate::Undetermined => Self::Undetermined(),
            GraphIrStereoConfigurationUpdate::Kinded { kind, coset } => Self::Kinded(
                StereoKind::from_rust(*kind),
                coset
                    .as_ref()
                    .map(|coset| {
                        StereoCoset::from_rust(py, coset)
                            .and_then(|coset| into_py_variant(py, coset))
                    })
                    .transpose()?,
            ),
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrStereoConfigurationUpdate {
        match self {
            Self::Unchanged() => GraphIrStereoConfigurationUpdate::Unchanged,
            Self::Undetermined() => GraphIrStereoConfigurationUpdate::Undetermined,
            Self::Kinded(kind, coset) => GraphIrStereoConfigurationUpdate::Kinded {
                kind: kind.to_rust(),
                coset: coset
                    .as_ref()
                    .map(|coset| coset.bind(py).borrow().to_rust(py)),
            },
        }
    }
}

impl_py_lattice!(
    StereoConfigurationAst,
    GraphIrStereoConfigurationForm,
    |value: &StereoConfigurationAst, py: Python<'_>| -> PyResult<GraphIrStereoConfigurationForm> {
        Ok(value.to_rust(py))
    },
    |py: Python<'_>, value: GraphIrStereoConfigurationForm| -> PyResult<StereoConfigurationAst> {
        StereoConfigurationAst::from_rust(py, &value)
    }
);

impl StereoConfigurationAst {
    pub(crate) fn from_rust(
        py: Python<'_>,
        ast: &GraphIrStereoConfigurationForm,
    ) -> PyResult<Self> {
        Ok(match ast {
            GraphIrStereoConfigurationForm::Undetermined => Self::Undetermined(),
            GraphIrStereoConfigurationForm::Kinded(kind, coset) => Self::Kinded(
                StereoKind::from_rust(*kind),
                into_py_variant(py, StereoCoset::from_rust(py, coset)?)?,
            ),
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrStereoConfigurationForm {
        match self {
            Self::Undetermined() => GraphIrStereoConfigurationForm::Undetermined,
            Self::Kinded(kind, coset) => GraphIrStereoConfigurationForm::Kinded(
                kind.to_rust(),
                coset.bind(py).borrow().to_rust(py),
            ),
        }
    }
}

/// Setter coercion for a stereo `configuration` field: the `TetrahedralConfiguration`
/// (`Ccw`/`Cw`) or `CisTransConfiguration` (`Z`/`E`) per-kind coset shorthand, or a
/// `StereoConfigurationAst`
/// passthrough. Axial/square-planar/etc. have no shorthand — use the full `Kinded` form.
#[derive(FromPyObject)]
pub(crate) enum StereoConfigurationLike {
    Tetrahedral(TetrahedralConfiguration),
    CisTrans(CisTransConfiguration),
    Ast(Py<StereoConfigurationAst>),
}

impl StereoConfigurationLike {
    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrStereoConfigurationForm {
        match self {
            StereoConfigurationLike::Tetrahedral(t) => {
                let coset = match t.to_rust() {
                    GraphIrTetrahedralConfiguration::Ccw => 0,
                    GraphIrTetrahedralConfiguration::Cw => 1,
                };
                GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(coset),
                )
            }
            StereoConfigurationLike::CisTrans(c) => {
                let coset = match c.to_rust() {
                    GraphIrCisTransConfiguration::Z => 0,
                    GraphIrCisTransConfiguration::E => 1,
                };
                GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::CisTrans,
                    GraphIrStereoCoset::Lit(coset),
                )
            }
            StereoConfigurationLike::Ast(a) => a.bind(py).borrow().to_rust(py),
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
    pub(crate) fn from_rust(ast: GraphIrLigandPermutation) -> Self {
        LigandPermutation {
            permutation: Permutation::from_inner(ast.0),
        }
    }

    pub(crate) fn to_rust(self) -> GraphIrLigandPermutation {
        GraphIrLigandPermutation(self.permutation.inner())
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
    pub(crate) fn from_rust(ast: GraphIrOrientedLigandPermutation) -> Self {
        OrientedLigandPermutation {
            permutation: LigandPermutation::from_rust(ast.permutation),
            orientation: Orientation::from_rust(ast.orientation),
        }
    }

    pub(crate) fn to_rust(self) -> GraphIrOrientedLigandPermutation {
        GraphIrOrientedLigandPermutation {
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
        StereoLigandPair::from_rust(GraphIrStereoLigandPair::new(
            GraphIrStereoLigandPosition(a),
            GraphIrStereoLigandPosition(b),
        ))
    }

    pub(crate) fn __repr__(&self) -> String {
        format!("StereoLigandPair({}, {})", self.first, self.second)
    }
}

impl StereoLigandPair {
    pub(crate) fn from_rust(ast: GraphIrStereoLigandPair) -> Self {
        StereoLigandPair {
            first: ast.first().0,
            second: ast.second().0,
        }
    }

    pub(crate) fn to_rust(self) -> GraphIrStereoLigandPair {
        GraphIrStereoLigandPair::new(
            GraphIrStereoLigandPosition(self.first),
            GraphIrStereoLigandPosition(self.second),
        )
    }
}

#[cfg(test)]
use crate::constraint::stereo::{
    FluxionalityAst, LigandSymmetryAst, StereoAtomConstraintKey, StereoAtomConstraintsUpdate,
    StereogenicityAst, TopicityAst, TopicityRelationAst, TopicityRelationLike,
};
use crate::constraint::stereo::{
    StereoAtomConstraintAst, StereoAtomConstraintsAst, StereoAtomConstraintsBacking,
    StereoAtomConstraintsLike, StereoAtomConstraintsView, StereoBondConstraintAst,
    StereoBondConstraintsAst, StereoBondConstraintsBacking, StereoBondConstraintsLike,
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
        $value:ident, $ast_value:ident, $constraint:ident, $constraints:ident, $like:ident,
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
                configuration: StereoConfigurationLike,
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
            fn set_configuration(&mut self, py: Python<'_>, value: StereoConfigurationLike) {
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
            fn set_constraints(slf: Py<Self>, py: Python<'_>, value: $like) -> PyResult<()> {
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

        impl_py_lattice!(
            $value,
            $ast_value,
            |value: &$value, _py: Python<'_>| -> PyResult<$ast_value> {
                Ok(value.inner().clone())
            },
            |_py: Python<'_>, value: $ast_value| -> PyResult<$value> {
                Ok($value::from_inner(value))
            }
        );
    };
}

stereo_value! {
    StereoAtomAst, GraphIrStereoAtomForm, StereoAtomConstraintAst, StereoAtomConstraintsAst,
    StereoAtomConstraintsLike, StereoAtomConstraintsView, StereoAtomConstraintsBacking, production,
}

stereo_value! {
    StereoBondAst, GraphIrStereoBondForm, StereoBondConstraintAst, StereoBondConstraintsAst,
    StereoBondConstraintsLike, StereoBondConstraintsView, StereoBondConstraintsBacking, production,
}

/// Per-entity molecule-embedded stereo view — `StereoAtomView` / `StereoBondView` — a handle
/// to the molecule plus the entity's id. Field reads rebuild the transient Rust view; the
/// molecule is never copied. The site atom/bond and ligands are read-only topology; the
/// configuration and constraints are the mutable value.
macro_rules! stereo_view {
    (
        $view:ident, $ast_view:ident, $ast_id:ident, $namespace:ident, $entity_mut:ident,
        $id_error:literal, $constraint:ident, $constraints_view:ident, $constraints_backing:ident,
        $like:ident $(,)?
    ) => {
        #[pyclass]
        pub struct $view {
            owner: Py<MoleculeAst>,
            id: $ast_id,
        }

        impl $view {
            /// Rebuild the transient AST view for this entity, or `IndexError` if the id is
            /// no longer present.
            fn view<'a>(&self, molecule: &'a GraphIrMolecule) -> PyResult<$ast_view<'a>> {
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
            fn coset(&self, py: Python<'_>) -> PyResult<StereoCoset> {
                let molecule = self.owner.bind(py).borrow();
                StereoCoset::from_rust(py, self.view(molecule.inner())?.coset())
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
            fn set_configuration(&self, py: Python<'_>, value: StereoConfigurationLike) {
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
            fn set_constraints(&self, py: Python<'_>, value: $like) -> PyResult<()> {
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
    StereoAtomView, GraphIrStereoAtomView, GraphIrStereoAtomId, stereo_atoms, stereo_atom_mut,
    "stereo atom id out of range", StereoAtomConstraintAst, StereoAtomConstraintsView,
    StereoAtomConstraintsBacking, StereoAtomConstraintsLike,
}

stereo_view! {
    StereoBondView, GraphIrStereoBondView, GraphIrStereoBondId, stereo_bonds, stereo_bond_mut,
    "stereo bond id out of range", StereoBondConstraintAst, StereoBondConstraintsView,
    StereoBondConstraintsBacking, StereoBondConstraintsLike,
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
        fn $resolve_index(molecule: &GraphIrMolecule, index: isize) -> PyResult<$ast_id> {
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
                let ligands: Vec<GraphIrStereoLigand> =
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
    StereoAtomViews, StereoAtomView, StereoAtomViewIter, GraphIrStereoAtomId, GraphIrAtomId, stereo_atoms,
    stereo_atom_mut, StereoAtomAst, resolve_stereo_atom_index, "stereo atom id out of range",
}

stereo_views! {
    StereoBondViews, StereoBondView, StereoBondViewIter, GraphIrStereoBondId, GraphIrBondId, stereo_bonds,
    stereo_bond_mut, StereoBondAst, resolve_stereo_bond_index, "stereo bond id out of range",
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_chem::element::Element as ChemElement;
    use umol_graph_ir::ir::{
        AtomForm as GraphIrAtomForm, BondForm as GraphIrBondForm,
        FluxionalityForm as GraphIrFluxionalityForm,
        LigandSymmetryForm as GraphIrLigandSymmetryForm, MoleculeEntries,
        StereoAtomConstraintForm as GraphIrStereoAtomConstraintForm,
        StereoAtomConstraintKey as GraphIrStereoAtomConstraintKey,
        StereoAtomConstraintsForm as GraphIrStereoAtomConstraintsForm,
        StereoBondConstraintForm as GraphIrStereoBondConstraintForm,
        StereoBondConstraintsForm as GraphIrStereoBondConstraintsForm,
        StereoLigandPair as GraphIrStereoLigandPair,
        StereoLigandPosition as GraphIrStereoLigandPosition,
        StereogenicityForm as GraphIrStereogenicityForm, TopicityForm as GraphIrTopicityForm,
        TopicityRelationForm as GraphIrTopicityRelationForm,
    };

    use super::*;
    use crate::boolean::{BooleanAst, BooleanLike};

    #[rstest]
    #[case(vec![0, 1, 2, 3])]
    #[case(vec![1, 0, 2, 3])]
    #[case(vec![2, 0, 1])]
    fn test_permutation_image(#[case] image: Vec<usize>) {
        let permutation = Permutation::new(image.clone()).unwrap();
        assert_eq!(permutation.image(), image);
        assert_eq!(permutation.degree(), image.len());
    }

    #[rstest]
    fn test_permutation_identity() {
        assert_eq!(Permutation::identity(4).image(), vec![0, 1, 2, 3]);
    }

    #[rstest]
    #[case(GraphIrStereoTerm::Lit(1))]
    #[case(GraphIrStereoTerm::LitSet(BTreeSet::from([0, 2])))]
    #[case(GraphIrStereoTerm::Var(Box::new(("x".to_string(), None))))]
    #[case(GraphIrStereoTerm::Var(Box::new(("y".to_string(), Some(BTreeSet::from([0, 1]))))))]
    #[case(GraphIrStereoTerm::Swap(Box::new(GraphIrStereoTerm::Lit(0))))]
    #[case(GraphIrStereoTerm::Mirror(Box::new(GraphIrStereoTerm::Lit(0))))]
    #[case(GraphIrStereoTerm::Apply(Box::new(GraphIrStereoTerm::Lit(0)), PermPermutation::from_image(&[1, 0, 2, 3])))]
    fn test_stereo_term_roundtrip(#[case] ast: GraphIrStereoTerm) {
        Python::attach(|py| {
            assert_eq!(StereoTerm::from_rust(py, &ast).unwrap().to_rust(py), ast);
        });
    }

    #[rstest]
    #[case(GraphIrStereoCoset::Undetermined)]
    #[case(GraphIrStereoCoset::Lit(1))]
    #[case(GraphIrStereoCoset::LitSet(BTreeSet::from([0, 1])))]
    #[case(GraphIrStereoCoset::Term(Box::new(GraphIrStereoTerm::Lit(1))))]
    fn test_stereo_coset_roundtrip(#[case] ast: GraphIrStereoCoset) {
        Python::attach(|py| {
            assert_eq!(StereoCoset::from_rust(py, &ast).unwrap().to_rust(py), ast);
        });
    }

    #[rstest]
    #[case(GraphIrTetrahedralStereoForm::Undetermined)]
    #[case(GraphIrTetrahedralStereoForm::NotStereo)]
    #[case(GraphIrTetrahedralStereoForm::Stereo(GraphIrStereoCoset::Lit(1)))]
    #[case(
        GraphIrTetrahedralStereoForm::Stereo(GraphIrStereoCoset::Term(Box::new(
            GraphIrStereoTerm::Lit(0)
        )))
    )]
    fn test_tetrahedral_stereo_ast_roundtrip(#[case] ast: GraphIrTetrahedralStereoForm) {
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
    #[case(
        GraphIrTetrahedralStereoForm::NotStereo,
        Some(GraphIrTetrahedralStereo::NotStereo)
    )]
    #[case(
        GraphIrTetrahedralStereoForm::Stereo(GraphIrStereoCoset::Lit(2)),
        Some(GraphIrTetrahedralStereo::Stereo(2))
    )]
    #[case(GraphIrTetrahedralStereoForm::Undetermined, None)]
    #[case(GraphIrTetrahedralStereoForm::Stereo(GraphIrStereoCoset::LitSet(BTreeSet::from([0, 1]))), None)]
    fn test_tetrahedral_stereo_ast_as_lit(
        #[case] ast: GraphIrTetrahedralStereoForm,
        #[case] expected: Option<GraphIrTetrahedralStereo>,
    ) {
        Python::attach(|py| {
            assert_eq!(
                TetrahedralStereoAst::from_rust(py, &ast)
                    .unwrap()
                    .as_lit(py)
                    .map(TetrahedralStereo::to_rust),
                expected
            );
        });
    }

    #[rstest]
    #[case(TetrahedralConfiguration::Ccw, GraphIrTetrahedralConfiguration::Ccw)]
    #[case(TetrahedralConfiguration::Cw, GraphIrTetrahedralConfiguration::Cw)]
    fn test_tetrahedral_configuration_to_rust(
        #[case] config: TetrahedralConfiguration,
        #[case] expected: GraphIrTetrahedralConfiguration,
    ) {
        assert_eq!(config.to_rust(), expected);
    }

    #[rstest]
    #[case(GraphIrCisTransStereoForm::Undetermined)]
    #[case(GraphIrCisTransStereoForm::NotStereo)]
    #[case(GraphIrCisTransStereoForm::Stereo(GraphIrStereoCoset::Lit(1)))]
    #[case(GraphIrCisTransStereoForm::Stereo(GraphIrStereoCoset::Term(Box::new(
        GraphIrStereoTerm::Lit(0)
    ))))]
    fn test_cis_trans_stereo_ast_roundtrip(#[case] ast: GraphIrCisTransStereoForm) {
        Python::attach(|py| {
            assert_eq!(
                CisTransStereoAst::from_rust(py, &ast).unwrap().to_rust(py),
                ast
            );
        });
    }

    #[rstest]
    #[case(
        GraphIrCisTransStereoForm::NotStereo,
        Some(GraphIrCisTransStereo::NotStereo)
    )]
    #[case(
        GraphIrCisTransStereoForm::Stereo(GraphIrStereoCoset::Lit(1)),
        Some(GraphIrCisTransStereo::Stereo(1))
    )]
    #[case(GraphIrCisTransStereoForm::Undetermined, None)]
    #[case(GraphIrCisTransStereoForm::Stereo(GraphIrStereoCoset::LitSet(BTreeSet::from([0, 1]))), None)]
    fn test_cis_trans_stereo_ast_as_lit(
        #[case] ast: GraphIrCisTransStereoForm,
        #[case] expected: Option<GraphIrCisTransStereo>,
    ) {
        Python::attach(|py| {
            assert_eq!(
                CisTransStereoAst::from_rust(py, &ast)
                    .unwrap()
                    .as_lit(py)
                    .map(CisTransStereo::to_rust),
                expected
            );
        });
    }

    #[rstest]
    #[case(CisTransConfiguration::Z, GraphIrCisTransConfiguration::Z)]
    #[case(CisTransConfiguration::E, GraphIrCisTransConfiguration::E)]
    fn test_cis_trans_configuration_to_rust(
        #[case] config: CisTransConfiguration,
        #[case] expected: GraphIrCisTransConfiguration,
    ) {
        assert_eq!(config.to_rust(), expected);
    }

    #[rstest]
    #[case(GraphIrStereoKind::Tetrahedral)]
    #[case(GraphIrStereoKind::CisTrans)]
    #[case(GraphIrStereoKind::Axial)]
    #[case(GraphIrStereoKind::SquarePlanar)]
    #[case(GraphIrStereoKind::TrigonalBipyramidal)]
    #[case(GraphIrStereoKind::Octahedral)]
    fn test_stereo_kind_roundtrip(#[case] ast: GraphIrStereoKind) {
        assert_eq!(StereoKind::from_rust(ast).to_rust(), ast);
    }

    #[rstest]
    #[case(GraphIrStereoLigandKind::Atom)]
    #[case(GraphIrStereoLigandKind::ImplicitHydrogen)]
    #[case(GraphIrStereoLigandKind::LonePair)]
    fn test_stereo_ligand_kind_roundtrip(#[case] ast: GraphIrStereoLigandKind) {
        assert_eq!(StereoLigandKind::from_rust(ast).to_rust(), ast);
    }

    #[rstest]
    #[case(GraphIrTopicity::Homotopic)]
    #[case(GraphIrTopicity::Enantiotopic)]
    #[case(GraphIrTopicity::Diastereotopic)]
    fn test_topicity_roundtrip(#[case] ast: GraphIrTopicity) {
        assert_eq!(Topicity::from_rust(ast).to_rust(), ast);
    }

    #[rstest]
    #[case(GraphIrStereogenicity::Symmetric)]
    #[case(GraphIrStereogenicity::Prochiral)]
    #[case(GraphIrStereogenicity::Stereogenic)]
    fn test_stereogenicity_roundtrip(#[case] ast: GraphIrStereogenicity) {
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
    #[case(GraphIrStereoLigand::new(GraphIrAtomId(0), GraphIrStereoLigandKind::Atom))]
    #[case(GraphIrStereoLigand::new(GraphIrAtomId(5), GraphIrStereoLigandKind::LonePair))]
    fn test_stereo_ligand_roundtrip(#[case] ast: GraphIrStereoLigand) {
        assert_eq!(StereoLigand::from_rust(ast).to_rust(), ast);
    }

    #[rstest]
    fn test_stereo_configuration_ast_roundtrip() {
        Python::attach(|py| {
            for ast in [
                GraphIrStereoConfigurationForm::Undetermined,
                GraphIrStereoConfigurationForm::kinded(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(1),
                ),
                GraphIrStereoConfigurationForm::kinded(
                    GraphIrStereoKind::Octahedral,
                    GraphIrStereoCoset::Undetermined,
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
    #[case::unchanged(GraphIrStereoConfigurationUpdate::Unchanged)]
    #[case::undetermined(GraphIrStereoConfigurationUpdate::Undetermined)]
    #[case::kind_only(GraphIrStereoConfigurationUpdate::Kinded {
        kind: GraphIrStereoKind::Tetrahedral,
        coset: None,
    })]
    #[case::undetermined_coset(GraphIrStereoConfigurationUpdate::Kinded {
        kind: GraphIrStereoKind::Tetrahedral,
        coset: Some(GraphIrStereoCoset::Undetermined),
    })]
    #[case::absolute(GraphIrStereoConfigurationUpdate::Kinded {
        kind: GraphIrStereoKind::CisTrans,
        coset: Some(GraphIrStereoCoset::Lit(1)),
    })]
    fn test_stereo_configuration_update_roundtrip(
        #[case] update: GraphIrStereoConfigurationUpdate,
    ) {
        Python::attach(|py| {
            assert_eq!(
                StereoConfigurationUpdate::from_rust(py, &update)
                    .unwrap()
                    .to_rust(py),
                update
            );
        });
    }

    #[rstest]
    fn test_stereo_configuration_ast_kind_coset() {
        Python::attach(|py| {
            let coset = into_py_variant(py, StereoCoset::Lit(1)).unwrap();
            let config = StereoConfigurationAst::Kinded(StereoKind::Tetrahedral, coset);
            assert_eq!(config.kind(), Some(StereoKind::Tetrahedral));
            assert_eq!(
                config.coset(py).unwrap().bind(py).borrow().to_rust(py),
                GraphIrStereoCoset::Lit(1)
            );
            let undetermined = StereoConfigurationAst::Undetermined();
            assert_eq!(undetermined.kind(), None);
            assert!(undetermined.coset(py).is_none());
        });
    }

    #[rstest]
    fn test_stereo_configuration_like_to_rust() {
        Python::attach(|py| {
            // the Th shorthand → Kinded(Tetrahedral, coset)
            assert_eq!(
                StereoConfigurationLike::Tetrahedral(TetrahedralConfiguration::Cw).to_rust(py),
                GraphIrStereoConfigurationForm::kinded(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(1)
                )
            );
            // the Ct shorthand → Kinded(CisTrans, coset)
            assert_eq!(
                StereoConfigurationLike::CisTrans(CisTransConfiguration::E).to_rust(py),
                GraphIrStereoConfigurationForm::kinded(
                    GraphIrStereoKind::CisTrans,
                    GraphIrStereoCoset::Lit(1)
                )
            );
            // a StereoConfigurationAst passes through
            let config = Py::new(py, StereoConfigurationAst::Undetermined()).unwrap();
            assert_eq!(
                StereoConfigurationLike::Ast(config).to_rust(py),
                GraphIrStereoConfigurationForm::Undetermined
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
        let ligand_permutation =
            LigandPermutation::new(Permutation::new(vec![1, 0, 2, 3]).unwrap());
        assert_eq!(ligand_permutation.permutation.image(), vec![1, 0, 2, 3]);
        assert_eq!(
            ligand_permutation.__repr__(),
            "LigandPermutation([1, 0, 2, 3])"
        );
    }

    #[rstest]
    #[case(GraphIrLigandPermutation(PermPermutation::identity(4)))]
    #[case(GraphIrLigandPermutation(PermPermutation::from_image(&[1, 0, 2, 3])))]
    fn test_ligand_permutation_roundtrip(#[case] ast: GraphIrLigandPermutation) {
        assert_eq!(LigandPermutation::from_rust(ast).to_rust(), ast);
    }

    #[rstest]
    #[case::equal(vec![1, 0, 2, 3], vec![1, 0, 2, 3], true)]
    #[case::different(vec![1, 0, 2, 3], vec![0, 1, 2, 3], false)]
    fn test_ligand_permutation_matches(
        #[case] a: Vec<usize>,
        #[case] b: Vec<usize>,
        #[case] expected: bool,
    ) {
        let a = LigandPermutation::new(Permutation::new(a).unwrap());
        let b = LigandPermutation::new(Permutation::new(b).unwrap());
        assert_eq!(a.matches(&b), expected);
    }

    #[rstest]
    #[case::equal(vec![1, 0, 2, 3], Orientation::Proper, vec![1, 0, 2, 3], Orientation::Proper, true)]
    #[case::different_orientation(vec![1, 0, 2, 3], Orientation::Proper, vec![1, 0, 2, 3], Orientation::Improper, false)]
    #[case::different_permutation(vec![1, 0, 2, 3], Orientation::Proper, vec![0, 1, 2, 3], Orientation::Proper, false)]
    fn test_oriented_ligand_permutation_matches(
        #[case] a_permutation: Vec<usize>,
        #[case] a_orientation: Orientation,
        #[case] b_permutation: Vec<usize>,
        #[case] b_orientation: Orientation,
        #[case] expected: bool,
    ) {
        let a = OrientedLigandPermutation::new(
            LigandPermutation::new(Permutation::new(a_permutation).unwrap()),
            a_orientation,
        );
        let b = OrientedLigandPermutation::new(
            LigandPermutation::new(Permutation::new(b_permutation).unwrap()),
            b_orientation,
        );
        assert_eq!(a.matches(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case(GraphIrOrientedLigandPermutation { permutation: GraphIrLigandPermutation(PermPermutation::from_image(&[1, 0, 2, 3])), orientation: PermOrientation::Proper })]
    #[case(GraphIrOrientedLigandPermutation { permutation: GraphIrLigandPermutation(PermPermutation::identity(4)), orientation: PermOrientation::Improper })]
    fn test_oriented_ligand_permutation_roundtrip(#[case] ast: GraphIrOrientedLigandPermutation) {
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
    #[case(GraphIrStereoLigandPair::new(
        GraphIrStereoLigandPosition(0),
        GraphIrStereoLigandPosition(3)
    ))]
    #[case(GraphIrStereoLigandPair::new(
        GraphIrStereoLigandPosition(2),
        GraphIrStereoLigandPosition(1)
    ))]
    fn test_stereo_ligand_pair_roundtrip(#[case] ast: GraphIrStereoLigandPair) {
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
    #[case(GraphIrTopicityRelationForm::Undetermined)]
    #[case(GraphIrTopicityRelationForm::Lit(GraphIrTopicity::Homotopic))]
    #[case(GraphIrTopicityRelationForm::LitSet(BTreeSet::from([
        GraphIrTopicity::Homotopic,
        GraphIrTopicity::Enantiotopic,
    ])))]
    #[case(GraphIrTopicityRelationForm::NotSet(BTreeSet::from([GraphIrTopicity::Diastereotopic])))]
    fn test_topicity_relation_ast_roundtrip(#[case] ast: GraphIrTopicityRelationForm) {
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
    #[case(GraphIrStereogenicityForm::Undetermined)]
    #[case(GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic))]
    #[case(GraphIrStereogenicityForm::LitSet(BTreeSet::from([
        GraphIrStereogenicity::Symmetric,
        GraphIrStereogenicity::Prochiral,
    ])))]
    #[case(GraphIrStereogenicityForm::NotSet(BTreeSet::from([GraphIrStereogenicity::Stereogenic])))]
    fn test_stereogenicity_ast_roundtrip(#[case] ast: GraphIrStereogenicityForm) {
        assert_eq!(StereogenicityAst::from_rust(&ast).to_rust(), ast);
    }

    #[rstest]
    fn test_ligand_symmetry_ast_new() {
        Python::attach(|py| {
            let permutation = OrientedLigandPermutation::new(
                LigandPermutation::new(Permutation::new(vec![1, 0, 2, 3]).unwrap()),
                Orientation::Proper,
            );
            let value = LigandSymmetryAst::new(py, permutation, BooleanLike::Lit(true)).unwrap();
            assert!(value.permutation() == permutation);
            assert_eq!(
                value.invariant.bind(py).borrow().to_rust(),
                GraphIrBooleanForm::Lit(true)
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
                LigandPermutation::new(Permutation::new(vec![1, 0, 2, 3]).unwrap()),
                Orientation::Proper,
            );
            let other_permutation = OrientedLigandPermutation::new(
                LigandPermutation::new(Permutation::new(vec![0, 1, 2, 3]).unwrap()),
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
            assert!(wildcard.matches(py, &invariant_true).unwrap());
            assert!(!invariant_true.matches(py, &invariant_false).unwrap());
            assert!(!invariant_true.matches(py, &other).unwrap());
        });
    }

    #[rstest]
    fn test_ligand_symmetry_ast_roundtrip() {
        Python::attach(|py| {
            for ast in [
                GraphIrLigandSymmetryForm {
                    permutation: GraphIrOrientedLigandPermutation {
                        permutation: GraphIrLigandPermutation(PermPermutation::from_image(&[
                            1, 0, 2, 3,
                        ])),
                        orientation: PermOrientation::Proper,
                    },
                    invariant: GraphIrBooleanForm::Lit(true),
                },
                GraphIrLigandSymmetryForm {
                    permutation: GraphIrOrientedLigandPermutation {
                        permutation: GraphIrLigandPermutation(PermPermutation::identity(4)),
                        orientation: PermOrientation::Improper,
                    },
                    invariant: GraphIrBooleanForm::Undetermined,
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
            let permutation = LigandPermutation::new(Permutation::new(vec![1, 0, 2, 3]).unwrap());
            let value = FluxionalityAst::new(py, permutation, BooleanLike::Lit(false)).unwrap();
            assert!(value.permutation() == permutation);
            assert_eq!(
                value.active.bind(py).borrow().to_rust(),
                GraphIrBooleanForm::Lit(false)
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
            let permutation = LigandPermutation::new(Permutation::new(vec![1, 0, 2, 3]).unwrap());
            let other_permutation =
                LigandPermutation::new(Permutation::new(vec![0, 1, 2, 3]).unwrap());
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
            assert!(wildcard.matches(py, &active_true).unwrap());
            assert!(!active_true.matches(py, &active_false).unwrap());
            assert!(!active_true.matches(py, &other).unwrap());
        });
    }

    #[rstest]
    fn test_fluxionality_ast_roundtrip() {
        Python::attach(|py| {
            for ast in [
                GraphIrFluxionalityForm {
                    permutation: GraphIrLigandPermutation(PermPermutation::from_image(&[
                        1, 0, 2, 3,
                    ])),
                    active: GraphIrBooleanForm::Lit(false),
                },
                GraphIrFluxionalityForm {
                    permutation: GraphIrLigandPermutation(PermPermutation::identity(4)),
                    active: GraphIrBooleanForm::Undetermined,
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
                TopicityAst::new(py, pair, TopicityRelationLike::Lit(Topicity::Homotopic)).unwrap();
            assert!(value.pair() == pair);
            assert_eq!(
                value.relation.bind(py).borrow().to_rust(),
                GraphIrTopicityRelationForm::Lit(GraphIrTopicity::Homotopic)
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
            assert!(wildcard.matches(py, &homotopic).unwrap());
            assert!(!homotopic.matches(py, &enantiotopic).unwrap());
            assert!(!homotopic.matches(py, &other).unwrap());
        });
    }

    #[rstest]
    fn test_topicity_ast_roundtrip() {
        Python::attach(|py| {
            for ast in [
                GraphIrTopicityForm {
                    pair: GraphIrStereoLigandPair::new(
                        GraphIrStereoLigandPosition(0),
                        GraphIrStereoLigandPosition(2),
                    ),
                    relation: GraphIrTopicityRelationForm::Lit(GraphIrTopicity::Homotopic),
                },
                GraphIrTopicityForm {
                    pair: GraphIrStereoLigandPair::new(
                        GraphIrStereoLigandPosition(1),
                        GraphIrStereoLigandPosition(3),
                    ),
                    relation: GraphIrTopicityRelationForm::Undetermined,
                },
            ] {
                assert_eq!(TopicityAst::from_rust(py, &ast).unwrap().to_rust(py), ast);
            }
        });
    }

    #[rustfmt::skip]
    #[rstest]
    #[case(GraphIrStereoAtomConstraintForm::LigandSymmetry(GraphIrLigandSymmetryForm { permutation: GraphIrOrientedLigandPermutation { permutation: GraphIrLigandPermutation(PermPermutation::from_image(&[1, 0, 2, 3])), orientation: PermOrientation::Proper }, invariant: GraphIrBooleanForm::Lit(true) }))]
    #[case(GraphIrStereoAtomConstraintForm::Fluxionality(GraphIrFluxionalityForm { permutation: GraphIrLigandPermutation(PermPermutation::identity(4)), active: GraphIrBooleanForm::Lit(false) }))]
    #[case(GraphIrStereoAtomConstraintForm::Topicity(GraphIrTopicityForm { pair: GraphIrStereoLigandPair::new(GraphIrStereoLigandPosition(0), GraphIrStereoLigandPosition(1)), relation: GraphIrTopicityRelationForm::Lit(GraphIrTopicity::Homotopic) }))]
    #[case(GraphIrStereoAtomConstraintForm::Stereogenicity(GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic)))]
    fn test_stereo_atom_constraint_ast_roundtrip(#[case] ast: GraphIrStereoAtomConstraintForm) {
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
            let ast = GraphIrStereoAtomConstraintForm::Topicity(GraphIrTopicityForm {
                pair: GraphIrStereoLigandPair::new(
                    GraphIrStereoLigandPosition(0),
                    GraphIrStereoLigandPosition(1),
                ),
                relation: GraphIrTopicityRelationForm::Lit(GraphIrTopicity::Homotopic),
            });
            let key = StereoAtomConstraintAst::from_rust(py, &ast)
                .unwrap()
                .key(py)
                .unwrap();
            assert_eq!(
                key.to_rust(py),
                GraphIrStereoAtomConstraintKey::Topicity(GraphIrStereoLigandPair::new(
                    GraphIrStereoLigandPosition(0),
                    GraphIrStereoLigandPosition(1),
                ))
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_constraints_ast_get() {
        Python::attach(|py| {
            let stereogenicity = GraphIrStereoAtomConstraintForm::Stereogenicity(
                GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
            );
            let mut ast_cs = GraphIrStereoAtomConstraintsForm::new();
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
                    &GraphIrStereoAtomConstraintForm::Stereogenicity(
                        GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
                    ),
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
                GraphIrStereoAtomConstraintForm::Stereogenicity(GraphIrStereogenicityForm::Lit(
                    GraphIrStereogenicity::Stereogenic
                ))
            );
            assert_eq!(constraints.__len__(), 0);
            assert!(constraints.pop(py, key).unwrap().is_none());
        });
    }

    #[rstest]
    fn test_stereo_atom_constraints_ast_accessors() {
        Python::attach(|py| {
            let mut ast_cs = GraphIrStereoAtomConstraintsForm::new();
            ast_cs.extend([
                GraphIrStereoAtomConstraintForm::LigandSymmetry(GraphIrLigandSymmetryForm {
                    permutation: GraphIrOrientedLigandPermutation {
                        permutation: GraphIrLigandPermutation(PermPermutation::from_image(&[
                            1, 0, 2, 3,
                        ])),
                        orientation: PermOrientation::Proper,
                    },
                    invariant: GraphIrBooleanForm::Lit(true),
                }),
                GraphIrStereoAtomConstraintForm::Topicity(GraphIrTopicityForm {
                    pair: GraphIrStereoLigandPair::new(
                        GraphIrStereoLigandPosition(0),
                        GraphIrStereoLigandPosition(1),
                    ),
                    relation: GraphIrTopicityRelationForm::Lit(GraphIrTopicity::Homotopic),
                }),
                GraphIrStereoAtomConstraintForm::Stereogenicity(GraphIrStereogenicityForm::Lit(
                    GraphIrStereogenicity::Stereogenic,
                )),
            ]);
            let constraints = StereoAtomConstraintsAst::from_inner(ast_cs);

            assert_eq!(
                constraints.stereogenicity().to_rust(),
                GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic)
            );
            assert_eq!(
                constraints.topicity(StereoLigandPair::new(0, 1)).to_rust(),
                GraphIrTopicityRelationForm::Lit(GraphIrTopicity::Homotopic)
            );
            let ligand_symmetries = constraints.ligand_symmetries(py).unwrap();
            assert_eq!(ligand_symmetries.len(), 1);
            assert_eq!(
                ligand_symmetries[0].to_rust(py).invariant,
                GraphIrBooleanForm::Lit(true)
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_constraints_ast_iter() {
        Python::attach(|py| {
            let mut ast_cs = GraphIrStereoAtomConstraintsForm::new();
            ast_cs.extend([
                GraphIrStereoAtomConstraintForm::Topicity(GraphIrTopicityForm {
                    pair: GraphIrStereoLigandPair::new(
                        GraphIrStereoLigandPosition(0),
                        GraphIrStereoLigandPosition(1),
                    ),
                    relation: GraphIrTopicityRelationForm::Lit(GraphIrTopicity::Homotopic),
                }),
                GraphIrStereoAtomConstraintForm::Stereogenicity(GraphIrStereogenicityForm::Lit(
                    GraphIrStereogenicity::Stereogenic,
                )),
            ]);
            let constraints = StereoAtomConstraintsAst::from_inner(ast_cs);

            let keys: Vec<GraphIrStereoAtomConstraintKey> = constraints
                .keys(py)
                .unwrap()
                .keys
                .map(|k| k.bind(py).borrow().to_rust(py))
                .collect();
            assert_eq!(
                keys,
                vec![
                    GraphIrStereoAtomConstraintKey::Topicity(GraphIrStereoLigandPair::new(
                        GraphIrStereoLigandPosition(0),
                        GraphIrStereoLigandPosition(1),
                    )),
                    GraphIrStereoAtomConstraintKey::Stereogenicity,
                ]
            );
            let values: Vec<GraphIrStereoAtomConstraintForm> = constraints
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
                    &GraphIrStereoAtomConstraintForm::Stereogenicity(
                        GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
                    ),
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
    fn test_stereo_atom_constraints_like_to_rust() {
        Python::attach(|py| {
            let entry = into_py_variant(
                py,
                StereoAtomConstraintAst::from_rust(
                    py,
                    &GraphIrStereoAtomConstraintForm::Stereogenicity(
                        GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
                    ),
                )
                .unwrap(),
            )
            .unwrap();
            let container = Py::new(py, StereoAtomConstraintsAst::new(py, vec![entry])).unwrap();
            let arg = StereoAtomConstraintsLike::Container(container);
            let mut expected = GraphIrStereoAtomConstraintsForm::new();
            expected.extend([GraphIrStereoAtomConstraintForm::Stereogenicity(
                GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
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
            let stereogenicity = GraphIrStereoBondConstraintForm::Stereogenicity(
                GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
            );
            let mut ast_cs = GraphIrStereoBondConstraintsForm::new();
            ast_cs.extend([stereogenicity.clone()]);
            let constraints = StereoBondConstraintsAst::from_inner(ast_cs);
            assert_eq!(constraints.__len__(), 1);
            assert_eq!(
                constraints.stereogenicity().to_rust(),
                GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic)
            );

            let mut container_ast = GraphIrStereoBondConstraintsForm::new();
            container_ast.extend([stereogenicity.clone()]);
            let container =
                Py::new(py, StereoBondConstraintsAst::from_inner(container_ast)).unwrap();
            let arg = StereoBondConstraintsLike::Container(container);
            let mut expected = GraphIrStereoBondConstraintsForm::new();
            expected.extend([stereogenicity]);
            assert_eq!(arg.to_rust(py).unwrap(), expected);
        });
    }

    #[rstest]
    fn test_stereo_atom_constraints_view_set() {
        Python::attach(|py| {
            let value = Py::new(
                py,
                StereoAtomAst::from_inner(GraphIrStereoAtomForm::new(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(0),
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
                    &GraphIrStereoAtomConstraintForm::Stereogenicity(
                        GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
                    ),
                )
                .unwrap(),
            )
            .unwrap();
            view.set(py, stereogenicity);
            assert_eq!(
                value.borrow(py).inner().constraints.stereogenicity(),
                GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic)
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_constraints_view_pop() {
        Python::attach(|py| {
            let mut ast_cs = GraphIrStereoAtomConstraintsForm::new();
            ast_cs.extend([GraphIrStereoAtomConstraintForm::Stereogenicity(
                GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
            )]);
            let value = Py::new(
                py,
                StereoAtomAst::from_inner(GraphIrStereoAtomForm {
                    configuration: GraphIrStereoConfigurationForm::Kinded(
                        GraphIrStereoKind::Tetrahedral,
                        GraphIrStereoCoset::Lit(0),
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
                GraphIrStereoAtomConstraintForm::Stereogenicity(GraphIrStereogenicityForm::Lit(
                    GraphIrStereogenicity::Stereogenic
                ))
            );
            assert_eq!(value.borrow(py).inner().constraints.len(), 0);
        });
    }

    #[rstest]
    fn test_stereo_atom_constraints_view_getitem() {
        Python::attach(|py| {
            let mut ast_cs = GraphIrStereoAtomConstraintsForm::new();
            ast_cs.extend([GraphIrStereoAtomConstraintForm::Stereogenicity(
                GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
            )]);
            let value = Py::new(
                py,
                StereoAtomAst::from_inner(GraphIrStereoAtomForm {
                    configuration: GraphIrStereoConfigurationForm::Kinded(
                        GraphIrStereoKind::Tetrahedral,
                        GraphIrStereoCoset::Lit(0),
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
                GraphIrStereoAtomConstraintForm::Stereogenicity(GraphIrStereogenicityForm::Lit(
                    GraphIrStereogenicity::Stereogenic
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
            let mut ast_cs = GraphIrStereoAtomConstraintsForm::new();
            ast_cs.extend([
                GraphIrStereoAtomConstraintForm::Topicity(GraphIrTopicityForm {
                    pair: GraphIrStereoLigandPair::new(
                        GraphIrStereoLigandPosition(0),
                        GraphIrStereoLigandPosition(1),
                    ),
                    relation: GraphIrTopicityRelationForm::Lit(GraphIrTopicity::Homotopic),
                }),
                GraphIrStereoAtomConstraintForm::Stereogenicity(GraphIrStereogenicityForm::Lit(
                    GraphIrStereogenicity::Stereogenic,
                )),
            ]);
            let value = Py::new(
                py,
                StereoAtomAst::from_inner(GraphIrStereoAtomForm {
                    configuration: GraphIrStereoConfigurationForm::Kinded(
                        GraphIrStereoKind::Tetrahedral,
                        GraphIrStereoCoset::Lit(0),
                    ),
                    constraints: ast_cs,
                }),
            )
            .unwrap();
            let view = StereoAtomConstraintsView {
                backing: StereoAtomConstraintsBacking::Value(value.clone_ref(py)),
            };
            let keys: Vec<GraphIrStereoAtomConstraintKey> = view
                .keys(py)
                .unwrap()
                .keys
                .map(|k| k.bind(py).borrow().to_rust(py))
                .collect();
            assert_eq!(
                keys,
                vec![
                    GraphIrStereoAtomConstraintKey::Topicity(GraphIrStereoLigandPair::new(
                        GraphIrStereoLigandPosition(0),
                        GraphIrStereoLigandPosition(1),
                    )),
                    GraphIrStereoAtomConstraintKey::Stereogenicity,
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
                StereoAtomAst::from_inner(GraphIrStereoAtomForm::new(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(0),
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
                    &GraphIrStereoAtomConstraintForm::Stereogenicity(
                        GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
                    ),
                )
                .unwrap(),
            )
            .unwrap();
            view.update(py, StereoAtomConstraintsUpdate::Entries(vec![entry]))
                .unwrap();
            assert_eq!(
                value.borrow(py).inner().constraints.stereogenicity(),
                GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic)
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_constraints_view_accessors() {
        Python::attach(|py| {
            let mut ast_cs = GraphIrStereoAtomConstraintsForm::new();
            ast_cs.extend([
                GraphIrStereoAtomConstraintForm::Topicity(GraphIrTopicityForm {
                    pair: GraphIrStereoLigandPair::new(
                        GraphIrStereoLigandPosition(0),
                        GraphIrStereoLigandPosition(1),
                    ),
                    relation: GraphIrTopicityRelationForm::Lit(GraphIrTopicity::Homotopic),
                }),
                GraphIrStereoAtomConstraintForm::Stereogenicity(GraphIrStereogenicityForm::Lit(
                    GraphIrStereogenicity::Stereogenic,
                )),
            ]);
            let value = Py::new(
                py,
                StereoAtomAst::from_inner(GraphIrStereoAtomForm {
                    configuration: GraphIrStereoConfigurationForm::Kinded(
                        GraphIrStereoKind::Tetrahedral,
                        GraphIrStereoCoset::Lit(0),
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
                GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic)
            );
            assert_eq!(
                view.topicity(py, StereoLigandPair::new(0, 1))
                    .unwrap()
                    .to_rust(),
                GraphIrTopicityRelationForm::Lit(GraphIrTopicity::Homotopic)
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_ast_constraints() {
        Python::attach(|py| {
            let value = Py::new(
                py,
                StereoAtomAst::from_inner(GraphIrStereoAtomForm::new(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(0),
                )),
            )
            .unwrap();
            let view = StereoAtomAst::constraints(value.clone_ref(py));
            let stereogenicity = into_py_variant(
                py,
                StereoAtomConstraintAst::from_rust(
                    py,
                    &GraphIrStereoAtomConstraintForm::Stereogenicity(
                        GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
                    ),
                )
                .unwrap(),
            )
            .unwrap();
            view.set(py, stereogenicity);
            assert_eq!(
                value.borrow(py).inner().constraints.stereogenicity(),
                GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic)
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_ast_set_constraints_self() {
        Python::attach(|py| {
            let mut ast_cs = GraphIrStereoAtomConstraintsForm::new();
            ast_cs.extend([GraphIrStereoAtomConstraintForm::Stereogenicity(
                GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
            )]);
            let value = Py::new(
                py,
                StereoAtomAst::from_inner(GraphIrStereoAtomForm {
                    configuration: GraphIrStereoConfigurationForm::Kinded(
                        GraphIrStereoKind::Tetrahedral,
                        GraphIrStereoCoset::Lit(0),
                    ),
                    constraints: ast_cs,
                }),
            )
            .unwrap();
            let own_view = StereoAtomAst::constraints(value.clone_ref(py));
            StereoAtomAst::set_constraints(
                value.clone_ref(py),
                py,
                StereoAtomConstraintsLike::View(Py::new(py, own_view).unwrap()),
            )
            .unwrap();
            assert_eq!(
                value.borrow(py).inner().constraints.stereogenicity(),
                GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic)
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_constraints_view_update_self() {
        Python::attach(|py| {
            let mut ast_cs = GraphIrStereoAtomConstraintsForm::new();
            ast_cs.extend([GraphIrStereoAtomConstraintForm::Stereogenicity(
                GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
            )]);
            let value = Py::new(
                py,
                StereoAtomAst::from_inner(GraphIrStereoAtomForm {
                    configuration: GraphIrStereoConfigurationForm::Kinded(
                        GraphIrStereoKind::Tetrahedral,
                        GraphIrStereoCoset::Lit(0),
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
                StereoBondAst::from_inner(GraphIrStereoBondForm::new(
                    GraphIrStereoKind::CisTrans,
                    GraphIrStereoCoset::Lit(0),
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
                    &GraphIrStereoBondConstraintForm::Stereogenicity(
                        GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
                    ),
                )
                .unwrap(),
            )
            .unwrap();
            view.set(py, stereogenicity);
            assert_eq!(
                value.borrow(py).inner().constraints.stereogenicity(),
                GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic)
            );
        });
    }

    #[rstest]
    #[case::ccw(
        StereoConfigurationLike::Tetrahedral(TetrahedralConfiguration::Ccw),
        GraphIrStereoAtomForm::new(GraphIrStereoKind::Tetrahedral, GraphIrStereoCoset::Lit(0))
    )]
    #[case::cw(
        StereoConfigurationLike::Tetrahedral(TetrahedralConfiguration::Cw),
        GraphIrStereoAtomForm::new(GraphIrStereoKind::Tetrahedral, GraphIrStereoCoset::Lit(1))
    )]
    fn test_stereo_atom_ast_new(
        #[case] configuration: StereoConfigurationLike,
        #[case] expected: GraphIrStereoAtomForm,
    ) {
        Python::attach(|py| {
            let value = StereoAtomAst::new(py, configuration, None);
            assert_eq!(*value.inner(), expected);
        });
    }

    #[rstest]
    fn test_stereo_atom_ast_new_constraints() {
        Python::attach(|py| {
            let stereogenicity = GraphIrStereoAtomConstraintForm::Stereogenicity(
                GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
            );
            let mut ast_cs = GraphIrStereoAtomConstraintsForm::new();
            ast_cs.extend([stereogenicity.clone()]);
            let container = Py::new(py, StereoAtomConstraintsAst::from_inner(ast_cs)).unwrap();
            let value = StereoAtomAst::new(
                py,
                StereoConfigurationLike::Tetrahedral(TetrahedralConfiguration::Ccw),
                Some(container),
            );
            let mut expected_cs = GraphIrStereoAtomConstraintsForm::new();
            expected_cs.extend([stereogenicity]);
            assert_eq!(
                *value.inner(),
                GraphIrStereoAtomForm {
                    configuration: GraphIrStereoConfigurationForm::Kinded(
                        GraphIrStereoKind::Tetrahedral,
                        GraphIrStereoCoset::Lit(0)
                    ),
                    constraints: expected_cs,
                }
            );
        });
    }

    #[rstest]
    #[case::ccw(
        "Th0",
        GraphIrStereoConfigurationForm::Kinded(
            GraphIrStereoKind::Tetrahedral,
            GraphIrStereoCoset::Lit(0)
        )
    )]
    #[case::undetermined_coset(
        "Th*",
        GraphIrStereoConfigurationForm::Kinded(
            GraphIrStereoKind::Tetrahedral,
            GraphIrStereoCoset::Undetermined
        )
    )]
    fn test_stereo_atom_ast_parse(
        #[case] input: &str,
        #[case] expected: GraphIrStereoConfigurationForm,
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
        GraphIrStereoAtomForm::new(GraphIrStereoKind::Tetrahedral, GraphIrStereoCoset::Lit(0)),
        "Th0"
    )]
    #[case::square_planar(
        GraphIrStereoAtomForm::new(GraphIrStereoKind::SquarePlanar, GraphIrStereoCoset::Lit(2)),
        "Sp2"
    )]
    fn test_stereo_atom_ast_str(#[case] ast: GraphIrStereoAtomForm, #[case] expected: &str) {
        let value = StereoAtomAst::from_inner(ast);
        assert_eq!(value.__str__(), expected);
    }

    #[rstest]
    fn test_stereo_atom_ast_repr() {
        let value = StereoAtomAst::from_inner(GraphIrStereoAtomForm::new(
            GraphIrStereoKind::Tetrahedral,
            GraphIrStereoCoset::Lit(0),
        ));
        assert_eq!(value.__repr__(), "StereoAtomAst.parse('Th0')");
    }

    #[rstest]
    fn test_stereo_atom_ast_configuration() {
        Python::attach(|py| {
            let value = StereoAtomAst::from_inner(GraphIrStereoAtomForm::new(
                GraphIrStereoKind::Tetrahedral,
                GraphIrStereoCoset::Lit(0),
            ));
            assert_eq!(
                value.configuration(py).unwrap().to_rust(py),
                GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(0)
                )
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_ast_set_configuration() {
        Python::attach(|py| {
            let mut value = StereoAtomAst::from_inner(GraphIrStereoAtomForm::new(
                GraphIrStereoKind::Tetrahedral,
                GraphIrStereoCoset::Lit(0),
            ));
            value.set_configuration(
                py,
                StereoConfigurationLike::Tetrahedral(TetrahedralConfiguration::Cw),
            );
            assert_eq!(
                value.inner().configuration,
                GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(1)
                )
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_ast_set_constraints() {
        Python::attach(|py| {
            let value = Py::new(
                py,
                StereoAtomAst::from_inner(GraphIrStereoAtomForm::new(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(0),
                )),
            )
            .unwrap();
            let stereogenicity = GraphIrStereoAtomConstraintForm::Stereogenicity(
                GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
            );
            let mut ast_cs = GraphIrStereoAtomConstraintsForm::new();
            ast_cs.extend([stereogenicity.clone()]);
            let container = Py::new(py, StereoAtomConstraintsAst::from_inner(ast_cs)).unwrap();
            StereoAtomAst::set_constraints(
                value.clone_ref(py),
                py,
                StereoAtomConstraintsLike::Container(container),
            )
            .unwrap();
            let mut expected_cs = GraphIrStereoAtomConstraintsForm::new();
            expected_cs.extend([stereogenicity]);
            assert_eq!(value.borrow(py).inner().constraints, expected_cs);
        });
    }

    #[rstest]
    fn test_stereo_atom_ast_asdict() {
        Python::attach(|py| {
            let stereogenicity = GraphIrStereoAtomConstraintForm::Stereogenicity(
                GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
            );
            let mut ast_cs = GraphIrStereoAtomConstraintsForm::new();
            ast_cs.extend([stereogenicity]);
            let value = StereoAtomAst::from_inner(GraphIrStereoAtomForm {
                configuration: GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(0),
                ),
                constraints: ast_cs,
            });
            let dict = value.asdict(py).unwrap();
            let configuration = dict.get_item("configuration").unwrap().unwrap();
            let expected = into_py_variant(
                py,
                StereoConfigurationAst::from_rust(
                    py,
                    &GraphIrStereoConfigurationForm::Kinded(
                        GraphIrStereoKind::Tetrahedral,
                        GraphIrStereoCoset::Lit(0),
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
                StereoConfigurationLike::CisTrans(CisTransConfiguration::Z),
                None,
            );
            assert_eq!(
                *value.inner(),
                GraphIrStereoBondForm::new(GraphIrStereoKind::CisTrans, GraphIrStereoCoset::Lit(0))
            );
            assert_eq!(value.__str__(), "Ct0");
        });
    }

    #[rstest]
    #[case::z(
        "Ct0",
        GraphIrStereoBondForm::new(GraphIrStereoKind::CisTrans, GraphIrStereoCoset::Lit(0))
    )]
    #[case::e(
        "Ct1",
        GraphIrStereoBondForm::new(GraphIrStereoKind::CisTrans, GraphIrStereoCoset::Lit(1))
    )]
    fn test_stereo_bond_ast_parse(#[case] input: &str, #[case] expected: GraphIrStereoBondForm) {
        let value = StereoBondAst::parse(input).unwrap();
        assert_eq!(*value.inner(), expected);
    }

    #[rstest]
    fn test_stereo_bond_ast_str() {
        let value = StereoBondAst::from_inner(GraphIrStereoBondForm::new(
            GraphIrStereoKind::CisTrans,
            GraphIrStereoCoset::Lit(1),
        ));
        assert_eq!(value.__str__(), "Ct1");
    }

    fn stereo_atom_molecule(py: Python<'_>) -> Py<MoleculeAst> {
        let molecule = GraphIrMolecule::from_entries(MoleculeEntries {
            atoms: vec![GraphIrAtomForm::from_element(ChemElement::C); 5],
            stereo_atoms: vec![(
                GraphIrAtomId(0),
                vec![
                    GraphIrStereoLigand::new(GraphIrAtomId(1), GraphIrStereoLigandKind::Atom),
                    GraphIrStereoLigand::new(GraphIrAtomId(2), GraphIrStereoLigandKind::Atom),
                    GraphIrStereoLigand::new(GraphIrAtomId(3), GraphIrStereoLigandKind::Atom),
                    GraphIrStereoLigand::new(GraphIrAtomId(4), GraphIrStereoLigandKind::Atom),
                ],
                GraphIrStereoAtomForm::new(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(0),
                ),
            )],
            ..Default::default()
        });
        Py::new(py, MoleculeAst::from_rust(molecule)).unwrap()
    }

    fn stereo_bond_molecule(py: Python<'_>) -> Py<MoleculeAst> {
        let molecule = GraphIrMolecule::from_entries(MoleculeEntries {
            atoms: vec![GraphIrAtomForm::from_element(ChemElement::C); 4],
            bonds: vec![
                (
                    GraphIrAtomId(0),
                    GraphIrAtomId(1),
                    GraphIrBondForm::from_order(2),
                ),
                (
                    GraphIrAtomId(0),
                    GraphIrAtomId(2),
                    GraphIrBondForm::from_order(1),
                ),
                (
                    GraphIrAtomId(1),
                    GraphIrAtomId(3),
                    GraphIrBondForm::from_order(1),
                ),
            ],
            stereo_bonds: vec![(
                GraphIrBondId(0),
                vec![
                    GraphIrStereoLigand::new(GraphIrAtomId(2), GraphIrStereoLigandKind::Atom),
                    GraphIrStereoLigand::new(GraphIrAtomId(3), GraphIrStereoLigandKind::Atom),
                ],
                GraphIrStereoBondForm::new(GraphIrStereoKind::CisTrans, GraphIrStereoCoset::Lit(0)),
            )],
            ..Default::default()
        });
        Py::new(py, MoleculeAst::from_rust(molecule)).unwrap()
    }

    #[rstest]
    fn test_stereo_atom_view_id() {
        Python::attach(|py| {
            let view = StereoAtomView {
                owner: stereo_atom_molecule(py),
                id: GraphIrStereoAtomId(0),
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
                id: GraphIrStereoAtomId(0),
            };
            assert_eq!(view.site_id(py).unwrap(), 0);
        });
    }

    #[rstest]
    fn test_stereo_atom_view_ligands() {
        Python::attach(|py| {
            let view = StereoAtomView {
                owner: stereo_atom_molecule(py),
                id: GraphIrStereoAtomId(0),
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
                id: GraphIrStereoAtomId(0),
            };
            assert_eq!(view.kind(py).unwrap(), StereoKind::Tetrahedral);
        });
    }

    #[rstest]
    fn test_stereo_atom_view_coset() {
        Python::attach(|py| {
            let view = StereoAtomView {
                owner: stereo_atom_molecule(py),
                id: GraphIrStereoAtomId(0),
            };
            assert_eq!(
                view.coset(py).unwrap().to_rust(py),
                GraphIrStereoCoset::Lit(0)
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_view_configuration() {
        Python::attach(|py| {
            let view = StereoAtomView {
                owner: stereo_atom_molecule(py),
                id: GraphIrStereoAtomId(0),
            };
            assert_eq!(
                view.configuration(py).unwrap().to_rust(py),
                GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(0)
                )
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_view_set_configuration() {
        Python::attach(|py| {
            let view = StereoAtomView {
                owner: stereo_atom_molecule(py),
                id: GraphIrStereoAtomId(0),
            };
            view.set_configuration(
                py,
                StereoConfigurationLike::Tetrahedral(TetrahedralConfiguration::Cw),
            );
            assert_eq!(
                view.configuration(py).unwrap().to_rust(py),
                GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(1)
                )
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_view_constraints() {
        Python::attach(|py| {
            let view = StereoAtomView {
                owner: stereo_atom_molecule(py),
                id: GraphIrStereoAtomId(0),
            };
            let stereogenicity = into_py_variant(
                py,
                StereoAtomConstraintAst::from_rust(
                    py,
                    &GraphIrStereoAtomConstraintForm::Stereogenicity(
                        GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
                    ),
                )
                .unwrap(),
            )
            .unwrap();
            view.constraints(py).set(py, stereogenicity);
            // a fresh molecule-backed handle proves the write hit the molecule
            assert_eq!(
                view.constraints(py).stereogenicity(py).unwrap().to_rust(),
                GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic)
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_view_set_constraints() {
        Python::attach(|py| {
            let view = StereoAtomView {
                owner: stereo_atom_molecule(py),
                id: GraphIrStereoAtomId(0),
            };
            let mut ast_cs = GraphIrStereoAtomConstraintsForm::new();
            ast_cs.extend([GraphIrStereoAtomConstraintForm::Stereogenicity(
                GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic),
            )]);
            let container = Py::new(py, StereoAtomConstraintsAst::from_inner(ast_cs)).unwrap();
            view.set_constraints(py, StereoAtomConstraintsLike::Container(container))
                .unwrap();
            assert_eq!(
                view.constraints(py).stereogenicity(py).unwrap().to_rust(),
                GraphIrStereogenicityForm::Lit(GraphIrStereogenicity::Stereogenic)
            );
        });
    }

    #[rstest]
    fn test_stereo_atom_view_asdict() {
        Python::attach(|py| {
            let view = StereoAtomView {
                owner: stereo_atom_molecule(py),
                id: GraphIrStereoAtomId(0),
            };
            let dict = view.asdict(py).unwrap();
            let configuration = dict.get_item("configuration").unwrap().unwrap();
            let expected = into_py_variant(
                py,
                StereoConfigurationAst::from_rust(
                    py,
                    &GraphIrStereoConfigurationForm::Kinded(
                        GraphIrStereoKind::Tetrahedral,
                        GraphIrStereoCoset::Lit(0),
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
                id: GraphIrStereoAtomId(5),
            };
            assert!(view.site_id(py).is_err());
        });
    }

    #[rstest]
    fn test_stereo_bond_view() {
        Python::attach(|py| {
            let view = StereoBondView {
                owner: stereo_bond_molecule(py),
                id: GraphIrStereoBondId(0),
            };
            assert_eq!(view.id(), 0);
            assert_eq!(view.site_id(py).unwrap(), 0);
            assert_eq!(
                view.configuration(py).unwrap().to_rust(py),
                GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::CisTrans,
                    GraphIrStereoCoset::Lit(0)
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
                StereoAtomAst::from_inner(GraphIrStereoAtomForm::new(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(1),
                )),
            )
            .unwrap();
            views.__setitem__(py, 0, replacement.borrow(py)).unwrap();
            let view = views.__getitem__(py, 0).unwrap();
            // value replaced
            assert_eq!(
                view.configuration(py).unwrap().to_rust(py),
                GraphIrStereoConfigurationForm::Kinded(
                    GraphIrStereoKind::Tetrahedral,
                    GraphIrStereoCoset::Lit(1)
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
