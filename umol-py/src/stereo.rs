//! Stereo sub-ASTs mirroring `umol_ast::ast::stereo` and `umol_perm` (S5a): the
//! `Permutation` value, the recursive `StereoTerm` transformation algebra, the
//! `StereoCosetAst` coset, and the `TetrahedralStereoAst` atom-stereo state.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use std::collections::BTreeSet;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
// The stereo value-leaf AST mirrors are consumed only by the `#[cfg(test)]` `from_ast`/`to_ast`
// until S0c/S1b/S4a wire them into the config/constraint leaves and the entity view.
#[cfg(test)]
use umol_ast::ast::{
    AtomId as AstAtomId, StereoKind as AstStereoKind, StereoLigand as AstStereoLigand,
    StereoLigandKind as AstStereoLigandKind, Stereogenicity as AstStereogenicity,
    Topicity as AstTopicity,
};
use umol_ast::ast::{
    CisTransStereoAst as AstCisTransStereoAst, StereoCosetAst as AstStereoCosetAst,
    StereoTerm as AstStereoTerm, TetrahedralStereoAst as AstTetrahedralStereoAst,
};
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

    #[cfg(test)]
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
/// A fieldless, hashable value enum mirroring the Rust `Topicity`.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
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

    #[cfg(test)]
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
}
