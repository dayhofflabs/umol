//! Stereo sub-ASTs mirroring `umol_ast::ast::stereo` and `umol_perm` (S5a): the
//! `Permutation` value, the recursive `StereoTerm` transformation algebra, the
//! `StereoCosetAst` coset, and the `TetrahedralStereoAst` atom-stereo state.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use std::collections::BTreeSet;

use pyo3::prelude::*;
use umol_ast::ast::{
    StereoCosetAst as AstStereoCosetAst, StereoTerm as AstStereoTerm,
    TetrahedralStereoAst as AstTetrahedralStereoAst,
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
        variant_repr(slf.bind(py).as_any(), "TetrahedralStereoAst", variant, arity)
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
#[pyclass(eq, from_py_object)]
#[derive(Clone, Copy, PartialEq)]
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
}
