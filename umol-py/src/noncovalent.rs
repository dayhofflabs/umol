//! Noncovalent-bond bindings. This file grows over the B5 slice; it opens with the
//! kind leaf: `NoncovalentBondKind` (the interaction kind) and `NoncovalentBondKindAst`
//! (`Undetermined | Lit(kind)`), mirroring `umol_ast::ast::{NoncovalentBondKind,
//! NoncovalentBondKindAst}` — the noncovalent analog of `atom.element: ElementAst`
//! over the `Element` value enum.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use pyo3::prelude::*;
use umol_ast::ast::{
    AsLit, NoncovalentBondKind as AstNoncovalentBondKind,
    NoncovalentBondKindAst as AstNoncovalentBondKindAst,
};

use crate::convert::{hash_ast, variant_repr};

/// A noncovalent interaction kind. A fieldless, hashable value enum whose members
/// mirror the Rust `NoncovalentBondKind` exactly.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum NoncovalentBondKind {
    HydrogenBond,
    HalogenBond,
    ChalcogenBond,
    Ionic,
    VanDerWaals,
}

impl NoncovalentBondKind {
    pub(crate) fn from_ast(ast: AstNoncovalentBondKind) -> Self {
        match ast {
            AstNoncovalentBondKind::HydrogenBond => Self::HydrogenBond,
            AstNoncovalentBondKind::HalogenBond => Self::HalogenBond,
            AstNoncovalentBondKind::ChalcogenBond => Self::ChalcogenBond,
            AstNoncovalentBondKind::Ionic => Self::Ionic,
            AstNoncovalentBondKind::VanDerWaals => Self::VanDerWaals,
        }
    }

    pub(crate) fn to_ast(self) -> AstNoncovalentBondKind {
        match self {
            Self::HydrogenBond => AstNoncovalentBondKind::HydrogenBond,
            Self::HalogenBond => AstNoncovalentBondKind::HalogenBond,
            Self::ChalcogenBond => AstNoncovalentBondKind::ChalcogenBond,
            Self::Ionic => AstNoncovalentBondKind::Ionic,
            Self::VanDerWaals => AstNoncovalentBondKind::VanDerWaals,
        }
    }
}

/// A noncovalent bond's interaction kind: undetermined, or a concrete
/// `NoncovalentBondKind`. Mirrors `NoncovalentBondKindAst`.
#[pyclass]
pub enum NoncovalentBondKindAst {
    Undetermined(),
    Lit(NoncovalentBondKind),
}

#[pymethods]
impl NoncovalentBondKindAst {
    /// The concrete interaction kind, or `None` when undetermined.
    fn as_lit(&self) -> Option<NoncovalentBondKind> {
        self.to_ast().as_lit().map(NoncovalentBondKind::from_ast)
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.to_ast() == other.to_ast()
    }

    fn __hash__(&self) -> u64 {
        hash_ast(&self.to_ast())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            NoncovalentBondKindAst::Undetermined() => ("Undetermined", 0),
            NoncovalentBondKindAst::Lit(_) => ("Lit", 1),
        };
        variant_repr(
            slf.bind(py).as_any(),
            "NoncovalentBondKindAst",
            variant,
            arity,
        )
    }
}

impl NoncovalentBondKindAst {
    #[cfg(test)]
    pub(crate) fn from_ast(ast: &AstNoncovalentBondKindAst) -> Self {
        match ast {
            AstNoncovalentBondKindAst::Undetermined => Self::Undetermined(),
            AstNoncovalentBondKindAst::Lit(k) => Self::Lit(NoncovalentBondKind::from_ast(*k)),
        }
    }

    pub(crate) fn to_ast(&self) -> AstNoncovalentBondKindAst {
        match self {
            Self::Undetermined() => AstNoncovalentBondKindAst::Undetermined,
            Self::Lit(k) => AstNoncovalentBondKindAst::Lit(k.to_ast()),
        }
    }
}

/// Setter coercion for a noncovalent `kind` field: a bare `NoncovalentBondKind` →
/// `Lit`, or a `NoncovalentBondKindAst` passthrough (mirroring the `Undetermined |
/// Lit` structure).
#[cfg(test)]
#[derive(FromPyObject)]
pub(crate) enum NoncovalentBondKindArg {
    Kind(NoncovalentBondKind),
    Ast(Py<NoncovalentBondKindAst>),
}

#[cfg(test)]
impl NoncovalentBondKindArg {
    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstNoncovalentBondKindAst {
        match self {
            NoncovalentBondKindArg::Kind(k) => AstNoncovalentBondKindAst::Lit(k.to_ast()),
            NoncovalentBondKindArg::Ast(a) => a.bind(py).borrow().to_ast(),
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case(AstNoncovalentBondKindAst::Undetermined)]
    #[case(AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::HydrogenBond))]
    #[case(AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::VanDerWaals))]
    fn test_noncovalent_bond_kind_ast_roundtrip(#[case] ast: AstNoncovalentBondKindAst) {
        assert_eq!(NoncovalentBondKindAst::from_ast(&ast).to_ast(), ast);
    }

    #[rstest]
    #[case(AstNoncovalentBondKindAst::Undetermined, None)]
    #[case(
        AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::Ionic),
        Some(NoncovalentBondKind::Ionic)
    )]
    fn test_noncovalent_bond_kind_ast_as_lit(
        #[case] ast: AstNoncovalentBondKindAst,
        #[case] expected: Option<NoncovalentBondKind>,
    ) {
        assert_eq!(NoncovalentBondKindAst::from_ast(&ast).as_lit(), expected);
    }

    #[rstest]
    #[case(AstNoncovalentBondKind::HydrogenBond)]
    #[case(AstNoncovalentBondKind::ChalcogenBond)]
    fn test_noncovalent_bond_kind_roundtrip(#[case] ast: AstNoncovalentBondKind) {
        assert_eq!(NoncovalentBondKind::from_ast(ast).to_ast(), ast);
    }

    #[rstest]
    fn test_noncovalent_bond_kind_arg_to_ast() {
        Python::attach(|py| {
            // a bare kind coerces to Lit
            assert_eq!(
                NoncovalentBondKindArg::Kind(NoncovalentBondKind::HydrogenBond).to_ast(py),
                AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::HydrogenBond)
            );
            // a NoncovalentBondKindAst passes through
            let ast = Py::new(py, NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic)).unwrap();
            assert_eq!(
                NoncovalentBondKindArg::Ast(ast).to_ast(py),
                AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::Ionic)
            );
        });
    }
}
