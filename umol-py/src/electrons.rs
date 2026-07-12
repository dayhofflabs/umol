//! Electron-counts AST mirror: `Undetermined`, or a positional per-member-atom
//! count vector. A value leaf shared by the aromatic-system and multicenter-bond
//! bindings; the vector is positional (cell = member atom), aligned to the owning
//! entity's participant order.

use pyo3::prelude::*;
use umol_ast::ast::{AsLit, ElectronCountsAst as AstElectronCountsAst};

use crate::convert::{hash_ast, variant_repr};

/// A per-member-atom electron-count vector: undetermined, or a concrete list of
/// counts positionally aligned to the owning entity's atoms.
#[pyclass]
pub enum ElectronCountsAst {
    Undetermined(),
    Lit(Vec<i64>),
}

#[pymethods]
impl ElectronCountsAst {
    /// The concrete count vector, or `None` when undetermined.
    fn as_lit(&self) -> Option<Vec<i64>> {
        self.to_ast().as_lit()
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.to_ast() == other.to_ast()
    }

    fn __hash__(&self) -> u64 {
        hash_ast(&self.to_ast())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            ElectronCountsAst::Undetermined() => ("Undetermined", 0),
            ElectronCountsAst::Lit(_) => ("Lit", 1),
        };
        variant_repr(slf.bind(py).as_any(), "ElectronCountsAst", variant, arity)
    }
}

impl ElectronCountsAst {
    pub(crate) fn from_ast(ast: &AstElectronCountsAst) -> Self {
        match ast {
            AstElectronCountsAst::Undetermined => Self::Undetermined(),
            AstElectronCountsAst::Lit(counts) => Self::Lit(counts.clone()),
        }
    }

    pub(crate) fn to_ast(&self) -> AstElectronCountsAst {
        match self {
            Self::Undetermined() => AstElectronCountsAst::Undetermined,
            Self::Lit(counts) => AstElectronCountsAst::Lit(counts.clone()),
        }
    }
}

/// Setter coercion for an electron-counts field: a Python `list[int]` → `Lit`, or an
/// `ElectronCountsAst` passthrough (mirroring `impl From<Vec<i64>>`).
#[derive(FromPyObject)]
pub(crate) enum ElectronCountsArg {
    Lit(Vec<i64>),
    Ast(Py<ElectronCountsAst>),
}

impl ElectronCountsArg {
    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstElectronCountsAst {
        match self {
            ElectronCountsArg::Lit(counts) => AstElectronCountsAst::Lit(counts.clone()),
            ElectronCountsArg::Ast(a) => a.bind(py).borrow().to_ast(),
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case(AstElectronCountsAst::Undetermined)]
    #[case(AstElectronCountsAst::Lit(vec![1, 1, 1, 1, 1, 1]))]
    #[case(AstElectronCountsAst::Lit(vec![]))]
    fn test_electron_counts_ast_roundtrip(#[case] ast: AstElectronCountsAst) {
        assert_eq!(ElectronCountsAst::from_ast(&ast).to_ast(), ast);
    }

    #[rstest]
    #[case(AstElectronCountsAst::Undetermined, None)]
    #[case(AstElectronCountsAst::Lit(vec![2, 0, 2]), Some(vec![2, 0, 2]))]
    fn test_electron_counts_ast_as_lit(
        #[case] ast: AstElectronCountsAst,
        #[case] expected: Option<Vec<i64>>,
    ) {
        assert_eq!(ElectronCountsAst::from_ast(&ast).as_lit(), expected);
    }

    #[rstest]
    fn test_electron_counts_arg_to_ast() {
        Python::attach(|py| {
            // a bare list coerces to Lit
            assert_eq!(
                ElectronCountsArg::Lit(vec![1, 0, 1]).to_ast(py),
                AstElectronCountsAst::Lit(vec![1, 0, 1])
            );
            // an ElectronCountsAst passes through
            let ast = Py::new(py, ElectronCountsAst::Lit(vec![2, 2])).unwrap();
            assert_eq!(
                ElectronCountsArg::Ast(ast).to_ast(py),
                AstElectronCountsAst::Lit(vec![2, 2])
            );
        });
    }
}
