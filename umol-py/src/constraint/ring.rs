//! Ring scope and membership constraint payloads used by atom and bond families.

use pyo3::prelude::*;
use umol_ast::ast::{RingMembershipAst as AstRingMembershipAst, RingScope as AstRingScope};

use crate::convert::{hash_ast, into_py_variant, variant_repr};
use crate::value::{ValueArg, ValueAst};

#[pyclass]
pub enum RingScope {
    All(),
    Size(u8),
}

#[pymethods]
impl RingScope {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_ast() == other.to_ast()
    }

    fn __hash__(&self) -> u64 {
        hash_ast(&self.to_ast())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            Self::All() => ("All", 0),
            Self::Size(_) => ("Size", 1),
        };
        variant_repr(slf.bind(py).as_any(), "RingScope", variant, arity)
    }
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
        Ok(Self {
            scope,
            count: count.to_py(py)?,
        })
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "RingMembershipAst({}, {})",
            self.scope.bind(py).as_any().repr()?.extract::<String>()?,
            self.count.bind(py).as_any().repr()?.extract::<String>()?,
        ))
    }
}

impl RingMembershipAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstRingMembershipAst) -> PyResult<Self> {
        Ok(Self {
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

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

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
}
