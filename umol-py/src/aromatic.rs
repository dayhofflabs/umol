//! Aromatic system value type and aromatic-constraint surface mirroring
//! `umol_ast::ast`: `AromaticSystemAst`, the `AromaticSystemConstraintAst` enum, the
//! `AromaticSystemConstraintsAst` container, and the `AromaticSystemConstraintsView`
//! live handle. An aromatic system is a single unordered set of member atoms; the
//! value carries a positional per-atom `electrons` vector plus charge, spin, and
//! constraints. The member atoms are the participants of the owning molecule's
//! aromatic relation, so they are topology (the view half) rather than value.

use pyo3::prelude::*;
use umol_ast::ast::{
    AromaticSystemConstraintAst as AstAromaticSystemConstraintAst,
    AromaticSystemConstraintKey as AstAromaticSystemConstraintKey,
};

#[cfg(test)]
use crate::convert::into_py_variant;
use crate::convert::{hash_ast, variant_repr};
use crate::value::ValueAst;

/// The key (identity) of an aromatic-system constraint, for keyed lookup. The
/// single key `ElectronCount` is the bare discriminant (no sub-key).
#[pyclass]
pub enum AromaticSystemConstraintKey {
    ElectronCount(),
}

#[pymethods]
impl AromaticSystemConstraintKey {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_ast() == other.to_ast()
    }

    fn __hash__(&self) -> u64 {
        hash_ast(&self.to_ast())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            AromaticSystemConstraintKey::ElectronCount() => ("ElectronCount", 0),
        };
        variant_repr(
            slf.bind(py).as_any(),
            "AromaticSystemConstraintKey",
            variant,
            arity,
        )
    }
}

impl AromaticSystemConstraintKey {
    pub(crate) fn from_ast(ast: &AstAromaticSystemConstraintKey) -> Self {
        match ast {
            AstAromaticSystemConstraintKey::ElectronCount => Self::ElectronCount(),
        }
    }

    pub(crate) fn to_ast(&self) -> AstAromaticSystemConstraintKey {
        match self {
            Self::ElectronCount() => AstAromaticSystemConstraintKey::ElectronCount,
        }
    }
}

/// An aromatic-system-scope constraint: the asserted total π-electron count of the
/// system (cross-checked against `sum(AromaticSystemAst::electrons)`).
#[pyclass]
pub enum AromaticSystemConstraintAst {
    ElectronCount(Py<ValueAst>),
}

#[pymethods]
impl AromaticSystemConstraintAst {
    /// The constraint's key (identity).
    #[getter]
    fn key(&self, py: Python<'_>) -> AromaticSystemConstraintKey {
        AromaticSystemConstraintKey::from_ast(&self.to_ast(py).key())
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let variant = match &*slf.bind(py).borrow() {
            AromaticSystemConstraintAst::ElectronCount(_) => "ElectronCount",
        };
        variant_repr(
            slf.bind(py).as_any(),
            "AromaticSystemConstraintAst",
            variant,
            1,
        )
    }
}

impl AromaticSystemConstraintAst {
    /// The AST→Python bridge, paired with `to_ast`. Test-only until the constraints
    /// container (S1b) consumes it.
    #[cfg(test)]
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstAromaticSystemConstraintAst) -> PyResult<Self> {
        Ok(match ast {
            AstAromaticSystemConstraintAst::ElectronCount(v) => {
                Self::ElectronCount(into_py_variant(py, ValueAst::from_ast(py, v)?)?)
            }
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstAromaticSystemConstraintAst {
        match self {
            Self::ElectronCount(v) => {
                AstAromaticSystemConstraintAst::ElectronCount(v.bind(py).borrow().to_ast(py))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::ValueAst as AstValueAst;

    use super::*;

    #[rstest]
    fn test_aromatic_system_constraint_key_roundtrip() {
        let key =
            AromaticSystemConstraintKey::from_ast(&AstAromaticSystemConstraintKey::ElectronCount);
        assert_eq!(key.to_ast(), AstAromaticSystemConstraintKey::ElectronCount);
    }

    #[rstest]
    fn test_aromatic_system_constraint_ast_key() {
        Python::attach(|py| {
            let constraint = AstAromaticSystemConstraintAst::electron_count(6);
            let key = AromaticSystemConstraintAst::from_ast(py, &constraint)
                .unwrap()
                .key(py);
            assert_eq!(key.to_ast(), AstAromaticSystemConstraintKey::ElectronCount);
        });
    }

    #[rstest]
    #[case(AstAromaticSystemConstraintAst::electron_count(6))]
    #[case(AstAromaticSystemConstraintAst::electron_count(AstValueAst::Undetermined))]
    fn test_aromatic_system_constraint_ast_roundtrip(#[case] ast: AstAromaticSystemConstraintAst) {
        Python::attach(|py| {
            assert_eq!(
                AromaticSystemConstraintAst::from_ast(py, &ast)
                    .unwrap()
                    .to_ast(py),
                ast
            );
        });
    }
}
