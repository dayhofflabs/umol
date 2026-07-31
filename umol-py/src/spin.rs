//! Unpaired-electron AST values.

use pyo3::prelude::*;
use umol_ast::ast::UnpairedElectronsAst as AstUnpairedElectronsAst;

use crate::convert::{hash_rust, into_py_variant};
use crate::value::{ValueArg, ValueAst};

/// Unpaired-electron count and multiplicity as independent value fields.
#[pyclass]
pub struct UnpairedElectronsAst {
    #[pyo3(get)]
    count: Py<ValueAst>,
    #[pyo3(get)]
    multiplicity: Py<ValueAst>,
}

#[pymethods]
impl UnpairedElectronsAst {
    #[new]
    fn new(py: Python<'_>, count: ValueArg, multiplicity: ValueArg) -> PyResult<Self> {
        Ok(UnpairedElectronsAst {
            count: count.to_py(py)?,
            multiplicity: multiplicity.to_py(py)?,
        })
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "UnpairedElectronsAst({}, {})",
            self.count.bind(py).as_any().repr()?.extract::<String>()?,
            self.multiplicity
                .bind(py)
                .as_any()
                .repr()?
                .extract::<String>()?,
        ))
    }
}

impl UnpairedElectronsAst {
    pub(crate) fn from_rust(
        py: Python<'_>,
        ast: &AstUnpairedElectronsAst,
    ) -> PyResult<UnpairedElectronsAst> {
        Ok(UnpairedElectronsAst {
            count: into_py_variant(py, ValueAst::from_rust(py, &ast.count)?)?,
            multiplicity: into_py_variant(py, ValueAst::from_rust(py, &ast.multiplicity)?)?,
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> AstUnpairedElectronsAst {
        AstUnpairedElectronsAst {
            count: self.count.bind(py).borrow().to_rust(py),
            multiplicity: self.multiplicity.bind(py).borrow().to_rust(py),
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::{UnpairedElectronsAst as AstUnpairedElectronsAst, ValueAst as AstValueAst};

    use super::*;

    #[rstest]
    #[case::complete(AstUnpairedElectronsAst {
        count: AstValueAst::Lit(1),
        multiplicity: AstValueAst::Lit(2),
    })]
    #[case::physics_invalid(AstUnpairedElectronsAst {
        count: AstValueAst::Lit(2),
        multiplicity: AstValueAst::Lit(2),
    })]
    #[case::partial(AstUnpairedElectronsAst {
        count: AstValueAst::Undetermined,
        multiplicity: AstValueAst::Undetermined,
    })]
    fn test_unpaired_electrons_ast_roundtrip(#[case] ast: AstUnpairedElectronsAst) {
        Python::attach(|py| {
            assert_eq!(
                UnpairedElectronsAst::from_rust(py, &ast)
                    .unwrap()
                    .to_rust(py),
                ast
            );
        });
    }
}
