//! Binding exceptions.

use pyo3::exceptions::PyException;
use pyo3::{create_exception, PyErr};
use umol_ast::ast::Contradiction as AstContradiction;
use umol_ast::dsl::ParseError as AstParseError;

create_exception!(
    umol,
    ContradictionError,
    PyException,
    "Raised when a umol operation reaches a contradiction."
);

create_exception!(
    umol,
    ParseError,
    PyException,
    "Raised when a umol DSL string fails to parse."
);

/// Map an `umol_ast` parse error onto the catchable `umol.ParseError`.
pub(crate) fn parse_error(error: AstParseError) -> PyErr {
    ParseError::new_err(error.to_string())
}

/// Map an `umol_ast` contradiction onto the catchable `umol.ContradictionError`.
pub(crate) fn contradiction_error(error: AstContradiction) -> PyErr {
    ContradictionError::new_err(error.to_string())
}

#[cfg(test)]
mod tests {
    use pyo3::prelude::*;
    use rstest::rstest;

    use super::*;

    #[rstest]
    fn test_contradiction_error() {
        Python::attach(|py| {
            let error = contradiction_error(AstContradiction);
            assert!(error.is_instance_of::<ContradictionError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "reached a contradiction"
            );
        });
    }
}
