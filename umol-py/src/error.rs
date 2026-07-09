//! Binding exceptions.

use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::PyErr;
use umol_ast::dsl::ParseError as AstParseError;

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
